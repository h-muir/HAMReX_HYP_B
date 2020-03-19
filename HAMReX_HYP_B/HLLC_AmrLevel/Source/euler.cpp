#include "structdefs.H"
#include "funcdefs.H"
#include "classdefs.H"
#include <cmath>
#include <algorithm>

#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_TagBox.H>
#include <AmrLevelAdv.H>
#include <AMReX_AmrLevel.H>

#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBCellFlag.H>


using namespace amrex;
using namespace std;

/* ---------------------GLOBAL VARIABLES: ----------------------------*/
string EoS;		//system Equation of State --> set in initialiseStructs function.
Plasma19 AirPlasma("mixture19_cns.txt");
string scheme;
/*--------------------------------------------------------------------*/

/*----------------------------VECTOR DEFINITIONS--------------------------------------------------
 * 		
 * 		names:			   State vector (S):	     Flux vector (F):    	 Flux vector (G):
 * 
 * 	{	mass	}			{	 rho	}			{	  rho*u		}		{	  rho*v		}
 *  {	mom_x	} 			{	rho*u	}			{   rho*u^2 +p	}		{    rho*v*u	}
 *  {	mom_y	}			{	rho*v	}			{	 rho*u*v	}		{   rho*v^2 +p	}
 *  {	 en		}			{	  E		}			{	u*[E + p]	}		{	v*[E + p]	}
 *  {	rho		}			{	 rho	}			{  	    0		}		{  	    0		}
 *  {	 u		}			{	  u		}			{  	    0		}		{  	    0		}
 *  {	 v		}			{	  v		}			{	    0		}		{  	    0		}
 *  {	 p 		}			{	  p 	}			{	    0 		}		{  	    0		}
 * 
 *   			 			  where: E = rho*e + 1/2*rho*|(u,v)|^2
 * 						  				and: e = EoS(rho,p)
 -----------------------------------------------------------------------------------------------*/

//the following variable name sare assumed for the current system:
Vector<string> assumed_variables{"mass", "mom_x", "mom_y", "en", "rho", "u", "v", "p", "a", "T", "sigma", 
								"EoS19", "EB", "level_set", "nx", "ny", "B_r", "B_z", "B_theta", "B_abs", 
								 "J_r", "J_z", "J_theta", "J_abs", "phi"};

void check_assumed_variables(Vector<string> variable_names){
	//check that the variable names which are called in this system 
	// all exist in the inputted variables list
	
	for(int i = 0; i<assumed_variables.size(); ++i){
		if(std::find(variable_names.begin(), variable_names.end(), assumed_variables[i]) == variable_names.end())
		{
			cout << assumed_variables[i] << endl;
			amrex::Abort("above assumed variable was not found in the inputted list of variable names");
		}else{
			//assumed variable found
		}
	}
	cout << "all assumed variables have been inputted correctly from amr.derive_plot_vars" << endl;
}


void initialiseStructs(SettingsStruct& SimSettings, ParameterStruct& Parameters)
{
	/* -----------------------------------------------------
     * Global simulation parameters.
     * -----------------------------------------------------*/
     
    ParmParse pp;
	pp.get("nsteps_max", SimSettings.nsteps_max);
	pp.get("plot_int", SimSettings.plot_int);
    pp.get("ndmin", SimSettings.ndmin);
    pp.get("startT", SimSettings.startT);
    pp.get("chatter", SimSettings.chatter);
    pp.get("testcase", SimSettings.testcase);
    pp.get("solver", SimSettings.solver);
    pp.get("NCOMP", SimSettings.NCOMP);
    pp.get("NGROW", SimSettings.NGROW);
    pp.get("MUSCL", SimSettings.MUSCL);
    pp.get("scheme", SimSettings.scheme);
    pp.get("EoS", SimSettings.EoS);
    SimSettings.EBgeom = "none";
    pp.query("EBgeom", SimSettings.EBgeom);
    SimSettings.Bfield = "none";
    pp.query("Bfield", SimSettings.Bfield);

    
    /* -----------------------------------------------------
     * Test Case specific parameters.
     * -----------------------------------------------------*/   
    
    ParmParse ppGeom("geometry");
    ppGeom.getarr("prob_lo", Parameters.prob_lo); 
    ppGeom.getarr("prob_hi", Parameters.prob_hi);
    ppGeom.get("coord_sys", Parameters.coord_sys);
    ppGeom.getarr("is_periodic", Parameters.is_periodic);
     
    ParmParse ppTest(SimSettings.testcase);
    ppTest.get("finalT", SimSettings.finalT);
    ppTest.get("n_vars", Parameters.n_vars);
    ppTest.getarr("vars", Parameters.vars);
    ppTest.getarr("PL", Parameters.PL);
    ppTest.getarr("PR", Parameters.PR);
    ppTest.getarr("n_cells", Parameters.n_cells);
    ppTest.get("IC", Parameters.IC);
    ppTest.getarr("bc_conds", Parameters.bc_conds);
    if(Parameters.IC == "x0" || Parameters.IC == "y0" || Parameters.IC == "x0_perturbed"){ 
		ppTest.get("x0", Parameters.x0);
	}else if(Parameters.IC == "mach_theta"){
		ppTest.get("x0", Parameters.x0);
		ppTest.get("mach", Parameters.mach);
		ppTest.get("theta",Parameters.theta);	
	}
	
	ParmParse ppEB(SimSettings.EBgeom);
	if(SimSettings.EBgeom == "sphere"){
		ppEB.get("R", Parameters.EB_R);
		ppEB.getarr("centre", Parameters.EB_centre);
	}else if(SimSettings.EBgeom == "plane" || SimSettings.EBgeom == "planeX" || SimSettings.EBgeom == "triangle"){
		ppEB.getarr("point", Parameters.EB_point);
		ppEB.getarr("normal", Parameters.EB_normal);
		if(SimSettings.EBgeom == "planeX" || SimSettings.EBgeom == "triangle"){
			ppEB.getarr("point2", Parameters.EB_point2);
			ppEB.getarr("normal2", Parameters.EB_normal2);
			if(SimSettings.EBgeom == "triangle"){
				ppEB.getarr("point3", Parameters.EB_point3);
				ppEB.getarr("normal3", Parameters.EB_normal3);
			}	
		}
	}else if(SimSettings.EBgeom == "none"){
		//no additional settings
		cout << "EBgeom set to none. " << endl;
	}else{
		amrex::Abort("invalid EBgeom passed from inputs file");
	}
	
	if(SimSettings.Bfield != "none"){
		ParmParse ppB(SimSettings.Bfield);
		Parameters.Bfield = SimSettings.Bfield;
		if(SimSettings.Bfield == "dipole"){
			ppB.getarr("m", Parameters.dipole_m);
			ppB.getarr("centre", Parameters.dipole_centre);
		}else{
			amrex::Abort("invalid Bfield passed from inputs file");
		}
	}
    
    ParmParse ppAmr("amr");
    ppAmr.getarr("derive_plot_vars", Parameters.problem_variables);
    ppAmr.get("max_level", Parameters.max_level);
    ppAmr.get("max_grid_size", Parameters.max_grid_size);
    ppAmr.get("refinement_prop", Parameters.refinement_prop);
    ppAmr.get("refinement_condition", Parameters.refinement_condition);
    ppAmr.getarr("refinement_grad_fracs", Parameters.refinement_grad_fracs);
    
    //pp.getarr("problem_variables", Parameters.problem_variables);
      
    if( Parameters.problem_variables.size() != SimSettings.NCOMP){ 
		cout << "Parameters.problem_variables.size() : " << Parameters.problem_variables.size() << endl;
		cout << "SimSettings.NCOMP : " << SimSettings.NCOMP << endl;
		amrex::Abort("problem variables list must be of length NCOMP");
	}
	
	check_assumed_variables(Parameters.problem_variables);
    
    //set the global variable
    EoS = SimSettings.EoS;
    if(EoS != "ideal" && EoS != "plasma19"){
		amrex::Abort("invalid EoS input in test input file");
	}
    cout << "global variable EoS set to: " << EoS << endl;
    
    scheme = SimSettings.scheme;
    if(scheme != "SuperBee" && scheme != "MinBee" && scheme != "UltraBee" && scheme != "VanLeer"){
		amrex::Abort("invalid MUSCL limiter scheme in test input file");
	}
	
	//AirPlasma.testing_function();
		
}

void setBoundaryConditions(Vector<BCRec>& bc_vec, const ParameterStruct& p, const int& nc, AccessVariable& V){
	
	/*----------------------------------------------
					  hi_bc[1]
				 ________________
				|				 |
				|				 |
	  lo_bc[0]  |				 |  hi_bc[0]
				|				 |
				|________________|
				
					  lo_bc[1]
	
	
	order of bc_conds = [left, right, bottom, top]
	
	-----------------------------------------------*/
	
	int lo_bc[BL_SPACEDIM];
    int hi_bc[BL_SPACEDIM];
    
    //check BC's are all avild:
    for(int s = 0; s<4; ++s){
		if( !(p.bc_conds[s] == "transmissive") && !(p.bc_conds[s] == "reflective")){
			amrex::Abort("bc_conds from input file are invalid - must be 4 x transmissive/reflective");
		}
	}
    
    if(BL_SPACEDIM > 2){
		amrex::Abort("setBoundaryConditions not configured for BL_SPACEDIM > 2" );
	}
    
	for(int n=0; n<nc; ++n){
		for (int i = 0; i < BL_SPACEDIM; ++i) {
				
			if(i == 0){
				//lo_bc --> left:
				if(p.bc_conds[0] == "transmissive"){
					lo_bc[i] =  BCType::foextrap; //transmissive 
				}else if(p.bc_conds[0] == "reflective"){
					if(n == V["u"] || n == V["mom_x"]){
						lo_bc[i] = BCType::reflect_odd;  //reflect x velocities
					}else{
					lo_bc[i] = BCType::reflect_even; //transmit property
					}
				}
				//hi_bc --> right:
				if(p.bc_conds[1] == "transmissive"){
					hi_bc[i] = BCType::foextrap; //transmissive
				}else if(p.bc_conds[1] == "reflective"){
					if(n == V["u"] || n == V["mom_x"]){
						hi_bc[i] = BCType::reflect_odd;  //reflect x velocities
					}else{
						hi_bc[i] = BCType::reflect_even; //transmit property
					}
				}
			}
			if(i == 1){
				//lo_bc --> bottom:
				 if(p.bc_conds[2] == "transmissive"){
					lo_bc[i] =  BCType::foextrap; //transmissive
				}else if(p.bc_conds[2] == "reflective"){
					if(n == V["v"] || n == V["mom_y"]){
						lo_bc[i] = BCType::reflect_odd;  //reflect y velocities
					}else{
						lo_bc[i] = BCType::reflect_even; //transmit property
					}
				}
				//hi_bc --> top:
				if(p.bc_conds[3] == "transmissive"){
					hi_bc[i] =  BCType::foextrap; //transmissive
				}else if(p.bc_conds[3] == "reflective"){
					if(n == V["v"] || n == V["mom_y"]){
						hi_bc[i] = BCType::reflect_odd;  //reflect y velocities
					}else{
						hi_bc[i] = BCType::reflect_even; //transmit property
					}
				}
			}
			
		}//end i<dim 
		
		//store BC in bc_vec:
		BCRec bc(lo_bc, hi_bc);
		bc_vec[n] = bc;
		
	}//end n<ncomp		
	
}

//EoS function:
int EoS_properties(std::string eos, const double& rho, const double& p, const double& e, double& out_val, std::string prop){
	//EoS read from case: "ideal" or "EoS19"
	//properties: rho, p and e passed by reference
	//calculated property determined by strong prop: "e", "p", "a", "sigma_e"

	double gamma = 1.4;	
	int EoS19 = 0;
	
	if(rho < AirPlasma.rho_min){
		eos = "ideal";
	}else if(prop != "p" && p < AirPlasma.p_min){
		eos = "ideal";
	}
	
	if(eos == "ideal"){
		EoS19 = 0;
		if(prop == "e"){
			out_val = p/(rho*(gamma-1));
		}else if(prop == "p"){
			out_val = rho*(gamma - 1)*e;
		}else if(prop == "a"){
			out_val = pow(gamma*p/rho,0.5);
		}else if(prop == "sigma"){
			out_val = 200.0;
			//"can't produce conductivity vals from ideal EoS"
			//amrex::Abort("can't produce conductivity vals from ideal EoS");
		}else if(prop == "T"){
			out_val = p/(rho*287);
		}else if(prop == "gamma"){
			out_val = gamma;
		}else{
			amrex::Abort("property function not defined");
		}
	}else if(eos == "plasma19"){
		EoS19 = 1;
		if(prop == "e"){
			out_val = AirPlasma.getSpecificInternalEnergy(rho, p);
		}else if(prop == "p"){
			out_val = AirPlasma.getPressure(rho, e);
			if(out_val < AirPlasma.p_min){
				EoS19 = 0;
			}
		}else if(prop == "a"){
			out_val = AirPlasma.getSoundSpeed(rho, p);
		}else if(prop == "sigma"){
			out_val = AirPlasma.getConductivity(rho, p);
		}else if(prop == "T"){
			out_val = AirPlasma.getTemperature(rho, p);
		}else if(prop == "gamma"){
			out_val = AirPlasma.getAdiabaticIndex();
		}else{
			amrex::Abort("property function not defined");
		}
		//EoS19 = AirPlasma.EOS19;
	}else{
		amrex::Abort("EoS not defined");
	}
	
	return EoS19;
}

//conversion functions:

void S_conWtoU(Box const& box, Array4<Real> const& prop_arr, AccessVariable& V){
	
	const auto lo = lbound(box);
    const auto hi = ubound(box);
    
    double e;
    
    for(int k = lo.z; k <= hi.z; ++k) {
		for(int j = lo.y; j <= hi.y; ++j) {
			for(int i = lo.x; i <= hi.x; ++i) {
				prop_arr(i,j,k,V["mass"]) = prop_arr(i,j,k,V["rho"]);
				prop_arr(i,j,k,V["mom_x"]) = prop_arr(i,j,k,V["rho"])*prop_arr(i,j,k,V["u"]);
				prop_arr(i,j,k,V["mom_y"]) = prop_arr(i,j,k,V["rho"])*prop_arr(i,j,k,V["v"]);
				prop_arr(i,j,k,V["EoS19"]) = EoS_properties(EoS, prop_arr(i,j,k,V["rho"]), prop_arr(i,j,k,V["p"]), 0.0, e, "e"); //compute e
				prop_arr(i,j,k,V["en"]) = prop_arr(i,j,k,V["rho"])*e + 0.5*prop_arr(i,j,k,V["rho"])*\
											(pow(prop_arr(i,j,k,V["u"]),2) + pow(prop_arr(i,j,k,V["v"]),2));
				prop_arr(i,j,k,V["EoS19"]) = EoS_properties(EoS, prop_arr(i,j,k,V["rho"]), prop_arr(i,j,k,V["p"]), e, prop_arr(i,j,k,V["a"]), "a");
				prop_arr(i,j,k,V["EoS19"]) = EoS_properties(EoS, prop_arr(i,j,k,V["rho"]), prop_arr(i,j,k,V["p"]), e, prop_arr(i,j,k,V["T"]), "T");
				prop_arr(i,j,k,V["EoS19"]) = EoS_properties(EoS, prop_arr(i,j,k,V["rho"]), prop_arr(i,j,k,V["p"]), e, prop_arr(i,j,k,V["sigma"]), "sigma");
			}
		}
	}
	
}

void S_conUtoW(Box const& box, Array4<Real> const& prop_arr, AccessVariable& V){
	
	const auto lo = lbound(box);
    const auto hi = ubound(box);
    
    double e, p;
    
    for(int k = lo.z; k <= hi.z; ++k) {
		for(int j = lo.y; j <= hi.y; ++j) {
			for(int i = lo.x; i <= hi.x; ++i) {	
				prop_arr(i,j,k,V["rho"]) = prop_arr(i,j,k,V["mass"]);
				prop_arr(i,j,k,V["u"])	 = prop_arr(i,j,k,V["mom_x"])/prop_arr(i,j,k,V["rho"]);
				prop_arr(i,j,k,V["v"])	 = prop_arr(i,j,k,V["mom_y"])/prop_arr(i,j,k,V["rho"]);
				e = 1.0/prop_arr(i,j,k,V["rho"])*(prop_arr(i,j,k,V["en"]) - \
						0.5*prop_arr(i,j,k,V["rho"])*(pow(prop_arr(i,j,k,V["u"]),2) + pow(prop_arr(i,j,k,V["v"]),2))); 
				prop_arr(i,j,k,V["EoS19"]) = EoS_properties(EoS, prop_arr(i,j,k,V["rho"]), 0.0, e, p, "p"); //compute p
				prop_arr(i,j,k,V["p"]) = p;
				prop_arr(i,j,k,V["EoS19"]) = EoS_properties(EoS, prop_arr(i,j,k,V["rho"]), p, e, prop_arr(i,j,k,V["a"]), "a");
				prop_arr(i,j,k,V["EoS19"]) = EoS_properties(EoS, prop_arr(i,j,k,V["rho"]), p, e, prop_arr(i,j,k,V["T"]), "T");
				prop_arr(i,j,k,V["EoS19"]) = EoS_properties(EoS, prop_arr(i,j,k,V["rho"]), p, e, prop_arr(i,j,k,V["sigma"]), "sigma");
			}
		}
	}
   
}	

void F_conWtoF(Box const& box, Array4<Real> const& prop_arr, AccessVariable& V){
	
	const auto lo = lbound(box);
    const auto hi = ubound(box);
    
    double e, E;
    
    for(int k = lo.z; k <= hi.z; ++k) {
		for(int j = lo.y; j <= hi.y; ++j) {
			for(int i = lo.x; i <= hi.x; ++i) {	
				prop_arr(i,j,k,V["mass"]) = prop_arr(i,j,k,V["rho"])*prop_arr(i,j,k,V["u"]);
				prop_arr(i,j,k,V["mom_x"]) = prop_arr(i,j,k,V["rho"])*pow(prop_arr(i,j,k,V["u"]),2)+prop_arr(i,j,k,V["p"]);
				prop_arr(i,j,k,V["mom_y"]) = prop_arr(i,j,k,V["rho"])*prop_arr(i,j,k,V["u"])*prop_arr(i,j,k,V["v"]);
				prop_arr(i,j,k,V["EoS19"]) = EoS_properties(EoS, prop_arr(i,j,k,V["rho"]), prop_arr(i,j,k,V["p"]), 0.0, e, "e");
				E = prop_arr(i,j,k,V["rho"])*e + 0.5*prop_arr(i,j,k,V["rho"])*\
				    (pow(prop_arr(i,j,k,V["u"]),2) + pow(prop_arr(i,j,k,V["v"]),2));
				prop_arr(i,j,k,V["en"]) = prop_arr(i,j,k,V["u"])*(E+prop_arr(i,j,k,V["p"]));
			}
		}
	}
   
}			
	

//Conversion functions between classes:

ConsF conWtoF(Prim const& W){
    double mass = W.rho*W.u_vec[0];
    Vector<Real> mom_vec(3);
    mom_vec[0] = W.rho*pow(W.u_vec[0],2)+W.p;
    mom_vec[1] = W.rho*W.u_vec[0]*W.u_vec[1];
    mom_vec[2] = 0;
    double e;
    int EoS19 = EoS_properties(EoS, W.rho, W.p, 0.0, e, "e");
    double E =  W.rho*e + 0.5*W.rho*pow(mag(W.u_vec),2);
    double en = W.u_vec[0]*(E + W.p );
    ConsF F(mass, mom_vec, en);
    return F;
}

Prim conUtoW(ConsU const& U){
    double rho = U.rho;
    Vector<Real> u_vec(3);
    u_vec[0] = U.rhou_vec[0]/rho;
    u_vec[1] = U.rhou_vec[1]/rho;
    u_vec[2] = 0;
    double e = 1/rho*(U.E - 0.5*rho*pow(mag(u_vec),2));
    double p;
    int EoS19 = EoS_properties(EoS, rho, 0.0, e, p, "p");
    Prim W(rho, u_vec, p);
    W.r = U.r;
    return W;
}

ConsU conWtoU(Prim const& W){
    double rho = W.rho;
    Vector<Real> rhou_vec(3);
    rhou_vec[0] = W.rho*W.u_vec[0];
    rhou_vec[1] = W.rho*W.u_vec[1];
    rhou_vec[2] = W.rho*W.u_vec[2];
    double e;
    int EoS19 = EoS_properties(EoS, W.rho, W.p, 0.0, e, "e");
    double E = W.rho*e + 0.5*W.rho*pow(mag(W.u_vec),2);
    ConsU U(rho, rhou_vec, E);
    U.r = W.r;
    return U;
}


ConsF conWtoG(Prim const& W){
    double mass = W.rho*W.u_vec[1];
    Vector<Real> mom_vec(3);
    mom_vec[0] = W.rho*W.u_vec[1]*W.u_vec[0];
    mom_vec[1] = W.rho*W.u_vec[1]*W.u_vec[1]+W.p;
    mom_vec[2] = 0;
    double e;
    int EoS19 = EoS_properties(EoS, W.rho, W.p, 0.0, e, "e");
    double E =  W.rho*e + 0.5*W.rho*pow(mag(W.u_vec),2);
    double en = W.u_vec[1]*(E + W.p);
    ConsF G(mass, mom_vec, en);
    return G;
}

ConsF conWtoSE(Prim const& W){
    double mass = ((W.r > 0)? -1/W.r * W.rho*W.u_vec[0] : 0);
    Vector<Real> mom_vec(3);
    mom_vec[0] = ((W.r > 0)? -1/W.r * W.rho*pow(W.u_vec[0],2) : 0);
    mom_vec[1] = 0;
    double e;
    int EoS19 = EoS_properties(EoS, W.rho, W.p, 0.0, e, "e");
    double E =  W.rho*e + 0.5*W.rho*pow(mag(W.u_vec),2);
    double en = ((W.r > 0)? -1/W.r * W.u_vec[0]*(E+W.p) : 0);
    ConsF SE(mass, mom_vec, en);
    return SE;
}

//rankine-hugoniot calculation:
Prim rankine_hugoniot(Prim const& W, double const& MS, Vector<Real> const& normal_vec){
	//rankine hugoniot relation derived for ideal gas with gamma = 1.4
	//assume intitial conditions are in ideal regime when this function 
	//is called
	//also assumes a normal shock along normalised_vec:
	double gamma = 1.4;
	double eR;
	int EoS19 = EoS_properties(EoS, W.rho, W.p, 0.0, eR, "e");
	if(EoS19){ 
		cout << "WARNING: rankine hugoniot calc assumes ideal gas regime \
					but plasma validity range redetected" << endl;
		EoS_properties("ideal", W.rho, W.p, 0.0, eR, "e");
	}
	double aR;
	EoS19 = EoS_properties("ideal", W.rho, W.p, eR, aR, "a");
	Vector<Real> MR_vec{W.u_vec[0]/aR, W.u_vec[1]/aR, 0};
	double MR 	= mag(MR_vec);
	double rho 	= (((gamma+1)*pow(MR-MS,2))/((gamma-1)*pow(MR-MS,2)+2))*W.rho;
	double p   	= ((2*gamma*pow(MR-MS,2)-(gamma-1))/(gamma+1))*W.p;
	double S3 	= MS*aR;
	Vector<Real> S3_vec = S3*normal_vec;
	Vector<Real> u_vec = (1-W.rho/rho)*S3_vec + (W.rho/rho)*W.u_vec;
	
	cout << "computed properties from analytic rankine hugoniot relations: \n";
	cout << "gamma = " << gamma << endl;
	cout << "aR = " << aR << endl;
	cout << "MR = " << MR << endl;
	cout << "rho = " << rho << endl;
	cout << "p = " << p << endl;
	cout << "S3_vec : \n";
	printVector(S3_vec);
	cout << "u_vec : \n";
	printVector(u_vec);

	Prim WL(rho, u_vec, p);
    return WL;
	
}

//problem initialisation:
void initial_conditions(int const& nc, Box const& box, Array4<Real> const& prop_arr, const Real* dx, ParameterStruct const& p, AccessVariable& V){

    const auto lo = lbound(box);
    const auto hi = ubound(box);
	
	//amrex::Print() << "lo : " << lo << ", hi : " << hi << "\n";

	//will need changing to *dx when cell sizes change in AMR
    int int_x0 = (p.x0/p.prob_hi[0])*p.n_cells[0];
    int int_y0 = (p.x0/p.prob_hi[1])*p.n_cells[1];
	
	//AccessVariable V(p.problem_variables);
	
	//cout << "PL and PR size() " << p.PL.size() << " : " << p.PR.size() << endl;
	
	Real R = 0.2;
	Real h;
	
	//set all intitial data to zero - this seems to make visit happy:
	for(int k = lo.z; k <= hi.z; ++k) {
		for(int j = lo.y; j <= hi.y; ++j) {
			for(int i = lo.x; i <= hi.x; ++i) {
				for(int n=0; n<nc; ++n){
					prop_arr(i,j,k,n) = 0.0;
				}
			}
		}
	}
	
	
	if(p.IC == "source"){
		for(int n = 0; n<p.n_vars; ++n){
		for(int k = lo.z; k <= hi.z; ++k) {
			for(int j = lo.y; j <= hi.y; ++j) {
				for(int i = lo.x; i <= hi.x; ++i) {
					
					h = sqrt(pow((i+0.5)*dx[0]-0.5,2) + pow((j+0.5)*dx[1]-0.5,2));
					if(h <= R){
						prop_arr(i,j,k,V[p.vars[n]]) = p.PL[n];
					}else{
						prop_arr(i,j,k,V[p.vars[n]]) = p.PR[n];
					}
					
				}
			}
		}
	}		
	}else if(p.IC == "x0"){
		for(int n = 0; n<p.n_vars; ++n){
		for(int k = lo.z; k <= hi.z; ++k) {
			for(int j = lo.y; j <= hi.y; ++j) {
				for(int i = lo.x; i <= hi.x; ++i) {
					if(i*dx[0] <= p.x0){
						prop_arr(i,j,k,V[p.vars[n]]) = p.PL[n]; 
					}else{
						prop_arr(i,j,k,V[p.vars[n]]) = p.PR[n];
					}
					
				}
			}
		}
		}
	}else if(p.IC == "x0_perturbed"){
		for(int n = 0; n<p.n_vars; ++n){
		for(int k = lo.z; k <= hi.z; ++k) {
			for(int j = lo.y; j <= hi.y; ++j) {
				for(int i = lo.x; i <= hi.x; ++i) {
					if(i*dx[0] <= p.x0){
						prop_arr(i,j,k,V[p.vars[n]]) = p.PL[n]; 
					}else{
						prop_arr(i,j,k,V[p.vars[n]]) = p.PR[n];
					}
					//perturbation along centreline:
					if(i*dx[0] <=  p.x0){
						prop_arr(i,j,k,V["rho"]) += pow(-1,j)*0.05;
					}
					
				}
			}
		}
		}
	}else if(p.IC == "y0"){
		for(int n = 0; n<p.n_vars; ++n){
		for(int k = lo.z; k <= hi.z; ++k) {
			for(int j = lo.y; j <= hi.y; ++j) {
				for(int i = lo.x; i <= hi.x; ++i) {
					if(j*dx[1] <= p.x0){
						prop_arr(i,j,k,V[p.vars[n]]) = p.PL[n]; 
					}else{
						prop_arr(i,j,k,V[p.vars[n]]) = p.PR[n];
					}
					
				}
			}
		}
		}
	}else if(p.IC == "mach_theta"){
		for(int n = 0; n<p.n_vars; ++n){
		for(int k = lo.z; k <= hi.z; ++k) {
			for(int j = lo.y; j <= hi.y; ++j) {
				for(int i = lo.x; i <= hi.x; ++i) {
					if(i*dx[0] <= p.x0 + (j*dx[1])*tan(p.theta)){
						prop_arr(i,j,k,V[p.vars[n]]) = p.PL[n]; 
					}else{
						prop_arr(i,j,k,V[p.vars[n]]) = p.PR[n];
					}	
				}
			}
		}
		}
	}else if(p.IC == "log_continuous"){
		for(int k = lo.z; k <= hi.z; ++k) {
			for(int j = lo.y; j <= hi.y; ++j) {
				for(int i = lo.x; i <= hi.x; ++i) {
					prop_arr(i,j,k,V["rho"]) = pow(10, log10(p.PL[0]) + (i*dx[0])*(log10(p.PR[0]) - log10(p.PL[0]))); 
					prop_arr(i,j,k,V["u"]) = p.PL[1];
					prop_arr(i,j,k,V["v"]) = p.PL[2];
					prop_arr(i,j,k,V["p"]) = pow(10, log10(p.PL[3]) + (j*dx[1])*(log10(p.PR[3]) - log10(p.PL[3])));
				}
			}
		}
	}else{
		amrex::Abort("invalid intitial condition - check inputs file for adv.IC");
	}

	S_conWtoU(box, prop_arr, V);
	
}

void initialise_Bfield(Box const& box, Array4<Real> const& prop_arr, const Real* dx, ParameterStruct const& p, AccessVariable& V){

    const auto lo = lbound(box);
    const auto hi = ubound(box);
	
	
	if(p.Bfield == "dipole"){
		Vector<Real> r_vec(3, 0.0);
		Vector<Real> B_vec(3, 0.0);
		Real r_mag;
		for(int k = lo.z; k <= hi.z; ++k) {
			for(int j = lo.y; j <= hi.y; ++j) {
				r_vec[1] = j*dx[1] - p.dipole_centre[1];
				for(int i = lo.x; i <= hi.x; ++i) {
					if(prop_arr(i,j,k,V["EB"]) <= 0){
						prop_arr(i,j,k,V["B_r"]) = 0.0; 
						prop_arr(i,j,k,V["B_z"]) = 0.0;
						prop_arr(i,j,k,V["B_theta"]) = 0.0;
						prop_arr(i,j,k,V["B_abs"]) = 0.0;
					}else{
						r_vec[0] = i*dx[0] - p.dipole_centre[0];
						r_mag = mag(r_vec);
						if(r_mag == 0){
							amrex::Abort("dipole centre has been initialised outside of geometry - not allowed (singularity)");
						}
						B_vec = (1/(4*M_PI*pow(r_mag,5)))*(3*dot(r_vec, p.dipole_m)*r_vec - pow(r_mag, 2)*p.dipole_m);
						prop_arr(i,j,k,V["B_r"]) = B_vec[0]; 
						prop_arr(i,j,k,V["B_z"]) = B_vec[1];
						prop_arr(i,j,k,V["B_theta"]) = B_vec[2];
						prop_arr(i,j,k,V["B_abs"]) = mag(B_vec);
					}
				}
			}
		}//end i,j,k
	}
	
}

void initialise_problem(MultiFab& S_new, const FabArray<EBCellFlagFab>& flags, Geometry const& geom, SettingsStruct const& SimSettings, 
								ParameterStruct const& Parameters, AccessVariable& V){
	
	const Real* dx  = geom.CellSize();
	
	for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
	{
		const Box& box     = mfi.validbox();
		FArrayBox& fab 	   = S_new[mfi];
		const int* lo      = box.loVect();
		const int* hi      = box.hiVect();
		
		Array4<Real> const& S_new_arr = fab.array();	
		
		initial_conditions(S_new.nComp(), box, S_new_arr, dx, Parameters, V);
		
		if(SimSettings.EBgeom == "none"){
			// all cells set to fluid:
			S_new[mfi].setVal(1.0, box, V["EB"], 1);
		}else{
			
			FabType EB_fab_type		 	= 	flags[mfi].getType(box); //EB
			auto const& flags_bx		= 	flags[mfi];
			Array4<const EBCellFlag> flags_arr	= flags_bx.array();
			
			if(EB_fab_type == FabType::covered) {
				//cout << "box is covered \n";
				//rigid body
				S_new[mfi].setVal(-1.0, box, V["EB"], 1);
			}else if(EB_fab_type == FabType::regular){
				//cout << "box is regular \n";
				//fluid
				S_new[mfi].setVal(1.0, box, V["EB"], 1);
			}else if(EB_fab_type == FabType::singlevalued){
				//cout << "box has cut cells: \n";
				find_embedded_boundary(box, S_new_arr, flags_arr, S_new.nComp(), V);
			}
		}
		
		if(SimSettings.Bfield != "none"){
			initialise_Bfield(box, S_new_arr, dx, Parameters, V);
		} 

	}//end MFIter
			
}

void rigid_body_step1(MultiFab& S_new, MultiFab& S_old, const FabArray<EBCellFlagFab>& flags, const Real* dx, AccessVariable& V)
{
	//Rigid body step 1:
	for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
        const Box& box     = mfi.validbox();
        FArrayBox& fab 	   = S_old[mfi];
        const int* lo      = box.loVect();
        const int* hi      = box.hiVect();
        
        Array4<Real> const& S_old_arr = fab.array();
        
        FabType EB_fab_type		 	= 	flags[mfi].getType(box); //EB
	    auto const& flags_bx		= 	flags[mfi];
		Array4<const EBCellFlag> flags_arr	= flags_bx.array();
		
		if(EB_fab_type == FabType::singlevalued){
			//cout << "box has cut cells: \n";
			initialise_boundary_cells(box, S_old_arr, flags_arr, S_new.nComp(), V, dx);
		}else if(EB_fab_type == FabType::covered){
			//cout << "box is regular \n";
			//S_new[mfi].setVal(-1e5, box, V["p"], 1);
			initialise_boundary_cells(box, S_old_arr, flags_arr, S_new.nComp(), V, dx);
		}
		  
    }
		
}

void rigid_body_step2(MultiFab& S_new, MultiFab& S_old, const FabArray<EBCellFlagFab>& flags, const Real* dx, AccessVariable& V, 
							double &global_p_min)
{ 
        
	for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
	{
		const Box& box     = mfi.validbox();
		FArrayBox& fab 	   = S_old[mfi];
		const int* lo      = box.loVect();
		const int* hi      = box.hiVect();
		
		Array4<Real> const& S_old_arr = fab.array();
		
		FabType EB_fab_type		 	= 	flags[mfi].getType(box); //EB
		auto const& flags_bx		= 	flags[mfi];
		Array4<const EBCellFlag> flags_arr	= flags_bx.array();
		
		if(EB_fab_type == FabType::singlevalued){
			//cout << "box has cut cells: \n";
			sweep(box, S_old_arr, flags_arr, S_new.nComp(), V, dx, global_p_min);
		}else if(EB_fab_type == FabType::covered){
			//only operate on boundary cells
			sweep(box, S_old_arr, flags_arr, S_new.nComp(), V, dx, global_p_min);
		}
		  
	}
	
}


void special_boundary(Box const& dom, MultiFab& S_old, SettingsStruct const& sim, const Real* dx, AccessVariable& V){
    
    //S_new box passed here. box is whole physical domain as a box
    const auto prob_lo = lbound(dom); 
    const auto prob_hi = ubound(dom);
    
    //cout << "prob_lo.x = " << prob_lo.x << endl;
    //cout << "prob_hi.x = " << prob_hi.x << endl;
    
    ParmParse ppTest(sim.testcase);
    string specialBC = "none";
    ppTest.query("specialBC", specialBC);
    if(specialBC == "partial_right" || specialBC == "partial_lower"){
		Real dist = 0.2;
		ppTest.query("transmissive_dist", dist);
		
		for (MFIter mfi(S_old, true); mfi.isValid(); ++mfi){
		//H --------------------------------------------------------------------
	    const Box& bx = mfi.tilebox();

	    Array4<Real> const& S_old_arr = S_old[mfi].array();	//H
	    //--------------------------------------------------------------------H/		
		
		const auto lo = lbound(bx); //lo bound of real domain (excluding ghost cells)
		const auto hi = ubound(bx);
		
		if(specialBC == "partial_lower" && lo.y == prob_lo.y && lo.x*dx[0] <= dist){	
		//overwrite bottom reflective portion with transmissive
			for(int k = lo.z; k <= hi.z; ++k) {
				for(int j = prob_lo.y-S_old.nGrow(); j < prob_lo.y; ++j) {
					for(int i = lo.x-S_old.nGrow(); i*dx[0] <= dist; ++i) {
						S_old_arr(i,j,k,V["mass"]) 	= S_old_arr(i,prob_lo.y,k,V["mass"]);
						S_old_arr(i,j,k,V["mom"]) 	= S_old_arr(i,prob_lo.y,k,V["mom"]);
						S_old_arr(i,j,k,V["en"]) 	= S_old_arr(i,prob_lo.y,k,V["en"]);
						S_old_arr(i,j,k,V["rho"]) 	= S_old_arr(i,prob_lo.y,k,V["rho"]);
						S_old_arr(i,j,k,V["u"]) 	= S_old_arr(i,prob_lo.y,k,V["u"]);
						S_old_arr(i,j,k,V["v"]) 	= S_old_arr(i,prob_lo.y,k,V["v"]);
						S_old_arr(i,j,k,V["p"]) 	= S_old_arr(i,prob_lo.y,k,V["p"]);
					}
				}
			}//end i,j,k loop
		
		}//end partial_lower altered portion update
		
		if(specialBC == "partial_right" && hi.x == prob_hi.x && lo.y*dx[1] <= dist){
		//overwrite right reflective portion with transmissive
			for(int k = lo.z; k <= hi.z; ++k) {
				for(int j = lo.y; j*dx[1] <= dist; ++j) {
					for(int i = prob_hi.x; i < prob_hi.x + S_old.nGrow(); ++i) {
						S_old_arr(i,j,k,V["mass"]) 	= S_old_arr(prob_hi.x,j,k,V["mass"]);
						S_old_arr(i,j,k,V["mom"]) 	= S_old_arr(prob_hi.x,j,k,V["mom"]);
						S_old_arr(i,j,k,V["en"]) 	= S_old_arr(prob_hi.x,j,k,V["en"]);
						S_old_arr(i,j,k,V["rho"]) 	= S_old_arr(prob_hi.x,j,k,V["rho"]);
						S_old_arr(i,j,k,V["u"]) 	= S_old_arr(prob_hi.x,j,k,V["u"]);
						S_old_arr(i,j,k,V["v"]) 	= S_old_arr(prob_hi.x,j,k,V["v"]);
						S_old_arr(i,j,k,V["p"]) 	= S_old_arr(prob_hi.x,j,k,V["p"]);
					}
				}
			}//end i,j,k loop
		
		}//end partial_right altered portion update
		
		}//end MFIter
	}else if(specialBC == "cylindrical_Z"){ 
		
		Real delta_phi = 0; //property gradient
		
		for (MFIter mfi(S_old, true); mfi.isValid(); ++mfi){

			const Box& bx = mfi.tilebox();
			Array4<Real> const& S_old_arr = S_old[mfi].array();	//H	
			
			const auto lo = lbound(bx); //lo bound of real domain (excluding ghost cells)
			const auto hi = ubound(bx);
			double e;
			int i = 0;
			
			if(lo.x == prob_lo.x){
			//overwrite right reflective portion with transmissive
				for(int k = lo.z; k <= hi.z; ++k) {
					for(int j = lo.y; j <= hi.y; ++j) {
						for(int ng = 1; ng < S_old.nGrow()+1; ++ng) {
							for(int v = 0; v < sim.NCOMP; ++v){
								//delta_phi = S_old_arr(prob_lo.x+ng-1,j,k,v) - S_old_arr(prob_lo.x+ng,j,k,v);
								//S_old_arr(prob_lo.x-ng,j,k,v) = S_old_arr(prob_lo.x-(ng-1),j,k,v) + delta_phi;
								//S_old_arr(prob_lo.x-ng,j,k,v) = S_old_arr(prob_lo.x,j,k,v);
							}
							i = prob_lo.x-ng;
							S_old_arr(i,j,k,V["u"]) = 0.0;
							S_old_arr(i,j,k,V["mom_x"]) = 0.0;
							S_old_arr(i,j,k,V["EoS19"]) = EoS_properties(EoS, S_old_arr(i,j,k,V["rho"]), S_old_arr(i,j,k,V["p"]), 0.0, e, "e"); //compute e
							S_old_arr(i,j,k,V["en"]) = S_old_arr(i,j,k,V["rho"])*e + 0.5*S_old_arr(i,j,k,V["rho"])*\
											(pow(S_old_arr(i,j,k,V["u"]),2) + pow(S_old_arr(i,j,k,V["v"]),2)); 
							
						}
					}
				}//end i,j,k loop
			
			}//end partial_right altered portion update
			
		}//end MFIter
	}else{
		//no additional boundary treatements necessary
	}
}
 

double r_calc(double n_val, double l_val){
    //alternative conventions for dividing by zero RHS slope:
    // -3 encoded to represent -ve infinity
    // +3 encoded to represent +ve infinity 
    double r = l_val;
    if(n_val == 0){
        if(l_val== 0){
            r = 0;           //both left and right slopes are zero
        }else if(l_val < 0){
            r = -3;          //negative slope on left, zero slope on right
        }else{
            r = 3;           //positive slope on left, zero slope on right
        }
    }else{ r = r/n_val;}     //permissible division by delta_i_n

    return r;
}

double delta_bar_calc(string& scheme, double& r, double& delta_i, double& sigma_R, double& sigma_L){
    double delta_bar;
    if(scheme == "SuperBee"){
        if(r < 0){
            delta_bar = 0;
        }else if(0 <= r && r < 0.5){
            delta_bar = 2*r*delta_i;
        }else if(0.5 <= r && r < 1){
            delta_bar = delta_i;
        }else if(1 <= r && r < 2){
            delta_bar = ((r < sigma_R)? 
                    r*delta_i : sigma_R*delta_i);
        }else{
            delta_bar = ((sigma_R < 2)? 
                    sigma_R*delta_i : 2*delta_i);
        }
    }else if(scheme == "MinBee"){

        if(r <= 0){
            delta_bar = 0;
        }else if(0 <= r && r <= 1.0){
            delta_bar = r*delta_i;
        }else{
            delta_bar = ((1 < sigma_R)? delta_i : sigma_R*delta_i);
        }
    }else if(scheme == "VanLeer"){
        if(r < 0){
            delta_bar = 0;
        }else{
            delta_bar = ((2*r/(1+r) < sigma_R)? 
                    (2*r/(1+r))*delta_i : sigma_R*delta_i);
        }
    }else if(scheme == "UltraBee"){
        if(r <= 0){
            delta_bar = 0;
        }else{
            delta_bar = ((sigma_L < sigma_R)? 
                    sigma_L*delta_i : sigma_R*delta_i);
        }
    }else{
		cout << "no scheme name given to delta_bar_calc" << endl;
        return delta_i;
    }

    return delta_bar;
}

Prim delta(double& omega, Prim& Wl, Prim& Wi, Prim& Wr, string& scheme){
    //omega E [-1, 1]
    //Wl = W_{i-1}, Wi = Wi, Wr = W_{i+1}
    Prim delta_i_l(Wi.rho - Wl.rho, Wi.u_vec - Wi.u_vec, 
                    Wi.p - Wl.p);
    Prim delta_i_n(Wr.rho - Wi.rho, Wr.u_vec - Wi.u_vec, 
                    Wr.p - Wi.p);
                    
    //r calc for each W property
    double r_rho = r_calc(delta_i_n.rho, delta_i_l.rho);
    vector<double> r_u_vec = delta_i_l.u_vec;
    for(int i = 0; i<3; ++i){
		r_u_vec[i] = r_calc(delta_i_n.u_vec[i], delta_i_l.u_vec[i]);
	}
    double r_p = r_calc(delta_i_n.p, delta_i_l.p);
	
	//delta calc:
    Prim delta_i(0.5*(1+omega)*delta_i_l.rho+0.5*(1-omega)*delta_i_n.rho,
                  0.5*(1+omega)*delta_i_l.u_vec+0.5*(1-omega)*delta_i_n.u_vec,
                  0.5*(1+omega)*delta_i_l.p+0.5*(1-omega)*delta_i_n.p);
	
	//sigma calcs:
	Vector<Real> sigma_u_vec(3);
	
	for(int i=0; i<3; ++i){
		sigma_u_vec[i] = 2.0/(1-omega+(1+omega)*r_u_vec[i]);
	}
	
    Prim sigma_R(2/(1-omega+(1+omega)*r_rho), sigma_u_vec, 
                    2/(1-omega+(1+omega)*r_p));
                    
    for(int i=0; i<3; ++i){
		sigma_u_vec[i] = 2.0*r_u_vec[i]/(1-omega+(1+omega)*r_u_vec[i]);
	}

    Prim sigma_L(2*r_rho/(1-omega+(1+omega)*r_rho), sigma_u_vec, 
                 2*r_p/(1-omega+(1+omega)*r_p));

    //adjustments for infinity encoded cases:
    if(r_rho == -3 || r_rho == 3){
        sigma_R.rho = 0;
    }
    for(int i=0; i<3; ++i){
		if(r_u_vec[i] == -3 || r_u_vec[i] == 3){
			sigma_R.u_vec[i] = 0;
		}
	}
    if(r_p == -3 || r_p == 3){
        sigma_R.p = 0;
    }
    
    //delta bar calc for each W property
    Vector<Real> v0(3); //zero vec for initialisation 
    Prim delta_bar(0,v0,0);
    
    delta_bar.rho = delta_bar_calc(scheme, r_rho, delta_i.rho, sigma_R.rho, sigma_L.rho);
    for(int i=0; i<3; ++i){
		delta_bar.u_vec[i] = delta_bar_calc(scheme, r_u_vec[i], delta_i.u_vec[i], sigma_R.u_vec[i], sigma_L.u_vec[i]);
	}
	delta_bar.p = delta_bar_calc(scheme, r_p, delta_i.p, sigma_R.p, sigma_L.p);
    
    return delta_bar;

}

double compute_alpha(const ConsU& U){
	
	Prim W = conUtoW(U);
	double V = mag(W.u_vec);
	double a;
	int EoS19 = EoS_properties(EoS, W.rho, W.p, 0.0, a, "a");
	double mach = V/a;
	double alpha = (5.0/12.0)*mach;
	
	//cout << "alpha = " << alpha << endl;
	
	alpha = 4.0;
	
	return alpha;
	
}

double compute_eps_x(Array4<Real> const& S_old, const int& i, const int& j, const int& k, AccessVariable& V){
	
	/* 

	 * limiter for interface between 3-4:
	 
		|     |     |
		|  2  |  1  |
		|_____|_____|
		|     |     |
		|  4  |  3  |
		|_____|_____|
		|     |     |
		|  6  |  5  |
		|     |     |
	 	
	*/
	
	double v1,v2,v3,v4,v5,v6;
	
	v3 = S_old(i,j,k,V["v"]);
	v4 = S_old(i-1,j,k,V["v"]);
	v1 = S_old(i,j+1,k,V["v"]);
	v2 = S_old(i-1,j+1,k,V["v"]);
	v5 = S_old(i,j-1,k,V["v"]);
	v6 = S_old(i-1,j-1,k,V["v"]);
	
	double n1, n2, n3, n4;
		
	n1 = fabs(v1 - v3);
	n2 = fabs(v2 - v4);
	n3 = fabs(v4 - v6);
	n4 = fabs(v3 - v5);

	double n_max = max(n1,n2,n3,n4);

	return n_max;	
}

double compute_eps_y(Array4<Real> const& S_old, const int& i, const int& j, const int& k, AccessVariable& V){
	
	/* 
	   G-stencil:
	 * limiter for interface between 3-4:
		 _____ _____ _____
		|     |     |     |
		|  1  |  3  |  5  |
		|_____|_____|_____|
		|     |     |     |
		|  2  |  4  |  6  |
		|_____|_____|_____|
	 	
	*/
	double u1, u2, u3, u4, u5, u6;
	
	u3 = S_old(i,j,k,V["u"]);
	u4 = S_old(i,j-1,k,V["u"]);
	u1 = S_old(i-1,j,k,V["u"]);
	u2 = S_old(i-1,j-1,k,V["u"]);
	u5 = S_old(i+1,j,k,V["u"]);
	u6 = S_old(i+1,j-1,k,V["u"]);
	
	double n1, n2, n3, n4;	

	n1 = fabs(u1 - u3);
	n2 = fabs(u2 - u4);
	n3 = fabs(u4 - u6);
	n4 = fabs(u3 - u5);
	
	double n_max = max(n1,n2,n3,n4);

	return n_max;	
	
}
	

void calc_fluxes_HLLC_HS_x(Box const& box, Array4<Real> const& S_old, Array4<Real> const& F, 
					Real const dt, Real const dx, int const nc, AccessVariable& V, bool MUSCL)
{
	const auto lo = lbound(box);
    const auto hi = ubound(box);
    
    Vector<Real> v0(3);  

	//MUSCL-Hancock parameters (initialised):    
	Prim Wll(0,v0,0), Wl(0,v0,0); 
	Prim WbarL(0,v0,0), WbarR(0,v0,0);
	Prim delta_l(0,v0,0), delta_i(0,v0,0);
	//tuning parameters:
	double omega = 0;
	//string scheme = "UltraBee";
	
	
	Prim Wi(0,v0,0), Wr(0,v0,0);
	Prim WL(0,v0,0), WR(0,v0,0);
    
    //4 cons. states:
	ConsU UL(0,v0,0), UR(0,v0,0);
	ConsU ULstar(0,v0,0), URstar(0,v0,0);
	//4 cons. fluxes:   
	ConsF FL(0,v0,0), FR(0,v0,0);
	ConsF FLstar(0,v0,0), FRstar(0,v0,0);
	ConsF FLstarHLLC(0,v0,0), FRstarHLLC(0,v0,0);
	ConsF FstarHLL(0,v0,0);
	double epsF, epsG;
	double theta, alpha;
	double tol = 1e-6;
	//Determined flux: (F_i-1/2):                      
	ConsF Fi(0,v0,0);        
    
    double aL, aR;
    double SL, SR, Sstar, Splus;    //wave speeds used in HLLC
    
    int LL,L,I,R;
    
    /*careful with indexing:
	 
	F-2	 F-1   F0   ..   Fi  Fi+1  ..  Fn+1  Fn+2  Fn+3
	............ ____ ____ ____ ____ ____ ..... .....
	|  G  |  G  | S  | S  | S  | S  | S  |  G  |  G  |
	|..-2.|..-1.|__0_|_.._|__i_|_.._|__n_|..+1.|..+2.| 
	
	*/
    
	//assumes NUM_GROW >= 2
    
    for(int k = lo.z; k <= hi.z; ++k) {
        for(int j = lo.y; j <= hi.y; ++j) {
            for(int i = lo.x; i <= hi.x+1; ++i) { 
				LL = i-2;
				L = i-1;
				I = i;
				R = i+1;
				//-------------------convert S to W---------------------------
				if(MUSCL){
				v0[0] = S_old(LL,j,k,V["u"]);
				v0[1] = S_old(LL,j,k,V["v"]);
				Wll = Prim(S_old(LL,j,k,V["rho"]), v0, S_old(LL,j,k,V["p"]));
				v0[0] = S_old(R,j,k,V["u"]);
				v0[1] = S_old(R,j,k,V["v"]);
				Wr = Prim(S_old(R,j,k,V["rho"]), v0, S_old(R,j,k,V["p"]));
				}
				
				v0[0] = S_old(L,j,k,V["u"]);
				v0[1] = S_old(L,j,k,V["v"]);
				Wl = Prim(S_old(L,j,k,V["rho"]), v0, S_old(L,j,k,V["p"]));
				v0[0] = S_old(I,j,k,V["u"]);
				v0[1] = S_old(I,j,k,V["v"]);
				Wi = Prim(S_old(I,j,k,V["rho"]), v0, S_old(I,j,k,V["p"]));
				//-------------------convert S to W---------------------------
				
				
				//--------------------------------------MUSCL half time step----------------------------------------------
				if(MUSCL){
				delta_l = delta(omega, Wll, Wl, Wi, scheme);
				delta_i = delta(omega, Wl, Wi, Wr, scheme);
				
				//WbarL:
				WbarL.rho = Wl.rho + 0.5*((1-dt/dx*Wl.u_vec[0])*delta_l.rho - dt/dx*Wl.rho*delta_l.u_vec[0]);
				WbarL.u_vec[0] = Wl.u_vec[0] + 0.5*((1-dt/dx*Wl.u_vec[0])*delta_l.u_vec[0] - dt/dx*((1/Wl.rho)*delta_l.p));
				WbarL.u_vec[1] = Wl.u_vec[1] + 0.5*((1-dt/dx*Wl.u_vec[0])*delta_l.u_vec[1]);
				WbarL.u_vec[2] = Wl.u_vec[2] + 0.5*((1-dt/dx*Wl.u_vec[0])*delta_l.u_vec[2]);
				S_old(L,j,k,V["EoS19"]) = EoS_properties(EoS, Wl.rho, Wl.p, 0.0, aL, "a");
				//aL = pow(gamma*(Wl.p)/Wl.rho, 0.5);
				WbarL.p = Wl.p + 0.5*((1-dt/dx*Wl.u_vec[0])*delta_l.p - dt/dx*(Wl.rho*pow(aL,2))*delta_l.u_vec[0]);					
				
				//WbarR:
				WbarR.rho = Wi.rho - 0.5*((1+dt/dx*Wi.u_vec[0])*delta_i.rho + dt/dx*Wi.rho*delta_i.u_vec[0]);
				WbarR.u_vec[0] = Wi.u_vec[0] - 0.5*((1+dt/dx*Wi.u_vec[0])*delta_i.u_vec[0] + dt/dx*((1/Wi.rho)*delta_i.p));
				WbarR.u_vec[1] = Wi.u_vec[1] - 0.5*((1+dt/dx*Wi.u_vec[0])*delta_i.u_vec[1]);
				WbarR.u_vec[2] = Wi.u_vec[2] - 0.5*((1+dt/dx*Wi.u_vec[0])*delta_i.u_vec[2]);
				S_old(I,j,k,V["EoS19"]) = EoS_properties(EoS, Wi.rho, Wi.p, 0.0, aR, "a");
				//aR = pow(gamma*(Wi.p)/Wi.rho, 0.5);
				WbarR.p = Wi.p - 0.5*((1+dt/dx*Wi.u_vec[0])*delta_i.p + dt/dx*(Wi.rho*pow(aL,2))*delta_i.u_vec[0]);
				
				WL = WbarL; WR = WbarR;
				}else{
					WL = Wl; WR = Wi;
				}
					
				//-------------------------------------------------------------------------------------------------------
				
				//direct wave speed estimates:
				S_old(L,j,k,V["EoS19"]) = EoS_properties(EoS, WL.rho, WL.p, 0.0, aL, "a");
				//cout << "aL = " << aL << endl;
				S_old(I,j,k,V["EoS19"]) = EoS_properties(EoS, WR.rho, WR.p, 0.0, aR, "a");
				
				SL = ((WL.u_vec[0] - aL < WR.u_vec[0] - aR)? WL.u_vec[0] - aL : WR.u_vec[0] - aR);
				SR = ((WL.u_vec[0] + aL > WR.u_vec[0] + aR)? WL.u_vec[0] + aL : WR.u_vec[0] + aR);
				Sstar = (WR.u_vec[0]*WR.rho*(SR-WR.u_vec[0]) - WL.u_vec[0]*WL.rho*(SL-WL.u_vec[0]) + (WL.p-WR.p))\
						/(WR.rho*(SR-WR.u_vec[0])-WL.rho*(SL-WL.u_vec[0]));
				Splus = ((fabs(WL.u_vec[0])+aL > fabs(WR.u_vec[0])+aR)? \
						  fabs(WL.u_vec[0])+aL : fabs(WR.u_vec[0])+aR);
						  
						  
				//W-->U and W-->F conversions:
				UL = conWtoU(WL);
				UR = conWtoU(WR);
				FL = conWtoF(WL);
				FR = conWtoF(WR);
				
				//flux evaluation at (x/t)=0:	
				if( 0 < SL){
					Fi = FL;
				}else if(SL <= 0 && 0 < SR){
					
					FstarHLL.mass = (SR*FL.mass - SL*FR.mass + SL*SR*(UR.rho-UL.rho))/(SR-SL);
					FstarHLL.mom_vec[0] = (SR*FL.mom_vec[0] - SL*FR.mom_vec[0] + SL*SR*(UR.rhou_vec[0]-UL.rhou_vec[0]))/(SR-SL);
					FstarHLL.mom_vec[1] = (SR*FL.mom_vec[1] - SL*FR.mom_vec[1] + SL*SR*(UR.rhou_vec[1]-UL.rhou_vec[1]))/(SR-SL);
					FstarHLL.en = (SR*FL.en - SL*FR.en + SL*SR*(UR.E-UL.E))/(SR-SL);
					
					if(S_old(i,j,k,V["level_set"]) > 0){
						epsF = 0;
					}else{
						epsF = compute_eps_x(S_old, i, j, k, V);
					}
					
					if(S_old(i,j,k,V["level_set"]) > 0){
						epsG = 0;
					}else{
						epsG = compute_eps_y(S_old, i, j, k, V);
					}
									
					
					theta = atan(epsG/(epsF+tol));
					alpha = cos(theta)*((epsF/aR)/(1 + epsF/aR));
				
				
					if(SL <= 0 && 0 < Sstar){
						//cout << "left star state \n";
						//left star state terms:
						ULstar.rho = WL.rho*((SL-WL.u_vec[0])/(SL-Sstar));
						ULstar.rhou_vec[0] = WL.rho*Sstar*((SL-WL.u_vec[0])/(SL-Sstar));
						ULstar.rhou_vec[1] = WL.rho*WL.u_vec[1]*((SL-WL.u_vec[0])/(SL-Sstar)); 
						ULstar.E = WL.rho*((SL-WL.u_vec[0])/(SL-Sstar))*(UL.E/WL.rho+(Sstar-WL.u_vec[0])\
									*(Sstar+WL.p/(WL.rho*(SL-WL.u_vec[0]))));
						if(ULstar.E < 0){
							cout << "negative energy computed in cell" << endl;
							amrex::Abort("aborted in MUSCL_HLLC_x");
							break;
						}
						//left star state flux:
						FLstarHLLC.mass = FL.mass + SL*(ULstar.rho - UL.rho);
						FLstarHLLC.mom_vec = FL.mom_vec + SL*(ULstar.rhou_vec - UL.rhou_vec);  //check vector addition
						FLstarHLLC.en = FL.en + SL*(ULstar.E - UL.E);
					
					
						Fi.mass = (1-alpha)*FLstarHLLC.mass + alpha*FstarHLL.mass;
						Fi.mom_vec[0] = (1-alpha)*FLstarHLLC.mom_vec[0] + alpha*FstarHLL.mom_vec[0];
						Fi.mom_vec[1] = (1-alpha)*FLstarHLLC.mom_vec[1] + alpha*FstarHLL.mom_vec[1];
						Fi.en = (1-alpha)*FLstarHLLC.en + alpha*FstarHLL.en;
					
					}else if(Sstar <= 0 && 0 < SR){
						//cout << "right star state \n";
						//right star state terms:
						URstar.rho = WR.rho*((SR-WR.u_vec[0])/(SR-Sstar));
						URstar.rhou_vec[0] = WR.rho*Sstar*((SR-WR.u_vec[0])/(SR-Sstar));
						URstar.rhou_vec[1] = WR.rho*WR.u_vec[1]*((SR-WR.u_vec[0])/(SR-Sstar)); 
						URstar.E = WR.rho*((SR-WR.u_vec[0])/(SR-Sstar))*(UR.E/WR.rho+(Sstar-WR.u_vec[0])\
									*(Sstar+WR.p/(WR.rho*(SR-WR.u_vec[0]))));
						if(ULstar.E < 0){
							cout << "negative energy computed in cell" << endl;
							amrex::Abort("aborted in MUSCL_HLLC_x");
							break;
						}
						FRstarHLLC.mass = FR.mass + SR*(URstar.rho - UR.rho);
						FRstarHLLC.mom_vec = FR.mom_vec + SR*(URstar.rhou_vec - UR.rhou_vec); //check vector addition
						FRstarHLLC.en = FR.en + SR*(URstar.E - UR.E);
					
						Fi.mass = (1-alpha)*FRstarHLLC.mass + alpha*FstarHLL.mass;
						Fi.mom_vec[0] = (1-alpha)*FRstarHLLC.mom_vec[0] + alpha*FstarHLL.mom_vec[0];
						Fi.mom_vec[1] = (1-alpha)*FRstarHLLC.mom_vec[1] + alpha*FstarHLL.mom_vec[1];
						Fi.en = (1-alpha)*FRstarHLLC.en + alpha*FstarHLL.en;
					}
							
				}else if(SR <= 0){
					Fi = FR; 
				}
				
				//-------------------convert F to flux_arr------------------------
				for(int n = 0; n<nc; ++n){
					F(i,j,k,n) = 0; //default set all fluxes to zero
				}
				//overwrite calculated fluxes:
				F(i,j,k,V["mass"]) = Fi.mass;
				F(i,j,k,V["mom_x"]) = Fi.mom_vec[0];
				F(i,j,k,V["mom_y"]) = Fi.mom_vec[1];
				F(i,j,k,V["en"]) = Fi.en;
				//----------------------------------------------------------------										
			
			}
		}
	}//end i,j,k
	
}


void calc_fluxes_HLL_x(Box const& box, Array4<Real> const& S_old, Array4<Real> const& F, 
					Real const dt, Real const dx, int const nc, AccessVariable& V, bool MUSCL)
{
	const auto lo = lbound(box);
    const auto hi = ubound(box);
    
    Vector<Real> v0(3);  

	//MUSCL-Hancock parameters (initialised):    
	Prim Wll(0,v0,0), Wl(0,v0,0); 
	Prim WbarL(0,v0,0), WbarR(0,v0,0);
	Prim delta_l(0,v0,0), delta_i(0,v0,0);
	//tuning parameters:
	double omega = 0;
	//string scheme = "UltraBee";
	
	
	Prim Wi(0,v0,0), Wr(0,v0,0);
	Prim WL(0,v0,0), WR(0,v0,0);
    
    //4 cons. states:
	ConsU UL(0,v0,0), UR(0,v0,0);
	//4 cons. fluxes:   
	ConsF FL(0,v0,0), FR(0,v0,0);
	ConsF FstarHLL(0,v0,0);
	//Determined flux: (F_i-1/2):                      
	ConsF Fi(0,v0,0);        
    
    double aL, aR;
    double SL, SR, Sstar, Splus;    //wave speeds used in HLLC
    
    int LL,L,I,R;
    
    /*careful with indexing:
	 
	F-2	 F-1   F0   ..   Fi  Fi+1  ..  Fn+1  Fn+2  Fn+3
	............ ____ ____ ____ ____ ____ ..... .....
	|  G  |  G  | S  | S  | S  | S  | S  |  G  |  G  |
	|..-2.|..-1.|__0_|_.._|__i_|_.._|__n_|..+1.|..+2.| 
	
	*/
    
	//assumes NUM_GROW >= 2
    
    for(int k = lo.z; k <= hi.z; ++k) {
        for(int j = lo.y; j <= hi.y; ++j) {
            for(int i = lo.x; i <= hi.x+1; ++i) { 
				LL = i-2;
				L = i-1;
				I = i;
				R = i+1;
				//-------------------convert S to W---------------------------
				if(MUSCL){
				v0[0] = S_old(LL,j,k,V["u"]);
				v0[1] = S_old(LL,j,k,V["v"]);
				Wll = Prim(S_old(LL,j,k,V["rho"]), v0, S_old(LL,j,k,V["p"]));
				v0[0] = S_old(R,j,k,V["u"]);
				v0[1] = S_old(R,j,k,V["v"]);
				Wr = Prim(S_old(R,j,k,V["rho"]), v0, S_old(R,j,k,V["p"]));
				}
				
				v0[0] = S_old(L,j,k,V["u"]);
				v0[1] = S_old(L,j,k,V["v"]);
				Wl = Prim(S_old(L,j,k,V["rho"]), v0, S_old(L,j,k,V["p"]));
				v0[0] = S_old(I,j,k,V["u"]);
				v0[1] = S_old(I,j,k,V["v"]);
				Wi = Prim(S_old(I,j,k,V["rho"]), v0, S_old(I,j,k,V["p"]));
				//-------------------convert S to W---------------------------
				
				
				//--------------------------------------MUSCL half time step----------------------------------------------
				if(MUSCL){
				delta_l = delta(omega, Wll, Wl, Wi, scheme);
				delta_i = delta(omega, Wl, Wi, Wr, scheme);
				
				//WbarL:
				WbarL.rho = Wl.rho + 0.5*((1-dt/dx*Wl.u_vec[0])*delta_l.rho - dt/dx*Wl.rho*delta_l.u_vec[0]);
				WbarL.u_vec[0] = Wl.u_vec[0] + 0.5*((1-dt/dx*Wl.u_vec[0])*delta_l.u_vec[0] - dt/dx*((1/Wl.rho)*delta_l.p));
				WbarL.u_vec[1] = Wl.u_vec[1] + 0.5*((1-dt/dx*Wl.u_vec[0])*delta_l.u_vec[1]);
				WbarL.u_vec[2] = Wl.u_vec[2] + 0.5*((1-dt/dx*Wl.u_vec[0])*delta_l.u_vec[2]);
				S_old(L,j,k,V["EoS19"]) = EoS_properties(EoS, Wl.rho, Wl.p, 0.0, aL, "a");
				//aL = pow(gamma*(Wl.p)/Wl.rho, 0.5);
				WbarL.p = Wl.p + 0.5*((1-dt/dx*Wl.u_vec[0])*delta_l.p - dt/dx*(Wl.rho*pow(aL,2))*delta_l.u_vec[0]);					
				
				//WbarR:
				WbarR.rho = Wi.rho - 0.5*((1+dt/dx*Wi.u_vec[0])*delta_i.rho + dt/dx*Wi.rho*delta_i.u_vec[0]);
				WbarR.u_vec[0] = Wi.u_vec[0] - 0.5*((1+dt/dx*Wi.u_vec[0])*delta_i.u_vec[0] + dt/dx*((1/Wi.rho)*delta_i.p));
				WbarR.u_vec[1] = Wi.u_vec[1] - 0.5*((1+dt/dx*Wi.u_vec[0])*delta_i.u_vec[1]);
				WbarR.u_vec[2] = Wi.u_vec[2] - 0.5*((1+dt/dx*Wi.u_vec[0])*delta_i.u_vec[2]);
				S_old(I,j,k,V["EoS19"]) = EoS_properties(EoS, Wi.rho, Wi.p, 0.0, aR, "a");
				//aR = pow(gamma*(Wi.p)/Wi.rho, 0.5);
				WbarR.p = Wi.p - 0.5*((1+dt/dx*Wi.u_vec[0])*delta_i.p + dt/dx*(Wi.rho*pow(aL,2))*delta_i.u_vec[0]);
				
				WL = WbarL; WR = WbarR;
				}else{
					WL = Wl; WR = Wi;
				}
					
				//-------------------------------------------------------------------------------------------------------
				
				//direct wave speed estimates:
				S_old(L,j,k,V["EoS19"]) = EoS_properties(EoS, WL.rho, WL.p, 0.0, aL, "a");
				//cout << "aL = " << aL << endl;
				S_old(I,j,k,V["EoS19"]) = EoS_properties(EoS, WR.rho, WR.p, 0.0, aR, "a");
				
				SL = ((WL.u_vec[0] - aL < WR.u_vec[0] - aR)? WL.u_vec[0] - aL : WR.u_vec[0] - aR);
				SR = ((WL.u_vec[0] + aL > WR.u_vec[0] + aR)? WL.u_vec[0] + aL : WR.u_vec[0] + aR);
				Sstar = (WR.u_vec[0]*WR.rho*(SR-WR.u_vec[0]) - WL.u_vec[0]*WL.rho*(SL-WL.u_vec[0]) + (WL.p-WR.p))\
						/(WR.rho*(SR-WR.u_vec[0])-WL.rho*(SL-WL.u_vec[0]));
				Splus = ((fabs(WL.u_vec[0])+aL > fabs(WR.u_vec[0])+aR)? \
						  fabs(WL.u_vec[0])+aL : fabs(WR.u_vec[0])+aR);
						  
						  
				//W-->U and W-->F conversions:
				UL = conWtoU(WL);
				UR = conWtoU(WR);
				FL = conWtoF(WL);
				FR = conWtoF(WR);
				
				//flux evaluation at (x/t)=0:	
				if( 0 < SL){
					Fi = FL;
				}else if(SL <= 0 && 0 < Sstar){
					//cout << "left star state \n";
					//left star state terms:
	
					FstarHLL.mass = (SR*FL.mass - SL*FR.mass + SL*SR*(UR.rho-UL.rho))/(SR-SL);
					FstarHLL.mom_vec[0] = (SR*FL.mom_vec[0] - SL*FR.mom_vec[0] + SL*SR*(UR.rhou_vec[0]-UL.rhou_vec[0]))/(SR-SL);
					FstarHLL.mom_vec[1] = (SR*FL.mom_vec[1] - SL*FR.mom_vec[1] + SL*SR*(UR.rhou_vec[1]-UL.rhou_vec[1]))/(SR-SL);
					FstarHLL.en = (SR*FL.en - SL*FR.en + SL*SR*(UR.E-UL.E))/(SR-SL);
										
					Fi = FstarHLL;
					
				}else if(Sstar <= 0 && 0 < SR){
					//cout << "right star state \n";
					//right star state terms:
					
					FstarHLL.mass = (SR*FL.mass - SL*FR.mass + SL*SR*(UR.rho-UL.rho))/(SR-SL);
					FstarHLL.mom_vec[0] = (SR*FL.mom_vec[0] - SL*FR.mom_vec[0] + SL*SR*(UR.rhou_vec[0]-UL.rhou_vec[0]))/(SR-SL);
					FstarHLL.mom_vec[1] = (SR*FL.mom_vec[1] - SL*FR.mom_vec[1] + SL*SR*(UR.rhou_vec[1]-UL.rhou_vec[1]))/(SR-SL);
					FstarHLL.en = (SR*FL.en - SL*FR.en + SL*SR*(UR.E-UL.E))/(SR-SL);
					
					Fi = FstarHLL;
					
				}else if(SR <= 0){
					Fi = FR; 
				}
				
				//-------------------convert F to flux_arr------------------------
				for(int n = 0; n<nc; ++n){
					F(i,j,k,n) = 0; //default set all fluxes to zero
				}
				//overwrite calculated fluxes:
				F(i,j,k,V["mass"]) = Fi.mass;
				F(i,j,k,V["mom_x"]) = Fi.mom_vec[0];
				F(i,j,k,V["mom_y"]) = Fi.mom_vec[1];
				F(i,j,k,V["en"]) = Fi.en;
				//----------------------------------------------------------------										
			
			}
		}
	}//end i,j,k
	
}


    

void calc_fluxes_HLLC_x(Box const& box, Array4<Real> const& S_old, Array4<Real> const& F, 
					Real const dt, Real const dx, int const nc, AccessVariable& V, bool MUSCL)
{
	const auto lo = lbound(box);
    const auto hi = ubound(box);
    
    Vector<Real> v0(3);  

	//MUSCL-Hancock parameters (initialised):    
	Prim Wll(0,v0,0), Wl(0,v0,0); 
	Prim WbarL(0,v0,0), WbarR(0,v0,0);
	Prim delta_l(0,v0,0), delta_i(0,v0,0);
	//tuning parameters:
	double omega = 0;
	//string scheme = "UltraBee";
	
	
	Prim Wi(0,v0,0), Wr(0,v0,0);
	Prim WL(0,v0,0), WR(0,v0,0);
    
    //4 cons. states:
	ConsU UL(0,v0,0), UR(0,v0,0);
	ConsU ULstar(0,v0,0), URstar(0,v0,0);
	//4 cons. fluxes:   
	ConsF FL(0,v0,0), FR(0,v0,0);
	ConsF FLstar(0,v0,0), FRstar(0,v0,0);
	//Determined flux: (F_i-1/2):                      
	ConsF Fi(0,v0,0);        
    
    double aL, aR;
    double SL, SR, Sstar, Splus;    //wave speeds used in HLLC
    
    int LL,L,I,R;
    
    /*careful with indexing:
	 
	F-2	 F-1   F0   ..   Fi  Fi+1  ..  Fn+1  Fn+2  Fn+3
	............ ____ ____ ____ ____ ____ ..... .....
	|  G  |  G  | S  | S  | S  | S  | S  |  G  |  G  |
	|..-2.|..-1.|__0_|_.._|__i_|_.._|__n_|..+1.|..+2.| 
	
	*/
    
	//assumes NUM_GROW >= 2
    
    for(int k = lo.z; k <= hi.z; ++k) {
        for(int j = lo.y; j <= hi.y; ++j) {
            for(int i = lo.x; i <= hi.x+1; ++i) { 
				LL = i-2;
				L = i-1;
				I = i;
				R = i+1;
				//-------------------convert S to W---------------------------
				if(MUSCL){
				v0[0] = S_old(LL,j,k,V["u"]);
				v0[1] = S_old(LL,j,k,V["v"]);
				Wll = Prim(S_old(LL,j,k,V["rho"]), v0, S_old(LL,j,k,V["p"]));
				v0[0] = S_old(R,j,k,V["u"]);
				v0[1] = S_old(R,j,k,V["v"]);
				Wr = Prim(S_old(R,j,k,V["rho"]), v0, S_old(R,j,k,V["p"]));
				}
				
				v0[0] = S_old(L,j,k,V["u"]);
				v0[1] = S_old(L,j,k,V["v"]);
				Wl = Prim(S_old(L,j,k,V["rho"]), v0, S_old(L,j,k,V["p"]));
				v0[0] = S_old(I,j,k,V["u"]);
				v0[1] = S_old(I,j,k,V["v"]);
				Wi = Prim(S_old(I,j,k,V["rho"]), v0, S_old(I,j,k,V["p"]));
				//-------------------convert S to W---------------------------
				
				
				//--------------------------------------MUSCL half time step----------------------------------------------
				if(MUSCL){
				delta_l = delta(omega, Wll, Wl, Wi, scheme);
				delta_i = delta(omega, Wl, Wi, Wr, scheme);
				
				//WbarL:
				WbarL.rho = Wl.rho + 0.5*((1-dt/dx*Wl.u_vec[0])*delta_l.rho - dt/dx*Wl.rho*delta_l.u_vec[0]);
				WbarL.u_vec[0] = Wl.u_vec[0] + 0.5*((1-dt/dx*Wl.u_vec[0])*delta_l.u_vec[0] - dt/dx*((1/Wl.rho)*delta_l.p));
				WbarL.u_vec[1] = Wl.u_vec[1] + 0.5*((1-dt/dx*Wl.u_vec[0])*delta_l.u_vec[1]);
				WbarL.u_vec[2] = Wl.u_vec[2] + 0.5*((1-dt/dx*Wl.u_vec[0])*delta_l.u_vec[2]);
				S_old(L,j,k,V["EoS19"]) = EoS_properties(EoS, Wl.rho, Wl.p, 0.0, aL, "a");
				//aL = pow(gamma*(Wl.p)/Wl.rho, 0.5);
				WbarL.p = Wl.p + 0.5*((1-dt/dx*Wl.u_vec[0])*delta_l.p - dt/dx*(Wl.rho*pow(aL,2))*delta_l.u_vec[0]);					
				
				//WbarR:
				WbarR.rho = Wi.rho - 0.5*((1+dt/dx*Wi.u_vec[0])*delta_i.rho + dt/dx*Wi.rho*delta_i.u_vec[0]);
				WbarR.u_vec[0] = Wi.u_vec[0] - 0.5*((1+dt/dx*Wi.u_vec[0])*delta_i.u_vec[0] + dt/dx*((1/Wi.rho)*delta_i.p));
				WbarR.u_vec[1] = Wi.u_vec[1] - 0.5*((1+dt/dx*Wi.u_vec[0])*delta_i.u_vec[1]);
				WbarR.u_vec[2] = Wi.u_vec[2] - 0.5*((1+dt/dx*Wi.u_vec[0])*delta_i.u_vec[2]);
				S_old(I,j,k,V["EoS19"]) = EoS_properties(EoS, Wi.rho, Wi.p, 0.0, aR, "a");
				//aR = pow(gamma*(Wi.p)/Wi.rho, 0.5);
				WbarR.p = Wi.p - 0.5*((1+dt/dx*Wi.u_vec[0])*delta_i.p + dt/dx*(Wi.rho*pow(aL,2))*delta_i.u_vec[0]);
				
				WL = WbarL; WR = WbarR;
				}else{
					WL = Wl; WR = Wi;
				}
					
				//-------------------------------------------------------------------------------------------------------
				
				//direct wave speed estimates:
				S_old(L,j,k,V["EoS19"]) = EoS_properties(EoS, WL.rho, WL.p, 0.0, aL, "a");
				//cout << "aL = " << aL << endl;
				S_old(I,j,k,V["EoS19"]) = EoS_properties(EoS, WR.rho, WR.p, 0.0, aR, "a");
				
				SL = ((WL.u_vec[0] - aL < WR.u_vec[0] - aR)? WL.u_vec[0] - aL : WR.u_vec[0] - aR);
				SR = ((WL.u_vec[0] + aL > WR.u_vec[0] + aR)? WL.u_vec[0] + aL : WR.u_vec[0] + aR);
				Sstar = (WR.u_vec[0]*WR.rho*(SR-WR.u_vec[0]) - WL.u_vec[0]*WL.rho*(SL-WL.u_vec[0]) + (WL.p-WR.p))\
						/(WR.rho*(SR-WR.u_vec[0])-WL.rho*(SL-WL.u_vec[0]));
				Splus = ((fabs(WL.u_vec[0])+aL > fabs(WR.u_vec[0])+aR)? \
						  fabs(WL.u_vec[0])+aL : fabs(WR.u_vec[0])+aR);
						  
						  
				//W-->U and W-->F conversions:
				UL = conWtoU(WL);
				UR = conWtoU(WR);
				FL = conWtoF(WL);
				FR = conWtoF(WR);
				
				//flux evaluation at (x/t)=0:	
				if( 0 < SL){
					Fi = FL;
				}else if(SL <= 0 && 0 < Sstar){
					//cout << "left star state \n";
					//left star state terms:
					ULstar.rho = WL.rho*((SL-WL.u_vec[0])/(SL-Sstar));
					ULstar.rhou_vec[0] = WL.rho*Sstar*((SL-WL.u_vec[0])/(SL-Sstar));
					ULstar.rhou_vec[1] = WL.rho*WL.u_vec[1]*((SL-WL.u_vec[0])/(SL-Sstar)); 
					ULstar.E = WL.rho*((SL-WL.u_vec[0])/(SL-Sstar))*(UL.E/WL.rho+(Sstar-WL.u_vec[0])\
								*(Sstar+WL.p/(WL.rho*(SL-WL.u_vec[0]))));
					if(ULstar.E < 0){
						cout << "negative energy computed in cell" << endl;
						amrex::Abort("aborted in MUSCL_HLLC_x");
						break;
					}
					//left star state flux:
					FLstar.mass = FL.mass + SL*(ULstar.rho - UL.rho);
					FLstar.mom_vec = FL.mom_vec + SL*(ULstar.rhou_vec - UL.rhou_vec);  //check vector addition
					FLstar.en = FL.en + SL*(ULstar.E - UL.E);
					Fi = FLstar;
				}else if(Sstar <= 0 && 0 < SR){
					//cout << "right star state \n";
					//right star state terms:
					URstar.rho = WR.rho*((SR-WR.u_vec[0])/(SR-Sstar));
					URstar.rhou_vec[0] = WR.rho*Sstar*((SR-WR.u_vec[0])/(SR-Sstar));
					URstar.rhou_vec[1] = WR.rho*WR.u_vec[1]*((SR-WR.u_vec[0])/(SR-Sstar)); 
					URstar.E = WR.rho*((SR-WR.u_vec[0])/(SR-Sstar))*(UR.E/WR.rho+(Sstar-WR.u_vec[0])\
								*(Sstar+WR.p/(WR.rho*(SR-WR.u_vec[0]))));
					if(ULstar.E < 0){
						cout << "negative energy computed in cell" << endl;
						amrex::Abort("aborted in MUSCL_HLLC_x");
						break;
					}
					FRstar.mass = FR.mass + SR*(URstar.rho - UR.rho);
					FRstar.mom_vec = FR.mom_vec + SR*(URstar.rhou_vec - UR.rhou_vec); //check vector addition
					FRstar.en = FR.en + SR*(URstar.E - UR.E);
					Fi = FRstar;
				}else if(SR <= 0){
					Fi = FR; 
				}
				
				//-------------------convert F to flux_arr------------------------
				for(int n = 0; n<nc; ++n){
					F(i,j,k,n) = 0; //default set all fluxes to zero
				}
				//overwrite calculated fluxes:
				F(i,j,k,V["mass"]) = Fi.mass;
				F(i,j,k,V["mom_x"]) = Fi.mom_vec[0];
				F(i,j,k,V["mom_y"]) = Fi.mom_vec[1];
				F(i,j,k,V["en"]) = Fi.en;
				//----------------------------------------------------------------										
			
			}
		}
	}//end i,j,k
	
}

void calc_fluxes_HLLC_HS_y(Box const& box, Array4<Real> const& S_old, Array4<Real> const& G, 
					Real const dt, Real const dy, int const nc, AccessVariable& V, bool MUSCL)
{
	const auto lo = lbound(box);
    const auto hi = ubound(box);
    
    Vector<Real> v0(3);  
    
  
	//MUSCL-Hancock parameters (initialised):    
	Prim Wll(0,v0,0), Wl(0,v0,0); 
	Prim WbarL(0,v0,0), WbarR(0,v0,0);
	Prim delta_l(0,v0,0), delta_i(0,v0,0);
	//tuning parameters:
	double omega = 0;
	//string scheme = "UltraBee";
	
	
	Prim Wi(0,v0,0), Wr(0,v0,0);
	Prim WL(0,v0,0), WR(0,v0,0);
    
    //4 cons. states:
	ConsU UL(0,v0,0), UR(0,v0,0);
	ConsU ULstar(0,v0,0), URstar(0,v0,0);
	//4 cons. fluxes:   
	ConsF GL(0,v0,0), GR(0,v0,0);
	ConsF GLstar(0,v0,0), GRstar(0,v0,0);
	//SWM parameters:
	ConsF GLstarHLLC(0,v0,0), GRstarHLLC(0,v0,0);
	ConsF GstarHLL(0,v0,0);
	double epsF, epsG;
	double theta, alpha;
	double tol = 1e-6;
	//Determined flux: (F_i-1/2):                      
	ConsF Gi(0,v0,0);    
	 
    
    double aL, aR;
    double SL, SR, Sstar, Splus;    //wave speeds used in HLLC

    int LL,L,I,R;
    
    /*careful with indexing:
	 
	F-2	 F-1   F0   ..   Fi  Fi+1  ..  Fn+1  Fn+2  Fn+3
	............ ____ ____ ____ ____ ____ ..... .....
	|  G  |  G  | S  | S  | S  | S  | S  |  G  |  G  |
	|..-2.|..-1.|__0_|_.._|__i_|_.._|__n_|..+1.|..+2.| 
	
	*/
    
	//assumes NUM_GROW >= 2
    
    for(int k = lo.z; k <= hi.z; ++k) {
        for(int j = lo.y; j <= hi.y+1; ++j) {
            for(int i = lo.x; i <= hi.x; ++i) { 
				LL = j-2;
				L = j-1;
				I = j;
				R = j+1;
				//-------------------convert S to W---------------------------
				if(MUSCL){
				v0[0] = S_old(i,LL,k,V["u"]);
				v0[1] = S_old(i,LL,k,V["v"]);
				Wll = Prim(S_old(i,LL,k,V["rho"]), v0, S_old(i,LL,k,V["p"]));
				v0[0] = S_old(i,R,k,V["u"]);
				v0[1] = S_old(i,R,k,V["v"]);
				Wr = Prim(S_old(i,R,k,V["rho"]), v0, S_old(i,R,k,V["p"]));
				}
				
				v0[0] = S_old(i,L,k,V["u"]);
				v0[1] = S_old(i,L,k,V["v"]);
				Wl = Prim(S_old(i,L,k,V["rho"]), v0, S_old(i,L,k,V["p"]));
				v0[0] = S_old(i,I,k,V["u"]);
				v0[1] = S_old(i,I,k,V["v"]);
				Wi = Prim(S_old(i,I,k,V["rho"]), v0, S_old(i,I,k,V["p"]));
				//-------------------convert S to W---------------------------
				
				
				//--------------------------------------MUSCL half time step----------------------------------------------
				if(MUSCL){
				delta_l = delta(omega, Wll, Wl, Wi, scheme);
				delta_i = delta(omega, Wl, Wi, Wr, scheme);
				
				//WbarL:
				WbarL.rho = Wl.rho + 0.5*((1-dt/dy*Wl.u_vec[1])*delta_l.rho - dt/dy*Wl.rho*delta_l.u_vec[1]);
				WbarL.u_vec[0] = Wl.u_vec[0] + 0.5*((1-dt/dy*Wl.u_vec[1])*delta_l.u_vec[0]);
				WbarL.u_vec[1] = Wl.u_vec[1] + 0.5*((1-dt/dy*Wl.u_vec[1])*delta_l.u_vec[1] - dt/dy*((1/Wl.rho)*delta_l.p));
				WbarL.u_vec[2] = Wl.u_vec[2] + 0.5*((1-dt/dy*Wl.u_vec[1])*delta_l.u_vec[2]);
				S_old(i,L,k,V["EoS19"]) = EoS_properties(EoS, Wl.rho, Wl.p, 0.0, aL, "a");
				//aL = pow(gamma*(Wl.p)/Wl.rho, 0.5);
				WbarL.p = Wl.p + 0.5*((1-dt/dy*Wl.u_vec[1])*delta_l.p - dt/dy*(Wl.rho*pow(aL,2))*delta_l.u_vec[1]);					
				
				//WbarR:
				WbarR.rho = Wi.rho - 0.5*((1+dt/dy*Wi.u_vec[1])*delta_i.rho + dt/dy*Wi.rho*delta_i.u_vec[1]);
				WbarR.u_vec[0] = Wi.u_vec[0] - 0.5*((1+dt/dy*Wi.u_vec[1])*delta_i.u_vec[0]);
				WbarR.u_vec[1] = Wi.u_vec[1] - 0.5*((1+dt/dy*Wi.u_vec[1])*delta_i.u_vec[1] + dt/dy*((1/Wi.rho)*delta_i.p));
				WbarR.u_vec[2] = Wi.u_vec[2] - 0.5*((1+dt/dy*Wi.u_vec[1])*delta_i.u_vec[2]);
				S_old(i,I,k,V["EoS19"]) = EoS_properties(EoS, Wi.rho, Wi.p, 0.0, aR, "a");
				//aR = pow(gamma*(Wi.p)/Wi.rho, 0.5);
				WbarR.p = Wi.p - 0.5*((1+dt/dy*Wi.u_vec[1])*delta_i.p + dt/dy*(Wi.rho*pow(aL,2))*delta_i.u_vec[1]);
				
				WL = WbarL; WR = WbarR;
				}else{
				WL = Wl; WR = Wi;
				}
				//--------------------------------------MUSCL half time step----------------------------------------------
				//direct wave speed estimates:
				S_old(i,L,k,V["EoS19"]) = EoS_properties(EoS, WL.rho, WL.p, 0.0, aL, "a");
				//cout << "aL = " << aL << endl;
				S_old(i,I,k,V["EoS19"]) = EoS_properties(EoS, WR.rho, WR.p, 0.0, aR, "a");
			
				SL = ((WL.u_vec[1] - aL < WR.u_vec[1] - aR)? WL.u_vec[1] - aL : WR.u_vec[1] - aR);
				SR = ((WL.u_vec[1] + aL > WR.u_vec[1] + aR)? WL.u_vec[1] + aL : WR.u_vec[1] + aR);
				Sstar = (WR.u_vec[1]*WR.rho*(SR-WR.u_vec[1]) - WL.u_vec[1]*WL.rho*(SL-WL.u_vec[1]) + (WL.p-WR.p))\
						/(WR.rho*(SR-WR.u_vec[1])-WL.rho*(SL-WL.u_vec[1]));
				Splus = ((fabs(WL.u_vec[1])+aL > fabs(WR.u_vec[1])+aR)? \
						  fabs(WL.u_vec[1])+aL : fabs(WR.u_vec[1])+aR);
						  
						  
				//W-->U and W-->G conversions:
				UL = conWtoU(WL);
				UR = conWtoU(WR);
				GL = conWtoG(WL);
				GR = conWtoG(WR);
				
				//flux evaluation at (x/t)=0:	
				if( 0 < SL){
					Gi = GL;
				}else if(SL <= 0 && 0 < SR){
					
					GstarHLL.mass = (SR*GL.mass - SL*GR.mass + SL*SR*(UR.rho-UL.rho))/(SR-SL);
					GstarHLL.mom_vec[0] = (SR*GL.mom_vec[0] - SL*GR.mom_vec[0] + SL*SR*(UR.rhou_vec[0]-UL.rhou_vec[0]))/(SR-SL);
					GstarHLL.mom_vec[1] = (SR*GL.mom_vec[1] - SL*GR.mom_vec[1] + SL*SR*(UR.rhou_vec[1]-UL.rhou_vec[1]))/(SR-SL);
					GstarHLL.en = (SR*GL.en - SL*GR.en + SL*SR*(UR.E-UL.E))/(SR-SL);
					
					if(S_old(i,j,k,V["level_set"]) > 0){
						epsG = 0;
					}else{
						epsG = compute_eps_y(S_old, i, j, k, V);
					}
					
					if(S_old(i,j,k,V["level_set"]) > 0){
						epsF = 0;
					}else{
						epsF = compute_eps_x(S_old, i, j, k, V);
					}					
					
					theta = atan(epsF/(epsG+tol));
					alpha = cos(theta)*((epsG/aR)/(1 + epsG/aR));					
					
					
					if(SL <= 0 && 0 < Sstar){
						//cout << "left star state \n";
						//left star state terms:
						ULstar.rho = WL.rho*((SL-WL.u_vec[1])/(SL-Sstar));
						ULstar.rhou_vec[0] = WL.rho*WL.u_vec[0]*((SL-WL.u_vec[1])/(SL-Sstar));
						ULstar.rhou_vec[1] = ULstar.rho*Sstar; 
						ULstar.E = ((SL-WL.u_vec[1])/(SL-Sstar))*(UL.E+WL.rho*(Sstar-WL.u_vec[1])*\
										(Sstar+WL.p/(WL.rho*(SL-WL.u_vec[1]))));
						if(ULstar.E < 0){
							cout << "negative energy computed in cell" << endl;
							break;
						}
						//left star state flux:
						GLstarHLLC.mass = GL.mass + SL*(ULstar.rho - UL.rho);
						GLstarHLLC.mom_vec = GL.mom_vec + SL*(ULstar.rhou_vec - UL.rhou_vec);  //check vector addition
						GLstarHLLC.en = GL.en + SL*(ULstar.E - UL.E);
						
						Gi.mass = (1-alpha)*GLstarHLLC.mass + alpha*GstarHLL.mass;
						Gi.mom_vec[0] = (1-alpha)*GLstarHLLC.mom_vec[0] + alpha*GstarHLL.mom_vec[0];
						Gi.mom_vec[1] = (1-alpha)*GLstarHLLC.mom_vec[1] + alpha*GstarHLL.mom_vec[1];
						Gi.en = (1-alpha)*GLstarHLLC.en + alpha*GstarHLL.en;
					
					
					}else if(Sstar <= 0 && 0 < SR){
						//cout << "right star state \n";
						//right star state terms:
						URstar.rho = WR.rho*((SR-WR.u_vec[1])/(SR-Sstar));
						URstar.rhou_vec[0] = WR.rho*WR.u_vec[0]*((SR-WR.u_vec[1])/(SR-Sstar));
						URstar.rhou_vec[1] = URstar.rho*Sstar; 
						URstar.E = ((SR-WR.u_vec[1])/(SR-Sstar))*(UR.E + WR.rho*(Sstar-WR.u_vec[1])*\
										(Sstar+WR.p/(WR.rho*(SR-WR.u_vec[1]))));
						if(ULstar.E < 0){
							cout << "negative energy computed in cell" << endl;
							break;
						}
						GRstarHLLC.mass = GR.mass + SR*(URstar.rho - UR.rho);
						GRstarHLLC.mom_vec = GR.mom_vec + SR*(URstar.rhou_vec - UR.rhou_vec); //check vector addition
						GRstarHLLC.en = GR.en + SR*(URstar.E - UR.E);
					
						Gi.mass = (1-alpha)*GRstarHLLC.mass + alpha*GstarHLL.mass;
						Gi.mom_vec[0] = (1-alpha)*GRstarHLLC.mom_vec[0] + alpha*GstarHLL.mom_vec[0];
						Gi.mom_vec[1] = (1-alpha)*GRstarHLLC.mom_vec[1] + alpha*GstarHLL.mom_vec[1];
						Gi.en = (1-alpha)*GRstarHLLC.en + alpha*GstarHLL.en;
					}
					
				}else if(SR <= 0){
					Gi = GR; 
				}
				
				//-------------------convert G to flux_arr------------------------
				for(int n = 0; n<nc; ++n){
					G(i,j,k,n) = 0; //default set all fluxes to zero
				}
				//overwrite computed fluxes:
				G(i,j,k,V["mass"]) = Gi.mass;
				G(i,j,k,V["mom_x"]) = Gi.mom_vec[0];
				G(i,j,k,V["mom_y"]) = Gi.mom_vec[1];
				G(i,j,k,V["en"]) = Gi.en;
				//----------------------------------------------------------------		
				
				
			}
		}
	}//end i,j,k
	
}//end HLLC_y

void calc_fluxes_HLL_y(Box const& box, Array4<Real> const& S_old, Array4<Real> const& G, 
					Real const dt, Real const dy, int const nc, AccessVariable& V, bool MUSCL)
{
	const auto lo = lbound(box);
    const auto hi = ubound(box);
    
    Vector<Real> v0(3);  
    
  
	//MUSCL-Hancock parameters (initialised):    
	Prim Wll(0,v0,0), Wl(0,v0,0); 
	Prim WbarL(0,v0,0), WbarR(0,v0,0);
	Prim delta_l(0,v0,0), delta_i(0,v0,0);
	//tuning parameters:
	double omega = 0;
	//string scheme = "UltraBee";
	
	
	Prim Wi(0,v0,0), Wr(0,v0,0);
	Prim WL(0,v0,0), WR(0,v0,0);
    
    //4 cons. states:
	ConsU UL(0,v0,0), UR(0,v0,0);
	ConsU ULstar(0,v0,0), URstar(0,v0,0);
	//4 cons. fluxes:   
	ConsF GL(0,v0,0), GR(0,v0,0);
	ConsF GLstar(0,v0,0), GRstar(0,v0,0);
	//SWM parameters:
	ConsF GstarHLL(0,v0,0);
	//Determined flux: (F_i-1/2):                      
	ConsF Gi(0,v0,0);        
    
    double aL, aR;
    double SL, SR, Sstar, Splus;    //wave speeds used in HLLC

    int LL,L,I,R;
    
    /*careful with indexing:
	 
	F-2	 F-1   F0   ..   Fi  Fi+1  ..  Fn+1  Fn+2  Fn+3
	............ ____ ____ ____ ____ ____ ..... .....
	|  G  |  G  | S  | S  | S  | S  | S  |  G  |  G  |
	|..-2.|..-1.|__0_|_.._|__i_|_.._|__n_|..+1.|..+2.| 
	
	*/
    
	//assumes NUM_GROW >= 2
    
    for(int k = lo.z; k <= hi.z; ++k) {
        for(int j = lo.y; j <= hi.y+1; ++j) {
            for(int i = lo.x; i <= hi.x; ++i) { 
				LL = j-2;
				L = j-1;
				I = j;
				R = j+1;
				//-------------------convert S to W---------------------------
				if(MUSCL){
				v0[0] = S_old(i,LL,k,V["u"]);
				v0[1] = S_old(i,LL,k,V["v"]);
				Wll = Prim(S_old(i,LL,k,V["rho"]), v0, S_old(i,LL,k,V["p"]));
				v0[0] = S_old(i,R,k,V["u"]);
				v0[1] = S_old(i,R,k,V["v"]);
				Wr = Prim(S_old(i,R,k,V["rho"]), v0, S_old(i,R,k,V["p"]));
				}
				
				v0[0] = S_old(i,L,k,V["u"]);
				v0[1] = S_old(i,L,k,V["v"]);
				Wl = Prim(S_old(i,L,k,V["rho"]), v0, S_old(i,L,k,V["p"]));
				v0[0] = S_old(i,I,k,V["u"]);
				v0[1] = S_old(i,I,k,V["v"]);
				Wi = Prim(S_old(i,I,k,V["rho"]), v0, S_old(i,I,k,V["p"]));
				//-------------------convert S to W---------------------------
				
				
				//--------------------------------------MUSCL half time step----------------------------------------------
				if(MUSCL){
				delta_l = delta(omega, Wll, Wl, Wi, scheme);
				delta_i = delta(omega, Wl, Wi, Wr, scheme);
				
				//WbarL:
				WbarL.rho = Wl.rho + 0.5*((1-dt/dy*Wl.u_vec[1])*delta_l.rho - dt/dy*Wl.rho*delta_l.u_vec[1]);
				WbarL.u_vec[0] = Wl.u_vec[0] + 0.5*((1-dt/dy*Wl.u_vec[1])*delta_l.u_vec[0]);
				WbarL.u_vec[1] = Wl.u_vec[1] + 0.5*((1-dt/dy*Wl.u_vec[1])*delta_l.u_vec[1] - dt/dy*((1/Wl.rho)*delta_l.p));
				WbarL.u_vec[2] = Wl.u_vec[2] + 0.5*((1-dt/dy*Wl.u_vec[1])*delta_l.u_vec[2]);
				S_old(i,L,k,V["EoS19"]) = EoS_properties(EoS, Wl.rho, Wl.p, 0.0, aL, "a");
				//aL = pow(gamma*(Wl.p)/Wl.rho, 0.5);
				WbarL.p = Wl.p + 0.5*((1-dt/dy*Wl.u_vec[1])*delta_l.p - dt/dy*(Wl.rho*pow(aL,2))*delta_l.u_vec[1]);					
				
				//WbarR:
				WbarR.rho = Wi.rho - 0.5*((1+dt/dy*Wi.u_vec[1])*delta_i.rho + dt/dy*Wi.rho*delta_i.u_vec[1]);
				WbarR.u_vec[0] = Wi.u_vec[0] - 0.5*((1+dt/dy*Wi.u_vec[1])*delta_i.u_vec[0]);
				WbarR.u_vec[1] = Wi.u_vec[1] - 0.5*((1+dt/dy*Wi.u_vec[1])*delta_i.u_vec[1] + dt/dy*((1/Wi.rho)*delta_i.p));
				WbarR.u_vec[2] = Wi.u_vec[2] - 0.5*((1+dt/dy*Wi.u_vec[1])*delta_i.u_vec[2]);
				S_old(i,I,k,V["EoS19"]) = EoS_properties(EoS, Wi.rho, Wi.p, 0.0, aR, "a");
				//aR = pow(gamma*(Wi.p)/Wi.rho, 0.5);
				WbarR.p = Wi.p - 0.5*((1+dt/dy*Wi.u_vec[1])*delta_i.p + dt/dy*(Wi.rho*pow(aL,2))*delta_i.u_vec[1]);
				
				WL = WbarL; WR = WbarR;
				}else{
				WL = Wl; WR = Wi;
				}
				//--------------------------------------MUSCL half time step----------------------------------------------
				//direct wave speed estimates:
				S_old(i,L,k,V["EoS19"]) = EoS_properties(EoS, WL.rho, WL.p, 0.0, aL, "a");
				//cout << "aL = " << aL << endl;
				S_old(i,I,k,V["EoS19"]) = EoS_properties(EoS, WR.rho, WR.p, 0.0, aR, "a");
			
				SL = ((WL.u_vec[1] - aL < WR.u_vec[1] - aR)? WL.u_vec[1] - aL : WR.u_vec[1] - aR);
				SR = ((WL.u_vec[1] + aL > WR.u_vec[1] + aR)? WL.u_vec[1] + aL : WR.u_vec[1] + aR);
				Sstar = (WR.u_vec[1]*WR.rho*(SR-WR.u_vec[1]) - WL.u_vec[1]*WL.rho*(SL-WL.u_vec[1]) + (WL.p-WR.p))\
						/(WR.rho*(SR-WR.u_vec[1])-WL.rho*(SL-WL.u_vec[1]));
				Splus = ((fabs(WL.u_vec[1])+aL > fabs(WR.u_vec[1])+aR)? \
						  fabs(WL.u_vec[1])+aL : fabs(WR.u_vec[1])+aR);
						  
						  
				//W-->U and W-->G conversions:
				UL = conWtoU(WL);
				UR = conWtoU(WR);
				GL = conWtoG(WL);
				GR = conWtoG(WR);
				
				//flux evaluation at (x/t)=0:	
				if( 0 < SL){
					Gi = GL;
				}else if(SL <= 0 && 0 < Sstar){
					//cout << "left star state \n";
					//left star state terms:
					
					GstarHLL.mass = (SR*GL.mass - SL*GR.mass + SL*SR*(UR.rho-UL.rho))/(SR-SL);
					GstarHLL.mom_vec[0] = (SR*GL.mom_vec[0] - SL*GR.mom_vec[0] + SL*SR*(UR.rhou_vec[0]-UL.rhou_vec[0]))/(SR-SL);
					GstarHLL.mom_vec[1] = (SR*GL.mom_vec[1] - SL*GR.mom_vec[1] + SL*SR*(UR.rhou_vec[1]-UL.rhou_vec[1]))/(SR-SL);
					GstarHLL.en = (SR*GL.en - SL*GR.en + SL*SR*(UR.E-UL.E))/(SR-SL);

					Gi = GstarHLL;
					
					
				}else if(Sstar <= 0 && 0 < SR){
					//cout << "right star state \n";
					//right star state terms:
					
					GstarHLL.mass = (SR*GL.mass - SL*GR.mass + SL*SR*(UR.rho-UL.rho))/(SR-SL);
					GstarHLL.mom_vec[0] = (SR*GL.mom_vec[0] - SL*GR.mom_vec[0] + SL*SR*(UR.rhou_vec[0]-UL.rhou_vec[0]))/(SR-SL);
					GstarHLL.mom_vec[1] = (SR*GL.mom_vec[1] - SL*GR.mom_vec[1] + SL*SR*(UR.rhou_vec[1]-UL.rhou_vec[1]))/(SR-SL);
					GstarHLL.en = (SR*GL.en - SL*GR.en + SL*SR*(UR.E-UL.E))/(SR-SL);
					
					Gi = GstarHLL;
				}else if(SR <= 0){
					Gi = GR; 
				}
				
				//-------------------convert G to flux_arr------------------------
				for(int n = 0; n<nc; ++n){
					G(i,j,k,n) = 0; //default set all fluxes to zero
				}
				//overwrite computed fluxes:
				G(i,j,k,V["mass"]) = Gi.mass;
				G(i,j,k,V["mom_x"]) = Gi.mom_vec[0];
				G(i,j,k,V["mom_y"]) = Gi.mom_vec[1];
				G(i,j,k,V["en"]) = Gi.en;
				//----------------------------------------------------------------		
				
				
			}
		}
	}//end i,j,k
	
}//end HLL_y


void calc_fluxes_HLLC_y(Box const& box, Array4<Real> const& S_old, Array4<Real> const& G, 
					Real const dt, Real const dy, int const nc, AccessVariable& V, bool MUSCL)
{
	const auto lo = lbound(box);
    const auto hi = ubound(box);
    
    Vector<Real> v0(3);  
    
  
	//MUSCL-Hancock parameters (initialised):    
	Prim Wll(0,v0,0), Wl(0,v0,0); 
	Prim WbarL(0,v0,0), WbarR(0,v0,0);
	Prim delta_l(0,v0,0), delta_i(0,v0,0);
	//tuning parameters:
	double omega = 0;
	//string scheme = "UltraBee";
	
	
	Prim Wi(0,v0,0), Wr(0,v0,0);
	Prim WL(0,v0,0), WR(0,v0,0);
    
    //4 cons. states:
	ConsU UL(0,v0,0), UR(0,v0,0);
	ConsU ULstar(0,v0,0), URstar(0,v0,0);
	//4 cons. fluxes:   
	ConsF GL(0,v0,0), GR(0,v0,0);
	ConsF GLstar(0,v0,0), GRstar(0,v0,0);
	//Determined flux: (F_i-1/2):                      
	ConsF Gi(0,v0,0);        
    
    double aL, aR;
    double SL, SR, Sstar, Splus;    //wave speeds used in HLLC

    int LL,L,I,R;
    
    /*careful with indexing:
	 
	F-2	 F-1   F0   ..   Fi  Fi+1  ..  Fn+1  Fn+2  Fn+3
	............ ____ ____ ____ ____ ____ ..... .....
	|  G  |  G  | S  | S  | S  | S  | S  |  G  |  G  |
	|..-2.|..-1.|__0_|_.._|__i_|_.._|__n_|..+1.|..+2.| 
	
	*/
    
	//assumes NUM_GROW >= 2
    
    for(int k = lo.z; k <= hi.z; ++k) {
        for(int j = lo.y; j <= hi.y+1; ++j) {
            for(int i = lo.x; i <= hi.x; ++i) { 
				LL = j-2;
				L = j-1;
				I = j;
				R = j+1;
				//-------------------convert S to W---------------------------
				if(MUSCL){
				v0[0] = S_old(i,LL,k,V["u"]);
				v0[1] = S_old(i,LL,k,V["v"]);
				Wll = Prim(S_old(i,LL,k,V["rho"]), v0, S_old(i,LL,k,V["p"]));
				v0[0] = S_old(i,R,k,V["u"]);
				v0[1] = S_old(i,R,k,V["v"]);
				Wr = Prim(S_old(i,R,k,V["rho"]), v0, S_old(i,R,k,V["p"]));
				}
				
				v0[0] = S_old(i,L,k,V["u"]);
				v0[1] = S_old(i,L,k,V["v"]);
				Wl = Prim(S_old(i,L,k,V["rho"]), v0, S_old(i,L,k,V["p"]));
				v0[0] = S_old(i,I,k,V["u"]);
				v0[1] = S_old(i,I,k,V["v"]);
				Wi = Prim(S_old(i,I,k,V["rho"]), v0, S_old(i,I,k,V["p"]));
				//-------------------convert S to W---------------------------
				
				
				//--------------------------------------MUSCL half time step----------------------------------------------
				if(MUSCL){
				delta_l = delta(omega, Wll, Wl, Wi, scheme);
				delta_i = delta(omega, Wl, Wi, Wr, scheme);
				
				//WbarL:
				WbarL.rho = Wl.rho + 0.5*((1-dt/dy*Wl.u_vec[1])*delta_l.rho - dt/dy*Wl.rho*delta_l.u_vec[1]);
				WbarL.u_vec[0] = Wl.u_vec[0] + 0.5*((1-dt/dy*Wl.u_vec[1])*delta_l.u_vec[0]);
				WbarL.u_vec[1] = Wl.u_vec[1] + 0.5*((1-dt/dy*Wl.u_vec[1])*delta_l.u_vec[1] - dt/dy*((1/Wl.rho)*delta_l.p));
				WbarL.u_vec[2] = Wl.u_vec[2] + 0.5*((1-dt/dy*Wl.u_vec[1])*delta_l.u_vec[2]);
				S_old(i,L,k,V["EoS19"]) = EoS_properties(EoS, Wl.rho, Wl.p, 0.0, aL, "a");
				//aL = pow(gamma*(Wl.p)/Wl.rho, 0.5);
				WbarL.p = Wl.p + 0.5*((1-dt/dy*Wl.u_vec[1])*delta_l.p - dt/dy*(Wl.rho*pow(aL,2))*delta_l.u_vec[1]);					
				
				//WbarR:
				WbarR.rho = Wi.rho - 0.5*((1+dt/dy*Wi.u_vec[1])*delta_i.rho + dt/dy*Wi.rho*delta_i.u_vec[1]);
				WbarR.u_vec[0] = Wi.u_vec[0] - 0.5*((1+dt/dy*Wi.u_vec[1])*delta_i.u_vec[0]);
				WbarR.u_vec[1] = Wi.u_vec[1] - 0.5*((1+dt/dy*Wi.u_vec[1])*delta_i.u_vec[1] + dt/dy*((1/Wi.rho)*delta_i.p));
				WbarR.u_vec[2] = Wi.u_vec[2] - 0.5*((1+dt/dy*Wi.u_vec[1])*delta_i.u_vec[2]);
				S_old(i,I,k,V["EoS19"]) = EoS_properties(EoS, Wi.rho, Wi.p, 0.0, aR, "a");
				//aR = pow(gamma*(Wi.p)/Wi.rho, 0.5);
				WbarR.p = Wi.p - 0.5*((1+dt/dy*Wi.u_vec[1])*delta_i.p + dt/dy*(Wi.rho*pow(aL,2))*delta_i.u_vec[1]);
				
				WL = WbarL; WR = WbarR;
				}else{
				WL = Wl; WR = Wi;
				}
				//--------------------------------------MUSCL half time step----------------------------------------------
				//direct wave speed estimates:
				S_old(i,L,k,V["EoS19"]) = EoS_properties(EoS, WL.rho, WL.p, 0.0, aL, "a");
				//cout << "aL = " << aL << endl;
				S_old(i,I,k,V["EoS19"]) = EoS_properties(EoS, WR.rho, WR.p, 0.0, aR, "a");
			
				SL = ((WL.u_vec[1] - aL < WR.u_vec[1] - aR)? WL.u_vec[1] - aL : WR.u_vec[1] - aR);
				SR = ((WL.u_vec[1] + aL > WR.u_vec[1] + aR)? WL.u_vec[1] + aL : WR.u_vec[1] + aR);
				Sstar = (WR.u_vec[1]*WR.rho*(SR-WR.u_vec[1]) - WL.u_vec[1]*WL.rho*(SL-WL.u_vec[1]) + (WL.p-WR.p))\
						/(WR.rho*(SR-WR.u_vec[1])-WL.rho*(SL-WL.u_vec[1]));
				Splus = ((fabs(WL.u_vec[1])+aL > fabs(WR.u_vec[1])+aR)? \
						  fabs(WL.u_vec[1])+aL : fabs(WR.u_vec[1])+aR);
						  
						  
				//W-->U and W-->G conversions:
				UL = conWtoU(WL);
				UR = conWtoU(WR);
				GL = conWtoG(WL);
				GR = conWtoG(WR);
				
				//flux evaluation at (x/t)=0:	
				if( 0 < SL){
					Gi = GL;
				}else if(SL <= 0 && 0 < Sstar){
					//cout << "left star state \n";
					//left star state terms:
					ULstar.rho = WL.rho*((SL-WL.u_vec[1])/(SL-Sstar));
					ULstar.rhou_vec[0] = WL.rho*WL.u_vec[0]*((SL-WL.u_vec[1])/(SL-Sstar));
					ULstar.rhou_vec[1] = ULstar.rho*Sstar; 
					ULstar.E = ((SL-WL.u_vec[1])/(SL-Sstar))*(UL.E+WL.rho*(Sstar-WL.u_vec[1])*\
									(Sstar+WL.p/(WL.rho*(SL-WL.u_vec[1]))));
					if(ULstar.E < 0){
						cout << "negative energy computed in cell" << endl;
						break;
					}
					//left star state flux:
					GLstar.mass = GL.mass + SL*(ULstar.rho - UL.rho);
					GLstar.mom_vec = GL.mom_vec + SL*(ULstar.rhou_vec - UL.rhou_vec);  //check vector addition
					GLstar.en = GL.en + SL*(ULstar.E - UL.E);
					Gi = GLstar;
					
				}else if(Sstar <= 0 && 0 < SR){
					//cout << "right star state \n";
					//right star state terms:
					URstar.rho = WR.rho*((SR-WR.u_vec[1])/(SR-Sstar));
					URstar.rhou_vec[0] = WR.rho*WR.u_vec[0]*((SR-WR.u_vec[1])/(SR-Sstar));
					URstar.rhou_vec[1] = URstar.rho*Sstar; 
					URstar.E = ((SR-WR.u_vec[1])/(SR-Sstar))*(UR.E + WR.rho*(Sstar-WR.u_vec[1])*\
									(Sstar+WR.p/(WR.rho*(SR-WR.u_vec[1]))));
					if(ULstar.E < 0){
						cout << "negative energy computed in cell" << endl;
						break;
					}
					GRstar.mass = GR.mass + SR*(URstar.rho - UR.rho);
					GRstar.mom_vec = GR.mom_vec + SR*(URstar.rhou_vec - UR.rhou_vec); //check vector addition
					GRstar.en = GR.en + SR*(URstar.E - UR.E);
					Gi = GRstar;
					
				}else if(SR <= 0){
					Gi = GR; 
				}
				
				//-------------------convert G to flux_arr------------------------
				for(int n = 0; n<nc; ++n){
					G(i,j,k,n) = 0; //default set all fluxes to zero
				}
				//overwrite computed fluxes:
				G(i,j,k,V["mass"]) = Gi.mass;
				G(i,j,k,V["mom_x"]) = Gi.mom_vec[0];
				G(i,j,k,V["mom_y"]) = Gi.mom_vec[1];
				G(i,j,k,V["en"]) = Gi.en;
				//----------------------------------------------------------------		
				
				
			}
		}
	}//end i,j,k
	
}//end HLLC_y
			
	
void cons_update(const std::string sweep, Box const& box, Array4<Real> const& F, Array4<Real> const& S_old, 
        Array4<Real> const& S_new, Real const dt, Real const dx, int const nc, AccessVariable& V){
    
    const auto lo = lbound(box);
    const auto hi = ubound(box);
    
    Vector<Real> v0(3);
    Prim Wn(0,v0,0);
    Prim Wnew(0,v0,0);
    ConsU Un(0,v0,0);
    ConsU Unew(0,v0,0);
    
    double eps1, eps2, aR;
    double theta, alpha;
    int eos_temp;
    double tol = 1e-9;
    
    for(int k = lo.z; k <= hi.z; ++k) {
        for(int j = lo.y; j <= hi.y; ++j) {
            for(int i = lo.x; i <= hi.x; ++i) {
				
				//-------------------convert P to U---------------------------
				v0[0] = S_old(i,j,k,V["u"]);
				v0[1] = S_old(i,j,k,V["v"]);
				Wn = Prim(S_old(i,j,k,V["rho"]), v0, S_old(i,j,k,V["p"]));
				Un = conWtoU(Wn);
				
				//-------------conservative update on Un----------------------
				//where flux_prop_arr_x(rho*u, rho*u^2+p, rho*u*v, u(E+p)
                if(sweep == "x-sweep"){

					Unew.rho = Un.rho + dt/dx * (F(i,j,k,V["mass"]) - F(i+1,j,k,V["mass"])); 
					Unew.rhou_vec[0] = Un.rhou_vec[0] + dt/dx * (F(i,j,k,V["mom_x"]) - F(i+1,j,k,V["mom_x"])); 
					Unew.rhou_vec[1] = Un.rhou_vec[1] + dt/dx * (F(i,j,k,V["mom_y"]) - F(i+1,j,k,V["mom_y"])); 
					Unew.E = Un.E + dt/dx * (F(i,j,k,V["en"]) - F(i+1,j,k,V["en"]));
					
                    
                }else if(sweep == "y-sweep"){
                    Unew.rho = Un.rho + dt/dx * (F(i,j,k,V["mass"]) - F(i,j+1,k,V["mass"])); 
					Unew.rhou_vec[0] = Un.rhou_vec[0] + dt/dx * (F(i,j,k,V["mom_x"]) - F(i,j+1,k,V["mom_x"])); 
					Unew.rhou_vec[1] = Un.rhou_vec[1] + dt/dx * (F(i,j,k,V["mom_y"]) - F(i,j+1,k,V["mom_y"])); 
					Unew.E = Un.E + dt/dx * (F(i,j,k,V["en"]) - F(i,j+1,k,V["en"])); 
					
					/*
					//visualisation of shock delection method for HLLC-HS
					if(S_old(i,j,k,V["level_set"]) > 0){
						eps1 = 0;
					}else{
						eps1 = compute_eps_y(S_old, i, j, k, V);
					}
					
					
					if(S_old(i,j,k,V["level_set"]) > 0){
						eps2 = 0;
					}else{
						eps2 = compute_eps_x(S_old, i, j, k, V);
					}
					eos_temp = EoS_properties(EoS, S_old(i,j,k,V["rho"]), S_old(i,j,k,V["p"]), 0.0, aR, "a");
					
					theta = atan(eps2/(eps1+tol));
					alpha = cos(theta)*((eps1/aR)/(1 + eps1/aR));
					S_new(i,j,k,V["sigma"]) = alpha;
					
					theta = atan(eps1/(eps2+tol));
					alpha = cos(theta)*((eps2/aR)/(1 + eps2/aR));
					S_new(i,j,k,V["T"]) = alpha;
					*/
					
                }
				
				//--------------convert Unew to S_new-----------------------
				//Wnew = conUtoW(Unew);
				S_new(i,j,k,V["mass"]) = Unew.rho;
				S_new(i,j,k,V["mom_x"]) = Unew.rhou_vec[0];
				S_new(i,j,k,V["mom_y"]) = Unew.rhou_vec[1];
				S_new(i,j,k,V["en"]) = Unew.E;
								
            }
        }
    }//end i,j,k loop
    
    S_conUtoW(box, S_new, V);
    
}


void radial_sourceterm_integration(MultiFab& S_old, MultiFab& S_new, Real const dt, Real const dx, AccessVariable& V){
     
    //NOTE for half time step integration: assumes dt is passed as: dt/2
    
    Real R;
    
	for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
	{
		
		const Box& bx = mfi.tilebox();
		const auto lo = lbound(bx);
		const auto hi = ubound(bx);

		FArrayBox& stateout      =   S_new[mfi];
		FArrayBox& statein 		 = 	 S_old[mfi];
		
		Array4<Real> const& S_old_arr = statein.array();
		Array4<Real> const& S_new_arr = stateout.array();
		
		for(int k = lo.z; k <= hi.z; ++k) {
			for(int j = lo.y; j <= hi.y; ++j) {
				for(int i = lo.x; i <= hi.x; ++i) {
					
					//-------------------radial position---------------------------
					R = (i+0.5)*dx; //cell centre: should avoid div by R=0
					//-------------radial source term integration----------------------
					S_new_arr(i,j,k,V["mass"]) = S_old_arr(i,j,k,V["mass"]) + dt * (-1.0/R * S_old_arr(i,j,k,V["rho"]) * S_old_arr(i,j,k,V["u"]));
					S_new_arr(i,j,k,V["mom_x"]) = S_old_arr(i,j,k,V["mom_x"]) + dt * (-1.0/R * S_old_arr(i,j,k,V["rho"]) * pow(S_old_arr(i,j,k,V["u"]), 2));
					S_new_arr(i,j,k,V["en"]) = S_old_arr(i,j,k,V["en"]) + dt * (-1.0/R * S_old_arr(i,j,k,V["u"]) * (S_old_arr(i,j,k,V["en"]) + S_old_arr(i,j,k,V["p"])) );								
				}
			}
		}//end i,j,k loop
		
		//update primitive variables from new conservative variables: 
		S_conUtoW(bx, S_new_arr, V);	
	}    
}


Real max_dim_wave_speed(Box const& bx, Array4<Real> const& S_new, int const dim, AccessVariable& V){
	
	//convert from i dimension to n_index
	
	//this calculation should be sufficient to track EoS19 variable (?)
	
	string u_dim;
	if(dim == 0){
		u_dim = "u";
	}else if(dim==1){
		u_dim = "v";
	}
	
	
	Real max_wave_speed = 0.0;
	Real sound_speed, wave_speed;
	
	const auto lo = lbound(bx);
    const auto hi = ubound(bx);
    
    for(int k = lo.z; k <= hi.z; ++k) {
        for(int j = lo.y; j <= hi.y; ++j) {
            for(int i = lo.x; i <= hi.x; ++i) {
				
				S_new(i,j,k,V["EoS19"]) = EoS_properties(EoS, S_new(i,j,k,V["rho"]), S_new(i,j,k,V["p"]), 0.0, sound_speed, "a"); 
				wave_speed = sound_speed + fabs(S_new(i,j,k,V[u_dim]));
				
				if(wave_speed > max_wave_speed){
					max_wave_speed = wave_speed;
				}
				
			}
		}
	}//end i,j,k loop
	
	//cout << "max_wave_speed = " << max_wave_speed << endl;
	return max_wave_speed;
	
}



void rescale_fluxes(const std::string sweep, Box const& box, Array4<Real> const& F, 
					Real const dt, Real const dx, int const nc, AccessVariable& V){
	//nodal box is passed in 
	
	const auto lo = lbound(box);
    const auto hi = ubound(box);

	int ynode;
	int xnode;
	if(sweep == "x-sweep"){
		ynode = 0;
		xnode = 1;
	}else if(sweep == "y-sweep"){
		ynode = 1;
		xnode = 0;
	}
	
	//int const nc = 1;
	for(int n = 0; n<nc; ++n){
    for(int k = lo.z; k <= hi.z; ++k) {
        for(int j = lo.y; j <= hi.y + ynode; ++j) {
            for(int i = lo.x; i <= hi.x + xnode; ++i) {
				
                if(sweep == "x-sweep"){
                    F(i,j,k,n) = F(i,j,k,n) * (dt * dx); //dx = dy for Fx
                }else if(sweep == "y-sweep"){
                    F(i,j,k,n) = F(i,j,k,n) * (dt * dx); //dx = dx for Fy
                }
				

            }
        }
    }
	}
    
}

void determine_grad_max(const Box& box, Array4<Real> const& S_old_arr, const Real* dx, string prop, 
							string const refinement_condition, Real& grad_max, AccessVariable& V)
{
	
	const auto lo = lbound(box);
    const auto hi = ubound(box);
    
    Real local_grad, phi_x, phi_y;
    
    Real alpha = 1.0; 	//log scaling parameter
    
    //find max gradient
    for(int k = lo.z; k <= hi.z; ++k) {
		for(int j = lo.y; j <= hi.y; ++j) {
			for(int i = lo.x; i <= hi.x; ++i) {
				phi_x = fabs((S_old_arr(i+1,j,k,V[prop])-S_old_arr(i-1,j,k,V[prop]))/(2*dx[0]));
				phi_y = fabs((S_old_arr(i,j+1,k,V[prop])-S_old_arr(i,j-1,k,V[prop]))/(2*dx[1]));
				if(refinement_condition == "grad"){
					local_grad = sqrt(pow(phi_x,2) + pow(phi_y,2));
					//local_grad = (phi_x > phi_y ? phi_x : phi_y);
					if(prop == "mass"){
						local_grad = local_grad/fabs(S_old_arr(i,j,k,V[prop])); //scale by property magnitude 
					}															//for rho only as must be non-zero everywhere
				}else if(refinement_condition == "log_grad"){
					local_grad = (phi_x > phi_y ? phi_x : phi_y); 
					if(prop == "mass"){
						local_grad = local_grad/fabs(S_old_arr(i,j,k,V[prop])); //scale by property magnitude 
					}															//for rho only as must be non-zero everywhere
					local_grad = log(alpha*local_grad + 1);
				}else{
					amrex::Abort("refinement_condition invalid - check inputs file");
				}
				if(local_grad > grad_max){
					grad_max = local_grad;
				}
			}
		}
	}//end i,j,k
	
	//cout << "grad_max = " << grad_max << endl;
	
}
	

void amr_tagging(Array4<char> const& tagarr, const Box& box, Array4<Real> const& S_old_arr, const Real* dx, Real grad_max, Real grad_frac,
					string prop, string const refinement_condition, AccessVariable& V)
{
	//tagging function for refinement
	//based on gradient of rho
	//assumes rho is ncomp element 0
	//box has S_new dimensions, but S_old_arr has at least +1 ngrow
	//therefore, we can use central different function for gradient
	
	const auto lo = lbound(box);
    const auto hi = ubound(box);
    
    Real local_grad, phi_x, phi_y;
    //Real prop_max = 0.5;
        
    Real alpha = 1.0; 	//log scaling parameter
    
    
    
	//set tags based on gradient relative to max gradient: 
	//ie. prop_max is proportion of max gradient
	for(int k = lo.z; k <= hi.z; ++k) {
		for(int j = lo.y; j <= hi.y; ++j) {
			for(int i = lo.x; i <= hi.x; ++i) {
                if(S_old_arr(i,j,k,V["EB"]) < -0.99){
                    tagarr(i,j,k,V[prop]) = TagBox::CLEAR;   //never refine on rigid body cells
                }else{
                    phi_x = fabs((S_old_arr(i+1,j,k,V[prop])-S_old_arr(i-1,j,k,V[prop]))/(2*dx[0]));
                    phi_y = fabs((S_old_arr(i,j+1,k,V[prop])-S_old_arr(i,j-1,k,V[prop]))/(2*dx[1]));
                    if(refinement_condition == "grad"){
                        local_grad = sqrt(pow(phi_x,2) + pow(phi_y,2));
                        //local_grad = (phi_x > phi_y ? phi_x : phi_y);
                        if(prop == "mass"){
                            local_grad = local_grad/fabs(S_old_arr(i,j,k,V[prop])); //scale by property magnitude 
                        }															//for rho only as must be non-zero everywhere
                    }else if(refinement_condition == "log_grad"){
                        local_grad = sqrt(pow(phi_x,2) + pow(phi_y,2));
                        //local_grad = (phi_x > phi_y ? phi_x : phi_y); 
                        if(prop == "mass"){
                            local_grad = local_grad/fabs(S_old_arr(i,j,k,V[prop])); //scale by property magnitude 
                        }															//for rho only as must be non-zero everywhere
                        local_grad = log(alpha*local_grad + 1);
                    }else{
                        amrex::Abort("refinement_condition invalid - check inputs file");
                    }
                    /*
                    if(prop == "mass"){
                        S_old_arr(i,j,k,V["grad_rho"]) = local_grad;
                    }
                    */
                    if(local_grad > grad_frac*grad_max){
                        tagarr(i,j,k,V[prop]) = TagBox::SET;
                        //cout << "tagged" << endl;
					}else{
                        tagarr(i,j,k,V[prop]) = TagBox::CLEAR;
                    }
                    
                    
                    if((S_old_arr(i,j,k,V["EB"]) > -1 && S_old_arr(i,j,k,V["EB"]) < 1)){
						//BE CAREFUL OF POTENTIAL SEGFAULT HERE
						tagarr(i,j,k,V[prop]) = TagBox::SET;   //additionally refine on the boundary cells
					}
					
                }
			}
		}
	}//end i,j,k
		
}

void Qtransform(Array4<Real> const& S_old_arr, AccessVariable& V, const int& i, const int& j, const int& k, const int& p0_i, const int& p0_j, 
						const double& D00, const double& D10, const double& D01, const double& D11){
		
		double vx, vy;
		double nx, ny; 
							
		S_old_arr(i,j,k,V["rho"]) = (D00*S_old_arr(p0_i,p0_j,k,V["rho"]) + D10*S_old_arr(p0_i+1,p0_j,k,V["rho"]) \
										+ D01*S_old_arr(p0_i,p0_j+1,k,V["rho"]) + D11*S_old_arr(p0_i+1,p0_j+1,k,V["rho"])) / \
										(D00 + D10 + D01 + D11);
		S_old_arr(i,j,k,V["p"]) = (D00*S_old_arr(p0_i,p0_j,k,V["p"]) + D10*S_old_arr(p0_i+1,p0_j,k,V["p"]) \
										+ D01*S_old_arr(p0_i,p0_j+1,k,V["p"]) + D11*S_old_arr(p0_i+1,p0_j+1,k,V["p"])) / \
										(D00 + D10 + D01 + D11);
		vx = (D00*S_old_arr(p0_i,p0_j,k,V["u"]) + D10*S_old_arr(p0_i+1,p0_j,k,V["u"]) \
										+ D01*S_old_arr(p0_i,p0_j+1,k,V["u"]) + D11*S_old_arr(p0_i+1,p0_j+1,k,V["u"])) / \
										(D00 + D10 + D01 + D11);
		vy = (D00*S_old_arr(p0_i,p0_j,k,V["v"]) + D10*S_old_arr(p0_i+1,p0_j,k,V["v"]) \
										+ D01*S_old_arr(p0_i,p0_j+1,k,V["v"]) + D11*S_old_arr(p0_i+1,p0_j+1,k,V["v"])) / \
										(D00 + D10 + D01 + D11);
		nx = (D00*S_old_arr(p0_i,p0_j,k,V["nx"]) + D10*S_old_arr(p0_i+1,p0_j,k,V["nx"]) \
										+ D01*S_old_arr(p0_i,p0_j+1,k,V["nx"]) + D11*S_old_arr(p0_i+1,p0_j+1,k,V["nx"])) / \
										(D00 + D10 + D01 + D11);
		ny = (D00*S_old_arr(p0_i,p0_j,k,V["ny"]) + D10*S_old_arr(p0_i+1,p0_j,k,V["ny"]) \
										+ D01*S_old_arr(p0_i,p0_j+1,k,V["ny"]) + D11*S_old_arr(p0_i+1,p0_j+1,k,V["ny"])) / \
										(D00 + D10 + D01 + D11);
		 //check which nx, ny should be used here
		 //will rigid body point normals be the same as projected point? 
		S_old_arr(i,j,k,V["u"]) = vx - 2*(nx*vx + ny*vy)*nx;
		S_old_arr(i,j,k,V["v"]) = vy - 2*(nx*vx + ny*vy)*ny;
	
}

void Qtransform2(Array4<Real> const& S_old_arr, AccessVariable& V, const int& i, const int& j, const int& k, const int& p0_i, const int& p0_j, 
						const double& D00, const double& D10, const double& D01, const double& D11){
		
	
		S_old_arr(i,j,k,V["rho"]) = (D00*S_old_arr(p0_i,p0_j,k,V["rho"]) + D10*S_old_arr(p0_i+1,p0_j,k,V["rho"]) \
										+ D01*S_old_arr(p0_i,p0_j+1,k,V["rho"]) + D11*S_old_arr(p0_i+1,p0_j+1,k,V["rho"])) / \
										(D00 + D10 + D01 + D11);
		S_old_arr(i,j,k,V["p"]) = (D00*S_old_arr(p0_i,p0_j,k,V["p"]) + D10*S_old_arr(p0_i+1,p0_j,k,V["p"]) \
										+ D01*S_old_arr(p0_i,p0_j+1,k,V["p"]) + D11*S_old_arr(p0_i+1,p0_j+1,k,V["p"])) / \
										(D00 + D10 + D01 + D11);
		S_old_arr(i,j,k,V["u"]) = (D00*S_old_arr(p0_i,p0_j,k,V["u"]) + D10*S_old_arr(p0_i+1,p0_j,k,V["u"]) \
										+ D01*S_old_arr(p0_i,p0_j+1,k,V["u"]) + D11*S_old_arr(p0_i+1,p0_j+1,k,V["u"])) / \
										(D00 + D10 + D01 + D11);
		S_old_arr(i,j,k,V["v"]) = (D00*S_old_arr(p0_i,p0_j,k,V["v"]) + D10*S_old_arr(p0_i+1,p0_j,k,V["v"]) \
										+ D01*S_old_arr(p0_i,p0_j+1,k,V["v"]) + D11*S_old_arr(p0_i+1,p0_j+1,k,V["v"])) / \
										(D00 + D10 + D01 + D11);
	
}
					

void initialise_boundary_cells(const Box& box, Array4<Real> const& S_old_arr, Array4<const EBCellFlag>& flags_arr, const int& nc, 
									AccessVariable& V, const Real* dx){
	
	const auto lo = lbound(box);
    const auto hi = ubound(box);
    
    double ip, jp;
    int p0_i, p0_j;
    double delta_i, delta_j;
    
    double DX = pow(pow(dx[0],2) + pow(dx[1],2),0.5);
    double D00, D10, D01, D11;
    double tol = 1e-6;
    
    vector<int> i_list;
    vector<int> j_list;
    
    //find cut cells, covered cells and regular cells:
    for(int k = lo.z; k <= hi.z; ++k) {
		for(int j = lo.y; j <= hi.y; ++j) {
			for(int i = lo.x; i <= hi.x; ++i) {
				
				if(flags_arr(i,j,k).isSingleValued()){
					//cut cell
					//S_old_arr(i,j,k, V["EB"]) = 0;
					//projection point into fluid region
					ip = i + (DX/dx[0])*S_old_arr(i,j,k, V["nx"]); //x-position on integer scale
					jp = j + (DX/dx[1])*S_old_arr(i,j,k, V["ny"]); //y-position on integer scale
					p0_i = (int)(floor(ip));	//bottom left x-coord
					p0_j = (int)(floor(jp));	//bottom left y-coord
					if(p0_i < 0 || p0_j < 0){
						amrex::Abort("error in boundary cells function 1");	//make sure it doesn't exceed domain
					}
					//haven't done top-right domain bounds check
					delta_i = ip-p0_i;	//fraction of cell size into 4-point mixture zone
					delta_j = jp-p0_j;
					if(delta_i + delta_j < tol){ 
						//projected point hits the centre of another point
						//cout << "too small" << endl;
						D00 = 1.0; D10 = 0.0; D01 = 0.0; D11 = 0.0;
						if(!flags_arr(p0_i,p0_j,k).isRegular()){
							D00 = 0;
							i_list.push_back(i);
							j_list.push_back(j);
							//amrex::Abort("in boundary cells calc, projection lands exactly on a non-fluid cell");
						}
					}else{
								
						if(flags_arr(p0_i,p0_j,k).isRegular()){ //only include in mixture rule if cell is fluid (regular)
							D00 = pow(pow(delta_i*dx[0],2) + pow(delta_j*dx[1],2),-0.5);
						}else{
							D00 = 0;
						}
						if(flags_arr(p0_i+1,p0_j,k).isRegular()){
							D10 = pow(pow((1-delta_i)*dx[0],2) + pow(delta_j*dx[1],2),-0.5);
						}else{
							D10 = 0;
						}
						if(flags_arr(p0_i,p0_j+1,k).isRegular()){
							D01 = pow(pow(delta_i*dx[0],2) + pow((1-delta_j)*dx[1],2),-0.5);
						}else{
							D01 = 0;
						}
						if(flags_arr(p0_i+1,p0_j+1,k).isRegular()){
							D11 = pow(pow((1-delta_i)*dx[0],2) + pow((1-delta_j)*dx[1],2),-0.5);
						}else{
							D11 = 0;
						}
						
					}
					
					//maxD = pow(pow(tol + tol,-0.5);
					if(D00 + D10 + D01 + D11 == 0){
						i_list.push_back(i);
						j_list.push_back(j);
						//amrex::Abort("in boundary cells calc, projection does not reach fluid region");
					}else{
						Qtransform(S_old_arr, V, i, j, k, p0_i, p0_j, D00, D10, D01, D11);
					}
					
				}else if(flags_arr(i,j,k).isCovered()){
					//fully solid
					if(fabs(S_old_arr(i,j,k, V["nx"])) + fabs(S_old_arr(i,j,k, V["ny"])) != 0){ 
						//otherwise p will never be overwritten
						S_old_arr(i,j,k, V["p"]) = -1e5; 
					}
				}else if(flags_arr(i,j,k).isRegular()){
					//fully fluid
				}

			}
		}
	}//end i,j,k loop
	
	//handle special cell cases
	//for cells which don't reach valid fluid region
	//once all other boundary cells are filled
	//allow the calc to be a mixture of the valid boundary cells 
	for(int s = 0; s<i_list.size(); ++s){
		int i = i_list[s];
		int j = j_list[s];
		int k = 0;
		ip = i + (DX/dx[0])*S_old_arr(i,j,k, V["nx"]); //x-position on integer scale
		jp = j + (DX/dx[1])*S_old_arr(i,j,k, V["ny"]); //y-position on integer scale
		p0_i = (int)(floor(ip));	//bottom left x-coord
		p0_j = (int)(floor(jp));	//bottom left y-coord
		if(p0_i < 0 || p0_j < 0){
			amrex::Abort("error in boundary cells function 1");	//make sure it doesn't exceed domain
		}
		//haven't done top-right domain bounds check
		delta_i = ip-p0_i;	//fraction of cell size into 4-point mixture zone
		delta_j = jp-p0_j;
		if(delta_i + delta_j < tol){ 
			//projected point hits the centre of another point
			//cout << "too small" << endl;
			D00 = 1.0; D10 = 0.0; D01 = 0.0; D11 = 0.0;
		}else{	
			if(!flags_arr(p0_i,p0_j,k).isCovered()){ //only include in mixture rule if cell is fluid (regular)
				D00 = pow(pow(delta_i*dx[0],2) + pow(delta_j*dx[1],2),-0.5);
			}else{
				D00 = 0;
			}
			if(!flags_arr(p0_i+1,p0_j,k).isCovered()){
				D10 = pow(pow((1-delta_i)*dx[0],2) + pow(delta_j*dx[1],2),-0.5);
			}else{
				D10 = 0;
			}
			if(!flags_arr(p0_i,p0_j+1,k).isCovered()){
				D01 = pow(pow(delta_i*dx[0],2) + pow((1-delta_j)*dx[1],2),-0.5);
			}else{
				D01 = 0;
			}
			if(!flags_arr(p0_i+1,p0_j+1,k).isCovered()){
				D11 = pow(pow((1-delta_i)*dx[0],2) + pow((1-delta_j)*dx[1],2),-0.5);
			}else{
				D11 = 0;
			}
			
		}
		
		//maxD = pow(pow(tol + tol,-0.5);
		if(D00 + D10 + D01 + D11 == 0){
			amrex::Abort("special case projection goes entirely into solid region");
		}else{
			Qtransform2(S_old_arr, V, i, j, k, p0_i, p0_j, D00, D10, D01, D11);
		}
	}
		
	
	
	S_conWtoU(box, S_old_arr, V);
	
}

void Qcalc(Array4<Real> const& S_arr, const int& i, const int& j, const int& k, AccessVariable& V, const Real* dx){
	
	double qi, qj;
	
	double test_val;
	
	if(S_arr(i,j,k, V["nx"]) >= 0){
		qi = i+1;
	}else{
		qi = i-1;
		if(qi < 0){
			qi = 0;
			//cout << "FLAG" << endl;
		}
	}
	if(S_arr(i,j,k, V["ny"]) >= 0){
		qj = j+1;
	}else{
		qj = j-1;
	}
	
	if(fabs(S_arr(i,j,k, V["nx"])) + fabs(S_arr(i,j,k, V["ny"])) == 0){
		//S_arr(i,j,k) remains as at previous
	}else{
		S_arr(i,j,k, V["rho"]) = (fabs(S_arr(i,j,k, V["nx"]))/dx[0] * S_arr(qi,j,k, V["rho"]) + \
										fabs(S_arr(i,j,k, V["ny"]))/dx[1] * S_arr(i,qj,k, V["rho"]))/ \
											(fabs(S_arr(i,j,k, V["nx"]))/dx[0] + fabs(S_arr(i,j,k, V["ny"]))/dx[1]);
		S_arr(i,j,k, V["p"]) = (fabs(S_arr(i,j,k, V["nx"]))/dx[0] * S_arr(qi,j,k, V["p"]) + \
										fabs(S_arr(i,j,k, V["ny"]))/dx[1] * S_arr(i,qj,k, V["p"]))/ \
											(fabs(S_arr(i,j,k, V["nx"]))/dx[0] + fabs(S_arr(i,j,k, V["ny"]))/dx[1]);	
		S_arr(i,j,k, V["u"]) = (fabs(S_arr(i,j,k, V["nx"]))/dx[0] * S_arr(qi,j,k, V["u"]) + \
										fabs(S_arr(i,j,k, V["ny"]))/dx[1] * S_arr(i,qj,k, V["u"]))/ \
											(fabs(S_arr(i,j,k, V["nx"]))/dx[0] + fabs(S_arr(i,j,k, V["ny"]))/dx[1]);									
		S_arr(i,j,k, V["v"]) = (fabs(S_arr(i,j,k, V["nx"]))/dx[0] * S_arr(qi,j,k, V["v"]) + \
										fabs(S_arr(i,j,k, V["ny"]))/dx[1] * S_arr(i,qj,k, V["v"]))/ \
											(fabs(S_arr(i,j,k, V["nx"]))/dx[0] + fabs(S_arr(i,j,k, V["ny"]))/dx[1]);
	}
	
	
}

void sweep(const Box& box, Array4<Real> const& S_arr, Array4<const EBCellFlag>& flags_arr, const int& nc, 
									AccessVariable& V, const Real* dx, double& global_p_min){
	
	//boxes passed are either cut or regular (solid)
	//box in MFI is S_new and S_arr passed is S_old (with NGROW)
	//so finite diff stencils shoul not exceed bounds
	
	const auto lo = lbound(box);
    const auto hi = ubound(box);
    
    //conduct 4 sweeps and operate on cells which are regular only
    
    //sweep 1: y /\  x --> 
    for(int k = lo.z; k <= hi.z; ++k) {
		for(int j = lo.y; j <= hi.y; ++j) {
			for(int i = lo.x; i <= hi.x; ++i) {
				
				if(flags_arr(i,j,k).isCovered()){
					//rigid body cell
					Qcalc(S_arr, i,j,k, V, dx);
				}

			}
		}
	}//end i,j,k loop
	
	//sweep 2: y \/  x-->
    for(int k = lo.z; k <= hi.z; ++k) {
		for(int j = hi.y; j >= lo.y; --j) {
			for(int i = lo.x; i <= hi.x; ++i) {
				
				if(flags_arr(i,j,k).isCovered()){
					//rigid body cell
					Qcalc(S_arr, i,j,k, V, dx);
				}

			}
		}
	}//end i,j,k loop
	
	//sweep 3:  y /\  x <-- 
    for(int k = lo.z; k <= hi.z; ++k) {
		for(int j = lo.y; j <= hi.y; ++j) {
			for(int i = hi.x; i >= lo.x; --i) {
				
				if(flags_arr(i,j,k).isCovered()){
					//rigid body cell
					Qcalc(S_arr, i,j,k, V, dx);
				}

			}
		}
	}//end i,j,k loop
	
	//sweep 4:  y \/  x <-- 
    for(int k = lo.z; k <= hi.z; ++k) {
		for(int j = hi.y; j >= lo.y; --j) {
			for(int i = hi.x; i >= lo.x; --i) {
				
				if(flags_arr(i,j,k).isCovered()){
					//rigid body cell
					Qcalc(S_arr, i,j,k, V, dx);
				}

                if(S_arr(i,j,k, V["p"]) < global_p_min){
                    global_p_min = S_arr(i,j,k, V["p"]);
                }

			}
		}
	}//end i,j,k loop
	
	//S_conWtoU(box, S_arr, V);
	
}


