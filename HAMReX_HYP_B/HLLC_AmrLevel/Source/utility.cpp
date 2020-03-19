#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EBAmrUtil.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EB_levelset.H>
#include <AMReX_EB_LSCore.H>
#include <AMReX_EB_LSCoreBase.H>

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <assert.h>
#include <functional>
#include <algorithm>
#include <iomanip>
#include <limits>

#include "structdefs.H"
#include "classdefs.H"
#include "funcdefs.H"

using namespace amrex;
using namespace std;


//typedef vector<double> Vector;

//template<typename T>
void printVector(Vector<Real> vec){
    int len = vec.size();
    for(int i; i<len; ++i){
        cout << vec[i] << ' ';
    }
    cout << endl;
}//Vector-array print function 

//operator overloading for addition, subtraction and multiplication of vectors

//template<typename T>
Vector<Real> operator+(Vector<Real> const& first, Vector<Real> const& second)
{
    assert(first.size() == second.size());

    Vector<Real> result;
    result.reserve(first.size());

    transform(first.begin(), first.end(), second.begin(), back_inserter(result), plus<Real>());
    return result;
}

//template<typename T>
Vector<Real> operator-(Vector<Real> const& first, Vector<Real> const& second)
{
    assert(first.size() == second.size());

    Vector<Real> result;
    result.reserve(first.size());

    transform(first.begin(), first.end(), second.begin(), back_inserter(result), minus<Real>());
    return result;
}

//template<typename T>
Vector<Real> operator*(double const& scal, Vector<Real> const& first)
{
    
    vector<double> second(first.size(), scal);
    assert(first.size() == second.size());

    Vector<Real> result;
    result.reserve(first.size());

    transform(first.begin(), first.end(), second.begin(), back_inserter(result), multiplies<Real>());
    return result;
}



//vector functions:
double mag(Vector<Real> const& v){
    assert(v.size() == 3);
    double mag_v = pow( pow(v[0],2) + pow(v[1],2) + pow(v[2],2), 0.5);
    return mag_v;
}

double dot(Vector<Real> const& v1, Vector<Real> const& v2){
    assert(v1.size() == 3 && v2.size() == 3);
    double dot_v = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
    return dot_v;
}

Vector<Real> cross(Vector<Real> const& v1, Vector<Real> const& v2){
    assert(v1.size() == 3 && v2.size() == 3);
    Vector<Real> out_vec(3,0);
    out_vec[0] = v1[1]*v2[2] - v1[2]*v2[1];
    out_vec[1] = -(v1[0]*v2[2] - v1[2]*v2[0]);
    out_vec[2] = v1[0]*v2[1] - v1[1]*v2[0];
    return out_vec;
}

//accessvariable wrapper:

AccessVariable::AccessVariable(Vector<string> const &problem_variables){
	
	int num_vars = problem_variables.size();
	
	for(int n =0; n<num_vars; ++n){
		variables.insert(std::pair<string,int>(problem_variables[n],n));
	}
	
}

int& AccessVariable::operator[](string const &var)
{
	return variables[var];
}


//Embedded boundary stuff:


void makeEmbeddedBoundary(Geometry& geom, const int& max_level){
	
	ParmParse pp;
	string EBgeom;
	pp.get("EBgeom", EBgeom);
	ParmParse ppEB(EBgeom);
	int ls_pad = 4;		//number of cells away from EB in which level set is computed
	if(EBgeom == "sphere"){
		
		Real radius;
		Vector<Real> EB_centre(AMREX_SPACEDIM);
		ppEB.get("R", radius);
		ppEB.getarr("centre", EB_centre);
		Array<Real,AMREX_SPACEDIM> centre{0.0, 0.0}; //default Center of the sphere
		for(int i = 0; i<AMREX_SPACEDIM; ++i){
			centre[i] = EB_centre[i];
		}
		
		bool inside = false;  // Is the fluid inside the sphere --> false
		EB2::SphereIF sphere(radius, centre, inside);

		auto shop = EB2::makeShop(sphere);

		EB2::Build(shop, geom, max_level, max_level, ls_pad); 
		//arguments: (gshop, geom, required_coarsening_level, max_coarsening_level, ngrow)
	
	}else if(EBgeom == "plane"){
		Vector<Real> point, normal;
		ppEB.getarr("point", point);
		ppEB.getarr("normal", normal);
		
		EB2::PlaneIF plane({AMREX_D_DECL(point[0], point[1], 0.)},
                             {AMREX_D_DECL(normal[0], normal[1], 0.)}, false);
        
        auto shop = EB2::makeShop(plane);

		EB2::Build(shop, geom, max_level, max_level, ls_pad); 
		//arguments: (gshop, geom, required_coarsening_level, max_coarsening_level, ngrow)                     
                             
		
	}else if(EBgeom == "planeX"){
		Vector<Real> point, normal;
		Vector<Real> point2, normal2;
		ppEB.getarr("point", point);
		ppEB.getarr("normal", normal);
		ppEB.getarr("point2", point2);
		ppEB.getarr("normal2", normal2);
		
		EB2::PlaneIF plane1({AMREX_D_DECL(point[0], point[1], 0.)},
                             {AMREX_D_DECL(normal[0], normal[1], 0.)}, false);
        EB2::PlaneIF plane2({AMREX_D_DECL(point2[0], point2[1], 0.)},
                             {AMREX_D_DECL(normal2[0], normal2[1], 0.)}, false);
        
        auto planeX = EB2::makeIntersection(plane1, plane2);
        
        auto shop = EB2::makeShop(planeX);

		EB2::Build(shop, geom, max_level, max_level, ls_pad); 
		//arguments: (gshop, geom, required_coarsening_level, max_coarsening_level, ngrow)                     
                             
		
	}else if(EBgeom == "triangle"){
		Vector<Real> point, normal;
		Vector<Real> point2, normal2;
		Vector<Real> point3, normal3;
		ppEB.getarr("point", point);
		ppEB.getarr("normal", normal);
		ppEB.getarr("point2", point2);
		ppEB.getarr("normal2", normal2);
		ppEB.getarr("point3", point3);
		ppEB.getarr("normal3", normal3);
		
		EB2::PlaneIF plane1({AMREX_D_DECL(point[0], point[1], 0.)},
                             {AMREX_D_DECL(normal[0], normal[1], 0.)}, false);
        EB2::PlaneIF plane2({AMREX_D_DECL(point2[0], point2[1], 0.)},
                             {AMREX_D_DECL(normal2[0], normal2[1], 0.)}, false);
        EB2::PlaneIF plane3({AMREX_D_DECL(point3[0], point3[1], 0.)},
                             {AMREX_D_DECL(normal3[0], normal3[1], 0.)}, false);
        
        auto planeX = EB2::makeIntersection(plane1, plane2, plane3);
        
        auto shop = EB2::makeShop(planeX);

		EB2::Build(shop, geom, max_level, max_level, ls_pad); 
		//arguments: (gshop, geom, required_coarsening_level, max_coarsening_level, ngrow)                     
                             
		
	}else if(EBgeom == "none"){
		//no additional settings
		cout << "EBgeom is set to \"none\" \n";
		
		//EB2 build still required for valid intitialisation under EB = TRUE in makefile
		//make a random shop to satisfy EB = TRUE
		Real radius = 1.0;
		Array<Real,AMREX_SPACEDIM> centre{0.0, 0.0}; //default Center of the sphere
		bool inside = false;  // Is the fluid inside the sphere --> false
		EB2::SphereIF sphere(radius, centre, inside);

		auto shop = EB2::makeShop(sphere);
		
		EB2::Build(shop, geom, max_level, max_level, ls_pad);
		
		
	}else{
		amrex::Abort("invalid EBgeom passed from inputs file");
	}
	
}

void initialise_levelset_geometry(MultiFab& S_new, const int level, Geometry& geom, SettingsStruct const& sim, 
										ParameterStruct const& p, AccessVariable& V)
{
	int ls_pad = 4; //number of cells away from EB in which level set is computed
	int eb_pad = 2; //sim.NGROW;
	int ls_ref = 1; //1
	int eb_ref = 2; //2;
	int ebt_size = 32;
	
	const BoxArray ba = S_new.boxArray();
	const DistributionMapping dm = S_new.DistributionMap();
	
	if(sim.EBgeom == "sphere"){
		
		Real radius = p.EB_R;
		Array<Real,AMREX_SPACEDIM> centre{p.EB_centre[0], p.EB_centre[1]}; 
		
		//Define EB:
		EB2::SphereIF sphere(radius, centre, false); // Is the fluid inside the sphere --> false
        EB2::GeometryShop<EB2::SphereIF> sphere_gshop(sphere);
        
        // Build level-set factory
		LSFactory level_set(level, ls_ref, eb_ref, ls_pad, eb_pad, ba, geom, dm);
		
		// Build EB
		const Geometry& eb_geom = level_set.get_eb_geom();
		EB2::Build(sphere_gshop, eb_geom, p.max_level, p.max_level);
		//parameters: gshop, eb_geom, required_coarsening_level, max_coarsening_level
		
		const EB2::IndexSpace& sphere_ebis = EB2::IndexSpace::top();
		const EB2::Level&      sphere_lev  = sphere_ebis.getLevel(eb_geom);

		// Build EB factory
		EBFArrayBoxFactory eb_factory(sphere_lev, eb_geom, ba, dm, {eb_pad, eb_pad, eb_pad}, EBSupport::full); 
		
		// Fill level-set (factory)
		GShopLSFactory<EB2::SphereIF> sphere_lsgs(sphere_gshop, level_set);
		std::unique_ptr<MultiFab> sphere_mf_impfunc = sphere_lsgs.fill_impfunc();

		MultiFab::Copy(S_new, *sphere_mf_impfunc, 0, V["level_set"], 1, S_new.nGrow());
		
	}else if(sim.EBgeom == "plane"){
		
		//Define EB:
		EB2::PlaneIF plane({AMREX_D_DECL(p.EB_point[0], p.EB_point[1], 0.)},
                             {AMREX_D_DECL(p.EB_normal[0], p.EB_normal[1], 0.)}, false);
        EB2::GeometryShop<EB2::PlaneIF> plane_gshop(plane);
        
        // Build level-set factory
		LSFactory level_set(level, ls_ref, eb_ref, ls_pad, eb_pad, ba, geom, dm);
		
		// Build EB
		const Geometry& eb_geom = level_set.get_eb_geom();
		EB2::Build(plane_gshop, eb_geom, p.max_level, p.max_level);
		//parameters: gshop, eb_geom, required_coarsening_level, max_coarsening_level
		
		const EB2::IndexSpace& plane_ebis = EB2::IndexSpace::top();
		const EB2::Level&      plane_lev  = plane_ebis.getLevel(eb_geom);

		// Build EB factory
		EBFArrayBoxFactory eb_factory(plane_lev, eb_geom, ba, dm, {eb_pad, eb_pad, eb_pad}, EBSupport::full);
		
		// Fill level-set (factory)
		GShopLSFactory<EB2::PlaneIF> plane_lsgs(plane_gshop, level_set);
		std::unique_ptr<MultiFab> plane_mf_impfunc = plane_lsgs.fill_impfunc();

		MultiFab::Copy(S_new, *plane_mf_impfunc, 0, V["level_set"], 1, S_new.nGrow());
		
	}else if(sim.EBgeom == "planeX"){
		
		//Define EB:
		EB2::PlaneIF plane1({AMREX_D_DECL(p.EB_point[0], p.EB_point[1], 0.)},
                             {AMREX_D_DECL(p.EB_normal[0], p.EB_normal[1], 0.)}, false);
        
        //EB2::GeometryShop<EB2::PlaneIF> plane_gshop(plane1);
        
        EB2::PlaneIF plane2({AMREX_D_DECL(p.EB_point2[0], p.EB_point2[1], 0.)},
                             {AMREX_D_DECL(p.EB_normal2[0], p.EB_normal2[1], 0.)}, false);
        
        auto planeX = EB2::makeIntersection(plane1, plane2);
        
        auto planeX_gshop = EB2::makeShop(planeX);
        
        // Build level-set factory
		LSFactory level_set(level, ls_ref, eb_ref, ls_pad, eb_pad, ba, geom, dm);
		
		// Build EB
		const Geometry& eb_geom = level_set.get_eb_geom();
		EB2::Build(planeX_gshop, eb_geom, p.max_level, p.max_level);
		//parameters: gshop, eb_geom, required_coarsening_level, max_coarsening_level
		
		const EB2::IndexSpace& plane_ebis = EB2::IndexSpace::top();
		const EB2::Level&      plane_lev  = plane_ebis.getLevel(eb_geom);

		// Build EB factory
		EBFArrayBoxFactory eb_factory(plane_lev, eb_geom, ba, dm, {eb_pad, eb_pad, eb_pad}, EBSupport::full);
		
		// Fill level-set (factory)
		GShopLSFactory< EB2::IntersectionIF<EB2::PlaneIF, EB2::PlaneIF> > plane_lsgs(planeX_gshop, level_set);
		std::unique_ptr<MultiFab> plane_mf_impfunc = plane_lsgs.fill_impfunc();

		MultiFab::Copy(S_new, *plane_mf_impfunc, 0, V["level_set"], 1, S_new.nGrow());
		
	}else if(sim.EBgeom == "triangle"){
		
		//Define EB:
		EB2::PlaneIF plane1({AMREX_D_DECL(p.EB_point[0], p.EB_point[1], 0.)},
                             {AMREX_D_DECL(p.EB_normal[0], p.EB_normal[1], 0.)}, false);

        EB2::PlaneIF plane2({AMREX_D_DECL(p.EB_point2[0], p.EB_point2[1], 0.)},
                             {AMREX_D_DECL(p.EB_normal2[0], p.EB_normal2[1], 0.)}, false);
                             
        EB2::PlaneIF plane3({AMREX_D_DECL(p.EB_point3[0], p.EB_point3[1], 0.)},
                             {AMREX_D_DECL(p.EB_normal3[0], p.EB_normal3[1], 0.)}, false);
        
        auto planeX = EB2::makeIntersection(plane1, plane2, plane3);
        
        auto planeX_gshop = EB2::makeShop(planeX);
        
        // Build level-set factory
		LSFactory level_set(level, ls_ref, eb_ref, ls_pad, eb_pad, ba, geom, dm);
		
		// Build EB
		const Geometry& eb_geom = level_set.get_eb_geom();
		EB2::Build(planeX_gshop, eb_geom, p.max_level, p.max_level);
		//parameters: gshop, eb_geom, required_coarsening_level, max_coarsening_level
		
		const EB2::IndexSpace& plane_ebis = EB2::IndexSpace::top();
		const EB2::Level&      plane_lev  = plane_ebis.getLevel(eb_geom);

		// Build EB factory
		EBFArrayBoxFactory eb_factory(plane_lev, eb_geom, ba, dm, {eb_pad, eb_pad, eb_pad}, EBSupport::full);
		
		// Fill level-set (factory)
		GShopLSFactory< EB2::IntersectionIF<EB2::PlaneIF, EB2::PlaneIF, EB2::PlaneIF> > plane_lsgs(planeX_gshop, level_set);
		std::unique_ptr<MultiFab> plane_mf_impfunc = plane_lsgs.fill_impfunc();

		MultiFab::Copy(S_new, *plane_mf_impfunc, 0, V["level_set"], 1, S_new.nGrow());
		
	}else{
		amrex::Abort("invalid EBgeom setting");
	}
	
}

void initialise_levelset_normals(MultiFab& S_old, MultiFab& S_new, const Real* dx, AccessVariable& V)
{
    Real nx, ny;
    
    for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
        const Box& box     		= mfi.validbox();
        const auto lo 			= lbound(box);
		const auto hi 			= ubound(box);
		
        FArrayBox& fab_new 	  	= S_new[mfi];
	    FArrayBox& fab_old 		= S_old[mfi];	
		
        Array4<Real> const& S_new_arr = fab_new.array();
        Array4<Real> const& S_old_arr = fab_old.array();

		for(int k = lo.z; k <= hi.z; ++k) {
			for(int j = lo.y; j <= hi.y; ++j) {
				for(int i = lo.x; i <= hi.x; ++i) {
					if(i == lo.x){
						nx = (S_old_arr(i+1,j,k,V["level_set"]) - S_old_arr(i,j,k,V["level_set"]))/(dx[0]);
					}else if(i == hi.x){
						nx = (S_old_arr(i,j,k,V["level_set"]) - S_old_arr(i-1,j,k,V["level_set"]))/(dx[0]);
					}else{
						nx = (S_old_arr(i+1,j,k,V["level_set"]) - S_old_arr(i-1,j,k,V["level_set"]))/(2*dx[0]);
					}
					if(j==lo.y){
						ny = (S_old_arr(i,j+1,k,V["level_set"]) - S_old_arr(i,j,k,V["level_set"]))/(dx[1]);
					}else if(j==hi.y){
						ny = (S_old_arr(i,j,k,V["level_set"]) - S_old_arr(i,j-1,k,V["level_set"]))/(dx[1]);
					}else{
						ny = (S_old_arr(i,j+1,k,V["level_set"]) - S_old_arr(i,j-1,k,V["level_set"]))/(2*dx[1]);
					}
					
					if(pow(nx,2) + pow(ny,2) == 0){
						S_new_arr(i,j,k,V["nx"]) = 0.0;
						S_new_arr(i,j,k,V["ny"]) = 0.0;
					}else{
						S_new_arr(i,j,k,V["nx"]) = -nx/pow(pow(nx,2) + pow(ny,2),0.5);
						S_new_arr(i,j,k,V["ny"]) = -ny/pow(pow(nx,2) + pow(ny,2),0.5);
					}
				}
			}
		} // end i,j,k loop

	}//end: MFIter
    	
}

	
void find_embedded_boundary(const Box& box, Array4<Real> const& S_new_arr, Array4<const EBCellFlag>& flags_arr, const int& nc, 
								AccessVariable& V){
	
	/*
	 find the embeddded boundary location and set EB_celltype_arr vals:
	   *  cut cells     =  0 (EB passes through cell)
	   *  covered cells = -1 (inside geom)
	   *  regular cells =  1 (outside geom)
	 */
	
	const auto lo = lbound(box);
    const auto hi = ubound(box);
    
    //find cut cells, covered cells and regular cells:
    for(int k = lo.z; k <= hi.z; ++k) {
		for(int j = lo.y; j <= hi.y; ++j) {
			for(int i = lo.x; i <= hi.x; ++i) {
				
				if(flags_arr(i,j,k).isCovered()){
					S_new_arr(i,j,k, V["EB"]) = -1;
				}else if(flags_arr(i,j,k).isSingleValued()){
					S_new_arr(i,j,k, V["EB"]) = 0;
				}else if(flags_arr(i,j,k).isRegular()){
					S_new_arr(i,j,k, V["EB"]) = 1;
				}

			}
		}
	}//end i,j,k loop
	
}

Plasma19::Plasma19(string DataFile){
			
	ifstream data(DataFile);	//"mixture19_cns.txt"
	
	if(!data){		
		amrex::Abort("data empty - file not found");
	}
	
	data >> n_rho; 		//500
	data >> n_p;		//500
	data >> n_species; 	//19
	
	//cout << "first row of data \n";
	//cout << n_rho << ' ' << n_p << ' ' << n_species << endl;
	
	v_SpecificGasConstants.resize(n_species);
	v_HeatsofFormations.resize(n_species);
	v_VibrationalTemperatures.resize(n_species);
	v_Densities.resize(n_rho);
	v_Pressures.resize(n_p);
	
	m_SoundSpeeds.resize(n_p, Vector<Real>(n_rho));
	m_InternalEnergies.resize(n_p, Vector<Real>(n_rho));
	m_Temperatures.resize(n_p, Vector<Real>(n_rho));
	m_SigmaElectrical.resize(n_p, Vector<Real>(n_rho));
	m3_SpeciesMassFractions.resize(n_species, 
						Vector< Vector<Real> >(n_p, Vector<Real>(n_rho)));
	m_SigmaThermal.resize(n_p, Vector<Real>(n_rho));
	
	//populate property vectors and matrices with the data:
	
	for(int i=0; i<n_species; ++i){
		data >> v_SpecificGasConstants[i];
	}
	for(int i=0; i<n_species; ++i){
		data >> v_HeatsofFormations[i];
	}
	for(int i=0; i<n_species; ++i){
		data >> v_VibrationalTemperatures[i];
	}
	for(int i=0; i<n_rho; ++i){
		data >> v_Densities[i];
	}
	for(int i=0; i<n_p; ++i){
		data >> v_Pressures[i];
	}
	for(int j=0; j<n_p; ++j){
		for(int i=0; i<n_rho; ++i){
			data >> m_SoundSpeeds[j][i];
		}
	}
	for(int j=0; j<n_p; ++j){
		for(int i=0; i<n_rho; ++i){
			data >> m_InternalEnergies[j][i];
		}
	}
	for(int j=0; j<n_p; ++j){
		for(int i=0; i<n_rho; ++i){
			data >> m_Temperatures[j][i];
		}
	}
	for(int j=0; j<n_p; ++j){
		for(int i=0; i<n_rho; ++i){
			data >> m_SigmaElectrical[j][i];
		}
	}
	for(int k=0; k<n_species; ++k){
		for(int j=0; j<n_p; ++j){
			for(int i=0; i<n_rho; ++i){
				data >> m3_SpeciesMassFractions[k][j][i];
			}
		}
	}
	for(int j=0; j<n_p; ++j){
		for(int i=0; i<n_rho; ++i){
			data >> m_SigmaThermal[j][i];
		}
	}
	
	cout << "plasma data initialised." << endl;
}


Real interp1D(const Vector<Real>& vec_ref, const Real& val, int &val_i, bool &flagOOB){
	
	int index;
	int vec_length = vec_ref.size();
	Real val1, val2;
	Real delta_val;
	
	if(vec_ref[0] < vec_ref[vec_length-1]){
		//property ascending order:
		if( (val < vec_ref[0]) ){
			flagOOB = 1; //Flag Out Of Bounds (below plasma data range) - switch to ideal gas
			return 0;
		}else if( (val > vec_ref[vec_length-1]) ){
			cout << "val, last : " << val << ", " << vec_ref[vec_length-1] << endl;
			amrex::Abort("Plasma19 OOB: value passed to interp1D function is out of range (high)");
		}
			
		index = 0;
		while(val >= vec_ref[index]){
			++index;
			if(index >= vec_length-1){
				amrex::Abort("error in interp1D (1)");
			}
		}
		--index;
	}else if(vec_ref[0] > vec_ref[vec_length-1]){
		//property descending order:
		if( (val > vec_ref[0]) ){
			flagOOB = 1; //Flag Out Of Bounds (below plasma data range) - switch to ideal gas
			return 0;
		}else if( (val < vec_ref[vec_length-1]) ){
			amrex::Abort("Plasma19 OOB: value passed to interp1D function is out of range (high), descending");
		}
			
		index = 0;
		while(val <= vec_ref[index]){
			++index;
			if(index >= vec_length-1){
				amrex::Abort("error in interp1D (1)");
			}
		}
		--index;
	}
	val1 = vec_ref[index];
	val2 = vec_ref[index+1];
	
	if(val2 == val1){
		cout << "val1 == val2 in interp" << endl;
		amrex::Abort("this is a problem!");
	}
	
	if(index <0 || index >= vec_length-1){
		std::cout << "index = " << index << std::endl;
		std::cout << "first, last, val = " << vec_ref[0] << ", " << vec_ref[vec_length-1] << ", " << val << endl;
		amrex::Abort("this is the problem!");
	}
	
	val_i = index;
	
	if(val2 == val1){
		delta_val = 0;
	}else{
		delta_val = (val-val1)/(val2-val1);
	}
	
	return delta_val;
		
} 

Real bilinear_interp(const Vector< Vector<Real> >& m_property, const Real& rho, const Real& p, const int& rho_i, 
						const int& p_i, const Real& delta_rho, const Real& delta_p){
	
	Real val;
	Real prop1, prop2;
	
	prop1 = m_property[p_i][rho_i] + delta_rho*(m_property[p_i][rho_i+1] - m_property[p_i][rho_i]);
	prop2 = m_property[p_i+1][rho_i] + delta_rho*(m_property[p_i+1][rho_i+1] - m_property[p_i+1][rho_i]);
	
	val = prop1 + delta_p*(prop2-prop1);
	
	
	return val;
		
}


Real Plasma19::getSpecificInternalEnergy(const Real rho, const Real p){
	
	int rho_i, p_i;
	Real delta_rho, delta_p;
	bool flagOOB = 0;
	Real e;
	
	delta_rho = interp1D(v_Densities, rho, rho_i, flagOOB);
	delta_p = interp1D(v_Pressures, p, p_i, flagOOB);
	
	if(flagOOB == 0){
		e = bilinear_interp(m_InternalEnergies, rho, p, rho_i, p_i, delta_rho, delta_p);
		return e;
	}else{
		//state below low range of plasma data
		//evaluated as ideal gas instead
		return p/((m_adiabaticIndex - 1.0)*rho);
	}
	
}

Real Plasma19::getSoundSpeed(const Real rho, const Real p){
	
	int rho_i, p_i;
	Real delta_rho, delta_p;
	bool flagOOB = 0;
	Real s;
	
	delta_rho = interp1D(v_Densities, rho, rho_i, flagOOB);
	delta_p = interp1D(v_Pressures, p, p_i, flagOOB);
	
	if(flagOOB == 0){
		s = bilinear_interp(m_SoundSpeeds, rho, p, rho_i, p_i, delta_rho, delta_p);
		return s;
	}else{
		//state below low range of plasma data
		//evaluated as ideal gas instead
		return sqrt(m_adiabaticIndex * p/rho);
	}
	
}

Real Plasma19::getPressure(const Real rho, const Real e){
	
	int rho_i, e_i;
	Real delta_rho, delta_e;
	bool flagOOB = 0;
	Real p;
	double p1, p2;
	
	delta_rho = interp1D(v_Densities, rho, rho_i, flagOOB);
	
	if(flagOOB){
		p = (m_adiabaticIndex - 1.0)*rho*e;
		return p;
	}
	
	Vector<Real> e_vec1(n_p);
	Vector<Real> e_vec2(n_p);
	
	for(int j =0; j<n_p; ++j){
		e_vec1[j] = m_InternalEnergies[j][rho_i];
		e_vec2[j] = m_InternalEnergies[j][rho_i+1];
	}

	delta_e = interp1D(e_vec1, e, e_i, flagOOB);
	//std::cout << "e_i = " << e_i << std::endl;
	if(!flagOOB){
		p1 = v_Pressures[e_i] + delta_e*(v_Pressures[e_i+1]-v_Pressures[e_i]);
	}else{
		p = (m_adiabaticIndex - 1.0)*rho*e;
		return p;
	}
	delta_e = interp1D(e_vec2, e, e_i, flagOOB);
	if(!flagOOB){
		p2 = v_Pressures[e_i] + delta_e*(v_Pressures[e_i+1]-v_Pressures[e_i]);
	}else{
		p = (m_adiabaticIndex - 1.0)*rho*e;
		return p;
	}
	
	if(!flagOOB){
		p = p1 + delta_rho*(p2-p1);
		return p;
	}else{
		//state below low range of plasma data
		//evaluated as ideal gas instead
		p = (m_adiabaticIndex - 1.0)*rho*e;
		return p;
	}
	
}

Real Plasma19::getDensity(const Real p, const Real e){
	
	int p_i, e_i;
	Real delta_p, delta_e;
	bool flagOOB = 0;
	Real rho;
	
	delta_p = interp1D(v_Pressures, p, p_i, flagOOB);
	
	Vector<Real> e_vec1(n_rho);
	Vector<Real> e_vec2(n_rho);
	
	for(int i =0; i<n_rho; ++i){
		e_vec1[i] = m_InternalEnergies[p_i][i];
		e_vec2[i] = m_InternalEnergies[p_i+1][i];
	}

	delta_e = interp1D(e_vec1, e, e_i, flagOOB);
	double rho1 = v_Densities[e_i] + delta_e*(v_Densities[e_i+1]-v_Densities[e_i]);
	delta_e = interp1D(e_vec2, e, e_i, flagOOB);
	double rho2 = v_Densities[e_i] + delta_e*(v_Densities[e_i+1]-v_Densities[e_i]);
	
	if(flagOOB == 0){
		rho = rho1 + delta_p*(rho2-rho1);
		return rho;
	}else{
		//state below low range of plasma data
		//evaluated as ideal gas instead
		return p/((m_adiabaticIndex - 1.0)*e);
	}
	
}

Real Plasma19::getTemperature(const Real rho, const Real p){
	
	int rho_i, p_i;
	Real delta_rho, delta_p;
	bool flagOOB = 0;
	Real T;
	
	delta_rho = interp1D(v_Densities, rho, rho_i, flagOOB);
	delta_p = interp1D(v_Pressures, p, p_i, flagOOB);
	
	if(flagOOB == 0){
		T = bilinear_interp(m_Temperatures, rho, p, rho_i, p_i, delta_rho, delta_p);
		return T;
	}else{
		//state below low range of plasma data
		//evaluated as ideal gas instead
		return p/(GasConstantR*rho); 
	}
	
}

Real Plasma19::getConductivity(const Real rho, const Real p){
	
	int rho_i, p_i;
	Real delta_rho, delta_p;
	bool flagOOB = 0;
	Real C;
	
	delta_rho = interp1D(v_Densities, rho, rho_i, flagOOB);
	delta_p = interp1D(v_Pressures, p, p_i, flagOOB);
	
	if(flagOOB == 0){
		C = bilinear_interp(m_SigmaElectrical, rho, p, rho_i, p_i, delta_rho, delta_p);
		return C;
	}else{
		//state below low range of plasma data
		//evaluated as ideal gas instead
		return 1.5816e-58; //this is actually heavily temperature dependent
							// should be replaced by analytic relation
	}
	
}

Real Plasma19::getThermalConductivity(const Real rho, const Real p){
	
	int rho_i, p_i;
	Real delta_rho, delta_p;
	bool flagOOB = 0;
	Real TC;
	
	delta_rho = interp1D(v_Densities, rho, rho_i, flagOOB);
	delta_p = interp1D(v_Pressures, p, p_i, flagOOB);
	
	if(flagOOB == 0){
		TC = bilinear_interp(m_SigmaThermal, rho, p, rho_i, p_i, delta_rho, delta_p);
		return TC;
	}else{
		//state below low range of plasma data
		//evaluated as ideal gas instead
		return 0.02527;  //value currently taken as constant at ambirent air conditions
							// should be replaced by analytic relation
	}
	
}

Real Plasma19::getSpeciesMassFraction(const Real rho, const Real p, const int s){
	
	int rho_i, p_i;
	Real delta_rho, delta_p;
	bool flagOOB = 0;
	Real SMF;
	
	delta_rho = interp1D(v_Densities, rho, rho_i, flagOOB);
	delta_p = interp1D(v_Pressures, p, p_i, flagOOB);
	
	if(flagOOB == 0){
		SMF = bilinear_interp(m3_SpeciesMassFractions[s], rho, p, rho_i, p_i, delta_rho, delta_p);
		return SMF;
	}else{
		//state below low range of plasma data
		//evaluated as ideal gas instead
		return 0.0; //CHECK (???) 
	}
	
}

Real Plasma19::getAdiabaticIndex(){
	
	return m_adiabaticIndex;
	
}

void Plasma19::testing_function(){
	
	Real e = getSpecificInternalEnergy(1.225, 102000);
	std::cout << "test e should be ~ 207000 J/kg " << std::endl;
	std::cout << "computed e = " << e << std::endl;
	
	Real s = getSoundSpeed(1.225, 102000);
	std::cout << "test s should be ~ 342 m/s " << std::endl;
	std::cout << "computed s = " << s << std::endl;
	
	Real p = getPressure(1.225, e);
	std::cout << "test p should be ~ 102000 Pa " << std::endl;
	std::cout << "computed p = " << p << std::endl;
	
	Real rho = getDensity(p, e);
	std::cout << "test rho should be ~ 1.225kg/m^3 " << std::endl;
	std::cout << "computed rho = " << rho << std::endl;
	
	Real T = getTemperature(1.225, 102000);
	std::cout << "test T should be ~ 290 K " << std::endl;
	std::cout << "computed T = " << T << std::endl;
	
	Real sigma = getConductivity(1.225, 102000);
	std::cout << "test sigma should be v small " << std::endl;
	std::cout << "computed sigma = " << sigma << std::endl;
	
	Real sigmaT = getThermalConductivity(1.225, 102000);
	std::cout << "test sigma_thermal should be fairly small " << std::endl;
	std::cout << "computed sigma_thermal = " << sigmaT << std::endl;
	
	Real SMF = getSpeciesMassFraction(1.225, 102000, 0);
	std::cout << "test species mass fraction - no idea which species or what this should be: " << std::endl;
	std::cout << "computed species mass fraction = " << SMF << std::endl;
	
	Real gamma = getAdiabaticIndex();
	std::cout << "test adiabatic index should = 1.4" << std::endl;
	std::cout << "computed adiabatic indexn= " << gamma << std::endl;
	
	
}








