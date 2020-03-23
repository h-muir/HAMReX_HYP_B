
#include <AmrLevelAdv.H>
#include <Adv_F.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_FluxRegister.H>

#include <AMReX_BC_TYPES.H>
#include <AMReX_BCUtil.H>

#include <AMReX_EB2.H>
#include <AMReX_EB2_Level.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EB_levelset.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EB_LSCore.H>
#include <AMReX_EB_LSCoreBase.H>
#include <AMReX_EBAmrUtil.H>
#include <AMReX_EBMultiFabUtil.H>

#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>

#include "structdefs.H"
#include "funcdefs.H"
#include "magnetism.H"


using namespace amrex;
using namespace std;

int      AmrLevelAdv::verbose         = 0;
Real     AmrLevelAdv::cfl             = 0.9;
int      AmrLevelAdv::do_reflux       = 1;

int      AmrLevelAdv::NUM_STATE       = 25;  // One variable in the state
int      AmrLevelAdv::NUM_GROW        = 2;  // number of ghost cells

ParameterStruct  AmrLevelAdv::Parameters;
SettingsStruct  AmrLevelAdv::SimSettings;

//string 	 AmrLevelAdv::EoS; //specify EoS as either: "ideal" or "plasma19"
//Plasma19 AmrLevelAdv::AirPlasma("mixture19_cns.txt");


Vector<BCRec> AmrLevelAdv::bc_vec(NUM_STATE);

//
//Default constructor.  Builds invalid object.
//
AmrLevelAdv::AmrLevelAdv ()
{
	cout << "begin AmrLevel constructor 1" << endl;
    flux_reg = 0;
}

//
//The basic constructor.
//
AmrLevelAdv::AmrLevelAdv (Amr&            papa,
     	                  int             lev,
                          const Geometry& level_geom,
                          const BoxArray& bl,
                          const DistributionMapping& dm,
                          Real            time)
    :
    AmrLevel(papa,lev,level_geom,bl,dm,time) 
{
    flux_reg = 0;
    if (level > 0 && do_reflux)
        flux_reg = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);
}

//
//The destructor.
//
AmrLevelAdv::~AmrLevelAdv () 
{
    delete flux_reg;
}

//
//Restart from a checkpoint file.
//
void
AmrLevelAdv::restart (Amr&          papa,
	              std::istream& is,
                      bool          bReadSpecial)
{
    AmrLevel::restart(papa,is,bReadSpecial);

    BL_ASSERT(flux_reg == 0);
    if (level > 0 && do_reflux)
        flux_reg = new FluxRegister(grids,dmap,crse_ratio,level,NUM_STATE);
}

void 
AmrLevelAdv::checkPoint (const std::string& dir,
		         std::ostream&      os,
                         VisMF::How         how,
                         bool               dump_old) 
{
  AmrLevel::checkPoint(dir, os, how, dump_old);
}

//
//Write a plotfile to specified directory.
//

void
AmrLevelAdv::writePlotFile (const std::string& dir, std::ostream& os, VisMF::How how)
{
	//variable types specified in AmrLevelAdv::variableSetUp () below
    AmrLevel::writePlotFile(dir, os, how);
 
}

//
//Define data descriptors.
//
void
AmrLevelAdv::variableSetUp ()
{
	
    BL_ASSERT(desc_lst.size() == 0);

    // Get options, set phys_bc
    read_params();

    desc_lst.addDescriptor(Phi_Type,IndexType::TheCellType(),
                           StateDescriptor::Point,0,NUM_STATE,
			   &cell_cons_interp);
    
    /* H-----------------------------------------------------
     * SimSettings and Parameters are member variables of 
     * AmrLevelAdv, of struct types: SettingsStruct and 
     * ParameterStruct respectively, defined in StructDefs
     * -----------------------------------------------------*/ 
    
    initialiseStructs(SimSettings, Parameters);
    
    AccessVariable V(Parameters.problem_variables);
    
    setBoundaryConditions(bc_vec, Parameters, NUM_STATE, V);
    
    if(Parameters.IC == "mach_theta"){
		if(Parameters.vars[0] == "rho" && Parameters.vars[1] == "u" && \
		   Parameters.vars[2] == "v"   && Parameters.vars[3] == "p"){
			Vector<Real> u_vec{Parameters.PR[1], Parameters.PR[2], 0};
			Prim WR(Parameters.PR[0], u_vec, Parameters.PR[3]);
			Vector<Real> normal_vec{cos(Parameters.theta),-sin(Parameters.theta),0};
			Prim WL = rankine_hugoniot(WR, Parameters.mach, normal_vec);
			//automatically compute downstream shock state from upstream state
			//as a function of mach number and angle (rankine hugoniot relations)
			Parameters.PL[0] = WL.rho;
			Parameters.PL[1] = WL.u_vec[0]; 
			Parameters.PL[2] = WL.u_vec[1];
			Parameters.PL[3] = WL.p;
		}else{
			amrex::Abort("P vars input properties not input as assumed");
		}
	}
    
    if(NUM_STATE != SimSettings.NCOMP){
		amrex::Abort("NUM_STATE != NCOMP in input file");
	}
    
    Vector<std::string> variable_names = Parameters.problem_variables;
    
    StateDescriptor::BndryFunc 	func;
    func = StateDescriptor::BndryFunc(phifill);
    
    
	for(int comp = 0; comp < NUM_STATE; ++comp){
		amrex::Print() << "variable name " << comp << ": " << variable_names[comp] << "\n";
		desc_lst.setComponent(Phi_Type, comp, variable_names[comp], bc_vec[comp], 
			  func);
	}
	
	cout << "end of variableSetUp()" << endl;
	
}

//
//Cleanup data descriptors at end of run.
//
void
AmrLevelAdv::variableCleanUp () 
{
    desc_lst.clear();
}

//
//Initialize grid data at problem start-up.
//
void
AmrLevelAdv::initData ()
{
    //
    // Loop over grids, call FORTRAN function to init with data.
    //
    const Real* dx  = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();
    MultiFab& S_new = get_new_data(Phi_Type);
    Real cur_time   = state[Phi_Type].curTime();
    Real prev_time   = state[Phi_Type].prevTime();

    const Box& prob_domain = geom.Domain(); 

    amrex::Print() << "nsteps_max: " << Parameters.nsteps_max << "\n";
    
    AccessVariable V(Parameters.problem_variables);

    if (verbose) {
        amrex::Print() << "Initializing the data at level " << level << std::endl;
    }
    
    //create new multifab with ngrow:
	MultiFab S_old(S_new.boxArray(), dmap, NUM_STATE, NUM_GROW); //H
	FillPatch(*this, S_old, NUM_GROW, cur_time, Phi_Type, 0, NUM_STATE); //H
	special_boundary(prob_domain, S_old, SimSettings, dx, V);  
    
    //EBGeom:
	const Vector<int> ngrow_vec{0,0,0};   
	std::unique_ptr<amrex::EBFArrayBoxFactory> factory_ptr = makeEBFabFactory(geom, S_new.boxArray(), dmap, ngrow_vec, EBSupport::full);
    auto factory = dynamic_cast<EBFArrayBoxFactory const*>(&(S_new.Factory()));  
    const FabArray<EBCellFlagFab>& flags = factory->getMultiEBCellFlagFab();
    MultiCutFab const& centroid = factory->getCentroid();
    
	initialise_problem(S_new, flags, geom, SimSettings, Parameters, V);
	
	//set boundaries:
	S_new.FillBoundary(geom.periodicity());
	special_boundary(prob_domain, S_new, SimSettings, dx, V);  
	
	if(SimSettings.EBgeom != "none"){
		
		initialise_levelset_geometry(S_new, level, geom, SimSettings, Parameters, V);
		
		MultiFab::Copy(S_old, S_new, 0, 0, NUM_STATE, S_new.nGrow());
		FillPatch(*this, S_old, NUM_GROW, cur_time, Phi_Type, 0, NUM_STATE);
		special_boundary(prob_domain, S_old, SimSettings, dx, V); 
		
		initialise_levelset_normals(S_old, S_new, dx, V);
		
		MultiFab::Copy(S_old, S_new, V["nx"], V["nx"], 2, S_new.nGrow());
		FillPatch(*this, S_old, NUM_GROW, cur_time, Phi_Type, 0, NUM_STATE);
		special_boundary(prob_domain, S_old, SimSettings, dx, V);
		
		rigid_body_step1(S_new, S_old, flags, dx, V);
		
		MultiFab::Copy(S_new, S_old, 0, 0, NUM_STATE, S_new.nGrow());
		FillPatch(*this, S_old, NUM_GROW, cur_time, Phi_Type, 0, NUM_STATE);
		special_boundary(prob_domain, S_old, SimSettings, dx, V);	
		
		double global_p_min = -1e5;
		int counter = 0;
		
		//Rigid body step 2:
		while(global_p_min < 0){
			++counter;
			if(counter >= 30){break;} 
			global_p_min = 1.0;   
			
			rigid_body_step2(S_new, S_old, flags, dx, V, global_p_min);
			
			MultiFab::Copy(S_new, S_old, 0, 0, NUM_STATE, S_new.nGrow());
			FillPatch(*this, S_old, NUM_GROW, cur_time, Phi_Type, 0, NUM_STATE);
			special_boundary(prob_domain, S_old, SimSettings, dx, V);

			ParallelDescriptor::ReduceRealMax(global_p_min);
			//cout << "global p min = " << global_p_min << endl;
		}
		
		rigid_body_step2(S_new, S_old, flags, dx, V, global_p_min);
			
		MultiFab::Copy(S_new, S_old, 0, 0, NUM_STATE, S_new.nGrow());
		FillPatch(*this, S_old, NUM_GROW, cur_time, Phi_Type, 0, NUM_STATE);
		special_boundary(prob_domain, S_old, SimSettings, dx, V);
		
		for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
		{
				const Box& box     = mfi.validbox();
				FArrayBox& fab 	   = S_new[mfi];
				
				Array4<Real> const& S_new_arr = fab.array();

				S_conWtoU(box, S_new_arr, V);
		}
		
	
	}//end (if EBgeom != none)

	
    if (verbose) {
	amrex::Print() << "Done initializing the level " << level 
                       << " data " << std::endl; 
    }

}

//
//Initialize data on this level from another AmrLevelAdv (during regrid).
//
void
AmrLevelAdv::init (AmrLevel &old)
{  
    AmrLevelAdv* oldlev = (AmrLevelAdv*) &old;
    //
    // Create new grid data by fillpatching from old.
    //
    Real dt_new    = parent->dtLevel(level);
    Real cur_time  = oldlev->state[Phi_Type].curTime();
    Real prev_time = oldlev->state[Phi_Type].prevTime();
    Real dt_old    = cur_time - prev_time;
    setTimeLevel(cur_time,dt_old,dt_new);

    MultiFab& S_new = get_new_data(Phi_Type);

    FillPatch(old, S_new, 0, cur_time, Phi_Type, 0, NUM_STATE);
}

//
//Initialize data on this level after regridding if old level did not previously exist
//
void
AmrLevelAdv::init ()
{
    Real dt        = parent->dtLevel(level);
    Real cur_time  = getLevel(level-1).state[Phi_Type].curTime();
    Real prev_time = getLevel(level-1).state[Phi_Type].prevTime();

    Real dt_old = (cur_time - prev_time)/(Real)parent->MaxRefRatio(level-1);

    setTimeLevel(cur_time,dt_old,dt);
    MultiFab& S_new = get_new_data(Phi_Type);
    FillCoarsePatch(S_new, 0, cur_time, Phi_Type, 0, NUM_STATE);
}

//
//Advance grids at this level in time.
//
Real
AmrLevelAdv::advance (Real time,
                      Real dt,
                      int  iteration,
                      int  ncycle)
{

    MultiFab& S_mm = get_new_data(Phi_Type);
    Real maxval = S_mm.max(0);
    Real minval = S_mm.min(0);
    
    amrex::Print() << "phi max = " << maxval << ", min = " << minval  << std::endl;
    
    const Real prev_time = state[Phi_Type].prevTime();
    const Real cur_time = state[Phi_Type].curTime();
    const Real ctr_time = 0.5*(prev_time + cur_time);
    
    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();
    
    AccessVariable V(Parameters.problem_variables);

    MultiFab& S_new = get_new_data(Phi_Type);
    
    MultiFab S_old(S_new.boxArray(), dmap, S_new.nComp(), NUM_GROW); //H
    FillPatch(*this, S_old, NUM_GROW, time, Phi_Type, 0, NUM_STATE); //H
    const Box& prob_domain = geom.Domain();
    //S_old.FillBoundary(geom.periodicity());
   // FillDomainBoundary(S_old, geom, bc_vec);
	special_boundary(prob_domain, S_old, SimSettings, dx, V);

    //
    // Get pointers to Flux registers, or set pointer to zero if not there.
    //
    FluxRegister *fine    = 0;
    FluxRegister *current = 0;
    
    int finest_level = parent->finestLevel();

    if (do_reflux && level < finest_level) {
		fine = &getFluxReg(level+1);
		fine->setVal(0.0);
    }

    if (do_reflux && level > 0) {
		current = &getFluxReg(level);
    }

    MultiFab fluxes[BL_SPACEDIM]; //Q: does this make an array of MultiFabs, of length BL_SPACEDIM,
									//like python syntax? 

    if (do_reflux)
    {
	for (int j = 0; j < BL_SPACEDIM; j++)
	{
	    BoxArray ba = S_new.boxArray();
	    ba.surroundingNodes(j);
	    fluxes[j].define(ba, dmap, NUM_STATE, 0);
	}
    }
	
#ifdef _OPENMP
#pragma omp parallel
#endif
    
    
    //H------------------------------------------------------------------------- 
    //create a multifab array to store the fluxes:
	MultiFab flux_arr[BL_SPACEDIM];
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
        // flux(dir) has one component, zero ghost cells, and is nodal in direction dir
        BoxArray edge_ba = S_new.boxArray();
        edge_ba.surroundingNodes(dir); //switches box from type:cell to type:node in direction dir ? 
        flux_arr[dir].define(edge_ba, dmap, NUM_STATE, 0); //ngrow=0

    }
    //------------------------------------------------------------------------H/
    
    //-------------------------------------------------------------------------
    
    {
	//FArrayBox flux[BL_SPACEDIM];
	//H-------------------------------------------------------------------------
	Real dxr = dx[0];
	Real dyr = dx[1];
	
	//HALF TIME STEP RADIAL SOURCE TERM INTEGRATION (half step 1):
	if(Parameters.coord_sys == 1){
		radial_sourceterm_integration(S_old, S_new, dt/2.0, dxr, V);
		MultiFab::Copy(S_old, S_new, 0, 0, NUM_STATE, S_new.nGrow());
		FillPatch(*this, S_old, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
		special_boundary(prob_domain, S_old, SimSettings, dx, V);
	}

    //x-sweep:
	for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
	{
		//H --------------------------------------------------------------------
	    const Box& bx = mfi.tilebox();

	    FArrayBox& stateout      =   S_new[mfi];
	    FArrayBox& statein 		 = 	 S_old[mfi];	
	    FArrayBox& flux_fab_x 	 = 	 flux_arr[0][mfi];	
	    	    
	    //--------------------------------------------------------------------H/		

		//H --------------------------------------------------------------------
		//access the array from the array box:
		
        Array4<Real> const& S_old_arr 		= statein.array();
        Array4<Real> const& S_new_arr 		= stateout.array();
        Array4<Real> const& flux_prop_arr_x = flux_fab_x.array();
		
		if(SimSettings.solver == "HLL"){
			calc_fluxes_HLL_x(bx, S_old_arr, flux_prop_arr_x, dt, dxr, NUM_STATE, V, SimSettings.MUSCL);
		}else if(SimSettings.solver == "HLLC"){
			calc_fluxes_HLLC_x(bx, S_old_arr, flux_prop_arr_x, dt, dxr, NUM_STATE, V, SimSettings.MUSCL);
		}else if(SimSettings.solver == "HLLC-HS"){
			calc_fluxes_HLLC_HS_x(bx, S_old_arr, flux_prop_arr_x, dt, dxr, NUM_STATE, V, SimSettings.MUSCL);
		}else{
			amrex::Abort("invalid solver - check settings file");
		}
		cons_update("x-sweep", bx, flux_prop_arr_x, S_old_arr, S_new_arr, dt, dxr, NUM_STATE, V);
		
        
        //--------------------------------------------------------------------H/
		    		
	    if (do_reflux) {
			rescale_fluxes("x-sweep", bx, flux_prop_arr_x, dt, dyr, NUM_STATE, V); //rescale by dt*dy for Fx
			//for (int i = 0; i < BL_SPACEDIM ; i++)
			fluxes[0][mfi].copy(flux_fab_x,mfi.nodaltilebox(0));
			//fluxes[i][mfi].copy(flux[i],mfi.nodaltilebox(i));
	    }
	    
	} //end x-sweepMFIter	
	
    //copy data from updated S_new --> S_old between dimensional sweeps:
	MultiFab::Copy(S_old, S_new, 0, 0, NUM_STATE, S_new.nGrow());
    FillPatch(*this, S_old, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
    special_boundary(prob_domain, S_old, SimSettings, dx, V);
    //note: S_new.nGrow() = 0
    
	//y-sweep:--------------------------------------------------------------------------------------
	for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
	{
		//H --------------------------------------------------------------------
	    const Box& bx = mfi.tilebox();

	    FArrayBox& stateout      =   S_new[mfi];
	    FArrayBox& statein 		 = 	 S_old[mfi];	//H
	    FArrayBox& flux_fab_y 	 =   flux_arr[1][mfi];	//H
	    //H --------------------------------------------------------------------
		
		//H --------------------------------------------------------------------
		//access the array from the array box:
		
        Array4<Real> const& S_old_arr = statein.array();
        Array4<Real> const& S_new_arr = stateout.array();
        Array4<Real> const& flux_prop_arr_y = flux_fab_y.array();
		
		if(SimSettings.solver == "HLLC"){
			calc_fluxes_HLLC_y(bx, S_old_arr, flux_prop_arr_y, dt, dyr, NUM_STATE, V, SimSettings.MUSCL);
		}else if(SimSettings.solver == "HLL"){
			calc_fluxes_HLL_y(bx, S_old_arr, flux_prop_arr_y, dt, dyr, NUM_STATE, V, SimSettings.MUSCL);
		}else if(SimSettings.solver == "HLLC-HS"){
			calc_fluxes_HLLC_HS_y(bx, S_old_arr, flux_prop_arr_y, dt, dyr, NUM_STATE, V, SimSettings.MUSCL);
		}else{
			amrex::Abort("invalid solver defined in settings file");
		}
		
        cons_update("y-sweep", bx, flux_prop_arr_y, S_old_arr, S_new_arr, dt, dyr, NUM_STATE, V);
        
        //----------------------------------------------------------------------

		    		
	    if (do_reflux) {
			rescale_fluxes("y-sweep", bx, flux_prop_arr_y, dt, dxr, NUM_STATE, V); //rescale by dt*dx for Fy
			//for (int i = 0; i < BL_SPACEDIM ; i++)
			fluxes[1][mfi].copy(flux_fab_y,mfi.nodaltilebox(1));
			//fluxes[i][mfi].copy(flux[i],mfi.nodaltilebox(i));
	    }
	    
	} //end y-sweep MFIter	
	MultiFab::Copy(S_old, S_new, 0, 0, NUM_STATE, S_new.nGrow());
    FillPatch(*this, S_old, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
    special_boundary(prob_domain, S_old, SimSettings, dx, V);
	
	//HALF TIME STEP RADIAL SOURCE TERM INTEGRATION (half step 2):
	if(Parameters.coord_sys == 1){
		radial_sourceterm_integration(S_old, S_new, dt/2.0, dxr, V);
	}
	 	 
    } //end local scope
    
	
    if (do_reflux) {
		
		if (current) {
			//for (int i = 0; i < BL_SPACEDIM ; i++)
			current->FineAdd(fluxes[0],0,0,0,NUM_STATE,1.);
			current->FineAdd(fluxes[1],1,0,0,NUM_STATE,1.);
		}
		
		if (fine) {
			//for (int i = 0; i < BL_SPACEDIM ; i++)
			fine->CrseInit(fluxes[0],0,0,0,NUM_STATE,-1.);
			fine->CrseInit(fluxes[1],1,0,0,NUM_STATE,-1.);
		}
    }
	
    //------------------------EMBEDDED BOUNDARY--------------------------------
	if(SimSettings.EBgeom != "none"){
		const Vector<int> ngrow_vec{NUM_GROW,NUM_GROW,NUM_GROW}; //not sure if this should be {0,0,0} or {NUM_GROW,NUM_GROW,NUM_GROW}? 
		std::unique_ptr<amrex::EBFArrayBoxFactory> factory_ptr = makeEBFabFactory(geom, S_new.boxArray(), dmap, ngrow_vec, EBSupport::full);
		auto factory = dynamic_cast<EBFArrayBoxFactory const*>(&(S_new.Factory()));  
		const FabArray<EBCellFlagFab>& flags = factory->getMultiEBCellFlagFab();
		MultiCutFab const& centroid = factory->getCentroid();

		rigid_body_step1(S_new, S_old, flags, dx, V);
			
		MultiFab::Copy(S_new, S_old, 0, 0, NUM_STATE, S_new.nGrow());
		FillPatch(*this, S_old, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
		special_boundary(prob_domain, S_old, SimSettings, dx, V);	
		
		double global_p_min = -1e5;
		int counter = 0;
		
		//Rigid body step 2:
		while(global_p_min < 0 || counter < 10){
			++counter;
			if(counter >= 100){break;} 
			global_p_min = 1.0;   
			
			rigid_body_step2(S_new, S_old, flags, dx, V, global_p_min);
			
			MultiFab::Copy(S_new, S_old, 0, 0, NUM_STATE, S_new.nGrow());
			FillPatch(*this, S_old, NUM_GROW, cur_time, Phi_Type, 0, NUM_STATE);
			special_boundary(prob_domain, S_old, SimSettings, dx, V);

			ParallelDescriptor::ReduceRealMax(global_p_min);
			//cout << "global p min = " << global_p_min << endl;
		}
		
		rigid_body_step2(S_new, S_old, flags, dx, V, global_p_min);
	
		MultiFab::Copy(S_new, S_old, 0, 0, NUM_STATE, S_new.nGrow());
		FillPatch(*this, S_old, NUM_GROW, time, Phi_Type, 0, NUM_STATE);
		special_boundary(prob_domain, S_old, SimSettings, dx, V);
		
		for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
		{
				const Box& box     = mfi.validbox();
				FArrayBox& fab 	   = S_new[mfi];
				
				Array4<Real> const& S_new_arr = fab.array();

				S_conWtoU(box, S_new_arr, V);
		}
		
		special_boundary(prob_domain, S_new, SimSettings, dx, V); 
	}//end EBgeom section
	
	//------------------------EMField----------------------------------

    //EM field:
	if(SimSettings.Bfield == "dipole"){
		
		cout << "EMField Linear solve on level: " << level << endl;
		/*
		bool has_old_data = state[0].hasOldData();
  
		// allocOldData does old_data.reset if old_data is null, otherwise does nothing
		state[0].allocOldData();
		if(!has_old_data)
		{
		state[0].oldData().setVal(0.0);
		}
		// Moves current data into old data (and vice-versa)
		state[0].swapTimeLevels(dt);
		
		MultiFab& S_old_data = get_old_data(Phi_Type);
		
		MultiFab& S_new_data = get_new_data(Phi_Type);		
		*/
		EMField EM_lev(Parameters, geom, grids, dmap);
		
		LPInfo info;
		info.setAgglomeration(EM_lev.agglomeration);
		info.setConsolidation(EM_lev.consolidation);
		info.setMaxCoarseningLevel(EM_lev.max_coarsening_level);
		//info.setMetricTerm(false);
		
		MLABecLaplacian mlabec({geom}, {grids}, {dmap}, info);
		
		mlabec.setMaxOrder(EM_lev.linop_maxorder);
		
		mlabec.setDomainBC({AMREX_D_DECL(LinOpBCType::Neumann,   //x-lo
									 LinOpBCType::Neumann,   //y-lo
									 LinOpBCType::Neumann)}, //z-lo
							{AMREX_D_DECL(LinOpBCType::Dirichlet,   //x-hi
									 LinOpBCType::Dirichlet,   //y-hi
									 LinOpBCType::Neumann)});//z-hi
		
		int ng = 1;	
		MultiFab Sol_crse;
		
		if(level > 0){
			AmrLevelAdv& crse_level = getLevel(level-1);
			Real crse_time = crse_level.state[0].curTime();
			const BoxArray& crse_box = crse_level.boxArray();
			const DistributionMapping& crse_dmap = crse_level.DistributionMap();
			Sol_crse.define(crse_box,parent->DistributionMap(level-1),1,ng);
			FillPatch(crse_level, Sol_crse, ng, time, Phi_Type, V["phi"], 1);
			mlabec.setCoarseFineBC(&Sol_crse, EM_lev.ref_ratio);
		}
		
		MultiFab S_EM(grids, dmap, NUM_STATE, ng);
		FillPatch(*this, S_EM, ng, time, Phi_Type, 0, NUM_STATE);
		S_EM.FillBoundary();
		FillDomainBoundary(S_EM, geom, bc_vec);
		
		EM_lev.defineSolution(S_EM, V);
		mlabec.setLevelBC(0,&EM_lev.solution);
		
		EM_lev.computeAlphaFab();
		mlabec.setACoeffs(0, EM_lev.alpha_fab);
		
		EM_lev.defineScalars();
		mlabec.setScalars(EM_lev.A_scalar, EM_lev.B_scalar);
		
		
		EM_lev.computeFaceSigmaFabs(S_EM, V);
		mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(EM_lev.face_sigma_fabs));
		
		EM_lev.computeRHS(S_EM, V);
		
		MLMG mlmg(mlabec);
		mlmg.setMaxIter(EM_lev.max_iter);
		mlmg.setMaxFmgIter(EM_lev.max_fmg_iter);
		mlmg.setVerbose(EM_lev.verbose);
		mlmg.setBottomVerbose(EM_lev.bottom_verbose);
		
		mlmg.solve({&EM_lev.solution}, {&EM_lev.rhs}, EM_lev.tol_rel, EM_lev.tol_abs);
		
		
		MultiFab::Copy(S_EM, EM_lev.solution, 0, V["phi"], 1, S_EM.nGrow()); 		//fills S_EM with phi solution
		MultiFab gradPhi(grids, dmap, 3, 0);
		EM_lev.compute_gradPhi(geom, EM_lev.solution, gradPhi);
		EM_lev.computeJfield(S_EM, gradPhi, EM_lev.VcrossB, EM_lev.sigma_fab, V); 		//fills S_EM with J_R, J_z, J_theta, J_abs solution
		
		MultiFab::Copy(S_new, S_EM, 0, 0, NUM_STATE, S_new.nGrow());
		
	}
    //-----------------------------------------------------------------

    return dt;
}

//
//Estimate time step.
//
Real
AmrLevelAdv::estTimeStep (Real)
{
    // This is just a dummy value to start with 
    Real dt_est  = 1.0e+20;

    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();
    const Real cur_time = state[Phi_Type].curTime();
    MultiFab& S_new = get_new_data(Phi_Type);

#ifdef _OPENMP
#pragma omp parallel reduction(min:dt_est)
#endif
    {
	AccessVariable V(Parameters.problem_variables);
	for (MFIter mfi(S_new, true); mfi.isValid(); ++mfi)
	{
		
		const Box& bx = mfi.tilebox();
	    FArrayBox& stateout     	  =   		S_new[mfi];
	    Array4<Real> const& S_new_arr =   stateout.array();
		
	    for (int i = 0; i < 2; ++i) { //configured for x and y directions only
			Real umax = max_dim_wave_speed(bx, S_new_arr, i, V); //uface[i].norm(0);
			if (umax > 1.e-100) {
				dt_est = std::min(dt_est, dx[i] / umax);
			}
	    }
	}
    }

    ParallelDescriptor::ReduceRealMin(dt_est);
    dt_est *= cfl;

    if (verbose) {
	amrex::Print() << "AmrLevelAdv::estTimeStep at level " << level 
                       << ":  dt_est = " << dt_est << std::endl;
    }
    
    return dt_est;
}

//
//Compute initial time step.
//
Real
AmrLevelAdv::initialTimeStep ()
{
    return estTimeStep(0.0);
}

//
//Compute initial `dt'.
//
void
AmrLevelAdv::computeInitialDt (int                   finest_level,
							   int                   sub_cycle,
                               Vector<int>&           n_cycle,
                               const Vector<IntVect>& ref_ratio,
                               Vector<Real>&          dt_level,
                               Real                  stop_time)
{
    //
    // Grids have been constructed, compute dt for all levels.
    //
    if (level > 0)
        return;

    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        dt_level[i] = getLevel(i).initialTimeStep();
        n_factor   *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_level[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[Phi_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

//
//Compute new `dt'.
//
void
AmrLevelAdv::computeNewDt (int                   finest_level,
		           int                   sub_cycle,
                           Vector<int>&           n_cycle,
                           const Vector<IntVect>& ref_ratio,
                           Vector<Real>&          dt_min,
                           Vector<Real>&          dt_level,
                           Real                  stop_time,
                           int                   post_regrid_flag)
{
    //
    // We are at the end of a coarse grid timecycle.
    // Compute the timesteps for the next iteration.
    //
    if (level > 0)
        return;

    for (int i = 0; i <= finest_level; i++)
    {
        AmrLevelAdv& adv_level = getLevel(i);
        dt_min[i] = adv_level.estTimeStep(dt_level[i]);
    }

    if (post_regrid_flag == 1) 
    {
	//
	// Limit dt's by pre-regrid dt
	//
	for (int i = 0; i <= finest_level; i++)
	{
	    dt_min[i] = std::min(dt_min[i],dt_level[i]);
	}
    }
    else 
    {
	//
	// Limit dt's by change_max * old dt
	//
	static Real change_max = 1.1;
	for (int i = 0; i <= finest_level; i++)
	{
	    dt_min[i] = std::min(dt_min[i],change_max*dt_level[i]);
	}
    }
    
    //
    // Find the minimum over all levels
    //
    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_min[i]);
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001*dt_0;
    Real cur_time  = state[Phi_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

//
//Do work after timestep().
//
void
AmrLevelAdv::post_timestep (int iteration)
{
    //
    // Integration cycle on fine level grids is complete
    // do post_timestep stuff here.
    //
    int finest_level = parent->finestLevel();

    if (do_reflux && level < finest_level)
        reflux();

    if (level < finest_level)
        avgDown();
}

//
//Do work after regrid().
//
void
AmrLevelAdv::post_regrid (int lbase, int new_finest) {
//particle function (removed from here)
}

//
//Do work after a restart().
//
void
AmrLevelAdv::post_restart() 
{
//particle function (removed from here)
}

//
//Do work after init().
//
void
AmrLevelAdv::post_init (Real stop_time)
{
    if (level > 0)
        return;
    //
    // Average data down from finer levels
    // so that conserved data is consistent between levels.
    //
    int finest_level = parent->finestLevel();
    for (int k = finest_level-1; k>= 0; k--)
        getLevel(k).avgDown();
}

//
//Error estimation for regridding.
//
void
AmrLevelAdv::errorEst (TagBoxArray& tags,
	               int          clearval,
                       int          tagval,
                       Real         time,
                       int          n_error_buf,
                       int          ngrow)
{
    const Real* dx        = geom.CellSize();
    const Real* prob_lo   = geom.ProbLo();

    MultiFab& S_new = get_new_data(Phi_Type);
    
    //with bordering ghost cells:
    MultiFab S_old(S_new.boxArray(), dmap, S_new.nComp(), S_new.nGrow()+1); //H
    FillPatch(*this, S_old, S_new.nGrow()+1, time, Phi_Type, 0, NUM_STATE); //H
    
    //cout << "level: " << level << endl;

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
	AccessVariable V(Parameters.problem_variables);
	string	prop = Parameters.refinement_prop;		//refinement based off gradient of property: mass
	string refinement_condition = Parameters.refinement_condition;
	Real grad_frac = Parameters.refinement_grad_fracs[level];
	Real grad_max = 0;
	
	//determine max rho gradient:
	for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
	{
	    const Box&  bx  				= mfi.tilebox();
	    
	    FArrayBox& fab  				= S_old[mfi];
	    
	    Array4<Real> const& S_old_arr 	= fab.array();

        TagBox& tagfab  				= tags[mfi]; //tags is passed to ErrorEst as argument
        
        Array4<char> const& tagarr 		= tagfab.array();
        

		determine_grad_max(bx,S_old_arr,dx, prop, refinement_condition, grad_max, V);
		

	}
	ParallelDescriptor::ReduceRealMax(grad_max);	
	//cout << "grad_max = " << grad_max << endl;
	//tagging based on propertion of max_rho
	for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
	{
	    const Box&  bx  				= mfi.tilebox();
	    
	    FArrayBox& fab  				= S_old[mfi];
	    
	    Array4<Real> const& S_old_arr 	= fab.array();

        TagBox& tagfab  				= tags[mfi]; //creates tagbox with MFI dimensions ?
        
        Array4<char> const& tagarr 		= tagfab.array();
 

		amr_tagging(tagarr,bx,S_old_arr,dx,grad_max, grad_frac, prop, refinement_condition, V);

	    
	}
    }//end local scope
    MultiFab::Copy(S_new, S_old, 0, 0, NUM_STATE, S_new.nGrow()); //copy compute grad(rho) to S_new 
}

//can probably remove this whole function:
void
AmrLevelAdv::read_params ()
{
	
    static bool done = false;

    if (done) return;

    done = true;

    ParmParse pp("adv");   

    pp.query("v",verbose);
    pp.query("cfl",cfl);
    pp.query("do_reflux",do_reflux);
	
    //
    // read tagging Parameters from probin file
    //

    std::string probin_file("probin");

    ParmParse ppa("amr");
    ppa.query("probin_file",probin_file);

    int probin_file_length = probin_file.length();
    Vector<int> probin_file_name(probin_file_length);

    for (int i = 0; i < probin_file_length; i++)
	probin_file_name[i] = probin_file[i];

    // use a fortran routine to
    // read in tagging parameters from probin file
    get_tagging_params(probin_file_name.dataPtr(), &probin_file_length);
    
    //cout << "end of read_params() function" << endl;

}

void
AmrLevelAdv::reflux ()
{
    BL_ASSERT(level<parent->finestLevel());

    const Real strt = amrex::second();

    getFluxReg(level+1).Reflux(get_new_data(Phi_Type),1.0,0,0,NUM_STATE,geom);
    
    if (verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real      end    = amrex::second() - strt;
	
        ParallelDescriptor::ReduceRealMax(end,IOProc);
	
        amrex::Print() << "AmrLevelAdv::reflux() at level " << level 
                       << " : time = " << end << std::endl;
    }
}

void
AmrLevelAdv::avgDown ()
{
    if (level == parent->finestLevel()) return;
    avgDown(Phi_Type);
}

void
AmrLevelAdv::avgDown (int state_indx)
{
    if (level == parent->finestLevel()) return;

    AmrLevelAdv& fine_lev = getLevel(level+1);
    MultiFab&  S_fine   = fine_lev.get_new_data(state_indx);
    MultiFab&  S_crse   = get_new_data(state_indx);
    
    amrex::average_down(S_fine,S_crse,
                         fine_lev.geom,geom,
                         0,S_fine.nComp(),parent->refRatio(level));
}

