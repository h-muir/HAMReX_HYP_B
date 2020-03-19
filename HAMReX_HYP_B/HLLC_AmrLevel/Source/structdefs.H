#ifndef STRUCTDEFS_H_
#define STRUCTDEFS_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>
#include <string>

using namespace amrex;

struct SettingsStruct
{
	public:
		 int nsteps_max;
		 int plot_int;
		 int ndmin;
		 Real startT;
		 Real finalT;
		 bool chatter;
		 std::string testcase;
		 std::string solver;
		 std::string EBgeom;
		 int NCOMP;
		 int NGROW;
		 bool MUSCL;
		 std::string scheme;
		 std::string EoS;
		 std::string Bfield;
			 
};

struct ParameterStruct  : public SettingsStruct
{
	public:
		 Vector<Real> prob_lo;
		 Vector<Real> prob_hi;
		 Vector<int>  n_cells;
		 Vector<int> is_periodic; 	//boundary conditions
		 Vector<Real> dx;
		 int coord_sys; 			//cartesian, cylindrical or pherical
		
		 //Test specific conditions
		 Vector<std::string> problem_variables;
		 int n_vars;
		 Vector<std::string> vars;		
		 std::string IC;		 
		 Vector<Real> PL;
		 Vector<Real> PR;
		 Real x0;
		 Real mach;
		 Real theta;
		 
		 int max_level;
		 int max_grid_size;
		 std::string refinement_prop;
		 std::string refinement_condition;
		 Vector<Real> refinement_grad_fracs;
		 
		 Vector<std::string> bc_conds;
		 
		 //embedded boundary parameters:
		 Real EB_R;
		 Vector<Real> EB_centre;
		 Vector<Real> EB_point;
		 Vector<Real> EB_normal;
		 Vector<Real> EB_point2;
		 Vector<Real> EB_normal2;
		 Vector<Real> EB_point3;
		 Vector<Real> EB_normal3;
		 
		 //magnetic field parameters
		 std::string Bfield;
		 Vector<Real> dipole_m;
		 Vector<Real> dipole_centre;
		 

    
    ParameterStruct()
    {
		prob_lo 	= Vector<Real>(AMREX_SPACEDIM);
		prob_hi 	= Vector<Real>(AMREX_SPACEDIM);
		n_cells 	= Vector<int>(AMREX_SPACEDIM);
		is_periodic = Vector<int>(AMREX_SPACEDIM);
		dx 		 	= Vector<Real>(AMREX_SPACEDIM);
		vars		= Vector<std::string>(n_vars);
		PL		 	= Vector<Real>(n_vars);
		PR		 	= Vector<Real>(n_vars);
		
		refinement_grad_fracs = Vector<Real>(max_level);
		
		problem_variables = Vector<std::string>(NCOMP);
		bc_conds		  = Vector<std::string>(4);
		
		EB_centre 		 = Vector<Real>(BL_SPACEDIM);
		EB_point 		 = Vector<Real>(AMREX_SPACEDIM);
		EB_normal		 = Vector<Real>(AMREX_SPACEDIM);
		EB_point2 		 = Vector<Real>(AMREX_SPACEDIM);
		EB_normal2		 = Vector<Real>(AMREX_SPACEDIM);
		EB_point3 		 = Vector<Real>(AMREX_SPACEDIM);
		EB_normal3		 = Vector<Real>(AMREX_SPACEDIM);
		
		dipole_m		 = Vector<Real>(AMREX_SPACEDIM);
		dipole_centre	 = Vector<Real>(AMREX_SPACEDIM);
		
	}
};





#endif
