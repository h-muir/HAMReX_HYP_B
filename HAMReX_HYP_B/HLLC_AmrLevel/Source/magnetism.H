#ifndef EMFIELD_H_
#define EMFIELD_H_

#include "structdefs.H"
#include "classdefs.H"

#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>

using namespace amrex;

class EMField
{
	public:
	
	EMField(ParameterStruct& p, Geometry const& geom, AccessVariable& V);
	void initialise_params(ParameterStruct& p);
	void computeVcrossB(MultiFab& S_EM, MultiFab &S_out, AccessVariable V);
	void compute_divVcrossB(Geometry const& geom, MultiFab& sigma_fab, MultiFab& S_VB, MultiFab &S_out);
	void compute_gradPhi(Geometry const& geom, MultiFab& S_phi, MultiFab &S_out);
	void computeJfield(MultiFab& S_EM, MultiFab& S_gradPhi, MultiFab& S_VB, MultiFab& sigma_fab, AccessVariable V);
	void initialise_fabs(MultiFab& S, Geometry const& geom, AccessVariable& V);
	void set_Derichlet_BC(Geometry const& geom, MultiFab &S_border, int boundary, int ng, Real val);
	void solve_by_level(MLABecLaplacian &mlabec, MultiFab &S_EM, MultiFab &Sol_crse, int const& ilev, AccessVariable& V);
	
	
	
	int max_level;
    int ref_ratio;
    Vector<int> n_cells;
    int max_grid_size;


	// For MLMG solver
    int verbose = 2;
    int bottom_verbose = 0;
    int max_iter = 100;
    int max_fmg_iter = 0;
    //Real reltol = 1.e-11;
    int linop_maxorder = 2;
    bool agglomeration = true;
    bool consolidation = true;
    int max_coarsening_level = 30;
    Real tol_rel = 1.e-10;
    Real tol_abs = 0.0;
    
    private:
    
	/*
	Geometry S_geom;
	BoxArray S_ba;
	DistributionMapping S_dmap;
	
	MultiFab solution;
	MultiFab rhs;
	MultiFab alpha_fab;
	MultiFab sigma_fab;
	
	Real A_scalar;
	Real B_scalar;
	
	MultiFab VcrossB;
	*/
};


#endif
