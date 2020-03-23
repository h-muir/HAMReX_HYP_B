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
	
	EMField(ParameterStruct& p, Geometry const& geom, BoxArray const& grids, DistributionMapping const& dmap);
	void initialise_params(ParameterStruct& p);
	void defineSolution(MultiFab& S_EM, AccessVariable& V);
	void computeAlphaFab();
	void defineScalars();
	void computeFaceSigmaFabs(MultiFab& S_EM, AccessVariable& V);
	void computeRHS(MultiFab& S_EM, AccessVariable& V);
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
    
    int ng = 1;


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
    
	Geometry S_geom;
	BoxArray S_grids;
	DistributionMapping S_dmap;
	
	MultiFab solution;
	
	Real A_scalar;
	Real B_scalar;
	
	
	MultiFab alpha_fab;
	MultiFab sigma_fab;
	Array<MultiFab,AMREX_SPACEDIM> face_sigma_fabs;
	
	MultiFab rhs;
	MultiFab VcrossB;
	MultiFab divVcrossB;
	
};


#endif
