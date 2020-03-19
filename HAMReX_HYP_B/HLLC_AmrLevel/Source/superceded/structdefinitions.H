#ifndef STRUCTDEFINITIONS
#define STRUCTDEFINITIONS

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>
#include <AMReX_ParmParse.H>
#include <string>

using namespace amrex;

struct ParameterStruct
{
	Vector<Real> dimL;
    Vector<int>  n_cells;
    Vector<Real> a;
 
	int Ncomp;
	int Nghost;
	
	Real phiL;
	Real phiR;
    Real CFL;
    Real x0;
    Real dx, dy;
    
    std::string IC;
    
    ParameterStruct()
    {
		dimL 	= Vector<Real>(AMREX_SPACEDIM);
		n_cells = Vector<int>(AMREX_SPACEDIM);
		Vector<Real> a(AMREX_SPACEDIM);
	}
};

void main_main ();

void initial_conditions(Box const&, Array4<Real> const&, ParameterStruct const&);

void calc_fluxes(const std::string, Box const&, Array4<Real> const&, Array4<Real> const&, Real const, Real const);

void update(const std::string, Box const&, Array4<Real> const&, Array4<Real> const&, Array4<Real> const&, Real const);

void advance(MultiFab&, MultiFab&, Array<MultiFab, AMREX_SPACEDIM>& flux_arr, Geometry const& geom,
                Real const, ParameterStruct const&);



#endif
