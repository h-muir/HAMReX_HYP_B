#ifndef FUNCDEFS_H_
#define FUNCDEFS_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCRec.H>
#include <string>

#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_EB_LSCore.H>
#include <AMReX_EB_LSCoreBase.H>

#include "classdefs.H"

void main_main ();

// in utility.cpp:
//operator functions 
void printVector(Vector<Real>);
Vector<Real> operator+(Vector<Real> const&, Vector<Real> const&);
Vector<Real> operator-(Vector<Real> const&, Vector<Real> const&);
Vector<Real> operator*(double const&, Vector<Real> const&);

double mag(Vector<Real> const&); 
double dot(Vector<Real> const&, Vector<Real> const&);
Vector<Real> cross(Vector<Real> const&, Vector<Real> const&);

//embedded boundary functions:
void makeEmbeddedBoundary(Geometry&, const int&);
void find_embedded_boundary(const Box&, Array4<Real> const&, Array4<const EBCellFlag>&, const int&, AccessVariable&);
void constructLevelSet(const int&, Geometry&, const BoxArray&, const DistributionMapping&, 
						SettingsStruct const&, ParameterStruct const&, MultiFab&, iMultiFab&);
void initialise_levelset_geometry(MultiFab&, const int, Geometry&, SettingsStruct const&, 
							ParameterStruct const&, AccessVariable&);
void initialise_levelset_normals(MultiFab&, MultiFab&, const Real*, AccessVariable&);

//conversion functions:
void S_conWtoU(Box const&, Array4<Real> const&, AccessVariable&);
void S_conUtoW(Box const&, Array4<Real> const&, AccessVariable&);
void F_conWtoF(Box const&, Array4<Real> const&, AccessVariable&);

//class conversion functions:
ConsU conWtoU(Prim const&); ConsF conWtoF(Prim const&); ConsF conWtoG(Prim const&); 
ConsF conWtoSE(Prim const&); Prim conUtoW(ConsU const&);
Prim rankine_hugoniot(Prim const&, double const&, Vector<Real> const&);

//Plasma19 interpolation and look-up table functions:
Real bilinear_interp(const Vector< Vector<Real> >&, const Real&, const Real&, const int&, 
						const int&, const Real&, const Real&);
Real interp1D(const Vector<Real>&, const Real&, int&, bool&);
						
//euler.cpp functions:
void initialiseStructs(SettingsStruct&, ParameterStruct&);
void setBoundaryConditions(Vector<BCRec>&, const ParameterStruct&, const int&, AccessVariable&);
void special_boundary(Box const&, MultiFab&, SettingsStruct const&, const Real*, AccessVariable&);
void initial_conditions(int const&, Box const&, Array4<Real> const&, const Real* dx, ParameterStruct const&, AccessVariable&);
void initialise_Bfield(Box const&, Array4<Real> const&, const Real*, ParameterStruct const&, AccessVariable&);
void initialise_problem(MultiFab&, const FabArray<EBCellFlagFab>&, Geometry const&, SettingsStruct const&, ParameterStruct const&, AccessVariable&);
void rigid_body_step1(MultiFab&, MultiFab&, const FabArray<EBCellFlagFab>&, const Real*, AccessVariable&);
void rigid_body_step2(MultiFab&, MultiFab&, const FabArray<EBCellFlagFab>&, const Real*, AccessVariable&, double&);
void check_assumed_variables(Vector<std::string>);
void determine_grad_max(const Box&, Array4<Real> const&, const Real*, std::string, std::string const, Real&, AccessVariable&);
void amr_tagging(Array4<char> const&, const Box&, Array4<Real> const&, const Real*, 
					Real, Real, std::string, std::string const, AccessVariable&);
int EoS_properties(std::string,const double&, const double&, const double&, double&, std::string);
double delta_bar_calc(std::string&, double&, double&, double&, double&);
Prim delta(double&, Prim&, Prim&, Prim&, std::string&);

double compute_alpha(const ConsU&);
double compute_eps_x(Array4<Real> const&, const int&, const int&, const int&, AccessVariable&);
double compute_eps_y(Array4<Real> const&, const int&, const int&, const int&, AccessVariable&);
void calc_fluxes_HLLC_HS_x(Box const&, Array4<Real> const&, Array4<Real> const&, 
					Real const, Real const, int const, AccessVariable&, bool);
void calc_fluxes_HLL_x(Box const&, Array4<Real> const&, Array4<Real> const&, 
					Real const, Real const, int const, AccessVariable&, bool);
void calc_fluxes_HLLC_x(Box const&, Array4<Real> const&, Array4<Real> const&, 
					Real const, Real const, int const, AccessVariable&, bool);
void calc_fluxes_HLLC_HS_y(Box const&, Array4<Real> const&, Array4<Real> const&, 
					Real const, Real const, int const, AccessVariable&, bool);
void calc_fluxes_HLL_y(Box const&, Array4<Real> const&, Array4<Real> const&, 
					Real const, Real const, int const, AccessVariable&, bool);
void calc_fluxes_HLLC_y(Box const&, Array4<Real> const&, Array4<Real> const&, 
					Real const, Real const, int const, AccessVariable&, bool);
void cons_update(const std::string, Box const&, Array4<Real> const&, Array4<Real> const&, 
        Array4<Real> const&, Real const, Real const, int const, AccessVariable&);
void radial_sourceterm_integration(MultiFab&, MultiFab&, Real const, Real const, AccessVariable&);
Real max_dim_wave_speed(Box const&, Array4<Real> const&, int const, AccessVariable&);
void rescale_fluxes(const std::string, Box const&, Array4<Real> const&, Real const, 
						Real const, int const, AccessVariable&);					
void Qtransform(Array4<Real> const&, AccessVariable&, const int&, const int&, const int&, const int&, const int&, 
						const double&, const double&, const double&, const double&);
void Qtransform2(Array4<Real> const&, AccessVariable&, const int&, const int&, const int&, const int&, const int&, 
						const double&, const double&, const double&, const double&);
void initialise_boundary_cells(const Box&, Array4<Real> const&, Array4<const EBCellFlag>&, const int&, 
									AccessVariable&, const Real*);
void Qcalc(Array4<Real> const&, const int&, const int&, const int&, AccessVariable&, const Real*);
void sweep(const Box& box, Array4<Real> const&, Array4<const EBCellFlag>&, const int&, 
									AccessVariable&, const Real*, double&);


#endif 
