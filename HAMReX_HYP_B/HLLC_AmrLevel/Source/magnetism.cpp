#include "magnetism.H"
//#include "classdefs.H"
#include "funcdefs.H"

#include <AMReX.H>
//#include <AMReX_MLMG.H>
#include <AMReX_Amr.H>
//#include <AMReX_Interpolater.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_BCRec.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_BCUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MLABecLaplacian.H>
//MLNodeLaplacian solves `div dot (sigma grad phi) = rhs` using the nodal solver

using namespace amrex;



EMField::EMField(ParameterStruct& p, Geometry const& geom, AccessVariable& V){
	
	initialise_params(p);
	//initialise_fabs(S, geom, V);
	
}
	
void EMField::initialise_params(ParameterStruct& p){
	
	max_level = p.max_level;
    ref_ratio = 2;
    n_cells.resize(AMREX_SPACEDIM);
    n_cells = p.n_cells;
    max_grid_size = p.max_grid_size;
    
}

void EMField::initialise_fabs(MultiFab& S_EM, Geometry const& geom, AccessVariable& V){
    /*
    S_geom = geom;
    S_ba = S_EM.boxArray();
    S_dmap = S_EM.DistributionMap();
    
    rhs.define(S_EM.boxArray(), S_EM.DistributionMap(), 1, S_EM.nGrow());
    solution.define(S_EM.boxArray(), S_EM.DistributionMap(), 1, S_EM.nGrow());
    alpha_fab.define(S_EM.boxArray(), S_EM.DistributionMap(), 1, S_EM.nGrow());
    sigma_fab.define(S_EM.boxArray(), S_EM.DistributionMap(), 1, S_EM.nGrow());
    
    alpha_fab.setVal(1.0);
    MultiFab::Copy(sigma_fab, S_EM, V["sigma"], 0, 1, S_EM.nGrow());
    MultiFab::Copy(solution, S_EM, V["phi"], 0, 1, S_EM.nGrow());
    //solution.setVal(1.0);
    A_scalar = 0;
    B_scalar = 1;
    
    int lo_bc[BL_SPACEDIM];
    int hi_bc[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; ++i) {
		lo_bc[i] = BCType::foextrap;
		hi_bc[i] = BCType::foextrap;
	}
	Vector<BCRec> bc_vec(3);
	for(int n = 0; n<3; ++n){
		BCRec bc(lo_bc, hi_bc);
		bc_vec[n] = bc;
	}
    
    VcrossB.define(S_EM.boxArray(), S_EM.DistributionMap(), 3, S_EM.nGrow());
    computeVcrossB(S_EM, VcrossB, V);
    VcrossB.FillBoundary();
    FillDomainBoundary(VcrossB, geom, bc_vec);
    //MultiFab::Copy(S_EM, VcrossB, 0, V["J_r"], 3, S_EM.nGrow());
    MultiFab divVcrossB(S_EM.boxArray(), S_EM.DistributionMap(), 1, 0);
    compute_divVcrossB(geom, sigma_fab, VcrossB, divVcrossB);
    //MultiFab::Copy(S_EM, divVcrossB, 0, V["phi"], 1, 0);
    //MultiFab::Copy(rhs, VcrossB, 2, 0, 1, S_EM.nGrow());
    MultiFab::Copy(rhs, divVcrossB, 0, 0, 1, 0);
    rhs.FillBoundary();
    //FillDomainBoundary(rhs, geom, bc_vec[0]);
    //rhs.setVal(0.0);
    */
}

void EMField::computeVcrossB(MultiFab& S_EM, MultiFab &S_out, AccessVariable V){
	
	//S_out has nc = 3 (i,j,k)
	
	Vector<Real> out_vec(3,0.0);
	Vector<Real> V_vec(3,0.0);
	Vector<Real> B_vec(3,0.0);
	
	for (MFIter mfi(S_out, true); mfi.isValid(); ++mfi){
		
	    const Box& bx = mfi.tilebox();	
		const auto lo = lbound(bx); 
		const auto hi = ubound(bx);
	
		Array4<Real> const& S_out_arr = S_out[mfi].array();
		Array4<Real> const& S_EM_arr = S_EM[mfi].array();
		
		for(int k = lo.z; k <= hi.z; ++k) {
			for(int j = lo.y; j <= hi.y; ++j) {
				for(int i = lo.x; i <= hi.x; ++i) {
					
					V_vec[0] = S_EM_arr(i,j,k,V["u"]);
					V_vec[1] = S_EM_arr(i,j,k,V["v"]);
					B_vec[0] = S_EM_arr(i,j,k,V["B_r"]);
					B_vec[1] = S_EM_arr(i,j,k,V["B_z"]);
					B_vec[2] = S_EM_arr(i,j,k,V["B_theta"]);
					
					out_vec = cross(V_vec, B_vec);
					
					S_out_arr(i,j,k,0) = out_vec[0];
					S_out_arr(i,j,k,1) = out_vec[1];
					S_out_arr(i,j,k,2) = out_vec[2];
					
				}
			}
		}//end i,j,k
		
	}//end MFIter
		
}

void EMField::compute_divVcrossB(Geometry const& geom, MultiFab& sigma_fab, MultiFab& S_VB, MultiFab &S_out){
	
	//S_out has nc = 1, ng = 0
	//S_VB has nc = 3 (i,j,k), ng = 1
	
	const Real* dx  = geom.CellSize();
	
	for (MFIter mfi(S_out, true); mfi.isValid(); ++mfi){
		
	    const Box& bx = mfi.tilebox();	
		const auto lo = lbound(bx); 
		const auto hi = ubound(bx);
	
		Array4<Real> const& S_out_arr = S_out[mfi].array();
		Array4<Real> const& S_VB_arr = S_VB[mfi].array();
		Array4<Real> const& sigma_fab_arr = sigma_fab[mfi].array();
		
		for(int k = lo.z; k <= hi.z; ++k) {
			for(int j = lo.y; j <= hi.y; ++j) {
				for(int i = lo.x; i <= hi.x; ++i) {
						
					S_out_arr(i,j,k) = sigma_fab_arr(i,j,k)*((S_VB_arr(i+1,j,k,2)-S_VB_arr(i-1,j,k,2))/(2*dx[0]) + \
										(S_VB_arr(i,j+1,k,2)-S_VB_arr(i,j-1,k,2))/(2*dx[1]));
					
				}
			}
		}//end i,j,k
		
	}//end MFIter
		
}

void EMField::compute_gradPhi(Geometry const& geom, MultiFab& S_phi, MultiFab &S_out){
	
	//S_out has nc = 3, ng = 0
	//S_phi has nc = 1 (i,j,k), ng = 1
	
	const Real* dx  = geom.CellSize();
	
	for (MFIter mfi(S_out, true); mfi.isValid(); ++mfi){
		
	    const Box& bx = mfi.tilebox();	
		const auto lo = lbound(bx); 
		const auto hi = ubound(bx);
	
		Array4<Real> const& S_out_arr = S_out[mfi].array();
		Array4<Real> const& S_phi_arr = S_phi[mfi].array();
		
		for(int k = lo.z; k <= hi.z; ++k) {
			for(int j = lo.y; j <= hi.y; ++j) {
				for(int i = lo.x; i <= hi.x; ++i) {
						
					S_out_arr(i,j,k,0) = (S_phi_arr(i+1,j,k)-S_phi_arr(i-1,j,k))/(2*dx[0]); 
					S_out_arr(i,j,k,1) = (S_phi_arr(i,j+1,k)-S_phi_arr(i,j-1,k))/(2*dx[1]);
					S_out_arr(i,j,k,2) =  0.0;
					
				}
			}
		}//end i,j,k
		
	}//end MFIter
		
}

void EMField::computeJfield(MultiFab& S_EM, MultiFab& S_gradPhi, MultiFab& S_VB, MultiFab& sigma_fab, AccessVariable V){
	
	//S_gradPhi and S_VB have nc = 3, ng = 0
	//sigma_fab has nc = 1, ng = 1
	//S_EM has all components and ng = 1
	
	Vector<Real> J_vec(3,0.0);
	
	for (MFIter mfi(S_VB, true); mfi.isValid(); ++mfi){
		
	    const Box& bx = mfi.tilebox();	
		const auto lo = lbound(bx); 
		const auto hi = ubound(bx);
		
		Array4<Real> const& S_EM_arr = S_EM[mfi].array();
		Array4<Real> const& S_phi_arr = S_gradPhi[mfi].array();
		Array4<Real> const& S_VB_arr = S_VB[mfi].array();
		Array4<Real> const& sigma_fab_arr = sigma_fab[mfi].array();	
		
		for(int k = lo.z; k <= hi.z; ++k) {
			for(int j = lo.y; j <= hi.y; ++j) {
				for(int i = lo.x; i <= hi.x; ++i) {
					
					S_EM_arr(i,j,k,V["J_r"]) = sigma_fab_arr(i,j,k)*(-S_phi_arr(i,j,k,0) + S_VB_arr(i,j,k,0));
					S_EM_arr(i,j,k,V["J_z"]) = sigma_fab_arr(i,j,k)*(-S_phi_arr(i,j,k,1) + S_VB_arr(i,j,k,1));
					S_EM_arr(i,j,k,V["J_theta"]) = sigma_fab_arr(i,j,k)*(-S_phi_arr(i,j,k,2) + S_VB_arr(i,j,k,2));
					J_vec[0] = S_EM_arr(i,j,k,V["J_r"]);
					J_vec[1] = S_EM_arr(i,j,k,V["J_z"]);
					J_vec[2] = S_EM_arr(i,j,k,V["J_theta"]);
					S_EM_arr(i,j,k,V["J_abs"]) = mag(J_vec);
					
				}
			}
		}//end i,j,k
		
	}//end MFIter
		
}

void EMField::set_Derichlet_BC(Geometry const& geom, MultiFab &S_border, int boundary_ref, int ng, Real val){
	
	//dom is whole problem domain, as a box; 
	const Box& dom = geom.Domain();
    const auto prob_lo = lbound(dom); 
    const auto prob_hi = ubound(dom);
	
	for (MFIter mfi(S_border, true); mfi.isValid(); ++mfi){
		
	    const Box& bx = mfi.tilebox();	
		const auto lo = lbound(bx); //lo bound of current box in dom (including ghost cells)
		const auto hi = ubound(bx);
		
		if(ng < 1){
			amrex::Abort("ghost cells assumed >= 1 in EMField::set_Derichlet_BC function");
		}
		
		//only set up for 2D boundary domain, with int_ref's: [0,1,3,4] = [lo.x, lo,y, hi.x, hi.y]
		if(boundary_ref == 0){
			
			Array4<Real> const& S_arr = S_border[mfi].array();	
			
			for(int k = lo.z; k <= hi.z; ++k) {
				for(int j = lo.y; j <= hi.y; ++j) {
					for(int i = prob_hi.x; i <= prob_hi.x+(ng-1); ++i) {
						
						S_arr(i,j,k) = val;
						
					}
				}
			}
			
		}else if(boundary_ref == 1){
			
			Array4<Real> const& S_arr = S_border[mfi].array();	
			
			for(int k = lo.z; k <= hi.z; ++k) {
				for(int j = prob_lo.y; j <= prob_lo.y+(ng-1); ++j) {
					for(int i = lo.x; i <= hi.x; ++i) {
						
						S_arr(i,j,k) = val;
						
					}
				}
			}	
			
		}else if((boundary_ref == 3) && (hi.x == prob_hi.x)){
			
			Array4<Real> const& S_arr = S_border[mfi].array();	
			
			for(int k = lo.z; k <= hi.z; ++k) {
				for(int j = lo.y; j <= hi.y; ++j) {
					for(int i = prob_hi.x-(ng-1); i <= prob_hi.x; ++i) {
						
						S_arr(i,j,k) = val;
						
					}
				}
			}	
			
		}else if((boundary_ref == 4) && (hi.y == prob_hi.y)){
			
			Array4<Real> const& S_arr = S_border[mfi].array();	
			
			for(int k = lo.z; k <= hi.z; ++k) {
				for(int j = prob_hi.y-(ng-1); j <= prob_hi.y; ++j) {
					for(int i = lo.x; i <= hi.x; ++i) {
						
						S_arr(i,j,k) = val;
						
					}
				}
			}	
			
		}//end boundary cases
	}//end MFIter
			
}

void EMField::solve_by_level(MLABecLaplacian &mlabec, MultiFab &S_EM, MultiFab &Sol_crse, int const& ilev, AccessVariable& V){
	
	//ilev is current level for solve
	/*
	LPInfo info;
    info.setAgglomeration(agglomeration);
    info.setConsolidation(consolidation);
    info.setMaxCoarseningLevel(max_coarsening_level);
	
	MLABecLaplacian mlabec({S_geom}, {S_ba}, {S_dmap}, info);
	
	mlabec.setMaxOrder(linop_maxorder);
	
	// This is a 3d problem with homogeneous Neumann BC
	mlabec.setDomainBC({AMREX_D_DECL(LinOpBCType::Neumann,   //x-lo
									 LinOpBCType::Neumann,   //y-lo
									 LinOpBCType::Neumann)}, //z-lo
					   {AMREX_D_DECL(LinOpBCType::Dirichlet,   //x-hi
									 LinOpBCType::Dirichlet,   //y-hi
									 LinOpBCType::Neumann)});//z-hi
	
									 
	if (ilev > 0) {
		mlabec.setCoarseFineBC(&Sol_crse, ref_ratio); 
	}
	
	// for problem with pure homogeneous Neumann BC, we could pass a nullptr
    //mlabec.setLevelBC(0, nullptr);
    
    //boundary coreection for derichlet BC's: 
    //arguments: (multifab, boundary_ref, nGhost, derichlet val)
    //set_Derichlet_BC(solution, 3, 1, 0.0); //x-hi
    //set_Derichlet_BC(solution, 4, 1, 0.0); //y-hi
    
	mlabec.setLevelBC(0, &solution);
	
	mlabec.setScalars(A_scalar, B_scalar);

    mlabec.setACoeffs(0, alpha_fab);
	
	Array<MultiFab,AMREX_SPACEDIM> face_sigma_fab;
	
	for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
	{
		const BoxArray& ba = amrex::convert(sigma_fab.boxArray(),
											IntVect::TheDimensionVector(idim));
		face_sigma_fab[idim].define(ba, sigma_fab.DistributionMap(), 1, 0);
		
	}
	
	amrex::average_cellcenter_to_face(GetArrOfPtrs(face_sigma_fab), sigma_fab, S_geom);
	
	mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(face_sigma_fab));
	
	MLMG mlmg(mlabec);
	mlmg.setMaxIter(max_iter);
	mlmg.setMaxFmgIter(max_fmg_iter);
	mlmg.setVerbose(verbose);
	mlmg.setBottomVerbose(bottom_verbose);
	
	mlmg.solve({&solution}, {&rhs}, tol_rel, tol_abs);
	
	MultiFab::Copy(S_EM, solution, 0, V["phi"], 1, S_EM.nGrow()); 		//fills S_EM with phi solution
	MultiFab gradPhi(S_EM.boxArray(), S_EM.DistributionMap(), 3, 0);
	compute_gradPhi(solution, gradPhi);
	computeJfield(S_EM, gradPhi, VcrossB, V); 							//fills S_EM with J_R, J_z, J_theta, J_abs solution
	*/
}
