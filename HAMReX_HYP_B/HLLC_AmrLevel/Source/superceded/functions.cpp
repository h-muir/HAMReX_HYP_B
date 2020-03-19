#include "structdefinitions.H"
#include <cmath>

using namespace amrex;


//currently just set-up for x-direction discontinuity
void initial_conditions(Box const& box, Array4<Real> const& prop_arr, ParameterStruct const& p){

    const auto lo = lbound(box);
    const auto hi = ubound(box);
	
	//amrex::Print() << "lo : " << lo << ", hi : " << hi << "\n";

	
    int int_x0 = (p.x0/p.dimL[0])*p.n_cells[0];
    int int_y0 = (p.x0/p.dimL[1])*p.n_cells[1];
	
	
	if(p.IC == "source"){
		Real h, R, rp;
		rp = 0.1;
		R = rp*p.dimL[0];
		for(int k = lo.z; k <= hi.z; ++k) {
			for(int j = lo.y; j <= hi.y; ++j) {
				for(int i = lo.x; i <= hi.x; ++i) {
					h = 100.*(pow(i*p.dx-0.5,2) + pow(j*p.dy-0.5,2));
					//amrex::Print() << "h: " << h << "\n";
					prop_arr(i,j,k) = 1.0 + exp(-h);

				}
			}
		}		
	}else{
		for(int k = lo.z; k <= hi.z; ++k) {
			for(int j = lo.y; j <= hi.y; ++j) {
				for(int i = lo.x; i <= hi.x; ++i) {

					if(i <= int_x0){
						prop_arr(i,j,k) = p.phiL;
					}else{
						prop_arr(i,j,k) = p.phiR;
					}

				}
			}
		}
	}
}

//Lax-Friedriks scheme
void calc_fluxes(const std::string sweep, Box const& box, Array4<Real> const& flux_prop_arr, Array4<Real> const& prop_arr,
                 Real const a, Real const CFL){

    const auto lo = lbound(box);
    const auto hi = ubound(box);

    Real xflux_in, xflux_out;
    Real yflux_in, yflux_out;

    for(int k = lo.z; k <= hi.z; ++k) {
        for(int j = lo.y; j <= hi.y; ++j) {
            for(int i = lo.x; i <= hi.x; ++i) {
                if(sweep == "x-sweep"){
                    //xflux_in located on the node on the LHS of cell centre, for same index i
                    xflux_in = (1.0+CFL)/(2.0*CFL)*a*prop_arr(i-1,j,k) + (CFL-1.0)/(2.0*CFL)*a*prop_arr(i,j,k);
                    flux_prop_arr(i,j,k) = xflux_in;
                    //for the most RHS cell centre there is an additional flux node to the RHS at index i+1:
                    if(i == hi.x){
                        xflux_out = (1.0+CFL)/(2.0*CFL)*a*prop_arr(i,j,k) + (CFL-1.0)/(2.0*CFL)*a*prop_arr(i+1,j,k);
                        flux_prop_arr(i+1,j,k) = xflux_out;
                    }
                }else if(sweep == "y-sweep"){
                    //yflux_in located on the node below cell centre, for same index j
                    yflux_in = (1.0+CFL)/(2.0*CFL)*a*prop_arr(i,j-1,k) + (CFL-1.0)/(2.0*CFL)*a*prop_arr(i,j,k);
                    flux_prop_arr(i,j,k) = yflux_in;
                    //for the most upward cell centre there is an additional flux node to the top at index j+1:
                    if(j == hi.y){
                        yflux_out = (1.0+CFL)/(2.0*CFL)*a*prop_arr(i,j,k) + (CFL-1.0)/(2.0*CFL)*a*prop_arr(i,j+1,k);
                        flux_prop_arr(i,j+1,k) = yflux_out;
                    }


                }

                
            }
        }
    }

    

}


void update(const std::string sweep, Box const& box, Array4<Real> const& flux_prop_arr, Array4<Real> const& prop_arr, 
        Array4<Real> const& prop_arr_new, Real const dtdx){
    
    const auto lo = lbound(box);
    const auto hi = ubound(box);

    for(int k = lo.z; k <= hi.z; ++k) {
        for(int j = lo.y; j <= hi.y; ++j) {
            for(int i = lo.x; i <= hi.x; ++i) {
                if(sweep == "x-sweep"){
                    prop_arr_new(i,j,k) = prop_arr(i,j,k) + dtdx*(flux_prop_arr(i,j,k) - flux_prop_arr(i+1,j,k));
                }else if(sweep == "y-sweep"){
                    prop_arr_new(i,j,k) = prop_arr(i,j,k) + dtdx*(flux_prop_arr(i,j,k) - flux_prop_arr(i,j+1,k));
                }

            }
        }
    }


}



void advance(MultiFab& phi_old, MultiFab& phi_new, Array<MultiFab, AMREX_SPACEDIM>& flux_arr, Geometry const& geom, 
        Real const dtdx,  ParameterStruct const& p)
{

    //amrex::Print() << "advance function called \n";

    //compute each dimensional flux (x-sweep then y-sweep) 
    for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();
        
        //access current array box
        FArrayBox& fab_old = phi_old[mfi];
        FArrayBox& fab_new = phi_new[mfi];
        FArrayBox& flux_fab_x = flux_arr[0][mfi];
        
        //access the array from the array box:
        Array4<Real> const& prop_arr = fab_old.array();
        Array4<Real> const& prop_arr_new = fab_new.array();
        Array4<Real> const& flux_prop_arr_x = flux_fab_x.array();

        calc_fluxes("x-sweep", bx, flux_prop_arr_x, prop_arr, p.a[0], p.CFL);
        update("x-sweep", bx, flux_prop_arr_x, prop_arr, prop_arr_new, dtdx);
        
    }
    
    //Q: no way to avoid second MFIter with dimensional splitting ?

    phi_new.FillBoundary(); //fills ghost cells for neighbouring patches, not including periodic boundaries
    phi_new.FillBoundary(geom.periodicity()); //fills periodic boundary ghost cells
    
    //copy phi_new to phi_old, where phi_old is then operated on again for the y-sweep
    MultiFab::Copy(phi_old, phi_new, 0, 0, p.Ncomp, p.Nghost); // 0,0 --> starting integer in each multifab
    
    //y-sweep
    for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();
        
        //access current array box
        FArrayBox& fab_old = phi_old[mfi];
        FArrayBox& fab_new = phi_new[mfi];
        FArrayBox& flux_fab_y = flux_arr[1][mfi];
        
        //access the array from the array box:
        Array4<Real> const& prop_arr = fab_old.array();
        Array4<Real> const& prop_arr_new = fab_new.array();
        Array4<Real> const& flux_prop_arr_y = flux_fab_y.array();

        calc_fluxes("y-sweep", bx, flux_prop_arr_y, prop_arr, p.a[1], p.CFL);
        update("y-sweep", bx, flux_prop_arr_y, prop_arr, prop_arr_new, dtdx);

    }

}

