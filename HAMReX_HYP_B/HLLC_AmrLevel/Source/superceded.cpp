

void calc_fluxes_HLLC_SWM_original_x(Box const& box, Array4<Real> const& S_old, Array4<Real> const& F, 
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
	ConsF ALstar(0,v0,0), ARstar(0,v0,0);
	ConsF D_HLL(0,v0,0);
	double epsL, epsR, alpha;
	double SLbar, SRbar;
	//Determined flux: (F_i-1/2):                      
	ConsF Fi(0,v0,0);        
	
	double tol  = 1e-6;
    
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
					FLstarHLLC.mass = FL.mass + SL*(ULstar.rho - UL.rho);
					FLstarHLLC.mom_vec = FL.mom_vec + SL*(ULstar.rhou_vec - UL.rhou_vec);  //check vector addition
					FLstarHLLC.en = FL.en + SL*(ULstar.E - UL.E);
					
					FstarHLL.mass = (SR*FL.mass - SL*FR.mass + SL*SR*(UR.rho-UL.rho))/(SR-SL);
					FstarHLL.mom_vec[0] = (SR*FL.mom_vec[0] - SL*FR.mom_vec[0] + SL*SR*(UR.rhou_vec[0]-UL.rhou_vec[0]))/(SR-SL);
					FstarHLL.mom_vec[1] = (SR*FL.mom_vec[1] - SL*FR.mom_vec[1] + SL*SR*(UR.rhou_vec[1]-UL.rhou_vec[1]))/(SR-SL);
					FstarHLL.en = (SR*FL.en - SL*FR.en + SL*SR*(UR.E-UL.E))/(SR-SL);
					
					ALstar.mass = FLstarHLLC.mass - FstarHLL.mass;
					ALstar.mom_vec = FLstarHLLC.mom_vec - FstarHLL.mom_vec;
					ALstar.en = FLstarHLLC.en - FstarHLL.en;
					
					
					alpha = compute_alpha(UR);
					
					Prim W1(0,v0,0), W2(0,v0,0), W3(0,v0,0), W4(0,v0,0), W5(0,v0,0), W6(0,v0,0);
					v0[0] = S_old(i,j,k,V["u"]);
					v0[1] = S_old(i,j,k,V["v"]);
					W3 = Prim(S_old(i,j,k,V["rho"]), v0, S_old(i,j,k,V["p"]));
					v0[0] = S_old(i-1,j,k,V["u"]);
					v0[1] = S_old(i-1,j,k,V["v"]);
					W4 = Prim(S_old(i-1,j,k,V["rho"]), v0, S_old(i-1,j,k,V["p"]));
					v0[0] = S_old(i,j+1,k,V["u"]);
					v0[1] = S_old(i,j+1,k,V["v"]);
					W1 = Prim(S_old(i,j+1,k,V["rho"]), v0, S_old(i,j+1,k,V["p"]));
					v0[0] = S_old(i-1,j+1,k,V["u"]);
					v0[1] = S_old(i-1,j+1,k,V["v"]);
					W2 = Prim(S_old(i-1,j+1,k,V["rho"]), v0, S_old(i-1,j+1,k,V["p"]));
					v0[0] = S_old(i-1,j-1,k,V["u"]);
					v0[1] = S_old(i-1,j-1,k,V["v"]);
					W6 = Prim(S_old(i-1,j-1,k,V["rho"]), v0, S_old(i-1,j-1,k,V["p"]));
					v0[0] = S_old(i,j-1,k,V["u"]);
					v0[1] = S_old(i,j-1,k,V["v"]);
					W5 = Prim(S_old(i,j-1,k,V["rho"]), v0, S_old(i,j-1,k,V["p"]));
					
					if(S_old(i,j,k,V["level_set"]) > 0){
						epsL = 0;
					}else{
						epsL = compute_eps_x(W1, W2, W3, W4, W5, W6);
					}
					epsR = epsL;
					
					
					SLbar = SL-alpha*epsL;
					SRbar = SR+alpha*epsR;
					
					D_HLL.mass = (fabs(SRbar)-fabs(SLbar))/(2*(SR-SL))*(FL.mass - FR.mass)+\
										(fabs(SLbar)*SR-fabs(SRbar)*SL)/(2*(SR-SL))*(UL.rho-UR.rho);				//Dbar
					D_HLL.mom_vec[0] = (fabs(SR)-fabs(SL))/(2*(SR-SL))*(FL.mom_vec[0] - FR.mom_vec[0])+\
										(fabs(SL)*SR-fabs(SR)*SL)/(2*(SR-SL))*(UL.rhou_vec[0]-UR.rhou_vec[0]);		//D
					D_HLL.mom_vec[1] = (fabs(SRbar)-fabs(SLbar))/(2*(SR-SL))*(FL.mom_vec[1] - FR.mom_vec[1])+\
										(fabs(SLbar)*SR-fabs(SRbar)*SL)/(2*(SR-SL))*(UL.rhou_vec[1]-UR.rhou_vec[1]);//Dbar
					D_HLL.en = (fabs(SR)-fabs(SL))/(2*(SR-SL))*(FL.en - FR.en)+\
										(fabs(SL)*SR-fabs(SR)*SL)/(2*(SR-SL))*(UL.E-UR.E);							//D
					
					FLstar.mass = 0.5*(FL.mass + FR.mass) + D_HLL.mass + ALstar.mass;
					FLstar.mom_vec = 0.5*(FL.mom_vec + FR.mom_vec) + D_HLL.mom_vec + ALstar.mom_vec;
					FLstar.en = 0.5*(FL.en + FR.en) + D_HLL.en + ALstar.en;
					
					/*
					if(epsL > tol){
						Fi = FstarHLL;
					}else{
						Fi = FLstarHLLC;
					}
					*/
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
					FRstarHLLC.mass = FR.mass + SR*(URstar.rho - UR.rho);
					FRstarHLLC.mom_vec = FR.mom_vec + SR*(URstar.rhou_vec - UR.rhou_vec); //check vector addition
					FRstarHLLC.en = FR.en + SR*(URstar.E - UR.E);
					
					FstarHLL.mass = (SR*FL.mass - SL*FR.mass + SL*SR*(UR.rho-UL.rho))/(SR-SL);
					FstarHLL.mom_vec[0] = (SR*FL.mom_vec[0] - SL*FR.mom_vec[0] + SL*SR*(UR.rhou_vec[0]-UL.rhou_vec[0]))/(SR-SL);
					FstarHLL.mom_vec[1] = (SR*FL.mom_vec[1] - SL*FR.mom_vec[1] + SL*SR*(UR.rhou_vec[1]-UL.rhou_vec[1]))/(SR-SL);
					FstarHLL.en = (SR*FL.en - SL*FR.en + SL*SR*(UR.E-UL.E))/(SR-SL);
					
					ARstar.mass = FRstarHLLC.mass - FstarHLL.mass;
					ARstar.mom_vec = FRstarHLLC.mom_vec - FstarHLL.mom_vec;
					ARstar.en = FRstarHLLC.en - FstarHLL.en;
					
					alpha = compute_alpha(UR);
					S_old(i,j,k,V["sigma"]) = alpha;
					
					Prim W1(0,v0,0), W2(0,v0,0), W3(0,v0,0), W4(0,v0,0), W5(0,v0,0), W6(0,v0,0);
					v0[0] = S_old(i,j,k,V["u"]);
					v0[1] = S_old(i,j,k,V["v"]);
					W3 = Prim(S_old(i,j,k,V["rho"]), v0, S_old(i,j,k,V["p"]));
					v0[0] = S_old(i-1,j,k,V["u"]);
					v0[1] = S_old(i-1,j,k,V["v"]);
					W4 = Prim(S_old(i-1,j,k,V["rho"]), v0, S_old(i-1,j,k,V["p"]));
					v0[0] = S_old(i,j+1,k,V["u"]);
					v0[1] = S_old(i,j+1,k,V["v"]);
					W1 = Prim(S_old(i,j+1,k,V["rho"]), v0, S_old(i,j+1,k,V["p"]));
					v0[0] = S_old(i-1,j+1,k,V["u"]);
					v0[1] = S_old(i-1,j+1,k,V["v"]);
					W2 = Prim(S_old(i-1,j+1,k,V["rho"]), v0, S_old(i-1,j+1,k,V["p"]));
					v0[0] = S_old(i-1,j-1,k,V["u"]);
					v0[1] = S_old(i-1,j-1,k,V["v"]);
					W6 = Prim(S_old(i-1,j-1,k,V["rho"]), v0, S_old(i-1,j-1,k,V["p"]));
					v0[0] = S_old(i,j-1,k,V["u"]);
					v0[1] = S_old(i,j-1,k,V["v"]);
					W5 = Prim(S_old(i,j-1,k,V["rho"]), v0, S_old(i,j-1,k,V["p"]));
					
					if(S_old(i,j,k,V["level_set"]) > 0){
						epsL = 0;
					}else{
						epsL = compute_eps_x(W1, W2, W3, W4, W5, W6);
					}
					epsR = epsL;
					
					SLbar = SL-alpha*epsL;
					SRbar = SR+alpha*epsR;
					
					D_HLL.mass = (fabs(SRbar)-fabs(SLbar))/(2*(SR-SL))*(FL.mass - FR.mass)+\
										(fabs(SLbar)*SR-fabs(SRbar)*SL)/(2*(SR-SL))*(UL.rho-UR.rho);				//Dbar
					D_HLL.mom_vec[0] = (fabs(SR)-fabs(SL))/(2*(SR-SL))*(FL.mom_vec[0] - FR.mom_vec[0])+\
										(fabs(SL)*SR-fabs(SR)*SL)/(2*(SR-SL))*(UL.rhou_vec[0]-UR.rhou_vec[0]);		//D
					D_HLL.mom_vec[1] = (fabs(SRbar)-fabs(SLbar))/(2*(SR-SL))*(FL.mom_vec[1] - FR.mom_vec[1])+\
										(fabs(SLbar)*SR-fabs(SRbar)*SL)/(2*(SR-SL))*(UL.rhou_vec[1]-UR.rhou_vec[1]);//Dbar
					D_HLL.en = (fabs(SR)-fabs(SL))/(2*(SR-SL))*(FL.en - FR.en)+\
										(fabs(SL)*SR-fabs(SR)*SL)/(2*(SR-SL))*(UL.E-UR.E);							//D
					

					FRstar.mass = 0.5*(FL.mass + FR.mass) + D_HLL.mass + ARstar.mass;
					FRstar.mom_vec = 0.5*(FL.mom_vec + FR.mom_vec) + D_HLL.mom_vec + ARstar.mom_vec;
					FRstar.en = 0.5*(FL.en + FR.en) + D_HLL.en + ARstar.en;
					
					/*
					if(epsL > tol){
						Fi = FstarHLL;
					}else{
						Fi = FRstarHLLC;
					}
					*/					
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
