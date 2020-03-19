void calc_fluxes_HLLC_x(Box const& box, Array4<Real> const& F, Array4<Real> const& S_old, 
							int const nc, AccessVariable& V){
	
	const auto lo = lbound(box);
    const auto hi = ubound(box);
    
    double aL, aR;        			//left and right sound speeds
								
	double SL, SR, Sstar, Splus;    //wave speeds used in HLLC
	
	
	Vector<Real> v0(3);
	//4 cons. states:
	ConsU UL(0,v0,0), UR(0,v0,0);
	ConsU ULstar(0,v0,0), URstar(0,v0,0);
	//4 cons. fluxes:   
	ConsF FL(0,v0,0), FR(0,v0,0);
	ConsF FLstar(0,v0,0), FRstar(0,v0,0);
	//Determined flux: (F_i-1/2):                      
	ConsF Fi(0,v0,0);        
	
	//Primitive left/right states:
	
	Prim WL(0,v0,0);
	Prim WR(0,v0,0);   
	
	//note: all the copying here might be really slow ?
	
	/*careful with indexing:
	 
		Fi   Fi+1  ..  Fn+1
	...... ____ ____ ____ .....
	|  G  | S  | S  | S  |  G  |
	|..-1.|__i_|_.._|__n_|..+1.| 
	
	*/
	
	int L, R;
	
	for(int k = lo.z; k <= hi.z; ++k) {
        for(int j = lo.y; j <= hi.y; ++j) {
            for(int i = lo.x; i <= hi.x+1; ++i) { 
				L = i-1;
				R = i;
				//-------------------convert S to W---------------------------
				//assuming P(rho, u, v, p)
				v0[0] = S_old(L,j,k,V["u"]);
				v0[1] = S_old(L,j,k,V["v"]);
				WL = Prim(S_old(L,j,k,V["rho"]), v0, S_old(L,j,k,V["p"]));
				v0[0] = S_old(R,j,k,V["u"]);
				v0[1] = S_old(R,j,k,V["v"]);
				WR = Prim(S_old(R,j,k,V["rho"]), v0, S_old(R,j,k,V["p"]));
				//-------------------convert S to W---------------------------
				
				//direct wave speed estimates:
				aL = EoS_properties(EoS, WL.rho, WL.p, 0, "a", EoS19);
				//cout << "aL = " << aL << endl;
				aR = EoS_properties(EoS, WR.rho, WR.p, 0, "a", EoS19);
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
						amrex::Abort("aborted in calc_fluxes_HLLC_x");
						break;
					}
					//left star state flux:
					FLstar.mass = FL.mass + SL*(ULstar.rho - UL.rho);
					FLstar.mom_vec = FL.mom_vec + SL*(ULstar.rhou_vec - UL.rhou_vec);  //check vector addition
					FLstar.en = FL.en + SL*(ULstar.E - UL.E);
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
						amrex::Abort("aborted in calc_fluxes_HLLC_x");
						break;
					}
					FRstar.mass = FR.mass + SR*(URstar.rho - UR.rho);
					FRstar.mom_vec = FR.mom_vec + SR*(URstar.rhou_vec - UR.rhou_vec); //check vector addition
					FRstar.en = FR.en + SR*(URstar.E - UR.E);
					Fi = FRstar;
				}else if(SR <= 0){
					Fi = FR; 
				}
				
				//-------------------convert F to flux_arr------------------------
				F(i,j,k,V["mass"]) = Fi.mass;
				F(i,j,k,V["mom_x"]) = Fi.mom_vec[0];
				F(i,j,k,V["mom_y"]) = Fi.mom_vec[1];
				F(i,j,k,V["en"]) = Fi.en;
				F(i,j,k,V["rho"]) = 0.0;
				F(i,j,k,V["u"]) = 0.0;
				F(i,j,k,V["v"]) = 0.0;
				F(i,j,k,V["p"]) = 0.0;
				//----------------------------------------------------------------							
			}
		}
	}//end i,j,k loop
	
	//return max_wave_speed;
				         	
} //end HLLC_x 

void calc_fluxes_HLLC_y(Box const& box, Array4<Real> const& G, Array4<Real> const& S_old, 
							int const nc, AccessVariable& V){
	
	Real max_wave_speed = 0.0;
	
	const auto lo = lbound(box);
    const auto hi = ubound(box);
    
    double aL, aR;        			//left and right sound speeds
									//and fast megneto-acoustic speeds
	double SL, SR, Sstar, Splus;    //wave speeds used in HLLC
	
	
	Vector<Real> v0(3);
	//4 cons. states:
	ConsU UL(0,v0,0), UR(0,v0,0);
	ConsU ULstar(0,v0,0), URstar(0,v0,0);
	//4 cons. fluxes:   
	ConsF GL(0,v0,0), GR(0,v0,0);
	ConsF GLstar(0,v0,0), GRstar(0,v0,0);
	//Determined flux: (F_i-1/2):                      
	ConsF Gi(0,v0,0);        
	
	//Primitive left/right states:
	
	Prim WL(0,v0,0);
	Prim WR(0,v0,0);   
	
	//note: all the copying here might be really slow ?
	
	/*careful with indexing:
	 
		Fi   Fi+1  ..  Fn+1
	...... ____ ____ ____ .....
	|  G  | S  | S  | S  |  G  |
	|..-1.|__i_|_.._|__n_|..+1.| 
	
	*/
	
	int L, R;
	
	for(int k = lo.z; k <= hi.z; ++k) {
        for(int j = lo.y; j <= hi.y+1; ++j) {
            for(int i = lo.x; i <= hi.x; ++i) { 
				L = j-1;
				R = j;
				//-------------------convert S to W---------------------------
				//assuming P(rho, u, v, p)
				v0[0] = S_old(i,L,k,V["u"]);
				v0[1] = S_old(i,L,k,V["v"]);
				WL = Prim(S_old(i,L,k,V["rho"]), v0, S_old(i,L,k,V["p"]));
				v0[0] = S_old(i,R,k,V["u"]);
				v0[1] = S_old(i,R,k,V["v"]);
				WR = Prim(S_old(i,R,k,V["rho"]), v0, S_old(i,R,k,V["p"]));
				//-------------------convert S to W---------------------------
				
				//direct wave speed estimates:
				aL = EoS_properties(EoS, WL.rho, WL.p, 0, "a", EoS19);
				//cout << "aL = " << aL << endl;
				aR = EoS_properties(EoS, WR.rho, WR.p, 0, "a", EoS19);
				SL = ((WL.u_vec[1] - aL < WR.u_vec[1] - aR)? WL.u_vec[1] - aL : WR.u_vec[1] - aR);
				SR = ((WL.u_vec[1] + aL > WR.u_vec[1] + aR)? WL.u_vec[1] + aL : WR.u_vec[1] + aR);
				Sstar = (WR.u_vec[1]*WR.rho*(SR-WR.u_vec[1]) - WL.u_vec[1]*WL.rho*(SL-WL.u_vec[1]) + (WL.p-WR.p))\
						/(WR.rho*(SR-WR.u_vec[1])-WL.rho*(SL-WL.u_vec[1]));
				Splus = ((fabs(WL.u_vec[1])+aL > fabs(WR.u_vec[1])+aR)? \
						  fabs(WL.u_vec[1])+aL : fabs(WR.u_vec[1])+aR);
						  
				if(Splus > max_wave_speed){
					max_wave_speed = Splus;
				}
						  
				//W-->U and W-->G conversions:
				UL = conWtoU(WL);
				UR = conWtoU(WR);
				GL = conWtoG(WL);
				GR = conWtoG(WR);
				
				//flux evaluation at (x/t)=0:	
				if( 0 < SL){
					Gi = GL;
				}else if(SL <= 0 && 0 < Sstar){
					//cout << "left star state \n";
					//left star state terms:
					ULstar.rho = WL.rho*((SL-WL.u_vec[1])/(SL-Sstar));
					ULstar.rhou_vec[0] = WL.rho*WL.u_vec[0]*((SL-WL.u_vec[1])/(SL-Sstar));
					ULstar.rhou_vec[1] = ULstar.rho*Sstar; 
					ULstar.E = ((SL-WL.u_vec[1])/(SL-Sstar))*(UL.E+WL.rho*(Sstar-WL.u_vec[1])*\
									(Sstar+WL.p/(WL.rho*(SL-WL.u_vec[1]))));
					if(ULstar.E < 0){
						cout << "negative energy computed in cell" << endl;
						break;
					}
					//left star state flux:
					GLstar.mass = GL.mass + SL*(ULstar.rho - UL.rho);
					GLstar.mom_vec = GL.mom_vec + SL*(ULstar.rhou_vec - UL.rhou_vec);  //check vector addition
					GLstar.en = GL.en + SL*(ULstar.E - UL.E);
					Gi = GLstar;
				}else if(Sstar <= 0 && 0 < SR){
					//cout << "right star state \n";
					//right star state terms:
					URstar.rho = WR.rho*((SR-WR.u_vec[1])/(SR-Sstar));
					URstar.rhou_vec[0] = WR.rho*WR.u_vec[0]*((SR-WR.u_vec[1])/(SR-Sstar));
					URstar.rhou_vec[1] = URstar.rho*Sstar; 
					URstar.E = ((SR-WR.u_vec[1])/(SR-Sstar))*(UR.E + WR.rho*(Sstar-WR.u_vec[1])*\
									(Sstar+WR.p/(WR.rho*(SR-WR.u_vec[1]))));
					if(ULstar.E < 0){
						cout << "negative energy computed in cell" << endl;
						break;
					}
					GRstar.mass = GR.mass + SR*(URstar.rho - UR.rho);
					GRstar.mom_vec = GR.mom_vec + SR*(URstar.rhou_vec - UR.rhou_vec); //check vector addition
					GRstar.en = GR.en + SR*(URstar.E - UR.E);
					Gi = GRstar;
				}else if(SR <= 0){
					Gi = GR; 
				}
				
				//-------------------convert G to flux_arr------------------------
				G(i,j,k,V["mass"]) = Gi.mass;
				G(i,j,k,V["mom_x"]) = Gi.mom_vec[0];
				G(i,j,k,V["mom_y"]) = Gi.mom_vec[1];
				G(i,j,k,V["en"]) = Gi.en;
				G(i,j,k,V["rho"]) = 0.0;
				G(i,j,k,V["u"]) = 0.0;
				G(i,j,k,V["v"]) = 0.0;
				G(i,j,k,V["p"]) = 0.0;
				//----------------------------------------------------------------							
			}
		}
	}//end i,j,k loop
	
	//return max_wave_speed;
				         	
} //end HLLC_y 
