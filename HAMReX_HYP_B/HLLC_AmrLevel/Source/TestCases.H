//Module: Toro Test input cases
#ifndef TESTCASES_H_
#define TESTCASES_H_


//zero-default data types for initialisation purposes
Vector<Real> v0(3,0);
#include "classdefs.H"

class TestCase{
    public:
        int n;          //x spatial discretisation
        int m;			//y spatial discretisation
        double R, Z;    //spatial Length and height
        double Tf;      //final time
        double T0;      //initial time
        double CFL;     //Courant number
        double dr, dz;  //spatial step sizes
        double r0, z0;  //property discontinuity positions
        double rc, zc;	//additional parameters for geometry construction
        double rl, zl;
        int nG;         //number of ghost cells on each side of domain
        string construction;
        string geometry;  		//option for rigid body geometry
        bool rigid_bodies;
        Vector<Real> geom_vec; 		//assigned after geomtery is defined
        vector< vector<double> > geom_grads;
        vector< vector<double> > signed_grads;
        Prim WL;
        Prim WR;
        Vector<Real> radial_vals;
        Vector<Real> z_vals;
        string gpname, datfile, title, giffile;
        string exact_sol_file;
        string BCs;
        string EoS; 	//"ideal" or "EoS19" 
        double material_gamma;
        bool EM_source;
        bool cylindrical_coords;
        int nT;         //number of time steps,
                        //not assigned until after solution and number of time steps are known
        TestCase(string CASE):
            WL(0,v0,0), 
            WR(0,v0,0) {
                exact_sol_file = "none";
                material_gamma = 1.4;		//gamma = 1.4 unless otherwise specific in CASE
                EM_source = false;			//EM source terms are off unless otherwise specified
                rigid_bodies = false;
                cylindrical_coords = true;
                if(CASE == "Test1"){
                    title = "torotest1";
                    datfile = "torotest1.dat";
                    gpname = "./visualisation/torotest1.gp";
                    giffile = "./visualisation/torotest1.gif";
                    construction = "default_x";
                    BCs = "transmissive";
                    exact_sol_file = "./ExactSols/Toro411";
                    WL.rho = 1.0; WL.u_vec[0] = 0.0; WL.p = 1.0; 
                    WR.rho = 0.125; WR.u_vec[0] = 0.0; WR.p = 0.1; 
                    n = 80;
                    m = 80;
                    R = 1.0;    //space domain 0 < x < 1
                    Z = 1.0;
                    Tf = 0.25;
                    T0 = 0.0;
                    CFL = 0.9;                    
                    nG = 2;
                    dr = R/(n-2*nG-1);
                    dz = Z/(m-2*nG-1);
                    r0 = 0.5;
                    z0 = 0.5;
                    EoS = "ideal";
                    cylindrical_coords = false;
                }
                if(CASE == "Test2"){
                    title = "torotest2";
                    datfile = "torotest2.dat";
                    gpname = "./visualisation/torotest2.gp";
                    giffile = "./visualisation/torotest2.gif";
                    construction = "default_x";
                    exact_sol_file = "./ExactSols/Toro412";
                    WL.rho = 1.0; WL.u_vec[0] = -2.0; WL.p = 0.4;
                    WR.rho = 1.0; WR.u_vec[0] = 2.0; WR.p = 0.4;
                    n = 80;
                    m = 80;
                    R = 1.0;    //space domain 0 < x < 1
                    Z = 1.0;
                    Tf = 0.15;
                    T0 = 0;
                    CFL = 0.9;
                    nG = 2;
                    dr = R/(n-2*nG-1);
                    dz = Z/(m-2*nG-1);
                    r0 = 0.5;
                    z0 = 0.5;
                    BCs = "transmissive";
                    EoS = "ideal";
                    cylindrical_coords = false;
                }
                
                if(CASE == "Test3"){
                    title = "torotest3";
                    datfile = "torotest3.dat";
                    gpname = "./visualisation/torotest3.gp";
                    giffile = "./visualisation/torotest3.gif";
                    construction = "default_x";
                    exact_sol_file = "./ExactSols/Toro413";
                    WL.rho = 1.0; WL.u_vec[0] = 0.0; WL.p = 1000.0;
                    WR.rho = 1.0; WR.u_vec[0] = 0.0; WR.p = 0.01;
                    n = 100;
                    m = 100;
                    R = 1.0;    //space domain 0 < x < 1
                    Z = 1.0;
                    Tf = 0.012;
                    T0 = 0;
                    CFL = 0.9;
                    nG = 2;
                    dr = R/(n-2*nG-1);
                    dz = Z/(m-2*nG-1);
                    r0 = 0.5;
                    z0 = 0.5;
                    BCs = "transmissive";
                    EoS = "ideal";
                    cylindrical_coords = false;
                }
                /*
                if(CASE == "Test4"){
                    title = "torotest4";
                    datfile = "torotest4.dat";
                    gpname = "gifplot_torotest4.gp";
                    giffile = "./visualisation/torotest4.gif";
                    construction = "default";
                    exact_sol_file = "./ExactSols/Toro414";
                    WL.rho = 1.0; WL.u_vec[0] = 0.0; WL.p = 0.01;
                    WR.rho = 1.0; WR.u_vec[0] = 0.0; WR.p = 100.0;
                    n = 200;
                    R = 1.0;    //space domain 0 < x < 1
                    Z = 1.0;
                    Tf = 0.035;
                    T0 = 0;
                    CFL = 0.9;
                    nG = 2;
                    dr = R/(n-2*nG-1);
                    dz = Z/(m-2*nG-1);
                    r0 = 0.5;
                }
                if(CASE == "Test5"){
                    title = "torotest5";
                    datfile = "torotest5.dat";
                    gpname = "gifplot_torotest5.gp";
                    giffile = "./visualisation/torotest5.gif";
                    construction = "default";
                    exact_sol_file = "./ExactSols/Toro415";
                    WL.rho = 5.99924; WL.u_vec[0] = 19.59745; WL.p = 460.894;
                    WR.rho = 5.99242; WR.u_vec[0] = -6.19633; WR.p = 46.0950;
                    n = 200;
                    R = 1.0;    //space domain 0 < x < 1
                    Z = 1.0;
                    Tf = 0.035;
                    T0 = 0;
                    CFL = 0.9;
                    nG = 2;
                    dr = R/(n-2*nG-1);
                    dz = Z/(m-2*nG-1);
                    r0 = 0.5;
                }
				*/
				
				if(CASE == "cylindrical_explosion"){
                    title = "cylindrical_explosion";
                    datfile = "cylindrical_explosion.dat";
                    gpname = "./visualisation/cylindrical_explosion.gp";
                    giffile = "./visualisation/cylindrical_explosion.gif";
                    construction = "default_x";
                    BCs = "transmissive";
                    exact_sol_file = "./ExactSols/cylindrical_explosion";
                    WL.rho = 1.0; WL.u_vec[0] = 0.0; WL.p = 1.0; 
                    WR.rho = 0.125; WR.u_vec[0] = 0.0; WR.p = 0.1; 
                    n = 150;
                    m = 150;
                    R = 1.0;    //space domain 0 < x < 1
                    Z = 1.0;
                    Tf = 0.25;
                    T0 = 0.0;
                    CFL = 0.9;                    
                    nG = 2;
                    dr = R/(n-2*nG-1);
                    dz = Z/(m-2*nG-1);
                    r0 = 0.4;
                    z0 = 0.4;
                    EoS = "ideal";
                }
                if(CASE == "real_cylindrical_explosion"){
                    title = "real_cylindrical_explosion";
                    datfile = "real_cylindrical_explosion.dat";
                    gpname = "./visualisation/real_cylindrical_explosion.gp";
                    giffile = "./visualisation/real_cylindrical_explosion.gif";
                    construction = "default_x";
                    BCs = "transmissive";
                    exact_sol_file = "none";
                    WL.rho = 1.0; WL.u_vec[0] = 0.0; WL.p = 1.0*101.325e3; 
                    WR.rho = 0.125; WR.u_vec[0] = 0.0; WR.p = 0.1*101.325e3; 
                    n = 100;
                    m = 100;
                    R = 1.0;    //space domain 0 < x < 1
                    Z = 1.0;
                    Tf = 0.6e-3;
                    T0 = 0.0;
                    CFL = 0.9;                    
                    nG = 2;
                    dr = R/(n-2*nG-1);
                    dz = Z/(m-2*nG-1);
                    r0 = 0.4;
                    z0 = 0.4;
                    EoS = "ideal";
                }
                if(CASE == "real_cylindrical_explosion_EoS19"){
                    title = "real_cylindrical_explosion_EoS19";
                    datfile = "real_cylindrical_explosion_EoS19.dat";
                    gpname = "./visualisation/real_cylindrical_explosion_EoS19.gp";
                    giffile = "./visualisation/real_cylindrical_explosion_EoS19.gif";
                    construction = "default_x";
                    BCs = "transmissive";
                    exact_sol_file = "none";
                    WL.rho = 1.0; WL.u_vec[0] = 0.0; WL.p = 1.0*101.325e3; 
                    WR.rho = 0.125; WR.u_vec[0] = 0.0; WR.p = 0.1*101.325e3; 
                    n = 100;
                    m = 100;
                    R = 1.0;    //space domain 0 < x < 1
                    Z = 1.0;
                    Tf = 0.6e-3;
                    T0 = 0.0;
                    CFL = 0.9;                    
                    nG = 2;
                    dr = R/(n-2*nG-1);
                    dz = Z/(m-2*nG-1);
                    r0 = 0.4;
                    z0 = 0.4;
                    EoS = "EoS19";
                }
                if(CASE == "lightning"){
                    title = "lightning";
                    datfile = "lightning.dat";
                    gpname = "./visualisation/lightning.gp";
                    giffile = "./visualisation/lightning.gif";
                    construction = "2Delectrode"; //"default_x"; //
                    geometry = "quadrilateral";
                    BCs = "transmissive";
                    exact_sol_file = "none";
                    WL.rho = 0.9; WL.u_vec[0] = 100.0; WL.p = 40*101.325e3; 
                    WR.rho = 1.225; WR.u_vec[0] = 0.0; WR.p = 1*101.325e3; 
                    n = 104;
                    m = 104;
                    R = 0.20;    //space domain 0 < x < 1
                    Z = 0.20;
                    Tf = 10.0e-6; //48e-6;
                    T0 = 0.0;
                    CFL = 0.9;                    
                    nG = 2;
                    dr = R/(n-2*nG-1);
                    dz = Z/(m-2*nG-1);
                    r0 = 0.016;
                    z0 = 0.20;
                    rc = 0.01-0.5*nG*dr;
                    zc = 0.15+0.5*nG*dz;
                    rl = 0.02+nG*dr;		//radial length for rigid body
                    zl = 0.10+nG*dr;		//z-length
                    EoS = "EoS19";
                    EM_source = true;
                    rigid_bodies = true;
                    cylindrical_coords = true;
                }
                if(CASE == "electrode"){
                    title = "electrode";
                    datfile = "electrode.dat";
                    gpname = "./visualisation/electrode.gp";
                    giffile = "./visualisation/electrode.gif";
                    construction = "2Delectrode"; //"default_x"; //
                    geometry = "quadrilateral";
                    BCs = "transmissive";
                    exact_sol_file = "none";
                    WL.rho = 1.225; WL.u_vec[0] = 100.0; WL.p = 40*101.325e3; 
                    WR.rho = 1.225; WR.u_vec[0] = 0.0; WR.p = 1*101.325e3; 
                    n = 104;
                    m = 104;
                    R = 0.20;    //space domain 0 < x < 1
                    Z = 0.20;
                    Tf = 0.0e-6; //48e-6;
                    T0 = 0.0;
                    CFL = 0.9;                    
                    nG = 2;
                    dr = R/(n-2*nG-1);
                    dz = Z/(m-2*nG-1);
                    r0 = 0.016;
                    z0 = 0.20;
                    rc = 0.008-0.5*nG*dr;
                    zc = 0.15+0.5*nG*dz;
                    rl = 0.016+nG*dr;		//radial length for rigid body
                    zl = 0.10+nG*dr;		//z-length
                    EoS = "ideal";
                    EM_source = false;
                    rigid_bodies = true;
                    cylindrical_coords = true;
                }
                
                if(CASE == "spherical"){
                    title = "spherical";
                    datfile = "spherical.dat";
                    gpname = "./visualisation/spherical.gp";
                    giffile = "./visualisation/spherical.gif";
                    construction = "default_y"; //
                    geometry = "spherical";
                    BCs = "transmissive";
                    exact_sol_file = "none";
                    WL.rho = 0.195e-3; WL.u_vec[0] = 0; WL.u_vec[1] = -5000; WL.p = 23.6; //downstream of shock, v=-5562.0
                    WR.rho = 0.195e-3; WR.u_vec[0] = 0; WR.u_vec[1] = -5000; WR.p = 23.6;  //upstream of shock
                    n = 104;
                    m = 104;
                    R = 1.0;    //space domain 0 < x < 1
                    Z = 1.0;
                    Tf = 50.0e-5; //48e-6;
                    T0 = 0.0;
                    CFL = 0.9;                    
                    nG = 2;
                    dr = R/(n-2*nG-1);
                    dz = Z/(m-2*nG-1);
                    r0 = 0.40;		//spherical nose cone radius
                    z0 = 0.50; 		//height of initial shock wave
                    rc = 0.00;
                    zc = -0.10; 		//spherical centre
                    EoS = "ideal";
                    EM_source = false;
                    rigid_bodies = true;
                    cylindrical_coords = false;
              
                }
                if(CASE == "OtsuAbe"){
                    title = "OtsuAbe";
                    datfile = "OtsuAbe.dat";
                    gpname = "./visualisation/OtsuAbe.gp";
                    giffile = "./visualisation/OtsuAbe.gif";
                    construction = "default_y"; //
                    geometry = "spherical";
                    BCs = "transmissive";
                    exact_sol_file = "none";
                    //US standard atmosphere at 69km altitude
                    WL.rho = 0.8568e-4; WL.u_vec[0] = 0; WL.u_vec[1] = -6500; WL.p = 5.417;  //downstream of shock, v=-5562.0
                    WR.rho = 0.8568e-4; WR.u_vec[0] = 0; WR.u_vec[1] = -6500; WR.p = 5.417;  //upstream of shock
                    n = 104;
                    m = 104;
                    R = 2.0;    //space domain 0 < x < 1
                    Z = 2.0;
                    Tf = 50.0e-5; //48e-6;
                    T0 = 0.0;
                    CFL = 0.9;                    
                    nG = 2;
                    dr = R/(n-2*nG-1);
                    dz = Z/(m-2*nG-1);
                    r0 = 1.0;		//spherical nose cone radius
                    z0 = 0.50; 		//height of initial shock wave
                    rc = 0.00;
                    zc = -0.20; 		//spherical centre
                    EoS = "ideal";
                    EM_source = false;
                    rigid_bodies = true;
                    cylindrical_coords = true;
              
                }
                
                if(CASE == "OtsuAbe19"){
                    title = "OtsuAbe19";
                    datfile = title + ".dat";
                    gpname = "./visualisation/" + title + ".gp";
                    giffile = "./visualisation/" + title + ".gif";
                    construction = "default_y"; //
                    geometry = "spherical";
                    BCs = "transmissive";
                    exact_sol_file = "none";
                    //US standard atmosphere at 69km altitude
                    WL.rho = 0.86e-3; WL.u_vec[0] = 0; WL.u_vec[1] = -6500; WL.p = 67.0;  //downstream of shock, v=-5562.0
                    WR.rho = 0.86e-3; WR.u_vec[0] = 0; WR.u_vec[1] = -6500; WR.p = 67.0;  //upstream of shock
                    n = 104;
                    m = 104;
                    R = 2.0;    //space domain 0 < x < 1
                    Z = 2.0;
                    Tf = 10.0e-5; //48e-6;
                    T0 = 0.0;
                    CFL = 0.9;                    
                    nG = 2;
                    dr = R/(n-2*nG-1);
                    dz = Z/(m-2*nG-1);
                    r0 = 1.0;		//spherical nose cone radius
                    z0 = 0.50; 		//height of initial shock wave
                    rc = 0.00;
                    zc = -0.20; 		//spherical centre
                    EoS = "EoS19";
                    EM_source = false;
                    rigid_bodies = true;
                    cylindrical_coords = true;
              
                }
                
                if(CASE == "rigid_validation1"){
                    title = "rigid_validation1";
                    datfile = "rigid_validation1.dat";
                    gpname = "./visualisation/rigid_validation1.gp";
                    giffile = "./visualisation/rigid_validation1.gif";
                    construction = "default_x"; //
                    geometry = "spherical";
                    BCs = "transmissive";
                    exact_sol_file = "none";
                    WL.rho = 1.3764; WL.u_vec[0] = 0.394; WL.u_vec[1] = 0.0; WL.p = 1.5698; //down or left
                    WR.rho = 1.0; WR.u_vec[0] = 0.0; WR.u_vec[1] = 0.0; WR.p = 1.0;  //up or right
                    n = 204;
                    m = 204;
                    R = 1.0;    //space domain 0 < x < 1
                    Z = 1.0;
                    Tf = 0.40; 
                    T0 = 0.0;
                    CFL = 0.9;                    
                    nG = 2;
                    dr = R/(n-2*nG-1);
                    dz = Z/(m-2*nG-1);
                    r0 = 0.2;		//spherical nose cone radius
                    z0 = 0.2; 		//height of initial shock wave
                    rc = 0.6;
                    zc = 0.5; 		//spherical centre
                    EoS = "ideal";
                    EM_source = false;
                    rigid_bodies = true;
                    cylindrical_coords = false;
                }
                
                if(CASE == "rigid_validation2"){
                    title = "rigid_validation2";
                    datfile = "rigid_validation2.dat";
                    gpname = "./visualisation/rigid_validation2.gp";
                    giffile = "./visualisation/rigid_validation2.gif";
                    construction = "default_x"; //
                    geometry = "spherical";
                    BCs = "transmissive";
                    exact_sol_file = "none";
                    WL.rho = 1.0; WL.u_vec[0] = 1.0; WL.u_vec[1] = 0.0; WL.p = 1.0; //down or left
                    WR.rho = 1.0; WR.u_vec[0] = 1.0; WR.u_vec[1] = 0.0; WR.p = 1.0;  //up or right
                    n = 204;
                    m = 204;
                    R = 1.0;    //space domain 0 < x < 1
                    Z = 1.0;
                    Tf = 1.0; 
                    T0 = 0.0;
                    CFL = 0.9;                    
                    nG = 2;
                    dr = R/(n-2*nG-1);
                    dz = Z/(m-2*nG-1);
                    r0 = 0.2;		//spherical nose cone radius
                    z0 = 0.2; 		//height of initial shock wave
                    rc = 0.5;
                    zc = 0.5; 		//spherical centre
                    EoS = "ideal";
                    EM_source = false;
                    rigid_bodies = true;
                    cylindrical_coords = false;
                }
                
                if(CASE == "rigid_validation3"){
                    title = "rigid_validation3";
                    datfile = "rigid_validation3.dat";
                    gpname = "./visualisation/rigid_validation3.gp";
                    giffile = "./visualisation/rigid_validation3.gif";
                    construction = "default_x"; //
                    geometry = "quadrilateral";
                    BCs = "transmissive";
                    exact_sol_file = "none";
                    WL.rho = 1.3764; WL.u_vec[0] = 0.394; WL.u_vec[1] = 0.0; WL.p = 1.5698; //down or left
                    WR.rho = 1.0; WR.u_vec[0] = 0.0; WR.u_vec[1] = 0.0; WR.p = 1.0;  //up or right
                    n = 104;
                    m = 104;
                    R = 1.0;    //space domain 0 < x < 1
                    Z = 1.0;
                    Tf = 0.4; 
                    T0 = 0.0;
                    CFL = 0.9;                    
                    nG = 2;
                    dr = R/(n-2*nG-1);
                    dz = Z/(m-2*nG-1);
                    r0 = 0.2;		//spherical nose cone radius
                    z0 = 0.2; 		//height of initial shock wave
                    rc = 0.6;		//radial centre for rigid body
                    zc = 0.5; 		//z-centre
                    rl = 0.4;		//radial length for rigid body
                    zl = 0.4;		//z-length
                    EoS = "ideal";
                    EM_source = false;
                    rigid_bodies = true;
                    cylindrical_coords = false;
                }
								
                radial_vals.resize(n);
                for(int i = 0; i<n; ++i){
                    radial_vals[i] = (i-nG)*dr + dr/2.0; 
                }
                z_vals.resize(m);
                for(int j = 0; j<m; ++j){
                    z_vals[j] = (j-nG)*dz + dz/2.0; 
                }
                
                geom_vec.resize(n*m);  //fluid = positive region, rigid body = negative region
				vector<double> grads(2);
				geom_grads.resize(m*n, grads);
				signed_grads.resize(m*n, grads);
            }
};

#endif 
