//these classes exist in .cpp files, but their forward declaration is included here:
#ifndef CLASSDEFS_H_
#define CLASSDEFS_H_

using namespace std;

/*-------------------class definitions for: ---------------------------
 *   - Primitive variables (W) 
 *   - Conservative variables (U)
 *   - Conservative fluxes in x direction (F): 
 * 	 	- Conservative fluxes in y direction (G) 
 *   	- Radial source terms for cylindrical euler system (SE)
 * ------------------------------------------------------------------*/

class Prim{
    public:
        //Prim(double, Vector<Real>, double);
        double rho;
        Vector<Real> u_vec;
        double p;
        double r;
        
        Prim(double _rho, Vector<Real> _u_vec, double _p){
		rho = _rho; u_vec = _u_vec; p = _p;
		}
};

class ConsU{
    public:
        //ConsU(double, Vector<Real>, double);
        double rho;
        Vector<Real> rhou_vec;
        double E;
        double r;
        
        ConsU(double _rho, Vector<Real> _rhou_vec, double _E){
			rho = _rho; rhou_vec = _rhou_vec; E = _E;
		}
};

class ConsF{
    public:
        //ConsF(double, Vector<Real>, double);
        double mass;
        Vector<Real> mom_vec;
        double en;
        
        ConsF(double _mass, Vector<Real> _mom_vec, double _en){
			mass = _mass; mom_vec = _mom_vec; en = _en;
		}
};

class AccessVariable
{
public:

	AccessVariable(Vector<std::string> const &problem_variables);
	
	int& operator[](std::string const &var);
	
	std::map<std::string, int> variables;
	
};

//PLASMA19 class definition:
//"mixture19_cns.txt" file is structured as follows:
/*
 * 
 * [num_densities=500 num_pressures=500	num_species=19] 		1x3
 * 
 * [specific gas constant for each species] 					1x19
 * 
 * [heats of formation for each species]						1x19
 * 
 * [characteristic vibrational temperature for each species] 	1x19
 * 
 * [list of densities] 				1x500
 * 
 * [list of pressures]				1x500
 * 
 * ------------------> rho
 * |
 * |
 * |	sound speeds				500x500
 * |
 * v
 * p
 * 
 * ------------------> rho
 * |
 * |
 * |	internal energies			500x500
 * |
 * v
 * p
 * 
 * ... x24 properties
 * 
 * The 24 properties are listed in order within the enumerated 
 * properties variable below. MFi is the mass fraction of the 
 * i-th species 
 * 
 */
					
class Plasma19{
	
	public:

	//below functions all defined in utility.cpp
	Plasma19 (string);
	Real getSpecificInternalEnergy(const Real rho, const Real p);
	Real getSoundSpeed(const Real rho, const Real p);
	Real getPressure(const Real rho, const Real e);
	Real getDensity(const Real p, const Real e);
	Real getTemperature(const Real rho, const Real p);
	Real getConductivity(const Real rho, const Real p);
	Real getThermalConductivity(const Real rho, const Real p);
	Real getSpeciesMassFraction(const Real rho, const Real p, const int s);
	Real getAdiabaticIndex();
	
	void testing_function();
	
	int n_EoS19; 
	
	//look-up table data range for rho-p:
	
	const Real rho_min = 0.001; //kg/m^3
	const Real rho_max = 10.0;  //kg/m^3
	
	const Real p_min = 1013.25; //Pa
	const Real p_max = 18238500; //Pa
	
	private:
		
	//some default values for ideal gas:
	Real m_adiabaticIndex = 1.4;
	Real m_cV = 1.0;
	Real m_tempRef = 0.0;
	Real GasConstantR = 287.0;

	/* n_: integer number
	 * v_: 1xn vector of values
	 * m_: nxn matrix of values
	 * m3_: nxnxn matrix of values
	 */

	int n_rho; 
	int n_p; 
	int n_species;

	Vector<Real> v_SpecificGasConstants;
	Vector<Real> v_HeatsofFormations;
	Vector<Real> v_VibrationalTemperatures;
	Vector<Real> v_Densities;
	Vector<Real> v_Pressures;

	Vector< Vector<Real> > m_SoundSpeeds;
	Vector< Vector<Real> > m_InternalEnergies;
	Vector< Vector<Real> > m_Temperatures;
	Vector< Vector<Real> > m_SigmaElectrical;
	Vector< Vector< Vector<Real> > > m3_SpeciesMassFractions;
	Vector< Vector<Real> > m_SigmaThermal;
	
	const int num_properties = 24;
	enum properties{ 	SoundSpeed = 1, 
				InternalEnergy,
				Temperature,
				SigmaElectrical,
				MF1, MF2, MF3, MF4, MF5,
				MF6, MF7, MF8, MF9, MF10,
				MF11,MF12,MF13,MF14,MF15,
				MF16,MF17,MF18,MF19,
				SigmaThermal		};
				
};

#endif
