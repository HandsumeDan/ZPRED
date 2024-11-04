# include <iostream>
# include <string>
# include "zpredUtils.h"

using namespace std;

# ifndef solution_H
# define solution_H

	struct excessVolume
	{double R;double Q;double a_CH3;double a_CH2;double a_CH;double a_benzene;double a_cyclohexane;double a_toluene;double a_OH;double a_CH3CHOHCH3;
	double a_H2O;double a_phenol;double a_CH3CO;double a_CHO;double a_COOH;double a_CH3COOH;double a_CH2O;double a_CHCl2;double a_CCl3;double a_aniline;double a_CNH2;};

	// Solution Properties Class
	class Solution
		{public:
			// Constructor
			Solution(string solvent,string solventConc,string solventConcType,string solute,string soluteConc,double T,string delimiter);
			Solution();
			string interpretMoleculeName(string Molecule);
			double calcMoleculeMolecularWeight(string structure);
			double getBondLength(string atomicBond);
			double calcPureSolventDensity(string solventName,double temperature);
			double calcPureSolventViscosity(string solventName,double temperature);
			double calcPureSolventDielectric(string solventName,double temperature);
			double* getDensityCoefficients(string saltName);
			double* getViscosityCoefficients(string saltName);
			double getDielectricLinearCoefficient(string soluteName);
			double calcSolutionDensity(double T);
			double calcSolutionViscosity(double T);
			double calcSolutionDielectric(double T);
			excessVolume getFunctionGroupInformation(string group);
			double getFunctionGroupInteractionParameter(string group,string group2);
			string getIonData(string saltConc,string saltName,string delimiter);
			int getSolventFunctionGroups(string theSolvent,string &theFunctionGroups);
			string interpretSolute(string solute);
			void initializeFunctionGroupInformation();
			double calcDebyeLength(double T);
			//void generateLEaPFile(string Molecule,string outFile);
			// Ion Data: valence,molarConc,angstromRadius,
			string ionData;
			// Solution Debye Length [Angstroms]
			double debyeLength;
			// Largest Ion [Angstroms]
			string largestIon;
			// Computed Solution Viscosity [Pa s]
			double viscosity;
			// Computed Solution Density [kg/L]
			double density;
			// Computed Solution Relative Dielectric
			double dielectric;
			// Solution Temperature [Kelvin]
			double temperature;
			int numSolvents;
			int numSolutes;
			int numIon;
			// Structural Formula(e) of Ions(s) under Assessment
			string* ion;
			// Ion Valence
			int* ionVal;
			// Molar Concentration [mol/L] of Ion(s)
			double* ionConc;
			// Ionic Radii of Ion(s) [Angstroms]
			double* ionRad;
			// Structural Formula(e) of Solvent(s) under Assessment
			string* solvents;
			// Function Groups Composing the Molecule
			string **solventFunctionGroup;		// [numSolvents][numFunctionGroup[numSolvents]]
			int *numFunctionGroup;
			// Mole Fractions of Solvent(s)
			double* solventsConc;
			// Molecular Weight of Solvent(s) [g/mol]
			double* solventsMW;
			// Structural Formula(e) of Solute(s) under Assessment
			string* solutes;
			// Molar Concentration [mol/L] of Solute(s)
			double* solutesConc;
			// Molecular Weight of Solute(s) [g/mol]
			double* solutesMW;
			double* moleFrac2MassFrac(double* x,double* MW,int numComps);
			double* massFrac2MoleFrac(double* w,double* MW,int numComps);
			double* massFrac2MoleFrac(string* sW,double* MW,int numComps);
			double* volFrac2MoleFrac(double* v,double* MW,double* density,int numComps);
			double* volFrac2MoleFrac(string* sV,double* MW,double* density,int numComps);
			double* molarity2MassFrac(double* M,double* MW,int numComps,double solventDensity);
			// Number of Elements on Periodic Table
			static const int numElements=118;
			// Atomic Symbols of Elements on Periodic Table
			string atomicSymbol[numElements]={"H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","ut","Fl","up","Lv","us","uo"};
			// Atomic Names of ELements on Periodic Table
			string atomicName[numElements]={"Hydrogen","Helium","Lithium","Beryllium","Boron","Carbon","Nitrogen","Oxygen","Fluorine","Neon","Sodium","Magnesium","Aluminium","Silicon","Phosphorus","Sulfur","Chlorine","Argon","Potassium","Calcium","Scandium","Titanium","Vanadium","Chromium","Manganese","Iron","Cobalt","Nickel","Copper","Zinc","Gallium","Germanium","Arsenic","Selenium","Bromine","Krypton","Rubidium","Strontium","Yttrium","Zirconium","Niobium","Molybdenum","Technetium","Ruthenium","Rhodium","Palladium","Silver","Cadmium","Indium","Tin","Antimony","Tellurium","Iodine","Xenon","Cesium","Barium","Lanthanum","Cerium","Praseodymium","Neodymium","Promethium","Samarium","Europium","Gadolinium","Terbium","Dysprosium","Holmium","Erbium","Thulium","Ytterbium","Lutetium","Hafnium","Tantalum","Tungsten","Rhenium","Osmium","Iridium","Platinum","Gold","Mercury","Thallium","Lead","Bismuth","Polonium","Astatine","Radon","Francium","Radium","Actinium","Thorium","Protactinium","Uranium","Neptunium","Plutonium","Americium","Curium","Berkelium","Californium","Einsteinium","Fermium","Mendelevium","Nobelium","Lawrencium","Rutherfordium","Dubnium","Seaborgium","Bohrium","Hassium","Meitnerium","Darmstadtium","Roentgenium","Copernicium","Ununtrium","Flerovium","Ununpentium","Livermorium","Ununseptium","Ununoctium"};
			// Atomic Molecular Weight of Elements on Periodic Table [g/mol]
			double atomicMolecularWeight[numElements]={1.0079,4.0026,6.941,9.0122,10.811,12.0107,14.0067,15.9994,18.9984,20.1797,22.9897,24.305,26.9815,28.0855,30.9738,32.065,35.453,39.948,39.0983,40.078,44.9559,47.867,50.9415,51.9961,54.938,55.845,58.9332,58.6934,63.546,65.39,69.723,72.64,74.9216,78.96,79.904,83.8,85.4678,87.62,88.9059,91.224,92.9064,95.94,98,101.07,102.9055,106.42,107.8682,112.411,114.818,118.71,121.76,127.6,126.9045,131.293,132.9055,137.327,138.9055,140.116,140.9077,144.24,145,150.36,151.964,157.25,158.9253,162.5,164.9303,167.259,168.9342,173.04,174.967,178.49,180.9479,183.84,186.207,190.23,192.217,195.078,196.9665,200.59,204.38,207.2,208.9804,209,210,222,223,226,227,232.0381,231.0359,238.0289,237,244,243,247,247,251,252,257,258,259,262,261,262,263,262,265,266,281,272,285,284,289,288,292,0,294};
			// Atomic Charge of Elements on Periodic Table
			string atomicCharge[numElements]={"1+","0","1+","2+","3+","4+","3-","2-","1-","0","1+","2+","3+","4+","3-","2-","1-","0","1+","2+","3+","4+|3+","5+|4+","3+|2+","2+|4+","3+|2+","2+|3+","2+|3+","2+|1+","2+","3+","4+","3-","2","1-","0","1+","2+","3+","4+","5+|3+","6+","7+","3+|4+","3+","2+|4+","1+","2+","3+","4+|2+","3+|5+","2-","1-","0","1+","2+","3+","3+","3+","3+","3+","3+|2+","3+|2+","3+","3+","3+","3+","3+","3+","3+|2+","3+","4+","2","6+","7+","4+","4+","4+|2+","3+|1+","2+|1+","2","2+|4+","3+|5+","2+|4+","1-","0","1+","2+","3+","4+","5+|4+","6+|4+","5+","4+|6+","3+|4+","3+","3+|4+","3+","3+","3+","2+|3+","2+|3+","3+","NAN","NAN","NAN","NAN","NAN","NAN","NAN","NAN","NAN","NAN","NAN","NAN","NAN","NAN","NAN"};
			// Ion Properties Data
			static const int numIons=79;
			// Ion Names that are used to identify input ions
			string modeledIons[numIons]={"Ag","Al","AsO4","B(OH)4","Ba","Be","BO3","Br","BrO","BrO3","Ca","Cd","CH3COO","CH3NH3","Cl","ClO","ClO2","ClO3","ClO4","CN","Co","Co","CO3","HCO2","Cr","Cr","CrO4","Cr2O7","Cs","Cu","Cu","C2O4","HC2O4","F","Fe","Fe","H","Hg","In","IO","IO3","IO4","K","Li","Mg","Mn","Mn","MnO4","Na","NH4","Ni","Ni","NO3","NO2","O2","OH","HCO3","HSO3","HSO4","H2PO4","HPO4","PO4","HPO3","PO3","Rb","Rh","SeO4","Sn","Sn","SO4","SO3","S2O3","Sr","Th","Tl","Tl","Y","Zn","Citrate"};
			// Ion Structural Formulae
			string modeledIonsStructures[numIons]={"Ag","Al","As[::O]O3","B[OH]4","Ba","Be","BO3","Br","BrO","Br[::O]2O","Ca","Cd","CH3C[::O]O","CH3NH3","Cl","ClO","Cl[::O]O","Cl[::O]2O","Cl[::O]3O","C:::N","Co{2}","Co{3}","C[::O]O2","HC[::O]O","Cr{2}","Cr{3}","Cr[::O]2O2","OCr[::O]2OCr[::O]2O","Cs","Cu{1}","Cu{2}","OC[::O]C[::O]O","OC[::O]C[::O]OH","F","Fe{2}","Fe{3}","H","Hg{2}","In","IO","I[::O]2O","I[::O]3O","K","Li","Mg","Mn{2}","Mn{3}","Mn[::O]3O","Na","NH4","Ni{2}","Ni{4}","N[::O]2O","N[::O]O","O2","OH","C[::O][O]OH","S[::O][O]OH","S[::O]2[O]OH","P[::O][O][OH]2","P[::O][O]2OH","P[::O][O]3","P[O]2OH","PO3","Rb","Rh{3}","Se[::O]2[O]2","Sn{2}","Sn{4}","S[::O]2[O]2","S[::O][O]2","S::S[::O][O]2","Sr","Th","Tl{1}","Tl{3}","Y","Zn","OC[::O]C[H]2C[O][C[::O]O]C[H]2C[::O]O"};
			// Ion Chemical Names
			string modeledIonsNames[numIons]={"Silver","Aluminum","Arsenate","Tetrahydroxyborate","Barium","Beryllium","Borate","Bromide","Hyprobromite","Bromate","Calcium","Cadmium","Acetate","Methylammonium","Chloride","Hypochlorite","Chlorite","Chlorate","Perchlorate","Cyanide","Cobalt (II)","Cobalt (III)","Carbonate","Formate","Chromium (II)","Chromium (III)","Chromate","Dichromate","Cesium","Copper (I)","Copper (II)","Oxalate","Hydrogenoxalate","Floride","Iron (II)","Iron (III)","Hydrogen","Mercury (II)","Indium","Hypoiodite","Iodate","Periodate","Potassium","Lithium","Magnesium","Manganese (II)","Manganese (III)","Permanganate","Sodium","Ammonium","Nickel (II)","Nickel (IV)","Nitrate","Nitrite","Peroxide","Hydroxide","Bicarbonate","Bisulfite","Bisulfate","Dihydrogen Phosphate","Hydrogen Phosphate","Phosphate","Hydrogen Phosphite","Phosphite","Rubidium","Rhodium","Selenate","Tin (II)","Tin (IV)","Sulfate","Sulfite","Thiosulfate","Strontium","Thorium","Thallium (I)","Thallium (III)","Yttrium","Zinc","Citrate"};
			// Ionic Radii (A); Reference: 1988. Ionic radii in aqueous solutions
			string modeledIonsRadius[numIons]={"0.99","0.46","2.44","2.44","1.48","0.33","1.91","1.95","1.49","1.54","1.00","0.88","2.18","2.28","1.76","2.16","2.18","1.71","2.28","1.91","0.68","0.68","1.78","1.28","0.54","0.54","2.56","2.98","1.71","0.98","0.98","2.32","2.32","1.45","0.69","0.61","0.25","1.00","0.73","1.60","1.22","2.49","1.37","0.66","0.67","0.77","0.77","2.29","0.93","1.18","0.64","0.64","2.03","1.92","1.42","1.33","1.56","2.17","2.28","2.38","2.30","2.23","2.22","2.22","1.50","0.61","2.53","1.20","1.20","2.39","2.39","2.44","1.22","1.11","0.81","0.81","0.93","0.67","3.47"};
			// Ion Charge Valence
			string modeledIonsCharge[numIons]={"1","3","-3","-1","2","2","-3","-1","-1","-1","2","2","-1","1","-1","-1","-1","-1","-1","-1","2","3","-2","-1","2","3","-2","-2","1","1","2","-2","-1","-1","2","3","1","2","3","-1","-1","-1","1","1","2","2","3","-1","1","1","2","4","-1","-1","-2","-1","-1","-1","-1","-1","-2","-3","-2","-3","1","3","-2","2","4","-2","-2","-2","2","4","1","3","3","2","-3"};

			// Functional Groups Modeled
			static const int numFunctionalGroups=19;
			string functionGroups[numFunctionalGroups]={"CH3","CH2","CH","*1CH(CH)4*1CH","*1CH2(CH2)4*1CH2","*1CH(CH)4*1C[CH3]","OH","CH3C[OH]HCH3","H2O","*1CH(CH)4*1C[OH]","CH3C::O","C[::O]H","C[::O]OH",\
												"CH3C[::O]OH","CH2O","CHCl2","CCl3","*1CH(CH)4*1C[NH2]","CNH2"};
			excessVolume CH3;
			excessVolume CH2;
			excessVolume CH;
			excessVolume benzene;
			excessVolume cyclohexane;
			excessVolume toluene;
			excessVolume OH;
			excessVolume CH3CHOHCH3;
			excessVolume H2O;
			excessVolume phenol;
			excessVolume CH3CO;
			excessVolume CHO;
			excessVolume COOH;
			excessVolume CH3COOH;
			excessVolume CH2O;
			excessVolume CHCl2;
			excessVolume CCl3;
			excessVolume aniline;
			excessVolume CNH2;
		private:
			int countDelimiter(string Data,string delimiter);
			int* fillIntArray(string Data,int numPnts,string delimiter);
			double* fillDoubleArray(string Data,int numPnts,string delimiter);
			string* fillStringArray(string Data,int numPnts,string delimiter);
			static const int numModeledSolvents=39;
			// List of Modeled Solvents Structures
			string modeledSolvents[numModeledSolvents]={"CH3C::[O]OH","CH3C[CH3]::O","CH3C:::N","CH3(CH2)3OH",\
												"CH3(CH2)3C[::O]CH3","*1CHCHC[Cl]CHCH*1CH","CHCl3","*1CH2(CH2)4*1CH2","CH3(CH2)3O(CH2)3CH3","C[Cl]H2C[Cl]H2","ClC[Cl]H2","CH3OCH2CH2OCH3",\
												"CH3C[::O]N[CH3]CH3","CH3N[CH3]C[::O]H","*1CH2OCH2CH2O*1CH2","CH3S[::O]CH3","CH3CH2OH","CH3C[::O]OCH2CH3","CH3CH2OCH2CH3",\
												"HOCH2CH2OH","CH3(CH2)3C[CH2CH3]HCH2OH","NH2C[::O]H","HOCH2C[OH]HCH2OH","CH3(CH2)4CH3","CH3C[CH3]HCH2OH","CH3OH","CH3O(CH2)2OH",\
												"CH3CH2C[::O]CH3","CH3C[CH3]HC[::O]CH3","CH3N[::O]O","CH3C[OH]HCH2OH","CH3(CH2)2OH","CH3C[OH]HCH3","*1CH2(CH2)3*1O",\
												"*1CH(CH)2C[CH2(CH2)2*2CH2]*2C*1CH","CH3C[CH3]HCH2C[::O]CH3","*1CH(CH)4*1C[CH3]","ClC[Cl]::C[Cl]H","H2O"};
			// Structural Formula Interpretation:
			// Paranthesis () should always be followed by a number indicating the number of repeats
			// Square brackets [] indicate bonding with previous atom but not next atom, i.e. a branch
			// Square brackets can be followed by a number indicating the number of repeats of the branch connected to the central atom
			// Any Square brackets following square brackets means atom bonds to previous atom with no square brackets
			// for cyclic and polycyclic compounds, connections are specified by *Num preceding atom
			// Squiggly Brackets {} designates the oxidation state of the preceding atom
			// (single bond is assumed so must indicate double bond as :: and triple bond as :::)
			string modeledSolventsNames[numModeledSolvents]={"acetic acid","acetone","acetonitrile","butanol","butyl methyl ketone","chlorobenzene","chloroform","cyclohexane","dibutyl ether",\
													"1,2-dichloroethane","dichloromethane","1,2-dimethoxyethane","n,n-dimethylacetamide","n,n-dimethylformamide","1,4-dioxane",\
													"dimethyl sulfoxide","ethanol","ethyl acetate","diethyl ether","ethylene glycol","2-ethylhexanol","formamide","glycerol",\
													"hexane","isobutanol","methanol","2-methoxyethanol","methyl ethyl ketone","methyl isopropyl ketone","nitromethane",\
													"propylene glycol","propanol","isopropanol","tetrahydrofuran","tetralin","methyl isobutyl ketone","toluene","trichloroethene","water"};
		};

# endif
