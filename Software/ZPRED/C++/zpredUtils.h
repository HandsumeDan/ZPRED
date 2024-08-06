# include <algorithm>
# include <cmath>
# include <cstdio>					// defines remove
# include <cstdlib>
# include <cstring>
# include <ctime>					// defines clock_t, clock, CLOCKS_PER_SEC
# include <exception>
# include <fstream>
# include <iostream>
# include <new>
# include <sstream>
# include <string>
# include <thread>
# include <utility>
# include <vector>
# include <fcntl.h>					// Defines O_CREAT flag
# include <ftw.h>					// for remove and nftw
# include <sys/reg.h>
# include <sys/syscall.h>			// Defines OS System Calls
# include <sys/ptrace.h>
# include <sys/types.h>
# include <sys/user.h>
# include <sys/wait.h>
# include <sys/stat.h>				// Defines chmod function
# include <boost/filesystem.hpp>		// create_directory
# include <boost/numeric/odeint.hpp>
# include <openssl/conf.h>
# include <openssl/evp.h>				// ciphers and message digests
# include <openssl/err.h>
# include <openssl/ssl.h>				// for libssl
# include <unistd.h>
# include "ap.h"
# include "solvers.h"
# include "stdafx.h"

using namespace std;
using namespace boost::numeric::odeint;

# ifndef zpredUtils_H
# define zpredUtils_H

	// Acceptable Deviation During Newton-Ralphson Iteration
	# define ACCEPTED_ERROR_PERCENT 5e-4
	// Default Length of a Number to consider for cnvrtNumToStrng
	# define BUFFER_SIZE 12
	// Number of digits behind  point for ion radii and inflation distance/slip plane position
	# define DISTANCE_SIG_FIGS 2
	// Number of digits behind  point for solvent dielectric constant
	# define ER_SIG_FIGS 2
	// Default Number of Threads to Execute ZPRED on (if not user specified)
	# define MAX_THREADS 3
	// Maximum Number of Rows for Terminal Display
	# define MAX_DISPLAY_ROWS 25
	// Maximum Number of Columns for Terminal Display
	# define MAX_DISPLAY_COLS 6
	// Maximum Number of Iterations to Solve for PBE Coefficient
	# define MAX_PBE_ITERATIONS 1000
	// Maximu Number of File Open Attempts
	# define MAX_FILE_OPEN_ATTEMPTS 100000

	// Physical Constants
	// Avogadro's Number (Particles per mole)
	# define Na 6.022140857e23
	// Universal Gas Constant (Joule/(Kelvin*mol))
	# define Rg 8.3144598
	//Boltzmann Constant (Joule/Kelvin)
	# define kb Rg/Na
	//Electron (or Elementary) Charge (Coulomb)
	# define ec 1.6021766208e-19
	// PI
	# define pi 3.14159265358979323846
	// Permeability of Free Space (Tesla*meter/Ampere)
	# define uo pi*0.0000004
	// Speed of Light in Vacuum (meter/second)
	# define sol 299792458
	// Permittivity of Free Space (Farad/meter)
	# define eo 1/(uo*sol*sol)

	// Executables
	# define PYTHON_EXE "python2.7"	

	struct Reference
	{double** E;double** y;double* ka;int numPlots;int* numPnts;};

	// ZPRED Component Software Executables
	struct softwareComponents
	{string apbsExecutable;string multivalueExecutable;string hydroproExecutable;string msmsExecutable;string msmsAtmTypeNumbers;string msmsPdbToXyzr;string msmsPdbToXyzrn;string pdb2pqr;};

	// MSMS Input Parameter Data Class
	struct msmsInput
	{string PDB;string pdbFile;string msmsFldr;string pdb_msmsFldr;string pdb_surfFldr;string surfPointDensity;string surfPointHiDensity;string probeRadius;string PDB_TO_XYZRN;string MSMS;bool SUCCESS;};

	// ZPRED Input Parameter Data Structure Array (for taking data from GUI)
	struct zInput
	{int numSolCond;
	string* id;
	double* shapeFactor;
	double* molecularWeight;
	string* files;
	int numFiles;
	string name;
	string* pH;
	string** solvent;
	int* numSolvent;
	double** solventConc;
	string* solventConcType;
	string** solute;
	int* numSolute;
	double** soluteConc;
	string* soluteConcType;
	double* temperature;
	double* proteinConc;
	string outputTitle;
	string* apbsWriteTypes;
	int numApbsWriteTypes;
	softwareComponents SC;
	msmsInput msms;
	string hydroproCalcType;
	string forceField;
	bool saveTemps;
	string pdie;
	string dime_x;string dime_y;string dime_z;
	bool RUN_APBS;bool RUN_PDB2PQR;bool RUN_HYDROPRO;bool RUN_HULLRAD;bool RUN_ZPRED;bool RUN_COLLIDE;
	int maxThreads;
	// User Specified Values (-1 if not specified)
	double* dielectric;double* viscosity;double* density;double* inflationDistance;
	// Generate Pot Profile
	bool* GENERATE_POT_PROFILE;
	string plot;};

	// ZPRED Input Parameter Structure Array for Parallel Execution of the Problem (each represents only one file and solution condition)
	struct zExeInput
	{string id;
	string file;
	string pH;
	string* solvent;
	int numSolvent;
	double* solventConc;
	string solventConcType;
	string* solute;
	int numSolute;
	double* soluteConc;
	string soluteConcType;
	double temperature;
	double proteinConc;
	string outputTitle;
	string* apbsWriteTypes;
	int numApbsWriteTypes;
	softwareComponents SC;
	msmsInput msms;
	string hydroproCalcType;
	string forceField;
	bool saveTemps;
	string pdie;
	string dime_x;string dime_y;string dime_z;
	bool RUN_APBS;bool RUN_PDB2PQR;bool RUN_HYDROPRO;bool RUN_HULLRAD;bool RUN_ZPRED;bool RUN_COLLIDE;
	int maxThreads;
	// User Specified Values (-1 if not specified)
	double dielectric;double viscosity;double density;double inflationDistance;
	// Generate Pot Profile
	bool GENERATE_POT_PROFILE;
	};

	// Single Solution Properties Data Structure Array
	struct solProp
	{string pH;string temperature;string dielectric;string viscosity;string density;int numSolute;string* solute;string* soluteConc;string debyeLength;int numSolvent;string* solvent;string* solventConc;};

	// Total Solution Properties Data Structure Array
	struct allSolProp
	{string* pH;string* temperature;string* dielectric;string* viscosity;string* density;int numSolute;string** solute;string** soluteConc;string* debyeLength;int numSolvent;string* solvent;string* solventConc;};

	struct electricProfile
	{double* position;double* potential;int numPnts;};

	// Single ZPRED Output Data Structure Array
	struct zOutput
	{string proteinName;string solCondNum;double molecularWeight;double proteinRadius;double solvatedRadius;double zetaPotential;string Xsp;double charge;double semMobility;double henryMobility;double kuwabaraMobility;double diffusivity;solProp solution;electricProfile EP;};

	// COLLIDE Input Data Structure Array
	struct colInput
	{int numSols;string* pRadius;double* proteinRadius;double* proteinRadiusStd;string* sRadius;double* solvatedRadius;double* solvatedRadiusStd;string* zPotential;double* zetaPotential;double* zetaPotentialStd;string* chrg;double* charge;double* chargeStd;string* eMob;double* mobility;double* mobilityStd;string* diff;double* diffusivity;double* diffusivityStd;allSolProp solution;double* collisionEfficiency;double* coagulationTime;};

	// General Vector
	struct Vec
	{double x;double y;double z;};

	// 2D Vector
	struct Vec2d
	{double x;double y;};

	// String Vector
	struct strVec
	{string x;string y;string z;};	

	// PDB File Data Type
	struct pdbData
	{string* recordType;int* atomNum;string* atomLabel;string* residue;string* chainLabel;int* atmChainNum;double* x;double* y;double* z;double* atmOccupancy;double* bFactor;string* terminus;};

	// Vertices Data Structure Array [output of MSMS]
	struct vertData
	{double* x;double* y;double* z;double* nx;double* ny;double* nz;int* analyticSurf;int* sphereNum;int* vertType;string* closestAtom;int numVertices;};

	// Structure Array for Ion Data
	struct ion_Data
	{double* conc;int* valence;double* mobility;double* radius;int numIons;};

	// Collision Computation Parameters
	struct collisionParameters
	{double zeta;double Rp;double Rh;double Xsp;};

	void perform_ZPRED_computation(zInput Z);
	double calc_coagulation_time(double protConc,double protMW,double viscosity,double T);
	double calc_collision_efficiency(double Rp,double Rh,double zetaPot,double er,double debyeLength,double T);
	string getBaseFolder(string f);
	string checkWhiteSpaceInFilePath(string iFile);
	double interpolate(double x,double x1,double x2,double y1,double y2);
	string array2String(int* Data,int numPnts,string delimiter);
	double calcHydrodynamicRadius(double diffusivity,double temperature,double viscosity);
	string checkFilePathParantheses(string f);
	void collect_all_files_in_dir(string dir2Search,string &fileLst,string delimiter,int &N);
	bool directory_exist(string fPath);
	bool file_exist(string fPath);
	string* get_all_files_in_dir(string dir2Search,int &N);
	string get_containing_folder(string fPath);
	string* get_dirs_in_dir(string fPath,int &N);
	string get_file_name(string fPath);
	string* get_files_in_dir(string fPath,int &N);
	bool IN_STRING_ARRAY(string X,string* L,int N);
	bool IN_LIST(string X,string list,int N,string delimiter);
	string nameBuffer(string soluteName,string soluteConc);
	void sort_array(string *theArr,int N);
	void chmod_x(string binFile);
	//void copyFile(string inFile,string outFile);
	void convert_pqr_to_pdb(string saveAsPdbFile,string pqrFile,string pdbFile);
	bool copyFile(string srcFile,string destFile);
	int countNumLinesMSMSVertFile(string fileName);
	void determine_dx_file_name_components(string fNm,string &pHVal,string &tVal,string &srVal);
	void make_directory(string fPath);
	void make_folder(string Fldr);
	string makeUpperCase(string X);
	string makeWhiteSpace(int Sz);
	string makeHTMLWhiteSpace(int Sz);
	int remove_directory(string fPath);
	void remove_file(string fPath);
	int unlink_cb(const char *fPath,const struct stat *sb, int typeFlag, struct FTW *ftwbuf);
	void resizeTerminal(int width,int height);
	zInput read_zpred_input_file(string inFile);
	zExeInput cnvrtStr2Input(string cmd,string delimiter,string& Fldr,string& calcNum,string& display);
	string* cnvrtzInput2Str(zInput In,string delimiter,string zFldr);
	string cnvrtNumToStrng(int Num,int numberAfterDecimalpoint);
	string cnvrtNumToStrng(double Num,int numberAfterDecimalpoint);
	string cnvrtNumToStrng(long double Num,int numberAfterDecimalpoint);
	string cnvrtNumToStrng(float Num,int numberAfterDecimalpoint);
	string cnvrtScientificNumToStrng(double Num);
	string cnvrtScientificNumToStrng(long double Num);
	string cnvrtScientificNumToStrng(float Num);
	int count_delimiter(string Data,string delimiter);
	void define_msmsInput(string PDB,string pdbFile,string msmsFldr,string pdb_msmsFldr,string pdb_surfFldr,string surfPointDensity,string surfPointHiDensity,string probeRad,string msmsBinFile,msmsInput& In);
	void erase_display();
	void handleErrors(void);
	static inline bool is_base64(unsigned char c);
	void decryptString(string encStr,const string salt,const string pass,string outFile);
	string base64_decode(string const& encoded_string);
	string base64_encode(unsigned char const* bytes_to_encode, unsigned int in_len);
	string readFile2Encrypt(string inFile);
	string encryptFile(string inFile,const string salt,const string pass);
	void extract_models_from_pqr_file(string pqrFile,string modelFldr);
	string getFileName(string filePath);
	double* fill_double_array(string Data,int numPnts,string delimiter);
	int* fill_int_array(string Data,int numPnts,string delimiter);
	string* fill_string_array(string Data,int numPnts,string delimiter);
	string formatNumberString(string Num);
	long getFileSize(string inFile);
	string initialize_ZPRED_display(int totalCalc,bool& MULTI_DISPLAY);
	string initialize_ZPRED_gui_display(int totalCalc,bool& MULTI_DISPLAY);
	string repStr(string Input,int Num);
	string repeatString(string Input,int Num);
	void rewind_display(string display);
	void* runZPRED(void* Ptr);
	void updateDisplayString(string dispStr,string calcNum,int updatePos,int nRows,string& dispLn);
	void updateGuiDisplayString(string dispStr,string calcNum,int updatePos,int nRows,string dispLn,string guiDisplayFile);
	string update_display(int totalCalc,int currCalc,int threadOffset);
	void writeAPBSInputGen(string oFile);
	void writeHullRad(string oFile);
	void writePsize(string oFile);
	void* runMSMS(msmsInput& In,string calcNum);
	void readMSMSVertFile(string fileName,vertData& Data);
	void writeHydroproInput(string calcType,string density,string molecularWeight,string proteinName,string oFile,string PDB,double specificVolume,double temperature,double viscosity);
	void runHYDROPRO(string PDB,string pdb_hydroFldr,string inFile,string logFile,string calcNum);
	double getDiffusivity(string hydroproFile);
	void runPROPKA_PDB2PQR(string pdbFile,string pqrFile,string calcNum,string forceField,string pdb_pqrFldr,string propkaOutFile,string pH,string PDB2PQR,string logFile);
	void runAPBS(string exeFile,string inFile,string logFile);
	int countNumLines(string file);
	void readPDBFile(string file,pdbData& Data);
	double calc_protein_radius_std(vertData Data,int Sz,Vec centerCoord,double Radius);
	double calc_protein_radius(vertData Data,int Sz,Vec centerCoord);
	double calc_distance(Vec Vec1,Vec Vec2);
	Vec calc_center_of_mass(pdbData Data,int Sz);
	double getSpecificVolume(string msmsOutputFile,double molecularWeight);
	void definePDBData(pdbData& Data,int Sz);
	double calcProteinMolecularWeight(string pdbFile);
	void runMULTIVALUE(string exeFile,string csvFile,string potFile,string mvFile,string logFile);
	void inFileEditer(string* writeTypes,int N,string pdie,strVec dim,string pqrFile,string apbsInputgen,string ionData,string temperature,string largestIon,string er,string fldrExt);
	void vert2csv(string inVertFile,string outCsvFile,bool VERBOSE);
	void inflateVert(double inflationDistance,string inVertFile,string outVertFile,bool VERBOSE);
	void writeMSMSVertFile(string fileName,vertData Data);
	float readPQRFile(string pqrFile);
	double readPROPKAFile(string proFile);
	void clear_ion_Data(ion_Data& In);
	double calc_OF_for_Sphere(ion_Data Ion,double er,double temperature,double zetaPot,double ka);
	double calc_HF_for_Cylinder(double ka);
	double calc_HF_for_Sphere(double ka);
	double calc_OHWF_for_Sphere(ion_Data Ion,double er,double temperature,double zetaPot,double ka);
	double getDiffusivity_hullrad(string hullradFile);
	//
	double calc_Ion_Mobility(int ion_val,double ion_radius,double viscosity);
	double calc_Ion_Friction_Coefficient(int ion_val,double ion_mobility);
	zOutput read_zpred_output_file(string zOutFile);
	double calc_std_dev(double* X,int Sz,double Avg);
	double calc_average(double* X,int Sz);
	double calc_Kuwabara_Correction_for_Sphere(double ka,double volFrac,double Zeta,double a);
	double shoelace(Vec2d p1,Vec2d p2,Vec2d p3);
	double shoelace(Vec2d p1,Vec2d p2,Vec2d p3,Vec2d p4);
	double calcAreaUnderCurve(double* x,double* y,int N);
	void write_saveAsPdb(string fPath);
	void remove_solvent(string removeSolventFile,string pdbFile);

# endif
