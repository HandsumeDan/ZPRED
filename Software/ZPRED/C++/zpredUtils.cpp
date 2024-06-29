# include <dirent.h>
# include "mobility.h"
# include "solution.h"
# include "zpredUtils.h"

double calc_average(double* X,int Sz)
	{double Output=0;
	int Counter=0;
	for(int i=0;i<Sz;i++)
		{if(!isnan(X[i]))
			{Output+=X[i];
			Counter++;}
		}
	if(Counter==0){return 0;}
	return Output/Counter;}

double calc_coagulation_time(double protConc,double protMW,double viscosity,double T)
	{// Input
	// protConc		protein concentration (grams/Liter)
	// protMW			molecular weight of protein structure (grams/mol)
	// viscosity		dynamic viscosity (Pa s)
	// T				absolute tempertaure (Kelvin)
	double Output;		// [seconds]
	
	double trm1=3*viscosity;
	double trm2=4*kb*T;
	trm2*=Na;
	trm2*=protConc;
	trm2/=protMW;
	trm2*=1000;

	Output=trm1/trm2;
	return Output;}

double calc_collision_efficiency(double Rp,double Rh,double zetaPot,double er,double debyeLength,double T)
	{// Input
	// Rp			anhydrous protein radius (Angstroms)
	// Rh			solvated protein radius (Angstroms)
	// zetaPot		zeta potential (volts)
	// er			solution relative diecltric constant (unitless)
	// debyeLength		Debye/EDL Length (Angstroms)
	// T				absolute temperature (Kelvin)
	double Output;		// fraction of particles successfully colliding and sticking
/*	
	// Attraction Energy
	double hamakerConstant=1e-21;		// Joules
	double Ea=-hamakerConstant*Rp/(12.0*d);
	
	// Repulsion Energy
	double trm1=exp(-(d-2*(Rh-Rp))/debyeLength);
	double trm2=log(1+trm1);
	double Er=2*pi*eo*er*Rh*zetaPot*zetaPot*trm2;


	double incrmnt=0.1;		// tenth an Angstrom
	double minDist=2*Rp;
	double maxDist=minDist+10*debyeLength;
	int minNumPnts=250;
	int numPnts=(maxDist-minDist)/incrmnt + 1;
	if(numPnts<minNumPnts)
		{numPnts=minNumPnts;
		incrmnt=(maxDist-minDist)/(numPnts-1);}
	double *x=new double[numPnts];
	double *y=new double[numPnts];
	for(int i=0;i<numPnts;i++)
		{x[i]=minDist+i*incrmnt;
		// Attraction Energy
		Ea=-hamakerConstant*Rp/(12.0*x[i]*1e-10);
		// Repulsion Energy
		trm1=exp(-(x[i]-2*(Rh-Rp))/debyeLength);
		trm2=log(1+trm1);
		Er=2*pi*eo*er*Rh*1e-10*zetaPot*zetaPot*trm2;
		// Total Interaction Energy
		Et=Ea+Er;
		y=exp(Et/(kb*T))/(x[i]*x[i]*1e-20);
		}
	double W=calcAreaUnderCurve(x,y,numPnts);
	Output=1/(2*Rp*1e-10*W);
*/
	// Attraction Energy
	double hamakerConstant=4e-21;		// Joules
	double Ea,Er,Et,trm1,trm2,Xsp;

	Xsp=Rh-Rp;
	if(Xsp<0){Xsp=0;}

	double incrmnt=0.1;		// tenth an Angstrom
	double minDist=2*Rp;
	double maxDist=minDist+100*debyeLength;
	int minNumPnts=250;
	int numPnts=(maxDist-minDist)/incrmnt + 1;
	if(numPnts<minNumPnts)
		{numPnts=minNumPnts;
		incrmnt=(maxDist-minDist)/(numPnts-1);}
	double *x=new double[numPnts];
	double *y=new double[numPnts];
	for(int i=0;i<numPnts;i++)
		{// Inter-Particle Distance
		x[i]=minDist+i*incrmnt;
		// Attraction Energy (J)
		Ea=-hamakerConstant*Rp/(12.0*x[i]);
		// Repulsion Energy (J)
		trm1=exp(-(x[i]-2*Xsp)/debyeLength);
		trm2=log(1+trm1);
		Er=2*pi*eo*er*Rh*1e-10*zetaPot*zetaPot*trm2;
		// Total Interaction Energy
		Et=Ea+Er;
/*		if(i%10==0)
			{cout<<x[i]<<"|"<<i<<"/"<<numPnts<<endl;
			cout<<"Attraction: "<<Ea<<endl;
			cout<<"Repulsion:  "<<Er<<endl;
			cout<<"Total:      "<<Et<<endl;}
*/
		trm1=Et/(kb*T);
		trm2=exp(trm1);
		//cout<<trm1<<"|"<<trm2<<endl;
		y[i]=trm1/(x[i]*x[i]*1e-20);
		}
	
	
	double W=calcAreaUnderCurve(x,y,numPnts);

	Output=abs(1/(2*Rp*1e-10*W));
	
	return Output;}


double calc_std_dev(double* X,int Sz,double Avg)
	{double Output=0;
	int Counter=0;
	for(int i=0;i<Sz;i++)
		{if(!isnan(X[i]))
			{Output+=(X[i]-Avg)*(X[i]-Avg);
			Counter++;}
		}
	if(Counter==1){return sqrt(Output/Counter);}
	else if(Counter==0){return 0;}
	else{return sqrt(Output/(Counter-1));}}

// 1994. A Simple Expression for Henry's Funtion for the Retardation Effect in Electrophoresis of Spherical Colloidal Particles
double calc_HF_for_Sphere(double ka)
    {double trm1=2.5/(1+2*exp(-ka));
	double trm2=1+trm1/ka;
    double f=1 + 1/(2*trm2*trm2*trm2);
	return 2*f/3.0;}

// 1996. Henry's Function for Electrophoresis of a Cylindrical Colloidal Particle
double calc_HF_for_Cylinder(double ka)
    {double trm1=1;
     double trm2=ka*(1+exp(-ka));
     double trm3=1+(2.55/trm2);
     double f=trm1+1/(trm3*trm3);
     return (1+2*f)/3.0;}

// 2001. Approximate Analytic Expression for the Electrophroteic Mobility of Spherical Colloidal Particle
double calc_OF_for_Sphere(ion_Data Ion,double er,double temperature,double zetaPot,double ka)
	{if(Ion.numIons>2){/*cerr<<"Can not use Ohshima approximation, because too many ions\n";*/return 0;}
	double mCo,mCtr,zCtr,zCo;
	// Dimensionless Ionic Drag Coefficient
	if(zetaPot>0)
		{if(Ion.valence[0]>0)
			{// First Ion is Co-Ion
			mCo=4*pi*er*eo*kb*temperature*Ion.radius[0]*1e-10;
			mCo/=Ion.valence[0]*Ion.valence[0]*ec*ec;
			zCo=Ion.valence[0];
			// Second Ion is Counter-Ion
			mCtr=4*pi*er*eo*kb*temperature*Ion.radius[1]*1e-10;
			mCtr/=Ion.valence[1]*Ion.valence[1]*ec*ec;
			zCtr=Ion.valence[1];}
		else if(Ion.valence[0]<0)
			{// First Ion is Counter-Ion
			mCtr=4*pi*er*eo*kb*temperature*Ion.radius[0]*1e-10;
			mCtr/=Ion.valence[0]*Ion.valence[0]*ec*ec;
			zCtr=Ion.valence[0];
			// Second Ion is Co-Ion
			mCo=4*pi*er*eo*kb*temperature*Ion.radius[1]*1e-10;
			mCo/=Ion.valence[1]*Ion.valence[1]*ec*ec;
			zCo=Ion.valence[1];}}
	else if(zetaPot<0)
		{if(Ion.valence[0]>0)
			{// First Ion is Counter-Ion
			mCtr=4*pi*er*eo*kb*temperature*Ion.radius[0]*1e-10;
			mCtr/=Ion.valence[0]*Ion.valence[0]*ec*ec;
			zCtr=Ion.valence[0];
			// Second Ion is Co-Ion
			mCo=4*pi*er*eo*kb*temperature*Ion.radius[1]*1e-10;
			mCo/=Ion.valence[1]*Ion.valence[1]*ec*ec;
			zCo=Ion.valence[1];}
		else if(Ion.valence[0]<0)
			{// First Ion is Co-Ion
			mCo=4*pi*er*eo*kb*temperature*Ion.radius[0]*1e-10;
			mCo/=Ion.valence[0]*Ion.valence[0]*ec*ec;
			zCo=Ion.valence[0];
			// Second Ion is Counter-Ion
			mCtr=4*pi*er*eo*kb*temperature*Ion.radius[1]*1e-10;
			mCtr/=Ion.valence[1]*Ion.valence[1]*ec*ec;
			zCtr=Ion.valence[1];}}
	
	//if(abs(Ion.valence[0])>abs(Ion.valence[1]))
	//	{zCtr=Ion.valence[0];}
	//else
	//	{zCtr=Ion.valence[1];}

	double fHenry=calc_HF_for_Sphere(ka);
	double f3=ka*(ka+1.3*exp(-0.18*ka)+2.5);
	f3/=2*(ka+1.2*exp(-7.4*ka)+4.8)*(ka+1.2*exp(-7.4*ka)+4.8)*(ka+1.2*exp(-7.4*ka)+4.8);
	double f4=9*ka*(ka+5.2*exp(-3.9*ka)+5.6);
	f4/=8*(ka-1.55*exp(-0.32*ka)+6.02)*(ka-1.55*exp(-0.32*ka)+6.02)*(ka-1.55*exp(-0.32*ka)+6.02);

	double output;
	double trm1,trm2;
	double dZ=ec*zetaPot/(kb*temperature);
	//cout<<dZ<<endl;
	if(abs(Ion.valence[0])!=abs(Ion.valence[1]))
		{if(zetaPot>=0)
			{//output=3*fHenry/2-dZ*(abs(zCtr)-zCo)*f4-dZ*dZ*((zCo*mCo+abs(zCtr)*mCtr)/(abs(zCtr)+zCo))*f4;
			trm1=dZ*(abs(zCtr)-zCo)*f4;
			//cout<<trm1<<endl;
			trm2=dZ*dZ*((zCo*mCo+abs(zCtr)*mCtr)/(abs(zCtr)+zCo))*f4;
			//cout<<trm2<<endl;
			output=3*fHenry/2-trm1-trm2;}
			//output=fHenry-trm1-trm2;}
		else if(zetaPot<0)
			{//output=3*fHenry/2-dZ*(zCtr-abs(zCo))*f4-dZ*dZ*((abs(zCo)*mCo+zCtr*mCtr)/(zCtr+abs(zCo)))*f4;
			trm1=dZ*(zCtr-abs(zCo))*f4;
			//cout<<trm1<<endl;
			trm2=dZ*dZ*((abs(zCo)*mCo+zCtr*mCtr)/(zCtr+abs(zCo)))*f4;
			//cout<<trm2<<endl;
			output=3*fHenry/2-trm1-trm2;}}
			//output=fHenry-trm1-trm2;}}
	else		
		{output=3*fHenry/2-(zCtr*ec*zetaPot/(kb*temperature))*(zCtr*ec*zetaPot/(kb*temperature))*(f3+(mCo+mCtr)*f4/2);}
		//output=fHenry-(zCtr*ec*zetaPot/(kb*temperature))*(zCtr*ec*zetaPot/(kb*temperature))*(f3+(mCo+mCtr)*f4/2);}
	return 2*output/3;}

// 1983. Approximate Analytical Expressions for the Electrophoretic Mobility of Spherical Colloidal Particles and the Conductivity of their Dilute Suspensions
double calc_OHWF_for_Sphere(ion_Data Ion,double er,double temperature,double zetaPot,double ka)
	{if(Ion.numIons>2){/*cerr<<"Can not use Ohshima-Healy-White approximation, because too many ions"<<endl;*/return 0;}
	// Define Dimensionless Zeta Potential
	double zCtr,zCo,rCo,rCtr;
	if(abs(Ion.valence[0])>abs(Ion.valence[1]))
		{zCtr=Ion.valence[0];}
	else
		{zCtr=Ion.valence[1];}	
	double dZ=zCtr*ec*abs(zetaPot)/(kb*temperature);

	if(zetaPot>0)
		{if(Ion.valence[0]>0)
			{// First Ion is Co-Ion
			zCo=Ion.valence[0];
			rCo=Ion.radius[0]*1e-10;
			// Second Ion is Counter-Ion
			zCtr=Ion.valence[1];
			rCtr=Ion.radius[1]*1e-10;}
		else if(Ion.valence[0]<0)
			{// First Ion is Counter-Ion
			zCtr=Ion.valence[0];
			rCtr=Ion.radius[0]*1e-10;
			// Second Ion is Co-Ion
			zCo=Ion.valence[1];
			rCo=Ion.radius[1]*1e-10;}}
	else if(zetaPot<0)
		{if(Ion.valence[0]>0)
			{// First Ion is Counter-Ion
			zCtr=Ion.valence[0];
			rCtr=Ion.radius[0]*1e-10;
			// Second Ion is Co-Ion
			zCo=Ion.valence[1];
			rCo=Ion.radius[1]*1e-10;}
		else if(Ion.valence[0]<0)
			{// First Ion is Co-Ion
			zCo=Ion.valence[0];
			rCo=Ion.radius[0]*1e-10;
			// Second Ion is Counter-Ion
			zCtr=Ion.valence[1];
			rCtr=Ion.radius[1]*1e-10;}}

	double A=2.0*(1.0+(12.0*pi*er*eo*kb*temperature*rCtr/(ec*ec*zCtr*zCtr*abs(zCtr))))*(exp(dZ/2.0)-1.0)/ka;
	double B=log(1.0+exp(dZ/2.0))-log(2.0);
	double C=1.0-25.0*exp(-ka*dZ/(6.0*(ka+6.0)))/(3.0*(ka+10.0));
	double D=log(1.0+exp(-dZ/2.0))-log(2.0);
	double t=tanh(dZ/4);
	double W=10.0*A*(t+7.0*t*t/20.0+t*t*t/9.0)/(1.0+A)-12.0*C*(t+t*t*t/9.0);
	double X=4.0*D*(1.0+12.0*pi*er*eo*kb*temperature*rCo/(ec*ec*zCo*zCo*abs(zCo)))*(1.0-exp(-dZ/2.0));
	double Y=8.0*A*B/((1.0+A)*(1.0+A))+6.0*dZ*(4.0*pi*er*eo*kb*temperature*rCo*D/(ec*ec*zCo*zCo*abs(zCo))+4.0*pi*er*eo*kb*temperature*rCtr*B/(ec*ec*zCtr*zCtr*abs(zCtr)))/(1.0+A);
	double Z=24.0*A*(4.0*pi*er*eo*kb*temperature*rCo*D*D/(ec*ec*zCo*zCo*abs(zCo))+4.0*pi*er*eo*kb*temperature*rCtr*B*B/(ec*ec*zCtr*zCtr*abs(zCtr)*(1.0+A)))/(1.0+A);

	double output=1-2*A*B/(dZ*(1+A))+(W-X+Y-Z)/(dZ*ka);
	return 2*output/3.0;}

// Input: ion_val = sign charge of ion
// Input: ion_radius = ion radius [meter]
// Input: viscosity [Newton*second/meter^2]
// Output: mobility [meter^2/(Volt*second)]
double calc_Ion_Mobility(int ion_val,double ion_radius,double viscosity){return ion_val*ec/(6*pi*viscosity*ion_radius);}

// Input: ion_val = sign charge of ion
// Input: ion_mobility = electrophoretic mobility of ion in m^2/Vs (calculated by calc_Ion_Mobility)
// Output: friction (drag force) coefficient in Ns/m
// Reference: Micro- and Nanoscale Fliud Mechanics by Kirby p. 253
double calc_Ion_Friction_Coefficient(int ion_val,double ion_mobility){return ion_val*ec/ion_mobility;}

void convert_pqr_to_pdb(string saveAsPdbFile,string pqrFile,string pdbFile)
	{//string savePdb=CURRENT_FOLDER; savePdb+="saveAsPdb.py"; //get_current_dir_name(); 
	string cmd="pymol -d \"run "+saveAsPdbFile+"; saveAsPdb "+pqrFile+", "+pdbFile+"; quit;\"";
	int statusCode=system(cmd.c_str());
	}

void collect_all_files_in_dir(string dir2Search,string &fileLst,string delimiter,int &N)
	{int numFldrs=0,numFiles=0;
	string *theFldrs=get_dirs_in_dir(dir2Search,numFldrs);
	string *theFiles=get_files_in_dir(dir2Search,numFiles);
	if(numFiles>0)
		{// Add Files Found in Directory to Output File List
		for(int i=0;i<numFiles;i++){fileLst+=dir2Search+theFiles[i]+delimiter; N++;}
		}
	if(numFldrs>0)
		{// Recursively loop through sub-folders
		for(int i=0;i<numFldrs;i++){collect_all_files_in_dir(dir2Search+theFldrs[i]+"/",fileLst,delimiter,N);}
		}
	}

void determine_dx_file_name_components(string fNm,string &pHVal,string &tVal,string &srVal)
	{int pos,pos2;
	pos=fNm.find("TIGERTOWN",0);
	if(pos!=string::npos){fNm.replace(pos,9,"");}

	pos=fNm.find("pH",0);
	if(pos==string::npos)
		{cerr<<"Error in determine_dx_file_name_components!Could not determine pH value\n";}
	else
		{pos2=fNm.find("T",pos);
		pHVal=fNm.substr(pos+2,pos2-pos-2);
		}
	pos=fNm.find("T",pos);
	if(pos==string::npos)
		{cerr<<"Error in determine_dx_file_name_components!Could not determine temperature\n";}
	else
		{pos2=fNm.find("SR",pos);
		if(pos2!=string::npos)
			{tVal=fNm.substr(pos+1,pos2-pos-1);
			pos=fNm.find("SR",pos);
			if(pos==string::npos)
				{cerr<<"Error in determine_dx_file_name_components!Could not determine ion radius\n";}
			else
				{srVal=fNm.substr(pos+2,fNm.find(".",pos)-pos+1);
				}
			}
		else
			{pos2=fNm.find("S",pos);
			if(pos2!=string::npos)
				{tVal=fNm.substr(pos+1,pos2-pos-1);}
			else
				{tVal=fNm.substr(pos+1,fNm.length()-pos-1);}
			}
		}

	//cout<<"|"<<pHVal<<"|"<<tVal<<"|"<<sRVal<<"|\n";
	}

string* get_all_files_in_dir(string dir2Search,int &N)
	{string theLst="",delimiter="\\";
	string tmp=dir2Search.substr(dir2Search.length()-1,1);
	if(tmp.compare("/")!=0){dir2Search+="/";}
	//cout<<dir2Search<<endl;
	collect_all_files_in_dir(dir2Search,theLst,delimiter,N);
	string *Output=fill_string_array(theLst,N,delimiter);
	sort_array(Output,N);
	return Output;}

string get_containing_folder(string fPath)
	{return fPath.substr(0,fPath.rfind("/")+1);}

string* get_dirs_in_dir(string fPath,int &N)
	{string lst="",*Output,tmp,dirPath,delimiter="\\";
	N=0;
	struct dirent *entry;
	struct stat sa;
	DIR *dir;
	int status;
	if(directory_exist(fPath))
		{dir=opendir(fPath.c_str());
		if(dir==NULL){cerr<<"Error in get_dirs_in_dir!\nCould not open directory"<<fPath<<"\n"; exit(EXIT_FAILURE);}
		while((entry=readdir(dir))!=NULL)
			{tmp=entry->d_name;
			dirPath=fPath+tmp;
			status=stat(dirPath.c_str(),&sa);
			if(status==0)
				{if(S_ISDIR(sa.st_mode))
					{tmp=entry->d_name;
					if(tmp.compare(".")!=0 && tmp.compare("..")!=0)
						{// Build List
						lst+=tmp+delimiter;
						}
					}
				}
			}
		closedir(dir);
		N=count_delimiter(lst,delimiter);
		Output=fill_string_array(lst,N,delimiter);
		}
	return Output;}

string get_file_name(string fPath)
	{return fPath.substr(fPath.rfind("/")+1,fPath.length()-fPath.rfind("/")-1);}

string* get_files_in_dir(string fPath,int &N)
	{string lst="",*Output,tmp,filePath,delimiter="\\";
	N=0;
	struct dirent *entry;
	struct stat sa;
	DIR *dir;
	int status;
	if(directory_exist(fPath))
		{dir=opendir(fPath.c_str());
		if(dir==NULL){cerr<<"Error in get_files_in_dir!\nCould not open directory:\n"<<fPath<<"\n"; exit(EXIT_FAILURE);}
		while((entry=readdir(dir))!=NULL)
			{tmp=entry->d_name;
			filePath=fPath+tmp;
			status=stat(filePath.c_str(),&sa);
			if(status==0)
				{if(S_ISREG(sa.st_mode))
					{tmp=entry->d_name;
					lst+=tmp+delimiter;
					}
				}
			}
		closedir(dir);
		N=count_delimiter(lst,delimiter);
		Output=fill_string_array(lst,N,delimiter);
		}
	return Output;}

double interpolate(double x,double x1,double x2,double y1,double y2){return (y2-y1)*(x-x1)/(x2-x1)+y1;}

string array2String(int* Data,int numPnts,string delimiter){string Output="";for(int i=0;i<numPnts;i++){Output+=cnvrtNumToStrng(Data[i],0)+delimiter;}return Output;}

bool IN_LIST(string X,string list,int N,string delimiter)
	{// Is X in list?
	bool Output=false;
	string* L,tmp;
	if(N==0){Output=false;}
	else{L=fill_string_array(list,N,delimiter);
		for(int i=0;i<N;i++)
			{tmp=L[i];
			if(tmp.compare(X)==0){Output=true;break;}}}
	return Output;}

bool IN_STRING_ARRAY(string X,string* L,int N)
	{// Is X in array?
	bool Output=false;
	string tmp;
	if(N==0){Output=false;}
	else{for(int i=0;i<N;i++)
			{tmp=L[i];
			if(tmp.compare(X)==0){Output=true;break;}}}
	return Output;}

/*void copyFile(string inFile,string outFile)
	{long fileSz=getFileSize(inFile);
	char* buf=new char[fileSz];
	ifstream fIn;
	bool ATTEMPTING=true;
	int numTries=100000,Counter=0;
	int Err;char errMsg[256];
	while(ATTEMPTING)
		{fIn.open(inFile.c_str(),ios::in|ios::binary);
		if(!fIn.fail()){fIn.read(buf,fileSz);fIn.close();ATTEMPTING=false;break;}
		else{Counter++;}
		if(Counter>=numTries)
			{Err=errno;erase_display();
			cerr<<"ERROR in copyFile!\nInput file could not be opened.\n"<<inFile<<endl;
			cerr<<strerror_r(Err,errMsg,256)<<"\n";
			exit(EXIT_FAILURE);}
		}
	
	ATTEMPTING=true;Counter=0;
	ofstream fOut;
	while(ATTEMPTING)
		{fOut.open(outFile.c_str(),ios::out|ios::binary);
		if(!fOut.fail()){fOut.write(buf,fileSz);fOut.close();ATTEMPTING=false;break;}
		else{Counter++;}
		if(Counter>=numTries)
			{Err=errno;erase_display();
			cerr<<"ERROR in copyFile!\nOutput file could not be opened.\n"<<outFile<<endl;
			cerr<<strerror_r(Err,errMsg,256)<<"\n";
			exit(EXIT_FAILURE);}
		}
	}
*/
bool copyFile(string srcFile,string destFile)
	{ifstream theSource(srcFile.c_str(),ifstream::in|ios::binary);
	ofstream theDestination(destFile.c_str(),ofstream::out|ofstream::trunc|ios::binary);
	theDestination<<theSource.rdbuf();
	return theSource && theDestination;}

string checkFilePathParantheses(string f)
	{string Output=f;
	int pos=Output.find("(",0);
	bool SEARCHING=true;
	while(SEARCHING)
		{if(pos!=string::npos)
			{Output.replace(pos,1,"\\(");
			pos=Output.find(")",pos+1);
			if(pos!=string::npos){Output.replace(pos,1,"\\)");}
			else{cerr<<"Error in checkFilePathParantheses!\nUnbalanced parantheses present.\n";exit(EXIT_FAILURE);}
			pos=Output.find("(",pos);}
		else{SEARCHING=false;break;}}
	return Output;}

string getBaseFolder(string f)
	{int pos=f.rfind("/",f.length()-1);
	return f.substr(0,pos)+"/";}

void chmod_x(string binFile)
	{int Res=chmod(binFile.c_str(),S_IRWXU);
	int Err;
	char errMsg[256];
	if(Res!=0)
		{// An Error has occured, capture error number with erno
		Err=errno;erase_display();
		cerr<<"ERROR in chmod_x while opening|"<<binFile<<"|END\n"<<strerror_r(Err,errMsg,256)<<"\n";
		exit(EXIT_FAILURE);}}

string nameBuffer(string soluteName,string soluteConc)
	{string Output="";
	int N=count_delimiter(soluteName,";");
	string* name=fill_string_array(soluteName,N,";");
	string* conc=fill_string_array(soluteConc,N,";");
	int pos;
	for(int i=0;i<N;i++)
		{soluteConc=formatNumberString(conc[i]);
		if(i==N-1){Output+=soluteConc+"M_"+name[i];}
		else{Output+=soluteConc+"M_"+name[i]+"_";}}
	delete [] name; delete [] conc;
	return Output;}

void sort_array(string *theArr,int N)
	{bool RUNNING=true,SORTED;
	string tmp;
	while(RUNNING)
		{SORTED=true;
		for(int i=0;i<N-1;i++)
			{if(theArr[i+1]<theArr[i])
				{SORTED=false;
				tmp=theArr[i+1];
				theArr[i+1]=theArr[i];
				theArr[i]=tmp;
				}
			}
		if(SORTED){RUNNING=false; break;}
		}
	}

int countNumLinesMSMSVertFile(string fileName)
	{float x,y,z,nx,ny,nz;
	int analyticSurf,sphereNum,vertType,Counter=0;
	char closestAtom[15000];

	string msms_vertices_file_format="%f %f %f %f %f %f %d %d %d %s";

	FILE *fIn;
    fIn=fopen(fileName.c_str(),"r");
    if(fIn==NULL)
        {cerr<<"ERROR in countNumLinesMSMSVertFile!\nInput Vertices File could not be found.\n\nFilepath:\t"<<fileName<<endl;
         exit(EXIT_FAILURE);}

	fseek(fIn,0,SEEK_END);
	long fileSize=ftell(fIn);
	rewind(fIn);

	int lineSz;
	while(!feof(fIn))
        {fscanf(fIn,msms_vertices_file_format.c_str(),&x,&y,&z,&nx,&ny,&nz,&analyticSurf,&sphereNum,&vertType,&closestAtom);		
		lineSz=ftell(fIn);
		if(lineSz==fileSize){break;}
         else{Counter++;}}
	fclose(fIn);
	return Counter;}

string cnvrtNumToStrng(int Num,int numberAfterDecimalpoint){stringstream ss;ss.setf(ios::fixed);if(numberAfterDecimalpoint>0){ss.setf(ios::showpoint);}ss.precision(numberAfterDecimalpoint);ss<<Num;return ss.str();}

string cnvrtNumToStrng(double Num,int numberAfterDecimalpoint)
	{stringstream ss;
	//ss.setf(ios::scientific); //ss.setf(ios::fixed);
	double testVal=log10(Num);
	double tmpVal=BUFFER_SIZE;
	tmpVal-=2;
	tmpVal*=-1;
	//if(numberAfterDecimalpoint>0){ss.setf(ios::showpoint);}
	//cout<<testVal<<"|"<<tmpVal<<endl;
	if(testVal>tmpVal)
		{ss.setf(ios::fixed);
		if(numberAfterDecimalpoint>0){ss.setf(ios::showpoint);}
		ss.precision(numberAfterDecimalpoint);}
	else
		{ss.setf(ios::scientific);}
	ss<<Num;
	return ss.str();}

string cnvrtNumToStrng(long double Num,int numberAfterDecimalpoint)
	{stringstream ss;
	double testVal=log10(Num);
	double tmpVal=BUFFER_SIZE;
	tmpVal-=2;
	tmpVal*=-1;
	if(testVal>tmpVal)
		{ss.setf(ios::fixed);
		if(numberAfterDecimalpoint>0){ss.setf(ios::showpoint);}
		ss.precision(numberAfterDecimalpoint);}
	else
		{ss.setf(ios::scientific);}
	ss<<Num;
	return ss.str();}

string cnvrtNumToStrng(float Num,int numberAfterDecimalpoint)
	{stringstream ss;
	double testVal=log10(Num);
	double tmpVal=BUFFER_SIZE;
	tmpVal-=2;
	tmpVal*=-1;
	if(testVal>tmpVal)
		{ss.setf(ios::fixed);
		if(numberAfterDecimalpoint>0){ss.setf(ios::showpoint);}
		ss.precision(numberAfterDecimalpoint);}
	else
		{ss.setf(ios::scientific);}
	ss<<Num;
	return ss.str();}

string* cnvrtzInput2Str(zInput In,string delimiter,string zFldr)
	{// Determine Number of Input Strings to Create
	int totalCalc=In.numFiles*In.numSolCond;
	string* Output=new string[totalCalc];
	string solvStr="",solvConcStr="",solvConcTypeStr="",soluStr="",soluConcStr="";
	int Counter;
	for(int p=0;p<In.numFiles;p++)
		{for(int j=0;j<In.numSolCond;j++)
			{Counter=p*In.numSolCond+j;
			// Define ZPRED Input String
			// Solution Condition Calculation Number (1)
			//Output[Counter]+=In.id[j]+delimiter;
			Output[Counter]=cnvrtNumToStrng(j,0)+delimiter;
			// PDB File Path (2)
			Output[Counter]+=In.files[p]+delimiter;
			// Job Name/Output Tile (3)
			Output[Counter]+=In.outputTitle+delimiter;
			// pH Value (4)
			Output[Counter]+=In.pH[j]+delimiter;
			// Temperature Value (5)
			Output[Counter]+=formatNumberString(cnvrtNumToStrng(In.temperature[j],BUFFER_SIZE))+delimiter;
			// Number of Solvents (6)
			Output[Counter]+=cnvrtNumToStrng(In.numSolvent[j],0)+delimiter;
			solvStr="";solvConcStr="";
			for(int i=0;i<In.numSolvent[j];i++)
				{solvStr+=In.solvent[j][i]+delimiter;
				solvConcStr+=formatNumberString(cnvrtNumToStrng(In.solventConc[j][i],BUFFER_SIZE))+delimiter;}
			// Solvent Name String	(already delimited) (7)
			Output[Counter]+=solvStr;
			// Solvent Concentration String (already delimited) (8)
			Output[Counter]+=solvConcStr;
			// Solvent Concentration Type String (9)
			Output[Counter]+=In.solventConcType[j]+delimiter;
			// Number of Solutes (10)
			Output[Counter]+=cnvrtNumToStrng(In.numSolute[j],0)+delimiter;
			soluStr="";soluConcStr="";
			for(int i=0;i<In.numSolute[j];i++)
				{soluStr+=In.solute[j][i]+delimiter;
				soluConcStr+=formatNumberString(cnvrtNumToStrng(In.soluteConc[j][i],BUFFER_SIZE))+delimiter;}
			// Solute Name String	(already delimited) (11)
			Output[Counter]+=soluStr;
			// Solute Concentration String (already delimited) (12)
			Output[Counter]+=soluConcStr;
			// Protein Concentration (13)
			Output[Counter]+=formatNumberString(cnvrtNumToStrng(In.proteinConc[j],BUFFER_SIZE))+delimiter;
			// Number of APBS Write Types (14)
			Output[Counter]+=cnvrtNumToStrng(In.numApbsWriteTypes,0)+delimiter;
			// APBS Write Types (15)
			for(int f=0;f<In.numApbsWriteTypes;f++){Output[Counter]+=In.apbsWriteTypes[f]+delimiter;}
			// Executable & Program File Paths (16)
			Output[Counter]+=In.SC.apbsExecutable+delimiter;
			Output[Counter]+=In.SC.multivalueExecutable+delimiter;	// (17)
			Output[Counter]+=In.SC.hydroproExecutable+delimiter;	// (18)
			Output[Counter]+=In.SC.msmsExecutable+delimiter;		// (19)
			Output[Counter]+=In.SC.msmsAtmTypeNumbers+delimiter;	// (20)
			Output[Counter]+=In.SC.msmsPdbToXyzr+delimiter;			// (21)
			Output[Counter]+=In.SC.msmsPdbToXyzrn+delimiter;		// (22)
			Output[Counter]+=In.SC.pdb2pqr+delimiter;				// (23)
			// MSMS Parameter: Low Surface Density (24)
			Output[Counter]+=In.msms.surfPointDensity+delimiter;
			// MSMS Parameter: High Surface Density (25)
			Output[Counter]+=In.msms.surfPointHiDensity+delimiter;
			// HYDROPRO Parameter: Calculation Type (26)
			Output[Counter]+=In.hydroproCalcType+delimiter;
			// PDB2PQR Parameter: Force Field (27)
			Output[Counter]+=In.forceField+delimiter;
			// ZPRED Parameter: Flag for Saving Calculation Files (28)
			if(In.saveTemps){Output[Counter]+="TRUE"+delimiter;}
			else{Output[Counter]+="FALSE"+delimiter;}
			// Protein Dielectric Constant (29)
			Output[Counter]+=In.pdie+delimiter;
			// APBS Parameter: X-Axis Grid Dimensions (30)
			Output[Counter]+=In.dime_x+delimiter;
			// APBS Parameter: Y-Axis Grid Dimensions (31)
			Output[Counter]+=In.dime_y+delimiter;
			// APBS Parameter: Z-Axis Grid Dimensions (32)
			Output[Counter]+=In.dime_z+delimiter;
			// ZPRED Parameter: APBS Control (33)
			if(In.RUN_APBS){Output[Counter]+="true"+delimiter;}
			else{Output[Counter]+="false"+delimiter;}
			// ZPRED Parameter: PDB2QPR Control (34)
			if(In.RUN_PDB2PQR){Output[Counter]+="true"+delimiter;}
			else{Output[Counter]+="false"+delimiter;}
			// ZPRED Parameter: HYDROPRO Control (35)
			if(In.RUN_HYDROPRO){Output[Counter]+="true"+delimiter;}
			else{Output[Counter]+="false"+delimiter;}
			// ZPRED Parameter: HULLRAD Control (36)
			if(In.RUN_HULLRAD){Output[Counter]+="true"+delimiter;}
			else{Output[Counter]+="false"+delimiter;}
			// ZPRED Parameter: Final Output Control (37)
			if(In.RUN_ZPRED){Output[Counter]+="true"+delimiter;}
			else{Output[Counter]+="false"+delimiter;}
			// ZPRED Parameter: COLLIDE Control (38)
			if(In.RUN_COLLIDE){Output[Counter]+="true"+delimiter;}
			else{Output[Counter]+="false"+delimiter;}					
			// ZPRED Base Folder (39)
			Output[Counter]+=zFldr+delimiter;
			// User-Specified Solution Density [kg/L] (40)
			Output[Counter]+=cnvrtNumToStrng(In.density[j],BUFFER_SIZE)+delimiter;
			// User-Specified Solution Relative Dielectric (41)
			Output[Counter]+=cnvrtNumToStrng(In.dielectric[j],ER_SIG_FIGS)+delimiter;
			// User-Specified Solution Inflation Distance/Slip Plane Position [A] (42)
			Output[Counter]+=cnvrtNumToStrng(In.inflationDistance[j],DISTANCE_SIG_FIGS)+delimiter;
			// User-Specified Solution Viscosity [Pa s] (43)
			Output[Counter]+=cnvrtNumToStrng(In.viscosity[j],BUFFER_SIZE)+delimiter;
			// Generate Potential Profile (44)
			if(In.GENERATE_POT_PROFILE[j]){Output[Counter]+="true"+delimiter;}
			else{Output[Counter]+="false"+delimiter;}
			
			Counter++;
			}
		}
	if(Counter!=totalCalc){cerr<<"ERROR in cnvrtzInput2Str\nNumber of Calculations does not match the number of strings written"<<endl;exit(EXIT_FAILURE);}
	return Output;}

zExeInput cnvrtStr2Input(string cmd,string delimiter,string& Fldr,string& calcNum,string& display)
	{zExeInput Out;	
	string tmp="";
	int pos,oldPos;
	pos=cmd.find(delimiter,0); Out.id=cmd.substr(0,pos); oldPos=pos+1; oldPos=pos+1;									// Solution Condition Calculation Number (1)
	pos=cmd.find(delimiter,pos+1); Out.file=cmd.substr(oldPos,pos-oldPos); oldPos=pos+1;									// PDB File Path (2)
	pos=cmd.find(delimiter,pos+1); Out.outputTitle=cmd.substr(oldPos,pos-oldPos); oldPos=pos+1;							// Job Name/Output Tile (3)
	pos=cmd.find(delimiter,pos+1); Out.pH=cmd.substr(oldPos,pos-oldPos); oldPos=pos+1;									// pH Value (4)
	pos=cmd.find(delimiter,pos+1); tmp=cmd.substr(oldPos,pos-oldPos); Out.temperature=strtod(tmp.c_str(),NULL); oldPos=pos+1;	// Temperature Value (5)
	pos=cmd.find(delimiter,pos+1); tmp=cmd.substr(oldPos,pos-oldPos); Out.numSolvent=atoi(tmp.c_str());						// Number of Solvents (6)
	// Solvent Name String	(already delimited) (7)
	// Solvent Concentration String (already delimited) (8)
	Out.solvent=new string[Out.numSolvent];
	Out.solventConc=new double[Out.numSolvent];
	// Extract Solvent(s) Name(s) (7)
	for(int i=0;i<Out.numSolvent;i++)
		{oldPos=pos+1;
		pos=cmd.find(delimiter,pos+1);
		Out.solvent[i]=cmd.substr(oldPos,pos-oldPos);}
	// Extract Solvent(s) Concentration(s) (8)
	for(int i=0;i<Out.numSolvent;i++)
		{oldPos=pos+1;
		pos=cmd.find(delimiter,pos+1);
		tmp=cmd.substr(oldPos,pos-oldPos);
		Out.solventConc[i]=strtod(tmp.c_str(),NULL);}
	oldPos=pos+1;
	pos=cmd.find(delimiter,pos+1); Out.solventConcType=cmd.substr(oldPos,pos-oldPos); oldPos=pos+1;						// Solvent Concentration Type (9)	
	pos=cmd.find(delimiter,pos+1); tmp=cmd.substr(oldPos,pos-oldPos); Out.numSolute=atoi(tmp.c_str());						// Number of Solutes (10)
	// Solute Name String	(already delimited) (11)
	// Solute Concentration String (already delimited) (12)
	Out.solute=new string[Out.numSolute];
	Out.soluteConc=new double[Out.numSolute];
	// Extract Solute Name(s) (11)
	for(int i=0;i<Out.numSolute;i++)
		{oldPos=pos+1;
		pos=cmd.find(delimiter,pos+1);
		Out.solute[i]=cmd.substr(oldPos,pos-oldPos);}
	// Extract Solute(s) Concentration(s) (12)
	for(int i=0;i<Out.numSolute;i++)
		{oldPos=pos+1;
		pos=cmd.find(delimiter,pos+1);
		tmp=cmd.substr(oldPos,pos-oldPos);
		Out.soluteConc[i]=strtod(tmp.c_str(),NULL);}
	oldPos=pos+1;
	pos=cmd.find(delimiter,pos+1); tmp=cmd.substr(oldPos,pos-oldPos); Out.proteinConc=strtod(tmp.c_str(),NULL); oldPos=pos+1;	// Protein Concentration (13)
	pos=cmd.find(delimiter,pos+1); tmp=cmd.substr(oldPos,pos-oldPos); Out.numApbsWriteTypes=atoi(tmp.c_str());				// Number of APBS Write Types (14)
	// APBS Write Types (15)
	Out.apbsWriteTypes=new string[Out.numApbsWriteTypes];
	for(int i=0;i<Out.numApbsWriteTypes;i++)
		{oldPos=pos+1;
		pos=cmd.find(delimiter,pos+1);
		Out.apbsWriteTypes[i]=cmd.substr(oldPos,pos-oldPos);}
	oldPos=pos+1;	
	pos=cmd.find(delimiter,pos+1); Out.SC.apbsExecutable=cmd.substr(oldPos,pos-oldPos); oldPos=pos+1;						// APBS EXECUTABLE (16)
	pos=cmd.find(delimiter,pos+1); Out.SC.multivalueExecutable=cmd.substr(oldPos,pos-oldPos); oldPos=pos+1;					// MULTIVALUE Executable (17)
	pos=cmd.find(delimiter,pos+1); Out.SC.hydroproExecutable=cmd.substr(oldPos,pos-oldPos); oldPos=pos+1;					// HYDROPRO Executable (18)
	pos=cmd.find(delimiter,pos+1); Out.SC.msmsExecutable=cmd.substr(oldPos,pos-oldPos); oldPos=pos+1;						// MSMS Executable (19)
	pos=cmd.find(delimiter,pos+1); Out.SC.msmsAtmTypeNumbers=cmd.substr(oldPos,pos-oldPos); oldPos=pos+1;					// MSMS AtmTypeNumbers File (20)
	pos=cmd.find(delimiter,pos+1); Out.SC.msmsPdbToXyzr=cmd.substr(oldPos,pos-oldPos); oldPos=pos+1;						// MSMS pdb_to_xyzr File (21)
	pos=cmd.find(delimiter,pos+1); Out.SC.msmsPdbToXyzrn=cmd.substr(oldPos,pos-oldPos); oldPos=pos+1;						// MSMS pdb_to_xyzrn File (22)
	pos=cmd.find(delimiter,pos+1); Out.SC.pdb2pqr=cmd.substr(oldPos,pos-oldPos); oldPos=pos+1;							// PDB2PQR main.py File (23)
	pos=cmd.find(delimiter,pos+1); Out.msms.surfPointDensity=cmd.substr(oldPos,pos-oldPos); oldPos=pos+1;					// MSMS Parameter: Low Surface Density (24)
	pos=cmd.find(delimiter,pos+1); Out.msms.surfPointHiDensity=cmd.substr(oldPos,pos-oldPos); oldPos=pos+1;					// MSMS Parameter: High Surface Density (25)
	pos=cmd.find(delimiter,pos+1); Out.hydroproCalcType=cmd.substr(oldPos,pos-oldPos); oldPos=pos+1;						// HYDROPRO Parameter: Calculation Type (26)
	pos=cmd.find(delimiter,pos+1); Out.forceField=cmd.substr(oldPos,pos-oldPos); oldPos=pos+1;							// PDB2PQR Parameter: Force Field (27)
	pos=cmd.find(delimiter,pos+1); tmp=cmd.substr(oldPos,pos-oldPos);												// ZPRED Parameter: saveTemps Flag for Saving Calculation Files (28)
	if(tmp.compare("FALSE")==0){Out.saveTemps=false;}
	else{Out.saveTemps=true;}
	oldPos=pos+1;
	pos=cmd.find(delimiter,pos+1); Out.pdie=cmd.substr(oldPos,pos-oldPos); oldPos=pos+1;									// Protein Dielectric Constant (29)
	pos=cmd.find(delimiter,pos+1); Out.dime_x=cmd.substr(oldPos,pos-oldPos); oldPos=pos+1;								// APBS Parameter: X-Axis Grid Dimensions (30)
	pos=cmd.find(delimiter,pos+1); Out.dime_y=cmd.substr(oldPos,pos-oldPos); oldPos=pos+1;								// APBS Parameter: Y-Axis Grid Dimensions (31)
	pos=cmd.find(delimiter,pos+1); Out.dime_z=cmd.substr(oldPos,pos-oldPos); oldPos=pos+1;								// APBS Parameter: Z-Axis Grid Dimensions (32)
	pos=cmd.find(delimiter,pos+1); tmp=cmd.substr(oldPos,pos-oldPos);												// ZPRED Parameter: APBS Control (33)
	if(tmp.compare("false")==0){Out.RUN_APBS=false;}
	else{Out.RUN_APBS=true;}
	oldPos=pos+1;
	pos=cmd.find(delimiter,pos+1); tmp=cmd.substr(oldPos,pos-oldPos);												// ZPRED Parameter: PDB2QPR Control (34)
	if(tmp.compare("false")==0){Out.RUN_PDB2PQR=false;}
	else{Out.RUN_PDB2PQR=true;}
	oldPos=pos+1;
	pos=cmd.find(delimiter,pos+1); tmp=cmd.substr(oldPos,pos-oldPos);												// ZPRED Parameter: HYDROPRO Control (35)
	if(tmp.compare("false")==0){Out.RUN_HYDROPRO=false;}
	else{Out.RUN_HYDROPRO=true;}
	oldPos=pos+1;
	pos=cmd.find(delimiter,pos+1); tmp=cmd.substr(oldPos,pos-oldPos);												// ZPRED Parameter: HULLRAD Control (36)
	if(tmp.compare("false")==0){Out.RUN_HULLRAD=false;}
	else{Out.RUN_HULLRAD=true;}
	oldPos=pos+1;
	pos=cmd.find(delimiter,pos+1); tmp=cmd.substr(oldPos,pos-oldPos);												// ZPRED Parameter: Final Output Control (37)
	if(tmp.compare("false")==0){Out.RUN_ZPRED=false;}
	else{Out.RUN_ZPRED=true;}
	oldPos=pos+1;
	pos=cmd.find(delimiter,pos+1); tmp=cmd.substr(oldPos,pos-oldPos);												// ZPRED Parameter: COLLIDE Control (38)
	if(tmp.compare("false")==0){Out.RUN_COLLIDE=false;}
	else{Out.RUN_COLLIDE=true;}
	oldPos=pos+1;
	pos=cmd.find(delimiter,pos+1); Fldr=cmd.substr(oldPos,pos-oldPos); oldPos=pos+1;									// ZPRED Base Folder (39)
	pos=cmd.find(delimiter,pos+1); tmp=cmd.substr(oldPos,pos-oldPos); Out.density=strtod(tmp.c_str(),NULL); oldPos=pos+1;		// User-Specified Solution Density [kg/L] (40)
	pos=cmd.find(delimiter,pos+1); tmp=cmd.substr(oldPos,pos-oldPos); Out.dielectric=strtod(tmp.c_str(),NULL); oldPos=pos+1;	// User-Specified Solution Relative Dielectric (41)
	pos=cmd.find(delimiter,pos+1); tmp=cmd.substr(oldPos,pos-oldPos); Out.inflationDistance=strtod(tmp.c_str(),NULL); oldPos=pos+1;// User-Specified Solution Inflation Distance/Slip Plane Position [A] (42)
	pos=cmd.find(delimiter,pos+1); tmp=cmd.substr(oldPos,pos-oldPos); Out.viscosity=strtod(tmp.c_str(),NULL); oldPos=pos+1;		// User-Specified Solution Viscosity [Pa s] (43)
	// Generate Electric Potential Profile? (44)	
	pos=cmd.find(delimiter,pos+1); tmp=cmd.substr(oldPos,pos-oldPos);
	if(tmp.compare("true")==0){Out.GENERATE_POT_PROFILE=true;}
	else if(tmp.compare("false")==0){Out.GENERATE_POT_PROFILE=false;}
	oldPos=pos+1;	
	pos=cmd.find(delimiter,pos+1); calcNum=cmd.substr(oldPos,pos-oldPos); oldPos=pos+1;									// ZPRED Calculation Number (45)
	pos=cmd.find(delimiter,pos+1); display=cmd.substr(oldPos,pos-oldPos);												// Computation Progress Display String (46)
	return Out;}

int count_delimiter(string Data,string delimiter)
	{int Counter=0,rng=delimiter.length()-1;
	string tmp="";
	for(int i=0;i<Data.length()-rng;i++)
		{tmp=Data[i];
		for(int j=1;j<rng+1;j++){tmp+=Data[i+j];}
		if(tmp.compare(delimiter)==0)
			{Counter++;}}
	return Counter;}

string checkWhiteSpaceInFilePath(string iFile)
	{int pos=iFile.find(" ",0),pos2,pos3;
	bool CHECKING=true;
	string tmp=iFile,tmp2;
	while(CHECKING)
		{if(pos!=string::npos)
			{// White Space Found, Add ' ' to folder Name
			// 1st /
			pos2=tmp.rfind("/",pos);
			// 2nd /
			pos3=tmp.find("/",pos);
			// Add Single Quotes to tmp
			tmp2=tmp.substr(0,pos2+1)+"\'"+tmp.substr(pos2+1,pos3-pos2-1)+"\'"+tmp.substr(pos3,tmp.length()-pos3);
			tmp=tmp2;
			//inputVal=inputVal.substr(1,inputVal.length()-1);pos=inputVal.find(" ",0);
			pos=tmp.find(" ",pos+3);
			}
		else if(pos==0)
			{// Whitespace Preceeding filepath, remove it
			tmp=tmp.substr(1,tmp.length()-1);
			pos=tmp.find(" ",0);
			}
		else{CHECKING=false;break;}
		}
	return tmp;}

void define_msmsInput(string PDB,string pdbFile,string msmsFldr,string pdb_msmsFldr,string pdb_surfFldr,string surfPointDensity,string surfPointHiDensity,string probeRad,string msmsBinFile,msmsInput& In)
	{// Current PDB File Name under assessment
	In.PDB=PDB;
	// File Path of Current PDB File Name under assessment
	In.pdbFile=checkWhiteSpaceInFilePath(pdbFile);
	// MSMS Computation Folder
	In.msmsFldr=checkWhiteSpaceInFilePath(msmsFldr);
	// File Path of MSMS Computation Folder
	In.pdb_msmsFldr=checkWhiteSpaceInFilePath(pdb_msmsFldr);
	// File Path of MSMS Output Folder	
	In.pdb_surfFldr=checkWhiteSpaceInFilePath(pdb_surfFldr);
	// MSMS Low Surface Density
	In.surfPointDensity=surfPointDensity;
	// MSMS High Surface Density
	In.surfPointHiDensity=surfPointHiDensity;
	// MSMS Spherical Probe Radius
	In.probeRadius=probeRad;
	// MSMS PDB_TO_XYZRN
	In.PDB_TO_XYZRN="pdb_to_xyzrn";
	// MSMS Executable File (file path is not needed)
	In.MSMS=checkWhiteSpaceInFilePath(getFileName(msmsBinFile));}

bool directory_exist(string fPath)
	{struct stat sb;
	if(stat(fPath.c_str(),&sb)==0 && S_ISDIR(sb.st_mode)){return true;}
	else{return false;}}

void erase_display(){cerr<<"\033[2J\033[1;H";}

string encryptFile(string inFile,const string salt,const string pass)
	{string s=readFile2Encrypt(inFile);	
	// Generate Key and IV from Password and Salt
	//const string pass="Password";
	unsigned char key[32], iv[32];
	// Allocate Space for Encrypted (Ciphered) text
	unsigned char ciphertext[s.length()];
  	// Buffer for the decrypted text
	unsigned char decryptedtext[s.length()];
	int keyLength=EVP_BytesToKey(EVP_aes_256_cbc(),EVP_sha1(),(unsigned char*)salt.c_str(),(unsigned char*)pass.c_str(),pass.length(),1,key,iv);
	//cout<<"Key Length:\n"<<keyLength<<"\nKey:\n"<<key<<endl;
	// Define Cipher Object
	EVP_CIPHER_CTX *ctx;
	int len,ret;
	// Create and initialize the context
	if(!(ctx=EVP_CIPHER_CTX_new())){handleErrors();}
	// Start Encryption
	if(1 != EVP_EncryptInit_ex(ctx,EVP_aes_256_cbc(),NULL,key,iv)){handleErrors();}	
	// Encrypt message and obtain output	
	if(1 != EVP_EncryptUpdate(ctx,ciphertext,&len,(unsigned char*)s.c_str(),s.length())){handleErrors();}
	int ciphertext_len=len;
	//cout<<ciphertext<<endl;
	// Finalize Encryption
	if(1 != EVP_EncryptFinal_ex(ctx,ciphertext+len,&len)){handleErrors();}
	ciphertext_len+=len;
	// Clean Up
	EVP_CIPHER_CTX_free(ctx);
	// Check ciphertext
	//cout<<"Ciphered text is:\n";
	//BIO_dump_fp(stdout,(const char*)ciphertext,ciphertext_len);
	// Base64 Encoded Encryption
	//string b64enc=base64_encode(reinterpret_cast<const unsigned char*>(cVal.c_str()), cVal.length());
	string Output=base64_encode(ciphertext,ciphertext_len);
	return Output;}

bool file_exist(string fPath)
	{struct stat sb;
	if(stat(fPath.c_str(),&sb)==0 && S_ISREG(sb.st_mode)){return true;}
	else{return false;}}

string readFile2Encrypt(string inFile)
	{long fileSz=getFileSize(inFile);
	char* buf=new char[fileSz];
	ifstream fIn;
	fIn.open(inFile.c_str(),ios::in|ios::binary);
	if(fIn.fail()){cerr<<"ERROR in readFile2Encrypt!\ninput file could not be opened.\n"<<inFile<<endl;exit(EXIT_FAILURE);}
	fIn.read(buf,fileSz);
	fIn.close();
	string val="";
	for(int i=0;i<fileSz;i++){val+=buf[i];}
	delete [] buf;
	string Output=val;
	return Output;}

string getFileName(string filePath)
	{string Output;
	int pos=filePath.rfind("/");
	Output=filePath.substr(pos+1,filePath.length()-pos-1);
	return Output;}

string base64_encode(unsigned char const* bytes_to_encode, unsigned int in_len)
	{static const string base64_chars = 
		 "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
		 "abcdefghijklmnopqrstuvwxyz"
		 "0123456789+/";
	string ret;
	int i=0,j=0;
	unsigned char char_array_3[3];
	unsigned char char_array_4[4];

	while(in_len--)
		{char_array_3[i++] = *(bytes_to_encode++);
		if(i==3)
			{char_array_4[0]=(char_array_3[0] & 0xfc) >> 2;
			char_array_4[1]=((char_array_3[0] & 0x03) << 4) + ((char_array_3[1] & 0xf0) >> 4);
			char_array_4[2]=((char_array_3[1] & 0x0f) << 2) + ((char_array_3[2] & 0xc0) >> 6);
			char_array_4[3]=char_array_3[2] & 0x3f;

			for(i=0;(i<4);i++){ret+=base64_chars[char_array_4[i]];}
			i=0;}}

	if(i)
		{for(j=i;j<3;j++){char_array_3[j]='\0';}

		char_array_4[0]=( char_array_3[0] & 0xfc) >> 2;
		char_array_4[1]=((char_array_3[0] & 0x03) << 4) + ((char_array_3[1] & 0xf0) >> 4);
		char_array_4[2]=((char_array_3[1] & 0x0f) << 2) + ((char_array_3[2] & 0xc0) >> 6);

		for(j=0;(j<i+1);j++){ret+=base64_chars[char_array_4[j]];}

		while((i++ < 3)){ret+='=';}}

	return ret;}

string base64_decode(string const& encoded_string)
	{static const string base64_chars = 
             "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
             "abcdefghijklmnopqrstuvwxyz"
             "0123456789+/";
	int in_len=encoded_string.size();
	int i=0,j=0,in_=0;
	unsigned char char_array_4[4], char_array_3[3];
	string ret;

	while(in_len-- && (encoded_string[in_]!='=') && is_base64(encoded_string[in_]))
		{char_array_4[i++]=encoded_string[in_];
		in_++;
		if(i==4)
			{for(i=0;i<4;i++){char_array_4[i]=base64_chars.find(char_array_4[i]);}

			char_array_3[0]=( char_array_4[0] << 2       ) + ((char_array_4[1] & 0x30) >> 4);
			char_array_3[1]=((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);
			char_array_3[2]=((char_array_4[2] & 0x3) << 6) +   char_array_4[3];

			for(i=0;(i<3);i++){ret += char_array_3[i];}
			i=0;}}

	if(i)
		{for(j=0;j<i;j++){char_array_4[j]=base64_chars.find(char_array_4[j]);}

		char_array_3[0]=(char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
		char_array_3[1]=((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);

		for(j=0;(j<i-1);j++){ret+=char_array_3[j];}}

	return ret;}

void decryptString(string encStr,const string salt,const string pass,string outFile)
	{// Base 64 Decode String\n";
	string decoded=base64_decode(encStr);
	// Generate Key and IV from Password and Salt\n";
	unsigned char key[32], iv[32];
	// Allocate Space for decrypted text\n";
	unsigned char decryptedtext[decoded.length()];
	int keyLength=EVP_BytesToKey(EVP_aes_256_cbc(),EVP_sha1(),(unsigned char*)salt.c_str(),(unsigned char*)pass.c_str(),pass.length(),1,key,iv);
	// Define Cipher Object
	EVP_CIPHER_CTX *ctx;
	int len,ret;
	// Create and initialize the context
	if(!(ctx=EVP_CIPHER_CTX_new())){handleErrors();}
	// Start Decryption
	if(1 != EVP_DecryptInit_ex(ctx,EVP_aes_256_cbc(),NULL,key,iv)){handleErrors();}
	// Decrypt message and obtain output
	if(1 != EVP_DecryptUpdate(ctx,decryptedtext,&len,(unsigned char*)decoded.c_str(),decoded.length())){handleErrors();}
	int decryptedtext_len=len;
	// Finalize Decryption
	if(1 != EVP_DecryptFinal_ex(ctx,decryptedtext+len,&len)){handleErrors();}
	decryptedtext_len+=len;
	// Clean Up
	EVP_CIPHER_CTX_free(ctx);
	// Add a NULL terminator. We are expecting printable text
	decryptedtext[decryptedtext_len]='\0';
	ofstream fOut;
	fOut.open(outFile.c_str(),ios::out|ios::binary);
	if(fOut.fail())
		{cerr<<"ERROR:\n input file could not be opened.\n"<<outFile<<endl;
		exit(1);}
	for(int i=0;i<decryptedtext_len;i++){fOut<<decryptedtext[i];}
	fOut.close();}

static inline bool is_base64(unsigned char c){return (isalnum(c) || (c == '+') || (c == '/'));}

void handleErrors(void){ERR_print_errors_fp(stderr);abort();}

double* fill_double_array(string Data,int numPnts,string delimiter)
	{double* Output=new double[numPnts];
	string bld="",tmp="";
	int Counter=0,rng=delimiter.length()-1;
	for(int i=0;i<Data.length();i++)
		{tmp=Data[i];
		for(int j=1;j<rng+1;j++){tmp+=Data[i+j];}
		if(tmp.compare(delimiter)==0 && Counter<numPnts)
			{Output[Counter]=strtod(bld.c_str(),NULL);
			Counter++;
			bld="";
			i=i+rng;}
		else{bld+=Data[i];}
		}
	return Output;}

int* fill_int_array(string Data,int numPnts,string delimiter)
	{int* Output=new int[numPnts];
	string bld="",tmp="";
	int Counter=0,rng=delimiter.length()-1;
	for(int i=0;i<Data.length();i++)
		{tmp=Data[i];
		for(int j=1;j<rng+1;j++){tmp+=Data[i+j];}
		if(tmp.compare(delimiter)==0 && Counter<numPnts)
			{Output[Counter]=atoi(bld.c_str());
			Counter++;
			bld="";
			i=i+rng;}
		else{bld+=Data[i];}}
	return Output;}

string* fill_string_array(string Data,int numPnts,string delimiter)
	{string* Output=new string[numPnts];
	string bld="",tmp="";
	int Counter=0,rng=delimiter.length()-1;
	for(int i=0;i<Data.length();i++)
		{tmp=Data[i];
		for(int j=1;j<rng+1;j++){tmp+=Data[i+j];}
		if(tmp.compare(delimiter)==0 && Counter<numPnts)
			{Output[Counter]=bld;
			Counter++;
			bld="";
			i=i+rng;}
		else{bld+=Data[i];}
		}
	return Output;}

string formatNumberString(string Num)
	{string Output,tmp,tmp2;
	int i,pos;
	/* Check if Scientific Notation */
	pos=Num.find("e",0);
	if(pos==string::npos)
		{pos=Num.find("E",0);
		if(pos!=string::npos){return Num;}
		}
	else if(pos!=string::npos){return Num;}
	/* Check if contains decimal*/
	pos=Num.find(".",0);
	if(pos==string::npos)
		{/* Number Value is an integer with no decimal*/
		for(i=0;i<Num.length();i++)
			{tmp=Num[i];
			/* Find 1st Non-Zero Number*/
			if(tmp.compare("0")!=0){break;}}
		Output=Num.substr(i,Num.length()-i);}
	else if(pos==0)
		{/* Precede decimal with zero*/
		Output="0";
		for(i=Num.length();i>0;i--)
			{tmp=Num[i-1];
			/* Search in reverse for 1st Non-Zero Number*/
			if(tmp.compare("0")!=0){break;}}
		Output+=Num.substr(0,i);}
	else{/* Number holds form XXX.XXX*/
		for(i=0;i<pos;i++)
			{tmp=Num[i];
			/* Find 1st Non-Zero Number*/
			if(tmp.compare("0")!=0){break;}}
		/* Define value preceding decimal*/
		if(i==pos){/* Only Zeros found up to the decimal*/Output="0";}
		else{Output=Num.substr(i,pos-i);}
		for(i=Num.length();i>pos;i--)
			{tmp=Num[i-1];
			/* Search in reverse for 1st Non-Zero Number*/
			if(tmp.compare("0")!=0){break;}}
		tmp2=Num.substr(pos,i-pos);
		/* will contain .XXX*/
		if(tmp2.length()>1){Output+=tmp2;}}
	return Output;}

long getFileSize(string inFile)
	{FILE *fI;
	int Err;
	char errMsg[256];
	fI=fopen(inFile.c_str(),"r");
	if(fI==NULL)
		{Err=errno;erase_display();
		cerr<<"ERROR in getFileSize!\nCould not open input file: "<<inFile<<endl;		
		cerr<<strerror_r(Err,errMsg,256)<<"\n";
		exit(EXIT_FAILURE);}
	fseek(fI,0,SEEK_END);
	long Output=ftell(fI);
	fclose(fI);
	return Output;}

string initialize_ZPRED_display(int totalCalc,bool& MULTI_DISPLAY)
	{int heightBuffer=5,widthBuffer=10,tWidth=0,tHeight=0;
	string Output="";
	if(totalCalc>MAX_DISPLAY_ROWS){tHeight=MAX_DISPLAY_ROWS;}
	else{tHeight=totalCalc;}
	// Determine Number of Columns to Display
	int numCols;
	if(totalCalc>=MAX_DISPLAY_ROWS*MAX_DISPLAY_COLS){numCols=MAX_DISPLAY_COLS;MULTI_DISPLAY=true;}	
	else{numCols=totalCalc/MAX_DISPLAY_ROWS+1;}
	int threadLabelWidth=log10(totalCalc)+1+3;	// [Num]+whitespace+........+| [Num]+whitespace+........+|
	int progressLabelWidth=8;
	int colWidth=threadLabelWidth+progressLabelWidth+2;
	tWidth=numCols*colWidth;
	//resizeTerminal(tWidth+widthBuffer,tHeight+heightBuffer);
	int df=0,val=log10(totalCalc)+1,index=0;
	string tNum="";
	for(int i=0;i<tHeight;i++)
		{for(int j=0;j<numCols;j++)
			{index=i+j*tHeight+1;
			if(index<=totalCalc)
				{tNum=cnvrtNumToStrng(index,0);
				df=val-tNum.length();
				if(j==numCols-1){Output+=makeWhiteSpace(df)+"["+tNum+"] ........|";}
				else{Output+=makeWhiteSpace(df)+"["+tNum+"] ........| ";}}}
		Output+="\n";}
	return Output;}

string initialize_ZPRED_gui_display(int totalCalc,bool& MULTI_DISPLAY)
	{int heightBuffer=5,widthBuffer=10,tWidth=0,tHeight=0;
	string Output="";
	if(totalCalc>MAX_DISPLAY_ROWS){tHeight=MAX_DISPLAY_ROWS;}
	else{tHeight=totalCalc;}
	// Determine Number of Columns to Display
	int numCols;
	if(totalCalc>=MAX_DISPLAY_ROWS*MAX_DISPLAY_COLS){numCols=MAX_DISPLAY_COLS;}//MULTI_DISPLAY=true;}	
	else{numCols=totalCalc/MAX_DISPLAY_ROWS+1;}
	int threadLabelWidth=log10(totalCalc)+1+3;	// [Num]+whitespace+........+| [Num]+whitespace+........+|
	int progressLabelWidth=8;
	int colWidth=threadLabelWidth+progressLabelWidth+2;
	tWidth=numCols*colWidth;
	//resizeTerminal(tWidth+widthBuffer,tHeight+heightBuffer);
	int df=0,val=log10(totalCalc)+1,index=0;
	string tNum="";
	for(int i=0;i<tHeight;i++)
		{for(int j=0;j<numCols;j++)
			{index=i+j*tHeight+1;
			if(index<=totalCalc)
				{tNum=cnvrtNumToStrng(index,0);
				df=val-tNum.length();
				if(j==numCols-1){Output+=makeHTMLWhiteSpace(2*df)+"<b>["+tNum+"] .&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;|</b>";}
				else{Output+=makeHTMLWhiteSpace(2*df)+"<b>["+tNum+"] .&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;|</b> ";}}}
		Output+="<br>";}
	return Output;}

void make_directory(string fPath)
	{int Err;
	char errMsg[256];
	int status=mkdir(fPath.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
	if(status==-1)
		{Err=errno;
		if(Err!=EEXIST)
			{cerr<<"Error in make_directory!\nError: "<<strerror_r(Err,errMsg,256)<<"\n"; exit(EXIT_FAILURE);}
		}
	}

void make_folder(string Fldr)
	{boost::filesystem::path dir(Fldr.c_str());
	boost::filesystem::create_directory(dir);}

string makeUpperCase(string X){string Output="";char letter,c;for(int i=0;i<X.length();i++){letter=X[i];Output+=toupper(letter);}return Output;}

string makeWhiteSpace(int Sz){string output="";for(int i=0;i<Sz;i++){output+=" ";}return output;}

string makeHTMLWhiteSpace(int Sz){string output="";for(int i=0;i<Sz;i++){output+="&nbsp;";}return output;}

int remove_directory(string fPath)
	{const char *path={fPath.c_str()};
	if(directory_exist(fPath))
		{return nftw(path, unlink_cb, 64, FTW_DEPTH | FTW_PHYS);}
	}

void remove_file(string fPath)
	{int result;
	if(file_exist(fPath))
		{result=remove(fPath.c_str());
		if(result!=0){perror(fPath.c_str());}
		}
	}

int unlink_cb(const char *fPath,const struct stat *sb, int typeFlag, struct FTW *ftwbuf)
	{int rv=remove(fPath);
	if(rv)
		{perror(fPath);}
	return rv;}

void resizeTerminal(int width,int height)
	{/* Terminal Dimension Values are in Characters*/
	string resizeTerminal="\e[8;"+cnvrtNumToStrng(height,0)+";"+cnvrtNumToStrng(width,0)+"t";printf("%s",resizeTerminal.c_str());}

zInput read_zpred_input_file(string inFile)
	{// , = internal delimiter| ; = external delimiter between solution conditinos
	zInput D;
	// Initialize Default Values
	//D.id="NoId";
	//D.name="molecule";
	D.numSolCond=0;
	D.numFiles=0;
	D.SC.apbsExecutable="";
	D.SC.multivalueExecutable="";
	D.SC.hydroproExecutable="";
	D.SC.msmsExecutable="";
	D.SC.msmsAtmTypeNumbers="";
	D.SC.msmsPdbToXyzr="";
	D.SC.msmsPdbToXyzrn="";
	D.SC.pdb2pqr="";			// main.py of PDB2PQR
	D.msms.PDB="";
	D.msms.pdbFile="";
	D.msms.pdb_msmsFldr="";
	D.msms.pdb_surfFldr="";
	D.msms.surfPointDensity="1.000000";
	D.msms.surfPointHiDensity="3.000000";
	D.hydroproCalcType="1";
	D.saveTemps=false; // ?3?
	D.pdie="4.0";
	D.dime_x="161"; D.dime_y="161"; D.dime_z="161";
	D.maxThreads=MAX_THREADS;
	D.RUN_APBS=true;
	D.RUN_PDB2PQR=true;
	D.RUN_HYDROPRO=false;
	D.RUN_HULLRAD=true;
	D.RUN_ZPRED=true;
	D.RUN_COLLIDE=true;
	D.outputTitle="Output/";
	D.numApbsWriteTypes=1;
	D.apbsWriteTypes=fill_string_array("pot;",1,";");
	D.forceField="PARSE";
	D.plot="";

	ifstream fIn;
	fIn.open(inFile.c_str());
	if(fIn.fail()){cerr<<"ERROR in read_zpred_input_file!\nInput file could not be opened.\n"<<inFile<<endl;exit(EXIT_FAILURE);}
	int Sz=10000,pos,oldPos,N;
	char Val[Sz];
	string bld,tmp,inputName,inputVal,delimiter="",*sArr;
	fIn.getline(Val,Sz);
	bool CHECKING;
	while(!fIn.eof())
		{tmp=Val;
		pos=tmp.find(".",0);
		if(pos!=string::npos && pos==0)
			{// Start of Line contains a period, now check if line defines an input parameter
			pos=tmp.find(":",pos);	// defines end of input parameter name
			if(pos!=string::npos)
				{// Extract Input Parameter Name
				inputName=tmp.substr(tmp.find(".",0)+1,pos-tmp.find(".",0)-1);
				inputName=makeUpperCase(inputName);
				// Extract Input Parameter Value(s)
				inputVal=tmp.substr(pos+1,tmp.length()-pos-1);
				// Check for & Remove Preceding WhiteSpace(s)
				CHECKING=true;
				pos=inputVal.find(" ",0);
				if(inputName.compare("FILES")==0||inputName.compare("APBSEXECUTABLE")==0||inputName.compare("MULTIVALUEEXECUTABLE")==0||inputName.compare("HYDROPROEXECUTABLE")==0||inputName.compare("MSMSEXECUTABLE")==0||inputName.compare("MSMSATMTYPENUMBERS")==0||inputName.compare("MSMSPDBTOXYZR")==0||inputName.compare("MSMSPDBTOXYZRN")==0||inputName.compare("PDB2PQR")==0)
					{// Remove preceeding white space if present
					if(pos==0){inputVal=inputVal.substr(1,inputVal.length()-1);}
					}
				else
					{
					while(CHECKING)
						{if(pos!=string::npos){inputVal=inputVal.substr(1,inputVal.length()-1);pos=inputVal.find(" ",0);}
						else{CHECKING=false;break;}}
					}
				if(inputVal.length()!=0)
					{// Identify Parameter and Load Values
					if(inputName.compare("ID")==0)
						{// Id for Computation (for organization purposes)
						D.numSolCond=count_delimiter(inputVal,";");
						D.id=fill_string_array(inputVal,D.numSolCond,";");}
					else if(inputName.compare("FILES")==0)
						{// Acquire PDB File Names
						D.numFiles=count_delimiter(inputVal,";");
						D.files=fill_string_array(inputVal,D.numFiles,";");}
					else if(inputName.compare("PH")==0)
						{// Array of pH Values
						D.pH=fill_string_array(inputVal,D.numSolCond,";");}
					else if(inputName.compare("SOLUTE")==0)
						{// Array of Solute Names ,;,;,;
						N=count_delimiter(inputVal,";");
						if(N!=D.numSolCond){cerr<<"ERROR in read_zpred_input_file!\nNumber of \'solute\' names ("<<N<<") does not match number of solution conditions.\nInput:"<<inputVal<<"\n";exit(EXIT_FAILURE);}
						sArr=fill_string_array(inputVal,N,";");
						D.solute=new string*[N];
						D.numSolute=new int[N];
						for(int i=0;i<N;i++)
							{tmp=sArr[i];
							D.numSolute[i]=count_delimiter(tmp,",");
							D.solute[i]=fill_string_array(tmp,D.numSolute[i],",");}
						delete [] sArr;}
					else if(inputName.compare("SOLUTECONC")==0)
						{// Array of Solute Concentrations soluteConc[numSolute][numSoluteConc]  (Ion#1Conc,Ion#1Conc,;Ion#2Conc,;)
						N=count_delimiter(inputVal,";");
						if(N!=D.numSolCond){cerr<<"ERROR in read_zpred_input_file!\nNumber of \'soluteConc\' ("<<N<<") does not match number of solution conditions.\nInput:"<<inputVal<<"\n";exit(EXIT_FAILURE);}
						sArr=fill_string_array(inputVal,N,";");
						D.soluteConc=new double*[N];
						for(int i=0;i<N;i++)
							{tmp=sArr[i];
							D.soluteConc[i]=fill_double_array(tmp,D.numSolute[i],",");}
						delete [] sArr;}
					else if(inputName.compare("SOLUTECONCTYPE")==0)
						{// Specifies type of concentration of soluteConc (always one value)
						N=count_delimiter(inputVal,";");
						if(N!=D.numSolCond){cerr<<"ERROR in read_zpred_input_file!\nNumber of \'soluteConcType\' ("<<N<<") does not match number of solution conditions.\nInput:"<<inputVal<<"\n";exit(EXIT_FAILURE);}
						/*sArr=fill_string_array(inputVal,N,";");
						D.soluteConcType=new string*[N];
						for(int i=0;i<N;i++)
							{tmp=sArr[i];
							D.soluteConcType[i]=fill_string_array(tmp,D.numSolute[i],",");}
						delete [] sArr;
						*/
						D.soluteConcType=fill_string_array(inputVal,N,";");
						}
					else if(inputName.compare("SOLVENT")==0)
						{// Array of Solvent Names ,;,;,;
						N=count_delimiter(inputVal,";");
						if(N!=D.numSolCond){cerr<<"ERROR in read_zpred_input_file!\nNumber of \'solvent\' names ("<<N<<") does not match number of solution conditions.\nInput:"<<inputVal<<"\n";exit(EXIT_FAILURE);}
						sArr=fill_string_array(inputVal,N,";");
						D.solvent=new string*[N];
						D.numSolvent=new int[N];
						for(int i=0;i<N;i++)
							{tmp=sArr[i];
							D.numSolvent[i]=count_delimiter(tmp,",");
							D.solvent[i]=fill_string_array(tmp,D.numSolvent[i],",");}
						delete [] sArr;}
					else if(inputName.compare("SOLVENTCONC")==0)
						{// Array of Solvent Concentrations
						N=count_delimiter(inputVal,";");
						if(N!=D.numSolCond){cerr<<"ERROR in read_zpred_input_file!\nNumber of \'solventConc\' ("<<N<<") does not match number of solution conditions.\nInput:"<<inputVal<<"\n";exit(EXIT_FAILURE);}
						sArr=fill_string_array(inputVal,N,";");
						D.solventConc=new double*[N];
						for(int i=0;i<N;i++)
							{tmp=sArr[i];
							D.solventConc[i]=fill_double_array(tmp,D.numSolvent[i],",");}
						delete [] sArr;}
					else if(inputName.compare("SOLVENTCONCTYPE")==0)
						{// Specifies type of concentration of solventConc
						N=count_delimiter(inputVal,";");
						if(N!=D.numSolCond){cerr<<"ERROR in read_zpred_input_file!\nNumber of \'solventConcType\' ("<<N<<") does not match number of solution conditions.\nInput:"<<inputVal<<"\n";exit(EXIT_FAILURE);}
						/*sArr=fill_string_array(inputVal,N,";");
						D.solventConcType=new string*[N];
						for(int i=0;i<N;i++)
							{tmp=sArr[i];
							D.solventConcType[i]=fill_string_array(tmp,D.numSolvent[i],",");}
						delete [] sArr;*/
						D.solventConcType=fill_string_array(inputVal,N,";");
						}
					else if(inputName.compare("TEMPERATURE")==0)
						{// Array of Temperatures (Kelvin)
						N=count_delimiter(inputVal,";");
						if(N!=D.numSolCond){cerr<<"ERROR in read_zpred_input_file!\nNumber of \'temperature\' ("<<N<<") does not match number of solution conditions.\nInput:"<<inputVal<<"\n";exit(EXIT_FAILURE);}
						D.temperature=fill_double_array(inputVal,N,";");}
					else if(inputName.compare("PROTEINCONC")==0)
						{// Protein Concentration (mg/mL)
						N=count_delimiter(inputVal,";");
						if(N!=D.numSolCond){cerr<<"ERROR in read_zpred_input_file!\nNumber of \'proteinConc\' ("<<N<<") does not match number of solution conditions.\nInput:"<<inputVal<<"\n";exit(EXIT_FAILURE);}
						D.proteinConc=fill_double_array(inputVal,N,";");}
					else if(inputName.compare("APBSWRITETYPES")==0)
						{// APBS Output Write Types
						N=count_delimiter(inputVal,";");
						D.numApbsWriteTypes=N;
						delete [] D.apbsWriteTypes;
						D.apbsWriteTypes=fill_string_array(inputVal,N,";");}
					else if(inputName.compare("APBSEXECUTABLE")==0)
						{// APBS Executable File Path
						D.SC.apbsExecutable=inputVal;}
					else if(inputName.compare("MULTIVALUEEXECUTABLE")==0)
						{// APBS Multivalue Executable File Path
						D.SC.multivalueExecutable=inputVal;}
					else if(inputName.compare("HYDROPROEXECUTABLE")==0)
						{// HYDROPRO Executable File Path
						D.SC.hydroproExecutable=inputVal;}
					else if(inputName.compare("MSMSEXECUTABLE")==0)
						{// MSMS Executable File Path
						D.SC.msmsExecutable=inputVal;}
					else if(inputName.compare("MSMSATMTYPENUMBERS")==0)
						{// MSMS Executable File Path
						D.SC.msmsAtmTypeNumbers=inputVal;}
					else if(inputName.compare("MSMSPDBTOXYZR")==0)
						{// MSMS Executable File Path
						D.SC.msmsPdbToXyzr=inputVal;}
					else if(inputName.compare("MSMSPDBTOXYZRN")==0)
						{// MSMS Executable File Path
						D.SC.msmsPdbToXyzrn=inputVal;}
					else if(inputName.compare("PDB2PQR")==0)
						{// File Path to PDB2PQR python script
						D.SC.pdb2pqr=inputVal;}
					else if(inputName.compare("MSMSSURFPOINTDENSITY")==0)
						{// MSMS Parameter: Low Surface Density
						D.msms.surfPointDensity=inputVal;}
					else if(inputName.compare("MSMSSURFPOINTHIDENSITY")==0)
						{// MSMS Parameter: High Surface Density
						D.msms.surfPointHiDensity=inputVal;}
					else if(inputName.compare("HYDROPROCALCTYPE")==0)
						{// HYDROPRO Parameter: Calculation Type selects model
						D.hydroproCalcType=inputVal;}
					else if(inputName.compare("FORCEFIELD")==0)
						{// PDB2PQR Force Field
						D.forceField=inputVal;}
					else if(inputName.compare("SAVETEMPS")==0)
						{// Flag for Saving Files Used in the Computation
						inputVal=makeUpperCase(inputVal);
						if(inputVal.compare("FALSE")==0){D.saveTemps=false;}
						else{D.saveTemps=true;}}
					else if(inputName.compare("PDIE")==0)
						{// APBS Parameter: Protein Relative Dielectric
						D.pdie=inputVal;}
					else if(inputName.compare("DIME_X")==0)
						{// APBS Parameter: Number of X-Axis Grid Points
						D.dime_x=inputVal;}
					else if(inputName.compare("DIME_Y")==0)
						{// APBS Parameter: Number of Y-Axis Grid Points
						D.dime_y=inputVal;}
					else if(inputName.compare("DIME_Z")==0)
						{// APBS Parameter: Number of Z-Axis Grid Points
						D.dime_z=inputVal;}
					else if(inputName.compare("RUN_APBS")==0)
						{// ZPRED Control Paramater
						inputVal=makeUpperCase(inputVal);
						if(inputVal.compare("FALSE")==0){D.RUN_APBS=false;}
						else{D.RUN_APBS=true;}}
					else if(inputName.compare("RUN_PDB2PQR")==0)
						{// ZPRED Control Paramater
						inputVal=makeUpperCase(inputVal);
						if(inputVal.compare("FALSE")==0){D.RUN_PDB2PQR=false;}
						else{D.RUN_PDB2PQR=true;}}
					else if(inputName.compare("RUN_HYDROPRO")==0)
						{// ZPRED Control Paramater
						inputVal=makeUpperCase(inputVal);
						if(inputVal.compare("FALSE")==0){D.RUN_HYDROPRO=false;}
						else{D.RUN_HYDROPRO=true;}}
					else if(inputName.compare("RUN_HULLRAD")==0)
						{// ZPRED Control Paramater
						inputVal=makeUpperCase(inputVal);
						if(inputVal.compare("FALSE")==0){D.RUN_HULLRAD=false;}
						else{D.RUN_HULLRAD=true;}}
					else if(inputName.compare("RUN_ZPRED")==0)
						{// ZPRED Control Paramater
						inputVal=makeUpperCase(inputVal);
						if(inputVal.compare("FALSE")==0){D.RUN_ZPRED=false;}
						else{D.RUN_ZPRED=true;}}
					else if(inputName.compare("RUN_COLLIDE")==0)
						{// ZPRED Control Paramater
						inputVal=makeUpperCase(inputVal);
						if(inputVal.compare("FALSE")==0){D.RUN_COLLIDE=false;}
						else{D.RUN_COLLIDE=true;}}
					else if(inputName.compare("OUTPUTTITLE")==0)
						{// ZPRED Output Folder Title
						D.outputTitle=inputVal;}
					else if(inputName.compare("MAXTHREADS")==0)
						{// ZPRED Control Paramater
						D.maxThreads=atoi(inputVal.c_str());}
					else if(inputName.compare("DENSITY")==0)
						{// User Specified Solution Property
						N=count_delimiter(inputVal,";");
						if(N!=D.numSolCond){cerr<<"ERROR in read_zpred_input_file!\nNumber of \'density\' ("<<N<<") does not match number of solution conditions.\n";exit(EXIT_FAILURE);}
						D.density=fill_double_array(inputVal,N,";");}
					else if(inputName.compare("DIELECTRIC")==0)
						{// User Specified Solution Property
						N=count_delimiter(inputVal,";");
						if(N!=D.numSolCond){cerr<<"ERROR in read_zpred_input_file!\nNumber of \'density\' ("<<N<<") does not match number of solution conditions.\n";exit(EXIT_FAILURE);}
						D.dielectric=fill_double_array(inputVal,N,";");}
					else if(inputName.compare("INFLATIONDISTANCE")==0)
						{// User Specified Solution Property
						N=count_delimiter(inputVal,";");
						if(N!=D.numSolCond){cerr<<"ERROR in read_zpred_input_file!\nNumber of \'density\' ("<<N<<") does not match number of solution conditions.\n";exit(EXIT_FAILURE);}
						D.inflationDistance=fill_double_array(inputVal,N,";");}
					else if(inputName.compare("VISCOSITY")==0)
						{// User Specified Solution Property
						N=count_delimiter(inputVal,";");
						if(N!=D.numSolCond){cerr<<"ERROR in read_zpred_input_file!\nNumber of \'density\' ("<<N<<") does not match number of solution conditions.\n";exit(EXIT_FAILURE);}
						D.viscosity=fill_double_array(inputVal,N,";");}
					else if(inputName.compare("GENERATEPROFILE")==0)
						{N=count_delimiter(inputVal,";");
						if(N!=D.numSolCond){cerr<<"ERROR in read_zpred_input_file!\nNumber of \'generateProfile\' ("<<N<<") does not match number of solution conditions.\n";exit(EXIT_FAILURE);}
						D.GENERATE_POT_PROFILE=new bool[N];
						sArr=fill_string_array(inputVal,N,";");
						for(int i=0;i<N;i++)
							{tmp=sArr[i];
							if(tmp.compare("true")==0){D.GENERATE_POT_PROFILE[i]=true;}
							else if(tmp.compare("false")==0){D.GENERATE_POT_PROFILE[i]=false;}
							else{cerr<<"ERROR in read_zpred_input_file!\nNumber of \'generateProfile\' has incorrect values.\n";exit(EXIT_FAILURE);}
							}
						delete [] sArr;}
					else if(inputName.compare("PLOT")==0)
						{// Plotting Feature
						D.plot=inputVal;}
					else{cerr<<"ERROR in read_zpred_input_file!\nUnrecognized Parameter ("<<inputName<<")."<<endl;exit(EXIT_FAILURE);}
					}
				}
			}
		fIn.getline(Val,Sz);}
	fIn.close();
	return D;}

string repStr(string Input,int Num)
	{return repeatString(Input,Num);}

string repeatString(string Input,int Num)
	{string Output="";
	for(int i=0;i<Num;i++){Output+=Input;}
	return Output;}

void rewind_display(string display)
	{// Determine Number of Rows
	int nRows=0,pos=0;
	bool SEARCHING=true;
	while(SEARCHING)
		{pos=display.find("\n",pos);
		if(pos==string::npos){break;SEARCHING=false;}
		else{nRows++;pos++;}}
	// Rewind to Beginning of Line of 1st Thread
	string upLn="\033["+cnvrtNumToStrng(nRows,0)+"A"+"\r\033[s";
	cerr<<upLn;}

void* runZPRED(void* Ptr)
	{bool VERBOSE=false;
	string delimiter=";";
	// Convert Input Pointer to String
	string *strPtr=reinterpret_cast<string*>(Ptr);
	string cmdString=*strPtr;
	string Fldr="",encSalt="",encPass="",calcNum="",display="";
	zExeInput D=cnvrtStr2Input(cmdString,";",Fldr,calcNum,display);

	// Determine Number of Rows
	int nRows=0,pos=0;
	bool SEARCHING=true;
	while(SEARCHING)
		{pos=display.find("\n",pos);
		if(pos==string::npos){break;SEARCHING=false;}
		else{nRows++;pos++;}}
	// Find Display String Line with Current Computation
	string dispLabel="["+calcNum+"]";
	pos=display.find(dispLabel,0);
	// Determine Beginning of Line
	int tPos;
	if(pos!=string::npos){tPos=pos-1;}
	else{tPos=0;}
	SEARCHING=true;
	string tmp="";
	while(SEARCHING)
		{if(tPos>=0){tmp=display.substr(tPos,1);}
		else if(tPos<0){tPos=0;break;}
		if(tmp.compare("\n")==0)
			{// Line Break Found
			tPos++;break;SEARCHING=false;}
		else{tPos--;}}
	// Acquire Entire Display Line
	string dispLn=display.substr(tPos,display.find("\n",tPos)-tPos);
	// Check if On Second Column or Later
	int currColOffset=0,cN=atoi(calcNum.c_str());
	// When on first column, cN%(MAX_DISPLAY_ROWS*MAX_DISPLAY_COLS)=1-25
	int test=cN%(MAX_DISPLAY_ROWS*MAX_DISPLAY_COLS);
	
	string bld;
	if(test==0 || test>nRows)
		{// Determine Current Column Offset from First Column
		currColOffset=((cN-1)%(MAX_DISPLAY_ROWS*MAX_DISPLAY_COLS))/nRows;
		// Update Display Line String
		pos=0;bld="";
		for(int i=0;i<currColOffset;i++)
			{bld+=dispLn.substr(pos,dispLn.find("]",pos)-pos);
			pos=dispLn.find("]",pos);
			if(pos==string::npos){cerr<<"ERROR in runZPRED ["<<calcNum<<"]!\nDisplay Line could not be updated.\n"<<dispLn<<endl; exit(EXIT_FAILURE);}
			bld+=dispLn.substr(pos,2)+"COMPLETE";
			pos=dispLn.find("|",pos);}
		if(pos!=string::npos){bld+=dispLn.substr(pos,dispLn.length()-pos);}
		else{bld+=dispLn.substr(0,dispLn.length());}
		
		dispLn=bld;}
	
	// Create list of files to delete after computation (when saveTemps=false)
	string files2Delete="";

	// Define Base Folders
	// ZPRED Input
	//string pdbFldr=Fldr+"Input/";
	// Define GUI Display File
	string guiDisplayFile=Fldr+".gui_display.txt";
	// Formatted Data Folders Folder (Also Holds ZPRED Configuration File)
	string fFldr=Fldr+"Files/";
	// COLLIDE Input File Containing List of ZPRED Output Files to Run
	string collideFile=fFldr+"zpredOutFilesList.txt";
	ofstream colIn;
	if(calcNum.compare("1")==0)
		{// At First Process Generate Empty File for Appending
		colIn.open(collideFile.c_str(),ifstream::trunc);
		if(colIn.fail()){cerr<<"ERROR in runZPRED!!!\nZPRED Output File List could not be opened.\n"<<collideFile<<endl; exit(EXIT_FAILURE);}
		colIn.close();}
	// Single ZPRED Computation Output Folder
	string zpredFldr=fFldr+"Z/";
	// Program Stdout/Stderr Log Folder
	string logFldr=fFldr+"LOG/";
	// APBS Output Folder (in DX Matrix Format)
	string dxFldr=fFldr+"DX/";
	// HYDROPRO Output Folder (HP stands for hydrodynamic properties)
	string hydroFldr=fFldr+"HP/";
	// APBS Input Command File Folder
	string inFldr=fFldr+"IN/";	
	// PROPKA Output Folder
	string propkaFldr=fFldr+"PROPKA/";	
	// PQR Formatted PDB Files Folder
	string pqrFldr=fFldr+"PQR/";
	// MSMS Surface Output Folder
	string surfFldr=fFldr+"S/";
	string surfFldrContents[6]={"CGO/","colorCGO/","CSV/","fromMSMS/","SAS/","fromMV/"};
	// Handle Special Operations Folders
	string msmsFldr=fFldr+"msms/";
	if(!directory_exist(msmsFldr)){make_folder(msmsFldr);}
	// APBS Inputgen Python Script
	string apbsInputgen=fFldr+"apbsInputgen.py"; //writeAPBSInputGen(apbsInputgen);
	// APBS psize.py supports apbsInputgen.py
	//string psize=fFldr+"psize.py"; writePsize(psize);
	// HullRad (Grisham Edit)
	string theHullRad=fFldr+"editedHullRad.py"; //writeHullRad(theHullRad);

	int PID,rChk=0,status=0;
	//chmod_x(D.SC.apbsExecutable);
	// APBS Multivalue (Provide executable permission)
	//chmod_x(D.SC.multivalueExecutable);
	// PDB2PQR
	string PDB2PQR=D.SC.pdb2pqr;
	// Acquire Specific File Name
	tmp=D.file;
	// Check if file is .pqr
	pos=tmp.rfind(".pqr");
	string PDB,pdbFile,pdb_pqrFldr,pqrFile,pHVal=D.pH,tmp2,saveAsPdbFile;
	if(pos!=string::npos)
		{pdbFile=tmp.substr(0,pos)+".pdb";
		pos=tmp.rfind("/",tmp.length()-1);
		PDB=tmp.substr(pos+1,tmp.find(".pqr",pos)-pos-1);
		pdb_pqrFldr=pqrFldr+PDB+"/";
		if(!directory_exist(pdb_pqrFldr)){make_folder(pdb_pqrFldr);}
		pqrFile=pdb_pqrFldr+"pH"+pHVal+"_"+calcNum+".pqr";
		if(copyFile(tmp,pqrFile)){}
		tmp2=Fldr.substr(0,Fldr.length()-2);
		saveAsPdbFile=get_containing_folder(tmp2)+"saveAsPdb.py";
		convert_pqr_to_pdb(saveAsPdbFile,pqrFile,pdbFile);
		}
	else
		{pos=tmp.rfind("/",tmp.length()-1);
		PDB=tmp.substr(pos+1,tmp.find(".pdb",pos)-pos-1);
		pdb_pqrFldr=pqrFldr+PDB+"/";
		if(!directory_exist(pdb_pqrFldr)){make_folder(pdb_pqrFldr);}
		// Define PDB File under Assessment
		pdbFile=D.file;
		}

	// Need to add shape assessment
	string shape="sphere",value;	
	double T=D.temperature;
	// Define String of Solvent Names
	string solvStr="";
	for(int i=0;i<D.numSolvent;i++){solvStr+=D.solvent[i]+delimiter;}
	// Define String of Solvent Concentrations
	string solvConcStr="";
	for(int i=0;i<D.numSolvent;i++){solvConcStr+=cnvrtNumToStrng(D.solventConc[i],BUFFER_SIZE)+delimiter;}
	// Define String of Solute Names
	string soluStr="";
	for(int i=0;i<D.numSolute;i++){soluStr+=D.solute[i]+delimiter;}
	// Define String of Solute Concentrations
	string soluConcStr="";
	for(int i=0;i<D.numSolute;i++){soluConcStr+=cnvrtNumToStrng(D.soluteConc[i],BUFFER_SIZE)+delimiter;}
	// Define Solution Name	
	string buffer=nameBuffer(soluStr,soluConcStr);
	// Define Folder Extension
	string fldrExt;
	if(buffer.length()<15){fldrExt="_"+buffer+"/";}
	else{fldrExt="_"+calcNum+"/";}
	//fldrExt="_"+buffer+"/";
	// Create PDB Specific Fldrs (use PDB file name not PDB id)
	// APBS Output Folder (dxFldr) Contents
	string pdb_dxFldr=dxFldr+PDB+fldrExt;
	if(!directory_exist(pdb_dxFldr)){make_folder(pdb_dxFldr);}
	string wtFldr="";
	for(int d=0;d<D.numApbsWriteTypes;d++)
		{wtFldr=pdb_dxFldr+D.apbsWriteTypes[d]+"/";
		if(!directory_exist(wtFldr)){make_folder(wtFldr);}
		}
	// APBS Input File Folder
	string pdb_inFldr=inFldr+PDB+fldrExt;
	if(!directory_exist(pdb_inFldr)){make_folder(pdb_inFldr);}
	// PROPKA Format Folder
	string pdb_propkaFldr=propkaFldr+PDB+"/";
	if(!directory_exist(pdb_propkaFldr)){make_folder(pdb_propkaFldr);}
	// ZPRED Output Folder for Single Calculation
	string pdb_zpredFldr=zpredFldr+PDB+fldrExt; 
	if(!directory_exist(pdb_zpredFldr)){make_folder(pdb_zpredFldr);}

	// Job Already Finished Check
	string zpredOutFile=pdb_zpredFldr+"pH"+pHVal+"T"+formatNumberString(cnvrtNumToStrng(T,BUFFER_SIZE))+"PC"+cnvrtNumToStrng(D.proteinConc,ER_SIG_FIGS)+"_"+calcNum+".txt";
	if(file_exist(zpredOutFile))
		{// Add ZPRED Output File to List
		colIn.open(collideFile.c_str(),ifstream::app);
		if(colIn.fail()){cerr<<"ERROR in runZPRED!!!\nZPRED Output File List could not be updated.\n"<<collideFile<<endl; exit(EXIT_FAILURE);}
		colIn<<zpredOutFile+"\n";
		colIn.close();
		return 0;}

	// PDB Specific Log Folder
	string pdb_logFldr=logFldr+PDB+"/"; make_folder(pdb_logFldr);
	string logFile=pdb_logFldr+calcNum+".log"; //files2Delete+=logFile+delimiter;
	// MSMS Computation Folder
	string pdb_msmsFldr=msmsFldr+PDB+"_"+calcNum+"_pH"+pHVal+"_T"+formatNumberString(cnvrtNumToStrng(T,BUFFER_SIZE))+fldrExt; make_folder(pdb_msmsFldr);
	string msmsBinFile=pdb_msmsFldr+getFileName(D.SC.msmsExecutable);
	if(copyFile(D.SC.msmsExecutable,msmsBinFile)){chmod_x(msmsBinFile);}
	files2Delete+=msmsBinFile+delimiter;
	string pdbToXyzrFile=pdb_msmsFldr+getFileName(D.SC.msmsPdbToXyzr);
	if(copyFile(D.SC.msmsPdbToXyzr,pdbToXyzrFile)){chmod_x(pdbToXyzrFile);}
	files2Delete+=pdbToXyzrFile+delimiter;
	string pdbToXyzrnFile=pdb_msmsFldr+getFileName(D.SC.msmsPdbToXyzrn);
	if(copyFile(D.SC.msmsPdbToXyzrn,pdbToXyzrnFile)){chmod_x(pdbToXyzrnFile);}
	files2Delete+=pdbToXyzrnFile+delimiter;
	string atmTypeNumbersFile=pdb_msmsFldr+getFileName(D.SC.msmsAtmTypeNumbers);
	if(copyFile(D.SC.msmsAtmTypeNumbers,atmTypeNumbersFile)){chmod_x(atmTypeNumbersFile);}
	files2Delete+=atmTypeNumbersFile+delimiter;
	// MSMS xyzrn file generated for input structure file
	string msmsXyzrnFile=pdb_msmsFldr+PDB+"_"+calcNum+".xyzrn"; files2Delete+=msmsXyzrnFile+delimiter;
	// HYDROPRO Output Folder (must be different for temperature,pH,and solute
	string pdb_hydroFldr=hydroFldr+PDB+"_"+calcNum+"_pH"+pHVal+"_T"+formatNumberString(cnvrtNumToStrng(T,BUFFER_SIZE))+fldrExt; make_folder(pdb_hydroFldr);
	// Copied PDB File
	string newPdbFile=pdb_hydroFldr+PDB+"_"+calcNum+".pdb"; files2Delete+=newPdbFile+delimiter;
	// Prepare to delete additional HYDROPRO files
	// HYDROPRO Atomic BEAD Model File
	string hydroBeaFile=pdb_hydroFldr+PDB+"_"+calcNum+"-pri.bea"; files2Delete+=hydroBeaFile+delimiter;
	// HYDROPRO VRML File
	string hydroVrmlFile=pdb_hydroFldr+PDB+"_"+calcNum+"-pri.vrml"; files2Delete+=hydroVrmlFile+delimiter;
	// HYDROPRO Sol File
	string hydroSolFile=pdb_hydroFldr+PDB+"_"+calcNum+"-sol.txt"; files2Delete+=hydroSolFile+delimiter;
	// HYDROPRO Fit File
	string hydroFitFile=pdb_hydroFldr+"hydropro-fit.txt"; files2Delete+=hydroFitFile+delimiter;
	// HYDROPRO Summary File
	string hydroSumFile=pdb_hydroFldr+"hydropro-sum.txt"; files2Delete+=hydroSumFile+delimiter;
	string transferString="";
	// Copy PDB File into HYDROPRO Computation Folder
	ifstream fIn;
	bool ATTEMPTING=true;
	int numTries=MAX_FILE_OPEN_ATTEMPTS,Counter=0;
	while(ATTEMPTING)
		{fIn.open(pdbFile.c_str());
		if(!fIn.fail()){break;ATTEMPTING=false;}
		else{Counter++;}
		if(Counter>=numTries)
			{cerr<<"ERROR in runZPRED("<<calcNum<<")!!!\nOriginal PDB File could not be copied.\n"<<pdbFile<<endl;
			exit(EXIT_FAILURE);}}
	ofstream fOut;
	fOut.open(newPdbFile.c_str());
	if(fOut.fail()){cerr<<"ERROR in runZPRED!!!\nPDB File could not be copied.\n"<<newPdbFile<<endl;exit(EXIT_FAILURE);}
	 // Copy Input File
    int Sz=15000;                      // overshoot actual line size, C++ automatically fixes c-string lengths
    char Val[Sz];
    fIn.getline(Val,Sz);
    while(!fIn.eof()){tmp=Val;transferString+=tmp+"\n";fIn.getline(Val,Sz);}
	fIn.close();fOut<<transferString;fOut.close();
	// Generate Temporary Executable
	string tempXFile=pdb_hydroFldr+"hydropro10-lnx"+calcNum+".exe";
	tmp=D.SC.hydroproExecutable;
	if(tmp.length()!=0)
		{if(copyFile(D.SC.hydroproExecutable,tempXFile))
			{chmod_x(tempXFile);}
		}
	// MSMS Surface (Output) Folder
	string pdb_surfFldr=surfFldr+PDB+fldrExt;make_folder(pdb_surfFldr);
	string sFldr="";
	for(int d=0;d<6;d++){sFldr=pdb_surfFldr+surfFldrContents[d];make_folder(sFldr);}	
	// Define Files (Need to Ensure File Path is less than MAX_FILE_PATH)
	// van der Waals Surface File (from MSMS)
	string vdwFile=pdb_surfFldr+"fromMSMS/SR0.00"+"_"+calcNum+".vert"; files2Delete+=vdwFile+delimiter;
	// MSMS van der Waals (VDW) Surface Output File
	string msmsOutputFile=pdb_surfFldr+"fromMSMS/msmsOutput_SR0.00"+"_"+calcNum+".txt"; files2Delete+=msmsOutputFile+delimiter;
	// MSMS VDW Area Output File
	string msmsAreaFile=pdb_surfFldr+"fromMSMS/SR0.00"+"_"+calcNum+".area"; files2Delete+=msmsAreaFile+delimiter;
	// MSMS VDW Face Output File
	string msmsFaceFile=pdb_surfFldr+"fromMSMS/SR0.00"+"_"+calcNum+".face"; files2Delete+=msmsFaceFile+delimiter;
	// HYDROPRO Input Command File (Must be named hydropro.dat)
	string hydroproInputFile=pdb_hydroFldr+"hydropro.dat"; files2Delete+=hydroproInputFile+delimiter;
	// HYDROPRO Output Computation File
	string hydroproOutputFile=pdb_hydroFldr+PDB+"_"+calcNum+"-res.txt"; files2Delete+=hydroproOutputFile+delimiter;
	// PQR File
	pqrFile=pdb_pqrFldr+"pH"+pHVal+"_"+calcNum+".pqr"; files2Delete+=pqrFile+delimiter;
	// PROPKA Output Computation File
	string propkaOutFile=pdb_propkaFldr+"pH"+pHVal+"_"+calcNum+".propka"; files2Delete+=propkaOutFile+delimiter;
	// APBS Input Command File
	string apbsInputFile=pdb_inFldr+"pH"+pHVal+"_"+calcNum+"T"+formatNumberString(cnvrtNumToStrng(T,BUFFER_SIZE))+"SR";
	// ZPRED Output File for Current Calculation (Not Total Calculation Output File!)
	zpredOutFile=pdb_zpredFldr+"pH"+pHVal+"T"+formatNumberString(cnvrtNumToStrng(T,BUFFER_SIZE));
	
	// Compute Solution Properties
	Solution S(solvStr,solvConcStr,D.solventConcType,soluStr,soluConcStr,T,delimiter);
	msmsInput msmsIn;
	vertData V;
	pdbData Data;
	strVec apbsGrdDim; apbsGrdDim.x=D.dime_x; apbsGrdDim.y=D.dime_y; apbsGrdDim.z=D.dime_z;
	ion_Data In;
	// Values of Mapping the Electrophoretic Space to define Electrokinetic Model
	double kaVals[132]={0.0100467,0.0102119,0.0105752,0.0110796,0.0116896,0.0123908,0.013134,0.0139218,0.0147913,0.0157151,0.0166577,0.0176569,0.0188034,0.0199777,0.0212255,0.0225511,0.0240154,0.0255153,0.0271721,0.0289364,0.0307437,0.0328163,0.0351104,0.0375649,0.040191,0.043101,0.0463294,0.0501491,0.0545372,0.0590336,0.0634555,0.0682086,0.0734888,0.0793625,0.0855062,0.0923404,0.10042,0.109207,0.118487,0.128854,0.140129,0.152391,0.166888,0.18319,0.201554,0.222795,0.247424,0.27606,0.310171,0.348496,0.393386,0.448216,0.515471,0.595586,0.691367,0.8063,0.946936,1.1147,1.31524,1.55186,1.83532,2.17056,2.56703,3.02886,3.56545,4.18735,4.91771,5.74861,6.68866,7.74625,8.90856,10.1977,11.6733,13.3003,15.0485,16.908,18.9529,21.1465,23.4842,26.0803,28.7617,31.5712,34.6553,37.952,41.2729,44.8843,48.8118,53.0829,57.3256,61.7633,66.5446,71.696,77.0664,82.4539,88.0128,93.9464,100.28,107.041,113.991,120.829,128.375,136.393,144.574,152.889,161.682,170.583,179.556,189,198.941,209.405,219.907,230.397,241.951,254.085,266.206,277.608,290.175,302.604,315.565,329.082,343.178,356.213,369.744,383.789,398.367,413.499,428.207,442.406,458.142,476.654,490.169,490.169};
	double EmVals[132]={5.1634,5.07903,4.99584,4.91384,4.83304,4.75223,4.67142,4.59061,4.51099,4.43137,4.35175,4.27094,4.19133,4.11171,4.03209,3.95247,3.87285,3.79441,3.7148,3.63518,3.55556,3.47594,3.39869,3.32026,3.24302,3.16578,3.08853,3.01248,2.9388,2.86275,2.7855,2.70945,2.63339,2.55734,2.4801,2.40404,2.32917,2.25431,2.17944,2.10458,2.0309,1.95722,1.88473,1.81462,1.7445,1.67558,1.60665,1.5413,1.4795,1.41889,1.36067,1.30838,1.26203,1.21925,1.18004,1.14795,1.12181,1.1016,1.08734,1.07665,1.07071,1.0719,1.07903,1.08972,1.10755,1.13012,1.15627,1.18717,1.22282,1.26322,1.30957,1.35829,1.40939,1.46405,1.52228,1.58408,1.64706,1.71242,1.78015,1.84789,1.91682,1.98812,2.05942,2.13191,2.20559,2.28045,2.35532,2.43018,2.50743,2.58348,2.66072,2.73797,2.81402,2.89245,2.97089,3.04932,3.12775,3.20737,3.28699,3.36661,3.44623,3.52585,3.60665,3.68746,3.76827,3.84908,3.92989,4.0107,4.09269,4.1735,4.25431,4.3363,4.4183,4.5003,4.58229,4.66548,4.74747,4.83066,4.91384,4.99703,5.08021,5.1634,5.24658,5.32977,5.41295,5.49614,5.57932,5.6637,5.74688,5.83125,5.91444,6};
	//double shapeFactor=D.shapeFactor[0];
	double molecularWeight=calcProteinMolecularWeight(pdbFile);
	// Solvent Viscosity (Solute(s)+Solvent(s)) (Pa s)
	double solventViscosity;
	if(D.viscosity==-1){solventViscosity=S.viscosity;}
	else{solventViscosity=D.viscosity;}
	// Solvent Density (kg/L)
	double solventDensity;
	if(D.density==-1){solventDensity=S.density;}
	else{solventDensity=D.density;}
	updateDisplayString("M",calcNum,1,nRows,dispLn);
	updateGuiDisplayString("M&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;",calcNum,1,nRows,dispLn,guiDisplayFile);
	// Initialize MSMS Input Parameters
	define_msmsInput(PDB,newPdbFile,msmsFldr,pdb_msmsFldr,pdb_surfFldr,D.msms.surfPointDensity,D.msms.surfPointHiDensity,"0.00",msmsBinFile,msmsIn);
	// Generate van der Waals Surface to Calculate Protein Radius
	runMSMS(msmsIn,calcNum);
	updateDisplayString(".",calcNum,1,nRows,dispLn);
	updateGuiDisplayString(".&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;",calcNum,1,nRows,dispLn,guiDisplayFile);
	// Calculate Protein Radius (A)
	// Read MSMS Vertices File
	readMSMSVertFile(vdwFile,V);
	// Read PDB File					
	int numLines=countNumLines(newPdbFile);
	definePDBData(Data,numLines);
	readPDBFile(newPdbFile,Data);
	// Calculate Center of Mass from PDB File
	Vec centerCoord=calc_center_of_mass(Data,numLines);
	// Calculate Average Protein Radius Using MSMS produced surface
	double proteinRadius=calc_protein_radius(V,V.numVertices,centerCoord);
	// Calculate Standard Deviation
	//radStd=calc_protein_radius_std(V,V.numVertices,centerCoord,Radius);
	// Calculate Specific Volume (L/g)
	double specificVolume=getSpecificVolume(msmsOutputFile,molecularWeight);
	// Calculate Solution Relative Dielectric (keep 2 significant digits for display purposes)
	double er;
	if(D.dielectric==-1){er=S.dielectric;}
	else{er=D.dielectric;}
	// Get Ion Data
	string largestIon=S.largestIon;	// should only keep 2 significant digits for display purposes
	string lastDigit=largestIon.substr(largestIon.length()-1,1);
	if(lastDigit.compare("0")==0){largestIon=formatNumberString(largestIon);}
	//largestIon=formatNumberString(largestIon);
	// iVal=strtod(largestIon.c_str(),NULL);
	//largestIon=cnvrtNumToStrng(iVal,2);
	//int iVal=largestIon.length();
	//pos=largestIon.find(".",0);
	//if(iVal-pos<3){largestIon+="0";}
	//<<"|"<<largestIon<<"|\n";
	string ionData=S.ionData;
	// Calculate Debye Length (A)
	double debyeLength=S.debyeLength;
	// Calculate Diffusivity
	long Result,fDesc,Result2;	
	double diffusivity,hydrodynamicRadius,slipPlanePosition;
	string hullradOutputFile=pdb_hydroFldr+"hullrad_output_"+calcNum+".txt";
	//string hullRadCMD="python2.7 "+theHullRad+" "+newPdbFile+" "+cnvrtNumToStrng(T,2)+" "+formatNumberString(cnvrtNumToStrng(solventViscosity,BUFFER_SIZE))+" "+formatNumberString(cnvrtNumToStrng(solventDensity,BUFFER_SIZE))+" > "+hullradOutputFile;
	string hullRadCMD=PYTHON_EXE;
	hullRadCMD+=" "+theHullRad+" "+newPdbFile+" "+cnvrtNumToStrng(T,2)+" "+formatNumberString(cnvrtNumToStrng(solventViscosity,BUFFER_SIZE))+" "+formatNumberString(cnvrtNumToStrng(solventDensity,BUFFER_SIZE))+" > "+hullradOutputFile;
	string theTxt;
	int sCode;
	if(D.inflationDistance==-1)
		{ATTEMPTING=true;
		if(D.RUN_HYDROPRO)
			{// Run HYDROPRO
			updateDisplayString("H",calcNum,2,nRows,dispLn);
			updateGuiDisplayString(".&nbsp;H&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;",calcNum,1,nRows,dispLn,guiDisplayFile);
			writeHydroproInput(D.hydroproCalcType,formatNumberString(cnvrtNumToStrng(solventDensity,BUFFER_SIZE)),formatNumberString(cnvrtNumToStrng(molecularWeight,BUFFER_SIZE)),D.outputTitle,hydroproInputFile,PDB+"_"+calcNum,specificVolume,T,solventViscosity);		
			PID=fork();
			if(PID==0)
				{// Run in Child Process, next time fork is executed it gives child process PID
				runHYDROPRO(PDB,pdb_hydroFldr,hydroproInputFile,logFile,calcNum);}
			else if(PID>0)
				{// Return to parent processs, fork should hold child process PID
				waitpid(PID,&status,0);}
			else if(PID==-1)
				{// fork() failed
				cerr<<"Fork failed."<<endl;exit(EXIT_FAILURE);}
			else{// Unpredicted Output
				cerr<<"Output of fork() unpredicted:\n"<<PID<<endl;}
			updateDisplayString(".",calcNum,2,nRows,dispLn);
			updateGuiDisplayString(".&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;",calcNum,1,nRows,dispLn,guiDisplayFile);
			// Remove Copied Temporary Executable File
			rChk=remove(tempXFile.c_str());
			if(rChk!=0){cerr<<"ERROR in runHYDROPRO\nCould not delete copied HYDROPRO Executable file."<<endl;exit(EXIT_FAILURE);}
			diffusivity=getDiffusivity(hydroproOutputFile);}
		else if(D.RUN_HULLRAD)
			{// Run HullRad
			updateDisplayString("H",calcNum,2,nRows,dispLn);
			//updateGuiDisplayString(".&nbsp;H&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;",calcNum,1,nRows,dispLn,guiDisplayFile);
			PID=fork();
			if(PID==0)
				{// Run in Child Process, next time fork is executed it gives child process PID
				sCode=system(hullRadCMD.c_str()); exit(EXIT_SUCCESS);
				//theTxt=encryptFile("/home/handsomedan/Desktop/ZPRED_III/editedHullRad.py","783e007c783e007c","Password");
				//fOut.open("/home/handsomedan/Desktop/ZPRED_III/theText.txt");
				//fOut<<theTxt;
				//fOut.close();
				}
			else if(PID>0)
				{// Return to parent processs, fork should hold child process PID
				waitpid(PID,&status,0);}
			else if(PID==-1)
				{// fork() failed
				cerr<<"Fork failed.\n";exit(EXIT_FAILURE);}
			else{// Unpredicted Output
				cerr<<"Output of fork() unpredicted:\n"<<PID<<endl;}
			updateDisplayString(".",calcNum,2,nRows,dispLn);
			//updateGuiDisplayString(".&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;",calcNum,1,nRows,dispLn,guiDisplayFile);
			diffusivity=getDiffusivity_hullrad(hullradOutputFile);
			}
		hydrodynamicRadius=calcHydrodynamicRadius(diffusivity,T,solventViscosity);
		slipPlanePosition=hydrodynamicRadius-proteinRadius;
		}

	// Determine Inflation Distance (keep 2 significant digits for display purposes)
	double inflationDistance;
	if(D.inflationDistance==-1)
		{if(slipPlanePosition<0){inflationDistance=0.00;}
		else{inflationDistance=slipPlanePosition;}
		if(isnan(inflationDistance)){cerr<<"ERROR defining inflationDistance!\nhydrodynamicRadius: "<<hydrodynamicRadius<<"\nproteinRadius: "<<proteinRadius<<endl;exit(EXIT_FAILURE);}}
	else{inflationDistance=D.inflationDistance;}	
	// Run PROPKA to hydrogenate & PDB2PQR to assign charge to PDB file	
	if(D.RUN_PDB2PQR)
		{updateDisplayString("P",calcNum,3,nRows,dispLn);
		updateGuiDisplayString(".&nbsp;.&nbsp;P&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;",calcNum,1,nRows,dispLn,guiDisplayFile);
		runPROPKA_PDB2PQR(pdbFile,pqrFile,calcNum,D.forceField,pdb_pqrFldr,propkaOutFile,pHVal,PDB2PQR,logFile);
		updateDisplayString(".",calcNum,3,nRows,dispLn);
		updateGuiDisplayString(".&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;",calcNum,1,nRows,dispLn,guiDisplayFile);}
	// Electrostatics with APBS
	apbsInputFile+=largestIon+"SD"+cnvrtNumToStrng(er,ER_SIG_FIGS)+".in";files2Delete+=apbsInputFile+delimiter;
	if(D.RUN_APBS)
		{updateDisplayString("A",calcNum,4,nRows,dispLn);
		updateGuiDisplayString(".&nbsp;.&nbsp;.&nbsp;A&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;",calcNum,1,nRows,dispLn,guiDisplayFile);
		//# Prepare IN Files for APBS with inputgen.py
		inFileEditer(D.apbsWriteTypes,D.numApbsWriteTypes,D.pdie,apbsGrdDim,pqrFile,apbsInputgen,ionData,formatNumberString(cnvrtNumToStrng(T,BUFFER_SIZE)),largestIon,cnvrtNumToStrng(er,ER_SIG_FIGS),fldrExt);
		// Run APBS to generate electrostatic fields
		PID=fork();
		if(PID==0){runAPBS(D.SC.apbsExecutable,apbsInputFile,logFile);}
		else if(PID>0){waitpid(PID,&status,0);}
		else if(PID==-1){cerr<<"Fork failed.\n";exit(EXIT_FAILURE);}
		else{cerr<<"Output of fork() unpredicted:\n"<<PID<<endl;}
		updateDisplayString(".",calcNum,4,nRows,dispLn);
		updateGuiDisplayString(".&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;",calcNum,1,nRows,dispLn,guiDisplayFile);}
	
	string sesFile="",sasFile="",csvFile="",potFile="",mvFile="",zpredSingleOutput="",shapeMsg="";
	ofstream zOut;
	ifstream zIn;
	int oldPos=0,ind=0;
	double Rh=0,ka=0,correction=0,pot=0,potSum=0,vol_frac=0,mobility=0,henryMobility=0,EmCompare=0,Em=0;
	float Q=0;
	bool RUNNING;
	//zpredOutFile+="SD"+cnvrtNumToStrng(er,ER_SIG_FIGS)+".txt";
	zpredOutFile+="PC"+cnvrtNumToStrng(D.proteinConc,ER_SIG_FIGS)+"_"+calcNum+".txt";

	// Add ZPRED Output File to List
	colIn.open(collideFile.c_str(),ifstream::app);
	if(colIn.fail()){cerr<<"ERROR in runZPRED!!!\nZPRED Output File List could not be updated.\n"<<collideFile<<endl;exit(EXIT_FAILURE);}
	colIn<<zpredOutFile+"\n";
	colIn.close();
	
	// Initialize MSMS Input Parameters
	updateDisplayString("M",calcNum,5,nRows,dispLn);
	define_msmsInput(PDB,newPdbFile,msmsFldr,pdb_msmsFldr,pdb_surfFldr,D.msms.surfPointDensity,D.msms.surfPointHiDensity,cnvrtNumToStrng(inflationDistance,DISTANCE_SIG_FIGS),msmsBinFile,msmsIn);
	// Generate Surface
	runMSMS(msmsIn,calcNum);
	// Define Error Values
	Error Er;
	// Determine if More Spherical or Cylindrical
	// Define Gyration Tensor gS
	double Rm,Rn,Value,Lx,Ly,Lz,Lsum,L,asphericity,sP,gyrationRadius,owMob,concCorrection;
	double** gS=new double*[3];
	for(int i=0;i<3;i++){gS[i]=new double[3];}
	// Solvated Radius
	Rh=proteinRadius+inflationDistance;
	// Define Electrokinetic Radius
	ka=Rh/debyeLength;
	shapeMsg="// Physical Shape Descriptor\n";
	//
	for(int i=0;i<3;i++)
		{for(int j=0;j<3;j++)
			{Value=0;
			for(int k=0;k<numLines;k++)
				{if(i==0)
					{// X Coordinate
					Rm=Data.x[k]-centerCoord.x;}
				else if(i==1)
					{// Y Coordinate
					Rm=Data.y[k]-centerCoord.y;}
				else if(i==2)
					{// Z Coordinate
					Rm=Data.z[k]-centerCoord.z;}
				if(j==0)
					{// X Coordinate
					Rn=Data.x[k]-centerCoord.x;}
				else if(j==1)
					{// Y Coordinate
					Rn=Data.y[k]-centerCoord.y;}
				else if(j==2)
					{// Z Coordinate
					Rn=Data.z[k]-centerCoord.z;}
				Value+=Rm*Rn;}
			Value/=numLines;
			gS[i][j]=Value;}}
	//for(int i=0;i<3;i++)
	//	{for(int j=0;j<3;j++)
	//		{cout<<S[i][j]<<" ";}
	//	cout<<"\n";}
	Lx=gS[0][0]; Ly=gS[1][1]; Lz=gS[2][2];
	Lsum=Lx+Ly+Lz;
	L=Lsum/3.0;
	// Define Asphericity
	asphericity=3*((Lx-L)*(Lx-L) + (Ly-L)*(Ly-L) + (Lz-L)*(Lz-L))/(2*Lsum*Lsum);	
	//cout<<"Asphericity:\t\t"<<asphericity<<endl;
	if(asphericity>=0.5){shape="cylinder";}
	// Define Shape Parameter
	sP=27*((Lx-L)*(Ly-L)*(Lz-L))/(Lsum*Lsum*Lsum);
	//cout<<"Shape Parameter:\t"<<sP<<endl;
	if(sP>=1){shape="cylinder";}
	// Define Gyration Radius
	gyrationRadius=sqrt(Lsum);
	//cout<<"Gyration Radius (Rg):\t"<<Rg<<" A"<<endl;
	double Q2;
	if(msmsIn.SUCCESS)
		{// Update Display String
		updateDisplayString("T",calcNum,5,nRows,dispLn);
		// MSMS Inflated Surface Output File
		msmsOutputFile=pdb_surfFldr+"fromMSMS/msmsOutput_SR"+cnvrtNumToStrng(inflationDistance,DISTANCE_SIG_FIGS)+"_"+calcNum+".txt"; files2Delete+=msmsOutputFile+delimiter;
		// MSMS VDW Area Output File
		msmsAreaFile=pdb_surfFldr+"fromMSMS/SR"+cnvrtNumToStrng(inflationDistance,DISTANCE_SIG_FIGS)+"_"+calcNum+".area"; files2Delete+=msmsAreaFile+delimiter;
		// MSMS VDW Face Output File
		msmsFaceFile=pdb_surfFldr+"fromMSMS/SR"+cnvrtNumToStrng(inflationDistance,DISTANCE_SIG_FIGS)+"_"+calcNum+".face"; files2Delete+=msmsFaceFile+delimiter;
		// Define Solvent Excluded Surface Vertices File
		sesFile=pdb_surfFldr+"fromMSMS/SR"+cnvrtNumToStrng(inflationDistance,DISTANCE_SIG_FIGS)+"_"+calcNum+".vert"; files2Delete+=sesFile+delimiter;
		sasFile=pdb_surfFldr+"SAS/SR"+cnvrtNumToStrng(inflationDistance,DISTANCE_SIG_FIGS)+"_"+calcNum+".vert"; files2Delete+=sasFile+delimiter;
		// Convert MSMS vert file to CSV
		csvFile=pdb_surfFldr+"CSV/SR"+cnvrtNumToStrng(inflationDistance,DISTANCE_SIG_FIGS)+"_"+calcNum+".csv"; files2Delete+=csvFile+delimiter;
		// Inflate SES to SAS
		inflateVert(inflationDistance,sesFile,sasFile,VERBOSE);
	
		vert2csv(sasFile,csvFile,VERBOSE);
		// Use APBS Multivalue Tool	
		// Define DX File Containing Computed Electric Potentials
		potFile=pdb_dxFldr+"pot/pH"+pHVal+"T"+formatNumberString(cnvrtNumToStrng(T,BUFFER_SIZE))+"SR"+largestIon+".dx"; files2Delete+=potFile+delimiter;
		// Define Output Multivalue File
		mvFile=pdb_surfFldr+"fromMV/pH"+pHVal+"T"+formatNumberString(cnvrtNumToStrng(T,BUFFER_SIZE))+"SR"+cnvrtNumToStrng(inflationDistance,DISTANCE_SIG_FIGS)+"_"+calcNum+".csv"; files2Delete+=mvFile+delimiter;
		updateDisplayString(".MV",calcNum,5,nRows,dispLn);		
		PID=fork();
		if(PID==0){runMULTIVALUE(D.SC.multivalueExecutable,csvFile,potFile,mvFile,logFile);}
		else if(PID>0){waitpid(PID,&status,0);}
		else if(PID==-1){cerr<<"Fork failed."<<endl;exit(EXIT_FAILURE);}
		else{cerr<<"Output of fork() unpredicted:\n"<<PID<<endl;}
		updateDisplayString(".Z",calcNum,6,nRows,dispLn);
		if(shape=="sphere")
			{correction=calc_HF_for_Sphere(ka);
			//shapeMsg+=".shape:\t\t\t\tGlobular\n";
			shapeMsg+=".shape: sphere\n";}
		else if(shape=="cylinder")
			{correction=calc_HF_for_Cylinder(ka);
			//shapeMsg+=".shape:\t\t\t\tFibrillar\n";
			shapeMsg+=".shape: cylinder\n";}
		zpredSingleOutput="//////////////////////////////////////////////////////\n";
		zpredSingleOutput+=shapeMsg;
		zpredSingleOutput+="// Asphericity\n";
		zpredSingleOutput+=".asphericity: "+formatNumberString(cnvrtNumToStrng(asphericity,BUFFER_SIZE))+"\n";
		zpredSingleOutput+="// Shape Parameter\n";
		zpredSingleOutput+=".shapeParameter: "+formatNumberString(cnvrtNumToStrng(sP,BUFFER_SIZE))+"\n";
		zpredSingleOutput+="// Gyration Radius [A]\n";
		zpredSingleOutput+=".gyrationRadius: "+formatNumberString(cnvrtNumToStrng(gyrationRadius,BUFFER_SIZE))+"\n";
		zpredSingleOutput+="// Solution Condition Number\n";
		zpredSingleOutput+=".solCondNum: "+D.id+"\n";
		zpredSingleOutput+="// Name of Protein Conformation\n";
		zpredSingleOutput+=".proteinName: "+PDB+"\n";
		zpredSingleOutput+="// Molecular Weight of Protein\n";
		zpredSingleOutput+=".molecularWeight: "+formatNumberString(cnvrtNumToStrng(calcProteinMolecularWeight(pdbFile),BUFFER_SIZE))+"\n";
		zpredSingleOutput+="// Anhydrous Protein Radius [A]\n";
		zpredSingleOutput+=".proteinRadius: "+formatNumberString(cnvrtNumToStrng(proteinRadius,BUFFER_SIZE))+"\n";
		zpredSingleOutput+="// Solvate (Hydrated) Protein Radius [A]\n";
		zpredSingleOutput+=".solvatedRadius: "+formatNumberString(cnvrtNumToStrng(Rh,BUFFER_SIZE))+"\n";
		//zpredSingleOutput+="// Shape Correction Factor\n";
		//zpredSingleOutput+=".shapeFactor: "+formatNumberString(cnvrtNumToStrng(shapeFactor,BUFFER_SIZE))+"\n";
		zpredSingleOutput+="// Protein Diffusivity [m^2/s]\n";
		zpredSingleOutput+=".diffusivity: "+formatNumberString(cnvrtNumToStrng(diffusivity,BUFFER_SIZE))+"\n";
		zpredSingleOutput+="// Solution Debye Length [A]\n";
		zpredSingleOutput+=".debyeLength: "+formatNumberString(cnvrtNumToStrng(debyeLength,BUFFER_SIZE))+"\n";
		zpredSingleOutput+="// Dimensionless Electrokinetic Radius\n";
		zpredSingleOutput+=".ka: "+formatNumberString(cnvrtNumToStrng(ka,BUFFER_SIZE))+"\n";
		zpredSingleOutput+="// Henry Correction Factor for Electrophoretic Retardation\n";
		zpredSingleOutput+=".henryCorrection: "+formatNumberString(cnvrtNumToStrng(correction,BUFFER_SIZE))+"\n";
		// Read Multivalue Vertices Output File
		zIn.open(mvFile.c_str());
		if(zIn.fail()){cerr<<"ERROR in runZPRED!!!\nMultivalue output file could not be opened.\n"<<mvFile<<endl;exit(EXIT_FAILURE);}
		// Skip Read Statement
		zIn.getline(Val,Sz);
		Counter=0;potSum=0;
		while(!zIn.eof())
			{value=Val;
			pos=value.rfind(",");
			tmp=value.substr(pos+1,value.length()-pos);
			pot=strtod(tmp.c_str(),NULL);
			if(!isnan(pot) && abs(pot)<1e9)
				{potSum+=pot;
				Counter++;}
			zIn.getline(Val,Sz);}
		zIn.close();
		// Acquire Protein Net Valence
		//Q=readPQRFile(pqrFile);
		Q2=readPROPKAFile(propkaOutFile);
		zpredSingleOutput+="// Protein Net Valence [e]\n";
		zpredSingleOutput+=".charge: "+formatNumberString(cnvrtNumToStrng(Q2,BUFFER_SIZE))+"\n";
		zpredSingleOutput+="// Hydration Layer Thickness [A]\n";
		zpredSingleOutput+=".Xsp: "+formatNumberString(cnvrtNumToStrng(inflationDistance,DISTANCE_SIG_FIGS))+"\n";
		// Undo Thermal Voltage and Take Average
		pot=potSum*kb*T/(Counter*ec);
		//cerr<<calcNum<<"|"<<pot<<" V\n";
		zpredSingleOutput+="// Electric Potential ("+formatNumberString(cnvrtNumToStrng(inflationDistance,DISTANCE_SIG_FIGS))+" A from surface) [V]\n";
		zpredSingleOutput+=".zetaPotential: "+formatNumberString(cnvrtNumToStrng(pot,BUFFER_SIZE))+"\n";
		zpredSingleOutput+="// Absolute Temperature [K]\n";
		zpredSingleOutput+=".temperature: "+formatNumberString(cnvrtNumToStrng(T,BUFFER_SIZE))+"\n";
		zpredSingleOutput+="// Solution pH\n";
		zpredSingleOutput+=".pH: "+pHVal+"\n";
		zpredSingleOutput+="// Solution Viscosity [Ns/m^2]\n";
		zpredSingleOutput+=".viscosity: "+formatNumberString(cnvrtNumToStrng(solventViscosity,BUFFER_SIZE))+"\n";
		zpredSingleOutput+="// Solution Relative Dielectric\n";
		zpredSingleOutput+=".dielectric: "+formatNumberString(cnvrtNumToStrng(er,BUFFER_SIZE))+"\n";
		zpredSingleOutput+="// Solution Density [kg/L]\n";
		zpredSingleOutput+=".density: "+formatNumberString(cnvrtNumToStrng(solventDensity,BUFFER_SIZE))+"\n";
		zpredSingleOutput+="// Solvent Components (Solvent Names)\n";
		zpredSingleOutput+=".solvent: "+solvStr+"\n";
		zpredSingleOutput+="// Solvent Components (Solvent Concentrations) [Mole Fraction]\n";
		zpredSingleOutput+=".solventConc: "+solvConcStr+"\n";
		zpredSingleOutput+="// Solvent Components (Solute Names)\n";
		zpredSingleOutput+=".solute: "+soluStr+"\n";
		zpredSingleOutput+="// Solvent Components (Solute Concentrations) [Molar]\n";
		zpredSingleOutput+=".soluteConc: "+soluConcStr+"\n";
		// Handle Ion Input
		In.numIons=count_delimiter(ionData,",")/3;
		In.conc=new double[In.numIons];
		In.valence=new int[In.numIons];
		In.mobility=new double[In.numIons];
		In.radius=new double[In.numIons];

		oldPos=0;
		for(int i=0;i<In.numIons;i++)
			{// Acquire valence
			pos=ionData.find(",",oldPos);tmp=ionData.substr(oldPos,pos-oldPos);In.valence[i]=atoi(tmp.c_str());oldPos=pos+1;
			// Acquire Concentration (Molar)		
			pos=ionData.find(",",oldPos);tmp=ionData.substr(oldPos,pos-oldPos);In.conc[i]=strtod(tmp.c_str(),NULL);oldPos=pos+1;
			// Acquire Ion radius (Angstroms)		
			pos=ionData.find(",",oldPos);tmp=ionData.substr(oldPos,pos-oldPos);In.radius[i]=strtod(tmp.c_str(),NULL);oldPos=pos+1;}

		// Electrophoretic Mobility (umcm/(Vs)) from O'Brien & White Algorithm
		// Absolute error is a threshold error value representing the acceptable error as the value of the measured state approaches
		// Relative error represents a percentage of the state's value (default: 1e-3 means the computed value will be accurate to within 0.1%
		Er.absolute=1e-12;Er.relative=1e-7;Er.a_x=1;Er.a_dxdt=1;
		freopen(logFile.c_str(),"a",stdout);
		switch(S.numIon)
			{case 1:
				break;
			case 2:
				owMob=two_ion_problem(pot,Rh*1e-10,S,Er,true);
				break;		
			case 3:
				owMob=three_ion_problem(pot,Rh*1e-10,S,Er,true);
				break;
			case 4:
				owMob=four_ion_problem(pot,Rh*1e-10,S,Er,true);
				break;
			case 5:
				owMob=five_ion_problem(pot,Rh*1e-10,S,Er,true);
				break;
			case 6:
				owMob=six_ion_problem(pot,Rh*1e-10,S,Er,true);
				break;
			case 7:
				owMob=seven_ion_problem(pot,Rh*1e-10,S,Er,true);
				break;
			case 8:
				owMob=eight_ion_problem(pot,Rh*1e-10,S,Er,true);
				break;
			case 9:
				owMob=nine_ion_problem(pot,Rh*1e-10,S,Er,true);
				break;
			default:
				break;
			}
	
		zpredSingleOutput+="// Specific Volume [L/g]\n";
		zpredSingleOutput+=".specificVolume: "+formatNumberString(cnvrtNumToStrng(specificVolume,BUFFER_SIZE))+"\n";
		vol_frac=specificVolume*D.proteinConc;
		zpredSingleOutput+="// Volume Fraction\n";
		zpredSingleOutput+=".volumeFraction: "+formatNumberString(cnvrtNumToStrng(vol_frac,BUFFER_SIZE))+"\n";
		// Calculate Electrophoretic Mobility umcm/(Vs) [still need to apply correction factor to mobility conversion]
		mobility=eo*er*pot*correction/solventViscosity;
		//
		zpredSingleOutput+="// Electrophoretic Mobility [umcm/(Vs)] Defined from ZPRED computed zeta potential and Henry Correction Factor\n";
		zpredSingleOutput+=".henryMobility: "+formatNumberString(cnvrtNumToStrng(mobility*1e8,BUFFER_SIZE))+"\n";
		//
		// Calculate Protein Concentration Correction Factor (in case of high concentrations)
		concCorrection=calc_Kuwabara_Correction_for_Sphere(ka,vol_frac,pot,Rh*1e-10);
		zpredSingleOutput+="// Electrophoretic Mobility [umcm/(Vs)] Defined from ZPRED computed zeta potential and Kuwabara Correction Factor\n";
		zpredSingleOutput+=".kuwabaraMobility: "+formatNumberString(cnvrtNumToStrng(mobility*1e8*concCorrection/correction,BUFFER_SIZE))+"\n";
		//
		zpredSingleOutput+="// Kuwabara Correction Factor for Particle Concentration\n";
		zpredSingleOutput+=".kuwabaraCorrection: "+formatNumberString(cnvrtNumToStrng(concCorrection,BUFFER_SIZE))+"\n";
		// Calculate Electrophoretic Mobility from Henry Eqn for comparison
		//henryMobility=Q*ec*correction/(4*pi*solventViscosity*Rh*1e-10*(1+ka)*shapeFactor);
		henryMobility=Q*ec*correction/(4*pi*solventViscosity*Rh*1e-10*(1+ka));
		zpredSingleOutput+="// Electrophoretic Mobility from Henry Equation [umcm/(Vs)]\n";
		zpredSingleOutput+=".henryModelMobility: "+formatNumberString(cnvrtNumToStrng(henryMobility*1e8,BUFFER_SIZE))+"\n";
		// from OW Algorithm
		if(isnan(owMob))
			{// OW Algorithm Failed
			zpredSingleOutput+="// Electrophoretic Mobility [umcm/(Vs)] Defined from ZPRED computed zeta potential and Henry Correction Factor\n";
			zpredSingleOutput+=".semMobility: NaN\n";
			// Concentration Corrected OW Algorithm Output
			zpredSingleOutput+="// Electrophoretic Mobility [umcm/(Vs)] Defined from ZPRED computed zeta potential and Kuwabara Correction Factor\n";
			zpredSingleOutput+=".semMobility2: NaN\n";
			}
		else
			{// OW Algorithm Successfully Completed
			zpredSingleOutput+="// Electrophoretic Mobility from O'Brien & White Algorithm (solves standard electrokinetic model (SEM)) [umcm/(Vs)]\n";
			zpredSingleOutput+=".semMobility: "+formatNumberString(cnvrtNumToStrng(owMob*1e8,BUFFER_SIZE))+"\n";
			// Concentration Corrected OW Algorithm Output
			zpredSingleOutput+="// Electrophoretic Mobility from O'Brien & White Algorithm with Kuwabara Correction Factor [umcm/(Vs)]\n";
			zpredSingleOutput+=".semMobility2: "+formatNumberString(cnvrtNumToStrng(owMob*1e8*concCorrection/correction,BUFFER_SIZE))+"\n";}
		zpredSingleOutput+="//////////////////////////////////////////////////////\n";
		// Define Dimensionless Electrophoretic Mobility
		Em=3*solventViscosity*ec*abs(mobility);
		Em/=2*er*eo*kb*T;
		zpredSingleOutput+="ka:\t"+formatNumberString(cnvrtNumToStrng(ka,BUFFER_SIZE))+"\n";			
		zpredSingleOutput+="Em:\t"+formatNumberString(cnvrtNumToStrng(Em,BUFFER_SIZE))+"\n";
		// Determine Appropriate Electrokinetic Model
		RUNNING=true;ind=0;
		while(RUNNING)
			{if(ka<kaVals[ind]){RUNNING=false;}
			else{ind++;}
			if(ind>131){ind--;RUNNING=false;}}
	
		if(ind==0){EmCompare=interpolate(ka,kaVals[ind],kaVals[ind+1],EmVals[ind],EmVals[ind+1]);}
		else{EmCompare=interpolate(ka,kaVals[ind-1],kaVals[ind],EmVals[ind-1],EmVals[ind]);}
		zpredSingleOutput+="Compare Em to:\t"+formatNumberString(cnvrtNumToStrng(EmCompare,BUFFER_SIZE))+"\n";
		if(Em<EmCompare){zpredSingleOutput+="Henry Equation Appropriate!!!\n";}
		else{if(ka<10)
				{zpredSingleOutput+="Ohshima Approximation Appropriate!!!\n";
				correction=calc_OF_for_Sphere(In,er,T,pot,ka);
				zpredSingleOutput+="Ohshima Correction: "+formatNumberString(cnvrtNumToStrng(correction,BUFFER_SIZE))+"\n";
				mobility=eo*er*pot*correction/solventViscosity;
				zpredSingleOutput+="Electrophoretic Mobility: "+formatNumberString(cnvrtNumToStrng(mobility*1e8,BUFFER_SIZE))+" umcm/(Vs)\n";}
			else{zpredSingleOutput+="Ohshima-Healy-White Approximation Appropriate!!!\n";				
				correction=calc_OHWF_for_Sphere(In,er,T,pot,ka);
				zpredSingleOutput+="Ohshima-Healy-White Correction: "+formatNumberString(cnvrtNumToStrng(correction,BUFFER_SIZE))+"\n";
				mobility=eo*er*pot*correction/solventViscosity;	
				zpredSingleOutput+="Electrophoretic Mobility: "+formatNumberString(cnvrtNumToStrng(mobility*1e8,BUFFER_SIZE))+" umcm/(Vs)\n";
				}}
		// Generate Electric Potential Profile?
		int profileLen=35;
		double* position=new double[profileLen];
		position[0]=0.1;
		position[1]=0.2;
		position[2]=0.3;
		position[3]=0.4;
		position[4]=0.5;
		position[5]=0.6;
		position[6]=0.7;
		position[7]=0.8;
		position[8]=0.9;
		position[9]=1.0;
		position[10]=1.25;
		position[11]=1.5;
		position[12]=1.75;
		position[13]=2.0;
		position[14]=2.25;
		position[15]=2.5;
		position[16]=2.75;
		position[17]=3.0;
		position[18]=3.25;
		position[19]=3.5;
		position[20]=3.75;
		position[21]=4.0;
		position[22]=4.25;
		position[23]=4.5;
		position[24]=4.75;
		position[25]=5.0;
		position[26]=5.5;
		position[27]=6.0;
		position[28]=6.5;
		position[29]=7.0;
		position[30]=7.5;
		position[31]=8.0;
		position[32]=8.5;
		position[33]=9.0;
		position[34]=9.5;
		if(D.GENERATE_POT_PROFILE)
			{zpredSingleOutput+="// Electric Potential Profile from Protein Surface\n";
			zpredSingleOutput+="Distance [A] | Potential [V]\n";
			zpredSingleOutput+=".generateProfile: ";
			for(int i=0;i<profileLen;i++)
				{define_msmsInput(PDB,newPdbFile,msmsFldr,pdb_msmsFldr,pdb_surfFldr,D.msms.surfPointDensity,D.msms.surfPointHiDensity,cnvrtNumToStrng(position[i],DISTANCE_SIG_FIGS),msmsBinFile,msmsIn);
				// Generate Surface
				runMSMS(msmsIn,calcNum);
				if(msmsIn.SUCCESS)
					{// Surface Generated Successfully
					// MSMS Inflated Surface Output File
					msmsOutputFile=pdb_surfFldr+"fromMSMS/msmsOutput_SR"+cnvrtNumToStrng(position[i],DISTANCE_SIG_FIGS)+"_"+calcNum+".txt"; files2Delete+=msmsOutputFile+delimiter;
					// MSMS VDW Area Output File
					msmsAreaFile=pdb_surfFldr+"fromMSMS/SR"+cnvrtNumToStrng(position[i],DISTANCE_SIG_FIGS)+"_"+calcNum+".area"; files2Delete+=msmsAreaFile+delimiter;
					// MSMS VDW Face Output File
					msmsFaceFile=pdb_surfFldr+"fromMSMS/SR"+cnvrtNumToStrng(position[i],DISTANCE_SIG_FIGS)+"_"+calcNum+".face"; files2Delete+=msmsFaceFile+delimiter;
					// Define Solvent Excluded Surface Vertices File
					sesFile=pdb_surfFldr+"fromMSMS/SR"+cnvrtNumToStrng(position[i],DISTANCE_SIG_FIGS)+"_"+calcNum+".vert"; files2Delete+=sesFile+delimiter;
					sasFile=pdb_surfFldr+"SAS/SR"+cnvrtNumToStrng(position[i],DISTANCE_SIG_FIGS)+"_"+calcNum+".vert"; files2Delete+=sasFile+delimiter;
					// Inflate SES to SAS
					inflateVert(position[i],sesFile,sasFile,VERBOSE);
					// Convert MSMS vert file to CSV
					csvFile=pdb_surfFldr+"CSV/SR"+cnvrtNumToStrng(position[i],DISTANCE_SIG_FIGS)+"_"+calcNum+".csv"; files2Delete+=csvFile+delimiter;
					vert2csv(sasFile,csvFile,VERBOSE);
					// Use APBS Multivalue Tool	
					// Define DX File Containing Computed Electric Potentials
					//potFile=pdb_dxFldr+"pot/pH"+pHVal+"T"+formatNumberString(cnvrtNumToStrng(T,BUFFER_SIZE))+"SR"+largestIon+".dx";
					// Define Output Multivalue File
					mvFile=pdb_surfFldr+"fromMV/pH"+pHVal+"T"+formatNumberString(cnvrtNumToStrng(T,BUFFER_SIZE))+"SR"+cnvrtNumToStrng(position[i],DISTANCE_SIG_FIGS)+"_"+calcNum+".csv"; files2Delete+=mvFile+delimiter;
					//updateDisplayString("MV",calcNum,6,nRows,dispLn);
					PID=fork();
					if(PID==0){runMULTIVALUE(D.SC.multivalueExecutable,csvFile,potFile,mvFile,logFile);}
					else if(PID>0){waitpid(PID,&status,0);}
					else if(PID==-1){cerr<<"Fork failed."<<endl;exit(EXIT_FAILURE);}
					else{cerr<<"Output of fork() unpredicted:\n"<<PID<<endl;}
					// Read Multivalue Vertices Output File
					zIn.open(mvFile.c_str());
					if(zIn.fail()){cerr<<"ERROR in runZPRED!!!\nMultivalue output file could not be opened.\n"<<mvFile<<endl;exit(EXIT_FAILURE);}
					// Skip Read Statement
					zIn.getline(Val,Sz);
					Counter=0;potSum=0;
					while(!zIn.eof())
						{value=Val;
						pos=value.rfind(",");
						tmp=value.substr(pos+1,value.length()-pos);
						pot=strtod(tmp.c_str(),NULL);
						if(!isnan(pot) && abs(pot)<1e9)
							{potSum+=pot;
							Counter++;}
						zIn.getline(Val,Sz);}
					zIn.close();
					// Undo Thermal Voltage and Take Average
					pot=potSum*kb*T/(Counter*ec);
					// Electric Potential Profile
					zpredSingleOutput+=cnvrtNumToStrng(position[i],DISTANCE_SIG_FIGS)+","+formatNumberString(cnvrtNumToStrng(pot,BUFFER_SIZE))+",;";
					}
				else
					{// Surface Generation Failed, provide NaN values
					// Electric Potential Profile
					zpredSingleOutput+=cnvrtNumToStrng(position[i],DISTANCE_SIG_FIGS)+",NaN,;";
					}
				}
			zpredSingleOutput+="\n";
			}
		else
			{zpredSingleOutput+="// Electric Potential Profile from Protein Surface\n";
			zpredSingleOutput+=".generateProfile: false\n";}
		// Output Zeta Potential Computation results
		zOut.open(zpredOutFile.c_str());
		if(zOut.fail()){cerr<<"ERROR in runZPRED!!!\nZPRED output file could not be opened.\n"<<zpredOutFile<<endl;exit(EXIT_FAILURE);}
		zOut<<zpredSingleOutput;
		zOut.close();
		clear_ion_Data(In);
		updateDisplayString(".",calcNum,7,nRows,dispLn);
		//updateGuiDisplayString(".&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;.&nbsp;",calcNum,1,nRows,dispLn,guiDisplayFile);
		updateDisplayString("COMPLETE",calcNum,1,nRows,dispLn);
		//updateGuiDisplayString("COMPLETE",calcNum,1,nRows,dispLn,guiDisplayFile);

		// Delete Computation Files (when saveTemps=false)
		string *deletedFiles;
		int N,Err;
		char errMsg[256];
		if(!D.saveTemps)
			{N=count_delimiter(files2Delete,delimiter);
			deletedFiles=fill_string_array(files2Delete,N,delimiter);
			for(int i=0;i<N;i++)
				{tmp=deletedFiles[i];
				rChk=remove(tmp.c_str());
				if(rChk!=0)
					{//Err=errno;erase_display();
					//cerr<<"ERROR in runZPRED\nCould not delete temp file.\n"<<tmp;
					//cerr<<strerror_r(Err,errMsg,256)<<"\n";
					/*exit(EXIT_FAILURE);*/}}
			}
		}
	else
		{// Update Display String
		updateDisplayString("F",calcNum,5,nRows,dispLn);
		if(shape=="sphere")
			{correction=calc_HF_for_Sphere(ka);
			//shapeMsg+=".shape:\t\t\t\tGlobular\n";
			shapeMsg+=".shape: sphere\n";}
		else if(shape=="cylinder")
			{correction=calc_HF_for_Cylinder(ka);
			//shapeMsg+=".shape:\t\t\t\tFibrillar\n";
			shapeMsg+=".shape: cylinder\n";}
		zpredSingleOutput="//////////////////////////////////////////////////////\n";
		zpredSingleOutput+=shapeMsg;
		zpredSingleOutput+="// Asphericity\n";
		zpredSingleOutput+=".asphericity: "+formatNumberString(cnvrtNumToStrng(asphericity,BUFFER_SIZE))+"\n";
		zpredSingleOutput+="// Shape Parameter\n";
		zpredSingleOutput+=".shapeParameter: "+formatNumberString(cnvrtNumToStrng(sP,BUFFER_SIZE))+"\n";
		zpredSingleOutput+="// Gyration Radius [A]\n";
		zpredSingleOutput+=".gyrationRadius: "+formatNumberString(cnvrtNumToStrng(gyrationRadius,BUFFER_SIZE))+"\n";
		zpredSingleOutput+="// Solution Condition Number\n";
		zpredSingleOutput+=".solCondNum: "+D.id+"\n";
		zpredSingleOutput+="// Name of Protein Conformation\n";
		zpredSingleOutput+=".proteinName: "+PDB+"\n";
		zpredSingleOutput+="// Anhydrous Protein Radius [A]\n";
		zpredSingleOutput+=".proteinRadius: "+formatNumberString(cnvrtNumToStrng(proteinRadius,BUFFER_SIZE))+"\n";
		zpredSingleOutput+="// Solvate (Hydrated) Protein Radius [A]\n";
		zpredSingleOutput+=".solvatedRadius: "+formatNumberString(cnvrtNumToStrng(Rh,BUFFER_SIZE))+"\n";
		//zpredSingleOutput+="// Shape Correction Factor\n";
		//zpredSingleOutput+=".shapeFactor: "+formatNumberString(cnvrtNumToStrng(shapeFactor,BUFFER_SIZE))+"\n";
		zpredSingleOutput+="// Protein Diffusivity [m^2/s]\n";
		zpredSingleOutput+=".diffusivity: "+formatNumberString(cnvrtNumToStrng(diffusivity,BUFFER_SIZE))+"\n";
		zpredSingleOutput+="// Solution Debye Length [A]\n";
		zpredSingleOutput+=".debyeLength: "+formatNumberString(cnvrtNumToStrng(debyeLength,BUFFER_SIZE))+"\n";
		zpredSingleOutput+="// Dimensionless Electrokinetic Radius\n";
		zpredSingleOutput+=".ka: "+formatNumberString(cnvrtNumToStrng(ka,BUFFER_SIZE))+"\n";
		zpredSingleOutput+="// Henry Correction Factor for Electrophoretic Retardation\n";
		zpredSingleOutput+=".henryCorrection: "+formatNumberString(cnvrtNumToStrng(correction,BUFFER_SIZE))+"\n";
		// Acquire Protein Net Valence
		//Q=readPQRFile(pqrFile);
		Q2=readPROPKAFile(propkaOutFile);
		zpredSingleOutput+="// Protein Net Valence [e]\n";
		zpredSingleOutput+=".charge: "+formatNumberString(cnvrtNumToStrng(Q2,BUFFER_SIZE))+"\n";
		zpredSingleOutput+="// Hydration Layer Thickness [A]\n";
		zpredSingleOutput+=".Xsp: "+formatNumberString(cnvrtNumToStrng(inflationDistance,DISTANCE_SIG_FIGS))+"\n";
		zpredSingleOutput+="// Electric Potential ("+formatNumberString(cnvrtNumToStrng(inflationDistance,DISTANCE_SIG_FIGS))+" A from surface) [V]\n";
		zpredSingleOutput+=".zetaPotential: NaN\n";
		zpredSingleOutput+="// Absolute Temperature [K]\n";
		zpredSingleOutput+=".temperature: "+formatNumberString(cnvrtNumToStrng(T,BUFFER_SIZE))+"\n";
		zpredSingleOutput+="// Solution pH\n";
		zpredSingleOutput+=".pH: "+pHVal+"\n";
		zpredSingleOutput+="// Solution Viscosity [Ns/m^2]\n";
		zpredSingleOutput+=".viscosity: "+formatNumberString(cnvrtNumToStrng(solventViscosity,BUFFER_SIZE))+"\n";
		zpredSingleOutput+="// Solution Relative Dielectric\n";
		zpredSingleOutput+=".dielectric: "+formatNumberString(cnvrtNumToStrng(er,BUFFER_SIZE))+"\n";
		zpredSingleOutput+="// Solution Density [kg/L]\n";
		zpredSingleOutput+=".density: "+formatNumberString(cnvrtNumToStrng(solventDensity,BUFFER_SIZE))+"\n";
		zpredSingleOutput+="// Solvent Components (Solvent Names)\n";
		zpredSingleOutput+=".solvent: "+solvStr+"\n";
		zpredSingleOutput+="// Solvent Components (Solvent Concentrations) [Mole Fraction]\n";
		zpredSingleOutput+=".solventConc: "+solvConcStr+"\n";
		zpredSingleOutput+="// Solvent Components (Solute Names)\n";
		zpredSingleOutput+=".solute: "+soluStr+"\n";
		zpredSingleOutput+="// Solvent Components (Solute Concentrations) [Molar]\n";
		zpredSingleOutput+=".soluteConc: "+soluConcStr+"\n";
		// Handle Ion Input
		In.numIons=count_delimiter(ionData,",")/3;
		In.conc=new double[In.numIons];
		In.valence=new int[In.numIons];
		In.mobility=new double[In.numIons];
		In.radius=new double[In.numIons];

		oldPos=0;
		for(int i=0;i<In.numIons;i++)
			{// Acquire valence
			pos=ionData.find(",",oldPos);tmp=ionData.substr(oldPos,pos-oldPos);In.valence[i]=atoi(tmp.c_str());oldPos=pos+1;
			// Acquire Concentration (Molar)		
			pos=ionData.find(",",oldPos);tmp=ionData.substr(oldPos,pos-oldPos);In.conc[i]=strtod(tmp.c_str(),NULL);oldPos=pos+1;
			// Acquire Ion radius (Angstroms)		
			pos=ionData.find(",",oldPos);tmp=ionData.substr(oldPos,pos-oldPos);In.radius[i]=strtod(tmp.c_str(),NULL);oldPos=pos+1;}
		zpredSingleOutput+="// Specific Volume [L/g]\n";
		zpredSingleOutput+=".specificVolume: "+formatNumberString(cnvrtNumToStrng(specificVolume,BUFFER_SIZE))+"\n";
		vol_frac=specificVolume*D.proteinConc;
		zpredSingleOutput+="// Volume Fraction\n";
		zpredSingleOutput+=".volumeFraction: "+formatNumberString(cnvrtNumToStrng(vol_frac,BUFFER_SIZE))+"\n";
		// Calculate Electrophoretic Mobility umcm/(Vs) [still need to apply correction factor to mobility conversion]
		//mobility=eo*er*pot*correction/solventViscosity;
		zpredSingleOutput+="// Electrophoretic Mobility [umcm/(Vs)] Defined from ZPRED computed zeta potential and Henry Correction Factor\n";
		zpredSingleOutput+=".henryMobility: NaN\n";
		// Calculate Protein Concentration Correction Factor (in case of high concentrations)
		//concCorrection=calc_Kuwabara_Correction_for_Sphere(ka,vol_frac,pot,Rh*1e-10);
		zpredSingleOutput+="// Electrophoretic Mobility [umcm/(Vs)] Defined from ZPRED computed zeta potential and Kuwabara Correction Factor\n";
		zpredSingleOutput+=".kuwabaraMobility: NaN\n";
		//
		zpredSingleOutput+="// Kuwabara Correction Factor for Particle Concentration\n";
		zpredSingleOutput+=".kuwabaraCorrection: NaN\n";
		// Calculate Electrophoretic Mobility from Henry Eqn for comparison
		//henryMobility=Q*ec*correction/(4*pi*solventViscosity*Rh*1e-10*(1+ka)*shapeFactor);
		henryMobility=Q*ec*correction/(4*pi*solventViscosity*Rh*1e-10*(1+ka));
		zpredSingleOutput+="// Electrophoretic Mobility from Henry Equation [umcm/(Vs)]\n";
		zpredSingleOutput+=".henryModelMobility: "+formatNumberString(cnvrtNumToStrng(henryMobility*1e8,BUFFER_SIZE))+"\n";
		zpredSingleOutput+="// Electrophoretic Mobility [umcm/(Vs)] Defined from ZPRED computed zeta potential and Henry Correction Factor\n";
		zpredSingleOutput+=".semMobility: NaN\n";
		// Concentration Corrected OW Algorithm Output
		zpredSingleOutput+="// Electrophoretic Mobility [umcm/(Vs)] Defined from ZPRED computed zeta potential and Kuwabara Correction Factor\n";
		zpredSingleOutput+=".semMobility2: NaN\n";
		zpredSingleOutput+="//////////////////////////////////////////////////////\n";
		// Generate Electric Potential Profile?
		zpredSingleOutput+="// Electric Potential Profile from Protein Surface\n";
		zpredSingleOutput+=".generateProfile: false\n";
		// Output Zeta Potential Computation results
		zOut.open(zpredOutFile.c_str());
		if(zOut.fail()){cerr<<"ERROR in runZPRED!!!\nZPRED output file could not be opened.\n"<<zpredOutFile<<endl;exit(EXIT_FAILURE);}
		zOut<<zpredSingleOutput;
		zOut.close();
		clear_ion_Data(In);
		updateDisplayString(".",calcNum,7,nRows,dispLn);
		updateDisplayString("COMPLETE",calcNum,1,nRows,dispLn);
		// Delete Computation Files (when saveTemps=false)
		string *deletedFiles;
		int N,Err;
		char errMsg[256];
		if(!D.saveTemps)
			{N=count_delimiter(files2Delete,delimiter);
			deletedFiles=fill_string_array(files2Delete,N,delimiter);
			for(int i=0;i<N;i++)
				{tmp=deletedFiles[i];
				rChk=remove(tmp.c_str());
				}
			}
		}
	return 0;}

// Ref: 1974. The Prediction of Electrokinetic Phenomena within Multiparticle Systems I. Electrophoresis and Electroosmosis
double calc_Kuwabara_Correction_for_Sphere(double ka,double volFrac,double Zeta,double a)
	{double Output;
	// Solvated Radius: a [m]
	// Zeta: [V]
	double exponent=1.0/3.0;
	//cerr<<"exp: "<<exponent<<"\n";
	double y=pow(volFrac,exponent);
	double y3=volFrac;		//y*y*y;
	//cerr<<"y: "<<y<<"\n";
	double b=a/y;
	// Extract Debye Length [m]
	double k=a/ka;
	double k_b=b/k;
	double ka2=ka*ka;
	double ka3=ka2*ka;
	double ka4=ka3*ka;
	double ka5=ka4*ka;
	double ka6=ka5*ka;

	double R=1/(sinh(k_b-ka)-k_b*cosh(k_b-ka));
	//cerr<<"R: "<<R<<"\n";
	double Q=(1-k_b*tanh(k_b-ka))/(tanh(k_b-ka)-k_b);
	//cerr<<"Q: "<<Q<<"\n";
	double A=R*(sinh(k_b)-k_b*cosh(k_b));
	//cerr<<"A: "<<A<<"\n";
	double B=R*(cosh(k_b)-k_b*sinh(k_b));
	//cerr<<"B: "<<B<<"\n";
	// Integration Setup for computing I variable
	int numPnts=500;
	double dt=(k_b-ka)/(double(numPnts)-1);
	//cerr<<"dt: "<<dt<<"\n";
	string tList="",iList="";
	double t=ka,iVal;
	for(int i=0;i<numPnts;i++)
		{iVal=(A*cosh(t)-B*sinh(t))/t;
		//cerr<<t<<" "<<iVal<<endl;
		tList+=formatNumberString(cnvrtNumToStrng(t,BUFFER_SIZE))+";";
		iList+=formatNumberString(cnvrtNumToStrng(iVal,BUFFER_SIZE))+";";
		t=t+dt;}
	double *t1=fill_double_array(tList,numPnts,";");
	double *i1=fill_double_array(iList,numPnts,";");
	double iArea=calcAreaUnderCurve(t1,i1,numPnts);
	//cerr<<"I: "<<iArea<<"\n";
	double trm1=1 + R*ka + ka2/16.0 - ka3*(5*Q/48.0 - R*y/4.0 + R*y3/12.0) - ka4/96.0;
	double trm2=ka5*(Q/96.0 - R*y/48.0) + iArea*(ka4/8.0 - ka6/96.0);
	double trm3=-Zeta*(trm1+trm2);
	//cerr<<"1st term: "<<trm3<<"\n";
	trm1=1 + 3/ka2 + 3*Q/ka + ka*(R/y3 - R/10.0 + Q/10.0) - ka2/40.0;
	trm2=ka3*(Q/120.0 - R*y3/30.0) - ka4/240.0 + ka5*(Q/240.0 - R*y/120.0) - iArea*ka6/240.0;
	double trm4=-Zeta*(trm1+trm2);
	//cerr<<"2nd term: "<<trm4<<"\n";
	double trm5=3*(1-y3)*Zeta/2.0;

	Output=(-trm3+y3*trm4)/trm5;

	return Output;}

double calcAreaUnderCurve(double* x,double* y,int N)
	{double Output=0,area=0,xIntercept=0;
	Vec2d p1,p2,p3,p4;
	for(int i=0;i<N-1;i++)
		{// Determine whether to add or subtract
		if(!isnan(y[i])&&!isnan(y[i+1]))
			{if(y[i]>=0 && y[i+1]>=0)
				{// Both Data Points Above X-Axis, so ADD area to total sum
				// Start at X-Axis
				p1.x=x[i];p1.y=0;
				// Define First Point
				p2.x=x[i];p2.y=y[i];
				// Define Next Point
				p3.x=x[i+1];p3.y=y[i+1];
				// Reconnect at X-Axis
				p4.x=x[i+1];p4.y=0;
				area=shoelace(p1,p2,p3,p4);
				if(!isnan(area)){Output+=area;}}
			else if(y[i]<0 && y[i+1]>=0)
				{// First Point is Negative, while Second is Positive
				// Define negative area
				p1.x=x[i];p1.y=0;
				p2.x=x[i];p2.y=-y[i];
				// Find X-Intercept
				xIntercept=interpolate(0,y[i],y[i+1],x[i],x[i+1]);
				p3.x=xIntercept;p3.y=0;
				area=-shoelace(p1,p2,p3);
				// Define positive area
				p1.x=xIntercept;p1.y=0;
				p2.x=x[i+1];p2.y=y[i+1];
				p3.x=x[i+1];p3.y=0;
				area+=shoelace(p1,p2,p3);
				if(!isnan(area)){Output+=area;}}
			else if(y[i]>=0 && y[i+1]<0)
				{// First Point is Positive, while Second is Negative
				// Define positive area
				p1.x=x[i];p1.y=0;
				p2.x=x[i];p2.y=y[i];
				// Find X-Intercept
				xIntercept=interpolate(0,y[i],y[i+1],x[i],x[i+1]);
				p3.x=xIntercept;p3.y=0;
				area=shoelace(p1,p2,p3);
				// Define negative area
				p1.x=xIntercept;p1.y=0;
				p2.x=x[i+1];p2.y=-y[i+1];
				p3.x=x[i+1];p3.y=0;
				area-=shoelace(p1,p2,p3);
				if(!isnan(area)){Output+=area;}}
			else if(y[i]<0 && y[i+1]<0)
				{// Both Data Points Below X-Axis, so SUBTRACT area to total sum
				// Start at X-Axis
				p1.x=x[i];p1.y=0;
				// Define First Point
				p2.x=x[i];p2.y=-y[i];
				// Define Next Point
				p3.x=x[i+1];p3.y=-y[i+1];
				// Reconnect at X-Axis
				p4.x=x[i+1];p4.y=0;
				area=shoelace(p1,p2,p3,p4);
				if(!isnan(area)){Output-=area;}}
			}
		}
	return Output;}

double shoelace(Vec2d p1,Vec2d p2,Vec2d p3,Vec2d p4)
	{// Computes area by shoelace formula for 4 point polygon
	double area=p1.x*p2.y+p2.x*p3.y+p3.x*p4.y+ p4.x*p1.y - p2.x*p1.y-p3.x*p2.y-p4.x*p3.y - p1.x*p4.y;
	area/=2;
	area=abs(area);
	return area;}

double shoelace(Vec2d p1,Vec2d p2,Vec2d p3)
	{// Computes area by shoelace formula for 3 point polygon
	double area=p1.x*p2.y+p2.x*p3.y+ p3.x*p1.y - p2.x*p1.y-p3.x*p2.y - p1.x*p3.y;
	area/=2;
	area=abs(area);
	return area;}

void clear_ion_Data(ion_Data& In)
	{In.numIons=0;
	delete [] In.conc;
	delete [] In.valence;
	delete [] In.mobility;
	delete [] In.radius;}

void updateDisplayString(string dispStr,string calcNum,int updatePos,int nRows,string& dispLn)
	{// Restore Cursor Position to top of display
	string setPos="\033[u";
	string clearLn="\033[K";
	// Move Cursor to Line Corresponding to Current Computation
	int cN=atoi(calcNum.c_str());cN--;
	int lnPos=cN%nRows;
	// Move to Line with Current Computation
	string dwnLn="\033["+cnvrtNumToStrng(lnPos,0)+"B";
	// Find the Current Computation in the Display String Line
	string dispLabel="["+calcNum+"]";
	int pos=dispLn.find(dispLabel,0);
	string before="",after="";
	if(pos==string::npos){cerr<<"Error in updateDisplayString!\nCould not find current computation."<<endl;exit(EXIT_FAILURE);}
	if(pos!=0){before=dispLn.substr(0,pos);}	
	string currDispLn=dispLn.substr(pos,dispLn.find("|",pos)-pos);
	after=dispLn.substr(dispLn.find("|",pos),dispLn.length()-dispLn.find("|",pos));
	// Update Current Computation at updatePos with dispStr
	string bld=currDispLn.substr(0,1+dispLabel.length());
	pos=1;
	for(int i=1+dispLabel.length();i<currDispLn.length();i++)
		{if(pos==updatePos){bld+=dispStr;i=i+dispStr.length()-1;pos=pos+dispStr.length()-1;}
		else{bld+=currDispLn[i];}
		if(pos>=8){break;}
		else{pos++;}}
	// Replace Current Computation Display with bld
	dispLn=before+bld+after;
	if(lnPos!=0){cerr<<setPos+dwnLn+clearLn+"\r"+dispLn+"\n";}
	else{cerr<<setPos+clearLn+"\r"+dispLn+"\n";}}

void updateGuiDisplayString(string dispStr,string calcNum,int updatePos,int nRows,string dispLn,string guiDisplayFile)
	{// Find the Current Computation in the Display String Line
	string dispLabel="["+calcNum+"]";
	int pos;//=dispLn.find(dispLabel,0);
	string before="",after="";
	string currDispLn;
	// Update Current Computation at updatePos with dispStr
	string bld;
	// Update GUI Display File
	int Sz=150000,Counter=0;
	char Val[Sz];
	string tmp="",output="";
	ifstream fIn; ofstream fOut;
	bool ATTEMPTING=true;
	while(ATTEMPTING)
		{fIn.open(guiDisplayFile.c_str());
		if(!fIn.fail())
			{fIn.getline(Val,Sz);
			tmp=Val;fIn.close();
			break;ATTEMPTING=false;}
		else if(Counter>MAX_FILE_OPEN_ATTEMPTS)
			{tmp="";break;ATTEMPTING=false;
			}
		else
			{Counter++;}
		}
	if(tmp.length()!=0 && tmp.find(dispLabel,0)!=string::npos)
		{// Find old Line and replace with new line
		pos=tmp.find(dispLabel,0);
		//if(pos==string::npos){cerr<<"Error in updateGuiDisplayString!\nCould not find current computation."<<endl;exit(EXIT_FAILURE);}
		if(pos!=0){before=tmp.substr(0,pos);}
		currDispLn=tmp.substr(pos,tmp.find("|",pos)-pos);
		after=tmp.substr(tmp.find("|",pos),tmp.length()-tmp.find("|",pos));
		bld=currDispLn.substr(0,1+dispLabel.length());
		pos=1;
		for(int i=1+dispLabel.length();i<currDispLn.length();i++)
			{if(pos==updatePos){bld+=dispStr;i=i+dispStr.length()-1;pos=pos+dispStr.length()-1;}
			else{bld+=currDispLn[i];}
			if(pos>=8){break;}
			else{pos++;}}
		// Replace Current Computation Display with bld
		output=before+bld+after;
		//	
		ATTEMPTING=true;Counter=0;
		while(ATTEMPTING)
			{fOut.open(guiDisplayFile.c_str(),ofstream::out|ofstream::trunc);
			if(!fOut.fail())
				{fOut<<output; fOut.close();
				break;ATTEMPTING=false;}
			else if(Counter>MAX_FILE_OPEN_ATTEMPTS)
				{break;ATTEMPTING=false;
				}
			else
				{Counter++;}
			}
		//if(fOut.fail()){cerr<<"Error in updateGuiDisplayString!\nGUI display file could not be opened.\n"<<guiDisplayFile;exit(EXIT_FAILURE);}
		}
	
	}

double calcHydrodynamicRadius(double diffusivity,double temperature,double viscosity)
	{double Rh=kb*temperature/(6*pi*viscosity*diffusivity);
	Rh*=1e10;	// Convert from meters to Angstroms
	return Rh;}

string update_display(int totalCalc,int currCalc,int threadOffset)
	{// Restore Cursor Position to top of display
	string setPos="\033[u";
	string clearLn="\033[K";
	// Move Cursor to End of Current Display
	string dwnLn="\033["+cnvrtNumToStrng(MAX_DISPLAY_ROWS,0)+"B";
	cerr<<setPos+dwnLn+"\n\r";
	int heightBuffer=5,widthBuffer=10,tWidth=0,tHeight=0;
	string Output="";
	if(currCalc>MAX_DISPLAY_ROWS){tHeight=MAX_DISPLAY_ROWS;}
	else{tHeight=currCalc;}
	// Determine Number of Columns to Display
	int numCols;
	if(currCalc>=MAX_DISPLAY_ROWS*MAX_DISPLAY_COLS){numCols=MAX_DISPLAY_COLS;}	
	else{numCols=currCalc/MAX_DISPLAY_ROWS+1;}
	int threadLabelWidth=log10(totalCalc)+1+3;
	int progressLabelWidth=8;
	int colWidth=threadLabelWidth+progressLabelWidth+2;
	tWidth=numCols*colWidth;
	//resizeTerminal(tWidth+widthBuffer,tHeight+heightBuffer);
	int df=0,val=log10(totalCalc)+1,index=0;
	string tNum="";
	for(int i=0;i<tHeight;i++)
		{for(int j=0;j<numCols;j++)
			{index=i+j*tHeight+1+threadOffset;
			if(index<=totalCalc)
				{tNum=cnvrtNumToStrng(index,0);
				df=val-tNum.length();
				if(j==numCols-1){Output+=makeWhiteSpace(df)+"["+tNum+"] ........|";}
				else{Output+=makeWhiteSpace(df)+"["+tNum+"] ........| ";}}}
		Output+="\n";}
	return Output;}

void* runMSMS(msmsInput& In,string calcNum)
	{string sphereInFile=In.pdb_msmsFldr+In.PDB+"_"+calcNum+".xyzrn";
	string triSurfOutFile=In.pdb_msmsFldr+In.PDB+"_"+calcNum+"_surface";
	string areaFile=In.pdb_msmsFldr+In.PDB+"_"+calcNum+"_area";
	string outputFile=In.pdb_surfFldr+"fromMSMS/msmsOutput_SR"+In.probeRadius+"_"+calcNum+".txt";
	// Programs must be executed in MSMS Computation Folder
	//chdir(In.msmsFldr.c_str());
	long Res;
	// Run PDB_TO_XYZRN
	string cmd="./"+In.PDB_TO_XYZRN+" "+checkFilePathParantheses(In.pdbFile)+" > "+checkFilePathParantheses(sphereInFile);
	bool ATTEMPTING=true;
	int statusCode;	
	while(ATTEMPTING)
		{Res=syscall(SYS_chdir,In.pdb_msmsFldr.c_str());
		statusCode=system(cmd.c_str());
		if(Res==0 && statusCode==0){break;ATTEMPTING=false;}}
	// Run MSMS
	cmd="./"+checkFilePathParantheses(In.MSMS)+" -if "+checkFilePathParantheses(sphereInFile)+" -of "+checkFilePathParantheses(triSurfOutFile)+" -af "+checkFilePathParantheses(areaFile)+" -probe_radius "+In.probeRadius+" -density "+In.surfPointDensity+" -hdensity "+In.surfPointHiDensity+" -no_header > "+checkFilePathParantheses(outputFile);
	//cmd="msms -if "+checkFilePathParantheses(sphereInFile)+" -of "+checkFilePathParantheses(triSurfOutFile)+" -af "+checkFilePathParantheses(areaFile)+" -probe_radius "+In.probeRadius+" -density "+In.surfPointDensity+" -hdensity "+In.surfPointHiDensity+" -no_header > "+checkFilePathParantheses(outputFile);
	//cmd="bash msms -if "+sphereInFile+" -of "+triSurfOutFile+" -af "+areaFile+" -probe_radius "+In.probeRadius+" -density "+In.surfPointDensity+" -hdensity "+In.surfPointHiDensity+" -no_header > "+outputFile;
	ATTEMPTING=true;
	int Err;
	char errMsg[256];

	while(ATTEMPTING)
		{try
			{Res=syscall(SYS_chdir,In.pdb_msmsFldr.c_str());
			statusCode=system(cmd.c_str());
			throw statusCode;
			}
		//catch(exception& Exc)
		catch(int Exc)
			{//cerr<<Exc.what()<<endl;
			//cerr<<"Exception occurred.\n"<<Exc<<endl;
			if(Exc==35584)
				{//cerr<<"HERE"<<endl;
				In.SUCCESS=false;
				return 0;
				//exit(EXIT_FAILURE);
				}
			}
		//cerr<<statusCode<<endl;
		if(Res==0 && statusCode==0)
			{In.SUCCESS=true;
			break;
			ATTEMPTING=false;}
		else if(statusCode==32512)
			{cmd="msms -if "+checkFilePathParantheses(sphereInFile)+" -of "+checkFilePathParantheses(triSurfOutFile)+" -af "+checkFilePathParantheses(areaFile)+" -probe_radius "+In.probeRadius+" -density "+In.surfPointDensity+" -hdensity "+In.surfPointHiDensity+" -no_header > "+checkFilePathParantheses(outputFile);}
		else if(statusCode==35584)
			{// Terminate Thread
			//cerr<<"Terminating Calculation Number "<<calcNum<<"...\n";
			//terminate();
			//exit(EXIT_FAILURE);
			//~thread();
			//return 0;
			}
		else if(statusCode!=0 && statusCode!=32512 && statusCode!=256)
			{// An Error has occured, capture error number with erno
			// 35584 = segmentation fault
			// 256 = file in use
			Err=errno;//erase_display();
			cerr<<"ERROR in runMSMS!\nStatus Code: "<<statusCode<<"\nError Number: "<<Err<<"\n"<<strerror_r(Err,errMsg,256)<<"\n";
			exit(EXIT_FAILURE);}}

/*
	int PID=fork(),status;
	if(PID==0)
		{// Run in Child Process, next time fork is executed it gives child process PID
		Res=syscall(SYS_chdir,In.pdb_msmsFldr.c_str());
		statusCode=system(cmd.c_str());
		exit(EXIT_SUCCESS);}
	else if(PID>0)
		{// Return to parent processs, fork should hold child process PID
		waitpid(PID,&status,0);}
	else if(PID==-1)
		{// fork() failed
		cerr<<"Fork failed.\n";exit(EXIT_FAILURE);}
	else{// Unpredicted Output
		cerr<<"Output of fork() unpredicted:\n"<<PID<<endl;}
*/

	// Move Files to Appropriate Folders
	// Face File
	string oldF=In.pdb_msmsFldr+In.PDB+"_"+calcNum+"_surface.face";
	string newF=In.pdb_surfFldr+"fromMSMS/SR"+In.probeRadius+"_"+calcNum+".face";
	ATTEMPTING=true;
	while(ATTEMPTING)
		{Res=syscall(SYS_rename,oldF.c_str(),newF.c_str());
		if(Res==0){break;ATTEMPTING=false;}}
	// Vertices File
	//pdb_surfFldr+"fromMSMS/SR0.00"+"_"+calcNum+".vert";
	oldF=In.pdb_msmsFldr+In.PDB+"_"+calcNum+"_surface.vert";
	newF=In.pdb_surfFldr+"fromMSMS/SR"+In.probeRadius+"_"+calcNum+".vert";
	ATTEMPTING=true;
	while(ATTEMPTING)
		{Res=syscall(SYS_rename,oldF.c_str(),newF.c_str());
		if(Res==0){break;ATTEMPTING=false;}}
	// Surface Area File
	oldF=In.pdb_msmsFldr+In.PDB+"_"+calcNum+"_area.area";
	newF=In.pdb_surfFldr+"fromMSMS/SR"+In.probeRadius+"_"+calcNum+".area";
	ATTEMPTING=true;
	while(ATTEMPTING)
		{Res=syscall(SYS_rename,oldF.c_str(),newF.c_str());
		if(Res==0){break;ATTEMPTING=false;}}
	}

void readMSMSVertFile(string fileName,vertData& Data)
	{int Num=countNumLinesMSMSVertFile(fileName);

	Data.x=new double[Num];
	Data.y=new double[Num];
	Data.z=new double[Num];
	Data.nx=new double[Num];
	Data.ny=new double[Num];
	Data.nz=new double[Num];
	Data.analyticSurf=new int[Num];
	Data.sphereNum=new int[Num];
	Data.vertType=new int[Num];
	Data.closestAtom=new string[Num];
	Data.numVertices=Num;

	float x,y,z,nx,ny,nz;
	int analyticSurf,sphereNum,vertType,Counter=0;
	char closestAtom[15000];
	//string msms_vertices_file_format="%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f %7d %7d %2d %s";	
	string msms_vertices_file_format="%f %f %f %f %f %f %d %d %d %s";

	FILE *fIn;
    fIn=fopen(fileName.c_str(),"r");
    if(fIn==NULL)
        {cerr<<"ERROR in readMSMSVertFile!\nInput Vertices File could not be found.\n\nFilepath:\t"<<fileName<<endl;exit(EXIT_FAILURE);}

	fseek(fIn,0,SEEK_END);
	long fileSize=ftell(fIn);
	rewind(fIn);	

	int lineSz;
	while(!feof(fIn))
        {fscanf(fIn,msms_vertices_file_format.c_str(),&x,&y,&z,&nx,&ny,&nz,&analyticSurf,&sphereNum,&vertType,&closestAtom);		
		lineSz=ftell(fIn);
		if(lineSz==fileSize){break;}
         else
            {Data.x[Counter]=x;
			Data.y[Counter]=y;
			Data.z[Counter]=z;
			Data.nx[Counter]=nx;
			Data.ny[Counter]=ny;
			Data.nz[Counter]=nz;
			Data.analyticSurf[Counter]=analyticSurf;
			Data.sphereNum[Counter]=sphereNum;
			Data.vertType[Counter]=vertType;
			Data.closestAtom[Counter]=closestAtom;
			Counter++;}}
	fclose(fIn);}

void writeHydroproInput(string calcType,string density,string molecularWeight,string proteinName,string oFile,string PDB,double specificVolume,double temperature,double viscosity)
	{string elementRadius,Output;
	ofstream fOut;
	fOut.open(oFile.c_str());
	if(fOut.fail()){cerr<<"ERROR in writeHydroproInput!!!\nFile could not be opened.\n"<<oFile<<endl;exit(EXIT_FAILURE);}
	Output=proteinName+"\t!Name of molecule\n";
	Output+=PDB+"\t\t\t\t!Name for output file\n";
	Output+=PDB+".pdb\t\t\t!Structural (PDB) file\n";
	if(calcType=="1"){elementRadius="2.9";}
	else if(calcType=="2"){elementRadius="4.8";}
	else if(calcType=="4"){elementRadius="6.1";}
	else{cerr<<"ERROR in writeHydroproInput!!!\nCalculation Type not recognized. Must be 1, 2 or 4."<<endl;exit(EXIT_FAILURE);}
	Output+=calcType+"\t\t\t\t\t!Type of calculation\n";
	Output+=elementRadius+",\t\t\t\t!AER, radius of primary elements\n";
	Output+="6,\t\t\t\t\t!NSIG\n";
	Output+="1.0,\t\t\t\t!Minimum radius of beads in the shell (SIGMIN)\n";
	Output+="2.0,\t\t\t\t!Maximum radius of beads in the shell (SIGMAX)\n";
	temperature-=273.15;	// Convert from Kelvin to Centigrade
	//Output+=cnvrtNumToStrng(temperature,2)+".,\t\t\t\t!T (temperature, centigrade)\n";
	Output+=formatNumberString(cnvrtNumToStrng(temperature,BUFFER_SIZE))+",\t\t\t\t!T (temperature, centigrade)\n";
	viscosity*=10;	// Convert from Pa s to poise	
	Output+=formatNumberString(cnvrtNumToStrng(viscosity,BUFFER_SIZE))+",\t\t!ETA (Viscosity of the solvent in poises)\n";
	Output+=molecularWeight+",\t\t\t!RM (Molecular weight)\n";
	specificVolume*=1000;	// Convert from L/g to mL/g
	Output+=formatNumberString(cnvrtNumToStrng(specificVolume,BUFFER_SIZE))+",\t\t\t!Partial specific volume, cm3/g\n";
	Output+=density+",\t\t!Solvent density, g/cm3\n";
	Output+="-1\t\t\t\t\t!Number of values of Q\n";
	Output+="-1\t\t\t\t\t!Number of intervals\n";
	Output+="0,\t\t\t\t\t!Number of trials for MC calculation of covolume\n";
	Output+="1\t\t\t\t\t!IDIF=1 (yes) for full diffusion tensors\n";
	Output+="*\t\t\t\t\t!End of file\n";
	fOut<<Output;
	fOut.close();}

void runHYDROPRO(string PDB,string pdb_hydroFldr,string inFile,string logFile,string calcNum)
	{// Define PDB File for Hydrodynamics Calculation
	string pdbFile=pdb_hydroFldr+PDB+"_"+calcNum+".pdb";
	// Generate Temporary Executable
	string tempXFile=pdb_hydroFldr+"hydropro10-lnx"+calcNum+".exe";
	char first[tempXFile.length()];
	for(int i=0;i<tempXFile.length();i++){first[i]=tempXFile[i];}
	char second[inFile.length()];
	for(int i=0;i<inFile.length();i++){second[i]=inFile[i];}
	char third[pdbFile.length()];
	for(int i=0;i<pdbFile.length();i++){third[i]=pdbFile[i];}
	char* argv[]={first,second,third,NULL};
	char* envp[]={NULL};
	bool ATTEMPTING=true;
	long Res,Res2;
	freopen(logFile.c_str(),"a",stdout);
	while(ATTEMPTING)
		{Res=syscall(SYS_chdir,pdb_hydroFldr.c_str());
		Res2=syscall(SYS_execve,tempXFile.c_str(),argv,envp);		
		if(Res==0 && Res2==0){break;ATTEMPTING=false;}}
	exit(3);}

double getDiffusivity(string hydroproFile)
	{ifstream fIn;
	fIn.open(hydroproFile.c_str());
	if(fIn.fail())
		{cerr<<"ERROR in getDiffusivity\nHYDROPRO output file could not be opened.\n"<<hydroproFile<<endl;
		 exit(EXIT_FAILURE);}
	// Copy Input File
	int Counter=0,Sz=15000;                      // overshoot actual line size, C++ automatically fixes c-string lengths
	char Val[Sz];
	string tmp,tmp1;
	// Skip HYDROPRO Header (7 lines)
	// Skip HYDROPRO data until translational diffusivity (cm^2/s)
	for(int i=0;i<25;i++){fIn.getline(Val,Sz);}
	tmp=Val;
	tmp1=tmp.substr(tmp.length()-16,9);
	double diffusivity=strtod(tmp1.c_str(),NULL);
	diffusivity/=10000;	// Convert from cm^2/s to m^2/s
	fIn.close();
	return diffusivity;}

double getDiffusivity_hullrad(string hullradFile)
	{ifstream fIn;
	fIn.open(hullradFile.c_str());
	if(fIn.fail())
		{cerr<<"ERROR in getDiffusivity_hullrad\nHullrad output file could not be opened.\n"<<hullradFile<<endl;
		 exit(EXIT_FAILURE);}
	// Copy Input File
	int Counter=0,Sz=15000;                      // overshoot actual line size, C++ automatically fixes c-string lengths
	char Val[Sz];
	string tmp,tmp1;
	// Skip 8 lines
	for(int i=0;i<9;i++){fIn.getline(Val,Sz);}
	// Acquire translational diffusivity (cm^2/s)
	tmp=Val;
	int pos=tmp.find(":",0);
	tmp1=tmp.substr(pos+1,tmp.find("cm^2/s",0)-pos-1);
	double diffusivity=strtod(tmp1.c_str(),NULL);
	diffusivity/=10000;	// Convert from cm^2/s to m^2/s
	fIn.close();
	return diffusivity;}

void runPROPKA_PDB2PQR(string pdbFile,string pqrFile,string calcNum,string forceField,string pdb_pqrFldr,string propkaOutFile,string pH,string PDB2PQR,string logFile)
	{//string cmd="python3.8 "+checkFilePathParantheses(PDB2PQR)+" --ff="+forceField+" --with-ph="+pH+" --whitespace -v "+checkFilePathParantheses(pdbFile)+" "+checkFilePathParantheses(pqrFile);
	string cmd="pdb2pqr --ff="+forceField+" --titration-state-method=propka --with-ph="+pH+" --whitespace "+checkFilePathParantheses(pdbFile)+" "+checkFilePathParantheses(pqrFile);
	bool ATTEMPTING=true;
	int statusCode;
	freopen(logFile.c_str(),"a",stdout);
	//freopen(logFile.c_str(),"a",stderr);
	while(ATTEMPTING)
		{statusCode=system(cmd.c_str());
		if(statusCode==0){break;ATTEMPTING=false;}}
	// Move PROPKA Output File to Appropriate Folder
	//string oldF=pdb_pqrFldr+"pH"+pH+"_"+calcNum+".propka";
	string oldF=pdb_pqrFldr+"pH"+pH+"_"+calcNum+".log";
	//long Res=syscall(SYS_rename,logFile.c_str(),oldF.c_str());
	long Res=syscall(SYS_rename,oldF.c_str(),propkaOutFile.c_str());	
	//
	//fflush(stderr);
	}

/*void runAPBS(string exeFile,string inFile,string logFile)
	{// Define Input Arguments
	char first[exeFile.length()+1];
	for(int i=0;i<exeFile.length();i++){first[i]=exeFile[i];}
	first[exeFile.length()]='\0';
	char second[inFile.length()+1];
	for(int i=0;i<inFile.length();i++){second[i]=inFile[i];} 
	int df=0,charLen=strlen(second),strLen=inFile.length();
	if(charLen>strLen)
		{df=charLen-strLen;
		second[charLen-df]=0;}
	char* argv[]={first,second,NULL};
	char* envp[]={NULL};
	long Res;
	bool ATTEMPTING=true;
	freopen(logFile.c_str(),"a",stdout);
	while(ATTEMPTING)
		{Res=syscall(SYS_execve,exeFile.c_str(),argv,envp);
		if(Res==0){break;ATTEMPTING=false;}}
	exit(EXIT_SUCCESS);}
*/
void runAPBS(string exeFile,string inFile,string logFile)
	{/*char** argv=new char*[3];
	argv[0]=new char[exeFile.length()+1];
	for(int i=0;i<exeFile.length();i++){argv[0][i]=exeFile[i];}
	argv[0][exeFile.length()]='\0';
	argv[1]=new char[inFile.length()+1];
	for(int i=0;i<inFile.length();i++){argv[1][i]=inFile[i];}
	argv[1][inFile.length()]='\0';
	argv[2]={NULL};
	char* envp[]={NULL};
	*/
	freopen(logFile.c_str(),"a",stdout);
	//freopen(logFile.c_str(),"a",stderr);
	//long Res=syscall(SYS_execve,exeFile.c_str(),argv,envp);

	string inFldr=get_containing_folder(inFile);
	string inFileName=get_file_name(inFile);
	string cmd="apbs "+inFileName;
	int PID=fork(),status,res;
	long Result;
	bool ATTEMPTING;
	if(PID==0)
		{ATTEMPTING=true;
		while(ATTEMPTING)
			{Result=syscall(SYS_chdir,inFldr.c_str());
			res=system(cmd.c_str());

			if(Result==0 && res==0){ATTEMPTING=false; break;}
			}
		}
	else if(PID>0)
		{waitpid(PID,&status,0);}
	exit(EXIT_SUCCESS);
	}

int countNumLines(string file)
    {ifstream fIn;
     fIn.open(file.c_str());
     int Sz=200;
     char Val[Sz];
     if(fIn.fail()){cerr<<"ERROR in countNumLines!\nPDB input file could not be opened.\n"<<file<<endl;exit(EXIT_FAILURE);}
     fIn.getline(Val,Sz);
     string chainLabel;
     string bld;
     int Count=0,bCount=0,cCount=0;
     int col1=6,col2=5+col1,col3=4+col2,col4=5+col3,col5=2+col4,col6=4+col5,col7=12+col6,col8=8+col7,col9=8+col8,col10=6+col9,col11=6+col10,col12=12+col11,col13=2+col12;
     while(!fIn.eof())
        {bld="";
         for(int i=0;i<col1;i++)
            {if(isalpha(Val[i]))
                {bld=bld+Val[i];}}
         if(bld=="ATOM")
            {Count++;}
         else if(bld=="HETATM")
            {Count++;}
         else if(bld=="TER")
            {Count++;}
         fIn.getline(Val,Sz);}
    // Account for counting of end
    //Count--;
    fIn.close();
    return Count;}

void definePDBData(pdbData& Data,int Sz)
    {string* recordType=new string[Sz];
     int* atomNum=new int[Sz];
     string* atomLabel=new string[Sz];
     string* residue=new string[Sz];
     string* chainLabel=new string[Sz];
     int* atmChnNum=new int[Sz];
     double* x=new double[Sz];
     double* y=new double[Sz];
     double* z=new double[Sz];
     double* atmOccupancy=new double[Sz];
     double* bFactor=new double[Sz];
     string* terminus=new string[Sz];

     Data.recordType=recordType;
     Data.atomNum=atomNum;
     Data.atomLabel=atomLabel;
     Data.residue=residue;
     Data.chainLabel=chainLabel;
     Data.atmChainNum=atmChnNum;
     Data.x=x;
     Data.y=y;
     Data.z=z;
     Data.atmOccupancy=atmOccupancy;
     Data.bFactor=bFactor;
     Data.terminus=terminus;}

void readPDBFile(string file,pdbData& Data)
    {ifstream fIn;
     fIn.open(file.c_str());
     if(fIn.fail())
        {cerr<<"ERROR in readPDBFile of ioProt!!!\n PDB file could not be opened.\n"<<file<<endl;
         exit(EXIT_FAILURE);}
     // Copy Input File
    int Counter=0,Sz=15000;                      // overshoot actual line size, C++ automatically fixes c-string lengths
    char Val[Sz];
    fIn.getline(Val,Sz);
    string recordType,atomNumber,atomLabel,residue,chainLabel,atmChainNumber,xPos,yPos,zPos,atomOccupancy,bFactor,terminus;
    int col1=6,col2=5+col1,col3=5+col2,col4=5+col3,col5=2+col4,col6=4+col5,col7=12+col6,col8=8+col7,col9=8+col8,col10=6+col9,col11=6+col10,col12=12+col11,col13=2+col12;
    while(!fIn.eof())
        {recordType="",atomNumber="",atomLabel="",residue="",chainLabel="",atmChainNumber="",xPos="",yPos="",zPos="",atomOccupancy="",bFactor="",terminus="";
         // Read First Column
         for(int i=0;i<col1;i++)
            {if(isalpha(Val[i]))
                {recordType=recordType+Val[i];}}
         if(recordType=="ATOM"||recordType=="HETATM"||recordType=="TER")
		 //if(recordType=="ATOM"||recordType=="HETATM")
            {// Read Second Column
             for(int i=col1;i<col2;i++)
                {if(Val[i]!=' ')
                    {atomNumber=atomNumber+Val[i];}}
             // Read Third Column
             for(int i=col2;i<col3;i++)
                {if(Val[i]!=' ')
                    {atomLabel=atomLabel+Val[i];}}
             // Read Fourth Column
             for(int i=col3;i<col4;i++)
                {if(Val[i]!=' ')
                    {residue=residue+Val[i];}}
             // Read Fifth Column
             for(int i=col4;i<col5;i++)
                {if(Val[i]!=' ')
                    {chainLabel=chainLabel+Val[i];}}
             // Read Sixth Column
             for(int i=col5;i<col6;i++)
                {if(Val[i]!=' ')
                    {atmChainNumber=atmChainNumber+Val[i];}}
             // Read Seventh Column
             for(int i=col6;i<col7;i++)
                {if(Val[i]!=' ')
                    {xPos=xPos+Val[i];}}
             // Read Eighth Column
             for(int i=col7;i<col8;i++)
                {if(Val[i]!=' ')
                    {yPos=yPos+Val[i];}}
             // Read Ninth Column
             for(int i=col8;i<col9;i++)
                {if(Val[i]!=' ')
                    {zPos=zPos+Val[i];}}
             // Read Tenth Column
             for(int i=col9;i<col10;i++)
                {if(Val[i]!=' ')
                    {atomOccupancy=atomOccupancy+Val[i];}}
             // Read Eleventh Column
             for(int i=col10;i<col11;i++)
                {if(Val[i]!=' ')
                    {bFactor=bFactor+Val[i];}}
             // Read Twelveth Column
             for(int i=col11;i<col12;i++)
                {if(Val[i]!=' ')
                    {terminus=terminus+Val[i];}}
             Data.recordType[Counter]=recordType;
             Data.atomNum[Counter]=atoi(atomNumber.c_str());
             Data.atomLabel[Counter]=atomLabel;
             Data.residue[Counter]=residue;
             Data.chainLabel[Counter]=chainLabel;
             Data.atmChainNum[Counter]=atoi(atmChainNumber.c_str());
             Data.x[Counter]=atof(xPos.c_str());
             Data.y[Counter]=atof(yPos.c_str());
             Data.z[Counter]=atof(zPos.c_str());
             Data.atmOccupancy[Counter]=atof(atomOccupancy.c_str());
             Data.bFactor[Counter]=atof(bFactor.c_str());
             Data.terminus[Counter]=terminus;
             Counter++;}
         fIn.getline(Val,Sz);}
         fIn.close();}

double getSpecificVolume(string msmsOutputFile,double molecularWeight)
	{ifstream fIn;
	fIn.open(msmsOutputFile.c_str());
     if(fIn.fail())
        {cerr<<"ERROR in read_hydropro_output_diffusivity!!!\nHYDROPRO output file could not be opened.\n"<<msmsOutputFile<<endl;
         exit(EXIT_FAILURE);}

	// Copy Input File
	int Counter=1,Sz=15000;                      // overshoot actual line size, C++ automatically fixes c-string lengths
	char Val[Sz];
	string tmp,tmp1;
	// Find Line Containing Volume
	string term="    Comp. probe_radius SES_volume SES_area)";
	bool RUNNING=true;
	while(RUNNING)
		{fIn.getline(Val,Sz);
		tmp=Val;
		if(tmp==term){RUNNING=false;}}
	// Capture line containing volume
	fIn.getline(Val,Sz);
	tmp=Val;
	// Find Space between volume and surfaec area
	RUNNING=true;
	while(RUNNING)
		{tmp1=tmp.substr(tmp.length()-Counter,1);
		if(tmp1==" "){RUNNING=false;}
		else{Counter++;}}
	// Find volume value
	RUNNING=true;
	while(RUNNING)
		{tmp1=tmp.substr(tmp.length()-Counter,1);
		if(tmp1!=" ")
			{RUNNING=false;}
		else
			{Counter++;}}
	// Record Volume
	string bld="";
	RUNNING=true;
	while(RUNNING)
		{tmp1=tmp.substr(tmp.length()-Counter,1);
		bld+=tmp1;
		if(tmp1==" ")
			{RUNNING=false;}
		else
			{Counter++;}}
	tmp="";
	for(int i=bld.length()-2;i>-1;i--)
		{tmp+=bld.substr(i,1);}	
	fIn.close();
	double Volume=strtod(tmp.c_str(),NULL);
	Volume/=1e27;	// Convert from A^3 to L
	double Mass=molecularWeight/Na;	// Determine mass in grams
	double specificVolume=Volume/Mass;	
	return specificVolume;}

Vec calc_center_of_mass(pdbData Data,int Sz)
    {Vec value;
	double xVal=0,yVal=0,zVal=0;
    for(int i=0;i<Sz;i++)
        {xVal+=Data.x[i];
         yVal+=Data.y[i];
         zVal+=Data.z[i];}
    xVal/=Sz;
    yVal/=Sz;
    zVal/=Sz;
    value.x=xVal;
    value.y=yVal;
    value.z=zVal;
	return value;}

double calc_distance(Vec Vec1,Vec Vec2)
    {return sqrt((Vec2.x-Vec1.x)*(Vec2.x-Vec1.x)+(Vec2.y-Vec1.y)*(Vec2.y-Vec1.y)+(Vec2.z-Vec1.z)*(Vec2.z-Vec1.z));}

double calc_protein_radius(vertData Data,int Sz,Vec centerCoord)
    {Vec pos;
	double Radius=0;
	for(int i=0;i<Sz;i++)
		{pos.x=Data.x[i];
		pos.y=Data.y[i];
		pos.z=Data.z[i];
		Radius+=calc_distance(pos,centerCoord);}
	Radius/=Sz;
	return Radius;}

double calc_protein_radius_std(vertData Data,int Sz,Vec centerCoord,double Radius)
     {double radStd=0,dst=0;
	Vec pos;
	for(int i=0;i<Sz;i++)
		{pos.x=Data.x[i];
		pos.y=Data.y[i];
		pos.z=Data.z[i];
		dst=calc_distance(pos,centerCoord);
		radStd+=(dst-Radius)*(dst-Radius);}
	return sqrt(radStd/(Sz-1));}

double calcProteinMolecularWeight(string pdbFile)
	{double Output=0;
	pdbData Data;
	int numLines=countNumLines(pdbFile);
	definePDBData(Data,numLines);
	readPDBFile(pdbFile,Data);
	string numAtomList="",atomList="",atomNm="",rpt="",branch="",delimiter=";",tmp,atm;
	char letter,second_letter;
	bool PASS,SEARCHING;
	string *Atom;
	int *numAtom;
	int pos,aCount=0,multiplier,Counter,numRepeats,i;

	for(int a=0;a<numLines;a++)
		{atm=Data.atomLabel[a];		
		// Get Atom Name
		PASS=false;atomNm="";multiplier=1;i=0;
		letter=atm[i];
		if(isupper(letter))
			{atomNm=letter;
			// Check if second letter completes atom name or gives number of atoms
			if(i+1<=atm.length()-1)
				{second_letter=atm[i+1];
				if(islower(second_letter)){atomNm+=second_letter;}}
			// Check if number behind Atom Name
			if(i+atomNm.length()<=atm.length()-1)
				{SEARCHING=true;Counter=0;tmp="";
				while(SEARCHING)
					{letter=atm[i+atomNm.length()+Counter];
					if(isdigit(letter)){tmp+=letter;Counter++;}
					else{SEARCHING=false;break;}}
				if(tmp.compare("")==0){multiplier=1;}
				else{multiplier=atoi(tmp.c_str());}
				i+=Counter;}			
			// Check is atom is already in Atom List, if not add
			if(atomList.compare("")==0)
				{// Initialize Atoms Present Lists
				atomList=atomNm+delimiter;
				numAtomList=cnvrtNumToStrng(multiplier,0)+delimiter;
				aCount++;
				Atom=fill_string_array(atomList,aCount,delimiter);
				numAtom=fill_int_array(numAtomList,aCount,delimiter);}
			else{// Check is atom is already in Atom List, if not add
				PASS=false;
				for(pos=0;pos<aCount;pos++){tmp=Atom[pos];if(tmp.compare(atomNm)==0){PASS=true;break;}}
				if(!PASS)
					{// Atom not found in list, so add
					atomList+=atomNm+delimiter;
					numAtomList+=cnvrtNumToStrng(multiplier,0)+delimiter;
					delete [] Atom;delete [] numAtom;
					aCount++;
					Atom=fill_string_array(atomList,aCount,delimiter);
					numAtom=fill_int_array(numAtomList,aCount,delimiter);}
				else{// Atom Found, so update numAtomList
					numAtom[pos]=numAtom[pos]+multiplier;
					numAtomList=array2String(numAtom,aCount,delimiter);}}
			}
		}
	Solution S;
	for(int i=0;i<aCount;i++)
		{// Identify Atom
		tmp=Atom[i];
		for(pos=0;pos<S.numElements;pos++){if(tmp.compare(S.atomicSymbol[pos])==0){break;}}
		// Get Atomic Molecular Weight
		Output+=S.atomicMolecularWeight[pos]*numAtom[i];}
	return Output;}

void inFileEditer(string* writeTypes,int N,string pdie,strVec dim,string pqrFile,string apbsInputgen,string ionData,string temperature,string largestIon,string er,string fldrExt)
	{// APBS IN File Parameters
	//int NumProcessors=12;
	//int NumMem=11700;
	string memMax="1024"; //11700/12
	string cfac="5";
	string fadd="80";	//80
	string Frmt="dx";
	string method="auto";
	// psize.py inputs
	string flags="--cfac="+cfac;
	flags+=" --fadd="+fadd;
	flags+=" --gmemceil="+memMax;
	flags+=" --method="+method;
	// APBS inputs
	flags+=" --pdie="+pdie;
	flags+=" --Frmt="+Frmt;
	flags+=" --istrng=1";
	flags+=" --dime-x="+dim.x;
	flags+=" --dime-y="+dim.y;
	flags+=" --dime-z="+dim.z;
	flags+=" --unfrmGrdSpcng";
	// APBS Flag write types
	string wt="";
	for(int i=0;i<N;i++)
		{wt=writeTypes[i];
		if(wt.compare("atompot")==0){flags+=" --atompot";}
		else if(wt.compare("charge")==0){flags+=" --charge";}
		else if(wt.compare("dielx")==0){flags+=" --dielx";}
		else if(wt.compare("diely")==0){flags+=" --diely";}
		else if(wt.compare("dielz")==0){flags+=" --dielz";}
		else if(wt.compare("edens")==0){flags+=" --edens";}		
		else if(wt.compare("ivdw")==0){flags+=" --ivdw";}
		else if(wt.compare("lap")==0){flags+=" --lap";}
		else if(wt.compare("ndens")==0){flags+=" --ndens";}
		else if(wt.compare("kappa")==0){flags+=" --kappa";}
		else if(wt.compare("pot")==0){flags+=" --pot";}
		else if(wt.compare("qdens")==0){flags+=" --qdens";}
		else if(wt.compare("smol")==0){flags+=" --smol";}
		else if(wt.compare("sspl")==0){flags+=" --sspl";}
		else if(wt.compare("vdw")==0){flags+=" --vdw";}
		else{cerr<<"ERROR in inFileEditer!\nUnrecognized APBS Write Type ("<<wt<<")"<<endl;}}
	// Non-Linear Poisson Boltzmann Calculation Flag
	flags+=" --npbe";
	// Ion Input Flag
	flags+=" --ion="+ionData;
	flags+=" --temp="+temperature;
	flags+=" --srad="+largestIon;
	flags+=" --sdie="+er;
	//string cmd="python3.8 "+checkFilePathParantheses(apbsInputgen)+" "+flags+" "+checkFilePathParantheses(pqrFile)+" "+checkFilePathParantheses(fldrExt);
	string cmd=PYTHON_EXE;
	cmd+=" "+checkFilePathParantheses(apbsInputgen)+" "+flags+" "+checkFilePathParantheses(pqrFile)+" "+checkFilePathParantheses(fldrExt);
	int statusCode;
	bool ATTEMPTING=true;
	while(ATTEMPTING)
		{statusCode=system(cmd.c_str());
		if(statusCode==0){break;ATTEMPTING=false;}}}

void runMULTIVALUE(string exeFile,string csvFile,string potFile,string mvFile,string logFile)
	{bool ATTEMPTING=true;
	int statusCode,pN=0;
	string potFldr=get_containing_folder(potFile),tmp,tmp2;	
	string *potFiles=get_all_files_in_dir(potFldr,pN);
	int *score=new int[pN],index,maxScore=0;
	for(int i=0;i<pN;i++){score[i]=0;}
	//cout<<"|"<<get_containing_folder(theHullRad)<<"|"<<get_file_name(theHullRad)<<"|"<<endl;
	string pHVal,tVal,sRVal,pHVal2,tVal2,sRVal2;
	determine_dx_file_name_components(get_file_name(potFile),pHVal,tVal,sRVal);
	long Res;
	if(pN==1)
		{// Only One File
		Res=syscall(SYS_rename,potFiles[0].c_str(),potFile.c_str());}
	else
		{for(int j=0;j<pN;j++)
			{determine_dx_file_name_components(get_file_name(potFiles[j]),pHVal2,tVal2,sRVal2);
			if(pHVal.compare(pHVal2)==0)
				{// Rename
				if(tVal.compare(tVal2)==0){Res=syscall(SYS_rename,potFiles[j].c_str(),potFile.c_str());}
				else
					{if(tVal2.length()<tVal.length())
						{for(int k=0;k<tVal2.length();k++)
							{tmp=tVal[k];
							tmp2=tVal2[k];
							if(tmp.compare(tmp2)==0){score[j]++;}
							}
						}
					}
				}
			}
		for(int j=0;j<pN;j++){if(score[j]>=maxScore){index=j; maxScore=score[j];}}
		Res=syscall(SYS_rename,potFiles[index].c_str(),potFile.c_str());
		}	

	string cmd=checkFilePathParantheses(exeFile)+" "+checkFilePathParantheses(csvFile)+" "+checkFilePathParantheses(potFile)+" "+checkFilePathParantheses(mvFile);
	freopen(logFile.c_str(),"a",stdout);
	while(ATTEMPTING)
		{statusCode=system(cmd.c_str());
		if(statusCode==0){fclose(stdout);break;ATTEMPTING=false;}}
	exit(EXIT_SUCCESS);}

void inflateVert(double inflationDistance,string inVertFile,string outVertFile,bool VERBOSE)
	{// Read MSMS Vertices File
	if(VERBOSE){fprintf(stdout,"Reading MSMS vertices file (%s)...",inVertFile.c_str());}
	vertData V;
	readMSMSVertFile(inVertFile,V);
	if(VERBOSE){fprintf(stdout,"%i vertices\n",V.numVertices);}	

	if(VERBOSE){fprintf(stdout,"Inflating Vertices by %f Angstroms...",inflationDistance);}
	for(int i=0;i<V.numVertices;i++)
		{V.x[i]+=inflationDistance*V.nx[i];
		V.y[i]+=inflationDistance*V.ny[i];
		V.z[i]+=inflationDistance*V.nz[i];}
	if(VERBOSE){fprintf(stdout,"Done.\n");}

	if(VERBOSE){fprintf(stdout,"Writing MSMS vertices file (%s)...",outVertFile.c_str());}
	writeMSMSVertFile(outVertFile,V);
	if(VERBOSE){fprintf(stdout,"Done.\n");}}

void vert2csv(string inVertFile,string outCsvFile,bool VERBOSE)
	{// Read MSMS Vertices File
	if(VERBOSE){fprintf(stdout,"Reading MSMS vertices file (%s)...",inVertFile.c_str());}
	vertData V;
	readMSMSVertFile(inVertFile,V);
	if(VERBOSE){fprintf(stdout,"%i vertices\n",V.numVertices);}	

	if(VERBOSE){fprintf(stdout,"Writing Comma Separarted Vector (CSV) file (%s)...",outCsvFile.c_str());}
	char Frmt[]="%f%s%f%s%f\n",delimiter[]=",";
	string Output="";
	//FILE *Out;
	//Out=fopen(outCsvFile.c_str(),"w");
	//for(int i=0;i<V.numVertices;i++){fprintf(Out,Frmt,V.x[i],delimiter,V.y[i],delimiter,V.z[i]);}
	//fclose(Out);
	for(int i=0;i<V.numVertices;i++)
		{Output+=cnvrtNumToStrng(V.x[i],BUFFER_SIZE)+delimiter;
		Output+=cnvrtNumToStrng(V.y[i],BUFFER_SIZE)+delimiter;
		Output+=cnvrtNumToStrng(V.z[i],BUFFER_SIZE)+"\n";}
	ofstream fOut;
	fOut.open(outCsvFile.c_str());
	if(fOut.fail()){cerr<<"ERROR in vert2csv!!!\nCSV output file could not be opened.\n"<<outCsvFile<<endl;exit(EXIT_FAILURE);}
	fOut<<Output;
	fOut.close();
	if(VERBOSE){fprintf(stdout,"Done.\n");}}

void writeMSMSVertFile(string fileName,vertData Data)
	{string spc1,spc2,spc3,spc4,spc5,spc6,spc7,spc8,spc9,spc10;
	string x,y,z,nx,ny,nz,analyticSurf,sphereNum,vertType,closestAtom;
	int frontSz,spcSz,col1=5,col2=6,col3=6,col4=6,col5=6,col6=6,col7=8,col8=8,col9=3;
	string Output="";
	for(int i=0;i<Data.numVertices;i++)	
		{// X Coordinate
		x=formatNumberString(cnvrtNumToStrng(Data.x[i],BUFFER_SIZE)); frontSz=x.find("."); spcSz=col1-frontSz; spc1=makeWhiteSpace(spcSz);
		// Y Coordinate
		y=formatNumberString(cnvrtNumToStrng(Data.y[i],BUFFER_SIZE)); frontSz=y.find("."); spcSz=col2-frontSz; spc2=makeWhiteSpace(spcSz);
		// Z Coordinate
		z=formatNumberString(cnvrtNumToStrng(Data.z[i],BUFFER_SIZE)); frontSz=z.find("."); spcSz=col3-frontSz; spc3=makeWhiteSpace(spcSz);
		// X Component of Normal Unit Vector
		nx=formatNumberString(cnvrtNumToStrng(Data.nx[i],BUFFER_SIZE)); frontSz=nx.find("."); spcSz=col4-frontSz; spc4=makeWhiteSpace(spcSz);
		// Y Component of Normal Unit Vector
		ny=formatNumberString(cnvrtNumToStrng(Data.ny[i],BUFFER_SIZE)); frontSz=ny.find("."); spcSz=col5-frontSz; spc5=makeWhiteSpace(spcSz);
		// Z Component of Normal Unit Vector
		nz=formatNumberString(cnvrtNumToStrng(Data.nz[i],BUFFER_SIZE)); frontSz=nz.find("."); spcSz=col6-frontSz; spc6=makeWhiteSpace(spcSz);
		// Analytic Surface Number
		analyticSurf=cnvrtNumToStrng(Data.analyticSurf[i],0); spcSz=col7-analyticSurf.length(); spc7=makeWhiteSpace(spcSz);
		// Sphere Number
		sphereNum=cnvrtNumToStrng(Data.sphereNum[i],0); spcSz=col8-sphereNum.length(); spc8=makeWhiteSpace(spcSz);
		// Vertice Type		
		vertType=cnvrtNumToStrng(Data.vertType[i],0); spcSz=col9-vertType.length(); spc9=makeWhiteSpace(spcSz);
		// Closest Atom
		closestAtom=Data.closestAtom[i]; 
		Output+=spc1+x+spc2+y+spc3+z;
		Output+=spc4+nx+spc5+ny+spc6+nz;
		Output+=spc7+analyticSurf+spc8+sphereNum;
		Output+=spc9+vertType+" "+closestAtom+"\n";}
	ofstream fOut;
	fOut.open(fileName.c_str());
	if(fOut.fail()){cerr<<"ERROR in writeMSMSVertFile\nOutput MSMS .vert file could not be opened.\n"<<fileName<<endl;exit(EXIT_FAILURE);}
	fOut<<Output;
	fOut.close();}

zOutput read_zpred_output_file(string zOutFile)
	{// Define and Initialize Values
	zOutput D;
	D.proteinName="";
	D.molecularWeight=-1;
	D.proteinRadius=-1;
	D.solvatedRadius=-1;
	D.zetaPotential=NAN;
	D.charge=NAN;
	D.diffusivity=NAN;
	D.semMobility=NAN;
	D.Xsp="";
	D.solution.pH="";
	D.solution.temperature="";
	D.solution.dielectric="";
	D.solution.viscosity="";
	D.solution.density="";
	D.solution.numSolute=-1;
	D.solution.numSolvent=-1;
	D.solution.debyeLength="";

	ifstream fIn;
	fIn.open(zOutFile.c_str());
	if(fIn.fail()){cerr<<"ERROR in read_zpred_output_file!\nZPRED output file could not be opened.\n"<<zOutFile<<endl;exit(EXIT_FAILURE);}
	int Sz=15000,pos;
	char Val[Sz];
	bool CHECKING,SPHERE=false,CYLINDER=false;
	string bld="",tmp="",inputName="",inputVal="",delimiter=";";
	string *sArr;
	double *dArr;
	fIn.getline(Val,Sz);
	while(!fIn.eof())
		{tmp=Val;
		pos=tmp.find(".",0);
		if(pos!=string::npos && pos==0)
			{// Line contains a period, now check if line defines an input parameter
			pos=tmp.find(":",pos);	// defines end of input parameter name
			if(pos!=string::npos)
				{// Extract Input Parameter Name
				inputName=tmp.substr(tmp.find(".",0)+1,pos-tmp.find(".",0)-1);
				// Extract Input Parameter Value(s)
				inputVal=tmp.substr(pos+1,tmp.length()-pos-1);
				// Check for Preceding WhiteSpace(s)
				CHECKING=true;
				pos=inputVal.find(" ",0);
				while(CHECKING)
					{if(pos!=string::npos){inputVal=inputVal.substr(1,inputVal.length()-1);pos=inputVal.find(" ",pos+1);}
					else{CHECKING=false;break;}}
				if(inputVal.length()!=0)
					{// Identify Parameter and Load Values
					if(inputName.compare("shape")==0)
						{// Shape Descriptor
						if(inputVal.compare("sphere")==0){SPHERE=true;}
						else if(inputVal.compare("cylinder")==0){CYLINDER=true;}
						}
					else if(inputName.compare("proteinRadius")==0)
						{// Anhydrous Protein Radius (A)
						D.proteinRadius=strtod(inputVal.c_str(),NULL);}
					else if(inputName.compare("proteinName")==0)
						{// Name of Specific PDB File under Assessment
						D.proteinName=inputVal;}
					else if(inputName.compare("solCondNum")==0)
						{// Name of Specific PDB File under Assessment
						D.solCondNum=inputVal;}
					else if(inputName.compare("molecularWeight")==0)
						{// Molecular Weight
						D.molecularWeight=strtod(inputVal.c_str(),NULL);}
					else if(inputName.compare("solvatedRadius")==0)
						{// Solvated Protein Radius (A)
						D.solvatedRadius=strtod(inputVal.c_str(),NULL);}
					else if(inputName.compare("charge")==0)
						{// Protein Net Charge Valence [e]
						D.charge=strtod(inputVal.c_str(),NULL);}
					else if(inputName.compare("diffusivity")==0)
						{// Protein Diffusivity [m^2/s]
						D.diffusivity=strtod(inputVal.c_str(),NULL);}
					else if(inputName.compare("henryMobility")==0)
						{// Ideal Single-Particle Electrophoretic Mobility [umcm/Vs]
						if(inputVal.compare("NaN")==0)
							{D.henryMobility=nan("");}
						else
							{D.henryMobility=strtod(inputVal.c_str(),NULL);}
						}
					else if(inputName.compare("kuwabaraMobility")==0)
						{// Kuwabara Concentration-Corrected Electrophoretic Mobility [umcm/Vs]
						if(inputVal.compare("NaN")==0)
							{D.kuwabaraMobility=nan("");}
						else
							{D.kuwabaraMobility=strtod(inputVal.c_str(),NULL);}
						}
					else if(inputName.compare("semMobility")==0)
						{// Ideal Single-Particle Electrophoretic Mobility [umcm/Vs]
						if(inputVal.compare("NaN")==0)
							{D.semMobility=nan("");}
						else
							{D.semMobility=strtod(inputVal.c_str(),NULL);}
						}
					//else if(SPHERE && inputName.compare("semMobility2")==0)
					//	{// Kuwabara Concentration-Corrected Electrophoretic Mobility [umcm/Vs]
					//	D.mobility=strtod(inputVal.c_str(),NULL);}
					else if(inputName.compare("Xsp")==0)
						{// Hydration Layer Thickness (A)
						D.Xsp=inputVal;}
					else if(inputName.compare("zetaPotential")==0)
						{// ZPRED Computed Zeta Potential (V)
						D.zetaPotential=strtod(inputVal.c_str(),NULL);}
					else if(inputName.compare("debyeLength")==0)
						{// Solution Debye Length (A)
						D.solution.debyeLength=inputVal;}
					else if(inputName.compare("pH")==0)
						{// Solution pH
						D.solution.pH=inputVal;}
					else if(inputName.compare("temperature")==0)
						{// Solution Temperature [K]
						D.solution.temperature=inputVal;}
					else if(inputName.compare("dielectric")==0)
						{// Solution Relative Dielectric
						D.solution.dielectric=inputVal;}
					else if(inputName.compare("viscosity")==0)
						{// Solution Viscosity [Pa s]
						D.solution.viscosity=inputVal;}
					else if(inputName.compare("density")==0)
						{// Solution Density [kg/L]
						D.solution.density=inputVal;}
					else if(inputName.compare("solvent")==0)
						{// Solution Components (Solvent Names)
						D.solution.numSolvent=count_delimiter(inputVal,delimiter);
						D.solution.solvent=fill_string_array(inputVal,D.solution.numSolvent,delimiter);}
					else if(inputName.compare("solventConc")==0)
						{// Solution Components (Solvent Concentrations) [Mole Fraction]
						D.solution.numSolvent=count_delimiter(inputVal,delimiter);
						D.solution.solventConc=fill_string_array(inputVal,D.solution.numSolvent,delimiter);}
					else if(inputName.compare("solute")==0)
						{// Solution Components (Solute Names)
						D.solution.numSolute=count_delimiter(inputVal,delimiter);
						D.solution.solute=fill_string_array(inputVal,D.solution.numSolute,delimiter);}
					else if(inputName.compare("soluteConc")==0)
						{// Solution Components (Solute Concentrations) [Molar]
						D.solution.numSolute=count_delimiter(inputVal,delimiter);
						D.solution.soluteConc=fill_string_array(inputVal,D.solution.numSolute,delimiter);}
					else if(inputName.compare("generateProfile")==0)
						{if(inputVal.compare("false")==0){D.EP.numPnts=0;}
						else
							{D.EP.numPnts=count_delimiter(inputVal,";");
							D.EP.position=new double[D.EP.numPnts];
							D.EP.potential=new double[D.EP.numPnts];
							sArr=fill_string_array(inputVal,D.EP.numPnts,";");
							for(int i=0;i<D.EP.numPnts;i++)
								{tmp=sArr[i];
								dArr=fill_double_array(tmp,2,",");
								D.EP.position[i]=dArr[0];
								D.EP.potential[i]=dArr[1];
								delete [] dArr;}
							delete [] sArr;}
						}
					else{//cerr<<"ERROR in read_zpred_output_file!\nUnrecognized Parameter ("<<inputName<<").\n";exit(EXIT_FAILURE);
						}
					}
				}
			}
		fIn.getline(Val,Sz);}
	fIn.close();
	return D;}

float readPQRFile(string pqrFile)
    {ifstream fIn;
    fIn.open(pqrFile.c_str());
    if(fIn.fail())
        {cerr<<"ERROR:\n PQR file could not be opened.\n"<<pqrFile<<endl;
         exit(EXIT_FAILURE);}

     // Copy Input File
    int Sz=15000;                      // overshoot actual line size, C++ automatically fixes c-string lengths
    char Val[Sz];
    string valBld,recordType,recordNum,tmp;
    int col1=6,col2=4+col1;
    fIn.getline(Val,Sz);
    while(!fIn.eof())
        {recordType="";
         recordNum="";
         if(!fIn.eof())
            {for(int i=0;i<col1;i++)
               {if(isalpha(Val[i]))
                   {recordType+=Val[i];}}
            if(recordType=="REMARK")
               {// Read Second Column
                for(int i=col1;i<col2;i++)
                   {if(Val[i]!=' ')
                       {recordNum+=Val[i];}}
                if(recordNum=="6")
                   {for(int i=col2;i<strlen(Val);i++)
                       {if(isdigit(Val[i]) || Val[i]=='.' || Val[i]=='-')
                           {valBld+=Val[i];}}}}
                else
                    {break;}}

         fIn.getline(Val,Sz);}
     fIn.close();
     return atof(valBld.c_str());}

double readPROPKAFile(string proFile)
	{string fNm=get_file_name(proFile),delimiter="\\";
	int pos=fNm.find("pH",0);
	string pHVal=fNm.substr(pos+2,fNm.find("_",pos)-pos-2);
	double pH=strtod(pHVal.c_str(),NULL);	//cout<<pH<<endl;
	ifstream fIn;
	fIn.open(proFile.c_str());
	if(fIn.fail()){cerr<<"ERROR:\nPROPKA file could not be opened.\n"<<proFile<<endl; return nan("");}
	int Sz=15000;                      // overshoot actual line size, C++ automatically fixes c-string lengths
	char Val[Sz],theChar;
	string tmp,pHLst="",q1Lst="",q2Lst="",lst="";
	fIn.getline(Val,Sz);
	bool RUNNING,ADD_DELIMITER=false;
	while(!fIn.eof())
		{tmp=Val;
		pos=tmp.find("Protein charge",0);
		if(pos!=string::npos)
			{// SKip Header
			fIn.getline(Val,Sz);
			fIn.getline(Val,Sz);
			RUNNING=true;
			while(RUNNING)
				{tmp=Val;
				//cout<<tmp<<endl;
				ADD_DELIMITER=true;
				for(int i=0;i<tmp.length();i++)
					{theChar=tmp[i];
					if(theChar==' ' && ADD_DELIMITER)
						{lst+=delimiter;
						ADD_DELIMITER=false;}
					else if(theChar==' ' && !ADD_DELIMITER)
						{}
					else if(isdigit(theChar) || theChar=='.' || theChar=='-')
						{lst+=theChar;
						ADD_DELIMITER=true;}
					}
				//cout<<lst<<endl;
				pos=tmp.find("The pI is",0);
				if(pos!=string::npos){RUNNING=false; break;}
				fIn.getline(Val,Sz);}
			}
	    fIn.getline(Val,Sz);}
	fIn.close();
	lst=lst.substr(1,lst.length()-1);
	//cout<<lst<<endl;
	int N=count_delimiter(lst,delimiter);
	N=N-2;
	string *tmpArr=fill_string_array(lst,N,delimiter);
	for(int i=0;i<N;i++)
		{if(i%3==0)
			{pHLst+=tmpArr[i]+delimiter;}
		else if(i%3==1)
			{q1Lst+=tmpArr[i]+delimiter;}
		else if(i%3==2)
			{q2Lst+=tmpArr[i]+delimiter;}
		}

	N=count_delimiter(pHLst,delimiter);
	double *thePH=fill_double_array(pHLst,N,delimiter);
	double *theQ1=fill_double_array(q1Lst,N,delimiter);		// Unfolded
	double *theQ2=fill_double_array(q2Lst,N,delimiter);		// Folded

	//cout<<thePH[0]<<"|"<<theQ1[0]<<endl;
	int index;
	for(int i=0;i<N;i++){if(pH<thePH[i]){index=i; break;}}

	return interpolate(pH,thePH[index-1],thePH[index],theQ2[index-1],theQ2[index]);
	}

void extract_models_from_pqr_file(string pqrFile,string modelFldr)
	{//string fNm=get_file_name(proFile),delimiter="\\";
	//int pos=fNm.find("pH",0);
	//string pHVal=fNm.substr(pos+2,fNm.find("_",pos)-pos-2);
	//double pH=strtod(pHVal.c_str(),NULL);	//cout<<pH<<endl;
	ifstream fIn;
	fIn.open(pqrFile.c_str());
	if(fIn.fail()){cerr<<"ERROR:\nPROPKA file could not be opened.\n"<<pqrFile<<endl; exit(EXIT_FAILURE);}
	
	ofstream fOut;

	int Sz=15000,pos;                      // overshoot actual line size, C++ automatically fixes c-string lengths
	char Val[Sz],theChar;
	string tmp,pHLst="",q1Lst="",q2Lst="",lst="",delimiter="\\",modelValue,modelOutput,fNm;
	fIn.getline(Val,Sz);
	bool RUNNING,ADD_DELIMITER=false;
	while(!fIn.eof())
		{tmp=Val;
		pos=tmp.find("MODEL",0);
		if(pos!=string::npos)
			{// Get Model Value
			ADD_DELIMITER=true;
			modelValue="";
			for(int i=0;i<tmp.length();i++)
				{theChar=tmp[i];
				if(theChar==' ' && ADD_DELIMITER)
					{lst+=delimiter;
					ADD_DELIMITER=false;}
				else if(theChar==' ' && !ADD_DELIMITER)
					{}
				else if(isdigit(theChar))
					{modelValue+=theChar;
					lst+=theChar;
					ADD_DELIMITER=true;}
				}
			fIn.getline(Val,Sz);
			RUNNING=true;
			modelOutput="";
			fNm=modelFldr+modelValue+".pqr";
			while(RUNNING)
				{tmp=Val;		//cout<<tmp<<endl;
				pos=tmp.find("ENDMDL",0);
				if(pos!=string::npos)
					{RUNNING=false; break;
					}
				else
					{modelOutput+=tmp+"\n";
					fIn.getline(Val,Sz);
					}
				}
			// Write Out Model
			fOut.open(fNm.c_str(),ofstream::out|ofstream::trunc);
			if(!fOut.fail()){fOut<<modelOutput; fOut.close();}
			}
		fIn.getline(Val,Sz);
		}
	fIn.close();
	lst=lst.substr(1,lst.length()-1);
	lst+=delimiter;
	//cout<<lst<<endl;
	}

void write_saveAsPdb(string fPath)
	{string text="from pymol import cmd\n\n";
	text+="def saveAsPdb(arg1,arg2):\n";
	text+="\t#print arg1\n";
	text+="\tcmd.load(arg1);\n";
	text+="\tfNm=arg1.split(\".\")\n";
	text+="\t#print fNm\n";
	text+="\tpdbFile=fNm[0] + \".pdb\"\n";
	text+="\tcmd.save(arg2)\n\n";
	text+="cmd.extend (\'saveAsPdb\',saveAsPdb)";
	
	ofstream fOut;
	fOut.open(fPath.c_str(),ofstream::out|ofstream::trunc);
	if(fOut.fail()){cerr<<"Error in write_saveAsPdb!\nCould not open python file to write: "<<fPath<<"\n"; exit(EXIT_FAILURE);}
	fOut<<text;
	fOut.close();}

void writeAPBSInputGen(string oFile)
	{const string A="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";
	decryptString(A,"783e007c783e007c","Password",oFile);}

void writeHullRad(string oFile)
	{const string A="Sz4n+IvIwQpFq4KqnxiqT4ZXzaqBg2O1eeZz4mfFGrJt7MKCubBhfFnV7dcOKUScMjyJSBaqEgNr+YEcgca0rTbvktqSIMNA893DeRX+5FwwvGgL9qqbEXyLJPMpI8EO+v8m/DwDpvY6OdBNCrNWeGyGYAJzYy/9biEE94vf+ScJgeg6L3CIuW3KYti1nZQQpeFEEdelb5gIRlw0j20OvporPaQ1FvnLUPJwGUrsi1Lyaisq+GNeDaUuCjrgcHsIb9uj/+D90nQIsvOln3KHnR6dv2c2lH+m75bV6NiGnXqKk9ykbRxMhWeJGjGDasqvq8lAh1QbO/Yrm5GoDeLBaFpcfU61aflnibSqPfuLVwB+FabDuWjU+8xB3Dk5mUTXICOk2KehANXa+tyBkSqeSygTvdDNRdidufwS3YxOIjV7An4Ni8o/QRLpL/dUbqYAi4cxCgsfUkVDAxQBMdg3mr4hhIGAMmtO3dV9dMDzdfQyP8ju3EA3DXu4wzLQaWmysMH6DoS3r7mS1vvqcPcKh3JY6PbEOqLIyEg/P2XUzOwknPD+b0WY9da6islhzsSP2SrcpZrxdcb63swOd9ydZ5/WU/ErZ73yct0W9lUWPxn/X2+k31uz7yGfKp5rP2ZGRGMT5yh+fjV+DUrRTQCnbZzarvMHJDMh34/QLCGYdMZbMfT0E412TtwuZPwT7ZnqJ4zrrQWeFV81KxEAkyj2yEA42pW42p2551qZVPSJPts4AcExF6UynSPGNHi34GyQabxv6E/RYQ/MiFwmuUOA3aCqUiJ5uqy4f4FsX015JttjPyKsJo1HstU91kF+rP/72rQXqZfPZmOYXS/wJLEawOpM5lLi48Kp5rxd5LrikbVyepR8UkbIJb4spAexYLaI9cuU6gYmrX41NXa9AACX4RJZdGW4oTbqbT4NRq4O4guEysHq4FKoeBheYg17PTsH0vDiwK2HzHsWYnSmpifFizzY7wLXVykFqcvDGpWB70eoTkloGM9ZmVbSoop9ll6F4LHnwHU8yNky+9WoiuJX4ryqxT0ranHBlCyb7k6G6MJnAkaQqPX/fiDQ6tmKxD9L6AaPvZwJ+o9DqH8hKOKZ5m7fsrRzIcSMHDLfzo9JoeGVACqvaeDTdBXu5pSXls0xf0ZtYe9ka4EI1FnbVdyvpL41zv5DkUA22xaLMLwpDCl8gnegykWL1lTm2XxLJzp7KNTwL4iA86uJG+LBOF2Kd0K83jeaMi1Dsq/iZ5L4WAW502iHezVQMSWfmRbRm98AH7h9Agus8no234NkA1aAEhb+FszSv5cw3Wn06z+6N+rOlR7ZBcbMQt6Xj8kgPwDMn4fvJcgubCWbL8pvSuiZp9tSEeWmR/tS4BZcdEnlnDtimkrJXM0B0qCFVNKLL+6XiMAgxx88A77X9/cO4Jt1p90fcRlNJkkorcpCwwK+pXmrUFNoLLYGFZ6FupzqQ6JENDs8mgyE+nRClLokuhzFJBS4TTti1VSaeWHG/n+XaYmcj9r5CFrqVwDZaCopQ6VJck6+dqYcfGuOdV6SDBN1pMUL7caUbW3HVqDGTWFuUYzGNCp7wY4/86FjUPkcQdi9nILbzb2rMwx9wXwPz5CWNnDLOjXxhrdRyfj78ZjmR7+eaUGcz+wb+shjGpNh0Fd2FBm55g67J4+HHe3VegWpMCpWLNS+MMSM4NJzF08w87tDrzq+nr46HPbQt6oKPuMJXDNTjpsylDSMmhvk4UzFA+e07k0ZwgA4GYDxhDqM2ozpGr9Qh7PF64MYmfPPh/lHQ6KchdvAxVT49QzbRCC4DF/yDTx7/x6QJ3/DdlcWNcfvY1cssIzX/V9kP3FYriCsdIr5Diuo6DXnMGyMsM0AP4nG20KfFOkwhBlgWy/db7BTQkjE8VnA61Z7UhrMcxqyMpbTUr13Gb3sD7sKwQZfg8uWSCP9klsGFaBQ39Expu8ooZme13x0kBlR0jJn0U69W5zEPzAYD9Z/jVw1MVF2cKChYiTUvJAx1kILj2vGDGhwCeYdzDth/4czzzkFWPdPg+kMJ/eicNsaB9ClTYjCvPXv3fP3sESnq+f5D+c4VkcS/rXFrNV1JfbhxktuK5AlwatGDniRtRsTLB1lKEYv+KiRBybD9ObTzEWLd1sna8Rk09e6lupIm3NZ+Pdr9hBE/XcDsUZ+cW1r/VmskiH3G6oX/f0ufxhVSQEZ/lHphKT3XSXPXN1A7WaeWxHpLMJbLMQd/Pn6TqEW71N+kmepUphAhU3PG7Bo4/kudDK0OD1uFgQVF/4gM73mJbFPev8P7DZnYkPH3NWkjzOum20yV2GHPkSYDw9MK1Mc8n0nY6TlWDMcX7SSaCS/3iy+SYjxmmdaUrjBmaoxfxgu4dsVQxTEduVzywb9kfWI0KzFWkcWtISLBSZhmPbG32RVTAc0TtiYk0pOMwdMRWr+d8jjo5O+UTU9SqjWYWsf9n0H8Qi97KlA0vzYY4uR2Wa5bsmaZfoFcxiXdy/YmxekPsw8vcKv+zO2Je7A+FMSCi42fcj/3G8Uw/VPQE+zE1d9Y/EknpejiUYgDf+T1y3yrXTpwbph1TmRxUBNnbwtlqGa/w7Rf4cZ2WftOrStzp5cVgIVvCI9NS1kxSt+3x4fdfCkAwkHBiFjBmcviV4vBMMnXkUOHptJWW5nncdIP4OxFHk/KcSCD/tLXzCilhfFN4gdfJoCxPM93TjiidJ0xPTMxjMMJk0BjU1jGJVLqB4IrXPc7/a0zXavGMNhLjMqRwgxo9APGqe/S6i3umfi0womGdfanBnYwjucumTsPwRo9Sxc8MuiSDmt2XA9gMcLEuAzYZGg5IK8HniHkhsZHM/EU9GR0mHIUAkjvQFuvSr3/Ci+lTK3UcEdHDsx30pjf52PBPB+GvhJEaQoy6V6pmjwu5fY0k3lZLv30VfQBE+4OmzO+RC8zKPaz6o2h+voY32HSG2+8DQl1S2mItU5RIbX8cyVzPTV2to7TWzhbyj3oBIQytAbJ2a2GFdv8GaOTFLF4fhKGw78+48TSQMZcF8iKy+epqWfHK7FTPS7pT7N1N8+w3TaUXCQj1o+AyDcfIyWXvt7AJUcnvT/C/xIo0OOvtC5Ur3phPkxKUST7EkFHR5ltl45f8KqCICIVm0cxxgxkN5u+bAY0UV2GXH9qRkY1QExpF70u4Fg48BUQr/1vGj/+qBWmUmBY6AD1oLpfc0P3fZKyxGZ6Wnd1WTyeRmtp3HYCFYZ6au/9j0fAS7i4/yQY9DlnpBRZo6N4loerXmYiqelQ4KPKwX72GDdyihycafMNVPZKwSO9ZqaZntYL/vQKsfvPFqZoTkMWq+X6N/96egbbSIdQrDtJuf9mMVYRrvjbh+bbs0sP/7R72YVvz+wmgHUOhYC/P/jr/YSx+QSTLJBy8y8tJtFk8KB+2eK+ibQtJbfq0SHFwKHEzv2gpAqWjdEDCixX9rN3bCZNYx3+mdqUIyA73aM9c0/iP2zQMb25MBw86P6jpVUIboU16Ydye4uOTKxDu4GgrIM2NLzl3pZxRB5bJzXq8gVdhK7qfXeYCfHo8GfvrU5DlDpLof8zS09Nsx4n/mt8ja2xOlwiOjmnt82OU0QAj1IOmGBfKuo4a7aLVbGRgZUCvpWu/81Ixs2+xfuWcO1ozEBPV5++lH5foV5agm9OfKZfXo2lowgRq12jaBHb1ZBi980OKoZwe4njxDV/yDrGVjkV4W68G6p9BIjw1s6QrQJqgRU2LABoN563pNlNvnJGgJ3aU9v71jLfb9EeE1sje3Eq/u5qPc5owwTdxcpJp8TQydAKiTbIfKNuTtmPLbpGrtXoKWQixvj72YZtOadO7lRzL10+RpR2C5X2Gnb9nLPEhiUUPerlnViag2t4Gs/jhOZyCA8uZsMFnW24k7br2+jl36eH0QhJQAOwzBHGzUOyxspDg12KjvgInPTwVZp4uy74ZBegEHDrIZs7k+4pQKy9L7+fTXKwRh6ZSno6TQn8fJQF6gojiUCGSKwVRUQjIHpatOWpYoPZvE8S/w1Dvkx/BnNk7LYUtja+KHluZ7JrYRn0zXxIHAQJkB/w1xdWAO2EnxGEv69DVkDwtBNURFro3c4fJkc4FLi7ahxJCFt3EhbkfWfQ34kCGcs5/r8gT3KOmEz7RqWuszZmvpkXwHeVWqSEol1u6GtWHbwUjFRFDvS/BesINqUy808vrRR0+mTATNC3sebJkQ6mCj6RAQgi0VedpenQxtTUvPG5aOQYNV+y9ZlOpPcFbH0aeyTGyVKMMghMSZMLeciHFjNe5Cm4UAumR1njGeBa7OHpY3IQw2FpgWaT94rs0Lv2vd+Se5WnOhIVcUAa0VZjRBJ83td0DcS6/Gqagu1pEAgJOHuQFPYaVny1l5kdQNSSZbD0NPwje/Pu/UpLf37E/WD6BO9mRzVM+XjvXS7yT561WM0KYDtHOwAgzM7bAYRaDLOo6lu4Qyuf+O/zqIEk1HrNko33fbqCSndYyhUTOTmcsNlm/ZSknRwe0W4tlqQv1DbUxl80JQ1MWFfTRegu5MAFr7XDGzFncPLloW5zn7lpnQs1uAS1BlVmg3lWkUchdJ0ngic+o+XKCfdJJih6ubeKSoNF4b/Y0n5OOWZkz1rlSCYvhWem02Sxqjia5ToGYZpMbmy85K6GVPUT2ZEzCRNhXnogZCpYkk3uy6sOarnUqc9XvFwERlbI/0bSMAVQ791GaoopDs7wml6N7bafv0p7KfkZnkxyAoYcul3GFOcUm+nYRMu/3p2BvlUeENLugP+rhHg4TZMQpTvaM+8++a3CPjT7RdW4oD6iFV8UZjOSIHGlGKPDXCx7vveabQgl6Q1ksEfmMT+THcJWwSCqFeNmtzm05R8/7InMnZKq4CSpTiVHzbcjUbzisKAXK+Cpz3Abwlofc7YTTG8A0dkqA4uB7lclgNho2gZPG9IcNc/74d08AMgxlBKMe1CjzdamP2I/0b1f99niS4qTwaZlUnYR2COElR0T207pMEvvDRbK0p0RGVzITb0jww6qdQq9b7k92k+l7CWnWniAfu7/3bm+t8XkWjzRTTMZsChg/KfekgDXpkdRlpDcAWeh/GW+1Z8PbsLhEh2Y06Xdb6cabB3h1Oc/KLOuuokixUh0Wt93QT25RE2eskcQpTRav2zse4ymEavYhQvGF+KF6/Vwp28nLFXOkRy1/uhiMo8i+gVukMS1LfWs3nmswn99RH9ZUh76mSTZ+bf40lDekL2DqDA2m8boRbYDRLfoGoVMO9S5tUoDH3QcLSSwJ5+ShQcSsLVRQijWxg4GEceyGkyzLg/EUYCMtuTD+FDXeOPWQZR4iDTkJSZpEzZFO7QKjNs8SxdywDrHep2BD6qmsTaV0p4bWv48vmXHIswlyzGzIJJkr/TAXJgfRiJERbDfpRrMNSLioZmcajiXSGQlwj7786ras6p4VVzmHcoiCKWkfXRlzrbznapJYT8PJrv1zwnWmsqhDrtfcrVWVf3Iz9ktrRFoeh1D9gBJP2aNLsrIc1AF87oiZvo0g0jqzNzxK1zbnCZj4UgPpHcbwRJ6dggwZq4vhxEBUHiionSk32ocQE1pT/O+Jh/dqQoVrEhI4lWOncIOt43FIh90zTlLd4YIZV5RXyjECePEBGyTp8Em/QAgwQCAgNBw/Nvo+m5VLTZ82lhbmU/LYTqnaJ1yIIcjkA/it+8kNuCMtuCY+wF6y5FNb/ZAaSeyuTwUOI2t+qp2Zq4VpZBA+xyuhPPCfUJcf+E1ep+QnF0eVly9iJKD8FSM86jgOqh9d1qm7b9SeYPopTpbLXnhMh9+NYuo7sgZz52CZ4E8TsKcll9EoPTFn7JL3sIy+iywX4L7q6wdqJnfoDNKh9n/UxGL8YFHDJuh79hqfrSu9UWED16VjQ8w6Xh88LYXe+PR8nzoj3itLtLmd0Nyrh4t8E3QRNZPOCS9IlkHo1jGbtvDmqUf9B52ukp/fpMpe+TdLE5Wy5KdeBDN7UnQXmetjXuXpjqnomRLwDIVhbr0RI0yacL9Z3GpJ6uv/taAk2xceGatXvgZnu0LPlIvRaz6JteY2JDiSWdj0fvmEZ1NDV8cnWOuk4FrD9opDmL52wA8yrM6UdigzKqV6quGZEDHu7rcy5c+VMTE3rQbwcrBT/B/y6MAZN2a5Z5BDNfHbwVpfYTd8koZVyTUS3HbLr0MdI48zwQNasvZS6I5UcnlifKDso4cq99gmCkSxhEbVV8pU1/bIJacOtpEFelDkKgCKU/Iqtq1NyHRhmqqAMEg1HMTeqVKl7T9R0YmjfAPwZFbdrwVAq4KjgfYiTPJa1lunHdr4ElstwVhAxWYozj26VlKwGnGSIlT40wINXKZ2jvy8jxNrkNn7bsvUyKJ6azLOkLGJ5cicHFdgAuKQxAzFJXhPs1nbn4Y9mCVQ4CQvXiVY9+rnx023iI2Mp9RF7AOb8I8d/SQEc7TabKJPq/6Dq35wsZqPU+cp2EyckOkVCgif37EPPZPeuuphpjPAWyPQHiHb37iHW1XUjgbNJ+mb3Qi1Lsu9rKI6VU+2YuXjqwdcpLS7rqOQwyJwKCSDZhoTY8bSnPPK6PhESKk1mvDeym6diC5Is94WXSpap7bjO49HkCkCcMmhKTez4jItWTmOPcK7RHmwnrlQZozQcyquXvpzCdN3ms5yjx0PzgEfyhXxwrDGZTxQyWwlzPUSY7MajbLMSPg/x23uV5J6z/MMwQ7gRSqhzS5v3mx6JFPtkdNe1LuyNy6zr314ADpXEcAqqksOvwzT5L8UUuqxTEy+uQXLIwUWVVvXMHYYaCxP9wpnc5AxajW1joxZZPyVSloFs3Srf218EayaPEfTjApvSMkIxyr+G9RtvkI3+sDcMSRBWoDYUHt5KQ/c+Wvw9wA2+gPrXtgY+Nhn0ghTZz5U9CP+9JY/zr/gc4P+fuRROjVC2xnVrPQA4Mzk7ECv/RUIZRooeGlCNOEpyKJ+uOQtcLV9A8Ql+WkVZkEq+eGtewKk98v40Er3q5KHVJQxk3g0D28jQpkh5HUKnuiA3OMGc/8LaREagMN6ZTOg5m7LJ8rOiBdlRG0movtAx6co6jqmAkq3TQtwKFp2bMICdYzZkhCFwxfR0vZFHbUmKdcG//vFMDVV9npHnKqZdel+ml1Hh4NLjklNYSduzHPnbDdVnYhtno4nVkSM66aRALor2KJq4sg+mPiKXoFaxjAhvZYtIw78FGWEKOcMt0OgBTnhJMY5DnvdyVKjK9jp5EZLjQqu6OVYO0+gunMMNrc/tVdcjdggYwprO2nkZaH011Ii2eqwHHfRnF9d3YF9X1/pNuP53Ol60HFSDfgXdNsfRXWqFES3fPCbXYaG/T5LTVaT7zFFbsnYFZoU7iUTG5L3fv8Tx1E5VL89xnSDwTvi+GA8WFyz5zdRdDgtNMn7R+5EBkxV/35j1UeEypoNFDz8BpXYFBNEuTxjzE3c8j1bYfGRtg75ATy0WHpVJ3cyBjkpzvW/S6yiIJmLnYw7nQKytOk0CfQi3DcjVNgg988Fh4HKONEw3p9x6dutMOSdVzqJJy0TdxZg6H3TpGswL5YPNfsZ3z6uRgAExYWQ3K5+wFMsLiGh5m273hp7sDSuLo+/NT5dYI1elGccRIy12y4I6ukDW8IblIG+4GOdzvjbLSOPROs0KYu5ow+R5JpRa/MmplCKW4AUrW2ouGM3yylJosxrE5QvvtGxpzh+5zPsm8Qj1jBl1VY0xKQoCZzOUjbJreYkATOQPVFxkgctP1mSriOKaNx6KIZVdK6oG/kMFc56G2n659nxDtAvSvA6tBsK/I2Wg8t4kGpxDz5eDsdA1dOH1itbG5lPzR6Ljv59IiNYyvVfI66QQwFyyhQTal+FFCof/PJNb0vhWil+igy07sC07oCjfVoIroUU2X0SJtqz6zeNf3E5JI33GK0ix+DOhW76BC5mFaHgEuVm5/dbp2X6MdJm3yyMBpByrWinTgJy+gFg5svmWTH8MbUy3AcILa7iJ/ncY5Nfd1vgdBIl5w4dHp0+2gB+rFjN6tSD/aKBXT1mxTgpfG3+uNgJFQzHw2IZ7HaGijlfGb0f5F8CMbPZKZ++M3jRkP1llggwr7jJMKf4fTmCw27BhK1JagpKSkvW17jbuWSrg8T5FDqZXcbPetM3TFZT7HvXOVHG0ceULPf2Kvf3Of2WiJ31Bid+8O/Wrq5PMpf3LkAbzJ0JuAFCaFApcLwDGJhe8OEtzZC3X7cxPYQZJs6jGEYD3aPH009PAMlivWxm10Ieh8mnNv8lYI6v74H/nzv3fFT+IhauTgPHIgB6ubdC4jWZhU2hHcaWtKVmxNNjlgfiiBmmwlttz71Q/zZkmp0hYjKZIhT/9/deGWpZ7b/fqEwwZKFjsDfgtA3Slfy+PsM9ZDZjp+lWPgWuel9IsW5PpWpzB8kI4kS/TihiVjIER3WqWEIqOhl7H2VDkqQ3wz1PoX1it2PUSpFIo5mlwe9pR90HqXsZw68po20FeHYWPDhaU7MaCU0w+fsksOaMRuJfrR5oRZnaaMC0mtT54dmm3HaNhCU2XmFN2vq+VB5bbsCkKtD1JiYaEphMGsngiRxyx9mB715l+2Hsvk6Px9NtTFNh/3gXrX2SaTyjwFNeH4dVEZ92Z24TYYvY2ZX9km3VSdOTxAI5ZHC95ORXOSCL2hs1+62KRxSAsOpjWcrLOdDqbiaClKd/4S99Cm5athrRnnIOJHzDAcMeR7p1fZAhyKE/VoeG7YrAcW9DYo3KCTbamSObzDW/JMNiYgcIldFjLWevsFofAyfacpw90k3E3J1mRB936XrUseEve5aNykmpYR28COfd5d0Ey5ZyiI0tv6F7i7LlJ2+JB9HUH4HMeOUnsfRdh0UjPIIpfxgJPMh7jK0kV1kbrn+5HApMEM7GDC5H7tNSr3KNT+6+wu4Jla+ZNOh3tFFbq8sRVdyJqeG9zddO++Yz/cbQwJ8Hpg9T4MbJu9p1tCagyvFz8uN3Zearx6SBmPMp/C7scPfP8VrG1t6xrVfBUjZ7m09S24JgeYidWLbDzoaswmljGhYmvOf91KND9SdOf+y+NwyxWmOnkMY+fB6rHSYGihcRoUKox2sefNEWA3+/IS4BoVnV+dj5macNyD4O+Pb0P7xCnrbcts2BlaMo7zpb3o71fSTzk/KB1tgszRyQSRzZQSVFYrfvUXjj/2l7oQwdBopdcF8sl4bw05RnYUtwCVNKjY+UglgGAzFmnXsDYQCFvhHgnZgx+FU0jqEsVqa5cRwfP0Hw+mBWRJNFrBCQcLBn03FkupZZmZ8hxunXyfaYE4K28iLL4NzBGMbid9AwWHvQPOYj1ajNu3MYsvSHAmFvcl4NRIkKQSKH+5alF2/Jd7l2WO0d3eCPABhq/2K82ANQaw/qgr52LbTPPbdSfg5GvkZ8ZFFc/ZKfZ+zENj0M7/m0Q0W1FbVaISUZ6XroloJUmL7/laFagnmrqBh0pGmydF/QQH0WTR49+J3Nce8s6HNUUwJ1zNcvWdZrd02BC7L9EwagfxtfJqXEXZ6xIIMw5xPFlDG9Udso5As+k69uAbkUs/0wcmnHs4ESBlv26CsSiIroWX4PFaKfexdbWIdllUQ8gdZLgl1ub3rzn+XQ2q1ICWl3G4TpcrXF4OYpMmBTcEdiNYMb8Zbpmheoc0nhUJbG8UUWMdsOqPZbBZ/4ZtQ6Q11fKteK9gcwTy7sV2fuEbyko3J3fNaeY4wfBBTAAp9YmJRbGey8zt7NNUBQcCRcEtmoIX0Mve2J/VcXXFucDljKiDjFdyNGqhlODhobLr7lSreoH6gNemFASGW30tEgUX8j4UDTTqJe0JXd2y+JtNeeXrTZMIZGkmSbK4o6ch23e1qlhG21IDPzkV/nbIX9Llixq6xsaSo4XQ+eQdB9w4JXmasIZ5XVytu7S+b2S1h+3RKED1dgQNGWpx99zcMGuWcPwI6h6QfOEZ0bi0CXvuO2Fot5jgDIV5NmHcA4Lz8biNlFYMSRQX9NSWhjWblVaxEGgSLv06DkjmwOsCU9QhSdKhuKrSWSyliLD7/SDFC76vM5WnddU+Wn/A08Do+IMzDeXUtjWZ6Tfr/EzS/foMPe/MXJox001sS4j3myMrHbWT89fvN1iOnBqdt/1RkJQNjUloCyDDnDjX3ILdtejanUYe5xmgfpAPMO1cbACjGydWSWmohOat/qxCc8eX1rnY4iDg7x4a0qplgVSclRhk7A9q+PFlZgMOIK+Ho3S486L3cWMYrnY0mqmjdjqJUX+wWkijJzjemYSiOgGd0VXFFsW89SBQAr2ue5iYTu/g4gIWhIW3F30Og5TUbiFyoOcR2l3OSlEnZdIy1WZ8DAugLSATOpVtwk0WaXt2mV/OSkGNhKje6ZxfqkyJ5lq6Uwk7FjtaOrNv8C9pCfD1+zh4nWfJ/XPXDemgSj0jQB+5VeFTripu0triRKHyk1vTbjWXwIP9oqqN49kDgEqd82pHxem3WKspWDg+bKrjvKEahljs0KlpFpui+Mk4cJtxzbR2rOK1nBM6lvhTSStBonVJJNLWRfTBNVAHKAVTv/H7nvWahs5qVRQbjDZPaz6Et4tLhRiEZ7NE+fhQWcpsudxmeP9SZbUfI+2seCTHL36WuRMMYBtqI0PlO5/8u3nDKMcfOrkwNaGPoqFVtzi+aSJRXzx1FguO3MeNiceplx43z4pGnPadbCBeKw7rW3supRWLbn+YNQw2bX6w3uT3cWoF40xUllMctR904oC3nLHWKG2k2XqgauN1X6zmOxT/Q8WvqRrt3XMPrVWB8ohQnYQfRFB2wuCVknIyg9DQ7QNHDonCT7CPrETQk1uSyw8lKvGJwywjZf0SN4Raj+X7KaHPAKjAbjLx/d09ITbR8paFQVVti1NM3WZr2jBkqvgUlGuSeahajnONxW6wXvCOWxxffvZvGbqXrO9ykJeqjDtRtk2r7Liz88kq50eJ/WMdmYlDyRTnFYroLEaycBIoqaJB/rbfZAEW/2+I3H4LA7juNnFtIScTKsPLUpvsIX4yYI8uBlv3PnBPvvB5eE3T24vPMqftz5hKm58t5hE63sQzSrpexKV24N+RGN8NZUU7dA98lPb/VK1wC9K+2sW76wjwGBL+q0qYVNbUKJroeQDAzQaW9ovYva5UOPkVile6NK82qhZ0dWKehvFjyApgIYlCOhKvxEism7UQ3w8w4+w58EZ+Z+t+a+3c5oGi7QRJJKjxpy7UvH2WZRamEomjmV6xqwc6rmCfs+7PGZowWZcwUw22DxZ6/b3mnEMT8ZLz9qFwjFwl1PR5rQONyuoP9OHa/mR3PR0/RbQDGC8LhlLdeao9QszHy3F3924LArUx2lREm/2dRbt+U1eBdQdwYBlsh/CiH7UB2DmADPR0b7JVSSnI6EHm7+G60WDP72wciP+NvsJs7/KDvXouWhDs4yzjQXRDz+04v093QVqK1YA+Kj7LJlrXjGF5WxepNZw47tZaNDERkpill/6BpoPP9PgLFUjFmCe/OQiz3fmPykioStryytKdUEkeZQEJdkoQlNg4Ejq1Gz1Wr5ouvqa1qRiPlcYXq1mf1w3YUiV8ktMGPOC/vH+JvJg2dPQ0tn9shb1A5wYZ1TInON8hPjLR5YpUtiFCNyiJE2V/m5olwRCdvmzjBj4w/0JEIwFmnAQupx80zj7SBNo5O/aJmeyS16Z4Ubt3XXqVRHO+eHzUFwp8KGW7AnXFW4cllsXHlAdirQf3u/MrY9PmxOAZqet3f8pBzZaXQrPG7to2h24Ea6USRS0XO3MUybrnqO8mHSczqE8ioifUUh6ZDsMh543s9QYrZ2jw/tcr3YxN73LlrvfoKmBP5Yj0MZbnfkOxYPTNafrMurz4UCltdNF66YxHlhcTDwpTChH88GQJo3Tsx/crik1VgHZ51A8bFwcp8U1LrgbdXw+5zIhsISnq19hrVqbxZCz4fwBH8NroVO2AWq6jCW1s0b0tzFpa9iqo+6wpvQeKEeTEuj1z1TtD3dd8C6F/D7RoMy2450OVBymx9l/LSFz7VvZiWJZYGVxwYxcD5EK67NSCvzFo0Ux8ThGhSVm7F+/APJUM91nHxgpMN2WYhWJ2hCPZevuGWON7Jt6fVRcABBBFgsQ3tFHskY3R4hlIEIpCWL3bdGMmRZjOBphIFwTb+u7xQSDSHnthr6RboG5MGCSuTWj28JbGXC+hZPkBvNjBMS5Bnnm85lQTVymN7VlPA1+1p9FrTzp1PqWTQT/4xVhCUef/7fVsl7SaytHrZXwNhJEF3EUxWpFi4c9bTGfZ/nVpy6eorLy8ABOS1lZiQKr7sFVoSICg/gS2A7dveMB1VUJOchZrk3BPnSc6iOH0ougqx8Y7AyTrMUrdgga5amZNsBssPo08PA6g0tb5nFBamehOgWIKkeELESwSGoBKaFodmNqKK2nWFd/2Bu+3d3aYoZJhF9EWPOyLnDzzDVMsWIsKt9alZIqhtbbasM1g5dugA4rmuu/aJa2Hy8zU52KLEPnwyaIQ+b0aBhsdwrijNysNR66quNXNj5bWuGQCGa6ArZMgqqbNBTf4aOfqqDaogvbKtN5qdu632d+XhuobNOnTL2uzOQfXcKa4w+OCbyhGcu/Kj8YLmqomiavjkzyRBmlUfGl9xnpCEDDPlwVsOSnBW7LlnlVYzL8eIayGCcLmTeeYSVv7W5/EfLTy/qgotWmJH84zPATIiScA2GAokkdHFoE05lLsDByK8lLibCpXSlI61W5bHTnsZOqZpHl4vVaZcZVguDPQsIoYbv5bCDWVVtIYhRuDe5LkaUAtu2Xi+zIKxFQd6jN7hhafwL+IBl8ntxwBzABojFm7KD8x7MacebJbpXQqxzRRaELPChDYycHgf4D9Wm3GfjF9E4IgceO3y4TizCr41YGJO6BNtmew1EazFcONHW3CLX3XasRRPYELiad0ulf6eBb/MFaCKprvSzScdyOLDPjW4CoKuOh2NyMQSvHvnDqNB2vGtg5YYEBqKKAeSeZZQ1ro6caQG0gpiFgcuLJ8v+L0q64rx38zrb+LO7Kn2uSGoQcsLKGj/d7NaKqg7kOCS4S+DeFU/bW0iAwmePnbTXFC2FUa39ZwTOxvylZOUK8u6Mp2jF/LgT8bkALnOOparMgWGSNN9taTF0s0nxxGOS7k0zgCZLyW5Bil5sxRIHGYgIT1MqDJsZarShedc5NqUbG6l62scxq2GEWihvdO/wZn1n6h6+94C24PP+vEDHEvjgdzu4UMjEFKoUtEzRsbx1xr+G4UIp/785L/Prz23I1+34ljwqO6dD2bjA6wK2gEVT+LBD9KmO5inIiy0tSxYHcaWFCoMqrDoPWdY+rZLRui0lPqEtTSBjqes7cmVr9iJZSSuuB43+bEU0fPRVOwD9q4Koar4MAsm3WSpaqYH9Rkmvj2j0AQhRHNqpSoHREti8fE+aySiNvi+dhsiyDwoi5J2OocIGHgJIMOrWI0jjVket/MkSeGZH8QF59DoxxiNLwwMKe5SPH94a0gTo58cavrz3m0O7OD0IA+RVMNknX7jvwnqOGUlccwHgQNelpaoxBgYIjozYo8B+j0YoEbsJDVnWICvARMJY/UvWk1XeQ+SYOaUTpTqzQwwikEE4W6Tc34w0HMl6EzeGOygy2KGxuQCBYSDzWDVYtpKJxWpnvF3hO42nKqJTlhcAwAdT+xfuHBgX5qKrMWrx+sTm2r5i5l0rCGw0KWiAhxqqOiAYy5G/yZcu7fGxPlXP+lgCviTStXbrp+TRAOJjUqdMnQ6zuxMyrUgD0mf56Q43WrhT30ipsezTlpuqT5AhT3VlezD9uscBOOChHzJnPgR9aNWJvNbuxe05O5PHmHbiGgvHyZsyQd8c7JhvU+HgRwoRpQNfPnQRGbGOt9S5cdACENZTmyKm2cwfNmoo8WJ6uEPqA7MTPge14sayOylFQkZwFbiDdW+uJznmLIdlR7nT/J30BbfLmMy6cWCkKcq33tQl14C51sSrKm2Z9Xe59vkU/F42WNd4Q+OKxWH2tTWrpmdF42r4kEWmd1R/yhH2buT0hOIF24HcKRy4+4jsT7ny8dyDLGyF6jEb+nx6PcznYLHPJi1TtUJN4bT5Y0Be+8z4zuWFtexfSwhZ/F1vnObBezBkB8+CzC6toF90PX8bOixb24whhTwkrb5Uq662eY1ZtI7PqZec3tGCYpRTiY12M8xdUeJNxjBWMboF4Bzz3J2LWXHeonDUfVW37rHGv6jIJyAHW2gBOdMMUXxzvacXUw9cPqpZ9Ev3dVv7Ys67xP4rlSC+5iBzDsJMV6dRb0BvNPR7gvqoh0YHqwQytBkzkSKwoeWbMxHNse/IxvYBYw8eDWc2mjgtfE/jVvr6Ol2UPVBkdFEm+4hU7cQ9jx5yBJqCC8ds2KuMQf+exGPlxUJK+tqlLQNU6P4QxCSd+BAE60Esr74tkUhNCF/fIXr3cANKIaaQctbav+Wpq9nVBOZlfC0fgnjRV3E+y3/grf+jhVKLO1++0CkQlyA02MZBsPVa3D0sljVuOYRMuwZxCHaFuKsFOQVYtXRXwB5Jtu6l5voLFTCDiVswpeROwZPuUJrgX20wEvVB+cv3EbqCLSpS+9xp08/FEPpPNaek+rIn2ieb9f7OtiLQSr74AwAsFNtJGTZwAwRtdxnapfaWDyJmjpDh4rTBF3EiRDSFA8w3EYFOrCpkIyGhi0FFjwVpcn/a4r1OFY3C2ljAYym7jOLbaBGjrmYEXHrkgEsSG/9pqbjS75lHm86i3tZFKW/bG9ZpGrjts8TGQimAITnjPvMVQ4PZYNM7rCX9J1gmftffzdUCZyGZdmhOmT8erSVJDd+MaA25M40ee2ROECZIc3BgMHT5cuXYQrL9LWwLfaY1IHvraVX4W5qmHrHPYq8lSFOrT/a3WYUqfggN89kiHPenIN5ZnqWDIawLddvbniK5hkxf5JJlECbK7XjbChtHYT9+4irH0b5KG/KZDgW2AHQma/7mlAEZhnj3idxUyVWzrwf/ZzK2WZFaYUV1kcRydBDEWYpsN++ZeH2RTIK0WV2h9Rr0xUEokGiZGO5hbX4lqY3ycJn3hfoFQbigmPqPeMNr8UYPdB/ja91/dNtsSKAY9zWG+n1+YCTEOP/AO0XbNMlUXEF0r/v6hcM0Gwfo63KIZmipGsvlMmMl19aou7KVKLXRzdISbbR0QtFoe0mG0mwa0AfvfzjRSSoTkMv6V+c0RSiUopoOnzjERd2uo8CKzGefcLS6zDvs2OWr1rwER+2VQf7Dy6HnA6YF41hUPDHBifnJGfHECVZ5Ce0XfV4Fau1O2SJJmbD7KWVGE3JRHlsJle2XN97D+CUCtzbfjX55XmS7dAVPob/qRpD/uqF9DE5S9vdfIzbb92t34d7HG49PPY9Hmn3vzhv6uubzH2kT3qD3QUgtCWxw9ufmJL1Dy2X5VuQ2qPXPGA23Ysxw8S0or9xRpujHNTwGcMUkKsk/PnO2sfXeTg9k0taSBQ/4LiJSrXjVPFKPi9g1Saqx+ZaxzugDckQ59cbz6c/QCgcTJvl1VhSRKdmPkg1OJ7d0OUEA5JBLIxNCOnJ3v/YYb2cFd6nNM8/M0BqUM938GeHkJBAB4SOlDr/MoPHXezt6ofwbLsy563jSKQ+umNhRq0lpUJ5knNQ20PhbeLZQZOFPaLA+aryS/BBYq92GgH0fntDsM64r65ntXOzLNUxlI8bNqTv24EZ2atLPENCEMNjMNLCVH2NmwymVTt26V4RXaMyA5lHJV73iKzbxZVhsBVG0LfHhCUf77TUwo6nRRf71KLqAiYO11QQIU02bSAr8Gz5/6nYjiE1XWFTa08D7nqitcFk8cxidgV5rYo7GhhvWsvjV82uhmO9DirOk20nw0bKfGsn7M1EJ5IeA5lmv3J/zChomqxzOWxN5+D7twa7YzqqED8bMxSsexfJihXsRQFMNfWY723LbJG79V3q4XRllnIX4IfiBVlaBYcEn9lBLyQeTZ5S3a7Ss7MCIwJtm2+lUADjF3MUXM4CJRIohOTKwqcFf5aNP8P8ygtgq7alKVM65e4AHmanDVd4xdxZQPbtD37pUnmBkknqbUu+XBZ7WOVkZzW7e6GNLJ9gU3PpV70CRtNqj3//nRAisPyf4ysfy944ssdtEK8lVhwxB2vSWqh+L4WZzLeo4wQ8Hq1HpF+YQaCGYkyksC0TcK2+XgOywqUqoqT7vUee6dWGWaVVIzroXMWqDISpWjP5ndgZoX9EgF67E+Rwl7buHxaTUgzRigQ3V3ZjQh9h1NqG5CblSJO2N9zfj5cWm6ObzOJJTKibTsCzEzCiwSmvmSh0bh6ycdngbzFZ/8q2sksBrZuJe4JqKEPtWTjQzpoTpNf7UN+OO0NwfX4dMEJV+MNPJCgicwCGpO4Qbi1LKm7srjOlLWVKJf5F1WRnwvt9zvKgNi2psXug1Yw+Scu2mVUUlUWsb8k+hcA2O6Y1tLNe9MmaVIAA3oiSOC+vAmq26JYHAwP4UbkFBzoCdzzxgH1gRwJvHW8/l/FzO2gcNZFqVzsU2NfMTVWEdCEsTPGRXOe4No/EbXxxIuACx7wcx/j4RAE7nzDMGOnYPOOqihcwLLGW3s1di/TfF9eYGx8dgcwX82EMlWNdNRfM+HkUwSFi16J0QE0iVqaliXGXkPzaMG9bwsK8/gw/b8WF6KzON9YoXXyzvmEJQ+dAQ6DTM1fAWLP3UhP4N6MUj/j8y7W7r2GesshEpUM4J7VviiRRCiiYg5ZhqNwiOsaWDryVHorY4V58MY0Wa07Bw+bXZGoleVmo1qL8GVI1seC9rPjOrD5w2kfi1ebR+TajjxnF4Ghnc1FXro7Ly+NmWcb6ZQuSB0VFbpov0/bQSOjnRrLpidDf8KQDqWdhTJXgb/KwzdVwBYRMVv6Hf8Ax+d6G/ugA8Z2thC5tc2uHvwgT1fpM8uFN01CKjEvcAChqJKZunFxzN6uGxXvxRZslIdQjXBcmvPwyYgM4fMnc8hnq9Yyg2nvt5JUsDgkPzWxUmf4lW8iDd/e5t2T2nxUNIBPPfssyGAMlGzHDX9iN1MX1G6VL9Wm1wJkllEmceKLkwum8yIUmEafeamYpy7RvpQ5sQsvoLDiZq4bUW6E56dwAjUBWEASZG3Ezp7g0AVsQochM/GHEs0h3goXzqH4dsrfasotAE2qbeXkyg16bjUwEqsS3EDTQBXZn3cfqVGJGJkwPjXGWeGfDQH9xh4PCseShNUGBHNEt+Fp6NUg/St99kQMNUd1uESyW05mlo87KtWmxGQERyPoDCZdL47uk0qGz05lP+S6IhkQotj1I5t5ML3GpB2krXxP8DA9tRS0BbfqnWep11cLixBgfjrcoWcYMcSm6uYlhkE88vrxiZyIeHFra1M4PXSHtI5MO6PKkMARJU/wTlZy+51F7kb7BGLKroQsubS9FqHdbnzgfCUHlwH1exUu3px/DxQjDkjpGrS2BYzPXhO9fhYZPypiQOOpb1rGumBT1ML+NOBphH1VyBlq273q3crkF7KdnzOxuUDBPGGkyHlvIw+f4pM19cUAqIB32OKwMFc8GbOd3iCdUnHgB2T/4cXDKsMoQerdxOT8KecQzVnk0UcD2CFM9ICGF3aBMWgnP7GZfpfXDIXqVEAmKBkPBinhL+cEbptjJhy+iIEheyjTHlBoiWI6/wKpdHszKP9OCWYSuNyhw7Icb2k2t8TVsY5zyIT+pfxKsWwEOb5rXsGbK72Zyrn6N+mJUjFVbJGgl8PF0h5XuHJxUnEhbqMZs+JK4XcL4lBmYV0ROe2RtCXOP5fPjnpFB2mdHE6d2LSd0jZHd2CfKAGGRUmIW8R++M3d9ugqRxN/XA8SonSrPlJZV9Mna1QL84bvWpHVx6qwbR93TcGsRdn+Sm4VvuW2m/RQJ7LKZnyJcs+6i7EdmworDyeWWhEp0zXYg/SU4erOPyEWW7qpQfxiDbrFC4Ss2FBMWhUfYYOtYxBFJ8wZbU84kaZF7M6qscit8RwgAd5TfCkZBBP+z/J8G4GHhaFsyh5krfGm8lMuthb/BVkwJEZlRUGm+JC+cNltDfjhMSRcVDHgOF1XLH8fvXqW/4zRd57/Jv0UV144Zw/WTF7jfJLPGegs0381rxFm5SlrTxYtqgEnIraL12O9kHf6hHgELXJnytRq37bzfgvfvGe1Nj1/7IHkFzPpFoDPqqS7e3U3UR2KOt0c+yP0J/aON6JYeh34T/AoUj/JPb5gVyj0rTykOSb4Vpl17eXIjCCpu8i8Xcr9/FuiFdfTPTRHCznYJw5UFpZdZbsBsdnkiEfMBB7JrENVMftQbzWJrZkA9Dw1nvoV4kTOqGRKN15zjiotrvbwkGOMUJGwiexj9PZ+WxGmiWnHooJks3sWmovqRoJ+ab2dCeAvKeE6Z+aixdOPCg5Se765FSWBIEprNQ8JbWxKrGGEGGHjopzkdW6ihcFmntJ053NHmb9CVQulzdQsRv7HbGZS8d5r1nVXcs9a/BV2VeIsewIA9GIwHfA2dsM2RwcvQmn7RP2uBiEgPvikmeMTouiqSf9kZ/ibcGwEfUikEsmxSE3HwGTwDgyKvzH8JhIG9grGG6U+GHpkJQuBmkn/7CTfB1ljiYn1AsqnQ9BtNLR1bKvmvFYYbUhLY9N6T+SYosRBOkOCCwwaGsy4hnlLMHaQTRIcKN8DfPK9JD05HfikmaT4pDLIcu+NhBINcZhArVWPw9CRJ5TZKMX+8YPy8ln+3+iFhp1Q86DRKVtT04bn7p5LcHkh4gyXtdZyIa528UYXhP1tkLQ92H3vyL1dZhFGogLqLnwau1Po0y3mbo+Rk8tX9QMSPSiaNDsWWFvNmYUC1uekQhUV785Klv3N/rBK6Yo19dFqwJt90SyEOQ9bYw/QvpfyU4ePpYzWSVrs5ka8SYBz7DqytY4WXqsEMMG5A4y5MR49FDTNf7Rpbv2wI5mEX+Cdu5mtGCtIWwLwzjmhXmntRQjauf1DZBfNrC9qouHavrth+NNBlxmQoIgNrCYnpSZjqVjw4qKiQQaiqgFhUcrCsR/CipfEEFMHxg6r4ydz4QJaMNDMkR3Fo9QSGZn33+VFDSqV7QF74OFq66F/8pltfdwG1EK7TE59rqY21PWo8Fi3NzNDTJhqcVh0goVlRHG5+uVwny4SVcHGD0R0E6/nGIFXmsXyO1hTI29TASDqnZg6akD5d7YxgFzlWS2NMU3ZTNao8CarLS5TDrDv/J5GgphA1xT/iZBcoCd2fX7lAtsUvXd7+QIgO6Vl4Y5av+gQyWVn5fmfwCy6xepBsOnfcAeFYS3zyuO744ihaU3X2j/Fk1GsBXRmz6I0Oc96VgqQ8n1O0hbccrJFcMzAYRfv/d756l4c0HeXoDpw5Bi65pkEY5cZzfCpWrQXtX752iwU634+IpeIOP9B5UuyjLS16sEUXVLvliHqoZCeoCiMGm3OtYz1t2bLRQvHDR6x1OIypAk0PwBF1E23eFQeiLuxMb5mE9PuuqcmYEdF8cjW7N/fkwSwWvgPIl64OXoY0LDLVTwvAv8KEIeW13XRkf2Hw44O+EUCXZg4hbTBeFTdmIjZB91elGIi/C15YIyAbY+G7vJtxZBX+lFsBj9GWkV6mnchpv416A2n6Br+DZYtt4B21BGa2wP0tKpqhbyceVv7av4BO5TUQWt4nzg5aLpC5o1IwCmyYEoNJPEd6K1Kf7ibwXKtkMGLETDoDkpyJsDt1vrx9yhk2ubnZ6QsL2PcAtXr2GBDEdsGW8/nrCLD4jIU1i5Tkt3TNVnrx7POUZS8X3kqhaxzkSw1RFdAUxWfmuD53YYPSNY0kKfHDZTpWz1901qCRC8I01G1qSV55kDqaATr/xngAYmQFFYiiriyjsFDDqGEDT5GjwJLHpv6xu/kUaH5o3LaNg1rY9xr24irAuMVUBpQodlRxcQnKcYu6L+0FqGAbVbo+VVMxF5fcDkY74eMXyMjwBt3t0ePgdKXTjaRLade/psU8TfrG3YBBndg53n6RdnUybW+kRnNFd9qoPgWng6Qvw7haRjC2IRf1s4NKAl2M2dYX1PazAC6l99fqikh2Tgrw8ByIT+FP6/LSEZg3aOrphdWE70jibnvHXJpggkFs7KuVrbM101gy5A4DQ4qXjO3YqFRM/0c//XVK3zyEvYZYtR0XR8wOQcPY4vHJMzepoQ+EUUJXHSRIso1ciSicvwjwNK5b863ldtgEfFbxZQxeE2gMQE4OXnakgA7oG0Di1rNmfAYgSaFuZsnMVr/aDQq5fGLB48Xc3Nfgc84ZvcnKwoYVz5tM2Vk9NnLahp/pDTenO66JC3fpfYahutTe+we58yG1IwL8ZYyKBE2C+gfzdCimqGgi5FM+j9ihbMJg0cQ003KCYUIyxyN9DavHcH0iK6+pRDhB+qrSQvzU5VZyRbkOqVqgYKNULVLkw1Uuh60tarf65f0jbYkbRbruN4hBMZ4P2nxY0uWDawWT4IkWcT49sMfmsgbQnvggoMBXQGgQD2eOugyIa5EaJgmOTkdShaRKVqQWUO8D68WCwwEjOVur0QegYE5aaeKuCj5KNJPAUke/aLS8AO1NuUdJu4+bFmDLmdybcxI6Rxjh7JXL4kFo5y7Ep0oyAP+qY3PjhUsQ6BWrKBB148LOJr+Oo6XjITH7VkKrlzZKCBLJuC3gocfGcp60PDoAo/2LTgZUyn0X6PpudpoCfCFgNbguzrzAMx2/dE/m/Bmnl/BDfGXk58c2biTXE/ZYYQ8bFvHIxGHS7YLAjawkefhAchPCrqE1Tq5diMd/sYH31h35F5cOJXId4i+cKn+tAjh20akbA4p//p7TPkw0BiDbFOttqO+DTKnm4TkYmvC/sH0yrdAqfmA+Xv2G3Ae9P/U2qA7G1g6nDMFCAuKtDn44jhdUDaczs7/ded/Agt0fwLXUN9M+1kUKKQA9ys/lVbww0skNvn7RB3t/KDh2RHKuca6JwcamBizuxzZRqmHL16fr8Hi2PjrVVqVtybbK8bWaUHrdazd3la5OSwNmBZIMIoKGs14H0rZWn8SXSRUyjtDrw9acT+pam5OiMpLZTdA4P1Zg4CQd1cx0SP+3JLC9Ujtubf4+60L/5uwBskXrUT+KM8lSE04hIHcQMRDLCoiS3NX9XQ2SOkyFu4kGUHjFVbQiV3nAxFnvp/fUkFS064LDJdhyH/jTNlRBDyYG7rV5bYzhtKxrUATQw56VHGTTtJTUJ9pv8loJIJ+3aLZvkIGwznbk+RJDfrCIFo15T0WulMyt/WeOUT2e2G5Opx2Vm9BaTcrVYuM/EyBXne4CjR53A+Wd/3Zg3Ghg2nnBrWDdYh5bbf1GgoPt2RagpBuMKUKv/iXEgMajFSQ0oWB+AinrqwwvU1OQvxgKHhLsHy6lHQ9Jfm91fl1Ur/ywTk/xIfJRXzv6Jhg+9f8Cc0NC4/danu3HD/A+1w1elyOWW+4yYXlbG7VRWJNnXxkma48Qn6RY7dR4gp3DeDDOXDrvI9JKfbgAgWAjwQ0o3RVOkfRFziLLANM5m9XHDOinu0xV7IAfODUF2ooGalG9pLI2p0ekz51sG4KO+/q+YGN9PBBIddaOWBsIJlIAaORYCTInfxiKmKHyeVR4lG7/rSo/n/2L/3ITXM8QxErGxzM9QyWR/EwOHxgVgxnutzU+j+VVjtpg6DASFoVDLD4An9lAjM8LJsQG6OflHKaLK5mII5P/lWSLpVEXickaPh3Pnw3GoqtrBT61qCPqXziRtqWZWO99F8G+eTrr/FFlpg3qAozvapVaTCM5fM7TO9w35OUKRoci0Kpas7Uwo5e6USFvbBhfdi+FfzP784azvjYiHFLyCgbfn2Cvl40/hS01xUQC6+pz3jqJoi1NGU8l/YZVLo6q9uN6io/eaOkspO9eQi1mUZwzxxk9r60lPN6JVwk42njDmSQ0w18NoLzO1uVyBC8Vb9kBlGWsIZBjxhojcOsMP+3gU5kSeO0iH9F777kPPcLeudyC5/TL+TrgOSMDyt6cQahF3WqN/+4mgjrcRtUzxBfOpDxcyF2MbkYY9GBPBH21Jt8GdATSF04+8G7FpahyxoSlp4HroJIyGppy+jPnVzY7AOdbkQnvASjwYgZmj8UhJPeOY6liydDmGmiHz8caE710Zn3XXEyiYxSQcggrJevDL5WKMtsmVYxuzQik9frc4zrRXni/rt1ZAb9qnG8jrducdtA3EnatFf7/SVgAKKGlXPlxJdV9ReFPPZ+NRyFWjSMcWSPfOM1t5Wrg+g18xf5vpr781ylW3xYw0MGMPw0hdIa+pID92Jl+lQGCQmmQLIc/NN24ld4344IXO447fAptKtn26E1pyNjR6mOx6nsRFSmYyLoAlddsM5ewlrtXPF/0rBkee2JawT4PRCJJ5f8HksIYJpD7d+vuXtE+1/lPoC5lAzkzQEDRPReXACbZrH3uGm+srX7zzVW8kuEHSDARIAjVrRbfTObM8HJB9w3c+VA+8ZVZmp3goCafgLj02ySlDdpc2aSAruENrlIgH4ZVrkFkrm1w1R3FdaAlSqWQKRADJ4mqEzwv+RvIJCRh7CHbfisnH7k/FtLPmj1JHE8OSmVPHJpPldXVmdadQBal2rn+yVLsU8DyjxipnXXcvSAdonSvqQKE44b2pksbe5VrrFc4Hkje7FE2zppfgxgcQdIelQdsIbD2FmdsNbU0O+j9W6VtXt/G4tJmItiTZgqPcpQdHz9YE0Co5Jo9OnrKXX/46goje+PjJEycb0WiGMDt5sCNs68TzE0N9xKf86ySClINqGeViG7iKIPuhrAZdn3bIW96sMFy7lqtPdR531GMIL8FhPzpoD6quEu/BDIBpCD52C2wADRdwOaB4uftsmJOvTrvNVsp7ubIS2nzejY57xAuiG/9jpr4RaFT4rfVJ2Nmbpy6VQ/ibjLbd3fF2nEKyw2+vZnCJtRkL1LxX1GYlpDDy7xcNX1aJBY8ZlTFWK0NnnCxyZqZI8WhfQdeT3rRuOqyy19m2LNlw1UDRYzBcohRk6y81ESyy0EPzLjSa2n2qLWthfKs3OnO3NmOOjx/18s1Z1STtXXxMrZql8NHJsQiQI9tEjLZazTSSmc48pRiJKUIpXsRyCnvNhKtTY9bVJCa/VG22XHIB/xicDrRGGBUPBQUJY2BA6qbodgbgDZ8pNhKaRD581vB3NYBDWnhcdWp0NulX7f6vEf0Wj0jbkxUPaTIuGb7K+uA5WtHsSFZ7efvcVSaU9yLNgtfgLBrXD39n+BY5ZjqhvXvAR6kGw95+CR5GCvqo/ceaJbmBHK295FtxK1YLEyYDG5o3qVwAwWOf7IPg42NFwnbGgJJHI/v4fZjnwbSnc06fr+xOfP1I1f3nZtNkDgUGT9U5Q0IWNOGrVRE6HYk2sfL/YQ+y2fwJCvHZXyfbJePZLTqMWVZBy8/JNHTJlhPCSzmZ0AWpNhhDbo/J3S5NULXuPDdCGXhJ+kWv1y2f2iI2AiU5F72Fbsfhy01lCFVhvjaXwVc9yoz8mhSlYNPRn1hgjnUBTtbvklYGvYvt/4bVd4xzxINc9J6dAy4EjSOduLTmY+c9Rz3Fakaoc8KvXEWVlj2EjyoGk+vcclmWy4mZhFHUKoEj7kCb/BKzTDPgyDmkzvHmhxqUqcD5DcinlnlcR81r9zwkiqGdzlWJl0DNiZSsr8rNo9GSH6vWA6IgGgGRCjY6C3oB9b/idx8dfPPP1gIyTAwrixkPvioNRh+OYofVVh6uj1sNaOE8w3DVuKwQx6T10debFzG2+OpbKUZO9B+eL6O7Zm8bNCfJOLHBgZ8W9yR13iJcs1lYva8ODhDnQgETndOazAuHSmYFbC8KLXwq4PSb1AZIvw/9LCVRVKPzn9ZOwdEWfCXKx3nO1ZNKtwl/VJcEy7INKb+EN6ZcOUJ3WAbYSJmSEkB0xgmeHG2y5QF/godykJXE4CM4zopMZ5gRqQ3oSgKJ24CbeCYbAad7JY3438OPLX8JD5eIsd1QIkRC+fEtb2b9gRcSOAPbnkyURj4yd9K1OCPOq8PWZ8v4eGFI1V2s6hsaTw3nG42BYCRVpLphD/P15MH7PhTxGVAd6VAAJ68xKyjZKgrE6HRJ/HGvi0v/Ax7AswZN6f0rhOrTj/wpj3txIgTsez1xSvJIXORhmtG7frACml84rJeQ+7J2/Blp53egN3t9l2uea1IRyvKQaJp75hzGGDbxyRyNOh73vy8N3bhQyxHTr4TwZa0H/B9mrogFni5/rznxV6/zfzS7reegb0KGWPGRYscAr5afLQfIsfqHDurFP8BLCrFN1sEdRDLHTXe4uj8orpqMf+U3gEMtHhrJK3RrRtlyymq8I/dvxiHzqU+QH4YE4zimeJQzunzzlfBYsiaIyteJVQonsEQ4OEFB/PSEF2ze5AwVEOrDjLKUpmCpC3Dhgt4G5GWXY7bfrVI5M3BWdqtxun2FK/ACsV+genTpJC6zzJMmVGbluqvGhnZCcGPZ5R39B9J2vP3JXJ/fGZpQVVlop0WvqCuU4jhpIOpcKAD25qJDp8spKHaU5gvBrxPccpzMrEMKlWDxnALHtNXfAGWAnDV+hXtaK6uqC/DwHSpqM7VWyr+DRc0+w+Bxeldxvqib5zvT3YiBe+YGrv2ymuV6wzsfo3FQZUL2yvjgKicOfXhL4EWEbw0/LBzBZA9nK7QDnSXaSRGIhJs4iq9ZhpWENN09zrCio2q3hhYIqotDuj7BSYqBh+KYsSZFqlyDL8ytBiQKg75cM27+9/7A7Hfy7aVddpiJnBHgZLAYLZ36eYwpeBgmiTl0hxnt7sxlV7CzGKrlQoUDwJOtMu4hzxMe3LtSlzuTjfRgX+reBsNburYVvdbDLRtf3ZBUDAQBHjXc5fXUeLL1RMbsrF66G6bMSL6j2l+WtNogLGuTiJe5XmUl8NL64h9ju3RNlLlss0eH42bRp/60LA/mlctnl+T1zijBH9B85o6V8yGTWlF55mgfqBG74DWsoflm4el/DEiJH73P4kw7CjQdGqgJ4yKo5ofZfQF9AUeko5Avc2ENAVjo0Ml+OB2xT09fmMMU9SAb8B7y2BsqxFKNI6FJje+xVbMxAt0MrtMNZreEIDdxKBPUViLmET4E7pDkyssHOPUjDFxnPheGTK98b+5LQZWVFqQJh4s+cy2CDBiprp9nLMYUMLNdPOjSiIwEogr4l4Z4UZ4gEuycD3GrCUZTBjw+SrvynUgFS16eD0VexYjQehSdY07L3fAU+S3B9YMEgAJIP7KpchYwSjqkRHs52Xy/flAMYoem6bCkoedellpQIDvwyt1qdgVrKtWaN3lRJhJ2Nb5/Fas5Ppgc0aGBDbpp2uqPDuWQSJCBcbfo53ZCElVs61lzC4y41ipOgFk1XC+2H5VOFzgz5wCyhCGeJkUrVzAY7njjgszgIQm0Rkiaii97N/k9uLXBHH1WA/M86fE7fPGeqK8pUA6nlKPeK/OhtairuFtqS0F0wMyOoRvpRnKZh4CIbqUyOf6lsqEiB+EYyi0hkrZnNxu0VlwZUIMR1ZrAkd+fEmepYaPUVmjRUXSddBIvFcQ5bHKnfvajRuEQwF0zQzZ5SyQFPunXDwfTXNlf/b/KSEDlGUTcPAGA16AAqN/eR23kJqsFcsVDPUjF4y/3ELppS8Pu7jvXHHuRasuIyMW++hBABdrCqwUXiCoHgbeZDAZgL2uOt+anknDhNSegaTMp00xDHXkexarVYysi3SfKnf8JUWH1j/fLqEYsZhlubaZJ2ZSzQPcin8jIeE27zsJQ8bM/0dLNWADHY9cMr+X0gHRmEydOJbEaGfvbtgm1aeaaci8IeKoXlyKHQdqeBtdvRBTVjZ3lxo2JQZZuYL6Wq+7FuAEdZ2QMbwPrVPz61w/MCJPQS3gseVYghCgBqBUK+9FZskvHXYI1V8sqoJOMTQXJegVRCfWKHWiFT+OIZXLlaRLc9WB2diurs0NGK5BnCk6gpzwD0CYREJ6M+vD6xbWTqFHMECug81rfi6aHwV+IZQ+RURmC183V1od8RY1XRZUCuvDP60dG5hGRCggeGcQtGFR9tl/OJHSm6wBX7Iymyq8qO2aBM2lwFCz6+OdJ0xbXyn0Q4GniGjPvYRkMPRGR/uV9TnGGtNy8L/3cDb4Pi1F/WIhE10pOuZ5cxli8KX2Brhd8JMH7Y3o4M98qyrmRdZSeXvHarSbhOQTFdC+cKICtL+zelTdQxHrsTLT2W/AE29uqaqLHdFy6Wn0qsxoT3qTuNMUWdawAvRUSURCyEAQM2CRbijyZJ3McXw9vXuziZB7pdONel7cVGrdJ2HJ1A4O4QKeUr+i6s3Bk4lhzIF4+jBDQcHjLsNIL7Pi5SikUtfNmmAdMYVk51bHwGQqbV051+1bqSRoLQUl+2RViuSjNtFAwOx4QQefaEyk+Cc6QgHV5YV2ThKOay1k+8QcGxswD5waTK22ncFHwuRENp1dPw7Qy8xr+kWBez26emFdXHpnS6f5FcgkuioYpyTB3zh1nrUyMVfaP1PkQy4iM5Y6NTYrN7kiaEI3V5KZnbl28lRjlqhOdBKbCOKfz/yS1ap5Q+l1aGz8a9myiyWdRaTOUerE+3KaoAR/K/4haWZO9sb35QFvRmgdukz27nHQ5QWOUbBn223tfXQjUspsEofP9NaYalo8fdNpqMTdLr+EeCfbh4bwWLloDU7Oy5d+GNSmaCaawvN/3YsqW53oomFwGiGJqMtsj2jwdHNOgZZJpcz76ww7hrYlt62ArULyezl21Ftnwyy6VL25l1Uv/kklxmTO5EB7tw6XM/7Hth7ggeh6ov8Zr5EminsW8CxWEVwesnMZAfSGINSVx2PBXZftfVdttbla5ORU6QifNWyWw2EKAjW96eHIkHqIJhTe6YBDH+oK/zR1QMH4C//80Ud6pHi3ks3y0i/Z31NNAc3IwfYrBR7ExQ86ITOOvXrFhTn4BlBhIhkjAgr8s4MfMulyZbQqKCNnoPLCFqlPiIOECnB3OUY7dYRwBfhlgmz6McF06VrxGF3qiYjSRpn1CuMGZ2l/o1FhGZNupCsQwIi5Ulpz9jfIAzTX2RxrBQDQqZbEyyHa5j0hB1qa9a7pklQzqm91/0I49PVqHqbRdl7ifBxbIcnrKQLtkBM8whWz+XUjKvdadqs3ohyOZQFLbcUvL8qcCsb5EVoAnCo26gJIA1kO4a/sDSatVRkovFWdoy1vhCO5qkWEiEPFcMq0zX3rMT+scc7Qgispu58Sk2UHzJqFGhD8IXn7xk8bi0eDvdBMtv7spbB5F95MJkK3tIxWO8S6DMz2PfIxV9B+tfOhCyCxSFgmKIcX9hfqShHrp9LL7JELOwLCPdDOFk8faq8EJFMgsDKN2RVgultO/wxkSYxVrLIioYiaw5VCstShkGWeH/3wVRthnjbfVH+ffrMk4RBVD+VSZ1ZeTD1Q0KjiObhlnx6IcaKZDWyMGezkUwTKE0qDk6zsNsdxWc73ksVknd4dqapHowq76u5XyPqXOkMHUilwixsmB+kvkk4lZvgrGVf8tTAsgvSrZ0CweVfOLwZKYf0sRHi7poTqyMcTcfRHHirP4bxVdByv2K58fBJIFlhPI9Bj2DaVwI85SGyGGyRPCSd/7kVGvYCpxiS/uvNrXplsi7R3B/wLK5EfVQ5dah8xrnvNe6e4outNxhgPWwyAfSbxavECF9clK5TojD8nHm3x8Q3RHbtlIKnitXQ7FF1Ee1DQn97dE19eh4bavnYPsfK42AUUB3nbuPkiGdplFgZAIeD8b+rHVirCkWh7J4UpHz/HQfMRjTgv+5gf0qwao3RUpKd724ONK8hxDusUbg25r+ohgoBDS+1BWN1gJ2dPgJnUGbnoxDKE2seArV67t483uqY4LBKFTFjHyD3FHYfhfDCt+ZqndvsfGGxniISkm2oaZ0qY1WZ7TSFssbGbjas3QQFdTTGxlHrKKGHzvkhHg++iEtlIAoulX5wHJ5EYIzTr62o2rIE0ylyEL1L44vBBOzgoJ3XB3C3Y6OBT8FmntHD7yzHu/RXS3HrRecmjgzApmhhsDUZUZluJI9zZcWGE1w+BRK7V6wDhzHgRvS/Mtk79NGbKJ5tLzg0uDSa/vjeDw0OA6HzzQV6/mTEYWuhBCLydVlEGehdMvj7qhcUaziPCpwN9CNnre3miiZuN38k1+bxolX5s9Br47GRU4yqGkAkz8e8OHoVv7bqB8PWpg7vSXPEjNr6Z2HvS1SVH5C3a/nqaRsOVOQ1i/XgLHP/3CqCcwcGi5eFBBf40MgT94DI1JwZR2/25oNVkcOknhu+KJq3omugz8ezEsHQWwKOy4e6qZCSTHFwsNO8frjwt+7wt31ytyHEPmh4fW7HW2y1O+dS0qYS1wq//cGtLUyVOjN3rfEaK6cCzalnLoPfsIRYj64JK/AcrSF9mDqhZSMJBnHO0ztnxILqOmwBxh6oAiZq6TGGisXSHUIaI86XFpIKZu02LcF2ln8W6mf3hHXZoGd5QmNKns/0n6sItvbk5ZAX04EOxgOFEeZnmbimFi4/FIJ9kzjKonIca9F69BPqVV5w0dT1NoNHs1Hp18Ql+LyzJSGUMYADbBz+qy/gTWx9bFlZWF+aoQgucB9k1UnJy/AztAO3WAlGdsM2prmbpq7vD/DW6YeDktbFkN4kghOQYTRZoebOKZjAI3IJEG9mckmMB+fXuezv1bppxxtdNBZ/k2aFNePsEMgtKusSBHLIBmB8Lwj9/2U8GmMqWjRDiW4GKiaTnDdDchth7Ca5hjdUsifzgZWVZOtQsZJePxSglcbn9vqVDCQau9zT8ek+1UJvMBAOpUj4aQnG+NCfQL89iVkJFjyF68QpzKkUDEYVNlhFKNeZKqv+Kk8cYrSmmKM19VeZfk5922NVppGqloPsvV8feO3SNOCjMj5Z5jYgrxMCxn81AmId0aKM7Vlgd58kdc0XTV3oUI8BuORF66rFqBcF/CMkGAyDWGMbfaPv/Ui86nmxb3+OGc1ULco6XU8Ph+9oyD5GuNWeVkQkqbmwVvJp8n7XQBPDSCRHKBsM7f9rTYDDTnbEhg316IXRPj/hWEvl5OhgNzCpqRG8uOuQkgNJ7t6YRYX5LaLIF1Lo7Yq1b9OYBnfBvPtW0rE8Zv3CgtfhSLNbl9TplgOZpkwlG5j8YGzE2zmgbmTrr3PV456ID9MEXjAAXxsHs83EGwQArbcCdtAF0IRdMwwvjzmNsuoJer6IlArIBTGbUnly9A200P9pTYLltJc/KNTZlR/jQ40gSuga6KX07E7WL5kclgdKCshkIlIfczk6UuaC6upskdLRj2AHOGAIYQXUdzD55jgNGDAf54249wL+c6HzMfoVmD4zquDOqZ8V7nVkFx6hgBE3XR8U/74lxE3FKQYOFj2F46X/LEPaFXI+djCpGp0Zh2VrkwOU8ke3W+jCflnSXBPwvXrpxLXOTSunbBUj1sP8nlriyiBi1HHiF5R9RiyRQ3B7J5o0nKxcXd/BS52KzfoQID5LE9JZk+3xZsrGE3ZFaDZiAssHjkiNFCCwmPQ1WqOfJEhE9ovQn4xzFPrW6MlVM1g/NFagNYmfTxmYIfrUwLSa+vJP3rlTEkfzL3sm2YyAkp8o1AF3J/Ur35b2VbXjsO6Scg3rOamCN+MdNVDo8lp6hr02oG2mX4LksTZzzhAn8r83S6y5XRufPEm9X49se+4SKlMHh29tGHKyGliAfFhOJi0JP2crmD+soryHUdfHd4Fe+agat1/ZWqSg4IW3URNIVPxrNQ85af+jNq/rvQbbXAdjZmxJxeTBfsYM+MvnriN50Yar9r23GVDfGJpsxBaRrivdk4B3I7tWrOL+idV20+TyqefiaHfOGazHGgEvx5AHpN1SKyX9IRFlw4fxVoU/JFGi25yrvjDKJHgP4ZdoVmiRbItioEbmGaUvFOl2KtdKUbOHUdcxrcJqbJVzBCPHbfFEyUr2nhpQOY8dElCTsU3xdRYENyy4Rxo7hsRJZyCFHYy59bpZVy/xVCqG8jcI+nzi7Q69+DCr9k6XLofrUXVoNXlWvJRNCMcPrsTONwp8p8XZEfjmK2lF8DM04K14q7prlxs2UpQxBtTS3IEujL3A+fAFXf0ARmc+wJGqImZ0I7I5erpz1tlUcV4hRPeWEhKU+OHfMN0qg5SYAZMwfE8TQCefTSRWBOMCxvx1FP8pb1ugxSECXetysb0cZKTILyMEXMQ71vCwZk2NFBZumnncJEb6a6Lvp1kAZCvMj1BuRLKS6gGAGuyaKeH8QrlO5tg70dEnR2IGnZ5a8Z5qg5yY3hxxoQTaEF4rbBsWJRjJzVpS0WXhNy32/eYqaNK4ax29NMJ60wWdbZBCJuYbnz6RIGhwz/UwNwPvEQD1mRZBUbQaXegj2o55dj1K1sgWlfmZUCJcHhVVXB1Zq14jet4+JbICpW0Dg7IM3rxVClcMMzzpv2L0Q/vgsEewk9B7TZwDAjmtQfJC9PiPDq9AXQo0NHYUQ1FGBDS+n6cg+Pq2phq4DMRHHSpvBTxSeDNtr58pNMcmbPSBjUSdhLlgB2t2xCrhZ6chnxWD8ppxqmuo4PuqQBRYMe8sUvoUjyej0ikdr06yCwoWo/z1yJ/2a3+uL/sPFT10/RZUgHVk3z7UFK80wcHJyHi3uoXtYRYR4y9skuAG+lUcpAXSkI5plSa+48jMl6kOlpBmBlD1uo6pq4Z0f3t8Iy1QnsRsAexPXGKM8UcgGMe7lV3WnCeWRIE/+5HowUPVpnfeYcYUbVKND/ev+VHb8txtQaQc0i21At8x/hbRbUnF5RwO9DBzW+icw0gw7CJQWa8auTlsVorfIoerzRaPS9SgM0MwK9sy8vkWa+JUC6f3DAZNDEBX144hIyex0clw0qAKEQnY2Hybna6yvJdNyPll3R+4ITsuSzsImhlQLryLprAMGie/0xe9Y6Dqn0MmxDRentHZ8p1IMSVF9TI0IXDH2hn0nJw5aaPO2hhI5PmOGIAp6QbE9oHMB27i5Ra/ap3eDnYtejMB/QG7CP6B9WsPpr8L1VsFyqPKHGbOtbXgJCpmie6F8It8iZ11x9GL9yzBNQwjcypafPN+eoqLF6RCmTbH2JaqiO+wpTwD/aMTwkDt+s7ynHF4PQA2OoV/wTVMs/1N1kI+87sRpX4e/ARDg1x4w6vQj2CNj6U7jwUwF2aSHqxeiM1lJ0GgQaUlBdhSlYlca5GsXWXM+pwsGyIvLiNBWnvf7xkf+q4S70aqPaSFBDup4D0jZljcWCWEEoPE5PgJt8InrAxL4f1/83cuoNhIKgNthuJ+HCUlaHXKAfMrkJa20Cl5z/xU+L6ZX6hG/K18PwLcH0aJkPEJ+oU8JI23+6QhISKmYraMW5A2nG7tml/lKoA0v5P6leYUTx/L7niHdbS4vINeB1QwVZ+2KGhTYe1sp8voplpIrj831CxKhQEXFgU1E0dhv3VKSXJuQg26PpW4LXhCuDbyH/lI3yxYL7OZZUPJvabQ7z2kAMcX8Wl2UW7K3yleoVOlQhBJ9YI/E6DG1LRWoRDT8nVPkANwYYurcVwrHymJN1Ft/Dxy6s7b7uC9YddPKuMBR0K6fPhYKpB+gXbyGai2Fvs5nS7hEwJeViO51HYljJ3fAwIv/xGwwFVdsbTm0A9AeG4qIYL0FxIO6eugVXGrRNKQguAx9Ux9eyUFedDufWpIfp6WiCDVH733whR7ZBueLS5nK41oLf62oWfGi8TS0zWA4L+kTKsl7Sv9CCF6sk3fKbj03TC488K0nzm/F1gq20jvPUZTQdOpkoAp7aFZ46eDYuJ3LGy1r1dtjXA9shZAxIVtFjXP17edvlkmIUDyVVc09RA0FoETWK2uyFoyw6njJqwCmtETv98+YqKNGrgPg5JgXLNy7AxqbrYocl5iWG1nIiTWwPozcALERSVibhiqzyz5LvOkHxdSNXqfzQ0Je0obxREeJisOlPARFVnCGcPYzvBE7hex2+BVPfG3WS+/Etr2PJCOnEBZVI+CVzZNxbKHth3A83tMJ25ynalb44zthXvRUu6Bqf2ypH2YTIhDwe351ncPw9t6oEWPhf0ayA2VpHCN341jBlvt4bhs02F4xUAneeIwtAOSx00YCkv5aY9juwKyrBd37Tfx79zq0Vaj3Auvs8nDd/EeeOokS01m+6+ZurVgxdah9wZUzN8uce2Eop5EtQTjitBB4c1qdN/phJIn9O2V4xjKKmCHmYwXHyGTRx9QGqA8upBMVlcyX+NPLnSiRxg4+P0A5HG0/2i+3pETp7Eoz2YSaDsayImPhAFnueBhGMAMooOieS3QoSQJvpsafELx76Ao1xnaSgIq5igZ9XKLZ8F19YKSdWgtc6o7E2ahgt6Pcn2Njf9T7+kj6hf2Q/qClQ3vEAUohCQiot0ukrwbAGhI8EmA/WMeM+ScsO3TrWQyzqKydxUb2DFjEj8O1F/epDi9ALVydYydl1mu4uiYwAjZ10nFN645FdaEY2zg/34rbfsKN4jmfl1ZyIvAN913os/RqZJItc4hyK1KI5S8yUaNTttmMsf4gIcHC4FR0oWbn6mSPfRDy9riXiSGLbwAR/mQXQVnliUv76T/JAwLW/gwdTrIao6tN7Etmye1ToCm7qDdaxQtS1bFnT/2ZS9Eufnw2aZhYST4pVglguGO7ZAykyfRk8MA3m/opuXfy5kVVJPu0BjLV4S/PAcPbNFWGj4im63Woq6ygoMTEfJJrOUr0yZ+CK6ZbrBeNI6ZZKxVSCMOfgWh/fVRD8cuO5TWRna+jqC69QxL3Y1SMOA4uh9XrbdoSpYrLOaI98wdnql8g2TDjf1BPPzZw/887vZUzm/6QyFSEvdx6ukg8LU+j3IgAI8g38W74KkH3WbBWs9n2z7Hii5sckMBk0VPU+2/YaOV4/bz97AZ7ZMqPcMnajpG5E14gw7BobTQY1NRbJQQg5fzekscnoNL1h2pt7ZS1I2CzMkaWsJNtyo9M9n4ueHs9e2KEu5r4XnWMzIbI/Yx+FEa1bgZINCuMgn6m/y2eE0Ji4kFJjVc9zTCeFAo0HF8T0c9cJpEBw8XjSgVwrAQjMVaV0FhLjEF8DowMVLaJKCfvayBkl/CZKsvlHSSpN3yOLgwTpJA2iFUiGMMk67FXARUyvOMHLcGJbr3uuNYkSgoFd3t7GfjD+aAzJtV92Du/3G86Wt/6pLWCt+BYZG0YGiZB50vahR5JlJWdEJz/yvmRO23eeUoy5BzDFwDGu+dxmZGVJwOrHLtzDMUhTf+ZqTisKvY8wely1nzm0X2T3/aj/snHh0pZkfFPTzKYafbtXNAcRPNJd2wdzqAhZkyrSbgyc5yOK82c6RZEzzRCXVmMSneV3AV5tYUxqRbhT2kSSASCEZaBtKoanUsr1SCKpCPet6GPeL986loyC0KIQhnhOrDbSh85Ld6CSsHb6+IhY0EOTcGur5K1co0gmmNNu3Cw89D+jXbPT2U20uCJkhKuXHVVAf7kfUvAs3mWU9FWwRgkEfshkbWXvZTQFPSTMuoPJ0NaDWcIMGYXj0N6THJXDvZ0q5SOtkr7nIUcyZL5CJV3NsYMDBBV4EY0/zc0jDAadW9XvoAASKfelpqHTqooq/MRae5IDqXTXfAH0+nkCoYW7tXCJLoUtaF+evacAPYpqvjpkFXZMHz9oAstG/r3xQQQaCPJMshpRAmmGB6EITfWpFmdSmNCNEvwuiim5wFD5JnlWpO2LNbFA9yFDUQpxT0ZwiiUYHaFS/5TZ1zwmB12NODtX8KfSlKRgOsNJm9WaTFIavJsrRVjbCH0TNu1OYP/UjAEeB/Cs21NeOV1ok69gugjXcuzo980ak3v3OSCFuQRR/iqYSBpSDLJP6nOIXuadCt6CcrB1n+AcsTPuQj8Y0PryGW7fHWOM/DY2+F4lr3ZXCUABL77iHVxfLo76zLCKf3L9HO9TTzsNv/EqG5bDyKLBydq9gzNOOilqnkSrT5f7cBLv49mriHOwyhEkbARgNX1ZwQ/6SDzXE2NYf6sNZrT/I1TCIJvL4JSLyi10czDyHeBdpehYVXq1YREh7eH8LvUok1dqSjFpwyhfwdHlb8Qew8ogHk0jPYqMMRW2TRvgIozYjIJFHXLkfP7ATjJjqqZ0tEHnSfbBSSUsiNLj+qNBSHriOdHqOK/a1Fpxw8gXNPdlbpZkRchsIXq9/tNx079iMRiGm0O6ldIyfnE/4qIqy2MHAkozgzZEQoUr43DbMnVMkGHw5k0wJAjdHetyw2bvSazh1LYC6S1vq4asUmCq9Oqo9IY8JeCMnn5vs+WVvRSwjr7AHV4zVa9HPOB0nNleM6d35Yn79/1rm+pgcK+ysN252Dud27qWDF2p2bTHfjk5WhoHmh9JZ0SXwFttSNSRWKM+7InmcxXeOSytYtbwKIrpKllnVi1MIkgDDxJvH/yv5jVGRfwtdW4tm4hgYVz7veIoN1quGGUTrvyt18QoqG8PpQ6U9uVu8pMatEP7fuP5zrfDUtuVy2wtSH/4H2l7+6nsoOl62oz3UFEq1Kgv/bzwXUe+7hRqNSct6BFrJjjRmUvEoo2LqAHZpujacc+BDmuV+r2llhp2mB5qKXgcMtBu8SC9STw41U0sQ4eVwEZlc6W2BkNErhwvDZg3aPFdDCvPMpRXg79FemC7536JPfymrHRNZJKnDEDEIj0TrziB3Cym7xTW1MDwPl8Y1j+ObamuKp77HUY/oayXAuEd+21LdX6QHr/VQd/3xdp7quyYn4wgMVsT9TC08au0wWcOW9K1dIi41JMlx+pYyNhLUOpmJ7ZB2n4tYnqTHouXqkCyRRpMBhlkVbd43ZV5vPMOpSjM/JT/zkyQBhQM5zzYqmIXOZayQXBv9gIKE9EHpc2WggENc4n8U0v0mF7LYLyF3lNyWsdpkKvUyMdgC9BEDtTLG5p0XK4g8o0iBa1zTEJO5SKBTyeMjB7y1yVpZ2jnXnLjoB+bZl5AgNL0f2rECJd2EYQVTypTdyM8CVsqK5JjtP0TtY6s3GqgqyXKF+cBXCghthDzvmTW+7vtsuTaj/2u/4/hVDWIX8r5XMTdzEwLYJRY3WWmCr7IgI0kJO9gHAisUAy4RIwYFVUWeMLrZh9G27YMALQ46vZji0CfyeclPAgaYHvX8sVPGa37CwO5g7qeHHqa+3ZIa/o5PrY1qb1b1ya99jLhEjKSWd44s9x442hapDhX38w8Un+vNNGdew+rGKlXfZjT6AI38pCkwAylZOwI+7bbwfeBNGWg4JV7B5mNjwuhr3+dRnWtYqolTe3oh70UgNDyQ7ti7LvI86AK+ESFsiVnTkB/3vacZvdawkQGb0Z7HrHDZ9hYDvcraKp0c9ng4ZAbJWVKtBjMzva6ad8tlfRkU2P1V2qlUWar23nqqLOM4EU12LtJ/RshpOVcLTHnRhAy03WCmLFsga4naPqQIFJgDeV0ni1+M6jRc2Ty1u9f6CYG/+3YvDxM05GhavUjbwafQNMUY+Ow6nxIb28VqSwqtgNBhCJVC1mtkN4S/Jsiwoh3iu3OIClx6W1kyyBzqVNUh0cgFO90X+iGyWZ/ZPCIh+XapX2Nkr9Lz0pRbINy7ksPv33WaLizummDPdePVRpg4CWVxly5UZHj5SyXMxRqwoEKsfaNlcxP4ay2SMUnnzbliv+7Bw8AlvNd0AR6iiKACQFwlcuZWSFwGSx25qjuQVI2QKvTxyN+acZhcEdhpKu3sRywlceohcjSNuWEDz30z91S33SKIaRARL5RKTFcJlMZ2GUWnJ4sTedzy/urp3JOafLyfnt5KFJDfDkDeK6m17jRd8oI4YURtZfGf6mlwHme0wklm3WVB/4xyH6er1lpU1wH3DcHhVheSV4I++7d2dA+gujQbKBTxWeOtx7WwYZrSt7dpavsH2mZtc7OjyZGSaEUPbLlVUC9IkWkvLQuRUJE+uhKpSX80uvoVX8oA6ecKDCBUUMsuIKrzKwLoNOVTXnP+P0mmZj7SyhGZsOnFA42weejgLStcpk3wAoZ6YrhHptm4bwu1R/U9vYHvsN4EkfGbhr4Du9BLZGBys1At0gVUdfSBoVfhAwdQmZlj+unotJq5s44+fUduePMtXBqcxcBELjcH8u6Tp8jlcY8aVQuoY1tRvJPL/p3UJZHthDRFCkbvXN5+k4qRVZVMJYp/y0MCgZ42hBSaOjacbUt0kPQUHiyTAHt6R6SAwgX2eHSQl1B3+fScgqob1Pm8vIyADMk2iLktktg2jSmX+LTOUkcuEut6RDxk1eamD2nQp0OSOMO1J8NdbPmmaeWIVjmJy3ooaTx3iZqt48/7T8cuKJ5ps4gynZk28GlCR5lkbhUmti+giVRV6Wk1XKnIp731E95VcsBt2dFnuK6PYRNqiJlD/3kJ3odR6Hi+yoS+jpIcfZIuiz5jyvWD2eowjwBxWXJRkEivst1Tm9gKMufgD5a+h4PKCWGpZ3Wbw+fUR7CiLGTxHWWAAAhDb8BVqWdyBNXF5LYt6K8N+FngS63jCELhEmOaGcy0kMr2FgrPh3GyMtZWBmDCi8P30BMz8EdkJZeONwZFEkclSR2LoeI1IvXK1raxhugMSqzYeQCWCfYSmcksDOE81LU/JgNC1FL056IrN58aUthUgcLIQmuzfLAVR5zmeZf1bfAqO5D8QW1N+T9Qh5Fcqzr3MEsL0u2BMx0p8HWEMpd8j4qQUYpsH2FG4zHShnLdv8OvBeHTulFA+s9fYd1TjItA/YDuFU0GaqgAJKOrq4KozGUZkit76ay2Vdxg0vjv4/Unl+UstziZK3lKQEipmrku/zCyAI1QT++q1wbkBDij/bD2oZ3xzVsYWbd7pVW34xwdy08hPSBOgSjNUQzya2IQBcoQBtw7TzQ2275TP9KW927o9I2bS7Oy6onAQm/RCY8vBI/6BLvNhMxo4uKZwNGf7rVJOdjbedgJjzukiwInE58DLyyqgPRmd8b4SCZhGEf180MbdOWQ+sYzQoeSzC7b3RtDjF/MwGdZO8AQjAbJSm2GHgDxprVBlnkamAfTd9sVGcZsZsLkOc6OUorOA+LViZykBRZFgwK2BVx10lTUwNtjBbndzooimLBzsIq/a7oHHB52kHgBdd98zx66Zt1e3A2W8hjruYlDroTMY2EBS26Uir0XUordlcXU9tnbf66GTLPaLSscNnRFO/aiKU+t7oOfQ8z2591xj2EHd3/pxQY6bGwjcQ3RobNF1Q6q9j4Xw/10ZC/4QTUrxRwYWI0k4JnPfTjgetElDoTzn8UFVJqCByReydpGgsh+NSfceWkuLguT8wJsaERITcsGtxZajT4iaoy1BqNRpJJu6+OTANX7yZXGkzZAvwFTvvapxrJkd2gMVwdRpSKG/uxnT/hve0vzV9XzChCA/x29dgaFAwIh57gmubYfth1itfq6DbJYelqjDMDs9LoiNpoBYlW8hBtXWov6b20EoRZYgx5Z9QTEeZCamDrZyJZ/DG95xWhz4BNaSXTOQxxmMWXXBpqF8nKT1qHY5BPfN84R1YZRLWTu4WZ1miNpAzGhkpKLA8zTC1A6LaF5yIAJoCJlfpeHjwIAfGiLGsFvx9zbmHA2WYefIWkR6UoEbJrPivWMbuddlN0qAvCzVPZOlaTAz/btPppKYGInmmtczDfOmbpFo3ZWl6EdSmz/YHKDYFAOllu4mG4v01WtlDpD4SGauJxOf4JgjvNLjHz3t9QdHd6r+e4xTBrwHlNRVIecTgfj9AuD2CCMZlNDX0mjk81rEy9i4E+aTxwjLD19Jovncpx/w5JFdAsjvqVeJNzx7ZBG8Azlj3aBim5g/nKTW/2ideUMNbvYv+WhRLcCJiSnqgs/6w/oPgQyffH/oRsZW8yV53guNcW3+lfyUT5QRrck0Q+k2Y1Afn60suQ1Vjb1w4i4/uGkx9VL4B+f7U9e8RaZya7ik2ztL2RGa2i2y/6h8CR4D/VvQHNx55gAOY1l+5n+ORhGRWkxsFzXEEgG/I30eXUBgqekUdhCbwlUMCpKnagLWfYvcP/YJp4A+qU2Tsr/PWBWW/TZ2E0A+20xH+98K0BWpOrT+1sWGy3CCalc59XLN0PAmB6Md9jJ8iTSRjItQvpJgobXK4wZIHv9QG+W0SRggdw3N2VxshySKNENtjEj7L9GeFCD9n4X32DoigYa3Tx/kKzvbynY2nnu3Ousn8smt2EtEaFc1Wra8JrKBQ3ClT0FVprHyiRJXESoR8d4nIWL5oJWpvE8WGsS/2tktVqMCQDcx+p488cD3qGmJS/P7GeGW8LAlWC5ACP0lPBiHl0R0lI5xxhKOhnLne+EDL+NnNSQ6IpFAZ6Nbm4GQIH2ee2E8Jn33db7dYHeINvdLTTnGe4EwhByVs3RwTuklqKYI6Li/ibupDszwMPrGaEvyl/d2sm0duOrOjrUa0QgMCsI4NIqXqIoB8Cn5AYgD8dkeA/XPob/Pi0/6bTYf0CiPtzT1Wq+QYk4obd6eaKBiNM506AWjGNnPndUmfji277Gf/3K5R/z5uCj8eLV8dqFfesbsvZRmYmYlXWRB4bxAAqef9MXmGt95q5+WOARW0/Ng7+75qtuBn83ZgEd9Y9cxgxIWqrIY0XVUYBVr2V0ybbm2pbMJ+DQz9XGNboLe1VIp6irO9I+UY9xN3OfkA7bjNLiNM5XmXSKl3V91rnvm/TyZYae5u/XwWMiicrDUf+g7scOUfLtlt2qrs3mshjm2w2TymkLS+DzDtxT3ACEBdRV6mZeFM9p1VnBGWvT4aJyJe1nsEgcPN6khWeYFasBSJbLR28PsTpVa9tl7yWuBTLfvkpVX/HmWi6ezx9h/qz5Cos0yyCp71ieVdTkgBcPNu5X3F6fxyT4gF69ulPOVSu/6E0XrtEixwHHnoxlai+AQBTIvrw8j23mh5fNt+dYdONXsIrsQvHRbcFQQdI4oJBqBnPLPL9+096GAh8PEhmSs/OTWUphkKNelDkS+AOsWCTBe0NUOXQuS0dkzrKioTFJYbttlIZsud5pgacTD99HvfN0numAo/mewDdFx+rJdpy+ZcxWhwP5rVjuHq59h7AaskX3nKM/T2nzVwvq/qQUFFQVgkH2CqL2THLpzySRRqu47DQObsBfqDRnoa5/17xCWY7EW8/J9yUY1+misBDQFmiupnQSGA4cOL9rdCj73NedVtBJkv6JMz2MTYMQftpMm/XhKqJ6nmK2C89S6sat9rKfjrLKWD5LU+7ryG3dRuKIX7ihciv6bKSS6yCxzlX4NiFb3EVVDhxlvdvlN9zplmIo7ErThdXINd4xhSoSZnsoQbpWmNnnmIqvdeNvEIE4NY4rGrOrZilwcWQrFEKpzrxsePR9DiIR2kARG8wYk6zOeRxn0WxsxOLhXlEDGDv3sRdj1duNWiWZzAlfHNPUHd9J8aLHu0rsgMZtM6sb04hC1Al/sicL1Pw3WL0l1JjqTdPz/yeN8btJ9ldezBTjEjm/82rDmeva4Zwj+XMBE1IUSgzNt82emt/M5tfiHkOvwq1XP+Nj2o1b/rHVEdOzLQUXXHytj+9z0xs5/TwEs5rRLoAUOF2MXtlkCLnnxWGa2Ftbo3h6GlfxkHe4JFE70pKNSu+Yyg0r6t1/Qx3MtMWqIemjrWXViy92jgrKwDSpLLDpJhIDpxdEumthROXPDRD/NmUBfp5lHV4Z1R1rCaWsxELCO2dFpCAVVux5srgv9E+VVOu0pQKINVTt1rkg505h4iH8RbEtWZZTHX4pNeos9pNvSnYv+lGNDJ18n/K5zG9egH7stRnwL1nkaPzUtvm+AZkg6LB6hpl+iHnVvX5DarIlg9lGWEvg73eDxfQ9I4vNuN7XNVy2GTAwXhAUmsnJLEFetoB3dVxsgX2dZVzcjhd7qM+GzOmJOsJGdH8UyuLrUDoIpsEbI+iouyY+07ZahQ7eZusm1CwObJ69G9Fb/a8Lxe4i78G+pYaklcePKjpalvoVL5SNu7X4IDxIWB4ljZH/jASowhKphIgEfbLkk/NQ0c6rq9XD/y5dHsMbgKZ4ZPLXrDeuBfIZuj7bCs9SKjN4npEcBuagsEg3nK9fkpfvNntZ2LlJXPDy506KcwwQ0I0NNCNSpSC4K/PJODqLSf5+PsWpnGhevqv/hx+qhbfV5ob407E4Tlb8/rusYhRe1Gwh4bMB7dExaFMxqasfZkJ9H1Us93jc0bmdp83MNg1hntTWpf1J9ErDdUOiA+dF7A3HlXrbLVxIKjP8GjClpa4y+5QJNEOJZv+a1wraA7KLnCn4/X+Q+zpoSUmzOjCq91A3kfDoPJDEHeYEGwJEd1Gowfx343WTR9QIER1wnRog8P5Kaqlpj5UAoNLz9BxbLiiq8KJdd4eos7xhVSzeL7yGpqeDTu5RkHpnrX+L1HRYU5kfwUojlVAAVHPtscCArtY4D1qkhXNc4wVjRyRWDYqPE8Qw8OfRdvPbKlcuzDqnD5OoQFvDQ83MVkr0BYUXxJJUsRO1Hc1baUyjgB9oE+UzLvTWf3qcEJMrcuAhZR7F3mjGLi5a9QjDiTRWDwtzelfXhoSAITRPRZIX05qTKzvAauZEAyNiwmRTrLTujZkG2+8z+R0DIorb2Tig/uC5WUdJf+bprK/mFA7n+YOWi+sJyUo8knGIzSjEDPfF1oIwXVN7Xj66SHUnDYvNd6b2fZrfVyMPmJwUTRJRxfSl8pEKhWfyB+aue+BGigVCaSrvDR1mrHiPm7gldRZKiWQ566hswIxRqhFXcbxpiCVHnHK2daU4Ws7oPBh/1PAT1IJXQ6N98vzi7MZoInPmd7Z8mjzdYdtuUcXic4OyJ9XagvaU5QzTLHaSkEz8OnB2Mpm6oSKLxS+meQtHO29BQ0/hTRWBeirLaOHcDawws5fCoQhH1kdJO/kyzMofpsIW9fRWDRJJsPihCjXQOIQfV3DQ56Ep7xdWCDb7grvlin5XMOUXPPNmcdbg5Flv3KN01fWS40eqHQFlsomP5AoXctiqymm08TJIs6qvV2yCQVowzh0BBehFuPI/t5K/W3e6SLbRMFUE2qYZtv2Altu3KsaXyxKad6KUNw2F4UNh9t4wsjiJMUj5bqHZAsBpRUMD0oTmzymi4A6dBfQEVzffuym+/KeF3PkRdUFQX/AazxwW+c7LNSafHgS9oBchunLwiCPBFVTb5YQrBplnpR5msq9XrhzemKJAUuQ4aRcErySfB+ftEKtw0V/t5HiTNP5/YUJCjgo/6UKiVVoiL7gdORlB/gW0itZHxDlesdXxW94mNuKCDD1OhsF0ITCzuX3CecwPWGXlVj9c29bOnsbQU5QSNk4XoJ58E+uxrP9iCR1ZjQQq4EYA5/sy/tfdNVlyZsM1LgYZG4DsMZ8C2xrHeKBVGxhAY/O4zFOD47gZ+1MKu6Z8TIkHM/EUUQpYP7FCTceN2BKCxp0mxFYaK1TF2aKK5SxK3PagxA5ZZ9rOTdIT4x/QFzGNYjTDpxqA/dYLdHYdVWKZ3Vz48h+kXa8clwK3OVF8/sr8d9BOpkrvPZYb6Xv5takbtj384xU/GwIv/9pEAzjii+ta1x510yAHnMd9B/e+7Oe/LRlDvyJNiF5X9HqCVah0j8Q8v8iwTKOA67UcBYnXkMc0j1SKMzNVdDqp6O1JfXQ/2d/3KoPjP4XGuK+QSZsHS35fg/1OoIDDV7OMQfq8lODq56YR/m1saHyvWgK8oszsWbC8qkCW1hsiz/+yoHha2yVIfLl9lzrGF9KwpzWSbqabVj5UQWgx2L2U5JBLk0DXJPbcFVmnxHmbGHvwXD1nEbM3aMzVi6/0EaiTiwblBlVVTh6fx1M7jyCuY30USBaTKU9bfIJP4flY1";
	decryptString(A,"783e007c783e007c","Password",oFile);}

void writePsize(string oFile)
	{const string P="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";
	decryptString(P,"783e007c783e007c","Password",oFile);}


