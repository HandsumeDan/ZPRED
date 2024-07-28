# include "zpredUtils.h"


using namespace std;
using namespace boost::numeric::odeint;


int main(int argc,char *argv[])
{	//string salt="783e007c783e007c";
	//string passWord="Password";
	string inFile,passWord,salt,tmp;
	bool VERBOSE=false,MULTI_DISPLAY=false;// boolean for controlling display updates beyond 150 threads
	int ch;
	while((ch=getopt(argc,argv,"i:v"))!=EOF)
		{switch(ch)
			{case 'i':
				// ZPRED Input File
				inFile=optarg;
				break;
			case 'v':
				// Turn Comments on or off
				VERBOSE=true;
				break;
			default:
				exit(EXIT_FAILURE);}}

	string aPqrFldr="/home/jenbunnyaqua/Desktop/ZPRED/forKen/MP654/";
	string proFile="/home/jenbunnyaqua/Desktop/ZPRED/forKen/MP654.pqr";
	extract_models_from_pqr_file(proFile,aPqrFldr);
	exit(1);

	// Get Current Working Directory
	string Fldr=get_current_dir_name(); Fldr+="/";

	zInput Z=read_zpred_input_file(inFile);
	// Check should be performed by GUI
	//cout<<display_zpred_input(Z)<<endl;
	
	string outputFldr=Fldr+Z.outputTitle+"/"; if(!directory_exist(outputFldr)){make_folder(outputFldr);}				// Main (ZPRED) Output Folder
	Fldr=outputFldr;
	string fFldr=Fldr+"Files/"; if(!directory_exist(fFldr)){make_folder(fFldr);}									// Formatted Data Folders Folder (Also Holds ZPRED Configuration File)	
	string zpredFldr=fFldr+"Z/"; if(!directory_exist(zpredFldr)){make_folder(zpredFldr);}							// ZPRED Output Folder
	string logFldr=fFldr+"LOG/"; if(!directory_exist(logFldr)){make_folder(logFldr);}								// Program Stdout/Stderr Log Folder	
	string dxFldr=fFldr+"DX/"; if(!directory_exist(dxFldr)){make_folder(dxFldr);}									// APBS Output Folder (in DX Matrix Format)
	string hydroFldr=fFldr+"HP/"; if(!directory_exist(hydroFldr)){make_folder(hydroFldr);}							// HYDROPRO Output Folder (HP stands for hydrodynamic properties)
	string inFldr=fFldr+"IN/"; if(!directory_exist(inFldr)){make_folder(inFldr);}									// APBS Input Command File Folder
	string propkaFldr=fFldr+"PROPKA/"; if(!directory_exist(propkaFldr)){make_folder(propkaFldr);}						// PROPKA Output Folder
	string pqrFldr=fFldr+"PQR/"; if(!directory_exist(pqrFldr)){make_folder(pqrFldr);}								// PQR Formatted PDB Files Folder
	string surfFldr=fFldr+"S/"; if(!directory_exist(surfFldr)){make_folder(surfFldr);}								// MSMS Surface Output Folder
	string surfFldrContents[6]={"CGO/","colorCGO/","CSV/","fromMSMS/","SAS/","fromMV/"};
	//
	string pqr2pdbFile=Fldr+"saveAsPdb.py"; write_saveAsPdb(pqr2pdbFile);											// PyMol PQR2PDB Converter
	string apbsInputgen=fFldr+"apbsInputgen.py"; writeAPBSInputGen(apbsInputgen);									// APBS Inputgen Python Script
	string psize=fFldr+"psize.py"; writePsize(psize);															// APBS psize.py supports apbsInputgen.py
	string theHullRad=fFldr+"editedHullRad.py"; writeHullRad(theHullRad);											// HullRad (Grisham Edit)

	//cout<<"|"<<get_containing_folder(theHullRad)<<"|"<<get_file_name(theHullRad)<<"|"<<endl;

	//string proFile="/home/jenbunnyaqua/ZPRED_II/Test/Files/PROPKA/6LYZ_1/pH7.00_1.propka";
	//cout<<readPROPKAFile(proFile)<<endl;

	//string pqFile="/home/jenbunnyaqua/ZPRED_II/fromKen/PQR/3.pqr";
	//string pdFile="/home/jenbunnyaqua/ZPRED_II/fromKen/PQR/3.pdb";
	//convert_pqr_to_pdb(pqFile,pdFile);
	//exit(1);

	//double dryRadius=16.507124680407;
	//double solvatedRadius=18.772025948976;
	//double zetaPot=0.037844942737; //0.018709780875;  //0.068863581929;
	//double kVal=30.359499874551; //9.544989291311;
	//double eVal=78.18514375; //77.28334375;
	//double vVal=0.000890166282;
	//double protConc=3;
	//double protMW=14316.9989;
	//double alpha=calc_collision_efficiency(dryRadius,solvatedRadius,zetaPot,eVal,kVal,298.15);
	//double coagTime=calc_coagulation_time(protConc,protMW,vVal,298.15);
	//cout<<alpha<<"|"<<coagTime<<endl;

	//cout<<cnvrtNumToStrng(1.2345e-20,BUFFER_SIZE)<<endl;
	//exit(1);

/*	try
		{
		perform_ZPRED_computation(Z);
		}
	catch (const exception& e) 
		{
		// print the exception
		cerr<<"Exception:\n"<<e.what()<<endl;
		}
*/

	perform_ZPRED_computation(Z);

	ofstream fOut;
	// After all ZPRED computations complete, store in user-specified folder (outputTitle parameter in input File)
	string zpredOutputFile=outputFldr+Z.outputTitle+".txt";
	// Read List of ZPRED Output Files
	string zOutListFile=fFldr+"zpredOutFilesList.txt",zFile="";	
	ifstream fIn;
	fIn.open(zOutListFile.c_str());
	if(fIn.fail()){cerr<<"ERROR in ZPRED.exe!\nZPRED Output File List file could not be opened.\n"<<zOutListFile<<endl; exit(EXIT_FAILURE);}
	int Sz=15000;
	char Val[Sz];
	zOutput zRes;
	// Collect and Average Zeta Potentials of the ensemble of structures for each solution condition
	int* fCounter=new int[Z.numSolCond];
	string *debyeLength=new string[Z.numSolCond];
	double *molecularWeight=new double[Z.numFiles];
	double** srVal=new double*[Z.numSolCond];	
	double** prVal=new double*[Z.numSolCond];	
	double** zpVal=new double*[Z.numSolCond];
	double** qVal=new double*[Z.numSolCond];
	double** kMobVal=new double*[Z.numSolCond];	
	double** hMobVal=new double*[Z.numSolCond];
	double** semMobVal=new double*[Z.numSolCond];
	double** difVal=new double*[Z.numSolCond];
	double*** electricPotential=new double**[Z.numSolCond];
	double* Sum;
	string*** position=new string**[Z.numSolCond];
	int** numProfilePnts=new int*[Z.numSolCond];
	string* sER=new string[Z.numSolCond];
	string* sDens=new string[Z.numSolCond];
	string* sVisc=new string[Z.numSolCond];
	string* sXsp=new string[Z.numSolCond];
	string* fNm=new string[Z.numFiles];
	for(int i=0;i<Z.numSolCond;i++)
		{fCounter[i]=0;
		debyeLength[i]="";
		srVal[i]=new double[Z.numFiles];
		prVal[i]=new double[Z.numFiles];
		zpVal[i]=new double[Z.numFiles];
		qVal[i]=new double[Z.numFiles];
		semMobVal[i]=new double[Z.numFiles];
		hMobVal[i]=new double[Z.numFiles];
		kMobVal[i]=new double[Z.numFiles];
		difVal[i]=new double[Z.numFiles];
		electricPotential[i]=new double*[Z.numFiles];
		position[i]=new string*[Z.numFiles];
		numProfilePnts[i]=new int[Z.numFiles];}
	int tabIndex;
	string s;
	fIn.getline(Val,Sz);
	while(!fIn.eof())
		{// Read Single ZPRED Output Results
		zFile=Val;
		zRes=read_zpred_output_file(zFile);
		s=zRes.solCondNum;tabIndex=atoi(s.c_str());
		debyeLength[tabIndex]=zRes.solution.debyeLength;
		molecularWeight[fCounter[tabIndex]]=zRes.molecularWeight;
		srVal[tabIndex][fCounter[tabIndex]]=zRes.solvatedRadius;
		prVal[tabIndex][fCounter[tabIndex]]=zRes.proteinRadius;
		zpVal[tabIndex][fCounter[tabIndex]]=zRes.zetaPotential;
		qVal[tabIndex][fCounter[tabIndex]]=zRes.charge;
		semMobVal[tabIndex][fCounter[tabIndex]]=zRes.semMobility;
		hMobVal[tabIndex][fCounter[tabIndex]]=zRes.henryMobility;
		kMobVal[tabIndex][fCounter[tabIndex]]=zRes.kuwabaraMobility;
		difVal[tabIndex][fCounter[tabIndex]]=zRes.diffusivity;
		sER[tabIndex]=zRes.solution.dielectric;
		sDens[tabIndex]=zRes.solution.density;
		sVisc[tabIndex]=zRes.solution.viscosity;
		sXsp[tabIndex]=zRes.Xsp;
		if(zRes.EP.numPnts>0)
			{// Electric Potential Profile Generated
			numProfilePnts[tabIndex][fCounter[tabIndex]]=zRes.EP.numPnts;
			electricPotential[tabIndex][fCounter[tabIndex]]=new double[zRes.EP.numPnts];
			position[tabIndex][fCounter[tabIndex]]=new string[zRes.EP.numPnts];			
			for(int i=0;i<zRes.EP.numPnts;i++)
				{position[tabIndex][fCounter[tabIndex]][i]=cnvrtNumToStrng(zRes.EP.position[i],DISTANCE_SIG_FIGS);
				electricPotential[tabIndex][fCounter[tabIndex]][i]=zRes.EP.potential[i];}
			}
		else
			{numProfilePnts[tabIndex][fCounter[tabIndex]]=0;}
		fNm[fCounter[tabIndex]]=zRes.proteinName;
		fCounter[tabIndex]++;
		fIn.getline(Val,Sz);}
	fIn.close();

	string Output="";
	double avgSR,stdSR,avgPR,stdPR,avgZP,stdZP,avgQ,stdQ,avgkMOB,stdkMOB,avghMOB,stdhMOB,avgsemMOB,stdsemMOB,avgDIF,stdDIF,avgPot,stdPot;
	double collisionEfficiency,coagulationTime,kValue,erValue,vValue,avgMW;
	double* dArr;
	// 9|9|8|9|6|8|8|
	// Determine header Width
	int lnPR=0,lnSR=0,lnZP=0,lnQ=0,lnMOB=0,lnDIF=0;
	string tmp2,tmp3,tmp4,tmp5,tmp6;
	for(int j=0;j<Z.numFiles;j++)
		{//fNm=cnvrtNumToStrng(j,0);
		//Output+="file #"+fNm+repeatString(" ",3-fNm.length())+"|";
		tmp=formatNumberString(cnvrtNumToStrng(prVal[0][j],BUFFER_SIZE)); if(tmp.length()>lnPR){lnPR=tmp.length();}						// Protein Radius
		tmp2=formatNumberString(cnvrtNumToStrng(srVal[0][j],BUFFER_SIZE)); if(tmp2.length()>lnSR){lnSR=tmp2.length();}					// Solvated Radius
		tmp3=formatNumberString(cnvrtNumToStrng(zpVal[0][j],BUFFER_SIZE)); if(tmp3.length()>lnZP){lnZP=tmp3.length();}					// Zeta Potential
		tmp4=formatNumberString(cnvrtNumToStrng(qVal[0][j],BUFFER_SIZE)); if(tmp4.length()>lnQ){lnQ=BUFFER_SIZE; lnQ+=2;}					// Net Charge Valence
		tmp5=formatNumberString(cnvrtNumToStrng(semMobVal[0][j],BUFFER_SIZE)); if(tmp5.length()>lnMOB){lnMOB=tmp5.length();}				// Electrophoretic Mobility
		tmp6=formatNumberString(cnvrtNumToStrng(difVal[0][j],BUFFER_SIZE)); if(tmp6.length()>lnDIF){lnDIF=tmp6.length();}					// Diffusivity
		}

	int theLineWidth=131;
	avgMW=calc_average(molecularWeight,Z.numFiles);
	for(int i=0;i<Z.numSolCond;i++)
		{Output+="Solution Condition #"+cnvrtNumToStrng(i,0)+"\n";
		Output+="pH "+Z.pH[i]+" ";
		Output+=formatNumberString(cnvrtNumToStrng(Z.temperature[i],BUFFER_SIZE))+" K\n";
		tmp="";
		for(int a=0;a<Z.numSolvent[i];a++)
			{tmp+=formatNumberString(cnvrtNumToStrng(Z.solventConc[i][a]*100,BUFFER_SIZE));//+"mol% "+Z.solvent[i][a]+"\n";
			tmp2=Z.solventConcType[i];
			if(tmp2.compare("massFraction")==0){tmp+=" mass% ";}
			else if(tmp2.compare("moleFraction")==0){tmp+=" mol% ";}
			else if(tmp2.compare("volumeFraction")==0){tmp+=" vol% ";}
			tmp+=Z.solvent[i][a]+"\n";
			}
		Output+=tmp;
		tmp="";
		for(int a=0;a<Z.numSolute[i];a++){tmp+=formatNumberString(cnvrtNumToStrng(Z.soluteConc[i][a],BUFFER_SIZE))+" M "+Z.solute[i][a]+"\n";}
		Output+=tmp;
		Output+="Relative Dielectric: "+sER[i]+"\n";
		Output+="Density: "+sDens[i]+" kg/L\n";
		Output+="Viscosity: "+sVisc[i]+" Ns/m^2\n";
		Output+="Hydration Layer Thickness: "+sXsp[i]+" A\n";
		Output+="Debye Length: "+debyeLength[i]+" A\n";
		Output+="Protein Concentration: "+formatNumberString(cnvrtNumToStrng(Z.proteinConc[i],BUFFER_SIZE))+" g/L\n";
		Output+="Protein Molecular Weight: "+formatNumberString(cnvrtNumToStrng(avgMW,BUFFER_SIZE))+" g/mol\n";
		Output+=repStr("-",theLineWidth)+"|\n";
		Output+=repStr("|",10)+repStr("-",lnPR)+"|"+repStr("-",lnSR)+"|"+repStr("-",lnZP)+"|"+repStr("-",lnQ)+"|SEM"+repStr("-",lnMOB-3)+"|Henry"+repStr("-",lnMOB-5)+"|Kuwabara"+repStr("-",lnMOB-8)+"|"+repStr("-",lnDIF)+"|\n";
		Output+=repStr("|",10)+"Anhydrous"+repStr("-",lnPR-9)+"|Solvated"+repStr("-",lnSR-8)+"|Zeta"+repStr("-",lnZP-4)+"|"+repStr("-",lnQ)+"|Electro."+repStr("-",lnMOB-8)+"|Electro."+repStr("-",lnMOB-8)+"|Electro."+repStr("-",lnMOB-8)+"|"+repStr("-",lnDIF)+"|\n";
		Output+=repStr("|",10)+"Radius"+repStr("-",lnPR-6)+"|Radius"+repStr("-",lnSR-6)+"|Potential"+repStr("-",lnZP-9)+"|Q"+repStr("-",lnQ-1)+"|Mobility"+repStr("-",lnMOB-8)+"|Mobility"+repStr("-",lnMOB-8)+"|Mobility"+repStr("-",lnMOB-8)+"|Diffusivity"+repStr("-",lnDIF-11)+"|\n";
		Output+=repStr("|",10)+"[A]"+repStr("-",lnPR-3)+"|[A]"+repStr("-",lnSR-3)+"|[Volts]"+repStr("-",lnZP-7)+"|"+repStr("-",lnQ)+"|[umcm/Vs]"+repStr("-",lnMOB-9)+"|[umcm/Vs]"+repStr("-",lnMOB-9)+"|[umcm/Vs]"+repStr("-",lnMOB-9)+"|[m^2/s]"+repStr("-",lnDIF-7)+"|\n";
		for(int j=0;j<Z.numFiles;j++)
			{tmp=fNm[j]; Output+=tmp+repeatString(" ",9-tmp.length())+"|";														// File Name
			tmp=formatNumberString(cnvrtNumToStrng(prVal[i][j],BUFFER_SIZE)); Output+=tmp+repeatString(" ",lnPR-tmp.length())+"|";			// Protein Radius
			tmp=formatNumberString(cnvrtNumToStrng(srVal[i][j],BUFFER_SIZE)); Output+=tmp+repeatString(" ",lnSR-tmp.length())+"|";			// Solvated Radius
			tmp=formatNumberString(cnvrtNumToStrng(zpVal[i][j],BUFFER_SIZE)); Output+=tmp+repeatString(" ",lnZP-tmp.length())+"|";			// Zeta Potential
			tmp=formatNumberString(cnvrtNumToStrng(qVal[i][j],BUFFER_SIZE)); Output+=tmp+repeatString(" ",lnQ-tmp.length())+"|";			// Charge Valence
			tmp=formatNumberString(cnvrtNumToStrng(semMobVal[i][j],BUFFER_SIZE)); Output+=tmp+repeatString(" ",lnMOB-tmp.length())+"|";		// Electrophoretic Mobility
			tmp=formatNumberString(cnvrtNumToStrng(hMobVal[i][j],BUFFER_SIZE)); Output+=tmp+repeatString(" ",lnMOB-tmp.length())+"|";		// Electrophoretic Mobility
			tmp=formatNumberString(cnvrtNumToStrng(kMobVal[i][j],BUFFER_SIZE)); Output+=tmp+repeatString(" ",lnMOB-tmp.length())+"|";		// Electrophoretic Mobility
			tmp=formatNumberString(cnvrtNumToStrng(difVal[i][j],BUFFER_SIZE)); Output+=tmp+repeatString(" ",lnDIF-tmp.length())+"|\n";		// Diffusivity
			}
		Output+=repeatString("-",theLineWidth)+"\n";
		// Compute Averages and Standard Deviations
		avgPR=calc_average(prVal[i],Z.numFiles); stdPR=calc_std_dev(prVal[i],Z.numFiles,avgPR);										// Anhydrous Protein Radius
		avgSR=calc_average(srVal[i],Z.numFiles); stdSR=calc_std_dev(srVal[i],Z.numFiles,avgSR);										// Solvated Protein Radius
		avgZP=calc_average(zpVal[i],Z.numFiles); stdZP=calc_std_dev(zpVal[i],Z.numFiles,avgZP);										// Zeta Potential
		avgQ=calc_average(qVal[i],Z.numFiles); stdQ=calc_std_dev(qVal[i],Z.numFiles,avgQ);											// Net Charge Valence
		avgkMOB=calc_average(kMobVal[i],Z.numFiles); stdkMOB=calc_std_dev(kMobVal[i],Z.numFiles,avgkMOB);								// Kuwabara Electrophoretic Mobility
		avghMOB=calc_average(hMobVal[i],Z.numFiles); stdhMOB=calc_std_dev(hMobVal[i],Z.numFiles,avghMOB);								// Henry Electrophoretic Mobility
		avgsemMOB=calc_average(semMobVal[i],Z.numFiles); stdsemMOB=calc_std_dev(semMobVal[i],Z.numFiles,avgsemMOB);						// Standard Electrokinetic Model Mobility
		avgDIF=calc_average(difVal[i],Z.numFiles); stdDIF=calc_std_dev(difVal[i],Z.numFiles,avgDIF);									// Single Protein Diffusivity
		// Write Output
		tmp=formatNumberString(cnvrtNumToStrng(avgPR,BUFFER_SIZE)); Output+="Average  |"+tmp+repeatString(" ",lnPR-tmp.length())+"|";		// Protein Radius
		tmp=formatNumberString(cnvrtNumToStrng(avgSR,BUFFER_SIZE)); Output+=tmp+repeatString(" ",lnSR-tmp.length())+"|";					// Solvated Radius
		tmp=formatNumberString(cnvrtNumToStrng(avgZP,BUFFER_SIZE)); Output+=tmp+repeatString(" ",lnZP-tmp.length())+"|";					// Zeta Potential
		tmp=formatNumberString(cnvrtNumToStrng(avgQ,BUFFER_SIZE)); Output+=tmp+repeatString(" ",lnQ-tmp.length())+"|";					// Charge Valence
		tmp=formatNumberString(cnvrtNumToStrng(avgsemMOB,BUFFER_SIZE)); Output+=tmp+repeatString(" ",lnMOB-tmp.length())+"|";				// Electrophoretic Mobility
		tmp=formatNumberString(cnvrtNumToStrng(avghMOB,BUFFER_SIZE)); Output+=tmp+repeatString(" ",lnMOB-tmp.length())+"|";				// Electrophoretic Mobility
		tmp=formatNumberString(cnvrtNumToStrng(avgkMOB,BUFFER_SIZE)); Output+=tmp+repeatString(" ",lnMOB-tmp.length())+"|";				// Electrophoretic Mobility
		tmp=formatNumberString(cnvrtNumToStrng(avgDIF,BUFFER_SIZE)); Output+=tmp+repeatString(" ",lnDIF-tmp.length())+"|\n";				// Diffusivity
		//
		tmp=formatNumberString(cnvrtNumToStrng(stdPR,BUFFER_SIZE)); Output+="Std.Dev. |"+tmp+repeatString(" ",lnPR-tmp.length())+"|";
		tmp=formatNumberString(cnvrtNumToStrng(stdSR,BUFFER_SIZE)); Output+=tmp+repeatString(" ",lnSR-tmp.length())+"|";
		tmp=formatNumberString(cnvrtNumToStrng(stdZP,BUFFER_SIZE)); Output+=tmp+repeatString(" ",lnZP-tmp.length())+"|";
		tmp=formatNumberString(cnvrtNumToStrng(stdQ,BUFFER_SIZE)); Output+=tmp+repeatString(" ",lnQ-tmp.length())+"|";
		tmp=formatNumberString(cnvrtNumToStrng(stdsemMOB,BUFFER_SIZE)); Output+=tmp+repeatString(" ",lnMOB-tmp.length())+"|";		
		tmp=formatNumberString(cnvrtNumToStrng(stdhMOB,BUFFER_SIZE)); Output+=tmp+repeatString(" ",lnMOB-tmp.length())+"|";
		tmp=formatNumberString(cnvrtNumToStrng(stdkMOB,BUFFER_SIZE)); Output+=tmp+repeatString(" ",lnMOB-tmp.length())+"|";
		tmp=formatNumberString(cnvrtNumToStrng(stdDIF,BUFFER_SIZE)); Output+=tmp+repeatString(" ",lnDIF-tmp.length())+"|\n";
		// Compute Collision Efficiency Based on Average Zeta Potential of Structures in a Single Solution Condition
		kValue=strtod(debyeLength[i].c_str(),NULL);		// Debye Length (A)
		erValue=strtod(sER[i].c_str(),NULL);			// Relative Dielectric Constant
		vValue=strtod(sVisc[i].c_str(),NULL);			// Viscosity
		collisionEfficiency=calc_collision_efficiency(avgPR,avgSR,avgZP,erValue,kValue,Z.temperature[i]);
		Output+="Collision Efficiency: "+formatNumberString(cnvrtScientificNumToStrng(collisionEfficiency))+"\n";
		// Compute Coagulation Time Depicting the Time Interval Between Collisions
		coagulationTime=calc_coagulation_time(Z.proteinConc[i],avgMW,vValue,Z.temperature[i]);
		Output+="Coagulation Time: "+formatNumberString(cnvrtScientificNumToStrng(coagulationTime))+" s\n";
		Output+=repeatString("-",theLineWidth)+"|\n\n";
		// Electric Potential Profile (if generated)
		if(numProfilePnts[i][0]>0)
			{Sum=new double[numProfilePnts[i][0]];
			for(int j=0;j<numProfilePnts[i][0];j++){Sum[j]=0;}
			Output+="Electric Potential Profile\n";
			Output+="Distance_[A] Electric_Potential_[V] Standard_Deviation\n";
			for(int j=0;j<Z.numFiles;j++)
				{if(j==0)
					{for(int k=0;k<numProfilePnts[i][j];k++){Output+=position[i][j][k]+"|";}
					Output+="\n";}
				tmp=fNm[j];
				Output+=tmp+repeatString(" ",9-tmp.length())+"|";
				for(int k=0;k<numProfilePnts[i][j];k++)
					{tmp=formatNumberString(cnvrtNumToStrng(electricPotential[i][j][k],BUFFER_SIZE));
					Output+=tmp+"|";
					if(!isnan(electricPotential[i][j][k])){Sum[k]+=electricPotential[i][j][k];}
					}
				Output+="\n";}
			//
			for(int j=0;j<numProfilePnts[i][0];j++)
				{avgPot=Sum[j]/Z.numFiles;
				dArr=new double[Z.numFiles];
				for(int k=0;k<Z.numFiles;k++){dArr[k]=electricPotential[i][k][j];}
				stdPot=calc_std_dev(dArr,Z.numFiles,avgPot);
				delete [] dArr;
				Output+=position[i][0][j]+"|"+formatNumberString(cnvrtNumToStrng(avgPot,BUFFER_SIZE))+"|"+formatNumberString(cnvrtNumToStrng(stdPot,BUFFER_SIZE))+"\n";
				}
			delete [] Sum;
			Output+="\n";}
		}
	// Write Output to File
	fOut.open(zpredOutputFile.c_str(),ofstream::out|ofstream::trunc);
	if(fOut.fail()){cerr<<"ERROR in ZPRED.exe\nZPRED Output file could not be opened.\n"<<zpredOutputFile; exit(EXIT_FAILURE);}
	fOut<<Output; fOut.close();
	cerr<<"Output stored in:\n"<<zpredOutputFile<<"\n";
	// Delete Folders (Exclude Z and LOG folders)
	string msmsFldr=fFldr+"msms/";
	if(!Z.saveTemps)
		{boost::filesystem::remove_all(dxFldr);
		boost::filesystem::remove_all(hydroFldr);
		boost::filesystem::remove_all(inFldr);
		boost::filesystem::remove_all(propkaFldr);
		boost::filesystem::remove_all(pqrFldr);
		boost::filesystem::remove_all(surfFldr);
		boost::filesystem::remove_all(msmsFldr);
		}
}

