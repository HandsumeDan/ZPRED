# include "zpredUtils.h"


using namespace std;
//using namespace boost::numeric::odeint;


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

	// Get Current Working Directory
	string Fldr=get_current_dir_name();Fldr+="/";

	zInput Z=read_zpred_input_file(inFile);
	// Check should be performed by GUI
	//cout<<display_zpred_input(Z)<<endl;
	// Main (ZPRED) Output Folder
	string outputFldr=Fldr+Z.outputTitle+"/"; make_folder(outputFldr);
	Fldr=outputFldr;
	// Formatted Data Folders Folder (Also Holds ZPRED Configuration File)	
	string fFldr=Fldr+"Files/"; make_folder(fFldr);
	//string fFldr=outputFldr+"Files/";make_folder(fFldr);
	// ZPRED Output Folder
	string zpredFldr=fFldr+"Z/"; make_folder(zpredFldr);
	// Program Stdout/Stderr Log Folder
	string logFldr=fFldr+"LOG/"; make_folder(logFldr);
	// APBS Output Folder (in DX Matrix Format)
	string dxFldr=fFldr+"DX/"; make_folder(dxFldr);
	// HYDROPRO Output Folder (HP stands for hydrodynamic properties)
	string hydroFldr=fFldr+"HP/"; make_folder(hydroFldr);
	// APBS Input Command File Folder
	string inFldr=fFldr+"IN/"; make_folder(inFldr);
	// PROPKA Output Folder
	string propkaFldr=fFldr+"PROPKA/"; make_folder(propkaFldr);
	// PQR Formatted PDB Files Folder
	string pqrFldr=fFldr+"PQR/"; make_folder(pqrFldr);
	// MSMS Surface Output Folder
	string surfFldr=fFldr+"S/"; make_folder(surfFldr);
	string surfFldrContents[6]={"CGO/","colorCGO/","CSV/","fromMSMS/","SAS/","fromMV/"};

	// PyMol PQR2PDB Converter
	string pqr2pdbFile=Fldr+"saveAsPdb.py"; write_saveAsPdb(pqr2pdbFile);
	// APBS Inputgen Python Script
	string apbsInputgen=fFldr+"apbsInputgen.py"; writeAPBSInputGen(apbsInputgen);
	// APBS psize.py supports apbsInputgen.py
	string psize=fFldr+"psize.py"; writePsize(psize);
	// HullRad (Grisham Edit)
	string theHullRad=fFldr+"editedHullRad.py"; writeHullRad(theHullRad);

	//cout<<"|"<<get_containing_folder(theHullRad)<<"|"<<get_file_name(theHullRad)<<"|"<<endl;

	//string proFile="/home/jenbunnyaqua/ZPRED_II/Test/Files/PROPKA/6LYZ_1/pH7.00_1.propka";
	//cout<<readPROPKAFile(proFile)<<endl;

	//string aPqrFldr="/home/jenbunnyaqua/ZPRED_II/fromKen/PQR/";
	//string proFile="/home/jenbunnyaqua/ZPRED_II/fromKen/mp90_frame.pqr";
	//extract_models_from_pqr_file(proFile,aPqrFldr);

	//string pqFile="/home/jenbunnyaqua/ZPRED_II/fromKen/PQR/3.pqr";
	//string pdFile="/home/jenbunnyaqua/ZPRED_II/fromKen/PQR/3.pdb";
	//convert_pqr_to_pdb(pqFile,pdFile);
	//exit(1);

	// String Array Delimiter
	string delimiter=";";
	// Determine Total Number of Zeta Potential Calculations
	int totalCalc=Z.numFiles*Z.numSolCond;
	string display,*cmdStrng,*cmd;
	vector<thread> t;
	int Counter=0,start=0,end=0,remainingCalc=totalCalc,numCalc=0;
	bool RUNNING_MULTI_DISPLAY;
	clock_t Tme;
	// Run ZPRED
	cerr<<"Starting ZPRED timer (";
	time_t timeResult=time(nullptr);
	char *startTimeValue=asctime(localtime(&timeResult));
	cerr<<timeResult<<")...\n";
	//tmp=startTimeValue;
	int startTime=timeResult;//atoi(startTimeValue);
	//Tme=clock();
	// Initialize GUI Display
	string guiHdr="<b>Executing "+cnvrtNumToStrng(totalCalc,0)+" ZPRED calculations ("+cnvrtNumToStrng(Z.numFiles,0)+" structures in "+cnvrtNumToStrng(totalCalc/Z.numFiles,0)+" solution conditions)</b><br>";
	guiHdr+="<b>&nbsp;---------------&nbsp;-----------------------&nbsp;&nbsp;--------------&nbsp;&nbsp;-----------------&nbsp;&nbsp;-------------&nbsp;&nbsp;----------- </b><br>";
	guiHdr+="<b>| A = APBS | H = HYDROPRO | M = MSMS | MV = MULTIVALUE | P = PDB2PQR | Z = ZPRED |</b><br>";
	guiHdr+="<b>&nbsp;---------------&nbsp;-----------------------&nbsp;&nbsp;--------------&nbsp;&nbsp;-----------------&nbsp;&nbsp;-------------&nbsp;&nbsp;----------- </b><br>";
	cerr<<"Executing "<<totalCalc<<" ZPRED calculations ("<<Z.numFiles<<" structures in "<<totalCalc/Z.numFiles<<" solution conditions)\n";
	cerr<<" ---------- -------------- ---------- ----------------- ------------- ----------- \n";
	if(Z.RUN_HULLRAD)
		{cerr<<"| A = APBS | H = HULLRAD  | M = MSMS | MV = MULTIVALUE | P = PDB2PQR | Z = ZPRED |\n";}
	else
		{cerr<<"| A = APBS | H = HYDROPRO | M = MSMS | MV = MULTIVALUE | P = PDB2PQR | Z = ZPRED |\n";}
	cerr<<" ---------- -------------- ---------- ----------------- ------------- ----------- \n";
	display=initialize_ZPRED_display(totalCalc,MULTI_DISPLAY); cerr<<display; rewind_display(display);
	string guiDisplayFile=Fldr+".gui_display.txt";
	ofstream fOut;
	fOut.open(guiDisplayFile.c_str(),ofstream::out|ofstream::trunc);
	if(fOut.fail()){cerr<<"ERROR in ZPRED\nZPRED Output file could not be opened.\n"<<guiDisplayFile;exit(EXIT_FAILURE);}
	fOut<<guiHdr+initialize_ZPRED_gui_display(totalCalc,MULTI_DISPLAY); fOut.close();
	// Convert ZPRED Input Parameters into string
	cmdStrng=cnvrtzInput2Str(Z,delimiter,Fldr);

	if(MULTI_DISPLAY)
		{// Only Perform 150 threads MAX at a time
		RUNNING_MULTI_DISPLAY=true;
		start=0;
		end=MAX_DISPLAY_ROWS*MAX_DISPLAY_COLS;
		while(RUNNING_MULTI_DISPLAY)
			{Counter=0;
			for(int i=start;i<end;i++)
				{// Execute ZPRED on Separate Thread for Each Structure and Solution conditions
				cmd=new string(cmdStrng[i]);
				// Append Calculation Number and initialized display string
				*cmd+=cnvrtNumToStrng(i+1,0)+";"+display+";";
				if(Counter<Z.maxThreads)
					{t.push_back(thread(runZPRED,(void*)cmd));
					Counter++;}
				else{// Counter Reached Maximum Threads, wait for their completion
					for(int j=i-Z.maxThreads;j<i;j++){t[j].join();}
					t.push_back(thread(runZPRED,(void*)cmd));
					Counter=1;}}
			for(int i=0;i<Counter;i++){t[end+i-Counter].join();}
			// Once Finished with 150 threads start again with remainder
			start+=MAX_DISPLAY_ROWS*MAX_DISPLAY_COLS;
			remainingCalc-=MAX_DISPLAY_ROWS*MAX_DISPLAY_COLS;
			if(remainingCalc>=MAX_DISPLAY_ROWS*MAX_DISPLAY_COLS){end+=MAX_DISPLAY_ROWS*MAX_DISPLAY_COLS;}
			else if(remainingCalc<=0){RUNNING_MULTI_DISPLAY=false; break;}
			else{end+=remainingCalc;}
			// Update Display
			numCalc=end-start;
			display=update_display(totalCalc,numCalc,start);
			cerr<<display;
			rewind_display(display);
			}
		}
	else{
		for(int i=0;i<totalCalc;i++)
			{cmd=new string(cmdStrng[i]);
			*cmd+=cnvrtNumToStrng(i+1,0)+";"+display+";";
			if(Counter<Z.maxThreads){t.push_back(thread(runZPRED,(void*)cmd));Counter++;}
			else{for(int j=i-Z.maxThreads;j<i;j++){t[j].join();}
				t.push_back(thread(runZPRED,(void*)cmd));
				Counter=1;}}
		for(int i=0;i<Counter;i++){t[totalCalc+i-Counter].join();}}
	erase_display();
	//cerr<<totalCalc<<" ZPRED calculations complete.\nOutput stored in:\n"<<outputFldr+Z.outputTitle.substr(0,Z.outputTitle.length()-1)+".txt"<<"\n";
	cerr<<totalCalc<<" ZPRED calculations complete.\n";
	//Tme=clock()-Tme;
	timeResult=time(nullptr);
	char *endTimeValue=asctime(localtime(&timeResult));
	//cerr<<endTimeValue<<endl;
	//tmp=endTimeValue;
	int endTime=timeResult;	//atoi(endTimeValue);
	//int endTime=(int)asctime(localtime(&timeResult));
	//cerr<<"ZPRED Timer = "<<((float)Tme)/CLOCKS_PER_SEC<<" s\n";
	cerr<<"ZPRED Timer = "<<endTime-startTime<<" s\n";
	// After all ZPRED computations complete, store in user-specified folder (outputTitle parameter in input File)
	string zpredOutputFile=outputFldr+Z.outputTitle+".txt";
	// Read List of ZPRED Output Files
	string zOutListFile=fFldr+"zpredOutFilesList.txt",zFile="";	
	ifstream fIn;
	fIn.open(zOutListFile.c_str());
	if(fIn.fail()){cerr<<"ERROR in theZPRED!\nZPRED Output File List file could not be opened.\n"<<zOutListFile<<endl;exit(EXIT_FAILURE);}
	int Sz=15000;
	char Val[Sz];
	zOutput zRes;
	// Collect and Average Zeta Potentials of the ensemble of structures for each solution condition
	int* fCounter=new int[Z.numSolCond];
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
	double* dArr;
	// 9|9|8|9|6|8|8|
	// Determine header Width
	int lnPR=0,lnSR=0,lnZP=0,lnQ=0,lnMOB=0,lnDIF=0;
	string tmp2,tmp3,tmp4,tmp5,tmp6;
	for(int j=0;j<Z.numFiles;j++)
		{//fNm=cnvrtNumToStrng(j,0);
		//Output+="file #"+fNm+repeatString(" ",3-fNm.length())+"|";
		tmp=formatNumberString(cnvrtNumToStrng(prVal[0][j],BUFFER_SIZE));
		if(tmp.length()>lnPR){lnPR=tmp.length();}
		tmp2=formatNumberString(cnvrtNumToStrng(srVal[0][j],BUFFER_SIZE));
		if(tmp2.length()>lnSR){lnSR=tmp2.length();}
		tmp3=formatNumberString(cnvrtNumToStrng(zpVal[0][j],BUFFER_SIZE));
		if(tmp3.length()>lnZP){lnZP=tmp3.length();}
		tmp4=formatNumberString(cnvrtNumToStrng(qVal[0][j],BUFFER_SIZE));
		if(tmp4.length()>lnQ)
			{//lnQ=tmp4.length();
			lnQ=BUFFER_SIZE;
			lnQ+=2;
			}
		tmp5=formatNumberString(cnvrtNumToStrng(semMobVal[0][j],BUFFER_SIZE));
		if(tmp5.length()>lnMOB){lnMOB=tmp5.length();}
		tmp6=formatNumberString(cnvrtNumToStrng(difVal[0][j],BUFFER_SIZE));
		if(tmp6.length()>lnDIF){lnDIF=tmp6.length();}
		}

	int theLineWidth=131;
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
		Output+="Protein Concentration: "+formatNumberString(cnvrtNumToStrng(Z.proteinConc[i],BUFFER_SIZE))+" g/L\n";
		Output+=repeatString("-",theLineWidth)+"|\n";
		Output+=repeatString("|",10)+repeatString("-",lnPR)+"|"+repeatString("-",lnSR)+"|"+repeatString("-",lnZP)+"|"+repeatString("-",lnQ)+"|SEM"+repeatString("-",lnMOB-3)+"|Henry"+repeatString("-",lnMOB-5)+"|Kuwabara"+repeatString("-",lnMOB-8)+"|"+repeatString("-",lnDIF)+"|\n";

		Output+=repeatString("|",10)+"Anhydrous"+repeatString("-",lnPR-9)+"|Solvated"+repeatString("-",lnSR-8)+"|Zeta"+repeatString("-",lnZP-4)+"|"+repeatString("-",lnQ)+"|Electro."+repeatString("-",lnMOB-8)+"|Electro."+repeatString("-",lnMOB-8)+"|Electro."+repeatString("-",lnMOB-8)+"|"+repeatString("-",lnDIF)+"|\n";
		Output+=repeatString("|",10)+"Radius"+repeatString("-",lnPR-6)+"|Radius"+repeatString("-",lnSR-6)+"|Potential"+repeatString("-",lnZP-9)+"|Q"+repeatString("-",lnQ-1)+"|Mobility"+repeatString("-",lnMOB-8)+"|Mobility"+repeatString("-",lnMOB-8)+"|Mobility"+repeatString("-",lnMOB-8)+"|Diffusivity"+repeatString("-",lnDIF-11)+"|\n";
		Output+=repeatString("|",10)+"[A]"+repeatString("-",lnPR-3)+"|[A]"+repeatString("-",lnSR-3)+"|[Volts]"+repeatString("-",lnZP-7)+"|"+repeatString("-",lnQ)+"|[umcm/Vs]"+repeatString("-",lnMOB-9)+"|[umcm/Vs]"+repeatString("-",lnMOB-9)+"|[umcm/Vs]"+repeatString("-",lnMOB-9)+"|[m^2/s]"+repeatString("-",lnDIF-7)+"|\n";
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
		// Anhydrous Protein Radius
		avgPR=calc_average(prVal[i],Z.numFiles); stdPR=calc_std_dev(prVal[i],Z.numFiles,avgPR);
		// Solvated Protein Radius
		avgSR=calc_average(srVal[i],Z.numFiles); stdSR=calc_std_dev(srVal[i],Z.numFiles,avgSR);
		// Zeta Potential
		avgZP=calc_average(zpVal[i],Z.numFiles); stdZP=calc_std_dev(zpVal[i],Z.numFiles,avgZP);
		// Net Valence
		avgQ=calc_average(qVal[i],Z.numFiles); stdQ=calc_std_dev(qVal[i],Z.numFiles,avgQ);
		// Kuwabara Electrophoretic Mobility
		avgkMOB=calc_average(kMobVal[i],Z.numFiles); stdkMOB=calc_std_dev(kMobVal[i],Z.numFiles,avgkMOB);
		// Henry Electrophoretic Mobility
		avghMOB=calc_average(hMobVal[i],Z.numFiles); stdhMOB=calc_std_dev(hMobVal[i],Z.numFiles,avghMOB);
		// Standard Electrokinetic Model Mobility
		avgsemMOB=calc_average(semMobVal[i],Z.numFiles); stdsemMOB=calc_std_dev(semMobVal[i],Z.numFiles,avgsemMOB);
		// Single Protein Diffusivity
		avgDIF=calc_average(difVal[i],Z.numFiles); stdDIF=calc_std_dev(difVal[i],Z.numFiles,avgDIF);
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
		//
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
	if(fOut.fail()){cerr<<"ERROR in ZPRED\nZPRED Output file could not be opened.\n"<<zpredOutputFile; exit(EXIT_FAILURE);}
	fOut<<Output; fOut.close();
	cerr<<"Output stored in:\n"<<zpredOutputFile<<"\n";
	// Delete Folders (Exclude Z and LOG folders)
	string msmsFldr=fFldr+"msms/";
	if(!Z.saveTemps)
		{
		boost::filesystem::remove_all(dxFldr);
		boost::filesystem::remove_all(hydroFldr);
		boost::filesystem::remove_all(inFldr);
		boost::filesystem::remove_all(propkaFldr);
		boost::filesystem::remove_all(pqrFldr);
		boost::filesystem::remove_all(surfFldr);
		boost::filesystem::remove_all(msmsFldr);
		}
}

