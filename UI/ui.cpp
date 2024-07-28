# include <cmath>
# include <cstdio>
# include <cstring>
# include <fstream>
# include <iostream>
# include <sstream>
# include <string>
# include <fcntl.h>					// Defines O_CREAT flag
# include <csignal>
# include <thread>
# include <sys/reg.h>
# include <sys/syscall.h>			// Defines OS System Calls
# include <sys/ptrace.h>
# include <sys/types.h>
# include <sys/user.h>
# include <sys/wait.h>
# include <sys/stat.h>				// Defines chmod function
# include <unistd.h>
# include <boost/filesystem.hpp>		// create_directory

# include <QtWidgets>
# include <QProcess>

# include "ui.h"

using namespace std;

UI::UI(QWidget *parent) : QWidget(parent)
	{string inputFileMsg="ZPRED uses the Protein Data Bank (PDB) format for its input molecular structures.  Ideally, they should be generated from a molecular dynamics simulation accounting for the structural motion of the folded structure.";
    // Initialize Values
	numSolvent=new int[1];numSolvent[0]=1;
	numSolute=new int[1];numSolute[0]=1;
	
	QString Fldr=QDir::currentPath();
	string sFldr=Fldr.toStdString();
	sFldr=checkFinalBackSlash(sFldr);
	// File Reference for When Downloads already Performed (makes .history.txt)
	refFile=getBaseFolder(sFldr)+".history.txt";
	// File Reference for When Optional Parameters are set (makes .optionalParameters.txt)
	optRefFile=getBaseFolder(sFldr)+".optionalParameters.txt";
	
	QFont font;
	font.setBold(true);

	ifstream fIn;
	ofstream fOut;
	fIn.open(refFile.c_str());
	string theTxt=".apbsFldr:\n.filesFldr:\n.hydroproFldr:\n.msmsFldr:\n.pdb2pqrFldr:\n.zpredFldr:\n";
	QPushButton* autoFillButton=new QPushButton(tr("Auto-Fill"));
	autoFillButton->setFixedWidth(100);
	connect(autoFillButton,SIGNAL(clicked()),this,SLOT(autoFillFilePaths()));
	if(!fIn.fail())
		{// Define Previous Download Folders
		fIn.close();defineDownloadFolders(refFile);
		}
	else
		{apbsFldr=sFldr;
		filesFldr=sFldr;
		hydroproFldr=sFldr;
		msmsFldr=sFldr;
		pdb2pqrFldr=sFldr;
		zpredFldr=sFldr;
		// Write New Ref
		fOut.open(refFile.c_str(),ofstream::out|ofstream::trunc);fOut<<theTxt;fOut.close();
		autoFillButton->hide();}

	fIn.open(optRefFile.c_str());
	string optTxt="// APBS Parameters\n.pdie:4.0\n.dx:161\n.dy:161\n.dz:161\n.atompot:false\n.charge:false\n.dielx:false\n.diely:false\n.dielz:false\n.edens:false\n.ivdw:false\n.kappa:false\n.lap:false\n.ndens:false\n.pot:true\n.qdens:false\n.smol:false\n.sspl:false\n.vdw:false\n";
	optTxt+="// HYDROPRO Parameters\n.calcType:1\n";
	optTxt+="// MSMS Parameters\n.lowPntDensity:1.00\n.highPntDensity:3.00\n";
	optTxt+="// PDB2PQR Parameters\n.forceField:PARSE\n";
	if(!fIn.fail()){fIn.close();}
	else
		{// Write New Ooptional Parameters Ref File
		fOut.open(optRefFile.c_str(),ofstream::out|ofstream::trunc);fOut<<optTxt;fOut.close();
		}

	// Solute(s) Label Array
	for(int i=0;i<MAX_SOLUTION_CONDITIONS;i++)
		{for(int j=0;j<MAX_SOLUTES;j++)
			{soluteLabel[i][j]=new QLabel(this);
			fSoluteLabel[i][j]=new QLabel(this);
			soluteConcLine[i][j]=new QLineEdit(this);
			soluteConcLabel[i][j]=new QLabel(this);
			//soluteConcTypeLabel[i][j]=new QLabel(this);
			}
		soluteConcTypeLabel[i]=new QLabel(this);}
	for(int i=0;i<MAX_SOLUTION_CONDITIONS;i++)
		{for(int j=0;j<MAX_SOLUTES;j++)
			{soluteLabel[i][j]->setText(tr(NO_VALUE));soluteLabel[i][j]->setAlignment(Qt::AlignCenter);
			soluteConcLine[i][j]->setText(tr(NO_VALUE));soluteConcLine[i][j]->setAlignment(Qt::AlignCenter);
			soluteConcLabel[i][j]->setText(tr(NO_VALUE));soluteConcLabel[i][j]->setAlignment(Qt::AlignCenter);
			//soluteConcTypeLabel[i][j]->setText(tr(NO_VALUE));soluteConcTypeLabel[i][j]->setAlignment(Qt::AlignCenter);
			}
		soluteConcTypeLabel[i]->setText(tr(NO_VALUE));soluteConcTypeLabel[i]->setAlignment(Qt::AlignCenter);}
	// Solvent(s) Label Array
	for(int i=0;i<MAX_SOLUTION_CONDITIONS;i++)
		{for(int j=0;j<MAX_SOLVENTS;j++)
			{solventLabel[i][j]=new QLabel(this);
			solventConcLine[i][j]=new QLineEdit(this);
			solventConcLabel[i][j]=new QLabel(this);
			}
		solventConcTypeLabel[i]=new QLabel(this);}
	for(int i=0;i<MAX_SOLUTION_CONDITIONS;i++)
		{for(int j=0;j<MAX_SOLVENTS;j++)
			{solventLabel[i][j]->setText(tr(NO_VALUE));solventLabel[i][j]->setAlignment(Qt::AlignCenter);
			solventConcLine[i][j]->setText(tr(NO_VALUE));solventConcLine[i][j]->setAlignment(Qt::AlignCenter);
			solventConcLabel[i][j]->setText(tr(NO_VALUE));solventConcLabel[i][j]->setAlignment(Qt::AlignCenter);
			}
		solventConcTypeLabel[i]->setText(tr(NO_VALUE));solventConcTypeLabel[i]->setAlignment(Qt::AlignCenter);}	
	// GUI Features
	for(int i=0;i<MAX_SOLUTION_CONDITIONS;i++)
		{// pH
		pH[i]=new QLineEdit;pH[i]->setAlignment(Qt::AlignCenter);pH[i]->setText(tr(NO_VALUE));
		connect(pH[i],SIGNAL(returnPressed()),this,SLOT(pH_defined()));
		pHLabel[i]=new QLabel(tr("<b><span style=\"color:red;\">pH</span></b>:"));pHLabel[i]->setBuddy(pH[i]);pHLabel[i]->setAlignment(Qt::AlignCenter|Qt::AlignRight);
		// Temperature
		T[i]=new QLineEdit;T[i]->setAlignment(Qt::AlignCenter);T[i]->setText(tr(NO_VALUE));
		connect(T[i],SIGNAL(returnPressed()),this,SLOT(T_defined()));
		TLabel[i]=new QLabel(tr("<b><span style=\"color:red;\">Temperature</span></b>:"));TLabel[i]->setBuddy(T[i]);TLabel[i]->setAlignment(Qt::AlignCenter|Qt::AlignRight);
		// Protein Concentration
		pC[i]=new QLineEdit;pC[i]->setText(tr(NO_VALUE));pC[i]->setAlignment(Qt::AlignCenter);
		connect(pC[i],SIGNAL(returnPressed()),this,SLOT(pC_defined()));
		pCLabel[i]=new QLabel(tr("<b><span style=\"color:red;\">Protein Concentration</span></b>:"));pCLabel[i]->setAlignment(Qt::AlignRight);
		// 'Select a Solvent' Labels
		seleSolvLabel[i]=new QLabel(tr("<b><span style=\"color:red;\">Select a Solvent</span></b>"));seleSolvLabel[i]->setAlignment(Qt::AlignCenter);
		seleSolvConcLabel[i]=new QLabel(tr("<b><span style=\"color:red;\">Fraction</span></b>"));seleSolvConcLabel[i]->setAlignment(Qt::AlignCenter);
		seleSolvConcTypeLabel[i]=new QLabel(tr("<b><span style=\"color:red;\">Units</span></b>"));seleSolvConcTypeLabel[i]->setAlignment(Qt::AlignCenter);
		// 'Select a Solute' Labels
		seleSoluLabel[i]=new QLabel(tr("<b><span style=\"color:red;\">Select a Solute</span></b>"));seleSoluLabel[i]->setAlignment(Qt::AlignCenter);
		seleSoluConcLabel[i]=new QLabel(tr("<b><span style=\"color:red;\">Concentration</span></b>"));seleSoluConcLabel[i]->setAlignment(Qt::AlignCenter);
		seleSoluConcTypeLabel[i]=new QLabel(tr("<b><span style=\"color:red;\">Units</span></b>"));seleSoluConcTypeLabel[i]->setAlignment(Qt::AlignCenter);
		// User Specified Density
		fDens[i]=new QLineEdit;fDens[i]->setAlignment(Qt::AlignCenter);fDens[i]->setText(tr(NO_VALUE));
		// User Specified Dielectric
		fDiel[i]=new QLineEdit;fDiel[i]->setAlignment(Qt::AlignCenter);fDiel[i]->setText(tr(NO_VALUE));
		// User Specified Viscosity
		fVisc[i]=new QLineEdit;fVisc[i]->setAlignment(Qt::AlignCenter);fVisc[i]->setText(tr(NO_VALUE));
		// User Specified Hydration Thickness
		fXsp[i]=new QLineEdit;fXsp[i]->setAlignment(Qt::AlignCenter);fXsp[i]->setText(tr(NO_VALUE));
		// User Specified Parameters Check Box
		forceParmBox[i]=new QCheckBox(tr("Specify Solution Parameters Instead of Computation?"));forceParmBox[i]->setChecked(false);forceParmBox[i]->setFont(font);
		connect(forceParmBox[i],SIGNAL(stateChanged(int)),this,SLOT(show_fBox()));
		// Generate Electric Potential Profile CHeckBox
		potProfileBox[i]=new QCheckBox(tr("Generate Electric Potential Profile?"));potProfileBox[i]->setChecked(false);potProfileBox[i]->setFont(font);
		}
	jobNameLabel=new QLabel(tr("<b><span style=\"color:red;\">Job Name</span></b>:"));jobNameLabel->setAlignment(Qt::AlignCenter|Qt::AlignRight);
	jobName=new QLineEdit;jobName->setAlignment(Qt::AlignCenter);jobName->setText(tr(NO_VALUE));jobName->setMaxLength(16);
	
	//connect(jobName,SIGNAL(editingFinished()),this,SLOT(jobName_defined()));
	connect(jobName,SIGNAL(returnPressed()),this,SLOT(jobName_defined()));

	QGroupBox* gBox=new QGroupBox(tr("|-INPUT PARAMETERS-|"));
	//gBox->setStyleSheet("QGroupBox { font-size:12pt; font-weight:bold; color:blue; } ");//gBox->setStyleSheet("QGroupBox { font-weight:bold; } ");	
	gBox->setFont(font);gBox->setAlignment(Qt::AlignRight);

	QPushButton* openButton=new QPushButton(tr("Select Structure Files"));openButton->setToolTip(tr(inputFileMsg.c_str()));
	openButton->setFixedWidth(200);openButton->setFont(font);
	connect(openButton,SIGNAL(clicked()),this,SLOT(openPdbFiles()));

	fileLabel=new QLabel(tr("<span style=\"color:red;\">No files selected</span>"));fileLabel->setAlignment(Qt::AlignRight);
	
	QPushButton* addSolCondButton=new QPushButton(tr("Add Solution Condition"));
	addSolCondButton->setFixedWidth(200);
	connect(addSolCondButton,SIGNAL(clicked()),this,SLOT(addSolCondTab()));
	
	remSolCondButton=new QPushButton(tr("Remove Solution Condition"));remSolCondButton->hide();
	remSolCondButton->setFixedWidth(200);
	connect(remSolCondButton,SIGNAL(clicked()),this,SLOT(remSolCondTab()));	
	// Solution Condition Tab
	scTab = new QTabWidget;
	
	QLabel* TUnitsLabel=new QLabel(tr("[Kelvin]"));TUnitsLabel->setAlignment(Qt::AlignLeft);
	QLabel* pCUnitsLabel=new QLabel(tr("[g/L]"));pCUnitsLabel->setAlignment(Qt::AlignLeft);
	QLabel* hiddenLabel1=new QLabel;
	
	// Define Solution Conditions Tab
	QWidget *tab1 = new QWidget;
	//
	TLabel[0]->setFixedWidth(150);
	T[0]->setFixedWidth(LINE_INPUT_SIZE);
	pH[0]->setFixedWidth(LINE_INPUT_SIZE);
	pC[0]->setFixedWidth(LINE_INPUT_SIZE);
	//
	QGridLayout *tab1hbox = new QGridLayout;
	tab1hbox->setSizeConstraint(QLayout::SetFixedSize);
	tab1hbox->addWidget(TLabel[0],0,0,1,1);tab1hbox->addWidget(T[0],0,1,1,1);tab1hbox->addWidget(TUnitsLabel,0,2,1,1);tab1hbox->addWidget(pHLabel[0],0,4,1,1);tab1hbox->addWidget(pH[0],0,5,1,1);
	tab1hbox->addWidget(create_oneSolventBox(0),1,0,1,3);tab1hbox->addWidget(create_oneSoluteBox(0),1,3,1,3);
	tab1hbox->addWidget(forceParmBox[0],2,0,1,3);tab1hbox->addWidget(pCLabel[0],2,3,1,1);tab1hbox->addWidget(pC[0],2,4,1,1);tab1hbox->addWidget(pCUnitsLabel,2,5,1,1);
	tab1hbox->addWidget(potProfileBox[0],3,0,1,3);
	tab1->setLayout(tab1hbox);

	scTab->addTab(tab1, tr("&0"));

	QGridLayout *ilayout = new QGridLayout;
	ilayout->setSizeConstraint(QLayout::SetFixedSize);
	ilayout->addWidget(jobNameLabel,0,0,1,1);ilayout->addWidget(jobName,0,1,1,1);ilayout->addWidget(hiddenLabel1,0,2,1,1);ilayout->addWidget(openButton,0,3,1,1);
	ilayout->addWidget(fileLabel,1,0,1,4);
	ilayout->addWidget(addSolCondButton,2,0,1,1);ilayout->addWidget(remSolCondButton,2,1,1,1);
	ilayout->addWidget(scTab,3,0,1,4);
	gBox->setLayout(ilayout);

	// Specify External Executables
	QGroupBox* eBox=new QGroupBox;
	eBox->setFont(font);
	// APBS
	QGroupBox* apbsBox=new QGroupBox(tr("APBS"));
	//apbsBox->setFont(font);
	apbsFilePath=new QLineEdit;apbsFilePath->setFixedWidth(100);
	QPushButton* aInfoButton=new QPushButton(tr("Installation"));
	aInfoButton->setFixedWidth(100);
	connect(aInfoButton,SIGNAL(clicked()),this,SLOT(launchApbsInstallInfo()));
	QPushButton* apbsButton=new QPushButton(tr("..."));apbsButton->setToolTip(tr("Specify the file path of the Adaptive Poisson-Boltzmann Solver (APBS) executable file."));
	apbsButton->setFixedWidth(40);
	connect(apbsButton,SIGNAL(clicked()),this,SLOT(openAPBS()));
	apbsLabel=new QLabel;
	apbsLabel->setText("<span style=\"color:red;\">apbs</span>:");
	apbsLabel->setAlignment(Qt::AlignCenter|Qt::AlignRight);apbsLabel->setToolTip(tr("APBS Description"));
	// MULTIVALUE
	mvFilePath=new QLineEdit;mvFilePath->setFixedWidth(100);
	QPushButton* mvButton=new QPushButton(tr("..."));mvButton->setToolTip(tr("Specify the file path of the \'multivalue\' executable file.  It is an APBS tool and will most likely be buried somewhere in the APBS files after its installation (e.g. APBS-1.5-linux64/share/apbs/tools/bin/ )."));
	mvButton->setFixedWidth(40);
	connect(mvButton,SIGNAL(clicked()),this,SLOT(openMULTIVALUE()));
	mvLabel=new QLabel(tr("<span style=\"color:red;\">multivalue</span>:"));mvLabel->setAlignment(Qt::AlignCenter|Qt::AlignRight);mvLabel->setToolTip(tr("[Brief Description]"));
	// APBS Box Layout
	QGridLayout *apbslayout=new QGridLayout;
	//apbslayout->setSizeConstraint(QLayout::SetFixedSize);
	apbslayout->addWidget(aInfoButton,0,1,1,2);
	apbslayout->addWidget(apbsLabel,1,0,1,1);apbslayout->addWidget(apbsFilePath,1,1,1,1);apbslayout->addWidget(apbsButton,1,2,1,1);
	apbslayout->addWidget(mvLabel,2,0,1,1);apbslayout->addWidget(mvFilePath,2,1,1,1);apbslayout->addWidget(mvButton,2,2,1,1);
	apbsBox->setLayout(apbslayout);
	//
	// HYDROPRO
	hydroBox=new QGroupBox(tr("HYDROPRO"));
	//
	QPushButton* hInfoButton=new QPushButton(tr("Installation"));
	hInfoButton->setFixedWidth(100);
	connect(hInfoButton,SIGNAL(clicked()),this,SLOT(launchHydroproInstallInfo()));
	hFilePath=new QLineEdit;hFilePath->setFixedWidth(100);
	QPushButton* hButton=new QPushButton(tr("..."));hButton->setToolTip(tr("Specify the file path of the HYDROPRO executable file."));
	hButton->setFixedWidth(40);
	connect(hButton,SIGNAL(clicked()),this,SLOT(openHYDROPRO()));
	hLabel=new QLabel;
	hLabel->setText("<span style=\"color:red;\">hydropro-lnx.exe</span>:");
	hLabel->setAlignment(Qt::AlignCenter|Qt::AlignRight);hLabel->setToolTip(tr("HYDROPRO Description"));
	// HYDROPRO Box Layout
	QGridLayout *hydrolayout=new QGridLayout;
	//hydrolayout->addWidget(dLabel2,0,0,1,3);
	hydrolayout->addWidget(hInfoButton,0,1,1,2);
	hydrolayout->addWidget(hLabel,1,0,1,1);hydrolayout->addWidget(hFilePath,1,1,1,1);hydrolayout->addWidget(hButton,1,2,1,1);
	hydroBox->setLayout(hydrolayout);
	hydroBox->hide();
	//
	// MSMS
	QGroupBox* msmsBox=new QGroupBox(tr("MSMS"));
	QPushButton* mInfoButton=new QPushButton(tr("Installation"));
	mInfoButton->setFixedWidth(100);
	connect(mInfoButton,SIGNAL(clicked()),this,SLOT(launchMsmsInstallInfo()));
	//
	msFilePath=new QLineEdit;msFilePath->setFixedWidth(100);
	QPushButton* msButton=new QPushButton(tr("..."));msButton->setToolTip(tr("Specify the file path of the MSMS executable file."));
	msButton->setFixedWidth(40);
	connect(msButton,SIGNAL(clicked()),this,SLOT(openMSMS()));
	msLabel=new QLabel;
	msLabel->setText("<span style=\"color:red;\">msms</span>:");
	msLabel->setAlignment(Qt::AlignCenter|Qt::AlignRight);msLabel->setToolTip(tr("MSMS Description"));
	// MSMS Atom Type File
	msaFilePath=new QLineEdit;msaFilePath->setFixedWidth(100);
	QPushButton* msaButton=new QPushButton(tr("..."));msaButton->setToolTip(tr("Specify the file path of the MSMS \'atmtypenumbers\' supplementary file."));
	msaButton->setFixedWidth(40);
	connect(msaButton,SIGNAL(clicked()),this,SLOT(openMSMS_A()));
	msaLabel=new QLabel(tr("<span style=\"color:red;\">atmTypeNumbers</span>:"));msaLabel->setBuddy(msaButton);msaLabel->setAlignment(Qt::AlignCenter|Qt::AlignRight);msaLabel->setToolTip(tr("[Brief Description]"));
	// MSMS XYZR File
	msxFilePath=new QLineEdit;msxFilePath->setFixedWidth(100);
	QPushButton* msxButton=new QPushButton(tr("..."));msxButton->setToolTip(tr("Specify the file path of the MSMS \'pdb_to_xyzr\' supplementary file."));
	msxButton->setFixedWidth(40);
	connect(msxButton,SIGNAL(clicked()),this,SLOT(openMSMS_X()));
	msxLabel=new QLabel(tr("<span style=\"color:red;\">pbd_to_xyzr</span>:"));msxLabel->setBuddy(msxButton);msxLabel->setAlignment(Qt::AlignCenter|Qt::AlignRight);msxLabel->setToolTip(tr("[Brief Description]"));
	// MSMS XYZRN File
	msxnFilePath=new QLineEdit;msxnFilePath->setFixedWidth(100);
	QPushButton* msxnButton=new QPushButton(tr("..."));msxnButton->setToolTip(tr("Specify the file path of the MSMS \'pdb_to_xyzrn\' supplementary file."));
	msxnButton->setFixedWidth(40);
	connect(msxnButton,SIGNAL(clicked()),this,SLOT(openMSMS_XN()));
	msxnLabel=new QLabel(tr("<span style=\"color:red;\">pdb_to_xyzrn</span>:"));msxnLabel->setBuddy(msxnButton);msxnLabel->setAlignment(Qt::AlignCenter|Qt::AlignRight);msxnLabel->setToolTip(tr("[Brief Description]"));
	// MSMS Box Layout
	QGridLayout *msmslayout=new QGridLayout;
	msmslayout->addWidget(mInfoButton,0,1,1,2);
	msmslayout->addWidget(msLabel,1,0,1,1);msmslayout->addWidget(msFilePath,1,1,1,1);msmslayout->addWidget(msButton,1,2,1,1);
	msmslayout->addWidget(msaLabel,2,0,1,1);msmslayout->addWidget(msaFilePath,2,1,1,1);msmslayout->addWidget(msaButton,2,2,1,1);
	msmslayout->addWidget(msxLabel,3,0,1,1);msmslayout->addWidget(msxFilePath,3,1,1,1);msmslayout->addWidget(msxButton,3,2,1,1);
	msmslayout->addWidget(msxnLabel,4,0,1,1);msmslayout->addWidget(msxnFilePath,4,1,1,1);msmslayout->addWidget(msxnButton,4,2,1,1);
	msmsBox->setLayout(msmslayout);

	// PDB2PQR
	QGroupBox* pBox=new QGroupBox(tr("PDB2PQR"));
	//
	QPushButton* pInfoButton=new QPushButton(tr("Installation"));
	pInfoButton->setFixedWidth(100);
	connect(pInfoButton,SIGNAL(clicked()),this,SLOT(launchPdb2pqrInstallInfo()));
	pFilePath=new QLineEdit;pFilePath->setFixedWidth(100);
	QPushButton* pButton=new QPushButton(tr("..."));pButton->setToolTip(tr("Specify the file path of the PDB2PQR 'main.py' python file."));
	pButton->setFixedWidth(40);
	connect(pButton,SIGNAL(clicked()),this,SLOT(openPDB2PQR()));
	pLabel=new QLabel;
	pLabel->setText("<span style=\"color:red;\">main.py</span>:");
	pLabel->setAlignment(Qt::AlignCenter|Qt::AlignRight);pLabel->setToolTip(tr("PDB2PQR Description"));
	// PDB2PQR Box Layout
	QGridLayout *playout=new QGridLayout;
	playout->addWidget(pInfoButton,0,1,1,2);
	playout->addWidget(pLabel,1,0,1,1);playout->addWidget(pFilePath,1,1,1,1);playout->addWidget(pButton,1,2,1,1);
	pBox->setLayout(playout);
	// ZPRED
	QGroupBox* zpredBox=new QGroupBox(tr("ZPRED"));
	QPushButton* zIInfoButton=new QPushButton(tr("Installation"));
	zIInfoButton->setFixedWidth(100);
	connect(zIInfoButton,SIGNAL(clicked()),this,SLOT(launchZpredInstallInfo()));

	zFilePath=new QLineEdit;zFilePath->setFixedWidth(100);
	QPushButton* zButton=new QPushButton(tr("..."));zButton->setToolTip(tr("Specify the file path of the ZPRED executable file."));
	zButton->setFixedWidth(40);
	connect(zButton,SIGNAL(clicked()),this,SLOT(openZPRED()));
	zLabel=new QLabel;
	zLabel->setText("<span style=\"color:red;\">ZPRED.exe</span>:");
	zLabel->setAlignment(Qt::AlignCenter|Qt::AlignRight);zLabel->setToolTip(tr("ZPRED Description"));
	// ZPRED Box Layout
	QGridLayout *zpredlayout=new QGridLayout;
	zpredlayout->addWidget(zIInfoButton,0,1,1,2);
	zpredlayout->addWidget(zLabel,1,0,1,1);zpredlayout->addWidget(zFilePath,1,1,1,1);zpredlayout->addWidget(zButton,1,2,1,1);
	zpredBox->setLayout(zpredlayout);

	connect(apbsButton,SIGNAL(clicked()),this,SLOT(apbsFilePath_defined()));
	connect(mvButton,SIGNAL(clicked()),this,SLOT(mvFilePath_defined()));
	connect(hButton,SIGNAL(clicked()),this,SLOT(hFilePath_defined()));
	connect(msButton,SIGNAL(clicked()),this,SLOT(msFilePath_defined()));
	connect(msaButton,SIGNAL(clicked()),this,SLOT(msaFilePath_defined()));
	connect(msxButton,SIGNAL(clicked()),this,SLOT(msxFilePath_defined()));
	connect(msxnButton,SIGNAL(clicked()),this,SLOT(msxnFilePath_defined()));
	connect(pButton,SIGNAL(clicked()),this,SLOT(pFilePath_defined()));
	connect(zButton,SIGNAL(clicked()),this,SLOT(zFilePath_defined()));

	QGridLayout *elayout = new QGridLayout;
	elayout->addWidget(autoFillButton,0,0,1,3);
	elayout->addWidget(apbsBox,1,0,1,3);
	elayout->addWidget(hydroBox,2,0,1,3);
	elayout->addWidget(msmsBox,3,0,1,3);
	elayout->addWidget(pBox,4,0,1,3);
	elayout->addWidget(zpredBox,5,0,1,3);
	eBox->setLayout(elayout);

	// Improved Executable Group Box
	QScrollArea* scrollArea=new QScrollArea();
	scrollArea->setWidget(eBox);  
	scrollArea->setWidgetResizable(true);  
	QGridLayout *exelayout=new QGridLayout;
	exelayout->addWidget(scrollArea);
	QGroupBox* exeBox=new QGroupBox(tr("|-PROGRAM FILE PATHS-|"));
	exeBox->setFont(font);
	exeBox->setLayout(exelayout);exeBox->setAlignment(Qt::AlignLeft);

	unsigned int number_of_threads = thread::hardware_concurrency();
	//std::cout << n << " concurrent threads are supported.\n";
	string numThreadValue;
	if(number_of_threads!=0)
		{numThreadValue=cnvrtNumToStrng(number_of_threads,0);}
	else
		{numThreadValue="4";}

	saveTempsBox=new QCheckBox(tr("Save Temporary Computational Files?"));saveTempsBox->setChecked(false);saveTempsBox->setFont(font);
	hydroproBox=new QCheckBox(tr("Use HYDROPRO instead of HULLRAD?"));hydroproBox->setChecked(false);hydroproBox->setFont(font);
	connect(hydroproBox,SIGNAL(stateChanged(int)),this,SLOT(show_hydroBox()));
	maxThreadsLine=new QLineEdit();
	maxThreadsLine->setAlignment(Qt::AlignCenter);
	maxThreadsLine->setText(tr(numThreadValue.c_str()));
	maxThreadsLine->setMaxLength(3);
	maxThreadsLine->setFixedWidth(80);
	maxThreadsLineLabel=new QLabel(tr(":<b><span style=\"color:red;\">Number of Simultaneous Computations</span></b>"));maxThreadsLineLabel->setAlignment(Qt::AlignCenter|Qt::AlignLeft);
	QPushButton *aOptButton=new QPushButton(tr("APBS\nOptions"));aOptButton->setFont(font);
	aOptButton->setFixedWidth(80);
	connect(aOptButton,SIGNAL(clicked()),this,SLOT(launchApbsOptions()));
	QPushButton *hOptButton=new QPushButton(tr("HYDROPRO\nOptions"));hOptButton->setFont(font);
	hOptButton->setFixedWidth(100);
	connect(hOptButton,SIGNAL(clicked()),this,SLOT(launchHydroproOptions()));
	QPushButton *mOptButton=new QPushButton(tr("MSMS\nOptions"));mOptButton->setFont(font);
	mOptButton->setFixedWidth(80);
	connect(mOptButton,SIGNAL(clicked()),this,SLOT(launchMsmsOptions()));
	QPushButton *pOptButton=new QPushButton(tr("PDB2PQR\nOptions"));pOptButton->setFont(font);
	pOptButton->setFixedWidth(100);
	connect(pOptButton,SIGNAL(clicked()),this,SLOT(launchPdb2pqrOptions()));
	connect(maxThreadsLine,SIGNAL(returnPressed()),this,SLOT(maxThreads_defined()));
	
	QPushButton *sInfoButton=new QPushButton(tr("Solution\nProperty\nComputation"));sInfoButton->setFont(font);
	connect(sInfoButton,SIGNAL(clicked()),this,SLOT(launchSolventInfo()));
	QPushButton *zInfoButton=new QPushButton(tr("How To\nUse\nZPRED"));zInfoButton->setFont(font);
	connect(zInfoButton,SIGNAL(clicked()),this,SLOT(launchZpredInfo()));
	QPushButton *vInfoButton=new QPushButton(tr("Experimental\nMethods and\nValidation"));vInfoButton->setFont(font);
	connect(vInfoButton,SIGNAL(clicked()),this,SLOT(launchValidationInfo()));

	QGroupBox* botLeft=new QGroupBox(tr("|-ADDITIONAL PARAMETERS-|"));
	QGridLayout *blayout=new QGridLayout;
	blayout->setSizeConstraint(QLayout::SetFixedSize);
	blayout->addWidget(saveTempsBox,0,0,1,2);blayout->addWidget(aOptButton,0,2,1,1);blayout->addWidget(hOptButton,0,3,1,1);
	blayout->addWidget(maxThreadsLine,1,0,1,1);blayout->addWidget(maxThreadsLineLabel,1,1,1,1);blayout->addWidget(mOptButton,1,2,1,1);blayout->addWidget(pOptButton,1,3,1,1);
	blayout->addWidget(hydroproBox,2,0,1,3);
	botLeft->setFont(font);
	botLeft->setLayout(blayout);botLeft->setAlignment(Qt::AlignCenter);

	QGroupBox* botCenter=new QGroupBox(tr("|-INFORMATION-|"));
	QGridLayout *bclayout=new QGridLayout;
	bclayout->setSizeConstraint(QLayout::SetFixedSize);
	bclayout->addWidget(sInfoButton,0,0);bclayout->addWidget(zInfoButton,0,1);bclayout->addWidget(vInfoButton,0,2);
	botCenter->setFont(font);
	botCenter->setLayout(bclayout);botCenter->setAlignment(Qt::AlignCenter);

	viewInputButton=new QPushButton(tr("View Input\nStructures")); viewInputButton->setFont(font); viewInputButton->hide();
	connect(viewInputButton,SIGNAL(clicked()),this,SLOT(launchPyMol()));
	QPushButton *calcButton=new QPushButton(tr("Calculate\nζ")); calcButton->setFont(font);
	connect(calcButton,SIGNAL(clicked()),this,SLOT(calcZetaPotential()));
	viewOutputButton=new QPushButton(tr("View Generated\nOutput")); viewOutputButton->setFont(font); viewOutputButton->hide();

	QGroupBox* botRight=new QGroupBox;
	QGridLayout *brlayout=new QGridLayout;
	//brlayout->setSizeConstraint(QLayout::SetFixedSize);
	brlayout->addWidget(viewInputButton,0,0);
	brlayout->addWidget(calcButton,1,0);
	brlayout->addWidget(viewOutputButton,2,0);
	botRight->setLayout(brlayout);

	QGridLayout *mainLayout = new QGridLayout;
	mainLayout->setSizeConstraint(QLayout::SetFixedSize);
	mainLayout->addWidget(gBox,0,0,1,2);
	mainLayout->addWidget(exeBox,0,2,1,1);
	mainLayout->addWidget(botLeft,1,0,1,1);
	mainLayout->addWidget(botCenter,1,1,1,1);
	mainLayout->addWidget(botRight,1,2,1,1);
	setLayout(mainLayout);
	setWindowTitle(tr("ZPRED: Zeta Potential (ζ) Prediction in a General Aqueous Electrolyte Solution"));
	resize(QDesktopWidget().availableGeometry(this).size());
}

void UI::show_fBox()
	{int tabIndex=scTab->currentIndex();
	QLabel* TUnitsLabel=new QLabel(tr("[Kelvin]"));TUnitsLabel->setAlignment(Qt::AlignLeft);
	QLabel* pCUnitsLabel=new QLabel(tr("[g/L]"));pCUnitsLabel->setAlignment(Qt::AlignLeft);
	//
	TLabel[tabIndex]->setFixedWidth(150);
	T[tabIndex]->setFixedWidth(LINE_INPUT_SIZE);
	pH[tabIndex]->setFixedWidth(LINE_INPUT_SIZE);
	pC[tabIndex]->setFixedWidth(LINE_INPUT_SIZE);
	// Solvent(s) Selection Group
	switch(numSolvent[tabIndex])
		{case 1: solventBox=create_oneSolventBox(tabIndex);break;
		case 2: solventBox=create_twoSolventBox(tabIndex);break;
		case 3: solventBox=create_threeSolventBox(tabIndex);break;
		case 4: solventBox=create_fourSolventBox(tabIndex);break;
		case 5: solventBox=create_fiveSolventBox(tabIndex);break;
		case 6: solventBox=create_sixSolventBox(tabIndex);break;
		case 7: solventBox=create_sevenSolventBox(tabIndex);break;
		case 8: solventBox=create_eightSolventBox(tabIndex);break;
		default: break;
		}

	// Solute(s) Selection Group
	switch(numSolute[tabIndex])
		{case 1: soluteBox=create_oneSoluteBox(tabIndex);break;
		case 2: soluteBox=create_twoSoluteBox(tabIndex);break;
		case 3: soluteBox=create_threeSoluteBox(tabIndex);break;
		case 4: soluteBox=create_fourSoluteBox(tabIndex);break;
		case 5: soluteBox=create_fiveSoluteBox(tabIndex);break;
		case 6: soluteBox=create_sixSoluteBox(tabIndex);break;
		case 7: soluteBox=create_sevenSoluteBox(tabIndex);break;
		case 8: soluteBox=create_eightSoluteBox(tabIndex);break;
		default: break;
		}

	// Define Solution Conditions Tab
	QWidget *tab1 = new QWidget;
	
	QGridLayout *tab1hbox = new QGridLayout;
	tab1hbox->setSizeConstraint(QLayout::SetFixedSize);
	tab1hbox->addWidget(TLabel[tabIndex],0,0,1,1);tab1hbox->addWidget(T[tabIndex],0,1,1,1);tab1hbox->addWidget(TUnitsLabel,0,2,1,1);tab1hbox->addWidget(pHLabel[tabIndex],0,4,1,1);tab1hbox->addWidget(pH[tabIndex],0,5,1,1);
	tab1hbox->addWidget(solventBox,1,0,1,3);tab1hbox->addWidget(soluteBox,1,3,1,3);
	if(forceParmBox[tabIndex]->isChecked())
		{tab1hbox->addWidget(forceParmBox[tabIndex],2,0,1,3);
		tab1hbox->addWidget(create_fBox(tabIndex),3,0,2,3);tab1hbox->addWidget(pCLabel[tabIndex],3,3,1,1);tab1hbox->addWidget(pC[tabIndex],3,4,1,1);tab1hbox->addWidget(pCUnitsLabel,3,5,1,1);
		tab1hbox->addWidget(potProfileBox[tabIndex],5,0,1,3);}
	else
		{tab1hbox->addWidget(forceParmBox[tabIndex],2,0,1,3);tab1hbox->addWidget(pCLabel[tabIndex],2,3,1,1);tab1hbox->addWidget(pC[tabIndex],2,4,1,1);tab1hbox->addWidget(pCUnitsLabel,2,5,1,1);
		tab1hbox->addWidget(potProfileBox[tabIndex],3,0,1,3);}
	tab1->setLayout(tab1hbox);
	// Rewrite tab
	string tmp="&"+cnvrtNumToStrng(tabIndex,0);
	// Clear Previous
	scTab->removeTab(tabIndex);
	scTab->insertTab(tabIndex,tab1,tmp.c_str());
	scTab->setCurrentIndex(tabIndex);
	}

void UI::show_hydroBox()
	{if(hydroproBox->isChecked()){hydroBox->show();}
	else{hydroBox->hide();}}


void UI::calcZetaPotential()
	{double *solvFrac,Sum,value;
	bool RUN_ZPRED=true;
	int pos;
	// Check for Sufficient Input While Gathering
	string delimiter=";",errMsg="<b><u>Error in User Interface Input!</u></b><br>";
	uiInput.delimiter=delimiter;
	uiInput.id=new string[numSolCond];
	string fLst="",sLst="",scLst="",sctLst="",tmp,bld,c;
	for(int i=0;i<numPDBs;i++){fLst+=pdbFiles[i]+delimiter;}
	uiInput.files=fLst;
	uiInput.outputTitle=jobName->text().toStdString();
	uiInput.pH=new string[numSolCond];
	uiInput.proteinConc=new string[numSolCond];
	uiInput.solvent=new string[numSolCond];
	uiInput.solventConc=new string[numSolCond];
	uiInput.solventConcType=new string[numSolCond];
	uiInput.solute=new string[numSolCond];
	uiInput.soluteConc=new string[numSolCond];
	uiInput.soluteConcType=new string[numSolCond];
	uiInput.temperature=new string[numSolCond];
	uiInput.fdensity=new string[numSolCond];
	uiInput.fdielectric=new string[numSolCond];
	uiInput.fviscosity=new string[numSolCond];
	uiInput.fXsp=new string[numSolCond];
	uiInput.genPotProfile=new string[numSolCond];
	uiInput.useHydropro=new string[numSolCond];
	for(int i=0;i<numSolCond;i++)
		{uiInput.id[i]=cnvrtNumToStrng(i,0);
		// Check pH
		tmp=pH[i]->text().toStdString();
		if(tmp.compare(NO_VALUE)==0){errMsg+="<b>Solution Condition #"+cnvrtNumToStrng(i,0)+"</b>: <b>pH</b> value is undefined.<br>";RUN_ZPRED=false;}
		else{uiInput.pH[i]=tmp;}
		// Check Protein Concentration
		tmp=pC[i]->text().toStdString();
		if(tmp.compare(NO_VALUE)==0){errMsg+="<b>Solution Condition #"+cnvrtNumToStrng(i,0)+"</b>: <b>Protien Concentration</b> is undefined.<br>";RUN_ZPRED=false;}
		else{uiInput.proteinConc[i]=tmp;}
		// Solvent
		sLst="";scLst="";sctLst="";
		for(int j=0;j<numSolvent[i];j++)
			{// Check Solvent Name
			tmp=solventLabel[i][j]->text().toStdString();
			pos=tmp.find(" ",0);
			if(pos!=string::npos){tmp.replace(pos,1,"_");}
			if(tmp.compare(NO_VALUE)==0){errMsg+="<b>Solution Condition #"+cnvrtNumToStrng(i,0)+"</b>: <b>Solvent Name</b> is undefined.<br>";RUN_ZPRED=false;}
			else{sLst+=tmp+",";}// use internal delimiter
			// Check Solvent Concentration/Fraction
			tmp=solventConcLabel[i][j]->text().toStdString();
			if(tmp.compare(NO_VALUE)==0){errMsg+="<b>Solution Condition #"+cnvrtNumToStrng(i,0)+"</b>: <b>Solvent Concentration</b> is undefined.<br>";RUN_ZPRED=false;}
			else{scLst+=tmp+",";}
			}
		uiInput.solvent[i]=sLst;
		// Check Solvent Fraction = 1
		solvFrac=fill_double_array(scLst,numSolvent[i],",");
		Sum=0;
		for(int j=0;j<numSolvent[i];j++){Sum+=solvFrac[j];}
		if(Sum!=1){errMsg+="<b>Solution Condition #"+cnvrtNumToStrng(i,0)+"</b>: <b>Solvent Fractions</b> do NOT add to one.<br>";RUN_ZPRED=false;}
		delete [] solvFrac;
		uiInput.solventConc[i]=scLst;
		tmp=solventConcTypeLabel[i]->text().toStdString();
		if(tmp.compare("Mass Fraction")==0){tmp="massFraction";}
		else if(tmp.compare("Mole Fraction")==0){tmp="moleFraction";}
		else if(tmp.compare("Volume Fraction")==0){tmp="volumeFraction";}
		else if(tmp.compare(NO_VALUE)==0){errMsg+="<b>Solution Condition #"+cnvrtNumToStrng(i,0)+"</b>: <b>Solvent Fraction Units</b> are undefined.<br>";RUN_ZPRED=false;}
		//uiInput.solventConcType[i]=sctLst;
		uiInput.solventConcType[i]=tmp;
		// Solute
		sLst="";scLst="";sctLst="";
		for(int j=0;j<numSolute[i];j++)
			{// Check Solute Name
			tmp=fSoluteLabel[i][j]->text().toStdString();
			if(tmp.compare(NO_VALUE)==0){errMsg+="<b>Solution Condition #"+cnvrtNumToStrng(i,0)+"</b>: <b>Solute Name</b> is undefined.<br>";RUN_ZPRED=false;}
			else{sLst+=tmp+",";}// use internal delimiter
			// Check Solute Concentration
			tmp=soluteConcLabel[i][j]->text().toStdString();
			if(tmp.compare(NO_VALUE)==0){errMsg+="<b>Solution Condition #"+cnvrtNumToStrng(i,0)+"</b>: <b>Solute Concentration</b> is undefined.<br>";RUN_ZPRED=false;}
			else{scLst+=tmp+",";}
			}
		uiInput.solute[i]=sLst;
		uiInput.soluteConc[i]=scLst;
		tmp=soluteConcTypeLabel[i]->text().toStdString();
		//if(tmp.compare("Mass Fraction")==0){tmp="massFraction";}
		//else if(tmp.compare("Mole Fraction")==0){tmp="moleFraction";}
		//else if(tmp.compare("Volume Fraction")==0){tmp="volumeFraction";}
		//else if(tmp.compare(NO_VALUE)==0){cerr<<"Error in Interface User Input!\nSolute Concentration Type for Solution Condition #"<<i<<" is undefined.\n";RUN_ZPRED=false;return;}
		//uiInput.soluteConcType[i]=sctLst;
		uiInput.soluteConcType[i]=tmp;
		// Check Temperature
		tmp=T[i]->text().toStdString();
		if(tmp.compare(NO_VALUE)==0){errMsg+="<b>Solution Condition #"+cnvrtNumToStrng(i,0)+"</b>: <b>Temperature</b> is undefined.<br>";RUN_ZPRED=false;}
		else{value=strtod(tmp.c_str(),NULL);
			if(value<273.15)
				{errMsg+="<b>Solution Condition #"+cnvrtNumToStrng(i,0)+"</b>: <b>Temperature</b> is below freezing.<br>";RUN_ZPRED=false;}
			else if(value>373.15)
				{errMsg+="<b>Solution Condition #"+cnvrtNumToStrng(i,0)+"</b>: <b>Temperature</b> is above boiling.<br>";RUN_ZPRED=false;}
			else
				{uiInput.temperature[i]=tmp;}
			}
		// Forced Values (if specified)
		// Density
		tmp=fDens[i]->text().toStdString();
		if(tmp.compare(NO_VALUE)==0){uiInput.fdensity[i]="-1";}
		else{uiInput.fdensity[i]=tmp;}
		// Relative Dielectric
		tmp=fDiel[i]->text().toStdString();
		if(tmp.compare(NO_VALUE)==0){uiInput.fdielectric[i]="-1";}
		else{uiInput.fdielectric[i]=tmp;}
		// Viscosity
		tmp=fVisc[i]->text().toStdString();
		if(tmp.compare(NO_VALUE)==0){uiInput.fviscosity[i]="-1";}
		else{uiInput.fviscosity[i]=tmp;}
		// Hydration Layer Thickness/Slip Plane Position
		tmp=fXsp[i]->text().toStdString();
		if(tmp.compare(NO_VALUE)==0){uiInput.fXsp[i]="-1";}
		else{uiInput.fXsp[i]=tmp;}
		// Electric Potential Profile Generation CheckBox
		if(potProfileBox[i]->isChecked()){uiInput.genPotProfile[i]="true";}
		else{uiInput.genPotProfile[i]="false";}
		}
	// Check APBS Executable
	tmp=apbsFilePath->text().toStdString();
	if(tmp.compare(NO_VALUE)==0||tmp.compare("")==0){errMsg+="<b>APBS Executable File Path</b> is undefined.<br>";RUN_ZPRED=false;}
	else{uiInput.apbs=tmp;}
	// Check APBS Multivalue
	tmp=mvFilePath->text().toStdString();
	if(tmp.compare(NO_VALUE)==0||tmp.compare("")==0){errMsg+="<b>APBS Tool (multivalue) File Path</b> is undefined.<br>";RUN_ZPRED=false;}
	else{uiInput.multivalue=tmp;}
	// Check HYDROPRO Executable
	if(hydroproBox->isChecked())
		{tmp=hFilePath->text().toStdString();
		if(tmp.compare(NO_VALUE)==0||tmp.compare("")==0){errMsg+="<b>HYDROPRO Executable File Path</b> is undefined.<br>";RUN_ZPRED=false;}
		else{uiInput.hydropro=tmp;}
		}
	// Check MSMS Executable
	tmp=msFilePath->text().toStdString();
	if(tmp.compare(NO_VALUE)==0||tmp.compare("")==0){errMsg+="<b>MSMS Executable File Path</b> is undefined.<br>";RUN_ZPRED=false;}
	else{uiInput.msms=tmp;}
	// Check MSMS atmtypenumbers File
	tmp=msaFilePath->text().toStdString();
	if(tmp.compare(NO_VALUE)==0||tmp.compare("")==0){errMsg+="<b>MSMS atmtypenumbers File Path</b> is undefined.<br>";RUN_ZPRED=false;}
	else{uiInput.msmsAtmTypeNumbers=tmp;}
	// Check MSMS pdb_to_xyzr File
	tmp=msxFilePath->text().toStdString();
	if(tmp.compare(NO_VALUE)==0||tmp.compare("")==0){errMsg+="<b>MSMS pdb_to_xyzr File Path</b> is undefined.<br>";RUN_ZPRED=false;}
	else{uiInput.msmsPdbToXyzr=tmp;}
	// Check MSMS pdb_to_xyzrn File
	tmp=msxnFilePath->text().toStdString();
	if(tmp.compare(NO_VALUE)==0||tmp.compare("")==0){errMsg+="<b>MSMS pdb_to_xyzrn File Path</b> is undefined.<br>";RUN_ZPRED=false;}
	else{uiInput.msmsPdbToXyzrn=tmp;}
	// Check PDB2PQR main.py File
	tmp=pFilePath->text().toStdString();
	if(tmp.compare(NO_VALUE)==0||tmp.compare("")==0){errMsg+="<b>PDB2PQR main.py File Path</b> is undefined.<br>";RUN_ZPRED=false;}
	else{uiInput.pdb2pqr=tmp;}
	// Check ZRPED Executable
	tmp=zFilePath->text().toStdString();
	if(tmp.compare(NO_VALUE)==0||tmp.compare("")==0){errMsg+="<b>ZPRED Executable File Path</b> is undefined.<br>";RUN_ZPRED=false;}
	else{uiInput.zpred=tmp;}
	// Check Maximum Threads
	tmp=maxThreadsLine->text().toStdString();
	if(tmp.compare(NO_VALUE)==0){errMsg+="<b>Number of Simultaneous Computations</b> is undefined.<br>";RUN_ZPRED=false;}
	else{uiInput.maxThreads=tmp;}	
	// Check Optional Parameters File Reference (.optionalParameters.txt)	
	ifstream fIn;
	fIn.open(optRefFile.c_str());
	if(fIn.fail()){cerr<<"ERROR in calcZetaPotential!\nOptional Parameters file could not be opened.\n"<<optRefFile<<endl;exit(EXIT_FAILURE);}
	uint Sz=15000;
	char Val[Sz];
	string inputName,inputVal;
	uiInput.apbsWriteTypes="";
	fIn.getline(Val,Sz);
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
				if(inputVal.length()!=0)
					{// Identify Parameter and Load Values
					if(inputName.compare("PDIE")==0)
						{// Protein Dielectric Constant
						uiInput.pdie=inputVal;}
					else if(inputName.compare("DX")==0)
						{// Num X Grid Pnts
						uiInput.dim_x=inputVal;}
					else if(inputName.compare("DY")==0)
						{// Num Y Grid Pnts
						uiInput.dim_y=inputVal;}
					else if(inputName.compare("DZ")==0)
						{// Num Z Grid Pnts
						uiInput.dim_z=inputVal;}
					else if(inputName.compare("ATOMPOT")==0)
						{// Write Type
						if(inputVal.compare("true")==0){uiInput.apbsWriteTypes+="atompot;";}
						else if(inputVal.compare("false")==0){}
						else{cerr<<"ERROR in calcZetaPotential!\nUnrecognized Optional Parameter ("<<inputName<<") Input:\n"<<inputVal<<endl;exit(EXIT_FAILURE);}
						}
					else if(inputName.compare("CHARGE")==0)
						{// Write Type
						if(inputVal.compare("true")==0){uiInput.apbsWriteTypes+="charge;";}
						else if(inputVal.compare("false")==0){}
						else{cerr<<"ERROR in calcZetaPotential!\nUnrecognized Parameter ("<<inputName<<") Input:\n"<<inputVal<<endl;exit(EXIT_FAILURE);}
						}
					else if(inputName.compare("DIELX")==0)
						{// Write Type
						if(inputVal.compare("true")==0){uiInput.apbsWriteTypes+="dielx;";}
						else if(inputVal.compare("false")==0){}
						else{cerr<<"ERROR in calcZetaPotential!\nUnrecognized Parameter ("<<inputName<<") Input:\n"<<inputVal<<endl;exit(EXIT_FAILURE);}
						}
					else if(inputName.compare("DIELY")==0)
						{// Write Type
						if(inputVal.compare("true")==0){uiInput.apbsWriteTypes+="diely;";}
						else if(inputVal.compare("false")==0){}
						else{cerr<<"ERROR in calcZetaPotential!\nUnrecognized Parameter ("<<inputName<<") Input:\n"<<inputVal<<endl;exit(EXIT_FAILURE);}
						}
					else if(inputName.compare("DIELZ")==0)
						{// Write Type
						if(inputVal.compare("true")==0){uiInput.apbsWriteTypes+="dielz;";}
						else if(inputVal.compare("false")==0){}
						else{cerr<<"ERROR in calcZetaPotential!\nUnrecognized Parameter ("<<inputName<<") Input:\n"<<inputVal<<endl;exit(EXIT_FAILURE);}
						}
					else if(inputName.compare("EDENS")==0)
						{// Write Type
						if(inputVal.compare("true")==0){uiInput.apbsWriteTypes+="edens;";}
						else if(inputVal.compare("false")==0){}
						else{cerr<<"ERROR in calcZetaPotential!\nUnrecognized Parameter ("<<inputName<<") Input:\n"<<inputVal<<endl;exit(EXIT_FAILURE);}
						}
					else if(inputName.compare("IVDW")==0)
						{// Write Type
						if(inputVal.compare("true")==0){uiInput.apbsWriteTypes+="ivdw;";}
						else if(inputVal.compare("false")==0){}
						else{cerr<<"ERROR in calcZetaPotential!\nUnrecognized Parameter ("<<inputName<<") Input:\n"<<inputVal<<endl;exit(EXIT_FAILURE);}
						}
					else if(inputName.compare("KAPPA")==0)
						{// Write Type
						if(inputVal.compare("true")==0){uiInput.apbsWriteTypes+="kappa;";}
						else if(inputVal.compare("false")==0){}
						else{cerr<<"ERROR in calcZetaPotential!\nUnrecognized Parameter ("<<inputName<<") Input:\n"<<inputVal<<endl;exit(EXIT_FAILURE);}
						}
					else if(inputName.compare("LAP")==0)
						{// Write Type
						if(inputVal.compare("true")==0){uiInput.apbsWriteTypes+="lap;";}
						else if(inputVal.compare("false")==0){}
						else{cerr<<"ERROR in calcZetaPotential!\nUnrecognized Parameter ("<<inputName<<") Input:\n"<<inputVal<<endl;exit(EXIT_FAILURE);}
						}
					else if(inputName.compare("NDENS")==0)
						{// Write Type
						if(inputVal.compare("true")==0){uiInput.apbsWriteTypes+="ndens;";}
						else if(inputVal.compare("false")==0){}
						else{cerr<<"ERROR in calcZetaPotential!\nUnrecognized Parameter ("<<inputName<<") Input:\n"<<inputVal<<endl;exit(EXIT_FAILURE);}
						}
					else if(inputName.compare("POT")==0)
						{// Write Type
						if(inputVal.compare("true")==0){uiInput.apbsWriteTypes+="pot;";}
						else if(inputVal.compare("false")==0){}
						else{cerr<<"ERROR in calcZetaPotential!\nUnrecognized Parameter ("<<inputName<<") Input:\n"<<inputVal<<endl;exit(EXIT_FAILURE);}
						}
					else if(inputName.compare("QDENS")==0)
						{// Write Type
						if(inputVal.compare("true")==0){uiInput.apbsWriteTypes+="qdens;";}
						else if(inputVal.compare("false")==0){}
						else{cerr<<"ERROR in calcZetaPotential!\nUnrecognized Parameter ("<<inputName<<") Input:\n"<<inputVal<<endl;exit(EXIT_FAILURE);}
						}
					else if(inputName.compare("SMOL")==0)
						{// Write Type
						if(inputVal.compare("true")==0){uiInput.apbsWriteTypes+="smol;";}
						else if(inputVal.compare("false")==0){}
						else{cerr<<"ERROR in calcZetaPotential!\nUnrecognized Parameter ("<<inputName<<") Input:\n"<<inputVal<<endl;exit(EXIT_FAILURE);}
						}
					else if(inputName.compare("SSPL")==0)
						{// Write Type
						if(inputVal.compare("true")==0){uiInput.apbsWriteTypes+="sspl;";}
						else if(inputVal.compare("false")==0){}
						else{cerr<<"ERROR in calcZetaPotential!\nUnrecognized Parameter ("<<inputName<<") Input:\n"<<inputVal<<endl;exit(EXIT_FAILURE);}
						}
					else if(inputName.compare("VDW")==0)
						{// Write Type
						if(inputVal.compare("true")==0){uiInput.apbsWriteTypes+="vdw;";}
						else if(inputVal.compare("false")==0){}
						else{cerr<<"ERROR in calcZetaPotential!\nUnrecognized Parameter ("<<inputName<<") Input:\n"<<inputVal<<endl;exit(EXIT_FAILURE);}
						}
					else if(inputName.compare("CALCTYPE")==0)
						{// HYDROPRO Calculation Type
						uiInput.calcType=inputVal;
						}
					else if(inputName.compare("LOWPNTDENSITY")==0)
						{// MSMS Low Point Surface Density
						uiInput.lowPntDensity=inputVal;
						}
					else if(inputName.compare("HIGHPNTDENSITY")==0)
						{// MSMS HI Point Surface Density
						uiInput.highPntDensity=inputVal;
						}
					else if(inputName.compare("FORCEFIELD")==0)
						{// PROPKA Force Field
						uiInput.forceField=inputVal;
						}
					//else{cerr<<"ERROR in apbsOptions::defineParameters!\nUnrecognized Parameter ("<<inputName<<")."<<endl;exit(EXIT_FAILURE);}
					}
				}
			}
		fIn.getline(Val,Sz);}
	fIn.close();
	// Pass Check Means Celebrate with Chariots of Fire
	QString theCurFldr=QDir::currentPath();
	string sFldr=theCurFldr.toStdString()+"/";
	string audioFile=sFldr+"UI/Images/chariots_of_fire.wav";
	QProcess qprocess;
	QStringList sL,zL;
	sL.append("-q");
	sL.append(audioFile.c_str());
	// Write ZPRED Input File per Solution Condition
	// Determine ZPRED Executable Folder
	string s=uiInput.zpred;
	pos=s.rfind("/",s.length()-1);
	string Fldr=s.substr(0,pos)+"/";
	string zInputFile=Fldr+"zpredInput.txt";
	//string cmd="/usr/bin/gnome-terminal -- \'"+uiInput.zpred+"\'";
	//string cmd=uiInput.zpred;//
	//zL.append("-i");
	//zL.append(zInputFile.c_str());
	string cmd=uiInput.zpred+" -i "+zInputFile;
	int statusCode,status;
	bool UPDATING=true;
	if(RUN_ZPRED)
		{writeZPREDInputFile(zInputFile);
		qprocess.startDetached("play",sL);
		//zpredPID=fork();
		//if(zpredPID==0)
		//	{//statusCode=system(cmd.c_str());
			//qprocess.startDetached(cmd.c_str(),zL);
		//	exit(EXIT_SUCCESS);}
		//else if(zpredPID>0)
		//	{launchZpredProgressDisplay();
		//	sleep(5);
			//for(int i=0;i<100000;i++){ZPD->updateDisplay();}
			//while(UPDATING)
			//	{ZPD->updateDisplay();
			//	sleep(10);}
		//	ZPD->updateDisplay();
		//	ZPD->show();
			//exit(EXIT_FAILURE);
			//waitpid(zpredPID,&status,0);
		//	}
		//else if(zpredPID==-1){cerr<<"Fork failed.\n";exit(EXIT_FAILURE);}
		//else{cerr<<"Output of fork() unpredicted:\n"<<zpredPID<<endl;}
		statusCode=system(cmd.c_str());
		//cerr<<"Status: "<<statusCode<<endl;
		viewOutputButton->show();
		}
	else
		{launchMsgBox(errMsg);}
	}

void UI::writeZPREDInputFile(string oFile)
	{string zIn="------------------------ZPRED Input----------------------------\n";
	string idLst="",pHLst="",pCLst="",solvLst="",solvCLst="",solvCTLst="",soluLst="",soluCLst="",soluCTLst="",tLst="",densLst="",dielLst="",viscLst="",xspLst="",gppLst="";
	// , = internal delimiter| ; = external delimiter of conditions
	for(int i=0;i<numSolCond;i++)
		{idLst+=uiInput.id[i]+";";
		pHLst+=uiInput.pH[i]+";";
		pCLst+=uiInput.proteinConc[i]+";";
		solvLst+=uiInput.solvent[i]+";";				// internal delimiter already included
		solvCLst+=uiInput.solventConc[i]+";";		// internal delimiter already included
		solvCTLst+=uiInput.solventConcType[i]+";";// internal delimiter already included
		soluLst+=uiInput.solute[i]+";";				// internal delimiter already included
		soluCLst+=uiInput.soluteConc[i]+";";		// internal delimiter already included
		soluCTLst+=uiInput.soluteConcType[i]+";";	// internal delimiter already included
		tLst+=uiInput.temperature[i]+";";
		densLst+=uiInput.fdensity[i]+";";
		dielLst+=uiInput.fdielectric[i]+";";
		viscLst+=uiInput.fviscosity[i]+";";
		xspLst+=uiInput.fXsp[i]+";";
		gppLst+=uiInput.genPotProfile[i]+";";}
	zIn+="// Solution Condition Number(s)\n";
	zIn+=".id: "+idLst+"\n";
	zIn+="// PDB Structure File(s)\n";
	zIn+=".files: "+uiInput.files+"\n";
	zIn+="// Job Name\n";
	zIn+=".outputTitle: "+uiInput.outputTitle+"\n";
	zIn+="// SOLUTION PARAMETERS\n";
	zIn+=".pH: "+pHLst+"\n";
	zIn+=".proteinConc: "+pCLst+"\n";
	zIn+=".solvent: "+solvLst+"\n";
	zIn+=".solventConc: "+solvCLst+"\n";
	zIn+=".solventConcType: "+solvCTLst+"\n";
	zIn+=".solute: "+soluLst+"\n";
	zIn+=".soluteConc: "+soluCLst+"\n";
	zIn+=".soluteConcType: "+soluCTLst+"\n";
	zIn+=".temperature: "+tLst+"\n";
	zIn+="// Program File Paths\n";
	zIn+=".apbsExecutable: "+uiInput.apbs+"\n";
	zIn+=".multivalueExecutable: "+uiInput.multivalue+"\n";
	if(hydroproBox->isChecked()){zIn+=".hydroproExecutable: "+uiInput.hydropro+"\n";}
	zIn+=".RUN_HYDROPRO: ";
	if(hydroproBox->isChecked()){zIn+="true\n";}
	else{zIn+="false\n";}
	zIn+=".RUN_HULLRAD: ";
	if(hydroproBox->isChecked()){zIn+="false\n";}
	else{zIn+="true\n";}
	zIn+=".msmsExecutable: "+uiInput.msms+"\n";
	zIn+=".msmsAtmTypeNumbers: "+uiInput.msmsAtmTypeNumbers+"\n";
	zIn+=".msmsPdbToXyzr: "+uiInput.msmsPdbToXyzr+"\n";
	zIn+=".msmsPdbToXyzrn: "+uiInput.msmsPdbToXyzrn+"\n";
	zIn+=".pdb2pqr: "+uiInput.pdb2pqr+"\n";
	zIn+="// Number of Computations to be Performed in Parallel\n";
	zIn+=".maxThreads: "+uiInput.maxThreads+"\n";
	zIn+="// Save Temporary Computational Files?\n";
	zIn+=".saveTemps: ";
	if(saveTempsBox->isChecked()){zIn+="true\n";}
	else{zIn+="false\n";}
	zIn+="// User-Specified Solution Parameters\n";
	zIn+=".density: "+densLst+"\n";
	zIn+=".dielectric: "+dielLst+"\n";
	zIn+=".viscosity: "+viscLst+"\n";
	zIn+=".inflationDistance: "+xspLst+"\n";
	zIn+="// Generate Electric Potential Profile?\n";
	zIn+=".generateProfile: "+gppLst+"\n";
	zIn+="// APBS OPTIONAL PARAMETERS\n";
	zIn+="// Protein Relative Dielectric\n";
	zIn+=".pdie: "+uiInput.pdie+"\n";
	zIn+="// Number of X Axis Grid Points\n";
	zIn+=".dime_x: "+uiInput.dim_x+"\n";
	zIn+="// Number of Y Axis Grid Points\n";
	zIn+=".dime_y: "+uiInput.dim_y+"\n";
	zIn+="// Number of Z Axis Grid Points\n";
	zIn+=".dime_z: "+uiInput.dim_z+"\n";
	zIn+="// Output Write Types\n";
	zIn+=".apbsWriteTypes: "+uiInput.apbsWriteTypes+"\n";
	zIn+="// HYDROPRO OPTIONAL PARAMETERS\n";
	zIn+=".hydroproCalcType: "+uiInput.calcType+"\n";
	zIn+="// MSMS OPTIONAL PARAMETERS\n";
	zIn+="// Low Point Surface Density\n";
	zIn+=".msmsSurfPointDensity: "+uiInput.lowPntDensity+"\n";
	zIn+="// High Point Surface Density\n";
	zIn+=".msmsSurfPointHiDensity: "+uiInput.highPntDensity+"\n";
	zIn+="// PDB2PQR OPTIONAL PARAMETERS\n";
	zIn+="// PROPKA Force Field\n";
	zIn+=".forceField: "+uiInput.forceField+"\n";
	zIn+="---------------------------------------------------------------\n";
	ofstream fOut;
	fOut.open(oFile.c_str(),ofstream::out|ofstream::trunc);
	if(fOut.fail()){cerr<<"ERROR in UI::writeZPREDInputFile(string oFile)!\nZPRED input file could not be written.\n"<<oFile<<endl;exit(EXIT_FAILURE);}
	fOut<<zIn;fOut.close();}

void UI::updateHistoryFile(string prgmNm,string Fldr)
	{// prgmNm is any: A=APBS | H=HYDROPRO | M=MSMS | P=PDB2PQR | Z=ZPRED
	// File Reference for When Downloads already Performed (makes .history.txt)
	//refFile=getBaseFolder(Fldr.toStdString())+".history.txt";
	ifstream fIn;
	ofstream fOut;
	fIn.open(refFile.c_str());
	uint Sz=15000,pos;
	char Val[Sz];
	string cpy="",tmp,inputName,inputVal;
	if(fIn.fail())
		{cout<<"Failed to open .history.txt\n";
		cpy=".apbsFldr:\n.filesFldr:\n.hydroproFldr:\n.msmsFldr:\n.pdb2pqrFldr:\n.zpredFldr:\n";
		fOut.open(refFile.c_str(),ofstream::out|ofstream::trunc);
		if(fOut.fail()){cerr<<"ERROR in updateHistoryFile!!!\nFile could not be opened.\n"<<refFile<<endl;exit(EXIT_FAILURE);}
		fOut<<cpy;fOut.close();
		// Reattempt
		fIn.open(refFile.c_str());
		// Copy & Update File
		fIn.getline(Val,Sz);
		cpy="";
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
					// Identify Parameter and Load Values
					if(inputName.compare("APBSFLDR")==0)
						{// APBS Download File Location
						if(prgmNm.compare("A")==0){cpy+=".apbsFldr:"+Fldr+"\n";}
						else{cpy+=".apbsFldr:"+inputVal+"\n";}
						}
					else if(inputName.compare("FILESFLDR")==0)
						{// HYDROPRO Download File Location
						if(prgmNm.compare("F")==0){cpy+=".filesFldr:"+Fldr+"\n";}
						else{cpy+=".filesFldr:"+inputVal+"\n";}
						}
					else if(inputName.compare("HYDROPROFLDR")==0)
						{// HYDROPRO Download File Location
						if(prgmNm.compare("H")==0){cpy+=".hydroproFldr:"+Fldr+"\n";}
						else{cpy+=".hydroproFldr:"+inputVal+"\n";}
						}
					else if(inputName.compare("MSMSFLDR")==0)
						{// MSMS Download File Location
						if(prgmNm.compare("M")==0){cpy+=".msmsFldr:"+Fldr+"\n";}
						else{cpy+=".msmsFldr:"+inputVal+"\n";}
						}
					else if(inputName.compare("PDB2PQRFLDR")==0)
						{// HYDROPRO Download File Location
						if(prgmNm.compare("P")==0){cpy+=".pdb2pqrFldr:"+Fldr+"\n";}
						else{cpy+=".pdb2pqrFldr:"+inputVal+"\n";}
						}
					else if(inputName.compare("ZPREDFLDR")==0)
						{// HYDROPRO Download File Location
						if(prgmNm.compare("Z")==0){cpy+=".zpredFldr:"+Fldr+"\n";}
						else{cpy+=".zpredFldr:"+inputVal+"\n";}
						}
					else{cerr<<"ERROR in updateHistoryFile!\nUnrecognized Parameter ("<<inputName<<")."<<endl;/*exit(EXIT_FAILURE);*/}
					}
				}
			fIn.getline(Val,Sz);}
		fIn.close();
		}
	else
		{// Copy & Update File
		//cout<<"Opened .history.txt\n";
		fIn.getline(Val,Sz);
		//cout<<Val<<endl;
		while(!fIn.eof())
			{tmp=Val;
			//cout<<tmp<<endl;
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
					//cout<<inputName<<" "<<inputVal<<endl;
					// Identify Parameter and Load Values
					if(inputName.compare("APBSFLDR")==0)
						{// APBS Download File Location
						if(prgmNm.compare("A")==0){cpy+=".apbsFldr:"+Fldr+"\n";}
						else{cpy+=".apbsFldr:"+inputVal+"\n";}
						}
					else if(inputName.compare("FILESFLDR")==0)
						{// HYDROPRO Download File Location
						if(prgmNm.compare("F")==0){cpy+=".filesFldr:"+Fldr+"\n";}
						else{cpy+=".filesFldr:"+inputVal+"\n";}
						}
					else if(inputName.compare("HYDROPROFLDR")==0)
						{// HYDROPRO Download File Location
						if(prgmNm.compare("H")==0){cpy+=".hydroproFldr:"+Fldr+"\n";}
						else{cpy+=".hydroproFldr:"+inputVal+"\n";}
						}
					else if(inputName.compare("MSMSFLDR")==0)
						{// MSMS Download File Location
						if(prgmNm.compare("M")==0){cpy+=".msmsFldr:"+Fldr+"\n";}
						else{cpy+=".msmsFldr:"+inputVal+"\n";}
						}
					else if(inputName.compare("PDB2PQRFLDR")==0)
						{// HYDROPRO Download File Location
						if(prgmNm.compare("P")==0){cpy+=".pdb2pqrFldr:"+Fldr+"\n";}
						else{cpy+=".pdb2pqrFldr:"+inputVal+"\n";}
						}
					else if(inputName.compare("ZPREDFLDR")==0)
						{// HYDROPRO Download File Location
						if(prgmNm.compare("Z")==0){cpy+=".zpredFldr:"+Fldr+"\n";}
						else{cpy+=".zpredFldr:"+inputVal+"\n";}
						}
					else{cerr<<"ERROR in updateHistoryFile!\nUnrecognized Parameter ("<<inputName<<")."<<endl;exit(EXIT_FAILURE);}
					}
				}
			fIn.getline(Val,Sz);}
		fIn.close();
		}
	//cout<<cpy<<endl;
	fOut.open(refFile.c_str(),ofstream::out|ofstream::trunc);
	if(fOut.fail()){cerr<<"ERROR in updateHistoryFile!!!\nFile could not be opened.\n"<<refFile<<endl;/*exit(EXIT_FAILURE);*/}
	else{fOut<<cpy;fOut.close();}
	//cout<<"Output complete.\n";
	}

void UI::defineDownloadFolders(string iFile)
	{ifstream fIn;
	fIn.open(iFile.c_str());
	if(fIn.fail()){cerr<<"ERROR in defineDownloadFolders!\nInput file could not be opened.\n"<<iFile<<endl;exit(EXIT_FAILURE);}
	uint Sz=15000,pos;
	char Val[Sz];
	string bld,tmp,inputName,inputVal;
	fIn.getline(Val,Sz);
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
				if(inputVal.length()!=0)
					{// Identify Parameter and Load Values
					if(inputName.compare("APBSFLDR")==0)
						{// APBS Download File Location
						apbsFldr=inputVal;}
					else if(inputName.compare("FILESFLDR")==0)
						{// PDB Files File Location
						filesFldr=inputVal;}
					else if(inputName.compare("HYDROPROFLDR")==0)
						{// HYDROPRO Download File Location
						hydroproFldr=inputVal;}
					else if(inputName.compare("MSMSFLDR")==0)
						{// MSMS Download File Location
						msmsFldr=inputVal;}
					else if(inputName.compare("PDB2PQRFLDR")==0)
						{// HYDROPRO Download File Location
						pdb2pqrFldr=inputVal;}
					else if(inputName.compare("ZPREDFLDR")==0)
						{// HYDROPRO Download File Location
						zpredFldr=inputVal;}
					else{cerr<<"ERROR in defineDownloadFolders!\nUnrecognized Parameter ("<<inputName<<")."<<endl;exit(EXIT_FAILURE);}
					}
				}
			}
		fIn.getline(Val,Sz);}
	fIn.close();
	}

void UI::launchMsgBox(string errTxt)
	{QString Fldr=QDir::currentPath();
	string sFldr=Fldr.toStdString()+"/";
	// File Reference for When Optional Parameters are set (makes .optionalParameters.txt)
	string msgFile=sFldr+".message.txt";
	ofstream fOut;
	fOut.open(msgFile.c_str(),ofstream::out|ofstream::trunc);fOut<<errTxt;fOut.close();
	MB=new msgBox();
	MB->show();}

void UI::launchMethodsInfo()
	{MI=new methodsInfo();
	MI->show();}

void UI::launchPyMol()
	{QString theCurFldr=QDir::currentPath();
	string sFldr=theCurFldr.toStdString()+"/";
	string pythonFunctionName="viewInput";
	string pythonFunctionFile=sFldr+"viewInput.py";
	string cmd="pymol -d run "+pythonFunctionFile+" -d "+pythonFunctionName;
	//writeZPREDInputFile(zInputFile);
	write_pymolViewInputFunction_File(pythonFunctionName,pythonFunctionFile);	
	int statusCode=system(cmd.c_str());
	//cerr<<"Status: "<<statusCode<<endl;
	}

void UI::launchSolventInfo()
	{SI=new solventInfo();
	SI->show();}

void UI::launchValidationInfo()
	{MI=new methodsInfo();
	MI->show();}

void UI::launchZpredInfo()
	{ZI=new zpredInfo();
	ZI->show();}

void UI::launchZpredInstallInfo()
	{ZII=new zpredInstallInfo();
	ZII->show();}

void UI::launchZpredProgressDisplay()
	{ZPD=new zpredProgressDisplay();
	QString Fldr=QDir::currentPath();
	string sFldr=Fldr.toStdString()+"/";
	ZPD->displayFldr=sFldr+jobName->text().toStdString()+"/";
	ZPD->PID=zpredPID;
	ZPD->show();}

void UI::launchApbsOptions()
	{AO=new apbsOptions();
	AO->show();}

void UI::launchApbsInstallInfo()
	{AII=new apbsInstallInfo();
	AII->show();}

void UI::launchHydroproOptions()
	{HO=new hydroproOptions();
	HO->show();}

void UI::launchHydroproInstallInfo()
	{HII=new hydroproInstallInfo();
	HII->show();}

void UI::launchMsmsOptions()
	{MO=new msmsOptions();
	MO->show();}

void UI::launchMsmsInstallInfo()
	{MII=new msmsInstallInfo();
	MII->show();}

void UI::launchPdb2pqrOptions()
	{PO=new pdb2pqrOptions();
	PO->show();}

void UI::launchPdb2pqrInstallInfo()
	{PII=new pInstallInfo();
	PII->show();}

void UI::maxThreads_defined(){maxThreadsLineLabel->setText(tr(":<b><span style=\"color:black;\">Number of Simultaneous Computations</span></b>"));}
void UI::pH_defined(){int tabIndex=scTab->currentIndex();pHLabel[tabIndex]->setText(tr("<b><span style=\"color:black;\">pH</span></b>:"));}
void UI::T_defined(){int tabIndex=scTab->currentIndex();TLabel[tabIndex]->setText(tr("<b><span style=\"color:black;\">Temperature</span></b>:"));}
void UI::pC_defined(){int tabIndex=scTab->currentIndex();pCLabel[tabIndex]->setText(tr("<b><span style=\"color:black;\">Protein Concentration</span></b>:"));}
void UI::solvent_defined(const QString&){int tabIndex=scTab->currentIndex();seleSolvLabel[tabIndex]->setText(tr("<b><span style=\"color:black;\">Select a Solvent</span></b>"));}
void UI::solventConc_defined(){int tabIndex=scTab->currentIndex();seleSolvConcLabel[tabIndex]->setText(tr("<b><span style=\"color:black;\">Fraction</span></b>"));}
void UI::solventConcType_defined(const QString&){int tabIndex=scTab->currentIndex();seleSolvConcTypeLabel[tabIndex]->setText(tr("<b><span style=\"color:black;\">Units</span></b>"));}
void UI::solute_defined(const QString&){int tabIndex=scTab->currentIndex();seleSoluLabel[tabIndex]->setText(tr("<b><span style=\"color:black;\">Select a Solute</span></b>"));}
void UI::soluteConc_defined(){int tabIndex=scTab->currentIndex();seleSoluConcLabel[tabIndex]->setText(tr("<b><span style=\"color:black;\">Concentration</span></b>"));}
void UI::soluteConcType_defined(const QString&){int tabIndex=scTab->currentIndex();seleSoluConcTypeLabel[tabIndex]->setText(tr("<b><span style=\"color:black;\">Units</span></b>"));}

void UI::apbsFilePath_defined(){apbsLabel->setText(tr("<span style=\"color:black;\">apbs</span>:"));}
void UI::mvFilePath_defined(){mvLabel->setText(tr("<span style=\"color:black;\">multivalue</span>:"));}
void UI::hFilePath_defined(){hLabel->setText("<span style=\"color:black;\">hydropro-lnx.exe</span>:");}
void UI::msFilePath_defined(){msLabel->setText(tr("<span style=\"color:black;\">msms</span>:"));}
void UI::msaFilePath_defined(){msaLabel->setText(tr("<span style=\"color:black;\">atmtypenumbers</span>:"));}
void UI::msxFilePath_defined(){msxLabel->setText(tr("<span style=\"color:black;\">pdb_to_xyzr</span>:"));}
void UI::msxnFilePath_defined(){msxnLabel->setText(tr("<span style=\"color:black;\">pdb_to_xyzrn</span>:"));}
void UI::pFilePath_defined(){pLabel->setText(tr("<span style=\"color:black;\">main.py</span>:"));}
void UI::zFilePath_defined(){zLabel->setText(tr("<span style=\"color:black;\">ZPRED.exe</span>:"));}

void UI::autoFillFilePaths()
	{// Define File Paths using .history.txt values
	// APBS Executable
	string tmp=apbsFldr+"bin/apbs";
	apbsFilePath->setText(tr(tmp.c_str()));
	apbsLabel->setText(tr("<span style=\"color:black;\">apbs</span>:"));
	// APBS Tool: Multivalue
	tmp=apbsFldr+"share/apbs/tools/bin/multivalue";
	mvFilePath->setText(tr(tmp.c_str()));
	mvLabel->setText(tr("<span style=\"color:black;\">multivalue</span>:"));
	// HYDROPRO Linux Executable
	tmp=hydroproFldr+"hydropro10-lnx.exe";
	hFilePath->setText(tr(tmp.c_str()));
	hLabel->setText("<span style=\"color:black;\">hydropro-lnx.exe</span>:");
	// MSMS Executable
	tmp=msmsFldr+"msms.i86Linux2.2.6.1";
	msFilePath->setText(tr(tmp.c_str()));
	msLabel->setText(tr("<span style=\"color:black;\">msms</span>:"));
	// MSMS: atmTypeNumbers File
	tmp=msmsFldr+"atmtypenumbers";
	msaFilePath->setText(tr(tmp.c_str()));
	msaLabel->setText(tr("<span style=\"color:black;\">atmtypenumbers</span>:"));
	// MSMS: XYZR File
	tmp=msmsFldr+"pdb_to_xyzr";
	msxFilePath->setText(tr(tmp.c_str()));
	msxLabel->setText(tr("<span style=\"color:black;\">pdb_to_xyzr</span>:"));
	// MSMS: XYZRN File
	tmp=msmsFldr+"pdb_to_xyzrn";
	msxnFilePath->setText(tr(tmp.c_str()));
	msxnLabel->setText(tr("<span style=\"color:black;\">pdb_to_xyzrn</span>:"));
	// PDB2PQR: main.py
	tmp=pdb2pqrFldr+"main.py";
	pFilePath->setText(tr(tmp.c_str()));
	pLabel->setText(tr("<span style=\"color:black;\">main.py</span>:"));
	// ZPRED Executable File Path
	tmp=zpredFldr+"ZPRED.exe";
	zFilePath->setText(tr(tmp.c_str()));
	zLabel->setText(tr("<span style=\"color:black;\">ZPRED.exe</span>:"));
	}

void UI::label1solute(int index)
	{int tabIndex=scTab->currentIndex();
	string tmp,tmp2,c,bld;
	tmp=modeledSolutes[index];
	// Get Temperature if Specified
	tmp2=T[tabIndex]->text().toStdString();
	double Temp=298.15;
	if(tmp2.compare(NO_VALUE)!=0){Temp=strtod(tmp2.c_str(),NULL);}
	cerr<<tmp<<" Aqueous Solubility at "<<Temp<<" K: "<<getSoluteAqueousSolubility(tmp,Temp)<<" M\n";
	fSoluteLabel[tabIndex][0]->setText(tmp.c_str());	// For Output
	bld="";
	for(uint j=0;j<tmp.length();j++)
		{c=tmp[j];
		bld+=fixSubscriptNumbers(c);}
	// Get Current File Path and Find Images
	//QString theCurFldr=QDir::currentPath();
	//string Fldr=theCurFldr.toStdString()+"/";
	//string iFile=Fldr+"UI/Images/"+tmp+".png";
	//QPixmap fig1(iFile.c_str());
	//int Sz=200;
	//soluteLabel[tabIndex][0]->setFixedHeight(Sz);
	//soluteLabel[tabIndex][0]->setFixedWidth(Sz);
	//soluteLabel[tabIndex][0]->setAlignment(Qt::AlignCenter);
	//soluteLabel[tabIndex][0]->setPixmap(fig1.scaled(Sz,Sz,Qt::KeepAspectRatio));soluteLabel[tabIndex][0]->setMask(fig1.mask());
	soluteLabel[tabIndex][0]->setText(bld.c_str());
	}

void UI::label2solute(int index)
	{int tabIndex=scTab->currentIndex();
	string tmp,tmp2,c,bld;
	tmp=modeledSolutes[index];
	// Get Temperature if Specified
	tmp2=T[tabIndex]->text().toStdString();
	double Temp=298.15;
	if(tmp2.compare(NO_VALUE)!=0){Temp=strtod(tmp2.c_str(),NULL);}
	cerr<<tmp<<" Aqueous Solubility at "<<Temp<<" K: "<<getSoluteAqueousSolubility(tmp,Temp)<<" M\n";
	fSoluteLabel[tabIndex][1]->setText(tmp.c_str());	// For Output
	bld="";
	for(uint j=0;j<tmp.length();j++)
		{c=tmp[j];
		bld+=fixSubscriptNumbers(c);}
	soluteLabel[tabIndex][1]->setText(bld.c_str());}

void UI::label3solute(int index)
	{int tabIndex=scTab->currentIndex();
	string tmp,tmp2,c,bld;
	tmp=modeledSolutes[index];
	// Get Temperature if Specified
	tmp2=T[tabIndex]->text().toStdString();
	double Temp=298.15;
	if(tmp2.compare(NO_VALUE)!=0){Temp=strtod(tmp2.c_str(),NULL);}
	cerr<<tmp<<" Aqueous Solubility at "<<Temp<<" K: "<<getSoluteAqueousSolubility(tmp,Temp)<<" M\n";
	fSoluteLabel[tabIndex][2]->setText(tmp.c_str());	// For Output
	bld="";
	for(uint j=0;j<tmp.length();j++)
		{c=tmp[j];
		bld+=fixSubscriptNumbers(c);}
	soluteLabel[tabIndex][2]->setText(bld.c_str());}

void UI::label4solute(int index)
	{int tabIndex=scTab->currentIndex();
	string tmp,tmp2,c,bld;
	tmp=modeledSolutes[index];
	// Get Temperature if Specified
	tmp2=T[tabIndex]->text().toStdString();
	double Temp=298.15;
	if(tmp2.compare(NO_VALUE)!=0){Temp=strtod(tmp2.c_str(),NULL);}
	cerr<<tmp<<" Aqueous Solubility at "<<Temp<<" K: "<<getSoluteAqueousSolubility(tmp,Temp)<<" M\n";
	fSoluteLabel[tabIndex][3]->setText(tmp.c_str());	// For Output
	bld="";
	for(uint j=0;j<tmp.length();j++)
		{c=tmp[j];
		bld+=fixSubscriptNumbers(c);}
	soluteLabel[tabIndex][3]->setText(bld.c_str());}

void UI::label5solute(int index)
	{int tabIndex=scTab->currentIndex();
	string tmp,tmp2,c,bld;
	tmp=modeledSolutes[index];
	// Get Temperature if Specified
	tmp2=T[tabIndex]->text().toStdString();
	double Temp=298.15;
	if(tmp2.compare(NO_VALUE)!=0){Temp=strtod(tmp2.c_str(),NULL);}
	cerr<<tmp<<" Aqueous Solubility at "<<Temp<<" K: "<<getSoluteAqueousSolubility(tmp,Temp)<<" M\n";
	fSoluteLabel[tabIndex][4]->setText(tmp.c_str());	// For Output
	bld="";
	for(uint j=0;j<tmp.length();j++)
		{c=tmp[j];
		bld+=fixSubscriptNumbers(c);}
	soluteLabel[tabIndex][4]->setText(bld.c_str());}

void UI::label6solute(int index)
	{int tabIndex=scTab->currentIndex();
	string tmp,tmp2,c,bld;
	tmp=modeledSolutes[index];
	// Get Temperature if Specified
	tmp2=T[tabIndex]->text().toStdString();
	double Temp=298.15;
	if(tmp2.compare(NO_VALUE)!=0){Temp=strtod(tmp2.c_str(),NULL);}
	cerr<<tmp<<" Aqueous Solubility at "<<Temp<<" K: "<<getSoluteAqueousSolubility(tmp,Temp)<<" M\n";
	fSoluteLabel[tabIndex][5]->setText(tmp.c_str());	// For Output
	bld="";
	for(uint j=0;j<tmp.length();j++)
		{c=tmp[j];
		bld+=fixSubscriptNumbers(c);}
	soluteLabel[tabIndex][5]->setText(bld.c_str());}

void UI::label7solute(int index)
	{int tabIndex=scTab->currentIndex();
	string tmp,tmp2,c,bld;
	tmp=modeledSolutes[index];
	// Get Temperature if Specified
	tmp2=T[tabIndex]->text().toStdString();
	double Temp=298.15;
	if(tmp2.compare(NO_VALUE)!=0){Temp=strtod(tmp2.c_str(),NULL);}
	cerr<<tmp<<" Aqueous Solubility at "<<Temp<<" K: "<<getSoluteAqueousSolubility(tmp,Temp)<<" M\n";
	fSoluteLabel[tabIndex][6]->setText(tmp.c_str());	// For Output
	bld="";
	for(uint j=0;j<tmp.length();j++)
		{c=tmp[j];
		bld+=fixSubscriptNumbers(c);}
	soluteLabel[tabIndex][6]->setText(bld.c_str());}

void UI::label8solute(int index)
	{int tabIndex=scTab->currentIndex();
	string tmp,tmp2,c,bld;
	tmp=modeledSolutes[index];
	// Get Temperature if Specified
	tmp2=T[tabIndex]->text().toStdString();
	double Temp=298.15;
	if(tmp2.compare(NO_VALUE)!=0){Temp=strtod(tmp2.c_str(),NULL);}
	cerr<<tmp<<" Aqueous Solubility at "<<Temp<<" K: "<<getSoluteAqueousSolubility(tmp,Temp)<<" M\n";
	fSoluteLabel[tabIndex][7]->setText(tmp.c_str());	// For Output
	bld="";
	for(uint j=0;j<tmp.length();j++)
		{c=tmp[j];
		bld+=fixSubscriptNumbers(c);}
	soluteLabel[tabIndex][7]->setText(bld.c_str());}

void UI::label1soluteConc(){soluteConcLabel[scTab->currentIndex()][0]->setText(soluteConcLine[scTab->currentIndex()][0]->text());}
void UI::label2soluteConc(){soluteConcLabel[scTab->currentIndex()][1]->setText(soluteConcLine[scTab->currentIndex()][1]->text());}
void UI::label3soluteConc(){soluteConcLabel[scTab->currentIndex()][2]->setText(soluteConcLine[scTab->currentIndex()][2]->text());}
void UI::label4soluteConc(){soluteConcLabel[scTab->currentIndex()][3]->setText(soluteConcLine[scTab->currentIndex()][3]->text());}
void UI::label5soluteConc(){soluteConcLabel[scTab->currentIndex()][4]->setText(soluteConcLine[scTab->currentIndex()][4]->text());}
void UI::label6soluteConc(){soluteConcLabel[scTab->currentIndex()][5]->setText(soluteConcLine[scTab->currentIndex()][5]->text());}
void UI::label7soluteConc(){soluteConcLabel[scTab->currentIndex()][6]->setText(soluteConcLine[scTab->currentIndex()][6]->text());}
void UI::label8soluteConc(){soluteConcLabel[scTab->currentIndex()][7]->setText(soluteConcLine[scTab->currentIndex()][7]->text());}

void UI::label1soluteConcType(const QString& s){int tabIndex=scTab->currentIndex();soluteConcTypeLabel[tabIndex]->setText(s);}

void UI::label1solvent(const QString& s){int tabIndex=scTab->currentIndex();solventLabel[tabIndex][0]->setText(s);}
void UI::label2solvent(const QString& s){int tabIndex=scTab->currentIndex();solventLabel[tabIndex][1]->setText(s);}
void UI::label3solvent(const QString& s){int tabIndex=scTab->currentIndex();solventLabel[tabIndex][2]->setText(s);}
void UI::label4solvent(const QString& s){int tabIndex=scTab->currentIndex();solventLabel[tabIndex][3]->setText(s);}
void UI::label5solvent(const QString& s){int tabIndex=scTab->currentIndex();solventLabel[tabIndex][4]->setText(s);}
void UI::label6solvent(const QString& s){int tabIndex=scTab->currentIndex();solventLabel[tabIndex][5]->setText(s);}
void UI::label7solvent(const QString& s){int tabIndex=scTab->currentIndex();solventLabel[tabIndex][6]->setText(s);}
void UI::label8solvent(const QString& s){int tabIndex=scTab->currentIndex();solventLabel[tabIndex][7]->setText(s);}

void UI::label1solventConc(){solventConcLabel[scTab->currentIndex()][0]->setText(solventConcLine[scTab->currentIndex()][0]->text());}
void UI::label2solventConc(){solventConcLabel[scTab->currentIndex()][1]->setText(solventConcLine[scTab->currentIndex()][1]->text());}
void UI::label3solventConc(){solventConcLabel[scTab->currentIndex()][2]->setText(solventConcLine[scTab->currentIndex()][2]->text());}
void UI::label4solventConc(){solventConcLabel[scTab->currentIndex()][3]->setText(solventConcLine[scTab->currentIndex()][3]->text());}
void UI::label5solventConc(){solventConcLabel[scTab->currentIndex()][4]->setText(solventConcLine[scTab->currentIndex()][4]->text());}
void UI::label6solventConc(){solventConcLabel[scTab->currentIndex()][5]->setText(solventConcLine[scTab->currentIndex()][5]->text());}
void UI::label7solventConc(){solventConcLabel[scTab->currentIndex()][6]->setText(solventConcLine[scTab->currentIndex()][6]->text());}
void UI::label8solventConc(){solventConcLabel[scTab->currentIndex()][7]->setText(solventConcLine[scTab->currentIndex()][7]->text());}

void UI::label1solventConcType(const QString& s){int tabIndex=scTab->currentIndex();solventConcTypeLabel[tabIndex]->setText(s);}

void UI::addSolute()
	{int tabIndex=scTab->currentIndex();
	//cout<<"i:"<<tabIndex<<"\np:"<<scTab->tabPosition()<<endl;
	// Update Current Number of Solutes
	if(numSolute[tabIndex]<8){numSolute[tabIndex]++;}	
	// Check if pH Defined
	QString val;//=pH[tabIndex]->text();
	string tmp;//=val.toStdString();
	QLabel* TUnitsLabel=new QLabel(tr("[Kelvin]"));TUnitsLabel->setAlignment(Qt::AlignLeft);
	QLabel* pCUnitsLabel=new QLabel(tr("[g/L]"));pCUnitsLabel->setAlignment(Qt::AlignLeft);
	//
	TLabel[tabIndex]->setFixedWidth(150);
	T[tabIndex]->setFixedWidth(LINE_INPUT_SIZE);
	pH[tabIndex]->setFixedWidth(LINE_INPUT_SIZE);
	pC[tabIndex]->setFixedWidth(LINE_INPUT_SIZE);
	// Solvent(s) Selection Group
	switch(numSolvent[tabIndex])
		{case 1: solventBox=create_oneSolventBox(tabIndex);break;
		case 2: solventBox=create_twoSolventBox(tabIndex);break;
		case 3: solventBox=create_threeSolventBox(tabIndex);break;
		case 4: solventBox=create_fourSolventBox(tabIndex);break;
		case 5: solventBox=create_fiveSolventBox(tabIndex);break;
		case 6: solventBox=create_sixSolventBox(tabIndex);break;
		case 7: solventBox=create_sevenSolventBox(tabIndex);break;
		case 8: solventBox=create_eightSolventBox(tabIndex);break;
		default: break;
		}
	// Solute(s) Selection Group
	switch(numSolute[tabIndex])
		{case 1: soluteBox=create_oneSoluteBox(tabIndex);break;
		case 2: soluteBox=create_twoSoluteBox(tabIndex);break;
		case 3: soluteBox=create_threeSoluteBox(tabIndex);break;
		case 4: soluteBox=create_fourSoluteBox(tabIndex);break;
		case 5: soluteBox=create_fiveSoluteBox(tabIndex);break;
		case 6: soluteBox=create_sixSoluteBox(tabIndex);break;
		case 7: soluteBox=create_sevenSoluteBox(tabIndex);break;
		case 8: soluteBox=create_eightSoluteBox(tabIndex);break;
		default: break;
		}
	// Define Solution Conditions Tab
	QWidget *tab1 = new QWidget;

	QGridLayout *tab1hbox = new QGridLayout;
	tab1hbox->setSizeConstraint(QLayout::SetFixedSize);
	tab1hbox->addWidget(TLabel[tabIndex],0,0,1,1);tab1hbox->addWidget(T[tabIndex],0,1,1,1);tab1hbox->addWidget(TUnitsLabel,0,2,1,1);tab1hbox->addWidget(pHLabel[tabIndex],0,4,1,1);tab1hbox->addWidget(pH[tabIndex],0,5,1,1);
	tab1hbox->addWidget(solventBox,1,0,1,3);tab1hbox->addWidget(soluteBox,1,3,1,3);	
	//tab1hbox->addWidget(create_fBox(tabIndex),2,0,2,3);tab1hbox->addWidget(pCLabel[tabIndex],2,3,1,1);tab1hbox->addWidget(pC[tabIndex],2,4,1,1);tab1hbox->addWidget(pCUnitsLabel,2,5,1,1);
	if(forceParmBox[tabIndex]->isChecked())
		{tab1hbox->addWidget(forceParmBox[tabIndex],2,0,1,3);
		tab1hbox->addWidget(create_fBox(tabIndex),3,0,2,3);tab1hbox->addWidget(pCLabel[tabIndex],3,3,1,1);tab1hbox->addWidget(pC[tabIndex],3,4,1,1);tab1hbox->addWidget(pCUnitsLabel,3,5,1,1);
		tab1hbox->addWidget(potProfileBox[tabIndex],5,0,1,3);}
	else
		{tab1hbox->addWidget(forceParmBox[tabIndex],2,0,1,3);tab1hbox->addWidget(pCLabel[tabIndex],2,3,1,1);tab1hbox->addWidget(pC[tabIndex],2,4,1,1);tab1hbox->addWidget(pCUnitsLabel,2,5,1,1);
		tab1hbox->addWidget(potProfileBox[tabIndex],3,0,1,3);}
	tab1->setLayout(tab1hbox);
	// Rewrite tab
	string ind=cnvrtNumToStrng(tabIndex,0);
	tmp="&"+ind;
	// Clear Previous
	scTab->removeTab(tabIndex);
	scTab->insertTab(tabIndex,tab1,tmp.c_str());
	scTab->setCurrentIndex(tabIndex);
	}

void UI::addSolvent()
	{int tabIndex=scTab->currentIndex();	
	// Update Current Number of Solvents
	if(numSolvent[tabIndex]<8){numSolvent[tabIndex]++;}
	// Check if pH Defined
	QString val;//=pH[tabIndex]->text();
	string tmp;//=val.toStdString();
	QLabel* TUnitsLabel=new QLabel(tr("[Kelvin]"));TUnitsLabel->setAlignment(Qt::AlignLeft);
	QLabel* pCUnitsLabel=new QLabel(tr("[g/L]"));pCUnitsLabel->setAlignment(Qt::AlignLeft);
	//
	TLabel[tabIndex]->setFixedWidth(150);
	T[tabIndex]->setFixedWidth(LINE_INPUT_SIZE);
	pH[tabIndex]->setFixedWidth(LINE_INPUT_SIZE);
	pC[tabIndex]->setFixedWidth(LINE_INPUT_SIZE);
	// Solvent(s) Selection Group
	switch(numSolvent[tabIndex])
		{case 1: solventBox=create_oneSolventBox(tabIndex);break;
		case 2: solventBox=create_twoSolventBox(tabIndex);break;
		case 3: solventBox=create_threeSolventBox(tabIndex);break;
		case 4: solventBox=create_fourSolventBox(tabIndex);break;
		case 5: solventBox=create_fiveSolventBox(tabIndex);break;
		case 6: solventBox=create_sixSolventBox(tabIndex);break;
		case 7: solventBox=create_sevenSolventBox(tabIndex);break;
		case 8: solventBox=create_eightSolventBox(tabIndex);break;
		default: break;
		}

	// Solute(s) Selection Group
	switch(numSolute[tabIndex])
		{case 1: soluteBox=create_oneSoluteBox(tabIndex);break;
		case 2: soluteBox=create_twoSoluteBox(tabIndex);break;
		case 3: soluteBox=create_threeSoluteBox(tabIndex);break;
		case 4: soluteBox=create_fourSoluteBox(tabIndex);break;
		case 5: soluteBox=create_fiveSoluteBox(tabIndex);break;
		case 6: soluteBox=create_sixSoluteBox(tabIndex);break;
		case 7: soluteBox=create_sevenSoluteBox(tabIndex);break;
		case 8: soluteBox=create_eightSoluteBox(tabIndex);break;
		default: break;
		}

	// Define Solution Conditions Tab
	QWidget *tab1 = new QWidget;
	
	QGridLayout *tab1hbox = new QGridLayout;
	tab1hbox->setSizeConstraint(QLayout::SetFixedSize);
	tab1hbox->addWidget(TLabel[tabIndex],0,0,1,1);tab1hbox->addWidget(T[tabIndex],0,1,1,1);tab1hbox->addWidget(TUnitsLabel,0,2,1,1);tab1hbox->addWidget(pHLabel[tabIndex],0,4,1,1);tab1hbox->addWidget(pH[tabIndex],0,5,1,1);
	tab1hbox->addWidget(solventBox,1,0,1,3);tab1hbox->addWidget(soluteBox,1,3,1,3);	
	//tab1hbox->addWidget(create_fBox(tabIndex),2,0,2,3);tab1hbox->addWidget(pCLabel[tabIndex],2,3,1,1);tab1hbox->addWidget(pC[tabIndex],2,4,1,1);tab1hbox->addWidget(pCUnitsLabel,2,5,1,1);
	if(forceParmBox[tabIndex]->isChecked())
		{tab1hbox->addWidget(forceParmBox[tabIndex],2,0,1,3);
		tab1hbox->addWidget(create_fBox(tabIndex),3,0,2,3);tab1hbox->addWidget(pCLabel[tabIndex],3,3,1,1);tab1hbox->addWidget(pC[tabIndex],3,4,1,1);tab1hbox->addWidget(pCUnitsLabel,3,5,1,1);
		tab1hbox->addWidget(potProfileBox[tabIndex],5,0,1,3);}
	else
		{tab1hbox->addWidget(forceParmBox[tabIndex],2,0,1,3);tab1hbox->addWidget(pCLabel[tabIndex],2,3,1,1);tab1hbox->addWidget(pC[tabIndex],2,4,1,1);tab1hbox->addWidget(pCUnitsLabel,2,5,1,1);
		tab1hbox->addWidget(potProfileBox[tabIndex],3,0,1,3);}
	tab1->setLayout(tab1hbox);
	// Rewrite tab
	tmp="&"+cnvrtNumToStrng(tabIndex,0);
	// Clear Previous
	scTab->removeTab(tabIndex);
	scTab->insertTab(tabIndex,tab1,tmp.c_str());
	scTab->setCurrentIndex(tabIndex);
	}

void UI::addSolCondTab()
	{//int tabIndex=scTab->currentIndex();
	// Copy Previous Data
	string solvLst="",soluLst="",delimiter=";";
	for(int i=0;i<numSolCond;i++)
		{solvLst+=cnvrtNumToStrng(numSolvent[i],0)+delimiter;
		soluLst+=cnvrtNumToStrng(numSolute[i],0)+delimiter;}
	// Increase NUmber of Solution Conditions
	numSolCond++;
	if(numSolCond>1)
		{// Show Removal Button
		remSolCondButton->show();}
	solvLst+="1"+delimiter;
	soluLst+="1"+delimiter;
	// Update Parameters
	delete [] numSolvent;delete [] numSolute;
	numSolvent=fill_int_array(solvLst,numSolCond,delimiter);
	numSolute=fill_int_array(soluLst,numSolCond,delimiter);
	// Define New Tab Index
	int N=scTab->count();
	//QLabel* pHLabel=new QLabel("<b><span style=\"color:red;\">pH</span></b>:");pHLabel->setBuddy(pH[N]);pHLabel->setAlignment(Qt::AlignRight);
	pHLabel[N]->setText(tr("<b><span style=\"color:red;\">pH</span></b>:"));
	//QLabel* TLabel=new QLabel("<b><span style=\"color:red;\">Temperature</span></b>:");TLabel->setBuddy(T[N]);TLabel->setAlignment(Qt::AlignRight);
	TLabel[N]->setText(tr("<b><span style=\"color:red;\">Temperature</span></b>:"));
	QLabel* TUnitsLabel=new QLabel(tr("[Kelvin]"));TUnitsLabel->setAlignment(Qt::AlignLeft);
	//QLabel* pCLabel=new QLabel("<b><span style=\"color:red;\">Protein Concentration</span></b>:");pCLabel->setBuddy(pC[N]);pCLabel->setAlignment(Qt::AlignRight);
	pCLabel[N]=new QLabel(tr("<b><span style=\"color:red;\">Protein Concentration</span></b>:"));pCLabel[N]->setAlignment(Qt::AlignRight);
	QLabel* pCUnitsLabel=new QLabel(tr("[g/L]"));pCUnitsLabel->setAlignment(Qt::AlignLeft);
	//
	TLabel[N]->setFixedWidth(150);
	T[N]->setFixedWidth(LINE_INPUT_SIZE);
	pH[N]->setFixedWidth(LINE_INPUT_SIZE);
	pC[N]->setFixedWidth(LINE_INPUT_SIZE);
	// Define Solution Conditions Tab
	QWidget *tab1 = new QWidget;
	QGridLayout *tab1hbox = new QGridLayout;
	tab1hbox->setSizeConstraint(QLayout::SetFixedSize);
	tab1hbox->addWidget(TLabel[N],0,0,1,1);tab1hbox->addWidget(T[N],0,1,1,1);tab1hbox->addWidget(TUnitsLabel,0,2,1,1);tab1hbox->addWidget(pHLabel[N],0,4,1,1);tab1hbox->addWidget(pH[N],0,5,1,1);
	tab1hbox->addWidget(create_oneSolventBox(N),1,0,1,3);tab1hbox->addWidget(create_oneSoluteBox(N),1,3,1,3);
	//tab1hbox->addWidget(create_fBox(N),2,0,2,3);tab1hbox->addWidget(pCLabel[N],2,3,1,1);tab1hbox->addWidget(pC[N],2,4,1,1);tab1hbox->addWidget(pCUnitsLabel,2,5,1,1);
	if(forceParmBox[N]->isChecked())
		{tab1hbox->addWidget(forceParmBox[N],2,0,1,3);
		tab1hbox->addWidget(create_fBox(N),3,0,2,3);tab1hbox->addWidget(pCLabel[N],3,3,1,1);tab1hbox->addWidget(pC[N],3,4,1,1);tab1hbox->addWidget(pCUnitsLabel,3,5,1,1);
		tab1hbox->addWidget(potProfileBox[N],5,0,1,3);}
	else
		{tab1hbox->addWidget(forceParmBox[N],2,0,1,3);tab1hbox->addWidget(pCLabel[N],2,3,1,1);tab1hbox->addWidget(pC[N],2,4,1,1);tab1hbox->addWidget(pCUnitsLabel,2,5,1,1);
		tab1hbox->addWidget(potProfileBox[N],3,0,1,3);}
	tab1->setLayout(tab1hbox);

	string tmp="&"+cnvrtNumToStrng(N,0);
	scTab->addTab(tab1, tr(tmp.c_str()));
	}

QGroupBox* UI::create_fBox(int tabIndex)
	{// Force Solution Conditions Group
	QGroupBox* Output = new QGroupBox;//(tr("Specify Solution Parameters Instead of Computation"));
	//Output->setCheckable(true);
	//Output->setChecked(false);
	//
	QLabel* fViscLabel=new QLabel(tr("User-Specified Viscosity (η):"));fViscLabel->setAlignment(Qt::AlignCenter|Qt::AlignRight);
	QLabel* fViscUnitsLabel=new QLabel(tr("[Ns/m²]"));fViscUnitsLabel->setAlignment(Qt::AlignLeft);
	//
	QLabel* fDensLabel=new QLabel(tr("User-Specified Density (ρ):"));fDensLabel->setAlignment(Qt::AlignCenter|Qt::AlignRight);
	QLabel* fDensUnitsLabel=new QLabel(tr("[kg/L]"));fDensUnitsLabel->setAlignment(Qt::AlignLeft);
	//
	QLabel* fDielLabel=new QLabel(tr("User-Specified Relative Dielectric (ε<sub>r</sub>):"));fDielLabel->setAlignment(Qt::AlignCenter|Qt::AlignRight);
	//
	QLabel* fXspLabel=new QLabel(tr("User-Specified Hydration Thickness (X<sub>SP</sub>):"));fXspLabel->setAlignment(Qt::AlignCenter|Qt::AlignRight);
	QLabel* fXspUnitsLabel=new QLabel(tr("[Å]"));fXspUnitsLabel->setAlignment(Qt::AlignLeft);
	//
	fDiel[tabIndex]->setFixedWidth(150);
	fDens[tabIndex]->setFixedWidth(150);
	fVisc[tabIndex]->setFixedWidth(150);
	fXsp[tabIndex]->setFixedWidth(150);
	//
	QGridLayout *layout = new QGridLayout;
	layout->addWidget(fDielLabel,0,0,1,1);layout->addWidget(fDiel[tabIndex],0,1,1,1);
	layout->addWidget(fDensLabel,1,0,1,1);layout->addWidget(fDens[tabIndex],1,1,1,1);layout->addWidget(fDensUnitsLabel,1,2,1,1);
	layout->addWidget(fViscLabel,2,0,1,1);layout->addWidget(fVisc[tabIndex],2,1,1,1);layout->addWidget(fViscUnitsLabel,2,2,1,1);
	layout->addWidget(fXspLabel,3,0,1,1);layout->addWidget(fXsp[tabIndex],3,1,1,1);layout->addWidget(fXspUnitsLabel,3,2,1,1);
	Output->setLayout(layout);
	return Output;}

QGroupBox* UI::create_oneSoluteBox(int tabIndex)
	{QGroupBox* Output=new QGroupBox(tr("Solute(s) Selection"));
	QPushButton* addSoluteButton=new QPushButton(tr("Add Solute"));
	addSoluteButton->setFixedWidth(150);
	connect(addSoluteButton,SIGNAL(clicked()),this,SLOT(addSolute()));

	string tmp,c,bld;
	QComboBox *selectSoluteButton_1 = new QComboBox;
	for(int i=0;i<numModeledSolutes;i++)
		{tmp=modeledSolutes[i];
		bld="";
		for(uint j=0;j<tmp.length();j++)
			{c=tmp[j];
			bld+=fixSubscriptNumbers(c);}
		selectSoluteButton_1->addItem(tr(bld.c_str()));}
	selectSoluteButton_1->setFixedWidth(150);	
	connect(selectSoluteButton_1, SIGNAL(activated(int)),this, SLOT(label1solute(int)));
	connect(selectSoluteButton_1,SIGNAL(activated(const QString&)),this,SLOT(solute_defined(const QString&)));

	QComboBox *selectSoluteConcTypeButton_1=new QComboBox;
	selectSoluteConcTypeButton_1->addItem(tr("Molar"));
	selectSoluteConcTypeButton_1->setFixedWidth(100);
	connect(selectSoluteConcTypeButton_1,SIGNAL(activated(const QString&)),this,SLOT(label1soluteConcType(const QString&)));
	connect(selectSoluteConcTypeButton_1,SIGNAL(activated(const QString&)),this,SLOT(soluteConcType_defined(const QString&)));

	soluteConcLine[tabIndex][0]->setFixedWidth(CONC_SIZE);
	soluteConcLabel[tabIndex][0]->setFixedWidth(CONC_SIZE);
	connect(soluteConcLine[tabIndex][0],SIGNAL(returnPressed()),this,SLOT(soluteConc_defined()));
	connect(soluteConcLine[tabIndex][0],SIGNAL(returnPressed()),this,SLOT(label1soluteConc()));

	QGridLayout *layout = new QGridLayout;
	layout->setSizeConstraint(QLayout::SetFixedSize);
	layout->addWidget(addSoluteButton,0,0,1,1);
	layout->addWidget(seleSoluLabel[tabIndex],1,0,1,1);layout->addWidget(seleSoluConcLabel[tabIndex],1,1,1,1);layout->addWidget(seleSoluConcTypeLabel[tabIndex],1,2,1,1);
	layout->addWidget(selectSoluteButton_1,2,0,1,1);layout->addWidget(soluteConcLine[tabIndex][0],2,1,1,1);layout->addWidget(selectSoluteConcTypeButton_1,2,2,1,1);
	layout->addWidget(soluteLabel[tabIndex][0],3,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][0],3,1,1,1);layout->addWidget(soluteConcTypeLabel[tabIndex],3,2,1,1);
	Output->setLayout(layout);
	return Output;}

QGroupBox* UI::create_oneSolventBox(int tabIndex)
	{// Solvent(s) Selection Group
	QGroupBox* Output=new QGroupBox(tr("Solvent(s) Selection"));
	QPushButton* addSolventButton=new QPushButton(tr("Add Solvent"));
	addSolventButton->setFixedWidth(150);
	connect(addSolventButton,SIGNAL(clicked()),this,SLOT(addSolvent()));

	string tmp;
	QComboBox *selectSolventButton_1=new QComboBox;
	for(int i=0;i<numModeledSolvents;i++){tmp=modeledSolvents[i];selectSolventButton_1->addItem(tr(tmp.c_str()));}
	selectSolventButton_1->setFixedWidth(180);
	selectSolventButton_1->setCurrentIndex(numModeledSolvents-1);
	connect(selectSolventButton_1,SIGNAL(activated(const QString&)),this,SLOT(label1solvent(const QString&)));
	connect(selectSolventButton_1,SIGNAL(activated(const QString&)),this,SLOT(solvent_defined(const QString&)));

	QComboBox *selectSolventConcTypeButton_1=new QComboBox;
	for(int i=0;i<numModeledSolventConcTypes;i++){tmp=modeledSolventConcTypes[i];selectSolventConcTypeButton_1->addItem(tr(tmp.c_str()));}
	selectSolventConcTypeButton_1->setFixedWidth(150);
	connect(selectSolventConcTypeButton_1,SIGNAL(activated(const QString&)),this,SLOT(label1solventConcType(const QString&)));
	connect(selectSolventConcTypeButton_1,SIGNAL(activated(const QString&)),this,SLOT(solventConcType_defined(const QString&)));

	solventConcLine[tabIndex][0]->setFixedWidth(CONC_SIZE);
	solventConcLabel[tabIndex][0]->setFixedWidth(CONC_SIZE);
	connect(solventConcLine[tabIndex][0],SIGNAL(returnPressed()),this,SLOT(solventConc_defined()));
	connect(solventConcLine[tabIndex][0],SIGNAL(returnPressed()),this,SLOT(label1solventConc()));

	QGridLayout *layout = new QGridLayout;
	layout->setSizeConstraint(QLayout::SetFixedSize);
	layout->addWidget(addSolventButton,0,0,1,1);
	layout->addWidget(seleSolvLabel[tabIndex],1,0,1,1);layout->addWidget(seleSolvConcLabel[tabIndex],1,1,1,1);layout->addWidget(seleSolvConcTypeLabel[tabIndex],1,2,1,1);
	layout->addWidget(selectSolventButton_1,2,0,1,1);layout->addWidget(solventConcLine[tabIndex][0],2,1,1,1);layout->addWidget(selectSolventConcTypeButton_1,2,2,1,1);
	layout->addWidget(solventLabel[tabIndex][0],3,0,1,1);layout->addWidget(solventConcLabel[tabIndex][0],3,1,1,1);layout->addWidget(solventConcTypeLabel[tabIndex],3,2,1,1);
	Output->setLayout(layout);
	return Output;}

QGroupBox* UI::create_twoSoluteBox(int tabIndex)
	{QGroupBox* Output=new QGroupBox(tr("Solute(s) Selection"));
	QPushButton* addSoluteButton=new QPushButton(tr("Add Solute"));
	addSoluteButton->setFixedWidth(150);
	connect(addSoluteButton,SIGNAL(clicked()),this,SLOT(addSolute()));
	QPushButton* remSoluteButton=new QPushButton(tr("Remove Solute"));
	remSoluteButton->setFixedWidth(150);
	connect(remSoluteButton,SIGNAL(clicked()),this,SLOT(remSolute()));
	//
	string tmp,c,bld;
	QComboBox *selectSoluteButton_1 = new QComboBox;
	QComboBox *selectSoluteButton_2 = new QComboBox;
	for(int i=0;i<numModeledSolutes;i++)
		{tmp=modeledSolutes[i];
		bld="";
		for(uint j=0;j<tmp.length();j++)
			{c=tmp[j];
			bld+=fixSubscriptNumbers(c);}
		selectSoluteButton_1->addItem(tr(bld.c_str()));selectSoluteButton_2->addItem(tr(bld.c_str()));}
	selectSoluteButton_1->setFixedWidth(150);
	connect(selectSoluteButton_1, SIGNAL(activated(int)),this, SLOT(label1solute(int)));
	connect(selectSoluteButton_2, SIGNAL(activated(int)),this, SLOT(label2solute(int)));
	connect(selectSoluteButton_1,SIGNAL(activated(const QString&)),this,SLOT(solute_defined(const QString&)));

	QComboBox *selectSoluteConcTypeButton_1=new QComboBox;
	selectSoluteConcTypeButton_1->addItem(tr("Molar"));
	selectSoluteConcTypeButton_1->setFixedWidth(100);
	connect(selectSoluteConcTypeButton_1,SIGNAL(activated(const QString&)),this,SLOT(label1soluteConcType(const QString&)));
	connect(selectSoluteConcTypeButton_1,SIGNAL(activated(const QString&)),this,SLOT(soluteConcType_defined(const QString&)));
	//
	for(int i=0;i<2;i++){soluteConcLine[tabIndex][i]->setFixedWidth(CONC_SIZE);soluteConcLabel[tabIndex][i]->setFixedWidth(CONC_SIZE);}
	//
	connect(soluteConcLine[tabIndex][0],SIGNAL(returnPressed()),this,SLOT(soluteConc_defined()));
	connect(soluteConcLine[tabIndex][0],SIGNAL(returnPressed()),this,SLOT(label1soluteConc()));
	connect(soluteConcLine[tabIndex][1],SIGNAL(returnPressed()),this,SLOT(label2soluteConc()));

	QGridLayout *layout = new QGridLayout;
	layout->setSizeConstraint(QLayout::SetFixedSize);
	layout->addWidget(addSoluteButton,0,0,1,1);layout->addWidget(remSoluteButton,0,1,1,2);
	layout->addWidget(seleSoluLabel[tabIndex],1,0,1,1);layout->addWidget(seleSoluConcLabel[tabIndex],1,1,1,1);layout->addWidget(seleSoluConcTypeLabel[tabIndex],1,2,1,1);
	layout->addWidget(selectSoluteButton_1,2,0,1,1);layout->addWidget(soluteConcLine[tabIndex][0],2,1,1,1);layout->addWidget(selectSoluteConcTypeButton_1,2,2,1,1);
	layout->addWidget(soluteLabel[tabIndex][0],3,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][0],3,1,1,1);layout->addWidget(soluteConcTypeLabel[tabIndex],3,2,1,1);
	layout->addWidget(selectSoluteButton_2,4,0,1,1);layout->addWidget(soluteConcLine[tabIndex][1],4,1,1,1);
	layout->addWidget(soluteLabel[tabIndex][1],5,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][1],5,1,1,1);
	Output->setLayout(layout);
	return Output;}

QGroupBox* UI::create_twoSolventBox(int tabIndex)
	{// Solvent(s) Selection Group
	QGroupBox* Output=new QGroupBox(tr("Solvent(s) Selection"));
	QPushButton* addSolventButton=new QPushButton(tr("Add Solvent"));
	addSolventButton->setFixedWidth(150);
	connect(addSolventButton,SIGNAL(clicked()),this,SLOT(addSolvent()));
	QPushButton* remSolventButton=new QPushButton(tr("Remove Solvent"));
	remSolventButton->setFixedWidth(150);
	connect(remSolventButton,SIGNAL(clicked()),this,SLOT(remSolvent()));
	//
	string tmp;
	QComboBox *selectSolventButton_1=new QComboBox;
	QComboBox *selectSolventButton_2=new QComboBox;
	for(int i=0;i<numModeledSolvents;i++)
		{tmp=modeledSolvents[i];selectSolventButton_1->addItem(tr(tmp.c_str()));selectSolventButton_2->addItem(tr(tmp.c_str()));}
	selectSolventButton_1->setFixedWidth(180);
	selectSolventButton_1->setCurrentIndex(numModeledSolvents-1);
	connect(selectSolventButton_1,SIGNAL(activated(const QString&)),this,SLOT(label1solvent(const QString&)));
	connect(selectSolventButton_2,SIGNAL(activated(const QString&)),this,SLOT(label2solvent(const QString&)));
	connect(selectSolventButton_1,SIGNAL(activated(const QString&)),this,SLOT(solvent_defined(const QString&)));

	QComboBox *selectSolventConcTypeButton_1=new QComboBox;
	for(int i=0;i<numModeledSolventConcTypes;i++){tmp=modeledSolventConcTypes[i];selectSolventConcTypeButton_1->addItem(tr(tmp.c_str()));}
	selectSolventConcTypeButton_1->setFixedWidth(150);
	connect(selectSolventConcTypeButton_1,SIGNAL(activated(const QString&)),this,SLOT(label1solventConcType(const QString&)));
	connect(selectSolventConcTypeButton_1,SIGNAL(activated(const QString&)),this,SLOT(solventConcType_defined(const QString&)));
	//
	for(int i=0;i<2;i++){solventConcLine[tabIndex][i]->setFixedWidth(CONC_SIZE);solventConcLabel[tabIndex][i]->setFixedWidth(CONC_SIZE);}
	//
	connect(solventConcLine[tabIndex][0],SIGNAL(returnPressed()),this,SLOT(solventConc_defined()));
	connect(solventConcLine[tabIndex][0],SIGNAL(returnPressed()),this,SLOT(label1solventConc()));
	connect(solventConcLine[tabIndex][1],SIGNAL(returnPressed()),this,SLOT(label2solventConc()));

	QGridLayout *layout = new QGridLayout;
	layout->setSizeConstraint(QLayout::SetFixedSize);
	layout->addWidget(addSolventButton,0,0,1,1);layout->addWidget(remSolventButton,0,1,1,2);
	layout->addWidget(seleSolvLabel[tabIndex],1,0,1,1);layout->addWidget(seleSolvConcLabel[tabIndex],1,1,1,1);layout->addWidget(seleSolvConcTypeLabel[tabIndex],1,2,1,1);
	layout->addWidget(selectSolventButton_1,2,0,1,1);layout->addWidget(solventConcLine[tabIndex][0],2,1,1,1);layout->addWidget(selectSolventConcTypeButton_1,2,2,1,1);
	layout->addWidget(solventLabel[tabIndex][0],3,0,1,1);layout->addWidget(solventConcLabel[tabIndex][0],3,1,1,1);layout->addWidget(solventConcTypeLabel[tabIndex],3,2,1,1);
	layout->addWidget(selectSolventButton_2,4,0,1,1);layout->addWidget(solventConcLine[tabIndex][1],4,1,1,1);
	layout->addWidget(solventLabel[tabIndex][1],5,0,1,1);layout->addWidget(solventConcLabel[tabIndex][1],5,1,1,1);
	Output->setLayout(layout);
	return Output;}

QGroupBox* UI::create_threeSoluteBox(int tabIndex)
	{QGroupBox* Output=new QGroupBox(tr("Solute(s) Selection"));
	QPushButton* addSoluteButton=new QPushButton(tr("Add Solute"));
	addSoluteButton->setFixedWidth(150);
	connect(addSoluteButton,SIGNAL(clicked()),this,SLOT(addSolute()));
	QPushButton* remSoluteButton=new QPushButton(tr("Remove Solute"));
	remSoluteButton->setFixedWidth(150);
	connect(remSoluteButton,SIGNAL(clicked()),this,SLOT(remSolute()));
	//
	string tmp,c,bld;
	QComboBox *selectSoluteButton_1 = new QComboBox;
	QComboBox *selectSoluteButton_2 = new QComboBox;
	QComboBox *selectSoluteButton_3 = new QComboBox;
	for(int i=0;i<numModeledSolutes;i++)
		{tmp=modeledSolutes[i];
		bld="";
		for(uint j=0;j<tmp.length();j++)
			{c=tmp[j];
			bld+=fixSubscriptNumbers(c);}
		selectSoluteButton_1->addItem(tr(bld.c_str()));selectSoluteButton_2->addItem(tr(bld.c_str()));selectSoluteButton_3->addItem(tr(bld.c_str()));}
	selectSoluteButton_1->setFixedWidth(150);
	connect(selectSoluteButton_1, SIGNAL(activated(int)),this, SLOT(label1solute(int)));
	connect(selectSoluteButton_2, SIGNAL(activated(int)),this, SLOT(label2solute(int)));
	connect(selectSoluteButton_3, SIGNAL(activated(int)),this, SLOT(label3solute(int)));
	connect(selectSoluteButton_1,SIGNAL(activated(const QString&)),this,SLOT(solute_defined(const QString&)));

	QComboBox *selectSoluteConcTypeButton_1=new QComboBox;
	selectSoluteConcTypeButton_1->addItem(tr("Molar"));
	selectSoluteConcTypeButton_1->setFixedWidth(100);
	connect(selectSoluteConcTypeButton_1,SIGNAL(activated(const QString&)),this,SLOT(label1soluteConcType(const QString&)));
	connect(selectSoluteConcTypeButton_1,SIGNAL(activated(const QString&)),this,SLOT(soluteConcType_defined(const QString&)));
	//
	for(int i=0;i<3;i++){soluteConcLine[tabIndex][i]->setFixedWidth(CONC_SIZE);soluteConcLabel[tabIndex][i]->setFixedWidth(CONC_SIZE);}
	connect(soluteConcLine[tabIndex][0],SIGNAL(returnPressed()),this,SLOT(soluteConc_defined()));
	connect(soluteConcLine[tabIndex][0],SIGNAL(returnPressed()),this,SLOT(label1soluteConc()));
	connect(soluteConcLine[tabIndex][1],SIGNAL(returnPressed()),this,SLOT(label2soluteConc()));
	connect(soluteConcLine[tabIndex][2],SIGNAL(returnPressed()),this,SLOT(label3soluteConc()));

	QGridLayout *layout = new QGridLayout;
	layout->setSizeConstraint(QLayout::SetFixedSize);
	layout->addWidget(addSoluteButton,0,0,1,1);layout->addWidget(remSoluteButton,0,1,1,2);
	layout->addWidget(seleSoluLabel[tabIndex],1,0,1,1);layout->addWidget(seleSoluConcLabel[tabIndex],1,1,1,1);layout->addWidget(seleSoluConcTypeLabel[tabIndex],1,2,1,1);
	layout->addWidget(selectSoluteButton_1,2,0,1,1);layout->addWidget(soluteConcLine[tabIndex][0],2,1,1,1);layout->addWidget(selectSoluteConcTypeButton_1,2,2,1,1);
	layout->addWidget(soluteLabel[tabIndex][0],3,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][0],3,1,1,1);layout->addWidget(soluteConcTypeLabel[tabIndex],3,2,1,1);
	layout->addWidget(selectSoluteButton_2,4,0,1,1);layout->addWidget(soluteConcLine[tabIndex][1],4,1,1,1);
	layout->addWidget(soluteLabel[tabIndex][1],5,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][1],5,1,1,1);
	layout->addWidget(selectSoluteButton_3,6,0,1,1);layout->addWidget(soluteConcLine[tabIndex][2],6,1,1,1);
	layout->addWidget(soluteLabel[tabIndex][2],7,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][2],7,1,1,1);
	Output->setLayout(layout);
	return Output;}

QGroupBox* UI::create_threeSolventBox(int tabIndex)
	{// Solvent(s) Selection Group
	QGroupBox* Output=new QGroupBox(tr("Solvent(s) Selection"));
	QPushButton* addSolventButton=new QPushButton(tr("Add Solvent"));
	addSolventButton->setFixedWidth(150);
	connect(addSolventButton,SIGNAL(clicked()),this,SLOT(addSolvent()));
	QPushButton* remSolventButton=new QPushButton(tr("Remove Solvent"));
	remSolventButton->setFixedWidth(150);
	connect(remSolventButton,SIGNAL(clicked()),this,SLOT(remSolvent()));
	//QPushButton *selectSolventButton_1 = new QPushButton(tr("Select a Solvent..."));
	string tmp;
	QComboBox *selectSolventButton_1=new QComboBox;
	QComboBox *selectSolventButton_2=new QComboBox;
	QComboBox *selectSolventButton_3=new QComboBox;
	for(int i=0;i<numModeledSolvents;i++)
		{tmp=modeledSolvents[i];selectSolventButton_1->addItem(tr(tmp.c_str()));selectSolventButton_2->addItem(tr(tmp.c_str()));selectSolventButton_3->addItem(tr(tmp.c_str()));}
	selectSolventButton_1->setFixedWidth(180);
	selectSolventButton_1->setCurrentIndex(numModeledSolvents-1);
	connect(selectSolventButton_1,SIGNAL(activated(const QString&)),this,SLOT(label1solvent(const QString&)));
	connect(selectSolventButton_2,SIGNAL(activated(const QString&)),this,SLOT(label2solvent(const QString&)));
	connect(selectSolventButton_3,SIGNAL(activated(const QString&)),this,SLOT(label3solvent(const QString&)));
	connect(selectSolventButton_1,SIGNAL(activated(const QString&)),this,SLOT(solvent_defined(const QString&)));

	QComboBox *selectSolventConcTypeButton_1=new QComboBox;
	for(int i=0;i<numModeledSolventConcTypes;i++){tmp=modeledSolventConcTypes[i];selectSolventConcTypeButton_1->addItem(tr(tmp.c_str()));}
	selectSolventConcTypeButton_1->setFixedWidth(150);
	connect(selectSolventConcTypeButton_1,SIGNAL(activated(const QString&)),this,SLOT(label1solventConcType(const QString&)));
	connect(selectSolventConcTypeButton_1,SIGNAL(activated(const QString&)),this,SLOT(solventConcType_defined(const QString&)));
	//
	for(int i=0;i<3;i++){solventConcLine[tabIndex][i]->setFixedWidth(CONC_SIZE);solventConcLabel[tabIndex][i]->setFixedWidth(CONC_SIZE);}
	connect(solventConcLine[tabIndex][0],SIGNAL(returnPressed()),this,SLOT(solventConc_defined()));
	connect(solventConcLine[tabIndex][0],SIGNAL(returnPressed()),this,SLOT(label1solventConc()));
	connect(solventConcLine[tabIndex][1],SIGNAL(returnPressed()),this,SLOT(label2solventConc()));
	connect(solventConcLine[tabIndex][2],SIGNAL(returnPressed()),this,SLOT(label3solventConc()));

	QGridLayout *layout = new QGridLayout;
	layout->setSizeConstraint(QLayout::SetFixedSize);
	layout->addWidget(addSolventButton,0,0,1,1);layout->addWidget(remSolventButton,0,1,1,2);
	layout->addWidget(seleSolvLabel[tabIndex],1,0,1,1);layout->addWidget(seleSolvConcLabel[tabIndex],1,1,1,1);layout->addWidget(seleSolvConcTypeLabel[tabIndex],1,2,1,1);
	layout->addWidget(selectSolventButton_1,2,0,1,1);layout->addWidget(solventConcLine[tabIndex][0],2,1,1,1);layout->addWidget(selectSolventConcTypeButton_1,2,2,1,1);
	layout->addWidget(solventLabel[tabIndex][0],3,0,1,1);layout->addWidget(solventConcLabel[tabIndex][0],3,1,1,1);layout->addWidget(solventConcTypeLabel[tabIndex],3,2,1,1);
	layout->addWidget(selectSolventButton_2,4,0,1,1);layout->addWidget(solventConcLine[tabIndex][1],4,1,1,1);
	layout->addWidget(solventLabel[tabIndex][1],5,0,1,1);layout->addWidget(solventConcLabel[tabIndex][1],5,1,1,1);
	layout->addWidget(selectSolventButton_3,6,0,1,1);layout->addWidget(solventConcLine[tabIndex][2],6,1,1,1);
	layout->addWidget(solventLabel[tabIndex][2],7,0,1,1);layout->addWidget(solventConcLabel[tabIndex][2],7,1,1,1);
	Output->setLayout(layout);
	return Output;}

QGroupBox* UI::create_fourSoluteBox(int tabIndex)
	{QGroupBox* Output=new QGroupBox(tr("Solute(s) Selection"));
	QPushButton* addSoluteButton=new QPushButton(tr("Add Solute"));
	addSoluteButton->setFixedWidth(150);
	connect(addSoluteButton,SIGNAL(clicked()),this,SLOT(addSolute()));
	QPushButton* remSoluteButton=new QPushButton(tr("Remove Solute"));
	remSoluteButton->setFixedWidth(150);
	connect(remSoluteButton,SIGNAL(clicked()),this,SLOT(remSolute()));
	//
	string tmp,c,bld;
	QComboBox *selectSoluteButton_1 = new QComboBox;
	QComboBox *selectSoluteButton_2 = new QComboBox;
	QComboBox *selectSoluteButton_3 = new QComboBox;
	QComboBox *selectSoluteButton_4 = new QComboBox;
	for(int i=0;i<numModeledSolutes;i++)
		{tmp=modeledSolutes[i];
		bld="";
		for(uint j=0;j<tmp.length();j++)
			{c=tmp[j];
			bld+=fixSubscriptNumbers(c);}
		selectSoluteButton_1->addItem(tr(bld.c_str()));selectSoluteButton_2->addItem(tr(bld.c_str()));
		selectSoluteButton_3->addItem(tr(bld.c_str()));selectSoluteButton_4->addItem(tr(bld.c_str()));}
	selectSoluteButton_1->setFixedWidth(150);
	connect(selectSoluteButton_1, SIGNAL(activated(int)),this, SLOT(label1solute(int)));
	connect(selectSoluteButton_2, SIGNAL(activated(int)),this, SLOT(label2solute(int)));
	connect(selectSoluteButton_3, SIGNAL(activated(int)),this, SLOT(label3solute(int)));
	connect(selectSoluteButton_4, SIGNAL(activated(int)),this, SLOT(label4solute(int)));
	connect(selectSoluteButton_1,SIGNAL(activated(const QString&)),this,SLOT(solute_defined(const QString&)));

	QComboBox *selectSoluteConcTypeButton_1=new QComboBox;
	selectSoluteConcTypeButton_1->addItem(tr("Molar"));
	selectSoluteConcTypeButton_1->setFixedWidth(100);
	connect(selectSoluteConcTypeButton_1,SIGNAL(activated(const QString&)),this,SLOT(label1soluteConcType(const QString&)));
	connect(selectSoluteConcTypeButton_1,SIGNAL(activated(const QString&)),this,SLOT(soluteConcType_defined(const QString&)));
	//
	for(int i=0;i<4;i++){soluteConcLine[tabIndex][i]->setFixedWidth(CONC_SIZE);soluteConcLabel[tabIndex][i]->setFixedWidth(CONC_SIZE);}
	connect(soluteConcLine[tabIndex][0],SIGNAL(returnPressed()),this,SLOT(soluteConc_defined()));
	connect(soluteConcLine[tabIndex][0],SIGNAL(returnPressed()),this,SLOT(label1soluteConc()));
	connect(soluteConcLine[tabIndex][1],SIGNAL(returnPressed()),this,SLOT(label2soluteConc()));
	connect(soluteConcLine[tabIndex][2],SIGNAL(returnPressed()),this,SLOT(label3soluteConc()));
	connect(soluteConcLine[tabIndex][3],SIGNAL(returnPressed()),this,SLOT(label4soluteConc()));

	QGridLayout *layout = new QGridLayout;
	layout->setSizeConstraint(QLayout::SetFixedSize);
	layout->addWidget(addSoluteButton,0,0,1,1);layout->addWidget(remSoluteButton,0,1,1,2);
	layout->addWidget(seleSoluLabel[tabIndex],1,0,1,1);layout->addWidget(seleSoluConcLabel[tabIndex],1,1,1,1);layout->addWidget(seleSoluConcTypeLabel[tabIndex],1,2,1,1);
	layout->addWidget(selectSoluteButton_1,2,0,1,1);layout->addWidget(soluteConcLine[tabIndex][0],2,1,1,1);layout->addWidget(selectSoluteConcTypeButton_1,2,2,1,1);
	layout->addWidget(soluteLabel[tabIndex][0],3,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][0],3,1,1,1);layout->addWidget(soluteConcTypeLabel[tabIndex],3,2,1,1);
	layout->addWidget(selectSoluteButton_2,4,0,1,1);layout->addWidget(soluteConcLine[tabIndex][1],4,1,1,1);
	layout->addWidget(soluteLabel[tabIndex][1],5,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][1],5,1,1,1);
	layout->addWidget(selectSoluteButton_3,6,0,1,1);layout->addWidget(soluteConcLine[tabIndex][2],6,1,1,1);
	layout->addWidget(soluteLabel[tabIndex][2],7,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][2],7,1,1,1);
	layout->addWidget(selectSoluteButton_4,8,0,1,1);layout->addWidget(soluteConcLine[tabIndex][3],8,1,1,1);
	layout->addWidget(soluteLabel[tabIndex][3],9,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][3],9,1,1,1);
	Output->setLayout(layout);
	return Output;}

QGroupBox* UI::create_fourSolventBox(int tabIndex)
	{// Solvent(s) Selection Group
	QGroupBox* Output=new QGroupBox(tr("Solvent(s) Selection"));
	QPushButton* addSolventButton=new QPushButton(tr("Add Solvent"));
	addSolventButton->setFixedWidth(150);
	connect(addSolventButton,SIGNAL(clicked()),this,SLOT(addSolvent()));
	QPushButton* remSolventButton=new QPushButton(tr("Remove Solvent"));
	remSolventButton->setFixedWidth(150);
	connect(remSolventButton,SIGNAL(clicked()),this,SLOT(remSolvent()));
	//QPushButton *selectSolventButton_1 = new QPushButton(tr("Select a Solvent..."));
	string tmp;
	QComboBox *selectSolventButton_1=new QComboBox;
	QComboBox *selectSolventButton_2=new QComboBox;
	QComboBox *selectSolventButton_3=new QComboBox;
	QComboBox *selectSolventButton_4=new QComboBox;
	for(int i=0;i<numModeledSolvents;i++)
		{tmp=modeledSolvents[i];selectSolventButton_1->addItem(tr(tmp.c_str()));selectSolventButton_2->addItem(tr(tmp.c_str()));selectSolventButton_3->addItem(tr(tmp.c_str()));
		selectSolventButton_4->addItem(tr(tmp.c_str()));}
	selectSolventButton_1->setFixedWidth(180);
	selectSolventButton_1->setCurrentIndex(numModeledSolvents-1);
	connect(selectSolventButton_1,SIGNAL(activated(const QString&)),this,SLOT(label1solvent(const QString&)));
	connect(selectSolventButton_2,SIGNAL(activated(const QString&)),this,SLOT(label2solvent(const QString&)));
	connect(selectSolventButton_3,SIGNAL(activated(const QString&)),this,SLOT(label3solvent(const QString&)));
	connect(selectSolventButton_4,SIGNAL(activated(const QString&)),this,SLOT(label4solvent(const QString&)));
	connect(selectSolventButton_1,SIGNAL(activated(const QString&)),this,SLOT(solvent_defined(const QString&)));

	QComboBox *selectSolventConcTypeButton_1=new QComboBox;
	for(int i=0;i<numModeledSolventConcTypes;i++){tmp=modeledSolventConcTypes[i];selectSolventConcTypeButton_1->addItem(tr(tmp.c_str()));}
	selectSolventConcTypeButton_1->setFixedWidth(150);
	connect(selectSolventConcTypeButton_1,SIGNAL(activated(const QString&)),this,SLOT(label1solventConcType(const QString&)));
	connect(selectSolventConcTypeButton_1,SIGNAL(activated(const QString&)),this,SLOT(solventConcType_defined(const QString&)));
	//
	for(int i=0;i<4;i++){solventConcLine[tabIndex][i]->setFixedWidth(CONC_SIZE);solventConcLabel[tabIndex][i]->setFixedWidth(CONC_SIZE);}
	connect(solventConcLine[tabIndex][0],SIGNAL(returnPressed()),this,SLOT(solventConc_defined()));
	connect(solventConcLine[tabIndex][0],SIGNAL(returnPressed()),this,SLOT(label1solventConc()));
	connect(solventConcLine[tabIndex][1],SIGNAL(returnPressed()),this,SLOT(label2solventConc()));
	connect(solventConcLine[tabIndex][2],SIGNAL(returnPressed()),this,SLOT(label3solventConc()));
	connect(solventConcLine[tabIndex][3],SIGNAL(returnPressed()),this,SLOT(label4solventConc()));

	QGridLayout *layout = new QGridLayout;
	layout->setSizeConstraint(QLayout::SetFixedSize);
	layout->addWidget(addSolventButton,0,0,1,1);layout->addWidget(remSolventButton,0,1,1,2);
	layout->addWidget(seleSolvLabel[tabIndex],1,0,1,1);layout->addWidget(seleSolvConcLabel[tabIndex],1,1,1,1);layout->addWidget(seleSolvConcTypeLabel[tabIndex],1,2,1,1);
	layout->addWidget(selectSolventButton_1,2,0,1,1);layout->addWidget(solventConcLine[tabIndex][0],2,1,1,1);layout->addWidget(selectSolventConcTypeButton_1,2,2,1,1);
	layout->addWidget(solventLabel[tabIndex][0],3,0,1,1);layout->addWidget(solventConcLabel[tabIndex][0],3,1,1,1);layout->addWidget(solventConcTypeLabel[tabIndex],3,2,1,1);
	layout->addWidget(selectSolventButton_2,4,0,1,1);layout->addWidget(solventConcLine[tabIndex][1],4,1,1,1);
	layout->addWidget(solventLabel[tabIndex][1],5,0,1,1);layout->addWidget(solventConcLabel[tabIndex][1],5,1,1,1);
	layout->addWidget(selectSolventButton_3,6,0,1,1);layout->addWidget(solventConcLine[tabIndex][2],6,1,1,1);
	layout->addWidget(solventLabel[tabIndex][2],7,0,1,1);layout->addWidget(solventConcLabel[tabIndex][2],7,1,1,1);
	layout->addWidget(selectSolventButton_4,8,0,1,1);layout->addWidget(solventConcLine[tabIndex][3],8,1,1,1);
	layout->addWidget(solventLabel[tabIndex][3],9,0,1,1);layout->addWidget(solventConcLabel[tabIndex][3],9,1,1,1);
	Output->setLayout(layout);
	return Output;}

QGroupBox* UI::create_fiveSoluteBox(int tabIndex)
	{QGroupBox* Output=new QGroupBox(tr("Solute(s) Selection"));
	QPushButton* addSoluteButton=new QPushButton(tr("Add Solute"));
	addSoluteButton->setFixedWidth(150);
	connect(addSoluteButton,SIGNAL(clicked()),this,SLOT(addSolute()));
	QPushButton* remSoluteButton=new QPushButton(tr("Remove Solute"));
	remSoluteButton->setFixedWidth(150);
	connect(remSoluteButton,SIGNAL(clicked()),this,SLOT(remSolute()));
	//
	string tmp,c,bld;
	QComboBox *selectSoluteButton_1 = new QComboBox;
	QComboBox *selectSoluteButton_2 = new QComboBox;
	QComboBox *selectSoluteButton_3 = new QComboBox;
	QComboBox *selectSoluteButton_4 = new QComboBox;
	QComboBox *selectSoluteButton_5 = new QComboBox;
	for(int i=0;i<numModeledSolutes;i++)
		{tmp=modeledSolutes[i];
		bld="";
		for(uint j=0;j<tmp.length();j++)
			{c=tmp[j];
			bld+=fixSubscriptNumbers(c);}
		selectSoluteButton_1->addItem(tr(bld.c_str()));selectSoluteButton_2->addItem(tr(bld.c_str()));
		selectSoluteButton_3->addItem(tr(bld.c_str()));selectSoluteButton_4->addItem(tr(bld.c_str()));selectSoluteButton_5->addItem(tr(bld.c_str()));}
	selectSoluteButton_1->setFixedWidth(150);
	connect(selectSoluteButton_1, SIGNAL(activated(int)),this, SLOT(label1solute(int)));
	connect(selectSoluteButton_2, SIGNAL(activated(int)),this, SLOT(label2solute(int)));
	connect(selectSoluteButton_3, SIGNAL(activated(int)),this, SLOT(label3solute(int)));
	connect(selectSoluteButton_4, SIGNAL(activated(int)),this, SLOT(label4solute(int)));
	connect(selectSoluteButton_5, SIGNAL(activated(int)),this, SLOT(label5solute(int)));
	connect(selectSoluteButton_1,SIGNAL(activated(const QString&)),this,SLOT(solute_defined(const QString&)));

	QComboBox *selectSoluteConcTypeButton_1=new QComboBox;
	selectSoluteConcTypeButton_1->addItem(tr("Molar"));
	selectSoluteConcTypeButton_1->setFixedWidth(100);
	connect(selectSoluteConcTypeButton_1,SIGNAL(activated(const QString&)),this,SLOT(label1soluteConcType(const QString&)));
	connect(selectSoluteConcTypeButton_1,SIGNAL(activated(const QString&)),this,SLOT(soluteConcType_defined(const QString&)));
	//
	for(int i=0;i<5;i++){soluteConcLine[tabIndex][i]->setFixedWidth(CONC_SIZE);soluteConcLabel[tabIndex][i]->setFixedWidth(CONC_SIZE);}
	connect(soluteConcLine[tabIndex][0],SIGNAL(returnPressed()),this,SLOT(soluteConc_defined()));
	connect(soluteConcLine[tabIndex][0],SIGNAL(returnPressed()),this,SLOT(label1soluteConc()));
	connect(soluteConcLine[tabIndex][1],SIGNAL(returnPressed()),this,SLOT(label2soluteConc()));
	connect(soluteConcLine[tabIndex][2],SIGNAL(returnPressed()),this,SLOT(label3soluteConc()));
	connect(soluteConcLine[tabIndex][3],SIGNAL(returnPressed()),this,SLOT(label4soluteConc()));
	connect(soluteConcLine[tabIndex][4],SIGNAL(returnPressed()),this,SLOT(label5soluteConc()));

	QGridLayout *layout = new QGridLayout;
	layout->setSizeConstraint(QLayout::SetFixedSize);
	layout->addWidget(addSoluteButton,0,0,1,1);layout->addWidget(remSoluteButton,0,1,1,2);
	layout->addWidget(seleSoluLabel[tabIndex],1,0,1,1);layout->addWidget(seleSoluConcLabel[tabIndex],1,1,1,1);layout->addWidget(seleSoluConcTypeLabel[tabIndex],1,2,1,1);
	layout->addWidget(selectSoluteButton_1,2,0,1,1);layout->addWidget(soluteConcLine[tabIndex][0],2,1,1,1);layout->addWidget(selectSoluteConcTypeButton_1,2,2,1,1);
	layout->addWidget(soluteLabel[tabIndex][0],3,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][0],3,1,1,1);layout->addWidget(soluteConcTypeLabel[tabIndex],3,2,1,1);
	layout->addWidget(selectSoluteButton_2,4,0,1,1);layout->addWidget(soluteConcLine[tabIndex][1],4,1,1,1);
	layout->addWidget(soluteLabel[tabIndex][1],5,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][1],5,1,1,1);
	layout->addWidget(selectSoluteButton_3,6,0,1,1);layout->addWidget(soluteConcLine[tabIndex][2],6,1,1,1);
	layout->addWidget(soluteLabel[tabIndex][2],7,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][2],7,1,1,1);
	layout->addWidget(selectSoluteButton_4,8,0,1,1);layout->addWidget(soluteConcLine[tabIndex][3],8,1,1,1);
	layout->addWidget(soluteLabel[tabIndex][3],9,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][3],9,1,1,1);
	layout->addWidget(selectSoluteButton_5,10,0,1,1);layout->addWidget(soluteConcLine[tabIndex][4],10,1,1,1);
	layout->addWidget(soluteLabel[tabIndex][4],11,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][4],11,1,1,1);
	Output->setLayout(layout);
	return Output;}

QGroupBox* UI::create_fiveSolventBox(int tabIndex)
	{// Solvent(s) Selection Group
	QGroupBox* Output=new QGroupBox(tr("Solvent(s) Selection"));
	QPushButton* addSolventButton=new QPushButton(tr("Add Solvent"));
	addSolventButton->setFixedWidth(150);
	connect(addSolventButton,SIGNAL(clicked()),this,SLOT(addSolvent()));
	QPushButton* remSolventButton=new QPushButton(tr("Remove Solvent"));
	remSolventButton->setFixedWidth(150);
	connect(remSolventButton,SIGNAL(clicked()),this,SLOT(remSolvent()));
	//QPushButton *selectSolventButton_1 = new QPushButton(tr("Select a Solvent..."));
	string tmp;
	QComboBox *selectSolventButton_1=new QComboBox;
	QComboBox *selectSolventButton_2=new QComboBox;
	QComboBox *selectSolventButton_3=new QComboBox;
	QComboBox *selectSolventButton_4=new QComboBox;
	QComboBox *selectSolventButton_5=new QComboBox;
	for(int i=0;i<numModeledSolvents;i++)
		{tmp=modeledSolvents[i];selectSolventButton_1->addItem(tr(tmp.c_str()));selectSolventButton_2->addItem(tr(tmp.c_str()));selectSolventButton_3->addItem(tr(tmp.c_str()));
		selectSolventButton_4->addItem(tr(tmp.c_str()));selectSolventButton_5->addItem(tr(tmp.c_str()));}
	selectSolventButton_1->setFixedWidth(180);
	selectSolventButton_1->setCurrentIndex(numModeledSolvents-1);
	connect(selectSolventButton_1,SIGNAL(activated(const QString&)),this,SLOT(label1solvent(const QString&)));
	connect(selectSolventButton_2,SIGNAL(activated(const QString&)),this,SLOT(label2solvent(const QString&)));
	connect(selectSolventButton_3,SIGNAL(activated(const QString&)),this,SLOT(label3solvent(const QString&)));
	connect(selectSolventButton_4,SIGNAL(activated(const QString&)),this,SLOT(label4solvent(const QString&)));
	connect(selectSolventButton_5,SIGNAL(activated(const QString&)),this,SLOT(label5solvent(const QString&)));
	connect(selectSolventButton_1,SIGNAL(activated(const QString&)),this,SLOT(solvent_defined(const QString&)));

	QComboBox *selectSolventConcTypeButton_1=new QComboBox;
	for(int i=0;i<numModeledSolventConcTypes;i++){tmp=modeledSolventConcTypes[i];selectSolventConcTypeButton_1->addItem(tr(tmp.c_str()));}
	selectSolventConcTypeButton_1->setFixedWidth(150);
	connect(selectSolventConcTypeButton_1,SIGNAL(activated(const QString&)),this,SLOT(label1solventConcType(const QString&)));
	connect(selectSolventConcTypeButton_1,SIGNAL(activated(const QString&)),this,SLOT(solventConcType_defined(const QString&)));
	//
	for(int i=0;i<5;i++){solventConcLine[tabIndex][i]->setFixedWidth(CONC_SIZE);solventConcLabel[tabIndex][i]->setFixedWidth(CONC_SIZE);}
	connect(solventConcLine[tabIndex][0],SIGNAL(returnPressed()),this,SLOT(solventConc_defined()));
	connect(solventConcLine[tabIndex][0],SIGNAL(returnPressed()),this,SLOT(label1solventConc()));
	connect(solventConcLine[tabIndex][1],SIGNAL(returnPressed()),this,SLOT(label2solventConc()));
	connect(solventConcLine[tabIndex][2],SIGNAL(returnPressed()),this,SLOT(label3solventConc()));
	connect(solventConcLine[tabIndex][3],SIGNAL(returnPressed()),this,SLOT(label4solventConc()));
	connect(solventConcLine[tabIndex][4],SIGNAL(returnPressed()),this,SLOT(label5solventConc()));

	QGridLayout *layout = new QGridLayout;
	layout->setSizeConstraint(QLayout::SetFixedSize);
	layout->addWidget(addSolventButton,0,0,1,1);layout->addWidget(remSolventButton,0,1,1,2);
	layout->addWidget(seleSolvLabel[tabIndex],1,0,1,1);layout->addWidget(seleSolvConcLabel[tabIndex],1,1,1,1);layout->addWidget(seleSolvConcTypeLabel[tabIndex],1,2,1,1);
	layout->addWidget(selectSolventButton_1,2,0,1,1);layout->addWidget(solventConcLine[tabIndex][0],2,1,1,1);layout->addWidget(selectSolventConcTypeButton_1,2,2,1,1);
	layout->addWidget(solventLabel[tabIndex][0],3,0,1,1);layout->addWidget(solventConcLabel[tabIndex][0],3,1,1,1);layout->addWidget(solventConcTypeLabel[tabIndex],3,2,1,1);
	layout->addWidget(selectSolventButton_2,4,0,1,1);layout->addWidget(solventConcLine[tabIndex][1],4,1,1,1);
	layout->addWidget(solventLabel[tabIndex][1],5,0,1,1);layout->addWidget(solventConcLabel[tabIndex][1],5,1,1,1);
	layout->addWidget(selectSolventButton_3,6,0,1,1);layout->addWidget(solventConcLine[tabIndex][2],6,1,1,1);
	layout->addWidget(solventLabel[tabIndex][2],7,0,1,1);layout->addWidget(solventConcLabel[tabIndex][2],7,1,1,1);
	layout->addWidget(selectSolventButton_4,8,0,1,1);layout->addWidget(solventConcLine[tabIndex][3],8,1,1,1);
	layout->addWidget(solventLabel[tabIndex][3],9,0,1,1);layout->addWidget(solventConcLabel[tabIndex][3],9,1,1,1);
	layout->addWidget(selectSolventButton_5,10,0,1,1);layout->addWidget(solventConcLine[tabIndex][4],10,1,1,1);
	layout->addWidget(solventLabel[tabIndex][4],11,0,1,1);layout->addWidget(solventConcLabel[tabIndex][4],11,1,1,1);
	Output->setLayout(layout);
	return Output;}

QGroupBox* UI::create_sixSoluteBox(int tabIndex)
	{QGroupBox* Output=new QGroupBox(tr("Solute(s) Selection"));
	QPushButton* addSoluteButton=new QPushButton(tr("Add Solute"));
	addSoluteButton->setFixedWidth(150);
	connect(addSoluteButton,SIGNAL(clicked()),this,SLOT(addSolute()));
	QPushButton* remSoluteButton=new QPushButton(tr("Remove Solute"));
	remSoluteButton->setFixedWidth(150);
	connect(remSoluteButton,SIGNAL(clicked()),this,SLOT(remSolute()));
	//
	string tmp,c,bld;
	QComboBox *selectSoluteButton_1 = new QComboBox;
	QComboBox *selectSoluteButton_2 = new QComboBox;
	QComboBox *selectSoluteButton_3 = new QComboBox;
	QComboBox *selectSoluteButton_4 = new QComboBox;
	QComboBox *selectSoluteButton_5 = new QComboBox;
	QComboBox *selectSoluteButton_6 = new QComboBox;
	for(int i=0;i<numModeledSolutes;i++)
		{tmp=modeledSolutes[i];
		bld="";
		for(uint j=0;j<tmp.length();j++)
			{c=tmp[j];
			bld+=fixSubscriptNumbers(c);}
		selectSoluteButton_1->addItem(tr(bld.c_str()));selectSoluteButton_2->addItem(tr(bld.c_str()));
		selectSoluteButton_3->addItem(tr(bld.c_str()));selectSoluteButton_4->addItem(tr(bld.c_str()));selectSoluteButton_5->addItem(tr(bld.c_str()));
		selectSoluteButton_6->addItem(tr(bld.c_str()));}
	selectSoluteButton_1->setFixedWidth(150);
	connect(selectSoluteButton_1, SIGNAL(activated(int)),this, SLOT(label1solute(int)));
	connect(selectSoluteButton_2, SIGNAL(activated(int)),this, SLOT(label2solute(int)));
	connect(selectSoluteButton_3, SIGNAL(activated(int)),this, SLOT(label3solute(int)));
	connect(selectSoluteButton_4, SIGNAL(activated(int)),this, SLOT(label4solute(int)));
	connect(selectSoluteButton_5, SIGNAL(activated(int)),this, SLOT(label5solute(int)));
	connect(selectSoluteButton_6, SIGNAL(activated(int)),this, SLOT(label6solute(int)));
	connect(selectSoluteButton_1,SIGNAL(activated(const QString&)),this,SLOT(solute_defined(const QString&)));

	QComboBox *selectSoluteConcTypeButton_1=new QComboBox;
	selectSoluteConcTypeButton_1->addItem(tr("Molar"));
	selectSoluteConcTypeButton_1->setFixedWidth(100);
	connect(selectSoluteConcTypeButton_1,SIGNAL(activated(const QString&)),this,SLOT(label1soluteConcType(const QString&)));
	connect(selectSoluteConcTypeButton_1,SIGNAL(activated(const QString&)),this,SLOT(solute_defined(const QString&)));
	//
	for(int i=0;i<6;i++){soluteConcLine[tabIndex][i]->setFixedWidth(CONC_SIZE);soluteConcLabel[tabIndex][i]->setFixedWidth(CONC_SIZE);}
	connect(soluteConcLine[tabIndex][0],SIGNAL(returnPressed()),this,SLOT(soluteConc_defined()));
	connect(soluteConcLine[tabIndex][0],SIGNAL(returnPressed()),this,SLOT(label1soluteConc()));
	connect(soluteConcLine[tabIndex][1],SIGNAL(returnPressed()),this,SLOT(label2soluteConc()));
	connect(soluteConcLine[tabIndex][2],SIGNAL(returnPressed()),this,SLOT(label3soluteConc()));
	connect(soluteConcLine[tabIndex][3],SIGNAL(returnPressed()),this,SLOT(label4soluteConc()));
	connect(soluteConcLine[tabIndex][4],SIGNAL(returnPressed()),this,SLOT(label5soluteConc()));
	connect(soluteConcLine[tabIndex][5],SIGNAL(returnPressed()),this,SLOT(label6soluteConc()));

	QGridLayout *layout = new QGridLayout;
	layout->setSizeConstraint(QLayout::SetFixedSize);
	layout->addWidget(addSoluteButton,0,0,1,1);layout->addWidget(remSoluteButton,0,1,1,2);
	layout->addWidget(seleSoluLabel[tabIndex],1,0,1,1);layout->addWidget(seleSoluConcLabel[tabIndex],1,1,1,1);layout->addWidget(seleSoluConcTypeLabel[tabIndex],1,2,1,1);
	layout->addWidget(selectSoluteButton_1,2,0,1,1);layout->addWidget(soluteConcLine[tabIndex][0],2,1,1,1);layout->addWidget(selectSoluteConcTypeButton_1,2,2,1,1);
	layout->addWidget(soluteLabel[tabIndex][0],3,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][0],3,1,1,1);layout->addWidget(soluteConcTypeLabel[tabIndex],3,2,1,1);
	layout->addWidget(selectSoluteButton_2,4,0,1,1);layout->addWidget(soluteConcLine[tabIndex][1],4,1,1,1);
	layout->addWidget(soluteLabel[tabIndex][1],5,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][1],5,1,1,1);
	layout->addWidget(selectSoluteButton_3,6,0,1,1);layout->addWidget(soluteConcLine[tabIndex][2],6,1,1,1);
	layout->addWidget(soluteLabel[tabIndex][2],7,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][2],7,1,1,1);
	layout->addWidget(selectSoluteButton_4,8,0,1,1);layout->addWidget(soluteConcLine[tabIndex][3],8,1,1,1);
	layout->addWidget(soluteLabel[tabIndex][3],9,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][3],9,1,1,1);
	layout->addWidget(selectSoluteButton_5,10,0,1,1);layout->addWidget(soluteConcLine[tabIndex][4],10,1,1,1);
	layout->addWidget(soluteLabel[tabIndex][4],11,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][4],11,1,1,1);
	layout->addWidget(selectSoluteButton_6,12,0,1,1);layout->addWidget(soluteConcLine[tabIndex][5],12,1,1,1);
	layout->addWidget(soluteLabel[tabIndex][5],13,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][5],13,1,1,1);
	Output->setLayout(layout);
	return Output;}

QGroupBox* UI::create_sixSolventBox(int tabIndex)
	{// Solvent(s) Selection Group
	QGroupBox* Output=new QGroupBox(tr("Solvent(s) Selection"));
	QPushButton* addSolventButton=new QPushButton(tr("Add Solvent"));
	addSolventButton->setFixedWidth(150);
	connect(addSolventButton,SIGNAL(clicked()),this,SLOT(addSolvent()));
	QPushButton* remSolventButton=new QPushButton(tr("Remove Solvent"));
	remSolventButton->setFixedWidth(150);
	connect(remSolventButton,SIGNAL(clicked()),this,SLOT(remSolvent()));
	//QPushButton *selectSolventButton_1 = new QPushButton(tr("Select a Solvent..."));
	string tmp;
	QComboBox *selectSolventButton_1=new QComboBox;
	QComboBox *selectSolventButton_2=new QComboBox;
	QComboBox *selectSolventButton_3=new QComboBox;
	QComboBox *selectSolventButton_4=new QComboBox;
	QComboBox *selectSolventButton_5=new QComboBox;
	QComboBox *selectSolventButton_6=new QComboBox;
	for(int i=0;i<numModeledSolvents;i++)
		{tmp=modeledSolvents[i];selectSolventButton_1->addItem(tr(tmp.c_str()));selectSolventButton_2->addItem(tr(tmp.c_str()));selectSolventButton_3->addItem(tr(tmp.c_str()));
		selectSolventButton_4->addItem(tr(tmp.c_str()));selectSolventButton_5->addItem(tr(tmp.c_str()));selectSolventButton_6->addItem(tr(tmp.c_str()));}
	selectSolventButton_1->setFixedWidth(180);
	selectSolventButton_1->setCurrentIndex(numModeledSolvents-1);
	connect(selectSolventButton_1,SIGNAL(activated(const QString&)),this,SLOT(label1solvent(const QString&)));
	connect(selectSolventButton_2,SIGNAL(activated(const QString&)),this,SLOT(label2solvent(const QString&)));
	connect(selectSolventButton_3,SIGNAL(activated(const QString&)),this,SLOT(label3solvent(const QString&)));
	connect(selectSolventButton_4,SIGNAL(activated(const QString&)),this,SLOT(label4solvent(const QString&)));
	connect(selectSolventButton_5,SIGNAL(activated(const QString&)),this,SLOT(label5solvent(const QString&)));
	connect(selectSolventButton_6,SIGNAL(activated(const QString&)),this,SLOT(label6solvent(const QString&)));
	connect(selectSolventButton_1,SIGNAL(activated(const QString&)),this,SLOT(solvent_defined(const QString&)));

	QComboBox *selectSolventConcTypeButton_1=new QComboBox;
	for(int i=0;i<numModeledSolventConcTypes;i++){tmp=modeledSolventConcTypes[i];selectSolventConcTypeButton_1->addItem(tr(tmp.c_str()));}
	selectSolventConcTypeButton_1->setFixedWidth(150);
	connect(selectSolventConcTypeButton_1,SIGNAL(activated(const QString&)),this,SLOT(label1solventConcType(const QString&)));
	connect(selectSolventConcTypeButton_1,SIGNAL(activated(const QString&)),this,SLOT(solventConcType_defined(const QString&)));
	//
	for(int i=0;i<6;i++){solventConcLine[tabIndex][i]->setFixedWidth(CONC_SIZE);solventConcLabel[tabIndex][i]->setFixedWidth(CONC_SIZE);}
	//
	connect(solventConcLine[tabIndex][0],SIGNAL(returnPressed()),this,SLOT(solventConc_defined()));
	connect(solventConcLine[tabIndex][0],SIGNAL(returnPressed()),this,SLOT(label1solventConc()));
	connect(solventConcLine[tabIndex][1],SIGNAL(returnPressed()),this,SLOT(label2solventConc()));
	connect(solventConcLine[tabIndex][2],SIGNAL(returnPressed()),this,SLOT(label3solventConc()));
	connect(solventConcLine[tabIndex][3],SIGNAL(returnPressed()),this,SLOT(label4solventConc()));
	connect(solventConcLine[tabIndex][4],SIGNAL(returnPressed()),this,SLOT(label5solventConc()));
	connect(solventConcLine[tabIndex][5],SIGNAL(returnPressed()),this,SLOT(label6solventConc()));

	QGridLayout *layout = new QGridLayout;
	layout->setSizeConstraint(QLayout::SetFixedSize);
	layout->addWidget(addSolventButton,0,0,1,1);layout->addWidget(remSolventButton,0,1,1,2);
	layout->addWidget(seleSolvLabel[tabIndex],1,0,1,1);layout->addWidget(seleSolvConcLabel[tabIndex],1,1,1,1);layout->addWidget(seleSolvConcTypeLabel[tabIndex],1,2,1,1);
	layout->addWidget(selectSolventButton_1,2,0,1,1);layout->addWidget(solventConcLine[tabIndex][0],2,1,1,1);layout->addWidget(selectSolventConcTypeButton_1,2,2,1,1);
	layout->addWidget(solventLabel[tabIndex][0],3,0,1,1);layout->addWidget(solventConcLabel[tabIndex][0],3,1,1,1);layout->addWidget(solventConcTypeLabel[tabIndex],3,2,1,1);
	layout->addWidget(selectSolventButton_2,4,0,1,1);layout->addWidget(solventConcLine[tabIndex][1],4,1,1,1);
	layout->addWidget(solventLabel[tabIndex][1],5,0,1,1);layout->addWidget(solventConcLabel[tabIndex][1],5,1,1,1);
	layout->addWidget(selectSolventButton_3,6,0,1,1);layout->addWidget(solventConcLine[tabIndex][2],6,1,1,1);
	layout->addWidget(solventLabel[tabIndex][2],7,0,1,1);layout->addWidget(solventConcLabel[tabIndex][2],7,1,1,1);
	layout->addWidget(selectSolventButton_4,8,0,1,1);layout->addWidget(solventConcLine[tabIndex][3],8,1,1,1);
	layout->addWidget(solventLabel[tabIndex][3],9,0,1,1);layout->addWidget(solventConcLabel[tabIndex][3],9,1,1,1);
	layout->addWidget(selectSolventButton_5,10,0,1,1);layout->addWidget(solventConcLine[tabIndex][4],10,1,1,1);
	layout->addWidget(solventLabel[tabIndex][4],11,0,1,1);layout->addWidget(solventConcLabel[tabIndex][4],11,1,1,1);
	layout->addWidget(selectSolventButton_6,12,0,1,1);layout->addWidget(solventConcLine[tabIndex][5],12,1,1,1);
	layout->addWidget(solventLabel[tabIndex][5],13,0,1,1);layout->addWidget(solventConcLabel[tabIndex][5],13,1,1,1);
	Output->setLayout(layout);
	return Output;}

QGroupBox* UI::create_sevenSoluteBox(int tabIndex)
	{QGroupBox* Output=new QGroupBox(tr("Solute(s) Selection"));
	QPushButton* addSoluteButton=new QPushButton(tr("Add Solute"));
	addSoluteButton->setFixedWidth(150);
	connect(addSoluteButton,SIGNAL(clicked()),this,SLOT(addSolute()));
	QPushButton* remSoluteButton=new QPushButton(tr("Remove Solute"));
	remSoluteButton->setFixedWidth(150);
	connect(remSoluteButton,SIGNAL(clicked()),this,SLOT(remSolute()));
	//
	string tmp,c,bld;
	QComboBox *selectSoluteButton_1 = new QComboBox;
	QComboBox *selectSoluteButton_2 = new QComboBox;
	QComboBox *selectSoluteButton_3 = new QComboBox;
	QComboBox *selectSoluteButton_4 = new QComboBox;
	QComboBox *selectSoluteButton_5 = new QComboBox;
	QComboBox *selectSoluteButton_6 = new QComboBox;
	QComboBox *selectSoluteButton_7 = new QComboBox;
	for(int i=0;i<numModeledSolutes;i++)
		{tmp=modeledSolutes[i];
		bld="";
		for(uint j=0;j<tmp.length();j++)
			{c=tmp[j];
			bld+=fixSubscriptNumbers(c);}
		selectSoluteButton_1->addItem(tr(bld.c_str()));selectSoluteButton_2->addItem(tr(bld.c_str()));
		selectSoluteButton_3->addItem(tr(bld.c_str()));selectSoluteButton_4->addItem(tr(bld.c_str()));selectSoluteButton_5->addItem(tr(bld.c_str()));
		selectSoluteButton_6->addItem(tr(bld.c_str()));selectSoluteButton_7->addItem(tr(bld.c_str()));}
	selectSoluteButton_1->setFixedWidth(150);
	connect(selectSoluteButton_1, SIGNAL(activated(int)),this, SLOT(label1solute(int)));
	connect(selectSoluteButton_2, SIGNAL(activated(int)),this, SLOT(label2solute(int)));
	connect(selectSoluteButton_3, SIGNAL(activated(int)),this, SLOT(label3solute(int)));
	connect(selectSoluteButton_4, SIGNAL(activated(int)),this, SLOT(label4solute(int)));
	connect(selectSoluteButton_5, SIGNAL(activated(int)),this, SLOT(label5solute(int)));
	connect(selectSoluteButton_6, SIGNAL(activated(int)),this, SLOT(label6solute(int)));
	connect(selectSoluteButton_7, SIGNAL(activated(int)),this, SLOT(label7solute(int)));

	connect(selectSoluteButton_1,SIGNAL(activated(const QString&)),this,SLOT(solute_defined(const QString&)));

	QComboBox *selectSoluteConcTypeButton_1=new QComboBox;
	selectSoluteConcTypeButton_1->addItem(tr("Molar"));
	selectSoluteConcTypeButton_1->setFixedWidth(100);
	connect(selectSoluteConcTypeButton_1,SIGNAL(activated(const QString&)),this,SLOT(label1soluteConcType(const QString&)));
	connect(selectSoluteConcTypeButton_1,SIGNAL(activated(const QString&)),this,SLOT(soluteConcType_defined(const QString&)));
	//
	for(int i=0;i<7;i++){soluteConcLine[tabIndex][i]->setFixedWidth(CONC_SIZE);soluteConcLabel[tabIndex][i]->setFixedWidth(CONC_SIZE);}
	connect(soluteConcLine[tabIndex][0],SIGNAL(returnPressed()),this,SLOT(soluteConc_defined()));
	connect(soluteConcLine[tabIndex][0],SIGNAL(returnPressed()),this,SLOT(label1soluteConc()));
	connect(soluteConcLine[tabIndex][1],SIGNAL(returnPressed()),this,SLOT(label2soluteConc()));
	connect(soluteConcLine[tabIndex][2],SIGNAL(returnPressed()),this,SLOT(label3soluteConc()));
	connect(soluteConcLine[tabIndex][3],SIGNAL(returnPressed()),this,SLOT(label4soluteConc()));
	connect(soluteConcLine[tabIndex][4],SIGNAL(returnPressed()),this,SLOT(label5soluteConc()));
	connect(soluteConcLine[tabIndex][5],SIGNAL(returnPressed()),this,SLOT(label6soluteConc()));
	connect(soluteConcLine[tabIndex][6],SIGNAL(returnPressed()),this,SLOT(label7soluteConc()));

	QGridLayout *layout = new QGridLayout;
	layout->setSizeConstraint(QLayout::SetFixedSize);
	layout->addWidget(addSoluteButton,0,0,1,1);layout->addWidget(remSoluteButton,0,1,1,2);
	layout->addWidget(seleSoluLabel[tabIndex],1,0,1,1);layout->addWidget(seleSoluConcLabel[tabIndex],1,1,1,1);layout->addWidget(seleSoluConcTypeLabel[tabIndex],1,2,1,1);
	layout->addWidget(selectSoluteButton_1,2,0,1,1);layout->addWidget(soluteConcLine[tabIndex][0],2,1,1,1);layout->addWidget(selectSoluteConcTypeButton_1,2,2,1,1);
	layout->addWidget(soluteLabel[tabIndex][0],3,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][0],3,1,1,1);layout->addWidget(soluteConcTypeLabel[tabIndex],3,2,1,1);
	layout->addWidget(selectSoluteButton_2,4,0,1,1);layout->addWidget(soluteConcLine[tabIndex][1],4,1,1,1);
	layout->addWidget(soluteLabel[tabIndex][1],5,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][1],5,1,1,1);
	layout->addWidget(selectSoluteButton_3,6,0,1,1);layout->addWidget(soluteConcLine[tabIndex][2],6,1,1,1);
	layout->addWidget(soluteLabel[tabIndex][2],7,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][2],7,1,1,1);
	layout->addWidget(selectSoluteButton_4,8,0,1,1);layout->addWidget(soluteConcLine[tabIndex][3],8,1,1,1);
	layout->addWidget(soluteLabel[tabIndex][3],9,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][3],9,1,1,1);
	layout->addWidget(selectSoluteButton_5,10,0,1,1);layout->addWidget(soluteConcLine[tabIndex][4],10,1,1,1);
	layout->addWidget(soluteLabel[tabIndex][4],11,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][4],11,1,1,1);
	layout->addWidget(selectSoluteButton_6,12,0,1,1);layout->addWidget(soluteConcLine[tabIndex][5],12,1,1,1);
	layout->addWidget(soluteLabel[tabIndex][5],13,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][5],13,1,1,1);
	layout->addWidget(selectSoluteButton_7,14,0,1,1);layout->addWidget(soluteConcLine[tabIndex][6],14,1,1,1);
	layout->addWidget(soluteLabel[tabIndex][6],15,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][6],15,1,1,1);
	Output->setLayout(layout);
	return Output;}

QGroupBox* UI::create_sevenSolventBox(int tabIndex)
	{// Solvent(s) Selection Group
	QGroupBox* Output=new QGroupBox(tr("Solvent(s) Selection"));
	QPushButton* addSolventButton=new QPushButton(tr("Add Solvent"));
	addSolventButton->setFixedWidth(150);
	connect(addSolventButton,SIGNAL(clicked()),this,SLOT(addSolvent()));
	QPushButton* remSolventButton=new QPushButton(tr("Remove Solvent"));
	remSolventButton->setFixedWidth(150);
	connect(remSolventButton,SIGNAL(clicked()),this,SLOT(remSolvent()));
	//QPushButton *selectSolventButton_1 = new QPushButton(tr("Select a Solvent..."));
	string tmp;
	QComboBox *selectSolventButton_1=new QComboBox;
	QComboBox *selectSolventButton_2=new QComboBox;
	QComboBox *selectSolventButton_3=new QComboBox;
	QComboBox *selectSolventButton_4=new QComboBox;
	QComboBox *selectSolventButton_5=new QComboBox;
	QComboBox *selectSolventButton_6=new QComboBox;
	QComboBox *selectSolventButton_7=new QComboBox;
	for(int i=0;i<numModeledSolvents;i++)
		{tmp=modeledSolvents[i];selectSolventButton_1->addItem(tr(tmp.c_str()));selectSolventButton_2->addItem(tr(tmp.c_str()));selectSolventButton_3->addItem(tr(tmp.c_str()));
		selectSolventButton_4->addItem(tr(tmp.c_str()));selectSolventButton_5->addItem(tr(tmp.c_str()));selectSolventButton_6->addItem(tr(tmp.c_str()));selectSolventButton_7->addItem(tr(tmp.c_str()));}
	selectSolventButton_1->setFixedWidth(180);
	selectSolventButton_1->setCurrentIndex(numModeledSolvents-1);
	connect(selectSolventButton_1,SIGNAL(activated(const QString&)),this,SLOT(label1solvent(const QString&)));
	connect(selectSolventButton_2,SIGNAL(activated(const QString&)),this,SLOT(label2solvent(const QString&)));
	connect(selectSolventButton_3,SIGNAL(activated(const QString&)),this,SLOT(label3solvent(const QString&)));
	connect(selectSolventButton_4,SIGNAL(activated(const QString&)),this,SLOT(label4solvent(const QString&)));
	connect(selectSolventButton_5,SIGNAL(activated(const QString&)),this,SLOT(label5solvent(const QString&)));
	connect(selectSolventButton_6,SIGNAL(activated(const QString&)),this,SLOT(label6solvent(const QString&)));
	connect(selectSolventButton_7,SIGNAL(activated(const QString&)),this,SLOT(label7solvent(const QString&)));
	connect(selectSolventButton_1,SIGNAL(activated(const QString&)),this,SLOT(solvent_defined(const QString&)));
	
	QComboBox *selectSolventConcTypeButton_1=new QComboBox;
	for(int i=0;i<numModeledSolventConcTypes;i++){tmp=modeledSolventConcTypes[i];selectSolventConcTypeButton_1->addItem(tr(tmp.c_str()));}
	selectSolventConcTypeButton_1->setFixedWidth(150);
	connect(selectSolventConcTypeButton_1,SIGNAL(activated(const QString&)),this,SLOT(label1solventConcType(const QString&)));
	connect(selectSolventConcTypeButton_1,SIGNAL(activated(const QString&)),this,SLOT(solventConcType_defined(const QString&)));
	//
	for(int i=0;i<7;i++){solventConcLine[tabIndex][i]->setFixedWidth(CONC_SIZE);solventConcLabel[tabIndex][i]->setFixedWidth(CONC_SIZE);}
	//
	connect(solventConcLine[tabIndex][0],SIGNAL(returnPressed()),this,SLOT(solventConc_defined()));
	connect(solventConcLine[tabIndex][0],SIGNAL(returnPressed()),this,SLOT(label1solventConc()));
	connect(solventConcLine[tabIndex][1],SIGNAL(returnPressed()),this,SLOT(label2solventConc()));
	connect(solventConcLine[tabIndex][2],SIGNAL(returnPressed()),this,SLOT(label3solventConc()));
	connect(solventConcLine[tabIndex][3],SIGNAL(returnPressed()),this,SLOT(label4solventConc()));
	connect(solventConcLine[tabIndex][4],SIGNAL(returnPressed()),this,SLOT(label5solventConc()));
	connect(solventConcLine[tabIndex][5],SIGNAL(returnPressed()),this,SLOT(label6solventConc()));
	connect(solventConcLine[tabIndex][6],SIGNAL(returnPressed()),this,SLOT(label7solventConc()));

	QGridLayout *layout = new QGridLayout;
	layout->setSizeConstraint(QLayout::SetFixedSize);
	layout->addWidget(addSolventButton,0,0,1,1);layout->addWidget(remSolventButton,0,1,1,2);
	layout->addWidget(seleSolvLabel[tabIndex],1,0,1,1);layout->addWidget(seleSolvConcLabel[tabIndex],1,1,1,1);layout->addWidget(seleSolvConcTypeLabel[tabIndex],1,2,1,1);
	layout->addWidget(selectSolventButton_1,2,0,1,1);layout->addWidget(solventConcLine[tabIndex][0],2,1,1,1);layout->addWidget(selectSolventConcTypeButton_1,2,2,1,1);
	layout->addWidget(solventLabel[tabIndex][0],3,0,1,1);layout->addWidget(solventConcLabel[tabIndex][0],3,1,1,1);layout->addWidget(solventConcTypeLabel[tabIndex],3,2,1,1);
	layout->addWidget(selectSolventButton_2,4,0,1,1);layout->addWidget(solventConcLine[tabIndex][1],4,1,1,1);
	layout->addWidget(solventLabel[tabIndex][1],5,0,1,1);layout->addWidget(solventConcLabel[tabIndex][1],5,1,1,1);
	layout->addWidget(selectSolventButton_3,6,0,1,1);layout->addWidget(solventConcLine[tabIndex][2],6,1,1,1);
	layout->addWidget(solventLabel[tabIndex][2],7,0,1,1);layout->addWidget(solventConcLabel[tabIndex][2],7,1,1,1);
	layout->addWidget(selectSolventButton_4,8,0,1,1);layout->addWidget(solventConcLine[tabIndex][3],8,1,1,1);
	layout->addWidget(solventLabel[tabIndex][3],9,0,1,1);layout->addWidget(solventConcLabel[tabIndex][3],9,1,1,1);
	layout->addWidget(selectSolventButton_5,10,0,1,1);layout->addWidget(solventConcLine[tabIndex][4],10,1,1,1);
	layout->addWidget(solventLabel[tabIndex][4],11,0,1,1);layout->addWidget(solventConcLabel[tabIndex][4],11,1,1,1);
	layout->addWidget(selectSolventButton_6,12,0,1,1);layout->addWidget(solventConcLine[tabIndex][5],12,1,1,1);
	layout->addWidget(solventLabel[tabIndex][5],13,0,1,1);layout->addWidget(solventConcLabel[tabIndex][5],13,1,1,1);
	layout->addWidget(selectSolventButton_7,14,0,1,1);layout->addWidget(solventConcLine[tabIndex][6],14,1,1,1);
	layout->addWidget(solventLabel[tabIndex][6],15,0,1,1);layout->addWidget(solventConcLabel[tabIndex][6],15,1,1,1);
	Output->setLayout(layout);
	return Output;}

QGroupBox* UI::create_eightSoluteBox(int tabIndex)
	{QGroupBox* Output=new QGroupBox(tr("Solute(s) Selection"));
	QPushButton* addSoluteButton=new QPushButton(tr("Add Solute"));
	addSoluteButton->setFixedWidth(150);
	connect(addSoluteButton,SIGNAL(clicked()),this,SLOT(addSolute()));addSoluteButton->hide();
	QPushButton* remSoluteButton=new QPushButton(tr("Remove Solute"));
	remSoluteButton->setFixedWidth(150);
	connect(remSoluteButton,SIGNAL(clicked()),this,SLOT(remSolute()));

	string tmp,c,bld;
	QComboBox *selectSoluteButton_1 = new QComboBox;
	QComboBox *selectSoluteButton_2 = new QComboBox;
	QComboBox *selectSoluteButton_3 = new QComboBox;
	QComboBox *selectSoluteButton_4 = new QComboBox;
	QComboBox *selectSoluteButton_5 = new QComboBox;
	QComboBox *selectSoluteButton_6 = new QComboBox;
	QComboBox *selectSoluteButton_7 = new QComboBox;
	QComboBox *selectSoluteButton_8 = new QComboBox;
	for(int i=0;i<numModeledSolutes;i++)
		{tmp=modeledSolutes[i];
		bld="";
		for(uint j=0;j<tmp.length();j++)
			{c=tmp[j];
			bld+=fixSubscriptNumbers(c);}
		selectSoluteButton_1->addItem(tr(bld.c_str()));selectSoluteButton_2->addItem(tr(bld.c_str()));
		selectSoluteButton_3->addItem(tr(bld.c_str()));selectSoluteButton_4->addItem(tr(bld.c_str()));selectSoluteButton_5->addItem(tr(bld.c_str()));
		selectSoluteButton_6->addItem(tr(bld.c_str()));selectSoluteButton_7->addItem(tr(bld.c_str()));selectSoluteButton_8->addItem(tr(bld.c_str()));}
	selectSoluteButton_1->setFixedWidth(150);
	// Label Selection
	connect(selectSoluteButton_1, SIGNAL(activated(int)),this, SLOT(label1solute(int)));
	connect(selectSoluteButton_2, SIGNAL(activated(int)),this, SLOT(label2solute(int)));
	connect(selectSoluteButton_3, SIGNAL(activated(int)),this, SLOT(label3solute(int)));
	connect(selectSoluteButton_4, SIGNAL(activated(int)),this, SLOT(label4solute(int)));
	connect(selectSoluteButton_5, SIGNAL(activated(int)),this, SLOT(label5solute(int)));
	connect(selectSoluteButton_6, SIGNAL(activated(int)),this, SLOT(label6solute(int)));
	connect(selectSoluteButton_7, SIGNAL(activated(int)),this, SLOT(label7solute(int)));
	connect(selectSoluteButton_8, SIGNAL(activated(int)),this, SLOT(label8solute(int)));
	// Go Black (indicate a value defined)
	connect(selectSoluteButton_1,SIGNAL(activated(const QString&)),this,SLOT(solute_defined(const QString&)));

	QComboBox *selectSoluteConcTypeButton_1=new QComboBox;
	selectSoluteConcTypeButton_1->addItem(tr("Molar"));
	selectSoluteConcTypeButton_1->setFixedWidth(100);
	connect(selectSoluteConcTypeButton_1,SIGNAL(activated(const QString&)),this,SLOT(label1soluteConcType(const QString&)));
	connect(selectSoluteConcTypeButton_1,SIGNAL(activated(const QString&)),this,SLOT(soluteConcType_defined(const QString&)));
	//
	for(int i=0;i<8;i++){soluteConcLine[tabIndex][i]->setFixedWidth(CONC_SIZE);soluteConcLabel[tabIndex][i]->setFixedWidth(CONC_SIZE);}
	connect(soluteConcLine[tabIndex][0],SIGNAL(returnPressed()),this,SLOT(soluteConc_defined()));
	connect(soluteConcLine[tabIndex][0],SIGNAL(returnPressed()),this,SLOT(label1soluteConc()));
	connect(soluteConcLine[tabIndex][1],SIGNAL(returnPressed()),this,SLOT(label2soluteConc()));
	connect(soluteConcLine[tabIndex][2],SIGNAL(returnPressed()),this,SLOT(label3soluteConc()));
	connect(soluteConcLine[tabIndex][3],SIGNAL(returnPressed()),this,SLOT(label4soluteConc()));
	connect(soluteConcLine[tabIndex][4],SIGNAL(returnPressed()),this,SLOT(label5soluteConc()));
	connect(soluteConcLine[tabIndex][5],SIGNAL(returnPressed()),this,SLOT(label6soluteConc()));
	connect(soluteConcLine[tabIndex][6],SIGNAL(returnPressed()),this,SLOT(label7soluteConc()));
	connect(soluteConcLine[tabIndex][7],SIGNAL(returnPressed()),this,SLOT(label8soluteConc()));
	//
	QGridLayout *layout = new QGridLayout;
	layout->setSizeConstraint(QLayout::SetFixedSize);
	layout->addWidget(addSoluteButton,0,0,1,1);layout->addWidget(remSoluteButton,0,1,1,2);
	layout->addWidget(seleSoluLabel[tabIndex],1,0,1,1);layout->addWidget(seleSoluConcLabel[tabIndex],1,1,1,1);layout->addWidget(seleSoluConcTypeLabel[tabIndex],1,2,1,1);
	layout->addWidget(selectSoluteButton_1,2,0,1,1);layout->addWidget(soluteConcLine[tabIndex][0],2,1,1,1);layout->addWidget(selectSoluteConcTypeButton_1,2,2,1,1);
	layout->addWidget(soluteLabel[tabIndex][0],3,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][0],3,1,1,1);layout->addWidget(soluteConcTypeLabel[tabIndex],3,2,1,1);
	layout->addWidget(selectSoluteButton_2,4,0,1,1);layout->addWidget(soluteConcLine[tabIndex][1],4,1,1,1);
	layout->addWidget(soluteLabel[tabIndex][1],5,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][1],5,1,1,1);
	layout->addWidget(selectSoluteButton_3,6,0,1,1);layout->addWidget(soluteConcLine[tabIndex][2],6,1,1,1);
	layout->addWidget(soluteLabel[tabIndex][2],7,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][2],7,1,1,1);
	layout->addWidget(selectSoluteButton_4,8,0,1,1);layout->addWidget(soluteConcLine[tabIndex][3],8,1,1,1);
	layout->addWidget(soluteLabel[tabIndex][3],9,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][3],9,1,1,1);
	layout->addWidget(selectSoluteButton_5,10,0,1,1);layout->addWidget(soluteConcLine[tabIndex][4],10,1,1,1);
	layout->addWidget(soluteLabel[tabIndex][4],11,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][4],11,1,1,1);
	layout->addWidget(selectSoluteButton_6,12,0,1,1);layout->addWidget(soluteConcLine[tabIndex][5],12,1,1,1);
	layout->addWidget(soluteLabel[tabIndex][5],13,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][5],13,1,1,1);
	layout->addWidget(selectSoluteButton_7,14,0,1,1);layout->addWidget(soluteConcLine[tabIndex][6],14,1,1,1);
	layout->addWidget(soluteLabel[tabIndex][6],15,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][6],15,1,1,1);
	layout->addWidget(selectSoluteButton_8,16,0,1,1);layout->addWidget(soluteConcLine[tabIndex][7],16,1,1,1);
	layout->addWidget(soluteLabel[tabIndex][7],17,0,1,1);layout->addWidget(soluteConcLabel[tabIndex][7],17,1,1,1);
	Output->setLayout(layout);
	return Output;}

QGroupBox* UI::create_eightSolventBox(int tabIndex)
	{// Solvent(s) Selection Group
	QGroupBox* Output=new QGroupBox(tr("Solvent(s) Selection"));
	QPushButton* addSolventButton=new QPushButton(tr("Add Solvent"));
	addSolventButton->setFixedWidth(150);
	connect(addSolventButton,SIGNAL(clicked()),this,SLOT(addSolvent()));addSolventButton->hide();
	QPushButton* remSolventButton=new QPushButton(tr("Remove Solvent"));
	remSolventButton->setFixedWidth(150);
	connect(remSolventButton,SIGNAL(clicked()),this,SLOT(remSolvent()));
	//QPushButton *selectSolventButton_1 = new QPushButton(tr("Select a Solvent..."));
	string tmp;
	QComboBox *selectSolventButton_1=new QComboBox;
	QComboBox *selectSolventButton_2=new QComboBox;
	QComboBox *selectSolventButton_3=new QComboBox;
	QComboBox *selectSolventButton_4=new QComboBox;
	QComboBox *selectSolventButton_5=new QComboBox;
	QComboBox *selectSolventButton_6=new QComboBox;
	QComboBox *selectSolventButton_7=new QComboBox;
	QComboBox *selectSolventButton_8=new QComboBox;
	for(int i=0;i<numModeledSolvents;i++)
		{tmp=modeledSolvents[i];
		selectSolventButton_1->addItem(tr(tmp.c_str()));selectSolventButton_2->addItem(tr(tmp.c_str()));selectSolventButton_3->addItem(tr(tmp.c_str()));
		selectSolventButton_4->addItem(tr(tmp.c_str()));selectSolventButton_5->addItem(tr(tmp.c_str()));selectSolventButton_6->addItem(tr(tmp.c_str()));selectSolventButton_7->addItem(tr(tmp.c_str()));
		selectSolventButton_8->addItem(tr(tmp.c_str()));}
	selectSolventButton_1->setFixedWidth(180);
	selectSolventButton_1->setCurrentIndex(numModeledSolvents-1);
	connect(selectSolventButton_1,SIGNAL(activated(const QString&)),this,SLOT(label1solvent(const QString&)));
	connect(selectSolventButton_2,SIGNAL(activated(const QString&)),this,SLOT(label2solvent(const QString&)));
	connect(selectSolventButton_3,SIGNAL(activated(const QString&)),this,SLOT(label3solvent(const QString&)));
	connect(selectSolventButton_4,SIGNAL(activated(const QString&)),this,SLOT(label4solvent(const QString&)));
	connect(selectSolventButton_5,SIGNAL(activated(const QString&)),this,SLOT(label5solvent(const QString&)));
	connect(selectSolventButton_6,SIGNAL(activated(const QString&)),this,SLOT(label6solvent(const QString&)));
	connect(selectSolventButton_7,SIGNAL(activated(const QString&)),this,SLOT(label7solvent(const QString&)));
	connect(selectSolventButton_8,SIGNAL(activated(const QString&)),this,SLOT(label8solvent(const QString&)));
	//
	QComboBox *selectSolventConcTypeButton_1=new QComboBox;
	for(int i=0;i<numModeledSolventConcTypes;i++){tmp=modeledSolventConcTypes[i];selectSolventConcTypeButton_1->addItem(tr(tmp.c_str()));}
	selectSolventConcTypeButton_1->setFixedWidth(150);
	connect(selectSolventConcTypeButton_1,SIGNAL(activated(const QString&)),this,SLOT(label1solventConcType(const QString&)));
	connect(selectSolventConcTypeButton_1,SIGNAL(activated(const QString&)),this,SLOT(solventConcType_defined(const QString&)));
	//
	for(int i=0;i<8;i++){solventConcLine[tabIndex][i]->setFixedWidth(CONC_SIZE);solventConcLabel[tabIndex][i]->setFixedWidth(CONC_SIZE);}
	//
	connect(solventConcLine[tabIndex][0],SIGNAL(returnPressed()),this,SLOT(solventConc_defined()));
	connect(solventConcLine[tabIndex][0],SIGNAL(returnPressed()),this,SLOT(label1solventConc()));
	connect(solventConcLine[tabIndex][1],SIGNAL(returnPressed()),this,SLOT(label2solventConc()));
	connect(solventConcLine[tabIndex][2],SIGNAL(returnPressed()),this,SLOT(label3solventConc()));
	connect(solventConcLine[tabIndex][3],SIGNAL(returnPressed()),this,SLOT(label4solventConc()));
	connect(solventConcLine[tabIndex][4],SIGNAL(returnPressed()),this,SLOT(label5solventConc()));
	connect(solventConcLine[tabIndex][5],SIGNAL(returnPressed()),this,SLOT(label6solventConc()));
	connect(solventConcLine[tabIndex][6],SIGNAL(returnPressed()),this,SLOT(label7solventConc()));
	connect(solventConcLine[tabIndex][7],SIGNAL(returnPressed()),this,SLOT(label8solventConc()));
	
	QGridLayout *layout = new QGridLayout;
	layout->setSizeConstraint(QLayout::SetFixedSize);
	layout->addWidget(addSolventButton,0,0,1,1);layout->addWidget(remSolventButton,0,1,1,2);
	layout->addWidget(seleSolvLabel[tabIndex],1,0,1,1);layout->addWidget(seleSolvConcLabel[tabIndex],1,1,1,1);layout->addWidget(seleSolvConcTypeLabel[tabIndex],1,2,1,1);
	layout->addWidget(selectSolventButton_1,2,0,1,1);layout->addWidget(solventConcLine[tabIndex][0],2,1,1,1);layout->addWidget(selectSolventConcTypeButton_1,2,2,1,1);
	layout->addWidget(solventLabel[tabIndex][0],3,0,1,1);layout->addWidget(solventConcLabel[tabIndex][0],3,1,1,1);layout->addWidget(solventConcTypeLabel[tabIndex],3,2,1,1);
	layout->addWidget(selectSolventButton_2,4,0,1,1);layout->addWidget(solventConcLine[tabIndex][1],4,1,1,1);
	layout->addWidget(solventLabel[tabIndex][1],5,0,1,1);layout->addWidget(solventConcLabel[tabIndex][1],5,1,1,1);
	layout->addWidget(selectSolventButton_3,6,0,1,1);layout->addWidget(solventConcLine[tabIndex][2],6,1,1,1);
	layout->addWidget(solventLabel[tabIndex][2],7,0,1,1);layout->addWidget(solventConcLabel[tabIndex][2],7,1,1,1);
	layout->addWidget(selectSolventButton_4,8,0,1,1);layout->addWidget(solventConcLine[tabIndex][3],8,1,1,1);
	layout->addWidget(solventLabel[tabIndex][3],9,0,1,1);layout->addWidget(solventConcLabel[tabIndex][3],9,1,1,1);
	layout->addWidget(selectSolventButton_5,10,0,1,1);layout->addWidget(solventConcLine[tabIndex][4],10,1,1,1);
	layout->addWidget(solventLabel[tabIndex][4],11,0,1,1);layout->addWidget(solventConcLabel[tabIndex][4],11,1,1,1);
	layout->addWidget(selectSolventButton_6,12,0,1,1);layout->addWidget(solventConcLine[tabIndex][5],12,1,1,1);
	layout->addWidget(solventLabel[tabIndex][5],13,0,1,1);layout->addWidget(solventConcLabel[tabIndex][5],13,1,1,1);
	layout->addWidget(selectSolventButton_7,14,0,1,1);layout->addWidget(solventConcLine[tabIndex][6],14,1,1,1);
	layout->addWidget(solventLabel[tabIndex][6],15,0,1,1);layout->addWidget(solventConcLabel[tabIndex][6],15,1,1,1);
	layout->addWidget(selectSolventButton_8,16,0,1,1);layout->addWidget(solventConcLine[tabIndex][7],16,1,1,1);
	layout->addWidget(solventLabel[tabIndex][7],17,0,1,1);layout->addWidget(solventConcLabel[tabIndex][7],17,1,1,1);
	Output->setLayout(layout);
	return Output;}

void UI::jobName_defined(){jobNameLabel->setText(tr("<b><span style=\"color:black;\">Job Name</span></b>:"));}

void UI::openAPBS()
	{//QStringList filenames = QFileDialog::getOpenFileNames(this,tr("files"),QDir::currentPath(),tr("All files (*)") );
	string subFldr=apbsFldr+"bin";
	QStringList filenames = QFileDialog::getOpenFileNames(this,tr("files"),tr(subFldr.c_str()),tr("All files (*)") );
	//for(int i=0;i<filenames.count();i++){cout<<i<<" "<<filenames.at(i).toLocal8Bit().constData()<<endl;}
	string tmp;
	if(filenames.count()!=0)
		{tmp=filenames.at(0).toLocal8Bit().constData();
		apbsFilePath->setText(tr(tmp.c_str()));}
	// Get Base Folder
	tmp=apbsFilePath->text().toStdString();
	apbsFldr=getBaseFolder(tmp);
	// Remove / at end
	apbsFldr=apbsFldr.substr(0,apbsFldr.length()-1);
	apbsFldr=getBaseFolder(apbsFldr);	
	updateHistoryFile("A",apbsFldr);
	}

void UI::openHYDROPRO()
	{QStringList filenames = QFileDialog::getOpenFileNames(this,tr("files"),tr(hydroproFldr.c_str()),tr("All files (*)") );
	string tmp;
	if(filenames.count()!=0)
		{tmp=filenames.at(0).toLocal8Bit().constData();
		hFilePath->setText(tr(tmp.c_str()));}
	// Get Base Folder
	tmp=hFilePath->text().toStdString();
	hydroproFldr=getBaseFolder(tmp);
	updateHistoryFile("H",hydroproFldr);
	}

void UI::openMULTIVALUE()
	{string subFldr=apbsFldr+"share/apbs/tools/bin";
	QStringList filenames = QFileDialog::getOpenFileNames(this,tr("files"),tr(subFldr.c_str()),tr("All files (*)") );
	string tmp;
	if(filenames.count()!=0)
		{tmp=filenames.at(0).toLocal8Bit().constData();
		mvFilePath->setText(tr(tmp.c_str()));}
	}

void UI::openMSMS()
	{QStringList filenames = QFileDialog::getOpenFileNames(this,tr("files"),tr(msmsFldr.c_str()),tr("All files (*)") );
	string tmp;
	if(filenames.count()!=0)
		{tmp=filenames.at(0).toLocal8Bit().constData();
		msFilePath->setText(tr(tmp.c_str()));}
	tmp=msFilePath->text().toStdString();
	msmsFldr=getBaseFolder(tmp);
	updateHistoryFile("M",msmsFldr);
	}

void UI::openMSMS_A()
	{QStringList filenames = QFileDialog::getOpenFileNames(this,tr("files"),tr(msmsFldr.c_str()),tr("All files (*)") );
	string tmp;
	if(filenames.count()!=0)
		{tmp=filenames.at(0).toLocal8Bit().constData();
		msaFilePath->setText(tr(tmp.c_str()));}
	tmp=msaFilePath->text().toStdString();
	msmsFldr=getBaseFolder(tmp);
	updateHistoryFile("M",msmsFldr);
	}

void UI::openMSMS_X()
	{QStringList filenames = QFileDialog::getOpenFileNames(this,tr("files"),tr(msmsFldr.c_str()),tr("All files (*)") );
	string tmp;
	if(filenames.count()!=0)
		{tmp=filenames.at(0).toLocal8Bit().constData();
		msxFilePath->setText(tr(tmp.c_str()));}
	tmp=msxFilePath->text().toStdString();
	msmsFldr=getBaseFolder(tmp);
	updateHistoryFile("M",msmsFldr);
	}

void UI::openMSMS_XN()
	{QStringList filenames = QFileDialog::getOpenFileNames(this,tr("files"),tr(msmsFldr.c_str()),tr("All files (*)") );
	string tmp;
	if(filenames.count()!=0)
		{tmp=filenames.at(0).toLocal8Bit().constData();
		msxnFilePath->setText(tr(tmp.c_str()));}
	tmp=msxnFilePath->text().toStdString();
	msmsFldr=getBaseFolder(tmp);
	updateHistoryFile("M",msmsFldr);
	}

void UI::openPDB2PQR()
	{QStringList filenames = QFileDialog::getOpenFileNames(this,tr("files"),tr(pdb2pqrFldr.c_str()),tr("All files (*)") );
	string tmp;
	if(filenames.count()!=0)
		{tmp=filenames.at(0).toLocal8Bit().constData();
		pFilePath->setText(tr(tmp.c_str()));}
	tmp=pFilePath->text().toStdString();
	pdb2pqrFldr=getBaseFolder(tmp);
	updateHistoryFile("P",pdb2pqrFldr);
	}

void UI::openZPRED()
	{QStringList filenames = QFileDialog::getOpenFileNames(this,tr("files"),tr(zpredFldr.c_str()),tr("All files (*)") );
	string tmp;
	if(filenames.count()!=0)
		{tmp=filenames.at(0).toLocal8Bit().constData();
		zFilePath->setText(tr(tmp.c_str()));}
	tmp=zFilePath->text().toStdString();
	zpredFldr=getBaseFolder(tmp);
	updateHistoryFile("Z",zpredFldr);
	}

void UI::openPdbFiles()
	{QStringList files=QFileDialog::getOpenFileNames(this,tr("files"),tr(filesFldr.c_str()),tr("All files (*)") );
	//for(int i=0;i<files.count();i++){cout<<i<<" "<<files.at(i).toLocal8Bit().constData()<<endl;}
	// Number of PDB Files Input by USer
	numPDBs=files.count();
	// Array of PDB File Names
	pdbFiles=new string[numPDBs];
	for(int i=0;i<numPDBs;i++){pdbFiles[i]=files.at(i).toLocal8Bit().constData();}
	string tmp="<span style=\"color:black;\">";
	tmp+=cnvrtNumToStrng(numPDBs,0);
	tmp+=" structure files will be used by ZPRED.</span>";
	fileLabel->setText(tmp.c_str());
	viewInputButton->show();
	tmp=pdbFiles[0];
	filesFldr=getBaseFolder(tmp);
	updateHistoryFile("F",filesFldr);
	}

void UI::remSolCondTab()
	{// Reduce Number of Solution Conditions
	numSolCond--;
	// Copy Previous Data
	string solvLst="",soluLst="",delimiter=";";
	for(int i=0;i<numSolCond;i++)
		{solvLst+=cnvrtNumToStrng(numSolvent[i],0)+delimiter;
		soluLst+=cnvrtNumToStrng(numSolute[i],0)+delimiter;}
	numSolvent=fill_int_array(solvLst,numSolCond,delimiter);
	numSolute=fill_int_array(soluLst,numSolCond,delimiter);
	if(numSolCond<=1)
		{// Hide Removal Button
		remSolCondButton->hide();}
	int N=scTab->count();
	scTab->removeTab(N-1);
	// Erase Data from deleted tab
	pH[N-1]->setText(tr(NO_VALUE));pHLabel[N-1]->setText(tr("<b><span style=\"color:red;\">pH</span></b>:"));
	T[N-1]->setText(tr(NO_VALUE));TLabel[N-1]->setText(tr("<b><span style=\"color:red;\">Temperature</span></b>:"));
	pC[N-1]->setText(tr(NO_VALUE));pCLabel[N-1]->setText(tr("<b><span style=\"color:red;\">Protein Concentration</span></b>:"));
	}

void UI::remSolute()
	{int tabIndex=scTab->currentIndex();	
	// Update Current Number of Solutes
	if(numSolute[tabIndex]>1){numSolute[tabIndex]--;}
	// Check if pH Defined
	QString val;//=pH[tabIndex]->text();
	string tmp;//=val.toStdString();
	QLabel* TUnitsLabel=new QLabel(tr("[Kelvin]"));TUnitsLabel->setAlignment(Qt::AlignLeft);
	QLabel* pCUnitsLabel=new QLabel(tr("[g/L]"));pCUnitsLabel->setAlignment(Qt::AlignLeft);
	//
	TLabel[tabIndex]->setFixedWidth(150);
	T[tabIndex]->setFixedWidth(LINE_INPUT_SIZE);
	pH[tabIndex]->setFixedWidth(LINE_INPUT_SIZE);
	pC[tabIndex]->setFixedWidth(LINE_INPUT_SIZE);
	// Solvent(s) Selection Group
	switch(numSolvent[tabIndex])
		{case 1: solventBox=create_oneSolventBox(tabIndex);break;
		case 2: solventBox=create_twoSolventBox(tabIndex);break;
		case 3: solventBox=create_threeSolventBox(tabIndex);break;
		case 4: solventBox=create_fourSolventBox(tabIndex);break;
		case 5: solventBox=create_fiveSolventBox(tabIndex);break;
		case 6: solventBox=create_sixSolventBox(tabIndex);break;
		case 7: solventBox=create_sevenSolventBox(tabIndex);break;
		case 8: solventBox=create_eightSolventBox(tabIndex);break;
		default: break;
		}

	// Solute(s) Selection Group
	switch(numSolute[tabIndex])
		{case 1: soluteBox=create_oneSoluteBox(tabIndex);break;
		case 2: soluteBox=create_twoSoluteBox(tabIndex);break;
		case 3: soluteBox=create_threeSoluteBox(tabIndex);break;
		case 4: soluteBox=create_fourSoluteBox(tabIndex);break;
		case 5: soluteBox=create_fiveSoluteBox(tabIndex);break;
		case 6: soluteBox=create_sixSoluteBox(tabIndex);break;
		case 7: soluteBox=create_sevenSoluteBox(tabIndex);break;
		case 8: soluteBox=create_eightSoluteBox(tabIndex);break;
		default: break;
		}
	// Define Solution Conditions Tab
	QWidget *tab1 = new QWidget;

	QGridLayout *tab1hbox = new QGridLayout;
	tab1hbox->setSizeConstraint(QLayout::SetFixedSize);
	tab1hbox->addWidget(TLabel[tabIndex],0,0,1,1);tab1hbox->addWidget(T[tabIndex],0,1,1,1);tab1hbox->addWidget(TUnitsLabel,0,2,1,1);tab1hbox->addWidget(pHLabel[tabIndex],0,4,1,1);tab1hbox->addWidget(pH[tabIndex],0,5,1,1);
	tab1hbox->addWidget(solventBox,1,0,1,3);tab1hbox->addWidget(soluteBox,1,3,1,3);	
	//tab1hbox->addWidget(create_fBox(tabIndex),2,0,2,3);tab1hbox->addWidget(pCLabel[tabIndex],2,3,1,1);tab1hbox->addWidget(pC[tabIndex],2,4,1,1);tab1hbox->addWidget(pCUnitsLabel,2,5,1,1);
		if(forceParmBox[tabIndex]->isChecked())
		{tab1hbox->addWidget(forceParmBox[tabIndex],2,0,1,3);
		tab1hbox->addWidget(create_fBox(tabIndex),3,0,2,3);tab1hbox->addWidget(pCLabel[tabIndex],3,3,1,1);tab1hbox->addWidget(pC[tabIndex],3,4,1,1);tab1hbox->addWidget(pCUnitsLabel,3,5,1,1);
		tab1hbox->addWidget(potProfileBox[tabIndex],5,0,1,3);}
	else
		{tab1hbox->addWidget(forceParmBox[tabIndex],2,0,1,3);tab1hbox->addWidget(pCLabel[tabIndex],2,3,1,1);tab1hbox->addWidget(pC[tabIndex],2,4,1,1);tab1hbox->addWidget(pCUnitsLabel,2,5,1,1);
		tab1hbox->addWidget(potProfileBox[tabIndex],3,0,1,3);}
	tab1->setLayout(tab1hbox);
	// Rewrite tab
	tmp="&"+cnvrtNumToStrng(tabIndex,0);
	// Clear Previous
	scTab->removeTab(tabIndex);
	scTab->insertTab(tabIndex,tab1,tmp.c_str());
	scTab->setCurrentIndex(tabIndex);
	}

void UI::remSolvent()
	{int tabIndex=scTab->currentIndex();	
	// Update Current Number of Solvents
	if(numSolvent[tabIndex]>1){numSolvent[tabIndex]--;}
	// Check if pH Defined
	QString val;//=pH[tabIndex]->text();
	string tmp;//=val.toStdString();
	QLabel* TUnitsLabel=new QLabel(tr("[Kelvin]"));TUnitsLabel->setAlignment(Qt::AlignLeft);
	QLabel* pCUnitsLabel=new QLabel(tr("[g/L]"));pCUnitsLabel->setAlignment(Qt::AlignLeft);
	//
	TLabel[tabIndex]->setFixedWidth(150);
	T[tabIndex]->setFixedWidth(LINE_INPUT_SIZE);
	pH[tabIndex]->setFixedWidth(LINE_INPUT_SIZE);
	pC[tabIndex]->setFixedWidth(LINE_INPUT_SIZE);
	// Solvent(s) Selection Group
	switch(numSolvent[tabIndex])
		{case 1: solventBox=create_oneSolventBox(tabIndex);break;
		case 2: solventBox=create_twoSolventBox(tabIndex);break;
		case 3: solventBox=create_threeSolventBox(tabIndex);break;
		case 4: solventBox=create_fourSolventBox(tabIndex);break;
		case 5: solventBox=create_fiveSolventBox(tabIndex);break;
		case 6: solventBox=create_sixSolventBox(tabIndex);break;
		case 7: solventBox=create_sevenSolventBox(tabIndex);break;
		case 8: solventBox=create_eightSolventBox(tabIndex);break;
		default: break;
		}

	// Solute(s) Selection Group
	switch(numSolute[tabIndex])
		{case 1: soluteBox=create_oneSoluteBox(tabIndex);break;
		case 2: soluteBox=create_twoSoluteBox(tabIndex);break;
		case 3: soluteBox=create_threeSoluteBox(tabIndex);break;
		case 4: soluteBox=create_fourSoluteBox(tabIndex);break;
		case 5: soluteBox=create_fiveSoluteBox(tabIndex);break;
		case 6: soluteBox=create_sixSoluteBox(tabIndex);break;
		case 7: soluteBox=create_sevenSoluteBox(tabIndex);break;
		case 8: soluteBox=create_eightSoluteBox(tabIndex);break;
		default: break;
		}
	// Define Solution Conditions Tab
	QWidget *tab1 = new QWidget;

	QGridLayout *tab1hbox = new QGridLayout;
	tab1hbox->setSizeConstraint(QLayout::SetFixedSize);
	tab1hbox->addWidget(TLabel[tabIndex],0,0,1,1);tab1hbox->addWidget(T[tabIndex],0,1,1,1);tab1hbox->addWidget(TUnitsLabel,0,2,1,1);tab1hbox->addWidget(pHLabel[tabIndex],0,4,1,1);tab1hbox->addWidget(pH[tabIndex],0,5,1,1);
	tab1hbox->addWidget(solventBox,1,0,1,3);tab1hbox->addWidget(soluteBox,1,3,1,3);	
	//tab1hbox->addWidget(create_fBox(tabIndex),2,0,2,3);tab1hbox->addWidget(pCLabel[tabIndex],2,3,1,1);tab1hbox->addWidget(pC[tabIndex],2,4,1,1);tab1hbox->addWidget(pCUnitsLabel,2,5,1,1);
		if(forceParmBox[tabIndex]->isChecked())
		{tab1hbox->addWidget(forceParmBox[tabIndex],2,0,1,3);
		tab1hbox->addWidget(create_fBox(tabIndex),3,0,2,3);tab1hbox->addWidget(pCLabel[tabIndex],3,3,1,1);tab1hbox->addWidget(pC[tabIndex],3,4,1,1);tab1hbox->addWidget(pCUnitsLabel,3,5,1,1);
		tab1hbox->addWidget(potProfileBox[tabIndex],5,0,1,3);}
	else
		{tab1hbox->addWidget(forceParmBox[tabIndex],2,0,1,3);tab1hbox->addWidget(pCLabel[tabIndex],2,3,1,1);tab1hbox->addWidget(pC[tabIndex],2,4,1,1);tab1hbox->addWidget(pCUnitsLabel,2,5,1,1);
		tab1hbox->addWidget(potProfileBox[tabIndex],3,0,1,3);}
	tab1->setLayout(tab1hbox);
	// Rewrite tab
	tmp="&"+cnvrtNumToStrng(tabIndex,0);
	// Clear Previous
	scTab->removeTab(tabIndex);
	scTab->insertTab(tabIndex,tab1,tmp.c_str());
	scTab->setCurrentIndex(tabIndex);
	}

void UI::write_pymolViewInputFunction_File(string pyFnNm,string outFile)
	{string Output="";
	// Write PyMol Function File
	Output+="#!/usr/bin/python\nfrom pymol import stored\nfrom time import sleep\n\n# Function Definition\n\n";
	Output+="def "+pyFnNm+"():\n\n";
	// Define PDB File(s)	
	for(int i=0;i<numPDBs;i++){Output+="\tpFile"+cnvrtNumToStrng(i+1,0)+"=\""+pdbFiles[i]+"\"\n";}
	Output+="\n";
	// Load PDB File(s)
	for(int i=0;i<numPDBs;i++){Output+="\tcmd.load(pFile"+cnvrtNumToStrng(i+1,0)+")\n";}
	Output+="\n";
	// More Functions
	Output+="\tcmd.do(\"util.cbag\")\n";
	Output+="\tcmd.zoom(\"all\",10)\n";
	Output+="\n";
	// End File
	Output+="\treturn\n\ncmd.extend(\""+pyFnNm+"\","+pyFnNm+")";

	ofstream fOut;
	fOut.open(outFile.c_str());
	if(fOut.fail()){cerr<<"Error in write_pymolViewInputFunction_File!\nCould not open file ("<<outFile<<")\n";}//exit(EXIT_FAILURE);}
	else{fOut<<Output;fOut.close();}
	}

string cnvrtNumToStrng(int Num,int numberAfterDecimalpoint)
	{stringstream ss;
	ss.setf(ios::fixed);
	if(numberAfterDecimalpoint>0)
		{ss.setf(ios::showpoint);}
	ss.precision(numberAfterDecimalpoint);
	ss<<Num;
	return ss.str();}

string cnvrtNumToStrng(unsigned int Num,int numberAfterDecimalpoint)
	{stringstream ss;
	ss.setf(ios::fixed);
	if(numberAfterDecimalpoint>0)
		{ss.setf(ios::showpoint);}
	ss.precision(numberAfterDecimalpoint);
	ss<<Num;
	return ss.str();}

string cnvrtNumToStrng(double Num,int numberAfterDecimalpoint){stringstream ss;ss.setf(ios::fixed);if(numberAfterDecimalpoint>0){ss.setf(ios::showpoint);}ss.precision(numberAfterDecimalpoint);ss<<Num;return ss.str();}

string cnvrtNumToStrng(long double Num,int numberAfterDecimalpoint){stringstream ss;ss.setf(ios::fixed);if(numberAfterDecimalpoint>0){ss.setf(ios::showpoint);}ss.precision(numberAfterDecimalpoint);ss<<Num;return ss.str();}

string cnvrtNumToStrng(float Num,int numberAfterDecimalpoint){stringstream ss;ss.setf(ios::fixed);if(numberAfterDecimalpoint>0){ss.setf(ios::showpoint);}ss.precision(numberAfterDecimalpoint);ss<<Num;return ss.str();}

double* fill_double_array(string Data,int numPnts,string delimiter){double* Output=new double[numPnts];string bld="",tmp="";int Counter=0;for(uint i=0;i<Data.length();i++){tmp=Data[i];if(tmp.compare(delimiter)==0 && Counter<numPnts){Output[Counter]=strtod(bld.c_str(),NULL);Counter++;bld="";}else{bld+=Data[i];}}return Output;}

int* fill_int_array(string Data,int numPnts,string delimiter){int* Output=new int[numPnts];string bld="",tmp="";int Counter=0;for(uint i=0;i<Data.length();i++){tmp=Data[i];if(tmp.compare(delimiter)==0 && Counter<numPnts){Output[Counter]=atoi(bld.c_str());Counter++;bld="";}else{bld+=Data[i];}}return Output;}

string* fill_string_array(string Data,int numPnts,string delimiter){string* Output=new string[numPnts];string bld="",tmp="";int Counter=0;for(uint i=0;i<Data.length();i++){tmp=Data[i];if(tmp.compare(delimiter)==0 && Counter<numPnts){Output[Counter]=bld;Counter++;bld="";}else{bld+=Data[i];}}return Output;}

string fixSubscriptNumbers(string s)
	{string Out;
	if(s.compare("0")==0){Out="₀";}	
	else if(s.compare("1")==0){Out="₁";}
	else if(s.compare("2")==0){Out="₂";}
	else if(s.compare("3")==0){Out="₃";}
	else if(s.compare("4")==0){Out="₄";}
	else if(s.compare("5")==0){Out="₅";}
	else if(s.compare("6")==0){Out="₆";}
	else if(s.compare("7")==0){Out="₇";}
	else if(s.compare("8")==0){Out="₈";}
	else if(s.compare("9")==0){Out="₉";}
	else{Out=s;}
	return Out;}

string getBaseFolder(string f)
	{int pos=f.rfind("/",f.length()-1);
	return f.substr(0,pos)+"/";}

string makeUpperCase(string X){string Output="";char letter;for(uint i=0;i<X.length();i++){letter=X[i];Output+=toupper(letter);}return Output;}

string checkFinalBackSlash(string s)
	{string tmp=s.substr(s.length()-1,1);
	if(tmp.compare("/")!=0){tmp=s+"/";}
	else{tmp=s;}
	return tmp;}

string setStringWidth(string In,int width)
	{int Counter=0;
	string Output="";
	for(uint i=0;i<In.length();i++)
		{if(Counter>=width)
			{Output+="\n";
			Output+=In[i];
			Counter=0;
			}
		else
			{Output+=In[i];
			Counter++;}
		}
	return Output;}

void apbsOptions::atompotChecked()
	{if(atompot->isChecked()){updateParameterFile("atompot","true");}
	else{updateParameterFile("atompot","false");}}

void apbsOptions::chargeChecked()
	{if(charge->isChecked()){updateParameterFile("charge","true");}
	else{updateParameterFile("charge","false");}}

void apbsOptions::dielxChecked()
	{if(dielx->isChecked()){updateParameterFile("dielx","true");}
	else{updateParameterFile("dielx","false");}}

void apbsOptions::dielyChecked()
	{if(diely->isChecked()){updateParameterFile("diely","true");}
	else{updateParameterFile("diely","false");}}

void apbsOptions::dielzChecked()
	{if(dielz->isChecked()){updateParameterFile("dielz","true");}
	else{updateParameterFile("dielz","false");}}

void apbsOptions::edensChecked()
	{if(edens->isChecked()){updateParameterFile("edens","true");}
	else{updateParameterFile("edens","false");}}

void apbsOptions::ivdwChecked()
	{if(ivdw->isChecked()){updateParameterFile("ivdw","true");}
	else{updateParameterFile("ivdw","false");}}

void apbsOptions::kappaChecked()
	{if(kappa->isChecked()){updateParameterFile("kappa","true");}
	else{updateParameterFile("kappa","false");}}

void apbsOptions::lapChecked()
	{if(lap->isChecked()){updateParameterFile("lap","true");}
	else{updateParameterFile("lap","false");}}

void apbsOptions::ndensChecked()
	{if(ndens->isChecked()){updateParameterFile("ndens","true");}
	else{updateParameterFile("ndens","false");}}

void apbsOptions::qdensChecked()
	{if(qdens->isChecked()){updateParameterFile("qdens","true");}
	else{updateParameterFile("qdens","false");}}

void apbsOptions::smolChecked()
	{if(smol->isChecked()){updateParameterFile("smol","true");}
	else{updateParameterFile("smol","false");}}

void apbsOptions::ssplChecked()
	{if(sspl->isChecked()){updateParameterFile("sspl","true");}
	else{updateParameterFile("sspl","false");}}

void apbsOptions::vdwChecked()
	{if(vdw->isChecked()){updateParameterFile("vdw","true");}
	else{updateParameterFile("vdw","false");}}

void apbsOptions::labelDX(const QString qS){string s=qS.toStdString();dx->setText(s.c_str());updateParameterFile("dx",s);}

void apbsOptions::labelDY(const QString qS){string s=qS.toStdString();dy->setText(s.c_str());updateParameterFile("dy",s);}

void apbsOptions::labelDZ(const QString qS){string s=qS.toStdString();dz->setText(s.c_str());updateParameterFile("dz",s);}

void apbsOptions::proteinDielectricDefined()
	{string s=proteinDielectricInput->text().toStdString();
	proteinDielectric->setText(s.c_str());
	updateParameterFile("pdie",s);}

void apbsOptions::defineParameters(string iFile)
	{ifstream fIn;
	fIn.open(iFile.c_str());
	if(fIn.fail()){cerr<<"ERROR in defineParameters!\nInput file could not be opened.\n"<<iFile<<endl;exit(EXIT_FAILURE);}
	uint Sz=15000,pos;
	char Val[Sz];
	string bld,tmp,inputName,inputVal;
	fIn.getline(Val,Sz);
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
				if(inputVal.length()!=0)
					{// Identify Parameter and Load Values
					if(inputName.compare("PDIE")==0)
						{// Protein Dielectric Constant
						proteinDielectric->setText(inputVal.c_str());}
					else if(inputName.compare("DX")==0)
						{// Num X Grid Pnts
						dx->setText(inputVal.c_str());}
					else if(inputName.compare("DY")==0)
						{// Num Y Grid Pnts
						dy->setText(inputVal.c_str());}
					else if(inputName.compare("DZ")==0)
						{// Num Z Grid Pnts
						dz->setText(inputVal.c_str());}
					else if(inputName.compare("ATOMPOT")==0)
						{// Write Type
						if(inputVal.compare("true")==0){atompot->setChecked(true);}
						else if(inputVal.compare("false")==0){atompot->setChecked(false);}
						else{cerr<<"ERROR in defineParameters!\nUnrecognized Parameter ("<<inputName<<") Input:\n"<<inputVal<<endl;exit(EXIT_FAILURE);}
						}
					else if(inputName.compare("CHARGE")==0)
						{// Write Type
						if(inputVal.compare("true")==0){charge->setChecked(true);}
						else if(inputVal.compare("false")==0){charge->setChecked(false);}
						else{cerr<<"ERROR in defineParameters!\nUnrecognized Parameter ("<<inputName<<") Input:\n"<<inputVal<<endl;exit(EXIT_FAILURE);}
						}
					else if(inputName.compare("DIELX")==0)
						{// Write Type
						if(inputVal.compare("true")==0){dielx->setChecked(true);}
						else if(inputVal.compare("false")==0){dielx->setChecked(false);}
						else{cerr<<"ERROR in defineParameters!\nUnrecognized Parameter ("<<inputName<<") Input:\n"<<inputVal<<endl;exit(EXIT_FAILURE);}
						}
					else if(inputName.compare("DIELY")==0)
						{// Write Type
						if(inputVal.compare("true")==0){diely->setChecked(true);}
						else if(inputVal.compare("false")==0){diely->setChecked(false);}
						else{cerr<<"ERROR in defineParameters!\nUnrecognized Parameter ("<<inputName<<") Input:\n"<<inputVal<<endl;exit(EXIT_FAILURE);}
						}
					else if(inputName.compare("DIELZ")==0)
						{// Write Type
						if(inputVal.compare("true")==0){dielz->setChecked(true);}
						else if(inputVal.compare("false")==0){dielz->setChecked(false);}
						else{cerr<<"ERROR in defineParameters!\nUnrecognized Parameter ("<<inputName<<") Input:\n"<<inputVal<<endl;exit(EXIT_FAILURE);}
						}
					else if(inputName.compare("EDENS")==0)
						{// Write Type
						if(inputVal.compare("true")==0){edens->setChecked(true);}
						else if(inputVal.compare("false")==0){edens->setChecked(false);}
						else{cerr<<"ERROR in defineParameters!\nUnrecognized Parameter ("<<inputName<<") Input:\n"<<inputVal<<endl;exit(EXIT_FAILURE);}
						}
					else if(inputName.compare("IVDW")==0)
						{// Write Type
						if(inputVal.compare("true")==0){ivdw->setChecked(true);}
						else if(inputVal.compare("false")==0){ivdw->setChecked(false);}
						else{cerr<<"ERROR in defineParameters!\nUnrecognized Parameter ("<<inputName<<") Input:\n"<<inputVal<<endl;exit(EXIT_FAILURE);}
						}
					else if(inputName.compare("KAPPA")==0)
						{// Write Type
						if(inputVal.compare("true")==0){kappa->setChecked(true);}
						else if(inputVal.compare("false")==0){kappa->setChecked(false);}
						else{cerr<<"ERROR in defineParameters!\nUnrecognized Parameter ("<<inputName<<") Input:\n"<<inputVal<<endl;exit(EXIT_FAILURE);}
						}
					else if(inputName.compare("LAP")==0)
						{// Write Type
						if(inputVal.compare("true")==0){lap->setChecked(true);}
						else if(inputVal.compare("false")==0){lap->setChecked(false);}
						else{cerr<<"ERROR in defineParameters!\nUnrecognized Parameter ("<<inputName<<") Input:\n"<<inputVal<<endl;exit(EXIT_FAILURE);}
						}
					else if(inputName.compare("NDENS")==0)
						{// Write Type
						if(inputVal.compare("true")==0){ndens->setChecked(true);}
						else if(inputVal.compare("false")==0){ndens->setChecked(false);}
						else{cerr<<"ERROR in defineParameters!\nUnrecognized Parameter ("<<inputName<<") Input:\n"<<inputVal<<endl;exit(EXIT_FAILURE);}
						}
					else if(inputName.compare("POT")==0)
						{// Write Type
						if(inputVal.compare("true")==0){pot->setChecked(true);}
						else if(inputVal.compare("false")==0){pot->setChecked(false);}
						else{cerr<<"ERROR in defineParameters!\nUnrecognized Parameter ("<<inputName<<") Input:\n"<<inputVal<<endl;exit(EXIT_FAILURE);}
						}
					else if(inputName.compare("QDENS")==0)
						{// Write Type
						if(inputVal.compare("true")==0){qdens->setChecked(true);}
						else if(inputVal.compare("false")==0){qdens->setChecked(false);}
						else{cerr<<"ERROR in defineParameters!\nUnrecognized Parameter ("<<inputName<<") Input:\n"<<inputVal<<endl;exit(EXIT_FAILURE);}
						}
					else if(inputName.compare("SMOL")==0)
						{// Write Type
						if(inputVal.compare("true")==0){smol->setChecked(true);}
						else if(inputVal.compare("false")==0){smol->setChecked(false);}
						else{cerr<<"ERROR in defineParameters!\nUnrecognized Parameter ("<<inputName<<") Input:\n"<<inputVal<<endl;exit(EXIT_FAILURE);}
						}
					else if(inputName.compare("SSPL")==0)
						{// Write Type
						if(inputVal.compare("true")==0){sspl->setChecked(true);}
						else if(inputVal.compare("false")==0){sspl->setChecked(false);}
						else{cerr<<"ERROR in defineParameters!\nUnrecognized Parameter ("<<inputName<<") Input:\n"<<inputVal<<endl;exit(EXIT_FAILURE);}
						}
					else if(inputName.compare("VDW")==0)
						{// Write Type
						if(inputVal.compare("true")==0){vdw->setChecked(true);}
						else if(inputVal.compare("false")==0){vdw->setChecked(false);}
						else{cerr<<"ERROR in defineParameters!\nUnrecognized Parameter ("<<inputName<<") Input:\n"<<inputVal<<endl;exit(EXIT_FAILURE);}
						}
					//else{cerr<<"ERROR in apbsOptions::defineParameters!\nUnrecognized Parameter ("<<inputName<<")."<<endl;exit(EXIT_FAILURE);}
					}
				}
			}
		fIn.getline(Val,Sz);}
	fIn.close();
	}

apbsOptions::apbsOptions(QWidget *parent) : QWidget(parent)
{	
	QString Fldr=QDir::currentPath();
	string sFldr=Fldr.toStdString();
	sFldr=checkFinalBackSlash(sFldr);
	// File Reference for When Optional Parameters are set (makes .optionalParameters.txt)
	string optRefFile=getBaseFolder(sFldr)+".optionalParameters.txt";
	ifstream fIn;
	ofstream fOut;
	fIn.open(optRefFile.c_str());
	string theTxt="// APBS Parameters\n.pdie:4.0\n.dx:161\n.dy:161\n.dz:161\n.atompot:false\n.charge:false\n.dielx:false\n.diely:false\n.dielz:false\n.edens:false\n.ivdw:false\n.kappa:false\n.lap:false\n.ndens:false\n.pot:true\n.qdens:false\n.smol:false\n.sspl:false\n.vdw:false\n";
	theTxt+="// HYDROPRO Parameters\n.calcType:1\n";
	theTxt+="// MSMS Parameters\n.lowPntDensity:1.00\n.highPntDensity:3.00\n";
	theTxt+="// PDB2PQR Parameters\n.forceField:PARSE\n";
	proteinDielectric=new QLabel;
	proteinDielectricInput=new QLineEdit;
	dx=new QLabel;
	dy=new QLabel;
	dz=new QLabel;
	atompot=new QCheckBox(tr("atompot"));
	charge=new QCheckBox(tr("charge"));
	dielx=new QCheckBox(tr("dielx"));
	diely=new QCheckBox(tr("diely"));
	dielz=new QCheckBox(tr("dielz"));
	edens=new QCheckBox(tr("edens"));
	ivdw=new QCheckBox(tr("ivdw"));
	kappa=new QCheckBox(tr("kappa"));
	lap=new QCheckBox(tr("lap"));
	ndens=new QCheckBox(tr("ndens"));
	pot=new QCheckBox(tr("pot"));
	qdens=new QCheckBox(tr("qdens"));
	smol=new QCheckBox(tr("smol"));
	sspl=new QCheckBox(tr("sspl"));
	vdw=new QCheckBox(tr("vdw"));
	if(!fIn.fail())
		{// Define Previous Download Folders
		defineParameters(optRefFile);
		}
	else
		{proteinDielectric->setText("4.0");
		dx->setText("161");dy->setText("161");dz->setText("161");
		atompot->setChecked(false);
		charge->setChecked(false);
		dielx->setChecked(false);
		diely->setChecked(false);
		dielz->setChecked(false);
		edens->setChecked(false);
		ivdw->setChecked(false);
		kappa->setChecked(false);
		lap->setChecked(false);
		ndens->setChecked(false);
		pot->setChecked(true);
		qdens->setChecked(false);
		smol->setChecked(false);
		sspl->setChecked(false);
		vdw->setChecked(false);
		// Write New Ref
		fOut.open(optRefFile.c_str(),ofstream::out|ofstream::trunc);fOut<<theTxt;fOut.close();
		}
	proteinDielectricInput->setText(tr(NO_VALUE));proteinDielectricInput->setAlignment(Qt::AlignCenter);
	connect(proteinDielectricInput,SIGNAL(returnPressed()),this,SLOT(proteinDielectricDefined()));
	
	QComboBox *dxB=new QComboBox;
	QComboBox *dyB=new QComboBox;
	QComboBox *dzB=new QComboBox;
	string dim[4]={"65","97","129","161"};
	string tmp;
	for(int i=0;i<4;i++)
		{tmp=dim[i];
		dxB->addItem(tr(tmp.c_str()));
		dyB->addItem(tr(tmp.c_str()));
		dzB->addItem(tr(tmp.c_str()));}
	connect(dxB,SIGNAL(activated(const QString&)),this,SLOT(labelDX(const QString&)));
	connect(dyB,SIGNAL(activated(const QString&)),this,SLOT(labelDY(const QString&)));
	connect(dzB,SIGNAL(activated(const QString&)),this,SLOT(labelDZ(const QString&)));

	QLabel* L1=new QLabel;QLabel* L2=new QLabel;QLabel* L3=new QLabel;QLabel* L4=new QLabel;
	L1->setText("<b>Protein Relative Dielectric</b>:");
	L2->setText("<b>Number of X Axis Grid Points (dx)</b>:");
	L3->setText("<b>Number of Y Axis Grid Points (dy)</b>:");
	L4->setText("<b>Number of Z Axis Grid Points (dz)</b>:");

	// APBS Parameters
	QGroupBox* pBox=new QGroupBox;
	QGridLayout *lot=new QGridLayout;
	lot->setSizeConstraint(QLayout::SetFixedSize);
	lot->addWidget(L1,0,0,1,1);lot->addWidget(proteinDielectricInput,0,1,1,1);lot->addWidget(proteinDielectric,0,2,1,1);
	lot->addWidget(L2,1,0,1,1);lot->addWidget(dxB,1,1,1,1);lot->addWidget(dx,1,2,1,1);
	lot->addWidget(L3,2,0,1,1);lot->addWidget(dyB,2,1,1,1);lot->addWidget(dy,2,2,1,1);
	lot->addWidget(L4,3,0,1,1);lot->addWidget(dzB,3,1,1,1);lot->addWidget(dz,3,2,1,1);
	pBox->setLayout(lot);
	pBox->setStyleSheet("background-color:lightgrey;");

	connect(atompot,SIGNAL(stateChanged(int)),this,SLOT(atompotChecked()));
	connect(charge,SIGNAL(stateChanged(int)),this,SLOT(chargeChecked()));
	connect(dielx,SIGNAL(stateChanged(int)),this,SLOT(dielxChecked()));
	connect(diely,SIGNAL(stateChanged(int)),this,SLOT(dielyChecked()));
	connect(dielz,SIGNAL(stateChanged(int)),this,SLOT(dielzChecked()));
	connect(edens,SIGNAL(stateChanged(int)),this,SLOT(edensChecked()));
	connect(ivdw,SIGNAL(stateChanged(int)),this,SLOT(ivdwChecked()));
	connect(kappa,SIGNAL(stateChanged(int)),this,SLOT(kappaChecked()));
	connect(lap,SIGNAL(stateChanged(int)),this,SLOT(lapChecked()));
	connect(ndens,SIGNAL(stateChanged(int)),this,SLOT(ndensChecked()));
	connect(qdens,SIGNAL(stateChanged(int)),this,SLOT(qdensChecked()));
	connect(smol,SIGNAL(stateChanged(int)),this,SLOT(smolChecked()));
	connect(sspl,SIGNAL(stateChanged(int)),this,SLOT(ssplChecked()));
	connect(vdw,SIGNAL(stateChanged(int)),this,SLOT(vdwChecked()));
	// APBS Write Types Box
	QGroupBox* wtBox=new QGroupBox(tr("Output Write Types"));
	QGridLayout *lo=new QGridLayout;
	lo->setSizeConstraint(QLayout::SetFixedSize);
	lo->addWidget(atompot,0,0,1,1);lo->addWidget(charge,0,1,1,1);lo->addWidget(dielx,0,2,1,1);lo->addWidget(diely,0,3,1,1);lo->addWidget(dielz,0,4,1,1);
	lo->addWidget(edens,1,0,1,1);lo->addWidget(ivdw,1,1,1,1);lo->addWidget(kappa,1,2,1,1);lo->addWidget(lap,1,3,1,1);lo->addWidget(ndens,1,4,1,1);
	lo->addWidget(pot,2,0,1,1);lo->addWidget(qdens,2,1,1,1);lo->addWidget(smol,2,2,1,1);lo->addWidget(sspl,2,3,1,1);lo->addWidget(vdw,2,4,1,1);
	wtBox->setLayout(lo);
	wtBox->setStyleSheet("background-color:lightgrey;");

	QGridLayout *mLayout = new QGridLayout;
	mLayout->addWidget(pBox);
	mLayout->addWidget(wtBox);
	setLayout(mLayout);
	setWindowTitle(tr("Additional APBS Parameters"));
}
///////////////////////////////////////////////////
void hydroproOptions::defineParameters(string iFile)
	{ifstream fIn;
	fIn.open(iFile.c_str());
	if(fIn.fail()){cerr<<"ERROR in hydroproOptions::defineParameters!\nInput file could not be opened.\n"<<iFile<<endl;exit(EXIT_FAILURE);}
	uint Sz=15000,pos;
	char Val[Sz];
	string bld,tmp,inputName,inputVal;
	fIn.getline(Val,Sz);
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
				if(inputVal.length()!=0)
					{// Identify Parameter and Load Values
					if(inputName.compare("CALCTYPE")==0)
						{// Protein Dielectric Constant
						if(inputVal.compare("1")==0){calcType->setText("Bead per Atom (2.9 Å)");}
						else if(inputVal.compare("2")==0){calcType->setText("Bead per Residue (4.8 Å)");}
						else if(inputVal.compare("4")==0){calcType->setText("Bead per Residue (6.1 Å)");}
						else{cerr<<"ERROR in hydroproOptions::defineParameters!\nUnrecognized calcType Value ("<<inputVal<<")."<<endl;exit(EXIT_FAILURE);}
						}
					//else{cerr<<"ERROR in hydroproOptions::defineParameters!\nUnrecognized Parameter ("<<inputName<<")."<<endl;exit(EXIT_FAILURE);}
					}
				}
			}
		fIn.getline(Val,Sz);}
	fIn.close();
	}

void hydroproOptions::labelCalcType(const QString qS)
	{string s=qS.toStdString();
	calcType->setText(s.c_str());
	if(s.compare("Bead per Atom (2.9 Å)")==0){updateParameterFile("calcType","1");}
	else if(s.compare("Bead per Residue (4.8 Å)")==0){updateParameterFile("calcType","2");}
	else if(s.compare("Bead per Residue (6.1 Å)")==0){updateParameterFile("calcType","4");}
	else{cerr<<"Error in hydroproOptions::labelCalcType!\nUnrecognized calcType Label ("<<s<<")\n";exit(EXIT_FAILURE);}
	}

hydroproOptions::hydroproOptions(QWidget *parent) : QWidget(parent)
{	
	QString Fldr=QDir::currentPath();
	string sFldr=Fldr.toStdString();
	sFldr=checkFinalBackSlash(sFldr);
	// File Reference for When Optional Parameters are set (makes .optionalParameters.txt)
	string optRefFile=getBaseFolder(sFldr)+".optionalParameters.txt";
	ifstream fIn;
	ofstream fOut;
	fIn.open(optRefFile.c_str());
	string theTxt="// APBS Parameters\n.pdie:4.0\n.dx:161\n.dy:161\n.dz:161\n.atompot:false\n.charge:false\n.dielx:false\n.diely:false\n.dielz:false\n.edens:false\n.ivdw:false\n.kappa:false\n.lap:false\n.ndens:false\n.pot:true\n.qdens:false\n.smol:false\n.sspl:false\n.vdw:false\n";
	theTxt+="// HYDROPRO Parameters\n.calcType:1\n";
	theTxt+="// MSMS Parameters\n.lowPntDensity:1.00\n.highPntDensity:3.00\n";
	theTxt+="// PDB2PQR Parameters\n.forceField:PARSE\n";
	calcType=new QLabel;
	if(!fIn.fail())
		{// Define Previous Download Folders
		defineParameters(optRefFile);
		}
	else
		{calcType->setText("Bead per Atom (2.9 Å)");
		// Write New Ref
		fOut.open(optRefFile.c_str(),ofstream::out|ofstream::trunc);fOut<<theTxt;fOut.close();
		}
	//Type of Calculation (1 = bead per atom (2.9A), 2 and 4 = bead per residue (4.8A and 6.1A))

	QComboBox* cB=new QComboBox;
	string dim[3]={"Bead per Atom (2.9 Å)","Bead per Residue (4.8 Å)","Bead per Residue (6.1 Å)"};
	string tmp;
	for(int i=0;i<3;i++)
		{tmp=dim[i];
		cB->addItem(tr(tmp.c_str()));}
	connect(cB,SIGNAL(activated(const QString&)),this,SLOT(labelCalcType(const QString&)));

	QLabel* L1=new QLabel;L1->setText("<b>Calculation Type</b>:");

	// APBS Parameters
	QGroupBox* pBox=new QGroupBox;
	QGridLayout *lot=new QGridLayout;
	lot->setSizeConstraint(QLayout::SetFixedSize);
	lot->addWidget(L1,0,0,1,1);lot->addWidget(cB,0,1,1,1);lot->addWidget(calcType,0,2,1,1);
	pBox->setLayout(lot);
	pBox->setStyleSheet("background-color:lightgrey;");

	QGridLayout *mLayout = new QGridLayout;
	mLayout->addWidget(pBox);
	setLayout(mLayout);
	setWindowTitle(tr("Additional HYDROPRO Parameters"));
}
//////////////////////////////////////////////////////
void msmsOptions::defineParameters(string iFile)
	{ifstream fIn;
	fIn.open(iFile.c_str());
	if(fIn.fail()){cerr<<"ERROR in msmsOptions::defineParameters!\nInput file could not be opened.\n"<<iFile<<endl;exit(EXIT_FAILURE);}
	uint Sz=15000,pos;
	char Val[Sz];
	string bld,tmp,inputName,inputVal;
	fIn.getline(Val,Sz);
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
				if(inputVal.length()!=0)
					{// Identify Parameter and Load Values
					if(inputName.compare("LOWPNTDENSITY")==0)
						{lowPntDensity->setText(inputVal.c_str());}
					else if(inputName.compare("HIGHPNTDENSITY")==0)
						{highPntDensity->setText(inputVal.c_str());}
					//else{cerr<<"ERROR in msmsOptions::defineParameters!\nUnrecognized Parameter ("<<inputName<<")."<<endl;exit(EXIT_FAILURE);}
					}
				}
			}
		fIn.getline(Val,Sz);}
	fIn.close();
	}

void msmsOptions::lowPntDensityDefined()
	{string s=lowPntDensityInput->text().toStdString();
	lowPntDensity->setText(tr(s.c_str()));
	updateParameterFile("lowPntDensity",s);}

void msmsOptions::highPntDensityDefined()
	{string s=highPntDensityInput->text().toStdString();
	highPntDensity->setText(tr(s.c_str()));
	updateParameterFile("highPntDensity",s);}

msmsOptions::msmsOptions(QWidget *parent) : QWidget(parent)
{	QString Fldr=QDir::currentPath();
	string sFldr=Fldr.toStdString();
	sFldr=checkFinalBackSlash(sFldr);
	// File Reference for When Optional Parameters are set (makes .optionalParameters.txt)
	string optRefFile=getBaseFolder(sFldr)+".optionalParameters.txt";
	ifstream fIn;
	ofstream fOut;
	fIn.open(optRefFile.c_str());
	string theTxt="// APBS Parameters\n.pdie:4.0\n.dx:161\n.dy:161\n.dz:161\n.atompot:false\n.charge:false\n.dielx:false\n.diely:false\n.dielz:false\n.edens:false\n.ivdw:false\n.kappa:false\n.lap:false\n.ndens:false\n.pot:true\n.qdens:false\n.smol:false\n.sspl:false\n.vdw:false\n";
	theTxt+="// HYDROPRO Parameters\n.calcType:1\n";
	theTxt+="// MSMS Parameters\n.lowPntDensity:1.00\n.highPntDensity:3.00\n";
	theTxt+="// PDB2PQR Parameters\n.forceField:PARSE\n";
	lowPntDensity=new QLabel;
	lowPntDensityInput=new QLineEdit;
	highPntDensity=new QLabel;
	highPntDensityInput=new QLineEdit;
	if(!fIn.fail())
		{// Define Previous Download Folders
		defineParameters(optRefFile);
		}
	else
		{// msms Low Surface Point Density
		lowPntDensity->setText("1.00");
		// msms High Surface Point Density
		highPntDensity->setText("3.00");		
		// Write New Ref
		fOut.open(optRefFile.c_str(),ofstream::out|ofstream::trunc);fOut<<theTxt;fOut.close();
		}
	lowPntDensityInput->setText(tr(NO_VALUE));lowPntDensityInput->setAlignment(Qt::AlignCenter);
	highPntDensityInput->setText(tr(NO_VALUE));highPntDensityInput->setAlignment(Qt::AlignCenter);

	connect(lowPntDensityInput,SIGNAL(returnPressed()),this,SLOT(lowPntDensityDefined()));
	connect(highPntDensityInput,SIGNAL(returnPressed()),this,SLOT(highPntDensityDefined()));

	QLabel* L1=new QLabel;QLabel* L2=new QLabel;
	L1->setText("<b>Low Surface Point Density</b>:");
	L2->setText("<b>High Surface Point Density</b>:");

	QGroupBox* pBox=new QGroupBox;
	QGridLayout *lot=new QGridLayout;
	lot->setSizeConstraint(QLayout::SetFixedSize);
	lot->addWidget(L1,0,0,1,1);lot->addWidget(lowPntDensityInput,0,1,1,1);lot->addWidget(lowPntDensity,0,2,1,1);
	lot->addWidget(L2,1,0,1,1);lot->addWidget(highPntDensityInput,1,1,1,1);lot->addWidget(highPntDensity,1,2,1,1);
	pBox->setLayout(lot);
	pBox->setStyleSheet("background-color:lightgrey;");

	QGridLayout *mLayout = new QGridLayout;
	mLayout->addWidget(pBox);
	setLayout(mLayout);
	setWindowTitle(tr("Additional MSMS Parameters"));
}
///////////////////////////////////////////////////
void pdb2pqrOptions::defineParameters(string iFile)
	{ifstream fIn;
	fIn.open(iFile.c_str());
	if(fIn.fail()){cerr<<"ERROR in pdb2pqrOptions::defineParameters!\nInput file could not be opened.\n"<<iFile<<endl;exit(EXIT_FAILURE);}
	uint Sz=15000,pos;
	char Val[Sz];
	string bld,tmp,inputName,inputVal;
	fIn.getline(Val,Sz);
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
				if(inputVal.length()!=0)
					{// Identify Parameter and Load Values
					if(inputName.compare("FORCEFIELD")==0)
						{forceField->setText(inputVal.c_str());}
					//else{cerr<<"ERROR in pdb2pqrOptions::defineParameters!\nUnrecognized Parameter ("<<inputName<<")."<<endl;exit(EXIT_FAILURE);}
					}
				}
			}
		fIn.getline(Val,Sz);}
	fIn.close();
	}

void pdb2pqrOptions::labelForceField(const QString qS){string s=qS.toStdString();forceField->setText(s.c_str());updateParameterFile("forceField",s);}

pdb2pqrOptions::pdb2pqrOptions(QWidget *parent) : QWidget(parent)
{	
	QString Fldr=QDir::currentPath();
	string sFldr=Fldr.toStdString();
	sFldr=checkFinalBackSlash(sFldr);
	// File Reference for When Optional Parameters are set (makes .optionalParameters.txt)
	string optRefFile=getBaseFolder(sFldr)+".optionalParameters.txt";
	ifstream fIn;
	ofstream fOut;
	fIn.open(optRefFile.c_str());
	string theTxt="// APBS Parameters\n.pdie:4.0\n.dx:161\n.dy:161\n.dz:161\n.atompot:false\n.charge:false\n.dielx:false\n.diely:false\n.dielz:false\n.edens:false\n.ivdw:false\n.kappa:false\n.lap:false\n.ndens:false\n.pot:true\n.qdens:false\n.smol:false\n.sspl:false\n.vdw:false\n";
	theTxt+="// HYDROPRO Parameters\n.calcType:1\n";
	theTxt+="// MSMS Parameters\n.lowPntDensity:1.00\n.highPntDensity:3.00\n";
	theTxt+="// PDB2PQR Parameters\n.forceField:PARSE\n";
	forceField=new QLabel;
	if(!fIn.fail())
		{// Define Previous Download Folders
		defineParameters(optRefFile);
		}
	else
		{forceField->setText("PARSE");
		// Write New Ref
		fOut.open(optRefFile.c_str(),ofstream::out|ofstream::trunc);fOut<<theTxt;fOut.close();
		}
	//PROPKA Input: Force Field (AMBER, CHARMM, PARSE, PEOEPB, SWANSON, TYL06)

	QComboBox* cB=new QComboBox;
	string dim[6]={"AMBER","CHARMM","PARSE","PEOEPB","SWANSON","TYL06"};
	string tmp;
	for(int i=0;i<6;i++)
		{tmp=dim[i];
		cB->addItem(tr(tmp.c_str()));}
	connect(cB,SIGNAL(activated(const QString&)),this,SLOT(labelForceField(const QString&)));

	QLabel* L1=new QLabel;L1->setText("<b>PROPKA Force Field</b>:");

	// APBS Parameters
	QGroupBox* pBox=new QGroupBox;
	QGridLayout *lot=new QGridLayout;
	lot->setSizeConstraint(QLayout::SetFixedSize);
	lot->addWidget(L1,0,0,1,1);lot->addWidget(cB,0,1,1,1);lot->addWidget(forceField,0,2,1,1);
	pBox->setLayout(lot);
	pBox->setStyleSheet("background-color:lightgrey;");

	QGridLayout *mLayout = new QGridLayout;
	mLayout->addWidget(pBox);
	setLayout(mLayout);
	setWindowTitle(tr("Additional PDB2PQR Parameters"));
}

void updateParameterFile(string Name,string Value)
	{// Name refers to
	Name=makeUpperCase(Name);
	QString Fldr=QDir::currentPath();
	string sFldr=Fldr.toStdString();
	sFldr=checkFinalBackSlash(sFldr);
	// File Reference for When Optional Parameters are set (makes .optionalParameters.txt)
	string optRefFile=getBaseFolder(sFldr)+".optionalParameters.txt";
	ifstream fIn;
	ofstream fOut;
	fIn.open(optRefFile.c_str());
	uint Sz=15000,pos;
	char Val[Sz];
	string cpy="",tmp,inputName,inputVal;
	if(fIn.fail())
		{cout<<"Failed to open .optionalParameters.txt\n";
		cpy="// APBS Parameters\n.pdie:4.0\n.dx:161\n.dy:161\n.dz:161\n.atompot:false\n.charge:false\n.dielx:false\n.diely:false\n.dielz:false\n.edens:false\n.ivdw:false\n.kappa:false\n.lap:false\n.ndens:false\n.pot:true\n.qdens:false\n.smol:false\n.sspl:false\n.vdw:false\n";
		cpy+="// HYDROPRO Parameters\n.calcType:1\n";
		cpy+="// MSMS Parameters\n.lowPntDensity:1.00\n.highPntDensity:3.00\n";
		cpy+="// PDB2PQR Parameters\n.forceField:PARSE\n";
		fOut.open(optRefFile.c_str(),ofstream::out|ofstream::trunc);
		if(fOut.fail()){cerr<<"ERROR in updateParameterFile!!!\nFile could not be opened.\n"<<optRefFile<<endl;exit(EXIT_FAILURE);}
		fOut<<cpy;fOut.close();
		// Reattempt
		fIn.open(optRefFile.c_str());
		// Copy & Update File
		fIn.getline(Val,Sz);
		cpy="";
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
					// Identify Parameter and Load Values
					if(inputName.compare("PDIE")==0)
						{// Protein Dielectric
						if(Name.compare("PDIE")==0){cpy+=".pdie:"+Value+"\n";}
						else{cpy+=".pdie:"+inputVal+"\n";}
						}
					else if(inputName.compare("DX")==0)
						{// APBS Num X Grid Points
						if(Name.compare("DX")==0){cpy+=".dx:"+Value+"\n";}
						else{cpy+=".dx:"+inputVal+"\n";}
						}
					else if(inputName.compare("DY")==0)
						{// APBS Num y Grid Points
						if(Name.compare("DY")==0){cpy+=".dy:"+Value+"\n";}
						else{cpy+=".dy:"+inputVal+"\n";}
						}
					else if(inputName.compare("DZ")==0)
						{// APBS Num z Grid Points
						if(Name.compare("DZ")==0){cpy+=".dz:"+Value+"\n";}
						else{cpy+=".dz:"+inputVal+"\n";}
						}
					else if(inputName.compare("ATOMPOT")==0)
						{// Write Type
						if(Name.compare("ATOMPOT")==0){cpy+=".atompot:"+Value+"\n";}
						else{cpy+=".atompot:"+inputVal+"\n";}
						}
					else if(inputName.compare("CHARGE")==0)
						{// Write Type
						if(Name.compare("CHARGE")==0){cpy+=".charge:"+Value+"\n";}
						else{cpy+=".charge:"+inputVal+"\n";}
						}
					else if(inputName.compare("DIELX")==0)
						{// Write Type
						if(Name.compare("DIELX")==0){cpy+=".dielx:"+Value+"\n";}
						else{cpy+=".dielx:"+inputVal+"\n";}
						}
					else if(inputName.compare("DIELY")==0)
						{// Write Type
						if(Name.compare("DIELY")==0){cpy+=".diely:"+Value+"\n";}
						else{cpy+=".diely:"+inputVal+"\n";}
						}
					else if(inputName.compare("DIELZ")==0)
						{// Write Type
						if(Name.compare("DIELZ")==0){cpy+=".dielz:"+Value+"\n";}
						else{cpy+=".dielz:"+inputVal+"\n";}
						}
					else if(inputName.compare("EDENS")==0)
						{// Write Type
						if(Name.compare("EDENS")==0){cpy+=".edens:"+Value+"\n";}
						else{cpy+=".edens:"+inputVal+"\n";}
						}
					else if(inputName.compare("IVDW")==0)
						{// Write Type
						if(Name.compare("IVDW")==0){cpy+=".ivdw:"+Value+"\n";}
						else{cpy+=".ivdw:"+inputVal+"\n";}
						}
					else if(inputName.compare("KAPPA")==0)
						{// Write Type
						if(Name.compare("KAPPA")==0){cpy+=".kappa:"+Value+"\n";}
						else{cpy+=".kappa:"+inputVal+"\n";}
						}
					else if(inputName.compare("LAP")==0)
						{// Write Type
						if(Name.compare("LAP")==0){cpy+=".lap:"+Value+"\n";}
						else{cpy+=".lap:"+inputVal+"\n";}
						}
					else if(inputName.compare("NDENS")==0)
						{// Write Type
						if(Name.compare("NDENS")==0){cpy+=".ndens:"+Value+"\n";}
						else{cpy+=".ndens:"+inputVal+"\n";}
						}
					else if(inputName.compare("POT")==0)
						{// Write Type
						if(Name.compare("POT")==0){cpy+=".pot:"+Value+"\n";}
						else{cpy+=".pot:"+inputVal+"\n";}
						}
					else if(inputName.compare("QDENS")==0)
						{// Write Type
						if(Name.compare("QDENS")==0){cpy+=".qdens:"+Value+"\n";}
						else{cpy+=".qdens:"+inputVal+"\n";}
						}
					else if(inputName.compare("SMOL")==0)
						{// Write Type
						if(Name.compare("SMOL")==0){cpy+=".smol:"+Value+"\n";}
						else{cpy+=".smol:"+inputVal+"\n";}
						}
					else if(inputName.compare("SSPL")==0)
						{// Write Type
						if(Name.compare("SSPL")==0){cpy+=".sspl:"+Value+"\n";}
						else{cpy+=".sspl:"+inputVal+"\n";}
						}
					else if(inputName.compare("VDW")==0)
						{// Write Type
						if(Name.compare("VDW")==0){cpy+=".vdw:"+Value+"\n";}
						else{cpy+=".vdw:"+inputVal+"\n";}
						}
					else if(inputName.compare("CALCTYPE")==0)
						{// HYDROPRO Calculation Type
						if(Name.compare("CALCTYPE")==0){cpy+=".calcType:"+Value+"\n";}
						else{cpy+=".calcType:"+inputVal+"\n";}
						}
					else if(inputName.compare("LOWPNTDENSITY")==0)
						{// MSMS Low Point Surface Density
						if(Name.compare("LOWPNTDENSITY")==0){cpy+=".lowPntDensity:"+Value+"\n";}
						else{cpy+=".lowPntDensity:"+inputVal+"\n";}
						}
					else if(inputName.compare("HIGHPNTDENSITY")==0)
						{// MSMS HI Point Surface Density
						if(Name.compare("HIGHPNTDENSITY")==0){cpy+=".highPntDensity:"+Value+"\n";}
						else{cpy+=".highPntDensity:"+inputVal+"\n";}
						}
					else if(inputName.compare("FORCEFIELD")==0)
						{// PROPKA Force Field
						if(Name.compare("FORCEFIELD")==0){cpy+=".forceField:"+Value+"\n";}
						else{cpy+=".forceField:"+inputVal+"\n";}
						}
					else{cerr<<"ERROR in updateParameterFile!\nUnrecognized Parameter ("<<inputName<<")."<<endl;exit(EXIT_FAILURE);}
					}
				}
			fIn.getline(Val,Sz);}
		fIn.close();
		}
	else
		{// Copy & Update File
		//cout<<"Opened .history.txt\n";
		fIn.getline(Val,Sz);
		//cout<<Val<<endl;
		while(!fIn.eof())
			{tmp=Val;
			//cout<<tmp<<endl;
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
					//cout<<inputName<<" "<<inputVal<<endl;
					// Identify Parameter and Load Values
					if(inputName.compare("PDIE")==0)
						{// Protein Dielectric
						if(Name.compare("PDIE")==0){cpy+=".pdie:"+Value+"\n";}
						else{cpy+=".pdie:"+inputVal+"\n";}
						}
					else if(inputName.compare("DX")==0)
						{// APBS Num X Grid Points
						if(Name.compare("DX")==0){cpy+=".dx:"+Value+"\n";}
						else{cpy+=".dx:"+inputVal+"\n";}
						}
					else if(inputName.compare("DY")==0)
						{// APBS Num y Grid Points
						if(Name.compare("DY")==0){cpy+=".dy:"+Value+"\n";}
						else{cpy+=".dy:"+inputVal+"\n";}
						}
					else if(inputName.compare("DZ")==0)
						{// APBS Num z Grid Points
						if(Name.compare("DZ")==0){cpy+=".dz:"+Value+"\n";}
						else{cpy+=".dz:"+inputVal+"\n";}
						}
					else if(inputName.compare("ATOMPOT")==0)
						{// Write Type
						if(Name.compare("ATOMPOT")==0){cpy+=".atompot:"+Value+"\n";}
						else{cpy+=".atompot:"+inputVal+"\n";}
						}
					else if(inputName.compare("CHARGE")==0)
						{// Write Type
						if(Name.compare("CHARGE")==0){cpy+=".charge:"+Value+"\n";}
						else{cpy+=".charge:"+inputVal+"\n";}
						}
					else if(inputName.compare("DIELX")==0)
						{// Write Type
						if(Name.compare("DIELX")==0){cpy+=".dielx:"+Value+"\n";}
						else{cpy+=".dielx:"+inputVal+"\n";}
						}
					else if(inputName.compare("DIELY")==0)
						{// Write Type
						if(Name.compare("DIELY")==0){cpy+=".diely:"+Value+"\n";}
						else{cpy+=".diely:"+inputVal+"\n";}
						}
					else if(inputName.compare("DIELZ")==0)
						{// Write Type
						if(Name.compare("DIELZ")==0){cpy+=".dielz:"+Value+"\n";}
						else{cpy+=".dielz:"+inputVal+"\n";}
						}
					else if(inputName.compare("EDENS")==0)
						{// Write Type
						if(Name.compare("EDENS")==0){cpy+=".edens:"+Value+"\n";}
						else{cpy+=".edens:"+inputVal+"\n";}
						}
					else if(inputName.compare("IVDW")==0)
						{// Write Type
						if(Name.compare("IVDW")==0){cpy+=".ivdw:"+Value+"\n";}
						else{cpy+=".ivdw:"+inputVal+"\n";}
						}
					else if(inputName.compare("KAPPA")==0)
						{// Write Type
						if(Name.compare("KAPPA")==0){cpy+=".kappa:"+Value+"\n";}
						else{cpy+=".kappa:"+inputVal+"\n";}
						}
					else if(inputName.compare("LAP")==0)
						{// Write Type
						if(Name.compare("LAP")==0){cpy+=".lap:"+Value+"\n";}
						else{cpy+=".lap:"+inputVal+"\n";}
						}
					else if(inputName.compare("NDENS")==0)
						{// Write Type
						if(Name.compare("NDENS")==0){cpy+=".ndens:"+Value+"\n";}
						else{cpy+=".ndens:"+inputVal+"\n";}
						}
					else if(inputName.compare("POT")==0)
						{// Write Type
						if(Name.compare("POT")==0){cpy+=".pot:"+Value+"\n";}
						else{cpy+=".pot:"+inputVal+"\n";}
						}
					else if(inputName.compare("QDENS")==0)
						{// Write Type
						if(Name.compare("QDENS")==0){cpy+=".qdens:"+Value+"\n";}
						else{cpy+=".qdens:"+inputVal+"\n";}
						}
					else if(inputName.compare("SMOL")==0)
						{// Write Type
						if(Name.compare("SMOL")==0){cpy+=".smol:"+Value+"\n";}
						else{cpy+=".smol:"+inputVal+"\n";}
						}
					else if(inputName.compare("SSPL")==0)
						{// Write Type
						if(Name.compare("SSPL")==0){cpy+=".sspl:"+Value+"\n";}
						else{cpy+=".sspl:"+inputVal+"\n";}
						}
					else if(inputName.compare("VDW")==0)
						{// Write Type
						if(Name.compare("VDW")==0){cpy+=".vdw:"+Value+"\n";}
						else{cpy+=".vdw:"+inputVal+"\n";}
						}
					else if(inputName.compare("CALCTYPE")==0)
						{// HYDROPRO Calculation Type
						if(Name.compare("CALCTYPE")==0){cpy+=".calcType:"+Value+"\n";}
						else{cpy+=".calcType:"+inputVal+"\n";}
						}
					else if(inputName.compare("LOWPNTDENSITY")==0)
						{// MSMS Low Point Surface Density
						if(Name.compare("LOWPNTDENSITY")==0){cpy+=".lowPntDensity:"+Value+"\n";}
						else{cpy+=".lowPntDensity:"+inputVal+"\n";}
						}
					else if(inputName.compare("HIGHPNTDENSITY")==0)
						{// MSMS HI Point Surface Density
						if(Name.compare("HIGHPNTDENSITY")==0){cpy+=".highPntDensity:"+Value+"\n";}
						else{cpy+=".highPntDensity:"+inputVal+"\n";}
						}
					else if(inputName.compare("FORCEFIELD")==0)
						{// PROPKA Force Field
						if(Name.compare("FORCEFIELD")==0){cpy+=".forceField:"+Value+"\n";}
						else{cpy+=".forceField:"+inputVal+"\n";}
						}
					else{cerr<<"ERROR in updateParameterFile!\nUnrecognized Parameter ("<<inputName<<")."<<endl;exit(EXIT_FAILURE);}
					}
				}
			fIn.getline(Val,Sz);}
		fIn.close();
		}
	//cout<<cpy<<endl;
	fOut.open(optRefFile.c_str(),ofstream::out|ofstream::trunc);
	if(fOut.fail()){cerr<<"ERROR in updateParameterFile!!!\nFile could not be opened.\n"<<optRefFile<<endl;exit(EXIT_FAILURE);}
	fOut<<cpy;fOut.close();
	//cout<<"Output complete.\n";
	}

double interpolate(double x,double x1,double x2,double y1,double y2){double output=(y2-y1)*(x-x1)/(x2-x1) + y1;return output;}

double getSoluteAqueousSolubility(string solute,double T)
	{double Output,MW;
	// Number of Experimental Solubilities at Different Temperature
	int N,index;
	// Convert Absolute temperature to C
	T-=273.15;
	// Water Empirical Constants for density as function of temperature
	double Aw=-2.8054253e-10,Bw=1.0556302e-7,Cw=4.6170461e-5,Dw=0.0079870401,Ew=16.945176,Fw=999.83952,Gw=0.01687985;
	// Pure Water Density (kg/m^3) (or g/L)
	double waterD=(((((Aw*T+Bw)*T-Cw)*T-Dw)*T+Ew)*T+Fw)/(1+Gw*T);
	bool RUNNING;
	if(solute.compare("AlCl3")==0)
		{// Aluminum Chloride
		MW=133.341;
		N=8;
		double s_AlCl3[N]={439,449,458,466,473,481,486,490};	// g/L
		double t_AlCl3[N]={0,10,20,30,40,60,80,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_AlCl3[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_AlCl3[index],t_AlCl3[index+1],s_AlCl3[index],s_AlCl3[index+1]);}
		else{Output=interpolate(T,t_AlCl3[index-1],t_AlCl3[index],s_AlCl3[index-1],s_AlCl3[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("Al2(SO4)3")==0)
		{// Aluminum Sulfate
		MW=342.15;
		N=11;
		double s_Al2SO43[N]={312,335,364,404,458,522,592,662,730,808,890};	// g/L
		double t_Al2SO43[N]={0,10,20,30,40,50,60,70,80,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_Al2SO43[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_Al2SO43[index],t_Al2SO43[index+1],s_Al2SO43[index],s_Al2SO43[index+1]);}
		else{Output=interpolate(T,t_Al2SO43[index-1],t_Al2SO43[index],s_Al2SO43[index-1],s_Al2SO43[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("BaCl2")==0)
		{// Barium Chloride
		MW=208.23;
		N=9;
		double s_BaCl2[N]={312,335,358,381,408,462,525,558,594};	// g/L
		double t_BaCl2[N]={0,10,20,30,40,60,80,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_BaCl2[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_BaCl2[index],t_BaCl2[index+1],s_BaCl2[index],s_BaCl2[index+1]);}
		else{Output=interpolate(T,t_BaCl2[index-1],t_BaCl2[index],s_BaCl2[index-1],s_BaCl2[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("CaCl2")==0)
		{// Calcium Chloride
		MW=110.98;
		N=9;
		double s_CaCl2[N]={595,647,745,1000,1280,1370,1470,1540,1590};	// g/L
		double t_CaCl2[N]={0,10,20,30,40,60,80,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_CaCl2[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_CaCl2[index],t_CaCl2[index+1],s_CaCl2[index],s_CaCl2[index+1]);}
		else{Output=interpolate(T,t_CaCl2[index-1],t_CaCl2[index],s_CaCl2[index-1],s_CaCl2[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("CdCl2")==0)
		{// Cadmium Chloride
		MW=183.31;
		N=8;
		double s_CdCl2[N]={1000,1350,1350,1350,1350,1360,1400,1470};	// g/L
		double t_CdCl2[N]={0,10,20,30,40,60,80,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_CdCl2[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_CdCl2[index],t_CdCl2[index+1],s_CdCl2[index],s_CdCl2[index+1]);}
		else{Output=interpolate(T,t_CdCl2[index-1],t_CdCl2[index],s_CdCl2[index-1],s_CdCl2[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("CdSO4")==0)
		{// Cadmium Sulfate
		MW=208.47;
		N=8;
		double s_CdSO4[N]={754,760,766,785,818,667,631,608};	// g/L
		double t_CdSO4[N]={0,10,20,40,60,80,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_CdSO4[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_CdSO4[index],t_CdSO4[index+1],s_CdSO4[index],s_CdSO4[index+1]);}
		else{Output=interpolate(T,t_CdSO4[index-1],t_CdSO4[index],s_CdSO4[index-1],s_CdSO4[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("CoCl2")==0)
		{// Cobalt (II) Chloride
		MW=129.839;
		N=9;
		double s_CoCl2[N]={435,477,529,597,695,938,976,1010,1060};	// g/L
		double t_CoCl2[N]={0,10,20,30,40,60,80,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_CoCl2[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_CoCl2[index],t_CoCl2[index+1],s_CoCl2[index],s_CoCl2[index+1]);}
		else{Output=interpolate(T,t_CoCl2[index-1],t_CoCl2[index],s_CoCl2[index-1],s_CoCl2[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("CoSO4")==0)
		{// Cobalt (II) Sulfate
		MW=154.996;
		N=9;
		double s_CoSO4[N]={255,305,361,420,488,550,538,453,389};	// g/L
		double t_CoSO4[N]={0,10,20,30,40,60,80,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_CoSO4[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_CoSO4[index],t_CoSO4[index+1],s_CoSO4[index],s_CoSO4[index+1]);}
		else{Output=interpolate(T,t_CoSO4[index-1],t_CoSO4[index],s_CoSO4[index-1],s_CoSO4[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
/*
	else if(solute.compare("CrCl3")==0)
		{// Chormium (III) Chloride
		MW=;
		N=;
		double s_CrCl3[N]={};	// g/L
		double t_CrCl3[N]={0,10,20,30,40,60,80,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_CrCl3[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_CrCl3[index],t_CrCl3[index+1],s_CrCl3[index],s_CrCl3[index+1]);}
		else{Output=interpolate(T,t_CrCl3[index-1],t_CrCl3[index],s_CrCl3[index-1],s_CrCl3[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("Cr2(SO4)3")==0)
		{// Chormium (III) Sulfate
		MW=392.16;
		N=;
		double s_Cr2SO43[N]={};	// g/L
		double t_Cr2SO43[N]={0,10,20,30,40,60,80,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_Cr2SO43[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_Cr2SO43[index],t_Cr2SO43[index+1],s_Cr2SO43[index],s_Cr2SO43[index+1]);}
		else{Output=interpolate(T,t_Cr2SO43[index-1],t_Cr2SO43[index],s_Cr2SO43[index-1],s_Cr2SO43[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
*/
	else if(solute.compare("CuCl2")==0)
		{// Copper (II) Chloride
		MW=134.45;
		N=9;
		double s_CuCl2[N]={686,709,730,773,876,965,1040,1080,1200};	// g/L
		double t_CuCl2[N]={0,10,20,30,40,60,80,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_CuCl2[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_CuCl2[index],t_CuCl2[index+1],s_CuCl2[index],s_CuCl2[index+1]);}
		else{Output=interpolate(T,t_CuCl2[index-1],t_CuCl2[index],s_CuCl2[index-1],s_CuCl2[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("CuSO4")==0)
		{// Copper (II) Sulfate
		MW=159.609;
		N=8;
		double s_CuSO4[N]={231,275,320,378,446,618,838,1140};	// g/L
		double t_CuSO4[N]={0,10,20,30,40,60,80,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_CuSO4[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_CuSO4[index],t_CuSO4[index+1],s_CuSO4[index],s_CuSO4[index+1]);}
		else{Output=interpolate(T,t_CuSO4[index-1],t_CuSO4[index],s_CuSO4[index-1],s_CuSO4[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("FeCl2")==0)
		{// Iron (II) Chloride
		MW=126.751;
		N=9;
		double s_FeCl2[N]={497,590,625,667,700,783,887,923,949};	// g/L
		double t_FeCl2[N]={0,10,20,30,40,60,80,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_FeCl2[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_FeCl2[index],t_FeCl2[index+1],s_FeCl2[index],s_FeCl2[index+1]);}
		else{Output=interpolate(T,t_FeCl2[index-1],t_FeCl2[index],s_FeCl2[index-1],s_FeCl2[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("FeSO4")==0)
		{// Iron (II) Sulfate
		MW=151.91;
		N=7;
		double s_FeSO4[N]={288,400,480,600,733,1010,799};	// g/L
		double t_FeSO4[N]={20,40,50,60,70,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_FeSO4[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_FeSO4[index],t_FeSO4[index+1],s_FeSO4[index],s_FeSO4[index+1]);}
		else{Output=interpolate(T,t_FeSO4[index-1],t_FeSO4[index],s_FeSO4[index-1],s_FeSO4[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("FeCl3")==0)
		{// Iron (III) Chloride
		MW=162.204;
		N=3;
		double s_FeCl3[N]={744,918,1070};	// g/L
		double t_FeCl3[N]={0,20,30};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_FeCl3[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_FeCl3[index],t_FeCl3[index+1],s_FeCl3[index],s_FeCl3[index+1]);}
		else{Output=interpolate(T,t_FeCl3[index-1],t_FeCl3[index],s_FeCl3[index-1],s_FeCl3[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
/*
	else if(solute.compare("Fe2(SO4)3")==0)
		{// Iron (III) Sulfate
		MW=399.88;
		N=;
		double s_Fe2SO43[N]={};	// g/L
		double t_Fe2SO43[N]={0,20,30};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_Fe2SO43[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_Fe2SO43[index],t_Fe2SO43[index+1],s_Fe2SO43[index],s_Fe2SO43[index+1]);}
		else{Output=interpolate(T,t_Fe2SO43[index-1],t_Fe2SO43[index],s_Fe2SO43[index-1],s_Fe2SO43[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
*/
	else if(solute.compare("HCl")==0)
		{// Hydrochloric Acid
		MW=36.46;
		N=11;
		double s_HCl[N]={810,750,700,655,610,575,530,500,470,430,400};	// g/L
		double t_HCl[N]={0,10,20,30,40,50,60,70,80,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_HCl[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_HCl[index],t_HCl[index+1],s_HCl[index],s_HCl[index+1]);}
		else{Output=interpolate(T,t_HCl[index-1],t_HCl[index],s_HCl[index-1],s_HCl[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
/*
	else if(solute.compare("HCN")==0)
		{// Hydrogen Cyanide
		MW=26.018;
		N=;
		double s_HCN[N]={};	// g/L
		double t_HCN[N]={0,10,20,30,40,50,60,70,80,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_HCN[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_HCN[index],t_HCN[index+1],s_HCN[index],s_HCN[index+1]);}
		else{Output=interpolate(T,t_HCN[index-1],t_HCN[index],s_HCN[index-1],s_HCN[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
*/
	else if(solute.compare("HNO3")==0)
		{// Nitric Acid
		MW=63.012;
		N=2;
		double s_HNO3[N]={0.95/(0.05/waterD),0.95/(0.05/waterD)};	// g/L
		double t_HNO3[N]={0,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_HNO3[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_HNO3[index],t_HNO3[index+1],s_HNO3[index],s_HNO3[index+1]);}
		else{Output=interpolate(T,t_HNO3[index-1],t_HNO3[index],s_HNO3[index-1],s_HNO3[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("H3PO4")==0)
		{// Phosphorous Acid
		MW=97.994;
		N=3;
		double s_H3PO4[N]={3694,4460,5480};	// g/L
		double t_H3PO4[N]={0.5,15,20};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_H3PO4[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_H3PO4[index],t_H3PO4[index+1],s_H3PO4[index],s_H3PO4[index+1]);}
		else{Output=interpolate(T,t_H3PO4[index-1],t_H3PO4[index],s_H3PO4[index-1],s_H3PO4[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("H2SO4")==0)
		{// Sulfric Acid (quite soluble in H2O)
		MW=98.079;
		N=2;
		double s_H2SO4[N]={18.4*MW,18.4*MW};	// g/L
		double t_H2SO4[N]={0,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_H2SO4[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_H2SO4[index],t_H2SO4[index+1],s_H2SO4[index],s_H2SO4[index+1]);}
		else{Output=interpolate(T,t_H2SO4[index-1],t_H2SO4[index],s_H2SO4[index-1],s_H2SO4[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("KCl")==0)
		{// Potassium Chloride
		MW=74.5513;
		N=10;
		double s_KCl[N]={280,312,342,372,401,426,458,513,539,563};	// g/L
		double t_KCl[N]={0,10,20,30,40,50,60,80,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_KCl[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_KCl[index],t_KCl[index+1],s_KCl[index],s_KCl[index+1]);}
		else{Output=interpolate(T,t_KCl[index-1],t_KCl[index],s_KCl[index-1],s_KCl[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("KClO4")==0)
		{// Potassium Perchlorate
		MW=138.55;
		N=9;
		double s_KClO4[N]={7.6,10.6,16.8,25.6,37.3,73,134,177,223};	// g/L
		double t_KClO4[N]={0,10,20,30,40,60,80,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_KClO4[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_KClO4[index],t_KClO4[index+1],s_KClO4[index],s_KClO4[index+1]);}
		else{Output=interpolate(T,t_KClO4[index-1],t_KClO4[index],s_KClO4[index-1],s_KClO4[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("K2CO3")==0)
		{// Potassium Carbonate
		MW=138.205;
		N=10;
		double s_K2CO3[N]={1050,1090,1110,1140,1170,1212,1270,1400,1480,1560};	// g/L
		double t_K2CO3[N]={0,10,20,30,40,50,60,80,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_K2CO3[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_K2CO3[index],t_K2CO3[index+1],s_K2CO3[index],s_K2CO3[index+1]);}
		else{Output=interpolate(T,t_K2CO3[index-1],t_K2CO3[index],s_K2CO3[index-1],s_K2CO3[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("KH2PO4")==0)
		{// Potassium Dihydrogen Phosphate
		MW=136.086;
		N=9;
		double s_KH2PO4[N]={148,183,226,280,355,410,502,704,835};	// g/L
		double t_KH2PO4[N]={0,10,20,30,40,50,60,80,90};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_KH2PO4[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_KH2PO4[index],t_KH2PO4[index+1],s_KH2PO4[index],s_KH2PO4[index+1]);}
		else{Output=interpolate(T,t_KH2PO4[index-1],t_KH2PO4[index],s_KH2PO4[index-1],s_KH2PO4[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("K2HPO4")==0)
		{// Potassium Hydrogen Phosphate
		MW=174.2;
		N=2;
		double s_K2HPO4[N]={1492.5,1500};	// g/L
		double t_K2HPO4[N]={10,20};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_K2HPO4[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_K2HPO4[index],t_K2HPO4[index+1],s_K2HPO4[index],s_K2HPO4[index+1]);}
		else{Output=interpolate(T,t_K2HPO4[index-1],t_K2HPO4[index],s_K2HPO4[index-1],s_K2HPO4[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("KNO3")==0)
		{// Potassium Nitrate
		MW=101.1032;
		N=11;
		double s_KNO3[N]={133,209,316,458,639,855,1100,1380,1690,2020,2460};	// g/L
		double t_KNO3[N]={0,10,20,30,40,50,60,70,80,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_KNO3[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_KNO3[index],t_KNO3[index+1],s_KNO3[index],s_KNO3[index+1]);}
		else{Output=interpolate(T,t_KNO3[index-1],t_KNO3[index],s_KNO3[index-1],s_KNO3[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("KOH")==0)
		{// Potassium Hydroxide
		MW=56.11;
		N=7;
		double s_KOH[N]={957,1030,1120,1260,1340,1540,1780};	// g/L
		double t_KOH[N]={0,10,20,30,40,60,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_KOH[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_KOH[index],t_KOH[index+1],s_KOH[index],s_KOH[index+1]);}
		else{Output=interpolate(T,t_KOH[index-1],t_KOH[index],s_KOH[index-1],s_KOH[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("K2SO4")==0)
		{// Potassium Sulfate
		MW=174.259;
		N=9;
		double s_K2SO4[N]={74,93,111,130,148,182,214,229,241};	// g/L
		double t_K2SO4[N]={0,10,20,30,40,60,80,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_K2SO4[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_K2SO4[index],t_K2SO4[index+1],s_K2SO4[index],s_K2SO4[index+1]);}
		else{Output=interpolate(T,t_K2SO4[index-1],t_K2SO4[index],s_K2SO4[index-1],s_K2SO4[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("LiCl")==0)
		{// Lithium Chloride
		MW=42.39;
		N=9;
		double s_LiCl[N]={692,745,835,862,898,984,1120,1210,1280};	// g/L
		double t_LiCl[N]={0,10,20,30,40,60,80,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_LiCl[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_LiCl[index],t_LiCl[index+1],s_LiCl[index],s_LiCl[index+1]);}
		else{Output=interpolate(T,t_LiCl[index-1],t_LiCl[index],s_LiCl[index-1],s_LiCl[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("Li2SO4")==0)
		{// Lithium Sulfate
		MW=109.94;
		N=8;
		double s_Li2SO4[N]={361,355,348,342,337,326,314,309};	// g/L
		double t_Li2SO4[N]={0,10,20,30,40,60,80,90};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_Li2SO4[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_Li2SO4[index],t_Li2SO4[index+1],s_Li2SO4[index],s_Li2SO4[index+1]);}
		else{Output=interpolate(T,t_Li2SO4[index-1],t_Li2SO4[index],s_Li2SO4[index-1],s_Li2SO4[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("MgCl2")==0)
		{// Magnesium Chloride
		MW=95.211;
		N=9;
		double s_MgCl2[N]={529,536,546,558,575,610,661,695,733};	// g/L
		double t_MgCl2[N]={0,10,20,30,40,60,80,90};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_MgCl2[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_MgCl2[index],t_MgCl2[index+1],s_MgCl2[index],s_MgCl2[index+1]);}
		else{Output=interpolate(T,t_MgCl2[index-1],t_MgCl2[index],s_MgCl2[index-1],s_MgCl2[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("MgSO4")==0)
		{// Magnesium Sulfate
		MW=120.366;
		N=11;
		double s_MgSO4[N]={255,304,351,397,447,504,548,592,548,529,502};	// g/L
		double t_MgSO4[N]={0,10,20,30,40,50,60,70,80,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_MgSO4[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_MgSO4[index],t_MgSO4[index+1],s_MgSO4[index],s_MgSO4[index+1]);}
		else{Output=interpolate(T,t_MgSO4[index-1],t_MgSO4[index],s_MgSO4[index-1],s_MgSO4[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("MnCl2")==0)
		{// Manganese (II) Chloride
		MW=125.844;
		N=9;
		double s_MnCl2[N]={634,681,739,808,885,1090,1130,1140,1150};	// g/L
		double t_MnCl2[N]={0,10,20,30,40,60,80,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_MnCl2[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_MnCl2[index],t_MnCl2[index+1],s_MnCl2[index],s_MnCl2[index+1]);}
		else{Output=interpolate(T,t_MnCl2[index-1],t_MnCl2[index],s_MnCl2[index-1],s_MnCl2[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("MnSO4")==0)
		{// Manganese (II) Sulfate
		MW=151.001;
		N=9;
		double s_MnSO4[N]={529,597,629,629,600,536,456,409,353};	// g/L
		double t_MnSO4[N]={0,10,20,30,40,60,80,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_MnSO4[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_MnSO4[index],t_MnSO4[index+1],s_MnSO4[index],s_MnSO4[index+1]);}
		else{Output=interpolate(T,t_MnSO4[index-1],t_MnSO4[index],s_MnSO4[index-1],s_MnSO4[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("NaBr")==0)
		{// Sodium Bromide
		MW=102.894;
		N=9;
		double s_NaBr[N]={802,852,908,984,1070,1180,1200,1210,1210};	// g/L
		double t_NaBr[N]={0,10,20,30,40,60,80,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_NaBr[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_NaBr[index],t_NaBr[index+1],s_NaBr[index],s_NaBr[index+1]);}
		else{Output=interpolate(T,t_NaBr[index-1],t_NaBr[index],s_NaBr[index-1],s_NaBr[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("NaCl")==0)
		{// Sodium Chloride
		MW=58.443;
		N=11;
		double s_NaCl[N]={356.5,357.2,358.9,360.9,363.7,366.9,370.4,374.6,379.3,384.7,389.9};	// g/L
		double t_NaCl[N]={0,10,20,30,40,50,60,70,80,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_NaCl[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_NaCl[index],t_NaCl[index+1],s_NaCl[index],s_NaCl[index+1]);}
		else{Output=interpolate(T,t_NaCl[index-1],t_NaCl[index],s_NaCl[index-1],s_NaCl[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("NaClO3")==0)
		{// Sodium Chlorate
		MW=106.44;
		N=9;
		double s_NaClO3[N]={796,876,959,1050,1150,1370,1670,1840,2040};	// g/L
		double t_NaClO3[N]={0,10,20,30,40,60,80,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_NaClO3[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_NaClO3[index],t_NaClO3[index+1],s_NaClO3[index],s_NaClO3[index+1]);}
		else{Output=interpolate(T,t_NaClO3[index-1],t_NaClO3[index],s_NaClO3[index-1],s_NaClO3[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("Na2CO3")==0)
		{// Sodium Carbonate
		MW=105.9888;
		N=9;
		double s_Na2CO3[N]={70,125,215,397,490,460,439,439,455};	// g/L
		double t_Na2CO3[N]={0,10,20,30,40,60,80,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_Na2CO3[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_Na2CO3[index],t_Na2CO3[index+1],s_Na2CO3[index],s_Na2CO3[index+1]);}
		else{Output=interpolate(T,t_Na2CO3[index-1],t_Na2CO3[index],s_Na2CO3[index-1],s_Na2CO3[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("NaF")==0)
		{// Sodium Fluoride
		MW=41.9881;
		N=7;
		double s_NaF[N]={36.6,40.6,42.2,44,46.8,48.9,50.8};	// g/L
		double t_NaF[N]={0,20,30,40,60,80,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_NaF[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_NaF[index],t_NaF[index+1],s_NaF[index],s_NaF[index+1]);}
		else{Output=interpolate(T,t_NaF[index-1],t_NaF[index],s_NaF[index-1],s_NaF[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("NaHCO3")==0)
		{// Sodium Hydrogen Carbonate
		MW=84.0066;
		N=6;
		double s_NaHCO3[N]={70,81,96,111,127,160};	// g/L
		double t_NaHCO3[N]={0,10,20,30,40,60};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_NaHCO3[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_NaHCO3[index],t_NaHCO3[index+1],s_NaHCO3[index],s_NaHCO3[index+1]);}
		else{Output=interpolate(T,t_NaHCO3[index-1],t_NaHCO3[index],s_NaHCO3[index-1],s_NaHCO3[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("NaH2PO4")==0)
		{// Mono-Sodium Phosphate
		MW=119.98;
		N=8;
		double s_NaH2PO4[N]={565,698,869,1070,1330,1720,2110,2340};	// g/L
		double t_NaH2PO4[N]={0,10,20,30,40,60};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_NaH2PO4[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_NaH2PO4[index],t_NaH2PO4[index+1],s_NaH2PO4[index],s_NaH2PO4[index+1]);}
		else{Output=interpolate(T,t_NaH2PO4[index-1],t_NaH2PO4[index],s_NaH2PO4[index-1],s_NaH2PO4[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("Na2HPO4")==0)
		{// Di-Sodium Phosphate Reference: Lange's Handbook of Chemistry by John A Dean Table 5.2
		MW=141.96;
		N=9;
		double s_Na2HPO4[N]={0.0168*waterD,0.0353*waterD,0.0783*waterD,0.22*waterD,0.553*waterD,0.828*waterD,0.923*waterD,1.02*waterD,1.04*waterD};	// g Na2HPO4/g Water
		double t_Na2HPO4[N]={0,10,20,30,40,60,80,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_Na2HPO4[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_Na2HPO4[index],t_Na2HPO4[index+1],s_Na2HPO4[index],s_Na2HPO4[index+1]);}
		else{Output=interpolate(T,t_Na2HPO4[index-1],t_Na2HPO4[index],s_Na2HPO4[index-1],s_Na2HPO4[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("NaHSO3")==0)
		{// Sodium Bisulfite
		MW=104.061;
		N=2;
		double s_NaHSO3[N]={420,420};	// g/L
		double t_NaHSO3[N]={20,25};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_NaHSO3[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_NaHSO3[index],t_NaHSO3[index+1],s_NaHSO3[index],s_NaHSO3[index+1]);}
		else{Output=interpolate(T,t_NaHSO3[index-1],t_NaHSO3[index],s_NaHSO3[index-1],s_NaHSO3[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("NaI")==0)
		{// Sodium Iodide
		MW=149.89;
		N=8;
		double s_NaI[N]={1590,1670,1780,1910,2050,2570,2950,3020};	// g/L
		double t_NaI[N]={0,10,20,30,40,60};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_NaI[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_NaI[index],t_NaI[index+1],s_NaI[index],s_NaI[index+1]);}
		else{Output=interpolate(T,t_NaI[index-1],t_NaI[index],s_NaI[index-1],s_NaI[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("Na2MoO4")==0)
		{// Sodium Molybdate
		MW=205.92;
		N=6;
		double s_Na2MoO4[N]={441,647,653,669,686,718};	// g/L
		double t_Na2MoO4[N]={0,10,20,30,40,60};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_Na2MoO4[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_Na2MoO4[index],t_Na2MoO4[index+1],s_Na2MoO4[index],s_Na2MoO4[index+1]);}
		else{Output=interpolate(T,t_Na2MoO4[index-1],t_Na2MoO4[index],s_Na2MoO4[index-1],s_Na2MoO4[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("NaNO2")==0)
		{// Sodium Nitrite
		MW=68.9953;
		N=8;
		double s_NaNO2[N]={712,751,808,876,949,1110,1330,1600};	// g/L
		double t_NaNO2[N]={0,10,20,30,40,60,80,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_NaNO2[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_NaNO2[index],t_NaNO2[index+1],s_NaNO2[index],s_NaNO2[index+1]);}
		else{Output=interpolate(T,t_NaNO2[index-1],t_NaNO2[index],s_NaNO2[index-1],s_NaNO2[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("NaNO3")==0)
		{// Sodium Nitrate
		MW=84.9947;
		N=8;
		double s_NaNO3[N]={730,808,876,949,1020,1220,1480,1800};	// g/L
		double t_NaNO3[N]={0,10,20,30,40,60,80,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_NaNO3[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_NaNO3[index],t_NaNO3[index+1],s_NaNO3[index],s_NaNO3[index+1]);}
		else{Output=interpolate(T,t_NaNO3[index-1],t_NaNO3[index],s_NaNO3[index-1],s_NaNO3[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("NaOH")==0)
		{// Sodium Hydroxide
		MW=39.9971;
		N=5;
		double s_NaOH[N]={980,1090,1190,1290,1740};	// g/L
		double t_NaOH[N]={10,20,30,40,60};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_NaOH[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_NaOH[index],t_NaOH[index+1],s_NaOH[index],s_NaOH[index+1]);}
		else{Output=interpolate(T,t_NaOH[index-1],t_NaOH[index],s_NaOH[index-1],s_NaOH[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("Na3PO4")==0)
		{// Sodium Phosphate
		MW=163.94;
		N=9;
		double s_Na3PO4[N]={45,82,121,163,202,209,600,681,770};	// g/L
		double t_Na3PO4[N]={0,10,20,30,40,60,80,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_Na3PO4[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_Na3PO4[index],t_Na3PO4[index+1],s_Na3PO4[index],s_Na3PO4[index+1]);}
		else{Output=interpolate(T,t_Na3PO4[index-1],t_Na3PO4[index],s_Na3PO4[index-1],s_Na3PO4[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("Na2SO3")==0)
		{// Sodium Sulfite
		MW=126.043;
		N=2;
		double s_Na2SO3[N]={270,270};	// g/L
		double t_Na2SO3[N]={10,20};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_Na2SO3[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_Na2SO3[index],t_Na2SO3[index+1],s_Na2SO3[index],s_Na2SO3[index+1]);}
		else{Output=interpolate(T,t_Na2SO3[index-1],t_Na2SO3[index],s_Na2SO3[index-1],s_Na2SO3[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("Na2S2O3")==0)
		{// Sodium Thiosulfate
		MW=158.11;
		N=5;
		double s_Na2S2O3[N]={715,730,776,908,972};	// g/L
		double t_Na2S2O3[N]={0,20,40,80,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_Na2S2O3[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_Na2S2O3[index],t_Na2S2O3[index+1],s_Na2S2O3[index],s_Na2S2O3[index+1]);}
		else{Output=interpolate(T,t_Na2S2O3[index-1],t_Na2S2O3[index],s_Na2S2O3[index-1],s_Na2S2O3[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("Na2SO4")==0)
		{// Sodium Sulfate
		MW=142.04;
		N=9;
		double s_Na2SO4[N]={49,91,195,408,488,453,437,427,425};	// g/L
		double t_Na2SO4[N]={0,10,20,30,40,60,80,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_Na2SO4[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_Na2SO4[index],t_Na2SO4[index+1],s_Na2SO4[index],s_Na2SO4[index+1]);}
		else{Output=interpolate(T,t_Na2SO4[index-1],t_Na2SO4[index],s_Na2SO4[index-1],s_Na2SO4[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("NH3")==0)
		{// Ammonia
		MW=17.031;
		N=3;
		//double s_NH3[N]={1176,900,702,565,428,333,252,188,138,100,88};	// L/L ??
		double s_NH3[N]={0.47/(0.53/waterD),0.31/(0.69/waterD),0.18/(0.82/waterD)};	// g/L
		double t_NH3[N]={0,25,50};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_NH3[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_NH3[index],t_NH3[index+1],s_NH3[index],s_NH3[index+1]);}
		else{Output=interpolate(T,t_NH3[index-1],t_NH3[index],s_NH3[index-1],s_NH3[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("NH4Cl")==0)
		{// Ammonium Chloride
		MW=53.49;
		N=11;
		double s_NH4Cl[N]={294,332,372,414,458,504,553,602,656,712,773};	// g/L
		double t_NH4Cl[N]={0,10,20,30,40,50,60,70,80,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_NH4Cl[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_NH4Cl[index],t_NH4Cl[index+1],s_NH4Cl[index],s_NH4Cl[index+1]);}
		else{Output=interpolate(T,t_NH4Cl[index-1],t_NH4Cl[index],s_NH4Cl[index-1],s_NH4Cl[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("NH4NO3")==0)
		{// Ammonium Nitrate
		MW=80.043;
		N=11;
		double s_NH4NO3[N]={118,150,192,242,297,344,421,499,580,740,871};	// g/L
		double t_NH4NO3[N]={0,10,20,30,40,50,60,70,80,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_NH4NO3[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_NH4NO3[index],t_NH4NO3[index+1],s_NH4NO3[index],s_NH4NO3[index+1]);}
		else{Output=interpolate(T,t_NH4NO3[index-1],t_NH4NO3[index],s_NH4NO3[index-1],s_NH4NO3[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("(NH4)2SO4")==0)
		{// Ammonium Sulfate
		MW=132.14;
		N=9;
		double s_NH42SO4[N]={706,730,754,781,812,843,874,941,1030};	// g/L
		double t_NH42SO4[N]={0,10,20,30,40,50,60,80,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_NH42SO4[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_NH42SO4[index],t_NH42SO4[index+1],s_NH42SO4[index],s_NH42SO4[index+1]);}
		else{Output=interpolate(T,t_NH42SO4[index-1],t_NH42SO4[index],s_NH42SO4[index-1],s_NH42SO4[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("NiCl2")==0)
		{// Nickel (II) Chloride
		MW=129.5994;
		N=8;
		double s_NiCl2[N]={534,563,668,706,732,812,866,876};	// g/L
		double t_NiCl2[N]={0,10,20,30,40,60,80,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_NiCl2[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_NiCl2[index],t_NiCl2[index+1],s_NiCl2[index],s_NiCl2[index+1]);}
		else{Output=interpolate(T,t_NiCl2[index-1],t_NiCl2[index],s_NiCl2[index-1],s_NiCl2[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("NiSO4")==0)
		{// Nickel (II) Sulfate
		MW=154.75;
		N=7;
		double s_NiSO4[N]={444,466,492,556,645,701,767};	// g/L
		double t_NiSO4[N]={20,30,40,60,80,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_NiSO4[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_NiSO4[index],t_NiSO4[index+1],s_NiSO4[index],s_NiSO4[index+1]);}
		else{Output=interpolate(T,t_NiSO4[index-1],t_NiSO4[index],s_NiSO4[index-1],s_NiSO4[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("SrCl2")==0)
		{// Strontium Chloride
		MW=158.53;
		N=8;
		double s_SrCl2[N]={435,477,529,587,653,818,905,1010};	// g/L
		double t_SrCl2[N]={0,10,20,30,40,60,80,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_SrCl2[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_SrCl2[index],t_SrCl2[index+1],s_SrCl2[index],s_SrCl2[index+1]);}
		else{Output=interpolate(T,t_SrCl2[index-1],t_SrCl2[index],s_SrCl2[index-1],s_SrCl2[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("ZnCl2")==0)
		{// Zinc Chloride
		MW=136.315;
		N=8;
		double s_ZnCl2[N]={3420,3630,3950,4370,4520,4880,5410,6140};	// g/L
		double t_ZnCl2[N]={0,10,20,30,40,60,80,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_ZnCl2[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_ZnCl2[index],t_ZnCl2[index+1],s_ZnCl2[index],s_ZnCl2[index+1]);}
		else{Output=interpolate(T,t_ZnCl2[index-1],t_ZnCl2[index],s_ZnCl2[index-1],s_ZnCl2[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("ZnSO4")==0)
		{// Zinc Sulfate
		MW=161.47;
		N=8;
		double s_ZnSO4[N]={416,472,538,613,705,754,711,605};	// g/L
		double t_ZnSO4[N]={0,10,20,30,40,60,80,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_ZnSO4[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_ZnSO4[index],t_ZnSO4[index+1],s_ZnSO4[index],s_ZnSO4[index+1]);}
		else{Output=interpolate(T,t_ZnSO4[index-1],t_ZnSO4[index],s_ZnSO4[index-1],s_ZnSO4[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else if(solute.compare("H3Citrate")==0)
		{// Citric Acid
		MW=192.123;
		N=10;
		double s_H3Citrate[N]={540,592,643,686,709,735,762,788,814,840};	// g/L
		double t_H3Citrate[N]={10,20,30,40,50,60,70,80,90,100};
		RUNNING=true;index=0;
		while(RUNNING)
			{if(T<=t_H3Citrate[index]){RUNNING=false;}
			else{index++;}
			if(index>N){index--;RUNNING=false;}
			}
		if(index==0){Output=interpolate(T,t_H3Citrate[index],t_H3Citrate[index+1],s_H3Citrate[index],s_H3Citrate[index+1]);}
		else{Output=interpolate(T,t_H3Citrate[index-1],t_H3Citrate[index],s_H3Citrate[index-1],s_H3Citrate[index]);}
		// Convert from g/L to Molar Concentration
		Output/=MW;
		}
	else
		{cerr<<"Error in getSoluteAqueousSolubility!\nUnrecognized solute ("<<solute<<")\n";exit(EXIT_FAILURE);}

	return Output;}

methodsInfo::methodsInfo(QWidget *parent) : QWidget(parent)
{	
	int imgWidth=400;
	int imgHeight=400;
	int figWidth=800;
	int figHeight=800;
	int lHeight=100;
	int lWidth=100;
	int propertyBoxSz=1150;
	// Get Current File Path and Find Images
	QString theCurFldr=QDir::currentPath();
	string Fldr=theCurFldr.toStdString()+"/";
	// Images
	// Lysozyme Monomer Animation File
	string lyzMonoFile=Fldr+"UI/Images/6LYZ.gif";
	// Lysozyme Dimer Animation File
	string lyzDiFile=Fldr+"UI/Images/4R0F.gif";
	// Lysozyme Charge Valence vs H Titration Curve
	string lyzChargeFile=Fldr+"UI/Images/LYZ_Charge_Valence.png";
	// Lysozyme Diffusivity
	string lyzDiffFile=Fldr+"UI/Images/LYZ_Diffusivity.png";
	// Bovine Serum Albumin Monomer Animation File
	string bsaMonoFile=Fldr+"UI/Images/4F5SA.gif";
	// Bovine Serum Albumin Dimer Animation File
	string bsaDiFile=Fldr+"UI/Images/4F5S.gif";
	// BSA Diffusivity
	string bsaDiffFile=Fldr+"UI/Images/BSA_Diffusivity.png";
	// Bovine Serum Albumin Charge Valence vs H Titration Curve
	string bsaChargeFile=Fldr+"UI/Images/BSA_Charge_Valence.png";
	// Bovine B-Lactoglobulin Monomer Animation File
	string blgMonoFile=Fldr+"UI/Images/3BLG.gif";
	// Bovine B-Lactoglobulin Dimer Animation File
	string blgDiFile=Fldr+"UI/Images/1BEB.gif";
	// Bovine B-Lactoglobulin Charge Valence vs H Titration Curve
	string blgChargeFile=Fldr+"UI/Images/BLG_Charge_Valence.png";
	// BLG Diffusivity
	string blgDiffFile=Fldr+"UI/Images/BLG_Diffusivity.png";

	QMovie *s11Movie = new QMovie(lyzMonoFile.c_str());
	QMovie *s12Movie = new QMovie(lyzDiFile.c_str());
	QMovie *s21Movie = new QMovie(blgMonoFile.c_str());
	QMovie *s22Movie = new QMovie(blgDiFile.c_str());
	QMovie *s31Movie = new QMovie(bsaMonoFile.c_str());
	QMovie *s32Movie = new QMovie(bsaDiFile.c_str());
	
	QPixmap fig1(lyzChargeFile.c_str());
	QPixmap fig2(bsaChargeFile.c_str());
	QPixmap fig3(blgChargeFile.c_str());
	//
	QPixmap fig4(lyzDiffFile.c_str());
	QPixmap fig5(bsaDiffFile.c_str());
	QPixmap fig6(blgDiffFile.c_str());
	//

	string s="<b>This interface exists to show the high accuracy that can be achieved by this software. A concise description of the experimental validation performed is presented with data. ZPRED accomplishes its computation through different software components, which determine the following experimentally verifiable parameters of the input structure(s):<ol><li>charge (<b>Q</b>)</li><li>single particle diffusivity (<b>D<sub>o</sub></b>)</li><li>electrophoretic mobility (<b>u<sub>e</sub></b>)</li></ol>Experimental validation for each of these three computed parameters is presented in numerical order below the images. Each software (described below) was applied to an ensemble of structures sampled from molecular dynamics and then averaged to provide the values presented. The actual structures of both the monomeric and dimeric forms are presented below in a gif image, which has exaggerated rotation to allow for a better view of the molecules.</b>";
	QLabel* desc=new QLabel;
	desc->setText(s.c_str());
	desc->setFixedWidth(propertyBoxSz);
	desc->setWordWrap(true);
	
	QLabel* L00=new QLabel;QLabel* L01=new QLabel;QLabel* L02=new QLabel;
	QLabel* L10=new QLabel;QLabel* L11=new QLabel;QLabel* L12=new QLabel;
	QLabel* L20=new QLabel;QLabel* L21=new QLabel;QLabel* L22=new QLabel;
	QLabel* L30=new QLabel;QLabel* L31=new QLabel;QLabel* L32=new QLabel;
	L00->setFixedHeight(lHeight);
	L00->setFixedWidth(lWidth);
	L00->setAlignment(Qt::AlignCenter);
	L01->setText("<b>Monomer</b>");
	L01->setFixedHeight(lHeight);
	L01->setFixedWidth(lWidth);
	L01->setAlignment(Qt::AlignCenter);	
	L02->setFixedHeight(lHeight);
	L02->setFixedWidth(lWidth);
	L02->setText("<b>Dimer</b>");L02->setAlignment(Qt::AlignBottom|Qt::AlignCenter);
	//
	L10->setText("<b>Hen Egg White Lysozyme</b>");
	L10->setFixedHeight(lHeight);
	L10->setFixedWidth(lWidth);
	L10->setAlignment(Qt::AlignRight);
	L10->setWordWrap(true);
	L11->setFixedHeight(imgHeight);
	L11->setFixedWidth(imgWidth);
	L11->setAlignment(Qt::AlignCenter);L11->setMovie(s11Movie);L11->setScaledContents(true);s11Movie->start();
	L12->setFixedHeight(imgHeight);
	L12->setFixedWidth(imgWidth);
	L12->setAlignment(Qt::AlignCenter);L12->setMovie(s12Movie);L12->setScaledContents(true);s12Movie->start();
	//
	L20->setText("<b>Bovine &beta;-Lactoglobulin</b>");
	L20->setFixedHeight(lHeight);
	L20->setFixedWidth(lWidth);
	L20->setAlignment(Qt::AlignRight);
	L20->setWordWrap(true);
	L21->setFixedHeight(imgHeight);
	L21->setFixedWidth(imgWidth);
	L21->setAlignment(Qt::AlignCenter);L21->setMovie(s21Movie);L21->setScaledContents(true);s21Movie->start();
	L22->setFixedHeight(imgHeight);
	L22->setFixedWidth(imgWidth);
	L22->setAlignment(Qt::AlignCenter);L22->setMovie(s22Movie);L22->setScaledContents(true);s22Movie->start();
	//
	L30->setText("<b>Bovine Serum Albumin</b>");
	L30->setFixedHeight(lHeight);
	L30->setFixedWidth(lWidth);
	L30->setAlignment(Qt::AlignRight);
	L30->setWordWrap(true);
	L31->setFixedHeight(imgHeight);
	L31->setFixedWidth(imgWidth);
	L31->setAlignment(Qt::AlignCenter);L31->setMovie(s31Movie);L31->setScaledContents(true);s31Movie->start();
	L32->setFixedHeight(imgHeight);
	L32->setFixedWidth(imgWidth);
	L32->setAlignment(Qt::AlignCenter);L32->setMovie(s32Movie);L32->setScaledContents(true);s32Movie->start();

	QGroupBox* inBox=new QGroupBox;
	QGridLayout *lo1=new QGridLayout;
	lo1->setSizeConstraint(QLayout::SetFixedSize);
	lo1->addWidget(L00,0,0,1,1);lo1->addWidget(L01,0,1,1,1);lo1->addWidget(L02,0,2,1,1);
	lo1->addWidget(L10,1,0,1,1);lo1->addWidget(L11,1,1,1,1);lo1->addWidget(L12,1,2,1,1);
	lo1->addWidget(L20,2,0,1,1);lo1->addWidget(L21,2,1,1,1);lo1->addWidget(L22,2,2,1,1);
	lo1->addWidget(L30,3,0,1,1);lo1->addWidget(L31,3,1,1,1);lo1->addWidget(L32,3,2,1,1);
	inBox->setLayout(lo1);
	// Scrolling Group Box
	QScrollArea* scrollArea=new QScrollArea();
	scrollArea->setWidget(inBox);  
	scrollArea->setWidgetResizable(true);  
	QGridLayout *outlo=new QGridLayout;
	outlo->addWidget(scrollArea);
	QGroupBox* topBox=new QGroupBox;
	topBox->setLayout(outlo);
	topBox->setMinimumHeight(imgHeight);
	// Charge Validation
	//1. Validation of Charge (Q) Computation.
	s="<b>PROPKA 3.0</b>[24] is the software component responsible for the assignment of charge to a given molecular structure. Experimental validation based on hydrogen ion titrations of three different protein structures (<b>hen egg white lysozyme</b>, <b>bovine serum albumin</b>, and <b>bovine &beta;-lactoglobulin</b>) are presented below. The charge of a structure depends on the pH of its solution, and thus validation involves comparison of experimentally and computationally determined charges across a broad range of pH values. To account for possible dimerization during experimental titrations, <b>PROPKA 3.0</b> was applied to both monomeric and dimeric forms of each of the proteins assessed.";
	QLabel* chargeDesc=new QLabel;
	chargeDesc->setText(s.c_str());
	chargeDesc->setFixedWidth(propertyBoxSz);
	chargeDesc->setWordWrap(true);

	QLabel* cL00=new QLabel;QLabel* cL10=new QLabel;QLabel* cL20=new QLabel;QLabel* cL30=new QLabel;QLabel* cL40=new QLabel;QLabel* cL50=new QLabel;QLabel* cL60=new QLabel;
	cL00->setFixedHeight(figHeight);
	cL00->setFixedWidth(figWidth);
	cL00->setAlignment(Qt::AlignCenter);cL00->setPixmap(fig1.scaled(figWidth,figHeight,Qt::KeepAspectRatio));cL00->setMask(fig1.mask());
	s="<b>Figure 1. Structures of the Hen Egg White Lysozyme Monomer </b>(<b>A</b>; PDB id: <a href=\"http://www.rcsb.org/structure/6lyz\">6lyz</a>) <b>and Dimer</b> (<b>B</b>; PDB id: <a href=\"http://www.rcsb.org/structure/4r0f\">4r0f</a>) <b>and Their Hydrogen Ion Titration Curves (C).</b> Experimental values came from[84]. Structures were obtained from the Protein Data Bank (PDB)[85] and their charges were calculated using PROPKA 3.0[24] on an ensemble of structures sampled from molecular dynamics.<br>";
	cL10->setText(s.c_str());
	cL10->setFixedHeight(lHeight);
	cL10->setFixedWidth(figWidth);
	cL10->setWordWrap(true);
	//
	cL20->setFixedHeight(figHeight);
	cL20->setFixedWidth(figWidth);
	cL20->setAlignment(Qt::AlignCenter);cL20->setPixmap(fig2.scaled(figWidth,figHeight,Qt::KeepAspectRatio));cL20->setMask(fig2.mask());
	s="<b>Figure 2. Structures of the Bovine Serum Albumin Monomer </b>(<b>A</b>; PDB id: <a href=\"http://www.rcsb.org/structure/4f5s\">4f5s (Chain A)</a>) <b>and Dimer</b> (<b>B</b>; PDB id: <a href=\"http://www.rcsb.org/structure/4f5s\">4f5s</a>) <b>and Their Hydrogen Ion Titration Curves (C).</b> Experimental values came from [86]. Structures were obtained from the Protein Data Bank (PDB)[85] and their charges were calculated using PROPKA 3.0[24] on an ensemble of structures sampled from molecular dynamics.<br>";
	cL30->setText(s.c_str());
	cL30->setFixedHeight(lHeight);
	cL30->setFixedWidth(figWidth);
	cL30->setWordWrap(true);
	//
	cL40->setFixedHeight(figHeight);
	cL40->setFixedWidth(figWidth);
	cL40->setAlignment(Qt::AlignCenter);cL40->setPixmap(fig3.scaled(figWidth,figHeight,Qt::KeepAspectRatio));cL40->setMask(fig3.mask());
	s="<b>Figure 3. Structures of the Bovine &#946-Lactoglobulin Monomer </b>(<b>A</b>; PDB id: <a href=\"http://www.rcsb.org/structure/3blg\">3blg</a>) <b>and Dimer</b> (<b>B</b>; PDB id: <a href=\"http://www.rcsb.org/structure/1beb\">1beb</a>) <b>and Their Hydrogen Ion Titration Curves (C).</b> Experimental values came from [87]. Structures were obtained from the Protein Data Bank (PDB)[85] and their charges were calculated using PROPKA 3.0[24] on an ensemble of structures sampled from molecular dynamics.<br>";
	cL50->setText(s.c_str());
	cL50->setFixedWidth(figWidth);
	cL50->setWordWrap(true);
	//
	s="As shown above in <b>Figs. 1-3</b>, excellent agreement between experimental and computational values was achieved by <b>PROPKA 3.0</b>[24]. In <b>Fig. 3</b>, experimental values for &beta;-Lactoglobulin deviate from computational values around the isoelectric point (pH 5.4[87]), which is expected as the structure can aggregate into an octamer in these conditions. There is currently no solved structure for the octamer oligomerization state. Additional experimental validation of PROPKA can be found elsewhere[25].<br><br>";
	cL60->setText(s.c_str());
	cL60->setFixedWidth(propertyBoxSz);
	cL60->setWordWrap(true);
	//
	QGroupBox* inCBox=new QGroupBox;
	QGridLayout *lo2=new QGridLayout;
	lo2->setSizeConstraint(QLayout::SetFixedSize);
	lo2->addWidget(cL00,0,0);
	lo2->addWidget(cL10,1,0);
	lo2->addWidget(cL20,2,0);
	lo2->addWidget(cL30,3,0);
	lo2->addWidget(cL40,4,0);
	lo2->addWidget(cL50,5,0);
	inCBox->setLayout(lo2);
	// Scrolling Group Box
	QScrollArea* sA=new QScrollArea();
	sA->setWidget(inCBox);  
	sA->setWidgetResizable(true);  
	QGridLayout *outlo2=new QGridLayout;
	//outlo2->setSizeConstraint(QLayout::SetFixedSize);
	outlo2->addWidget(sA);
	QGroupBox* outCBox=new QGroupBox;
	outCBox->setLayout(outlo2);
	outCBox->setMinimumHeight(figHeight);
	
	QGroupBox* cBox=new QGroupBox(tr("1. Validation of Charge (Q) Computation"));
	QGridLayout *cLo=new QGridLayout;
	cLo->addWidget(chargeDesc,0,0);
	cLo->addWidget(outCBox,1,0);
	cLo->addWidget(cL60,2,0);
	cBox->setLayout(cLo);
	//
	// Final Scrolling Group Box
	QScrollArea* lastsA=new QScrollArea();
	lastsA->setWidget(cBox);
	lastsA->setWidgetResizable(true);  
	QGridLayout *lastoutlo=new QGridLayout;
	lastoutlo->addWidget(lastsA);
	QGroupBox* scrollBox=new QGroupBox;
	scrollBox->setLayout(lastoutlo);
	//scrollBox->setMinimumHeight(figHeight);
	//
	QGridLayout *mLayout = new QGridLayout;
	//mLayout->setSizeConstraint(QLayout::SetFixedSize);
	mLayout->addWidget(desc,0,0);
	mLayout->addWidget(topBox,1,0);
	mLayout->addWidget(scrollBox,2,0);
	//mLayout->addWidget(cBox,2,0);

	setLayout(mLayout);
	setWindowTitle(tr("Experimental Methods and Validation"));
	resize(QDesktopWidget().availableGeometry(this).size());
}

zpredInfo::zpredInfo(QWidget *parent) : QWidget(parent)
{	
	int imgWidth=400;
	int imgHeight=400;
	int arwWidth=100;
	int arwHeight=100;
	int propertyBoxSz=1150;
	// Get Current File Path and Find Images
	QString theCurFldr=QDir::currentPath();
	string Fldr=theCurFldr.toStdString()+"/";
	// Images
	// Step 1 Animation File
	string s1File=Fldr+"UI/Images/ZPRED_Step_1.gif";
	// Step 2 Animation File
	string s2File=Fldr+"UI/Images/ZPRED_Step_2.png";
	// Step 3 Animation File
	string s3File=Fldr+"UI/Images/ZPRED_Step_3.png";
	// Step 4 Animation File
	string s4File=Fldr+"UI/Images/ZPRED_Step_4.gif";
	// Step 5 Animation File
	string s5File=Fldr+"UI/Images/ZPRED_Step_5.gif";
	// Step 6 Animation File
	string s6File=Fldr+"UI/Images/ZPRED_Step_6.gif";
	// Right Arrow
	string raFile=Fldr+"UI/Images/right_arrow.png";
	// Left Arrow
	string laFile=Fldr+"UI/Images/left_arrow.png";
	// Dpwn Arrow
	string daFile=Fldr+"UI/Images/down_arrow.png";
	// Long Left Arrow
	string llaFile=Fldr+"UI/Images/long_left_arrow.png";

	//QPixmap s1Pixmap(s1File.c_str());
	QMovie *s1Movie = new QMovie(s1File.c_str());
	QPixmap s2Pixmap(s2File.c_str());
	QPixmap s3Pixmap(s3File.c_str());
	QMovie *s4Movie = new QMovie(s4File.c_str());
	QMovie *s5Movie = new QMovie(s5File.c_str());
	QMovie *s6Movie = new QMovie(s6File.c_str());
	QPixmap raPixmap(raFile.c_str());
	QPixmap laPixmap(laFile.c_str());
	QPixmap daPixmap(daFile.c_str());
	QPixmap llaPixmap(llaFile.c_str());

	// Interface Description
	string s="<b>This interface shows how to properly use ZPRED to compute the zeta potential of a protein crystal structure.  As shown below, accurate computation is achieved through six primary steps.  Users of ZPRED are expected to already have completed the first two steps, which basically involves acquiring realistic conformations of the protein structure in solution.  This can be done using molecular dynamics software as described below.  Once structures are acquired and put into .pdb format, you are ready to use ZPRED.  For convenience, all required user inputs are highlighted in <span style=\"color:red;\">red</span> and go black once successfully input. First, use the \'Select Structure Files\' button to choose the ensemble of structures representing your protien in solution.  Then design your protein solution by specifying one or more solution conditions, which as mentioned previously will be initially highlighted in <span style=\"color:red;\">red</span>.</b><br>";
	QLabel* desc=new QLabel;
	desc->setText(s.c_str());
	desc->setFixedWidth(propertyBoxSz);
	desc->setWordWrap(true);

	// Specify External Executables
	QGroupBox* s1Box=new QGroupBox;QGroupBox* s1Box4=new QGroupBox;
	QGroupBox* s2Box=new QGroupBox;QGroupBox* s2Box4=new QGroupBox;
	QGroupBox* s3Box=new QGroupBox;QGroupBox* s3Box4=new QGroupBox;
	QGroupBox* s4Box=new QGroupBox;QGroupBox* s4Box4=new QGroupBox;
	QGroupBox* s5Box=new QGroupBox;QGroupBox* s5Box4=new QGroupBox;
	QGroupBox* s6Box=new QGroupBox;QGroupBox* s6Box4=new QGroupBox;
	QGroupBox* sBox=new QGroupBox;
	QLabel* L11=new QLabel;QLabel* L12=new QLabel;QLabel* L13=new QLabel;QLabel* L14=new QLabel;
	L11->setText("<b><span style=\"background-color:grey; color:black; font-size:12pt;\">(#1) Simulate Structural Motion of Molecule in Solution</span></b>");L11->setAlignment(Qt::AlignCenter);	
	L12->setFixedHeight(imgHeight);
	L12->setFixedWidth(imgWidth);
	L12->setAlignment(Qt::AlignCenter);L12->setMovie(s1Movie);L12->setScaledContents(true);s1Movie->start();
	L13->setText("<span style=\"background-color:grey; color:black; font-size:12pt;\"><b>Molecular Dynamics Simulation of Hen Egg White Lysozyme.</b> Solvent water molecules were removed for a better view of the protein structure. Rotational motion is exaggerated from the code making the gif formatted image.</span>");
	L13->setFixedWidth(imgWidth);
	L13->setWordWrap(true);
	s="The first step of ZPRED uses molecular dynamics software (e.g. Amber 2015) to simulate the structural motions of the molecule in solvent. In general, this step is comprised of four parts listed below:<ol><li>prepare the molecular structure for a molecular dynamics simulation</li><li>energetically minimize the molecular structure</li><li>thermally excite the molecular structure</li><li>simulate the molecular structure in explicit solvent</li></ol>";
	s+="In the <b>first part</b>, atomic coordinates of a molecular structure (e.g. a crystal structure from the <a href=\"https://www.rcsb.org/\">protein data bank (PDB)</a>) are prepared by removing any water molecules (e.g. using the Amber tool, pdb4amber) and protonating the structure (e.g. using the Amber tool, reduce). Prepared structures are loaded into a molecular dynamics simulation (e.g. as an UNIT object manipulated by the Amber program, teLeap) specifying force field parameters and generating initial topology and coordinates of the atoms of the prepared structure in a specified volume of solvent molecules. In general, globular (spherical) molecules are housed in water boxes extending 20 Å from the molecular surface and fibrillar (cylindrical) molecules are housed in water boxes extending 30 Å away. The <b>second part</b> takes the generated topology and coordinate files and performs a molecular dynamics simulation (e.g. using either sander or pmemd of Amber) to energetically minimize the structure, which optimizes its geometry in solution. The coordinates of the optimized structure provide a starting point for the simulation of the part 3. The <b>third part</b> gradually heats the structure from 0 K to a specified temperature, inducing thermal motion of the solvent and the molecule.<br><br>";
	s+="The final <b>fourth part</b> is the main molecular dynamics simulation and uses the coordinates of the prepared heated structure as input. This simulation is run until the structure reaches a steady-state (about 100 nanoseconds) based on the root mean squared displacement of the primary structural chain (e.g. a protein backbone). Structural steady-state should be assessed after centering the entire trajectory of the solvent and atomic coordinates around the molecule’s center of mass (e.g. using the Amber tool, cpptraj). Once a structural steady-state is reached, the molecule should switch between a limited number of molecular conformations, which are sampled based on the variation in the root mean squared displacement.";
	L14->setText(s.c_str());
	L14->setWordWrap(true);
	L14->setFixedWidth(imgWidth);
	QGridLayout *slo1=new QGridLayout;
	slo1->addWidget(L14,0,0);
	s1Box4->setLayout(slo1);

	// Scrolling Group Box
	QScrollArea* sA1=new QScrollArea();
	sA1->setWidget(s1Box4);  
	sA1->setWidgetResizable(true);  
	QGridLayout *soutlo1=new QGridLayout;
	soutlo1->setSizeConstraint(QLayout::SetFixedSize);
	soutlo1->addWidget(sA1);
	QGroupBox* s1OutBox=new QGroupBox;
	s1OutBox->setLayout(soutlo1);
	//
	QGridLayout *lo1=new QGridLayout;
	lo1->setSizeConstraint(QLayout::SetFixedSize);
	lo1->addWidget(L11,0,0);
	lo1->addWidget(L12,1,0);
	lo1->addWidget(L13,2,0);
	lo1->addWidget(s1OutBox,3,0);
	s1Box->setLayout(lo1);
	s1Box->setStyleSheet("background-color:grey;");
	// Step 2
	QLabel* L21=new QLabel;QLabel* L22=new QLabel;QLabel* L23=new QLabel;QLabel* L24=new QLabel;
	L21->setText("<b><span style=\"background-color:grey; color:black; font-size:12pt;\">(#2) Sample Molecular Conformations from Simulation</span></b>");L21->setAlignment(Qt::AlignCenter);	
	L22->setFixedHeight(imgHeight);
	L22->setFixedWidth(imgWidth);
	L22->setAlignment(Qt::AlignCenter);L22->setPixmap(s2Pixmap.scaled(imgWidth,imgHeight,Qt::KeepAspectRatio));L22->setMask(s2Pixmap.mask());
	L23->setText("<span style=\"background-color:grey; color:black; font-size:12pt;\"><b>Sampling of Hen Egg White Lysozyme Conformations Based on Root Mean Squared Deviation</b></span>");
	L23->setFixedWidth(imgWidth);
	L23->setWordWrap(true);
	s="There exist many clustering algorithms[72-80] for identifying and sampling molecular conformations in the most appropriate folded state from molecular dynamics simulations. For its experimental validation, ZPRED applied a simple structural sampling protocol during molecular dynamics, which allowed the simulation to run until the structure reached a steady-state based on the root mean squared displacement of the primary structural chain (the protein backbone). Structural steady-state was assessed after centering the entire trajectory of the solvent and atomic coordinates around the molecule’s center of mass (e.g. using the Amber tool, cpptraj). Once a structural steady-state was reached (as shown above), the molecule switched between a limited number of molecular conformations, which were sampled based on the variation in the root mean squared displacement. For consistency in this demonstration, simulations were run for 100 nanoseconds and twenty conformations were sampled from the last 20 nanoseconds, one per each nanosecond.<br>";
	L24->setText(s.c_str());
	L24->setWordWrap(true);
	L24->setFixedWidth(imgWidth);
	QGridLayout *slo2=new QGridLayout;
	slo2->addWidget(L24,0,0);
	s2Box4->setLayout(slo2);
	// Scrolling Group Box
	QScrollArea* sA2=new QScrollArea();
	sA2->setWidget(s2Box4);  
	sA2->setWidgetResizable(true);  
	QGridLayout *soutlo2=new QGridLayout;
	soutlo2->setSizeConstraint(QLayout::SetFixedSize);
	soutlo2->addWidget(sA2);
	QGroupBox* s2OutBox=new QGroupBox;
	s2OutBox->setLayout(soutlo2);
	//
	QGridLayout *lo2=new QGridLayout;
	lo2->setSizeConstraint(QLayout::SetFixedSize);
	lo2->addWidget(L21,0,0);
	lo2->addWidget(L22,1,0);
	lo2->addWidget(L23,2,0);
	lo2->addWidget(s2OutBox,3,0);
	s2Box->setLayout(lo2);
	s2Box->setStyleSheet("background-color:grey;");
	////////// Step 3
	QLabel* L31=new QLabel;QLabel* L32=new QLabel;QLabel* L33=new QLabel;QLabel* L34=new QLabel;
	L31->setText("<b><span style=\"background-color:grey; color:black; font-size:12pt;\">(#3) Assign Atomic Charge & Radii</span></b>");L31->setAlignment(Qt::AlignCenter);	
	L32->setFixedHeight(imgHeight);
	L32->setFixedWidth(imgWidth);
	L32->setAlignment(Qt::AlignCenter);L32->setPixmap(s3Pixmap.scaled(imgWidth,imgHeight,Qt::KeepAspectRatio));L32->setMask(s2Pixmap.mask());
	L33->setText("<span style=\"background-color:grey; color:black; font-size:12pt;\"><b>Hen Egg White Lysozyme PQR Structure with Defined Atomic Charges and Radii.</b></span>");
	L33->setFixedWidth(imgWidth);
	L33->setWordWrap(true);
	s="The third step of ZPRED assigns a charge and radius to each atomic coordinate of the molecular structure using available software. The software, PDB2PQR[22][23], converts a PDB formatted coordinates file into a PQR formatted file for use in electrostatics software. This involves checking the integrity of the structure (e.g. whether heavy atoms are missing or not) and then protonating it based on a pKa predictor (PROPKA 3.0[24]) at a specified pH. It’s worth noting, PROPKA 3.0 has proven to be a more accurate method for pKa prediction among other pKa predictors (MCCE, MEAD and UHBD)[25]. Following protonation, the position of hydrogens are determined by Monte Carlo optimization based on the global H-bonding network of the structure considering charge residue side chains and water-molecule interactions. Once atomic charges and radii have been assigned, the structure is ready for electrostatic calculations.";
	L34->setText(s.c_str());
	L34->setWordWrap(true);
	L34->setFixedWidth(imgWidth);
	QGridLayout *slo3=new QGridLayout;
	slo3->addWidget(L34,0,0);
	s3Box4->setLayout(slo3);
	// Scrolling Group Box
	QScrollArea* sA3=new QScrollArea();
	sA3->setWidget(s3Box4);  
	sA3->setWidgetResizable(true);  
	QGridLayout *soutlo3=new QGridLayout;
	soutlo3->setSizeConstraint(QLayout::SetFixedSize);
	soutlo3->addWidget(sA3);
	QGroupBox* s3OutBox=new QGroupBox;
	s3OutBox->setLayout(soutlo3);
	//
	QGridLayout *lo3=new QGridLayout;
	lo3->setSizeConstraint(QLayout::SetFixedSize);
	lo3->addWidget(L31,0,0);
	lo3->addWidget(L32,1,0);
	lo3->addWidget(L33,2,0);
	lo3->addWidget(s3OutBox,3,0);
	s3Box->setLayout(lo3);
	s3Box->setStyleSheet("background-color:grey;");
	/////////// Step 4
	QLabel* L41=new QLabel;QLabel* L42=new QLabel;QLabel* L43=new QLabel;QLabel* L44=new QLabel;
	L41->setText("<b><span style=\"background-color:grey; color:black; font-size:12pt;\">(#4) Compute Electric Potential Distribution</span></b>");L41->setAlignment(Qt::AlignCenter);	
	L42->setFixedHeight(imgHeight);
	L42->setFixedWidth(imgWidth);
	L42->setAlignment(Qt::AlignCenter);L42->setMovie(s4Movie);L42->setScaledContents(true);s4Movie->start();
	L43->setText("<span style=\"background-color:grey; color:black; font-size:12pt;\"><b>Computed Electric Potential Distribution of Hen Egg White Lysozyme.</b> Electric potential values (above) have been scaled by thermal energy (zeφ/k<sub>b</sub>T) and are colored based on the scale</span>");
	L43->setFixedWidth(imgWidth);
	L43->setWordWrap(true);
	s="The fourth step of ZPRED solves the Poisson-Boltzmann equation (PBE) over the structure, theoretically modeling the diffuse region of a Gouy-Chapman EDL around the structure. This is accomplished with the adaptive Poisson-Boltzmann solver (APBS)[26]. To solve the PBE, the problem domain is divided into two regions of different dielectrics: the molecule and the solvent. As shown below, the two regions are separated by a solvent-accessible surface generated over the molecular structure using the largest ion in the solvent. Thus, ionic radii are needed as an input.";
	L44->setText(s.c_str());
	L44->setWordWrap(true);
	L44->setFixedWidth(imgWidth);
	QGridLayout *slo4=new QGridLayout;
	slo4->addWidget(L44,0,0);
	s4Box4->setLayout(slo4);
	// Scrolling Group Box
	QScrollArea* sA4=new QScrollArea();
	sA4->setWidget(s4Box4);  
	sA4->setWidgetResizable(true);  
	QGridLayout *soutlo4=new QGridLayout;
	soutlo4->setSizeConstraint(QLayout::SetFixedSize);
	soutlo4->addWidget(sA4);
	QGroupBox* s4OutBox=new QGroupBox;
	s4OutBox->setLayout(soutlo4);
	//
	QGridLayout *lo4=new QGridLayout;
	lo4->setSizeConstraint(QLayout::SetFixedSize);
	lo4->addWidget(L41,0,0);
	lo4->addWidget(L42,1,0);
	lo4->addWidget(L43,2,0);
	lo4->addWidget(s4OutBox,3,0);
	s4Box->setLayout(lo4);
	s4Box->setStyleSheet("background-color:grey;");
	// Step 5
	QLabel* L51=new QLabel;QLabel* L52=new QLabel;QLabel* L53=new QLabel;QLabel* L54=new QLabel;
	L51->setText("<b><span style=\"background-color:grey; color:black; font-size:12pt;\">(#5) Estimate Slip Plane Position (X<sub>SP</sub>)</span></b>");L51->setAlignment(Qt::AlignCenter);	
	L52->setFixedHeight(imgHeight);
	L52->setFixedWidth(imgWidth);
	L52->setAlignment(Qt::AlignCenter);L52->setMovie(s5Movie);L52->setScaledContents(true);s5Movie->start();
	L53->setText("<span style=\"background-color:grey; color:black; font-size:12pt;\"><b>Estimating the Slip Plane Position (X<sub>SP</sub>) with Hydrodynamic Software.</b></span>");
	L53->setFixedWidth(imgWidth);
	L53->setWordWrap(true);
	s="The fifth step of ZPRED estimates the thickness of the hydration layer around the protein molecule. The position of the slip plane (i.e. edge of hydration) relative to the molecular surface must be either determined from experimental data or estimated computationally. The Stokes-Einstein hydrodynamic radius[30] (<b>R<sub>h</sub></b>) determined from measured diffusivities, and the electrophoretic radius[44] (<b>R<sub>e</sub></b>) determined from measured electrophoretic mobilities provide reasonable representations of the slip plane position. ZPRED assumes hydration is distributed based on the “uniform hydration layer model”[4]. Thus, the slip plane position (<b>X<sub>SP</sub>=R<sub>h</sub>-R<sub>P</sub></b>) can be estimated experimentally by subtracting the anhydrous molecular radius (<b>R<sub>P</sub></b>) from a measured solvated radius (a represents a general solvated radius and can be either <b>R<sub>h</sub></b> or <b>R<sub>e</sub></b>), which should always be greater than or equal to the molecular radius.<br>";
	s+="<b><span style=\"color:black;font-size:12pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;k<sub>b</sub>T<br></span></b>";
	s+="<b><span style=\"color:black;font-size:12pt;\">R<sub>h</sub> = </span><span style=\"color:black;font-size:4pt;\">----------------------------------<br></span></b>";
	s+="<b><span style=\"color:black;font-size:12pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;6&pi;&eta;D<sub>o</sub><br></span></b>";
	s+="where <b>k<sub>b</sub></b> is the Boltzmann constant, <b>T</b> is absolute temperature, <b>&eta;</b> is the protein-free solution viscosity and <b>D<sub>o</sub></b> is the single particle/protein diffusivity<br><br>";
	s+="It is important to note the Stokes-Einstein equation[30] (above) is limited to relatively high salt concentrations, and thus, other methods for determining molecular size must be used, such as the <b>R<sub>e</sub></b> determined from electrophoretic mobility measurements (below).<br>";
	s+="<b><span style=\"color:black;font-size:12pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Qef(&kappa;R<sub>e</sub>)<br></span></b>";
	s+="<b><span style=\"color:black;font-size:12pt;\">R<sub>e</sub> = </span><span style=\"color:black;font-size:4pt;\">--------------------------------------------------------------<br></span></b>";
	s+="<b><span style=\"color:black;font-size:12pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;6&pi;&eta;u<sub>e</sub>(1+&kappa;R<sub>e</sub>)f<sub>s</sub><br></span></b>";
	//
	s+="<b><span style=\"color:black;font-size:12pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1<br></span></b>";
	s+="<b><span style=\"color:black;font-size:12pt;\">f(&kappa;R<sub>e</sub>) = 1 + </span><span style=\"color:black;font-size:4pt;\">----------------------------<br></span></b>";
	s+="<b><span style=\"color:black;font-size:12pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2X<sup>3</sup><br></span></b>";
	//
	s+="<b><span style=\"color:black;font-size:12pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Y<br></span></b>";
	s+="<b><span style=\"color:black;font-size:12pt;\">X = 1 + </span><span style=\"color:black;font-size:4pt;\">----------------<br></span></b>";
	s+="<b><span style=\"color:black;font-size:12pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&kappa;R<sub>e</sub><br></span></b>";
	//
	s+="<b><span style=\"color:black;font-size:12pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2.5<br></span></b>";
	s+="<b><span style=\"color:black;font-size:12pt;\">Y = </span><span style=\"color:black;font-size:4pt;\">-------------------------------------------------------<br></span></b>";
	s+="<b><span style=\"color:black;font-size:12pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1+2exp(-&kappa;R<sub>e</sub>)<br></span></b>";
	//
	s+="where <b>Q</b> is the protein net charge valence, <b>e</b> is the electron charge, <b>f(&kappa;R<sub>e</sub>)</b> is the Henry correction factor for electrophoretic retardation (defined above)[12], <b>&eta;</b> is the protein-free solution viscosity, <b>u<sub>e</sub></b> is the protein electrophoretic mobility, <b>&kappa;</b> is the inverse Debye length and <b>f<sub>s</sub></b> is a frictional shape factor[28]<br><br>";
	//
	s+="Estimating the <b>X<sub>SP</sub></b> computationally requires estimating the molecular anhydrous radius (<b>R<sub>P</sub></b>) and <b>a</b>, the solvated radius. For globular molecules, <b>R<sub>P</sub></b> can be calculated as the average distance between the center of mass and the solvent-excluded surface of the structure under assessment.  The solvated radius (<b>a</b>) can be determined computationally as <b>R<sub>h</sub></b> or <b>R<sub>e</sub></b> depending entirely on the user's preference. Here a discussion of the computation of <b>R<sub>h</sub></b> is presented.<br><br>";
	s+="<b>R<sub>h</sub></b> depends on temperature, which is controlled by the user; leaving protein-free solution viscosity and single particle/protein diffusivity to be defined. ZPRED uses a number of empirical relationships for the pure solvent viscosity. If values cannot be determined, the viscosity of pure solvent (e.g. water) can be used as an estimate since added salt only affects viscosity at higher salt concentrations. Single particle diffusivity can be computed with the hydrodynamics software, HYDROPRO[27], which requires the following inputs:<ul><li>protein structure</li><li>its specific volume</li><li>its molecular weight</li><li>temperature</li><li>protein-free solution density</li><li>protein-free solution viscosity</li></ul>";
	s+="Molecular weight is computed from the summation of the molecular weights of atoms present. Specific volume of each structure requires calculating the volume of the structure (using <b>MSMS</b>[29]) and then dividing by the mass of the structure. Protein-free solution density and viscosity are computed from empirical relationships described in the Solution Property Computation Widget. If values can not be determined, the density and viscosity of pure solvent can be used as an estimate. Once viscosity and diffusivity values are obtained, the <b>R<sub>h</sub></b> is calculated with the Stokes-Einstein Equation. Alternatively, <b>R<sub>e</sub></b> can be used instead if the electrophoretic mobility is easier to acquire than the diffusivity. Just like experimental values, the estimated <b>X<sub>SP</sub></b> is calculated by subtracting the molecular radius from the estimated solvated raius (<b>a</b>) (e.g. <b>R<sub>h</sub></b> or <b>R<sub>e</sub></b>)<br>";
	s+="<b><span style=\"color:black;font-size:12pt;\">X<sub>SP</sub> = a - R<sub>P</sub><br></span></b>";
	s+="where <b>X<sub>SP</sub></b> is the slip plane position (i.e. hydration layer thickness), <b>a</b> is a solvated radius and <b>R<sub>P</sub></b> is the anhydrous protein radius";
	L54->setText(s.c_str());
	L54->setWordWrap(true);
	L54->setFixedWidth(imgWidth);
	QGridLayout *slo5=new QGridLayout;
	slo5->addWidget(L54,0,0);
	s5Box4->setLayout(slo5);
	// Scrolling Group Box
	QScrollArea* sA5=new QScrollArea();
	sA5->setWidget(s5Box4);  
	sA5->setWidgetResizable(true);  
	QGridLayout *soutlo5=new QGridLayout;
	soutlo5->setSizeConstraint(QLayout::SetFixedSize);
	soutlo5->addWidget(sA5);
	QGroupBox* s5OutBox=new QGroupBox;
	s5OutBox->setLayout(soutlo5);
	//
	QGridLayout *lo5=new QGridLayout;
	lo5->setSizeConstraint(QLayout::SetFixedSize);
	lo5->addWidget(L51,0,0);
	lo5->addWidget(L52,1,0);
	lo5->addWidget(L53,2,0);
	lo5->addWidget(s5OutBox,3,0);
	s5Box->setLayout(lo5);
	s5Box->setStyleSheet("background-color:grey;");
	// Step 6
	QLabel* L61=new QLabel;QLabel* L62=new QLabel;QLabel* L63=new QLabel;QLabel* L64=new QLabel;
	L61->setText("<b><span style=\"background-color:grey; color:black; font-size:12pt;\">(#6) Compute Zeta Potential of Molecular Structure</span></b>");L61->setAlignment(Qt::AlignCenter);	
	L62->setFixedHeight(imgHeight);
	L62->setFixedWidth(imgWidth);
	L62->setAlignment(Qt::AlignCenter);L62->setMovie(s6Movie);L62->setScaledContents(true);s6Movie->start();
	L63->setText("<span style=\"background-color:grey; color:black; font-size:12pt;\"><b>Capturing the Zeta Potential on Hen Egg White Lysozyme.<br></b></span>");
	L63->setFixedWidth(imgWidth);
	L63->setWordWrap(true);
	s="The final step of ZPRED averages the electric potentials computed by APBS at the estimated slip plane position. This involves generating a solvent-excluded surface (SES) over the molecular structure (using <b>MSMS</b>[29]). The SES generated is composed of Cartesian coordinates and their normal vectors directed away from the molecular surface. The SES is inflated to the <b>X<sub>SP</sub></b> by translating its initial coordinates along their respective normal vector by the estimated slip plane distance. Then the calculated electric potentials at the inflated coordinates are captured (using APBS[26] tool, <b>multivalue</b>). Using <b>multivalue</b> requires converting the inflated coordinates into a comma separated vector (CSV) file format, which is simply done by writing each coordinate on its own line and delimiting by commas in a text file. A <b>&#950;</b> value for each conformation is computed by averaging the captured potentials at the inflated SES. A final zeta potential output value is obtained by averaging the zeta potentials determined from each conformation. The final value represents what would be expected from the molecular structure in solution.";
	L64->setText(s.c_str());
	L64->setFixedWidth(imgWidth);
	L64->setWordWrap(true);
	QGridLayout *slo6=new QGridLayout;
	slo6->addWidget(L64,0,0);
	s6Box4->setLayout(slo6);
	// Scrolling Group Box
	QScrollArea* sA6=new QScrollArea();
	sA6->setWidget(s6Box4);  
	sA6->setWidgetResizable(true);  
	QGridLayout *soutlo6=new QGridLayout;
	soutlo6->setSizeConstraint(QLayout::SetFixedSize);
	soutlo6->addWidget(sA6);
	QGroupBox* s6OutBox=new QGroupBox;
	s6OutBox->setLayout(soutlo6);
	//
	QGridLayout *lo6=new QGridLayout;
	lo6->setSizeConstraint(QLayout::SetFixedSize);
	lo6->addWidget(L61,0,0);
	lo6->addWidget(L62,1,0);
	lo6->addWidget(L63,2,0);
	lo6->addWidget(s6OutBox,3,0);
	s6Box->setLayout(lo6);
	s6Box->setStyleSheet("background-color:grey;");
	// Right Arrow
	QLabel* Lra=new QLabel;QLabel* Lra2=new QLabel;
	Lra->setFixedHeight(arwHeight);
	Lra->setFixedWidth(arwWidth);
	Lra->setAlignment(Qt::AlignCenter);Lra->setPixmap(raPixmap.scaled(arwWidth,arwHeight,Qt::KeepAspectRatio));Lra->setMask(raPixmap.mask());
	Lra2->setFixedHeight(arwHeight);
	Lra2->setFixedWidth(arwWidth);
	Lra2->setAlignment(Qt::AlignCenter);Lra2->setPixmap(raPixmap.scaled(arwWidth,arwHeight,Qt::KeepAspectRatio));Lra2->setMask(raPixmap.mask());
	// Leftt Arrow
	QLabel* Lla=new QLabel;
	Lla->setFixedHeight(arwHeight);
	Lla->setFixedWidth(arwWidth);
	Lla->setAlignment(Qt::AlignCenter);Lla->setPixmap(laPixmap.scaled(arwWidth,arwHeight,Qt::KeepAspectRatio));Lla->setMask(laPixmap.mask());
	// Down Arrow
	QLabel* Lda=new QLabel;QLabel* Lda2=new QLabel;
	Lda->setFixedHeight(arwHeight);
	Lda->setFixedWidth(arwWidth);
	Lda->setAlignment(Qt::AlignCenter);Lda->setPixmap(daPixmap.scaled(arwWidth,arwHeight,Qt::KeepAspectRatio));Lda->setMask(daPixmap.mask());
	Lda2->setFixedHeight(arwHeight);
	Lda2->setFixedWidth(arwWidth);
	Lda2->setAlignment(Qt::AlignCenter);Lda2->setPixmap(daPixmap.scaled(arwWidth,arwHeight,Qt::KeepAspectRatio));Lda2->setMask(daPixmap.mask());
	// Long Left Arrow
	QLabel* Llla=new QLabel;
	Llla->setFixedHeight(2*arwHeight);
	Llla->setFixedWidth(2*imgWidth);
	Llla->setAlignment(Qt::AlignTop|Qt::AlignCenter);Llla->setPixmap(llaPixmap.scaled(2*imgWidth,2*arwHeight,Qt::KeepAspectRatio));Llla->setMask(llaPixmap.mask());
	// Hidden Label
	QLabel* Lh=new QLabel;QLabel* Lh2=new QLabel;QLabel* Lh3=new QLabel;QLabel* Lh4=new QLabel;QLabel* Lh5=new QLabel;
	Lh->setFixedHeight(arwHeight);
	Lh->setFixedWidth(imgWidth);
	Lh->setAlignment(Qt::AlignCenter|Qt::AlignCenter);//Lh->hide();
	Lh2->setFixedHeight(arwHeight);
	Lh2->setFixedWidth(arwWidth);
	Lh2->setAlignment(Qt::AlignCenter|Qt::AlignCenter);//Lh2->hide();
	Lh3->setFixedHeight(arwHeight);
	Lh3->setFixedWidth(arwWidth);
	Lh3->setAlignment(Qt::AlignCenter|Qt::AlignCenter);//Lh3->hide();
	Lh4->setFixedHeight(arwHeight);
	Lh4->setFixedWidth(arwWidth);
	Lh4->setAlignment(Qt::AlignCenter|Qt::AlignCenter);//Lh4->hide();
	Lh5->setFixedHeight(arwHeight);
	Lh5->setFixedWidth(imgWidth);
	Lh5->setAlignment(Qt::AlignTop|Qt::AlignCenter);//Lh5->hide();

	QGridLayout *inLo = new QGridLayout;
	//inLo->setSizeConstraint(QLayout::SetFixedSize);
	inLo->addWidget(s1Box,0,0,1,1);inLo->addWidget(Lra,0,1,1,1);inLo->addWidget(s2Box,0,2,1,1);inLo->addWidget(Lra2,0,3,1,1);inLo->addWidget(s3Box,0,4,1,1);
	inLo->addWidget(Lh,1,0,1,1);inLo->addWidget(Lh2,1,1,1,1);inLo->addWidget(Lda,1,2,1,1);inLo->addWidget(Lh3,1,3,1,1);inLo->addWidget(Lda2,1,4,1,1);
	inLo->addWidget(s6Box,2,0,1,1);inLo->addWidget(Lla,2,1,1,1);inLo->addWidget(s5Box,2,2,1,1);inLo->addWidget(Lh4,2,3,1,1);inLo->addWidget(s4Box,2,4,1,1);
	inLo->addWidget(Lh5,3,0,1,1);inLo->addWidget(Llla,3,1,1,3);
	sBox->setLayout(inLo);
	sBox->setStyleSheet("background-color:white;");
	//sBox->setMinimumWidth();

	// Scrolling Group Box
	QScrollArea* scrollArea=new QScrollArea();
	scrollArea->setWidget(sBox);  
	scrollArea->setWidgetResizable(true);  
	QGridLayout *outlo=new QGridLayout;
	//outlo->setSizeConstraint(QLayout::SetFixedSize);
	outlo->addWidget(scrollArea);
	QGroupBox* OutBox=new QGroupBox;
	OutBox->setLayout(outlo);
 
	QGridLayout *mLayout = new QGridLayout;
	//mLayout->setSizeConstraint(QLayout::SetFixedSize);
	//mLayout->addWidget(gOutBox);
	mLayout->addWidget(desc,0,0);
	mLayout->addWidget(OutBox,1,0);
	setLayout(mLayout);
	setWindowTitle(tr("ZPRED Computation Information"));
	resize(QDesktopWidget().availableGeometry(this).size());
}

solventInfo::solventInfo(QWidget *parent) : QWidget(parent)
{	
	// Get Current File Path and Find Images
	QString theCurFldr=QDir::currentPath();
	string Fldr=theCurFldr.toStdString()+"/";
	// Images
	// Acetic Acid
	string aaDielFile=Fldr+"UI/Images/Acetic_Acid_Dielectric.png";
	string aaDensFile=Fldr+"UI/Images/Acetic_Acid_Density.png";
	string aaViscFile=Fldr+"UI/Images/Acetic_Acid_Viscosity.png";
	QPixmap aaDielPixmap(aaDielFile.c_str());
	QPixmap aaDensPixmap(aaDensFile.c_str());
	QPixmap aaViscPixmap(aaViscFile.c_str());
	string aceticAcid="<b><span style=\"background-color:white; color:black; font-size:12pt;\">Acetic Acid<br><br></span></b>";
	aceticAcid+="<b><span style=\"background-color:white; color:black; font-size:12pt;\">&nbsp;&nbsp;H<br></span></b>";
	aceticAcid+="<b><span style=\"background-color:white; color:black; font-size:4pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\\<br></span></b>";
	aceticAcid+="<b><span style=\"background-color:white; color:black; font-size:12pt;\">H-C-C=O<br></span></b>";
	aceticAcid+="<b><span style=\"background-color:white; color:black; font-size:4pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;/&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\\<br></span></b>";
	aceticAcid+="<b><span style=\"background-color:white; color:black; font-size:12pt;\">&nbsp;&nbsp;H&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;O-H</span></b>";
	// Acetone
	string acDielFile=Fldr+"UI/Images/Acetone_Dielectric.png";
	string acDensFile=Fldr+"UI/Images/Acetone_Density.png";
	string acViscFile=Fldr+"UI/Images/Acetone_Viscosity.png";
	QPixmap acDielPixmap(acDielFile.c_str());
	QPixmap acDensPixmap(acDensFile.c_str());
	QPixmap acViscPixmap(acViscFile.c_str());
	string acetone="<b><span style=\"background-color:grey; color:black; font-size:12pt;\">Acetone<br><br></span></b>";
	acetone+="<b><span style=\"background-color:grey; color:black; font-size:12pt;\">&nbsp;&nbsp;H&nbsp;H<br></span></b>";
	acetone+="<b><span style=\"background-color:grey; color:black; font-size:4pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|/<br></span></b>";
	acetone+="<b><span style=\"background-color:grey; color:black; font-size:12pt;\">H-C<br></span></b>";
	acetone+="<b><span style=\"background-color:grey; color:black; font-size:4pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\\<br></span></b>";
	acetone+="<b><span style=\"background-color:grey; color:black; font-size:12pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;C=O<br></span></b>";
	acetone+="<b><span style=\"background-color:grey; color:black; font-size:4pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;/<br></span></b>";
	acetone+="<b><span style=\"background-color:grey; color:black; font-size:12pt;\">H-C<br></span></b>";
	acetone+="<b><span style=\"background-color:grey; color:black; font-size:4pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|\\<br></span></b>";
	acetone+="<b><span style=\"background-color:grey; color:black; font-size:12pt;\">&nbsp;&nbsp;H&nbsp;H</span></b>";
	// Butanol
	string bDielFile=Fldr+"UI/Images/Butanol_Dielectric.png";
	string bDensFile=Fldr+"UI/Images/Butanol_Density.png";
	string bViscFile=Fldr+"UI/Images/Butanol_Viscosity.png";
	QPixmap bDielPixmap(bDielFile.c_str());
	QPixmap bDensPixmap(bDensFile.c_str());
	QPixmap bViscPixmap(bViscFile.c_str());
	string butanol="<b><span style=\"background-color:white; color:black; font-size:12pt;\">Butanol<br><br></span></b>";
	butanol+="<b><span style=\"background-color:white; color:black; font-size:12pt;\">&nbsp;&nbsp;H&nbsp;&nbsp;&nbsp;H&nbsp;&nbsp;H&nbsp;H<br></span></b>";
	butanol+="<b><span style=\"background-color:white; color:black; font-size:4pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\\&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|<br></span></b>";
	butanol+="<b><span style=\"background-color:white; color:black; font-size:12pt;\">H-C-C-C-C-O-H<br></span></b>";
	butanol+="<b><span style=\"background-color:white; color:black; font-size:4pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;/&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|<br></span></b>";
	butanol+="<b><span style=\"background-color:white; color:black; font-size:12pt;\">&nbsp;&nbsp;H&nbsp;&nbsp;&nbsp;H&nbsp;&nbsp;H&nbsp;H</span></b>";
	// DMSO
	string dDielFile=Fldr+"UI/Images/DMSO_Dielectric.png";
	string dDensFile=Fldr+"UI/Images/DMSO_Density.png";
	string dViscFile=Fldr+"UI/Images/DMSO_Viscosity.png";
	QPixmap dDielPixmap(dDielFile.c_str());
	QPixmap dDensPixmap(dDensFile.c_str());
	QPixmap dViscPixmap(dViscFile.c_str());
	string dmso="<b><span style=\"background-color:grey; color:black; font-size:12pt;\">Dimethyl Sulfoxide<br><br></span></b>";
	dmso+="<b><span style=\"background-color:grey; color:black; font-size:12pt;\">&nbsp;&nbsp;H&nbsp;H<br></span></b>";
	dmso+="<b><span style=\"background-color:grey; color:black; font-size:4pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|/<br></span></b>";
	dmso+="<b><span style=\"background-color:grey; color:black; font-size:12pt;\">H-C<br></span></b>";
	dmso+="<b><span style=\"background-color:grey; color:black; font-size:4pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\\<br></span></b>";
	dmso+="<b><span style=\"background-color:grey; color:black; font-size:12pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;S=O<br></span></b>";
	dmso+="<b><span style=\"background-color:grey; color:black; font-size:4pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;/<br></span></b>";
	dmso+="<b><span style=\"background-color:grey; color:black; font-size:12pt;\">H-C<br></span></b>";
	dmso+="<b><span style=\"background-color:grey; color:black; font-size:4pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|\\<br></span></b>";
	dmso+="<b><span style=\"background-color:grey; color:black; font-size:12pt;\">&nbsp;&nbsp;H&nbsp;H</span></b>";
	// Ethanol
	string eDielFile=Fldr+"UI/Images/Ethanol_Dielectric.png";
	string eDensFile=Fldr+"UI/Images/Ethanol_Density.png";
	string eViscFile=Fldr+"UI/Images/Ethanol_Viscosity.png";
	QPixmap eDielPixmap(eDielFile.c_str());
	QPixmap eDensPixmap(eDensFile.c_str());
	QPixmap eViscPixmap(eViscFile.c_str());
	string ethanol="<b><span style=\"background-color:white; color:black; font-size:12pt;\">Ethanol<br><br></span></b>";
	ethanol+="<b><span style=\"background-color:white; color:black; font-size:12pt;\">&nbsp;&nbsp;H&nbsp;&nbsp;&nbsp;H<br></span></b>";
	ethanol+="<b><span style=\"background-color:white; color:black; font-size:4pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\\&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|<br></span></b>";
	ethanol+="<b><span style=\"background-color:white; color:black; font-size:12pt;\">H-C-C-O-H<br></span></b>";
	ethanol+="<b><span style=\"background-color:white; color:black; font-size:4pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;/&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|<br></span></b>";
	ethanol+="<b><span style=\"background-color:white; color:black; font-size:12pt;\">&nbsp;&nbsp;H&nbsp;&nbsp;&nbsp;H</span></b>";
	// Glycerol
	string gDielFile=Fldr+"UI/Images/Glycerol_Dielectric.png";
	string gDensFile=Fldr+"UI/Images/Glycerol_Density.png";
	string gViscFile=Fldr+"UI/Images/Glycerol_Viscosity.png";
	QPixmap gDielPixmap(gDielFile.c_str());
	QPixmap gDensPixmap(gDensFile.c_str());
	QPixmap gViscPixmap(gViscFile.c_str());
	string glycerol="<b><span style=\"background-color:grey; color:black; font-size:12pt;\">Glycerol<br><br></span></b>";
	glycerol+="<b><span style=\"background-color:grey; color:black; font-size:12pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;O-H<br></span></b>";
	glycerol+="<b><span style=\"background-color:grey; color:black; font-size:4pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|<br></span></b>";
	glycerol+="<b><span style=\"background-color:grey; color:black; font-size:12pt;\">H-O&nbsp;HC&nbsp;&nbsp;&nbsp;O-H<br></span></b>";
	glycerol+="<b><span style=\"background-color:grey; color:black; font-size:4pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\\&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;/&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\\&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;/<br></span></b>";
	glycerol+="<b><span style=\"background-color:grey; color:black; font-size:12pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;C&nbsp;&nbsp;&nbsp;&nbsp;C<br></span></b>";
	glycerol+="<b><span style=\"background-color:grey; color:black; font-size:4pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;/&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\\&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;/&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\\<br></span></b>";
	glycerol+="<b><span style=\"background-color:grey; color:black; font-size:12pt;\">&nbsp;&nbsp;&nbsp;&nbsp;H&nbsp;&nbsp;H&nbsp;H&nbsp;H<br></span></b>";
	// Isopropanol
	string iDielFile=Fldr+"UI/Images/Isopropanol_Dielectric.png";
	string iDensFile=Fldr+"UI/Images/Isopropanol_Density.png";
	string iViscFile=Fldr+"UI/Images/Isopropanol_Viscosity.png";
	QPixmap iDielPixmap(iDielFile.c_str());
	QPixmap iDensPixmap(iDensFile.c_str());
	QPixmap iViscPixmap(iViscFile.c_str());
	string isopropanol="<b><span style=\"background-color:white; color:black; font-size:12pt;\">Isopropanol<br><br></span></b>";
	isopropanol+="<b><span style=\"background-color:white; color:black; font-size:12pt;\">&nbsp;&nbsp;H&nbsp;H<br></span></b>";
	isopropanol+="<b><span style=\"background-color:white; color:black; font-size:4pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|/<br></span></b>";
	isopropanol+="<b><span style=\"background-color:white; color:black; font-size:12pt;\">H-C<br></span></b>";
	isopropanol+="<b><span style=\"background-color:white; color:black; font-size:4pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\\<br></span></b>";
	isopropanol+="<b><span style=\"background-color:white; color:black; font-size:12pt;\">H-C-O-H<br></span></b>";
	isopropanol+="<b><span style=\"background-color:white; color:black; font-size:4pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;/<br></span></b>";
	isopropanol+="<b><span style=\"background-color:white; color:black; font-size:12pt;\">H-C<br></span></b>";
	isopropanol+="<b><span style=\"background-color:white; color:black; font-size:4pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|\\<br></span></b>";
	isopropanol+="<b><span style=\"background-color:white; color:black; font-size:12pt;\">&nbsp;&nbsp;H&nbsp;H</span></b>";
	// Methanol
	string mDielFile=Fldr+"UI/Images/Methanol_Dielectric.png";
	string mDensFile=Fldr+"UI/Images/Methanol_Density.png";
	string mViscFile=Fldr+"UI/Images/Methanol_Viscosity.png";
	QPixmap mDielPixmap(mDielFile.c_str());
	QPixmap mDensPixmap(mDensFile.c_str());
	QPixmap mViscPixmap(mViscFile.c_str());
	string methanol="<b><span style=\"background-color:grey; color:black; font-size:12pt;\">Methanol<br><br></span></b>";
	methanol+="<b><span style=\"background-color:grey; color:black; font-size:12pt;\">&nbsp;&nbsp;H<br></span></b>";
	methanol+="<b><span style=\"background-color:grey; color:black; font-size:4pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\\<br></span></b>";
	methanol+="<b><span style=\"background-color:grey; color:black; font-size:12pt;\">H-C-O-H<br></span></b>";
	methanol+="<b><span style=\"background-color:grey; color:black; font-size:4pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;/<br></span></b>";
	methanol+="<b><span style=\"background-color:grey; color:black; font-size:12pt;\">&nbsp;&nbsp;H</span></b>";
	// Propanol
	string pDielFile=Fldr+"UI/Images/Propanol_Dielectric.png";
	string pDensFile=Fldr+"UI/Images/Propanol_Density.png";
	string pViscFile=Fldr+"UI/Images/Propanol_Viscosity.png";
	QPixmap pDielPixmap(pDielFile.c_str());
	QPixmap pDensPixmap(pDensFile.c_str());
	QPixmap pViscPixmap(pViscFile.c_str());
	string propanol="<b><span style=\"background-color:white; color:black; font-size:12pt;\">Propanol<br><br></span></b>";
	propanol+="<b><span style=\"background-color:white; color:black; font-size:12pt;\">&nbsp;&nbsp;H&nbsp;&nbsp;&nbsp;H&nbsp;&nbsp;H<br></span></b>";
	propanol+="<b><span style=\"background-color:white; color:black; font-size:4pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\\&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|<br></span></b>";
	propanol+="<b><span style=\"background-color:white; color:black; font-size:12pt;\">H-C-C-C-O-H<br></span></b>";
	propanol+="<b><span style=\"background-color:white; color:black; font-size:4pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;/&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|<br></span></b>";
	propanol+="<b><span style=\"background-color:white; color:black; font-size:12pt;\">&nbsp;&nbsp;H&nbsp;&nbsp;&nbsp;H&nbsp;&nbsp;H</span></b>";
	// Water
	string wDielFile=Fldr+"UI/Images/Water_Dielectric.png";
	string wDensFile=Fldr+"UI/Images/Water_Density.png";
	string wDensEqnFile=Fldr+"UI/Images/Water_Density_Eqn.png";
	string wViscFile=Fldr+"UI/Images/Water_Viscosity.png";
	QPixmap wDielPixmap(wDielFile.c_str());
	QPixmap wDensPixmap(wDensFile.c_str());
	QPixmap wDensEqnPixmap(wDensEqnFile.c_str());
	QPixmap wViscPixmap(wViscFile.c_str());
	string water="<b><span style=\"background-color:grey; color:black; font-size:12pt;\">Water<br><br></span></b>";
	water+="<b><span style=\"background-color:grey; color:black; font-size:12pt;\">&nbsp;&nbsp;H<br></span></b>";
	water+="<b><span style=\"background-color:grey; color:black; font-size:4pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\\<br></span></b>";
	water+="<b><span style=\"background-color:grey; color:black; font-size:12pt;\">&nbsp;&nbsp;&nbsp;&nbsp;O<br></span></b>";
	water+="<b><span style=\"background-color:grey; color:black; font-size:4pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;/<br></span></b>";
	water+="<b><span style=\"background-color:grey; color:black; font-size:12pt;\">&nbsp;&nbsp;H</span></b>";

	int propertyBoxSz=1150;
	QGroupBox* gBox = new QGroupBox;
	//gBox->setStyleSheet("QGroupBox { font-size:12pt; font-weight:bold; color:blue; } ");
	//gBox->setStyleSheet("QGroupBox { font-weight:bold; } ");
	QFont font;
	font.setBold(true);
	//gBox->setFont(font);//gBox->setAlignment(Qt::AlignRight);
	// Initial Description
	QLabel* description=new QLabel();
	string s="<b>This interface provides empirical relationships for calculating the solution properties (specifically, relative dielectric, density, Debye length and viscosity) used by ZPRED. In general, all properties are defined by first calculating the respective property of pure solvent (e.g. water) at a specified temperature and then correcting the value to consider the presence of solutes. As this work validating ZPRED focused on aqueous solutions, the discussion is primarily directed to accounting for solutes in water; however, the mixing rules presented may be useful for mixed solvents as well.</b><br>";
	//s=setStringWidth(s,width);
	description->setText(tr(s.c_str()));
	description->setWordWrap(true);
	description->setFixedWidth(propertyBoxSz);

	// Relative Dielectric Description
	QLabel* erDesc=new QLabel;
	s="The <b>relative dielectric</b> (also called <b>static relative permittivity</b>) (<b>ε<sub>r</sub></b>) represents a substance's ability to reduce the electric charge force between two approaching particles submerged in the substance. To compute the relative dielectric of a solution, it is first necessary to predict the relative dielectric of the pure solvent(s) (<b>ε<sub>rs</sub></b>) used. Pure solvent relative dielectric constants as a function of temperature are shown below for water and a few other common solvents. Plots in the table show the fit of the empirical relations to large sets of experimental data and provide a comparison among solvents in the relevant temperature range for water.";
	//s=setStringWidth(s,width);
	erDesc->setText(tr(s.c_str()));
	erDesc->setWordWrap(true);
	erDesc->setFixedWidth(propertyBoxSz);

	// Specify External Executables
	QGroupBox* erInBox=new QGroupBox;
	QGroupBox* erInBox1=new QGroupBox;
	QGroupBox* erInBox2=new QGroupBox;
	QGroupBox* erInBox3=new QGroupBox;
	QGroupBox* erInBox4=new QGroupBox;
	QGroupBox* erInBox5=new QGroupBox;
	QGroupBox* erInBox6=new QGroupBox;
	QGroupBox* erInBox7=new QGroupBox;
	QGroupBox* erInBox8=new QGroupBox;
	QGroupBox* erInBox9=new QGroupBox;
	QGroupBox* erInBox10=new QGroupBox;
	QGroupBox* erInBox11=new QGroupBox;

	//erInBox->setFont(font);	
	QLabel* L11=new QLabel;QLabel* L12=new QLabel;QLabel* L13=new QLabel;
	L11->setText("<b><span style=\"background-color:darkgrey; color:black; font-size:18pt;\">Solvent</span></b>");L11->setAlignment(Qt::AlignCenter);	
	L12->setText("<b><span style=\"background-color:darkgrey; color:black; font-size:18pt;\">ε<sub>rs</sub> as Function of Absolute Temperature(t)</span></b>");L12->setAlignment(Qt::AlignCenter);
	L13->setText("<span style=\"background-color:darkgrey; color:black; font-size:18pt;\"> </span>");
	//
	QGridLayout *lo1=new QGridLayout;
	//lo1->setSizeConstraint(QLayout::SetFixedSize);
	lo1->addWidget(L11,0,0,1,1);lo1->addWidget(L12,0,1,1,1);lo1->addWidget(L13,0,2,1,1);
	erInBox1->setLayout(lo1);
	erInBox1->setStyleSheet("background-color:darkgrey;");
	// Acetic Acid
	int imgWidth=400;
	int imgHeight=400;
	QLabel* L21=new QLabel;QLabel* L22=new QLabel;QLabel* L23=new QLabel;
	L21->setText(aceticAcid.c_str());//L21->setAlignment(Qt::AlignCenter);	
	L22->setText("<b><span style=\"background-color:white; color:black; font-size:12pt;\">-15.731 + 0.12662t − 0.00017738t<sup>2</sup></span></b>");L22->setAlignment(Qt::AlignCenter);
	//L23->setFixedHeight(imgHeight);
	L23->setFixedWidth(imgWidth);
	L23->setAlignment(Qt::AlignCenter);L23->setPixmap(aaDielPixmap.scaled(imgWidth,imgHeight,Qt::KeepAspectRatio));L23->setMask(aaDielPixmap.mask());
	QGridLayout *lo2=new QGridLayout;
	//lo2->setSizeConstraint(QLayout::SetFixedSize);
	lo2->addWidget(L21,0,0,1,1);lo2->addWidget(L22,0,1,1,1);lo2->addWidget(L23,0,2,1,1);
	erInBox2->setLayout(lo2);
	erInBox2->setStyleSheet("background-color:white;");
	// Acetone
	QLabel* L31=new QLabel;QLabel* L32=new QLabel;QLabel* L33=new QLabel;
	L31->setText(acetone.c_str());//L31->setAlignment(Qt::AlignCenter);	
	L32->setText("<b><span style=\"background-color:grey; color:black; font-size:12pt;\">88.157 − 0.343t + 0.00038925t<sup>2</sup></span></b>");L32->setAlignment(Qt::AlignCenter);
	//L33->setFixedHeight(imgHeight);
	L33->setFixedWidth(imgWidth);
	L33->setAlignment(Qt::AlignCenter);L33->setPixmap(acDielPixmap.scaled(imgWidth,imgHeight,Qt::KeepAspectRatio));L33->setMask(acDielPixmap.mask());
	QGridLayout *lo3=new QGridLayout;
	lo3->addWidget(L31,0,0,1,1);lo3->addWidget(L32,0,1,1,1);lo3->addWidget(L33,0,2,1,1);
	erInBox3->setLayout(lo3);
	erInBox3->setStyleSheet("background-color:grey;");
	// Butanol
	QLabel* L41=new QLabel;QLabel* L42=new QLabel;QLabel* L43=new QLabel;
	L41->setText(butanol.c_str());//L41->setAlignment(Qt::AlignCenter);	
	L42->setText("<b><span style=\"background-color:white; color:black; font-size:12pt;\">105.78 − 0.50587t + 0.00084733t<sup>2</sup> − 4.8841x10<sup>-7</sup>t<sup>3</sup></span></b>");L42->setAlignment(Qt::AlignCenter);
	//L43->setFixedHeight(imgHeight);
	L43->setFixedWidth(imgWidth);
	L43->setAlignment(Qt::AlignCenter);L43->setPixmap(bDielPixmap.scaled(imgWidth,imgHeight,Qt::KeepAspectRatio));L43->setMask(bDielPixmap.mask());
	QGridLayout *lo4=new QGridLayout;
	lo4->addWidget(L41,0,0,1,1);lo4->addWidget(L42,0,1,1,1);lo4->addWidget(L43,0,2,1,1);
	erInBox4->setLayout(lo4);
	erInBox4->setStyleSheet("background-color:white;");
	// DMSO
	QLabel* L51=new QLabel;QLabel* L52=new QLabel;QLabel* L53=new QLabel;
	L51->setText(dmso.c_str());//L51->setAlignment(Qt::AlignCenter);	
	L52->setText("<b><span style=\"background-color:grey; color:black; font-size:12pt;\">38.478 + 0.16939t − 0.00047423t<sup>2</sup></span></b>");L52->setAlignment(Qt::AlignCenter);
	//L53->setFixedHeight(imgHeight);
	L53->setFixedWidth(imgWidth);
	L53->setAlignment(Qt::AlignCenter);L53->setPixmap(dDielPixmap.scaled(imgWidth,imgHeight,Qt::KeepAspectRatio));L53->setMask(dDielPixmap.mask());
	QGridLayout *lo5=new QGridLayout;
	lo5->addWidget(L51,0,0,1,1);lo5->addWidget(L52,0,1,1,1);lo5->addWidget(L53,0,2,1,1);
	erInBox5->setLayout(lo5);
	erInBox5->setStyleSheet("background-color:grey;");
	// Ethanol
	QLabel* L61=new QLabel;QLabel* L62=new QLabel;QLabel* L63=new QLabel;	
	L61->setText(ethanol.c_str());//L61->setAlignment(Qt::AlignCenter);	
	L62->setText("<b><span style=\"background-color:white; color:black; font-size:12pt;\">151.45 − 0.87020t + 0.001957t<sup>2</sup> − 0.0000015512t<sup>3</sup></span></b>");L62->setAlignment(Qt::AlignCenter);
	L63->setFixedWidth(imgWidth);
	L63->setAlignment(Qt::AlignCenter);L63->setPixmap(eDielPixmap.scaled(imgWidth,imgHeight,Qt::KeepAspectRatio));L63->setMask(eDielPixmap.mask());
	QGridLayout *lo6=new QGridLayout;
	//lo6->setSizeConstraint(QLayout::SetFixedSize);
	lo6->addWidget(L61,0,0,1,1);lo6->addWidget(L62,0,1,1,1);lo6->addWidget(L63,0,2,1,1);
	erInBox6->setLayout(lo6);
	erInBox6->setStyleSheet("background-color:white;");
	// Glycerol
	QLabel* L71=new QLabel;QLabel* L72=new QLabel;QLabel* L73=new QLabel;
	L71->setText(glycerol.c_str());//L71->setAlignment(Qt::AlignCenter);	
	L72->setText("<b><span style=\"background-color:grey; color:black; font-size:12pt;\">77.503 − 0.037984t − 0.00023107t<sup>2</sup></span></b>");L72->setAlignment(Qt::AlignCenter);
	//L73->setFixedHeight(imgHeight);
	L73->setFixedWidth(imgWidth);
	L73->setAlignment(Qt::AlignCenter);L73->setPixmap(gDielPixmap.scaled(imgWidth,imgHeight,Qt::KeepAspectRatio));L73->setMask(gDielPixmap.mask());
	QGridLayout *lo7=new QGridLayout;
	//lo7->setSizeConstraint(QLayout::SetFixedSize);
	lo7->addWidget(L71,0,0,1,1);lo7->addWidget(L72,0,1,1,1);lo7->addWidget(L73,0,2,1,1);
	erInBox7->setLayout(lo7);
	erInBox7->setStyleSheet("background-color:grey;");
	// Isopropanol
	QLabel* L81=new QLabel;QLabel* L82=new QLabel;QLabel* L83=new QLabel;
	L81->setText(isopropanol.c_str());//L81->setAlignment(Qt::AlignCenter);	
	L82->setText("<b><span style=\"background-color:white; color:black; font-size:12pt;\">104.16 − 0.41011t + 0.00042049t<sup>2</sup></span></b>");L82->setAlignment(Qt::AlignCenter);
	//L83->setFixedHeight(imgHeight);
	L83->setFixedWidth(imgWidth);
	L83->setAlignment(Qt::AlignCenter);L83->setPixmap(iDielPixmap.scaled(imgWidth,imgHeight,Qt::KeepAspectRatio));L83->setMask(iDielPixmap.mask());
	QGridLayout *lo8=new QGridLayout;
	lo8->addWidget(L81,0,0,1,1);lo8->addWidget(L82,0,1,1,1);lo8->addWidget(L83,0,2,1,1);
	erInBox8->setLayout(lo8);
	erInBox8->setStyleSheet("background-color:white;");
	// Methanol
	QLabel* L91=new QLabel;QLabel* L92=new QLabel;QLabel* L93=new QLabel;
	L91->setText(methanol.c_str());
	L92->setText("<b><span style=\"background-color:grey; color:black; font-size:12pt;\">193.41 − 0.92211t + 0.0012839t<sup>2</sup></span></b>");L92->setAlignment(Qt::AlignCenter);
	//L93->setFixedHeight(imgHeight);
	L93->setFixedWidth(imgWidth);
	L93->setAlignment(Qt::AlignCenter);L93->setPixmap(mDielPixmap.scaled(imgWidth,imgHeight,Qt::KeepAspectRatio));L93->setMask(mDielPixmap.mask());
	QGridLayout *lo9=new QGridLayout;
	lo9->addWidget(L91,0,0,1,1);lo9->addWidget(L92,0,1,1,1);lo9->addWidget(L93,0,2,1,1);
	erInBox9->setLayout(lo9);
	erInBox9->setStyleSheet("background-color:grey;");
	// Propanol
	QLabel* L101=new QLabel;QLabel* L102=new QLabel;QLabel* L103=new QLabel;
	L101->setText(propanol.c_str());
	L102->setText("<b><span style=\"background-color:white; color:black; font-size:12pt;\">98.045 − 0.36860t + 0.00036422t<sup>2</sup></span></b>");L102->setAlignment(Qt::AlignCenter);
	//L103->setFixedHeight(imgHeight);
	L103->setFixedWidth(imgWidth);
	L103->setAlignment(Qt::AlignCenter);L103->setPixmap(pDielPixmap.scaled(imgWidth,imgHeight,Qt::KeepAspectRatio));L103->setMask(pDielPixmap.mask());
	QGridLayout *lo10=new QGridLayout;
	lo10->addWidget(L101,0,0,1,1);lo10->addWidget(L102,0,1,1,1);lo10->addWidget(L103,0,2,1,1);
	erInBox10->setLayout(lo10);
	erInBox10->setStyleSheet("background-color:white;");
	// Water
	QLabel* L111=new QLabel;QLabel* L112=new QLabel;QLabel* L113=new QLabel;
	L111->setText(water.c_str());
	L112->setText("<b><span style=\"background-color:grey; color:black; font-size:12pt;\">87.740 − 0.4008(t−273.15) + 0.0009398(t−273.15)<sup>2</sup> − 0.00000141(t−273.15)<sup>3</sup></span></b>");L112->setAlignment(Qt::AlignCenter);
	//L113->setFixedHeight(imgHeight);
	L113->setFixedWidth(imgWidth);
	L113->setAlignment(Qt::AlignCenter);L113->setPixmap(wDielPixmap.scaled(imgWidth,imgHeight,Qt::KeepAspectRatio));L113->setMask(wDielPixmap.mask());
	QGridLayout *lo11=new QGridLayout;
	lo11->addWidget(L111,0,0,1,1);lo11->addWidget(L112,0,1,1,1);lo11->addWidget(L113,0,2,1,1);
	erInBox11->setLayout(lo11);
	erInBox11->setStyleSheet("background-color:grey;");
	// Mixed Solvents
	QGroupBox* erInBox12=new QGroupBox(tr("Mixed Solvents"));//erInBox12->setFont(font);
	QLabel* L121=new QLabel;QLabel* L122=new QLabel;QLabel* L123=new QLabel;
	L121->setText("If multiple solvents are used, the following mixing rule (based on Oster's Rule)<b><sup>[5]</sup></b> can be applied to determine the relative dielectric of the pure solvent mixture (<b>ε<sub>rs,mix</sub></b>):<br>");
	L121->setWordWrap(true);
	L121->setFixedWidth(propertyBoxSz);
	s="<b><span style=\"color:black; font-size:12pt;\">ε<sub>rs,mix</sub> = <sup>1</sup>/<sub>4</sub>[1 + 9P<sub>m</sub> + 3(9P<sub>m</sub><sup>2</sup>+2P<sub>m</sub>+1)<sup>(1/2)</sup>]<br><br></span></b>";
	s+="<b><span style=\"color:black; font-size:12pt;\">P<sub>m</sub> = &Sigma;<sup>N</sup><sub>i=1</sub>x<sub>i</sub>V<sub>m,i</sub>P<sub>i</sub><br></span></b>";
	s+="<b><span style=\"color:black; font-size:4pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------------------------------------------------<br></span></b>";
	s+="<b><span style=\"color:black; font-size:12pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&Sigma;<sup>N</sup><sub>i=1</sub>x<sub>i</sub>V<sub>m,i</sub><br><br></span></b>";
	s+="<b><span style=\"color:black; font-size:12pt;\">P<sub>i</sub> = (ε<sub>rs,i</sub>-1)(2ε<sub>rs,i</sub>+1)<br></span></b>";
	s+="<b><span style=\"color:black; font-size:4pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--------------------------------------------------------------<br></span></b>";
	s+="<b><span style=\"color:black; font-size:12pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;9ε<sub>rs,i</sub><br></span></b>";
	L122->setText(s.c_str());
	s="where <b>P<sub>m</sub></b> is the polarization of the solvent mixture, <b>N</b> is the number of solvents, <b>x<sub>i</sub></b> is the mole fraction of the i-th solvent, <b>V<sub>m,i</sub></b> is the molar volume of the i-th solvent (determined by dividing the solvent molecular weight by its density), <b>P<sub>i</sub></b> is the polarization of the i-th solvent and <b>ε<sub>rs,i</sub></b> is the relative dielectric of i-th solvent (defined above)";
	L123->setText(s.c_str());
	L123->setWordWrap(true);
	L123->setFixedWidth(propertyBoxSz);
	QGridLayout *lo12=new QGridLayout;
	lo12->addWidget(L121,0,0);lo12->addWidget(L122,1,0);lo12->addWidget(L123,2,0);
	erInBox12->setLayout(lo12);	
	// Accounting for Solutes
	QGroupBox* erInBox13=new QGroupBox(tr("Accounting for Solutes"));
	QLabel* L131=new QLabel;QLabel* L132=new QLabel;QLabel* L133=new QLabel;
	L131->setText("The relative dielectric (<b>ε<sub>r</sub></b>) of a solution can either increase or decrease with addition of solutes depending on their interaction in the solvent. An increase in the dielectric can occur from ion pairing, which results in the ion pair holding a higher dipole moment than the solvent. On the other hand, the dielectric decreases when ions dissociate and this effect is referred to as the <b>dielectric decrement</b><b><sup>[6-9]</sup></b>. Under dilute ion concentrations (≪ 0.01M), the relative dielectric of the solution can be assumed to be equal to that of the pure solvent. However, with increasing ion concentration the relative dielectric can deviate significantly from that of the pure solvent and this must be accounted for.<br><br>Currently, ZPRED only accounts for the linear behavior of the dielectric decrement ([salt] < 1.5M) using the simple relation below<b><sup>[6][7]</sup></b>.<br><br>");
	L131->setWordWrap(true);
	L131->setFixedWidth(propertyBoxSz);
	s="<b><span style=\"color:black; font-size:12pt;\">ε<sub>r</sub> = ε<sub>rs</sub> − &alpha;C<br></span></b>";
	L132->setText(s.c_str());
	s="where <b>&alpha;</b> is an ion-specific parameter called the <b>total excess polarization<sup>[9]</sup></b> of the ionic species and <b>C</b> is the molar concentration of the ionic species of interest<br><br>This limited model did not affect validation of ZPRED as ion concentrations tested remained below the non-linear regime. This limitation of ZPRED has been presented with the intention to be corrected in future work. A good direction of advancement is referenced here<b><sup>[9]</sup></b>.";
	L133->setText(s.c_str());
	L133->setWordWrap(true);
	L133->setFixedWidth(propertyBoxSz);
	QGridLayout *lo13=new QGridLayout;
	lo13->addWidget(L131,0,0);lo13->addWidget(L132,1,0);lo13->addWidget(L133,2,0);
	erInBox13->setLayout(lo13);
	erInBox13->setStyleSheet("background-color:grey;");

	QGridLayout *inLo = new QGridLayout;
	//inLo->addWidget(erDesc,0,0,1,3);
	//inLo->addWidget(erInBox1,1,0,1,3);
	inLo->setSizeConstraint(QLayout::SetFixedSize);
	inLo->addWidget(erInBox2,0,0);
	inLo->addWidget(erInBox3,1,0);
	inLo->addWidget(erInBox4,2,0);
	inLo->addWidget(erInBox5,3,0);
	inLo->addWidget(erInBox6,4,0);
	inLo->addWidget(erInBox7,5,0);
	inLo->addWidget(erInBox8,6,0);
	inLo->addWidget(erInBox9,7,0);
	inLo->addWidget(erInBox10,8,0);
	inLo->addWidget(erInBox11,9,0);
	erInBox->setLayout(inLo);
	erInBox->setStyleSheet("background-color:black;");

	// Scrolling Group Box
	QScrollArea* scrollArea=new QScrollArea();
	scrollArea->setWidget(erInBox);  
	scrollArea->setWidgetResizable(true);  
	QGridLayout *outlo=new QGridLayout;
	outlo->setSizeConstraint(QLayout::SetFixedSize);
	outlo->addWidget(erDesc,0,0);
	outlo->addWidget(erInBox1,1,0);
	outlo->addWidget(scrollArea,2,0);
	outlo->addWidget(erInBox12,3,0);
	outlo->addWidget(erInBox13,4,0);
	QGroupBox* erOutBox=new QGroupBox(tr("|-RELATIVE DIELECTRIC-|"));
	//erOutBox->setFont(font);
	erOutBox->setLayout(outlo);

	// Density
	// Description
	QLabel* dDesc=new QLabel;
	s="The <b>density</b> (specifically <b>mass density</b>) of a solution is a measure of how condensed its mass is within the volume it occupies. To compute the density of a solution, it is first necessary to predict the density of the pure solvent(s) (<b>ρ<sub>s</sub></b>) used. Pure solvent densities as a function of temperature are shown below for water and a few other common solvents. Plots in the table show the fit of the empirical relations to large sets of experimental data and provide a comparison among solvents in the relevant temperature range for water.";
	//s=setStringWidth(s,width);
	dDesc->setText(tr(s.c_str()));
	dDesc->setWordWrap(true);
	dDesc->setFixedWidth(propertyBoxSz);

	// Specify External Executables
	QGroupBox* dInBox=new QGroupBox;
	QGroupBox* dInBox1=new QGroupBox;
	QGroupBox* dInBox2=new QGroupBox;
	QGroupBox* dInBox3=new QGroupBox;
	QGroupBox* dInBox4=new QGroupBox;
	QGroupBox* dInBox5=new QGroupBox;
	QGroupBox* dInBox6=new QGroupBox;
	QGroupBox* dInBox7=new QGroupBox;
	QGroupBox* dInBox8=new QGroupBox;
	QGroupBox* dInBox9=new QGroupBox;
	QGroupBox* dInBox10=new QGroupBox;
	QGroupBox* dInBox11=new QGroupBox;

	QLabel* dL11=new QLabel;QLabel* dL12=new QLabel;QLabel* dL13=new QLabel;
	dL11->setText("<b><span style=\"background-color:darkgrey; color:black; font-size:18pt;\">Solvent</span></b>");dL11->setAlignment(Qt::AlignCenter);	
	dL12->setText("<b><span style=\"background-color:darkgrey; color:black; font-size:18pt;\">ρ<sub>s</sub> [kg/L] as Function of Absolute Temperature(t)</span></b>");dL12->setAlignment(Qt::AlignCenter);
	dL13->setText("<span style=\"background-color:darkgrey; color:black; font-size:18pt;\"> </span>");
	//
	QGridLayout *dlo1=new QGridLayout;
	//dlo1->setSizeConstraint(QLayout::SetFixedSize);
	dlo1->addWidget(dL11,0,0,1,1);dlo1->addWidget(dL12,0,1,1,1);dlo1->addWidget(dL13,0,2,1,1);
	dInBox1->setLayout(dlo1);
	dInBox1->setStyleSheet("background-color:darkgrey;");
	// Acetic Acid
	QLabel* dL21=new QLabel;QLabel* dL22=new QLabel;QLabel* dL23=new QLabel;
	dL21->setText(aceticAcid.c_str());//dL21->setAlignment(Qt::AlignCenter);	
	dL22->setText("<b><span style=\"background-color:white; color:black; font-size:12pt;\">1.27334 − 4.77633x10<sup>-4</sup>t − 9.8132x10<sup>-7</sup>t<sup>2</sup></span></b>");dL22->setAlignment(Qt::AlignCenter);
	//dL23->setFixedHeight(imgHeight);
	dL23->setFixedWidth(imgWidth);
	dL23->setAlignment(Qt::AlignCenter);dL23->setPixmap(aaDensPixmap.scaled(imgWidth,imgHeight,Qt::KeepAspectRatio));dL23->setMask(aaDensPixmap.mask());
	QGridLayout *dlo2=new QGridLayout;
	dlo2->addWidget(dL21,0,0,1,1);dlo2->addWidget(dL22,0,1,1,1);dlo2->addWidget(dL23,0,2,1,1);
	dInBox2->setLayout(dlo2);
	dInBox2->setStyleSheet("background-color:white;");
	// Acetone
	QLabel* dL31=new QLabel;QLabel* dL32=new QLabel;QLabel* dL33=new QLabel;
	dL31->setText(acetone.c_str());//dL31->setAlignment(Qt::AlignCenter);	
	dL32->setText("<b><span style=\"background-color:grey; color:black; font-size:12pt;\">1.06011 − 6.68357x10<sup>-4</sup>t − 8.60506x10<sup>-7</sup>t<sup>2</sup></span></b>");dL32->setAlignment(Qt::AlignCenter);
	//dL33->setFixedHeight(imgHeight);
	dL33->setFixedWidth(imgWidth);
	dL33->setAlignment(Qt::AlignCenter);dL33->setPixmap(acDensPixmap.scaled(imgWidth,imgHeight,Qt::KeepAspectRatio));dL33->setMask(acDensPixmap.mask());
	QGridLayout *dlo3=new QGridLayout;
	dlo3->addWidget(dL31,0,0,1,1);dlo3->addWidget(dL32,0,1,1,1);dlo3->addWidget(dL33,0,2,1,1);
	dInBox3->setLayout(dlo3);
	dInBox3->setStyleSheet("background-color:grey;");
	// Butanol
	QLabel* dL41=new QLabel;QLabel* dL42=new QLabel;QLabel* dL43=new QLabel;
	dL41->setText(butanol.c_str());//dL41->setAlignment(Qt::AlignCenter);	
	dL42->setText("<b><span style=\"background-color:white; color:black; font-size:12pt;\">1.15309 − 2.13475x10<sup>-3</sup>t + 5.15573x10<sup>-6</sup>t<sup>2</sup> − 6.38112x10<sup>-9</sup>t<sup>3</sup></span></b>");dL42->setAlignment(Qt::AlignCenter);
	//dL43->setFixedHeight(imgHeight);
	dL43->setFixedWidth(imgWidth);
	dL43->setAlignment(Qt::AlignCenter);dL43->setPixmap(bDensPixmap.scaled(imgWidth,imgHeight,Qt::KeepAspectRatio));dL43->setMask(bDensPixmap.mask());
	QGridLayout *dlo4=new QGridLayout;
	dlo4->addWidget(dL41,0,0,1,1);dlo4->addWidget(dL42,0,1,1,1);dlo4->addWidget(dL43,0,2,1,1);
	dInBox4->setLayout(dlo4);
	dInBox4->setStyleSheet("background-color:white;");
	// DMSO
	QLabel* dL51=new QLabel;QLabel* dL52=new QLabel;QLabel* dL53=new QLabel;
	dL51->setText(dmso.c_str());//dL51->setAlignment(Qt::AlignCenter);	
	s="<b><span style=\"background-color:grey; color:black; font-size:12pt;\">0.3442[(0.2534)<sup>-y</sup>]<br><br>y = (1 - t/726)<sup>0.322</sup></span></b>";
	dL52->setText(s.c_str());dL52->setAlignment(Qt::AlignCenter);
	dL53->setFixedWidth(imgWidth);
	dL53->setAlignment(Qt::AlignCenter);dL53->setPixmap(dDensPixmap.scaled(imgWidth,imgHeight,Qt::KeepAspectRatio));dL53->setMask(dDensPixmap.mask());
	QGridLayout *dlo5=new QGridLayout;
	dlo5->addWidget(dL51,0,0,1,1);dlo5->addWidget(dL52,0,1,1,1);dlo5->addWidget(dL53,0,2,1,1);
	dInBox5->setLayout(dlo5);
	dInBox5->setStyleSheet("background-color:grey;");
	// Ethanol
	QLabel* dL61=new QLabel;QLabel* dL62=new QLabel;QLabel* dL63=new QLabel;	
	dL61->setText(ethanol.c_str());//dL61->setAlignment(Qt::AlignCenter);
	dL62->setText("<b><span style=\"background-color:white; color:black; font-size:12pt;\">1.16239 − 2.25788x10<sup>-3</sup>t + 5.30621x10<sup>-6</sup>t<sup>2</sup> − 6.6307x10<sup>-9</sup>t<sup>3</sup></span></b>");dL62->setAlignment(Qt::AlignCenter);
	//dL63->setFixedHeight(imgHeight);
	dL63->setFixedWidth(imgWidth);
	dL63->setAlignment(Qt::AlignCenter);dL63->setPixmap(eDensPixmap.scaled(imgWidth,imgHeight,Qt::KeepAspectRatio));dL63->setMask(eDensPixmap.mask());
	QGridLayout *dlo6=new QGridLayout;
	dlo6->addWidget(dL61,0,0,1,1);dlo6->addWidget(dL62,0,1,1,1);dlo6->addWidget(dL63,0,2,1,1);
	dInBox6->setLayout(dlo6);
	dInBox6->setStyleSheet("background-color:white;");
	// Glycerol
	QLabel* dL71=new QLabel;QLabel* dL72=new QLabel;QLabel* dL73=new QLabel;
	dL71->setText(glycerol.c_str());//dL71->setAlignment(Qt::AlignCenter);	
	dL72->setText("<b><span style=\"background-color:grey; color:black; font-size:12pt;\">1.43632 − 5.28889x10<sup>-4</sup>t − 2.29556x10<sup>-7</sup>t<sup>2</sup></span></b>");dL72->setAlignment(Qt::AlignCenter);
	//dL73->setFixedHeight(imgHeight);
	dL73->setFixedWidth(imgWidth);
	dL73->setAlignment(Qt::AlignCenter);dL73->setPixmap(gDensPixmap.scaled(imgWidth,imgHeight,Qt::KeepAspectRatio));dL73->setMask(gDensPixmap.mask());
	QGridLayout *dlo7=new QGridLayout;
	//dlo7->setSizeConstraint(QLayout::SetFixedSize);
	dlo7->addWidget(dL71,0,0,1,1);dlo7->addWidget(dL72,0,1,1,1);dlo7->addWidget(dL73,0,2,1,1);
	dInBox7->setLayout(dlo7);
	dInBox7->setStyleSheet("background-color:grey;");
	// Isopropanol
	QLabel* dL81=new QLabel;QLabel* dL82=new QLabel;QLabel* dL83=new QLabel;
	dL81->setText(isopropanol.c_str());//dL81->setAlignment(Qt::AlignCenter);	
	dL82->setText("<b><span style=\"background-color:white; color:black; font-size:12pt;\">0.993868 − 3.40464x10<sup>-8</sup>t − 6.95191x10<sup>-6</sup>t<sup>2</sup> + 2.40158x10<sup>-8</sup>t<sup>3</sup> − 2.92637x10<sup>-11</sup>t<sup>4</sup></span></b>");dL82->setAlignment(Qt::AlignCenter);
	//dL83->setFixedHeight(imgHeight);
	dL83->setFixedWidth(imgWidth);
	dL83->setAlignment(Qt::AlignCenter);dL83->setPixmap(iDensPixmap.scaled(imgWidth,imgHeight,Qt::KeepAspectRatio));dL83->setMask(iDensPixmap.mask());
	QGridLayout *dlo8=new QGridLayout;
	dlo8->addWidget(dL81,0,0,1,1);dlo8->addWidget(dL82,0,1,1,1);dlo8->addWidget(dL83,0,2,1,1);
	dInBox8->setLayout(dlo8);
	dInBox8->setStyleSheet("background-color:white;");
	// Methanol
	QLabel* dL91=new QLabel;QLabel* dL92=new QLabel;QLabel* dL93=new QLabel;
	dL91->setText(methanol.c_str());//dL91->setAlignment(Qt::AlignCenter);	
	dL92->setText("<b><span style=\"background-color:grey; color:black; font-size:12pt;\">1.14445 − 1.79908x10<sup>-3</sup>t + 3.1645x10<sup>-6</sup>t<sup>2</sup> − 3.87839x10<sup>-9</sup>t<sup>3</sup></span></b>");dL92->setAlignment(Qt::AlignCenter);
	//dL93->setFixedHeight(imgHeight);
	dL93->setFixedWidth(imgWidth);
	dL93->setAlignment(Qt::AlignCenter);dL93->setPixmap(mDensPixmap.scaled(imgWidth,imgHeight,Qt::KeepAspectRatio));dL93->setMask(mDensPixmap.mask());
	QGridLayout *dlo9=new QGridLayout;
	dlo9->addWidget(dL91,0,0,1,1);dlo9->addWidget(dL92,0,1,1,1);dlo9->addWidget(dL93,0,2,1,1);
	dInBox9->setLayout(dlo9);
	dInBox9->setStyleSheet("background-color:grey;");
	// Propanol
	QLabel* dL101=new QLabel;QLabel* dL102=new QLabel;QLabel* dL103=new QLabel;
	dL101->setText(propanol.c_str());//dL101->setAlignment(Qt::AlignCenter);	
	dL102->setText("<b><span style=\"background-color:white; color:black; font-size:12pt;\">1.01077 − 3.99649x10<sup>-8</sup>t − 6.64923x10<sup>-6</sup>t<sup>2</sup> + 2.16751x10<sup>-8</sup>t<sup>3</sup> − 2.46167x10<sup>-11</sup>t<sup>4</sup></span></b>");dL102->setAlignment(Qt::AlignCenter);
	//dL103->setFixedHeight(imgHeight);
	dL103->setFixedWidth(imgWidth);
	dL103->setAlignment(Qt::AlignCenter);dL103->setPixmap(pDensPixmap.scaled(imgWidth,imgHeight,Qt::KeepAspectRatio));dL103->setMask(pDensPixmap.mask());
	QGridLayout *dlo10=new QGridLayout;
	dlo10->addWidget(dL101,0,0,1,1);dlo10->addWidget(dL102,0,1,1,1);dlo10->addWidget(dL103,0,2,1,1);
	dInBox10->setLayout(dlo10);
	dInBox10->setStyleSheet("background-color:white;");
	// Water
	QLabel* dL111=new QLabel;QLabel* dL112=new QLabel;QLabel* dL113=new QLabel;
	dL111->setText(water.c_str());//dL111->setAlignment(Qt::AlignCenter);
	s="<b><span style=\"background-color:grey; color:black; font-size:12pt;\">((((−2.8054253x10<sup>-13</sup>y + 1.0556302x10<sup>-10</sup>)y − 4.6170461x10<sup>-8</sup>)y − 7.9870401x10<sup>-6</sup>)y + 0.016945176)y + 0.99983952<br></span></b>";
	s+="<b><span style=\"background-color:grey; color:black; font-size:4pt;\">---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------<br></span></b>";
	s+="<b><span style=\"background-color:grey; color:black; font-size:12pt;\">1 + 0.01687985y<br><br></span></b>";
	s+="<b><span style=\"background-color:grey; color:black; font-size:12pt;\">y = t − 273.15</span></b>";
	dL112->setText(tr(s.c_str()));dL112->setAlignment(Qt::AlignCenter);
	dL113->setFixedWidth(imgWidth);
	dL113->setAlignment(Qt::AlignCenter);dL113->setPixmap(wDensPixmap.scaled(imgWidth,imgHeight,Qt::KeepAspectRatio));dL113->setMask(wDensPixmap.mask());
	QGridLayout *dlo11=new QGridLayout;
	dlo11->addWidget(dL111,0,0,1,1);dlo11->addWidget(dL112,0,1,1,1);dlo11->addWidget(dL113,0,2,1,1);
	dInBox11->setLayout(dlo11);
	dInBox11->setStyleSheet("background-color:grey;");
	// Mixed Solvents
	QGroupBox* dInBox12=new QGroupBox(tr("Mixed Solvents"));
	QLabel* dL121=new QLabel;QLabel* dL122=new QLabel;QLabel* dL123=new QLabel;
	dL121->setText("ZPRED uses a simple mixing rule<b><sup>[15]</sup></b> for ideal solvent mixtures, which assumes the solvent components do not interact with each other when mixed (i.e. no excess volume). This is well known to be a bad assumption. However, since this work focuses on aqueous solutions involving a pure solvent of only water this issue did not affect validation of ZPRED and has been presented with the intention to be corrected in future work.<br><br>");
	dL121->setWordWrap(true);
	dL121->setFixedWidth(propertyBoxSz);
	s="<b><span style=\"color:black; font-size:12pt;\">&rho;<sub>s,mix</sub> = &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1<br></span></b>";
	s+="<b><span style=\"color:black; font-size:4pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;----------------------------------------<br></span></b>";
	s+="<b><span style=\"color:black; font-size:12pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&Sigma;<sup>N</sup><sub>i=1</sub>(w<sub>i</sub>/&rho;<sub>s,i</sub>)<br></span></b>";
	dL122->setText(s.c_str());
	s="where <b>N</b> is the number of solvents, <b>&rho;<sub>s,i</sub></b> is the i-th pure solvent density and <b>w<sub>i</sub></b> is the mass fraction of the i-th solvent";
	dL123->setText(s.c_str());
	dL123->setWordWrap(true);
	dL123->setFixedWidth(propertyBoxSz);
	QGridLayout *dlo12=new QGridLayout;
	dlo12->addWidget(dL121,0,0);dlo12->addWidget(dL122,1,0);dlo12->addWidget(dL123,2,0);
	dInBox12->setLayout(dlo12);	
	// Accounting for Solutes
	QGroupBox* dInBox13=new QGroupBox(tr("Accounting for Solutes"));
	QLabel* dL131=new QLabel;QLabel* dL132=new QLabel;QLabel* dL133=new QLabel;
	dL131->setText("The density of an aqueous solution typically increases in the presence of salts. ZPRED uses a generalized aqueous electrolyte solution density model<b><sup>[15]</sup></b>, which applies the following mixing rule to define the density (<b>&rho;</b>) of aqueous electrolyte solutions.<br><br>");
	dL131->setWordWrap(true);
	dL131->setFixedWidth(propertyBoxSz);
	s="<b><span style=\"color:black; font-size:12pt;\">&rho; = &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1<br></span></b>";
	s+="<b><span style=\"color:black; font-size:4pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-------------------------------------------------------------------------------<br></span></b>";
	s+="<b><span style=\"color:black; font-size:12pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(w<sub>s</sub>/&rho;<sub>s</sub>)+&Sigma;<sup>M</sup><sub>i=1</sub>(w<sub>i</sub>/&rho;<sub>i</sub>)<br><br></span></b>";
	s+="<b><span style=\"color:black; font-size:12pt;\">1&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;w<sub>i</sub>+C<sub>2</sub>+C<sub>3</sub>(t-273.15)<br></span></b>";
	s+="<b><span style=\"color:black; font-size:4pt;\">----&nbsp;&nbsp;</span><span style=\"color:black; font-size:12pt;\">&nbsp;&nbsp;=</span><span style=\"color:black; font-size:4pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;---------------------------------------------------------------------------------------------------------------------------------------------------------<br></span></b>";
	s+="<b><span style=\"color:black; font-size:12pt;\">&rho;<sub>i</sub>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(C<sub>0</sub>w<sub>i</sub>+C<sub>1</sub>)exp(0.000001(t+C<sub>4</sub>-273.15)<sup>2</sup>)<br></span></b>";
	dL132->setText(s.c_str());
	s="where <b>w<sub>s</sub></b> is the solvent mass fraction, <b>&rho;<sub>s</sub></b> is the pure solvent density, <b>M</b> is the number of solutes, <b>w<sub>i</sub></b> is the i-th solute mass fraction, <b>&rho;<sub>i</sub></b> is the i-th solute density defined by the empirical relationship (above), <b>C<sub>0-4</sub></b> are empirical constants fit to experimental data and <b>t</b> is absolute temperature";
	dL133->setText(s.c_str());
	dL133->setWordWrap(true);
	dL133->setFixedWidth(propertyBoxSz);
	QGridLayout *dlo13=new QGridLayout;
	dlo13->addWidget(dL131,0,0);dlo13->addWidget(dL132,1,0);dlo13->addWidget(dL133,2,0);
	dInBox13->setLayout(dlo13);
	dInBox13->setStyleSheet("background-color:grey;");

	QGridLayout *dinLo = new QGridLayout;
	dinLo->setSizeConstraint(QLayout::SetFixedSize);
	dinLo->addWidget(dInBox2,0,0);
	dinLo->addWidget(dInBox3,1,0);
	dinLo->addWidget(dInBox4,2,0);
	dinLo->addWidget(dInBox5,3,0);
	dinLo->addWidget(dInBox6,4,0);
	dinLo->addWidget(dInBox7,5,0);
	dinLo->addWidget(dInBox8,6,0);
	dinLo->addWidget(dInBox9,7,0);
	dinLo->addWidget(dInBox10,8,0);
	dinLo->addWidget(dInBox11,9,0);
	dInBox->setLayout(dinLo);
	dInBox->setStyleSheet("background-color:black;");

	// Scrolling Group Box
	QScrollArea* dscrollArea=new QScrollArea();
	dscrollArea->setWidget(dInBox);  
	dscrollArea->setWidgetResizable(true);  
	QGridLayout *doutlo=new QGridLayout;
	doutlo->setSizeConstraint(QLayout::SetFixedSize);
	doutlo->addWidget(dDesc,0,0);
	doutlo->addWidget(dInBox1,1,0);
	doutlo->addWidget(dscrollArea,2,0);
	doutlo->addWidget(dInBox12,3,0);
	doutlo->addWidget(dInBox13,4,0);
	QGroupBox* dOutBox=new QGroupBox(tr("|-DENSITY-|"));
	dOutBox->setLayout(doutlo);

	// Viscosity Description
	QLabel* vDesc=new QLabel;
	s="The <b>viscosity</b> (also called <b>dynamic viscosity</b> or <b>absolute viscosity</b>) (<b>η</b>) of a substance is a measure of its resistance to deformation under an applied shear stress. In other words, it expresses the magnitude of frictional force within the substance that resists motion while the substance is accelerated by a stress. To compute the viscosity of a solution, it is first necessary to predict the viscosity of the pure solvent(s) (<b>η<sub>s</sub></b>) used. Similar to before, pure solvent viscosities as a function of temperature are shown below for water and the other common solvents. Plots in the table show the fit of the empirical relations to large sets of experimental data and provide a comparison among solvents in the relevant temperature range for water.";
	//s=setStringWidth(s,width);
	vDesc->setText(tr(s.c_str()));
	vDesc->setWordWrap(true);
	vDesc->setFixedWidth(propertyBoxSz);

	// Specify External Executables
	QGroupBox* vInBox=new QGroupBox;
	QGroupBox* vInBox1=new QGroupBox;
	QGroupBox* vInBox2=new QGroupBox;
	QGroupBox* vInBox3=new QGroupBox;
	QGroupBox* vInBox4=new QGroupBox;
	QGroupBox* vInBox5=new QGroupBox;
	QGroupBox* vInBox6=new QGroupBox;
	QGroupBox* vInBox7=new QGroupBox;
	QGroupBox* vInBox8=new QGroupBox;
	QGroupBox* vInBox9=new QGroupBox;
	QGroupBox* vInBox10=new QGroupBox;
	QGroupBox* vInBox11=new QGroupBox;
	//erInBox->setFont(font);	
	QLabel* vL11=new QLabel;QLabel* vL12=new QLabel;QLabel* vL13=new QLabel;
	vL11->setText("<b><span style=\"background-color:darkgrey; color:black; font-size:18pt;\">Solvent</span></b>");vL11->setAlignment(Qt::AlignCenter);	
	vL12->setText("<b><span style=\"background-color:darkgrey; color:black; font-size:18pt;\">η<sub>s</sub> [mPas] as Function of Absolute Temperature(t)</span></b>");vL12->setAlignment(Qt::AlignCenter);
	vL13->setText("<b><span style=\"background-color:darkgrey; color:black; font-size:18pt;\"> </span></b>");
	//
	QGridLayout *vlo1=new QGridLayout;
	//vlo1->setSizeConstraint(QLayout::SetFixedSize);
	vlo1->addWidget(vL11,0,0,1,1);vlo1->addWidget(vL12,0,1,1,1);vlo1->addWidget(vL13,0,2,1,1);
	vInBox1->setLayout(vlo1);
	vInBox1->setStyleSheet("background-color:darkgrey;");
	// Acetic Acid
	QLabel* vL21=new QLabel;QLabel* vL22=new QLabel;QLabel* vL23=new QLabel;
	vL21->setText(aceticAcid.c_str());//vL21->setAlignment(Qt::AlignCenter);	
	vL22->setText("<b><span style=\"background-color:white; color:black; font-size:12pt;\">66.46849exp(−1.360185x10<sup>-2</sup>t)</span></b>");vL22->setAlignment(Qt::AlignCenter);
	//vL23->setFixedHeight(imgHeight);
	vL23->setFixedWidth(imgWidth);
	vL23->setAlignment(Qt::AlignCenter);vL23->setPixmap(aaViscPixmap.scaled(imgWidth,imgHeight,Qt::KeepAspectRatio));vL23->setMask(aaViscPixmap.mask());
	QGridLayout *vlo2=new QGridLayout;
	//vlo2->setSizeConstraint(QLayout::SetFixedSize);
	vlo2->addWidget(vL21,0,0,1,1);vlo2->addWidget(vL22,0,1,1,1);vlo2->addWidget(vL23,0,2,1,1);
	vInBox2->setLayout(vlo2);
	vInBox2->setStyleSheet("background-color:white;");
	// Acetone
	QLabel* vL31=new QLabel;QLabel* vL32=new QLabel;QLabel* vL33=new QLabel;
	vL31->setText(acetone.c_str());//vL31->setAlignment(Qt::AlignCenter);	
	vL32->setText("<b><span style=\"background-color:grey; color:black; font-size:12pt;\">4.62452exp(−9.015176x10<sup>-3</sup>t)</span></b>");vL32->setAlignment(Qt::AlignCenter);
	//vL33->setFixedHeight(imgHeight);
	vL33->setFixedWidth(imgWidth);
	vL33->setAlignment(Qt::AlignCenter);vL33->setPixmap(acViscPixmap.scaled(imgWidth,imgHeight,Qt::KeepAspectRatio));vL33->setMask(acViscPixmap.mask());
	QGridLayout *vlo3=new QGridLayout;
	vlo3->addWidget(vL31,0,0,1,1);vlo3->addWidget(vL32,0,1,1,1);vlo3->addWidget(vL33,0,2,1,1);
	vInBox3->setLayout(vlo3);
	vInBox3->setStyleSheet("background-color:grey;");
	// Butanol
	QLabel* vL41=new QLabel;QLabel* vL42=new QLabel;QLabel* vL43=new QLabel;
	vL41->setText(butanol.c_str());//vL41->setAlignment(Qt::AlignCenter);	
	vL42->setText("<b><span style=\"background-color:white; color:black; font-size:12pt;\">2146.186exp(−2.256471x10<sup>-2</sup>t)</span></b>");vL42->setAlignment(Qt::AlignCenter);
	//vL43->setFixedHeight(imgHeight);
	vL43->setFixedWidth(imgWidth);
	vL43->setAlignment(Qt::AlignCenter);vL43->setPixmap(bViscPixmap.scaled(imgWidth,imgHeight,Qt::KeepAspectRatio));vL43->setMask(bViscPixmap.mask());
	QGridLayout *vlo4=new QGridLayout;
	vlo4->addWidget(vL41,0,0,1,1);vlo4->addWidget(vL42,0,1,1,1);vlo4->addWidget(vL43,0,2,1,1);
	vInBox4->setLayout(vlo4);
	vInBox4->setStyleSheet("background-color:white;");
	// DMSO
	QLabel* vL51=new QLabel;QLabel* vL52=new QLabel;QLabel* vL53=new QLabel;
	vL51->setText(dmso.c_str());//vL51->setAlignment(Qt::AlignCenter);	
	vL52->setText("<b><span style=\"background-color:grey; color:black; font-size:12pt;\">126.2972exp(−1.406096x10<sup>-2</sup>t)</span></b>");vL52->setAlignment(Qt::AlignCenter);
	//vL53->setFixedHeight(imgHeight);
	vL53->setFixedWidth(imgWidth);
	vL53->setAlignment(Qt::AlignCenter);vL53->setPixmap(dViscPixmap.scaled(imgWidth,imgHeight,Qt::KeepAspectRatio));vL53->setMask(dViscPixmap.mask());
	QGridLayout *vlo5=new QGridLayout;
	vlo5->addWidget(vL51,0,0,1,1);vlo5->addWidget(vL52,0,1,1,1);vlo5->addWidget(vL53,0,2,1,1);
	vInBox5->setLayout(vlo5);
	vInBox5->setStyleSheet("background-color:grey;");
	// Ethanol
	QLabel* vL61=new QLabel;QLabel* vL62=new QLabel;QLabel* vL63=new QLabel;
	//vL61->setText("<span style=\"background-color:white; color:black; font-size:12pt;\">Ethanol</span>");//vL61->setAlignment(Qt::AlignCenter);	
	vL61->setText(ethanol.c_str());
	vL62->setText("<b><span style=\"background-color:white; color:black; font-size:12pt;\">10^(−6.4406 + 1120/t + 0.0137t − 1.55x10<sup>-5</sup>t<sup>2</sup>)</span></b>");vL62->setAlignment(Qt::AlignCenter);
	//vL63->setFixedHeight(imgHeight);
	vL63->setFixedWidth(imgWidth);
	vL63->setAlignment(Qt::AlignCenter);vL63->setPixmap(eViscPixmap.scaled(imgWidth,imgHeight,Qt::KeepAspectRatio));vL63->setMask(eViscPixmap.mask());
	QGridLayout *vlo6=new QGridLayout;
	vlo6->addWidget(vL61,0,0,1,1);vlo6->addWidget(vL62,0,1,1,1);vlo6->addWidget(vL63,0,2,1,1);
	vInBox6->setLayout(vlo6);
	vInBox6->setStyleSheet("background-color:white;");
	// Glycerol
	QLabel* vL71=new QLabel;QLabel* vL72=new QLabel;QLabel* vL73=new QLabel;
	vL71->setText(glycerol.c_str());//vL71->setAlignment(Qt::AlignCenter);	
	vL72->setText("<b><span style=\"background-color:grey; color:black; font-size:12pt;\">10^(−18.215 + 4230/t + 0.0287t − 1.86x10<sup>-5</sup>t<sup>2</sup>)</span></b>");vL72->setAlignment(Qt::AlignCenter);
	//vL73->setFixedHeight(imgHeight);
	vL73->setFixedWidth(imgWidth);
	vL73->setAlignment(Qt::AlignCenter);vL73->setPixmap(gViscPixmap.scaled(imgWidth,imgHeight,Qt::KeepAspectRatio));vL73->setMask(gViscPixmap.mask());
	QGridLayout *vlo7=new QGridLayout;
	vlo7->addWidget(vL71,0,0,1,1);vlo7->addWidget(vL72,0,1,1,1);vlo7->addWidget(vL73,0,2,1,1);
	vInBox7->setLayout(vlo7);
	vInBox7->setStyleSheet("background-color:grey;");
	// Isopropanol
	QLabel* vL81=new QLabel;QLabel* vL82=new QLabel;QLabel* vL83=new QLabel;
	vL81->setText(isopropanol.c_str());//vL81->setAlignment(Qt::AlignCenter);	
	vL82->setText("<b><span style=\"background-color:white; color:black; font-size:12pt;\">10^(−0.7009 + 842/t − 0.0086t + 8.3x10<sup>-6</sup>t<sup>2</sup>)</span></b>");vL82->setAlignment(Qt::AlignCenter);
	//vL83->setFixedHeight(imgHeight);
	vL83->setFixedWidth(imgWidth);
	vL83->setAlignment(Qt::AlignCenter);vL83->setPixmap(iViscPixmap.scaled(imgWidth,imgHeight,Qt::KeepAspectRatio));vL83->setMask(iViscPixmap.mask());
	QGridLayout *vlo8=new QGridLayout;
	vlo8->addWidget(vL81,0,0,1,1);vlo8->addWidget(vL82,0,1,1,1);vlo8->addWidget(vL83,0,2,1,1);
	vInBox8->setLayout(vlo8);
	vInBox8->setStyleSheet("background-color:white;");
	// Methanol
	QLabel* vL91=new QLabel;QLabel* vL92=new QLabel;QLabel* vL93=new QLabel;
	vL91->setText(methanol.c_str());//vL91->setAlignment(Qt::AlignCenter);	
	vL92->setText("<b><span style=\"background-color:grey; color:black; font-size:12pt;\">10^(−9.0562 + 1250/t + 0.0224t − 2.35x10<sup>-5</sup>t<sup>2</sup>)</span></b>");vL92->setAlignment(Qt::AlignCenter);
	//vL93->setFixedHeight(imgHeight);
	vL93->setFixedWidth(imgWidth);
	vL93->setAlignment(Qt::AlignCenter);vL93->setPixmap(mViscPixmap.scaled(imgWidth,imgHeight,Qt::KeepAspectRatio));vL93->setMask(mViscPixmap.mask());
	QGridLayout *vlo9=new QGridLayout;
	vlo9->addWidget(vL91,0,0,1,1);vlo9->addWidget(vL92,0,1,1,1);vlo9->addWidget(vL93,0,2,1,1);
	vInBox9->setLayout(vlo9);
	vInBox9->setStyleSheet("background-color:grey;");
	// Propanol
	QLabel* vL101=new QLabel;QLabel* vL102=new QLabel;QLabel* vL103=new QLabel;
	vL101->setText(propanol.c_str());//vL101->setAlignment(Qt::AlignCenter);	
	vL102->setText("<b><span style=\"background-color:white; color:black; font-size:12pt;\">10^(−3.7702 + 992/t + 0.0041t − 5.46x10<sup>-6</sup>t<sup>2</sup>)</span></b>");vL102->setAlignment(Qt::AlignCenter);
	//vL103->setFixedHeight(imgHeight);
	vL103->setFixedWidth(imgWidth);
	vL103->setAlignment(Qt::AlignCenter);vL103->setPixmap(pViscPixmap.scaled(imgWidth,imgHeight,Qt::KeepAspectRatio));vL103->setMask(pViscPixmap.mask());
	QGridLayout *vlo10=new QGridLayout;
	vlo10->addWidget(vL101,0,0,1,1);vlo10->addWidget(vL102,0,1,1,1);vlo10->addWidget(vL103,0,2,1,1);
	vInBox10->setLayout(vlo10);
	vInBox10->setStyleSheet("background-color:white;");
	// Water
	QLabel* vL111=new QLabel;QLabel* vL112=new QLabel;QLabel* vL113=new QLabel;
	vL111->setText(water.c_str());//vL111->setAlignment(Qt::AlignCenter);
	//s="<span style=\"background-color:grey; color:black; font-size:12pt;\">((((−2.8054253x10⁻¹³y + 1.0556302x10⁻¹⁰)y − 4.6170461x10⁻⁸)y − 7.9870401x10⁻⁶)y + 0.016945176)y + 0.99983952<br></span>";
	s="<b><span style=\"background-color:grey; color:black; font-size:12pt;\">(t − 273.15) + 246<br></span></b>";
	s+="<b><span style=\"background-color:grey; color:black; font-size:4pt;\">-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------<br></span></b>";
	s+="<b><span style=\"background-color:grey; color:black; font-size:12pt;\">137.37 + 5.2842(t − 273.15) + 0.05594(t − 273.15)<sup>2</sup></span></b>";
	vL112->setText(s.c_str());vL112->setAlignment(Qt::AlignCenter);
	//vL113->setFixedHeight(imgHeight);
	vL113->setFixedWidth(imgWidth);
	vL113->setAlignment(Qt::AlignCenter);vL113->setPixmap(wViscPixmap.scaled(imgWidth,imgHeight,Qt::KeepAspectRatio));vL113->setMask(wViscPixmap.mask());
	QGridLayout *vlo11=new QGridLayout;
	vlo11->addWidget(vL111,0,0,1,1);vlo11->addWidget(vL112,0,1,1,1);vlo11->addWidget(vL113,0,2,1,1);
	vInBox11->setLayout(vlo11);
	vInBox11->setStyleSheet("background-color:grey;");
	// Mixed Solvents
	QGroupBox* vInBox12=new QGroupBox(tr("Mixed Solvents"));
	QLabel* vL121=new QLabel;QLabel* vL122=new QLabel;QLabel* vL123=new QLabel;
	vL121->setText("ZPRED applies the Arrhenius mixing rule<b><sup>[18]</sup></b> for ideal solvent mixtures as shown:<br><br>");
	vL121->setWordWrap(true);
	vL121->setFixedWidth(propertyBoxSz);
	s="<b><span style=\"color:black; font-size:12pt;\">&eta;<sub>s,mix</sub> = exp(&Sigma;<sup>N</sup><sub>i=1</sub>w<sub>i</sub>ln(&eta;<sub>i</sub>))<br></span></b>";
	vL122->setText(s.c_str());
	s="where <b>N</b> is the number of solvents, <b>w<sub>i</sub></b> is the mass fraction occupied by the i-th solvent and <b>&eta;<sub>i</sub></b> is the i-th solvent viscosity";
	vL123->setText(s.c_str());
	vL123->setWordWrap(true);
	vL123->setFixedWidth(propertyBoxSz);
	QGridLayout *vlo12=new QGridLayout;
	vlo12->addWidget(vL121,0,0);vlo12->addWidget(vL122,1,0);vlo12->addWidget(vL123,2,0);
	vInBox12->setLayout(vlo12);	
	// Accounting for Solutes
	QGroupBox* vInBox13=new QGroupBox(tr("Accounting for Solutes"));
	QLabel* vL131=new QLabel;QLabel* vL132=new QLabel;QLabel* vL133=new QLabel;
	vL131->setText("The viscosity of an aqueous solution can either increase or decrease in the presences of salts<b><sup>[20]</sup></b>. ZPRED uses a generalized viscosity model<b><sup>[18]</sup></b> for aqueous electrolyte solutions, which applies the following mixing rule to define the viscosity (<b>&eta;</b>) of an aqueous electrolyte solution.<br><br>");
	vL131->setWordWrap(true);
	vL131->setFixedWidth(propertyBoxSz);
	s="<b><span style=\"color:black;font-size:12pt;\">&eta; = (&eta;<sub>s</sub>)^(w<sub>s</sub>)&Pi;<sup>M</sup><sub>i=1</sub>(&eta;<sub>i</sub>)^(w<sub>i</sub>)<br></span></b>";
	s+="<b><span style=\"color:black;font-size:12pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>_</b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>_</b><br></span></b>";
	s+="<b><span style=\"color:black;font-size:12pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>|</b>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;C<sub>1</sub>(1-w<sub>s</sub>)^(C<sub>2</sub>)+C<sub>3</sub>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>|</b><br></span></b>";
	s+="<b><span style=\"color:black;font-size:12pt;\">&eta;<sub>i</sub> = exp<b>|</b></span><span style=\"color:black; font-size:4pt;\">-----------------------------------------------------------------------------------------------------------------------------------------------------------</span><span style=\"color:black;font-size:12pt;\"><b>|</b><br></span></b>";
	s+="<b><span style=\"color:black;font-size:12pt;\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<b>|_</b>(1+C<sub>4</sub>(t-273.15))(1+C<sub>5</sub>(1-w<sub>s</sub>)^(C<sub>6</sub>))<b>_|</b><br></span></b>";
	vL132->setText(s.c_str());
	s="where <b>w<sub>s</sub></b> is the solvent mass fraction, <b>&eta;<sub>s</sub></b> is the pure solvent viscosity, <b>M</b> is the number of solutes, <b>w<sub>i</sub></b> is the i-th solute mass fraction, <b>&eta;<sub>i</sub></b> is the i-th solute viscosity defined by the empirical relationship (above), <b>C<sub>1-6</sub></b> are empirical constants fit to experimental data and <b>t</b> is absolute temperature";
	vL133->setText(s.c_str());
	vL133->setWordWrap(true);
	vL133->setFixedWidth(propertyBoxSz);
	QGridLayout *vlo13=new QGridLayout;
	vlo13->addWidget(vL131,0,0);vlo13->addWidget(vL132,1,0);vlo13->addWidget(vL133,2,0);
	vInBox13->setLayout(vlo13);
	vInBox13->setStyleSheet("background-color:grey;");

	QGridLayout *vinLo = new QGridLayout;
	vinLo->setSizeConstraint(QLayout::SetFixedSize);
	vinLo->addWidget(vInBox2,0,0);
	vinLo->addWidget(vInBox3,1,0);
	vinLo->addWidget(vInBox4,2,0);
	vinLo->addWidget(vInBox5,3,0);
	vinLo->addWidget(vInBox6,4,0);
	vinLo->addWidget(vInBox7,5,0);
	vinLo->addWidget(vInBox8,6,0);
	vinLo->addWidget(vInBox9,7,0);
	vinLo->addWidget(vInBox10,8,0);
	vinLo->addWidget(vInBox11,9,0);
	vInBox->setLayout(vinLo);
	vInBox->setStyleSheet("background-color:black;");

	// Scrolling Group Box
	QScrollArea* vscrollArea=new QScrollArea();
	vscrollArea->setWidget(vInBox);  
	vscrollArea->setWidgetResizable(true);  
	QGridLayout *voutlo=new QGridLayout;
	voutlo->setSizeConstraint(QLayout::SetFixedSize);
	voutlo->addWidget(vDesc,0,0);
	voutlo->addWidget(vInBox1,1,0);
	voutlo->addWidget(vscrollArea,2,0);
	voutlo->addWidget(vInBox12,3,0);
	voutlo->addWidget(vInBox13,4,0);
	QGroupBox* vOutBox=new QGroupBox(tr("|-VISCOSITY-|"));
	//vOutBox->setFont(font);
	vOutBox->setLayout(voutlo);

	// Reference List
	QGroupBox* rBox=new QGroupBox(tr("|-REFERENCES-|"));
	QLabel* rL11=new QLabel;QLabel* rL12=new QLabel;
	rL11->setText("<b>[1]</b>");rL11->setAlignment(Qt::AlignCenter);
	s="D. R. Lide, Ed., <u>CRC Handbook of Chemistry and Physics</u>, 86th ed., Boca Raton, FL: Taylor and Francis, 2006 [Online]. Available: <a href=\"http://hbcponline.com/faces/contents/ContentsSearch.xhtml\">http://hbcponline.com/faces/contents/ContentsSearch.xhtml</a>";	
	rL12->setText(s.c_str());rL12->setAlignment(Qt::AlignLeft);rL12->setTextInteractionFlags(Qt::TextBrowserInteraction);rL12->setOpenExternalLinks(true);
	QLabel* rL21=new QLabel;QLabel* rL22=new QLabel;
	rL21->setText("<b>[2]</b>");rL21->setAlignment(Qt::AlignCenter);
	s="C. Wohlfarth and O. Madelung, Eds., Static Dielectric Constants of Pure Liquids and Binary Liquid Mixtures, New series, vol. 6, Springer, 1991 [doi: <a href=\"https://doi.org/10.1007/978-3-662-48168-4\">10.1007/978-3-662-48168-4</a>]";
	rL22->setText(s.c_str());rL22->setAlignment(Qt::AlignLeft);rL22->setTextInteractionFlags(Qt::TextBrowserInteraction);rL22->setOpenExternalLinks(true);
	QLabel* rL31=new QLabel;QLabel* rL32=new QLabel;
	rL31->setText("<b>[3]</b>");rL31->setAlignment(Qt::AlignCenter);
	s="C. G. Malmberg and A. A. Maryott, \"Dielectric Constant of Water from 0 to 100 C,\" Journal of Research of the National Bureau of Standards, vol. 56, no. 1, pp. 1-8, 1956 [<a href=\"https://nvlpubs.nist.gov/nistpubs/jres/56/jresv56n1p1_a1b.pdf\">link</a>]";
	rL32->setText(s.c_str());rL32->setAlignment(Qt::AlignLeft);rL32->setTextInteractionFlags(Qt::TextBrowserInteraction);rL32->setOpenExternalLinks(true);
	QLabel* rL41=new QLabel;QLabel* rL42=new QLabel;
	rL41->setText("<b>[4]</b>");rL41->setAlignment(Qt::AlignCenter);
	s="C. Wohlfarth and M. D. Lechner, Eds., Landolt-Bornstein Numerical Data and Functional Relationships in Science and Technology Static Dielectric Constants of Pure Liquids and Binary Liquid Mixtures, New series, vol. 27, Springer, 2015 [doi: <a href=\"https://doi.org/10.1007/978-3-662-48168-4\">10.1007/978-3-662-48168-4</a>]";
	rL42->setText(s.c_str());rL42->setAlignment(Qt::AlignLeft);rL42->setTextInteractionFlags(Qt::TextBrowserInteraction);rL42->setOpenExternalLinks(true);
	QLabel* rL51=new QLabel;QLabel* rL52=new QLabel;
	rL51->setText("<b>[5]</b>");rL51->setAlignment(Qt::AlignCenter);
	s="P. Wang and A. Anderko, \"Computation of dielectric constants of solvent mixtures and electrolyte solutions,\" Fluid Phase Equilibria, vol. 186, pp.103-122, 2001 [doi: <a href=\"https://doi.org/10.1016/S0378-3812(01)00507-6\">10.1016/S0378-3812(01)00507-6</a>]";
	rL52->setText(s.c_str());rL52->setAlignment(Qt::AlignLeft);rL52->setTextInteractionFlags(Qt::TextBrowserInteraction);rL52->setOpenExternalLinks(true);
	QLabel* rL61=new QLabel;QLabel* rL62=new QLabel;
	rL61->setText("<b>[6]</b>");rL61->setAlignment(Qt::AlignCenter);
	s="J. B. Hasted, D. M. Ritson and C. H. Collie, \"Dielectric Properties of Aqueous Ionic Solutions. Part I and II,\" Journal of Chemical Physics, vol. 16, num. 1, pp. 1-21, 1948 [doi: <a href=\"https://doi.org/10.1063/1.1746645\">10.1063/1.1746645</a>]";
	rL62->setText(s.c_str());rL62->setAlignment(Qt::AlignLeft);rL62->setTextInteractionFlags(Qt::TextBrowserInteraction);rL62->setOpenExternalLinks(true);
	QLabel* rL71=new QLabel;QLabel* rL72=new QLabel;
	rL71->setText("<b>[7]</b>");rL71->setAlignment(Qt::AlignCenter);
	s="E. Glueckauf, \"Bulk dielectric constant of aqueous electrolyte solutions,\" Transactions of the Faraday Society, vol. 60, pp. 1637-1645, 1964 [doi: <a href=\"https://doi.org/10.1039/TF9646001637\">10.1039/TF9646001637</a>]";
	rL72->setText(s.c_str());rL72->setAlignment(Qt::AlignLeft);rL72->setTextInteractionFlags(Qt::TextBrowserInteraction);rL72->setOpenExternalLinks(true);
	QLabel* rL81=new QLabel;QLabel* rL82=new QLabel;
	rL81->setText("<b>[8]</b>");rL81->setAlignment(Qt::AlignCenter);
	s="A. Levy, D. Andelman and H. Orland, \"Dielectric Constant of Ionic Solutions: A Field-Theory Approach,\" Physical Review Letters, vol. 108, num. 22, pp. 227801, 2012 [doi: <a href=\"https://doi.org/10.1103/PhysRevLett.108.227801\">10.1103/PhysRevLett.108.227801</a>]";
	rL82->setText(s.c_str());rL82->setAlignment(Qt::AlignLeft);rL82->setTextInteractionFlags(Qt::TextBrowserInteraction);rL82->setOpenExternalLinks(true);
	QLabel* rL91=new QLabel;QLabel* rL92=new QLabel;
	rL91->setText("<b>[9]</b>");rL91->setAlignment(Qt::AlignCenter);
	s="N. Gavish and K. Promislow, \"Dependence of the dielectric constant of electrolyte solutions on ionic concentration: A microfield approach,\" Physical Review E, vol. 94, num. 1, pp. 012611, 2016 [doi: <a href=\"https://doi.org/10.1103/PhysRevE.94.012611\">10.1103/PhysRevE.94.012611</a>]";
	rL92->setText(s.c_str());rL92->setAlignment(Qt::AlignLeft);rL92->setTextInteractionFlags(Qt::TextBrowserInteraction);rL92->setOpenExternalLinks(true);
	QLabel* rL101=new QLabel;QLabel* rL102=new QLabel;
	rL101->setText("<b>[10]</b>");rL101->setAlignment(Qt::AlignCenter);
	s="R. H. Perry, D. W. Green and J. O. Maloney, <u>Perry's Chemical Engineers' Handbook</u>, 7th ed., McGraw-Hill, 1997";
	rL102->setText(s.c_str());rL102->setAlignment(Qt::AlignLeft);rL102->setTextInteractionFlags(Qt::TextBrowserInteraction);rL102->setOpenExternalLinks(true);
	QLabel* rL111=new QLabel;QLabel* rL112=new QLabel;
	rL111->setText("<b>[11]</b>");rL111->setAlignment(Qt::AlignCenter);
	s="M. Frenkel, X. Hong, Q. Dong, X. Yan and R. D. Chirico, <u>Landolt-Bornstein Numerical Data and Functional Relationships in Science and Technology Thermodynamic Properties of Organic Compounds and their Mixtures Subvolume I Densities of Phenols, Aldehydes, Ketones, Carboxylic Acids, Amines, Nitriles and Nitrohydrocarbons</u>, New series, Group IV: Physical Chemistry, vol. 8, Springer, 2002 [doi: <a href=\"https://doi.org/10.1007/b76771\">10.1007/b76771</a>]";
	rL112->setText(s.c_str());rL112->setAlignment(Qt::AlignLeft);rL112->setTextInteractionFlags(Qt::TextBrowserInteraction);rL112->setOpenExternalLinks(true);
	QLabel* rL121=new QLabel;QLabel* rL122=new QLabel;
	rL121->setText("<b>[12]</b>");rL121->setAlignment(Qt::AlignCenter);
	s="M. Frenkel, X. Hong, R. C. Willhoit and K. R. Hall, <u>Landolt-Bornstein Numerical Data and Functional Relationships in Science and Technology Thermodynamic Properties of Organic Compounds and their Mixtures Subvolume G Densities of Alcohols</u>, New series, Group IV: Physical Chemistry, vol. 8, Springer, 1998 [doi: <a href=\"https://doi.org/10.1007/b75928\">10.1007/b75928</a>]";
	rL122->setText(s.c_str());rL122->setAlignment(Qt::AlignLeft);rL122->setTextInteractionFlags(Qt::TextBrowserInteraction);rL122->setOpenExternalLinks(true);
	QLabel* rL131=new QLabel;QLabel* rL132=new QLabel;
	rL131->setText("<b>[13]</b>");rL131->setAlignment(Qt::AlignCenter);
	s="C. L. Yaws, <u>Yaw's Handbook of Thermodynamic and Physical Properties of Chemical Compounds</u>";
	rL132->setText(s.c_str());rL132->setAlignment(Qt::AlignLeft);rL132->setTextInteractionFlags(Qt::TextBrowserInteraction);rL132->setOpenExternalLinks(true);
	QLabel* rL141=new QLabel;QLabel* rL142=new QLabel;
	rL141->setText("<b>[14]</b>");rL141->setAlignment(Qt::AlignCenter);
	s="Gaylord Chemical, \"DMSO Physical Properties\" 2018. [Online]. Available: <a href=\"https://www.gaylordchemical.com/wp-content/uploads/2015/09/GC-Literature-101B-REV.pdf\">https://www.gaylordchemical.com/literature/dmso-physical-properties/</a>";
	rL142->setText(s.c_str());rL142->setAlignment(Qt::AlignLeft);rL142->setTextInteractionFlags(Qt::TextBrowserInteraction);rL142->setOpenExternalLinks(true);
	QLabel* rL151=new QLabel;QLabel* rL152=new QLabel;
	rL151->setText("<b>[15]</b>");rL151->setAlignment(Qt::AlignCenter);
	s="M. Laliberte and W. E. Cooper, \"Model for calculating the Density of Aqueous Electrolyte Solutions,\" Journal of Chemical Engineering Data, vol. 49, pp. 1141-1151, 2004 [doi: <a href=\"https://doi.org/10.1021/je0498659\">10.1021/je0498659</a>]";
	rL152->setText(s.c_str());rL152->setAlignment(Qt::AlignLeft);rL152->setTextInteractionFlags(Qt::TextBrowserInteraction);rL152->setOpenExternalLinks(true);
	QLabel* rL161=new QLabel;QLabel* rL162=new QLabel;
	rL161->setText("<b>[16]</b>");rL161->setAlignment(Qt::AlignCenter);
	s="C. T. Crowe, D. F. Elger, B. C. Williams and J. A. Roberson, <u>Engineering Fluid Mechanics</u>, 9th ed., Hoboken, NJ: John Wiley & Sons, Inc., 2009";
	rL162->setText(s.c_str());rL162->setAlignment(Qt::AlignLeft);rL162->setTextInteractionFlags(Qt::TextBrowserInteraction);rL162->setOpenExternalLinks(true);
	QLabel* rL171=new QLabel;QLabel* rL172=new QLabel;
	rL171->setText("<b>[17]</b>");rL171->setAlignment(Qt::AlignCenter);
	s="C. Wohlfarth, B. Wohlfarth and M. D. Lechner, <u>Landolt-Bornstein Numerical Data and Functional Relationships in Science and Technology Viscosity of Pure Organic Liquids and Binary Liquid Mixtures</u>, New series, Group IV: Physical Chemistry, vol. 2, Springer, 2002 [doi: <a href=\"https://doi.org/10.1007/b68239\">10.1007/b68239</a>]";
	rL172->setText(s.c_str());rL172->setAlignment(Qt::AlignLeft);rL172->setTextInteractionFlags(Qt::TextBrowserInteraction);rL172->setOpenExternalLinks(true);
	QLabel* rL181=new QLabel;QLabel* rL182=new QLabel;
	rL181->setText("<b>[18]</b>");rL181->setAlignment(Qt::AlignCenter);
	s="M. Laliberte, \"Model for Calculating the Viscosity of Aqueous Solutions,\" Journal of Chemical and Engineering Data, vol. 52, no. 2, pp. 321-335, 2007 [doi: <a href=\"https://doi.org/10.1021/je0604075\">10.1021/je0604075</a>]";
	rL182->setText(s.c_str());rL182->setAlignment(Qt::AlignLeft);rL182->setTextInteractionFlags(Qt::TextBrowserInteraction);rL182->setOpenExternalLinks(true);
	QLabel* rL191=new QLabel;QLabel* rL192=new QLabel;
	rL191->setText("<b>[19]</b>");rL191->setAlignment(Qt::AlignCenter);
	s="Y. A. Cengel and A. J. Ghajar, <u>Heat and Mass Transfer Fundamentals & Applications</u>, 4th ed., New York, NY: McGraw-Hill, 2007";
	rL192->setText(s.c_str());rL192->setAlignment(Qt::AlignLeft);rL192->setTextInteractionFlags(Qt::TextBrowserInteraction);rL192->setOpenExternalLinks(true);
	QLabel* rL201=new QLabel;QLabel* rL202=new QLabel;
	rL201->setText("<b>[20]</b>");rL201->setAlignment(Qt::AlignCenter);
	s="G. Jones and S. K. Talley, \"The Viscosity of Aqueous Solutions as a Function of the Concentration,\" Journal of the American Chemical Society, vol. 55, no. 2, pp. 624-642, 1933 [doi: <a href=\"https://doi.org/10.1021/ja01329a024\">10.1021/ja01329a024</a>]";
	rL202->setText(s.c_str());rL202->setAlignment(Qt::AlignLeft);rL202->setTextInteractionFlags(Qt::TextBrowserInteraction);rL202->setOpenExternalLinks(true);

	QGridLayout *routlo=new QGridLayout;
	routlo->setSizeConstraint(QLayout::SetFixedSize);
	routlo->addWidget(rL11,0,0);routlo->addWidget(rL12,0,1);
	routlo->addWidget(rL21,1,0);routlo->addWidget(rL22,1,1);
	routlo->addWidget(rL31,2,0);routlo->addWidget(rL32,2,1);
	routlo->addWidget(rL41,3,0);routlo->addWidget(rL42,3,1);
	routlo->addWidget(rL51,4,0);routlo->addWidget(rL52,4,1);
	routlo->addWidget(rL61,5,0);routlo->addWidget(rL62,5,1);
	routlo->addWidget(rL71,6,0);routlo->addWidget(rL72,6,1);
	routlo->addWidget(rL81,7,0);routlo->addWidget(rL82,7,1);
	routlo->addWidget(rL91,8,0);routlo->addWidget(rL92,8,1);
	routlo->addWidget(rL101,9,0);routlo->addWidget(rL102,9,1);
	routlo->addWidget(rL111,10,0);routlo->addWidget(rL112,10,1);
	routlo->addWidget(rL121,11,0);routlo->addWidget(rL122,11,1);
	routlo->addWidget(rL131,12,0);routlo->addWidget(rL132,12,1);
	routlo->addWidget(rL141,13,0);routlo->addWidget(rL142,13,1);
	routlo->addWidget(rL151,14,0);routlo->addWidget(rL152,14,1);
	routlo->addWidget(rL161,15,0);routlo->addWidget(rL162,15,1);
	routlo->addWidget(rL171,16,0);routlo->addWidget(rL172,16,1);
	routlo->addWidget(rL181,17,0);routlo->addWidget(rL182,17,1);
	routlo->addWidget(rL191,18,0);routlo->addWidget(rL192,18,1);
	routlo->addWidget(rL201,19,0);routlo->addWidget(rL202,19,1);
	rBox->setLayout(routlo);

	QGridLayout *ilayout = new QGridLayout;
	//ilayout->setSizeConstraint(QLayout::SetFixedSize);
	ilayout->addWidget(description,0,0);
	ilayout->addWidget(erOutBox,1,0);
	ilayout->addWidget(dOutBox,2,0);
	ilayout->addWidget(vOutBox,3,0);
	ilayout->addWidget(rBox,4,0);
	gBox->setLayout(ilayout);

	// Scrolling Group Box
	QScrollArea* gscrollArea=new QScrollArea();
	gscrollArea->setWidget(gBox);  
	gscrollArea->setWidgetResizable(true);  
	QGridLayout *goutlo=new QGridLayout;
	goutlo->addWidget(gscrollArea);
	QGroupBox* gOutBox=new QGroupBox;//(tr("|-SOLVENT PROPERTIES COMPUTATION INFORMATION-|"));
	//gOutBox->setFont(font);
	gOutBox->setLayout(goutlo);

	QGridLayout *mLayout = new QGridLayout;
	//mLayout->setSizeConstraint(QLayout::SetFixedSize);
	mLayout->addWidget(gOutBox);
	setLayout(mLayout);
	setWindowTitle(tr("Solvent Property Computation Information"));
	resize(QDesktopWidget().availableGeometry(this).size());
}

msgBox::msgBox(QWidget *parent) : QWidget(parent)
{	
	int width=800;
	QString Fldr=QDir::currentPath();
	string sFldr=Fldr.toStdString()+"/";
	// File Reference for When Optional Parameters are set (makes .optionalParameters.txt)
	string msgFile=sFldr+".message.txt",tmp="";
	ifstream fIn;
	fIn.open(msgFile.c_str());
	if(fIn.fail()){cerr<<"Error in msgBox!\nCould not open msg file ("<<msgFile<<")\n";exit(EXIT_FAILURE);}
	uint Sz=150000;
	char Val[Sz];
	fIn.getline(Val,Sz);
	tmp=Val;fIn.close();
	displayLabel=new QLabel;displayLabel->setText(tmp.c_str());
	displayLabel->setWordWrap(true);
	displayLabel->setFixedWidth(width);
	
	QGridLayout *mLayout = new QGridLayout;
	mLayout->setSizeConstraint(QLayout::SetFixedSize);
	mLayout->addWidget(displayLabel);
	setLayout(mLayout);
	setWindowTitle(tr("ZPRED Message"));
	//resize(QDesktopWidget().availableGeometry(this).size());
}

apbsInstallInfo::apbsInstallInfo(QWidget *parent) : QWidget(parent)
{	
	int width=500;
	string s="<b><u>Ubuntu</u></b><br>";
	s+="&#8226;&nbsp;&nbsp;Download APBS 1.5 at the following link: ";
	s+="<a href=\"https://sourceforge.net/projects/apbs/files/apbs/apbs-1.5/APBS-1.5-linux64.tar.gz\"><span style=\"color:blue;\">download</span></a><br><br>";
	s+="&#8226;&nbsp;&nbsp;After downloading the compressed APBS package, decompress it using the following command in a terminal:<br><br>";
	s+="<b>tar xzf APBS-1.5-linux64.tar.gz --strip-components 1 -C /usr/local</b><br><br>";
	s+="<b>Note</b>: It is necessary to enter the following commands as well:<br>";
	s+="&#8226;&nbsp;&nbsp;To install a necessary library:<br><br>";
	s+="<b>sudo apt-get install libgfortran3</b><br><br>";
	s+="&#8226;&nbsp;&nbsp;To correct a linking error with libreadline.so.6:<br><br>";
	s+="<b>cd /lib/x86_64-linux-gnu<br>sudo ln -s libreadline.so.7 libreadline.so.6</b><br><br>";
	s+="&#8226;&nbsp;&nbsp;Decompress the downloaded APBS package to a known location in order to locate the following executable binaries:<br><br>";
	s+="<b>APBS-1.5-linux64/bin/apbs</b><br>";
	s+="<b>APBS-1.5-linux64/share/apbs/tools/bin/multivalue</b><br>";	

	QLabel* tLabel=new QLabel;tLabel->setText(s.c_str());
	tLabel->setWordWrap(true);
	tLabel->setFixedWidth(width);
	tLabel->setTextInteractionFlags(Qt::TextBrowserInteraction);tLabel->setOpenExternalLinks(true);tLabel->setAlignment(Qt::AlignLeft);
	
	QGridLayout *mLayout = new QGridLayout;
	mLayout->setSizeConstraint(QLayout::SetFixedSize);
	mLayout->addWidget(tLabel);
	setLayout(mLayout);
	setWindowTitle(tr("APBS Installation Information"));
}

hydroproInstallInfo::hydroproInstallInfo(QWidget *parent) : QWidget(parent)
{	
	int width=500;
	string s="<b><u>Ubuntu</u></b><br>";
	s+="&#8226;&nbsp;&nbsp;Download HYDROPRO 10 at the following link: ";
	s+="<a href=\"http://leonardo.inf.um.es/macromol/programs/hydropro/hydropro10.zip\"><span style=\"color:blue;\">download</span></a><br><br>";
	s+="&#8226;&nbsp;&nbsp;After downloading the compressed HYDROPRO package, decompress it, and the HYDROPRO executable file (below) should be ready<br><br>";
	s+="<b>hydropro10-lnx.exe</b><br>";

	QLabel* tLabel=new QLabel;tLabel->setText(s.c_str());
	tLabel->setWordWrap(true);
	tLabel->setFixedWidth(width);
	tLabel->setTextInteractionFlags(Qt::TextBrowserInteraction);tLabel->setOpenExternalLinks(true);tLabel->setAlignment(Qt::AlignLeft);
	
	QGridLayout *mLayout = new QGridLayout;
	mLayout->setSizeConstraint(QLayout::SetFixedSize);
	mLayout->addWidget(tLabel);
	setLayout(mLayout);
	setWindowTitle(tr("HYDROPRO Installation Information"));
}

msmsInstallInfo::msmsInstallInfo(QWidget *parent)
    : QWidget(parent)
{	
	int width=500;
	string s="<b><u>Ubuntu</u></b><br>";
	s+="&#8226;&nbsp;&nbsp;Download MSMS 2.6.1 at the following link: ";
	s+="<a href=\"http://mgltools.scripps.edu/downloads/tars/releases/MSMSRELEASE/REL2.6.1/msms_i86Linux2_2.6.1.tar.gz\"><span style=\"color:blue;\">download</span></a><br><br>";
	s+="&#8226;&nbsp;&nbsp;After downloading the compressed MSMS package, decompress it, and the following files should be ready:<br><br>";
	s+="<b>msms.i86Linux2.2.6.1</b> (MSMS executable)<br>";
	s+="<b>atmtypenumbers</b> file<br>";
	s+="<b>pdb_to_xyzr</b> file<br>";
	s+="<b>pdb_to_xyzrn</b> file<br>";

	QLabel* tLabel=new QLabel;tLabel->setText(s.c_str());
	tLabel->setWordWrap(true);
	tLabel->setFixedWidth(width);
	tLabel->setTextInteractionFlags(Qt::TextBrowserInteraction);tLabel->setOpenExternalLinks(true);tLabel->setAlignment(Qt::AlignLeft);
	
	QGridLayout *mLayout = new QGridLayout;
	mLayout->setSizeConstraint(QLayout::SetFixedSize);
	mLayout->addWidget(tLabel);
	setLayout(mLayout);
	setWindowTitle(tr("MSMS Installation Information"));
}

pInstallInfo::pInstallInfo(QWidget *parent) : QWidget(parent)
{	
	int width=700;
	string s="<b><u>Ubuntu</u></b><br>";
	s+="&#8226;&nbsp;&nbsp;Download PDB2PQR 1.8 (and PROPKA 3.0) at the following link: ";
	s+="<a href=\"https://sourceforge.net/projects/pdb2pqr/files/pdb2pqr/pdb2pqr-1.8/pdb2pqr-1.8.tar.gz\"><span style=\"color:blue;\">download</span></a><br><br>";
	s+="&#8226;&nbsp;&nbsp;Python and Numpy must be installed before installation and can be obtained with the terminal commands below:<br><br>";
	s+="<b>sudo apt-get install python-dev<br>sudo apt-get install python-pip<br>pip install --upgrade pip<br>sudo pip install numpy scipy</b><br><br>";
	s+="&#8226;&nbsp;&nbsp;When configuring PDB2PQR, configure without pdb2pka using the following command:<br><br>";
	s+="<b>./configure --disable-pdb2pka</b><br>";
	s+="<b>sudo make</b><br>";
	s+="<b>sudo make install</b><br><br>";
	s+="&#8226;&nbsp;&nbsp;After installation, the following python file should be ready:<br><br>";
	s+="<b>main.py</b><br>";

	QLabel* tLabel=new QLabel;tLabel->setText(s.c_str());
	tLabel->setWordWrap(true);
	tLabel->setFixedWidth(width);
	tLabel->setTextInteractionFlags(Qt::TextBrowserInteraction);tLabel->setOpenExternalLinks(true);tLabel->setAlignment(Qt::AlignLeft);
	
	QGridLayout *mLayout = new QGridLayout;
	mLayout->setSizeConstraint(QLayout::SetFixedSize);
	mLayout->addWidget(tLabel);
	setLayout(mLayout);
	setWindowTitle(tr("PDB2PQR Installation Information"));
}

zpredInstallInfo::zpredInstallInfo(QWidget *parent) : QWidget(parent)
{	
	int width=700;
	QString Fldr=QDir::currentPath();
	string sFldr=Fldr.toStdString()+"/";
	// Write ZPRED.cpp file
	//string zpredFile=sFldr+"Software/ZPRED/C++/ZPRED.cpp";

	string s="<b><u>Ubuntu</u></b><br>";
	s+="&#8226;&nbsp;&nbsp;Download necessary ZPRED dependencies:<br>";
	s+="&#8226;&nbsp;&nbsp;<a href=\"https://www.boost.org/\">Boost</a> system and filesystem libraries can be obtained with the terminal command:<br><br>";
	s+="<b>sudo apt-get install libboost-all-dev</b><br><br>";
	s+="&#8226;&nbsp;&nbsp;<a href=\"https://gcc.gnu.org/\">g++</a> can be obtained with the terminal command:<br><br>";
	s+="<b>sudo apt-get install g++</b><br><br>";
	s+="&#8226;&nbsp;&nbsp;<a href=\"https://wiki.openssl.org/\">OpenSSL</a> Crypto library can be obtained with the terminal command:<br><br>";
	s+="<b>sudo apt-get install libssl-dev</b><br><br>";
	s+="&#8226;&nbsp;&nbsp;POSIX Pthread library can be obtained with the terminal command:<br><br>";
	s+="<b>sudo apt-get install libpthread-stubs0-dev</b><br><br>";
	//s+="&#8226;&nbsp;&nbsp;Download ZPRED source code: <a href=\""+zpredFile+"\"><span style=\"color:blue;\">download</span></a><br><br>";
	//s+="&#8226;&nbsp;&nbsp;Compile the ZPRED executable (ZPRED.exe) with the terminal command:<br><br>";
	//s+="<b>g++ -std=gnu++11 -o ZPRED.exe ZPRED.cpp -lboost_system -lboost_filesystem -lcrypto -lpthread</b><br><br>";
	s+="&#8226;&nbsp;&nbsp;ZPRED C++ source code is included with the user interface.  It can be located in Software/ZPRED/C++/ folder along wih the shell script, \"CompileZPRED.sh\". To compile the code, please first install the above dependencies, then compile the ZPRED.exe executable using the included shell script, \"CompileZPRED.sh\"<br><br>";

	QLabel* tLabel=new QLabel;tLabel->setText(s.c_str());
	tLabel->setWordWrap(true);
	tLabel->setFixedWidth(width);
	tLabel->setTextInteractionFlags(Qt::TextBrowserInteraction);tLabel->setOpenExternalLinks(true);tLabel->setAlignment(Qt::AlignLeft);
	
	QGridLayout *mLayout = new QGridLayout;
	mLayout->setSizeConstraint(QLayout::SetFixedSize);
	mLayout->addWidget(tLabel);
	setLayout(mLayout);
	setWindowTitle(tr("ZPRED Installation Information"));
}

zpredProgressDisplay::zpredProgressDisplay(QWidget *parent) : QWidget(parent)
{
	int width=1000;
	int height=500;
/*
	string msgFile=displayFldr+".gui_display.txt";
	ifstream fIn;
	fIn.open(msgFile.c_str());
	if(fIn.fail()){cerr<<"Error in zpredProgressDisplay!\nCould not open ZPRED gui_display file ("<<msgFile<<")\n";}//exit(EXIT_FAILURE);}
	int Sz=150000;
	char Val[Sz];
	fIn.getline(Val,Sz);
	string tmp=Val;fIn.close();
	displayLabel=new QLabel;displayLabel->setText(tmp.c_str());
	displayLabel->setWordWrap(true);
	displayLabel->setFixedWidth(width);
*/	
	displayLabel=new QLabel;displayLabel->setText("Initializing...");
	displayLabel->setFixedHeight(height);	
	displayLabel->setFixedWidth(width);

	QPushButton* uButton=new QPushButton(tr("UPDATE"));//sButton->setToolTip(tr("Specify the file path of the ZPRED executable file."));
	uButton->setFixedWidth(200);
	connect(uButton,SIGNAL(clicked()),this,SLOT(updateDisplay()));

	QPushButton* sButton=new QPushButton(tr("STOP"));//sButton->setToolTip(tr("Specify the file path of the ZPRED executable file."));
	sButton->setFixedWidth(200);
	connect(sButton,SIGNAL(clicked()),this,SLOT(killZPRED()));
	
	QGridLayout *mLayout = new QGridLayout;
	mLayout->setSizeConstraint(QLayout::SetFixedSize);
	mLayout->addWidget(displayLabel,0,0,1,3);
	mLayout->addWidget(uButton,1,0,1,1);mLayout->addWidget(sButton,1,2,1,1);
	setLayout(mLayout);
	setWindowTitle(tr("ZPRED Computations In Progress..."));
}

void zpredProgressDisplay::updateDisplay()
{
	int width=1000;
	int height=500;

	string msgFile=displayFldr+".gui_display.txt";
	ifstream fIn;		
	int Sz=150000;
	char Val[Sz];
	string tmp;
	fIn.open(msgFile.c_str());
	if(fIn.fail()){}//cerr<<"Error in zpredProgressDisplay!\nCould not open ZPRED gui_display file ("<<msgFile<<")\n";}//exit(EXIT_FAILURE);}
	else
		{fIn.getline(Val,Sz);
		tmp=Val;fIn.close();
		displayLabel->setText(tmp.c_str());
		displayLabel->setWordWrap(true);
		displayLabel->setFixedWidth(width);
		displayLabel->setFixedHeight(height);}
}

void zpredProgressDisplay::killZPRED()
{
	kill(PID,SIGTERM);
}

// makeGif UI
gUI::gUI(QWidget *parent) : QWidget(parent)
	{// Bold Font
	QFont font;
	font.setBold(true);

	QGroupBox* theBox=new QGroupBox;
	//gBox->setStyleSheet("QGroupBox { font-size:12pt; font-weight:bold; color:blue; } ");//gBox->setStyleSheet("QGroupBox { font-weight:bold; } ");	
	theBox->setFont(font);
	theBox->setAlignment(Qt::AlignRight);
	// Label
	pdbLabel=new QLabel;
	pdbLabel->setText("<b><span style=\"color:red;\">Select PDB file(s)</span>:</b>");
	pdbLabel->setAlignment(Qt::AlignLeft);
	// CheckBoxes
	xRotationBox=new QCheckBox(tr("Rotate Along X-Axis?")); xRotationBox->setChecked(false); xRotationBox->setFont(font);
	yRotationBox=new QCheckBox(tr("Rotate Along Y-Axis?")); yRotationBox->setChecked(true); yRotationBox->setFont(font);
	zRotationBox=new QCheckBox(tr("Rotate Along Z-Axis?")); zRotationBox->setChecked(false); zRotationBox->setFont(font);
	// Button
	QPushButton* pdbButton=new QPushButton(tr("..."));	
	pdbButton->setToolTip(tr("Tooltip"));
	pdbButton->setFixedWidth(40);
	pdbButton->setFont(font);
	connect(pdbButton,SIGNAL(clicked()),this,SLOT(open_pdbFiles()));

	QGridLayout *layout = new QGridLayout;
	layout->setSizeConstraint(QLayout::SetFixedSize);
	layout->addWidget(pdbLabel,0,0,1,1);layout->addWidget(pdbButton,0,1,1,1);
	layout->addWidget(xRotationBox,1,0,1,2);
	layout->addWidget(yRotationBox,2,0,1,2);
	layout->addWidget(zRotationBox,3,0,1,2);
	theBox->setLayout(layout);

	QPushButton *xButton=new QPushButton(tr("Make Gif"));
	xButton->setFont(font);
	
	connect(xButton,SIGNAL(clicked()),this,SLOT(makeGif()));

	QGroupBox* botRight=new QGroupBox;
	QGridLayout *brlayout=new QGridLayout;
	//brlayout->setSizeConstraint(QLayout::SetFixedSize);
	brlayout->addWidget(xButton,2,0,1,1);
	botRight->setLayout(brlayout);

	QGridLayout *mLayout = new QGridLayout;
	mLayout->setSizeConstraint(QLayout::SetFixedSize);
	mLayout->addWidget(theBox,0,0);
	mLayout->addWidget(botRight,0,1);
	setLayout(mLayout);
	setWindowTitle(tr("Make Gif"));
	resize(QDesktopWidget().availableGeometry(this).size());
}

void gUI::makeGif()
	{QString theCurFldr=QDir::currentPath();
	string sFldr=theCurFldr.toStdString()+"/PyMol_Images/";
	// Make PyMol Images Folder
	boost::filesystem::path dir(sFldr.c_str());
	boost::filesystem::create_directory(dir);

	string pythonFunctionName="Test";
	string pythonFunctionFile=sFldr+"Test.py";
	string cmd="pymol -d run "+pythonFunctionFile+" -d "+pythonFunctionName;
	//string cmd2="cd "+sFldr+"PyMol_Images/";
	string cmd2="convert PyMol_Images/*.png Output.gif";
	int statusCode,status;
	bool xROT,yROT,zROT;
	if(xRotationBox->isChecked()){xROT=true;}
	else{xROT=false;}
	if(yRotationBox->isChecked()){yROT=true;}
	else{yROT=false;}
	if(zRotationBox->isChecked()){zROT=true;}
	else{zROT=false;}
	//writeZPREDInputFile(zInputFile);
	write_pymolFunction_File(xROT,yROT,zROT,pythonFunctionName,pythonFunctionFile);		
	cerr<<"Creating stack of png files with PyMol...";
	statusCode=system(cmd.c_str());
	cerr<<"Done\nCreating gif file (Output.gif)...";
	statusCode=system(cmd2.c_str());
	cerr<<"Done\n";
	//cerr<<"Status: "<<statusCode<<endl;
	// Delete folder and contained files
	boost::filesystem::remove_all(sFldr);
	}

void gUI::open_pdbFiles()
	{QString Fldr=QDir::currentPath();
	string sFldr=Fldr.toStdString()+"/";
	QStringList filenames = QFileDialog::getOpenFileNames(this,tr("files"),tr(sFldr.c_str()),tr("All files (*)") );
	string tmp;
	if(filenames.count()!=0)
		{numPDBs=filenames.count();
		pdbFiles=new string[numPDBs];
		for(int i=0;i<numPDBs;i++)
			{tmp=filenames.at(i).toLocal8Bit().constData();
			pdbFiles[i]=tmp;}
		pdbLabel->setText("<b><span style=\"color:black;\">Select PDB file(s)</span>:</b>");}
	}

void gUI::write_pymolFunction_File(bool xROT,bool yROT,bool zROT,string pyFnNm,string outFile)
	{QString theCurFldr=QDir::currentPath();
	string sFldr=theCurFldr.toStdString()+"/PyMol_Images/";
	string Output="";
	// Write PyMol Function File
	Output+="#!/usr/bin/python\nfrom pymol import stored\nfrom time import sleep\n\n# Function Definition\n\n";
	Output+="def "+pyFnNm+"():\n\n";
	// Define PDB File(s)	
	for(int i=0;i<numPDBs;i++)
		{Output+="\tpFile"+cnvrtNumToStrng(i+1,0)+"=\""+pdbFiles[i]+"\"\n";}
	Output+="\n";
	// Load PDB File(s)
	for(int i=0;i<numPDBs;i++)
		{Output+="\tcmd.load(pFile"+cnvrtNumToStrng(i+1,0)+")\n";}
	Output+="\n";
	// More Functions
	Output+="\tcmd.do(\"util.cbag\")\n";
	Output+="\tcmd.zoom(\"all\",10)\n";
	
	Output+="\tfor a in range(0,360,1):\n";
	if(xROT){Output+="\t\tcmd.turn(\"x\",1)\n";}
	if(yROT){Output+="\t\tcmd.turn(\"y\",1)\n";}
	if(zROT){Output+="\t\tcmd.turn(\"z\",1)\n";}
	Output+="\t\tif a < 10:\n";
	Output+="\t\t\tpngFile=\""+sFldr+"\"+\"0000\"+str(a)+\".png\"\n";
	Output+="\t\telif a > 9 and a < 100:\n";
	Output+="\t\t\tpngFile=\""+sFldr+"\"+\"000\"+str(a)+\".png\"\n";
	Output+="\t\telif a > 99 and a < 1000:\n";
	Output+="\t\t\tpngFile=\""+sFldr+"\"+\"00\"+str(a)+\".png\"\n";
	Output+="\t\telif a > 999 and a < 10000:\n";
	Output+="\t\t\tpngFile=\""+sFldr+"\"+\"0\"+str(a)+\".png\"\n";
	Output+="\t\telse:\n";
	Output+="\t\t\tpngFile=\""+sFldr+"\"+str(a)+\".png\"\n";

	Output+="\t\tcmd.ray()\n";
	Output+="\t\tcmd.png(pngFile)\n";
	Output+="\n";
	Output+="\tcmd.quit()\n";
	Output+="\n";
	// End File
	Output+="\treturn\n\ncmd.extend(\""+pyFnNm+"\","+pyFnNm+")";

	ofstream fOut;
	fOut.open(outFile.c_str());
	if(fOut.fail()){cerr<<"Error in write_pymolFunction_File!\nCould not open file ("<<outFile<<")\n";}//exit(EXIT_FAILURE);}
	else
		{fOut<<Output;fOut.close();}
	}
