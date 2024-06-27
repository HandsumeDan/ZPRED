# include "mobility.h"


// INPUT: Zeta Potential (ZP) [Volts]
// INPUT: Solvated Radius (a) [meters]
// INPUT: Solution Structure Array
double two_ion_problem(double ZP,double a,Solution S,Error E,bool VERBOSE)
	{// Convert Debye Length from Angstroms to Meters
	double Debye=S.debyeLength*1e-10;
	// Final Position: Solvated Radius [meters]
	double rf=a;
	double r0;
	double dr=-1e-15;	// integrate from Inf to Protein Radius
	if(Debye>a){r0=5*Debye;}
	else{r0=5*a;}
		
	if(VERBOSE){cout<<"2 ION PROBLEM\nInitial Position:\t"<<r0<<" m"<<"\nFinal Position:\t\t"<<rf<<" m"<<"\nStep Size:\t\t"<<dr<<" m"<<endl;}

	// Calculate Poisson Boltzmann Coefficient and Determine Initial Position away from Protein
	double C=calc_PBE_coefficient(ZP,r0,rf,dr,S,E,VERBOSE);
	if(VERBOSE){cout<<"PBE Coeff: "<<C<<endl;}
	if(isnan(C)){return nan("");}
	int N=22;
	state_type y(N);
	// F Function
	// y[0] = f
	// y[1] = df/dr
	// y[2] = d2f/dr2
	// y[3] = d3f/dr3
	// y[4] = d4f/dr4
	// Equilibrium Electric Pontential
	// y[5] = potential
	// y[6] = d(potential)/dr
	// Induced Electrochemical Potential Function (2nd ion)
	// y[7] = IEP1
	// y[8] = d(IEP1)/dr
	// y[9] = IEP2
	// y[10] = d(IEP2)/dr	
	// Constants
	// y[11] = Temp					// [Kelvin]
	// y[12] = numIon
	// Ion Properties
	// y[13] = Anion.val;		
	// y[14] = Anion.conc;		// mM [mol/m^3]
	// y[15] = Anion.radius;		// meters
	// y[16] = Cation.val;
	// y[17] = Cation.conc;
	// y[18] = Cation.radius;
	// Solution Properties
	// y[19] = relative dielectric	[dimensionless]
	// y[20] = Debye Length [meters]
	// y[21] = viscosity  [Pa s]

	// Collect Error Data
	double abs_err=E.absolute;//1e-12;
	double rel_err=E.relative;//1e-7;
	double a_x=E.a_x;//1;
	double a_dxdt=E.a_dxdt;//1;
	controlled_stepper_type controlled_stepper(default_error_checker<double,range_algebra,default_operations>(abs_err,rel_err,a_x,a_dxdt));	

	// Calculate Mobility
	vector<state_type> y1,y2,y3,y4,y5,y6;
	vector<double> r1,r2,r3,r4,r5,r6;
	// Specify Boundary Conditions / Initial Values Infinitely Far from the Particle
	if(VERBOSE){cout<<"Solving Homogeneous Problem for 1st Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=1/(r0*r0);y[8]=-2/(r0*r0*r0);
	// Second Ion
	y[9]=0;y[10]=0;
	// Constants
	y[11]=S.temperature;y[12]=S.numIon;
	// Fist Ion Properties
	y[13]=S.ionVal[0];y[14]=S.ionConc[0]*1000;y[15]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[16]=S.ionVal[1];y[17]=S.ionConc[1]*1000;y[18]=S.ionRad[1]*1e-10;
	// Solution Properties
	y[19]=S.dielectric;y[20]=Debye;y[21]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps1=integrate_adaptive(controlled_stepper,homogeneous_form_2ion,y,r0,rf,dr,push_back_state_and_time(y1,r1));
	int Sz1=steps1+1;
	if(VERBOSE){cout<<"Done! ("<<Sz1<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 2nd Ion...";}	
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=1/(r0*r0);y[10]=-2/(r0*r0*r0);
	// Constants
	y[11]=S.temperature;y[12]=S.numIon;
	// Fist Ion Properties
	y[13]=S.ionVal[0];y[14]=S.ionConc[0]*1000;y[15]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[16]=S.ionVal[1];y[17]=S.ionConc[1]*1000;y[18]=S.ionRad[1]*1e-10;
	// Solution Properties
	y[19]=S.dielectric;y[20]=Debye;y[21]=S.viscosity;
	// Solve Homogeneous Form for Second Ion
	size_t steps2=integrate_adaptive(controlled_stepper,homogeneous_form_2ion,y,r0,rf,dr,push_back_state_and_time(y2,r2));
	int Sz2=steps2+1;
	if(VERBOSE){cout<<"Done! ("<<Sz2<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 1st f term...";}
	// F Function
	y[0]=r0;y[1]=1;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Constants
	y[11]=S.temperature;y[12]=S.numIon;
	// Fist Ion Properties
	y[13]=S.ionVal[0];y[14]=S.ionConc[0]*1000;y[15]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[16]=S.ionVal[1];y[17]=S.ionConc[1]*1000;y[18]=S.ionRad[1]*1e-10;
	// Solution Properties
	y[19]=S.dielectric;y[20]=Debye;y[21]=S.viscosity;
	// Solve Homogeneous Form for 1st term of F function
	size_t steps3=integrate_adaptive(controlled_stepper,homogeneous_form_2ion,y,r0,rf,dr,push_back_state_and_time(y3,r3));
	int Sz3=steps3+1;
	if(VERBOSE){cout<<"Done! ("<<Sz3<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 2nd f term...";}
	// F Function
	y[0]=1.0/r0;y[1]=-1.0/(r0*r0);y[2]=2.0/(r0*r0*r0);y[3]=-6.0/(r0*r0*r0*r0);y[4]=24.0/(r0*r0*r0*r0*r0);
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Constants
	y[11]=S.temperature;y[12]=S.numIon;
	// Fist Ion Properties
	y[13]=S.ionVal[0];y[14]=S.ionConc[0]*1000;y[15]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[16]=S.ionVal[1];y[17]=S.ionConc[1]*1000;y[18]=S.ionRad[1]*1e-10;
	// Solution Properties
	y[19]=S.dielectric;y[20]=Debye;y[21]=S.viscosity;
	// Solve Homogeneous Form for 2nd term of F function
	size_t steps4=integrate_adaptive(controlled_stepper,homogeneous_form_2ion,y,r0,rf,dr,push_back_state_and_time(y4,r4));
	int Sz4=steps4+1;
	if(VERBOSE){cout<<"Done! ("<<Sz4<<" points)"<<endl;}	

	if(VERBOSE){cout<<"Solving Flow Field Problem...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Constants
	y[11]=S.temperature;y[12]=S.numIon;
	// Fist Ion Properties
	y[13]=S.ionVal[0];y[14]=S.ionConc[0]*1000;y[15]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[16]=S.ionVal[1];y[17]=S.ionConc[1]*1000;y[18]=S.ionRad[1]*1e-10;
	// Solution Properties
	y[19]=S.dielectric;y[20]=Debye;y[21]=S.viscosity;
	// Solve Flow Field Perturbation 
	size_t steps5=integrate_adaptive(controlled_stepper,inhomogeneous_Prob1_2ion,y,r0,rf,dr,push_back_state_and_time(y5,r5));
	int Sz5=steps5+1;
	if(VERBOSE){cout<<"Done! ("<<Sz5<<" points)"<<endl;}

	if(VERBOSE){cout<<"Calculating Asymptotic Constants...";}
	alglib::real_2d_array A;
	alglib::real_1d_array C1,B1;
	// Define System Matrix containing coefficients of linear system
	A.setlength(4,4);
	for(int i=0;i<4;i++)
		{for(int j=0;j<4;j++)
			{A[i][j]=0;}}
	// d(IEP1)/dr
	A[0][0]=y1[Sz1-1][8];A[0][1]=y2[Sz2-1][8];A[0][2]=y3[Sz3-1][8];A[0][3]=y4[Sz4-1][8];
	// d(IEP2)/dr
	A[1][0]=y1[Sz1-1][10];A[1][1]=y2[Sz2-1][10];A[1][2]=y3[Sz3-1][10];A[1][3]=y4[Sz4-1][10];
	// f 1st (df/dr)
	A[2][0]=y1[Sz1-1][1];A[2][1]=y2[Sz2-1][1];A[2][2]=y3[Sz3-1][1];A[2][3]=y4[Sz4-1][1];
	// f 2nd (d2f/dr2)
	A[3][0]=y1[Sz1-1][2];A[3][1]=y2[Sz2-1][2];A[3][2]=y3[Sz3-1][2];A[3][3]=y4[Sz4-1][2];
	// Define Solution Matrix for Problem #1
	B1.setlength(4);
	// d(IEP1)/dr = 0
	B1[0]=0-y5[Sz5-1][8];
	// d(IEP2)/dr = 0
	B1[1]=0-y5[Sz5-1][10];
	// df/dr = -Radius/2
	B1[2]=-a/2 - y5[Sz5-1][1];
	// d2f/dr2 = -0.5
	B1[3]=-0.5 - y5[Sz5-1][2];
	// Create Output Parameters
	alglib::ae_int_t info;
	alglib::densesolverreport rep;
	alglib::rmatrixsolve(A,4,B1,info,rep,C1);
	if(VERBOSE)
		{cout<<"\n";
		cout<<C1[0]*A[0][0] + C1[1]*A[0][1] + C1[2]*A[0][2] + C1[3]*A[0][3] + y5[Sz5-1][8]<<endl;
		cout<<C1[0]*A[1][0] + C1[1]*A[1][1] + C1[2]*A[1][2] + C1[3]*A[1][3] + y5[Sz5-1][10]<<endl;
		cout<<C1[0]*A[2][0] + C1[1]*A[2][1] + C1[2]*A[2][2] + C1[3]*A[2][3] + y5[Sz5-1][1]<<endl;
		cout<<C1[0]*A[3][0] + C1[1]*A[3][1] + C1[2]*A[3][2] + C1[3]*A[3][3] + y5[Sz5-1][2]<<endl;
		cout<<"Done!"<<endl;
		cout<<"C1_1 = "<<C1[0]<<endl;
		cout<<"C2_1 = "<<C1[1]<<endl;
		cout<<"C3_1 = "<<C1[2]<<endl;
		cout<<"C4_1 = "<<C1[3]<<endl;}

	if(VERBOSE){cout<<"Solving Electric Field Problem...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Constants
	y[11]=S.temperature;y[12]=S.numIon;
	// Fist Ion Properties
	y[13]=S.ionVal[0];y[14]=S.ionConc[0]*1000;y[15]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[16]=S.ionVal[1];y[17]=S.ionConc[1]*1000;y[18]=S.ionRad[1]*1e-10;
	// Solution Properties
	y[19]=S.dielectric;y[20]=Debye;y[21]=S.viscosity;
	// Solve Electric Field Perturbation
	size_t steps6=integrate_adaptive(controlled_stepper,inhomogeneous_Prob2_2ion,y,r0,rf,dr,push_back_state_and_time(y6,r6));
	int Sz6=steps6+1;
	if(VERBOSE){cout<<"Done! ("<<Sz6<<" points)"<<endl;}
	
	if(VERBOSE){cout<<"Calculating Asymptotic Constants...";}
	alglib::real_1d_array C2,B2;
	// d(IEP1)/dr
	A[0][0]=y1[Sz1-1][8];A[0][1]=y2[Sz2-1][8];A[0][2]=y3[Sz3-1][8];A[0][3]=y4[Sz4-1][8];
	// d(IEP2)/dr
	A[1][0]=y1[Sz1-1][10];A[1][1]=y2[Sz2-1][10];A[1][2]=y3[Sz3-1][10];A[1][3]=y4[Sz4-1][10];
	// f 1st (df/dr)
	A[2][0]=y1[Sz1-1][1];A[2][1]=y2[Sz2-1][1];A[2][2]=y3[Sz3-1][1];A[2][3]=y4[Sz4-1][1];
	// f 2nd (d2f/dr2)
	A[3][0]=y1[Sz1-1][2];A[3][1]=y2[Sz2-1][2];A[3][2]=y3[Sz3-1][2];A[3][3]=y4[Sz4-1][2];
	// Define Solution Matrix for Problem #2
	B2.setlength(4);
	// d(IEP1)/dr = -1
	B2[0]=-1-y6[Sz6-1][8];
	// d(IEP2)/dr = -1
	B2[1]=-1-y6[Sz6-1][10];
	// df/dr = 0
	B2[2]=0-y6[Sz6-1][1];
	// d2f/dr2 = 0
	B2[3]=0-y6[Sz6-1][2];
	// Create Output Parameters
	alglib::rmatrixsolve(A,4,B2,info,rep,C2);

	if(VERBOSE)
		{cout<<"\n";
		cout<<C2[0]*A[0][0] + C2[1]*A[0][1] + C2[2]*A[0][2] + C2[3]*A[0][3] + y6[Sz6-1][8]<<endl;
		cout<<C2[0]*A[1][0] + C2[1]*A[1][1] + C2[2]*A[1][2] + C2[3]*A[1][3] + y6[Sz6-1][10]<<endl;
		cout<<C2[0]*A[2][0] + C2[1]*A[2][1] + C2[2]*A[2][2] + C2[3]*A[2][3] + y6[Sz6-1][1]<<endl;
		cout<<C2[0]*A[3][0] + C2[1]*A[3][1] + C2[2]*A[3][2] + C2[3]*A[3][3] + y6[Sz6-1][2]<<endl;
		cout<<"Done!"<<endl;
		cout<<"C1_2 = "<<C2[0]<<endl;
		cout<<"C2_2 = "<<C2[1]<<endl;
		cout<<"C3_2 = "<<C2[2]<<endl;
		cout<<"C4_2 = "<<C2[3]<<endl;}

	double mobility=-C2[2]/C1[2];
	
	if(VERBOSE)
		{cout<<"\n\nElectrophoretic Mobility:\t"<<mobility<<" m^2/Vs\n";
		cout<<"ka:\t\t"<<a/Debye<<endl;
		cout<<"Zeta Potential:\t"<<ZP<<" V\n";
		cout<<"E:\t\t"<<3*S.viscosity*ec*mobility/(2*eo*S.dielectric*kb*S.temperature)<<endl;
		cout<<"y:\t\t"<<ec*ZP/(kb*S.temperature)<<endl;}

	return mobility;}

// INPUT: Zeta Potential (ZP) [Volts]
// INPUT: Solvated Radius (a) [meters]
// INPUT: Solution Structure Array
double three_ion_problem(double ZP,double a,Solution S,Error E,bool VERBOSE)
	{// Convert Debye Length from Angstroms to Meters
	double Debye=S.debyeLength*1e-10;
	// Final Position: Solvated Radius [meters]
	double rf=a;
	double r0;
	double dr=-1e-15;	// integrate from Inf to Protein Radius
	if(Debye>a){r0=5*Debye;}
	else{r0=5*a;}
		
	if(VERBOSE){cout<<"3 ION PROBLEM\nInitial Position:\t"<<r0<<" m"<<"\nFinal Position:\t\t"<<rf<<" m"<<"\nStep Size:\t\t"<<dr<<" m"<<endl;}

	// Calculate Poisson Boltzmann Coefficient and Determine Initial Position away from Protein
	double C=calc_PBE_coefficient(ZP,r0,rf,dr,S,E,VERBOSE);
	if(VERBOSE){cout<<"PBE Coeff: "<<C<<endl;}
	if(isnan(C)){return nan("");}
	int N=27;
	state_type y(N);
	// F Function
	// y[0] = f
	// y[1] = df/dr
	// y[2] = d2f/dr2
	// y[3] = d3f/dr3
	// y[4] = d4f/dr4
	// Equilibrium Electric Pontential
	// y[5] = potential
	// y[6] = d(potential)/dr
	// Induced Electrochemical Potential Function (3 ion)
	// y[7] = IEP1
	// y[8] = d(IEP1)/dr
	// y[9] = IEP2
	// y[10] = d(IEP2)/dr
	// y[11] = IEP3
	// y[12] = d(IEP3)/dr
	// Constants
	// y[13] = Temp					// [Kelvin]
	// y[14] = numIon
	// Ion Properties
	// y[15] = Ion1.val;		
	// y[16] = Ion1.conc;		// mM [mol/m^3]
	// y[17] = Ion1.radius;		// meters
	// y[18] = Ion2.val;
	// y[19] = Ion2.conc;
	// y[20] = Ion2.radius;
	// y[21] = Ion3.val;
	// y[22] = Ion3.conc;
	// y[23] = Ion3.radius;
	// Solution Properties
	// y[24] = relative dielectric	[dimensionless]
	// y[25] = Debye Length [meters]
	// y[26] = viscosity  [Pa s]

	// Collect Error Data
	double abs_err=E.absolute;
	double rel_err=E.relative;
	double a_x=E.a_x;
	double a_dxdt=E.a_dxdt;
	controlled_stepper_type controlled_stepper(default_error_checker<double,range_algebra,default_operations>(abs_err,rel_err,a_x,a_dxdt));	

	// Calculate Mobility
	vector<state_type> y1,y2,y3,y4,y5,y6,y7;
	vector<double> r1,r2,r3,r4,r5,r6,r7;
	// Specify Boundary Conditions / Initial Values Infinitely Far from the Particle
	if(VERBOSE){cout<<"Solving Homogeneous Problem for 1st Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=1/(r0*r0);y[8]=-2/(r0*r0*r0);
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Constants
	y[13]=S.temperature;y[14]=S.numIon;
	// Fist Ion Properties
	y[15]=S.ionVal[0];y[16]=S.ionConc[0]*1000;y[17]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[18]=S.ionVal[1];y[19]=S.ionConc[1]*1000;y[20]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[21]=S.ionVal[2];y[22]=S.ionConc[2]*1000;y[23]=S.ionRad[2]*1e-10;
	// Solution Properties
	y[24]=S.dielectric;y[25]=Debye;y[26]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps1=integrate_adaptive(controlled_stepper,homogeneous_form_3ion,y,r0,rf,dr,push_back_state_and_time(y1,r1));
	int Sz1=steps1+1;
	if(VERBOSE){cout<<"Done! ("<<Sz1<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 2nd Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=1/(r0*r0);y[10]=-2/(r0*r0*r0);
	// Third Ion
	y[11]=0;y[12]=0;
	// Constants
	y[13]=S.temperature;y[14]=S.numIon;
	// Fist Ion Properties
	y[15]=S.ionVal[0];y[16]=S.ionConc[0]*1000;y[17]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[18]=S.ionVal[1];y[19]=S.ionConc[1]*1000;y[20]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[21]=S.ionVal[2];y[22]=S.ionConc[2]*1000;y[23]=S.ionRad[2]*1e-10;
	// Solution Properties
	y[24]=S.dielectric;y[25]=Debye;y[26]=S.viscosity;
	// Solve Homogeneous Form for Second Ion
	size_t steps2=integrate_adaptive(controlled_stepper,homogeneous_form_3ion,y,r0,rf,dr,push_back_state_and_time(y2,r2));
	int Sz2=steps2+1;
	if(VERBOSE){cout<<"Done! ("<<Sz2<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 3rd Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=1/(r0*r0);y[12]=-2/(r0*r0*r0);
	// Constants
	y[13]=S.temperature;y[14]=S.numIon;
	// Fist Ion Properties
	y[15]=S.ionVal[0];y[16]=S.ionConc[0]*1000;y[17]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[18]=S.ionVal[1];y[19]=S.ionConc[1]*1000;y[20]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[21]=S.ionVal[2];y[22]=S.ionConc[2]*1000;y[23]=S.ionRad[2]*1e-10;
	// Solution Properties
	y[24]=S.dielectric;y[25]=Debye;y[26]=S.viscosity;
	// Solve Homogeneous Form for 3rd Ion
	size_t steps3=integrate_adaptive(controlled_stepper,homogeneous_form_3ion,y,r0,rf,dr,push_back_state_and_time(y3,r3));
	int Sz3=steps3+1;
	if(VERBOSE){cout<<"Done! ("<<Sz3<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 1st f term...";}
	// F Function
	y[0]=r0;y[1]=1;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Constants
	y[13]=S.temperature;y[14]=S.numIon;
	// Fist Ion Properties
	y[15]=S.ionVal[0];y[16]=S.ionConc[0]*1000;y[17]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[18]=S.ionVal[1];y[19]=S.ionConc[1]*1000;y[20]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[21]=S.ionVal[2];y[22]=S.ionConc[2]*1000;y[23]=S.ionRad[2]*1e-10;
	// Solution Properties
	y[24]=S.dielectric;y[25]=Debye;y[26]=S.viscosity;
	// Solve Homogeneous Form
	size_t steps4=integrate_adaptive(controlled_stepper,homogeneous_form_3ion,y,r0,rf,dr,push_back_state_and_time(y4,r4));
	int Sz4=steps4+1;
	if(VERBOSE){cout<<"Done! ("<<Sz4<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 2nd f term...";}
	// F Function
	y[0]=1/r0;y[1]=-1/(r0*r0);y[2]=2/(r0*r0*r0);y[3]=-6/(r0*r0*r0*r0);y[4]=24/(r0*r0*r0*r0*r0);
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Constants
	y[13]=S.temperature;y[14]=S.numIon;
	// Fist Ion Properties
	y[15]=S.ionVal[0];y[16]=S.ionConc[0]*1000;y[17]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[18]=S.ionVal[1];y[19]=S.ionConc[1]*1000;y[20]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[21]=S.ionVal[2];y[22]=S.ionConc[2]*1000;y[23]=S.ionRad[2]*1e-10;
	// Solution Properties
	y[24]=S.dielectric;y[25]=Debye;y[26]=S.viscosity;
	// Solve Homogeneous Form
	size_t steps5=integrate_adaptive(controlled_stepper,homogeneous_form_3ion,y,r0,rf,dr,push_back_state_and_time(y5,r5));
	int Sz5=steps5+1;
	if(VERBOSE){cout<<"Done! ("<<Sz5<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Flow Field Problem...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Constants
	y[13]=S.temperature;y[14]=S.numIon;
	// Fist Ion Properties
	y[15]=S.ionVal[0];y[16]=S.ionConc[0]*1000;y[17]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[18]=S.ionVal[1];y[19]=S.ionConc[1]*1000;y[20]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[21]=S.ionVal[2];y[22]=S.ionConc[2]*1000;y[23]=S.ionRad[2]*1e-10;
	// Solution Properties
	y[24]=S.dielectric;y[25]=Debye;y[26]=S.viscosity;
	// Solve Flow Field Perturbation 
	size_t steps6=integrate_adaptive(controlled_stepper,inhomogeneous_Prob1_3ion,y,r0,rf,dr,push_back_state_and_time(y6,r6));
	int Sz6=steps6+1;
	if(VERBOSE){cout<<"Done! ("<<Sz6<<" points)"<<endl;}

	if(VERBOSE){cout<<"Calculating Asymptotic Constants...";}
	alglib::real_2d_array A;
	alglib::real_1d_array C1,B1;
	// Define System Matrix containing coefficients of linear system
	A.setlength(5,5);
	for(int i=0;i<5;i++)
		{for(int j=0;j<5;j++)
			{A[i][j]=0;}}
	// d(IEP1)/dr
	A[0][0]=y1[Sz1-1][8];A[0][1]=y2[Sz2-1][8];A[0][2]=y3[Sz3-1][8];A[0][3]=y4[Sz4-1][8];A[0][4]=y5[Sz5-1][8];
	// d(IEP2)/dr
	A[1][0]=y1[Sz1-1][10];A[1][1]=y2[Sz2-1][10];A[1][2]=y3[Sz3-1][10];A[1][3]=y4[Sz4-1][10];A[1][4]=y5[Sz5-1][10];
	// d(IEP3)/dr
	A[2][0]=y1[Sz1-1][12];A[2][1]=y2[Sz2-1][12];A[2][2]=y3[Sz3-1][12];A[2][3]=y4[Sz4-1][12];A[2][4]=y5[Sz5-1][12];
	// f 1st (df/dr)
	A[3][0]=y1[Sz1-1][1];A[3][1]=y2[Sz2-1][1];A[3][2]=y3[Sz3-1][1];A[3][3]=y4[Sz4-1][1];A[3][4]=y5[Sz5-1][1];
	// f 2nd (d2f/dr2)
	A[4][0]=y1[Sz1-1][2];A[4][1]=y2[Sz2-1][2];A[4][2]=y3[Sz3-1][2];A[4][3]=y4[Sz4-1][2];A[4][4]=y5[Sz5-1][2];
	// Define Solution Matrix for Problem #1
	B1.setlength(5);
	// d(IEP1)/dr = 0
	B1[0]=0-y6[Sz6-1][8];
	// d(IEP2)/dr = 0
	B1[1]=0-y6[Sz6-1][10];
	// d(IEP3)/dr = 0
	B1[2]=0-y6[Sz6-1][12];
	// df/dr = -Radius/2
	B1[3]=-a/2 - y6[Sz6-1][1];
	// d2f/dr2 = -0.5
	B1[4]=-0.5 - y6[Sz6-1][2];
	// Create Output Parameters
	alglib::ae_int_t info;
	alglib::densesolverreport rep;
	alglib::rmatrixsolve(A,5,B1,info,rep,C1);
	if(VERBOSE)
		{cout<<"\n";
		cout<<C1[0]*A[0][0] + C1[1]*A[0][1] + C1[2]*A[0][2] + C1[3]*A[0][3] + C1[4]*A[0][4] + y6[Sz6-1][8]<<endl;
		cout<<C1[0]*A[1][0] + C1[1]*A[1][1] + C1[2]*A[1][2] + C1[3]*A[1][3] + C1[4]*A[1][4] + y6[Sz6-1][10]<<endl;
		cout<<C1[0]*A[2][0] + C1[1]*A[2][1] + C1[2]*A[2][2] + C1[3]*A[2][3] + C1[4]*A[2][4] + y6[Sz6-1][12]<<endl;
		cout<<C1[0]*A[3][0] + C1[1]*A[3][1] + C1[2]*A[3][2] + C1[3]*A[3][3] + C1[4]*A[3][4] + y6[Sz6-1][1]<<endl;
		cout<<C1[0]*A[4][0] + C1[1]*A[4][1] + C1[2]*A[4][2] + C1[3]*A[4][3] + C1[4]*A[4][4] + y6[Sz6-1][2]<<endl;
		cout<<"Done!"<<endl;
		cout<<"C1_1 = "<<C1[0]<<endl;
		cout<<"C2_1 = "<<C1[1]<<endl;
		cout<<"C3_1 = "<<C1[2]<<endl;
		cout<<"C4_1 = "<<C1[3]<<endl;
		cout<<"C5_1 = "<<C1[4]<<endl;}
	
	if(VERBOSE){cout<<"Solving Electric Field Problem...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Constants
	y[13]=S.temperature;y[14]=S.numIon;
	// Fist Ion Properties
	y[15]=S.ionVal[0];y[16]=S.ionConc[0]*1000;y[17]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[18]=S.ionVal[1];y[19]=S.ionConc[1]*1000;y[20]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[21]=S.ionVal[2];y[22]=S.ionConc[2]*1000;y[23]=S.ionRad[2]*1e-10;
	// Solution Properties
	y[24]=S.dielectric;y[25]=Debye;y[26]=S.viscosity;
	// Solve Electric Field Perturbation
	size_t steps7=integrate_adaptive(controlled_stepper,inhomogeneous_Prob2_3ion,y,r0,rf,dr,push_back_state_and_time(y7,r7));
	int Sz7=steps7+1;
	if(VERBOSE){cout<<"Done! ("<<Sz7<<" points)"<<endl;}
	
	if(VERBOSE){cout<<"Calculating Asymptotic Constants...";}
	alglib::real_1d_array C2,B2;
	// d(IEP1)/dr
	A[0][0]=y1[Sz1-1][8];A[0][1]=y2[Sz2-1][8];A[0][2]=y3[Sz3-1][8];A[0][3]=y4[Sz4-1][8];A[0][4]=y5[Sz5-1][8];
	// d(IEP2)/dr
	A[1][0]=y1[Sz1-1][10];A[1][1]=y2[Sz2-1][10];A[1][2]=y3[Sz3-1][10];A[1][3]=y4[Sz4-1][10];A[1][4]=y5[Sz5-1][10];
	// d(IEP3)/dr
	A[2][0]=y1[Sz1-1][12];A[2][1]=y2[Sz2-1][12];A[2][2]=y3[Sz3-1][12];A[2][3]=y4[Sz4-1][12];A[2][4]=y5[Sz5-1][12];
	// f 1st (df/dr)
	A[3][0]=y1[Sz1-1][1];A[3][1]=y2[Sz2-1][1];A[3][2]=y3[Sz3-1][1];A[3][3]=y4[Sz4-1][1];A[3][4]=y5[Sz5-1][1];
	// f 2nd (d2f/dr2)
	A[4][0]=y1[Sz1-1][2];A[4][1]=y2[Sz2-1][2];A[4][2]=y3[Sz3-1][2];A[4][3]=y4[Sz4-1][2];A[4][4]=y5[Sz5-1][2];
	// Define Solution Matrix for Problem #1
	B2.setlength(5);
	// d(IEP1)/dr = -1
	B2[0]=-1-y7[Sz7-1][8];
	// d(IEP2)/dr = -1
	B2[1]=-1-y7[Sz7-1][10];
	// d(IEP3)/dr = -1
	B2[2]=-1-y7[Sz7-1][12];
	// df/dr = 0
	B2[3]=0-y7[Sz7-1][1];
	// d2f/dr2 = 0
	B2[4]=0-y7[Sz7-1][2];
	// Create Output Parameters
	alglib::rmatrixsolve(A,5,B2,info,rep,C2);
	if(VERBOSE)
		{cout<<"\n";
		cout<<C2[0]*A[0][0] + C2[1]*A[0][1] + C2[2]*A[0][2] + C2[3]*A[0][3] + C2[4]*A[0][4] + y7[Sz7-1][8]<<endl;
		cout<<C2[0]*A[1][0] + C2[1]*A[1][1] + C2[2]*A[1][2] + C2[3]*A[1][3] + C2[4]*A[1][4] + y7[Sz7-1][10]<<endl;
		cout<<C2[0]*A[2][0] + C2[1]*A[2][1] + C2[2]*A[2][2] + C2[3]*A[2][3] + C2[4]*A[2][4] + y7[Sz7-1][12]<<endl;
		cout<<C2[0]*A[3][0] + C2[1]*A[3][1] + C2[2]*A[3][2] + C2[3]*A[3][3] + C2[4]*A[3][4] + y7[Sz7-1][1]<<endl;
		cout<<C2[0]*A[4][0] + C2[1]*A[4][1] + C2[2]*A[4][2] + C2[3]*A[4][3] + C2[4]*A[4][4] + y7[Sz7-1][2]<<endl;
		cout<<"Done!"<<endl;
		cout<<"C1_2 = "<<C2[0]<<endl;
		cout<<"C2_2 = "<<C2[1]<<endl;
		cout<<"C3_2 = "<<C2[2]<<endl;
		cout<<"C4_2 = "<<C2[3]<<endl;
		cout<<"C5_2 = "<<C2[4]<<endl;}

	double mobility=-C2[3]/C1[3];
	
	if(VERBOSE)
		{cout<<"\n\nElectrophoretic Mobility:\t"<<mobility<<" m^2/Vs\n";
		cout<<"ka:\t\t"<<a/Debye<<endl;
		cout<<"Zeta Potential:\t"<<ZP<<" V\n";
		cout<<"E:\t\t"<<3*S.viscosity*ec*mobility/(2*eo*S.dielectric*kb*S.temperature)<<endl;
		cout<<"y:\t\t"<<ec*ZP/(kb*S.temperature)<<endl;}

	return mobility;}

double four_ion_problem(double ZP,double a,Solution S,Error E,bool VERBOSE)
	{// Convert Debye Length from Angstroms to Meters
	double Debye=S.debyeLength*1e-10;
	// Final Position: Solvated Radius [meters]
	double rf=a;
	double r0;
	double dr=-1e-15;	// integrate from Inf to Protein Radius
	if(Debye>a){r0=5*Debye;}
	else{r0=5*a;}
		
	if(VERBOSE){cout<<"4 ION PROBLEM\nInitial Position:\t"<<r0<<" m"<<"\nFinal Position:\t\t"<<rf<<" m"<<"\nStep Size:\t\t"<<dr<<" m"<<endl;}

	// Calculate Poisson Boltzmann Coefficient and Determine Initial Position away from Protein
	double C=calc_PBE_coefficient(ZP,r0,rf,dr,S,E,VERBOSE);
	if(VERBOSE){cout<<"PBE Coeff: "<<C<<endl;}
	if(isnan(C)){return nan("");}
	int N=7+2*S.numIon+2+3*S.numIon+3;//32;
	state_type y(N);
	// F Function
	// y[0] = f
	// y[1] = df/dr
	// y[2] = d2f/dr2
	// y[3] = d3f/dr3
	// y[4] = d4f/dr4
	// Equilibrium Electric Pontential
	// y[5] = potential
	// y[6] = d(potential)/dr
	// Induced Electrochemical Potential Function (4 ion)
	// y[7] = IEP1
	// y[8] = d(IEP1)/dr
	// y[9] = IEP2
	// y[10] = d(IEP2)/dr
	// y[11] = IEP3
	// y[12] = d(IEP3)/dr
	// y[13] = IEP4
	// y[14] = d(IEP4)/dr
	// Constants
	// y[15] = Temp					// [Kelvin]
	// y[16] = numIon
	// Ion Properties
	// y[17] = Ion1.val;		
	// y[18] = Ion1.conc;		// mM [mol/m^3]
	// y[19] = Ion1.radius;		// meters
	// y[20] = Ion2.val;
	// y[21] = Ion2.conc;
	// y[22] = Ion2.radius;
	// y[23] = Ion3.val;
	// y[24] = Ion3.conc;
	// y[25] = Ion3.radius;
	// y[26] = Ion4.val;
	// y[27] = Ion4.conc;
	// y[28] = Ion4.radius;
	// Solution Properties
	// y[29] = relative dielectric	[dimensionless]
	// y[30] = Debye Length [meters]
	// y[31] = viscosity  [Pa s]

	// Collect Error Data
	double abs_err=E.absolute;
	double rel_err=E.relative;
	double a_x=E.a_x;
	double a_dxdt=E.a_dxdt;
	controlled_stepper_type controlled_stepper(default_error_checker<double,range_algebra,default_operations>(abs_err,rel_err,a_x,a_dxdt));	

	// Calculate Mobility
	vector<state_type> y1,y2,y3,y4,y5,y6,y7,y8;
	vector<double> r1,r2,r3,r4,r5,r6,r7,r8;
	// Specify Boundary Conditions / Initial Values Infinitely Far from the Particle
	if(VERBOSE){cout<<"Solving Homogeneous Problem for 1st Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=1/(r0*r0);y[8]=-2/(r0*r0*r0);
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Constants
	y[15]=S.temperature;y[16]=S.numIon;
	// Fist Ion Properties
	y[17]=S.ionVal[0];y[18]=S.ionConc[0]*1000;y[19]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[20]=S.ionVal[1];y[21]=S.ionConc[1]*1000;y[22]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[23]=S.ionVal[2];y[24]=S.ionConc[2]*1000;y[25]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[26]=S.ionVal[3];y[27]=S.ionConc[3]*1000;y[28]=S.ionRad[3]*1e-10;
	// Solution Properties
	y[29]=S.dielectric;y[30]=Debye;y[31]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps1=integrate_adaptive(controlled_stepper,homogeneous_form_4ion,y,r0,rf,dr,push_back_state_and_time(y1,r1));
	int Sz1=steps1+1;
	if(VERBOSE){cout<<"Done! ("<<Sz1<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 2nd Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=1/(r0*r0);y[10]=-2/(r0*r0*r0);
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Constants
	y[15]=S.temperature;y[16]=S.numIon;
	// Fist Ion Properties
	y[17]=S.ionVal[0];y[18]=S.ionConc[0]*1000;y[19]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[20]=S.ionVal[1];y[21]=S.ionConc[1]*1000;y[22]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[23]=S.ionVal[2];y[24]=S.ionConc[2]*1000;y[25]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[26]=S.ionVal[3];y[27]=S.ionConc[3]*1000;y[28]=S.ionRad[3]*1e-10;
	// Solution Properties
	y[29]=S.dielectric;y[30]=Debye;y[31]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps2=integrate_adaptive(controlled_stepper,homogeneous_form_4ion,y,r0,rf,dr,push_back_state_and_time(y2,r2));
	int Sz2=steps2+1;
	if(VERBOSE){cout<<"Done! ("<<Sz2<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 3rd Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=1/(r0*r0);y[12]=-2/(r0*r0*r0);
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Constants
	y[15]=S.temperature;y[16]=S.numIon;
	// Fist Ion Properties
	y[17]=S.ionVal[0];y[18]=S.ionConc[0]*1000;y[19]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[20]=S.ionVal[1];y[21]=S.ionConc[1]*1000;y[22]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[23]=S.ionVal[2];y[24]=S.ionConc[2]*1000;y[25]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[26]=S.ionVal[3];y[27]=S.ionConc[3]*1000;y[28]=S.ionRad[3]*1e-10;
	// Solution Properties
	y[29]=S.dielectric;y[30]=Debye;y[31]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps3=integrate_adaptive(controlled_stepper,homogeneous_form_4ion,y,r0,rf,dr,push_back_state_and_time(y3,r3));
	int Sz3=steps3+1;
	if(VERBOSE){cout<<"Done! ("<<Sz3<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 4th Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=1/(r0*r0);y[14]=-2/(r0*r0*r0);
	// Constants
	y[15]=S.temperature;y[16]=S.numIon;
	// Fist Ion Properties
	y[17]=S.ionVal[0];y[18]=S.ionConc[0]*1000;y[19]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[20]=S.ionVal[1];y[21]=S.ionConc[1]*1000;y[22]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[23]=S.ionVal[2];y[24]=S.ionConc[2]*1000;y[25]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[26]=S.ionVal[3];y[27]=S.ionConc[3]*1000;y[28]=S.ionRad[3]*1e-10;
	// Solution Properties
	y[29]=S.dielectric;y[30]=Debye;y[31]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps4=integrate_adaptive(controlled_stepper,homogeneous_form_4ion,y,r0,rf,dr,push_back_state_and_time(y4,r4));
	int Sz4=steps4+1;
	if(VERBOSE){cout<<"Done! ("<<Sz4<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 1st f term...";}
	// F Function
	y[0]=r0;y[1]=1;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Constants
	y[15]=S.temperature;y[16]=S.numIon;
	// Fist Ion Properties
	y[17]=S.ionVal[0];y[18]=S.ionConc[0]*1000;y[19]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[20]=S.ionVal[1];y[21]=S.ionConc[1]*1000;y[22]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[23]=S.ionVal[2];y[24]=S.ionConc[2]*1000;y[25]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[26]=S.ionVal[3];y[27]=S.ionConc[3]*1000;y[28]=S.ionRad[3]*1e-10;
	// Solution Properties
	y[29]=S.dielectric;y[30]=Debye;y[31]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps5=integrate_adaptive(controlled_stepper,homogeneous_form_4ion,y,r0,rf,dr,push_back_state_and_time(y5,r5));
	int Sz5=steps5+1;
	if(VERBOSE){cout<<"Done! ("<<Sz5<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 2nd f term...";}
	// F Function
	y[0]=1/r0;y[1]=-1/(r0*r0);y[2]=2/(r0*r0*r0);y[3]=-6/(r0*r0*r0*r0);y[4]=24/(r0*r0*r0*r0*r0);
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Constants
	y[15]=S.temperature;y[16]=S.numIon;
	// Fist Ion Properties
	y[17]=S.ionVal[0];y[18]=S.ionConc[0]*1000;y[19]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[20]=S.ionVal[1];y[21]=S.ionConc[1]*1000;y[22]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[23]=S.ionVal[2];y[24]=S.ionConc[2]*1000;y[25]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[26]=S.ionVal[3];y[27]=S.ionConc[3]*1000;y[28]=S.ionRad[3]*1e-10;
	// Solution Properties
	y[29]=S.dielectric;y[30]=Debye;y[31]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps6=integrate_adaptive(controlled_stepper,homogeneous_form_4ion,y,r0,rf,dr,push_back_state_and_time(y6,r6));
	int Sz6=steps6+1;
	if(VERBOSE){cout<<"Done! ("<<Sz6<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Flow Field Problem...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Constants
	y[15]=S.temperature;y[16]=S.numIon;
	// Fist Ion Properties
	y[17]=S.ionVal[0];y[18]=S.ionConc[0]*1000;y[19]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[20]=S.ionVal[1];y[21]=S.ionConc[1]*1000;y[22]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[23]=S.ionVal[2];y[24]=S.ionConc[2]*1000;y[25]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[26]=S.ionVal[3];y[27]=S.ionConc[3]*1000;y[28]=S.ionRad[3]*1e-10;
	// Solution Properties
	y[29]=S.dielectric;y[30]=Debye;y[31]=S.viscosity;
	// Solve Flow Field Perturbation 
	size_t steps7=integrate_adaptive(controlled_stepper,inhomogeneous_Prob1_4ion,y,r0,rf,dr,push_back_state_and_time(y7,r7));
	int Sz7=steps7+1;
	if(VERBOSE){cout<<"Done! ("<<Sz7<<" points)"<<endl;}

	if(VERBOSE){cout<<"Calculating Asymptotic Constants...";}
	alglib::real_2d_array A;
	alglib::real_1d_array C1,B1;
	// Define System Matrix containing coefficients of linear system
	A.setlength(6,6);
	for(int i=0;i<6;i++)
		{for(int j=0;j<6;j++)
			{A[i][j]=0;}}
	// d(IEP1)/dr
	A[0][0]=y1[Sz1-1][8];A[0][1]=y2[Sz2-1][8];A[0][2]=y3[Sz3-1][8];A[0][3]=y4[Sz4-1][8];A[0][4]=y5[Sz5-1][8];A[0][5]=y6[Sz6-1][8];
	// d(IEP2)/dr
	A[1][0]=y1[Sz1-1][10];A[1][1]=y2[Sz2-1][10];A[1][2]=y3[Sz3-1][10];A[1][3]=y4[Sz4-1][10];A[1][4]=y5[Sz5-1][10];A[1][5]=y6[Sz6-1][10];
	// d(IEP3)/dr
	A[2][0]=y1[Sz1-1][12];A[2][1]=y2[Sz2-1][12];A[2][2]=y3[Sz3-1][12];A[2][3]=y4[Sz4-1][12];A[2][4]=y5[Sz5-1][12];A[2][5]=y6[Sz6-1][12];
	// d(IEP4)/dr
	A[3][0]=y1[Sz1-1][14];A[3][1]=y2[Sz2-1][14];A[3][2]=y3[Sz3-1][14];A[3][3]=y4[Sz4-1][14];A[3][4]=y5[Sz5-1][14];A[3][5]=y6[Sz6-1][14];
	// f 1st (df/dr)
	A[4][0]=y1[Sz1-1][1];A[4][1]=y2[Sz2-1][1];A[4][2]=y3[Sz3-1][1];A[4][3]=y4[Sz4-1][1];A[4][4]=y5[Sz5-1][1];A[4][5]=y6[Sz6-1][1];
	// f 2nd (d2f/dr2)
	A[5][0]=y1[Sz1-1][2];A[5][1]=y2[Sz2-1][2];A[5][2]=y3[Sz3-1][2];A[5][3]=y4[Sz4-1][2];A[5][4]=y5[Sz5-1][2];A[5][5]=y6[Sz6-1][2];
	// Define Solution Matrix for Problem #1
	B1.setlength(6);
	// d(IEP1)/dr = 0
	B1[0]=0-y7[Sz7-1][8];
	// d(IEP2)/dr = 0
	B1[1]=0-y7[Sz7-1][10];
	// d(IEP3)/dr = 0
	B1[2]=0-y7[Sz7-1][12];
	// d(IEP4)/dr = 0
	B1[3]=0-y7[Sz7-1][14];
	// df/dr = -Radius/2
	B1[4]=-a/2 - y7[Sz7-1][1];
	// d2f/dr2 = -0.5
	B1[5]=-0.5 - y7[Sz7-1][2];
	// Create Output Parameters
	alglib::ae_int_t info;
	alglib::densesolverreport rep;
	alglib::rmatrixsolve(A,6,B1,info,rep,C1);
	if(VERBOSE)
		{cout<<"\n";
		cout<<C1[0]*A[0][0] + C1[1]*A[0][1] + C1[2]*A[0][2] + C1[3]*A[0][3] + C1[4]*A[0][4] + C1[5]*A[0][5] + y7[Sz7-1][8]<<endl;
		cout<<C1[0]*A[1][0] + C1[1]*A[1][1] + C1[2]*A[1][2] + C1[3]*A[1][3] + C1[4]*A[1][4] + C1[5]*A[1][5] + y7[Sz7-1][10]<<endl;
		cout<<C1[0]*A[2][0] + C1[1]*A[2][1] + C1[2]*A[2][2] + C1[3]*A[2][3] + C1[4]*A[2][4] + C1[5]*A[2][5] + y7[Sz7-1][12]<<endl;
		cout<<C1[0]*A[3][0] + C1[1]*A[3][1] + C1[2]*A[3][2] + C1[3]*A[3][3] + C1[4]*A[3][4] + C1[5]*A[3][5] + y7[Sz7-1][14]<<endl;
		cout<<C1[0]*A[4][0] + C1[1]*A[4][1] + C1[2]*A[4][2] + C1[3]*A[4][3] + C1[4]*A[4][4] + C1[5]*A[4][5] + y7[Sz7-1][1]<<endl;
		cout<<C1[0]*A[5][0] + C1[1]*A[5][1] + C1[2]*A[5][2] + C1[3]*A[5][3] + C1[4]*A[5][4] + C1[5]*A[5][5] + y7[Sz7-1][2]<<endl;
		cout<<"Done!"<<endl;
		cout<<"C1_1 = "<<C1[0]<<endl;
		cout<<"C2_1 = "<<C1[1]<<endl;
		cout<<"C3_1 = "<<C1[2]<<endl;
		cout<<"C4_1 = "<<C1[3]<<endl;
		cout<<"C5_1 = "<<C1[4]<<endl;
		cout<<"C6_1 = "<<C1[5]<<endl;}

	if(VERBOSE){cout<<"Solving Electric Field Problem...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Constants
	y[15]=S.temperature;y[16]=S.numIon;
	// Fist Ion Properties
	y[17]=S.ionVal[0];y[18]=S.ionConc[0]*1000;y[19]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[20]=S.ionVal[1];y[21]=S.ionConc[1]*1000;y[22]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[23]=S.ionVal[2];y[24]=S.ionConc[2]*1000;y[25]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[26]=S.ionVal[3];y[27]=S.ionConc[3]*1000;y[28]=S.ionRad[3]*1e-10;
	// Solution Properties
	y[29]=S.dielectric;y[30]=Debye;y[31]=S.viscosity;
	// Solve Electric Field Perturbation
	size_t steps8=integrate_adaptive(controlled_stepper,inhomogeneous_Prob2_4ion,y,r0,rf,dr,push_back_state_and_time(y8,r8));
	int Sz8=steps8+1;
	if(VERBOSE){cout<<"Done! ("<<Sz8<<" points)"<<endl;}

	if(VERBOSE){cout<<"Calculating Asymptotic Constants...";}
	alglib::real_1d_array C2,B2;
	// d(IEP1)/dr
	A[0][0]=y1[Sz1-1][8];A[0][1]=y2[Sz2-1][8];A[0][2]=y3[Sz3-1][8];A[0][3]=y4[Sz4-1][8];A[0][4]=y5[Sz5-1][8];A[0][5]=y6[Sz6-1][8];
	// d(IEP2)/dr
	A[1][0]=y1[Sz1-1][10];A[1][1]=y2[Sz2-1][10];A[1][2]=y3[Sz3-1][10];A[1][3]=y4[Sz4-1][10];A[1][4]=y5[Sz5-1][10];A[1][5]=y6[Sz6-1][10];
	// d(IEP3)/dr
	A[2][0]=y1[Sz1-1][12];A[2][1]=y2[Sz2-1][12];A[2][2]=y3[Sz3-1][12];A[2][3]=y4[Sz4-1][12];A[2][4]=y5[Sz5-1][12];A[2][5]=y6[Sz6-1][12];
	// d(IEP4)/dr
	A[3][0]=y1[Sz1-1][14];A[3][1]=y2[Sz2-1][14];A[3][2]=y3[Sz3-1][14];A[3][3]=y4[Sz4-1][14];A[3][4]=y5[Sz5-1][14];A[3][5]=y6[Sz6-1][14];
	// f 1st (df/dr)
	A[4][0]=y1[Sz1-1][1];A[4][1]=y2[Sz2-1][1];A[4][2]=y3[Sz3-1][1];A[4][3]=y4[Sz4-1][1];A[4][4]=y5[Sz5-1][1];A[4][5]=y6[Sz6-1][1];
	// f 2nd (d2f/dr2)
	A[5][0]=y1[Sz1-1][2];A[5][1]=y2[Sz2-1][2];A[5][2]=y3[Sz3-1][2];A[5][3]=y4[Sz4-1][2];A[5][4]=y5[Sz5-1][2];A[5][5]=y6[Sz6-1][2];
	// Define Solution Matrix for Problem #2
	B2.setlength(6);
	// d(IEP1)/dr = -1
	B2[0]=-1-y8[Sz8-1][8];
	// d(IEP2)/dr = -1
	B2[1]=-1-y8[Sz8-1][10];
	// d(IEP3)/dr = -1
	B2[2]=-1-y8[Sz8-1][12];
	// d(IEP4)/dr = -1
	B2[3]=-1-y8[Sz8-1][14];
	// df/dr = 0
	B2[4]=0-y8[Sz8-1][1];
	// d2f/dr2 = 0
	B2[5]=0-y8[Sz8-1][2];
	// Create Output Parameters
	alglib::rmatrixsolve(A,6,B2,info,rep,C2);
	if(VERBOSE)
		{cout<<"\n";
		cout<<C2[0]*A[0][0] + C2[1]*A[0][1] + C2[2]*A[0][2] + C2[3]*A[0][3] + C2[4]*A[0][4] + C2[5]*A[0][5] + y8[Sz8-1][8]<<endl;
		cout<<C2[0]*A[1][0] + C2[1]*A[1][1] + C2[2]*A[1][2] + C2[3]*A[1][3] + C2[4]*A[1][4] + C2[5]*A[1][5] + y8[Sz8-1][10]<<endl;
		cout<<C2[0]*A[2][0] + C2[1]*A[2][1] + C2[2]*A[2][2] + C2[3]*A[2][3] + C2[4]*A[2][4] + C2[5]*A[2][5] + y8[Sz8-1][12]<<endl;
		cout<<C2[0]*A[3][0] + C2[1]*A[3][1] + C2[2]*A[3][2] + C2[3]*A[3][3] + C2[4]*A[3][4] + C2[5]*A[3][5] + y8[Sz8-1][14]<<endl;
		cout<<C2[0]*A[4][0] + C2[1]*A[4][1] + C2[2]*A[4][2] + C2[3]*A[4][3] + C2[4]*A[4][4] + C2[5]*A[4][5] + y8[Sz8-1][1]<<endl;
		cout<<C2[0]*A[5][0] + C2[1]*A[5][1] + C2[2]*A[5][2] + C2[3]*A[5][3] + C2[4]*A[5][4] + C2[5]*A[5][5] + y8[Sz8-1][2]<<endl;
		cout<<"Done!"<<endl;
		cout<<"C1_2 = "<<C2[0]<<endl;
		cout<<"C2_2 = "<<C2[1]<<endl;
		cout<<"C3_2 = "<<C2[2]<<endl;
		cout<<"C4_2 = "<<C2[3]<<endl;
		cout<<"C5_2 = "<<C2[4]<<endl;
		cout<<"C6_2 = "<<C2[5]<<endl;}

	double mobility=-C2[4]/C1[4];
	
	if(VERBOSE)
		{cout<<"\n\nElectrophoretic Mobility:\t"<<mobility<<" m^2/Vs\n";
		cout<<"ka:\t\t"<<a/Debye<<endl;
		cout<<"Zeta Potential:\t"<<ZP<<" V\n";
		cout<<"E:\t\t"<<3*S.viscosity*ec*mobility/(2*eo*S.dielectric*kb*S.temperature)<<endl;
		cout<<"y:\t\t"<<ec*ZP/(kb*S.temperature)<<endl;}

	return mobility;}

double five_ion_problem(double ZP,double a,Solution S,Error E,bool VERBOSE)
	{// Convert Debye Length from Angstroms to Meters
	double Debye=S.debyeLength*1e-10;
	// Final Position: Solvated Radius [meters]
	double rf=a;
	double r0;
	double dr=-1e-15;	// integrate from Inf to Protein Radius
	if(Debye>a){r0=5*Debye;}
	else{r0=5*a;}
		
	if(VERBOSE){cout<<"5 ION PROBLEM\nInitial Position:\t"<<r0<<" m"<<"\nFinal Position:\t\t"<<rf<<" m"<<"\nStep Size:\t\t"<<dr<<" m"<<endl;}

	// Calculate Poisson Boltzmann Coefficient and Determine Initial Position away from Protein
	double C=calc_PBE_coefficient(ZP,r0,rf,dr,S,E,VERBOSE);
	if(VERBOSE){cout<<"PBE Coeff: "<<C<<endl;}
	if(isnan(C)){return nan("");}
	int N=7+2*S.numIon+2+3*S.numIon+3;//37;
	state_type y(N);
	// F Function
	// y[0] = f
	// y[1] = df/dr
	// y[2] = d2f/dr2
	// y[3] = d3f/dr3
	// y[4] = d4f/dr4
	// Equilibrium Electric Pontential
	// y[5] = potential
	// y[6] = d(potential)/dr
	// Induced Electrochemical Potential Function (5 ion)
	// y[7] = IEP1
	// y[8] = d(IEP1)/dr
	// y[9] = IEP2
	// y[10] = d(IEP2)/dr
	// y[11] = IEP3
	// y[12] = d(IEP3)/dr
	// y[13] = IEP4
	// y[14] = d(IEP4)/dr
	// y[15] = IEP5
	// y[16] = d(IEP5)/dr
	// Constants
	// y[17] = Temp					// [Kelvin]
	// y[18] = numIon
	// Ion Properties
	// y[19] = Ion1.val;		
	// y[20] = Ion1.conc;		// mM [mol/m^3]
	// y[21] = Ion1.radius;		// meters
	// y[22] = Ion2.val;
	// y[23] = Ion2.conc;
	// y[24] = Ion2.radius;
	// y[25] = Ion3.val;
	// y[26] = Ion3.conc;
	// y[27] = Ion3.radius;
	// y[28] = Ion4.val;
	// y[29] = Ion4.conc;
	// y[30] = Ion4.radius;
	// y[31] = Ion5.val;
	// y[32] = Ion5.conc;
	// y[33] = Ion5.radius;
	// Solution Properties
	// y[34] = relative dielectric	[dimensionless]
	// y[35] = Debye Length [meters]
	// y[36] = viscosity  [Pa s]
	// Collect Error Data
	double abs_err=E.absolute;
	double rel_err=E.relative;
	double a_x=E.a_x;
	double a_dxdt=E.a_dxdt;
	controlled_stepper_type controlled_stepper(default_error_checker<double,range_algebra,default_operations>(abs_err,rel_err,a_x,a_dxdt));	

	// Calculate Mobility
	vector<state_type> y1,y2,y3,y4,y5,y6,y7,y8,y9;
	vector<double> r1,r2,r3,r4,r5,r6,r7,r8,r9;
	// Specify Boundary Conditions / Initial Values Infinitely Far from the Particle
	if(VERBOSE){cout<<"Solving Homogeneous Problem for 1st Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=1/(r0*r0);y[8]=-2/(r0*r0*r0);
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fiveth Ion
	y[15]=0;y[16]=0;
	// Constants
	y[17]=S.temperature;y[18]=S.numIon;
	// Fist Ion Properties
	y[19]=S.ionVal[0];y[20]=S.ionConc[0]*1000;y[21]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[22]=S.ionVal[1];y[23]=S.ionConc[1]*1000;y[24]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[25]=S.ionVal[2];y[26]=S.ionConc[2]*1000;y[27]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[28]=S.ionVal[3];y[29]=S.ionConc[3]*1000;y[30]=S.ionRad[3]*1e-10;
	// Fiveth Ion Properties
	y[31]=S.ionVal[4];y[32]=S.ionConc[4]*1000;y[33]=S.ionRad[4]*1e-10;
	// Solution Properties
	y[34]=S.dielectric;y[35]=Debye;y[36]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps1=integrate_adaptive(controlled_stepper,homogeneous_form_5ion,y,r0,rf,dr,push_back_state_and_time(y1,r1));
	int Sz1=steps1+1;
	if(VERBOSE){cout<<"Done! ("<<Sz1<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 2nd Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=1/(r0*r0);y[10]=-2/(r0*r0*r0);
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fiveth Ion
	y[15]=0;y[16]=0;
	// Constants
	y[17]=S.temperature;y[18]=S.numIon;
	// Fist Ion Properties
	y[19]=S.ionVal[0];y[20]=S.ionConc[0]*1000;y[21]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[22]=S.ionVal[1];y[23]=S.ionConc[1]*1000;y[24]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[25]=S.ionVal[2];y[26]=S.ionConc[2]*1000;y[27]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[28]=S.ionVal[3];y[29]=S.ionConc[3]*1000;y[30]=S.ionRad[3]*1e-10;
	// Fiveth Ion Properties
	y[31]=S.ionVal[4];y[32]=S.ionConc[4]*1000;y[33]=S.ionRad[4]*1e-10;
	// Solution Properties
	y[34]=S.dielectric;y[35]=Debye;y[36]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps2=integrate_adaptive(controlled_stepper,homogeneous_form_5ion,y,r0,rf,dr,push_back_state_and_time(y2,r2));
	int Sz2=steps2+1;
	if(VERBOSE){cout<<"Done! ("<<Sz2<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 3rd Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=1/(r0*r0);y[12]=-2/(r0*r0*r0);
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fiveth Ion
	y[15]=0;y[16]=0;
	// Constants
	y[17]=S.temperature;y[18]=S.numIon;
	// Fist Ion Properties
	y[19]=S.ionVal[0];y[20]=S.ionConc[0]*1000;y[21]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[22]=S.ionVal[1];y[23]=S.ionConc[1]*1000;y[24]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[25]=S.ionVal[2];y[26]=S.ionConc[2]*1000;y[27]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[28]=S.ionVal[3];y[29]=S.ionConc[3]*1000;y[30]=S.ionRad[3]*1e-10;
	// Fiveth Ion Properties
	y[31]=S.ionVal[4];y[32]=S.ionConc[4]*1000;y[33]=S.ionRad[4]*1e-10;
	// Solution Properties
	y[34]=S.dielectric;y[35]=Debye;y[36]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps3=integrate_adaptive(controlled_stepper,homogeneous_form_5ion,y,r0,rf,dr,push_back_state_and_time(y3,r3));
	int Sz3=steps3+1;
	if(VERBOSE){cout<<"Done! ("<<Sz3<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 4th Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=1/(r0*r0);y[14]=-2/(r0*r0*r0);
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Constants
	y[17]=S.temperature;y[18]=S.numIon;
	// Fist Ion Properties
	y[19]=S.ionVal[0];y[20]=S.ionConc[0]*1000;y[21]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[22]=S.ionVal[1];y[23]=S.ionConc[1]*1000;y[24]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[25]=S.ionVal[2];y[26]=S.ionConc[2]*1000;y[27]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[28]=S.ionVal[3];y[29]=S.ionConc[3]*1000;y[30]=S.ionRad[3]*1e-10;
	// Fiveth Ion Properties
	y[31]=S.ionVal[4];y[32]=S.ionConc[4]*1000;y[33]=S.ionRad[4]*1e-10;
	// Solution Properties
	y[34]=S.dielectric;y[35]=Debye;y[36]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps4=integrate_adaptive(controlled_stepper,homogeneous_form_5ion,y,r0,rf,dr,push_back_state_and_time(y4,r4));
	int Sz4=steps4+1;
	if(VERBOSE){cout<<"Done! ("<<Sz4<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 5th Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=1/(r0*r0);y[16]=-2/(r0*r0*r0);
	// Constants
	y[17]=S.temperature;y[18]=S.numIon;
	// Fist Ion Properties
	y[19]=S.ionVal[0];y[20]=S.ionConc[0]*1000;y[21]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[22]=S.ionVal[1];y[23]=S.ionConc[1]*1000;y[24]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[25]=S.ionVal[2];y[26]=S.ionConc[2]*1000;y[27]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[28]=S.ionVal[3];y[29]=S.ionConc[3]*1000;y[30]=S.ionRad[3]*1e-10;
	// Fiveth Ion Properties
	y[31]=S.ionVal[4];y[32]=S.ionConc[4]*1000;y[33]=S.ionRad[4]*1e-10;
	// Solution Properties
	y[34]=S.dielectric;y[35]=Debye;y[36]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps5=integrate_adaptive(controlled_stepper,homogeneous_form_5ion,y,r0,rf,dr,push_back_state_and_time(y5,r5));
	int Sz5=steps5+1;
	if(VERBOSE){cout<<"Done! ("<<Sz5<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 1st f term...";}
	// F Function
	y[0]=r0;y[1]=1;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Constants
	y[17]=S.temperature;y[18]=S.numIon;
	// Fist Ion Properties
	y[19]=S.ionVal[0];y[20]=S.ionConc[0]*1000;y[21]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[22]=S.ionVal[1];y[23]=S.ionConc[1]*1000;y[24]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[25]=S.ionVal[2];y[26]=S.ionConc[2]*1000;y[27]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[28]=S.ionVal[3];y[29]=S.ionConc[3]*1000;y[30]=S.ionRad[3]*1e-10;
	// Fiveth Ion Properties
	y[31]=S.ionVal[4];y[32]=S.ionConc[4]*1000;y[33]=S.ionRad[4]*1e-10;
	// Solution Properties
	y[34]=S.dielectric;y[35]=Debye;y[36]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps6=integrate_adaptive(controlled_stepper,homogeneous_form_5ion,y,r0,rf,dr,push_back_state_and_time(y6,r6));
	int Sz6=steps6+1;
	if(VERBOSE){cout<<"Done! ("<<Sz6<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 2nd f term...";}
	// F Function
	y[0]=1/r0;y[1]=-1/(r0*r0);y[2]=2/(r0*r0*r0);y[3]=-6/(r0*r0*r0*r0);y[4]=24/(r0*r0*r0*r0*r0);
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Constants
	y[17]=S.temperature;y[18]=S.numIon;
	// Fist Ion Properties
	y[19]=S.ionVal[0];y[20]=S.ionConc[0]*1000;y[21]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[22]=S.ionVal[1];y[23]=S.ionConc[1]*1000;y[24]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[25]=S.ionVal[2];y[26]=S.ionConc[2]*1000;y[27]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[28]=S.ionVal[3];y[29]=S.ionConc[3]*1000;y[30]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[31]=S.ionVal[4];y[32]=S.ionConc[4]*1000;y[33]=S.ionRad[4]*1e-10;
	// Solution Properties
	y[34]=S.dielectric;y[35]=Debye;y[36]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps7=integrate_adaptive(controlled_stepper,homogeneous_form_5ion,y,r0,rf,dr,push_back_state_and_time(y7,r7));
	int Sz7=steps7+1;
	if(VERBOSE){cout<<"Done! ("<<Sz7<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Flow Field Problem...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Constants
	y[17]=S.temperature;y[18]=S.numIon;
	// Fist Ion Properties
	y[19]=S.ionVal[0];y[20]=S.ionConc[0]*1000;y[21]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[22]=S.ionVal[1];y[23]=S.ionConc[1]*1000;y[24]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[25]=S.ionVal[2];y[26]=S.ionConc[2]*1000;y[27]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[28]=S.ionVal[3];y[29]=S.ionConc[3]*1000;y[30]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[31]=S.ionVal[4];y[32]=S.ionConc[4]*1000;y[33]=S.ionRad[4]*1e-10;
	// Solution Properties
	y[34]=S.dielectric;y[35]=Debye;y[36]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps8=integrate_adaptive(controlled_stepper,inhomogeneous_Prob1_5ion,y,r0,rf,dr,push_back_state_and_time(y8,r8));
	int Sz8=steps8+1;
	if(VERBOSE){cout<<"Done! ("<<Sz8<<" points)"<<endl;}

	if(VERBOSE){cout<<"Calculating Asymptotic Constants...";}
	alglib::real_2d_array A;
	alglib::real_1d_array C1,B1;
	// Define System Matrix containing coefficients of linear system
	A.setlength(7,7);
	for(int i=0;i<7;i++)
		{for(int j=0;j<7;j++)
			{A[i][j]=0;}}
	// d(IEP1)/dr
	A[0][0]=y1[Sz1-1][8];A[0][1]=y2[Sz2-1][8];A[0][2]=y3[Sz3-1][8];A[0][3]=y4[Sz4-1][8];A[0][4]=y5[Sz5-1][8];A[0][5]=y6[Sz6-1][8];A[0][6]=y7[Sz7-1][8];
	// d(IEP2)/dr
	A[1][0]=y1[Sz1-1][10];A[1][1]=y2[Sz2-1][10];A[1][2]=y3[Sz3-1][10];A[1][3]=y4[Sz4-1][10];A[1][4]=y5[Sz5-1][10];A[1][5]=y6[Sz6-1][10];A[1][6]=y7[Sz7-1][10];
	// d(IEP3)/dr
	A[2][0]=y1[Sz1-1][12];A[2][1]=y2[Sz2-1][12];A[2][2]=y3[Sz3-1][12];A[2][3]=y4[Sz4-1][12];A[2][4]=y5[Sz5-1][12];A[2][5]=y6[Sz6-1][12];A[2][6]=y7[Sz7-1][12];
	// d(IEP4)/dr
	A[3][0]=y1[Sz1-1][14];A[3][1]=y2[Sz2-1][14];A[3][2]=y3[Sz3-1][14];A[3][3]=y4[Sz4-1][14];A[3][4]=y5[Sz5-1][14];A[3][5]=y6[Sz6-1][14];A[3][6]=y7[Sz7-1][14];
	// d(IEP5)/dr
	A[4][0]=y1[Sz1-1][16];A[4][1]=y2[Sz2-1][16];A[4][2]=y3[Sz3-1][16];A[4][3]=y4[Sz4-1][16];A[4][4]=y5[Sz5-1][16];A[4][5]=y6[Sz6-1][16];A[4][6]=y7[Sz7-1][16];
	// f 1st (df/dr)
	A[5][0]=y1[Sz1-1][1];A[5][1]=y2[Sz2-1][1];A[5][2]=y3[Sz3-1][1];A[5][3]=y4[Sz4-1][1];A[5][4]=y5[Sz5-1][1];A[5][5]=y6[Sz6-1][1];A[5][6]=y7[Sz7-1][1];
	// f 2nd (d2f/dr2)
	A[6][0]=y1[Sz1-1][2];A[6][1]=y2[Sz2-1][2];A[6][2]=y3[Sz3-1][2];A[6][3]=y4[Sz4-1][2];A[6][4]=y5[Sz5-1][2];A[6][5]=y6[Sz6-1][2];A[6][6]=y7[Sz7-1][2];
	// Define Solution Matrix for Problem #1
	B1.setlength(7);
	// d(IEP1)/dr = 0
	B1[0]=0-y8[Sz8-1][8];
	// d(IEP2)/dr = 0
	B1[1]=0-y8[Sz8-1][10];
	// d(IEP3)/dr = 0
	B1[2]=0-y8[Sz8-1][12];
	// d(IEP4)/dr = 0
	B1[3]=0-y8[Sz8-1][14];
	// d(IEP5)/dr = 0
	B1[4]=0-y8[Sz8-1][16];
	// df/dr = -Radius/2
	B1[5]=-a/2 - y8[Sz8-1][1];
	// d2f/dr2 = -0.5
	B1[6]=-0.5 - y8[Sz8-1][2];
	// Create Output Parameters
	alglib::ae_int_t info;
	alglib::densesolverreport rep;
	alglib::rmatrixsolve(A,7,B1,info,rep,C1);
	if(VERBOSE)
		{cout<<"\n";
		cout<<C1[0]*A[0][0] + C1[1]*A[0][1] + C1[2]*A[0][2] + C1[3]*A[0][3] + C1[4]*A[0][4] + C1[5]*A[0][5] + C1[6]*A[0][6] + y8[Sz8-1][8]<<endl;
		cout<<C1[0]*A[1][0] + C1[1]*A[1][1] + C1[2]*A[1][2] + C1[3]*A[1][3] + C1[4]*A[1][4] + C1[5]*A[1][5] + C1[6]*A[1][6] + y8[Sz8-1][10]<<endl;
		cout<<C1[0]*A[2][0] + C1[1]*A[2][1] + C1[2]*A[2][2] + C1[3]*A[2][3] + C1[4]*A[2][4] + C1[5]*A[2][5] + C1[6]*A[2][6] + y8[Sz8-1][12]<<endl;
		cout<<C1[0]*A[3][0] + C1[1]*A[3][1] + C1[2]*A[3][2] + C1[3]*A[3][3] + C1[4]*A[3][4] + C1[5]*A[3][5] + C1[6]*A[3][6] + y8[Sz8-1][14]<<endl;
		cout<<C1[0]*A[4][0] + C1[1]*A[4][1] + C1[2]*A[4][2] + C1[3]*A[4][3] + C1[4]*A[4][4] + C1[5]*A[4][5] + C1[6]*A[4][6] + y8[Sz8-1][16]<<endl;
		cout<<C1[0]*A[5][0] + C1[1]*A[5][1] + C1[2]*A[5][2] + C1[3]*A[5][3] + C1[4]*A[5][4] + C1[5]*A[5][5] + C1[6]*A[5][6] + y8[Sz8-1][1]<<endl;
		cout<<C1[0]*A[6][0] + C1[1]*A[6][1] + C1[2]*A[6][2] + C1[3]*A[6][3] + C1[4]*A[6][4] + C1[5]*A[6][5] + C1[6]*A[6][6] + y8[Sz8-1][2]<<endl;
		cout<<"Done!"<<endl;
		cout<<"C1_1 = "<<C1[0]<<endl;
		cout<<"C2_1 = "<<C1[1]<<endl;
		cout<<"C3_1 = "<<C1[2]<<endl;
		cout<<"C4_1 = "<<C1[3]<<endl;
		cout<<"C5_1 = "<<C1[4]<<endl;
		cout<<"C6_1 = "<<C1[5]<<endl;
		cout<<"C7_1 = "<<C1[6]<<endl;}

	if(VERBOSE){cout<<"Solving Electric Field Problem...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Constants
	y[17]=S.temperature;y[18]=S.numIon;
	// Fist Ion Properties
	y[19]=S.ionVal[0];y[20]=S.ionConc[0]*1000;y[21]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[22]=S.ionVal[1];y[23]=S.ionConc[1]*1000;y[24]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[25]=S.ionVal[2];y[26]=S.ionConc[2]*1000;y[27]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[28]=S.ionVal[3];y[29]=S.ionConc[3]*1000;y[30]=S.ionRad[3]*1e-10;
	// Fiveth Ion Properties
	y[31]=S.ionVal[4];y[32]=S.ionConc[4]*1000;y[33]=S.ionRad[4]*1e-10;
	// Solution Properties
	y[34]=S.dielectric;y[35]=Debye;y[36]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps9=integrate_adaptive(controlled_stepper,inhomogeneous_Prob2_5ion,y,r0,rf,dr,push_back_state_and_time(y9,r9));
	int Sz9=steps9+1;
	if(VERBOSE){cout<<"Done! ("<<Sz9<<" points)"<<endl;}

	if(VERBOSE){cout<<"Calculating Asymptotic Constants...";}
	alglib::real_1d_array C2,B2;
	// d(IEP1)/dr
	A[0][0]=y1[Sz1-1][8];A[0][1]=y2[Sz2-1][8];A[0][2]=y3[Sz3-1][8];A[0][3]=y4[Sz4-1][8];A[0][4]=y5[Sz5-1][8];A[0][5]=y6[Sz6-1][8];A[0][6]=y7[Sz7-1][8];
	// d(IEP2)/dr
	A[1][0]=y1[Sz1-1][10];A[1][1]=y2[Sz2-1][10];A[1][2]=y3[Sz3-1][10];A[1][3]=y4[Sz4-1][10];A[1][4]=y5[Sz5-1][10];A[1][5]=y6[Sz6-1][10];A[1][6]=y7[Sz7-1][10];
	// d(IEP3)/dr
	A[2][0]=y1[Sz1-1][12];A[2][1]=y2[Sz2-1][12];A[2][2]=y3[Sz3-1][12];A[2][3]=y4[Sz4-1][12];A[2][4]=y5[Sz5-1][12];A[2][5]=y6[Sz6-1][12];A[2][6]=y7[Sz7-1][12];
	// d(IEP4)/dr
	A[3][0]=y1[Sz1-1][14];A[3][1]=y2[Sz2-1][14];A[3][2]=y3[Sz3-1][14];A[3][3]=y4[Sz4-1][14];A[3][4]=y5[Sz5-1][14];A[3][5]=y6[Sz6-1][14];A[3][6]=y7[Sz7-1][14];
	// d(IEP5)/dr
	A[4][0]=y1[Sz1-1][16];A[4][1]=y2[Sz2-1][16];A[4][2]=y3[Sz3-1][16];A[4][3]=y4[Sz4-1][16];A[4][4]=y5[Sz5-1][16];A[4][5]=y6[Sz6-1][16];A[4][6]=y7[Sz7-1][16];
	// f 1st (df/dr)
	A[5][0]=y1[Sz1-1][1];A[5][1]=y2[Sz2-1][1];A[5][2]=y3[Sz3-1][1];A[5][3]=y4[Sz4-1][1];A[5][4]=y5[Sz5-1][1];A[5][5]=y6[Sz6-1][1];A[5][6]=y7[Sz7-1][1];
	// f 2nd (d2f/dr2)
	A[6][0]=y1[Sz1-1][2];A[6][1]=y2[Sz2-1][2];A[6][2]=y3[Sz3-1][2];A[6][3]=y4[Sz4-1][2];A[6][4]=y5[Sz5-1][2];A[6][5]=y6[Sz6-1][2];A[6][6]=y7[Sz7-1][2];
	// Define Solution Matrix for Problem #2
	B2.setlength(7);
	// d(IEP1)/dr = -1
	B2[0]=-1-y9[Sz9-1][8];
	// d(IEP2)/dr = -1
	B2[1]=-1-y9[Sz9-1][10];
	// d(IEP3)/dr = -1
	B2[2]=-1-y9[Sz9-1][12];
	// d(IEP4)/dr = -1
	B2[3]=-1-y9[Sz9-1][14];
	// d(IEP5)/dr = -1
	B2[4]=-1-y9[Sz9-1][16];
	// df/dr = 0
	B2[5]=0-y9[Sz9-1][1];
	// d2f/dr2 = 0
	B2[6]=0-y9[Sz9-1][2];
	// Create Output Parameters
	alglib::rmatrixsolve(A,7,B2,info,rep,C2);
	if(VERBOSE)
		{cout<<"\n";
		cout<<C2[0]*A[0][0] + C2[1]*A[0][1] + C2[2]*A[0][2] + C2[3]*A[0][3] + C2[4]*A[0][4] + C2[5]*A[0][5] + C2[6]*A[0][6] + y9[Sz9-1][8]<<endl;
		cout<<C2[0]*A[1][0] + C2[1]*A[1][1] + C2[2]*A[1][2] + C2[3]*A[1][3] + C2[4]*A[1][4] + C2[5]*A[1][5] + C2[6]*A[1][6] + y9[Sz9-1][10]<<endl;
		cout<<C2[0]*A[2][0] + C2[1]*A[2][1] + C2[2]*A[2][2] + C2[3]*A[2][3] + C2[4]*A[2][4] + C2[5]*A[2][5] + C2[6]*A[2][6] + y9[Sz9-1][12]<<endl;
		cout<<C2[0]*A[3][0] + C2[1]*A[3][1] + C2[2]*A[3][2] + C2[3]*A[3][3] + C2[4]*A[3][4] + C2[5]*A[3][5] + C2[6]*A[3][6] + y9[Sz9-1][14]<<endl;
		cout<<C2[0]*A[4][0] + C2[1]*A[4][1] + C2[2]*A[4][2] + C2[3]*A[4][3] + C2[4]*A[4][4] + C2[5]*A[4][5] + C2[6]*A[4][6] + y9[Sz9-1][16]<<endl;
		cout<<C2[0]*A[5][0] + C2[1]*A[5][1] + C2[2]*A[5][2] + C2[3]*A[5][3] + C2[4]*A[5][4] + C2[5]*A[5][5] + C2[6]*A[5][6] + y9[Sz9-1][1]<<endl;
		cout<<C2[0]*A[6][0] + C2[1]*A[6][1] + C2[2]*A[6][2] + C2[3]*A[6][3] + C2[4]*A[6][4] + C2[5]*A[6][5] + C2[6]*A[6][6] + y9[Sz9-1][2]<<endl;
		cout<<"Done!"<<endl;
		cout<<"C1_2 = "<<C2[0]<<endl;
		cout<<"C2_2 = "<<C2[1]<<endl;
		cout<<"C3_2 = "<<C2[2]<<endl;
		cout<<"C4_2 = "<<C2[3]<<endl;
		cout<<"C5_2 = "<<C2[4]<<endl;
		cout<<"C6_2 = "<<C2[5]<<endl;
		cout<<"C7_2 = "<<C2[6]<<endl;}

	double mobility=-C2[5]/C1[5];
	
	if(VERBOSE)
		{cout<<"\n\nElectrophoretic Mobility:\t"<<mobility<<" m^2/Vs\n";
		cout<<"ka:\t\t"<<a/Debye<<endl;
		cout<<"Zeta Potential:\t"<<ZP<<" V\n";
		cout<<"E:\t\t"<<3*S.viscosity*ec*mobility/(2*eo*S.dielectric*kb*S.temperature)<<endl;
		cout<<"y:\t\t"<<ec*ZP/(kb*S.temperature)<<endl;}

	return mobility;}

double six_ion_problem(double ZP,double a,Solution S,Error E,bool VERBOSE)
	{// Convert Debye Length from Angstroms to Meters
	double Debye=S.debyeLength*1e-10;
	// Final Position: Solvated Radius [meters]
	double rf=a;
	double r0;
	double dr=-1e-15;	// integrate from Inf to Protein Radius
	if(Debye>a){r0=5*Debye;}
	else{r0=5*a;}
		
	if(VERBOSE){cout<<"6 ION PROBLEM\nInitial Position:\t"<<r0<<" m"<<"\nFinal Position:\t\t"<<rf<<" m"<<"\nStep Size:\t\t"<<dr<<" m"<<endl;}

	// Calculate Poisson Boltzmann Coefficient and Determine Initial Position away from Protein
	double C=calc_PBE_coefficient(ZP,r0,rf,dr,S,E,VERBOSE);
	if(VERBOSE){cout<<"PBE Coeff: "<<C<<endl;}
	if(isnan(C)){return nan("");}
	int N=7+2*S.numIon+2+3*S.numIon+3;//42;
	state_type y(N);
	// F Function
	// y[0] = f
	// y[1] = df/dr
	// y[2] = d2f/dr2
	// y[3] = d3f/dr3
	// y[4] = d4f/dr4
	// Equilibrium Electric Pontential
	// y[5] = potential
	// y[6] = d(potential)/dr
	// Induced Electrochemical Potential Function (6 ion)
	// y[7] = IEP1
	// y[8] = d(IEP1)/dr
	// y[9] = IEP2
	// y[10] = d(IEP2)/dr
	// y[11] = IEP3
	// y[12] = d(IEP3)/dr
	// y[13] = IEP4
	// y[14] = d(IEP4)/dr
	// y[15] = IEP5
	// y[16] = d(IEP5)/dr
	// y[17] = IEP6
	// y[18] = d(IEP6)/dr
	// Constants
	// y[19] = Temp					// [Kelvin]
	// y[20] = numIon
	// Ion Properties
	// y[21] = Ion1.val;		
	// y[22] = Ion1.conc;		// mM [mol/m^3]
	// y[23] = Ion1.radius;		// meters
	// y[24] = Ion2.val;
	// y[25] = Ion2.conc;
	// y[26] = Ion2.radius;
	// y[27] = Ion3.val;
	// y[28] = Ion3.conc;
	// y[29] = Ion3.radius;
	// y[30] = Ion4.val;
	// y[31] = Ion4.conc;
	// y[32] = Ion4.radius;
	// y[33] = Ion5.val;
	// y[34] = Ion5.conc;
	// y[35] = Ion5.radius;
	// y[36] = Ion6.val;
	// y[37] = Ion6.conc;
	// y[38] = Ion6.radius;
	// Solution Properties
	// y[39] = relative dielectric	[dimensionless]
	// y[40] = Debye Length [meters]
	// y[41] = viscosity  [Pa s]
	// Collect Error Data
	double abs_err=E.absolute;
	double rel_err=E.relative;
	double a_x=E.a_x;
	double a_dxdt=E.a_dxdt;
	controlled_stepper_type controlled_stepper(default_error_checker<double,range_algebra,default_operations>(abs_err,rel_err,a_x,a_dxdt));	

	// Calculate Mobility
	vector<state_type> y1,y2,y3,y4,y5,y6,y7,y8,y9,y10;
	vector<double> r1,r2,r3,r4,r5,r6,r7,r8,r9,r10;
	// Specify Boundary Conditions / Initial Values Infinitely Far from the Particle
	if(VERBOSE){cout<<"Solving Homogeneous Problem for 1st Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=1/(r0*r0);y[8]=-2/(r0*r0*r0);
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Constants
	y[19]=S.temperature;y[20]=S.numIon;
	// Fist Ion Properties
	y[21]=S.ionVal[0];y[22]=S.ionConc[0]*1000;y[23]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[24]=S.ionVal[1];y[25]=S.ionConc[1]*1000;y[26]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[27]=S.ionVal[2];y[28]=S.ionConc[2]*1000;y[29]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[30]=S.ionVal[3];y[31]=S.ionConc[3]*1000;y[32]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[33]=S.ionVal[4];y[34]=S.ionConc[4]*1000;y[35]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[36]=S.ionVal[5];y[37]=S.ionConc[5]*1000;y[38]=S.ionRad[5]*1e-10;
	// Solution Properties
	y[39]=S.dielectric;y[40]=Debye;y[41]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps1=integrate_adaptive(controlled_stepper,homogeneous_form_6ion,y,r0,rf,dr,push_back_state_and_time(y1,r1));
	int Sz1=steps1+1;
	if(VERBOSE){cout<<"Done! ("<<Sz1<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 2nd Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=1/(r0*r0);y[10]=-2/(r0*r0*r0);
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Constants
	y[19]=S.temperature;y[20]=S.numIon;
	// Fist Ion Properties
	y[21]=S.ionVal[0];y[22]=S.ionConc[0]*1000;y[23]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[24]=S.ionVal[1];y[25]=S.ionConc[1]*1000;y[26]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[27]=S.ionVal[2];y[28]=S.ionConc[2]*1000;y[29]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[30]=S.ionVal[3];y[31]=S.ionConc[3]*1000;y[32]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[33]=S.ionVal[4];y[34]=S.ionConc[4]*1000;y[35]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[36]=S.ionVal[5];y[37]=S.ionConc[5]*1000;y[38]=S.ionRad[5]*1e-10;
	// Solution Properties
	y[39]=S.dielectric;y[40]=Debye;y[41]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps2=integrate_adaptive(controlled_stepper,homogeneous_form_6ion,y,r0,rf,dr,push_back_state_and_time(y2,r2));
	int Sz2=steps2+1;
	if(VERBOSE){cout<<"Done! ("<<Sz2<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 3rd Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=1/(r0*r0);y[12]=-2/(r0*r0*r0);
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Constants
	y[19]=S.temperature;y[20]=S.numIon;
	// Fist Ion Properties
	y[21]=S.ionVal[0];y[22]=S.ionConc[0]*1000;y[23]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[24]=S.ionVal[1];y[25]=S.ionConc[1]*1000;y[26]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[27]=S.ionVal[2];y[28]=S.ionConc[2]*1000;y[29]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[30]=S.ionVal[3];y[31]=S.ionConc[3]*1000;y[32]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[33]=S.ionVal[4];y[34]=S.ionConc[4]*1000;y[35]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[36]=S.ionVal[5];y[37]=S.ionConc[5]*1000;y[38]=S.ionRad[5]*1e-10;
	// Solution Properties
	y[39]=S.dielectric;y[40]=Debye;y[41]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps3=integrate_adaptive(controlled_stepper,homogeneous_form_6ion,y,r0,rf,dr,push_back_state_and_time(y3,r3));
	int Sz3=steps3+1;
	if(VERBOSE){cout<<"Done! ("<<Sz3<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 4th Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=1/(r0*r0);y[14]=-2/(r0*r0*r0);
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Constants
	y[19]=S.temperature;y[20]=S.numIon;
	// Fist Ion Properties
	y[21]=S.ionVal[0];y[22]=S.ionConc[0]*1000;y[23]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[24]=S.ionVal[1];y[25]=S.ionConc[1]*1000;y[26]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[27]=S.ionVal[2];y[28]=S.ionConc[2]*1000;y[29]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[30]=S.ionVal[3];y[31]=S.ionConc[3]*1000;y[32]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[33]=S.ionVal[4];y[34]=S.ionConc[4]*1000;y[35]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[36]=S.ionVal[5];y[37]=S.ionConc[5]*1000;y[38]=S.ionRad[5]*1e-10;
	// Solution Properties
	y[39]=S.dielectric;y[40]=Debye;y[41]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps4=integrate_adaptive(controlled_stepper,homogeneous_form_6ion,y,r0,rf,dr,push_back_state_and_time(y4,r4));
	int Sz4=steps4+1;
	if(VERBOSE){cout<<"Done! ("<<Sz4<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 5th Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=1/(r0*r0);y[16]=-2/(r0*r0*r0);
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Constants
	y[19]=S.temperature;y[20]=S.numIon;
	// Fist Ion Properties
	y[21]=S.ionVal[0];y[22]=S.ionConc[0]*1000;y[23]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[24]=S.ionVal[1];y[25]=S.ionConc[1]*1000;y[26]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[27]=S.ionVal[2];y[28]=S.ionConc[2]*1000;y[29]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[30]=S.ionVal[3];y[31]=S.ionConc[3]*1000;y[32]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[33]=S.ionVal[4];y[34]=S.ionConc[4]*1000;y[35]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[36]=S.ionVal[5];y[37]=S.ionConc[5]*1000;y[38]=S.ionRad[5]*1e-10;
	// Solution Properties
	y[39]=S.dielectric;y[40]=Debye;y[41]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps5=integrate_adaptive(controlled_stepper,homogeneous_form_6ion,y,r0,rf,dr,push_back_state_and_time(y5,r5));
	int Sz5=steps5+1;
	if(VERBOSE){cout<<"Done! ("<<Sz5<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 6th Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=1/(r0*r0);y[18]=-2/(r0*r0*r0);
	// Constants
	y[19]=S.temperature;y[20]=S.numIon;
	// Fist Ion Properties
	y[21]=S.ionVal[0];y[22]=S.ionConc[0]*1000;y[23]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[24]=S.ionVal[1];y[25]=S.ionConc[1]*1000;y[26]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[27]=S.ionVal[2];y[28]=S.ionConc[2]*1000;y[29]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[30]=S.ionVal[3];y[31]=S.ionConc[3]*1000;y[32]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[33]=S.ionVal[4];y[34]=S.ionConc[4]*1000;y[35]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[36]=S.ionVal[5];y[37]=S.ionConc[5]*1000;y[38]=S.ionRad[5]*1e-10;
	// Solution Properties
	y[39]=S.dielectric;y[40]=Debye;y[41]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps6=integrate_adaptive(controlled_stepper,homogeneous_form_6ion,y,r0,rf,dr,push_back_state_and_time(y6,r6));
	int Sz6=steps6+1;
	if(VERBOSE){cout<<"Done! ("<<Sz6<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 1st f term...";}
	// F Function
	y[0]=r0;y[1]=1;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Constants
	y[19]=S.temperature;y[20]=S.numIon;
	// Fist Ion Properties
	y[21]=S.ionVal[0];y[22]=S.ionConc[0]*1000;y[23]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[24]=S.ionVal[1];y[25]=S.ionConc[1]*1000;y[26]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[27]=S.ionVal[2];y[28]=S.ionConc[2]*1000;y[29]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[30]=S.ionVal[3];y[31]=S.ionConc[3]*1000;y[32]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[33]=S.ionVal[4];y[34]=S.ionConc[4]*1000;y[35]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[36]=S.ionVal[5];y[37]=S.ionConc[5]*1000;y[38]=S.ionRad[5]*1e-10;
	// Solution Properties
	y[39]=S.dielectric;y[40]=Debye;y[41]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps7=integrate_adaptive(controlled_stepper,homogeneous_form_6ion,y,r0,rf,dr,push_back_state_and_time(y7,r7));
	int Sz7=steps7+1;
	if(VERBOSE){cout<<"Done! ("<<Sz7<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 2nd f term...";}
	// F Function
	y[0]=1/r0;y[1]=-1/(r0*r0);y[2]=2/(r0*r0*r0);y[3]=-6/(r0*r0*r0*r0);y[4]=24/(r0*r0*r0*r0*r0);
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Constants
	y[19]=S.temperature;y[20]=S.numIon;
	// Fist Ion Properties
	y[21]=S.ionVal[0];y[22]=S.ionConc[0]*1000;y[23]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[24]=S.ionVal[1];y[25]=S.ionConc[1]*1000;y[26]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[27]=S.ionVal[2];y[28]=S.ionConc[2]*1000;y[29]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[30]=S.ionVal[3];y[31]=S.ionConc[3]*1000;y[32]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[33]=S.ionVal[4];y[34]=S.ionConc[4]*1000;y[35]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[36]=S.ionVal[5];y[37]=S.ionConc[5]*1000;y[38]=S.ionRad[5]*1e-10;
	// Solution Properties
	y[39]=S.dielectric;y[40]=Debye;y[41]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps8=integrate_adaptive(controlled_stepper,homogeneous_form_6ion,y,r0,rf,dr,push_back_state_and_time(y8,r8));
	int Sz8=steps8+1;
	if(VERBOSE){cout<<"Done! ("<<Sz8<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Flow Field Problem...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Constants
	y[19]=S.temperature;y[20]=S.numIon;
	// Fist Ion Properties
	y[21]=S.ionVal[0];y[22]=S.ionConc[0]*1000;y[23]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[24]=S.ionVal[1];y[25]=S.ionConc[1]*1000;y[26]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[27]=S.ionVal[2];y[28]=S.ionConc[2]*1000;y[29]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[30]=S.ionVal[3];y[31]=S.ionConc[3]*1000;y[32]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[33]=S.ionVal[4];y[34]=S.ionConc[4]*1000;y[35]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[36]=S.ionVal[5];y[37]=S.ionConc[5]*1000;y[38]=S.ionRad[5]*1e-10;
	// Solution Properties
	y[39]=S.dielectric;y[40]=Debye;y[41]=S.viscosity;
	// Solve
	size_t steps9=integrate_adaptive(controlled_stepper,inhomogeneous_Prob1_6ion,y,r0,rf,dr,push_back_state_and_time(y9,r9));
	int Sz9=steps9+1;
	if(VERBOSE){cout<<"Done! ("<<Sz9<<" points)"<<endl;}

	if(VERBOSE){cout<<"Calculating Asymptotic Constants...";}
	alglib::real_2d_array A;
	alglib::real_1d_array C1,B1;
	// Define System Matrix containing coefficients of linear system
	A.setlength(8,8);
	for(int i=0;i<8;i++)
		{for(int j=0;j<8;j++)
			{A[i][j]=0;}}
	// d(IEP1)/dr
	A[0][0]=y1[Sz1-1][8];A[0][1]=y2[Sz2-1][8];A[0][2]=y3[Sz3-1][8];A[0][3]=y4[Sz4-1][8];A[0][4]=y5[Sz5-1][8];A[0][5]=y6[Sz6-1][8];A[0][6]=y7[Sz7-1][8];A[0][7]=y8[Sz8-1][8];
	// d(IEP2)/dr
	A[1][0]=y1[Sz1-1][10];A[1][1]=y2[Sz2-1][10];A[1][2]=y3[Sz3-1][10];A[1][3]=y4[Sz4-1][10];A[1][4]=y5[Sz5-1][10];A[1][5]=y6[Sz6-1][10];A[1][6]=y7[Sz7-1][10];A[1][7]=y8[Sz8-1][10];
	// d(IEP3)/dr
	A[2][0]=y1[Sz1-1][12];A[2][1]=y2[Sz2-1][12];A[2][2]=y3[Sz3-1][12];A[2][3]=y4[Sz4-1][12];A[2][4]=y5[Sz5-1][12];A[2][5]=y6[Sz6-1][12];A[2][6]=y7[Sz7-1][12];A[2][7]=y8[Sz8-1][12];
	// d(IEP4)/dr
	A[3][0]=y1[Sz1-1][14];A[3][1]=y2[Sz2-1][14];A[3][2]=y3[Sz3-1][14];A[3][3]=y4[Sz4-1][14];A[3][4]=y5[Sz5-1][14];A[3][5]=y6[Sz6-1][14];A[3][6]=y7[Sz7-1][14];A[3][7]=y8[Sz8-1][14];
	// d(IEP5)/dr
	A[4][0]=y1[Sz1-1][16];A[4][1]=y2[Sz2-1][16];A[4][2]=y3[Sz3-1][16];A[4][3]=y4[Sz4-1][16];A[4][4]=y5[Sz5-1][16];A[4][5]=y6[Sz6-1][16];A[4][6]=y7[Sz7-1][16];A[4][7]=y8[Sz8-1][16];
	// d(IEP6)/dr
	A[5][0]=y1[Sz1-1][18];A[5][1]=y2[Sz2-1][18];A[5][2]=y3[Sz3-1][18];A[5][3]=y4[Sz4-1][18];A[5][4]=y5[Sz5-1][18];A[5][5]=y6[Sz6-1][18];A[5][6]=y7[Sz7-1][18];A[5][7]=y8[Sz8-1][18];
	// f 1st (df/dr)
	A[6][0]=y1[Sz1-1][1];A[6][1]=y2[Sz2-1][1];A[6][2]=y3[Sz3-1][1];A[6][3]=y4[Sz4-1][1];A[6][4]=y5[Sz5-1][1];A[6][5]=y6[Sz6-1][1];A[6][6]=y7[Sz7-1][1];A[6][7]=y8[Sz8-1][1];
	// f 2nd (d2f/dr2)
	A[7][0]=y1[Sz1-1][2];A[7][1]=y2[Sz2-1][2];A[7][2]=y3[Sz3-1][2];A[7][3]=y4[Sz4-1][2];A[7][4]=y5[Sz5-1][2];A[7][5]=y6[Sz6-1][2];A[7][6]=y7[Sz7-1][2];A[7][7]=y8[Sz8-1][2];
	// Define Solution Matrix for Problem #1
	B1.setlength(8);
	// d(IEP1)/dr = 0
	B1[0]=0-y9[Sz9-1][8];
	// d(IEP2)/dr = 0
	B1[1]=0-y9[Sz9-1][10];
	// d(IEP3)/dr = 0
	B1[2]=0-y9[Sz9-1][12];
	// d(IEP4)/dr = 0
	B1[3]=0-y9[Sz9-1][14];
	// d(IEP5)/dr = 0
	B1[4]=0-y9[Sz9-1][16];
	// d(IEP6)/dr = 0
	B1[5]=0-y9[Sz9-1][18];
	// df/dr = -Radius/2
	B1[6]=-a/2 - y9[Sz9-1][1];
	// d2f/dr2 = -0.5
	B1[7]=-0.5 - y9[Sz9-1][2];
	// Create Output Parameters
	alglib::ae_int_t info;
	alglib::densesolverreport rep;
	alglib::rmatrixsolve(A,8,B1,info,rep,C1);
	if(VERBOSE)
		{cout<<"\n";
		cout<<C1[0]*A[0][0] + C1[1]*A[0][1] + C1[2]*A[0][2] + C1[3]*A[0][3] + C1[4]*A[0][4] + C1[5]*A[0][5] + C1[6]*A[0][6] + C1[7]*A[0][7] + y9[Sz9-1][8]<<endl;
		cout<<C1[0]*A[1][0] + C1[1]*A[1][1] + C1[2]*A[1][2] + C1[3]*A[1][3] + C1[4]*A[1][4] + C1[5]*A[1][5] + C1[6]*A[1][6] + C1[7]*A[1][7] + y9[Sz9-1][10]<<endl;
		cout<<C1[0]*A[2][0] + C1[1]*A[2][1] + C1[2]*A[2][2] + C1[3]*A[2][3] + C1[4]*A[2][4] + C1[5]*A[2][5] + C1[6]*A[2][6] + C1[7]*A[2][7] + y9[Sz9-1][12]<<endl;
		cout<<C1[0]*A[3][0] + C1[1]*A[3][1] + C1[2]*A[3][2] + C1[3]*A[3][3] + C1[4]*A[3][4] + C1[5]*A[3][5] + C1[6]*A[3][6] + C1[7]*A[3][7] + y9[Sz9-1][14]<<endl;
		cout<<C1[0]*A[4][0] + C1[1]*A[4][1] + C1[2]*A[4][2] + C1[3]*A[4][3] + C1[4]*A[4][4] + C1[5]*A[4][5] + C1[6]*A[4][6] + C1[7]*A[4][7] + y9[Sz9-1][16]<<endl;
		cout<<C1[0]*A[5][0] + C1[1]*A[5][1] + C1[2]*A[5][2] + C1[3]*A[5][3] + C1[4]*A[5][4] + C1[5]*A[5][5] + C1[6]*A[5][6] + C1[7]*A[5][7] + y9[Sz9-1][18]<<endl;
		cout<<C1[0]*A[6][0] + C1[1]*A[6][1] + C1[2]*A[6][2] + C1[3]*A[6][3] + C1[4]*A[6][4] + C1[5]*A[6][5] + C1[6]*A[6][6] + C1[7]*A[6][7] + y9[Sz9-1][1]<<endl;
		cout<<C1[0]*A[7][0] + C1[1]*A[7][1] + C1[2]*A[7][2] + C1[3]*A[7][3] + C1[4]*A[7][4] + C1[5]*A[7][5] + C1[6]*A[7][6] + C1[7]*A[7][7] + y9[Sz9-1][2]<<endl;
		cout<<"Done!"<<endl;
		cout<<"C1_1 = "<<C1[0]<<endl;
		cout<<"C2_1 = "<<C1[1]<<endl;
		cout<<"C3_1 = "<<C1[2]<<endl;
		cout<<"C4_1 = "<<C1[3]<<endl;
		cout<<"C5_1 = "<<C1[4]<<endl;
		cout<<"C6_1 = "<<C1[5]<<endl;
		cout<<"C7_1 = "<<C1[6]<<endl;
		cout<<"C8_1 = "<<C1[7]<<endl;}

	if(VERBOSE){cout<<"Solving Electric Field Problem...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Constants
	y[19]=S.temperature;y[20]=S.numIon;
	// Fist Ion Properties
	y[21]=S.ionVal[0];y[22]=S.ionConc[0]*1000;y[23]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[24]=S.ionVal[1];y[25]=S.ionConc[1]*1000;y[26]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[27]=S.ionVal[2];y[28]=S.ionConc[2]*1000;y[29]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[30]=S.ionVal[3];y[31]=S.ionConc[3]*1000;y[32]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[33]=S.ionVal[4];y[34]=S.ionConc[4]*1000;y[35]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[36]=S.ionVal[5];y[37]=S.ionConc[5]*1000;y[38]=S.ionRad[5]*1e-10;
	// Solution Properties
	y[39]=S.dielectric;y[40]=Debye;y[41]=S.viscosity;
	// Solve
	size_t steps10=integrate_adaptive(controlled_stepper,inhomogeneous_Prob2_6ion,y,r0,rf,dr,push_back_state_and_time(y10,r10));
	int Sz10=steps10+1;
	if(VERBOSE){cout<<"Done! ("<<Sz10<<" points)"<<endl;}

	if(VERBOSE){cout<<"Calculating Asymptotic Constants...";}
	alglib::real_1d_array C2,B2;
	// d(IEP1)/dr
	A[0][0]=y1[Sz1-1][8];A[0][1]=y2[Sz2-1][8];A[0][2]=y3[Sz3-1][8];A[0][3]=y4[Sz4-1][8];A[0][4]=y5[Sz5-1][8];A[0][5]=y6[Sz6-1][8];A[0][6]=y7[Sz7-1][8];A[0][7]=y8[Sz8-1][8];
	// d(IEP2)/dr
	A[1][0]=y1[Sz1-1][10];A[1][1]=y2[Sz2-1][10];A[1][2]=y3[Sz3-1][10];A[1][3]=y4[Sz4-1][10];A[1][4]=y5[Sz5-1][10];A[1][5]=y6[Sz6-1][10];A[1][6]=y7[Sz7-1][10];A[1][7]=y8[Sz8-1][10];
	// d(IEP3)/dr
	A[2][0]=y1[Sz1-1][12];A[2][1]=y2[Sz2-1][12];A[2][2]=y3[Sz3-1][12];A[2][3]=y4[Sz4-1][12];A[2][4]=y5[Sz5-1][12];A[2][5]=y6[Sz6-1][12];A[2][6]=y7[Sz7-1][12];A[2][7]=y8[Sz8-1][12];
	// d(IEP4)/dr
	A[3][0]=y1[Sz1-1][14];A[3][1]=y2[Sz2-1][14];A[3][2]=y3[Sz3-1][14];A[3][3]=y4[Sz4-1][14];A[3][4]=y5[Sz5-1][14];A[3][5]=y6[Sz6-1][14];A[3][6]=y7[Sz7-1][14];A[3][7]=y8[Sz8-1][14];
	// d(IEP5)/dr
	A[4][0]=y1[Sz1-1][16];A[4][1]=y2[Sz2-1][16];A[4][2]=y3[Sz3-1][16];A[4][3]=y4[Sz4-1][16];A[4][4]=y5[Sz5-1][16];A[4][5]=y6[Sz6-1][16];A[4][6]=y7[Sz7-1][16];A[4][7]=y8[Sz8-1][16];
	// d(IEP6)/dr
	A[5][0]=y1[Sz1-1][18];A[5][1]=y2[Sz2-1][18];A[5][2]=y3[Sz3-1][18];A[5][3]=y4[Sz4-1][18];A[5][4]=y5[Sz5-1][18];A[5][5]=y6[Sz6-1][18];A[5][6]=y7[Sz7-1][18];A[5][7]=y8[Sz8-1][18];
	// f 1st (df/dr)
	A[6][0]=y1[Sz1-1][1];A[6][1]=y2[Sz2-1][1];A[6][2]=y3[Sz3-1][1];A[6][3]=y4[Sz4-1][1];A[6][4]=y5[Sz5-1][1];A[6][5]=y6[Sz6-1][1];A[6][6]=y7[Sz7-1][1];A[6][7]=y8[Sz8-1][1];
	// f 2nd (d2f/dr2)
	A[7][0]=y1[Sz1-1][2];A[7][1]=y2[Sz2-1][2];A[7][2]=y3[Sz3-1][2];A[7][3]=y4[Sz4-1][2];A[7][4]=y5[Sz5-1][2];A[7][5]=y6[Sz6-1][2];A[7][6]=y7[Sz7-1][2];A[7][7]=y8[Sz8-1][2];
	// Define Solution Matrix for Problem #2
	B2.setlength(8);
	// d(IEP1)/dr = -1
	B2[0]=-1-y10[Sz10-1][8];
	// d(IEP2)/dr = -1
	B2[1]=-1-y10[Sz10-1][10];
	// d(IEP3)/dr = -1
	B2[2]=-1-y10[Sz10-1][12];
	// d(IEP4)/dr = -1
	B2[3]=-1-y10[Sz10-1][14];
	// d(IEP5)/dr = -1
	B2[4]=-1-y10[Sz10-1][16];
	// d(IEP6)/dr = -1
	B2[5]=-1-y10[Sz10-1][18];
	// df/dr = 0
	B2[6]=0-y10[Sz10-1][1];
	// d2f/dr2 = 0
	B2[7]=0-y10[Sz10-1][2];
	// Create Output Parameters
	alglib::rmatrixsolve(A,8,B2,info,rep,C2);
	if(VERBOSE)
		{cout<<"\n";
		cout<<C2[0]*A[0][0] + C2[1]*A[0][1] + C2[2]*A[0][2] + C2[3]*A[0][3] + C2[4]*A[0][4] + C2[5]*A[0][5] + C2[6]*A[0][6] + C2[7]*A[0][7] + y10[Sz10-1][8]<<endl;
		cout<<C2[0]*A[1][0] + C2[1]*A[1][1] + C2[2]*A[1][2] + C2[3]*A[1][3] + C2[4]*A[1][4] + C2[5]*A[1][5] + C2[6]*A[1][6] + C2[7]*A[1][7] + y10[Sz10-1][10]<<endl;
		cout<<C2[0]*A[2][0] + C2[1]*A[2][1] + C2[2]*A[2][2] + C2[3]*A[2][3] + C2[4]*A[2][4] + C2[5]*A[2][5] + C2[6]*A[2][6] + C2[7]*A[2][7] + y10[Sz10-1][12]<<endl;
		cout<<C2[0]*A[3][0] + C2[1]*A[3][1] + C2[2]*A[3][2] + C2[3]*A[3][3] + C2[4]*A[3][4] + C2[5]*A[3][5] + C2[6]*A[3][6] + C2[7]*A[3][7] + y10[Sz10-1][14]<<endl;
		cout<<C2[0]*A[4][0] + C2[1]*A[4][1] + C2[2]*A[4][2] + C2[3]*A[4][3] + C2[4]*A[4][4] + C2[5]*A[4][5] + C2[6]*A[4][6] + C2[7]*A[4][7] + y10[Sz10-1][16]<<endl;
		cout<<C2[0]*A[5][0] + C2[1]*A[5][1] + C2[2]*A[5][2] + C2[3]*A[5][3] + C2[4]*A[5][4] + C2[5]*A[5][5] + C2[6]*A[5][6] + C2[7]*A[5][7] + y10[Sz10-1][18]<<endl;
		cout<<C2[0]*A[6][0] + C2[1]*A[6][1] + C2[2]*A[6][2] + C2[3]*A[6][3] + C2[4]*A[6][4] + C2[5]*A[6][5] + C2[6]*A[6][6] + C2[7]*A[6][7] + y10[Sz10-1][1]<<endl;
		cout<<C2[0]*A[7][0] + C2[1]*A[7][1] + C2[2]*A[7][2] + C2[3]*A[7][3] + C2[4]*A[7][4] + C2[5]*A[7][5] + C2[6]*A[7][6] + C2[7]*A[7][7] + y10[Sz10-1][2]<<endl;
		cout<<"Done!"<<endl;
		cout<<"C1_2 = "<<C2[0]<<endl;
		cout<<"C2_2 = "<<C2[1]<<endl;
		cout<<"C3_2 = "<<C2[2]<<endl;
		cout<<"C4_2 = "<<C2[3]<<endl;
		cout<<"C5_2 = "<<C2[4]<<endl;
		cout<<"C6_2 = "<<C2[5]<<endl;
		cout<<"C7_2 = "<<C2[6]<<endl;
		cout<<"C8_2 = "<<C2[7]<<endl;}

	double mobility=-C2[6]/C1[6];
	
	if(VERBOSE)
		{cout<<"\n\nElectrophoretic Mobility:\t"<<mobility<<" m^2/Vs\n";
		cout<<"ka:\t\t"<<a/Debye<<endl;
		cout<<"Zeta Potential:\t"<<ZP<<" V\n";
		cout<<"E:\t\t"<<3*S.viscosity*ec*mobility/(2*eo*S.dielectric*kb*S.temperature)<<endl;
		cout<<"y:\t\t"<<ec*ZP/(kb*S.temperature)<<endl;}

	return mobility;}

double seven_ion_problem(double ZP,double a,Solution S,Error E,bool VERBOSE)
	{// Convert Debye Length from Angstroms to Meters
	double Debye=S.debyeLength*1e-10;
	// Final Position: Solvated Radius [meters]
	double rf=a;
	double r0;
	double dr=-1e-15;	// integrate from Inf to Protein Radius
	if(Debye>a){r0=5*Debye;}
	else{r0=5*a;}
		
	if(VERBOSE){cout<<"7 ION PROBLEM\nInitial Position:\t"<<r0<<" m"<<"\nFinal Position:\t\t"<<rf<<" m"<<"\nStep Size:\t\t"<<dr<<" m"<<endl;}

	// Calculate Poisson Boltzmann Coefficient and Determine Initial Position away from Protein
	double C=calc_PBE_coefficient(ZP,r0,rf,dr,S,E,VERBOSE);
	if(VERBOSE){cout<<"PBE Coeff: "<<C<<endl;}
	if(isnan(C)){return nan("");}
	int N=7+2*S.numIon+2+3*S.numIon+3;//47;
	state_type y(N);
	// F Function
	// y[0] = f
	// y[1] = df/dr
	// y[2] = d2f/dr2
	// y[3] = d3f/dr3
	// y[4] = d4f/dr4
	// Equilibrium Electric Pontential
	// y[5] = potential
	// y[6] = d(potential)/dr
	// Induced Electrochemical Potential Function (7 ion)
	// y[7] = IEP1
	// y[8] = d(IEP1)/dr
	// y[9] = IEP2
	// y[10] = d(IEP2)/dr
	// y[11] = IEP3
	// y[12] = d(IEP3)/dr
	// y[13] = IEP4
	// y[14] = d(IEP4)/dr
	// y[15] = IEP5
	// y[16] = d(IEP5)/dr
	// y[17] = IEP6
	// y[18] = d(IEP6)/dr
	// y[19] = IEP7
	// y[20] = d(IEP7)/dr
	// Constants
	// y[21] = Temp					// [Kelvin]
	// y[22] = numIon
	// Ion Properties
	// y[23] = Ion1.val;		
	// y[24] = Ion1.conc;		// mM [mol/m^3]
	// y[25] = Ion1.radius;		// meters
	// y[26] = Ion2.val;
	// y[27] = Ion2.conc;
	// y[28] = Ion2.radius;
	// y[29] = Ion3.val;
	// y[30] = Ion3.conc;
	// y[31] = Ion3.radius;
	// y[32] = Ion4.val;
	// y[33] = Ion4.conc;
	// y[34] = Ion4.radius;
	// y[35] = Ion5.val;
	// y[36] = Ion5.conc;
	// y[37] = Ion5.radius;
	// y[38] = Ion6.val;
	// y[39] = Ion6.conc;
	// y[40] = Ion6.radius;
	// y[41] = Ion7.val;
	// y[42] = Ion7.conc;
	// y[43] = Ion7.radius;
	// Solution Properties
	// y[44] = relative dielectric	[dimensionless]
	// y[45] = Debye Length [meters]
	// y[46] = viscosity  [Pa s]
	// Collect Error Data
	double abs_err=E.absolute,rel_err=E.relative,a_x=E.a_x,a_dxdt=E.a_dxdt;
	controlled_stepper_type controlled_stepper(default_error_checker<double,range_algebra,default_operations>(abs_err,rel_err,a_x,a_dxdt));	

	// Calculate Mobility
	vector<state_type> y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11;
	vector<double> r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11;
	// Specify Boundary Conditions / Initial Values Infinitely Far from the Particle
	if(VERBOSE){cout<<"Solving Homogeneous Problem for 1st Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=1/(r0*r0);y[8]=-2/(r0*r0*r0);
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Constants
	y[21]=S.temperature;y[22]=S.numIon;
	// Fist Ion Properties
	y[23]=S.ionVal[0];y[24]=S.ionConc[0]*1000;y[25]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[26]=S.ionVal[1];y[27]=S.ionConc[1]*1000;y[28]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[29]=S.ionVal[2];y[30]=S.ionConc[2]*1000;y[31]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[32]=S.ionVal[3];y[33]=S.ionConc[3]*1000;y[34]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[35]=S.ionVal[4];y[36]=S.ionConc[4]*1000;y[37]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[38]=S.ionVal[5];y[39]=S.ionConc[5]*1000;y[40]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[41]=S.ionVal[6];y[42]=S.ionConc[6]*1000;y[43]=S.ionRad[6]*1e-10;
	// Solution Properties
	y[44]=S.dielectric;y[45]=Debye;y[46]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps1=integrate_adaptive(controlled_stepper,homogeneous_form_7ion,y,r0,rf,dr,push_back_state_and_time(y1,r1));
	int Sz1=steps1+1;
	if(VERBOSE){cout<<"Done! ("<<Sz1<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 2nd Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=1/(r0*r0);y[10]=-2/(r0*r0*r0);
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Constants
	y[21]=S.temperature;y[22]=S.numIon;
	// Fist Ion Properties
	y[23]=S.ionVal[0];y[24]=S.ionConc[0]*1000;y[25]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[26]=S.ionVal[1];y[27]=S.ionConc[1]*1000;y[28]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[29]=S.ionVal[2];y[30]=S.ionConc[2]*1000;y[31]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[32]=S.ionVal[3];y[33]=S.ionConc[3]*1000;y[34]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[35]=S.ionVal[4];y[36]=S.ionConc[4]*1000;y[37]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[38]=S.ionVal[5];y[39]=S.ionConc[5]*1000;y[40]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[41]=S.ionVal[6];y[42]=S.ionConc[6]*1000;y[43]=S.ionRad[6]*1e-10;
	// Solution Properties
	y[44]=S.dielectric;y[45]=Debye;y[46]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps2=integrate_adaptive(controlled_stepper,homogeneous_form_7ion,y,r0,rf,dr,push_back_state_and_time(y2,r2));
	int Sz2=steps2+1;
	if(VERBOSE){cout<<"Done! ("<<Sz2<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 3rd Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=1/(r0*r0);y[12]=-2/(r0*r0*r0);
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Constants
	y[21]=S.temperature;y[22]=S.numIon;
	// Fist Ion Properties
	y[23]=S.ionVal[0];y[24]=S.ionConc[0]*1000;y[25]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[26]=S.ionVal[1];y[27]=S.ionConc[1]*1000;y[28]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[29]=S.ionVal[2];y[30]=S.ionConc[2]*1000;y[31]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[32]=S.ionVal[3];y[33]=S.ionConc[3]*1000;y[34]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[35]=S.ionVal[4];y[36]=S.ionConc[4]*1000;y[37]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[38]=S.ionVal[5];y[39]=S.ionConc[5]*1000;y[40]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[41]=S.ionVal[6];y[42]=S.ionConc[6]*1000;y[43]=S.ionRad[6]*1e-10;
	// Solution Properties
	y[44]=S.dielectric;y[45]=Debye;y[46]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps3=integrate_adaptive(controlled_stepper,homogeneous_form_7ion,y,r0,rf,dr,push_back_state_and_time(y3,r3));
	int Sz3=steps3+1;
	if(VERBOSE){cout<<"Done! ("<<Sz3<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 4th Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=1/(r0*r0);y[14]=-2/(r0*r0*r0);
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Constants
	y[21]=S.temperature;y[22]=S.numIon;
	// Fist Ion Properties
	y[23]=S.ionVal[0];y[24]=S.ionConc[0]*1000;y[25]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[26]=S.ionVal[1];y[27]=S.ionConc[1]*1000;y[28]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[29]=S.ionVal[2];y[30]=S.ionConc[2]*1000;y[31]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[32]=S.ionVal[3];y[33]=S.ionConc[3]*1000;y[34]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[35]=S.ionVal[4];y[36]=S.ionConc[4]*1000;y[37]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[38]=S.ionVal[5];y[39]=S.ionConc[5]*1000;y[40]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[41]=S.ionVal[6];y[42]=S.ionConc[6]*1000;y[43]=S.ionRad[6]*1e-10;
	// Solution Properties
	y[44]=S.dielectric;y[45]=Debye;y[46]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps4=integrate_adaptive(controlled_stepper,homogeneous_form_7ion,y,r0,rf,dr,push_back_state_and_time(y4,r4));
	int Sz4=steps4+1;
	if(VERBOSE){cout<<"Done! ("<<Sz4<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 5th Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=1/(r0*r0);y[16]=-2/(r0*r0*r0);
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Constants
	y[21]=S.temperature;y[22]=S.numIon;
	// Fist Ion Properties
	y[23]=S.ionVal[0];y[24]=S.ionConc[0]*1000;y[25]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[26]=S.ionVal[1];y[27]=S.ionConc[1]*1000;y[28]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[29]=S.ionVal[2];y[30]=S.ionConc[2]*1000;y[31]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[32]=S.ionVal[3];y[33]=S.ionConc[3]*1000;y[34]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[35]=S.ionVal[4];y[36]=S.ionConc[4]*1000;y[37]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[38]=S.ionVal[5];y[39]=S.ionConc[5]*1000;y[40]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[41]=S.ionVal[6];y[42]=S.ionConc[6]*1000;y[43]=S.ionRad[6]*1e-10;
	// Solution Properties
	y[44]=S.dielectric;y[45]=Debye;y[46]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps5=integrate_adaptive(controlled_stepper,homogeneous_form_7ion,y,r0,rf,dr,push_back_state_and_time(y5,r5));
	int Sz5=steps5+1;
	if(VERBOSE){cout<<"Done! ("<<Sz5<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 6th Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=1/(r0*r0);y[18]=-2/(r0*r0*r0);
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Constants
	y[21]=S.temperature;y[22]=S.numIon;
	// Fist Ion Properties
	y[23]=S.ionVal[0];y[24]=S.ionConc[0]*1000;y[25]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[26]=S.ionVal[1];y[27]=S.ionConc[1]*1000;y[28]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[29]=S.ionVal[2];y[30]=S.ionConc[2]*1000;y[31]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[32]=S.ionVal[3];y[33]=S.ionConc[3]*1000;y[34]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[35]=S.ionVal[4];y[36]=S.ionConc[4]*1000;y[37]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[38]=S.ionVal[5];y[39]=S.ionConc[5]*1000;y[40]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[41]=S.ionVal[6];y[42]=S.ionConc[6]*1000;y[43]=S.ionRad[6]*1e-10;
	// Solution Properties
	y[44]=S.dielectric;y[45]=Debye;y[46]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps6=integrate_adaptive(controlled_stepper,homogeneous_form_7ion,y,r0,rf,dr,push_back_state_and_time(y6,r6));
	int Sz6=steps6+1;
	if(VERBOSE){cout<<"Done! ("<<Sz6<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 7th Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=1/(r0*r0);y[20]=-2/(r0*r0*r0);
	// Constants
	y[21]=S.temperature;y[22]=S.numIon;
	// Fist Ion Properties
	y[23]=S.ionVal[0];y[24]=S.ionConc[0]*1000;y[25]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[26]=S.ionVal[1];y[27]=S.ionConc[1]*1000;y[28]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[29]=S.ionVal[2];y[30]=S.ionConc[2]*1000;y[31]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[32]=S.ionVal[3];y[33]=S.ionConc[3]*1000;y[34]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[35]=S.ionVal[4];y[36]=S.ionConc[4]*1000;y[37]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[38]=S.ionVal[5];y[39]=S.ionConc[5]*1000;y[40]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[41]=S.ionVal[6];y[42]=S.ionConc[6]*1000;y[43]=S.ionRad[6]*1e-10;
	// Solution Properties
	y[44]=S.dielectric;y[45]=Debye;y[46]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps7=integrate_adaptive(controlled_stepper,homogeneous_form_7ion,y,r0,rf,dr,push_back_state_and_time(y7,r7));
	int Sz7=steps7+1;
	if(VERBOSE){cout<<"Done! ("<<Sz7<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 1st f term...";}
	// F Function
	y[0]=r0;y[1]=1;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Constants
	y[21]=S.temperature;y[22]=S.numIon;
	// Fist Ion Properties
	y[23]=S.ionVal[0];y[24]=S.ionConc[0]*1000;y[25]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[26]=S.ionVal[1];y[27]=S.ionConc[1]*1000;y[28]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[29]=S.ionVal[2];y[30]=S.ionConc[2]*1000;y[31]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[32]=S.ionVal[3];y[33]=S.ionConc[3]*1000;y[34]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[35]=S.ionVal[4];y[36]=S.ionConc[4]*1000;y[37]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[38]=S.ionVal[5];y[39]=S.ionConc[5]*1000;y[40]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[41]=S.ionVal[6];y[42]=S.ionConc[6]*1000;y[43]=S.ionRad[6]*1e-10;
	// Solution Properties
	y[44]=S.dielectric;y[45]=Debye;y[46]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps8=integrate_adaptive(controlled_stepper,homogeneous_form_7ion,y,r0,rf,dr,push_back_state_and_time(y8,r8));
	int Sz8=steps8+1;
	if(VERBOSE){cout<<"Done! ("<<Sz8<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 2nd f term...";}
	// F Function
	y[0]=1/r0;y[1]=-1/(r0*r0);y[2]=2/(r0*r0*r0);y[3]=-6/(r0*r0*r0*r0);y[4]=24/(r0*r0*r0*r0*r0);
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Constants
	y[21]=S.temperature;y[22]=S.numIon;
	// Fist Ion Properties
	y[23]=S.ionVal[0];y[24]=S.ionConc[0]*1000;y[25]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[26]=S.ionVal[1];y[27]=S.ionConc[1]*1000;y[28]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[29]=S.ionVal[2];y[30]=S.ionConc[2]*1000;y[31]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[32]=S.ionVal[3];y[33]=S.ionConc[3]*1000;y[34]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[35]=S.ionVal[4];y[36]=S.ionConc[4]*1000;y[37]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[38]=S.ionVal[5];y[39]=S.ionConc[5]*1000;y[40]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[41]=S.ionVal[6];y[42]=S.ionConc[6]*1000;y[43]=S.ionRad[6]*1e-10;
	// Solution Properties
	y[44]=S.dielectric;y[45]=Debye;y[46]=S.viscosity;
	// Solve Homogeneous Form
	size_t steps9=integrate_adaptive(controlled_stepper,homogeneous_form_7ion,y,r0,rf,dr,push_back_state_and_time(y9,r9));
	int Sz9=steps9+1;
	if(VERBOSE){cout<<"Done! ("<<Sz9<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Flow Field Problem...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Constants
	y[21]=S.temperature;y[22]=S.numIon;
	// Fist Ion Properties
	y[23]=S.ionVal[0];y[24]=S.ionConc[0]*1000;y[25]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[26]=S.ionVal[1];y[27]=S.ionConc[1]*1000;y[28]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[29]=S.ionVal[2];y[30]=S.ionConc[2]*1000;y[31]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[32]=S.ionVal[3];y[33]=S.ionConc[3]*1000;y[34]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[35]=S.ionVal[4];y[36]=S.ionConc[4]*1000;y[37]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[38]=S.ionVal[5];y[39]=S.ionConc[5]*1000;y[40]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[41]=S.ionVal[6];y[42]=S.ionConc[6]*1000;y[43]=S.ionRad[6]*1e-10;
	// Solution Properties
	y[44]=S.dielectric;y[45]=Debye;y[46]=S.viscosity;
	// Solve 
	size_t steps10=integrate_adaptive(controlled_stepper,inhomogeneous_Prob1_7ion,y,r0,rf,dr,push_back_state_and_time(y10,r10));
	int Sz10=steps10+1;
	if(VERBOSE){cout<<"Done! ("<<Sz10<<" points)"<<endl;}
	
	if(VERBOSE){cout<<"Calculating Asymptotic Constants...";}
	alglib::real_2d_array A;
	alglib::real_1d_array C1,B1;
	// Define System Matrix containing coefficients of linear system
	A.setlength(9,9);
	for(int i=0;i<9;i++)
		{for(int j=0;j<9;j++)
			{A[i][j]=0;}}
	// d(IEP1)/dr
	A[0][0]=y1[Sz1-1][8];A[0][1]=y2[Sz2-1][8];A[0][2]=y3[Sz3-1][8];A[0][3]=y4[Sz4-1][8];A[0][4]=y5[Sz5-1][8];A[0][5]=y6[Sz6-1][8];A[0][6]=y7[Sz7-1][8];A[0][7]=y8[Sz8-1][8];A[0][8]=y9[Sz9-1][8];
	// d(IEP2)/dr
	A[1][0]=y1[Sz1-1][10];A[1][1]=y2[Sz2-1][10];A[1][2]=y3[Sz3-1][10];A[1][3]=y4[Sz4-1][10];A[1][4]=y5[Sz5-1][10];A[1][5]=y6[Sz6-1][10];A[1][6]=y7[Sz7-1][10];A[1][7]=y8[Sz8-1][10];A[1][8]=y9[Sz9-1][10];
	// d(IEP3)/dr
	A[2][0]=y1[Sz1-1][12];A[2][1]=y2[Sz2-1][12];A[2][2]=y3[Sz3-1][12];A[2][3]=y4[Sz4-1][12];A[2][4]=y5[Sz5-1][12];A[2][5]=y6[Sz6-1][12];A[2][6]=y7[Sz7-1][12];A[2][7]=y8[Sz8-1][12];A[2][8]=y9[Sz9-1][12];
	// d(IEP4)/dr
	A[3][0]=y1[Sz1-1][14];A[3][1]=y2[Sz2-1][14];A[3][2]=y3[Sz3-1][14];A[3][3]=y4[Sz4-1][14];A[3][4]=y5[Sz5-1][14];A[3][5]=y6[Sz6-1][14];A[3][6]=y7[Sz7-1][14];A[3][7]=y8[Sz8-1][14];A[3][8]=y9[Sz9-1][14];
	// d(IEP5)/dr
	A[4][0]=y1[Sz1-1][16];A[4][1]=y2[Sz2-1][16];A[4][2]=y3[Sz3-1][16];A[4][3]=y4[Sz4-1][16];A[4][4]=y5[Sz5-1][16];A[4][5]=y6[Sz6-1][16];A[4][6]=y7[Sz7-1][16];A[4][7]=y8[Sz8-1][16];A[4][8]=y9[Sz9-1][16];
	// d(IEP6)/dr
	A[5][0]=y1[Sz1-1][18];
	A[5][1]=y2[Sz2-1][18];
	A[5][2]=y3[Sz3-1][18];
	A[5][3]=y4[Sz4-1][18];
	A[5][4]=y5[Sz5-1][18];
	A[5][5]=y6[Sz6-1][18];
	A[5][6]=y7[Sz7-1][18];
	A[5][7]=y8[Sz8-1][18];
	A[5][8]=y9[Sz9-1][18];
	// d(IEP7)/dr
	A[6][0]=y1[Sz1-1][20];
	A[6][1]=y2[Sz2-1][20];
	A[6][2]=y3[Sz3-1][20];
	A[6][3]=y4[Sz4-1][20];
	A[6][4]=y5[Sz5-1][20];
	A[6][5]=y6[Sz6-1][20];
	A[6][6]=y7[Sz7-1][20];
	A[6][7]=y8[Sz8-1][20];
	A[6][8]=y9[Sz9-1][20];
	// f 1st (df/dr)
	A[7][0]=y1[Sz1-1][1];
	A[7][1]=y2[Sz2-1][1];
	A[7][2]=y3[Sz3-1][1];
	A[7][3]=y4[Sz4-1][1];
	A[7][4]=y5[Sz5-1][1];
	A[7][5]=y6[Sz6-1][1];
	A[7][6]=y7[Sz7-1][1];
	A[7][7]=y8[Sz8-1][1];
	A[7][8]=y9[Sz9-1][1];
	// f 2nd (d2f/dr2)
	A[8][0]=y1[Sz1-1][2];
	A[8][1]=y2[Sz2-1][2];
	A[8][2]=y3[Sz3-1][2];
	A[8][3]=y4[Sz4-1][2];
	A[8][4]=y5[Sz5-1][2];
	A[8][5]=y6[Sz6-1][2];
	A[8][6]=y7[Sz7-1][2];
	A[8][7]=y8[Sz8-1][2];
	A[8][8]=y9[Sz9-1][2];
	// Define Solution Matrix for Problem #1
	B1.setlength(9);
	// d(IEP1)/dr = 0
	B1[0]=0-y10[Sz10-1][8];
	// d(IEP2)/dr = 0
	B1[1]=0-y10[Sz10-1][10];
	// d(IEP3)/dr = 0
	B1[2]=0-y10[Sz10-1][12];
	// d(IEP4)/dr = 0
	B1[3]=0-y10[Sz10-1][14];
	// d(IEP5)/dr = 0
	B1[4]=0-y10[Sz10-1][16];
	// d(IEP6)/dr = 0
	B1[5]=0-y10[Sz10-1][18];
	// d(IEP7)/dr = 0
	B1[6]=0-y10[Sz10-1][20];
	// df/dr = -Radius/2
	B1[7]=-a/2 - y10[Sz10-1][1];
	// d2f/dr2 = -0.5
	B1[8]=-0.5 - y10[Sz10-1][2];
	// Create Output Parameters
	alglib::ae_int_t info;
	alglib::densesolverreport rep;
	alglib::rmatrixsolve(A,9,B1,info,rep,C1);
	if(VERBOSE)
		{cout<<"\n";
		cout<<C1[0]*A[0][0] + C1[1]*A[0][1] + C1[2]*A[0][2] + C1[3]*A[0][3] + C1[4]*A[0][4] + C1[5]*A[0][5] + C1[6]*A[0][6] + C1[7]*A[0][7] + C1[8]*A[0][8] + y10[Sz10-1][8]<<endl;
		cout<<C1[0]*A[1][0] + C1[1]*A[1][1] + C1[2]*A[1][2] + C1[3]*A[1][3] + C1[4]*A[1][4] + C1[5]*A[1][5] + C1[6]*A[1][6] + C1[7]*A[1][7] + C1[8]*A[1][8] + y10[Sz10-1][10]<<endl;
		cout<<C1[0]*A[2][0] + C1[1]*A[2][1] + C1[2]*A[2][2] + C1[3]*A[2][3] + C1[4]*A[2][4] + C1[5]*A[2][5] + C1[6]*A[2][6] + C1[7]*A[2][7] + C1[8]*A[2][8] + y10[Sz10-1][12]<<endl;
		cout<<C1[0]*A[3][0] + C1[1]*A[3][1] + C1[2]*A[3][2] + C1[3]*A[3][3] + C1[4]*A[3][4] + C1[5]*A[3][5] + C1[6]*A[3][6] + C1[7]*A[3][7] + C1[8]*A[3][8] + y10[Sz10-1][14]<<endl;
		cout<<C1[0]*A[4][0] + C1[1]*A[4][1] + C1[2]*A[4][2] + C1[3]*A[4][3] + C1[4]*A[4][4] + C1[5]*A[4][5] + C1[6]*A[4][6] + C1[7]*A[4][7] + C1[8]*A[4][8] + y10[Sz10-1][16]<<endl;
		cout<<C1[0]*A[5][0] + C1[1]*A[5][1] + C1[2]*A[5][2] + C1[3]*A[5][3] + C1[4]*A[5][4] + C1[5]*A[5][5] + C1[6]*A[5][6] + C1[7]*A[5][7] + C1[8]*A[5][8] + y10[Sz10-1][18]<<endl;
		cout<<C1[0]*A[6][0] + C1[1]*A[6][1] + C1[2]*A[6][2] + C1[3]*A[6][3] + C1[4]*A[6][4] + C1[5]*A[6][5] + C1[6]*A[6][6] + C1[7]*A[6][7] + C1[8]*A[6][8] + y10[Sz10-1][20]<<endl;
		cout<<C1[0]*A[7][0] + C1[1]*A[7][1] + C1[2]*A[7][2] + C1[3]*A[7][3] + C1[4]*A[7][4] + C1[5]*A[7][5] + C1[6]*A[7][6] + C1[7]*A[7][7] + C1[8]*A[7][8] + y10[Sz10-1][1]<<endl;
		cout<<C1[0]*A[8][0] + C1[1]*A[8][1] + C1[2]*A[8][2] + C1[3]*A[8][3] + C1[4]*A[8][4] + C1[5]*A[8][5] + C1[6]*A[8][6] + C1[7]*A[8][7] + C1[8]*A[8][8] + y10[Sz10-1][2]<<endl;
		cout<<"Done!"<<endl;
		cout<<"C1_1 = "<<C1[0]<<endl;
		cout<<"C2_1 = "<<C1[1]<<endl;
		cout<<"C3_1 = "<<C1[2]<<endl;
		cout<<"C4_1 = "<<C1[3]<<endl;
		cout<<"C5_1 = "<<C1[4]<<endl;
		cout<<"C6_1 = "<<C1[5]<<endl;
		cout<<"C7_1 = "<<C1[6]<<endl;
		cout<<"C8_1 = "<<C1[7]<<endl;
		cout<<"C9_1 = "<<C1[8]<<endl;}

	if(VERBOSE){cout<<"Solving Electric Field Problem...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Constants
	y[21]=S.temperature;y[22]=S.numIon;
	// Fist Ion Properties
	y[23]=S.ionVal[0];y[24]=S.ionConc[0]*1000;y[25]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[26]=S.ionVal[1];y[27]=S.ionConc[1]*1000;y[28]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[29]=S.ionVal[2];y[30]=S.ionConc[2]*1000;y[31]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[32]=S.ionVal[3];y[33]=S.ionConc[3]*1000;y[34]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[35]=S.ionVal[4];y[36]=S.ionConc[4]*1000;y[37]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[38]=S.ionVal[5];y[39]=S.ionConc[5]*1000;y[40]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[41]=S.ionVal[6];y[42]=S.ionConc[6]*1000;y[43]=S.ionRad[6]*1e-10;
	// Solution Properties
	y[44]=S.dielectric;y[45]=Debye;y[46]=S.viscosity;
	// Solve 
	size_t steps11=integrate_adaptive(controlled_stepper,inhomogeneous_Prob2_7ion,y,r0,rf,dr,push_back_state_and_time(y11,r11));
	int Sz11=steps11+1;
	if(VERBOSE){cout<<"Done! ("<<Sz11<<" points)"<<endl;}

	if(VERBOSE){cout<<"Calculating Asymptotic Constants...";}
	alglib::real_1d_array C2,B2;
	// d(IEP1)/dr
	A[0][0]=y1[Sz1-1][8];A[0][1]=y2[Sz2-1][8];A[0][2]=y3[Sz3-1][8];A[0][3]=y4[Sz4-1][8];A[0][4]=y5[Sz5-1][8];A[0][5]=y6[Sz6-1][8];A[0][6]=y7[Sz7-1][8];A[0][7]=y8[Sz8-1][8];A[0][8]=y9[Sz9-1][8];
	// d(IEP2)/dr
	A[1][0]=y1[Sz1-1][10];A[1][1]=y2[Sz2-1][10];A[1][2]=y3[Sz3-1][10];A[1][3]=y4[Sz4-1][10];A[1][4]=y5[Sz5-1][10];A[1][5]=y6[Sz6-1][10];A[1][6]=y7[Sz7-1][10];A[1][7]=y8[Sz8-1][10];A[1][8]=y9[Sz9-1][10];
	// d(IEP3)/dr
	A[2][0]=y1[Sz1-1][12];A[2][1]=y2[Sz2-1][12];A[2][2]=y3[Sz3-1][12];A[2][3]=y4[Sz4-1][12];A[2][4]=y5[Sz5-1][12];A[2][5]=y6[Sz6-1][12];A[2][6]=y7[Sz7-1][12];A[2][7]=y8[Sz8-1][12];A[2][8]=y9[Sz9-1][12];
	// d(IEP4)/dr
	A[3][0]=y1[Sz1-1][14];A[3][1]=y2[Sz2-1][14];A[3][2]=y3[Sz3-1][14];A[3][3]=y4[Sz4-1][14];A[3][4]=y5[Sz5-1][14];A[3][5]=y6[Sz6-1][14];A[3][6]=y7[Sz7-1][14];A[3][7]=y8[Sz8-1][14];A[3][8]=y9[Sz9-1][14];
	// d(IEP5)/dr
	A[4][0]=y1[Sz1-1][16];A[4][1]=y2[Sz2-1][16];A[4][2]=y3[Sz3-1][16];A[4][3]=y4[Sz4-1][16];A[4][4]=y5[Sz5-1][16];A[4][5]=y6[Sz6-1][16];A[4][6]=y7[Sz7-1][16];A[4][7]=y8[Sz8-1][16];A[4][8]=y9[Sz9-1][16];
	// d(IEP6)/dr
	A[5][0]=y1[Sz1-1][18];A[5][1]=y2[Sz2-1][18];A[5][2]=y3[Sz3-1][18];A[5][3]=y4[Sz4-1][18];A[5][4]=y5[Sz5-1][18];A[5][5]=y6[Sz6-1][18];A[5][6]=y7[Sz7-1][18];A[5][7]=y8[Sz8-1][18];A[5][8]=y9[Sz9-1][18];
	// d(IEP7)/dr
	A[6][0]=y1[Sz1-1][20];
	A[6][1]=y2[Sz2-1][20];
	A[6][2]=y3[Sz3-1][20];
	A[6][3]=y4[Sz4-1][20];
	A[6][4]=y5[Sz5-1][20];
	A[6][5]=y6[Sz6-1][20];
	A[6][6]=y7[Sz7-1][20];
	A[6][7]=y8[Sz8-1][20];
	A[6][8]=y9[Sz9-1][20];
	// f 1st (df/dr)
	A[7][0]=y1[Sz1-1][1];
	A[7][1]=y2[Sz2-1][1];
	A[7][2]=y3[Sz3-1][1];
	A[7][3]=y4[Sz4-1][1];
	A[7][4]=y5[Sz5-1][1];
	A[7][5]=y6[Sz6-1][1];
	A[7][6]=y7[Sz7-1][1];
	A[7][7]=y8[Sz8-1][1];
	A[7][8]=y9[Sz9-1][1];
	// f 2nd (d2f/dr2)
	A[8][0]=y1[Sz1-1][2];
	A[8][1]=y2[Sz2-1][2];
	A[8][2]=y3[Sz3-1][2];
	A[8][3]=y4[Sz4-1][2];
	A[8][4]=y5[Sz5-1][2];
	A[8][5]=y6[Sz6-1][2];
	A[8][6]=y7[Sz7-1][2];
	A[8][7]=y8[Sz8-1][2];
	A[8][8]=y9[Sz9-1][2];
	// Define Solution Matrix for Problem #2
	B2.setlength(9);
	// d(IEP1)/dr = -1
	B2[0]=-1-y11[Sz11-1][8];
	// d(IEP2)/dr = -1
	B2[1]=-1-y11[Sz11-1][10];
	// d(IEP3)/dr = -1
	B2[2]=-1-y11[Sz11-1][12];
	// d(IEP4)/dr = -1
	B2[3]=-1-y11[Sz11-1][14];
	// d(IEP5)/dr = -1
	B2[4]=-1-y11[Sz11-1][16];
	// d(IEP6)/dr = -1
	B2[5]=-1-y11[Sz11-1][18];
	// d(IEP7)/dr = -1
	B2[6]=-1-y11[Sz11-1][20];
	// df/dr = 0
	B2[7]=0-y11[Sz11-1][1];
	// d2f/dr2 = 0
	B2[8]=0-y11[Sz11-1][2];
	// Create Output Parameters
	alglib::rmatrixsolve(A,9,B2,info,rep,C2);
	if(VERBOSE)
		{cout<<"\n";
		cout<<C2[0]*A[0][0] + C2[1]*A[0][1] + C2[2]*A[0][2] + C2[3]*A[0][3] + C2[4]*A[0][4] + C2[5]*A[0][5] + C2[6]*A[0][6] + C2[7]*A[0][7] + C2[8]*A[0][8] + y11[Sz11-1][8]<<endl;
		cout<<C2[0]*A[1][0] + C2[1]*A[1][1] + C2[2]*A[1][2] + C2[3]*A[1][3] + C2[4]*A[1][4] + C2[5]*A[1][5] + C2[6]*A[1][6] + C2[7]*A[1][7] + C2[8]*A[1][8] + y11[Sz11-1][10]<<endl;
		cout<<C2[0]*A[2][0] + C2[1]*A[2][1] + C2[2]*A[2][2] + C2[3]*A[2][3] + C2[4]*A[2][4] + C2[5]*A[2][5] + C2[6]*A[2][6] + C2[7]*A[2][7] + C2[8]*A[2][8] + y11[Sz11-1][12]<<endl;
		cout<<C2[0]*A[3][0] + C2[1]*A[3][1] + C2[2]*A[3][2] + C2[3]*A[3][3] + C2[4]*A[3][4] + C2[5]*A[3][5] + C2[6]*A[3][6] + C2[7]*A[3][7] + C2[8]*A[3][8] + y11[Sz11-1][14]<<endl;
		cout<<C2[0]*A[4][0] + C2[1]*A[4][1] + C2[2]*A[4][2] + C2[3]*A[4][3] + C2[4]*A[4][4] + C2[5]*A[4][5] + C2[6]*A[4][6] + C2[7]*A[4][7] + C2[8]*A[4][8] + y11[Sz11-1][16]<<endl;
		cout<<C2[0]*A[5][0] + C2[1]*A[5][1] + C2[2]*A[5][2] + C2[3]*A[5][3] + C2[4]*A[5][4] + C2[5]*A[5][5] + C2[6]*A[5][6] + C2[7]*A[5][7] + C2[8]*A[5][8] + y11[Sz11-1][18]<<endl;
		cout<<C2[0]*A[6][0] + C2[1]*A[6][1] + C2[2]*A[6][2] + C2[3]*A[6][3] + C2[4]*A[6][4] + C2[5]*A[6][5] + C2[6]*A[6][6] + C2[7]*A[6][7] + C2[8]*A[6][8] + y11[Sz11-1][20]<<endl;
		cout<<C2[0]*A[7][0] + C2[1]*A[7][1] + C2[2]*A[7][2] + C2[3]*A[7][3] + C2[4]*A[7][4] + C2[5]*A[7][5] + C2[6]*A[7][6] + C2[7]*A[7][7] + C2[8]*A[7][8] + y11[Sz11-1][1]<<endl;
		cout<<C2[0]*A[8][0] + C2[1]*A[8][1] + C2[2]*A[8][2] + C2[3]*A[8][3] + C2[4]*A[8][4] + C2[5]*A[8][5] + C2[6]*A[8][6] + C2[7]*A[8][7] + C2[8]*A[8][8] + y11[Sz11-1][2]<<endl;
		cout<<"Done!"<<endl;
		cout<<"C1_2 = "<<C2[0]<<endl;
		cout<<"C2_2 = "<<C2[1]<<endl;
		cout<<"C3_2 = "<<C2[2]<<endl;
		cout<<"C4_2 = "<<C2[3]<<endl;
		cout<<"C5_2 = "<<C2[4]<<endl;
		cout<<"C6_2 = "<<C2[5]<<endl;
		cout<<"C7_2 = "<<C2[6]<<endl;
		cout<<"C8_2 = "<<C2[7]<<endl;
		cout<<"C9_2 = "<<C2[8]<<endl;}

	double mobility=-C2[7]/C1[7];
	
	if(VERBOSE)
		{cout<<"\n\nElectrophoretic Mobility:\t"<<mobility<<" m^2/Vs\n";
		cout<<"ka:\t\t"<<a/Debye<<endl;
		cout<<"Zeta Potential:\t"<<ZP<<" V\n";
		cout<<"E:\t\t"<<3*S.viscosity*ec*mobility/(2*eo*S.dielectric*kb*S.temperature)<<endl;
		cout<<"y:\t\t"<<ec*ZP/(kb*S.temperature)<<endl;}

	return mobility;}

double eight_ion_problem(double ZP,double a,Solution S,Error E,bool VERBOSE)
	{// Convert Debye Length from Angstroms to Meters
	double Debye=S.debyeLength*1e-10;
	// Final Position: Solvated Radius [meters]
	double rf=a;
	double r0;
	double dr=-1e-15;	// integrate from Inf to Protein Radius
	if(Debye>a){r0=5*Debye;}
	else{r0=5*a;}
		
	if(VERBOSE){cout<<"8 ION PROBLEM\nInitial Position:\t"<<r0<<" m"<<"\nFinal Position:\t\t"<<rf<<" m"<<"\nStep Size:\t\t"<<dr<<" m"<<endl;}

	// Calculate Poisson Boltzmann Coefficient and Determine Initial Position away from Protein
	double C=calc_PBE_coefficient(ZP,r0,rf,dr,S,E,VERBOSE);
	if(VERBOSE){cout<<"PBE Coeff: "<<C<<endl;}
	if(isnan(C)){return nan("");}
	int N=7+2*S.numIon+2+3*S.numIon+3;//52;
	state_type y(N);
	// F Function
	// y[0] = f
	// y[1] = df/dr
	// y[2] = d2f/dr2
	// y[3] = d3f/dr3
	// y[4] = d4f/dr4
	// Equilibrium Electric Pontential
	// y[5] = potential
	// y[6] = d(potential)/dr
	// Induced Electrochemical Potential Function (8 ion)
	// y[7] = IEP1
	// y[8] = d(IEP1)/dr
	// y[9] = IEP2
	// y[10] = d(IEP2)/dr
	// y[11] = IEP3
	// y[12] = d(IEP3)/dr
	// y[13] = IEP4
	// y[14] = d(IEP4)/dr
	// y[15] = IEP5
	// y[16] = d(IEP5)/dr
	// y[17] = IEP6
	// y[18] = d(IEP6)/dr
	// y[19] = IEP7
	// y[20] = d(IEP7)/dr
	// y[21] = IEP8
	// y[22] = d(IEP8)/dr
	// Constants
	// y[23] = Temp					// [Kelvin]
	// y[24] = numIon
	// Ion Properties
	// y[25] = Ion1.val;		
	// y[26] = Ion1.conc;		// mM [mol/m^3]
	// y[27] = Ion1.radius;		// meters
	// y[28] = Ion2.val;
	// y[29] = Ion2.conc;
	// y[30] = Ion2.radius;
	// y[31] = Ion3.val;
	// y[32] = Ion3.conc;
	// y[33] = Ion3.radius;
	// y[34] = Ion4.val;
	// y[35] = Ion4.conc;
	// y[36] = Ion4.radius;
	// y[37] = Ion5.val;
	// y[38] = Ion5.conc;
	// y[39] = Ion5.radius;
	// y[40] = Ion6.val;
	// y[41] = Ion6.conc;
	// y[42] = Ion6.radius;
	// y[43] = Ion7.val;
	// y[44] = Ion7.conc;
	// y[45] = Ion7.radius;
	// y[46] = Ion8.val;
	// y[47] = Ion8.conc;
	// y[48] = Ion8.radius;
	// Solution Properties
	// y[49] = relative dielectric	[dimensionless]
	// y[50] = Debye Length [meters]
	// y[51] = viscosity  [Pa s]
	// Collect Error Data
	double abs_err=E.absolute,rel_err=E.relative,a_x=E.a_x,a_dxdt=E.a_dxdt;
	controlled_stepper_type controlled_stepper(default_error_checker<double,range_algebra,default_operations>(abs_err,rel_err,a_x,a_dxdt));	

	// Calculate Mobility
	vector<state_type> y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12;
	vector<double> r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12;
	// Specify Boundary Conditions / Initial Values Infinitely Far from the Particle
	if(VERBOSE){cout<<"Solving Homogeneous Problem for 1st Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=1/(r0*r0);y[8]=-2/(r0*r0*r0);
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Eigth Ion
	y[21]=0;y[22]=0;
	// Constants
	y[23]=S.temperature;y[24]=S.numIon;
	// Fist Ion Properties
	y[25]=S.ionVal[0];y[26]=S.ionConc[0]*1000;y[27]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[28]=S.ionVal[1];y[29]=S.ionConc[1]*1000;y[30]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[31]=S.ionVal[2];y[32]=S.ionConc[2]*1000;y[33]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[34]=S.ionVal[3];y[35]=S.ionConc[3]*1000;y[36]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[37]=S.ionVal[4];y[38]=S.ionConc[4]*1000;y[39]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[40]=S.ionVal[5];y[41]=S.ionConc[5]*1000;y[42]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[43]=S.ionVal[6];y[44]=S.ionConc[6]*1000;y[45]=S.ionRad[6]*1e-10;
	// Eigth Ion Properties
	y[46]=S.ionVal[7];y[47]=S.ionConc[7]*1000;y[48]=S.ionRad[7]*1e-10;
	// Solution Properties
	y[49]=S.dielectric;y[50]=Debye;y[51]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps1=integrate_adaptive(controlled_stepper,homogeneous_form_8ion,y,r0,rf,dr,push_back_state_and_time(y1,r1));
	int Sz1=steps1+1;
	if(VERBOSE){cout<<"Done! ("<<Sz1<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 2nd Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=1/(r0*r0);y[10]=-2/(r0*r0*r0);
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Eigth Ion
	y[21]=0;y[22]=0;
	// Constants
	y[23]=S.temperature;y[24]=S.numIon;
	// Fist Ion Properties
	y[25]=S.ionVal[0];y[26]=S.ionConc[0]*1000;y[27]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[28]=S.ionVal[1];y[29]=S.ionConc[1]*1000;y[30]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[31]=S.ionVal[2];y[32]=S.ionConc[2]*1000;y[33]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[34]=S.ionVal[3];y[35]=S.ionConc[3]*1000;y[36]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[37]=S.ionVal[4];y[38]=S.ionConc[4]*1000;y[39]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[40]=S.ionVal[5];y[41]=S.ionConc[5]*1000;y[42]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[43]=S.ionVal[6];y[44]=S.ionConc[6]*1000;y[45]=S.ionRad[6]*1e-10;
	// Eigth Ion Properties
	y[46]=S.ionVal[7];y[47]=S.ionConc[7]*1000;y[48]=S.ionRad[7]*1e-10;
	// Solution Properties
	y[49]=S.dielectric;y[50]=Debye;y[51]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps2=integrate_adaptive(controlled_stepper,homogeneous_form_8ion,y,r0,rf,dr,push_back_state_and_time(y2,r2));
	int Sz2=steps2+1;
	if(VERBOSE){cout<<"Done! ("<<Sz2<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 3rd Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=1/(r0*r0);y[12]=-2/(r0*r0*r0);
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Eigth Ion
	y[21]=0;y[22]=0;
	// Constants
	y[23]=S.temperature;y[24]=S.numIon;
	// Fist Ion Properties
	y[25]=S.ionVal[0];y[26]=S.ionConc[0]*1000;y[27]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[28]=S.ionVal[1];y[29]=S.ionConc[1]*1000;y[30]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[31]=S.ionVal[2];y[32]=S.ionConc[2]*1000;y[33]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[34]=S.ionVal[3];y[35]=S.ionConc[3]*1000;y[36]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[37]=S.ionVal[4];y[38]=S.ionConc[4]*1000;y[39]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[40]=S.ionVal[5];y[41]=S.ionConc[5]*1000;y[42]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[43]=S.ionVal[6];y[44]=S.ionConc[6]*1000;y[45]=S.ionRad[6]*1e-10;
	// Eigth Ion Properties
	y[46]=S.ionVal[7];y[47]=S.ionConc[7]*1000;y[48]=S.ionRad[7]*1e-10;
	// Solution Properties
	y[49]=S.dielectric;y[50]=Debye;y[51]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps3=integrate_adaptive(controlled_stepper,homogeneous_form_8ion,y,r0,rf,dr,push_back_state_and_time(y3,r3));
	int Sz3=steps3+1;
	if(VERBOSE){cout<<"Done! ("<<Sz3<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 4th Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=1/(r0*r0);y[14]=-2/(r0*r0*r0);
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Eigth Ion
	y[21]=0;y[22]=0;
	// Constants
	y[23]=S.temperature;y[24]=S.numIon;
	// Fist Ion Properties
	y[25]=S.ionVal[0];y[26]=S.ionConc[0]*1000;y[27]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[28]=S.ionVal[1];y[29]=S.ionConc[1]*1000;y[30]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[31]=S.ionVal[2];y[32]=S.ionConc[2]*1000;y[33]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[34]=S.ionVal[3];y[35]=S.ionConc[3]*1000;y[36]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[37]=S.ionVal[4];y[38]=S.ionConc[4]*1000;y[39]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[40]=S.ionVal[5];y[41]=S.ionConc[5]*1000;y[42]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[43]=S.ionVal[6];y[44]=S.ionConc[6]*1000;y[45]=S.ionRad[6]*1e-10;
	// Eigth Ion Properties
	y[46]=S.ionVal[7];y[47]=S.ionConc[7]*1000;y[48]=S.ionRad[7]*1e-10;
	// Solution Properties
	y[49]=S.dielectric;y[50]=Debye;y[51]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps4=integrate_adaptive(controlled_stepper,homogeneous_form_8ion,y,r0,rf,dr,push_back_state_and_time(y4,r4));
	int Sz4=steps4+1;
	if(VERBOSE){cout<<"Done! ("<<Sz4<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 5th Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=1/(r0*r0);y[16]=-2/(r0*r0*r0);
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Eigth Ion
	y[21]=0;y[22]=0;
	// Constants
	y[23]=S.temperature;y[24]=S.numIon;
	// Fist Ion Properties
	y[25]=S.ionVal[0];y[26]=S.ionConc[0]*1000;y[27]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[28]=S.ionVal[1];y[29]=S.ionConc[1]*1000;y[30]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[31]=S.ionVal[2];y[32]=S.ionConc[2]*1000;y[33]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[34]=S.ionVal[3];y[35]=S.ionConc[3]*1000;y[36]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[37]=S.ionVal[4];y[38]=S.ionConc[4]*1000;y[39]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[40]=S.ionVal[5];y[41]=S.ionConc[5]*1000;y[42]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[43]=S.ionVal[6];y[44]=S.ionConc[6]*1000;y[45]=S.ionRad[6]*1e-10;
	// Eigth Ion Properties
	y[46]=S.ionVal[7];y[47]=S.ionConc[7]*1000;y[48]=S.ionRad[7]*1e-10;
	// Solution Properties
	y[49]=S.dielectric;y[50]=Debye;y[51]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps5=integrate_adaptive(controlled_stepper,homogeneous_form_8ion,y,r0,rf,dr,push_back_state_and_time(y5,r5));
	int Sz5=steps5+1;
	if(VERBOSE){cout<<"Done! ("<<Sz5<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 6th Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=1/(r0*r0);y[18]=-2/(r0*r0*r0);
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Eigth Ion
	y[21]=0;y[22]=0;
	// Constants
	y[23]=S.temperature;y[24]=S.numIon;
	// Fist Ion Properties
	y[25]=S.ionVal[0];y[26]=S.ionConc[0]*1000;y[27]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[28]=S.ionVal[1];y[29]=S.ionConc[1]*1000;y[30]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[31]=S.ionVal[2];y[32]=S.ionConc[2]*1000;y[33]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[34]=S.ionVal[3];y[35]=S.ionConc[3]*1000;y[36]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[37]=S.ionVal[4];y[38]=S.ionConc[4]*1000;y[39]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[40]=S.ionVal[5];y[41]=S.ionConc[5]*1000;y[42]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[43]=S.ionVal[6];y[44]=S.ionConc[6]*1000;y[45]=S.ionRad[6]*1e-10;
	// Eigth Ion Properties
	y[46]=S.ionVal[7];y[47]=S.ionConc[7]*1000;y[48]=S.ionRad[7]*1e-10;
	// Solution Properties
	y[49]=S.dielectric;y[50]=Debye;y[51]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps6=integrate_adaptive(controlled_stepper,homogeneous_form_8ion,y,r0,rf,dr,push_back_state_and_time(y6,r6));
	int Sz6=steps6+1;
	if(VERBOSE){cout<<"Done! ("<<Sz6<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 7th Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=1/(r0*r0);y[20]=-2/(r0*r0*r0);
	// Eigth Ion
	y[21]=0;y[22]=0;
	// Constants
	y[23]=S.temperature;y[24]=S.numIon;
	// Fist Ion Properties
	y[25]=S.ionVal[0];y[26]=S.ionConc[0]*1000;y[27]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[28]=S.ionVal[1];y[29]=S.ionConc[1]*1000;y[30]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[31]=S.ionVal[2];y[32]=S.ionConc[2]*1000;y[33]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[34]=S.ionVal[3];y[35]=S.ionConc[3]*1000;y[36]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[37]=S.ionVal[4];y[38]=S.ionConc[4]*1000;y[39]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[40]=S.ionVal[5];y[41]=S.ionConc[5]*1000;y[42]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[43]=S.ionVal[6];y[44]=S.ionConc[6]*1000;y[45]=S.ionRad[6]*1e-10;
	// Eigth Ion Properties
	y[46]=S.ionVal[7];y[47]=S.ionConc[7]*1000;y[48]=S.ionRad[7]*1e-10;
	// Solution Properties
	y[49]=S.dielectric;y[50]=Debye;y[51]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps7=integrate_adaptive(controlled_stepper,homogeneous_form_8ion,y,r0,rf,dr,push_back_state_and_time(y7,r7));
	int Sz7=steps7+1;
	if(VERBOSE){cout<<"Done! ("<<Sz7<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 8th Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Eigth Ion
	y[21]=1/(r0*r0);y[22]=-2/(r0*r0*r0);
	// Constants
	y[23]=S.temperature;y[24]=S.numIon;
	// Fist Ion Properties
	y[25]=S.ionVal[0];y[26]=S.ionConc[0]*1000;y[27]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[28]=S.ionVal[1];y[29]=S.ionConc[1]*1000;y[30]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[31]=S.ionVal[2];y[32]=S.ionConc[2]*1000;y[33]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[34]=S.ionVal[3];y[35]=S.ionConc[3]*1000;y[36]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[37]=S.ionVal[4];y[38]=S.ionConc[4]*1000;y[39]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[40]=S.ionVal[5];y[41]=S.ionConc[5]*1000;y[42]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[43]=S.ionVal[6];y[44]=S.ionConc[6]*1000;y[45]=S.ionRad[6]*1e-10;
	// Eigth Ion Properties
	y[46]=S.ionVal[7];y[47]=S.ionConc[7]*1000;y[48]=S.ionRad[7]*1e-10;
	// Solution Properties
	y[49]=S.dielectric;y[50]=Debye;y[51]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps8=integrate_adaptive(controlled_stepper,homogeneous_form_8ion,y,r0,rf,dr,push_back_state_and_time(y8,r8));
	int Sz8=steps8+1;
	if(VERBOSE){cout<<"Done! ("<<Sz8<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 1st f term...";}
	// F Function
	y[0]=r0;y[1]=1;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Eigth Ion
	y[21]=0;y[22]=0;
	// Constants
	y[23]=S.temperature;y[24]=S.numIon;
	// Fist Ion Properties
	y[25]=S.ionVal[0];y[26]=S.ionConc[0]*1000;y[27]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[28]=S.ionVal[1];y[29]=S.ionConc[1]*1000;y[30]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[31]=S.ionVal[2];y[32]=S.ionConc[2]*1000;y[33]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[34]=S.ionVal[3];y[35]=S.ionConc[3]*1000;y[36]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[37]=S.ionVal[4];y[38]=S.ionConc[4]*1000;y[39]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[40]=S.ionVal[5];y[41]=S.ionConc[5]*1000;y[42]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[43]=S.ionVal[6];y[44]=S.ionConc[6]*1000;y[45]=S.ionRad[6]*1e-10;
	// Eigth Ion Properties
	y[46]=S.ionVal[7];y[47]=S.ionConc[7]*1000;y[48]=S.ionRad[7]*1e-10;
	// Solution Properties
	y[49]=S.dielectric;y[50]=Debye;y[51]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps9=integrate_adaptive(controlled_stepper,homogeneous_form_8ion,y,r0,rf,dr,push_back_state_and_time(y9,r9));
	int Sz9=steps9+1;
	if(VERBOSE){cout<<"Done! ("<<Sz9<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 2nd f term...";}
	// F Function
	y[0]=1/r0;y[1]=-1/(r0*r0);y[2]=2/(r0*r0*r0);y[3]=-6/(r0*r0*r0*r0);y[4]=24/(r0*r0*r0*r0*r0);
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Eigth Ion
	y[21]=0;y[22]=0;
	// Constants
	y[23]=S.temperature;y[24]=S.numIon;
	// Fist Ion Properties
	y[25]=S.ionVal[0];y[26]=S.ionConc[0]*1000;y[27]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[28]=S.ionVal[1];y[29]=S.ionConc[1]*1000;y[30]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[31]=S.ionVal[2];y[32]=S.ionConc[2]*1000;y[33]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[34]=S.ionVal[3];y[35]=S.ionConc[3]*1000;y[36]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[37]=S.ionVal[4];y[38]=S.ionConc[4]*1000;y[39]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[40]=S.ionVal[5];y[41]=S.ionConc[5]*1000;y[42]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[43]=S.ionVal[6];y[44]=S.ionConc[6]*1000;y[45]=S.ionRad[6]*1e-10;
	// Eigth Ion Properties
	y[46]=S.ionVal[7];y[47]=S.ionConc[7]*1000;y[48]=S.ionRad[7]*1e-10;
	// Solution Properties
	y[49]=S.dielectric;y[50]=Debye;y[51]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps10=integrate_adaptive(controlled_stepper,homogeneous_form_8ion,y,r0,rf,dr,push_back_state_and_time(y10,r10));
	int Sz10=steps10+1;
	if(VERBOSE){cout<<"Done! ("<<Sz10<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Flow Field Problem...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Eigth Ion
	y[21]=0;y[22]=0;
	// Constants
	y[23]=S.temperature;y[24]=S.numIon;
	// Fist Ion Properties
	y[25]=S.ionVal[0];y[26]=S.ionConc[0]*1000;y[27]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[28]=S.ionVal[1];y[29]=S.ionConc[1]*1000;y[30]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[31]=S.ionVal[2];y[32]=S.ionConc[2]*1000;y[33]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[34]=S.ionVal[3];y[35]=S.ionConc[3]*1000;y[36]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[37]=S.ionVal[4];y[38]=S.ionConc[4]*1000;y[39]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[40]=S.ionVal[5];y[41]=S.ionConc[5]*1000;y[42]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[43]=S.ionVal[6];y[44]=S.ionConc[6]*1000;y[45]=S.ionRad[6]*1e-10;
	// Eigth Ion Properties
	y[46]=S.ionVal[7];y[47]=S.ionConc[7]*1000;y[48]=S.ionRad[7]*1e-10;
	// Solution Properties
	y[49]=S.dielectric;y[50]=Debye;y[51]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps11=integrate_adaptive(controlled_stepper,inhomogeneous_Prob1_8ion,y,r0,rf,dr,push_back_state_and_time(y11,r11));
	int Sz11=steps11+1;
	if(VERBOSE){cout<<"Done! ("<<Sz11<<" points)"<<endl;}
	
	if(VERBOSE){cout<<"Calculating Asymptotic Constants...";}
	alglib::real_2d_array A;
	alglib::real_1d_array C1,B1;
	// Define System Matrix containing coefficients of linear system
	A.setlength(10,10);
	for(int i=0;i<10;i++)
		{for(int j=0;j<10;j++)
			{A[i][j]=0;}}
	// d(IEP1)/dr
	A[0][0]=y1[Sz1-1][8];A[0][1]=y2[Sz2-1][8];A[0][2]=y3[Sz3-1][8];A[0][3]=y4[Sz4-1][8];A[0][4]=y5[Sz5-1][8];A[0][5]=y6[Sz6-1][8];A[0][6]=y7[Sz7-1][8];A[0][7]=y8[Sz8-1][8];
	A[0][8]=y9[Sz9-1][8];A[0][9]=y10[Sz10-1][8];
	// d(IEP2)/dr
	A[1][0]=y1[Sz1-1][10];A[1][1]=y2[Sz2-1][10];A[1][2]=y3[Sz3-1][10];A[1][3]=y4[Sz4-1][10];A[1][4]=y5[Sz5-1][10];A[1][5]=y6[Sz6-1][10];A[1][6]=y7[Sz7-1][10];A[1][7]=y8[Sz8-1][10];
	A[1][8]=y9[Sz9-1][10];A[1][9]=y10[Sz10-1][10];
	// d(IEP3)/dr
	A[2][0]=y1[Sz1-1][12];A[2][1]=y2[Sz2-1][12];A[2][2]=y3[Sz3-1][12];A[2][3]=y4[Sz4-1][12];A[2][4]=y5[Sz5-1][12];A[2][5]=y6[Sz6-1][12];A[2][6]=y7[Sz7-1][12];A[2][7]=y8[Sz8-1][12];
	A[2][8]=y9[Sz9-1][12];A[2][9]=y10[Sz10-1][12];
	// d(IEP4)/dr
	A[3][0]=y1[Sz1-1][14];A[3][1]=y2[Sz2-1][14];A[3][2]=y3[Sz3-1][14];A[3][3]=y4[Sz4-1][14];A[3][4]=y5[Sz5-1][14];A[3][5]=y6[Sz6-1][14];A[3][6]=y7[Sz7-1][14];A[3][7]=y8[Sz8-1][14];
	A[3][8]=y9[Sz9-1][14];A[3][9]=y10[Sz10-1][14];
	// d(IEP5)/dr
	A[4][0]=y1[Sz1-1][16];A[4][1]=y2[Sz2-1][16];A[4][2]=y3[Sz3-1][16];A[4][3]=y4[Sz4-1][16];A[4][4]=y5[Sz5-1][16];A[4][5]=y6[Sz6-1][16];A[4][6]=y7[Sz7-1][16];A[4][7]=y8[Sz8-1][16];
	A[4][8]=y9[Sz9-1][16];A[4][9]=y10[Sz10-1][16];
	// d(IEP6)/dr
	A[5][0]=y1[Sz1-1][18];A[5][1]=y2[Sz2-1][18];A[5][2]=y3[Sz3-1][18];A[5][3]=y4[Sz4-1][18];A[5][4]=y5[Sz5-1][18];A[5][5]=y6[Sz6-1][18];A[5][6]=y7[Sz7-1][18];A[5][7]=y8[Sz8-1][18];
	A[5][8]=y9[Sz9-1][18];A[5][9]=y10[Sz10-1][18];
	// d(IEP7)/dr
	A[6][0]=y1[Sz1-1][20];A[6][1]=y2[Sz2-1][20];A[6][2]=y3[Sz3-1][20];A[6][3]=y4[Sz4-1][20];A[6][4]=y5[Sz5-1][20];A[6][5]=y6[Sz6-1][20];A[6][6]=y7[Sz7-1][20];A[6][7]=y8[Sz8-1][20];
	A[6][8]=y9[Sz9-1][20];A[6][9]=y10[Sz10-1][20];
	// d(IEP8)/dr
	A[7][0]=y1[Sz1-1][22];
	A[7][1]=y2[Sz2-1][22];
	A[7][2]=y3[Sz3-1][22];
	A[7][3]=y4[Sz4-1][22];
	A[7][4]=y5[Sz5-1][22];
	A[7][5]=y6[Sz6-1][22];
	A[7][6]=y7[Sz7-1][22];
	A[7][7]=y8[Sz8-1][22];
	A[7][8]=y9[Sz9-1][22];
	A[7][9]=y10[Sz10-1][22];
	// f 1st (df/dr)
	A[8][0]=y1[Sz1-1][1];
	A[8][1]=y2[Sz2-1][1];
	A[8][2]=y3[Sz3-1][1];
	A[8][3]=y4[Sz4-1][1];
	A[8][4]=y5[Sz5-1][1];
	A[8][5]=y6[Sz6-1][1];
	A[8][6]=y7[Sz7-1][1];
	A[8][7]=y8[Sz8-1][1];
	A[8][8]=y9[Sz9-1][1];
	A[8][9]=y10[Sz10-1][1];
	// f 2nd (d2f/dr2)
	A[9][0]=y1[Sz1-1][2];
	A[9][1]=y2[Sz2-1][2];
	A[9][2]=y3[Sz3-1][2];
	A[9][3]=y4[Sz4-1][2];
	A[9][4]=y5[Sz5-1][2];
	A[9][5]=y6[Sz6-1][2];
	A[9][6]=y7[Sz7-1][2];
	A[9][7]=y8[Sz8-1][2];
	A[9][8]=y9[Sz9-1][2];
	A[9][9]=y10[Sz10-1][2];
	// Define Solution Matrix for Problem #1
	B1.setlength(10);
	// d(IEP1)/dr = 0
	B1[0]=0-y11[Sz11-1][8];
	// d(IEP2)/dr = 0
	B1[1]=0-y11[Sz11-1][10];
	// d(IEP3)/dr = 0
	B1[2]=0-y11[Sz11-1][12];
	// d(IEP4)/dr = 0
	B1[3]=0-y11[Sz11-1][14];
	// d(IEP5)/dr = 0
	B1[4]=0-y11[Sz11-1][16];
	// d(IEP6)/dr = 0
	B1[5]=0-y11[Sz11-1][18];
	// d(IEP7)/dr = 0
	B1[6]=0-y11[Sz11-1][20];
	// d(IEP8)/dr = 0
	B1[7]=0-y11[Sz11-1][22];
	// df/dr = -Radius/2
	B1[8]=-a/2 - y11[Sz11-1][1];
	// d2f/dr2 = -0.5
	B1[9]=-0.5 - y11[Sz11-1][2];
	// Create Output Parameters
	alglib::ae_int_t info;
	alglib::densesolverreport rep;
	alglib::rmatrixsolve(A,10,B1,info,rep,C1);
	if(VERBOSE)
	{cout<<"\n";
	cout<<C1[0]*A[0][0] + C1[1]*A[0][1] + C1[2]*A[0][2] + C1[3]*A[0][3] + C1[4]*A[0][4] + C1[5]*A[0][5] + C1[6]*A[0][6] + C1[7]*A[0][7] + C1[8]*A[0][8] + C1[9]*A[0][9] + y11[Sz11-1][8]<<endl;
	cout<<C1[0]*A[1][0] + C1[1]*A[1][1] + C1[2]*A[1][2] + C1[3]*A[1][3] + C1[4]*A[1][4] + C1[5]*A[1][5] + C1[6]*A[1][6] + C1[7]*A[1][7] + C1[8]*A[1][8] + C1[9]*A[1][9] + y11[Sz11-1][10]<<endl;
	cout<<C1[0]*A[2][0] + C1[1]*A[2][1] + C1[2]*A[2][2] + C1[3]*A[2][3] + C1[4]*A[2][4] + C1[5]*A[2][5] + C1[6]*A[2][6] + C1[7]*A[2][7] + C1[8]*A[2][8] + C1[9]*A[2][9] + y11[Sz11-1][12]<<endl;
	cout<<C1[0]*A[3][0] + C1[1]*A[3][1] + C1[2]*A[3][2] + C1[3]*A[3][3] + C1[4]*A[3][4] + C1[5]*A[3][5] + C1[6]*A[3][6] + C1[7]*A[3][7] + C1[8]*A[3][8] + C1[9]*A[3][9] + y11[Sz11-1][14]<<endl;
	cout<<C1[0]*A[4][0] + C1[1]*A[4][1] + C1[2]*A[4][2] + C1[3]*A[4][3] + C1[4]*A[4][4] + C1[5]*A[4][5] + C1[6]*A[4][6] + C1[7]*A[4][7] + C1[8]*A[4][8] + C1[9]*A[4][9] + y11[Sz11-1][16]<<endl;
	cout<<C1[0]*A[5][0] + C1[1]*A[5][1] + C1[2]*A[5][2] + C1[3]*A[5][3] + C1[4]*A[5][4] + C1[5]*A[5][5] + C1[6]*A[5][6] + C1[7]*A[5][7] + C1[8]*A[5][8] + C1[9]*A[5][9] + y11[Sz11-1][18]<<endl;
	cout<<C1[0]*A[6][0] + C1[1]*A[6][1] + C1[2]*A[6][2] + C1[3]*A[6][3] + C1[4]*A[6][4] + C1[5]*A[6][5] + C1[6]*A[6][6] + C1[7]*A[6][7] + C1[8]*A[6][8] + C1[9]*A[6][9] + y11[Sz11-1][20]<<endl;
	cout<<C1[0]*A[7][0] + C1[1]*A[7][1] + C1[2]*A[7][2] + C1[3]*A[7][3] + C1[4]*A[7][4] + C1[5]*A[7][5] + C1[6]*A[7][6] + C1[7]*A[7][7] + C1[8]*A[7][8] + C1[9]*A[7][9] + y11[Sz11-1][22]<<endl;
	cout<<C1[0]*A[8][0] + C1[1]*A[8][1] + C1[2]*A[8][2] + C1[3]*A[8][3] + C1[4]*A[8][4] + C1[5]*A[8][5] + C1[6]*A[8][6] + C1[7]*A[8][7] + C1[8]*A[8][8] + C1[9]*A[8][9] + y11[Sz11-1][1]<<endl;
	cout<<C1[0]*A[9][0] + C1[1]*A[9][1] + C1[2]*A[9][2] + C1[3]*A[9][3] + C1[4]*A[9][4] + C1[5]*A[9][5] + C1[6]*A[9][6] + C1[7]*A[9][7] + C1[8]*A[9][8] + C1[9]*A[9][9] + y11[Sz11-1][2]<<endl;
	cout<<"Done!"<<endl;
	cout<<"C1_1 = "<<C1[0]<<endl;
	cout<<"C2_1 = "<<C1[1]<<endl;
	cout<<"C3_1 = "<<C1[2]<<endl;
	cout<<"C4_1 = "<<C1[3]<<endl;
	cout<<"C5_1 = "<<C1[4]<<endl;
	cout<<"C6_1 = "<<C1[5]<<endl;
	cout<<"C7_1 = "<<C1[6]<<endl;
	cout<<"C8_1 = "<<C1[7]<<endl;
	cout<<"C9_1 = "<<C1[8]<<endl;
	cout<<"C10_1 = "<<C1[9]<<endl;}

	if(VERBOSE){cout<<"Solving Electric Field Problem...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Eigth Ion
	y[21]=0;y[22]=0;
	// Constants
	y[23]=S.temperature;y[24]=S.numIon;
	// Fist Ion Properties
	y[25]=S.ionVal[0];y[26]=S.ionConc[0]*1000;y[27]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[28]=S.ionVal[1];y[29]=S.ionConc[1]*1000;y[30]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[31]=S.ionVal[2];y[32]=S.ionConc[2]*1000;y[33]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[34]=S.ionVal[3];y[35]=S.ionConc[3]*1000;y[36]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[37]=S.ionVal[4];y[38]=S.ionConc[4]*1000;y[39]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[40]=S.ionVal[5];y[41]=S.ionConc[5]*1000;y[42]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[43]=S.ionVal[6];y[44]=S.ionConc[6]*1000;y[45]=S.ionRad[6]*1e-10;
	// Eigth Ion Properties
	y[46]=S.ionVal[7];y[47]=S.ionConc[7]*1000;y[48]=S.ionRad[7]*1e-10;
	// Solution Properties
	y[49]=S.dielectric;y[50]=Debye;y[51]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps12=integrate_adaptive(controlled_stepper,inhomogeneous_Prob2_8ion,y,r0,rf,dr,push_back_state_and_time(y12,r12));
	int Sz12=steps12+1;
	if(VERBOSE){cout<<"Done! ("<<Sz12<<" points)"<<endl;}

	if(VERBOSE){cout<<"Calculating Asymptotic Constants...";}
	alglib::real_1d_array C2,B2;
	// d(IEP1)/dr
	A[0][0]=y1[Sz1-1][8];A[0][1]=y2[Sz2-1][8];A[0][2]=y3[Sz3-1][8];A[0][3]=y4[Sz4-1][8];A[0][4]=y5[Sz5-1][8];A[0][5]=y6[Sz6-1][8];A[0][6]=y7[Sz7-1][8];A[0][7]=y8[Sz8-1][8];
	A[0][8]=y9[Sz9-1][8];A[0][9]=y10[Sz10-1][8];
	// d(IEP2)/dr
	A[1][0]=y1[Sz1-1][10];A[1][1]=y2[Sz2-1][10];A[1][2]=y3[Sz3-1][10];A[1][3]=y4[Sz4-1][10];A[1][4]=y5[Sz5-1][10];A[1][5]=y6[Sz6-1][10];A[1][6]=y7[Sz7-1][10];A[1][7]=y8[Sz8-1][10];
	A[1][8]=y9[Sz9-1][10];A[1][9]=y10[Sz10-1][10];
	// d(IEP3)/dr
	A[2][0]=y1[Sz1-1][12];A[2][1]=y2[Sz2-1][12];A[2][2]=y3[Sz3-1][12];A[2][3]=y4[Sz4-1][12];A[2][4]=y5[Sz5-1][12];A[2][5]=y6[Sz6-1][12];A[2][6]=y7[Sz7-1][12];A[2][7]=y8[Sz8-1][12];
	A[2][8]=y9[Sz9-1][12];A[2][9]=y10[Sz10-1][12];
	// d(IEP4)/dr
	A[3][0]=y1[Sz1-1][14];A[3][1]=y2[Sz2-1][14];A[3][2]=y3[Sz3-1][14];A[3][3]=y4[Sz4-1][14];A[3][4]=y5[Sz5-1][14];A[3][5]=y6[Sz6-1][14];A[3][6]=y7[Sz7-1][14];A[3][7]=y8[Sz8-1][14];
	A[3][8]=y9[Sz9-1][14];A[3][9]=y10[Sz10-1][14];
	// d(IEP5)/dr
	A[4][0]=y1[Sz1-1][16];A[4][1]=y2[Sz2-1][16];A[4][2]=y3[Sz3-1][16];A[4][3]=y4[Sz4-1][16];A[4][4]=y5[Sz5-1][16];A[4][5]=y6[Sz6-1][16];A[4][6]=y7[Sz7-1][16];A[4][7]=y8[Sz8-1][16];
	A[4][8]=y9[Sz9-1][16];A[4][9]=y10[Sz10-1][16];
	// d(IEP6)/dr
	A[5][0]=y1[Sz1-1][18];A[5][1]=y2[Sz2-1][18];A[5][2]=y3[Sz3-1][18];A[5][3]=y4[Sz4-1][18];A[5][4]=y5[Sz5-1][18];A[5][5]=y6[Sz6-1][18];A[5][6]=y7[Sz7-1][18];A[5][7]=y8[Sz8-1][18];
	A[5][8]=y9[Sz9-1][18];A[5][9]=y10[Sz10-1][18];
	// d(IEP7)/dr
	A[6][0]=y1[Sz1-1][20];
	A[6][1]=y2[Sz2-1][20];
	A[6][2]=y3[Sz3-1][20];
	A[6][3]=y4[Sz4-1][20];
	A[6][4]=y5[Sz5-1][20];
	A[6][5]=y6[Sz6-1][20];
	A[6][6]=y7[Sz7-1][20];
	A[6][7]=y8[Sz8-1][20];
	A[6][8]=y9[Sz9-1][20];
	A[6][9]=y10[Sz10-1][20];
	// d(IEP8)/dr
	A[7][0]=y1[Sz1-1][22];
	A[7][1]=y2[Sz2-1][22];
	A[7][2]=y3[Sz3-1][22];
	A[7][3]=y4[Sz4-1][22];
	A[7][4]=y5[Sz5-1][22];
	A[7][5]=y6[Sz6-1][22];
	A[7][6]=y7[Sz7-1][22];
	A[7][7]=y8[Sz8-1][22];
	A[7][8]=y9[Sz9-1][22];
	A[7][9]=y10[Sz10-1][22];
	// f 1st (df/dr)
	A[8][0]=y1[Sz1-1][1];
	A[8][1]=y2[Sz2-1][1];
	A[8][2]=y3[Sz3-1][1];
	A[8][3]=y4[Sz4-1][1];
	A[8][4]=y5[Sz5-1][1];
	A[8][5]=y6[Sz6-1][1];
	A[8][6]=y7[Sz7-1][1];
	A[8][7]=y8[Sz8-1][1];
	A[8][8]=y9[Sz9-1][1];
	A[8][9]=y10[Sz10-1][1];
	// f 2nd (d2f/dr2)
	A[9][0]=y1[Sz1-1][2];
	A[9][1]=y2[Sz2-1][2];
	A[9][2]=y3[Sz3-1][2];
	A[9][3]=y4[Sz4-1][2];
	A[9][4]=y5[Sz5-1][2];
	A[9][5]=y6[Sz6-1][2];
	A[9][6]=y7[Sz7-1][2];
	A[9][7]=y8[Sz8-1][2];
	A[9][8]=y9[Sz9-1][2];
	A[9][9]=y10[Sz10-1][2];
	// Define Solution Matrix for Problem #2
	B2.setlength(10);
	// d(IEP1)/dr = -1
	B2[0]=-1-y12[Sz12-1][8];
	// d(IEP2)/dr = -1
	B2[1]=-1-y12[Sz12-1][10];
	// d(IEP3)/dr = -1
	B2[2]=-1-y12[Sz12-1][12];
	// d(IEP4)/dr = -1
	B2[3]=-1-y12[Sz12-1][14];
	// d(IEP5)/dr = -1
	B2[4]=-1-y12[Sz12-1][16];
	// d(IEP6)/dr = -1
	B2[5]=-1-y12[Sz12-1][18];
	// d(IEP7)/dr = -1
	B2[6]=-1-y12[Sz12-1][20];
	// d(IEP8)/dr = -1
	B2[7]=-1-y12[Sz12-1][22];
	// df/dr = 0
	B2[8]=0-y12[Sz12-1][1];
	// d2f/dr2 = 0
	B2[9]=0-y12[Sz12-1][2];
	// Create Output Parameters
	alglib::rmatrixsolve(A,10,B2,info,rep,C2);
	if(VERBOSE)
	{cout<<"\n";
	cout<<C2[0]*A[0][0] + C2[1]*A[0][1] + C2[2]*A[0][2] + C2[3]*A[0][3] + C2[4]*A[0][4] + C2[5]*A[0][5] + C2[6]*A[0][6] + C2[7]*A[0][7] + C2[8]*A[0][8] + C2[9]*A[0][9] + y12[Sz12-1][8]<<endl;
	cout<<C2[0]*A[1][0] + C2[1]*A[1][1] + C2[2]*A[1][2] + C2[3]*A[1][3] + C2[4]*A[1][4] + C2[5]*A[1][5] + C2[6]*A[1][6] + C2[7]*A[1][7] + C2[8]*A[1][8] + C2[9]*A[1][9] + y12[Sz12-1][10]<<endl;
	cout<<C2[0]*A[2][0] + C2[1]*A[2][1] + C2[2]*A[2][2] + C2[3]*A[2][3] + C2[4]*A[2][4] + C2[5]*A[2][5] + C2[6]*A[2][6] + C2[7]*A[2][7] + C2[8]*A[2][8] + C2[9]*A[2][9] + y12[Sz12-1][12]<<endl;
	cout<<C2[0]*A[3][0] + C2[1]*A[3][1] + C2[2]*A[3][2] + C2[3]*A[3][3] + C2[4]*A[3][4] + C2[5]*A[3][5] + C2[6]*A[3][6] + C2[7]*A[3][7] + C2[8]*A[3][8] + C2[9]*A[3][9] + y12[Sz12-1][14]<<endl;
	cout<<C2[0]*A[4][0] + C2[1]*A[4][1] + C2[2]*A[4][2] + C2[3]*A[4][3] + C2[4]*A[4][4] + C2[5]*A[4][5] + C2[6]*A[4][6] + C2[7]*A[4][7] + C2[8]*A[4][8] + C2[9]*A[4][9] + y12[Sz12-1][16]<<endl;
	cout<<C2[0]*A[5][0] + C2[1]*A[5][1] + C2[2]*A[5][2] + C2[3]*A[5][3] + C2[4]*A[5][4] + C2[5]*A[5][5] + C2[6]*A[5][6] + C2[7]*A[5][7] + C2[8]*A[5][8] + C2[9]*A[5][9] + y12[Sz12-1][18]<<endl;
	cout<<C2[0]*A[6][0] + C2[1]*A[6][1] + C2[2]*A[6][2] + C2[3]*A[6][3] + C2[4]*A[6][4] + C2[5]*A[6][5] + C2[6]*A[6][6] + C2[7]*A[6][7] + C2[8]*A[6][8] + C2[9]*A[6][9] + y12[Sz12-1][20]<<endl;
	cout<<C2[0]*A[7][0] + C2[1]*A[7][1] + C2[2]*A[7][2] + C2[3]*A[7][3] + C2[4]*A[7][4] + C2[5]*A[7][5] + C2[6]*A[7][6] + C2[7]*A[7][7] + C2[8]*A[7][8] + C2[9]*A[7][9] + y12[Sz12-1][22]<<endl;
	cout<<C2[0]*A[8][0] + C2[1]*A[8][1] + C2[2]*A[8][2] + C2[3]*A[8][3] + C2[4]*A[8][4] + C2[5]*A[8][5] + C2[6]*A[8][6] + C2[7]*A[8][7] + C2[8]*A[8][8] + C2[9]*A[8][9] + y12[Sz12-1][1]<<endl;
	cout<<C2[0]*A[9][0] + C2[1]*A[9][1] + C2[2]*A[9][2] + C2[3]*A[9][3] + C2[4]*A[9][4] + C2[5]*A[9][5] + C2[6]*A[9][6] + C2[7]*A[9][7] + C2[8]*A[9][8] + C2[9]*A[9][9] + y12[Sz12-1][2]<<endl;
	cout<<"Done!"<<endl;
	cout<<"C1_2 = "<<C2[0]<<endl;
	cout<<"C2_2 = "<<C2[1]<<endl;
	cout<<"C3_2 = "<<C2[2]<<endl;
	cout<<"C4_2 = "<<C2[3]<<endl;
	cout<<"C5_2 = "<<C2[4]<<endl;
	cout<<"C6_2 = "<<C2[5]<<endl;
	cout<<"C7_2 = "<<C2[6]<<endl;
	cout<<"C8_2 = "<<C2[7]<<endl;
	cout<<"C9_2 = "<<C2[8]<<endl;
	cout<<"C10_2 = "<<C2[9]<<endl;}

	double mobility=-C2[8]/C1[8];
	
	if(VERBOSE)
		{cout<<"\n\nElectrophoretic Mobility:\t"<<mobility<<" m^2/Vs\n";
		cout<<"ka:\t\t"<<a/Debye<<endl;
		cout<<"Zeta Potential:\t"<<ZP<<" V\n";
		cout<<"E:\t\t"<<3*S.viscosity*ec*mobility/(2*eo*S.dielectric*kb*S.temperature)<<endl;
		cout<<"y:\t\t"<<ec*ZP/(kb*S.temperature)<<endl;}

	return mobility;}

double nine_ion_problem(double ZP,double a,Solution S,Error E,bool VERBOSE)
	{// Convert Debye Length from Angstroms to Meters
	double Debye=S.debyeLength*1e-10;
	// Final Position: Solvated Radius [meters]
	double rf=a;
	double r0;
	double dr=-1e-15;	// integrate from Inf to Protein Radius
	if(Debye>a){r0=5*Debye;}
	else{r0=5*a;}
		
	if(VERBOSE){cout<<"9 ION PROBLEM\nInitial Position:\t"<<r0<<" m"<<"\nFinal Position:\t\t"<<rf<<" m"<<"\nStep Size:\t\t"<<dr<<" m"<<endl;}

	// Calculate Poisson Boltzmann Coefficient and Determine Initial Position away from Protein
	double C=calc_PBE_coefficient(ZP,r0,rf,dr,S,E,VERBOSE);
	if(VERBOSE){cout<<"PBE Coeff: "<<C<<endl;}
	if(isnan(C)){return nan("");}
	int N=7+2*S.numIon+2+3*S.numIon+3;//57;
	state_type y(N);
	// F Function
	// y[0] = f
	// y[1] = df/dr
	// y[2] = d2f/dr2
	// y[3] = d3f/dr3
	// y[4] = d4f/dr4
	// Equilibrium Electric Pontential
	// y[5] = potential
	// y[6] = d(potential)/dr
	// Induced Electrochemical Potential Function (9 ion)
	// y[7] = IEP1
	// y[8] = d(IEP1)/dr
	// y[9] = IEP2
	// y[10] = d(IEP2)/dr
	// y[11] = IEP3
	// y[12] = d(IEP3)/dr
	// y[13] = IEP4
	// y[14] = d(IEP4)/dr
	// y[15] = IEP5
	// y[16] = d(IEP5)/dr
	// y[17] = IEP6
	// y[18] = d(IEP6)/dr
	// y[19] = IEP7
	// y[20] = d(IEP7)/dr
	// y[21] = IEP8
	// y[22] = d(IEP8)/dr
	// y[23] = IEP9
	// y[24] = d(IEP9)/dr
	// Constants
	// y[25] = Temp					// [Kelvin]
	// y[26] = numIon
	// Ion Properties
	// y[27] = Ion1.val;		
	// y[28] = Ion1.conc;		// mM [mol/m^3]
	// y[29] = Ion1.radius;		// meters
	// y[30] = Ion2.val;
	// y[31] = Ion2.conc;
	// y[32] = Ion2.radius;
	// y[33] = Ion3.val;
	// y[34] = Ion3.conc;
	// y[35] = Ion3.radius;
	// y[36] = Ion4.val;
	// y[37] = Ion4.conc;
	// y[38] = Ion4.radius;
	// y[39] = Ion5.val;
	// y[40] = Ion5.conc;
	// y[41] = Ion5.radius;
	// y[42] = Ion6.val;
	// y[43] = Ion6.conc;
	// y[44] = Ion6.radius;
	// y[45] = Ion7.val;
	// y[46] = Ion7.conc;
	// y[47] = Ion7.radius;
	// y[48] = Ion8.val;
	// y[49] = Ion8.conc;
	// y[50] = Ion8.radius;
	// y[51] = Ion9.val;
	// y[52] = Ion9.conc;
	// y[53] = Ion9.radius;
	// Solution Properties
	// y[54] = relative dielectric	[dimensionless]
	// y[55] = Debye Length [meters]
	// y[56] = viscosity  [Pa s]
	// Collect Error Data
	double abs_err=E.absolute,rel_err=E.relative,a_x=E.a_x,a_dxdt=E.a_dxdt;
	controlled_stepper_type controlled_stepper(default_error_checker<double,range_algebra,default_operations>(abs_err,rel_err,a_x,a_dxdt));	

	// Calculate Mobility
	vector<state_type> y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13;
	vector<double> r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13;
	// Specify Boundary Conditions / Initial Values Infinitely Far from the Particle
	if(VERBOSE){cout<<"Solving Homogeneous Problem for 1st Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=1/(r0*r0);y[8]=-2/(r0*r0*r0);
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Eigth Ion
	y[21]=0;y[22]=0;
	// Ninth Ion
	y[23]=0;y[24]=0;
	// Constants
	y[25]=S.temperature;y[26]=S.numIon;
	// Fist Ion Properties
	y[27]=S.ionVal[0];y[28]=S.ionConc[0]*1000;y[29]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[30]=S.ionVal[1];y[31]=S.ionConc[1]*1000;y[32]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[33]=S.ionVal[2];y[34]=S.ionConc[2]*1000;y[35]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[36]=S.ionVal[3];y[37]=S.ionConc[3]*1000;y[38]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[39]=S.ionVal[4];y[40]=S.ionConc[4]*1000;y[41]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[42]=S.ionVal[5];y[43]=S.ionConc[5]*1000;y[44]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[45]=S.ionVal[6];y[46]=S.ionConc[6]*1000;y[47]=S.ionRad[6]*1e-10;
	// Eigth Ion Properties
	y[48]=S.ionVal[7];y[49]=S.ionConc[7]*1000;y[50]=S.ionRad[7]*1e-10;
	// Ninth Ion Properties
	y[51]=S.ionVal[8];y[52]=S.ionConc[8]*1000;y[53]=S.ionRad[8]*1e-10;
	// Solution Properties
	y[54]=S.dielectric;y[55]=Debye;y[56]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps1=integrate_adaptive(controlled_stepper,homogeneous_form_9ion,y,r0,rf,dr,push_back_state_and_time(y1,r1));
	int Sz1=steps1+1;
	if(VERBOSE){cout<<"Done! ("<<Sz1<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 2nd Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=1/(r0*r0);y[10]=-2/(r0*r0*r0);
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Eigth Ion
	y[21]=0;y[22]=0;
	// Ninth Ion
	y[23]=0;y[24]=0;
	// Constants
	y[25]=S.temperature;y[26]=S.numIon;
	// Fist Ion Properties
	y[27]=S.ionVal[0];y[28]=S.ionConc[0]*1000;y[29]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[30]=S.ionVal[1];y[31]=S.ionConc[1]*1000;y[32]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[33]=S.ionVal[2];y[34]=S.ionConc[2]*1000;y[35]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[36]=S.ionVal[3];y[37]=S.ionConc[3]*1000;y[38]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[39]=S.ionVal[4];y[40]=S.ionConc[4]*1000;y[41]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[42]=S.ionVal[5];y[43]=S.ionConc[5]*1000;y[44]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[45]=S.ionVal[6];y[46]=S.ionConc[6]*1000;y[47]=S.ionRad[6]*1e-10;
	// Eigth Ion Properties
	y[48]=S.ionVal[7];y[49]=S.ionConc[7]*1000;y[50]=S.ionRad[7]*1e-10;
	// Ninth Ion Properties
	y[51]=S.ionVal[8];y[52]=S.ionConc[8]*1000;y[53]=S.ionRad[8]*1e-10;
	// Solution Properties
	y[54]=S.dielectric;y[55]=Debye;y[56]=S.viscosity;
	// Solve Homogeneous Form for First Ion
	size_t steps2=integrate_adaptive(controlled_stepper,homogeneous_form_9ion,y,r0,rf,dr,push_back_state_and_time(y2,r2));
	int Sz2=steps2+1;
	if(VERBOSE){cout<<"Done! ("<<Sz2<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 3rd Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=1/(r0*r0);y[12]=-2/(r0*r0*r0);
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Eigth Ion
	y[21]=0;y[22]=0;
	// Ninth Ion
	y[23]=0;y[24]=0;
	// Constants
	y[25]=S.temperature;y[26]=S.numIon;
	// Fist Ion Properties
	y[27]=S.ionVal[0];y[28]=S.ionConc[0]*1000;y[29]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[30]=S.ionVal[1];y[31]=S.ionConc[1]*1000;y[32]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[33]=S.ionVal[2];y[34]=S.ionConc[2]*1000;y[35]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[36]=S.ionVal[3];y[37]=S.ionConc[3]*1000;y[38]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[39]=S.ionVal[4];y[40]=S.ionConc[4]*1000;y[41]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[42]=S.ionVal[5];y[43]=S.ionConc[5]*1000;y[44]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[45]=S.ionVal[6];y[46]=S.ionConc[6]*1000;y[47]=S.ionRad[6]*1e-10;
	// Eigth Ion Properties
	y[48]=S.ionVal[7];y[49]=S.ionConc[7]*1000;y[50]=S.ionRad[7]*1e-10;
	// Ninth Ion Properties
	y[51]=S.ionVal[8];y[52]=S.ionConc[8]*1000;y[53]=S.ionRad[8]*1e-10;
	// Solution Properties
	y[54]=S.dielectric;y[55]=Debye;y[56]=S.viscosity;
	// Solve Homogeneous Form
	size_t steps3=integrate_adaptive(controlled_stepper,homogeneous_form_9ion,y,r0,rf,dr,push_back_state_and_time(y3,r3));
	int Sz3=steps3+1;
	if(VERBOSE){cout<<"Done! ("<<Sz3<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 4th Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=1/(r0*r0);y[14]=-2/(r0*r0*r0);
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Eigth Ion
	y[21]=0;y[22]=0;
	// Ninth Ion
	y[23]=0;y[24]=0;
	// Constants
	y[25]=S.temperature;y[26]=S.numIon;
	// Fist Ion Properties
	y[27]=S.ionVal[0];y[28]=S.ionConc[0]*1000;y[29]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[30]=S.ionVal[1];y[31]=S.ionConc[1]*1000;y[32]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[33]=S.ionVal[2];y[34]=S.ionConc[2]*1000;y[35]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[36]=S.ionVal[3];y[37]=S.ionConc[3]*1000;y[38]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[39]=S.ionVal[4];y[40]=S.ionConc[4]*1000;y[41]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[42]=S.ionVal[5];y[43]=S.ionConc[5]*1000;y[44]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[45]=S.ionVal[6];y[46]=S.ionConc[6]*1000;y[47]=S.ionRad[6]*1e-10;
	// Eigth Ion Properties
	y[48]=S.ionVal[7];y[49]=S.ionConc[7]*1000;y[50]=S.ionRad[7]*1e-10;
	// Ninth Ion Properties
	y[51]=S.ionVal[8];y[52]=S.ionConc[8]*1000;y[53]=S.ionRad[8]*1e-10;
	// Solution Properties
	y[54]=S.dielectric;y[55]=Debye;y[56]=S.viscosity;
	// Solve Homogeneous Form
	size_t steps4=integrate_adaptive(controlled_stepper,homogeneous_form_9ion,y,r0,rf,dr,push_back_state_and_time(y4,r4));
	int Sz4=steps4+1;
	if(VERBOSE){cout<<"Done! ("<<Sz4<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 5th Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=1/(r0*r0);y[16]=-2/(r0*r0*r0);
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Eigth Ion
	y[21]=0;y[22]=0;
	// Ninth Ion
	y[23]=0;y[24]=0;
	// Constants
	y[25]=S.temperature;y[26]=S.numIon;
	// Fist Ion Properties
	y[27]=S.ionVal[0];y[28]=S.ionConc[0]*1000;y[29]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[30]=S.ionVal[1];y[31]=S.ionConc[1]*1000;y[32]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[33]=S.ionVal[2];y[34]=S.ionConc[2]*1000;y[35]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[36]=S.ionVal[3];y[37]=S.ionConc[3]*1000;y[38]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[39]=S.ionVal[4];y[40]=S.ionConc[4]*1000;y[41]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[42]=S.ionVal[5];y[43]=S.ionConc[5]*1000;y[44]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[45]=S.ionVal[6];y[46]=S.ionConc[6]*1000;y[47]=S.ionRad[6]*1e-10;
	// Eigth Ion Properties
	y[48]=S.ionVal[7];y[49]=S.ionConc[7]*1000;y[50]=S.ionRad[7]*1e-10;
	// Ninth Ion Properties
	y[51]=S.ionVal[8];y[52]=S.ionConc[8]*1000;y[53]=S.ionRad[8]*1e-10;
	// Solution Properties
	y[54]=S.dielectric;y[55]=Debye;y[56]=S.viscosity;
	// Solve Homogeneous Form
	size_t steps5=integrate_adaptive(controlled_stepper,homogeneous_form_9ion,y,r0,rf,dr,push_back_state_and_time(y5,r5));
	int Sz5=steps5+1;
	if(VERBOSE){cout<<"Done! ("<<Sz5<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 6th Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=1/(r0*r0);y[18]=-2/(r0*r0*r0);
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Eigth Ion
	y[21]=0;y[22]=0;
	// Ninth Ion
	y[23]=0;y[24]=0;
	// Constants
	y[25]=S.temperature;y[26]=S.numIon;
	// Fist Ion Properties
	y[27]=S.ionVal[0];y[28]=S.ionConc[0]*1000;y[29]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[30]=S.ionVal[1];y[31]=S.ionConc[1]*1000;y[32]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[33]=S.ionVal[2];y[34]=S.ionConc[2]*1000;y[35]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[36]=S.ionVal[3];y[37]=S.ionConc[3]*1000;y[38]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[39]=S.ionVal[4];y[40]=S.ionConc[4]*1000;y[41]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[42]=S.ionVal[5];y[43]=S.ionConc[5]*1000;y[44]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[45]=S.ionVal[6];y[46]=S.ionConc[6]*1000;y[47]=S.ionRad[6]*1e-10;
	// Eigth Ion Properties
	y[48]=S.ionVal[7];y[49]=S.ionConc[7]*1000;y[50]=S.ionRad[7]*1e-10;
	// Ninth Ion Properties
	y[51]=S.ionVal[8];y[52]=S.ionConc[8]*1000;y[53]=S.ionRad[8]*1e-10;
	// Solution Properties
	y[54]=S.dielectric;y[55]=Debye;y[56]=S.viscosity;
	// Solve Homogeneous Form
	size_t steps6=integrate_adaptive(controlled_stepper,homogeneous_form_9ion,y,r0,rf,dr,push_back_state_and_time(y6,r6));
	int Sz6=steps6+1;
	if(VERBOSE){cout<<"Done! ("<<Sz6<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 7th Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=1/(r0*r0);y[20]=-2/(r0*r0*r0);
	// Eigth Ion
	y[21]=0;y[22]=0;
	// Ninth Ion
	y[23]=0;y[24]=0;
	// Constants
	y[25]=S.temperature;y[26]=S.numIon;
	// Fist Ion Properties
	y[27]=S.ionVal[0];y[28]=S.ionConc[0]*1000;y[29]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[30]=S.ionVal[1];y[31]=S.ionConc[1]*1000;y[32]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[33]=S.ionVal[2];y[34]=S.ionConc[2]*1000;y[35]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[36]=S.ionVal[3];y[37]=S.ionConc[3]*1000;y[38]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[39]=S.ionVal[4];y[40]=S.ionConc[4]*1000;y[41]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[42]=S.ionVal[5];y[43]=S.ionConc[5]*1000;y[44]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[45]=S.ionVal[6];y[46]=S.ionConc[6]*1000;y[47]=S.ionRad[6]*1e-10;
	// Eigth Ion Properties
	y[48]=S.ionVal[7];y[49]=S.ionConc[7]*1000;y[50]=S.ionRad[7]*1e-10;
	// Ninth Ion Properties
	y[51]=S.ionVal[8];y[52]=S.ionConc[8]*1000;y[53]=S.ionRad[8]*1e-10;
	// Solution Properties
	y[54]=S.dielectric;y[55]=Debye;y[56]=S.viscosity;
	// Solve Homogeneous Form
	size_t steps7=integrate_adaptive(controlled_stepper,homogeneous_form_9ion,y,r0,rf,dr,push_back_state_and_time(y7,r7));
	int Sz7=steps7+1;
	if(VERBOSE){cout<<"Done! ("<<Sz7<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 8th Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Eigth Ion
	y[21]=1/(r0*r0);y[22]=-2/(r0*r0*r0);
	// Ninth Ion
	y[23]=0;y[24]=0;
	// Constants
	y[25]=S.temperature;y[26]=S.numIon;
	// Fist Ion Properties
	y[27]=S.ionVal[0];y[28]=S.ionConc[0]*1000;y[29]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[30]=S.ionVal[1];y[31]=S.ionConc[1]*1000;y[32]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[33]=S.ionVal[2];y[34]=S.ionConc[2]*1000;y[35]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[36]=S.ionVal[3];y[37]=S.ionConc[3]*1000;y[38]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[39]=S.ionVal[4];y[40]=S.ionConc[4]*1000;y[41]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[42]=S.ionVal[5];y[43]=S.ionConc[5]*1000;y[44]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[45]=S.ionVal[6];y[46]=S.ionConc[6]*1000;y[47]=S.ionRad[6]*1e-10;
	// Eigth Ion Properties
	y[48]=S.ionVal[7];y[49]=S.ionConc[7]*1000;y[50]=S.ionRad[7]*1e-10;
	// Ninth Ion Properties
	y[51]=S.ionVal[8];y[52]=S.ionConc[8]*1000;y[53]=S.ionRad[8]*1e-10;
	// Solution Properties
	y[54]=S.dielectric;y[55]=Debye;y[56]=S.viscosity;
	// Solve Homogeneous Form
	size_t steps8=integrate_adaptive(controlled_stepper,homogeneous_form_9ion,y,r0,rf,dr,push_back_state_and_time(y8,r8));
	int Sz8=steps8+1;
	if(VERBOSE){cout<<"Done! ("<<Sz8<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 9th Ion...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Eigth Ion
	y[21]=0;y[22]=0;
	// Ninth Ion
	y[23]=1/(r0*r0);y[24]=-2/(r0*r0*r0);
	// Constants
	y[25]=S.temperature;y[26]=S.numIon;
	// Fist Ion Properties
	y[27]=S.ionVal[0];y[28]=S.ionConc[0]*1000;y[29]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[30]=S.ionVal[1];y[31]=S.ionConc[1]*1000;y[32]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[33]=S.ionVal[2];y[34]=S.ionConc[2]*1000;y[35]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[36]=S.ionVal[3];y[37]=S.ionConc[3]*1000;y[38]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[39]=S.ionVal[4];y[40]=S.ionConc[4]*1000;y[41]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[42]=S.ionVal[5];y[43]=S.ionConc[5]*1000;y[44]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[45]=S.ionVal[6];y[46]=S.ionConc[6]*1000;y[47]=S.ionRad[6]*1e-10;
	// Eigth Ion Properties
	y[48]=S.ionVal[7];y[49]=S.ionConc[7]*1000;y[50]=S.ionRad[7]*1e-10;
	// Ninth Ion Properties
	y[51]=S.ionVal[8];y[52]=S.ionConc[8]*1000;y[53]=S.ionRad[8]*1e-10;
	// Solution Properties
	y[54]=S.dielectric;y[55]=Debye;y[56]=S.viscosity;
	// Solve Homogeneous Form
	size_t steps9=integrate_adaptive(controlled_stepper,homogeneous_form_9ion,y,r0,rf,dr,push_back_state_and_time(y9,r9));
	int Sz9=steps9+1;
	if(VERBOSE){cout<<"Done! ("<<Sz9<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 1st f term...";}
	// F Function
	y[0]=r0;y[1]=1;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Eigth Ion
	y[21]=0;y[22]=0;
	// Ninth Ion
	y[23]=0;y[24]=0;
	// Constants
	y[25]=S.temperature;y[26]=S.numIon;
	// Fist Ion Properties
	y[27]=S.ionVal[0];y[28]=S.ionConc[0]*1000;y[29]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[30]=S.ionVal[1];y[31]=S.ionConc[1]*1000;y[32]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[33]=S.ionVal[2];y[34]=S.ionConc[2]*1000;y[35]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[36]=S.ionVal[3];y[37]=S.ionConc[3]*1000;y[38]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[39]=S.ionVal[4];y[40]=S.ionConc[4]*1000;y[41]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[42]=S.ionVal[5];y[43]=S.ionConc[5]*1000;y[44]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[45]=S.ionVal[6];y[46]=S.ionConc[6]*1000;y[47]=S.ionRad[6]*1e-10;
	// Eigth Ion Properties
	y[48]=S.ionVal[7];y[49]=S.ionConc[7]*1000;y[50]=S.ionRad[7]*1e-10;
	// Ninth Ion Properties
	y[51]=S.ionVal[8];y[52]=S.ionConc[8]*1000;y[53]=S.ionRad[8]*1e-10;
	// Solution Properties
	y[54]=S.dielectric;y[55]=Debye;y[56]=S.viscosity;
	// Solve Homogeneous Form
	size_t steps10=integrate_adaptive(controlled_stepper,homogeneous_form_9ion,y,r0,rf,dr,push_back_state_and_time(y10,r10));
	int Sz10=steps10+1;
	if(VERBOSE){cout<<"Done! ("<<Sz10<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Homogeneous Problem for 2nd f term...";}
	// F Function
	y[0]=1/r0;y[1]=-1/(r0*r0);y[2]=2/(r0*r0*r0);y[3]=-6/(r0*r0*r0*r0);y[4]=24/(r0*r0*r0*r0*r0);
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Eigth Ion
	y[21]=0;y[22]=0;
	// Ninth Ion
	y[23]=0;y[24]=0;
	// Constants
	y[25]=S.temperature;y[26]=S.numIon;
	// Fist Ion Properties
	y[27]=S.ionVal[0];y[28]=S.ionConc[0]*1000;y[29]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[30]=S.ionVal[1];y[31]=S.ionConc[1]*1000;y[32]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[33]=S.ionVal[2];y[34]=S.ionConc[2]*1000;y[35]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[36]=S.ionVal[3];y[37]=S.ionConc[3]*1000;y[38]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[39]=S.ionVal[4];y[40]=S.ionConc[4]*1000;y[41]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[42]=S.ionVal[5];y[43]=S.ionConc[5]*1000;y[44]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[45]=S.ionVal[6];y[46]=S.ionConc[6]*1000;y[47]=S.ionRad[6]*1e-10;
	// Eigth Ion Properties
	y[48]=S.ionVal[7];y[49]=S.ionConc[7]*1000;y[50]=S.ionRad[7]*1e-10;
	// Ninth Ion Properties
	y[51]=S.ionVal[8];y[52]=S.ionConc[8]*1000;y[53]=S.ionRad[8]*1e-10;
	// Solution Properties
	y[54]=S.dielectric;y[55]=Debye;y[56]=S.viscosity;
	// Solve Homogeneous Form
	size_t steps11=integrate_adaptive(controlled_stepper,homogeneous_form_9ion,y,r0,rf,dr,push_back_state_and_time(y11,r11));
	int Sz11=steps11+1;
	if(VERBOSE){cout<<"Done! ("<<Sz11<<" points)"<<endl;}

	if(VERBOSE){cout<<"Solving Flow Field Problem...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Eigth Ion
	y[21]=0;y[22]=0;
	// Ninth Ion
	y[23]=0;y[24]=0;
	// Constants
	y[25]=S.temperature;y[26]=S.numIon;
	// Fist Ion Properties
	y[27]=S.ionVal[0];y[28]=S.ionConc[0]*1000;y[29]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[30]=S.ionVal[1];y[31]=S.ionConc[1]*1000;y[32]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[33]=S.ionVal[2];y[34]=S.ionConc[2]*1000;y[35]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[36]=S.ionVal[3];y[37]=S.ionConc[3]*1000;y[38]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[39]=S.ionVal[4];y[40]=S.ionConc[4]*1000;y[41]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[42]=S.ionVal[5];y[43]=S.ionConc[5]*1000;y[44]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[45]=S.ionVal[6];y[46]=S.ionConc[6]*1000;y[47]=S.ionRad[6]*1e-10;
	// Eigth Ion Properties
	y[48]=S.ionVal[7];y[49]=S.ionConc[7]*1000;y[50]=S.ionRad[7]*1e-10;
	// Ninth Ion Properties
	y[51]=S.ionVal[8];y[52]=S.ionConc[8]*1000;y[53]=S.ionRad[8]*1e-10;
	// Solution Properties
	y[54]=S.dielectric;y[55]=Debye;y[56]=S.viscosity;
	// Solve InHomogeneous Form
	size_t steps12=integrate_adaptive(controlled_stepper,inhomogeneous_Prob1_9ion,y,r0,rf,dr,push_back_state_and_time(y12,r12));
	int Sz12=steps12+1;
	if(VERBOSE){cout<<"Done! ("<<Sz12<<" points)"<<endl;}
	
	if(VERBOSE){cout<<"Calculating Asymptotic Constants...";}
	alglib::real_2d_array A;
	alglib::real_1d_array C1,B1;
	// Define System Matrix containing coefficients of linear system
	A.setlength(11,11);
	for(int i=0;i<11;i++)
		{for(int j=0;j<11;j++)
			{A[i][j]=0;}}
	// d(IEP1)/dr
	A[0][0]=y1[Sz1-1][8];A[0][1]=y2[Sz2-1][8];A[0][2]=y3[Sz3-1][8];A[0][3]=y4[Sz4-1][8];A[0][4]=y5[Sz5-1][8];A[0][5]=y6[Sz6-1][8];A[0][6]=y7[Sz7-1][8];A[0][7]=y8[Sz8-1][8];
	A[0][8]=y9[Sz9-1][8];A[0][9]=y10[Sz10-1][8];A[0][10]=y11[Sz11-1][8];
	// d(IEP2)/dr
	A[1][0]=y1[Sz1-1][10];A[1][1]=y2[Sz2-1][10];A[1][2]=y3[Sz3-1][10];A[1][3]=y4[Sz4-1][10];A[1][4]=y5[Sz5-1][10];A[1][5]=y6[Sz6-1][10];A[1][6]=y7[Sz7-1][10];A[1][7]=y8[Sz8-1][10];
	A[1][8]=y9[Sz9-1][10];A[1][9]=y10[Sz10-1][10];A[1][10]=y11[Sz11-1][10];
	// d(IEP3)/dr
	A[2][0]=y1[Sz1-1][12];A[2][1]=y2[Sz2-1][12];A[2][2]=y3[Sz3-1][12];A[2][3]=y4[Sz4-1][12];A[2][4]=y5[Sz5-1][12];A[2][5]=y6[Sz6-1][12];A[2][6]=y7[Sz7-1][12];A[2][7]=y8[Sz8-1][12];
	A[2][8]=y9[Sz9-1][12];A[2][9]=y10[Sz10-1][12];A[2][10]=y11[Sz11-1][12];
	// d(IEP4)/dr
	A[3][0]=y1[Sz1-1][14];A[3][1]=y2[Sz2-1][14];A[3][2]=y3[Sz3-1][14];A[3][3]=y4[Sz4-1][14];A[3][4]=y5[Sz5-1][14];A[3][5]=y6[Sz6-1][14];A[3][6]=y7[Sz7-1][14];A[3][7]=y8[Sz8-1][14];
	A[3][8]=y9[Sz9-1][14];A[3][9]=y10[Sz10-1][14];A[3][10]=y11[Sz11-1][14];
	// d(IEP5)/dr
	A[4][0]=y1[Sz1-1][16];A[4][1]=y2[Sz2-1][16];A[4][2]=y3[Sz3-1][16];A[4][3]=y4[Sz4-1][16];A[4][4]=y5[Sz5-1][16];A[4][5]=y6[Sz6-1][16];A[4][6]=y7[Sz7-1][16];A[4][7]=y8[Sz8-1][16];
	A[4][8]=y9[Sz9-1][16];A[4][9]=y10[Sz10-1][16];A[4][10]=y11[Sz11-1][16];
	// d(IEP6)/dr
	A[5][0]=y1[Sz1-1][18];A[5][1]=y2[Sz2-1][18];A[5][2]=y3[Sz3-1][18];A[5][3]=y4[Sz4-1][18];A[5][4]=y5[Sz5-1][18];A[5][5]=y6[Sz6-1][18];A[5][6]=y7[Sz7-1][18];A[5][7]=y8[Sz8-1][18];
	A[5][8]=y9[Sz9-1][18];A[5][9]=y10[Sz10-1][18];A[5][10]=y11[Sz11-1][18];
	// d(IEP7)/dr
	A[6][0]=y1[Sz1-1][20];A[6][1]=y2[Sz2-1][20];A[6][2]=y3[Sz3-1][20];A[6][3]=y4[Sz4-1][20];A[6][4]=y5[Sz5-1][20];A[6][5]=y6[Sz6-1][20];A[6][6]=y7[Sz7-1][20];A[6][7]=y8[Sz8-1][20];
	A[6][8]=y9[Sz9-1][20];A[6][9]=y10[Sz10-1][20];A[6][10]=y11[Sz11-1][20];
	// d(IEP8)/dr
	A[7][0]=y1[Sz1-1][22];
	A[7][1]=y2[Sz2-1][22];
	A[7][2]=y3[Sz3-1][22];
	A[7][3]=y4[Sz4-1][22];
	A[7][4]=y5[Sz5-1][22];
	A[7][5]=y6[Sz6-1][22];
	A[7][6]=y7[Sz7-1][22];
	A[7][7]=y8[Sz8-1][22];
	A[7][8]=y9[Sz9-1][22];
	A[7][9]=y10[Sz10-1][22];
	A[7][10]=y11[Sz11-1][22];
	// d(IEP9)/dr
	A[8][0]=y1[Sz1-1][24];
	A[8][1]=y2[Sz2-1][24];
	A[8][2]=y3[Sz3-1][24];
	A[8][3]=y4[Sz4-1][24];
	A[8][4]=y5[Sz5-1][24];
	A[8][5]=y6[Sz6-1][24];
	A[8][6]=y7[Sz7-1][24];
	A[8][7]=y8[Sz8-1][24];
	A[8][8]=y9[Sz9-1][24];
	A[8][9]=y10[Sz10-1][24];
	A[8][10]=y11[Sz11-1][24];
	// f 1st (df/dr)
	A[9][0]=y1[Sz1-1][1];
	A[9][1]=y2[Sz2-1][1];
	A[9][2]=y3[Sz3-1][1];
	A[9][3]=y4[Sz4-1][1];
	A[9][4]=y5[Sz5-1][1];
	A[9][5]=y6[Sz6-1][1];
	A[9][6]=y7[Sz7-1][1];
	A[9][7]=y8[Sz8-1][1];
	A[9][8]=y9[Sz9-1][1];
	A[9][9]=y10[Sz10-1][1];
	A[9][10]=y11[Sz11-1][1];
	// f 2nd (d2f/dr2)
	A[10][0]=y1[Sz1-1][2];
	A[10][1]=y2[Sz2-1][2];
	A[10][2]=y3[Sz3-1][2];
	A[10][3]=y4[Sz4-1][2];
	A[10][4]=y5[Sz5-1][2];
	A[10][5]=y6[Sz6-1][2];
	A[10][6]=y7[Sz7-1][2];
	A[10][7]=y8[Sz8-1][2];
	A[10][8]=y9[Sz9-1][2];
	A[10][9]=y10[Sz10-1][2];
	A[10][10]=y11[Sz11-1][2];
	// Define Solution Matrix for Problem #1
	B1.setlength(11);
	// d(IEP1)/dr = 0
	B1[0]=0-y12[Sz12-1][8];
	// d(IEP2)/dr = 0
	B1[1]=0-y12[Sz12-1][10];
	// d(IEP3)/dr = 0
	B1[2]=0-y12[Sz12-1][12];
	// d(IEP4)/dr = 0
	B1[3]=0-y12[Sz12-1][14];
	// d(IEP5)/dr = 0
	B1[4]=0-y12[Sz12-1][16];
	// d(IEP6)/dr = 0
	B1[5]=0-y12[Sz12-1][18];
	// d(IEP7)/dr = 0
	B1[6]=0-y12[Sz12-1][20];
	// d(IEP8)/dr = 0
	B1[7]=0-y12[Sz12-1][22];
	// d(IEP9)/dr = 0
	B1[8]=0-y12[Sz12-1][24];
	// df/dr = -Radius/2
	B1[9]=-a/2 - y12[Sz12-1][1];
	// d2f/dr2 = -0.5
	B1[10]=-0.5 - y12[Sz12-1][2];
	// Create Output Parameters
	alglib::ae_int_t info;
	alglib::densesolverreport rep;
	alglib::rmatrixsolve(A,11,B1,info,rep,C1);
	if(VERBOSE)
	{cout<<"\n";
/*
	cout<<C1[0]*A[0][0] + C1[1]*A[0][1] + C1[2]*A[0][2] + C1[3]*A[0][3] + C1[4]*A[0][4] + C1[5]*A[0][5] + C1[6]*A[0][6] + C1[7]*A[0][7] + C1[8]*A[0][8] + C1[9]*A[0][9] + y11[Sz11-1][8]<<endl;
	cout<<C1[0]*A[1][0] + C1[1]*A[1][1] + C1[2]*A[1][2] + C1[3]*A[1][3] + C1[4]*A[1][4] + C1[5]*A[1][5] + C1[6]*A[1][6] + C1[7]*A[1][7] + C1[8]*A[1][8] + C1[9]*A[1][9] + y11[Sz11-1][10]<<endl;
	cout<<C1[0]*A[2][0] + C1[1]*A[2][1] + C1[2]*A[2][2] + C1[3]*A[2][3] + C1[4]*A[2][4] + C1[5]*A[2][5] + C1[6]*A[2][6] + C1[7]*A[2][7] + C1[8]*A[2][8] + C1[9]*A[2][9] + y11[Sz11-1][12]<<endl;
	cout<<C1[0]*A[3][0] + C1[1]*A[3][1] + C1[2]*A[3][2] + C1[3]*A[3][3] + C1[4]*A[3][4] + C1[5]*A[3][5] + C1[6]*A[3][6] + C1[7]*A[3][7] + C1[8]*A[3][8] + C1[9]*A[3][9] + y11[Sz11-1][14]<<endl;
	cout<<C1[0]*A[4][0] + C1[1]*A[4][1] + C1[2]*A[4][2] + C1[3]*A[4][3] + C1[4]*A[4][4] + C1[5]*A[4][5] + C1[6]*A[4][6] + C1[7]*A[4][7] + C1[8]*A[4][8] + C1[9]*A[4][9] + y11[Sz11-1][16]<<endl;
	cout<<C1[0]*A[5][0] + C1[1]*A[5][1] + C1[2]*A[5][2] + C1[3]*A[5][3] + C1[4]*A[5][4] + C1[5]*A[5][5] + C1[6]*A[5][6] + C1[7]*A[5][7] + C1[8]*A[5][8] + C1[9]*A[5][9] + y11[Sz11-1][18]<<endl;
	cout<<C1[0]*A[6][0] + C1[1]*A[6][1] + C1[2]*A[6][2] + C1[3]*A[6][3] + C1[4]*A[6][4] + C1[5]*A[6][5] + C1[6]*A[6][6] + C1[7]*A[6][7] + C1[8]*A[6][8] + C1[9]*A[6][9] + y11[Sz11-1][20]<<endl;
	cout<<C1[0]*A[7][0] + C1[1]*A[7][1] + C1[2]*A[7][2] + C1[3]*A[7][3] + C1[4]*A[7][4] + C1[5]*A[7][5] + C1[6]*A[7][6] + C1[7]*A[7][7] + C1[8]*A[7][8] + C1[9]*A[7][9] + y11[Sz11-1][22]<<endl;
	cout<<C1[0]*A[8][0] + C1[1]*A[8][1] + C1[2]*A[8][2] + C1[3]*A[8][3] + C1[4]*A[8][4] + C1[5]*A[8][5] + C1[6]*A[8][6] + C1[7]*A[8][7] + C1[8]*A[8][8] + C1[9]*A[8][9] + y11[Sz11-1][1]<<endl;
	cout<<C1[0]*A[9][0] + C1[1]*A[9][1] + C1[2]*A[9][2] + C1[3]*A[9][3] + C1[4]*A[9][4] + C1[5]*A[9][5] + C1[6]*A[9][6] + C1[7]*A[9][7] + C1[8]*A[9][8] + C1[9]*A[9][9] + y11[Sz11-1][2]<<endl;
*/
	cout<<"Done!"<<endl;
	cout<<"C1_1 = "<<C1[0]<<endl;
	cout<<"C2_1 = "<<C1[1]<<endl;
	cout<<"C3_1 = "<<C1[2]<<endl;
	cout<<"C4_1 = "<<C1[3]<<endl;
	cout<<"C5_1 = "<<C1[4]<<endl;
	cout<<"C6_1 = "<<C1[5]<<endl;
	cout<<"C7_1 = "<<C1[6]<<endl;
	cout<<"C8_1 = "<<C1[7]<<endl;
	cout<<"C9_1 = "<<C1[8]<<endl;
	cout<<"C10_1 = "<<C1[9]<<endl;}

	if(VERBOSE){cout<<"Solving Electric Field Problem...";}
	// F Function
	y[0]=0;y[1]=0;y[2]=0;y[3]=0;y[4]=0;
	// Poisson-Boltzmann Equation for Equilibrium Potential
	y[5]=C*exp(-r0/Debye)/r0;y[6]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// First Ion
	y[7]=0;y[8]=0;
	// Second Ion
	y[9]=0;y[10]=0;
	// Third Ion
	y[11]=0;y[12]=0;
	// Fourth Ion
	y[13]=0;y[14]=0;
	// Fifth Ion
	y[15]=0;y[16]=0;
	// Sixth Ion
	y[17]=0;y[18]=0;
	// Seventh Ion
	y[19]=0;y[20]=0;
	// Eigth Ion
	y[21]=0;y[22]=0;
	// Ninth Ion
	y[23]=0;y[24]=0;
	// Constants
	y[25]=S.temperature;y[26]=S.numIon;
	// Fist Ion Properties
	y[27]=S.ionVal[0];y[28]=S.ionConc[0]*1000;y[29]=S.ionRad[0]*1e-10;
	// Second Ion Properties
	y[30]=S.ionVal[1];y[31]=S.ionConc[1]*1000;y[32]=S.ionRad[1]*1e-10;
	// Third Ion Properties
	y[33]=S.ionVal[2];y[34]=S.ionConc[2]*1000;y[35]=S.ionRad[2]*1e-10;
	// Fourth Ion Properties
	y[36]=S.ionVal[3];y[37]=S.ionConc[3]*1000;y[38]=S.ionRad[3]*1e-10;
	// Fifth Ion Properties
	y[39]=S.ionVal[4];y[40]=S.ionConc[4]*1000;y[41]=S.ionRad[4]*1e-10;
	// Sixth Ion Properties
	y[42]=S.ionVal[5];y[43]=S.ionConc[5]*1000;y[44]=S.ionRad[5]*1e-10;
	// Seventh Ion Properties
	y[45]=S.ionVal[6];y[46]=S.ionConc[6]*1000;y[47]=S.ionRad[6]*1e-10;
	// Eigth Ion Properties
	y[48]=S.ionVal[7];y[49]=S.ionConc[7]*1000;y[50]=S.ionRad[7]*1e-10;
	// Ninth Ion Properties
	y[51]=S.ionVal[8];y[52]=S.ionConc[8]*1000;y[53]=S.ionRad[8]*1e-10;
	// Solution Properties
	y[54]=S.dielectric;y[55]=Debye;y[56]=S.viscosity;
	// Solve InHomogeneous Form
	size_t steps13=integrate_adaptive(controlled_stepper,inhomogeneous_Prob2_9ion,y,r0,rf,dr,push_back_state_and_time(y13,r13));
	int Sz13=steps13+1;
	if(VERBOSE){cout<<"Done! ("<<Sz13<<" points)"<<endl;}

	if(VERBOSE){cout<<"Calculating Asymptotic Constants...";}
	alglib::real_1d_array C2,B2;
	// d(IEP1)/dr
	A[0][0]=y1[Sz1-1][8];A[0][1]=y2[Sz2-1][8];A[0][2]=y3[Sz3-1][8];A[0][3]=y4[Sz4-1][8];A[0][4]=y5[Sz5-1][8];A[0][5]=y6[Sz6-1][8];A[0][6]=y7[Sz7-1][8];A[0][7]=y8[Sz8-1][8];
	A[0][8]=y9[Sz9-1][8];A[0][9]=y10[Sz10-1][8];A[0][10]=y11[Sz11-1][8];
	// d(IEP2)/dr
	A[1][0]=y1[Sz1-1][10];A[1][1]=y2[Sz2-1][10];A[1][2]=y3[Sz3-1][10];A[1][3]=y4[Sz4-1][10];A[1][4]=y5[Sz5-1][10];A[1][5]=y6[Sz6-1][10];A[1][6]=y7[Sz7-1][10];A[1][7]=y8[Sz8-1][10];
	A[1][8]=y9[Sz9-1][10];A[1][9]=y10[Sz10-1][10];A[1][10]=y11[Sz11-1][10];
	// d(IEP3)/dr
	A[2][0]=y1[Sz1-1][12];A[2][1]=y2[Sz2-1][12];A[2][2]=y3[Sz3-1][12];A[2][3]=y4[Sz4-1][12];A[2][4]=y5[Sz5-1][12];A[2][5]=y6[Sz6-1][12];A[2][6]=y7[Sz7-1][12];A[2][7]=y8[Sz8-1][12];
	A[2][8]=y9[Sz9-1][12];A[2][9]=y10[Sz10-1][12];A[2][10]=y11[Sz11-1][12];
	// d(IEP4)/dr
	A[3][0]=y1[Sz1-1][14];A[3][1]=y2[Sz2-1][14];A[3][2]=y3[Sz3-1][14];A[3][3]=y4[Sz4-1][14];A[3][4]=y5[Sz5-1][14];A[3][5]=y6[Sz6-1][14];A[3][6]=y7[Sz7-1][14];A[3][7]=y8[Sz8-1][14];
	A[3][8]=y9[Sz9-1][14];A[3][9]=y10[Sz10-1][14];A[3][10]=y11[Sz11-1][14];
	// d(IEP5)/dr
	A[4][0]=y1[Sz1-1][16];A[4][1]=y2[Sz2-1][16];A[4][2]=y3[Sz3-1][16];A[4][3]=y4[Sz4-1][16];A[4][4]=y5[Sz5-1][16];A[4][5]=y6[Sz6-1][16];A[4][6]=y7[Sz7-1][16];A[4][7]=y8[Sz8-1][16];
	A[4][8]=y9[Sz9-1][16];A[4][9]=y10[Sz10-1][16];A[4][10]=y11[Sz11-1][16];
	// d(IEP6)/dr
	A[5][0]=y1[Sz1-1][18];A[5][1]=y2[Sz2-1][18];A[5][2]=y3[Sz3-1][18];A[5][3]=y4[Sz4-1][18];A[5][4]=y5[Sz5-1][18];A[5][5]=y6[Sz6-1][18];A[5][6]=y7[Sz7-1][18];A[5][7]=y8[Sz8-1][18];
	A[5][8]=y9[Sz9-1][18];A[5][9]=y10[Sz10-1][18];A[5][10]=y11[Sz11-1][18];
	// d(IEP7)/dr
	A[6][0]=y1[Sz1-1][20];
	A[6][1]=y2[Sz2-1][20];
	A[6][2]=y3[Sz3-1][20];
	A[6][3]=y4[Sz4-1][20];
	A[6][4]=y5[Sz5-1][20];
	A[6][5]=y6[Sz6-1][20];
	A[6][6]=y7[Sz7-1][20];
	A[6][7]=y8[Sz8-1][20];
	A[6][8]=y9[Sz9-1][20];
	A[6][9]=y10[Sz10-1][20];
	A[6][10]=y11[Sz11-1][20];
	// d(IEP8)/dr
	A[7][0]=y1[Sz1-1][22];
	A[7][1]=y2[Sz2-1][22];
	A[7][2]=y3[Sz3-1][22];
	A[7][3]=y4[Sz4-1][22];
	A[7][4]=y5[Sz5-1][22];
	A[7][5]=y6[Sz6-1][22];
	A[7][6]=y7[Sz7-1][22];
	A[7][7]=y8[Sz8-1][22];
	A[7][8]=y9[Sz9-1][22];
	A[7][9]=y10[Sz10-1][22];
	A[7][10]=y11[Sz11-1][22];
	// d(IEP9)/dr
	A[8][0]=y1[Sz1-1][24];
	A[8][1]=y2[Sz2-1][24];
	A[8][2]=y3[Sz3-1][24];
	A[8][3]=y4[Sz4-1][24];
	A[8][4]=y5[Sz5-1][24];
	A[8][5]=y6[Sz6-1][24];
	A[8][6]=y7[Sz7-1][24];
	A[8][7]=y8[Sz8-1][24];
	A[8][8]=y9[Sz9-1][24];
	A[8][9]=y10[Sz10-1][24];
	A[8][10]=y11[Sz11-1][24];
	// f 1st (df/dr)
	A[9][0]=y1[Sz1-1][1];
	A[9][1]=y2[Sz2-1][1];
	A[9][2]=y3[Sz3-1][1];
	A[9][3]=y4[Sz4-1][1];
	A[9][4]=y5[Sz5-1][1];
	A[9][5]=y6[Sz6-1][1];
	A[9][6]=y7[Sz7-1][1];
	A[9][7]=y8[Sz8-1][1];
	A[9][8]=y9[Sz9-1][1];
	A[9][9]=y10[Sz10-1][1];
	A[9][10]=y11[Sz11-1][1];
	// f 2nd (d2f/dr2)
	A[10][0]=y1[Sz1-1][2];
	A[10][1]=y2[Sz2-1][2];
	A[10][2]=y3[Sz3-1][2];
	A[10][3]=y4[Sz4-1][2];
	A[10][4]=y5[Sz5-1][2];
	A[10][5]=y6[Sz6-1][2];
	A[10][6]=y7[Sz7-1][2];
	A[10][7]=y8[Sz8-1][2];
	A[10][8]=y9[Sz9-1][2];
	A[10][9]=y10[Sz10-1][2];
	A[10][10]=y11[Sz11-1][2];
	// Define Solution Matrix for Problem #2
	B2.setlength(11);
	// d(IEP1)/dr = -1
	B2[0]=-1-y13[Sz13-1][8];
	// d(IEP2)/dr = -1
	B2[1]=-1-y13[Sz13-1][10];
	// d(IEP3)/dr = -1
	B2[2]=-1-y13[Sz13-1][12];
	// d(IEP4)/dr = -1
	B2[3]=-1-y13[Sz13-1][14];
	// d(IEP5)/dr = -1
	B2[4]=-1-y13[Sz13-1][16];
	// d(IEP6)/dr = -1
	B2[5]=-1-y13[Sz13-1][18];
	// d(IEP7)/dr = -1
	B2[6]=-1-y13[Sz13-1][20];
	// d(IEP8)/dr = -1
	B2[7]=-1-y13[Sz13-1][22];
	// d(IEP9)/dr = -1
	B2[8]=-1-y13[Sz13-1][24];
	// df/dr = 0
	B2[9]=0-y13[Sz13-1][1];
	// d2f/dr2 = 0
	B2[10]=0-y13[Sz13-1][2];
	// Create Output Parameters
	alglib::rmatrixsolve(A,11,B2,info,rep,C2);
	if(VERBOSE)
	{cout<<"\n";
/*	cout<<C2[0]*A[0][0] + C2[1]*A[0][1] + C2[2]*A[0][2] + C2[3]*A[0][3] + C2[4]*A[0][4] + C2[5]*A[0][5] + C2[6]*A[0][6] + C2[7]*A[0][7] + C2[8]*A[0][8] + C2[9]*A[0][9] + y12[Sz12-1][8]<<endl;
	cout<<C2[0]*A[1][0] + C2[1]*A[1][1] + C2[2]*A[1][2] + C2[3]*A[1][3] + C2[4]*A[1][4] + C2[5]*A[1][5] + C2[6]*A[1][6] + C2[7]*A[1][7] + C2[8]*A[1][8] + C2[9]*A[1][9] + y12[Sz12-1][10]<<endl;
	cout<<C2[0]*A[2][0] + C2[1]*A[2][1] + C2[2]*A[2][2] + C2[3]*A[2][3] + C2[4]*A[2][4] + C2[5]*A[2][5] + C2[6]*A[2][6] + C2[7]*A[2][7] + C2[8]*A[2][8] + C2[9]*A[2][9] + y12[Sz12-1][12]<<endl;
	cout<<C2[0]*A[3][0] + C2[1]*A[3][1] + C2[2]*A[3][2] + C2[3]*A[3][3] + C2[4]*A[3][4] + C2[5]*A[3][5] + C2[6]*A[3][6] + C2[7]*A[3][7] + C2[8]*A[3][8] + C2[9]*A[3][9] + y12[Sz12-1][14]<<endl;
	cout<<C2[0]*A[4][0] + C2[1]*A[4][1] + C2[2]*A[4][2] + C2[3]*A[4][3] + C2[4]*A[4][4] + C2[5]*A[4][5] + C2[6]*A[4][6] + C2[7]*A[4][7] + C2[8]*A[4][8] + C2[9]*A[4][9] + y12[Sz12-1][16]<<endl;
	cout<<C2[0]*A[5][0] + C2[1]*A[5][1] + C2[2]*A[5][2] + C2[3]*A[5][3] + C2[4]*A[5][4] + C2[5]*A[5][5] + C2[6]*A[5][6] + C2[7]*A[5][7] + C2[8]*A[5][8] + C2[9]*A[5][9] + y12[Sz12-1][18]<<endl;
	cout<<C2[0]*A[6][0] + C2[1]*A[6][1] + C2[2]*A[6][2] + C2[3]*A[6][3] + C2[4]*A[6][4] + C2[5]*A[6][5] + C2[6]*A[6][6] + C2[7]*A[6][7] + C2[8]*A[6][8] + C2[9]*A[6][9] + y12[Sz12-1][20]<<endl;
	cout<<C2[0]*A[7][0] + C2[1]*A[7][1] + C2[2]*A[7][2] + C2[3]*A[7][3] + C2[4]*A[7][4] + C2[5]*A[7][5] + C2[6]*A[7][6] + C2[7]*A[7][7] + C2[8]*A[7][8] + C2[9]*A[7][9] + y12[Sz12-1][22]<<endl;
	cout<<C2[0]*A[8][0] + C2[1]*A[8][1] + C2[2]*A[8][2] + C2[3]*A[8][3] + C2[4]*A[8][4] + C2[5]*A[8][5] + C2[6]*A[8][6] + C2[7]*A[8][7] + C2[8]*A[8][8] + C2[9]*A[8][9] + y12[Sz12-1][1]<<endl;
	cout<<C2[0]*A[9][0] + C2[1]*A[9][1] + C2[2]*A[9][2] + C2[3]*A[9][3] + C2[4]*A[9][4] + C2[5]*A[9][5] + C2[6]*A[9][6] + C2[7]*A[9][7] + C2[8]*A[9][8] + C2[9]*A[9][9] + y12[Sz12-1][2]<<endl;
*/
	cout<<"Done!"<<endl;
	cout<<"C1_2 = "<<C2[0]<<endl;
	cout<<"C2_2 = "<<C2[1]<<endl;
	cout<<"C3_2 = "<<C2[2]<<endl;
	cout<<"C4_2 = "<<C2[3]<<endl;
	cout<<"C5_2 = "<<C2[4]<<endl;
	cout<<"C6_2 = "<<C2[5]<<endl;
	cout<<"C7_2 = "<<C2[6]<<endl;
	cout<<"C8_2 = "<<C2[7]<<endl;
	cout<<"C9_2 = "<<C2[8]<<endl;
	cout<<"C10_2 = "<<C2[9]<<endl;}

	double mobility=-C2[9]/C1[9];
	
	if(VERBOSE)
		{cout<<"\n\nElectrophoretic Mobility:\t"<<mobility<<" m^2/Vs\n";
		cout<<"ka:\t\t"<<a/Debye<<endl;
		cout<<"Zeta Potential:\t"<<ZP<<" V\n";
		cout<<"E:\t\t"<<3*S.viscosity*ec*mobility/(2*eo*S.dielectric*kb*S.temperature)<<endl;
		cout<<"y:\t\t"<<ec*ZP/(kb*S.temperature)<<endl;}

	return mobility;}

// INPUT: Zeta Potential (Zeta) [Volts]
// INPUT: Initial Position (r0) [meters]
// INPUT: Final Position (rf) [meters]
// INPUT: Integration Interval (dr) [meters]
// INPUT: Debye Length (k) [meters]
// INPUT: Temperature (Temp) [Kelvin]
// INPUT: Ion Concentrations (ion_conc) [Molar]
double calc_PBE_coefficient(double Zeta,double &r0,double rf,double &dr,Solution S,Error E,bool VERBOSE)
	{// Convert Debye Length to meters
	double k=S.debyeLength*1e-10;

	double abs_err=E.absolute;
	double rel_err=E.relative;
	double a_x=E.a_x;
	double a_dxdt=E.a_dxdt;
	
	controlled_stepper_type controlled_stepper(default_error_checker<double,range_algebra,default_operations>(abs_err,rel_err,a_x,a_dxdt));
	
	int N=5+2*S.numIon;	
	state_type y(N);
	// Potential [V]
	// y[0]=C*exp(-r0/Debye)/r0;
	// d(Potential)/dr
	// y[1]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// Constants
	// y[2]=S.temperature;		// Temperature [Kevlin]
	// y[3]=S.numIon;
	// y[4]=S.dielectric;
	// Ion Data
	//for(int i=0;i<S.numIon;i++)
	//	{// Ion Valence
	//	y[5+2*i]=S.ionVal[i];
	// Ion Concentration [milli-Molar] (or equalivalently [mol/meter^3])
	//	y[6+2*i]=S.ionConc[i]*1000;}}

	vector<state_type> y1,y2,y3;	// dependent variables (defined above)
	vector<double> r1,r2,r3;
	double C1,reductionFactor=1e12;
	bool RUNNING=true;
	size_t steps;
	int Sz1,Counter=0;
	while(RUNNING)
		{C1=(Zeta*r0/exp(-r0/k))/reductionFactor;
		if(VERBOSE){cout<<"Solving Poisson Boltzman Equation...\nFirst Guess:\tC = "<<C1<<"...";}
		initialize_PBE(y,r0,C1,S);
		steps=integrate_adaptive(controlled_stepper,poisson_boltzmann,y,r0,rf,dr,push_back_state_and_time(y1,r1));
		Sz1=steps+1;
		if(!isnan(y1[Sz1-1][0]))
			{// Good Choice of C1
			RUNNING=false;break;}
		else
			{Counter++;
			if(Counter>=MAX_PBE_ITERATIONS){return nan("");}
			// NaN Error
			if(VERBOSE){cout<<"Output NaN! Retrying...\n";}
			// Increase Reduction Factor
			reductionFactor*=10;
			y1.clear(); r1.clear();}
		}
	double Err1=Zeta-y1[Sz1-1][0];
	if(VERBOSE){cout<<"Done! ("<<Sz1<<" points)\t Calculated: "<<y1[Sz1-1][0]<<endl;}
	y1.clear(); r1.clear();

	//double C2=Zeta*r0/exp(-r0/k)/1e8;
	double C2=Zeta*r0/exp(-r0/k)/(reductionFactor*1e4);
	if(VERBOSE){cout<<"Second Guess:\tC = "<<C2<<"...";}
	initialize_PBE(y,r0,C2,S);
	steps=integrate_adaptive(controlled_stepper,poisson_boltzmann,y,r0,rf,dr,push_back_state_and_time(y2,r2));
	int Sz2=steps+1;
	double Err2=Zeta-y2[Sz2-1][0];
	if(VERBOSE){cout<<"Done! ("<<Sz2<<" points)\t Calculated: "<<y2[Sz2-1][0]<<endl;}
	y2.clear(); r2.clear();

	// Guess Coefficient Assuming Linear Relationship
	//double C=(C1-C2)*(Zeta-y2[Sz2-1][0])/(y1[Sz1-1][0]-y2[Sz2-1][0])+C2;
	double C;
	// Acceptable Error
	double Err3;
	RUNNING=true;
	bool nanCHECKING=true;
	int Sz3;
	Counter=0;
	while(RUNNING)
		{// Define New Coefficient Based on Error (C2=next to last,C1=last)
		C=(C2-C1)*(Err2/2-Err1)/(Err2-Err1)+C1;
		if(VERBOSE){cout<<"Next Guess:\tC = "<<C<<"...";}
		nanCHECKING=true;
		while(nanCHECKING)
			{initialize_PBE(y,r0,C,S);
			steps=integrate_adaptive(controlled_stepper,poisson_boltzmann,y,r0,rf,dr,push_back_state_and_time(y3,r3));
			Sz3=steps+1;
			if(!isnan(y3[Sz3-1][0]))
				{// Good Choice of C
				nanCHECKING=false;break;}
			else
				{Counter++;
				if(Counter>=MAX_PBE_ITERATIONS){return nan("");}
				// NaN Error
				if(VERBOSE){cout<<"Output NaN! Retrying...\n";}
				// Reduce C
				C/=2;
				y3.clear();r3.clear();}
			}
		Err3=Zeta-y3[Sz3-1][0];
		if(VERBOSE){cout<<"Done! ("<<Sz3<<" points)\t|Calculated: "<<y3[Sz3-1][0]<<"|Error %: "<<abs(Err3/Zeta)<<"|"<<ACCEPTED_ERROR_PERCENT<<endl;}
		
		if(abs(Err3/Zeta)>ACCEPTED_ERROR_PERCENT)
			{// Update Coefficient Before Last
			C1=C2;Err1=Err2;
			// Update Last Coefficient
			C2=C;Err2=Err3;			
			// Clear Memory
			y3.clear();r3.clear();}
		else{RUNNING=false;break;}
		}

	return C;}

	// y[0] = f
	// y[1] = df/dr
	// y[2] = d2f/dr2
	// y[3] = d3f/dr3
	// y[4] = d4f/dr4
	// Equilibrium Electric Pontential
	// y[5] = potential
	// y[6] = d(potential)/dr
	// Induced Electrochemical Potential Function (9 ion)
	// y[7] = IEP1
	// y[8] = d(IEP1)/dr
	// y[9] = IEP2
	// y[10] = d(IEP2)/dr
	// y[11] = IEP3
	// y[12] = d(IEP3)/dr
	// y[13] = IEP4
	// y[14] = d(IEP4)/dr
	// y[15] = IEP5
	// y[16] = d(IEP5)/dr
	// y[17] = IEP6
	// y[18] = d(IEP6)/dr
	// y[19] = IEP7
	// y[20] = d(IEP7)/dr
	// y[21] = IEP8
	// y[22] = d(IEP8)/dr
	// y[23] = IEP9
	// y[24] = d(IEP9)/dr
	// Constants
	// y[25] = Temp					// [Kelvin]
	// y[26] = numIon
	// Ion Properties
	// y[27] = Ion1.val;		
	// y[28] = Ion1.conc;		// mM [mol/m^3]
	// y[29] = Ion1.radius;		// meters
	// y[30] = Ion2.val;
	// y[31] = Ion2.conc;
	// y[32] = Ion2.radius;
	// y[33] = Ion3.val;
	// y[34] = Ion3.conc;
	// y[35] = Ion3.radius;
	// y[36] = Ion4.val;
	// y[37] = Ion4.conc;
	// y[38] = Ion4.radius;
	// y[39] = Ion5.val;
	// y[40] = Ion5.conc;
	// y[41] = Ion5.radius;
	// y[42] = Ion6.val;
	// y[43] = Ion6.conc;
	// y[44] = Ion6.radius;
	// y[45] = Ion7.val;
	// y[46] = Ion7.conc;
	// y[47] = Ion7.radius;
	// y[48] = Ion8.val;
	// y[49] = Ion8.conc;
	// y[50] = Ion8.radius;
	// y[51] = Ion9.val;
	// y[52] = Ion9.conc;
	// y[53] = Ion9.radius;
	// Solution Properties
	// y[54] = relative dielectric	[dimensionless]
	// y[55] = Debye Length [meters]
	// y[56] = viscosity  [Pa s]
void homogeneous_form_9ion(const state_type &y,state_type &dy,const double r)
	{double T=y[25];									// temperature [Kelvin]
	double viscosity=y[56];							// [Newton*second/meter^2]
	// Organize Ion Data
	int numIon=y[26];	
	double* ion_conc=new double[numIon];
	int* ion_val=new int[numIon];
	double* ion_rad=new double[numIon];
	double* ion_mob=new double[numIon];
	double* ion_fc=new double[numIon];
	// 1st Ion Properties
	ion_val[0]=y[27];
   ion_conc[0]=y[28];  		// mM or mol/m^3	    
	ion_rad[0]=y[29];		// meters
	ion_mob[0]=calc_Ion_Mobility(ion_val[0],ion_rad[0],viscosity);
	ion_fc[0]=calc_Ion_Friction_Coefficient(ion_val[0],ion_mob[0]);
	// 2nd Ion Properties
	ion_val[1]=y[30];
   ion_conc[1]=y[31];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[1]=y[32];													// [meter]
	ion_mob[1]=calc_Ion_Mobility(ion_val[1],ion_rad[1],viscosity);		// [meter^2/Volt*second]
	ion_fc[1]=calc_Ion_Friction_Coefficient(ion_val[1],ion_mob[1]);		// [Newton*second/meter]
	// 3rd Ion Properties
	ion_val[2]=y[33];
   ion_conc[2]=y[34];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[2]=y[35];													// [meter]
	ion_mob[2]=calc_Ion_Mobility(ion_val[2],ion_rad[2],viscosity);		// [meter^2/Volt*second]
	ion_fc[2]=calc_Ion_Friction_Coefficient(ion_val[2],ion_mob[2]);		// [Newton*second/meter]
	// 4th Ion Properties
	ion_val[3]=y[36];
   ion_conc[3]=y[37];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[3]=y[38];													// [meter]
	ion_mob[3]=calc_Ion_Mobility(ion_val[3],ion_rad[3],viscosity);		// [meter^2/Volt*second]
	ion_fc[3]=calc_Ion_Friction_Coefficient(ion_val[3],ion_mob[3]);		// [Newton*second/meter]
	// 5th Ion Properties
	ion_val[4]=y[39];
   ion_conc[4]=y[40];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[4]=y[41];													// [meter]
	ion_mob[4]=calc_Ion_Mobility(ion_val[4],ion_rad[4],viscosity);		// [meter^2/Volt*second]
	ion_fc[4]=calc_Ion_Friction_Coefficient(ion_val[4],ion_mob[4]);		// [Newton*second/meter]
	// 6th Ion Properties
	ion_val[5]=y[42];
   ion_conc[5]=y[43];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[5]=y[44];													// [meter]
	ion_mob[5]=calc_Ion_Mobility(ion_val[5],ion_rad[5],viscosity);		// [meter^2/Volt*second]
	ion_fc[5]=calc_Ion_Friction_Coefficient(ion_val[5],ion_mob[5]);		// [Newton*second/meter]
	// 7th Ion Properties
	ion_val[6]=y[45];
   ion_conc[6]=y[46];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[6]=y[47];													// [meter]
	ion_mob[6]=calc_Ion_Mobility(ion_val[6],ion_rad[6],viscosity);		// [meter^2/Volt*second]
	ion_fc[6]=calc_Ion_Friction_Coefficient(ion_val[6],ion_mob[6]);		// [Newton*second/meter]
	// 8th Ion Properties
	ion_val[7]=y[48];
   ion_conc[7]=y[49];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[7]=y[50];													// [meter]
	ion_mob[7]=calc_Ion_Mobility(ion_val[7],ion_rad[7],viscosity);		// [meter^2/Volt*second]
	ion_fc[7]=calc_Ion_Friction_Coefficient(ion_val[7],ion_mob[7]);		// [Newton*second/meter]
	// 9th Ion Properties
	ion_val[8]=y[51];
   ion_conc[8]=y[52];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[8]=y[53];													// [meter]
	ion_mob[8]=calc_Ion_Mobility(ion_val[8],ion_rad[8],viscosity);		// [meter^2/Volt*second]
	ion_fc[8]=calc_Ion_Friction_Coefficient(ion_val[8],ion_mob[8]);		// [Newton*second/meter]

	double Debye=y[55];				// meters	
	double dielectric=y[54];
	// Convert Concentration to Number Density (particles per volume (m^3))
	ion_conc[0]*=Na;
	ion_conc[1]*=Na;
	ion_conc[2]*=Na;
	ion_conc[3]*=Na;
	ion_conc[4]*=Na;
	ion_conc[5]*=Na;
	ion_conc[6]*=Na;
	ion_conc[7]*=Na;
	ion_conc[8]*=Na;
	// Calculation Constants
	double A;
	int ionIndex;
    // Linearized Governing Equations
	// F FUNCTION
	// df/dr = y[1]
	dy[0]=y[1];
	// d2f/dr2 = y[2]
	dy[1]=y[2];
	// d3f/dr3 = y[3]
	dy[2]=y[3];
	// d4f/dr4 = y[4]
	dy[3]=y[4];
	// d5f/dr5
	// 1st Ion Equilibrium Number Density
	double dN1_dr=(-ion_conc[0]*ion_val[0]*ec/(kb*T))*y[6]*exp(-ec*ion_val[0]*y[5]/(kb*T));
	double trm1=ion_val[0]*ec*dN1_dr*y[7];
	// 2nd Ion Equilibrium Number Density
	double dN2_dr=(-ion_conc[1]*ion_val[1]*ec/(kb*T))*y[6]*exp(-ec*ion_val[1]*y[5]/(kb*T));
	double trm2=ion_val[1]*ec*dN2_dr*y[9];
	// 3rd Ion Equilibrium Number Density
	double dN3_dr=(-ion_conc[2]*ion_val[2]*ec/(kb*T))*y[6]*exp(-ec*ion_val[2]*y[5]/(kb*T));
	double trm3=ion_val[2]*ec*dN3_dr*y[11];
	// 4th Ion Equilibrium Number Density
	double dN4_dr=(-ion_conc[3]*ion_val[3]*ec/(kb*T))*y[6]*exp(-ec*ion_val[3]*y[5]/(kb*T));
	double trm4=ion_val[3]*ec*dN4_dr*y[13];
	// 5th Ion Equilibrium Number Density
	double dN5_dr=(-ion_conc[4]*ion_val[4]*ec/(kb*T))*y[6]*exp(-ec*ion_val[4]*y[5]/(kb*T));
	double trm5=ion_val[4]*ec*dN5_dr*y[15];
	// 6th Ion Equilibrium Number Density
	double dN6_dr=(-ion_conc[5]*ion_val[5]*ec/(kb*T))*y[6]*exp(-ec*ion_val[5]*y[5]/(kb*T));
	double trm6=ion_val[5]*ec*dN6_dr*y[17];
	// 7th Ion Equilibrium Number Density
	double dN7_dr=(-ion_conc[6]*ion_val[6]*ec/(kb*T))*y[6]*exp(-ec*ion_val[6]*y[5]/(kb*T));
	double trm7=ion_val[6]*ec*dN7_dr*y[19];
	// 8th Ion Equilibrium Number Density
	double dN8_dr=(-ion_conc[7]*ion_val[7]*ec/(kb*T))*y[6]*exp(-ec*ion_val[7]*y[5]/(kb*T));
	double trm8=ion_val[7]*ec*dN8_dr*y[21];
	// 9th Ion Equilibrium Number Density
	double dN9_dr=(-ion_conc[8]*ion_val[8]*ec/(kb*T))*y[6]*exp(-ec*ion_val[8]*y[5]/(kb*T));
	double trm9=ion_val[8]*ec*dN8_dr*y[23];
	//
	dy[4]=-4*y[4]/r + 4*y[3]/(r*r) + (trm1+trm2+trm3+trm4+trm5+trm6+trm7+trm8)/(viscosity*r);
	// EQUILIBRIUM POTENTIAL
	// d(potential)/dr = y[6]
	dy[5]=y[6];
	// d2(potential)/dr2
	A=0;
	for(int i=0;i<numIon;i++){A+=ion_val[i]*ion_conc[i]*exp(-ec*ion_val[i]*y[5]/(kb*T));}
	dy[6]=-ec*A/(dielectric*eo) - 2*y[6]/r;
	// INDUCED ELECTROCHEMICAL POTENTIAL FUNCTION (2 ions)
	// 1st Ion
	// d(IEP)/dr
   dy[7]=y[8];
	// d2(IEP)/dr2
	A=(ion_val[0]*ec/(kb*T))*y[6];
   dy[8]=-2*y[8]/r + 2*y[7]/(r*r) + A*(y[8] - (2*ion_fc[0]/(r*ion_val[0]*ec))*y[1]);
	// 2nd Ion
	// d(IEP)/dr
   dy[9]=y[10];
	// d2(IEP)/dr2
	A=(ion_val[1]*ec/(kb*T))*y[6];
   dy[10]=-2*y[10]/r + 2*y[9]/(r*r) + A*(y[10] - (2*ion_fc[1]/(r*ion_val[1]*ec))*y[1]);
	// 3rd Ion
	// d(IEP)/dr
   dy[11]=y[12];
	// d2(IEP)/dr2
	A=(ion_val[2]*ec/(kb*T))*y[6];
   dy[12]=-2*y[12]/r + 2*y[11]/(r*r) + A*(y[12] - (2*ion_fc[2]/(r*ion_val[2]*ec))*y[1]);
	// 4th Ion
	// d(IEP)/dr
   dy[13]=y[14];
	// d2(IEP)/dr2
	A=(ion_val[3]*ec/(kb*T))*y[6];
   dy[14]=-2*y[14]/r + 2*y[13]/(r*r) + A*(y[14] - (2*ion_fc[3]/(r*ion_val[3]*ec))*y[1]);
	// 5th Ion
	// d(IEP)/dr
   dy[15]=y[16];
	// d2(IEP)/dr2
	A=(ion_val[4]*ec/(kb*T))*y[6];
   dy[16]=-2*y[16]/r + 2*y[15]/(r*r) + A*(y[16] - (2*ion_fc[4]/(r*ion_val[4]*ec))*y[1]);
	// 6th Ion
	// d(IEP)/dr
   dy[17]=y[18];
	// d2(IEP)/dr2
	A=(ion_val[5]*ec/(kb*T))*y[6];
   dy[18]=-2*y[18]/r + 2*y[17]/(r*r) + A*(y[18] - (2*ion_fc[5]/(r*ion_val[5]*ec))*y[1]);
	// 7th Ion
	// d(IEP)/dr
   dy[19]=y[20];
	// d2(IEP)/dr2
	A=(ion_val[6]*ec/(kb*T))*y[6];
   dy[20]=-2*y[20]/r + 2*y[19]/(r*r) + A*(y[20] - (2*ion_fc[6]/(r*ion_val[6]*ec))*y[1]);
	// 8th Ion
	// d(IEP)/dr
   dy[21]=y[22];
	// d2(IEP)/dr2
	A=(ion_val[7]*ec/(kb*T))*y[6];
   dy[22]=-2*y[22]/r + 2*y[21]/(r*r) + A*(y[22] - (2*ion_fc[7]/(r*ion_val[7]*ec))*y[1]);
	// 9th Ion
	// d(IEP)/dr
   dy[23]=y[24];
	// d2(IEP)/dr2
	A=(ion_val[8]*ec/(kb*T))*y[6];
   dy[24]=-2*y[24]/r + 2*y[23]/(r*r) + A*(y[24] - (2*ion_fc[8]/(r*ion_val[8]*ec))*y[1]);

	// Clear Memory
	delete [] ion_conc;delete [] ion_val;delete [] ion_rad;delete [] ion_mob;delete [] ion_fc;}

void inhomogeneous_Prob1_9ion(const state_type &y,state_type &dy,const double r)
	{double T=y[25];									// temperature [Kelvin]
	double viscosity=y[56];							// [Newton*second/meter^2]
	// Organize Ion Data
	int numIon=y[26];	
	double* ion_conc=new double[numIon];
	int* ion_val=new int[numIon];
	double* ion_rad=new double[numIon];
	double* ion_mob=new double[numIon];
	double* ion_fc=new double[numIon];
	// 1st Ion Properties
	ion_val[0]=y[27];
   ion_conc[0]=y[28];  		// mM or mol/m^3	    
	ion_rad[0]=y[29];		// meters
	ion_mob[0]=calc_Ion_Mobility(ion_val[0],ion_rad[0],viscosity);
	ion_fc[0]=calc_Ion_Friction_Coefficient(ion_val[0],ion_mob[0]);
	// 2nd Ion Properties
	ion_val[1]=y[30];
   ion_conc[1]=y[31];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[1]=y[32];													// [meter]
	ion_mob[1]=calc_Ion_Mobility(ion_val[1],ion_rad[1],viscosity);		// [meter^2/Volt*second]
	ion_fc[1]=calc_Ion_Friction_Coefficient(ion_val[1],ion_mob[1]);		// [Newton*second/meter]
	// 3rd Ion Properties
	ion_val[2]=y[33];
   ion_conc[2]=y[34];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[2]=y[35];													// [meter]
	ion_mob[2]=calc_Ion_Mobility(ion_val[2],ion_rad[2],viscosity);		// [meter^2/Volt*second]
	ion_fc[2]=calc_Ion_Friction_Coefficient(ion_val[2],ion_mob[2]);		// [Newton*second/meter]
	// 4th Ion Properties
	ion_val[3]=y[36];
   ion_conc[3]=y[37];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[3]=y[38];													// [meter]
	ion_mob[3]=calc_Ion_Mobility(ion_val[3],ion_rad[3],viscosity);		// [meter^2/Volt*second]
	ion_fc[3]=calc_Ion_Friction_Coefficient(ion_val[3],ion_mob[3]);		// [Newton*second/meter]
	// 5th Ion Properties
	ion_val[4]=y[39];
   ion_conc[4]=y[40];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[4]=y[41];													// [meter]
	ion_mob[4]=calc_Ion_Mobility(ion_val[4],ion_rad[4],viscosity);		// [meter^2/Volt*second]
	ion_fc[4]=calc_Ion_Friction_Coefficient(ion_val[4],ion_mob[4]);		// [Newton*second/meter]
	// 6th Ion Properties
	ion_val[5]=y[42];
   ion_conc[5]=y[43];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[5]=y[44];													// [meter]
	ion_mob[5]=calc_Ion_Mobility(ion_val[5],ion_rad[5],viscosity);		// [meter^2/Volt*second]
	ion_fc[5]=calc_Ion_Friction_Coefficient(ion_val[5],ion_mob[5]);		// [Newton*second/meter]
	// 7th Ion Properties
	ion_val[6]=y[45];
   ion_conc[6]=y[46];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[6]=y[47];													// [meter]
	ion_mob[6]=calc_Ion_Mobility(ion_val[6],ion_rad[6],viscosity);		// [meter^2/Volt*second]
	ion_fc[6]=calc_Ion_Friction_Coefficient(ion_val[6],ion_mob[6]);		// [Newton*second/meter]
	// 8th Ion Properties
	ion_val[7]=y[48];
   ion_conc[7]=y[49];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[7]=y[50];													// [meter]
	ion_mob[7]=calc_Ion_Mobility(ion_val[7],ion_rad[7],viscosity);		// [meter^2/Volt*second]
	ion_fc[7]=calc_Ion_Friction_Coefficient(ion_val[7],ion_mob[7]);		// [Newton*second/meter]
	// 9th Ion Properties
	ion_val[8]=y[51];
   ion_conc[8]=y[52];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[8]=y[53];													// [meter]
	ion_mob[8]=calc_Ion_Mobility(ion_val[8],ion_rad[8],viscosity);		// [meter^2/Volt*second]
	ion_fc[8]=calc_Ion_Friction_Coefficient(ion_val[8],ion_mob[8]);		// [Newton*second/meter]

	double Debye=y[55];				// meters	
	double dielectric=y[54];
	// Convert Concentration to Number Density (particles per volume (m^3))
	ion_conc[0]*=Na;
	ion_conc[1]*=Na;
	ion_conc[2]*=Na;
	ion_conc[3]*=Na;
	ion_conc[4]*=Na;
	ion_conc[5]*=Na;
	ion_conc[6]*=Na;
	ion_conc[7]*=Na;
	ion_conc[8]*=Na;
	// Calculation Constants
	double A,dN_dr;
	int ionIndex;
	// Linearized Governing Equations
	// F FUNCTION
	// df/dr = y[1]
	dy[0]=y[1];
	// d2f/dr2 = y[2]
	dy[1]=y[2];
	// d3f/dr3 = y[3]
	dy[2]=y[3];
	// d4f/dr4 = y[4]
	dy[3]=y[4];
	// d5f/dr5
	A=0;
	for(int i=0;i<numIon;i++)
		{// Ion Equilibrium Number Density
		dN_dr=(-ion_conc[i]*ion_val[i]*ec/(kb*T))*y[6]*exp(-ec*ion_val[i]*y[5]/(kb*T));
		// Ion Contribution
		ionIndex=7+2*i; 
		A+=ion_val[i]*ec*dN_dr*y[ionIndex];}
    dy[4]=-4*y[4]/r + 4*y[3]/(r*r) + A/(viscosity*r);
	// EQUILIBRIUM POTENTIAL
	// d(potential)/dr = y[6]
	dy[5]=y[6];
	// d2(potential)/dr2
	A=0;
	for(int i=0;i<numIon;i++){A+=ion_val[i]*ion_conc[i]*exp(-ec*ion_val[i]*y[5]/(kb*T));}
	dy[6]=-ec*A/(dielectric*eo) - 2*y[6]/r;
	// INDUCED ELECTROCHEMICAL POTENTIAL FUNCTION (2 ions)
	// 1st Ion
	// d(IEP)/dr = y[8]
   dy[7]=y[8];
	// d2(IEP)/dr2
	A=(ion_val[0]*ec/(kb*T))*y[6];
   dy[8]=-2*y[8]/r + 2*y[7]/(r*r) + A*(y[8] - (2*ion_fc[0]/(r*ion_val[0]*ec))*y[1] - ion_fc[0]/(ion_val[0]*ec));
	// 2nd Ion
	// d(IEP)/dr = y[10]
   dy[9]=y[10];
	// d2(IEP)/dr2
	A=(ion_val[1]*ec/(kb*T))*y[6];
   dy[10]=-2*y[10]/r + 2*y[9]/(r*r) + A*(y[10] - (2*ion_fc[1]/(r*ion_val[1]*ec))*y[1] - ion_fc[1]/(ion_val[1]*ec));
	// 3rd Ion
	// d(IEP)/dr
   dy[11]=y[12];
	// d2(IEP)/dr2
	A=(ion_val[2]*ec/(kb*T))*y[6];
   dy[12]=-2*y[12]/r + 2*y[11]/(r*r) + A*(y[12] - (2*ion_fc[2]/(r*ion_val[2]*ec))*y[1] - ion_fc[2]/(ion_val[2]*ec));
	// 4th Ion
	// d(IEP)/dr
   dy[13]=y[14];
	// d2(IEP)/dr2
	A=(ion_val[3]*ec/(kb*T))*y[6];
   dy[14]=-2*y[14]/r + 2*y[13]/(r*r) + A*(y[14] - (2*ion_fc[3]/(r*ion_val[3]*ec))*y[1] - ion_fc[3]/(ion_val[3]*ec));
	// 5th Ion
	// d(IEP)/dr
   dy[15]=y[16];
	// d2(IEP)/dr2
	A=(ion_val[4]*ec/(kb*T))*y[6];
   dy[16]=-2*y[16]/r + 2*y[15]/(r*r) + A*(y[16] - (2*ion_fc[4]/(r*ion_val[4]*ec))*y[1] - ion_fc[4]/(ion_val[4]*ec));
	// 6th Ion
	// d(IEP)/dr
   dy[17]=y[18];
	// d2(IEP)/dr2
	A=(ion_val[5]*ec/(kb*T))*y[6];
   dy[18]=-2*y[18]/r + 2*y[17]/(r*r) + A*(y[18] - (2*ion_fc[5]/(r*ion_val[5]*ec))*y[1] - ion_fc[5]/(ion_val[5]*ec));
	// 7th Ion
	// d(IEP)/dr
   dy[19]=y[20];
	// d2(IEP)/dr2
	A=(ion_val[6]*ec/(kb*T))*y[6];
   dy[20]=-2*y[20]/r + 2*y[19]/(r*r) + A*(y[20] - (2*ion_fc[6]/(r*ion_val[6]*ec))*y[1] - ion_fc[6]/(ion_val[6]*ec));
	// 8th Ion
	// d(IEP)/dr
   dy[21]=y[22];
	// d2(IEP)/dr2
	A=(ion_val[7]*ec/(kb*T))*y[6];
   dy[22]=-2*y[22]/r + 2*y[21]/(r*r) + A*(y[22] - (2*ion_fc[7]/(r*ion_val[7]*ec))*y[1] - ion_fc[7]/(ion_val[7]*ec));
	// 9th Ion
	// d(IEP)/dr
   dy[23]=y[24];
	// d2(IEP)/dr2
	A=(ion_val[8]*ec/(kb*T))*y[6];
   dy[24]=-2*y[24]/r + 2*y[23]/(r*r) + A*(y[24] - (2*ion_fc[8]/(r*ion_val[8]*ec))*y[1] - ion_fc[8]/(ion_val[8]*ec));

	// Clear Memory
	delete [] ion_conc;delete [] ion_val;delete [] ion_rad;delete [] ion_mob;delete [] ion_fc;}

void inhomogeneous_Prob2_9ion(const state_type &y,state_type &dy,const double r)
	{double T=y[25];									// temperature [Kelvin]
	double viscosity=y[56];							// [Newton*second/meter^2]
	// Organize Ion Data
	int numIon=y[26];	
	double* ion_conc=new double[numIon];
	int* ion_val=new int[numIon];
	double* ion_rad=new double[numIon];
	double* ion_mob=new double[numIon];
	double* ion_fc=new double[numIon];
	// 1st Ion Properties
	ion_val[0]=y[27];
   ion_conc[0]=y[28];  		// mM or mol/m^3	    
	ion_rad[0]=y[29];		// meters
	ion_mob[0]=calc_Ion_Mobility(ion_val[0],ion_rad[0],viscosity);
	ion_fc[0]=calc_Ion_Friction_Coefficient(ion_val[0],ion_mob[0]);
	// 2nd Ion Properties
	ion_val[1]=y[30];
   ion_conc[1]=y[31];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[1]=y[32];													// [meter]
	ion_mob[1]=calc_Ion_Mobility(ion_val[1],ion_rad[1],viscosity);		// [meter^2/Volt*second]
	ion_fc[1]=calc_Ion_Friction_Coefficient(ion_val[1],ion_mob[1]);		// [Newton*second/meter]
	// 3rd Ion Properties
	ion_val[2]=y[33];
   ion_conc[2]=y[34];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[2]=y[35];													// [meter]
	ion_mob[2]=calc_Ion_Mobility(ion_val[2],ion_rad[2],viscosity);		// [meter^2/Volt*second]
	ion_fc[2]=calc_Ion_Friction_Coefficient(ion_val[2],ion_mob[2]);		// [Newton*second/meter]
	// 4th Ion Properties
	ion_val[3]=y[36];
   ion_conc[3]=y[37];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[3]=y[38];													// [meter]
	ion_mob[3]=calc_Ion_Mobility(ion_val[3],ion_rad[3],viscosity);		// [meter^2/Volt*second]
	ion_fc[3]=calc_Ion_Friction_Coefficient(ion_val[3],ion_mob[3]);		// [Newton*second/meter]
	// 5th Ion Properties
	ion_val[4]=y[39];
   ion_conc[4]=y[40];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[4]=y[41];													// [meter]
	ion_mob[4]=calc_Ion_Mobility(ion_val[4],ion_rad[4],viscosity);		// [meter^2/Volt*second]
	ion_fc[4]=calc_Ion_Friction_Coefficient(ion_val[4],ion_mob[4]);		// [Newton*second/meter]
	// 6th Ion Properties
	ion_val[5]=y[42];
   ion_conc[5]=y[43];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[5]=y[44];													// [meter]
	ion_mob[5]=calc_Ion_Mobility(ion_val[5],ion_rad[5],viscosity);		// [meter^2/Volt*second]
	ion_fc[5]=calc_Ion_Friction_Coefficient(ion_val[5],ion_mob[5]);		// [Newton*second/meter]
	// 7th Ion Properties
	ion_val[6]=y[45];
   ion_conc[6]=y[46];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[6]=y[47];													// [meter]
	ion_mob[6]=calc_Ion_Mobility(ion_val[6],ion_rad[6],viscosity);		// [meter^2/Volt*second]
	ion_fc[6]=calc_Ion_Friction_Coefficient(ion_val[6],ion_mob[6]);		// [Newton*second/meter]
	// 8th Ion Properties
	ion_val[7]=y[48];
   ion_conc[7]=y[49];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[7]=y[50];													// [meter]
	ion_mob[7]=calc_Ion_Mobility(ion_val[7],ion_rad[7],viscosity);		// [meter^2/Volt*second]
	ion_fc[7]=calc_Ion_Friction_Coefficient(ion_val[7],ion_mob[7]);		// [Newton*second/meter]
	// 9th Ion Properties
	ion_val[8]=y[51];
   ion_conc[8]=y[52];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[8]=y[53];													// [meter]
	ion_mob[8]=calc_Ion_Mobility(ion_val[8],ion_rad[8],viscosity);		// [meter^2/Volt*second]
	ion_fc[8]=calc_Ion_Friction_Coefficient(ion_val[8],ion_mob[8]);		// [Newton*second/meter]

	double Debye=y[55];				// meters	
	double dielectric=y[54];
	// Convert Concentration to Number Density (particles per volume (m^3))
	ion_conc[0]*=Na;
	ion_conc[1]*=Na;
	ion_conc[2]*=Na;
	ion_conc[3]*=Na;
	ion_conc[4]*=Na;
	ion_conc[5]*=Na;
	ion_conc[6]*=Na;
	ion_conc[7]*=Na;
	ion_conc[8]*=Na;
    // Calculation Constants
	double A,dN_dr;
	int ionIndex;
    // Linearized Governing Equations
	// F FUNCTION
	// df/dr = y[1]
	dy[0]=y[1];
	// d2f/dr2 = y[2]
	dy[1]=y[2];
	// d3f/dr3 = y[3]
	dy[2]=y[3];
	// d4f/dr4 = y[4]
	dy[3]=y[4];
	// d5f/dr5
	A=0;
	for(int i=0;i<numIon;i++)
		{// Ion Equilibrium Number Density
		dN_dr=(-ion_conc[i]*ion_val[i]*ec/(kb*T))*y[6]*exp(-ec*ion_val[i]*y[5]/(kb*T));
		// Ion Contribution
		ionIndex=7+2*i; 
		A+=ion_val[i]*ec*dN_dr*(y[ionIndex]+r);}
   dy[4]=-4*y[4]/r + 4*y[3]/(r*r) + A/(viscosity*r);
	// EQUILIBRIUM POTENTIAL
	// d(potential)/dr = y[6]
	dy[5]=y[6];
	// d2(potential)/dr2
	A=0;
	for(int i=0;i<numIon;i++){A+=ion_val[i]*ion_conc[i]*exp(-ec*ion_val[i]*y[5]/(kb*T));}
	dy[6]=-ec*A/(dielectric*eo) - 2*y[6]/r;
	// INDUCED ELECTROCHEMICAL POTENTIAL FUNCTION (2 ions)
	// 1st Ion
	// d(IEP)/dr = y[8]
   dy[7]=y[8];
	// d2(IEP)/dr2
	A=(ion_val[0]*ec/(kb*T))*y[6];
   dy[8]=-2*y[8]/r + 2*y[7]/(r*r) + A*(y[8] - (2*ion_fc[0]/(r*ion_val[0]*ec))*y[1] + 1);
	// 2nd Ion
	// d(IEP)/dr = y[10]
   dy[9]=y[10];
	// d2(IEP)/dr2
	A=(ion_val[1]*ec/(kb*T))*y[6];
   dy[10]=-2*y[10]/r + 2*y[9]/(r*r) + A*(y[10] - (2*ion_fc[1]/(r*ion_val[1]*ec))*y[1] + 1);
	// 3rd Ion
	// d(IEP)/dr
   dy[11]=y[12];
	// d2(IEP)/dr2
	A=(ion_val[2]*ec/(kb*T))*y[6];
   dy[12]=-2*y[12]/r + 2*y[11]/(r*r) + A*(y[12] - (2*ion_fc[2]/(r*ion_val[2]*ec))*y[1] + 1);
	// 4th Ion
	// d(IEP)/dr
   dy[13]=y[14];
	// d2(IEP)/dr2
	A=(ion_val[3]*ec/(kb*T))*y[6];
   dy[14]=-2*y[14]/r + 2*y[13]/(r*r) + A*(y[14] - (2*ion_fc[3]/(r*ion_val[3]*ec))*y[1] + 1);
	// 5th Ion
	// d(IEP)/dr
   dy[15]=y[16];
	// d2(IEP)/dr2
	A=(ion_val[4]*ec/(kb*T))*y[6];
   dy[16]=-2*y[16]/r + 2*y[15]/(r*r) + A*(y[16] - (2*ion_fc[4]/(r*ion_val[4]*ec))*y[1] + 1);
	// 6th Ion
	// d(IEP)/dr
   dy[17]=y[18];
	// d2(IEP)/dr2
	A=(ion_val[5]*ec/(kb*T))*y[6];
   dy[18]=-2*y[18]/r + 2*y[17]/(r*r) + A*(y[18] - (2*ion_fc[5]/(r*ion_val[5]*ec))*y[1] + 1);
	// 7th Ion
	// d(IEP)/dr
   dy[19]=y[20];
	// d2(IEP)/dr2
	A=(ion_val[6]*ec/(kb*T))*y[6];
   dy[20]=-2*y[20]/r + 2*y[19]/(r*r) + A*(y[20] - (2*ion_fc[6]/(r*ion_val[6]*ec))*y[1] + 1);
	// 8th Ion
	// d(IEP)/dr
   dy[21]=y[22];
	// d2(IEP)/dr2
	A=(ion_val[7]*ec/(kb*T))*y[6];
   dy[22]=-2*y[22]/r + 2*y[21]/(r*r) + A*(y[22] - (2*ion_fc[7]/(r*ion_val[7]*ec))*y[1] + 1);
	// 9th Ion
	// d(IEP)/dr
   dy[23]=y[24];
	// d2(IEP)/dr2
	A=(ion_val[8]*ec/(kb*T))*y[6];
   dy[24]=-2*y[24]/r + 2*y[23]/(r*r) + A*(y[24] - (2*ion_fc[8]/(r*ion_val[8]*ec))*y[1] + 1);

	// Clear Memory
	delete [] ion_conc;delete [] ion_val;delete [] ion_rad;delete [] ion_mob;delete [] ion_fc;}


	// y[0] = f
	// y[1] = df/dr
	// y[2] = d2f/dr2
	// y[3] = d3f/dr3
	// y[4] = d4f/dr4
	// Equilibrium Electric Pontential
	// y[5] = potential
	// y[6] = d(potential)/dr
	// Induced Electrochemical Potential Function (3 ion)
	// y[7] = IEP1
	// y[8] = d(IEP1)/dr
	// y[9] = IEP2
	// y[10] = d(IEP2)/dr
	// y[11] = IEP3
	// y[12] = d(IEP3)/dr
	// y[13] = IEP4
	// y[14] = d(IEP4)/dr
	// y[15] = IEP5
	// y[16] = d(IEP5)/dr
	// y[17] = IEP6
	// y[18] = d(IEP6)/dr
	// y[19] = IEP7
	// y[20] = d(IEP7)/dr
	// y[21] = IEP8
	// y[22] = d(IEP8)/dr
	// Constants
	// y[23] = Temp					// [Kelvin]
	// y[24] = numIon
	// Ion Properties
	// y[25] = Ion1.val;		
	// y[26] = Ion1.conc;		// mM [mol/m^3]
	// y[27] = Ion1.radius;		// meters
	// y[28] = Ion2.val;
	// y[29] = Ion2.conc;
	// y[30] = Ion2.radius;
	// y[31] = Ion3.val;
	// y[32] = Ion3.conc;
	// y[33] = Ion3.radius;
	// y[34] = Ion4.val;
	// y[35] = Ion4.conc;
	// y[36] = Ion4.radius;
	// y[37] = Ion5.val;
	// y[38] = Ion5.conc;
	// y[39] = Ion5.radius;
	// y[40] = Ion6.val;
	// y[41] = Ion6.conc;
	// y[42] = Ion6.radius;
	// y[43] = Ion7.val;
	// y[44] = Ion7.conc;
	// y[45] = Ion7.radius;
	// y[46] = Ion8.val;
	// y[47] = Ion8.conc;
	// y[48] = Ion8.radius;
	// Solution Properties
	// y[49] = relative dielectric	[dimensionless]
	// y[50] = Debye Length [meters]
	// y[51] = viscosity  [Pa s]
void homogeneous_form_8ion(const state_type &y,state_type &dy,const double r)
	{double T=y[23];									// temperature [Kelvin]
	double viscosity=y[51];							// [Newton*second/meter^2]
	// Organize Ion Data
	int numIon=y[24];	
	double* ion_conc=new double[numIon];
	int* ion_val=new int[numIon];
	double* ion_rad=new double[numIon];
	double* ion_mob=new double[numIon];
	double* ion_fc=new double[numIon];
	// 1st Ion Properties
	ion_val[0]=y[25];
   ion_conc[0]=y[26];  		// mM or mol/m^3	    
	ion_rad[0]=y[27];		// meters
	ion_mob[0]=calc_Ion_Mobility(ion_val[0],ion_rad[0],viscosity);
	ion_fc[0]=calc_Ion_Friction_Coefficient(ion_val[0],ion_mob[0]);
	// 2nd Ion Properties
	ion_val[1]=y[28];
   ion_conc[1]=y[29];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[1]=y[30];													// [meter]
	ion_mob[1]=calc_Ion_Mobility(ion_val[1],ion_rad[1],viscosity);		// [meter^2/Volt*second]
	ion_fc[1]=calc_Ion_Friction_Coefficient(ion_val[1],ion_mob[1]);		// [Newton*second/meter]
	// 3rd Ion Properties
	ion_val[2]=y[31];
   ion_conc[2]=y[32];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[2]=y[33];													// [meter]
	ion_mob[2]=calc_Ion_Mobility(ion_val[2],ion_rad[2],viscosity);		// [meter^2/Volt*second]
	ion_fc[2]=calc_Ion_Friction_Coefficient(ion_val[2],ion_mob[2]);		// [Newton*second/meter]
	// 4th Ion Properties
	ion_val[3]=y[34];
   ion_conc[3]=y[35];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[3]=y[36];													// [meter]
	ion_mob[3]=calc_Ion_Mobility(ion_val[3],ion_rad[3],viscosity);		// [meter^2/Volt*second]
	ion_fc[3]=calc_Ion_Friction_Coefficient(ion_val[3],ion_mob[3]);		// [Newton*second/meter]
	// 5th Ion Properties
	ion_val[4]=y[37];
   ion_conc[4]=y[38];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[4]=y[39];													// [meter]
	ion_mob[4]=calc_Ion_Mobility(ion_val[4],ion_rad[4],viscosity);		// [meter^2/Volt*second]
	ion_fc[4]=calc_Ion_Friction_Coefficient(ion_val[4],ion_mob[4]);		// [Newton*second/meter]
	// 6th Ion Properties
	ion_val[5]=y[40];
   ion_conc[5]=y[41];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[5]=y[42];													// [meter]
	ion_mob[5]=calc_Ion_Mobility(ion_val[5],ion_rad[5],viscosity);		// [meter^2/Volt*second]
	ion_fc[5]=calc_Ion_Friction_Coefficient(ion_val[5],ion_mob[5]);		// [Newton*second/meter]
	// 7th Ion Properties
	ion_val[6]=y[43];
   ion_conc[6]=y[44];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[6]=y[45];													// [meter]
	ion_mob[6]=calc_Ion_Mobility(ion_val[6],ion_rad[6],viscosity);		// [meter^2/Volt*second]
	ion_fc[6]=calc_Ion_Friction_Coefficient(ion_val[6],ion_mob[6]);		// [Newton*second/meter]
	// 8th Ion Properties
	ion_val[7]=y[46];
   ion_conc[7]=y[47];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[7]=y[48];													// [meter]
	ion_mob[7]=calc_Ion_Mobility(ion_val[7],ion_rad[7],viscosity);		// [meter^2/Volt*second]
	ion_fc[7]=calc_Ion_Friction_Coefficient(ion_val[7],ion_mob[7]);		// [Newton*second/meter]

	double Debye=y[50];				// meters	
	double dielectric=y[49];
	// Convert Concentration to Number Density (particles per volume (m^3))
	ion_conc[0]*=Na;
	ion_conc[1]*=Na;
	ion_conc[2]*=Na;
	ion_conc[3]*=Na;
	ion_conc[4]*=Na;
	ion_conc[5]*=Na;
	ion_conc[6]*=Na;
	ion_conc[7]*=Na;
	// Calculation Constants
	double A;
	int ionIndex;
    // Linearized Governing Equations
	// F FUNCTION
	// df/dr = y[1]
	dy[0]=y[1];
	// d2f/dr2 = y[2]
	dy[1]=y[2];
	// d3f/dr3 = y[3]
	dy[2]=y[3];
	// d4f/dr4 = y[4]
	dy[3]=y[4];
	// d5f/dr5
	// 1st Ion Equilibrium Number Density
	double dN1_dr=(-ion_conc[0]*ion_val[0]*ec/(kb*T))*y[6]*exp(-ec*ion_val[0]*y[5]/(kb*T));
	double trm1=ion_val[0]*ec*dN1_dr*y[7];
	// 2nd Ion Equilibrium Number Density
	double dN2_dr=(-ion_conc[1]*ion_val[1]*ec/(kb*T))*y[6]*exp(-ec*ion_val[1]*y[5]/(kb*T));
	double trm2=ion_val[1]*ec*dN2_dr*y[9];
	// 3rd Ion Equilibrium Number Density
	double dN3_dr=(-ion_conc[2]*ion_val[2]*ec/(kb*T))*y[6]*exp(-ec*ion_val[2]*y[5]/(kb*T));
	double trm3=ion_val[2]*ec*dN3_dr*y[11];
	// 4th Ion Equilibrium Number Density
	double dN4_dr=(-ion_conc[3]*ion_val[3]*ec/(kb*T))*y[6]*exp(-ec*ion_val[3]*y[5]/(kb*T));
	double trm4=ion_val[3]*ec*dN4_dr*y[13];
	// 5th Ion Equilibrium Number Density
	double dN5_dr=(-ion_conc[4]*ion_val[4]*ec/(kb*T))*y[6]*exp(-ec*ion_val[4]*y[5]/(kb*T));
	double trm5=ion_val[4]*ec*dN5_dr*y[15];
	// 6th Ion Equilibrium Number Density
	double dN6_dr=(-ion_conc[5]*ion_val[5]*ec/(kb*T))*y[6]*exp(-ec*ion_val[5]*y[5]/(kb*T));
	double trm6=ion_val[5]*ec*dN6_dr*y[17];
	// 7th Ion Equilibrium Number Density
	double dN7_dr=(-ion_conc[6]*ion_val[6]*ec/(kb*T))*y[6]*exp(-ec*ion_val[6]*y[5]/(kb*T));
	double trm7=ion_val[6]*ec*dN7_dr*y[19];
	// 8th Ion Equilibrium Number Density
	double dN8_dr=(-ion_conc[7]*ion_val[7]*ec/(kb*T))*y[6]*exp(-ec*ion_val[7]*y[5]/(kb*T));
	double trm8=ion_val[7]*ec*dN8_dr*y[21];
	//
	dy[4]=-4*y[4]/r + 4*y[3]/(r*r) + (trm1+trm2+trm3+trm4+trm5+trm6+trm7+trm8)/(viscosity*r);
	// EQUILIBRIUM POTENTIAL
	// d(potential)/dr = y[6]
	dy[5]=y[6];
	// d2(potential)/dr2
	A=0;
	for(int i=0;i<numIon;i++){A+=ion_val[i]*ion_conc[i]*exp(-ec*ion_val[i]*y[5]/(kb*T));}
	dy[6]=-ec*A/(dielectric*eo) - 2*y[6]/r;
	// INDUCED ELECTROCHEMICAL POTENTIAL FUNCTION (2 ions)
	// 1st Ion
	// d(IEP)/dr
   dy[7]=y[8];
	// d2(IEP)/dr2
	A=(ion_val[0]*ec/(kb*T))*y[6];
   dy[8]=-2*y[8]/r + 2*y[7]/(r*r) + A*(y[8] - (2*ion_fc[0]/(r*ion_val[0]*ec))*y[1]);
	// 2nd Ion
	// d(IEP)/dr
   dy[9]=y[10];
	// d2(IEP)/dr2
	A=(ion_val[1]*ec/(kb*T))*y[6];
   dy[10]=-2*y[10]/r + 2*y[9]/(r*r) + A*(y[10] - (2*ion_fc[1]/(r*ion_val[1]*ec))*y[1]);
	// 3rd Ion
	// d(IEP)/dr
   dy[11]=y[12];
	// d2(IEP)/dr2
	A=(ion_val[2]*ec/(kb*T))*y[6];
   dy[12]=-2*y[12]/r + 2*y[11]/(r*r) + A*(y[12] - (2*ion_fc[2]/(r*ion_val[2]*ec))*y[1]);
	// 4th Ion
	// d(IEP)/dr
   dy[13]=y[14];
	// d2(IEP)/dr2
	A=(ion_val[3]*ec/(kb*T))*y[6];
   dy[14]=-2*y[14]/r + 2*y[13]/(r*r) + A*(y[14] - (2*ion_fc[3]/(r*ion_val[3]*ec))*y[1]);
	// 5th Ion
	// d(IEP)/dr
   dy[15]=y[16];
	// d2(IEP)/dr2
	A=(ion_val[4]*ec/(kb*T))*y[6];
   dy[16]=-2*y[16]/r + 2*y[15]/(r*r) + A*(y[16] - (2*ion_fc[4]/(r*ion_val[4]*ec))*y[1]);
	// 6th Ion
	// d(IEP)/dr
   dy[17]=y[18];
	// d2(IEP)/dr2
	A=(ion_val[5]*ec/(kb*T))*y[6];
   dy[18]=-2*y[18]/r + 2*y[17]/(r*r) + A*(y[18] - (2*ion_fc[5]/(r*ion_val[5]*ec))*y[1]);
	// 7th Ion
	// d(IEP)/dr
   dy[19]=y[20];
	// d2(IEP)/dr2
	A=(ion_val[6]*ec/(kb*T))*y[6];
   dy[20]=-2*y[20]/r + 2*y[19]/(r*r) + A*(y[20] - (2*ion_fc[6]/(r*ion_val[6]*ec))*y[1]);
	// 8th Ion
	// d(IEP)/dr
   dy[21]=y[22];
	// d2(IEP)/dr2
	A=(ion_val[7]*ec/(kb*T))*y[6];
   dy[22]=-2*y[22]/r + 2*y[21]/(r*r) + A*(y[22] - (2*ion_fc[7]/(r*ion_val[7]*ec))*y[1]);

	// Clear Memory
	delete [] ion_conc;delete [] ion_val;delete [] ion_rad;delete [] ion_mob;delete [] ion_fc;}

void inhomogeneous_Prob1_8ion(const state_type &y,state_type &dy,const double r)
	{double T=y[23];									// temperature [Kelvin]
	double viscosity=y[51];							// [Newton*second/meter^2]
	// Organize Ion Data
	int numIon=y[24];	
	double* ion_conc=new double[numIon];
	int* ion_val=new int[numIon];
	double* ion_rad=new double[numIon];
	double* ion_mob=new double[numIon];
	double* ion_fc=new double[numIon];
	// 1st Ion Properties
	ion_val[0]=y[25];
   ion_conc[0]=y[26];  		// mM or mol/m^3	    
	ion_rad[0]=y[27];		// meters
	ion_mob[0]=calc_Ion_Mobility(ion_val[0],ion_rad[0],viscosity);
	ion_fc[0]=calc_Ion_Friction_Coefficient(ion_val[0],ion_mob[0]);
	// 2nd Ion Properties
	ion_val[1]=y[28];
   ion_conc[1]=y[29];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[1]=y[30];													// [meter]
	ion_mob[1]=calc_Ion_Mobility(ion_val[1],ion_rad[1],viscosity);		// [meter^2/Volt*second]
	ion_fc[1]=calc_Ion_Friction_Coefficient(ion_val[1],ion_mob[1]);		// [Newton*second/meter]
	// 3rd Ion Properties
	ion_val[2]=y[31];
   ion_conc[2]=y[32];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[2]=y[33];													// [meter]
	ion_mob[2]=calc_Ion_Mobility(ion_val[2],ion_rad[2],viscosity);		// [meter^2/Volt*second]
	ion_fc[2]=calc_Ion_Friction_Coefficient(ion_val[2],ion_mob[2]);		// [Newton*second/meter]
	// 4th Ion Properties
	ion_val[3]=y[34];
   ion_conc[3]=y[35];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[3]=y[36];													// [meter]
	ion_mob[3]=calc_Ion_Mobility(ion_val[3],ion_rad[3],viscosity);		// [meter^2/Volt*second]
	ion_fc[3]=calc_Ion_Friction_Coefficient(ion_val[3],ion_mob[3]);		// [Newton*second/meter]
	// 5th Ion Properties
	ion_val[4]=y[37];
   ion_conc[4]=y[38];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[4]=y[39];													// [meter]
	ion_mob[4]=calc_Ion_Mobility(ion_val[4],ion_rad[4],viscosity);		// [meter^2/Volt*second]
	ion_fc[4]=calc_Ion_Friction_Coefficient(ion_val[4],ion_mob[4]);		// [Newton*second/meter]
	// 6th Ion Properties
	ion_val[5]=y[40];
   ion_conc[5]=y[41];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[5]=y[42];													// [meter]
	ion_mob[5]=calc_Ion_Mobility(ion_val[5],ion_rad[5],viscosity);		// [meter^2/Volt*second]
	ion_fc[5]=calc_Ion_Friction_Coefficient(ion_val[5],ion_mob[5]);		// [Newton*second/meter]
	// 7th Ion Properties
	ion_val[6]=y[43];
   ion_conc[6]=y[44];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[6]=y[45];													// [meter]
	ion_mob[6]=calc_Ion_Mobility(ion_val[6],ion_rad[6],viscosity);		// [meter^2/Volt*second]
	ion_fc[6]=calc_Ion_Friction_Coefficient(ion_val[6],ion_mob[6]);		// [Newton*second/meter]
	// 8th Ion Properties
	ion_val[7]=y[46];
   ion_conc[7]=y[47];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[7]=y[48];													// [meter]
	ion_mob[7]=calc_Ion_Mobility(ion_val[7],ion_rad[7],viscosity);		// [meter^2/Volt*second]
	ion_fc[7]=calc_Ion_Friction_Coefficient(ion_val[7],ion_mob[7]);		// [Newton*second/meter]

	double Debye=y[50];				// meters	
	double dielectric=y[49];
	// Convert Concentration to Number Density (particles per volume (m^3))
	ion_conc[0]*=Na;
	ion_conc[1]*=Na;
	ion_conc[2]*=Na;
	ion_conc[3]*=Na;
	ion_conc[4]*=Na;
	ion_conc[5]*=Na;
	ion_conc[6]*=Na;
	ion_conc[7]*=Na;
	// Calculation Constants
	double A,dN_dr;
	int ionIndex;
	// Linearized Governing Equations
	// F FUNCTION
	// df/dr = y[1]
	dy[0]=y[1];
	// d2f/dr2 = y[2]
	dy[1]=y[2];
	// d3f/dr3 = y[3]
	dy[2]=y[3];
	// d4f/dr4 = y[4]
	dy[3]=y[4];
	// d5f/dr5
	A=0;
	for(int i=0;i<numIon;i++)
		{// Ion Equilibrium Number Density
		dN_dr=(-ion_conc[i]*ion_val[i]*ec/(kb*T))*y[6]*exp(-ec*ion_val[i]*y[5]/(kb*T));
		// Ion Contribution
		ionIndex=7+2*i; 
		A+=ion_val[i]*ec*dN_dr*y[ionIndex];}
    dy[4]=-4*y[4]/r + 4*y[3]/(r*r) + A/(viscosity*r);
	// EQUILIBRIUM POTENTIAL
	// d(potential)/dr = y[6]
	dy[5]=y[6];
	// d2(potential)/dr2
	A=0;
	for(int i=0;i<numIon;i++){A+=ion_val[i]*ion_conc[i]*exp(-ec*ion_val[i]*y[5]/(kb*T));}
	dy[6]=-ec*A/(dielectric*eo) - 2*y[6]/r;
	// INDUCED ELECTROCHEMICAL POTENTIAL FUNCTION (2 ions)
	// 1st Ion
	// d(IEP)/dr = y[8]
   dy[7]=y[8];
	// d2(IEP)/dr2
	A=(ion_val[0]*ec/(kb*T))*y[6];
   dy[8]=-2*y[8]/r + 2*y[7]/(r*r) + A*(y[8] - (2*ion_fc[0]/(r*ion_val[0]*ec))*y[1] - ion_fc[0]/(ion_val[0]*ec));
	// 2nd Ion
	// d(IEP)/dr = y[10]
   dy[9]=y[10];
	// d2(IEP)/dr2
	A=(ion_val[1]*ec/(kb*T))*y[6];
   dy[10]=-2*y[10]/r + 2*y[9]/(r*r) + A*(y[10] - (2*ion_fc[1]/(r*ion_val[1]*ec))*y[1] - ion_fc[1]/(ion_val[1]*ec));
	// 3rd Ion
	// d(IEP)/dr
   dy[11]=y[12];
	// d2(IEP)/dr2
	A=(ion_val[2]*ec/(kb*T))*y[6];
   dy[12]=-2*y[12]/r + 2*y[11]/(r*r) + A*(y[12] - (2*ion_fc[2]/(r*ion_val[2]*ec))*y[1] - ion_fc[2]/(ion_val[2]*ec));
	// 4th Ion
	// d(IEP)/dr
   dy[13]=y[14];
	// d2(IEP)/dr2
	A=(ion_val[3]*ec/(kb*T))*y[6];
   dy[14]=-2*y[14]/r + 2*y[13]/(r*r) + A*(y[14] - (2*ion_fc[3]/(r*ion_val[3]*ec))*y[1] - ion_fc[3]/(ion_val[3]*ec));
	// 5th Ion
	// d(IEP)/dr
   dy[15]=y[16];
	// d2(IEP)/dr2
	A=(ion_val[4]*ec/(kb*T))*y[6];
   dy[16]=-2*y[16]/r + 2*y[15]/(r*r) + A*(y[16] - (2*ion_fc[4]/(r*ion_val[4]*ec))*y[1] - ion_fc[4]/(ion_val[4]*ec));
	// 6th Ion
	// d(IEP)/dr
   dy[17]=y[18];
	// d2(IEP)/dr2
	A=(ion_val[5]*ec/(kb*T))*y[6];
   dy[18]=-2*y[18]/r + 2*y[17]/(r*r) + A*(y[18] - (2*ion_fc[5]/(r*ion_val[5]*ec))*y[1] - ion_fc[5]/(ion_val[5]*ec));
	// 7th Ion
	// d(IEP)/dr
   dy[19]=y[20];
	// d2(IEP)/dr2
	A=(ion_val[6]*ec/(kb*T))*y[6];
   dy[20]=-2*y[20]/r + 2*y[19]/(r*r) + A*(y[20] - (2*ion_fc[6]/(r*ion_val[6]*ec))*y[1] - ion_fc[6]/(ion_val[6]*ec));
	// 8th Ion
	// d(IEP)/dr
   dy[21]=y[22];
	// d2(IEP)/dr2
	A=(ion_val[7]*ec/(kb*T))*y[6];
   dy[22]=-2*y[22]/r + 2*y[21]/(r*r) + A*(y[22] - (2*ion_fc[7]/(r*ion_val[7]*ec))*y[1] - ion_fc[7]/(ion_val[7]*ec));

	// Clear Memory
	delete [] ion_conc;delete [] ion_val;delete [] ion_rad;delete [] ion_mob;delete [] ion_fc;}

void inhomogeneous_Prob2_8ion(const state_type &y,state_type &dy,const double r)
	{double T=y[23];									// temperature [Kelvin]
	double viscosity=y[51];							// [Newton*second/meter^2]
	// Organize Ion Data
	int numIon=y[24];	
	double* ion_conc=new double[numIon];
	int* ion_val=new int[numIon];
	double* ion_rad=new double[numIon];
	double* ion_mob=new double[numIon];
	double* ion_fc=new double[numIon];
	// 1st Ion Properties
	ion_val[0]=y[25];
   ion_conc[0]=y[26];  		// mM or mol/m^3	    
	ion_rad[0]=y[27];		// meters
	ion_mob[0]=calc_Ion_Mobility(ion_val[0],ion_rad[0],viscosity);
	ion_fc[0]=calc_Ion_Friction_Coefficient(ion_val[0],ion_mob[0]);
	// 2nd Ion Properties
	ion_val[1]=y[28];
   ion_conc[1]=y[29];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[1]=y[30];													// [meter]
	ion_mob[1]=calc_Ion_Mobility(ion_val[1],ion_rad[1],viscosity);		// [meter^2/Volt*second]
	ion_fc[1]=calc_Ion_Friction_Coefficient(ion_val[1],ion_mob[1]);		// [Newton*second/meter]
	// 3rd Ion Properties
	ion_val[2]=y[31];
   ion_conc[2]=y[32];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[2]=y[33];													// [meter]
	ion_mob[2]=calc_Ion_Mobility(ion_val[2],ion_rad[2],viscosity);		// [meter^2/Volt*second]
	ion_fc[2]=calc_Ion_Friction_Coefficient(ion_val[2],ion_mob[2]);		// [Newton*second/meter]
	// 4th Ion Properties
	ion_val[3]=y[34];
   ion_conc[3]=y[35];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[3]=y[36];													// [meter]
	ion_mob[3]=calc_Ion_Mobility(ion_val[3],ion_rad[3],viscosity);		// [meter^2/Volt*second]
	ion_fc[3]=calc_Ion_Friction_Coefficient(ion_val[3],ion_mob[3]);		// [Newton*second/meter]
	// 5th Ion Properties
	ion_val[4]=y[37];
   ion_conc[4]=y[38];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[4]=y[39];													// [meter]
	ion_mob[4]=calc_Ion_Mobility(ion_val[4],ion_rad[4],viscosity);		// [meter^2/Volt*second]
	ion_fc[4]=calc_Ion_Friction_Coefficient(ion_val[4],ion_mob[4]);		// [Newton*second/meter]
	// 6th Ion Properties
	ion_val[5]=y[40];
   ion_conc[5]=y[41];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[5]=y[42];													// [meter]
	ion_mob[5]=calc_Ion_Mobility(ion_val[5],ion_rad[5],viscosity);		// [meter^2/Volt*second]
	ion_fc[5]=calc_Ion_Friction_Coefficient(ion_val[5],ion_mob[5]);		// [Newton*second/meter]
	// 7th Ion Properties
	ion_val[6]=y[43];
   ion_conc[6]=y[44];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[6]=y[45];													// [meter]
	ion_mob[6]=calc_Ion_Mobility(ion_val[6],ion_rad[6],viscosity);		// [meter^2/Volt*second]
	ion_fc[6]=calc_Ion_Friction_Coefficient(ion_val[6],ion_mob[6]);		// [Newton*second/meter]
	// 8th Ion Properties
	ion_val[7]=y[46];
   ion_conc[7]=y[47];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[7]=y[48];													// [meter]
	ion_mob[7]=calc_Ion_Mobility(ion_val[7],ion_rad[7],viscosity);		// [meter^2/Volt*second]
	ion_fc[7]=calc_Ion_Friction_Coefficient(ion_val[7],ion_mob[7]);		// [Newton*second/meter]

	double Debye=y[50];				// meters	
	double dielectric=y[49];
	// Convert Concentration to Number Density (particles per volume (m^3))
	ion_conc[0]*=Na;
	ion_conc[1]*=Na;
	ion_conc[2]*=Na;
	ion_conc[3]*=Na;
	ion_conc[4]*=Na;
	ion_conc[5]*=Na;
	ion_conc[6]*=Na;
	ion_conc[7]*=Na;
    // Calculation Constants
	double A,dN_dr;
	int ionIndex;
    // Linearized Governing Equations
	// F FUNCTION
	// df/dr = y[1]
	dy[0]=y[1];
	// d2f/dr2 = y[2]
	dy[1]=y[2];
	// d3f/dr3 = y[3]
	dy[2]=y[3];
	// d4f/dr4 = y[4]
	dy[3]=y[4];
	// d5f/dr5
	A=0;
	for(int i=0;i<numIon;i++)
		{// Ion Equilibrium Number Density
		dN_dr=(-ion_conc[i]*ion_val[i]*ec/(kb*T))*y[6]*exp(-ec*ion_val[i]*y[5]/(kb*T));
		// Ion Contribution
		ionIndex=7+2*i; 
		A+=ion_val[i]*ec*dN_dr*(y[ionIndex]+r);}
   dy[4]=-4*y[4]/r + 4*y[3]/(r*r) + A/(viscosity*r);
	// EQUILIBRIUM POTENTIAL
	// d(potential)/dr = y[6]
	dy[5]=y[6];
	// d2(potential)/dr2
	A=0;
	for(int i=0;i<numIon;i++){A+=ion_val[i]*ion_conc[i]*exp(-ec*ion_val[i]*y[5]/(kb*T));}
	dy[6]=-ec*A/(dielectric*eo) - 2*y[6]/r;
	// INDUCED ELECTROCHEMICAL POTENTIAL FUNCTION (2 ions)
	// 1st Ion
	// d(IEP)/dr = y[8]
   dy[7]=y[8];
	// d2(IEP)/dr2
	A=(ion_val[0]*ec/(kb*T))*y[6];
   dy[8]=-2*y[8]/r + 2*y[7]/(r*r) + A*(y[8] - (2*ion_fc[0]/(r*ion_val[0]*ec))*y[1] + 1);
	// 2nd Ion
	// d(IEP)/dr = y[10]
   dy[9]=y[10];
	// d2(IEP)/dr2
	A=(ion_val[1]*ec/(kb*T))*y[6];
   dy[10]=-2*y[10]/r + 2*y[9]/(r*r) + A*(y[10] - (2*ion_fc[1]/(r*ion_val[1]*ec))*y[1] + 1);
	// 3rd Ion
	// d(IEP)/dr
   dy[11]=y[12];
	// d2(IEP)/dr2
	A=(ion_val[2]*ec/(kb*T))*y[6];
   dy[12]=-2*y[12]/r + 2*y[11]/(r*r) + A*(y[12] - (2*ion_fc[2]/(r*ion_val[2]*ec))*y[1] + 1);
	// 4th Ion
	// d(IEP)/dr
   dy[13]=y[14];
	// d2(IEP)/dr2
	A=(ion_val[3]*ec/(kb*T))*y[6];
   dy[14]=-2*y[14]/r + 2*y[13]/(r*r) + A*(y[14] - (2*ion_fc[3]/(r*ion_val[3]*ec))*y[1] + 1);
	// 5th Ion
	// d(IEP)/dr
   dy[15]=y[16];
	// d2(IEP)/dr2
	A=(ion_val[4]*ec/(kb*T))*y[6];
   dy[16]=-2*y[16]/r + 2*y[15]/(r*r) + A*(y[16] - (2*ion_fc[4]/(r*ion_val[4]*ec))*y[1] + 1);
	// 6th Ion
	// d(IEP)/dr
   dy[17]=y[18];
	// d2(IEP)/dr2
	A=(ion_val[5]*ec/(kb*T))*y[6];
   dy[18]=-2*y[18]/r + 2*y[17]/(r*r) + A*(y[18] - (2*ion_fc[5]/(r*ion_val[5]*ec))*y[1] + 1);
	// 7th Ion
	// d(IEP)/dr
   dy[19]=y[20];
	// d2(IEP)/dr2
	A=(ion_val[6]*ec/(kb*T))*y[6];
   dy[20]=-2*y[20]/r + 2*y[19]/(r*r) + A*(y[20] - (2*ion_fc[6]/(r*ion_val[6]*ec))*y[1] + 1);
	// 8th Ion
	// d(IEP)/dr
   dy[21]=y[22];
	// d2(IEP)/dr2
	A=(ion_val[7]*ec/(kb*T))*y[6];
   dy[22]=-2*y[22]/r + 2*y[21]/(r*r) + A*(y[22] - (2*ion_fc[7]/(r*ion_val[7]*ec))*y[1] + 1);

	// Clear Memory
	delete [] ion_conc;delete [] ion_val;delete [] ion_rad;delete [] ion_mob;delete [] ion_fc;}

	// y[0] = f
	// y[1] = df/dr
	// y[2] = d2f/dr2
	// y[3] = d3f/dr3
	// y[4] = d4f/dr4
	// Equilibrium Electric Pontential
	// y[5] = potential
	// y[6] = d(potential)/dr
	// Induced Electrochemical Potential Function (3 ion)
	// y[7] = IEP1
	// y[8] = d(IEP1)/dr
	// y[9] = IEP2
	// y[10] = d(IEP2)/dr
	// y[11] = IEP3
	// y[12] = d(IEP3)/dr
	// y[13] = IEP4
	// y[14] = d(IEP4)/dr
	// y[15] = IEP5
	// y[16] = d(IEP5)/dr
	// y[17] = IEP6
	// y[18] = d(IEP6)/dr
	// y[19] = IEP7
	// y[20] = d(IEP7)/dr
	// Constants
	// y[21] = Temp					// [Kelvin]
	// y[22] = numIon
	// Ion Properties
	// y[23] = Ion1.val;		
	// y[24] = Ion1.conc;		// mM [mol/m^3]
	// y[25] = Ion1.radius;		// meters
	// y[26] = Ion2.val;
	// y[27] = Ion2.conc;
	// y[28] = Ion2.radius;
	// y[29] = Ion3.val;
	// y[30] = Ion3.conc;
	// y[31] = Ion3.radius;
	// y[32] = Ion4.val;
	// y[33] = Ion4.conc;
	// y[34] = Ion4.radius;
	// y[35] = Ion5.val;
	// y[36] = Ion5.conc;
	// y[37] = Ion5.radius;
	// y[38] = Ion6.val;
	// y[39] = Ion6.conc;
	// y[40] = Ion6.radius;
	// y[41] = Ion7.val;
	// y[42] = Ion7.conc;
	// y[43] = Ion7.radius;
	// Solution Properties
	// y[44] = relative dielectric	[dimensionless]
	// y[45] = Debye Length [meters]
	// y[46] = viscosity  [Pa s]
void homogeneous_form_7ion(const state_type &y,state_type &dy,const double r)
	{double T=y[21];									// temperature [Kelvin]
	double viscosity=y[46];							// [Newton*second/meter^2]
	// Organize Ion Data
	int numIon=y[22];	
	double* ion_conc=new double[numIon];
	int* ion_val=new int[numIon];
	double* ion_rad=new double[numIon];
	double* ion_mob=new double[numIon];
	double* ion_fc=new double[numIon];
	// 1st Ion Properties
	ion_val[0]=y[23];
   ion_conc[0]=y[24];  		// mM or mol/m^3	    
	ion_rad[0]=y[25];		// meters
	ion_mob[0]=calc_Ion_Mobility(ion_val[0],ion_rad[0],viscosity);
	ion_fc[0]=calc_Ion_Friction_Coefficient(ion_val[0],ion_mob[0]);
	// 2nd Ion Properties
	ion_val[1]=y[26];
   ion_conc[1]=y[27];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[1]=y[28];													// [meter]
	ion_mob[1]=calc_Ion_Mobility(ion_val[1],ion_rad[1],viscosity);		// [meter^2/Volt*second]
	ion_fc[1]=calc_Ion_Friction_Coefficient(ion_val[1],ion_mob[1]);		// [Newton*second/meter]
	// 3rd Ion Properties
	ion_val[2]=y[29];
   ion_conc[2]=y[30];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[2]=y[31];													// [meter]
	ion_mob[2]=calc_Ion_Mobility(ion_val[2],ion_rad[2],viscosity);		// [meter^2/Volt*second]
	ion_fc[2]=calc_Ion_Friction_Coefficient(ion_val[2],ion_mob[2]);		// [Newton*second/meter]
	// 4th Ion Properties
	ion_val[3]=y[32];
   ion_conc[3]=y[33];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[3]=y[34];													// [meter]
	ion_mob[3]=calc_Ion_Mobility(ion_val[3],ion_rad[3],viscosity);		// [meter^2/Volt*second]
	ion_fc[3]=calc_Ion_Friction_Coefficient(ion_val[3],ion_mob[3]);		// [Newton*second/meter]
	// 5th Ion Properties
	ion_val[4]=y[35];
   ion_conc[4]=y[36];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[4]=y[37];													// [meter]
	ion_mob[4]=calc_Ion_Mobility(ion_val[4],ion_rad[4],viscosity);		// [meter^2/Volt*second]
	ion_fc[4]=calc_Ion_Friction_Coefficient(ion_val[4],ion_mob[4]);		// [Newton*second/meter]
	// 6th Ion Properties
	ion_val[5]=y[38];
   ion_conc[5]=y[39];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[5]=y[40];													// [meter]
	ion_mob[5]=calc_Ion_Mobility(ion_val[5],ion_rad[5],viscosity);		// [meter^2/Volt*second]
	ion_fc[5]=calc_Ion_Friction_Coefficient(ion_val[5],ion_mob[5]);		// [Newton*second/meter]
	// 7th Ion Properties
	ion_val[6]=y[41];
   ion_conc[6]=y[42];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[6]=y[43];													// [meter]
	ion_mob[6]=calc_Ion_Mobility(ion_val[6],ion_rad[6],viscosity);		// [meter^2/Volt*second]
	ion_fc[6]=calc_Ion_Friction_Coefficient(ion_val[6],ion_mob[6]);		// [Newton*second/meter]

	double Debye=y[45];				// meters	
	double dielectric=y[44];
	// Convert Concentration to Number Density (particles per volume (m^3))
	ion_conc[0]*=Na;
	ion_conc[1]*=Na;
	ion_conc[2]*=Na;
	ion_conc[3]*=Na;
	ion_conc[4]*=Na;
	ion_conc[5]*=Na;
	ion_conc[6]*=Na;
	// Calculation Constants
	double A;
	int ionIndex;
    // Linearized Governing Equations
	// F FUNCTION
	// df/dr = y[1]
	dy[0]=y[1];
	// d2f/dr2 = y[2]
	dy[1]=y[2];
	// d3f/dr3 = y[3]
	dy[2]=y[3];
	// d4f/dr4 = y[4]
	dy[3]=y[4];
	// d5f/dr5
	// 1st Ion Equilibrium Number Density
	double dN1_dr=(-ion_conc[0]*ion_val[0]*ec/(kb*T))*y[6]*exp(-ec*ion_val[0]*y[5]/(kb*T));
	double trm1=ion_val[0]*ec*dN1_dr*y[7];
	// 2nd Ion Equilibrium Number Density
	double dN2_dr=(-ion_conc[1]*ion_val[1]*ec/(kb*T))*y[6]*exp(-ec*ion_val[1]*y[5]/(kb*T));
	double trm2=ion_val[1]*ec*dN2_dr*y[9];
	// 3rd Ion Equilibrium Number Density
	double dN3_dr=(-ion_conc[2]*ion_val[2]*ec/(kb*T))*y[6]*exp(-ec*ion_val[2]*y[5]/(kb*T));
	double trm3=ion_val[2]*ec*dN3_dr*y[11];
	// 4th Ion Equilibrium Number Density
	double dN4_dr=(-ion_conc[3]*ion_val[3]*ec/(kb*T))*y[6]*exp(-ec*ion_val[3]*y[5]/(kb*T));
	double trm4=ion_val[3]*ec*dN4_dr*y[13];
	// 5th Ion Equilibrium Number Density
	double dN5_dr=(-ion_conc[4]*ion_val[4]*ec/(kb*T))*y[6]*exp(-ec*ion_val[4]*y[5]/(kb*T));
	double trm5=ion_val[4]*ec*dN5_dr*y[15];
	// 6th Ion Equilibrium Number Density
	double dN6_dr=(-ion_conc[5]*ion_val[5]*ec/(kb*T))*y[6]*exp(-ec*ion_val[5]*y[5]/(kb*T));
	double trm6=ion_val[5]*ec*dN6_dr*y[17];
	// 7th Ion Equilibrium Number Density
	double dN7_dr=(-ion_conc[6]*ion_val[6]*ec/(kb*T))*y[6]*exp(-ec*ion_val[6]*y[5]/(kb*T));
	double trm7=ion_val[6]*ec*dN7_dr*y[19];
	//
	dy[4]=-4*y[4]/r + 4*y[3]/(r*r) + (trm1+trm2+trm3+trm4+trm5+trm6+trm7)/(viscosity*r);
	// EQUILIBRIUM POTENTIAL
	// d(potential)/dr = y[6]
	dy[5]=y[6];
	// d2(potential)/dr2
	A=0;
	for(int i=0;i<numIon;i++){A+=ion_val[i]*ion_conc[i]*exp(-ec*ion_val[i]*y[5]/(kb*T));}
	dy[6]=-ec*A/(dielectric*eo) - 2*y[6]/r;
	// INDUCED ELECTROCHEMICAL POTENTIAL FUNCTION (2 ions)
	// 1st Ion
	// d(IEP)/dr
   dy[7]=y[8];
	// d2(IEP)/dr2
	A=(ion_val[0]*ec/(kb*T))*y[6];
   dy[8]=-2*y[8]/r + 2*y[7]/(r*r) + A*(y[8] - (2*ion_fc[0]/(r*ion_val[0]*ec))*y[1]);
	// 2nd Ion
	// d(IEP)/dr
   dy[9]=y[10];
	// d2(IEP)/dr2
	A=(ion_val[1]*ec/(kb*T))*y[6];
   dy[10]=-2*y[10]/r + 2*y[9]/(r*r) + A*(y[10] - (2*ion_fc[1]/(r*ion_val[1]*ec))*y[1]);
	// 3rd Ion
	// d(IEP)/dr
   dy[11]=y[12];
	// d2(IEP)/dr2
	A=(ion_val[2]*ec/(kb*T))*y[6];
   dy[12]=-2*y[12]/r + 2*y[11]/(r*r) + A*(y[12] - (2*ion_fc[2]/(r*ion_val[2]*ec))*y[1]);
	// 4th Ion
	// d(IEP)/dr
   dy[13]=y[14];
	// d2(IEP)/dr2
	A=(ion_val[3]*ec/(kb*T))*y[6];
   dy[14]=-2*y[14]/r + 2*y[13]/(r*r) + A*(y[14] - (2*ion_fc[3]/(r*ion_val[3]*ec))*y[1]);
	// 5th Ion
	// d(IEP)/dr
   dy[15]=y[16];
	// d2(IEP)/dr2
	A=(ion_val[4]*ec/(kb*T))*y[6];
   dy[16]=-2*y[16]/r + 2*y[15]/(r*r) + A*(y[16] - (2*ion_fc[4]/(r*ion_val[4]*ec))*y[1]);
	// 6th Ion
	// d(IEP)/dr
   dy[17]=y[18];
	// d2(IEP)/dr2
	A=(ion_val[5]*ec/(kb*T))*y[6];
   dy[18]=-2*y[18]/r + 2*y[17]/(r*r) + A*(y[18] - (2*ion_fc[5]/(r*ion_val[5]*ec))*y[1]);
	// 7th Ion
	// d(IEP)/dr
   dy[19]=y[20];
	// d2(IEP)/dr2
	A=(ion_val[6]*ec/(kb*T))*y[6];
   dy[20]=-2*y[20]/r + 2*y[19]/(r*r) + A*(y[20] - (2*ion_fc[6]/(r*ion_val[6]*ec))*y[1]);

	// Clear Memory
	delete [] ion_conc;delete [] ion_val;delete [] ion_rad;delete [] ion_mob;delete [] ion_fc;}

void inhomogeneous_Prob1_7ion(const state_type &y,state_type &dy,const double r)
	{double T=y[21];									// temperature [Kelvin]
	double viscosity=y[46];							// [Newton*second/meter^2]
	// Organize Ion Data
	int numIon=y[22];	
	double* ion_conc=new double[numIon];
	int* ion_val=new int[numIon];
	double* ion_rad=new double[numIon];
	double* ion_mob=new double[numIon];
	double* ion_fc=new double[numIon];
	// 1st Ion Properties
	ion_val[0]=y[23];
   ion_conc[0]=y[24];  		// mM or mol/m^3	    
	ion_rad[0]=y[25];		// meters
	ion_mob[0]=calc_Ion_Mobility(ion_val[0],ion_rad[0],viscosity);
	ion_fc[0]=calc_Ion_Friction_Coefficient(ion_val[0],ion_mob[0]);
	// 2nd Ion Properties
	ion_val[1]=y[26];
   ion_conc[1]=y[27];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[1]=y[28];													// [meter]
	ion_mob[1]=calc_Ion_Mobility(ion_val[1],ion_rad[1],viscosity);		// [meter^2/Volt*second]
	ion_fc[1]=calc_Ion_Friction_Coefficient(ion_val[1],ion_mob[1]);		// [Newton*second/meter]
	// 3rd Ion Properties
	ion_val[2]=y[29];
   ion_conc[2]=y[30];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[2]=y[31];													// [meter]
	ion_mob[2]=calc_Ion_Mobility(ion_val[2],ion_rad[2],viscosity);		// [meter^2/Volt*second]
	ion_fc[2]=calc_Ion_Friction_Coefficient(ion_val[2],ion_mob[2]);		// [Newton*second/meter]
	// 4th Ion Properties
	ion_val[3]=y[32];
   ion_conc[3]=y[33];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[3]=y[34];													// [meter]
	ion_mob[3]=calc_Ion_Mobility(ion_val[3],ion_rad[3],viscosity);		// [meter^2/Volt*second]
	ion_fc[3]=calc_Ion_Friction_Coefficient(ion_val[3],ion_mob[3]);		// [Newton*second/meter]
	// 5th Ion Properties
	ion_val[4]=y[35];
   ion_conc[4]=y[36];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[4]=y[37];													// [meter]
	ion_mob[4]=calc_Ion_Mobility(ion_val[4],ion_rad[4],viscosity);		// [meter^2/Volt*second]
	ion_fc[4]=calc_Ion_Friction_Coefficient(ion_val[4],ion_mob[4]);		// [Newton*second/meter]
	// 6th Ion Properties
	ion_val[5]=y[38];
   ion_conc[5]=y[39];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[5]=y[40];													// [meter]
	ion_mob[5]=calc_Ion_Mobility(ion_val[5],ion_rad[5],viscosity);		// [meter^2/Volt*second]
	ion_fc[5]=calc_Ion_Friction_Coefficient(ion_val[5],ion_mob[5]);		// [Newton*second/meter]
	// 7th Ion Properties
	ion_val[6]=y[41];
   ion_conc[6]=y[42];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[6]=y[43];													// [meter]
	ion_mob[6]=calc_Ion_Mobility(ion_val[6],ion_rad[6],viscosity);		// [meter^2/Volt*second]
	ion_fc[6]=calc_Ion_Friction_Coefficient(ion_val[6],ion_mob[6]);		// [Newton*second/meter]

	double Debye=y[45];				// meters	
	double dielectric=y[44];
	// Convert Concentration to Number Density (particles per volume (m^3))
	ion_conc[0]*=Na;
	ion_conc[1]*=Na;
	ion_conc[2]*=Na;
	ion_conc[3]*=Na;
	ion_conc[4]*=Na;
	ion_conc[5]*=Na;
	ion_conc[6]*=Na;
	// Calculation Constants
	double A,dN_dr;
	int ionIndex;
	// Linearized Governing Equations
	// F FUNCTION
	// df/dr = y[1]
	dy[0]=y[1];
	// d2f/dr2 = y[2]
	dy[1]=y[2];
	// d3f/dr3 = y[3]
	dy[2]=y[3];
	// d4f/dr4 = y[4]
	dy[3]=y[4];
	// d5f/dr5
	A=0;
	for(int i=0;i<numIon;i++)
		{// Ion Equilibrium Number Density
		dN_dr=(-ion_conc[i]*ion_val[i]*ec/(kb*T))*y[6]*exp(-ec*ion_val[i]*y[5]/(kb*T));
		// Ion Contribution
		ionIndex=7+2*i; 
		A+=ion_val[i]*ec*dN_dr*y[ionIndex];}
    dy[4]=-4*y[4]/r + 4*y[3]/(r*r) + A/(viscosity*r);
	// EQUILIBRIUM POTENTIAL
	// d(potential)/dr = y[6]
	dy[5]=y[6];
	// d2(potential)/dr2
	A=0;
	for(int i=0;i<numIon;i++){A+=ion_val[i]*ion_conc[i]*exp(-ec*ion_val[i]*y[5]/(kb*T));}
	dy[6]=-ec*A/(dielectric*eo) - 2*y[6]/r;
	// INDUCED ELECTROCHEMICAL POTENTIAL FUNCTION (2 ions)
	// 1st Ion
	// d(IEP)/dr = y[8]
   dy[7]=y[8];
	// d2(IEP)/dr2
	A=(ion_val[0]*ec/(kb*T))*y[6];
   dy[8]=-2*y[8]/r + 2*y[7]/(r*r) + A*(y[8] - (2*ion_fc[0]/(r*ion_val[0]*ec))*y[1] - ion_fc[0]/(ion_val[0]*ec));
	// 2nd Ion
	// d(IEP)/dr = y[10]
   dy[9]=y[10];
	// d2(IEP)/dr2
	A=(ion_val[1]*ec/(kb*T))*y[6];
   dy[10]=-2*y[10]/r + 2*y[9]/(r*r) + A*(y[10] - (2*ion_fc[1]/(r*ion_val[1]*ec))*y[1] - ion_fc[1]/(ion_val[1]*ec));
	// 3rd Ion
	// d(IEP)/dr
   dy[11]=y[12];
	// d2(IEP)/dr2
	A=(ion_val[2]*ec/(kb*T))*y[6];
   dy[12]=-2*y[12]/r + 2*y[11]/(r*r) + A*(y[12] - (2*ion_fc[2]/(r*ion_val[2]*ec))*y[1] - ion_fc[2]/(ion_val[2]*ec));
	// 4th Ion
	// d(IEP)/dr
   dy[13]=y[14];
	// d2(IEP)/dr2
	A=(ion_val[3]*ec/(kb*T))*y[6];
   dy[14]=-2*y[14]/r + 2*y[13]/(r*r) + A*(y[14] - (2*ion_fc[3]/(r*ion_val[3]*ec))*y[1] - ion_fc[3]/(ion_val[3]*ec));
	// 5th Ion
	// d(IEP)/dr
   dy[15]=y[16];
	// d2(IEP)/dr2
	A=(ion_val[4]*ec/(kb*T))*y[6];
   dy[16]=-2*y[16]/r + 2*y[15]/(r*r) + A*(y[16] - (2*ion_fc[4]/(r*ion_val[4]*ec))*y[1] - ion_fc[4]/(ion_val[4]*ec));
	// 6th Ion
	// d(IEP)/dr
   dy[17]=y[18];
	// d2(IEP)/dr2
	A=(ion_val[5]*ec/(kb*T))*y[6];
   dy[18]=-2*y[18]/r + 2*y[17]/(r*r) + A*(y[18] - (2*ion_fc[5]/(r*ion_val[5]*ec))*y[1] - ion_fc[5]/(ion_val[5]*ec));
	// 7th Ion
	// d(IEP)/dr
   dy[19]=y[20];
	// d2(IEP)/dr2
	A=(ion_val[6]*ec/(kb*T))*y[6];
   dy[20]=-2*y[20]/r + 2*y[19]/(r*r) + A*(y[20] - (2*ion_fc[6]/(r*ion_val[6]*ec))*y[1] - ion_fc[6]/(ion_val[6]*ec));

	// Clear Memory
	delete [] ion_conc;delete [] ion_val;delete [] ion_rad;delete [] ion_mob;delete [] ion_fc;}

void inhomogeneous_Prob2_7ion(const state_type &y,state_type &dy,const double r)
	{double T=y[21];									// temperature [Kelvin]
	double viscosity=y[46];							// [Newton*second/meter^2]
	// Organize Ion Data
	int numIon=y[22];	
	double* ion_conc=new double[numIon];
	int* ion_val=new int[numIon];
	double* ion_rad=new double[numIon];
	double* ion_mob=new double[numIon];
	double* ion_fc=new double[numIon];
	// 1st Ion Properties
	ion_val[0]=y[23];
   ion_conc[0]=y[24];  		// mM or mol/m^3	    
	ion_rad[0]=y[25];		// meters
	ion_mob[0]=calc_Ion_Mobility(ion_val[0],ion_rad[0],viscosity);
	ion_fc[0]=calc_Ion_Friction_Coefficient(ion_val[0],ion_mob[0]);
	// 2nd Ion Properties
	ion_val[1]=y[26];
   ion_conc[1]=y[27];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[1]=y[28];													// [meter]
	ion_mob[1]=calc_Ion_Mobility(ion_val[1],ion_rad[1],viscosity);		// [meter^2/Volt*second]
	ion_fc[1]=calc_Ion_Friction_Coefficient(ion_val[1],ion_mob[1]);		// [Newton*second/meter]
	// 3rd Ion Properties
	ion_val[2]=y[29];
   ion_conc[2]=y[30];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[2]=y[31];													// [meter]
	ion_mob[2]=calc_Ion_Mobility(ion_val[2],ion_rad[2],viscosity);		// [meter^2/Volt*second]
	ion_fc[2]=calc_Ion_Friction_Coefficient(ion_val[2],ion_mob[2]);		// [Newton*second/meter]
	// 4th Ion Properties
	ion_val[3]=y[32];
   ion_conc[3]=y[33];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[3]=y[34];													// [meter]
	ion_mob[3]=calc_Ion_Mobility(ion_val[3],ion_rad[3],viscosity);		// [meter^2/Volt*second]
	ion_fc[3]=calc_Ion_Friction_Coefficient(ion_val[3],ion_mob[3]);		// [Newton*second/meter]
	// 5th Ion Properties
	ion_val[4]=y[35];
   ion_conc[4]=y[36];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[4]=y[37];													// [meter]
	ion_mob[4]=calc_Ion_Mobility(ion_val[4],ion_rad[4],viscosity);		// [meter^2/Volt*second]
	ion_fc[4]=calc_Ion_Friction_Coefficient(ion_val[4],ion_mob[4]);		// [Newton*second/meter]
	// 6th Ion Properties
	ion_val[5]=y[38];
   ion_conc[5]=y[39];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[5]=y[40];													// [meter]
	ion_mob[5]=calc_Ion_Mobility(ion_val[5],ion_rad[5],viscosity);		// [meter^2/Volt*second]
	ion_fc[5]=calc_Ion_Friction_Coefficient(ion_val[5],ion_mob[5]);		// [Newton*second/meter]
	// 7th Ion Properties
	ion_val[6]=y[41];
   ion_conc[6]=y[42];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[6]=y[43];													// [meter]
	ion_mob[6]=calc_Ion_Mobility(ion_val[6],ion_rad[6],viscosity);		// [meter^2/Volt*second]
	ion_fc[6]=calc_Ion_Friction_Coefficient(ion_val[6],ion_mob[6]);		// [Newton*second/meter]

	double Debye=y[45];				// meters	
	double dielectric=y[44];
	// Convert Concentration to Number Density (particles per volume (m^3))
	ion_conc[0]*=Na;
	ion_conc[1]*=Na;
	ion_conc[2]*=Na;
	ion_conc[3]*=Na;
	ion_conc[4]*=Na;
	ion_conc[5]*=Na;
	ion_conc[6]*=Na;
    // Calculation Constants
	double A,dN_dr;
	int ionIndex;
    // Linearized Governing Equations
	// F FUNCTION
	// df/dr = y[1]
	dy[0]=y[1];
	// d2f/dr2 = y[2]
	dy[1]=y[2];
	// d3f/dr3 = y[3]
	dy[2]=y[3];
	// d4f/dr4 = y[4]
	dy[3]=y[4];
	// d5f/dr5
	A=0;
	for(int i=0;i<numIon;i++)
		{// Ion Equilibrium Number Density
		dN_dr=(-ion_conc[i]*ion_val[i]*ec/(kb*T))*y[6]*exp(-ec*ion_val[i]*y[5]/(kb*T));
		// Ion Contribution
		ionIndex=7+2*i; 
		A+=ion_val[i]*ec*dN_dr*(y[ionIndex]+r);}
   dy[4]=-4*y[4]/r + 4*y[3]/(r*r) + A/(viscosity*r);
	// EQUILIBRIUM POTENTIAL
	// d(potential)/dr = y[6]
	dy[5]=y[6];
	// d2(potential)/dr2
	A=0;
	for(int i=0;i<numIon;i++){A+=ion_val[i]*ion_conc[i]*exp(-ec*ion_val[i]*y[5]/(kb*T));}
	dy[6]=-ec*A/(dielectric*eo) - 2*y[6]/r;
	// INDUCED ELECTROCHEMICAL POTENTIAL FUNCTION (2 ions)
	// 1st Ion
	// d(IEP)/dr = y[8]
   dy[7]=y[8];
	// d2(IEP)/dr2
	A=(ion_val[0]*ec/(kb*T))*y[6];
   dy[8]=-2*y[8]/r + 2*y[7]/(r*r) + A*(y[8] - (2*ion_fc[0]/(r*ion_val[0]*ec))*y[1] + 1);
	// 2nd Ion
	// d(IEP)/dr = y[10]
   dy[9]=y[10];
	// d2(IEP)/dr2
	A=(ion_val[1]*ec/(kb*T))*y[6];
   dy[10]=-2*y[10]/r + 2*y[9]/(r*r) + A*(y[10] - (2*ion_fc[1]/(r*ion_val[1]*ec))*y[1] + 1);
	// 3rd Ion
	// d(IEP)/dr
   dy[11]=y[12];
	// d2(IEP)/dr2
	A=(ion_val[2]*ec/(kb*T))*y[6];
   dy[12]=-2*y[12]/r + 2*y[11]/(r*r) + A*(y[12] - (2*ion_fc[2]/(r*ion_val[2]*ec))*y[1] + 1);
	// 4th Ion
	// d(IEP)/dr
   dy[13]=y[14];
	// d2(IEP)/dr2
	A=(ion_val[3]*ec/(kb*T))*y[6];
   dy[14]=-2*y[14]/r + 2*y[13]/(r*r) + A*(y[14] - (2*ion_fc[3]/(r*ion_val[3]*ec))*y[1] + 1);
	// 5th Ion
	// d(IEP)/dr
   dy[15]=y[16];
	// d2(IEP)/dr2
	A=(ion_val[4]*ec/(kb*T))*y[6];
   dy[16]=-2*y[16]/r + 2*y[15]/(r*r) + A*(y[16] - (2*ion_fc[4]/(r*ion_val[4]*ec))*y[1] + 1);
	// 6th Ion
	// d(IEP)/dr
   dy[17]=y[18];
	// d2(IEP)/dr2
	A=(ion_val[5]*ec/(kb*T))*y[6];
   dy[18]=-2*y[18]/r + 2*y[17]/(r*r) + A*(y[18] - (2*ion_fc[5]/(r*ion_val[5]*ec))*y[1] + 1);
	// 7th Ion
	// d(IEP)/dr
   dy[19]=y[20];
	// d2(IEP)/dr2
	A=(ion_val[6]*ec/(kb*T))*y[6];
   dy[20]=-2*y[20]/r + 2*y[19]/(r*r) + A*(y[20] - (2*ion_fc[6]/(r*ion_val[6]*ec))*y[1] + 1);

	// Clear Memory
	delete [] ion_conc;delete [] ion_val;delete [] ion_rad;delete [] ion_mob;delete [] ion_fc;}

	// y[0] = f
	// y[1] = df/dr
	// y[2] = d2f/dr2
	// y[3] = d3f/dr3
	// y[4] = d4f/dr4
	// Equilibrium Electric Pontential
	// y[5] = potential
	// y[6] = d(potential)/dr
	// Induced Electrochemical Potential Function (3 ion)
	// y[7] = IEP1
	// y[8] = d(IEP1)/dr
	// y[9] = IEP2
	// y[10] = d(IEP2)/dr
	// y[11] = IEP3
	// y[12] = d(IEP3)/dr
	// y[13] = IEP4
	// y[14] = d(IEP4)/dr
	// y[15] = IEP5
	// y[16] = d(IEP5)/dr
	// y[17] = IEP6
	// y[18] = d(IEP6)/dr
	// Constants
	// y[19] = Temp					// [Kelvin]
	// y[20] = numIon
	// Ion Properties
	// y[21] = Ion1.val;		
	// y[22] = Ion1.conc;		// mM [mol/m^3]
	// y[23] = Ion1.radius;		// meters
	// y[24] = Ion2.val;
	// y[25] = Ion2.conc;
	// y[26] = Ion2.radius;
	// y[27] = Ion3.val;
	// y[28] = Ion3.conc;
	// y[29] = Ion3.radius;
	// y[30] = Ion4.val;
	// y[31] = Ion4.conc;
	// y[32] = Ion4.radius;
	// y[33] = Ion5.val;
	// y[34] = Ion5.conc;
	// y[35] = Ion5.radius;
	// y[36] = Ion6.val;
	// y[37] = Ion6.conc;
	// y[38] = Ion6.radius;
	// Solution Properties
	// y[39] = relative dielectric	[dimensionless]
	// y[40] = Debye Length [meters]
	// y[41] = viscosity  [Pa s]
void homogeneous_form_6ion(const state_type &y,state_type &dy,const double r)
	{double T=y[19];									// temperature [Kelvin]
	double viscosity=y[41];							// [Newton*second/meter^2]
	// Organize Ion Data
	int numIon=y[20];	
	double* ion_conc=new double[numIon];
	int* ion_val=new int[numIon];
	double* ion_rad=new double[numIon];
	double* ion_mob=new double[numIon];
	double* ion_fc=new double[numIon];
	// 1st Ion Properties
	ion_val[0]=y[21];
   ion_conc[0]=y[22];  		// mM or mol/m^3	    
	ion_rad[0]=y[23];		// meters
	ion_mob[0]=calc_Ion_Mobility(ion_val[0],ion_rad[0],viscosity);
	ion_fc[0]=calc_Ion_Friction_Coefficient(ion_val[0],ion_mob[0]);
	// 2nd Ion Properties
	ion_val[1]=y[24];
   ion_conc[1]=y[25];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[1]=y[26];													// [meter]
	ion_mob[1]=calc_Ion_Mobility(ion_val[1],ion_rad[1],viscosity);		// [meter^2/Volt*second]
	ion_fc[1]=calc_Ion_Friction_Coefficient(ion_val[1],ion_mob[1]);		// [Newton*second/meter]
	// 3rd Ion Properties
	ion_val[2]=y[27];
   ion_conc[2]=y[28];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[2]=y[29];													// [meter]
	ion_mob[2]=calc_Ion_Mobility(ion_val[2],ion_rad[2],viscosity);		// [meter^2/Volt*second]
	ion_fc[2]=calc_Ion_Friction_Coefficient(ion_val[2],ion_mob[2]);		// [Newton*second/meter]
	// 4th Ion Properties
	ion_val[3]=y[30];
   ion_conc[3]=y[31];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[3]=y[32];													// [meter]
	ion_mob[3]=calc_Ion_Mobility(ion_val[3],ion_rad[3],viscosity);		// [meter^2/Volt*second]
	ion_fc[3]=calc_Ion_Friction_Coefficient(ion_val[3],ion_mob[3]);		// [Newton*second/meter]
	// 5th Ion Properties
	ion_val[4]=y[33];
   ion_conc[4]=y[34];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[4]=y[35];													// [meter]
	ion_mob[4]=calc_Ion_Mobility(ion_val[4],ion_rad[4],viscosity);		// [meter^2/Volt*second]
	ion_fc[4]=calc_Ion_Friction_Coefficient(ion_val[4],ion_mob[4]);		// [Newton*second/meter]
	// 6th Ion Properties
	ion_val[5]=y[36];
   ion_conc[5]=y[37];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[5]=y[38];													// [meter]
	ion_mob[5]=calc_Ion_Mobility(ion_val[5],ion_rad[5],viscosity);		// [meter^2/Volt*second]
	ion_fc[5]=calc_Ion_Friction_Coefficient(ion_val[5],ion_mob[5]);		// [Newton*second/meter]

	double Debye=y[40];				// meters	
	double dielectric=y[39];
	// Convert Concentration to Number Density (particles per volume (m^3))
	ion_conc[0]*=Na;
	ion_conc[1]*=Na;
	ion_conc[2]*=Na;
	ion_conc[3]*=Na;
	ion_conc[4]*=Na;
	ion_conc[5]*=Na;
	// Calculation Constants
	double A;
	int ionIndex;
    // Linearized Governing Equations
	// F FUNCTION
	// df/dr = y[1]
	dy[0]=y[1];
	// d2f/dr2 = y[2]
	dy[1]=y[2];
	// d3f/dr3 = y[3]
	dy[2]=y[3];
	// d4f/dr4 = y[4]
	dy[3]=y[4];
	// d5f/dr5
	// 1st Ion Equilibrium Number Density
	double dN1_dr=(-ion_conc[0]*ion_val[0]*ec/(kb*T))*y[6]*exp(-ec*ion_val[0]*y[5]/(kb*T));
	double trm1=ion_val[0]*ec*dN1_dr*y[7];
	// 2nd Ion Equilibrium Number Density
	double dN2_dr=(-ion_conc[1]*ion_val[1]*ec/(kb*T))*y[6]*exp(-ec*ion_val[1]*y[5]/(kb*T));
	double trm2=ion_val[1]*ec*dN2_dr*y[9];
	// 3rd Ion Equilibrium Number Density
	double dN3_dr=(-ion_conc[2]*ion_val[2]*ec/(kb*T))*y[6]*exp(-ec*ion_val[2]*y[5]/(kb*T));
	double trm3=ion_val[2]*ec*dN3_dr*y[11];
	// 4th Ion Equilibrium Number Density
	double dN4_dr=(-ion_conc[3]*ion_val[3]*ec/(kb*T))*y[6]*exp(-ec*ion_val[3]*y[5]/(kb*T));
	double trm4=ion_val[3]*ec*dN4_dr*y[13];
	// 5th Ion Equilibrium Number Density
	double dN5_dr=(-ion_conc[4]*ion_val[4]*ec/(kb*T))*y[6]*exp(-ec*ion_val[4]*y[5]/(kb*T));
	double trm5=ion_val[4]*ec*dN5_dr*y[15];
	// 6th Ion Equilibrium Number Density
	double dN6_dr=(-ion_conc[5]*ion_val[5]*ec/(kb*T))*y[6]*exp(-ec*ion_val[5]*y[5]/(kb*T));
	double trm6=ion_val[5]*ec*dN6_dr*y[17];
	//
	dy[4]=-4*y[4]/r + 4*y[3]/(r*r) + (trm1+trm2+trm3+trm4+trm5+trm6)/(viscosity*r);
	// EQUILIBRIUM POTENTIAL
	// d(potential)/dr = y[6]
	dy[5]=y[6];
	// d2(potential)/dr2
	A=0;
	for(int i=0;i<numIon;i++){A+=ion_val[i]*ion_conc[i]*exp(-ec*ion_val[i]*y[5]/(kb*T));}
	dy[6]=-ec*A/(dielectric*eo) - 2*y[6]/r;
	// INDUCED ELECTROCHEMICAL POTENTIAL FUNCTION (2 ions)
	// 1st Ion
	// d(IEP)/dr
   dy[7]=y[8];
	// d2(IEP)/dr2
	A=(ion_val[0]*ec/(kb*T))*y[6];
   dy[8]=-2*y[8]/r + 2*y[7]/(r*r) + A*(y[8] - (2*ion_fc[0]/(r*ion_val[0]*ec))*y[1]);
	// 2nd Ion
	// d(IEP)/dr
   dy[9]=y[10];
	// d2(IEP)/dr2
	A=(ion_val[1]*ec/(kb*T))*y[6];
   dy[10]=-2*y[10]/r + 2*y[9]/(r*r) + A*(y[10] - (2*ion_fc[1]/(r*ion_val[1]*ec))*y[1]);
	// 3rd Ion
	// d(IEP)/dr
   dy[11]=y[12];
	// d2(IEP)/dr2
	A=(ion_val[2]*ec/(kb*T))*y[6];
   dy[12]=-2*y[12]/r + 2*y[11]/(r*r) + A*(y[12] - (2*ion_fc[2]/(r*ion_val[2]*ec))*y[1]);
	// 4th Ion
	// d(IEP)/dr
   dy[13]=y[14];
	// d2(IEP)/dr2
	A=(ion_val[3]*ec/(kb*T))*y[6];
   dy[14]=-2*y[14]/r + 2*y[13]/(r*r) + A*(y[14] - (2*ion_fc[3]/(r*ion_val[3]*ec))*y[1]);
	// 5th Ion
	// d(IEP)/dr
   dy[15]=y[16];
	// d2(IEP)/dr2
	A=(ion_val[4]*ec/(kb*T))*y[6];
   dy[16]=-2*y[16]/r + 2*y[15]/(r*r) + A*(y[16] - (2*ion_fc[4]/(r*ion_val[4]*ec))*y[1]);
	// 6th Ion
	// d(IEP)/dr
   dy[17]=y[18];
	// d2(IEP)/dr2
	A=(ion_val[5]*ec/(kb*T))*y[6];
   dy[18]=-2*y[18]/r + 2*y[17]/(r*r) + A*(y[18] - (2*ion_fc[5]/(r*ion_val[5]*ec))*y[1]);

	// Clear Memory
	delete [] ion_conc;delete [] ion_val;delete [] ion_rad;delete [] ion_mob;delete [] ion_fc;}

void inhomogeneous_Prob1_6ion(const state_type &y,state_type &dy,const double r)
	{double T=y[19];									// temperature [Kelvin]
	double viscosity=y[41];							// [Newton*second/meter^2]
	// Organize Ion Data
	int numIon=y[20];	
	double* ion_conc=new double[numIon];
	int* ion_val=new int[numIon];
	double* ion_rad=new double[numIon];
	double* ion_mob=new double[numIon];
	double* ion_fc=new double[numIon];
	// 1st Ion Properties
	ion_val[0]=y[21];
   ion_conc[0]=y[22];  		// mM or mol/m^3	    
	ion_rad[0]=y[23];		// meters
	ion_mob[0]=calc_Ion_Mobility(ion_val[0],ion_rad[0],viscosity);
	ion_fc[0]=calc_Ion_Friction_Coefficient(ion_val[0],ion_mob[0]);
	// 2nd Ion Properties
	ion_val[1]=y[24];
   ion_conc[1]=y[25];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[1]=y[26];													// [meter]
	ion_mob[1]=calc_Ion_Mobility(ion_val[1],ion_rad[1],viscosity);		// [meter^2/Volt*second]
	ion_fc[1]=calc_Ion_Friction_Coefficient(ion_val[1],ion_mob[1]);		// [Newton*second/meter]
	// 3rd Ion Properties
	ion_val[2]=y[27];
   ion_conc[2]=y[28];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[2]=y[29];													// [meter]
	ion_mob[2]=calc_Ion_Mobility(ion_val[2],ion_rad[2],viscosity);		// [meter^2/Volt*second]
	ion_fc[2]=calc_Ion_Friction_Coefficient(ion_val[2],ion_mob[2]);		// [Newton*second/meter]
	// 4th Ion Properties
	ion_val[3]=y[30];
   ion_conc[3]=y[31];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[3]=y[32];													// [meter]
	ion_mob[3]=calc_Ion_Mobility(ion_val[3],ion_rad[3],viscosity);		// [meter^2/Volt*second]
	ion_fc[3]=calc_Ion_Friction_Coefficient(ion_val[3],ion_mob[3]);		// [Newton*second/meter]
	// 5th Ion Properties
	ion_val[4]=y[33];
   ion_conc[4]=y[34];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[4]=y[35];													// [meter]
	ion_mob[4]=calc_Ion_Mobility(ion_val[4],ion_rad[4],viscosity);		// [meter^2/Volt*second]
	ion_fc[4]=calc_Ion_Friction_Coefficient(ion_val[4],ion_mob[4]);		// [Newton*second/meter]
	// 6th Ion Properties
	ion_val[5]=y[36];
   ion_conc[5]=y[37];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[5]=y[38];													// [meter]
	ion_mob[5]=calc_Ion_Mobility(ion_val[5],ion_rad[5],viscosity);		// [meter^2/Volt*second]
	ion_fc[5]=calc_Ion_Friction_Coefficient(ion_val[5],ion_mob[5]);		// [Newton*second/meter]

	double Debye=y[40];	// meters	
	double dielectric=y[39];
	// Convert Concentration to Number Density (particles per volume (m^3))
	ion_conc[0]*=Na;
	ion_conc[1]*=Na;
   ion_conc[2]*=Na;
	ion_conc[3]*=Na;
	ion_conc[4]*=Na;
	ion_conc[5]*=Na;
	// Calculation Constants
	double A,dN_dr;
	int ionIndex;
	// Linearized Governing Equations
	// F FUNCTION
	// df/dr = y[1]
	dy[0]=y[1];
	// d2f/dr2 = y[2]
	dy[1]=y[2];
	// d3f/dr3 = y[3]
	dy[2]=y[3];
	// d4f/dr4 = y[4]
	dy[3]=y[4];
	// d5f/dr5
	A=0;
	for(int i=0;i<numIon;i++)
		{// Ion Equilibrium Number Density
		dN_dr=(-ion_conc[i]*ion_val[i]*ec/(kb*T))*y[6]*exp(-ec*ion_val[i]*y[5]/(kb*T));
		// Ion Contribution
		ionIndex=7+2*i; 
		A+=ion_val[i]*ec*dN_dr*y[ionIndex];}
    dy[4]=-4*y[4]/r + 4*y[3]/(r*r) + A/(viscosity*r);
	// EQUILIBRIUM POTENTIAL
	// d(potential)/dr = y[6]
	dy[5]=y[6];
	// d2(potential)/dr2
	A=0;
	for(int i=0;i<numIon;i++){A+=ion_val[i]*ion_conc[i]*exp(-ec*ion_val[i]*y[5]/(kb*T));}
	dy[6]=-ec*A/(dielectric*eo) - 2*y[6]/r;
	// INDUCED ELECTROCHEMICAL POTENTIAL FUNCTION (2 ions)
	// 1st Ion
	// d(IEP)/dr = y[8]
   dy[7]=y[8];
	// d2(IEP)/dr2
	A=(ion_val[0]*ec/(kb*T))*y[6];
   dy[8]=-2*y[8]/r + 2*y[7]/(r*r) + A*(y[8] - (2*ion_fc[0]/(r*ion_val[0]*ec))*y[1] - ion_fc[0]/(ion_val[0]*ec));
	// 2nd Ion
	// d(IEP)/dr = y[10]
   dy[9]=y[10];
	// d2(IEP)/dr2
	A=(ion_val[1]*ec/(kb*T))*y[6];
   dy[10]=-2*y[10]/r + 2*y[9]/(r*r) + A*(y[10] - (2*ion_fc[1]/(r*ion_val[1]*ec))*y[1] - ion_fc[1]/(ion_val[1]*ec));
	// 3rd Ion
	// d(IEP)/dr
   dy[11]=y[12];
	// d2(IEP)/dr2
	A=(ion_val[2]*ec/(kb*T))*y[6];
   dy[12]=-2*y[12]/r + 2*y[11]/(r*r) + A*(y[12] - (2*ion_fc[2]/(r*ion_val[2]*ec))*y[1] - ion_fc[2]/(ion_val[2]*ec));
	// 4th Ion
	// d(IEP)/dr
   dy[13]=y[14];
	// d2(IEP)/dr2
	A=(ion_val[3]*ec/(kb*T))*y[6];
   dy[14]=-2*y[14]/r + 2*y[13]/(r*r) + A*(y[14] - (2*ion_fc[3]/(r*ion_val[3]*ec))*y[1] - ion_fc[3]/(ion_val[3]*ec));
	// 5th Ion
	// d(IEP)/dr
   dy[15]=y[16];
	// d2(IEP)/dr2
	A=(ion_val[4]*ec/(kb*T))*y[6];
   dy[16]=-2*y[16]/r + 2*y[15]/(r*r) + A*(y[16] - (2*ion_fc[4]/(r*ion_val[4]*ec))*y[1] - ion_fc[4]/(ion_val[4]*ec));
	// 6th Ion
	// d(IEP)/dr
   dy[17]=y[18];
	// d2(IEP)/dr2
	A=(ion_val[5]*ec/(kb*T))*y[6];
   dy[18]=-2*y[18]/r + 2*y[17]/(r*r) + A*(y[18] - (2*ion_fc[5]/(r*ion_val[5]*ec))*y[1] - ion_fc[5]/(ion_val[5]*ec));

	// Clear Memory
	delete [] ion_conc;delete [] ion_val;delete [] ion_rad;delete [] ion_mob;delete [] ion_fc;}

void inhomogeneous_Prob2_6ion(const state_type &y,state_type &dy,const double r)
	{double T=y[19];								// temperature [Kelvin]
	double viscosity=y[41];						// [Newton*second/meter^2]
	// Organize Ion Data
	int numIon=y[20];	
	double* ion_conc=new double[numIon];
	int* ion_val=new int[numIon];
	double* ion_rad=new double[numIon];
	double* ion_mob=new double[numIon];
	double* ion_fc=new double[numIon];
	// 1st Ion Properties
	ion_val[0]=y[21];
   ion_conc[0]=y[22];  		// mM or mol/m^3	    
	ion_rad[0]=y[23];		// meters
	ion_mob[0]=calc_Ion_Mobility(ion_val[0],ion_rad[0],viscosity);
	ion_fc[0]=calc_Ion_Friction_Coefficient(ion_val[0],ion_mob[0]);
	// 2nd Ion Properties
	ion_val[1]=y[24];
   ion_conc[1]=y[25];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[1]=y[26];													// [meter]
	ion_mob[1]=calc_Ion_Mobility(ion_val[1],ion_rad[1],viscosity);		// [meter^2/Volt*second]
	ion_fc[1]=calc_Ion_Friction_Coefficient(ion_val[1],ion_mob[1]);		// [Newton*second/meter]
	// 3rd Ion Properties
	ion_val[2]=y[27];
   ion_conc[2]=y[28];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[2]=y[29];													// [meter]
	ion_mob[2]=calc_Ion_Mobility(ion_val[2],ion_rad[2],viscosity);		// [meter^2/Volt*second]
	ion_fc[2]=calc_Ion_Friction_Coefficient(ion_val[2],ion_mob[2]);		// [Newton*second/meter]
	// 4th Ion Properties
	ion_val[3]=y[30];
   ion_conc[3]=y[31];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[3]=y[32];													// [meter]
	ion_mob[3]=calc_Ion_Mobility(ion_val[3],ion_rad[3],viscosity);		// [meter^2/Volt*second]
	ion_fc[3]=calc_Ion_Friction_Coefficient(ion_val[3],ion_mob[3]);		// [Newton*second/meter]
	// 5th Ion Properties
	ion_val[4]=y[33];
   ion_conc[4]=y[34];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[4]=y[35];													// [meter]
	ion_mob[4]=calc_Ion_Mobility(ion_val[4],ion_rad[4],viscosity);		// [meter^2/Volt*second]
	ion_fc[4]=calc_Ion_Friction_Coefficient(ion_val[4],ion_mob[4]);		// [Newton*second/meter]
	// 6th Ion Properties
	ion_val[5]=y[36];
   ion_conc[5]=y[37];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[5]=y[38];													// [meter]
	ion_mob[5]=calc_Ion_Mobility(ion_val[5],ion_rad[5],viscosity);		// [meter^2/Volt*second]
	ion_fc[5]=calc_Ion_Friction_Coefficient(ion_val[5],ion_mob[5]);		// [Newton*second/meter]

	double Debye=y[40];	// meters	
	double dielectric=y[39];
	// Convert Concentration to Number Density (particles per volume (m^3))
	ion_conc[0]*=Na;
	ion_conc[1]*=Na;
   ion_conc[2]*=Na;
	ion_conc[3]*=Na;
	ion_conc[4]*=Na;
	ion_conc[5]*=Na;
    // Calculation Constants
	double A,dN_dr;
	int ionIndex;
    // Linearized Governing Equations
	// F FUNCTION
	// df/dr = y[1]
	dy[0]=y[1];
	// d2f/dr2 = y[2]
	dy[1]=y[2];
	// d3f/dr3 = y[3]
	dy[2]=y[3];
	// d4f/dr4 = y[4]
	dy[3]=y[4];
	// d5f/dr5
	A=0;
	for(int i=0;i<numIon;i++)
		{// Ion Equilibrium Number Density
		dN_dr=(-ion_conc[i]*ion_val[i]*ec/(kb*T))*y[6]*exp(-ec*ion_val[i]*y[5]/(kb*T));
		// Ion Contribution
		ionIndex=7+2*i; 
		A+=ion_val[i]*ec*dN_dr*(y[ionIndex]+r);}
   dy[4]=-4*y[4]/r + 4*y[3]/(r*r) + A/(viscosity*r);
	// EQUILIBRIUM POTENTIAL
	// d(potential)/dr = y[6]
	dy[5]=y[6];
	// d2(potential)/dr2
	A=0;
	for(int i=0;i<numIon;i++){A+=ion_val[i]*ion_conc[i]*exp(-ec*ion_val[i]*y[5]/(kb*T));}
	dy[6]=-ec*A/(dielectric*eo) - 2*y[6]/r;
	// INDUCED ELECTROCHEMICAL POTENTIAL FUNCTION (2 ions)
	// 1st Ion
	// d(IEP)/dr = y[8]
   dy[7]=y[8];
	// d2(IEP)/dr2
	A=(ion_val[0]*ec/(kb*T))*y[6];
   dy[8]=-2*y[8]/r + 2*y[7]/(r*r) + A*(y[8] - (2*ion_fc[0]/(r*ion_val[0]*ec))*y[1] + 1);
	// 2nd Ion
	// d(IEP)/dr = y[10]
   dy[9]=y[10];
	// d2(IEP)/dr2
	A=(ion_val[1]*ec/(kb*T))*y[6];
   dy[10]=-2*y[10]/r + 2*y[9]/(r*r) + A*(y[10] - (2*ion_fc[1]/(r*ion_val[1]*ec))*y[1] + 1);
	// 3rd Ion
	// d(IEP)/dr
   dy[11]=y[12];
	// d2(IEP)/dr2
	A=(ion_val[2]*ec/(kb*T))*y[6];
   dy[12]=-2*y[12]/r + 2*y[11]/(r*r) + A*(y[12] - (2*ion_fc[2]/(r*ion_val[2]*ec))*y[1] + 1);
	// 4th Ion
	// d(IEP)/dr
   dy[13]=y[14];
	// d2(IEP)/dr2
	A=(ion_val[3]*ec/(kb*T))*y[6];
   dy[14]=-2*y[14]/r + 2*y[13]/(r*r) + A*(y[14] - (2*ion_fc[3]/(r*ion_val[3]*ec))*y[1] + 1);
	// 5th Ion
	// d(IEP)/dr
   dy[15]=y[16];
	// d2(IEP)/dr2
	A=(ion_val[4]*ec/(kb*T))*y[6];
   dy[16]=-2*y[16]/r + 2*y[15]/(r*r) + A*(y[16] - (2*ion_fc[4]/(r*ion_val[4]*ec))*y[1] + 1);
	// 6th Ion
	// d(IEP)/dr
   dy[17]=y[18];
	// d2(IEP)/dr2
	A=(ion_val[5]*ec/(kb*T))*y[6];
   dy[18]=-2*y[18]/r + 2*y[17]/(r*r) + A*(y[18] - (2*ion_fc[5]/(r*ion_val[5]*ec))*y[1] + 1);

	// Clear Memory
	delete [] ion_conc;delete [] ion_val;delete [] ion_rad;delete [] ion_mob;delete [] ion_fc;}

	// y[0] = f
	// y[1] = df/dr
	// y[2] = d2f/dr2
	// y[3] = d3f/dr3
	// y[4] = d4f/dr4
	// Equilibrium Electric Pontential
	// y[5] = potential
	// y[6] = d(potential)/dr
	// Induced Electrochemical Potential Function (3 ion)
	// y[7] = IEP1
	// y[8] = d(IEP1)/dr
	// y[9] = IEP2
	// y[10] = d(IEP2)/dr
	// y[11] = IEP3
	// y[12] = d(IEP3)/dr
	// y[13] = IEP4
	// y[14] = d(IEP4)/dr
	// y[15] = IEP5
	// y[16] = d(IEP5)/dr
	// Constants
	// y[17] = Temp					// [Kelvin]
	// y[18] = numIon
	// Ion Properties
	// y[19] = Ion1.val;		
	// y[20] = Ion1.conc;		// mM [mol/m^3]
	// y[21] = Ion1.radius;		// meters
	// y[22] = Ion2.val;
	// y[23] = Ion2.conc;
	// y[24] = Ion2.radius;
	// y[25] = Ion3.val;
	// y[26] = Ion3.conc;
	// y[27] = Ion3.radius;
	// y[28] = Ion4.val;
	// y[29] = Ion4.conc;
	// y[30] = Ion4.radius;
	// y[31] = Ion5.val;
	// y[32] = Ion5.conc;
	// y[33] = Ion5.radius;
	// Solution Properties
	// y[34] = relative dielectric	[dimensionless]
	// y[35] = Debye Length [meters]
	// y[36] = viscosity  [Pa s]
void homogeneous_form_5ion(const state_type &y,state_type &dy,const double r)
	{double T=y[17];									// temperature [Kelvin]
	double viscosity=y[36];							// [Newton*second/meter^2]
	// Organize Ion Data
	int numIon=y[18];	
	double* ion_conc=new double[numIon];
	int* ion_val=new int[numIon];
	double* ion_rad=new double[numIon];
	double* ion_mob=new double[numIon];
	double* ion_fc=new double[numIon];
	// 1st Ion Properties
	ion_val[0]=y[19];
   ion_conc[0]=y[20];  		// mM or mol/m^3	    
	ion_rad[0]=y[21];		// meters
	ion_mob[0]=calc_Ion_Mobility(ion_val[0],ion_rad[0],viscosity);
	ion_fc[0]=calc_Ion_Friction_Coefficient(ion_val[0],ion_mob[0]);
	// 2nd Ion Properties
	ion_val[1]=y[22];
   ion_conc[1]=y[23];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[1]=y[24];													// [meter]
	ion_mob[1]=calc_Ion_Mobility(ion_val[1],ion_rad[1],viscosity);		// [meter^2/Volt*second]
	ion_fc[1]=calc_Ion_Friction_Coefficient(ion_val[1],ion_mob[1]);		// [Newton*second/meter]
	// 3rd Ion Properties
	ion_val[2]=y[25];
   ion_conc[2]=y[26];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[2]=y[27];													// [meter]
	ion_mob[2]=calc_Ion_Mobility(ion_val[2],ion_rad[2],viscosity);		// [meter^2/Volt*second]
	ion_fc[2]=calc_Ion_Friction_Coefficient(ion_val[2],ion_mob[2]);		// [Newton*second/meter]
	// 4th Ion Properties
	ion_val[3]=y[28];
   ion_conc[3]=y[29];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[3]=y[30];													// [meter]
	ion_mob[3]=calc_Ion_Mobility(ion_val[3],ion_rad[3],viscosity);		// [meter^2/Volt*second]
	ion_fc[3]=calc_Ion_Friction_Coefficient(ion_val[3],ion_mob[3]);		// [Newton*second/meter]
	// 5th Ion Properties
	ion_val[4]=y[31];
   ion_conc[4]=y[32];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[4]=y[33];													// [meter]
	ion_mob[4]=calc_Ion_Mobility(ion_val[4],ion_rad[4],viscosity);		// [meter^2/Volt*second]
	ion_fc[4]=calc_Ion_Friction_Coefficient(ion_val[4],ion_mob[4]);		// [Newton*second/meter]

	double Debye=y[35];				// meters	
	double dielectric=y[34];
	// Convert Concentration to Number Density (particles per volume (m^3))
	ion_conc[0]*=Na;
	ion_conc[1]*=Na;
	ion_conc[2]*=Na;
	ion_conc[3]*=Na;
	ion_conc[4]*=Na;
	// Calculation Constants
	double A;
	int ionIndex;
    // Linearized Governing Equations
	// F FUNCTION
	// df/dr = y[1]
	dy[0]=y[1];
	// d2f/dr2 = y[2]
	dy[1]=y[2];
	// d3f/dr3 = y[3]
	dy[2]=y[3];
	// d4f/dr4 = y[4]
	dy[3]=y[4];
	// d5f/dr5
	// 1st Ion Equilibrium Number Density
	double dN1_dr=(-ion_conc[0]*ion_val[0]*ec/(kb*T))*y[6]*exp(-ec*ion_val[0]*y[5]/(kb*T));
	double trm1=ion_val[0]*ec*dN1_dr*y[7];
	// 2nd Ion Equilibrium Number Density
	double dN2_dr=(-ion_conc[1]*ion_val[1]*ec/(kb*T))*y[6]*exp(-ec*ion_val[1]*y[5]/(kb*T));
	double trm2=ion_val[1]*ec*dN2_dr*y[9];
	// 3rd Ion Equilibrium Number Density
	double dN3_dr=(-ion_conc[2]*ion_val[2]*ec/(kb*T))*y[6]*exp(-ec*ion_val[2]*y[5]/(kb*T));
	double trm3=ion_val[2]*ec*dN3_dr*y[11];
	// 4th Ion Equilibrium Number Density
	double dN4_dr=(-ion_conc[3]*ion_val[3]*ec/(kb*T))*y[6]*exp(-ec*ion_val[3]*y[5]/(kb*T));
	double trm4=ion_val[3]*ec*dN4_dr*y[13];
	// 5th Ion Equilibrium Number Density
	double dN5_dr=(-ion_conc[4]*ion_val[4]*ec/(kb*T))*y[6]*exp(-ec*ion_val[4]*y[5]/(kb*T));
	double trm5=ion_val[4]*ec*dN5_dr*y[15];
	//
	dy[4]=-4*y[4]/r + 4*y[3]/(r*r) + (trm1+trm2+trm3+trm4+trm5)/(viscosity*r);
	// EQUILIBRIUM POTENTIAL
	// d(potential)/dr = y[6]
	dy[5]=y[6];
	// d2(potential)/dr2
	A=0;
	for(int i=0;i<numIon;i++){A+=ion_val[i]*ion_conc[i]*exp(-ec*ion_val[i]*y[5]/(kb*T));}
	dy[6]=-ec*A/(dielectric*eo) - 2*y[6]/r;
	// INDUCED ELECTROCHEMICAL POTENTIAL FUNCTION (2 ions)
	// 1st Ion
	// d(IEP)/dr
   dy[7]=y[8];
	// d2(IEP)/dr2
	A=(ion_val[0]*ec/(kb*T))*y[6];
   dy[8]=-2*y[8]/r + 2*y[7]/(r*r) + A*(y[8] - (2*ion_fc[0]/(r*ion_val[0]*ec))*y[1]);
	// 2nd Ion
	// d(IEP)/dr
   dy[9]=y[10];
	// d2(IEP)/dr2
	A=(ion_val[1]*ec/(kb*T))*y[6];
   dy[10]=-2*y[10]/r + 2*y[9]/(r*r) + A*(y[10] - (2*ion_fc[1]/(r*ion_val[1]*ec))*y[1]);
	// 3rd Ion
	// d(IEP)/dr
   dy[11]=y[12];
	// d2(IEP)/dr2
	A=(ion_val[2]*ec/(kb*T))*y[6];
   dy[12]=-2*y[12]/r + 2*y[11]/(r*r) + A*(y[12] - (2*ion_fc[2]/(r*ion_val[2]*ec))*y[1]);
	// 4th Ion
	// d(IEP)/dr
   dy[13]=y[14];
	// d2(IEP)/dr2
	A=(ion_val[3]*ec/(kb*T))*y[6];
   dy[14]=-2*y[14]/r + 2*y[13]/(r*r) + A*(y[14] - (2*ion_fc[3]/(r*ion_val[3]*ec))*y[1]);
	// 5th Ion
	// d(IEP)/dr
   dy[15]=y[16];
	// d2(IEP)/dr2
	A=(ion_val[4]*ec/(kb*T))*y[6];
   dy[16]=-2*y[16]/r + 2*y[15]/(r*r) + A*(y[16] - (2*ion_fc[4]/(r*ion_val[4]*ec))*y[1]);

	// Clear Memory
	delete [] ion_conc;delete [] ion_val;delete [] ion_rad;delete [] ion_mob;delete [] ion_fc;}

void inhomogeneous_Prob1_5ion(const state_type &y,state_type &dy,const double r)
	{double T=y[17];									// temperature [Kelvin]
	double viscosity=y[36];							// [Newton*second/meter^2]
	// Organize Ion Data
	int numIon=y[18];	
	double* ion_conc=new double[numIon];
	int* ion_val=new int[numIon];
	double* ion_rad=new double[numIon];
	double* ion_mob=new double[numIon];
	double* ion_fc=new double[numIon];
	// 1st Ion Properties
	ion_val[0]=y[19];
   ion_conc[0]=y[20];  		// mM or mol/m^3	    
	ion_rad[0]=y[21];		// meters
	ion_mob[0]=calc_Ion_Mobility(ion_val[0],ion_rad[0],viscosity);
	ion_fc[0]=calc_Ion_Friction_Coefficient(ion_val[0],ion_mob[0]);
	// 2nd Ion Properties
	ion_val[1]=y[22];
   ion_conc[1]=y[23];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[1]=y[24];													// [meter]
	ion_mob[1]=calc_Ion_Mobility(ion_val[1],ion_rad[1],viscosity);		// [meter^2/Volt*second]
	ion_fc[1]=calc_Ion_Friction_Coefficient(ion_val[1],ion_mob[1]);		// [Newton*second/meter]
	// 3rd Ion Properties
	ion_val[2]=y[25];
   ion_conc[2]=y[26];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[2]=y[27];													// [meter]
	ion_mob[2]=calc_Ion_Mobility(ion_val[2],ion_rad[2],viscosity);		// [meter^2/Volt*second]
	ion_fc[2]=calc_Ion_Friction_Coefficient(ion_val[2],ion_mob[2]);		// [Newton*second/meter]
	// 4th Ion Properties
	ion_val[3]=y[28];
   ion_conc[3]=y[29];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[3]=y[30];													// [meter]
	ion_mob[3]=calc_Ion_Mobility(ion_val[3],ion_rad[3],viscosity);		// [meter^2/Volt*second]
	ion_fc[3]=calc_Ion_Friction_Coefficient(ion_val[3],ion_mob[3]);		// [Newton*second/meter]
	// 5th Ion Properties
	ion_val[4]=y[31];
   ion_conc[4]=y[32];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[4]=y[33];													// [meter]
	ion_mob[4]=calc_Ion_Mobility(ion_val[4],ion_rad[4],viscosity);		// [meter^2/Volt*second]
	ion_fc[4]=calc_Ion_Friction_Coefficient(ion_val[4],ion_mob[4]);		// [Newton*second/meter]

	double Debye=y[35];	// meters	
	double dielectric=y[34];
	// Convert Concentration to Number Density (particles per volume (m^3))
	ion_conc[0]*=Na;
	ion_conc[1]*=Na;
   ion_conc[2]*=Na;
	ion_conc[3]*=Na;
	ion_conc[4]*=Na;
	// Calculation Constants
	double A,dN_dr;
	int ionIndex;

	// Linearized Governing Equations
	// F FUNCTION
	// df/dr = y[1]
	dy[0]=y[1];
	// d2f/dr2 = y[2]
	dy[1]=y[2];
	// d3f/dr3 = y[3]
	dy[2]=y[3];
	// d4f/dr4 = y[4]
	dy[3]=y[4];
	// d5f/dr5
	A=0;
	for(int i=0;i<numIon;i++)
		{// Ion Equilibrium Number Density
		dN_dr=(-ion_conc[i]*ion_val[i]*ec/(kb*T))*y[6]*exp(-ec*ion_val[i]*y[5]/(kb*T));
		// Ion Contribution
		ionIndex=7+2*i; 
		A+=ion_val[i]*ec*dN_dr*y[ionIndex];}
    dy[4]=-4*y[4]/r + 4*y[3]/(r*r) + A/(viscosity*r);
	// EQUILIBRIUM POTENTIAL
	// d(potential)/dr = y[6]
	dy[5]=y[6];
	// d2(potential)/dr2
	A=0;
	for(int i=0;i<numIon;i++){A+=ion_val[i]*ion_conc[i]*exp(-ec*ion_val[i]*y[5]/(kb*T));}
	dy[6]=-ec*A/(dielectric*eo) - 2*y[6]/r;
	// INDUCED ELECTROCHEMICAL POTENTIAL FUNCTION (2 ions)
	// 1st Ion
	// d(IEP)/dr = y[8]
   dy[7]=y[8];
	// d2(IEP)/dr2
	A=(ion_val[0]*ec/(kb*T))*y[6];
   dy[8]=-2*y[8]/r + 2*y[7]/(r*r) + A*(y[8] - (2*ion_fc[0]/(r*ion_val[0]*ec))*y[1] - ion_fc[0]/(ion_val[0]*ec));
	// 2nd Ion
	// d(IEP)/dr = y[10]
   dy[9]=y[10];
	// d2(IEP)/dr2
	A=(ion_val[1]*ec/(kb*T))*y[6];
   dy[10]=-2*y[10]/r + 2*y[9]/(r*r) + A*(y[10] - (2*ion_fc[1]/(r*ion_val[1]*ec))*y[1] - ion_fc[1]/(ion_val[1]*ec));
	// 3rd Ion
	// d(IEP)/dr
   dy[11]=y[12];
	// d2(IEP)/dr2
	A=(ion_val[2]*ec/(kb*T))*y[6];
   dy[12]=-2*y[12]/r + 2*y[11]/(r*r) + A*(y[12] - (2*ion_fc[2]/(r*ion_val[2]*ec))*y[1] - ion_fc[2]/(ion_val[2]*ec));
	// 4th Ion
	// d(IEP)/dr
   dy[13]=y[14];
	// d2(IEP)/dr2
	A=(ion_val[3]*ec/(kb*T))*y[6];
   dy[14]=-2*y[14]/r + 2*y[13]/(r*r) + A*(y[14] - (2*ion_fc[3]/(r*ion_val[3]*ec))*y[1] - ion_fc[3]/(ion_val[3]*ec));
	// 5th Ion
	// d(IEP)/dr
   dy[15]=y[16];
	// d2(IEP)/dr2
	A=(ion_val[4]*ec/(kb*T))*y[6];
   dy[16]=-2*y[16]/r + 2*y[15]/(r*r) + A*(y[16] - (2*ion_fc[4]/(r*ion_val[4]*ec))*y[1] - ion_fc[4]/(ion_val[4]*ec));

	// Clear Memory
	delete [] ion_conc;delete [] ion_val;delete [] ion_rad;delete [] ion_mob;delete [] ion_fc;}

void inhomogeneous_Prob2_5ion(const state_type &y,state_type &dy,const double r)
	{double T=y[17];								// temperature [Kelvin]
	double viscosity=y[36];						// [Newton*second/meter^2]
	// Organize Ion Data
	int numIon=y[18];	
	double* ion_conc=new double[numIon];
	int* ion_val=new int[numIon];
	double* ion_rad=new double[numIon];
	double* ion_mob=new double[numIon];
	double* ion_fc=new double[numIon];
	// 1st Ion Properties
	ion_val[0]=y[19];
   ion_conc[0]=y[20];  		// mM or mol/m^3	    
	ion_rad[0]=y[21];		// meters
	ion_mob[0]=calc_Ion_Mobility(ion_val[0],ion_rad[0],viscosity);
	ion_fc[0]=calc_Ion_Friction_Coefficient(ion_val[0],ion_mob[0]);
	// 2nd Ion Properties
	ion_val[1]=y[22];
   ion_conc[1]=y[23];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[1]=y[24];													// [meter]
	ion_mob[1]=calc_Ion_Mobility(ion_val[1],ion_rad[1],viscosity);		// [meter^2/Volt*second]
	ion_fc[1]=calc_Ion_Friction_Coefficient(ion_val[1],ion_mob[1]);		// [Newton*second/meter]
	// 3rd Ion Properties
	ion_val[2]=y[25];
   ion_conc[2]=y[26];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[2]=y[27];													// [meter]
	ion_mob[2]=calc_Ion_Mobility(ion_val[2],ion_rad[2],viscosity);		// [meter^2/Volt*second]
	ion_fc[2]=calc_Ion_Friction_Coefficient(ion_val[2],ion_mob[2]);		// [Newton*second/meter]
	// 4th Ion Properties
	ion_val[3]=y[28];
   ion_conc[3]=y[29];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[3]=y[30];													// [meter]
	ion_mob[3]=calc_Ion_Mobility(ion_val[3],ion_rad[3],viscosity);		// [meter^2/Volt*second]
	ion_fc[3]=calc_Ion_Friction_Coefficient(ion_val[3],ion_mob[3]);		// [Newton*second/meter]
	// 5th Ion Properties
	ion_val[4]=y[31];
   ion_conc[4]=y[32];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[4]=y[33];													// [meter]
	ion_mob[4]=calc_Ion_Mobility(ion_val[4],ion_rad[4],viscosity);		// [meter^2/Volt*second]
	ion_fc[4]=calc_Ion_Friction_Coefficient(ion_val[4],ion_mob[4]);		// [Newton*second/meter]

	double Debye=y[35];	// meters	
	double dielectric=y[34];
	// Convert Concentration to Number Density (particles per volume (m^3))
	ion_conc[0]*=Na;
	ion_conc[1]*=Na;
   ion_conc[2]*=Na;
	ion_conc[3]*=Na;
	ion_conc[4]*=Na;
    // Calculation Constants
	double A,dN_dr;
	int ionIndex;
    // Linearized Governing Equations
	// F FUNCTION
	// df/dr = y[1]
	dy[0]=y[1];
	// d2f/dr2 = y[2]
	dy[1]=y[2];
	// d3f/dr3 = y[3]
	dy[2]=y[3];
	// d4f/dr4 = y[4]
	dy[3]=y[4];
	// d5f/dr5
	A=0;
	for(int i=0;i<numIon;i++)
		{// Ion Equilibrium Number Density
		dN_dr=(-ion_conc[i]*ion_val[i]*ec/(kb*T))*y[6]*exp(-ec*ion_val[i]*y[5]/(kb*T));
		// Ion Contribution
		ionIndex=7+2*i; 
		A+=ion_val[i]*ec*dN_dr*(y[ionIndex]+r);}
   dy[4]=-4*y[4]/r + 4*y[3]/(r*r) + A/(viscosity*r);
	// EQUILIBRIUM POTENTIAL
	// d(potential)/dr = y[6]
	dy[5]=y[6];
	// d2(potential)/dr2
	A=0;
	for(int i=0;i<numIon;i++){A+=ion_val[i]*ion_conc[i]*exp(-ec*ion_val[i]*y[5]/(kb*T));}
	dy[6]=-ec*A/(dielectric*eo) - 2*y[6]/r;
	// INDUCED ELECTROCHEMICAL POTENTIAL FUNCTION (2 ions)
	// 1st Ion
	// d(IEP)/dr = y[8]
   dy[7]=y[8];
	// d2(IEP)/dr2
	A=(ion_val[0]*ec/(kb*T))*y[6];
   dy[8]=-2*y[8]/r + 2*y[7]/(r*r) + A*(y[8] - (2*ion_fc[0]/(r*ion_val[0]*ec))*y[1] + 1);
	// 2nd Ion
	// d(IEP)/dr = y[10]
   dy[9]=y[10];
	// d2(IEP)/dr2
	A=(ion_val[1]*ec/(kb*T))*y[6];
   dy[10]=-2*y[10]/r + 2*y[9]/(r*r) + A*(y[10] - (2*ion_fc[1]/(r*ion_val[1]*ec))*y[1] + 1);
	// 3rd Ion
	// d(IEP)/dr
   dy[11]=y[12];
	// d2(IEP)/dr2
	A=(ion_val[2]*ec/(kb*T))*y[6];
   dy[12]=-2*y[12]/r + 2*y[11]/(r*r) + A*(y[12] - (2*ion_fc[2]/(r*ion_val[2]*ec))*y[1] + 1);
	// 4th Ion
	// d(IEP)/dr
   dy[13]=y[14];
	// d2(IEP)/dr2
	A=(ion_val[3]*ec/(kb*T))*y[6];
   dy[14]=-2*y[14]/r + 2*y[13]/(r*r) + A*(y[14] - (2*ion_fc[3]/(r*ion_val[3]*ec))*y[1] + 1);
	// 5th Ion
	// d(IEP)/dr
   dy[15]=y[16];
	// d2(IEP)/dr2
	A=(ion_val[4]*ec/(kb*T))*y[6];
   dy[16]=-2*y[16]/r + 2*y[15]/(r*r) + A*(y[16] - (2*ion_fc[4]/(r*ion_val[4]*ec))*y[1] + 1);

	// Clear Memory
	delete [] ion_conc;delete [] ion_val;delete [] ion_rad;delete [] ion_mob;delete [] ion_fc;}

	// y[0] = f
	// y[1] = df/dr
	// y[2] = d2f/dr2
	// y[3] = d3f/dr3
	// y[4] = d4f/dr4
	// Equilibrium Electric Pontential
	// y[5] = potential
	// y[6] = d(potential)/dr
	// Induced Electrochemical Potential Function (3 ion)
	// y[7] = IEP1
	// y[8] = d(IEP1)/dr
	// y[9] = IEP2
	// y[10] = d(IEP2)/dr
	// y[11] = IEP3
	// y[12] = d(IEP3)/dr
	// y[13] = IEP4
	// y[14] = d(IEP4)/dr
	// Constants
	// y[15] = Temp					// [Kelvin]
	// y[16] = numIon
	// Ion Properties
	// y[17] = Ion1.val;		
	// y[18] = Ion1.conc;		// mM [mol/m^3]
	// y[19] = Ion1.radius;		// meters
	// y[20] = Ion2.val;
	// y[21] = Ion2.conc;
	// y[22] = Ion2.radius;
	// y[23] = Ion3.val;
	// y[24] = Ion3.conc;
	// y[25] = Ion3.radius;
	// y[26] = Ion4.val;
	// y[27] = Ion4.conc;
	// y[28] = Ion4.radius;
	// Solution Properties
	// y[29] = relative dielectric	[dimensionless]
	// y[30] = Debye Length [meters]
	// y[31] = viscosity  [Pa s]
void homogeneous_form_4ion(const state_type &y,state_type &dy,const double r)
	{double T=y[15];									// temperature [Kelvin]
	double viscosity=y[31];							// [Newton*second/meter^2]
	// Organize Ion Data
	int numIon=y[16];	
	double* ion_conc=new double[numIon];
	int* ion_val=new int[numIon];
	double* ion_rad=new double[numIon];
	double* ion_mob=new double[numIon];
	double* ion_fc=new double[numIon];
	// 1st Ion Properties
	ion_val[0]=y[17];
   ion_conc[0]=y[18];  		// mM or mol/m^3	    
	ion_rad[0]=y[19];		// meters
	ion_mob[0]=calc_Ion_Mobility(ion_val[0],ion_rad[0],viscosity);
	ion_fc[0]=calc_Ion_Friction_Coefficient(ion_val[0],ion_mob[0]);
	// 2nd Ion Properties
	ion_val[1]=y[20];
   ion_conc[1]=y[21];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[1]=y[22];													// [meter]
	ion_mob[1]=calc_Ion_Mobility(ion_val[1],ion_rad[1],viscosity);		// [meter^2/Volt*second]
	ion_fc[1]=calc_Ion_Friction_Coefficient(ion_val[1],ion_mob[1]);		// [Newton*second/meter]
	// 3rd Ion Properties
	ion_val[2]=y[23];
   ion_conc[2]=y[24];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[2]=y[25];													// [meter]
	ion_mob[2]=calc_Ion_Mobility(ion_val[2],ion_rad[2],viscosity);		// [meter^2/Volt*second]
	ion_fc[2]=calc_Ion_Friction_Coefficient(ion_val[2],ion_mob[2]);		// [Newton*second/meter]
	// 4th Ion Properties
	ion_val[3]=y[26];
   ion_conc[3]=y[27];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[3]=y[28];													// [meter]
	ion_mob[3]=calc_Ion_Mobility(ion_val[3],ion_rad[3],viscosity);		// [meter^2/Volt*second]
	ion_fc[3]=calc_Ion_Friction_Coefficient(ion_val[3],ion_mob[3]);		// [Newton*second/meter]

	double Debye=y[30];				// meters	
	double dielectric=y[29];
	// Convert Concentration to Number Density (particles per volume (m^3))
	ion_conc[0]*=Na;
	ion_conc[1]*=Na;
	ion_conc[2]*=Na;
	ion_conc[3]*=Na;
	// Calculation Constants
	double A;
	int ionIndex;
    // Linearized Governing Equations
	// F FUNCTION
	// df/dr = y[1]
	dy[0]=y[1];
	// d2f/dr2 = y[2]
	dy[1]=y[2];
	// d3f/dr3 = y[3]
	dy[2]=y[3];
	// d4f/dr4 = y[4]
	dy[3]=y[4];
	// d5f/dr5
	// 1st Ion Equilibrium Number Density
	double dN1_dr=(-ion_conc[0]*ion_val[0]*ec/(kb*T))*y[6]*exp(-ec*ion_val[0]*y[5]/(kb*T));
	double trm1=ion_val[0]*ec*dN1_dr*y[7];
	// 2nd Ion Equilibrium Number Density
	double dN2_dr=(-ion_conc[1]*ion_val[1]*ec/(kb*T))*y[6]*exp(-ec*ion_val[1]*y[5]/(kb*T));
	double trm2=ion_val[1]*ec*dN2_dr*y[9];
	// 3rd Ion Equilibrium Number Density
	double dN3_dr=(-ion_conc[2]*ion_val[2]*ec/(kb*T))*y[6]*exp(-ec*ion_val[2]*y[5]/(kb*T));
	double trm3=ion_val[2]*ec*dN3_dr*y[11];
	// 4th Ion Equilibrium Number Density
	double dN4_dr=(-ion_conc[3]*ion_val[3]*ec/(kb*T))*y[6]*exp(-ec*ion_val[3]*y[5]/(kb*T));
	double trm4=ion_val[3]*ec*dN4_dr*y[13];
	//
	dy[4]=-4*y[4]/r + 4*y[3]/(r*r) + (trm1+trm2+trm3+trm4)/(viscosity*r);
	// EQUILIBRIUM POTENTIAL
	// d(potential)/dr = y[6]
	dy[5]=y[6];
	// d2(potential)/dr2
	A=0;
	for(int i=0;i<numIon;i++){A+=ion_val[i]*ion_conc[i]*exp(-ec*ion_val[i]*y[5]/(kb*T));}
	dy[6]=-ec*A/(dielectric*eo) - 2*y[6]/r;
	// INDUCED ELECTROCHEMICAL POTENTIAL FUNCTION (2 ions)
	// 1st Ion
	// d(IEP)/dr
   dy[7]=y[8];
	// d2(IEP)/dr2
	A=(ion_val[0]*ec/(kb*T))*y[6];
   dy[8]=-2*y[8]/r + 2*y[7]/(r*r) + A*(y[8] - (2*ion_fc[0]/(r*ion_val[0]*ec))*y[1]);
	// 2nd Ion
	// d(IEP)/dr
   dy[9]=y[10];
	// d2(IEP)/dr2
	A=(ion_val[1]*ec/(kb*T))*y[6];
   dy[10]=-2*y[10]/r + 2*y[9]/(r*r) + A*(y[10] - (2*ion_fc[1]/(r*ion_val[1]*ec))*y[1]);
	// 3rd Ion
	// d(IEP)/dr
   dy[11]=y[12];
	// d2(IEP)/dr2
	A=(ion_val[2]*ec/(kb*T))*y[6];
   dy[12]=-2*y[12]/r + 2*y[11]/(r*r) + A*(y[12] - (2*ion_fc[2]/(r*ion_val[2]*ec))*y[1]);
	// 4th Ion
	// d(IEP)/dr
   dy[13]=y[14];
	// d2(IEP)/dr2
	A=(ion_val[3]*ec/(kb*T))*y[6];
   dy[14]=-2*y[14]/r + 2*y[13]/(r*r) + A*(y[14] - (2*ion_fc[3]/(r*ion_val[3]*ec))*y[1]);

	// Clear Memory
	delete [] ion_conc;delete [] ion_val;delete [] ion_rad;delete [] ion_mob;delete [] ion_fc;}

void inhomogeneous_Prob1_4ion(const state_type &y,state_type &dy,const double r)
	{double T=y[15];									// temperature [Kelvin]
	double viscosity=y[31];							// [Newton*second/meter^2]
	// Organize Ion Data
	int numIon=y[16];	
	double* ion_conc=new double[numIon];
	int* ion_val=new int[numIon];
	double* ion_rad=new double[numIon];
	double* ion_mob=new double[numIon];
	double* ion_fc=new double[numIon];
	// 1st Ion Properties
	ion_val[0]=y[17];
   ion_conc[0]=y[18];  		// mM or mol/m^3	    
	ion_rad[0]=y[19];		// meters
	ion_mob[0]=calc_Ion_Mobility(ion_val[0],ion_rad[0],viscosity);
	ion_fc[0]=calc_Ion_Friction_Coefficient(ion_val[0],ion_mob[0]);
	// 2nd Ion Properties
	ion_val[1]=y[20];
   ion_conc[1]=y[21];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[1]=y[22];													// [meter]
	ion_mob[1]=calc_Ion_Mobility(ion_val[1],ion_rad[1],viscosity);		// [meter^2/Volt*second]
	ion_fc[1]=calc_Ion_Friction_Coefficient(ion_val[1],ion_mob[1]);		// [Newton*second/meter]
	// 3rd Ion Properties
	ion_val[2]=y[23];
   ion_conc[2]=y[24];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[2]=y[25];													// [meter]
	ion_mob[2]=calc_Ion_Mobility(ion_val[2],ion_rad[2],viscosity);		// [meter^2/Volt*second]
	ion_fc[2]=calc_Ion_Friction_Coefficient(ion_val[2],ion_mob[2]);		// [Newton*second/meter]
	// 4th Ion Properties
	ion_val[3]=y[26];
   ion_conc[3]=y[27];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[3]=y[28];													// [meter]
	ion_mob[3]=calc_Ion_Mobility(ion_val[3],ion_rad[3],viscosity);		// [meter^2/Volt*second]
	ion_fc[3]=calc_Ion_Friction_Coefficient(ion_val[3],ion_mob[3]);		// [Newton*second/meter]

	double Debye=y[30];	// meters	
	double dielectric=y[29];
	// Convert Concentration to Number Density (particles per volume (m^3))
	ion_conc[0]*=Na;
	ion_conc[1]*=Na;
   ion_conc[2]*=Na;
	ion_conc[3]*=Na;
	// Calculation Constants
	double A,dN_dr;
	int ionIndex;

	// Linearized Governing Equations
	// F FUNCTION
	// df/dr = y[1]
	dy[0]=y[1];
	// d2f/dr2 = y[2]
	dy[1]=y[2];
	// d3f/dr3 = y[3]
	dy[2]=y[3];
	// d4f/dr4 = y[4]
	dy[3]=y[4];
	// d5f/dr5
	A=0;
	for(int i=0;i<numIon;i++)
		{// Ion Equilibrium Number Density
		dN_dr=(-ion_conc[i]*ion_val[i]*ec/(kb*T))*y[6]*exp(-ec*ion_val[i]*y[5]/(kb*T));
		// Ion Contribution
		ionIndex=7+2*i; 
		A+=ion_val[i]*ec*dN_dr*y[ionIndex];}
    dy[4]=-4*y[4]/r + 4*y[3]/(r*r) + A/(viscosity*r);
	// EQUILIBRIUM POTENTIAL
	// d(potential)/dr = y[6]
	dy[5]=y[6];
	// d2(potential)/dr2
	A=0;
	for(int i=0;i<numIon;i++){A+=ion_val[i]*ion_conc[i]*exp(-ec*ion_val[i]*y[5]/(kb*T));}
	dy[6]=-ec*A/(dielectric*eo) - 2*y[6]/r;
	// INDUCED ELECTROCHEMICAL POTENTIAL FUNCTION (2 ions)
	// 1st Ion
	// d(IEP)/dr = y[8]
   dy[7]=y[8];
	// d2(IEP)/dr2
	A=(ion_val[0]*ec/(kb*T))*y[6];
   dy[8]=-2*y[8]/r + 2*y[7]/(r*r) + A*(y[8] - (2*ion_fc[0]/(r*ion_val[0]*ec))*y[1] - ion_fc[0]/(ion_val[0]*ec));
	// 2nd Ion
	// d(IEP)/dr = y[10]
   dy[9]=y[10];
	// d2(IEP)/dr2
	A=(ion_val[1]*ec/(kb*T))*y[6];
   dy[10]=-2*y[10]/r + 2*y[9]/(r*r) + A*(y[10] - (2*ion_fc[1]/(r*ion_val[1]*ec))*y[1] - ion_fc[1]/(ion_val[1]*ec));
	// 3rd Ion
	// d(IEP)/dr
   dy[11]=y[12];
	// d2(IEP)/dr2
	A=(ion_val[2]*ec/(kb*T))*y[6];
   dy[12]=-2*y[12]/r + 2*y[11]/(r*r) + A*(y[12] - (2*ion_fc[2]/(r*ion_val[2]*ec))*y[1] - ion_fc[2]/(ion_val[2]*ec));
	// 4th Ion
	// d(IEP)/dr
   dy[13]=y[14];
	// d2(IEP)/dr2
	A=(ion_val[3]*ec/(kb*T))*y[6];
   dy[14]=-2*y[14]/r + 2*y[13]/(r*r) + A*(y[14] - (2*ion_fc[3]/(r*ion_val[3]*ec))*y[1] - ion_fc[3]/(ion_val[3]*ec));

	// Clear Memory
	delete [] ion_conc;delete [] ion_val;delete [] ion_rad;delete [] ion_mob;delete [] ion_fc;}

void inhomogeneous_Prob2_4ion(const state_type &y,state_type &dy,const double r)
	{double T=y[15];								// temperature [Kelvin]
	double viscosity=y[31];						// [Newton*second/meter^2]
	// Organize Ion Data
	int numIon=y[16];	
	double* ion_conc=new double[numIon];
	int* ion_val=new int[numIon];
	double* ion_rad=new double[numIon];
	double* ion_mob=new double[numIon];
	double* ion_fc=new double[numIon];
	// 1st Ion Properties
	ion_val[0]=y[17];
   ion_conc[0]=y[18];  		// mM or mol/m^3	    
	ion_rad[0]=y[19];		// meters
	ion_mob[0]=calc_Ion_Mobility(ion_val[0],ion_rad[0],viscosity);
	ion_fc[0]=calc_Ion_Friction_Coefficient(ion_val[0],ion_mob[0]);
	// 2nd Ion Properties
	ion_val[1]=y[20];
   ion_conc[1]=y[21];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[1]=y[22];													// [meter]
	ion_mob[1]=calc_Ion_Mobility(ion_val[1],ion_rad[1],viscosity);		// [meter^2/Volt*second]
	ion_fc[1]=calc_Ion_Friction_Coefficient(ion_val[1],ion_mob[1]);		// [Newton*second/meter]
	// 3rd Ion Properties
	ion_val[2]=y[23];
   ion_conc[2]=y[24];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[2]=y[25];													// [meter]
	ion_mob[2]=calc_Ion_Mobility(ion_val[2],ion_rad[2],viscosity);		// [meter^2/Volt*second]
	ion_fc[2]=calc_Ion_Friction_Coefficient(ion_val[2],ion_mob[2]);		// [Newton*second/meter]
	// 4th Ion Properties
	ion_val[3]=y[26];
   ion_conc[3]=y[27];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[3]=y[28];													// [meter]
	ion_mob[3]=calc_Ion_Mobility(ion_val[3],ion_rad[3],viscosity);		// [meter^2/Volt*second]
	ion_fc[3]=calc_Ion_Friction_Coefficient(ion_val[3],ion_mob[3]);		// [Newton*second/meter]

	double Debye=y[30];	// meters	
	double dielectric=y[29];
	// Convert Concentration to Number Density (particles per volume (m^3))
	ion_conc[0]*=Na;
	ion_conc[1]*=Na;
   ion_conc[2]*=Na;
	ion_conc[3]*=Na;
    // Calculation Constants
	double A,dN_dr;
	int ionIndex;
    // Linearized Governing Equations
	// F FUNCTION
	// df/dr = y[1]
	dy[0]=y[1];
	// d2f/dr2 = y[2]
	dy[1]=y[2];
	// d3f/dr3 = y[3]
	dy[2]=y[3];
	// d4f/dr4 = y[4]
	dy[3]=y[4];
	// d5f/dr5
	A=0;
	for(int i=0;i<numIon;i++)
		{// Ion Equilibrium Number Density
		dN_dr=(-ion_conc[i]*ion_val[i]*ec/(kb*T))*y[6]*exp(-ec*ion_val[i]*y[5]/(kb*T));
		// Ion Contribution
		ionIndex=7+2*i; 
		A+=ion_val[i]*ec*dN_dr*(y[ionIndex]+r);}
   dy[4]=-4*y[4]/r + 4*y[3]/(r*r) + A/(viscosity*r);
	// EQUILIBRIUM POTENTIAL
	// d(potential)/dr = y[6]
	dy[5]=y[6];
	// d2(potential)/dr2
	A=0;
	for(int i=0;i<numIon;i++){A+=ion_val[i]*ion_conc[i]*exp(-ec*ion_val[i]*y[5]/(kb*T));}
	dy[6]=-ec*A/(dielectric*eo) - 2*y[6]/r;
	// INDUCED ELECTROCHEMICAL POTENTIAL FUNCTION (2 ions)
	// 1st Ion
	// d(IEP)/dr = y[8]
   dy[7]=y[8];
	// d2(IEP)/dr2
	A=(ion_val[0]*ec/(kb*T))*y[6];
   dy[8]=-2*y[8]/r + 2*y[7]/(r*r) + A*(y[8] - (2*ion_fc[0]/(r*ion_val[0]*ec))*y[1] + 1);
	// 2nd Ion
	// d(IEP)/dr = y[10]
   dy[9]=y[10];
	// d2(IEP)/dr2
	A=(ion_val[1]*ec/(kb*T))*y[6];
   dy[10]=-2*y[10]/r + 2*y[9]/(r*r) + A*(y[10] - (2*ion_fc[1]/(r*ion_val[1]*ec))*y[1] + 1);
	// 3rd Ion
	// d(IEP)/dr
   dy[11]=y[12];
	// d2(IEP)/dr2
	A=(ion_val[2]*ec/(kb*T))*y[6];
   dy[12]=-2*y[12]/r + 2*y[11]/(r*r) + A*(y[12] - (2*ion_fc[2]/(r*ion_val[2]*ec))*y[1] + 1);
	// 4th Ion
	// d(IEP)/dr
   dy[13]=y[14];
	// d2(IEP)/dr2
	A=(ion_val[3]*ec/(kb*T))*y[6];
   dy[14]=-2*y[14]/r + 2*y[13]/(r*r) + A*(y[14] - (2*ion_fc[3]/(r*ion_val[3]*ec))*y[1] + 1);

	// Clear Memory
	delete [] ion_conc;delete [] ion_val;delete [] ion_rad;delete [] ion_mob;delete [] ion_fc;}

	// y[0] = f
	// y[1] = df/dr
	// y[2] = d2f/dr2
	// y[3] = d3f/dr3
	// y[4] = d4f/dr4
	// Equilibrium Electric Pontential
	// y[5] = potential
	// y[6] = d(potential)/dr
	// Induced Electrochemical Potential Function (3 ion)
	// y[7] = IEP1
	// y[8] = d(IEP1)/dr
	// y[9] = IEP2
	// y[10] = d(IEP2)/dr
	// y[11] = IEP3
	// y[12] = d(IEP3)/dr
	// Constants
	// y[13] = Temp					// [Kelvin]
	// y[14] = numIon
	// Ion Properties
	// y[15] = Ion1.val;		
	// y[16] = Ion1.conc;		// mM [mol/m^3]
	// y[17] = Ion1.radius;		// meters
	// y[18] = Ion2.val;
	// y[19] = Ion2.conc;
	// y[20] = Ion2.radius;
	// y[21] = Ion3.val;
	// y[22] = Ion3.conc;
	// y[23] = Ion3.radius;
	// Solution Properties
	// y[24] = relative dielectric	[dimensionless]
	// y[25] = Debye Length [meters]
	// y[26] = viscosity  [Pa s]
void homogeneous_form_3ion(const state_type &y,state_type &dy,const double r)
	{double T=y[13];									// temperature [Kelvin]
	double viscosity=y[26];							// [Newton*second/meter^2]
	// Organize Ion Data
	int numIon=y[14];	
	double* ion_conc=new double[numIon];
	int* ion_val=new int[numIon];
	double* ion_rad=new double[numIon];
	double* ion_mob=new double[numIon];
	double* ion_fc=new double[numIon];
	// 1st Ion Properties
	ion_val[0]=y[15];
   ion_conc[0]=y[16];  		// mM or mol/m^3	    
	ion_rad[0]=y[17];		// meters
	ion_mob[0]=calc_Ion_Mobility(ion_val[0],ion_rad[0],viscosity);
	ion_fc[0]=calc_Ion_Friction_Coefficient(ion_val[0],ion_mob[0]);
	// 2nd Ion Properties
	ion_val[1]=y[18];
   ion_conc[1]=y[19];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[1]=y[20];													// [meter]
	ion_mob[1]=calc_Ion_Mobility(ion_val[1],ion_rad[1],viscosity);		// [meter^2/Volt*second]
	ion_fc[1]=calc_Ion_Friction_Coefficient(ion_val[1],ion_mob[1]);		// [Newton*second/meter]
	// 3rd Ion Properties
	ion_val[2]=y[21];
   ion_conc[2]=y[22];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[2]=y[23];													// [meter]
	ion_mob[2]=calc_Ion_Mobility(ion_val[2],ion_rad[2],viscosity);		// [meter^2/Volt*second]
	ion_fc[2]=calc_Ion_Friction_Coefficient(ion_val[2],ion_mob[2]);		// [Newton*second/meter]

	double Debye=y[25];				// meters	
	double dielectric=y[24];
	// Convert Concentration to Number Density (particles per volume (m^3))
	ion_conc[0]*=Na;
	ion_conc[1]*=Na;
	ion_conc[2]*=Na;
	// Calculation Constants
	double A;
	int ionIndex;
    // Linearized Governing Equations
	// F FUNCTION
	// df/dr = y[1]
	dy[0]=y[1];
	// d2f/dr2 = y[2]
	dy[1]=y[2];
	// d3f/dr3 = y[3]
	dy[2]=y[3];
	// d4f/dr4 = y[4]
	dy[3]=y[4];
	// d5f/dr5
	// 1st Ion Equilibrium Number Density
	double dN1_dr=(-ion_conc[0]*ion_val[0]*ec/(kb*T))*y[6]*exp(-ec*ion_val[0]*y[5]/(kb*T));
	double trm1=ion_val[0]*ec*dN1_dr*y[7];
	// 2nd Ion Equilibrium Number Density
	double dN2_dr=(-ion_conc[1]*ion_val[1]*ec/(kb*T))*y[6]*exp(-ec*ion_val[1]*y[5]/(kb*T));
	double trm2=ion_val[1]*ec*dN2_dr*y[9];
	// 3rd Ion Equilibrium Number Density
	double dN3_dr=(-ion_conc[2]*ion_val[2]*ec/(kb*T))*y[6]*exp(-ec*ion_val[2]*y[5]/(kb*T));
	double trm3=ion_val[2]*ec*dN3_dr*y[11];
	//
	dy[4]=-4*y[4]/r + 4*y[3]/(r*r) + (trm1+trm2+trm3)/(viscosity*r);
	// EQUILIBRIUM POTENTIAL
	// d(potential)/dr = y[6]
	dy[5]=y[6];
	// d2(potential)/dr2
	A=0;
	for(int i=0;i<numIon;i++){A+=ion_val[i]*ion_conc[i]*exp(-ec*ion_val[i]*y[5]/(kb*T));}
	dy[6]=-ec*A/(dielectric*eo) - 2*y[6]/r;
	// INDUCED ELECTROCHEMICAL POTENTIAL FUNCTION (2 ions)
	// 1st Ion
	// d(IEP)/dr
   dy[7]=y[8];
	// d2(IEP)/dr2
	A=(ion_val[0]*ec/(kb*T))*y[6];
   dy[8]=-2*y[8]/r + 2*y[7]/(r*r) + A*(y[8] - (2*ion_fc[0]/(r*ion_val[0]*ec))*y[1]);
	// 2nd Ion
	// d(IEP)/dr
   dy[9]=y[10];
	// d2(IEP)/dr2
	A=(ion_val[1]*ec/(kb*T))*y[6];
   dy[10]=-2*y[10]/r + 2*y[9]/(r*r) + A*(y[10] - (2*ion_fc[1]/(r*ion_val[1]*ec))*y[1]);
	// 3rd Ion
	// d(IEP)/dr
   dy[11]=y[12];
	// d2(IEP)/dr2
	A=(ion_val[2]*ec/(kb*T))*y[6];
   dy[12]=-2*y[12]/r + 2*y[11]/(r*r) + A*(y[12] - (2*ion_fc[2]/(r*ion_val[2]*ec))*y[1]);

	// Clear Memory
	delete [] ion_conc;delete [] ion_val;delete [] ion_rad;delete [] ion_mob;delete [] ion_fc;}

void inhomogeneous_Prob1_3ion(const state_type &y,state_type &dy,const double r)
	{double T=y[13];									// temperature [Kelvin]
	double viscosity=y[26];							// [Newton*second/meter^2]
	// Organize Ion Data
	int numIon=y[14];	
	double* ion_conc=new double[numIon];
	int* ion_val=new int[numIon];
	double* ion_rad=new double[numIon];
	double* ion_mob=new double[numIon];
	double* ion_fc=new double[numIon];
	// 1st Ion Properties
	ion_val[0]=y[15];
   ion_conc[0]=y[16];  		// mM or mol/m^3	    
	ion_rad[0]=y[17];		// meters
	ion_mob[0]=calc_Ion_Mobility(ion_val[0],ion_rad[0],viscosity);
	ion_fc[0]=calc_Ion_Friction_Coefficient(ion_val[0],ion_mob[0]);
	// 2nd Ion Properties
	ion_val[1]=y[18];
   ion_conc[1]=y[19];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[1]=y[20];													// [meter]
	ion_mob[1]=calc_Ion_Mobility(ion_val[1],ion_rad[1],viscosity);		// [meter^2/Volt*second]
	ion_fc[1]=calc_Ion_Friction_Coefficient(ion_val[1],ion_mob[1]);		// [Newton*second/meter]
	// 3rd Ion Properties
	ion_val[2]=y[21];
   ion_conc[2]=y[22];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[2]=y[23];													// [meter]
	ion_mob[2]=calc_Ion_Mobility(ion_val[2],ion_rad[2],viscosity);		// [meter^2/Volt*second]
	ion_fc[2]=calc_Ion_Friction_Coefficient(ion_val[2],ion_mob[2]);		// [Newton*second/meter]

	double Debye=y[25];	// meters	
	double dielectric=y[24];
	// Convert Concentration to Number Density (particles per volume (m^3))
	ion_conc[0]*=Na;
	ion_conc[1]*=Na;
   ion_conc[2]*=Na;
	// Calculation Constants
	double A,dN_dr;
	int ionIndex;

	// Linearized Governing Equations
	// F FUNCTION
	// df/dr = y[1]
	dy[0]=y[1];
	// d2f/dr2 = y[2]
	dy[1]=y[2];
	// d3f/dr3 = y[3]
	dy[2]=y[3];
	// d4f/dr4 = y[4]
	dy[3]=y[4];
	// d5f/dr5
	A=0;
	for(int i=0;i<numIon;i++)
		{// Ion Equilibrium Number Density
		dN_dr=(-ion_conc[i]*ion_val[i]*ec/(kb*T))*y[6]*exp(-ec*ion_val[i]*y[5]/(kb*T));
		// Ion Contribution
		ionIndex=7+2*i; 
		A+=ion_val[i]*ec*dN_dr*y[ionIndex];}
    dy[4]=-4*y[4]/r + 4*y[3]/(r*r) + A/(viscosity*r);
	// EQUILIBRIUM POTENTIAL
	// d(potential)/dr = y[6]
	dy[5]=y[6];
	// d2(potential)/dr2
	A=0;
	for(int i=0;i<numIon;i++){A+=ion_val[i]*ion_conc[i]*exp(-ec*ion_val[i]*y[5]/(kb*T));}
	dy[6]=-ec*A/(dielectric*eo) - 2*y[6]/r;
	// INDUCED ELECTROCHEMICAL POTENTIAL FUNCTION (2 ions)
	// 1st Ion
	// d(IEP)/dr = y[8]
   dy[7]=y[8];
	// d2(IEP)/dr2
	A=(ion_val[0]*ec/(kb*T))*y[6];
   dy[8]=-2*y[8]/r + 2*y[7]/(r*r) + A*(y[8] - (2*ion_fc[0]/(r*ion_val[0]*ec))*y[1] - ion_fc[0]/(ion_val[0]*ec));
	// 2nd Ion
	// d(IEP)/dr = y[10]
   dy[9]=y[10];
	// d2(IEP)/dr2
	A=(ion_val[1]*ec/(kb*T))*y[6];
   dy[10]=-2*y[10]/r + 2*y[9]/(r*r) + A*(y[10] - (2*ion_fc[1]/(r*ion_val[1]*ec))*y[1] - ion_fc[1]/(ion_val[1]*ec));
	// 3rd Ion
	// d(IEP)/dr
   dy[11]=y[12];
	// d2(IEP)/dr2
	A=(ion_val[2]*ec/(kb*T))*y[6];
   dy[12]=-2*y[12]/r + 2*y[11]/(r*r) + A*(y[12] - (2*ion_fc[2]/(r*ion_val[2]*ec))*y[1] - ion_fc[2]/(ion_val[2]*ec));

	// Clear Memory
	delete [] ion_conc;delete [] ion_val;delete [] ion_rad;delete [] ion_mob;delete [] ion_fc;}

void inhomogeneous_Prob2_3ion(const state_type &y,state_type &dy,const double r)
	{double T=y[13];								// temperature [Kelvin]
	double viscosity=y[26];						// [Newton*second/meter^2]
	// Organize Ion Data
	int numIon=y[14];	
	double* ion_conc=new double[numIon];
	int* ion_val=new int[numIon];
	double* ion_rad=new double[numIon];
	double* ion_mob=new double[numIon];
	double* ion_fc=new double[numIon];
	// 1st Ion Properties
	ion_val[0]=y[15];
   ion_conc[0]=y[16];  		// mM or mol/m^3	    
	ion_rad[0]=y[17];		// meters
	ion_mob[0]=calc_Ion_Mobility(ion_val[0],ion_rad[0],viscosity);
	ion_fc[0]=calc_Ion_Friction_Coefficient(ion_val[0],ion_mob[0]);
	// 2nd Ion Properties
	ion_val[1]=y[18];
   ion_conc[1]=y[19];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[1]=y[20];													// [meter]
	ion_mob[1]=calc_Ion_Mobility(ion_val[1],ion_rad[1],viscosity);		// [meter^2/Volt*second]
	ion_fc[1]=calc_Ion_Friction_Coefficient(ion_val[1],ion_mob[1]);		// [Newton*second/meter]
	// 3rd Ion Properties
	ion_val[2]=y[21];
   ion_conc[2]=y[22];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[2]=y[23];													// [meter]
	ion_mob[2]=calc_Ion_Mobility(ion_val[2],ion_rad[2],viscosity);		// [meter^2/Volt*second]
	ion_fc[2]=calc_Ion_Friction_Coefficient(ion_val[2],ion_mob[2]);		// [Newton*second/meter]

	double Debye=y[25];	// meters	
	double dielectric=y[24];
	// Convert Concentration to Number Density (particles per volume (m^3))
	ion_conc[0]*=Na;
	ion_conc[1]*=Na;
   ion_conc[2]*=Na;
    // Calculation Constants
	double A,dN_dr;
	int ionIndex;
    // Linearized Governing Equations
	// F FUNCTION
	// df/dr = y[1]
	dy[0]=y[1];
	// d2f/dr2 = y[2]
	dy[1]=y[2];
	// d3f/dr3 = y[3]
	dy[2]=y[3];
	// d4f/dr4 = y[4]
	dy[3]=y[4];
	// d5f/dr5
	A=0;
	for(int i=0;i<numIon;i++)
		{// Ion Equilibrium Number Density
		dN_dr=(-ion_conc[i]*ion_val[i]*ec/(kb*T))*y[6]*exp(-ec*ion_val[i]*y[5]/(kb*T));
		// Ion Contribution
		ionIndex=7+2*i; 
		A+=ion_val[i]*ec*dN_dr*(y[ionIndex]+r);}
   dy[4]=-4*y[4]/r + 4*y[3]/(r*r) + A/(viscosity*r);
	// EQUILIBRIUM POTENTIAL
	// d(potential)/dr = y[6]
	dy[5]=y[6];
	// d2(potential)/dr2
	A=0;
	for(int i=0;i<numIon;i++){A+=ion_val[i]*ion_conc[i]*exp(-ec*ion_val[i]*y[5]/(kb*T));}
	dy[6]=-ec*A/(dielectric*eo) - 2*y[6]/r;
	// INDUCED ELECTROCHEMICAL POTENTIAL FUNCTION (2 ions)
	// 1st Ion
	// d(IEP)/dr = y[8]
   dy[7]=y[8];
	// d2(IEP)/dr2
	A=(ion_val[0]*ec/(kb*T))*y[6];
   dy[8]=-2*y[8]/r + 2*y[7]/(r*r) + A*(y[8] - (2*ion_fc[0]/(r*ion_val[0]*ec))*y[1] + 1);
	// 2nd Ion
	// d(IEP)/dr = y[10]
   dy[9]=y[10];
	// d2(IEP)/dr2
	A=(ion_val[1]*ec/(kb*T))*y[6];
   dy[10]=-2*y[10]/r + 2*y[9]/(r*r) + A*(y[10] - (2*ion_fc[1]/(r*ion_val[1]*ec))*y[1] + 1);
	// 3rd Ion
	// d(IEP)/dr
   dy[11]=y[12];
	// d2(IEP)/dr2
	A=(ion_val[2]*ec/(kb*T))*y[6];
   dy[12]=-2*y[12]/r + 2*y[11]/(r*r) + A*(y[12] - (2*ion_fc[2]/(r*ion_val[2]*ec))*y[1] + 1);

	// Clear Memory
	delete [] ion_conc;delete [] ion_val;delete [] ion_rad;delete [] ion_mob;delete [] ion_fc;}

void homogeneous_form_2ion(const state_type &y,state_type &dy,const double r)
	{double T=y[11];									// temperature [Kelvin]
	double viscosity=y[21];							// [Newton*second/meter^2]
	// Organize Ion Data
	int numIon=y[12];	
	double* ion_conc=new double[numIon];
	int* ion_val=new int[numIon];
	double* ion_rad=new double[numIon];
	double* ion_mob=new double[numIon];
	double* ion_fc=new double[numIon];
	// 1st Ion Properties
	ion_val[0]=y[13];
   ion_conc[0]=y[14];  		// mM or mol/m^3	    
	ion_rad[0]=y[15];		// meters
	ion_mob[0]=calc_Ion_Mobility(ion_val[0],ion_rad[0],viscosity);
	ion_fc[0]=calc_Ion_Friction_Coefficient(ion_val[0],ion_mob[0]);
	//ion_fc[0]=calc_Ion_Friction_Coefficient(ion_val[0],0.00675354);
	// 2nd Ion Properties
	ion_val[1]=y[16];
   ion_conc[1]=y[17];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[1]=y[18];													// [meter]
	ion_mob[1]=calc_Ion_Mobility(ion_val[1],ion_rad[1],viscosity);		// [meter^2/Volt*second]
	ion_fc[1]=calc_Ion_Friction_Coefficient(ion_val[1],ion_mob[1]);		// [Newton*second/meter]
	//ion_fc[1]=calc_Ion_Friction_Coefficient(ion_val[1],0.00648494);

	double Debye=y[20];				// meters	
	double dielectric=y[19];
	// Convert Concentration to Number Density (particles per volume (m^3))
	ion_conc[0]*=Na;
	ion_conc[1]*=Na;	

	//cout<<r<<" "<<y[0]<<" "<<y[5]<<" "<<y[7]<<" "<<y[9]<<endl;
	// Calculation Constants
	double A;
	int ionIndex;
    // Linearized Governing Equations
	// F FUNCTION
    // df/dr = y[1]
    dy[0]=y[1];
    // d2f/dr2 = y[2]
    dy[1]=y[2];
    // d3f/dr3 = y[3]
    dy[2]=y[3];
    // d4f/dr4 = y[4]
    dy[3]=y[4];
    // d5f/dr5
	// 1st Ion Equilibrium Number Density
	double dN1_dr=(-ion_conc[0]*ion_val[0]*ec/(kb*T))*y[6]*exp(-ec*ion_val[0]*y[5]/(kb*T));
	double trm1=ion_val[0]*ec*dN1_dr*y[7];
	// 2nd Ion Equilibrium Number Density
	double dN2_dr=(-ion_conc[1]*ion_val[1]*ec/(kb*T))*y[6]*exp(-ec*ion_val[1]*y[5]/(kb*T));
	double trm2=ion_val[1]*ec*dN2_dr*y[9];
	//A=0;
	//for(int i=0;i<numIon;i++)
	//	{// Ion Equilibrium Number Density
	//	dN_dr=(-ion_conc[i]*ion_val[i]*ec/(kb*T))*y[6]*exp(-ec*ion_val[i]*y[5]/(kb*T));
		// Ion Contribution
	//	ionIndex=7+2*i; 
	//	A+=ion_val[i]*ec*dN_dr*y[ionIndex];}
    //dy[4]=-4.0*y[4]/r + 4.0*y[3]/(r*r) + A/(viscosity*r);
	dy[4]=-4*y[4]/r + 4*y[3]/(r*r) + (trm1+trm2)/(viscosity*r);
	// EQUILIBRIUM POTENTIAL
	// d(potential)/dr = y[6]
	dy[5]=y[6];
	// d2(potential)/dr2
	A=0;
	for(int i=0;i<numIon;i++)
		{A+=ion_val[i]*ion_conc[i]*exp(-ec*ion_val[i]*y[5]/(kb*T));}
	dy[6]=-ec*A/(dielectric*eo) - 2*y[6]/r;
	// INDUCED ELECTROCHEMICAL POTENTIAL FUNCTION (2 ions)
	// 1st Ion
	// d(IEP)/dr = y[8]
   dy[7]=y[8];
	// d2(IEP)/dr2
	A=(ion_val[0]*ec/(kb*T))*y[6];
   dy[8]=-2*y[8]/r + 2*y[7]/(r*r) + A*(y[8] - (2*ion_fc[0]/(r*ion_val[0]*ec))*y[1]);
	// 2nd Ion
	// d(IEP)/dr = y[10]
   dy[9]=y[10];
	// d2(IEP)/dr2
	A=(ion_val[1]*ec/(kb*T))*y[6];
   dy[10]=-2*y[10]/r + 2*y[9]/(r*r) + A*(y[10] - (2*ion_fc[1]/(r*ion_val[1]*ec))*y[1]);

	// Clear Memory
	delete [] ion_conc;delete [] ion_val;delete [] ion_rad;delete [] ion_mob;delete [] ion_fc;}

void inhomogeneous_Prob1_2ion(const state_type &y,state_type &dy,const double r)
	{double T=y[11];									// temperature [Kelvin]
	double viscosity=y[21];							// [Newton*second/meter^2]
	// Organize Ion Data
	int numIon=y[12];	
	double* ion_conc=new double[numIon];
	int* ion_val=new int[numIon];
	double* ion_rad=new double[numIon];
	double* ion_mob=new double[numIon];
	double* ion_fc=new double[numIon];
	// 1st Ion Properties
	ion_val[0]=y[13];
   ion_conc[0]=y[14];  		// mM or mol/m^3	    
	ion_rad[0]=y[15];		// meters
	ion_mob[0]=calc_Ion_Mobility(ion_val[0],ion_rad[0],viscosity);
	ion_fc[0]=calc_Ion_Friction_Coefficient(ion_val[0],ion_mob[0]);

	// 2nd Ion Properties
	ion_val[1]=y[16];
   ion_conc[1]=y[17];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[1]=y[18];													// [meter]
	ion_mob[1]=calc_Ion_Mobility(ion_val[1],ion_rad[1],viscosity);		// [meter^2/Volt*second]
	ion_fc[1]=calc_Ion_Friction_Coefficient(ion_val[1],ion_mob[1]);		// [Newton*second/meter]

	double Debye=y[20];	// meters	
	double dielectric=y[19];
	// Convert Concentration to Number Density (particles per volume (m^3))
	ion_conc[0]*=Na;
	ion_conc[1]*=Na;
    
	// Calculation Constants
	double A,dN_dr;
	int ionIndex;

	// Linearized Governing Equations
	// F FUNCTION
	// df/dr = y[1]
	dy[0]=y[1];
	// d2f/dr2 = y[2]
	dy[1]=y[2];
	// d3f/dr3 = y[3]
	dy[2]=y[3];
	// d4f/dr4 = y[4]
	dy[3]=y[4];
	// d5f/dr5
	A=0;
	for(int i=0;i<numIon;i++)
		{// Ion Equilibrium Number Density
		dN_dr=(-ion_conc[i]*ion_val[i]*ec/(kb*T))*y[6]*exp(-ec*ion_val[i]*y[5]/(kb*T));
		// Ion Contribution
		ionIndex=7+2*i; 
		A+=ion_val[i]*ec*dN_dr*y[ionIndex];}
    dy[4]=-4*y[4]/r + 4*y[3]/(r*r) + A/(viscosity*r);
	// EQUILIBRIUM POTENTIAL
	// d(potential)/dr = y[6]
	dy[5]=y[6];
	// d2(potential)/dr2
	A=0;
	for(int i=0;i<numIon;i++)
		{A+=ion_val[i]*ion_conc[i]*exp(-ec*ion_val[i]*y[5]/(kb*T));}
	dy[6]=-ec*A/(dielectric*eo) - 2*y[6]/r;
	// INDUCED ELECTROCHEMICAL POTENTIAL FUNCTION (2 ions)
	// 1st Ion
	// d(IEP)/dr = y[8]
   dy[7]=y[8];
	// d2(IEP)/dr2
	A=(ion_val[0]*ec/(kb*T))*y[6];
   dy[8]=-2*y[8]/r + 2*y[7]/(r*r) + A*(y[8] - (2*ion_fc[0]/(r*ion_val[0]*ec))*y[1] - ion_fc[0]/(ion_val[0]*ec));
	// 2nd Ion
	// d(IEP)/dr = y[10]
   dy[9]=y[10];
	// d2(IEP)/dr2
	A=(ion_val[1]*ec/(kb*T))*y[6];
   dy[10]=-2*y[10]/r + 2*y[9]/(r*r) + A*(y[10] - (2*ion_fc[1]/(r*ion_val[1]*ec))*y[1] - ion_fc[1]/(ion_val[1]*ec));

	// Clear Memory
	delete [] ion_conc;delete [] ion_val;delete [] ion_rad;delete [] ion_mob;delete [] ion_fc;}

void inhomogeneous_Prob2_2ion(const state_type &y,state_type &dy,const double r)
	{double T=y[11];								// temperature [Kelvin]
	double viscosity=y[21];						// [Newton*second/meter^2]
	// Organize Ion Data
	int numIon=y[12];	
	double* ion_conc=new double[numIon];
	int* ion_val=new int[numIon];
	double* ion_rad=new double[numIon];
	double* ion_mob=new double[numIon];
	double* ion_fc=new double[numIon];
	// 1st Ion Properties
	ion_val[0]=y[13];
   ion_conc[0]=y[14];  		// mM or mol/m^3	    
	ion_rad[0]=y[15];		// meters
	ion_mob[0]=calc_Ion_Mobility(ion_val[0],ion_rad[0],viscosity);
	ion_fc[0]=calc_Ion_Friction_Coefficient(ion_val[0],ion_mob[0]);

	// 2nd Ion Properties
	ion_val[1]=y[16];
   ion_conc[1]=y[17];  												// milli-Molar or [mol/meter^3]	    
	ion_rad[1]=y[18];													// [meter]
	ion_mob[1]=calc_Ion_Mobility(ion_val[1],ion_rad[1],viscosity);		// [meter^2/Volt*second]
	ion_fc[1]=calc_Ion_Friction_Coefficient(ion_val[1],ion_mob[1]);		// [Newton*second/meter]

	double Debye=y[20];	// meters	
	double dielectric=y[19];
	// Convert Concentration to Number Density (particles per volume (m^3))
	ion_conc[0]*=Na;
	ion_conc[1]*=Na;
    
    // Calculation Constants
	double A,dN_dr;
	int ionIndex;
    // Linearized Governing Equations
	// F FUNCTION
	// df/dr = y[1]
	dy[0]=y[1];
	// d2f/dr2 = y[2]
	dy[1]=y[2];
	// d3f/dr3 = y[3]
	dy[2]=y[3];
	// d4f/dr4 = y[4]
	dy[3]=y[4];
	// d5f/dr5
	A=0;
	for(int i=0;i<numIon;i++)
		{// Ion Equilibrium Number Density
		dN_dr=(-ion_conc[i]*ion_val[i]*ec/(kb*T))*y[6]*exp(-ec*ion_val[i]*y[5]/(kb*T));
		// Ion Contribution
		ionIndex=7+2*i; 
		A+=ion_val[i]*ec*dN_dr*(y[ionIndex]+r);}
   dy[4]=-4*y[4]/r + 4*y[3]/(r*r) + A/(viscosity*r);
	// EQUILIBRIUM POTENTIAL
	// d(potential)/dr = y[6]
	dy[5]=y[6];
	// d2(potential)/dr2
	A=0;
	for(int i=0;i<numIon;i++)
		{A+=ion_val[i]*ion_conc[i]*exp(-ec*ion_val[i]*y[5]/(kb*T));}
	dy[6]=-ec*A/(dielectric*eo) - 2*y[6]/r;
	// INDUCED ELECTROCHEMICAL POTENTIAL FUNCTION (2 ions)
	// 1st Ion
	// d(IEP)/dr = y[8]
   dy[7]=y[8];
	// d2(IEP)/dr2
	A=(ion_val[0]*ec/(kb*T))*y[6];
   dy[8]=-2*y[8]/r + 2*y[7]/(r*r) + A*(y[8] - (2*ion_fc[0]/(r*ion_val[0]*ec))*y[1] + 1);
	// 2nd Ion
	// d(IEP)/dr = y[10]
   dy[9]=y[10];
	// d2(IEP)/dr2
	A=(ion_val[1]*ec/(kb*T))*y[6];
   dy[10]=-2*y[10]/r + 2*y[9]/(r*r) + A*(y[10] - (2*ion_fc[1]/(r*ion_val[1]*ec))*y[1] + 1);

	// Clear Memory
	delete [] ion_conc;delete [] ion_val;delete [] ion_rad;delete [] ion_mob;delete [] ion_fc;}

// INPUT:
// Potential
// y[0] = p;
// y[1] = dp/dr;
// Constants
// y[2] = Temperature [Kelvin]
// y[3] = numIons
// Ion Data (Example shown for 2 ion; numIon = 2)
// y[4] = 1st Ion Valence
// y[5] = 1st Ion Concentration [milli-Molar] (or equalivalently [mol/meter^3])
// y[6] = 2nd Ion Valence
// y[7] = 2nd Ion Concentration [milli-Molar] (or equalivalently [mol/meter^3])
void poisson_boltzmann(const state_type &y,state_type &dy,const double r)
	{// Input Handling
	double T=y[2];			// Kelvin
	double dielectric=y[4];
	int numIon=y[3];
	
	// Organize Ion Data
	double* ion_conc=new double[numIon];
	int* ion_val=new int[numIon];

	for(int i=0;i<numIon;i++)
		{ion_val[i]=y[5+2*i];
		ion_conc[i]=y[6+2*i]; // mM or mol/m^3
		// Convert Concentration to Number Density (particles per volume (m^3))
		ion_conc[i]*=Na;}
	
	// d(potential)/dr
	dy[0]=y[1];
	// d2(potential)/dr2
	double A=0;
	for(int i=0;i<numIon;i++){A+=ion_val[i]*ion_conc[i]*exp(-ec*ion_val[i]*y[0]/(kb*T));}
   dy[1]=-ec*A/(dielectric*eo) - 2*y[1]/r;
	//double trm1=ion_val[0]*ion_conc[0]*exp(-ec*ion_val[0]*y[0]/(kb*T));
	//double trm2=ion_val[1]*ion_conc[1]*exp(-ec*ion_val[1]*y[0]/(kb*T));
	//dy[1]=-ec*(trm1+trm2)/(dielectric*eo) - 2*y[1]/r;
	// Clear Memory
	delete [] ion_val;delete [] ion_conc;}

void initialize_PBE(state_type &y,double r0,double C,Solution S)
	{// Convert Debye Length to meters
	double Debye=S.debyeLength*1e-10;
	// Potential [V]
	y[0]=C*exp(-r0/Debye)/r0;
	// d(Potential)/dr
	y[1]=-C*(Debye+r0)*exp(-r0/Debye)/(Debye*r0*r0);
	// Constants
	y[2]=S.temperature;		// Temperature [Kevlin]
	y[3]=S.numIon;
	y[4]=S.dielectric;
	// Ion Data
	for(int i=0;i<S.numIon;i++)
		{// Ion Valence
		y[5+2*i]=S.ionVal[i];
		// Ion Concentration [milli-Molar] (or equalivalently [mol/meter^3])
		y[6+2*i]=S.ionConc[i]*1000;}}
