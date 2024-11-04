# include "solution.h"

Solution::Solution(string solvent,string solventConc,string solventConcType,string solute,string soluteConc,double T,string delimiter)
	{// Handle Input Solvent(s) and Solute(s)
	// Solvent Concentration = mole fraction
	// Solute Concentration = mole/L
	// Basic modeling principle: Calculate property of the pure (or mixed) solvent(s), then correct for presence of solutes

	initializeFunctionGroupInformation();		// for parsing solvent

	// Pure Solvent(s) Parameters
	numSolvents=countDelimiter(solvent,delimiter);
	solvents=new string[numSolvents];
	string* val=fillStringArray(solvent,numSolvents,delimiter);
	for(int i=0;i<numSolvents;i++){solvents[i]=interpretMoleculeName(val[i]);}
	delete [] val;
	// Calculate Molecular Weight of Solvent(s)
	solventsMW=new double[numSolvents];
	double* solventsDensity=new double[numSolvents];
	for(int i=0;i<numSolvents;i++){solventsMW[i]=calcMoleculeMolecularWeight(solvents[i]);}
	// Convert Solvent Concentration Units to mole fraction, if not already
	solventsConc=fillDoubleArray(solventConc,numSolvents,delimiter);
	if(solventConcType.compare("massFraction")==0)
		{solventsConc=massFrac2MoleFrac(solventsConc,solventsMW,numSolvents);}
	else if(solventConcType.compare("moleFraction")==0){/*Do Nothing*/}
	else if(solventConcType.compare("volumeFraction")==0)
		{// Get Pure Solvent Densities
		// Compute Pure Solvent(s) Density [kg/L] in g/L
		for(int i=0;i<numSolvents;i++){solventsDensity[i]=1000*calcPureSolventDensity(solvents[i],T);}
		solventsConc=volFrac2MoleFrac(solventsConc,solventsMW,solventsDensity,numSolvents);}
	else
		{cerr<<"Error in Solution::Solution\nUnrecognized solvent concentration units ("<<solventConcType<<")\n";exit(EXIT_FAILURE);}
	// Get Functional Groups Composing Solvent Molecule
	solventFunctionGroup=new string*[numSolvents];
	int tmpN;
	string tmpStr,nLst="",fgLst="",*tmpArr;	
	for(int i=0;i<numSolvents;i++)
		{tmpN=getSolventFunctionGroups(solvents[i],tmpStr);
		nLst+=cnvrtNumToStrng(tmpN,0)+delimiter;
		fgLst+=tmpStr+"\\";
		}
	numFunctionGroup=fillIntArray(nLst,numSolvents,delimiter);
	tmpArr=fillStringArray(fgLst,numSolvents,"\\");
	for(int i=0;i<numSolvents;i++){solventFunctionGroup[i]=fillStringArray(tmpArr[i],numFunctionGroup[i],GLOBAL_DELIMITER);}
	delete [] tmpArr;

	// Solute Parameters
	numSolutes=countDelimiter(solute,delimiter);
	solutes=new string[numSolutes];
	val=fillStringArray(solute,numSolutes,delimiter);
	for(int i=0;i<numSolutes;i++){solutes[i]=interpretMoleculeName(val[i]);}
	delete [] val;
	solutesConc=fillDoubleArray(soluteConc,numSolutes,delimiter);
	// Calculate Molecular Weight of Solute(s)
	solutesMW=new double[numSolutes];
	for(int i=0;i<numSolutes;i++){solutesMW[i]=calcMoleculeMolecularWeight(solutes[i]);}
	
	// Solution Temperature [K]
	temperature=T;
	// Calculate Solution (Solvent(s)+Solute(s)) Density [kg/L]
	density=calcSolutionDensity(T);
	//cout<<density<<" kg/L"<<endl;
	// Calculate Solution (Solvent(s)+Solute(s)) Viscosity [Pa s]
	viscosity=calcSolutionViscosity(T);
	//cout<<viscosity<<" Pa s"<<endl;
	// Calculate Solution (Solvent(s)+Solute(s)) Relative Dielectric
	dielectric=calcSolutionDielectric(T);
	//cout<<dielectric<<endl;
	// Get Ion Data
	ionData=getIonData(soluteConc,solute,delimiter);
	//cout<<ionData<<endl;
	int N=countDelimiter(ionData,",");
	string* All=fill_string_array(ionData,N,",");
	N/=3;
	// Identify Unique Ion(s) valence/size
	string* valence=new string[N];
	string* ionConcentration=new string[N];
	string* ionRadius=new string[N];
	// Define Number of Ionic Species from Solute(s)
	numIon=0;
	string tmp="",tmp2="",tmp3="",rLst="",vLst="",cLst="";
	int indexLst;
	for(int i=0;i<N;i++)
		{valence[i]=All[3*i];
		ionConcentration[i]=All[1+3*i];
		ionRadius[i]=All[2+3*i];}
	bool FOUND;
	// Count Unique valences and radii	
	for(int i=0;i<N;i++)
		{for(int j=0;j<N;j++)
			{if(i!=j)
				{// Determine Unique Radii
				if(!IN_LIST(ionRadius[i],rLst,numIon,delimiter))
					{numIon++;
					rLst+=ionRadius[i]+delimiter;
					vLst+=valence[i]+delimiter;
					cLst+=ionConcentration[i]+delimiter;
					break;}
				else
					{if(!IN_LIST(valence[i],vLst,numIon,delimiter))
						{// Unique Ion Found with same size but different valence
						numIon++;
						rLst+=ionRadius[i]+delimiter;
						vLst+=valence[i]+delimiter;
						cLst+=ionConcentration[i]+delimiter;
						break;}
					}
				}
			}
		}
	
	// Organize Ion Data
	ionVal=fill_int_array(vLst,numIon,delimiter);
	ionConc=fill_double_array(cLst,numIon,delimiter);
	ionRad=fill_double_array(rLst,numIon,delimiter);

	for(int i=0;i<numIon;i++)
		{ionConc[i]=0;
		for(int j=0;j<N;j++)
			{tmp=valence[j];tmp2=ionRadius[j];tmp3=ionConcentration[j];
			if(ionVal[i]==atoi(tmp.c_str()) && ionRad[i]==strtod(tmp2.c_str(),NULL))
				{// Same Size and Valence
				ionConc[i]+=strtod(tmp3.c_str(),NULL);
				}
			}
		}

	double Max=0,value;
	for(int i=0;i<numIon;i++)
		{// Record Ion Radii [Angstrom]
		tmp=All[2+3*i];
		value=strtod(tmp.c_str(),NULL);
		if(value>Max){Max=value;}}
	largestIon=cnvrtNumToStrng(Max,2);
	// Calculate Solution Debye Length (A)
	debyeLength=calcDebyeLength(T);
	//cout<<ionData<<"|"<<debyeLength<<endl;
	delete [] ionRadius;delete [] ionConcentration;delete [] valence;}

double Solution::calcSolutionDensity(double T)
	{double Output;
	// Convert Mole Fractions of Solvent(s) into Mass Fractions
	double* w=moleFrac2MassFrac(solventsConc,solventsMW,numSolvents);
	// Initialize array for pure solvent density values [kg/L]
	double* p=new double[numSolvents];
	// Compute Pure Solvent(s) Density [kg/L]
	for(int i=0;i<numSolvents;i++){p[i]=calcPureSolventDensity(solvents[i],T);}
	double *r=new double[numSolvents];
	double *q=new double[numSolvents];
	double *Q=new double[numSolvents];
	double *o=new double[numSolvents];
	double *RK=new double[numSolvents];
	double *Ar=new double[numSolvents];
	double *Ac=new double[numSolvents];
	double *A=new double[numSolvents];
	double *term=new double[numSolvents];
	double *X;
	double a,tmpDbl,value,total,Om,On,Pnm,Pkm,Pmk,trm1=0,trm2=0,trm3=0,trm4=0,trm5=0,Rk,Rk1,Rk2,Ar1,Ar2,A1,A2,Ac1,Ac2;
	int tmpInt,numTotalGroups,Counter,lastIndex=0;
	string tmp,nVal,group,group2,grp,*totalGroups;
	excessVolume eV,eV2;
	double G=0,Pressure=101325.0*1e3,mixedDens;
	// Define Pure Solvent(s) Mixture Density (Assuming ideal mixture) w/o solutes
	double Sum=0,solventDensity;
	if(numSolvents>=2)
		{// Binary Mixture or Ternary Mixture
		for(int i=0;i<numSolvents;i++)
			{tmpDbl=0;
			value=0;
			for(int j=0;j<numFunctionGroup[i];j++)
				{tmp=solventFunctionGroup[i][j];
				nVal=tmp.substr(0,1);
				tmpInt=atoi(nVal.c_str());
				group=tmp.substr(1,tmp.length()-1);
				//cout<<"|"<<tmpInt<<"|"<<tmp<<"|\n";
				eV=getFunctionGroupInformation(group);
				tmpDbl+=tmpInt*eV.R;
				value+=tmpInt*eV.Q;
				}
			r[i]=tmpDbl;
			q[i]=value;
			}
		//
		total=0;
		value=0;
		numTotalGroups=0;
		for(int i=0;i<numSolvents;i++)
			{total+=r[i]*solventsConc[i];
			value+=q[i]*solventsConc[i];
			numTotalGroups+=numFunctionGroup[i];}
		//
		for(int i=0;i<numSolvents;i++)
			{Q[i]=r[i]*solventsConc[i]/total;
			o[i]=q[i]*solventsConc[i]/value;}
		//
		X=new double[numTotalGroups];
		totalGroups=new string[numTotalGroups];
		total=0;
		Counter=0;
		for(int i=0;i<numSolvents;i++)
			{for(int j=0;j<numFunctionGroup[i];j++)
				{tmp=solventFunctionGroup[i][j];
				nVal=tmp.substr(0,1);
				tmpInt=atoi(nVal.c_str());
				total+=tmpInt*solventsConc[i];
				group=tmp.substr(1,tmp.length()-1);
				totalGroups[Counter]=group;
				Counter++;
				}
			}
		Counter=0;
		for(int i=0;i<numSolvents;i++)
			{for(int j=0;j<numFunctionGroup[i];j++)
				{tmp=solventFunctionGroup[i][j];
				nVal=tmp.substr(0,1);
				tmpInt=atoi(nVal.c_str());
				X[Counter]=tmpInt*solventsConc[i]/total;
				Counter++;
				}
			}
		//
		Sum=0;
		for(int k=0;k<numTotalGroups;k++)
			{group=totalGroups[k];
			eV=getFunctionGroupInformation(group);
			Sum+=eV.Q*X[k];
			}
		//
		for(int s=0;s<numSolvents;s++)
			{//
			term[s]=0;
			if(s==0)
				{for(int k=0;k<numTotalGroups;k++)
					{trm2=0;
					group=totalGroups[k];
					for(int i=0;i<numFunctionGroup[0];i++)
						{group2=totalGroups[i];
						a=getFunctionGroupInteractionParameter(group,group2);
						eV2=getFunctionGroupInformation(group2);
						Om=eV2.Q*X[i]/Sum;
						Pmk=exp(-a/temperature);
						trm2+=Om*Pmk;
						//
						a=getFunctionGroupInteractionParameter(group2,group);
						Pkm=exp(-a/temperature);
						value=Om*Pkm;
						total=0;
						for(int j=0;j<numTotalGroups;j++)
							{grp=totalGroups[j];
							eV=getFunctionGroupInformation(grp);
							On=eV.Q*X[j]/Sum;
							a=getFunctionGroupInteractionParameter(group2,grp);
							Pnm=exp(-a/temperature);
							total+=On+Pnm;
							}
						trm3+=value/total;
						}
					trm2=log(trm2);
					group=totalGroups[k];
					eV=getFunctionGroupInformation(group);
					Rk1=eV.Q*(1-trm2-trm3);
					//trm4+=Rk1;
					term[s]+=Rk1;
					}
				}
			else
				{lastIndex+=numFunctionGroup[s-1];
				for(int k=0;k<numTotalGroups;k++)
					{trm2=0;
					group=totalGroups[k];
					for(int i=lastIndex;i<lastIndex+numFunctionGroup[s];i++)
						{group2=totalGroups[i];
						a=getFunctionGroupInteractionParameter(group,group2);
						eV2=getFunctionGroupInformation(group2);
						Om=eV2.Q*X[i]/Sum;
						Pmk=exp(-a/temperature);
						trm2+=Om*Pmk;
						//
						a=getFunctionGroupInteractionParameter(group2,group);
						Pkm=exp(-a/temperature);
						value=Om*Pkm;
						total=0;
						for(int j=0;j<numTotalGroups;j++)
							{grp=totalGroups[j];
							eV=getFunctionGroupInformation(grp);
							On=eV.Q*X[j]/Sum;
							a=getFunctionGroupInteractionParameter(group2,grp);
							Pnm=exp(-a/temperature);
							total+=On+Pnm;
							}
						trm3+=value/total;
						}
					trm2=log(trm2);
					group=totalGroups[k];
					eV=getFunctionGroupInformation(group);
					Rk2=eV.Q*(1-trm2-trm3);
					//trm5+=Rk2;
					term[s]+=Rk2;
					}				
				}			
			}
		//
		for(int k=0;k<numTotalGroups;k++)
			{trm2=0;
			group=totalGroups[k];
			for(int i=0;i<numTotalGroups;i++)
				{group2=totalGroups[i];
				a=getFunctionGroupInteractionParameter(group,group2);
				Pmk=exp(-a/temperature);
				eV2=getFunctionGroupInformation(group2);
				Om=eV2.Q*X[i]/Sum;
				trm2+=Om*Pmk;
				//
				a=getFunctionGroupInteractionParameter(group2,group);
				Pkm=exp(-a/temperature);
				value=Om*Pkm;
				total=0;
				for(int j=0;j<numTotalGroups;j++)
					{grp=totalGroups[j];
					eV=getFunctionGroupInformation(grp);
					On=eV.Q*X[j]/Sum;
					a=getFunctionGroupInteractionParameter(group2,grp);
					Pnm=exp(-a/temperature);
					total+=On+Pnm;
					}
				trm3+=value/total;
				}
			trm2=log(trm2);
			group=totalGroups[k];
			eV=getFunctionGroupInformation(group);
			Rk=eV.Q*(1-trm2-trm3);
			trm1+=Rk;
			}
		//
		// Activity Coefficient R
		for(int i=0;i<numSolvents;i++){Ar[i]=trm1-term[i];}
		
		// Activity Coefficient C
		trm4=0;
		for(int i=0;i<numSolvents;i++){trm4+=solventsConc[i]*(5*(r[i]-q[i])-(r[i]-1));}
		for(int i=0;i<numSolvents;i++)
			{trm1=log(Q[i]/solventsConc[i]);
			trm2=5*q[i]*log(o[i]/Q[i]);
			trm3=5*(r[i]-q[i])-(r[i]-1);
			trm5=Q[i]*trm4/solventsConc[i];
			//Ac[i]=trm1+trm2+trm3-trm5;
			Ac[i]=1-(trm1+trm2+trm3-trm5);
			}
		// ACtivity Coefficient
		for(int i=0;i<numSolvents;i++){A[i]=Ac[i]+Ar[i];}

		//cout<<"Ethanol A: "<<A[0]<<"\n"<<Ac[0]<<"\n"<<Ar[0]<<"\n";
		//cout<<"Water A: "<<A[1]<<"\n"<<Ac[1]<<"\n"<<Ar[1]<<"\n";
		//G+=x_1*log(abs(A1));
		//G+=x_2*log(abs(A2));
		G=0;
		for(int i=0;i<numSolvents;i++){G+=solventsConc[i]*log(Ac[i]);}
		trm1=Rg;
		G*=trm1;
		G*=temperature;
/*		cout<<"Gibb's Free Energy Loss ("<<T<<" K): "<<G<<" J/mol\n";
		cout<<"Excess Volume: "<<G/Pressure<<" L/mol\n";
		for(int i=0;i<numSolvents;i++)
			{cout<<solvents[i]<<"\n";
			cout<<"A: "<<A[i]<<"|"<<Ac[i]<<"|"<<Ar[i]<<"\n";
			cout<<1000*p[i]<<" g/L\n";
			cout<<solventsMW[i]<<" g/mol\n";
			cout<<solventsConc[i]<<"\n";
			cout<<"--------------------\n";}		
*/
		mixedDens=0;
		for(int i=0;i<numSolvents;i++){mixedDens+=solventsConc[i]*solventsMW[i];}
		trm1=0;
		for(int i=0;i<numSolvents;i++){trm1+=solventsConc[i]*solventsMW[i]/(1000*p[i]);}
		trm3=G/Pressure+trm1;
		//trm3=trm1+trm2;
		mixedDens/=trm3;
		solventDensity=mixedDens/1000;		// convert to kg/L from g/L
		//cout<<"Mixture Density: "<<mixedDens<<" g/L\n";
		}
	else
		{for(int i=0;i<numSolvents;i++){Sum+=w[i]/p[i];}delete [] w;
		solventDensity=1/Sum;
		}
	
	// Convert from kg/L to kg/m^3
	solventDensity*=1000;
	double d0,d1,d2,d3,d4,vSolute,wTotal=0,wSolvent,MW,S;
	double a1,a2,a3,b1,b2,b3,c1,c2,c3;
	double *Coeff;
	tmp=solutes[0];
	if(numSolutes==1 && tmp.compare("KH2PO4")==0)
		{// 2007. Temperature and Concentration Dependence of Density of Model Liquid Foods
		// Convert Absolute Temperature to celsius
		T-=273.15;
		// KH2PO4 Molecular Weight (g/mol)
		MW=136.086;
		// Solute Concentration
		S=solutesConc[0]*MW/10;	// g KH2PO4/100 mL Water
		// Empirical Coefficients
		a1=1.00263;b1=0.00004211;c1=0.00000273;
		a2=-0.00013;b2=5.395e-8;c2=-6.891e-9;
		a3=-3.147e-6;b3=-1.079e-9;c3=6.952e-11;
		// Solution Density [kg/L]
		Output=a1+b1*S+c1*S*S + T*(a2+b2*S+c2*S*S) + T*T*(a3+b3*S+c3*S*S);
		}
	else
		{// Convert Molar Concentration of Solute(s) into Mass Fractions
		w=molarity2MassFrac(solutesConc,solutesMW,numSolutes,solventDensity);	
		Sum=0;
		// Convert Absolute Temperature to celsius
		T-=273.15;
		// Reference: 2004. Model for calculating the density of aqueous electrolyte solutions		
		for(int i=0;i<numSolutes;i++)
			{Coeff=getDensityCoefficients(solutes[i]);
			d0=Coeff[0];d1=Coeff[1];d2=Coeff[2];d3=Coeff[3];d4=Coeff[4];		
			// Calculate Apparent Specific Volume of Each Salt Component
			if(d0==0&&d1==0&&d2==0&&d3==0&&d4==0){vSolute=0;}
			else{vSolute=(w[i]+d2+d3*T)/((d0*w[i]+d1)*exp(0.000001*(T+d4)*(T+d4)));}
			Sum+=w[i]*vSolute;
			wTotal+=w[i];}
		// Calculate Pure Solvent Mix Mass Fraction
		wSolvent=1-wTotal;	
		// Calculate Solution (Solute(s)+Solvent(s)) Density
		Output=1/((wSolvent/solventDensity)+Sum);
		// Convert from kg/m^3 to kg/L
		Output/=1000;}
	return Output;}

double Solution::calcSolutionViscosity(double T)
	{double Output;
	// Convert Mole Fractions of Solvent(s) into Mass Fractions
	double* w=moleFrac2MassFrac(solventsConc,solventsMW,numSolvents);
	// Initialize array for pure solvent viscosity values
	double* v=new double[numSolvents];
	// Compute Pure Solvent(s) Viscosity [mPas]
	for(int i=0;i<numSolvents;i++){v[i]=calcPureSolventViscosity(solvents[i],T);}
	// Define Pure Solvent(s) Mixture Viscosity (Assuming ideal mixture) w/o solutes
	double solventViscosity=1;
	for(int i=0;i<numSolvents;i++){solventViscosity*=pow(v[i],w[i]);}delete [] w;	
	// Convert Molar Concentration of Solute(s) into Mass Fractions
	w=molarity2MassFrac(solutesConc,solutesMW,numSolutes,density*1000);
	// Determine Mass Fraction of Solvent
	double wSolvent=1;
	for(int i=0;i<numSolutes;i++){wSolvent-=w[i];}
	double* vCoeff;
	double A,B,C,D,E,F,trm,trm1,trm2,trm3,A1;
	double viscSolute,A2=1,A3;
	double u,vo,n,val,val2,w1,w2;
	string tmp=solutes[0],tmp2;
	if(numSolutes==1 && tmp.compare("???")==0)
		{
		}
	else if(numSolutes==2 && (tmp.compare("H3Citrate")==0 || tmp.compare("Na2HPO4")==0) && (tmp2.compare("H3Citrate")==0 && tmp2.compare("Na2HPO4")==0))
		{tmp2=solutes[1];
		// Citrate Phosphate Buffer
		if(tmp.compare("H3Citrate")==0 && tmp2.compare("Na2HPO4")==0)
			{// Calculate Citric Acid Solution Viscosity (mPas)
			w1=w[0];	// Citric Acid Mass Fraction
			w2=w[1];	// Na2HPO4 Mass Fraction
			// Ref. 2014. Modelling of the thermophysical properties of citric acid aqueous solutions. Density and viscosity
			A=-0.85692331;B=-0.00043388783;C=0.14271232;D=0.0004995215;E=-0.18401808;
			trm=A+B*w1*100+C*log(T);
			trm1=1+D*w1*100+E*log(T);
			val=trm/trm1;		// Citric Acid Solution Viscosity (mPas)
			//cerr<<"H3Citrate Viscosity: "<<val<<" mPas\n";
			// Calculate Na2HPO4 Solution Viscosity
			vCoeff=getViscosityCoefficients(tmp2);
			A=vCoeff[0];B=vCoeff[1];C=vCoeff[2];D=vCoeff[3];E=vCoeff[4];F=vCoeff[5];

			trm=pow(1-wSolvent,B);
			trm1=A*trm+C;
			trm2=D*T+1;
			trm=pow(1-wSolvent,F);
			trm3=E*trm+1;

			viscSolute=exp(trm1/(trm2*trm3));
			A2=pow(viscSolute,w2);
			A1=pow(solventViscosity,wSolvent);
			val2=A1*A2;			// Na2HPO4 Solution Viscosity (mPas)
			//cerr<<"Na2HPO4 Viscosity: "<<val2<<" mPas\n";
			A3=pow(val,w1);
			A2=pow(val2,w2);
			Output=A1*A2*A3;
			Output/=1000;	// Convert from mPas to Pas
			}
		else if(tmp.compare("Na2HPO4")==0 && tmp2.compare("H3Citrate")==0)
			{
			w1=w[0];	// Na2HPO4 Mass Fraction
			w2=w[1];	// Citric Acid Mass Fraction
			// Calculate Na2HPO4 Solution Viscosity
			vCoeff=getViscosityCoefficients(tmp);
			A=vCoeff[0];B=vCoeff[1];C=vCoeff[2];D=vCoeff[3];E=vCoeff[4];F=vCoeff[5];

			trm=pow(1-wSolvent,B);
			trm1=A*trm+C;
			trm2=D*T+1;
			trm=pow(1-wSolvent,F);
			trm3=E*trm+1;

			viscSolute=exp(trm1/(trm2*trm3));
			A2=pow(viscSolute,w1);
			A1=pow(solventViscosity,wSolvent);
			val2=A1*A2;			// Na2HPO4 Solution Viscosity (mPas)
			//cerr<<"Na2HPO4 Viscosity: "<<val2<<" mPas\n";
			// Calculate Citric Acid Solution Viscosity
			// Ref. 2014. Modelling of the thermophysical properties of citric acid aqueous solutions. Density and viscosity
			A=-0.85692331;B=-0.00043388783;C=0.14271232;D=0.0004995215;E=-0.18401808;
			trm=A+B*w2*100+C*log(T);
			trm1=1+D*w2*100+E*log(T);
			val=trm/trm1;		// Citric Acid Solution Viscosity (mPas)
			//cerr<<"H3Citrate Viscosity: "<<val<<" mPas\n";
			A3=pow(val,w2);
			A2=pow(val2,w1);
			Output=A1*A2*A3;
			Output/=1000;	// Convert from mPas to Pas
			}
		}
	else
		{// Convert Absolute Temperature to Celsius
		T-=273.15;
		for(int i=0;i<numSolutes;i++)
			{// Empirical Constants
			vCoeff=getViscosityCoefficients(solutes[i]);
			A=vCoeff[0];B=vCoeff[1];C=vCoeff[2];D=vCoeff[3];E=vCoeff[4];F=vCoeff[5];

			trm=pow(1-wSolvent,B);
			trm1=A*trm+C;
			trm2=D*T+1;
			trm=pow(1-wSolvent,F);
			trm3=E*trm+1;

			viscSolute=exp(trm1/(trm2*trm3));
			A2*=pow(viscSolute,w[i]);}
		A1=pow(solventViscosity,wSolvent);
		Output=A1*A2;
		// Convert from mPas to Pas
		Output/=1000;}
	return Output;}

double Solution::calcSolutionDielectric(double T)
	{// Reference: 2001. Computaton of dielectric constants of solvent mixtures and electrolyte solutions
	double Output;
	// Calculate Molar Volume of Solvent(s)
	double* Vm=new double[numSolvents];
	// Initialize array for pure solvent density values [kg/L]
	double* d=new double[numSolvents];
	// Compute Pure Solvent(s) Density [kg/L]
	for(int i=0;i<numSolvents;i++){d[i]=calcPureSolventDensity(solvents[i],T);}
	for(int i=0;i<numSolvents;i++){Vm[i]=solventsMW[i]/(d[i]*1000);}
	// Array for Pure Solvent Relative Dielectric Values
	double* pure_er=new double[numSolvents];
	// Compute Pure Solvent(s) Dielectric
	for(int i=0;i<numSolvents;i++){pure_er[i]=calcPureSolventDielectric(solvents[i],T);}
	// Array for Polarization Values
	double* p=new double[numSolvents];
	// Compute Solvent(s) Polarization
	for(int i=0;i<numSolvents;i++){p[i]=(pure_er[i]-1)*(2*pure_er[i]+1)/(9*pure_er[i]);}
	// Compute Polarization of Solvent Mixture (Oster's Rule)
	double top=0,bottom=0;
	for(int i=0;i<numSolvents;i++){top+=solventsConc[i]*Vm[i]*p[i];bottom+=solventsConc[i]*Vm[i];}
	double solventPolarization=top/bottom;
	// Compute Relative Dielectric of Solvent(s) Mixture
	double solventDielectric=(9*solventPolarization+1+3*sqrt(9*solventPolarization*solventPolarization+2*solventPolarization+1))/4;		
	// Correct for Presence of Solutes (association=increase dielectric,dissociation=decreased dielectric)
	double B,maxB;
	int pos;
	for(int i=0;i<numSolutes;i++)
		{B=getDielectricLinearCoefficient(solutes[i]);
		if(i==0){maxB=B;pos=i;}
		else if(B>maxB){maxB=B;pos=i;}}
	// Currently Only Modeling Dielectric Decrement
	// Applies Greatest Linear Coefficient
	Output=solventDielectric+maxB*solutesConc[pos];
	return Output;}

double Solution::getDielectricLinearCoefficient(string soluteName)
	{// Reference: 2011. Dielectric decrement as a source of ion-specific effects
	double Output;
	if(soluteName.compare("HCl")==0){Output=-19.05;}
	else if(soluteName.compare("LiCl")==0){Output=-12.72;}
	else if(soluteName.compare("NaCl")==0){Output=-11.59;}
	else if(soluteName.compare("KCl")==0){Output=-10.02;}
	else if(soluteName.compare("CsCl")==0){Output=-10.20;}
	else if(soluteName.compare("RbCl")==0){Output=-8.98;}
	else if(soluteName.compare("NaF")==0){Output=-11.90;}
	else if(soluteName.compare("KF")==0){Output=-12.00;}
	else if(soluteName.compare("CsF")==0){Output=-12.60;}
	else if(soluteName.compare("LiI")==0){Output=-14.95;}
	else if(soluteName.compare("NaI")==0){Output=-14.59;}
	else if(soluteName.compare("KI")==0){Output=-14.69;}
	else if(soluteName.compare("CsI")==0){Output=-14.17;}
	else if(soluteName.compare("NaOH")==0){Output=-21;}
	else{//cerr<<"Error in Solution::getDielectricLinearCoefficient!\nUnrecognized solute ("<<soluteName<<").\nUsing pure solvent dielectric as an estimate...\n";
		Output=0;}
	return Output;}

// Reference: 2004. Model for calculating the density of aqueous solutions
double* Solution::getDensityCoefficients(string saltName)
	{double d0,d1,d2,d3,d4;
	double* Output=new double[5];
	if(saltName=="AlCl3"){d0=221.68;d1=160.90;d2=0.15125;d3=0.002500;d4=1500.0;}
	else if(saltName=="Al2(SO4)3"){d0=-0.0017202;d1=0.0018967;d2=-0.030904;d3=0.004087;d4=3804.2;}
	else if(saltName=="BaCl2"){d0=-0.0030518;d1=0.00526670;d2=4.1785;d3=0.068274;d4=3971.9;}
	else if(saltName=="CaCl2"){d0=-0.63254;d1=0.93995;d2=4.2785;d3=0.048319;d4=3180.9;}
	else if(saltName=="CdCl2"){d0=-0.090879;d1=0.29116;d2=7.3827;d3=-0.031855;d4=-3477.5;}
	else if(saltName=="CdSO4"){d0=-1.0440e-7;d1=6.1070e-8;d2=-0.003761;d3=0.004108;d4=5007.7;}
	else if(saltName=="CoCl2"){d0=-8.0924e-8;d1=8.0261e-8;d2=410.24;d3=9.1808;d4=5619.8;}
	else if(saltName=="CoSO4"){d0=-118.36;d1=1368.1;d2=0.01304;d3=-0.000145;d4=-294.02;}
	else if(saltName=="CrCl3"){d0=3.1469;d1=232.160;d2=0.20191;d3=0.002500;d4=1500.0;}
	else if(saltName=="Cr2(SO4)3"){d0=1.0045;d1=1.7697;d2=-0.085017;d3=0.002500;d4=1500.0;}
	else if(saltName=="CuCl2"){d0=1868.5;d1=1137.20;d2=0.07185;d3=0.002565;d4=575.7;}
	else if(saltName=="CuSO4"){d0=-1.9827e-7;d1=1.0883e-7;d2=-0.12506;d3=0.003831;d4=4936.8;}
	else if(saltName=="FeCl2"){d0=98.654;d1=199.51;d2=0.33639;d3=0.0038444;d4=1650.1;}
	else if(saltName=="FeSO4"){d0=-3.6737e-6;d1=1.2465e-4;d2=-0.062861;d3=0.0015696;d4=3943.0;}
	else if(saltName=="FeCl3"){d0=-1333.8;d1=4369.2;d2=1.5298;d3=0.007099;d4=829.21;}
	else if(saltName=="Fe2(SO4)3"){d0=0.47444;d1=-0.64624;d2=-713.10;d3=-25.569;d4=4023.2;}
	else if(saltName=="HCl"){d0=-80.061;d1=255.42;d2=118.42;d3=1.0164;d4=2619.5;}
	else if(saltName=="HCN"){d0=255.82;d1=283.11;d2=0.66888;d3=0.0062057;d4=891.0;}
	else if(saltName=="HNO3"){d0=12.993;d1=-23.579;d2=-3.6070;d3=0.0079416;d4=-2427.1;}
	else if(saltName=="H3PO4"){d0=1358.3;d1=-4327.7;d2=-4.5950;d3=0.0043831;d4=-912.45;}
	else if(saltName=="H2SO4"){d0=89.891;d1=224.48;d2=0.82285;d3=0.0068422;d4=1571.5;}
	else if(saltName=="KCl"){d0=-0.46782;d1=4.30800;d2=2.3780;d3=0.022044;d4=2714.0;}
	else if(saltName=="K2CO3"){d0=-1.4313;d1=2.49170;d2=1.1028;d3=0.013116;d4=2836.0;}
	else if(saltName=="KNO3"){d0=7.5436;d1=26.38800;d2=1.2396;d3=0.011656;d4=2214.0;}
	else if(saltName=="KOH"){d0=194.85;d1=407.31;d2=0.14542;d3=0.002040;d4=1180.9;}
	else if(saltName=="K2SO4"){d0=-2.6681e-5;d1=3.0412e-5;d2=0.97118;d3=0.019816;d4=4366.1;}
	else if(saltName=="LiCl"){d0=17.807;d1=32.011;d2=1.3951;d3=-0.006234;d4=-2131.6;}
	else if(saltName=="Li2SO4"){d0=0.0014730;d1=0.0026934;d2=0.17699;d3=0.0041319;d4=3640.7;}
	else if(saltName=="MgCl2"){d0=-0.00051500;d1=0.0013444;d2=0.58358;d3=0.0085832;d4=3843.6;}
	else if(saltName=="MgSO4"){d0=3.9412e-7;d1=1.4425e-6;d2=-0.05372;d3=0.002062;d4=4563.3;}
	else if(saltName=="MnCl2"){d0=0.000001869;d1=0.00004545;d2=1.5758;d3=-0.010776;d4=-4369.9;}
	else if(saltName=="MnSO4"){d0=0.0032447;d1=0.057246;d2=0.05136;d3=0.002146;d4=3287.8;}
	else if(saltName=="NaBr"){d0=109.770;d1=513.04;d2=1.54540;d3=0.011019;d4=1618.1;}
	else if(saltName=="NaCl"){d0=-0.00433;d1=0.06471;d2=1.01660;d3=0.014624;d4=3315.6;}
	else if(saltName=="NaClO3"){d0=0.014763;d1=0.024913;d2=1.2924;d3=-0.0076175;d4=-3454.3;}
	else if(saltName=="Na2CO3"){d0=0.012755;d1=0.014217;d2=-0.091456;d3=0.0021342;d4=3342.4;}
	else if(saltName=="NaF"){d0=2.8191e-6;d1=2.1777e-7;d2=-0.041483;d3=0.00021765;d4=4586.9;}
	else if(saltName=="NaHCO3"){d0=-9.4794e-8;d1=1.5657e-7;d2=0.9912;d3=0.022644;d4=4900.2;}
	else if(saltName=="NaH2PO4"){d0=208.77;d1=641.05;d2=0.78893;d3=0.0045520;d4=1198.4;}
	else if(saltName=="Na2HPO4"){d0=1096.7;d1=937.57;d2=0.01424;d3=-0.0005595;d4=-860.20;}
	else if(saltName=="NaHSO3"){d0=6.1384e-6;d1=1.3029e-6;d2=0.13635;d3=-0.0014624;d4=-4472.5;}
	else if(saltName=="NaI"){d0=626.15;d1=1858.2;d2=1.7387;d3=0.010500;d4=1203.3;}
	else if(saltName=="Na2MoO4"){d0=-2.0813;d1=4.8446;d2=4.4342;d3=-0.020815;d4=-2942.3;}
	else if(saltName=="NaNO2"){d0=78.365;d1=298.00;d2=0.96246;d3=0.0021999;d4=1500.0;}
	else if(saltName=="NaNO3"){d0=49.209;d1=94.737;d2=0.77927;d3=0.0075451;d4=1819.2;}
	else if(saltName=="NaOH"){d0=385.55;d1=753.47;d2=-0.10938;d3=0.0006953;d4=542.88;}
	else if(saltName=="Na3PO4"){d0=1015.6;d1=1533.7;d2=-0.15180;d3=0.00013660;d4=173.71;}
	else if(saltName=="Na2SO3"){d0=1.5197e-5;d1=4.3766e-7;d2=0.10296;d3=-0.0015271;d4=-4500.9;}
	else if(saltName=="Na2S2O3"){d0=0.84462;d1=-1.5142;d2=-42.949;d3=0.19335;d4=-3425.9;}
	else if(saltName=="Na2SO4"){d0=-1.2095e-7;d1=4.3474e-7;d2=0.15364;d3=0.0072514;d4=4731.5;}
	else if(saltName=="NH3"){d0=0.12693;d1=0.10470;d2=1.0302;d3=-0.0050803;d4=-2973.7;}
	else if(saltName=="NH4Cl"){d0=6.56150;d1=89.772;d2=4.9024;d3=-0.016574;d4=-2089.3;}
	else if(saltName=="NH4NO3"){d0=1379.3;d1=1124.4;d2=0.65598;d3=0.0014106;d4=176.41;}
	else if(saltName=="(NH4)2SO4"){d0=-123.22;d1=452.59;d2=3.2898;d3=0.016292;d4=1692.4;}
	else if(saltName=="NiCl2"){d0=-1.3900e-6;d1=4.1879e-6;d2=0.77734;d3=-0.0066936;d4=-4638.1;}
	else if(saltName=="NiSO4"){d0=-0.03894;d1=0.22109;d2=-0.14443;d3=0.0009867;d4=3073.8;}
	else if(saltName=="SrCl2"){d0=1.3534e-6;d1=-7.4877e-7;d2=-1.9356;d3=0.010704;d4=-4882.1;}
	else if(saltName=="ZnCl2"){d0=1943.6;d1=304.34;d2=-0.013753;d3=0.0011543;d4=573.79;}
	else if(saltName=="ZnSO4"){d0=18.378;d1=35.927;d2=-0.089193;d3=0.0010773;d4=2066.3;}
	else{//cerr<<"Warning from Solution::getDensityCoefficients!\nThe salt ("<<saltName<<") does not have any data associated with it.\nUsing pure solvent density as an estimate...\n";
		d0=0; d1=0; d2=0; d3=0; d4=0;}
	Output[0]=d0; Output[1]=d1; Output[2]=d2; Output[3]=d3; Output[4]=d4;
	return Output;}

double Solution::calcPureSolventDensity(string solventName,double temperature)
	{// temperature is in K
	double Output,T=temperature,trm,trm1,trm2;
	double c1,c2,c3,c4;
	// Molecular Formulae are in a structural format
	// Water Empirical Constants
	double Aw=-2.8054253e-10,Bw=1.0556302e-7,Cw=4.6170461e-5,Dw=0.0079870401,Ew=16.945176,Fw=999.83952,Gw=0.01687985;	
	// Reference: NIST Thermophysical Property Data
	if(solventName.compare("CH3C[::O]OH")==0)
		{// Acetic Acid = CH3C[::O]OH
		Output=1.27334-4.77633e-4*T-9.8132e-7*T*T;}
	else if(solventName.compare("CH3C[CH3]::O")==0)
		{// Acetone = CH3C[CH3]::O
		Output=1.06011-6.68357e-4*T-8.60506e-7*T*T;}
	else if(solventName.compare("CH3C:::N")==0)
		{// Acetonitrile = CH3C:::N
		T-=273.15;c1=-0.0012148333;c2=0.8075233333;Output=c1*T+c2;}
	else if(solventName.compare("CH3(CH2)3OH")==0)
		{// Butanol = "CH3(CH2)3OH"
		Output=1.15309-2.13475e-3*T+5.15573e-6*T*T-6.38112e-9*T*T*T;}
	else if(solventName.compare("CH3(CH2)3C[::O]CH3")==0)
		{// Butyl Methyl Ketone = CH3(CH2)3C[::O]CH3
		// Reference: Yaws' Handbook of Thermodynamic and Physical Properties of Chemical Compounds
		c1=0.2714;c2=0.2607;c3=587.05;c4=0.2963;
		trm=1-T/c3;trm1=pow(trm,c4);trm2=pow(c2,-trm1);
		//Output=c1*(c2^-((1-T/c3)^c4)); // [kg/L]
		Output=c1*trm2;}
	else if(solventName.compare("*1CHCHC[Cl]CHCH*1CH")==0)
		{// Chlorobenzene = *1CHCHC[Cl]CHCH*1CH
		T-=273.15;c1=-0.0010754545;c2=1.1275909091;
		Output=c1*T+c2;}
	else if(solventName.compare("CHCl3")==0)
		{// Chloroform = CHCl3 [kg/L]
		T-=273.15;c1=-0.0018595238;c2=1.5255833333;
		Output=c1*T+c2;}
	else if(solventName.compare("*1CH2(CH2)4*1CH2")==0)
		{// Cyclohexane = *1CH2(CH2)4*1CH2
		T-=273.15;c1=-0.0009791515;c2=0.7984133333;
		Output=c1*T+c2;}
	else if(solventName.compare("CH3(CH2)3O(CH2)3CH3")==0)
		{// Dibutyl Ether = CH3(CH2)3O(CH2)3CH3 [kg/L]
		// Reference: Yaws' Handbook of Thermodynamic and Physical Properties of Chemical Compounds
		c1=0.2552;c2=0.2599;c3=581;c4=0.2857;
		trm=1-T/c3;trm1=pow(trm,c4);trm2=pow(c2,-trm1);
		//Output=c1*(c2^-((1-T/c3)^c4)); // [kg/L]
		Output=c1*trm2;}	
	else if(solventName.compare("C[Cl]H2C[Cl]H2")==0)
		{// 1,2-Dichloroethane = C[Cl]H2C[Cl]H2
		// Reference: Yaws' Handbook of Thermodynamic and Physical Properties of Chemical Compounds
		c1=0.465;c2=0.2874;c3=561;c4=0.3104;
		trm=1-T/c3;trm1=pow(trm,c4);trm2=pow(c2,-trm1);Output=c1*trm2;}
	else if(solventName.compare("ClC[Cl]H2")==0)
		{// Dichloromethane = ClC[Cl]H2
		T-=273.15;c1=-0.0020708182;c2=1.3647318182;
		Output=c1*T+c2;}
	else if(solventName.compare("CH3OCH2CH2OCH3")==0)
		{// 1,2-Dimethoxyethane = CH3OCH2CH2OCH3
		// Reference: Yaws' Handbook of Thermodynamic and Physical Properties of Chemical Compounds
		c1=0.2989;c2=0.2618;c3=536.15;c4=0.2857;
		trm=1-T/c3;trm1=pow(trm,c4);trm2=pow(c2,-trm1);Output=c1*trm2;}
	else if(solventName.compare("CH3C[::O]N[CH3]CH3")==0)
		{// N,N-Dimethylacetamide = CH3C[::O]N[CH3]CH3 [kg/L]
		// Reference: Yaws' Handbook of Thermodynamic and Physical Properties of Chemical Compounds
		c1=0.2714;c2=0.2327;c3=658;c4=0.2702;
		trm=1-T/c3;trm1=pow(trm,c4);trm2=pow(c2,-trm1);Output=c1*trm2;}
	else if(solventName.compare("CH3N[CH3]C[::O]H")==0)
		{// N,N-dimethylformamide = CH3N[CH3]C[::O]H
		// Reference: Yaws' Handbook of Thermodynamic and Physical Properties of Chemical Compounds
		c1=0.2738;c2=0.2301;c3=647;c4=0.2763;
		trm=1-T/c3;trm1=pow(trm,c4);trm2=pow(c2,-trm1);Output=c1*trm2;}
	else if(solventName.compare("*1CH2OCH2CH2O*1CH2")==0)
		{// 1,4-Dioxane = *1CH2OCH2CH2O*1CH2
		c1=0.3702;c2=10.2813;c3=587;c4=0.3047;
		trm=1-T/c3;trm1=pow(trm,c4);trm2=pow(c2,-trm1);Output=c1*trm2;}
	else if(solventName.compare("CH3S[::O]CH3")==0)
		{// Dimethyl sulfoxide = CH3S[::O]CH3 [kg/L]
		c1=0.3442;c2=0.2534;c3=726;c4=0.322;
		trm=1-T/c3;trm1=pow(trm,c4);trm2=pow(c2,-trm1);Output=c1*trm2;}
	else if(solventName.compare("CH3CH2OH")==0)
		{// Ethanol = CH3CH2OH (kg/L)
		Output=1.16239-2.25788e-3*T+5.30621e-6*T*T-6.6307e-9*T*T*T;}
	else if(solventName.compare("CH3C[::O]OCH2CH3")==0)
		{// Ethyl acetate = CH3C[::O]OCH2CH3 (kg/L)
		T-=273.15;c1=-0.0012899091;c2=0.9263681818;Output=c1*T+c2;}
	else if(solventName.compare("CH3CH2OCH2CH3")==0)
		{// Ethyl ether = CH3CH2OCH2CH3 (kg/L)
		T-=273.15;c1=-0.001259;c2=0.7386954545;Output=c1*T+c2;}
	else if(solventName.compare("HOCH2CH2OH")==0)
		{// Ethylene glycol = HOCH2CH2OH (kg/L)
		c1=0.325;c2=0.255;c3=645;c4=0.172;
		trm=1-T/c3;trm1=pow(trm,c4);trm2=pow(c2,-trm1);Output=c1*trm2;}		
	else if(solventName.compare("CH3(CH2)3C[CH2CH3]HCH2OH")==0)
		{// 2-Ethylhexanol = CH3(CH2)3C[CH2CH3]HCH2OH (kg/L)
		c1=0.2685;c2=0.2613;c3=640.25;c4=0.2773;
		trm=1-T/c3;trm1=pow(trm,c4);trm2=pow(c2,-trm1);Output=c1*trm2;}
	else if(solventName.compare("NH2C[::O]H")==0)
		{// Formamide = NH2C[::O]H (kg/L)
		c1=0.2763;c2=0.2035;c3=771;c4=0.2518;
		trm=1-T/c3;trm1=pow(trm,c4);trm2=pow(c2,-trm1);Output=c1*trm2;}
	else if(solventName.compare("HOCH2C[OH]HCH2OH")==0)
		{// Glycerol = HOCH2C[OH]HCH2OH [kg/L]
		Output=1.43632-5.28889e-4*T-2.29556e-7*T*T;}
	else if(solventName.compare("CH3(CH2)4CH3")==0)
		{// Hexane = CH3(CH2)4CH3 (kg/L)
		T-=273.15;c1=-0.0009644545;c2=0.6785772727;Output=c1*T+c2;}
	else if(solventName.compare("CH3C[CH3]HCH2OH")==0)
		{// Isobutanol = CH3C[CH3]HCH2OH kg/L
		c1=0.2698;c2=0.2721;c3=547.73;c4=0.2344;
		trm=1-T/c3;trm1=pow(trm,c4);trm2=pow(c2,-trm1);Output=c1*trm2;}
	else if(solventName.compare("CH3OH")==0)
		{// Methanol = CH3OH (kg/L)
		Output=1.14445-1.79908e-3*T+3.1645e-6*T*T-3.87839e-9*T*T*T;}
	else if(solventName.compare("CH3O(CH2)2OH")==0)
		{// 2-Methoxyethanol = CH3O(CH2)2OH kg/L
		c1=0.3188;c2=0.255;c3=564;c4=0.2857;
		trm=1-T/c3;trm1=pow(trm,c4);trm2=pow(c2,-trm1);Output=c1*trm2;}
	else if(solventName.compare("CH3CH2C[::O]CH3")==0)
		{// Methyl ethyl ketone = CH3CH2C[::O]CH3 [kg/L]
		c1=0.2676;c2=0.2514;c3=535.5;c4=0.2857;
		trm=1-T/c3;trm1=pow(trm,c4);trm2=pow(c2,-trm1);Output=c1*trm2;}
	else if(solventName.compare("CH3C[CH3]HC[::O]CH3")==0)
		{// Methyl isopropyl ketone = CH3C[CH3]HC[::O]CH3 [kg/L]
		c1=0.2753;c2=0.262;c3=553;c4=0.2857;
		trm=1-T/c3;trm1=pow(trm,c4);trm2=pow(c2,-trm1);Output=c1*trm2;}
	else if(solventName.compare("CH3N::[O]O")==0)
		{// Nitromethane = CH3N[::O]O (kg/L)
		T-=273.15;c1=-0.0014116667;c2=1.1674777778;Output=c1*T+c2;}
	else if(solventName.compare("CH3C[OH]HCH2OH")==0)
		{// Propylene glycol = CH3C[OH]HCH2OH (kg/L)
		T-=273.15;c1=-0.0008510824;c2=1.2898675019;Output=c1*T+c2;}
	else if(solventName.compare("CH3(CH2)2OH")==0)
		{// Propanol = CH3(CH2)2OH (kg/L)
		Output=1.01077-3.99649e-8*T-6.64923e-6*T*T+2.16751e-8*T*T*T-2.46167e-11*T*T*T*T;}
	else if(solventName.compare("CH3C[OH]HCH3")==0)
		{// Isopropanol = CH3C[OH]HCH3 (kg/L)
		Output=0.993868-3.40464e-8*T-6.95191e-6*T*T+2.40158e-8*T*T*T-2.92637e-11*T*T*T*T;}
	else if(solventName.compare("*1CH2(CH2)3*1O")==0)
		{// Tetrahydrofuran = *1CH2(CH2)3*1O [kg/L]
		c1=0.3221;c2=0.2808;c3=540.15;c4=0.2912;
		trm=1-T/c3;trm1=pow(trm,c4);trm2=pow(c2,-trm1);Output=c1*trm2;}
	else if(solventName.compare("*1CH(CH)2C[CH2(CH2)2*2CH2]*2C*1CH")==0)
		{// Tetralin = *1CH(CH)2C[CH2(CH2)2*2CH2]*2C*1CH [kg/L]
		c1=0.2984;c2=0.2575;c3=720.15;c4=0.2677;
		trm=1-T/c3;trm1=pow(trm,c4);trm2=pow(c2,-trm1);Output=c1*trm2;}
	else if(solventName.compare("CH3C[CH3]HCH2C[::O]CH3")==0)
		{// Methyl isobutyl ketone = CH3C[CH3]HCH2C[::O]CH3 kg/L
		c1=0.2665;c2=0.2589;c3=571.4;c4=0.2857;
		trm=1-T/c3;trm1=pow(trm,c4);trm2=pow(c2,-trm1);Output=c1*trm2;}
	else if(solventName.compare("*1CH(CH)4*1C[CH3]")==0)
		{// Toluene = *1CH(CH)4*1C[CH3] (kg/L)
		T-=273.15;c1=-0.000953;c2=0.8859045455;Output=c1*T+c2;}
	else if(solventName.compare("ClC[Cl]::C[Cl]H")==0)
		{// Trichloroethene = ClC[Cl]::C[Cl]H (kg/L)
		c1=0.5042;c2=0.2695;c3=571;c4=0.2857;
		trm=1-T/c3;trm1=pow(trm,c4);trm2=pow(c2,-trm1);Output=c1*trm2;}
	else if(solventName.compare("H2O")==0)
		{// Pure Water Density (kg/m^3) (or g/L)
		T-=273.15;	// Convert to centigrade
		// Reference: 2004. Model for calculating the density of aqueous electrolyte solutions
		Output=(((((Aw*T+Bw)*T-Cw)*T-Dw)*T+Ew)*T+Fw)/(1+Gw*T);
		// Convert to kg/L
		Output/=1000;}
	else{cerr<<"Error in Solution::calcPureSolventDensity!!!\nSolvent ("<<solventName<<") is not registered\n";exit(EXIT_FAILURE);}
	return Output;}

double Solution::calcPureSolventViscosity(string solventName,double temperature)
	{// temperature is in K
	double Output,T=temperature,trm;
	double c1,c2,c3,c4;
	// Molecular Formulae are in a structural format
	if(solventName.compare("CH3C[::O]OH")==0)
		{// Acetic Acid = CH3C[::O]OH (mPas) (298.15K-373.15K)
		Output=66.46849*exp(-1.360185e-2*T);}
	else if(solventName.compare("CH3C[CH3]::O")==0)
		{// Acetone = CH3C[CH3]::O (mPas) (248.15K-323.15K)
		Output=4.62452*exp(-9.015176e-3*T);}
	else if(solventName.compare("CH3C:::N")==0)
		{// Acetonitrile = CH3C:::N (mPas) (273.15K-348.15k)
		c1=0.0332045441;c2=679.809658852;Output=c1*exp(c2/T);}
	else if(solventName.compare("CH3(CH2)3OH")==0)
		{// Butanol = CH3(CH2)3OH (mPas) (248.15K-373.15K)
		Output=2146.186*exp(-2.256471e-2*T);}
	else if(solventName.compare("CH3(CH2)3C[::O]CH3")==0)
		{// Butyl Methyl Ketone = CH3(CH2)3C[::O]CH3 (2-Hexanone) (mPas) (248.15-373.15)
		c1=0.0108971435;c2=1186.5587172499;Output=c1*exp(c2/T);}
	else if(solventName.compare("*1CHCHC[Cl]CHCH*1CH")==0)
		{// Chlorobenzene = *1CHCHC[Cl]CHCH*1CH (mPas) (248.15K-373.15K)
		c1=0.0177206558;c2=1132.9079132319;Output=c1*exp(c2/T);}
	else if(solventName.compare("CHCl3")==0)
		{// Chloroform = CHCl3 (mPas) (248.15-323.15)
		c1=0.0266062965;c2=896.9468102971;Output=c1*exp(c2/T);}
	else if(solventName.compare("*1CH2(CH2)4*1CH2")==0)
		{// Cyclohexane = *1CH2(CH2)4*1CH2 (mPas) (298.15K-348.15K)
		c1=0.0071658161;c2=1438.9863354079;Output=c1*exp(c2/T);}
	else if(solventName.compare("CH3(CH2)3O(CH2)3CH3")==0)
		{// Dibutyl Ether = CH3(CH2)3O(CH2)3CH3 (mPas) (248.15-373.15)
		c1=0.0113182343;c2=1198.5352582595;Output=c1*exp(c2/T);}
	else if(solventName.compare("C[Cl]H2C[Cl]H2")==0)
		{// 1,2-Dichloroethane = C[Cl]H2C[Cl]H2 (mPas) (273.15-348.15)
		c1=0.0155035101;c2=1170.3034886562;Output=c1*exp(c2/T);}
	else if(solventName.compare("ClC[Cl]H2")==0)
		{// Dichloromethane = ClC[Cl]H2 (mPas) (248.15K-298.15K)
		c1=0.0249525467;c2=836.7495385873;Output=c1*exp(c2/T);}
	else if(solventName.compare("CH3OCH2CH2OCH3")==0)
		{// 1,2-Dimethoxyethane = CH3OCH2CH2OCH3 (mPas) (215-536)
		c1=-5.7926;c2=8.13e2;c3=0.0136;c4=-1.50e-5;trm=c1+c2/T+c3*T+c4*T*T;Output=pow(10,trm);}
	else if(solventName.compare("CH3C[::O]N[CH3]CH3")==0)
		{// N,N-Dimethylacetamide = CH3C[::O]N[CH3]CH3 (cP=mPas) (253K-658K)
		c1=-4.653;c2=8.36e2;c3=0.0083;c4=-7.83e-6;trm=c1+c2/T+c3*T+c4*T*T;Output=pow(10,trm);}
	else if(solventName.compare("CH3N[CH3]C[::O]H")==0)
		{// N,N-dimethylformamide = CH3N[CH3]C[::O]H (mPas) (273.15-323.15)
		c1=0.0195723424;c2=1118.7559320768;Output=c1*exp(c2/T);}
	else if(solventName.compare("*1CH2OCH2CH2O*1CH2")==0)
		{// 1,4-Dioxane = *1CH2OCH2CH2O*1CH2 (mPas) (298.15-348.15)
		c1=0.0074611653;c2=1508.9408756918;Output=c1*exp(c2/T);}
	else if(solventName.compare("CH3S[::O]CH3")==0)
		{// Dimethyl sulfoxide = CH3S[::O]CH3 (mPas) (298.15-323.15)
		Output=126.2972*exp(-1.406096e-2*T);}
	else if(solventName.compare("CH3CH2OH")==0)
		{// Ethanol = CH3CH2OH (mPas) (248.15-348.15)
		trm=-6.4406+1120/T+0.0137*T-1.55e-5*T*T;
		Output=pow(10,trm);}
	else if(solventName.compare("CH3C[::O]OCH2CH3")==0)
		{// Ethyl acetate = CH3C[::O]OCH2CH3 (mPas) (273.15-348.15)
		c1=0.0139186245;c2=1017.8514197246;Output=c1*exp(c2/T);}
	else if(solventName.compare("CH3CH2OCH2CH3")==0)
		{// Ethyl ether = CH3CH2OCH2CH3 (mPas) (273.15K-298.15k)
		c1=0.0174116237;c2=761.6265724486;Output=c1*exp(c2/T);}
	else if(solventName.compare("HOCH2CH2OH")==0)
		{// Ethylene glycol = HOCH2CH2OH (mPas) (298.15-373.15)
		c1=0.0004756685;c2=3108.8461427998;Output=c1*exp(c2/T);}		
	else if(solventName.compare("CH3(CH2)3C[CH2CH3]HCH2OH")==0)
		{// 2-Ethylhexanol = CH3(CH2)3C[CH2CH3]HCH2OH (mPas) (273.15-373.15)
		c1=0.0001158614;c2=3303.2710820549;Output=c1*exp(c2/T);}
	else if(solventName.compare("NH2C[::O]H")==0)
		{// Formamide = NH2C[::O]H (mPas) (273.15-323.15)
		c1=0.0011144366;c2=2393.042157062;Output=c1*exp(c2/T);}
	else if(solventName.compare("HOCH2C[OH]HCH2OH")==0)
		{// Glycerol = HOCH2C[OH]HCH2OH (mPas) (298.15-373.15)
		trm=-18.215+4230/T+0.0287*T-1.86e-5*T*T;
		Output=pow(10,trm);}
	else if(solventName.compare("CH3(CH2)4CH3")==0)
		{// Hexane = CH3(CH2)4CH3 (mPas) (273.15-323.15)
		c1=0.0137650712;c2=923.725764523;Output=c1*exp(c2/T);}
	else if(solventName.compare("CH3C[CH3]HCH2OH")==0)
		{// Isobutanol = CH3C[CH3]HCH2OH (mPas) (211-548)
		c1=-1.20e1;c2=2.18e3;c3=0.0238;c4=-2.14e-5;trm=c1+c2/T+c3*T+c4*T*T;Output=pow(10,trm);}
	else if(solventName.compare("CH3OH")==0)
		{// Methanol = CH3OH (mPas) (28.15-298.15)
		trm=-9.0562+1250/T+0.0224*T-2.35e-5*T*T;
		Output=pow(10,trm);}
	else if(solventName.compare("CH3O(CH2)2OH")==0)
		{// 2-Methoxyethanol = CH3O(CH2)2OH (cP=mPas) (293.15-535.8)
		c1=-31.096;c2=4.47e3;c3=0.0721;c4=-5.92e-5;trm=c1+c2/T+c3*T+c4*T*T;Output=pow(10,trm);}
	else if(solventName.compare("CH3CH2C[::O]CH3")==0)
		{// Methyl ethyl ketone = CH3CH2C[::O]CH3 (mPas) (248.15-348.15)
		c1=0.0178606573;c2=917.323904847;Output=c1*exp(c2/T);}
	else if(solventName.compare("CH3C[CH3]HC[::O]CH3")==0)
		{// Methyl isopropyl ketone = CH3C[CH3]HC[::O]CH3 (mPas) (293.15-525.35)
		c1=-9.3257;c2=1.28e3;c3=0.022;c4=-2.13e-5;trm=c1+c2/T+c3*T+c4*T*T;Output=pow(10,trm);}
	else if(solventName.compare("CH3N[::O]O")==0)
		{// Nitromethane = CH3N[::O]O (mPas) (248.15-373.15)
		c1=0.0189276739;c2=1051.6399479432;Output=c1*exp(c2/T);}
	else if(solventName.compare("CH3C[OH]HCH2OH")==0)
		{// Propylene glycol = CH3C[OH]HCH2OH (mPas) (273.15-373.15)
		c1=1.2558135966842e-5;c2=4588.5295496621;Output=c1*exp(c2/T);}
	else if(solventName.compare("CH3(CH2)2OH")==0)
		{// Propanol = CH3(CH2)2OH (mPas) (248.15-348.15)
		trm=-3.7702+992/T+0.0041*T-5.46e-6*T*T;
		Output=pow(10,trm);}
	else if(solventName.compare("CH3C[OH]HCH3")==0)
		{// Isopropanol = CH3C[OH]HCH3 (mPas) (273.15-348.15)
		trm=-0.7009+842/T-0.0086*T+8.3e-6*T*T;
		Output=pow(10,trm);}
	else if(solventName.compare("*1CH2(CH2)3*1O")==0)
		{// Tetrahydrofuran = *1CH2(CH2)3*1O (mPas) (248.15-323.15)
		c1=0.02080995;c2=920.2960256817;Output=c1*exp(c2/T);}
	else if(solventName.compare("*1CH(CH)2C[CH2(CH2)2*2CH2]*2C*1CH")==0)
		{// Tetralin = *1CH(CH)2C[CH2(CH2)2*2CH2]*2C*1CH (mPas) (273-720)
		c1=-6.371;c2=1.27e3;c3=0.0105;c4=-8.12e-6;trm=c1+c2/T+c3*T+c4*T*T;Output=pow(10,trm);}
	else if(solventName.compare("CH3C[CH3]HCH2C[::O]CH3")==0)
		{// Methyl isobutyl ketone = CH3C[CH3]HCH2C[::O]CH3 (mPas) (298.15-323.15)
		c1=0.012121122;c2=1134.710075196;Output=c1*exp(c2/T);}
	else if(solventName.compare("*1CH(CH)4*1C[CH3]")==0)
		{// Toluene = *1CH(CH)4*1C[CH3] (mPas) (248.15-373.15)
		c1=0.0148196137;c2=1083.056833036;Output=c1*exp(c2/T);}
	else if(solventName.compare("ClC[Cl]::C[Cl]H")==0)
		{// Trichloroethene = ClC[Cl]::C[Cl]H (mPas) (273.15-348.15)
		c1=0.0384955636;c2=793.4499135534;Output=c1*exp(c2/T);}
	else if(solventName.compare("H2O")==0)
		{// Pure Water Viscosity (mPa s) Reference: 2007. Model for Calculating the Viscosity of Aqueous Solutions		
		T-=273.15;	// Convert Absolute Temperature to Celsius
		Output=(T+246)/(137.37 + 5.2842*T + 0.05594*T*T);}
	else{cerr<<"Error in Solution::calcPureSolventViscosity\nSolvent ("<<solventName<<") is not registered\n";exit(EXIT_FAILURE);}
	return Output;}

double Solution::calcPureSolventDielectric(string solventName,double temperature)
	{// temperature is in K
	double Output,T=temperature;
	double c1,c2,c3,c4;
	// Molecular Formulae are in a structural format
	if(solventName.compare("CH3C[::O]OH")==0)
		{// Acetic Acid = CH3C[::O]OH (293K-363K) Ref; CRC Handbook of Chemistry and Physics
		c1=-0.15731e2;c2=0.12662;c3=-0.17738e-3;Output=c1+c2*T+c3*T*T;}
	else if(solventName.compare("CH3C[CH3]::O")==0)
		{// Acetone = CH3C[CH3]::O (273K-323K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.88157e2;c2=-0.34300;c3=0.38925e-3;Output=c1+c2*T+c3*T*T;}
	else if(solventName.compare("CH3C:::N")==0)
		{// Acetonitrile = CH3C:::N (288K-333K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.29724e3;c2=-0.15508e1;c3=0.22591e-2;Output=c1+c2*T+c3*T*T;}
	else if(solventName.compare("CH3(CH2)3OH")==0)
		{// Butanol = CH3(CH2)3OH (193K-553K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.10578e3;c2=-0.50587;c3=0.84733e-3;c4=-0.48841e-6;Output=c1+c2*T+c3*T*T+c4*T*T*T;}
	else if(solventName.compare("CH3(CH2)3C[::O]CH3")==0)
		{// Butyl Methyl Ketone = CH3(CH2)3C[::O]CH3 (243K-293K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.70378e2;c2=-0.29385;c3=0.35289e-3;Output=c1+c2*T+c3*T*T;}
	else if(solventName.compare("*1CHCHC[Cl]CHCH*1CH")==0)
		{// Chlorobenzene = *1CHCHC[Cl]CHCH*1CH (293K-430K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.19471e2;c2=-0.70786e-1;c3=0.82466e-4;Output=c1+c2*T+c3*T*T;}
	else if(solventName.compare("CHCl3")==0)
		{// Chloroform = CHCl3 (218K-323K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.15115e2;c2=-0.51830e-1;c3=0.56803e-4;Output=c1+c2*T+c3*T*T;}
	else if(solventName.compare("*1CH2(CH2)4*1CH2")==0)
		{// Cyclohexane = *1CH2(CH2)4*1CH2 (283K-333K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.24293e1;c2=-0.12095e-2;c3=-0.58741e-6;Output=c1+c2*T+c3*T*T;}
	else if(solventName.compare("CH3(CH2)3O(CH2)3CH3")==0)
		{// Dibutyl Ether = CH3(CH2)3O(CH2)3CH3 (293K-314K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.65383e1;c2=-0.16172e-1;c3=0.14969e-4;Output=c1+c2*T+c3*T*T;}
	else if(solventName.compare("C[Cl]H2C[Cl]H2")==0)
		{// 1,2-Dichloroethane = C[Cl]H2C[Cl]H2 (293K-343K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.24404e2;c2=-0.47892e-1;Output=c1+c2*T;}
	else if(solventName.compare("ClC[Cl]H2")==0)
		{// Dichloromethane = ClC[Cl]H2 (184K-338K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.40452e2;c2=-0.17748;c3=0.23942e-3;Output=c1+c2*T+c3*T*T;}
	else if(solventName.compare("CH3OCH2CH2OCH3")==0)
		{// 1,2-Dimethoxyethane = CH3OCH2CH2OCH3 (256K-318K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.48832e2;c2=-0.24218;c3=0.34413e-3;Output=c1+c2*T+c3*T*T;}
	else if(solventName.compare("CH3C[::O]N[CH3]CH3")==0)
		{// N,N-Dimethylacetamide = CH3C[::O]N[CH3]CH3 (294K-433K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.15420e3;c2=-0.57506;c3=0.61911e-3;Output=c1+c2*T+c3*T*T;}
	else if(solventName.compare("CH3N[CH3]C[::O]H")==0)
		{// N,N-dimethylformamide = CH3N[CH3]C[::O]H (213K-353K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.15364e3;c2=-0.60367;c3=0.71505e-3;Output=c1+c2*T+c3*T*T;}
	else if(solventName.compare("*1CH2OCH2CH2O*1CH2")==0)
		{// 1,4-Dioxane = *1CH2OCH2CH2O*1CH2 (293K-313K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.27299e1;c2=-0.17440e-2;Output=c1+c2*T;}
	else if(solventName.compare("CH3S[::O]CH3")==0)
		{// Dimethyl sulfoxide = CH3S[::O]CH3 (288K-343K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.38478e2;c2=0.16939;c3=-0.47423e-3;Output=c1+c2*T+c3*T*T;}
	else if(solventName.compare("CH3CH2OH")==0)
		{// Ethanol = CH3CH2OH (163K-523K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.15145e3;c2=-0.87020;c3=0.19570e-2;c4=-0.15512e-5;Output=c1+c2*T+c3*T*T+c4*T*T*T;}
	else if(solventName.compare("CH3C[::O]OCH2CH3")==0)
		{// Ethyl acetate = CH3C[::O]OCH2CH3 (293K-433K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.15646e2;c2=-0.44066e-1;c3=0.39137e-4;Output=c1+c2*T+c3*T*T;}
	else if(solventName.compare("CH3CH2OCH2CH3")==0)
		{// Ethyl ether = CH3CH2OCH2CH3 (283K-301K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.79725e1;c2=-0.12519e-1;Output=c1+c2*T;}
	else if(solventName.compare("HOCH2CH2OH")==0)
		{// Ethylene glycol = HOCH2CH2OH (293K-423K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.14355e3;c2=-0.48573;c3=0.46703e-3;Output=c1+c2*T+c3*T*T;}
	else if(solventName.compare("CH3(CH2)3C[CH2CH3]HCH2OH")==0)
		{// 2-Ethylhexanol = CH3(CH2)3C[CH2CH3]HCH2OH (208K-318K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.86074e2;c2=-0.42636;c3=0.55078e-3;Output=c1+c2*T+c3*T*T;}
	else if(solventName.compare("NH2C[::O]H")==0)
		{// Formamide = NH2C[::O]H (278K-333K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.26076e3;c2=-0.61145;c3=0.34296e-3;Output=c1+c2*T+c3*T*T;}
	else if(solventName.compare("HOCH2C[OH]HCH2OH")==0)
		{// Glycerol = HOCH2C[OH]HCH2OH (288K-343K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.77503e2;c2=-0.37984e-1;c3=-0.23107e-3;Output=c1+c2*T+c3*T*T;}
	else if(solventName.compare("CH3(CH2)4CH3")==0)
		{// Hexane = CH3(CH2)4CH3 (293K-473K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.19768e1;c2=0.70933e-3;c3=-0.34470e-5;Output=c1+c2*T+c3*T*T;}
	else if(solventName.compare("CH3C[CH3]HCH2OH")==0)
		{// Isobutanol = CH3C[CH3]HCH2OH (173K-533K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.10762e3;c2=-0.51398;c3=0.83702e-3;c4=-0.45299e-6;Output=c1+c2*T+c3*T*T+c4*T*T*T;}
	else if(solventName.compare("CH3OH")==0)
		{// Methanol = CH3OH (177K-293K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.19341e3;c2=-0.92211;c3=0.12839e-2;Output=c1+c2*T+c3*T*T;}
	else if(solventName.compare("CH3O(CH2)2OH")==0)
		{// 2-Methoxyethanol = CH3O(CH2)2OH (254K-318K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.11803e3;c2=-0.58000;c3=0.81001e-3;Output=c1+c2*T+c3*T*T;}
	else if(solventName.compare("CH3CH2C[::O]CH3")==0)
		{// Methyl ethyl ketone = CH3CH2C[::O]CH3 (293K-333K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.15457e2;c2=0.90152e-1;c3=-0.27100e-3;Output=c1+c2*T+c3*T*T;}
	else if(solventName.compare("CH3C[CH3]HC[::O]CH3")==0)
		{// Methyl isopropyl ketone = CH3C[CH3]HC[::O]CH3 (293K-328K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.30695e2;c2=-0.10962;c3=0.13810e-3;Output=c1+c2*T+c3*T*T;}
	else if(solventName.compare("CH3N[::O]O")==0)
		{// Nitromethane = CH3N[::O]O (288K-343K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.11227e3;c2=-0.35591;c3=0.34206e-3;Output=c1+c2*T+c3*T*T;}
	else if(solventName.compare("CH3C[OH]HCH2OH")==0)
		{// Propylene glycol = CH3C[OH]HCH2OH (193K-403K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.24546e3;c2=-0.15738e1;c3=0.38068e-2;c4=-0.32544e-5;Output=c1+c2*T+c3*T*T+c4*T*T*T;}
	else if(solventName.compare("CH3(CH2)2OH")==0)
		{// Propanol = CH3(CH2)2OH (193K-493K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.98045e2;c2=-0.36860;c3=0.36422e-3;Output=c1+c2*T+c3*T*T;}
	else if(solventName.compare("CH3C[OH]HCH3")==0)
		{// Isopropanol = CH3C[OH]HCH3 (193K-493K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.10416e3;c2=-0.41011;c3=0.42049e-3;Output=c1+c2*T+c3*T*T;}
	else if(solventName.compare("*1CH2(CH2)3*1O")==0)
		{// Tetrahydrofuran = *1CH2(CH2)3*1O (224K-295K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.30739e2;c2=-0.12946;c3=0.17195e-3;Output=c1+c2*T+c3*T*T;}
	else if(solventName.compare("*1CH(CH)2C[CH2(CH2)2*2CH2]*2C*1CH")==0)
		{// Tetralin = *1CH(CH)2C[CH2(CH2)2*2CH2]*2C*1CH (298K-343K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.29172e1;c2=0.12832e-2;c3=-0.59453e-5;Output=c1+c2*T+c3*T*T;}
	else if(solventName.compare("CH3C[CH3]HCH2C[::O]CH3")==0)
		{// Methyl isobutyl ketone = CH3C[CH3]HCH2C[::O]CH3 (204K-373K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.36341e2;c2=-0.97119e-1;c3=0.61896e-4;Output=c1+c2*T+c3*T*T;}
	else if(solventName.compare("*1CH(CH)4*1C[CH3]")==0)
		{// Toluene = *1CH(CH)4*1C[CH3] (207k-316K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.32584e1;c2=-0.34410e-2;c3=0.15937e-5;Output=c1+c2*T+c3*T*T;}
	else if(solventName.compare("ClC[Cl]::C[Cl]H")==0)
		{// Trichloroethene = ClC[Cl]::C[Cl]H (302K-338K) Ref; CRC Handbook of Chemistry and Physics
		c1=0.58319e1;c2=-0.80828e-2;Output=c1+c2*T;}
	else if(solventName.compare("H2O")==0)
		{// Reference: 1956. Dielectric constant of water from 0 to 100C
		T=temperature-273.15;c1=87.740;c2=-0.4008;c3=0.0009398;c4=-0.00000141;Output=c1+c2*T+c3*T*T+c4*T*T*T;}
	else{cerr<<"Error in Dielectric::get_pure_solvent_dielectric\nSolvent ("<<solventName<<") is not registered\n";exit(EXIT_FAILURE);}
	return Output;}

int Solution::countDelimiter(string Data,string delimiter){int Counter=0;string tmp="";for(int i=0;i<Data.length();i++){tmp=Data[i];if(tmp.compare(delimiter)==0){Counter++;}}return Counter;}

double* Solution::fillDoubleArray(string Data,int numPnts,string delimiter)
	{double* Output=new double[numPnts];
	string bld="",tmp="";
	int Counter=0;
	for(int i=0;i<Data.length();i++)
		{tmp=Data[i];
		if(tmp.compare(delimiter)==0 && Counter<numPnts)
			{Output[Counter]=strtod(bld.c_str(),NULL);
			Counter++;
			bld="";}
		else{bld+=Data[i];}
		}
	return Output;}

int* Solution::fillIntArray(string Data,int numPnts,string delimiter)
	{int* Output=new int[numPnts];
	string bld="",tmp="";
	int Counter=0;
	for(int i=0;i<Data.length();i++)
		{tmp=Data[i];
		if(tmp.compare(delimiter)==0 && Counter<numPnts)
			{Output[Counter]=atoi(bld.c_str());
			Counter++;
			bld="";}
		else{bld+=Data[i];}
		}
	return Output;}

string* Solution::fillStringArray(string Data,int numPnts,string delimiter){string* Output=new string[numPnts];string bld="",tmp="";int Counter=0;for(int i=0;i<Data.length();i++){tmp=Data[i];if(tmp.compare(delimiter)==0 && Counter<numPnts){Output[Counter]=bld;Counter++;bld="";}else{bld+=Data[i];}}return Output;}

// convert mole fraction (x) [mol comp/mol total] to mass fraction (w) [g comp/g total]
// Molecular Weight (MW) [g/mol]
double* Solution::moleFrac2MassFrac(double* x,double* MW,int numComps)
	{double* w=new double[numComps];
	double Sum=0;
	for(int i=0;i<numComps;i++){Sum+=x[i]*MW[i];}
	for(int i=0;i<numComps;i++){w[i]=x[i]*MW[i]/Sum;}
	return w;}

// convert mass fraction (x) [g comp/g total] to mole fraction (w) [mol comp/mol total]
// Molecular Weight (MW) [g/mol]
double* Solution::massFrac2MoleFrac(double* w,double* MW,int numComps)
	{double* x=new double[numComps];
	double Sum=0;
	for(int i=0;i<numComps;i++){Sum+=w[i]/MW[i];}
	for(int i=0;i<numComps;i++){x[i]=(w[i]/MW[i])/Sum;}
	return x;}

double* Solution::massFrac2MoleFrac(string* sW,double* MW,int numComps)
	{double* x=new double[numComps];
	// Convert String to Double
	string tmp;
	double* w=new double[numComps];
	for(int i=0;i<numComps;i++)
		{tmp=sW[i];
		w[i]=strtod(tmp.c_str(),NULL);}
	double Sum=0;
	for(int i=0;i<numComps;i++){Sum+=w[i]/MW[i];}
	for(int i=0;i<numComps;i++){x[i]=(w[i]/MW[i])/Sum;}
	return x;}

// Density [g/L]
// Molecular Weight (MW) [g/mol]
double* Solution::volFrac2MoleFrac(double* v,double* MW,double* density,int numComps)
	{double* x=new double[numComps];
	double Sum=0;
	for(int i=0;i<numComps;i++){Sum+=v[i]*density[i]/MW[i];}
	for(int i=0;i<numComps;i++){x[i]=(v[i]*density[i]/MW[i])/Sum;}
	return x;}

double* Solution::volFrac2MoleFrac(string* sV,double* MW,double* density,int numComps)
	{double* x=new double[numComps];
	// Convert String to Double
	string tmp;
	double* v=new double[numComps];
	for(int i=0;i<numComps;i++)
		{tmp=sV[i];
		v[i]=strtod(tmp.c_str(),NULL);}
	double Sum=0;
	for(int i=0;i<numComps;i++){Sum+=v[i]*density[i]/MW[i];}
	for(int i=0;i<numComps;i++){x[i]=(v[i]*density[i]/MW[i])/Sum;}
	return x;}

// convert molarity (M) [mol solute/total solvent volume] to mass fraction (w) [g comp/g total]
// Molecular Weight (MW) [g/mol]
// solventDensity [g/L]
double* Solution::molarity2MassFrac(double* M,double* MW,int numComps,double solventDensity)
	{double* w=new double[numComps];
	double Sum=solventDensity;
	/* Estimate Solution Density by sum of mass concentrations of each component*/
	for(int i=0;i<numComps;i++){Sum+=MW[i]*M[i];}
	/* Calculate mass fraction of each component*/
	for(int i=0;i<numComps;i++){w[i]=MW[i]*M[i]/Sum;}
	return w;}

double Solution::calcMoleculeMolecularWeight(string structure)
	{double Output=0;
	// Assess Moleculalur Structural Formula
	string numAtomList="",atomList="",atomNm="",rpt="",branch="",delimiter=GLOBAL_DELIMITER,tmp;
	char letter,second_letter;
	bool PASS,SEARCHING;
	string *Atom;
	int *numAtom;
	int pos,aCount=0,multiplier,Counter,numRepeats;
	// Find Number of Each Atom and Atoms Present
	for(int i=0;i<structure.length();i++)
		{// Get Atom Name
		PASS=false;atomNm="";multiplier=1;
		letter=structure[i];
		if(isupper(letter))
			{atomNm=letter;
			// Check if second letter completes atom name or gives number of atoms
			if(i+1<=structure.length()-1)
				{second_letter=structure[i+1];
				if(islower(second_letter)){atomNm+=second_letter;}}
			// Check if number behind Atom Name
			if(i+atomNm.length()<=structure.length()-1)
				{SEARCHING=true;Counter=0;tmp="";
				while(SEARCHING)
					{letter=structure[i+atomNm.length()+Counter];
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
		else if(strcmp(&letter,"(")==0)
			{// Chain Repeat Identified					345678
			// i=position of (, find end of repeat )Num e.g. (CH2)2
			pos=structure.find(")",i);
			// Acquire Repeating Atoms Unit
			rpt=structure.substr(i+1,pos-i-1);
			// Acquire Number of Repeats of Unit
			SEARCHING=true;Counter=0;tmp="";
			while(SEARCHING)
				{letter=structure[i+rpt.length()+2+Counter];
				if(isdigit(letter)){tmp+=letter;Counter++;}
				else{SEARCHING=false;break;}}
			if(tmp.compare("")==0){numRepeats=1;}
			else{numRepeats=atoi(tmp.c_str());}
			// Advance Loop Forward
			//cout<<"i="<<i<<endl;
			i+=1+rpt.length()+tmp.length();
			//cout<<"i="<<i<<endl;
			//cout<<numRepeats<<endl;
			// Determine and Count Atoms in Repeat Unit
			for(int j=0;j<rpt.length();j++)
				{// Get Atom Name
				PASS=false;atomNm="";
				letter=rpt[j];
				if(isupper(letter))
					{atomNm=letter;
					// Check if second letter completes atom name or gives number of atoms
					if(j+1<=rpt.length()-1){second_letter=rpt[j+1];if(islower(second_letter)){atomNm+=second_letter;j++;}}
					//cout<<j<<" "<<atomNm<<endl;
					// Check if number behind Atom Name
					if(j+atomNm.length()<=rpt.length()-1)
						{SEARCHING=true;Counter=0;tmp="";
						while(SEARCHING)
							{letter=rpt[j+atomNm.length()+Counter];
							if(isdigit(letter)){tmp+=letter;Counter++;}
							else{SEARCHING=false;break;}}
						if(tmp.compare("")==0){multiplier=1;}
						else{multiplier=atoi(tmp.c_str());}
						j+=Counter;}
					// Check is atom is already in Atom List, if not add
					if(atomList.compare("")==0)
						{// Initialize Atoms Present Lists
						atomList=atomNm+delimiter;
						numAtomList=cnvrtNumToStrng(multiplier*numRepeats,0)+delimiter;
						aCount++;
						Atom=fill_string_array(atomList,aCount,delimiter);
						numAtom=fill_int_array(numAtomList,aCount,delimiter);}
					else{// Check is atom is already in Atom List, if not add
						PASS=false;
						for(pos=0;pos<aCount;pos++){tmp=Atom[pos];if(tmp.compare(atomNm)==0){PASS=true;break;}}
						if(!PASS)
							{// Atom not found in list, so add
							atomList+=atomNm+delimiter;
							numAtomList+=cnvrtNumToStrng(multiplier*numRepeats,0)+delimiter;
							delete [] Atom;delete [] numAtom;
							aCount++;
							Atom=fill_string_array(atomList,aCount,delimiter);
							numAtom=fill_int_array(numAtomList,aCount,delimiter);}
						else{// Atom Found, so update numAtomList
							numAtom[pos]=numAtom[pos]+multiplier*numRepeats;
							numAtomList=array2String(numAtom,aCount,delimiter);}}}
				}
			}
		else if(strcmp(&letter,"[")==0)
			{// ?4? Molecular Branch Identified
			// i=position of [, find end of branch] e.g. [CH3]
			pos=structure.find("]",i);
			// Acquire Branching Atoms
			branch=structure.substr(i+1,pos-i-1);
			// Advance Loop Forward
			i+=1+branch.length();
			// Determine and Count Atoms in Branch
			for(int j=0;j<branch.length();j++)
				{// Check for Internal Repeating Unit in Branch
				letter=branch[j];
				if(strcmp(&letter,"(")==0)
					{// Chain Repeat Identified
					// j=position of (, find end of repeat )Num e.g. (CH2)2
					pos=branch.find(")",j);
					// Acquire Repeating Atoms Unit
					rpt=branch.substr(j+1,pos-j-1);
					// Acquire Number of Repeats of Unit
					SEARCHING=true;Counter=0;tmp="";
					while(SEARCHING)
						{letter=branch[j+rpt.length()+2+Counter];
						if(isdigit(letter)){tmp+=letter;Counter++;}
						else{SEARCHING=false;break;}}
					if(tmp.compare("")==0){numRepeats=1;}
					else{numRepeats=atoi(tmp.c_str());}
					// Advance Loop Forward
					j+=1+rpt.length()+tmp.length();
					// Determine and Count Atoms in Repeat Unit
					for(int k=0;k<rpt.length();k++)
						{// Get Atom Name
						PASS=false;atomNm="";
						letter=rpt[k];
						if(isupper(letter))
							{atomNm=letter;
							// Check if second letter completes atom name or gives number of atoms
							if(k+1<=rpt.length()-1){second_letter=rpt[k+1];if(islower(second_letter)){atomNm+=second_letter;k++;}}
							//cout<<j<<" "<<atomNm<<endl;
							// Check if number behind Atom Name
							if(k+atomNm.length()<=rpt.length()-1)
								{SEARCHING=true;Counter=0;tmp="";
								while(SEARCHING)
									{letter=rpt[k+atomNm.length()+Counter];
									if(isdigit(letter)){tmp+=letter;Counter++;}
									else{SEARCHING=false;break;}}
								if(tmp.compare("")==0){multiplier=1;}
								else{multiplier=atoi(tmp.c_str());}
								k+=Counter;}
							// Check is atom is already in Atom List, if not add
							if(atomList.compare("")==0)
								{// Initialize Atoms Present Lists
								atomList=atomNm+delimiter;
								numAtomList=cnvrtNumToStrng(multiplier*numRepeats,0)+delimiter;
								aCount++;
								Atom=fill_string_array(atomList,aCount,delimiter);
								numAtom=fill_int_array(numAtomList,aCount,delimiter);}
							else{// Check is atom is already in Atom List, if not add
								PASS=false;
								for(pos=0;pos<aCount;pos++){tmp=Atom[pos];if(tmp.compare(atomNm)==0){PASS=true;break;}}
								if(!PASS)
									{// Atom not found in list, so add
									atomList+=atomNm+delimiter;
									numAtomList+=cnvrtNumToStrng(multiplier*numRepeats,0)+delimiter;
									delete [] Atom;delete [] numAtom;
									aCount++;
									Atom=fill_string_array(atomList,aCount,delimiter);
									numAtom=fill_int_array(numAtomList,aCount,delimiter);}
								else{// Atom Found, so update numAtomList
									numAtom[pos]=numAtom[pos]+multiplier*numRepeats;
									numAtomList=array2String(numAtom,aCount,delimiter);}}}}
					}
				else if(isupper(letter))
					{atomNm=letter;
					// Check if second letter completes atom name or gives number of atoms
					if(j+1<=branch.length()-1)
						{second_letter=branch[j+1];
						if(islower(second_letter)){atomNm+=second_letter;j++;}}
					// Check if number behind Atom Name
					if(j+atomNm.length()<=branch.length()-1)
						{SEARCHING=true;Counter=0;tmp="";
						while(SEARCHING)
							{letter=branch[j+atomNm.length()+Counter];
							if(isdigit(letter)){tmp+=letter;Counter++;}
							else{SEARCHING=false;break;}}
						if(tmp.compare("")==0){multiplier=1;}
						else{multiplier=atoi(tmp.c_str());}
						j+=Counter;}
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
			}
		}
	//for(int i=0;i<aCount;i++){cout<<Atom[i]<<numAtom[i];}cout<<"\n";
	for(int i=0;i<aCount;i++)
		{// Identify Atom
		tmp=Atom[i];
		for(pos=0;pos<numElements;pos++){if(tmp.compare(atomicSymbol[pos])==0){break;}}
		// Get Atomic Molecular Weight
		Output+=atomicMolecularWeight[pos]*numAtom[i];}
	return Output;}

// Reference: 2007. Model for calculating the viscosity of aqueous solutions
double* Solution::getViscosityCoefficients(string saltName)
	{double v1,v2,v3,v4,v5,v6;
	double* Output=new double[6];
	if(saltName=="NH4NH4SO4" || saltName=="(NH4)2SO4")
		{// Valid at 15-60C for max mass fraction of 0.46
		v1=11.567;v2=2.5709;v3=1.9046;v4=0.0073357;v5=134.87;v6=6.7567;}
	else if(saltName=="AlCl3")
		{// Valid at 25C for max mass fraction of 0.04
		v1=-8335.5; v2=4823.7; v3=18.29; v4=0.0083713; v5=19265; v6=0.37058;}
	else if(saltName=="BaCl2")
		{// Valid at 10-70C for max mass fraction of 0.33
		v1=2.9763; v2=1.104; v3=4.0186; v4=0.0029981; v5=19.538; v6=0.095807;}
	else if(saltName=="CaNO3NO3" || saltName=="Ca(NO3)2")
		{// Valid at -10-100C for max mass fraction of 0.51
		v1=218.21; v2=5.4442; v3=5.0356; v4=0.010467; v5=176320; v6=15.86;}
	else if(saltName=="CaCl2")
		{// Valid at 0-100C for max mass fraction of 0.51
		v1=32.028; v2=0.78792; v3=-1.1495; v4=0.0026995; v5=780860; v6=5.8442;}
	else if(saltName=="CdNO3NO3" || saltName=="Cd(NO3)2")
		{// Valid at 15-55C for max mass fraction of 0.36
		v1=9.2559; v2=1.2793; v3=3.1227; v4=0.0077153; v5=-17.82; v6=24.838;}
	else if(saltName=="CdCl2")
		{// Valid at 25C for max mass fraction of 0.1592
		v1=6.5301; v2=4.0727; v3=0.044872; v4=-0.038783; v5=-0.778; v6=4.064;}
	else if(saltName=="CdSO4")
		{// Valid at 18-75C for max mass fraction of 0.3574
		v1=34.268; v2=2.9192; v3=4.375; v4=0.012318; v5=-27.966; v6=35.285;}
	else if(saltName=="CoCl2")
		{// Valid at 20-50C for max mass fraction of 0.3445
		v1=16.08; v2=2.3042; v3=7.9624; v4=0.0041625; v5=48.864; v6=-0.14434;}
	else if(saltName=="CoSO4")
		{// Valid at 25-75C for max mass fraction of 0.3305
		v1=26.343; v2=2.6995; v3=5.1901; v4=0.012926; v5=-0.099594; v6=1;}
	else if(saltName=="Cr2SO4SO4SO4" || saltName=="Cr2(SO4)3")
		{// Valid at 25C for max mass fraction of 0.0015
		v1=12.673; v2=2.4726; v3=12.007; v4=0.0083713; v5=7.3974; v6=2.0697;}
	else if(saltName=="CrCl3")
		{// Valid at 20-50C for max mass fraction of 0.2535
		v1=-120.2; v2=4.8139; v3=5.655; v4=0.013; v5=-11.341; v6=1.92;}
	else if(saltName=="CuNO3NO3" || saltName=="Cu(NO3)2")
		{// Valid at 25C for max mass fraction of 0.3331
		v1=45.382; v2=1.9172; v3=3.897; v4=0.0083713; v5=9.4796; v6=2.0697;}
	else if(saltName=="CuCl2")
		{// Valid at 20-50C for max mass fraction of 0.3751
		v1=6.9303; v2=1.8668; v3=7.4786; v4=0.0060045; v5=65.035; v6=0.15662;}
	else if(saltName=="CuSO4")
		{// Valid at 15-60C for max mass fraction of 0.4129
		v1=-130.52; v2=4.0392; v3=10.017; v4=0.060995; v5=-2324.8; v6=9.2724;}
	else if(saltName=="Fe2SO4SO4SO4" || saltName=="Fe2(SO4)3")
		{// Valid at 20-50C for max mass fraction of 0.3009
		v1=42.849; v2=2.5502; v3=7.1619; v4=0.011598; v5=1.239; v6=-0.57468;}
	else if(saltName=="FeCl2")
		{// Valid at 18-40C for max mass fraction of 0.0366
		v1=-0.29592; v2=18.533; v3=8.8165; v4=0.0021966; v5=385.52; v6=0.23964;}
	else if(saltName=="FeSO4")
		{// Valid at 25-75C for max mass fraction of 0.2109
		v1=27.074; v2=0.88474; v3=0.86143; v4=0.0051165; v5=7322.9; v6=3.9248;}
	else if(saltName=="H2O2")
		{// Valid at 0-20C for max mass fraction of 1.00
		v1=227.72; v2=0.0009402; v3=-226.89; v4=0.04474; v5=0.26782; v6=3.3224;}
	else if(saltName=="H2SO4")
		{// Valid at -10-75C for max mass fraction of 0.782
		v1=7.777; v2=1.2344; v3=1.9164; v4=0.0056363; v5=46.007; v6=2.8307;}
	else if(saltName=="H3PO4")
		{// Valid at 20-25C for max mass fraction of 0.8188
		v1=5044.2; v2=-0.017169; v3=-5029.7; v4=0.0083713; v5=1497.1; v6=-74.319;}
	else if(saltName=="HCH3COO" || saltName=="HCH3CO2" || saltName=="CH3COOH" || saltName=="CH3CO2H")
		{// Valid at 15-55C for max mass fraction of 1.00
		v1=-0.5661; v2=-0.21327; v3=3.4874; v4=0.074792; v5=1.4946; v6=27.123;}
	else if(saltName=="HCHOO" || saltName=="HCHO2" || saltName=="CHOOH" || saltName=="CHO2H")
		{// Valid at 15-55C for max mass fraction of 1.00
		v1=10.991; v2=0.039543; v3=-6.9252; v4=0.0064144; v5=19.113; v6=0.37167;}
	else if(saltName=="HCl")
		{// Valid at 10-42.5C for max mass fraction of 0.36
		v1=7.124; v2=1.1919; v3=1.6648; v4=0.00096271; v5=22.185; v6=1.479;}
	else if(saltName=="HCN")
		{// Valid at 0C for max mass fraction of 1.00
		v1=1.1509; v2=3.8268; v3=-0.38751; v4=0.0083713; v5=7.385; v6=2.0697;}
	else if(saltName=="HNO3")
		{// Valid at 4-25C for max mass fraction of 0.309
		v1=0.000082264; v2=-1.3972; v3=0.11962; v4=-0.01992; v5=-1.4448; v6=0.90694;}
	else if(saltName=="K2CO3")
		{// Valid at 19-89C for max mass fraction of 0.36
		v1=15.804; v2=1.8414; v3=2.4469; v4=0.0060894; v5=36.771; v6=3.1555;}
	else if(saltName=="K2Cr2O7")
		{// Valid at 0.2 89.3C for max mass fraction of 0.4
		v1=3.2419; v2=1.1035; v3=14.602; v4=0.00024669; v5=5404950; v6=0.26452;}
	else if(saltName=="K2HPO4")
		{// Valid at 20.95-49.95C for max mass fraction of 0.181
		v1=5873.1; v2=5.555; v3=4.3106; v4=0.096537; v5=-0.81881; v6=0.051171;}
	else if(saltName=="K2SO4")
		{// Valid at 0-89.5C for max mass fraction of 0.155
		v1=-983.76; v2=0.00019331; v3=984.52; v4=0.0038473; v5=-9.5001; v6=2.1916;}
	else if(saltName=="K3PO4")
		{// Valid at 19.95-49.95C for max mass fraction of 0.209
		v1=2.3507; v2=0.46935; v3=2.8206; v4=0.0097273; v5=-6.3284; v6=1.9539;}
	else if(saltName=="KBr")
		{// Valid at 0-95C for max mass fraction of 0.452
		v1=187.46; v2=-0.00068486; v3=-188.58; v4=-0.0027457; v5=-1.0079; v6=0.47054;}
	else if(saltName=="KCH3COO" || saltName=="KCH3CO2")
		{// Valid at 15-55C for max mass fraction of 0.625
		v1=6.4061; v2=3.0168; v3=3.6581; v4=0.013116; v5=0.36839; v6=-0.45637;}
	else if(saltName=="KCHOO" || saltName=="KCHO2")
		{// Valid at 15-55C for max mass fraction of 0.813
		v1=4.8875; v2=4.7802; v3=1.2964; v4=0.044035; v5=-1.0602; v6=1.202;}
	else if(saltName=="KCl")
		{// Valid at 5-150C for max mass fraction of 0.306
		v1=6.4883; v2=1.3175; v3=-0.77785; v4=0.092722; v5=-1.3; v6=2.0811;}
	else if(saltName=="KH2PO4")
		{// Valid at 19.95-44.95C for max mass fraction of 0.198
		v1=1358.1; v2=3.8539; v3=1.6617; v4=0.012129; v5=-1.0516; v6=4.1301;}
	else if(saltName=="KI")
		{// Valid at 5-95C for max mass fraction of 0.627
		v1=4.9108; v2=1.1146; v3=-1.2651; v4=0.060535; v5=0.52054; v6=0.09659;}
	else if(saltName=="KNO3")
		{// Valid at 15-60C for max mass fraction of 0.495
		v1=3.9621; v2=0.39973; v3=4.259; v4=0.00076988; v5=925.92; v6=0.33733;}
	else if(saltName=="KOH")
		{// Valid at -14.1-40C for max mass fraction of 0.519
		v1=57.433; v2=3.8049; v3=2.6226; v4=0.010312; v5=5219; v6=10.946;}
	else if(saltName=="Li2SO4")
		{// Valid at 5-128.58C for max mass fraction of 0.26
		v1=42.491; v2=2.2133; v3=6.064; v4=0.0081429; v5=24.167; v6=2.96;}
	else if(saltName=="LiCl")
		{// Valid at 0-95C for max mass fraction of 0.454
		v1=19.265; v2=1.9701; v3=1.6661; v4=0.01149; v5=-0.79691; v6=-0.017456;}
	else if(saltName=="LiNO3")
		{// Valid at 0-110C for max mass fraction of 0.671
		v1=16.905; v2=1.5469; v3=0.63451; v4=0.0083965; v5=156.27; v6=4.3685;}
	else if(saltName=="LiOH")
		{// Valid at 20-40C for max mass fraction of 0.113
		v1=4017.7; v2=3.1299; v3=14.397; v4=0.008346; v5=669.91; v6=2.1575;}
	else if(saltName=="MgNO3NO3" || saltName=="Mg(NO3)2")
		{// Valid at 0-50C for max mass fraction of 0.37
		v1=31.394; v2=2.4805; v3=3.7822; v4=0.0072105; v5=1195.6; v6=31.967;}
	else if(saltName=="MgCl2")
		{// Valid at 15-70C for max mass fraction of 0.386
		v1=24.032; v2=2.2694; v3=3.7108; v4=0.021853; v5=-1.1236; v6=0.14474;}
	else if(saltName=="MgSO4")
		{// Valid at 15-150C for max mass fraction of 0.303
		v1=23.286; v2=0.99267; v3=5.4094; v4=0.0063869; v5=63.574; v6=1.6689;}
	else if(saltName=="MnCl2")
		{// Valid at 25C for max mass fraction of 0.42
		v1=28.172; v2=2.4694; v3=3.8408; v4=0.0083713; v5=7.385; v6=2.0697;}
	else if(saltName=="MnSO4")
		{// Valid at 20-80C for max mass fraction of 0.364
		v1=24.038; v2=2.1341; v3=6.4724; v4=0.0074808; v5=0.67399; v6=-0.53884;}
	else if(saltName=="Na2CO3")
		{// Valid at 20-90C for max mass fraction of 0.308
		v1=16.179; v2=0.48606; v3=1.6033; v4=0.01761; v5=-6.9093; v6=2.2356;}
	else if(saltName=="Na2HPO4")
		{// Valid at 20-50C for max mass fraction of 0.099
		v1=99.53; v2=2.7075; v3=0.23625; v4=-0.017106; v5=-0.99535; v6=0.0060289;}
	else if(saltName=="Na2S2O3")
		{// Valid at 20-50C for max mass fraction of 0.575
		v1=19.509; v2=1.5231; v3=2.7793; v4=0.015159; v5=2.01403e12; v6=49.214;}
	else if(saltName=="Na2SO3")
		{// Valid at 25-40C for max mass fraction of 0.06
		v1=0.000044078; v2=-2.2823; v3=5.5871; v4=0.01463; v5=-0.2519; v6=4.7169;}
	else if(saltName=="Na2SO4")
		{// Valid at 15-150C for max mass fraction of 0.331
		v1=26.519; v2=1.5746; v3=3.4966; v4=0.010388; v5=106.23; v6=2.9738;}
	else if(saltName=="Na3PO4")
		{// Valid at 19.95-49.95C for max mass fraction of 0.076
		v1=24.444; v2=0.63525; v3=12.95; v4=0.0025154; v5=1639040; v6=1.5312;}
	else if(saltName=="NaBr")
		{// Valid at 5-60C for max mass fraction of 0.54
		v1=13.029; v2=1.7478; v3=0.60413; v4=0.010804; v5=17.681; v6=2.3831;}
	else if(saltName=="NaCH3COO" || saltName=="NaCH3CO2")
		{// Valid at 25-55C for max mass fraction of 0.456
		v1=13.239; v2=1.6335; v3=5.6914; v4=0.020441; v5=23994; v6=14.214;}
	else if(saltName=="NaCl")
		{// Valid at 5-154C for max mass fraction of 0.264
		v1=16.222; v2=1.3229; v3=1.4849; v4=0.0074691; v5=30.78; v6=2.0583;}
	else if(saltName=="NaClO3")
		{// Valid at 25-55C for max mass fraction of 0.583
		v1=13.697; v2=0.37882; v3=-1.1782; v4=0.0015768; v5=4065.9; v6=2.1278;}
	else if(saltName=="NaF")
		{// Valid at 5-55C for max mass fraction of 0.032
		v1=158.35; v2=1.6418; v3=5.1181; v4=0.0014974; v5=0.42269; v6=6.1895;}
	else if(saltName=="NaH2PO4")
		{// Valid at 20-50C for max mass fraction of 0.3
		v1=3.9294; v2=0.36692; v3=3.4152; v4=0.0097606; v5=3.2217; v6=0.53556;}
	else if(saltName=="NaI")
		{// Valid at 5-97.83C for max mass fraction of 0.629
		v1=12.353; v2=1.7906; v3=-0.22864; v4=0.0084003; v5=43.458; v6=3.1429;}
	else if(saltName=="NaNO3")
		{// Valid at 10-60C for max mass fraction of 0.552
		v1=5.5367; v2=1.3221; v3=1.4038; v4=0.014242; v5=0.43102; v6=-0.12645;}
	else if(saltName=="NaOH")
		{// Valid at 12.5-70C for max mass fraction of 0.56
		v1=440.2; v2=0.0089764; v3=-423.67; v4=0.015949; v5=107.6; v6=4.6489;}
	else if(saltName=="NH3")
		{// Valid at 19.85-39.85C for max mass fraction of 0.62
		v1=2.2287; v2=1.17; v3=3.4326; v4=0.0085807; v5=54.787; v6=0.84697;}
	else if(saltName=="NH4Cl")
		{// Valid at 10-73.5C for max mass fraction of 0.324
		v1=12.396; v2=1.5039; v3=-1.7756; v4=0.23471; v5=-2.7591; v6=2.8408;}
	else if(saltName=="NH4NO3")
		{// Valid at 15-60C for max mass fraction of 0.785
		v1=3.5287; v2=2.5639; v3=0.78627; v4=0.0095033; v5=0.65; v6=-0.48778;}
	else if(saltName=="NiCl2")
		{// Valid at 20-50C for max mass fraction of 0.379
		v1=272.76; v2=5.5487; v3=7.4913; v4=0.0062831; v5=5.4278; v6=-0.70825;}
	else if(saltName=="NiSO4")
		{// Valid at 15-60C for max mass fraction of 0.353
		v1=26.187; v2=1.887; v3=5.3108; v4=0.010252; v5=-0.31459; v6=0.41465;}
	else if(saltName=="PbNO3NO3" || saltName=="Pb(NO3)2")
		{// Valid at 25-50C for max mass fraction of 0.375
		v1=3939.2; v2=0.00036493; v3=-3932.2; v4=0.0042192; v5=437.27; v6=2.1463;}
	else if(saltName=="SrNO3NO3" || saltName=="Sr(NO3)2")
		{// Valid at 25C for max mass fraction of 0.25
		v1=18.338; v2=1.3569; v3=2.0865; v4=0.0083713; v5=7.385; v6=2.0697;}
	else if(saltName=="SrCl2")
		{// Valid at 10-72.9C for max mass fraction of 0.454
		v1=6.3341; v2=1.6291; v3=1.9197; v4=0.008781; v5=-4.265; v6=2.5302;}
	else if(saltName=="sucrose")
		{// Valid at 15-55C for max mass fraction of 0.25
		v1=59.669; v2=1.7701; v3=10.201; v4=0.010025; v5=2716.5; v6=4.8553;}
	else if(saltName=="ZnCl2")
		{// Valid at 25 C for max mass fraction of 0.52
		v1=12.697; v2=2.8245; v3=2.8195; v4=0.0083713; v5=7.385; v6=2.0697;}
	else if(saltName=="ZnSO4")
		{// Valid at 15-55C for max mass fraction of 0.318
		v1=13.593; v2=1.354; v3=10.556; v4=0.0036034; v5=613.22; v6=0.19583;}
	else{//cerr<<"The salt ("<<saltName<<") does not have any viscosity data associated with it.\nUsing pure solvent viscosity as an estimate...\n";
		v1=0; v2=0; v3=0; v4=0; v5=0; v6=0;}
	Output[0]=v1; Output[1]=v2; Output[2]=v3; Output[3]=v4; Output[4]=v5; Output[5]=v6;
	return Output;}

Solution::Solution(){}

string Solution::interpretMoleculeName(string Molecule)
	{string Output="";
	// Does Molecule Name Input Contain Only Atoms, if so pass
	char letter,second_letter;
	string atomNm="";
	bool PASS,PASS2=true;
	// Determine if molecule name is atomic symbol sequence, structural sequence or some name
	for(int i=0;i<Molecule.length();i++)
		{PASS=false;atomNm="";
		letter=Molecule[i];
		if(isupper(letter))
			{atomNm=letter;
			// Check if second letter completes atom name
			if(i+1<=Molecule.length()-1)
				{second_letter=Molecule[i+1];
				if(islower(second_letter)){atomNm+=second_letter;i++;}}
			// Check if possible atom name is indeed an atom
			for(int j=0;j<numElements;j++){if(atomNm.compare(atomicSymbol[j])==0){PASS=true;break;}}
			if(!PASS){/*Unrecognized atom symbol in molecule name*/PASS2=false;break;}
			Output+=atomNm;}
		else if(isdigit(letter)||strcmp(&letter,":")==0||strcmp(&letter,"[")==0||strcmp(&letter,"]")==0||strcmp(&letter,"(")==0||strcmp(&letter,")")==0)
			{// Number in name can be skipped
			// * <- Indicator of Cyclic compound, can be skipped
			// : <- double (::) or triple (:::) bond indicator, can be skipped
			// [] <- branching chain indicator, can be skipped
			// () <- chain repeat indicator, can be skipped
			Output+=letter;}
		else{/*Unrecognized letter in molecule name*/PASS2=false;break;}}
	if(!PASS2)
		{//Contains Something other than atomic symbols (e.g. molecule name like acetonitrile or ethyl acetate)
		if(Molecule.compare("acetic acid")==0||Molecule.compare("acetic_acid")==0||Molecule.compare("ethanoic acid")==0||Molecule.compare("ethanoic_acid")==0)
			{// Acetic Acid = CH3C[::O]OH
			Output="CH3C[::O]OH";}
		else if(Molecule.compare("acetone")==0||Molecule.compare("propanone")==0)
			{// Acetone = CH3C[CH3]::O
			Output="CH3C[CH3]::O";}
		else if(Molecule.compare("acetonitrile")==0)
			{// Acetonitrile = CH3C:::N
			Output="CH3C:::N";}
		else if(Molecule.compare("butanol")==0||Molecule.compare("butyl alcohol")==0||Molecule.compare("butyl_alcohol")==0)
			{// Butanol = "CH3(CH2)3OH"
			Output="CH3(CH2)3OH";}
		else if(Molecule.compare("butyl methyl ketone")==0||Molecule.compare("butyl_methyl_ketone")==0||Molecule.compare("2-hexanone")==0||Molecule.compare("methyl butyl ketone")==0||Molecule.compare("methyl_butyl_ketone")==0||Molecule.compare("propylacetone")==0)
			{// Butyl Methyl Ketone = CH3(CH2)3C[::O]CH3
			Output="CH3(CH2)3C[::O]CH3";}
		else if(Molecule.compare("chlorobenzene")==0||Molecule.compare("benzene chloride")==0||Molecule.compare("benzene_chloride")==0||Molecule.compare("phenylchloride")==0||Molecule.compare("chlorobenzol")==0)
			{// Chlorobenzene = *1CHCHC[Cl]CHCH*1CH
			Output="*1CHCHC[Cl]CHCH*1CH";}
		else if(Molecule.compare("chloroform")==0||Molecule.compare("trichloromethane")==0)
			{// Chloroform = CHCl3
			Output="CHCl3";}
		else if(Molecule.compare("cyclohexane")==0)
			{// Cyclohexane = *1CH2(CH2)4*1CH2
			Output="*1CH2(CH2)4*1CH2";}
		else if(Molecule.compare("dibutyl ether")==0||Molecule.compare("dibutyl_ether")==0||Molecule.compare("1-butoxybutane")==0)
			{// Dibutyl Ether = CH3(CH2)3O(CH2)3CH3
			Output="CH3(CH2)3O(CH2)3CH3";}
		else if(Molecule.compare("1,2-dichloroethane")==0||Molecule.compare("ethane dichloride")==0||Molecule.compare("ethane_dichloride")==0)
			{// 1,2-Dichloroethane = "C[Cl]H2C[Cl]H2"
			Output="C[Cl]H2C[Cl]H2";}
		else if(Molecule.compare("dichloromethane")==0)
			{// Dichloromethane = ClC[Cl]H2
			Output="ClC[Cl]H2";}
		else if(Molecule.compare("1,2-dimethoxyethane")==0)
			{// 1,2-Dimethoxyethane = CH3OCH2CH2OCH3
			Output="CH3OCH2CH2OCH3";}
		else if(Molecule.compare("n,n-dimethylacetamide")==0)
			{// N,N-Dimethylacetamide = CH3C[::O]N[CH3]CH3
			Output="CH3C[::O]N[CH3]CH3";}
		else if(Molecule.compare("n,n-dimethylformamide")==0)
			{// N,N-dimethylformamide = CH3N[CH3]C[::O]H
			Output="CH3N[CH3]C[::O]H";}
		else if(Molecule.compare("1,4-dioxane")==0)
			{// 1,4-Dioxane = *1CH2OCH2CH2O*1CH2
			Output="*1CH2OCH2CH2O*1CH2";}
		else if(Molecule.compare("dimethyl sulfoxide")==0||Molecule.compare("dimethyl_sulfoxide")==0)
			{// Dimethyl sulfoxide = CH3S[::O]CH3
			Output="CH3S[::O]CH3";}
		else if(Molecule.compare("ethanol")==0||Molecule.compare("ethyl alcohol")==0||Molecule.compare("ethyl_alcohol")==0)
			{// Ethanol = CH3CH2OH
			Output="CH3CH2OH";}
		else if(Molecule.compare("ethyl acetate")==0||Molecule.compare("ethyl_acetate")==0||Molecule.compare("ethyl ethanoate")==0||Molecule.compare("ethyl_ethanoate")==0)
			{// Ethyl acetate = CH3C[::O]OCH2CH3
			Output="CH3C[::O]OCH2CH3";}
		else if(Molecule.compare("diethyl ether")==0||Molecule.compare("diethyl_ether")==0||Molecule.compare("ethyl ether")==0||Molecule.compare("ethyl_ether")==0||Molecule.compare("ethoxyethane")==0)
			{// Ethyl ether = CH3CH2OCH2CH3
			Output="CH3CH2OCH2CH3";}
		else if(Molecule.compare("ethylene glycol")==0||Molecule.compare("ethylene_glycol")==0||Molecule.compare("ethane-1,2-diol")==0)
			{// Ethylene glycol = HOCH2CH2OH
			Output="HOCH2CH2OH";}
		else if(Molecule.compare("2-ethylhexanol")==0)
			{// 2-Ethylhexanol = CH3(CH2)3C[CH2CH3]HCH2OH
			Output="CH3(CH2)3C[CH2CH3]HCH2OH";}
		else if(Molecule.compare("formamide")==0||Molecule.compare("methanamide")==0||Molecule.compare("carbamaldehyde")==0)
			{// Formamide = NH2C[::O]H
			Output="NH2C[::O]H";}
		else if(Molecule.compare("glycerol")==0||Molecule.compare("glycerine")==0||Molecule.compare("glycerin")==0||Molecule.compare("propane-1,2,3-triol")==0)
			{// Glycerol = HOCH2C[OH]HCH2OH
			Output="HOCH2C[OH]HCH2OH";}
		else if(Molecule.compare("hexane")==0)
			{// Hexane = CH3(CH2)4CH3
			Output="CH3(CH2)4CH3";}
		else if(Molecule.compare("isobutanol")==0||Molecule.compare("2-methylpropan-1-ol")==0||Molecule.compare("isobutyl alcohol")==0||Molecule.compare("isobutyl_alcohol")==0)
			{// Isobutanol = CH3C[CH3]HCH2OH
			Output="CH3C[CH3]HCH2OH";}
		else if(Molecule.compare("methanol")==0||Molecule.compare("methyl alcohol")==0||Molecule.compare("methyl_alcohol")==0)
			{// Methanol = CH3OH
			Output="CH3OH";}
		else if(Molecule.compare("2-methoxyethanol")==0)
			{// 2-Methoxyethanol = CH3O(CH2)2OH
			Output="CH3O(CH2)2OH";}
		else if(Molecule.compare("methyl ethyl ketone")==0||Molecule.compare("methyl_ethyl_ketone")==0||Molecule.compare("ethyl methyl ketone")==0||Molecule.compare("ethyl_methyl_ketone")==0||Molecule.compare("butanone")==0)
			{// Methyl ethyl ketone = CH3CH2C[::O]CH3
			Output="CH3CH2C[::O]CH3";}
		else if(Molecule.compare("methyl isopropyl ketone")==0||Molecule.compare("methyl_isopropyl_ketone")==0||Molecule.compare("3-methylbutan-2-one")==0)
			{// Methyl isopropyl ketone = CH3C[CH3]HC[::O]CH3
			Output="CH3C[CH3]HC[::O]CH3";}
		else if(Molecule.compare("nitromethane")==0||Molecule.compare("nitrocarbol")==0)
			{// Nitromethane = CH3N[::O]O
			Output="CH3N[::O]O";}
		else if(Molecule.compare("propylene glycol")==0||Molecule.compare("propylene_glycol")==0||Molecule.compare("propane-1,2-diol")==0)
			{// Propylene glycol = CH3C[OH]HCH2OH
			Output="CH3C[OH]HCH2OH";}
		else if(Molecule.compare("propanol")==0||Molecule.compare("1-propanol")==0||Molecule.compare("propyl alcohol")==0||Molecule.compare("propyl_alcohol")==0)
			{// Propanol = CH3(CH2)2OH
			Output="CH3(CH2)2OH";}
		else if(Molecule.compare("isopropanol")==0||Molecule.compare("2-propanol")==0||Molecule.compare("isopropyl alcohol")==0||Molecule.compare("isopropyl_alcohol")==0)
			{// Isopropanol = CH3C[OH]HCH3
			Output="CH3C[OH]HCH3";}
		else if(Molecule.compare("tetrahydrofuran")==0||Molecule.compare("oxolane")==0)
			{// Tetrahydrofuran = *1CH2(CH2)3*1O
			Output="*1CH2(CH2)3*1O";}
		else if(Molecule.compare("tetralin")==0||Molecule.compare("1,2,3,4-tetrahydronaphthalene")==0)
			{// Tetralin = *1CH(CH)2C[CH2(CH2)2*2CH2]*2C*1CH
			Output="*1CH(CH)2C[CH2(CH2)2*2CH2]*2C*1CH";}
		else if(Molecule.compare("methyl isobutyl ketone")==0||Molecule.compare("methyl_isobutyl_ketone")==0||Molecule.compare("4-methylpentan-2-one")==0)
			{// Methyl isobutyl ketone = CH3C[CH3]HCH2C[::O]CH3
			Output="CH3C[CH3]HCH2C[::O]CH3";}
		else if(Molecule.compare("toluene")==0||Molecule.compare("methyl benzene")==0||Molecule.compare("methyl_benzene")==0)
			{// Toluene = *1CH(CH)4*1C[CH3]
			Output="*1CH(CH)4*1C[CH3]";}
		else if(Molecule.compare("trichloroethene")==0||Molecule.compare("trichloroethylene")==0)
			{// Trichloroethene = ClC[Cl]::C[Cl]H
			Output="ClC[Cl]::C[Cl]H";}
		else if(Molecule.compare("water")==0)
			{// Water = H2O
			Output="H2O";}
		else{cout<<"Error in interpretMoleculeName!!!\nUnrecognized molecule ("<<Molecule<<")"<<endl;
			//exit(EXIT_FAILURE);
			Output=Molecule;}
		}
	return Output;}

void Solution::initializeFunctionGroupInformation()
	{CH3.R=0.9011; CH3.Q=0.848; CH3.a_CH3=0; CH3.a_CH2=0; CH3.a_CH=0; CH3.a_benzene=-114.8; CH3.a_cyclohexane=-115.7; CH3.a_toluene=-115.7; CH3.a_OH=644.6; CH3.a_CH3CHOHCH3=310.7; CH3.a_H2O=1300; CH3.a_phenol=2255; CH3.a_CH3CO=472.6; CH3.a_CHO=158.1; CH3.a_COOH=139.4; CH3.a_CH3COOH=972.4; CH3.a_CH2O=662.1; CH3.a_CHCl2=-243.9; CH3.a_CCl3=7.5; CH3.a_aniline=902.2; CH3.a_CNH2=422.1;
	//
	CH2.R=0.6744; CH2.Q=0.54; CH2.a_CH3=0; CH2.a_CH2=0; CH2.a_CH=0; CH2.a_benzene=-114.8; CH2.a_cyclohexane=-115.7; CH2.a_toluene=-115.7; CH2.a_OH=644.6; CH2.a_CH3CHOHCH3=310.7; CH2.a_H2O=1300; CH2.a_phenol=2255; CH2.a_CH3CO=472.6; CH2.a_CHO=158.1; CH2.a_COOH=139.4; CH2.a_CH3COOH=972.4; CH2.a_CH2O=662.1; CH2.a_CHCl2=-243.9; CH2.a_CCl3=7.5; CH2.a_aniline=902.2; CH2.a_CNH2=422.1;
	//
	CH.R=0.4469; CH.Q=0.228; CH.a_CH3=0; CH.a_CH2=0; CH.a_CH=0; CH.a_benzene=-114.8; CH.a_cyclohexane=-115.7; CH.a_toluene=-115.7; CH.a_OH=644.6; CH.a_CH3CHOHCH3=310.7; CH.a_H2O=1300; CH.a_phenol=2255; CH.a_CH3CO=472.6; CH.a_CHO=158.1; CH.a_COOH=139.4; CH.a_CH3COOH=972.4; CH.a_CH2O=662.1; CH.a_CHCl2=-243.9; CH.a_CCl3=7.5; CH.a_aniline=902.2; CH.a_CNH2=422.1;
	//
	benzene.R=0.5313; benzene.Q=0.4; benzene.a_CH3=156.5; benzene.a_CH2=156.5; benzene.a_CH=156.5; benzene.a_benzene=0; benzene.a_cyclohexane=167; benzene.a_toluene=167; benzene.a_OH=703.9; benzene.a_CH3CHOHCH3=577.3; benzene.a_H2O=859.4; benzene.a_phenol=1649; benzene.a_CH3CO=593.7; benzene.a_CHO=362.3; benzene.a_COOH=461.8; benzene.a_CH3COOH=6; benzene.a_CH2O=32.14; benzene.a_CHCl2=0; benzene.a_CCl3=-231.9; benzene.a_aniline=1.64; benzene.a_CNH2=179.7;
	//
	cyclohexane.R=1.0396; cyclohexane.Q=0.66; cyclohexane.a_CH3=104.4; cyclohexane.a_CH2=104.4; cyclohexane.a_CH=104.4; cyclohexane.a_benzene=-146.8; cyclohexane.a_cyclohexane=0; cyclohexane.a_toluene=0; cyclohexane.a_OH=4000; cyclohexane.a_CH3CHOHCH3=906.8; cyclohexane.a_H2O=5695; cyclohexane.a_phenol=292.6; cyclohexane.a_CH3CO=916.7; cyclohexane.a_CHO=1218; cyclohexane.a_COOH=339.1; cyclohexane.a_CH3COOH=5688; cyclohexane.a_CH2O=213.1; cyclohexane.a_CHCl2=0; cyclohexane.a_CCl3=-12.14; cyclohexane.a_aniline=689.6; cyclohexane.a_CNH2=0;
	//
	toluene.R=1.2663; toluene.Q=0.968; toluene.a_CH3=104.4; toluene.a_CH2=104.4; toluene.a_CH=104.4; toluene.a_benzene=-146.8; toluene.a_cyclohexane=0; toluene.a_toluene=0; toluene.a_OH=4000; toluene.a_CH3CHOHCH3=906.8; toluene.a_H2O=5695; toluene.a_phenol=292.6; toluene.a_CH3CO=916.7; toluene.a_CHO=1218; toluene.a_COOH=339.1; toluene.a_CH3COOH=5688; toluene.a_CH2O=213.1; toluene.a_CHCl2=0; toluene.a_CCl3=-12.14; toluene.a_aniline=689.6; toluene.a_CNH2=0;
	//
	OH.R=1; OH.Q=1.2; OH.a_CH3=328.2; OH.a_CH2=328.2; OH.a_CH=328.2; OH.a_benzene=-9.21; OH.a_cyclohexane=1.27; OH.a_toluene=1.27; OH.a_OH=0; OH.a_CH3CHOHCH3=991.3; OH.a_H2O=28.73; OH.a_phenol=-195.5; OH.a_CH3CO=67.07; OH.a_CHO=1409; OH.a_COOH=-104; OH.a_CH3COOH=195.6; OH.a_CH2O=262.5; OH.a_CHCl2=272.2; OH.a_CCl3=-61.57; OH.a_aniline=-348.2; OH.a_CNH2=-166.8;							
	//
	CH3CHOHCH3.R=3.2491; CH3CHOHCH3.Q=3.124; CH3CHOHCH3.a_CH3=-131.9; CH3CHOHCH3.a_CH2=-131.9; CH3CHOHCH3.a_CH=-131.9; CH3CHOHCH3.a_benzene=-252; CH3CHOHCH3.a_cyclohexane=-273.6; CH3CHOHCH3.a_toluene=-273.6; CH3CHOHCH3.a_OH=-268.8; CH3CHOHCH3.a_CH3CHOHCH3=0; CH3CHOHCH3.a_H2O=5.89; CH3CHOHCH3.a_phenol=-153.2; CH3CHOHCH3.a_CH3CO=353.8; CH3CHOHCH3.a_CHO=-338.6; CH3CHOHCH3.a_COOH=-57.98; CH3CHOHCH3.a_CH3COOH=487.1; CH3CHOHCH3.a_CH2O=1970; CH3CHOHCH3.a_CHCl2=507.8; CH3CHOHCH3.a_CCl3=1544; CH3CHOHCH3.a_aniline=0; CH3CHOHCH3.a_CNH2=0;
	//
	H2O.R=0.92; H2O.Q=1.4; H2O.a_CH3=342.4; H2O.a_CH2=342.4; H2O.a_CH=342.4; H2O.a_benzene=372.8; H2O.a_cyclohexane=203.7; H2O.a_toluene=203.7; H2O.a_OH=-122.4; H2O.a_CH3CHOHCH3=104.9; H2O.a_H2O=0; H2O.a_phenol=344.5; H2O.a_CH3CO=-171.8; H2O.a_CHO=-349.9; H2O.a_COOH=-465.7; H2O.a_CH3COOH=-6.32; H2O.a_CH2O=64.42; H2O.a_CHCl2=370.7; H2O.a_CCl3=356.8; H2O.a_aniline=-109.8; H2O.a_CNH2=385.3;
	//
	phenol.R=0.8952; phenol.Q=0.68; phenol.a_CH3=-159.8; phenol.a_CH2=-159.8; phenol.a_CH=-159.8; phenol.a_benzene=-473.2; phenol.a_cyclohexane=-470.4; phenol.a_toluene=-470.4; phenol.a_OH=-63.15; phenol.a_CH3CHOHCH3=-547.2; phenol.a_H2O=-595.9; phenol.a_phenol=0; phenol.a_CH3CO=-825.7; phenol.a_CHO=0; phenol.a_COOH=0; phenol.a_CH3COOH=-898.3; phenol.a_CH2O=0; phenol.a_CHCl2=0; phenol.a_CCl3=0; phenol.a_aniline=-851.6; phenol.a_CNH2=0;
	//
	CH3CO.R=1.6724; CH3CO.Q=1.488; CH3CO.a_CH3=66.56; CH3CO.a_CH2=66.56; CH3CO.a_CH=66.56; CH3CO.a_benzene=-78.31; CH3CO.a_cyclohexane=-73.87; CH3CO.a_toluene=-73.87; CH3CO.a_OH=216; CH3CO.a_CH3CHOHCH3=-127.6; CH3CO.a_H2O=634.8; CH3CO.a_phenol=-568; CH3CO.a_CH3CO=0; CH3CO.a_CHO=-37.36; CH3CO.a_COOH=1247; CH3CO.a_CH3COOH=258.7; CH3CO.a_CH2O=5.202; CH3CO.a_CHCl2=-301; CH3CO.a_CCl3=12.01; CH3CO.a_aniline=1010; CH3CO.a_CNH2=0;
	//
	CHO.R=0.998; CHO.Q=0.948; CHO.a_CH3=146.1; CHO.a_CH2=146.1; CHO.a_CH=146.1; CHO.a_benzene=-75.3; CHO.a_cyclohexane=223.2; CHO.a_toluene=223.2; CHO.a_OH=-431.3; CHO.a_CH3CHOHCH3=231.4; CHO.a_H2O=623.7; CHO.a_phenol=0; CHO.a_CH3CO=128; CHO.a_CHO=0; CHO.a_COOH=0.75; CHO.a_CH3COOH=-245.8; CHO.a_CH2O=0; CHO.a_CHCl2=0; CHO.a_CCl3=0; CHO.a_aniline=0; CHO.a_CNH2=0;
	//
	COOH.R=1.3013; COOH.Q=1.224; COOH.a_CH3=1744; COOH.a_CH2=1744; COOH.a_CH=1744; COOH.a_benzene=75.49; COOH.a_cyclohexane=147.3; COOH.a_toluene=147.3; COOH.a_OH=118.4; COOH.a_CH3CHOHCH3=349.1; COOH.a_H2O=652.3; COOH.a_phenol=0; COOH.a_CH3CO=-101.3; COOH.a_CHO=1051; COOH.a_COOH=0; COOH.a_CH3COOH=-117.6; COOH.a_CH2O=-96.62; COOH.a_CHCl2=1670; COOH.a_CCl3=48.15; COOH.a_aniline=942.2; COOH.a_CNH2=0;
	//
	CH3COOH.R=1.9031; CH3COOH.Q=1.728; CH3COOH.a_CH3=-320.1; CH3COOH.a_CH2=-320.1; CH3COOH.a_CH-320.1; CH3COOH.a_benzene=114.8; CH3COOH.a_cyclohexane=-170; CH3COOH.a_toluene=-170; CH3COOH.a_OH=180.6; CH3COOH.a_CH3CHOHCH3=-152.8; CH3COOH.a_H2O=385.9; CH3COOH.a_phenol=-337.3; CH3COOH.a_CH3CO=58.84; CH3COOH.a_CHO=1090; CH3COOH.a_COOH=1417; CH3COOH.a_CH3COOH=0; CH3COOH.a_CH2O=-235.7; CH3COOH.a_CHCl2=108.9; CH3COOH.a_CCl3=-209.7; CH3COOH.a_aniline=0; CH3COOH.a_CNH2=0;
	//
	CH2O.R=0.9183; CH2O.Q=0.78; CH2O.a_CH3=1571; CH2O.a_CH2=1571; CH2O.a_CH=1571; CH2O.a_benzene=52.13; CH2O.a_cyclohexane=65.69; CH2O.a_toluene=65.69; CH2O.a_OH=137.1; CH2O.a_CH3CHOHCH3=-218.1; CH2O.a_H2O=212.8; CH2O.a_phenol=0; CH2O.a_CH3CO=52.38; CH2O.a_CHO=0; CH2O.a_COOH=1402; CH2O.a_CH3COOH=461.3; CH2O.a_CH2O=0; CH2O.a_CHCl2=137.8; CH2O.a_CCl3=-154.3; CH2O.a_aniline=0; CH2O.a_CNH2=0;
	//
	CHCl2.R=2.0606; CHCl2.Q=1.684; CHCl2.a_CH3=27.9; CHCl2.a_CH2=27.9; CHCl2.a_CH=27.9; CHCl2.a_benzene=0; CHCl2.a_cyclohexane=0; CHCl2.a_toluene=0; CHCl2.a_OH=669.2; CHCl2.a_CH3CHOHCH3=-401.6; CHCl2.a_H2O=740.4; CHCl2.a_phenol=0; CHCl2.a_CH3CO=550.6; CHCl2.a_CHO=0; CHCl2.a_COOH=437.7; CHCl2.a_CH3COOH=-132.9; CHCl2.a_CH2O=-197.7; CHCl2.a_CHCl2=0; CHCl2.a_CCl3=0; CHCl2.a_aniline=0; CHCl2.a_CNH2=0;
	//
	CCl3.R=2.6401; CCl3.Q=2.184; CCl3.a_CH3=21.23; CCl3.a_CH2=21.23; CCl3.a_CH=21.23; CCl3.a_benzene=288.5; CCl3.a_cyclohexane=33.61; CCl3.a_toluene=33.61; CCl3.a_OH=418.4; CCl3.a_CH3CHOHCH3=-465.7; CCl3.a_H2O=793.2; CCl3.a_phenol=0; CCl3.a_CH3CO=-825.7; CCl3.a_CHO=0; CCl3.a_COOH=370.4; CCl3.a_CH3COOH=-898.3; CCl3.a_CH2O=-20.93; CCl3.a_CHCl2=0; CCl3.a_CCl3=0; CCl3.a_aniline=-75.5; CCl3.a_CNH2=0;
	//
	aniline.R=1.06; aniline.Q=0.816; aniline.a_CH3=175.8; aniline.a_CH2=175.8; aniline.a_CH=175.8; aniline.a_benzene=-218.9; aniline.a_cyclohexane=-15.41; aniline.a_toluene=-15.41; aniline.a_OH=529; aniline.a_CH3CHOHCH3=0; aniline.a_H2O=-239.8; aniline.a_phenol=-860.3; aniline.a_CH3CO=857.7; aniline.a_CHO=0; aniline.a_COOH=681.4; aniline.a_CH3COOH=0; aniline.a_CH2O=0; aniline.a_CHCl2=0; aniline.a_CCl3=-216.3; aniline.a_aniline=0; aniline.a_CNH2=0;
	//
	CNH2.R=1.3692; CNH2.Q=1.236; CNH2.a_CH3=-16.74; CNH2.a_CH2=-16.74; CNH2.a_CH=-16.74; CNH2.a_benzene=-38.64; CNH2.a_cyclohexane=0; CNH2.a_toluene=0; CNH2.a_OH=0; CNH2.a_CH3CHOHCH3=0; CNH2.a_H2O=-527.7; CNH2.a_phenol=0; CNH2.a_CH3CO=0; CNH2.a_CHO=0; CNH2.a_COOH=0; CNH2.a_CH3COOH=0; CNH2.a_CH2O=0; CNH2.a_CHCl2=0; CNH2.a_CCl3=0; CNH2.a_aniline=0; CNH2.a_CNH2=0;
	}

// Reference: Lange's Handbook of Chemistry (Table >=4.9)
double Solution::getBondLength(string atomicBond)
	{// Output Bond Length in Angstroms (1A=0.1nm)
	double Output;
	if(atomicBond.compare("CC")==0){Output=1.54;}
	else if(atomicBond.compare("C::C")==0){Output=1.34;}
	else if(atomicBond.compare("C:::C")==0){Output=1.20;}
	else if(atomicBond.compare("CH")==0||atomicBond.compare("HC")==0){Output=1.09;}
	else if(atomicBond.compare("CF")==0||atomicBond.compare("FC")==0){Output=1.38;}
	else if(atomicBond.compare("CCl")==0||atomicBond.compare("ClC")==0){Output=1.76;}
	else if(atomicBond.compare("CBr")==0||atomicBond.compare("BrC")==0){Output=1.94;}
	else if(atomicBond.compare("CI")==0||atomicBond.compare("IC")==0){Output=2.14;}
	else if(atomicBond.compare("CN")==0||atomicBond.compare("NC")==0){Output=1.47;}
	else if(atomicBond.compare("C::N")==0||atomicBond.compare("N::C")==0){Output=1.32;}
	else if(atomicBond.compare("C:::N")==0||atomicBond.compare("N:::C")==0){Output=1.16;}
	else if(atomicBond.compare("CO")==0||atomicBond.compare("OC")==0){Output=1.43;}
	else if(atomicBond.compare("C::O")==0||atomicBond.compare("O::C")==0){Output=1.13;}
	else if(atomicBond.compare("CSe")==0||atomicBond.compare("SeC")==0){Output=1.98;}
	else if(atomicBond.compare("C::Se")==0||atomicBond.compare("Se::C")==0){Output=1.71;}
	else if(atomicBond.compare("CSi")==0||atomicBond.compare("SiC")==0){Output=1.87;}
	else if(atomicBond.compare("CS")==0||atomicBond.compare("SC")==0){Output=1.82;}
	else if(atomicBond.compare("C::S")==0||atomicBond.compare("S::C")==0){Output=1.71;}
	else if(atomicBond.compare("CAl")==0||atomicBond.compare("AlC")==0){Output=2.24;}
	else if(atomicBond.compare("CAs")==0||atomicBond.compare("AsC")==0){Output=1.98;}
	else if(atomicBond.compare("CB")==0||atomicBond.compare("BC")==0){Output=1.56;}
	else if(atomicBond.compare("CBe")==0||atomicBond.compare("BeC")==0){Output=1.93;}
	else if(atomicBond.compare("CBi")==0||atomicBond.compare("BiC")==0){Output=2.30;}
	else if(atomicBond.compare("CCo")==0||atomicBond.compare("CoC")==0){Output=1.83;}
	else if(atomicBond.compare("CCr")==0||atomicBond.compare("CrC")==0){Output=1.92;}
	else if(atomicBond.compare("CFe")==0||atomicBond.compare("FeC")==0){Output=1.84;}
	else if(atomicBond.compare("CGe")==0||atomicBond.compare("GeC")==0){Output=1.93;}
	else if(atomicBond.compare("CHg")==0||atomicBond.compare("HgC")==0){Output=2.07;}
	else if(atomicBond.compare("CIn")==0||atomicBond.compare("InC")==0){Output=2.16;}
	else if(atomicBond.compare("CMo")==0||atomicBond.compare("MoC")==0){Output=2.08;}
	else if(atomicBond.compare("CNi")==0||atomicBond.compare("NiC")==0){Output=2.11;}
	else if(atomicBond.compare("CPb")==0||atomicBond.compare("PbC")==0){Output=2.30;}
	else if(atomicBond.compare("CPd")==0||atomicBond.compare("PdC")==0){Output=2.27;}
	else if(atomicBond.compare("CSb")==0||atomicBond.compare("SbC")==0){Output=2.20;}
	else if(atomicBond.compare("CSn")==0||atomicBond.compare("SnC")==0){Output=2.14;}
	else if(atomicBond.compare("CTe")==0||atomicBond.compare("TeC")==0){Output=1.90;}
	else if(atomicBond.compare("CTl")==0||atomicBond.compare("TlC")==0){Output=2.70;}
	else if(atomicBond.compare("CW")==0||atomicBond.compare("WC")==0){Output=2.06;}
	else if(atomicBond.compare("BB")==0){Output=1.77;}
	else if(atomicBond.compare("BBr")==0||atomicBond.compare("BrB")==0){Output=1.87;}
	else if(atomicBond.compare("BCl")==0||atomicBond.compare("ClB")==0){Output=1.72;}
	else if(atomicBond.compare("BF")==0||atomicBond.compare("FB")==0){Output=1.29;}
	else if(atomicBond.compare("BH")==0||atomicBond.compare("HB")==0){Output=1.21;}
	else if(atomicBond.compare("BN")==0||atomicBond.compare("NB")==0){Output=1.42;}
	else if(atomicBond.compare("BO")==0||atomicBond.compare("OB")==0){Output=1.36;}
	else if(atomicBond.compare("HAl")==0||atomicBond.compare("AlH")==0){Output=1.65;}
	else if(atomicBond.compare("HAs")==0||atomicBond.compare("AsH")==0){Output=1.52;}
	else if(atomicBond.compare("HBe")==0||atomicBond.compare("BeH")==0){Output=1.34;}
	else if(atomicBond.compare("HBr")==0||atomicBond.compare("BrH")==0){Output=1.41;}
	else if(atomicBond.compare("HCa")==0||atomicBond.compare("CaH")==0){Output=2.00;}
	else if(atomicBond.compare("HCl")==0||atomicBond.compare("ClH")==0){Output=1.27;}
	else if(atomicBond.compare("HF")==0||atomicBond.compare("FH")==0){Output=0.92;}
	else if(atomicBond.compare("HGe")==0||atomicBond.compare("GeH")==0){Output=1.53;}
	else if(atomicBond.compare("HI")==0||atomicBond.compare("IH")==0){Output=1.61;}
	else if(atomicBond.compare("HK")==0||atomicBond.compare("KH")==0){Output=2.24;}
	else if(atomicBond.compare("HLi")==0||atomicBond.compare("LiH")==0){Output=1.60;}
	else if(atomicBond.compare("HMg")==0||atomicBond.compare("MgH")==0){Output=1.73;}
	else if(atomicBond.compare("HNa")==0||atomicBond.compare("NaH")==0){Output=1.89;}
	else if(atomicBond.compare("HSb")==0||atomicBond.compare("SbH")==0){Output=1.71;}
	else if(atomicBond.compare("HSe")==0||atomicBond.compare("SeH")==0){Output=1.46;}
	else if(atomicBond.compare("HSn")==0||atomicBond.compare("SnH")==0){Output=1.70;}
	else if(atomicBond.compare("NN")==0){Output=1.02;}
	else if(atomicBond.compare("NCl")==0||atomicBond.compare("ClN")==0){Output=1.79;}
	else if(atomicBond.compare("NF")==0||atomicBond.compare("FN")==0){Output=1.36;}
	else if(atomicBond.compare("NH")==0||atomicBond.compare("HN")==0){Output=1.02;}
	else if(atomicBond.compare("NO")==0||atomicBond.compare("ON")==0){Output=1.24;}
	else if(atomicBond.compare("N::O")==0||atomicBond.compare("O::N")==0){Output=1.19;}
	else if(atomicBond.compare("NSi")==0||atomicBond.compare("SiN")==0){Output=1.57;}
	else if(atomicBond.compare("OH")==0||atomicBond.compare("HO")==0){Output=0.96;}
	else if(atomicBond.compare("OO")==0){Output=1.48;}
	else if(atomicBond.compare("OAl")==0||atomicBond.compare("AlO")==0){Output=1.62;}
	else if(atomicBond.compare("OAs")==0||atomicBond.compare("AsO")==0){Output=1.79;}
	else if(atomicBond.compare("OBa")==0||atomicBond.compare("BaO")==0){Output=1.90;}
	else if(atomicBond.compare("OCl")==0||atomicBond.compare("ClO")==0){Output=1.48;}
	else if(atomicBond.compare("OMg")==0||atomicBond.compare("MgO")==0){Output=1.75;}
	else if(atomicBond.compare("OOs")==0||atomicBond.compare("OsO")==0){Output=1.66;}
	else if(atomicBond.compare("OPb")==0||atomicBond.compare("PbO")==0){Output=1.93;}
	else if(atomicBond.compare("PBr")==0||atomicBond.compare("BrP")==0){Output=2.23;}
	else if(atomicBond.compare("PCl")==0||atomicBond.compare("ClP")==0){Output=2.00;}
	else if(atomicBond.compare("PF")==0||atomicBond.compare("FP")==0){Output=1.55;}
	else if(atomicBond.compare("PH")==0||atomicBond.compare("HP")==0){Output=1.42;}
	else if(atomicBond.compare("PI")==0||atomicBond.compare("IP")==0){Output=2.52;}
	else if(atomicBond.compare("PN")==0||atomicBond.compare("NP")==0){Output=1.49;}
	else if(atomicBond.compare("PO")==0||atomicBond.compare("OP")==0){Output=1.45;}
	else if(atomicBond.compare("PS")==0||atomicBond.compare("SP")==0){Output=2.12;}
	else if(atomicBond.compare("PC")==0||atomicBond.compare("CP")==0){Output=1.56;}
	else if(atomicBond.compare("SiBr")==0||atomicBond.compare("BrSi")==0){Output=2.16;}
	else if(atomicBond.compare("SiCl")==0||atomicBond.compare("ClSi")==0){Output=2.02;}
	else if(atomicBond.compare("SiF")==0||atomicBond.compare("FSi")==0){Output=1.56;}
	else if(atomicBond.compare("SiH")==0||atomicBond.compare("HSi")==0){Output=1.48;}
	else if(atomicBond.compare("SiI")==0||atomicBond.compare("ISi")==0){Output=2.34;}
	else if(atomicBond.compare("SiO")==0||atomicBond.compare("OSi")==0){Output=1.53;}
	else if(atomicBond.compare("SiSi")==0){Output=2.30;}
	else if(atomicBond.compare("SBr")==0||atomicBond.compare("BrS")==0){Output=2.27;}
	else if(atomicBond.compare("SCl")==0||atomicBond.compare("ClS")==0){Output=1.59;}
	else if(atomicBond.compare("SF")==0||atomicBond.compare("FS")==0){Output=1.59;}
	else if(atomicBond.compare("SH")==0||atomicBond.compare("HS")==0){Output=1.33;}
	else if(atomicBond.compare("SO")==0||atomicBond.compare("OS")==0){Output=1.43;}
	else if(atomicBond.compare("SS")==0){Output=2.05;}
	else{cout<<"Error in Solution::getBondLength!!!\nUnrecognized Bond ("<<atomicBond<<")\n";exit(EXIT_FAILURE);}
	return Output;}

string Solution::interpretSolute(string solute)
	{string Output="",ionList="",chargeList="",radiusList="",delimiter=GLOBAL_DELIMITER;
	// Generate List of Possible Ion(s) contained in solute
	int Counter=0,pos,multiplier,pos2;	
	for(int i=0;i<numIons;i++)
		{pos=solute.find(modeledIons[i],0);
		if(pos!=string::npos)
			{// Add to list of culprits
			ionList+=modeledIons[i]+delimiter;Counter++;
			chargeList+=modeledIonsCharge[i]+delimiter;
			radiusList+=modeledIonsRadius[i]+delimiter;}}
	if(ionList.length()==0){cerr<<"Error in Solution::interpretSolute!\nNo ions identified in solute ("<<solute<<").\n";exit(EXIT_FAILURE);}
	string* iList=fill_string_array(ionList,Counter,delimiter);
	string* cList=fill_string_array(chargeList,Counter,delimiter);
	string* rList=fill_string_array(radiusList,Counter,delimiter);
	// Get Largest Culprit, which is guaranteed to be an ion w/o mistake
	int Max=iList[0].length(),maxPos=0;
	for(int i=1;i<Counter;i++){if(iList[i].length()>Max){maxPos=i;Max=iList[i].length();}}
	string largest=iList[maxPos];
	// If largest is transition metal, iList will hold a duplicate accounting for other possible oxidation states
	int Counter2=1;
	for(int i=0;i<Counter;i++){if(i!=maxPos && largest.compare(iList[i])==0){Counter2++;}}
	// Check if Largest Culprit is a Transition Metal
	if(Counter2>1)
		{// Largest Culprit is transition metal, reset maxPos to analyze next largest ion
		Max=0;
		for(int i=0;i<Counter;i++)
			{if(iList[i].length()>Max && largest.compare(iList[i])!=0){maxPos=i;Max=iList[i].length();}}
		}
	// Update Output charge valence of first ion
	Output=cList[maxPos]+",";
	// Check for any number(s) following largest ion
	pos=solute.find(largest.c_str(),0);
	bool SEARCHING=true;Counter=0;
	string tmp="";
	char letter;
	while(SEARCHING)
		{letter=solute[pos+largest.length()+Counter];
		if(isdigit(letter)){tmp+=letter;Counter++;}
		else{SEARCHING=false;break;}}
	if(tmp.compare("")==0){multiplier=1;}
	else{multiplier=atoi(tmp.c_str());}
	// Check for preceding and/or finishing parantheses
	int numRepeats=0;
	if(pos-1>=0)
		{letter=solute[pos-1];
		if(strcmp(&letter,"(")==0)
			{// Preceding Parantheses found, find end
			pos2=solute.find(")",pos-1);
			if(pos2==string::npos){cerr<<"Error in Solution::interpretSolute!\nUnbalanced parantheses in solute ("<<solute<<").\n";exit(EXIT_FAILURE);}
			// Acquire Number of Repeats of Unit
			SEARCHING=true;Counter=0;tmp="";
			while(SEARCHING)
				{letter=solute[pos2+1+Counter];
				if(isdigit(letter)){tmp+=letter;Counter++;}
				else{SEARCHING=false;break;}}
			if(tmp.compare("")==0){numRepeats=1;}
			else{numRepeats=atoi(tmp.c_str());}}
		}
	// Update Output Concentration Multiplier
	if(numRepeats==0){Output+=cnvrtNumToStrng(multiplier,0)+",";}
	else if(numRepeats>0){Output+=cnvrtNumToStrng(multiplier*numRepeats,0)+",";}
	else{cerr<<"Error in Solution::interpretSolute!\nBad parantheses in solute ("<<solute<<").\n";exit(EXIT_FAILURE);}
	// Update Output: Ionic Radius (A)
	Output+=rList[maxPos]+",";

	// Assess Parantheses to define remaining ion(s)
	string remainder="";
	int multiplierWidth;
	if(multiplier!=1){multiplierWidth=floor(log10(multiplier))+1;}
	else if(multiplier==1){multiplierWidth=0;}
	if(numRepeats==0)
		{// No Parantheses found, so simply remove largest ion
		if(pos==0)
			{// largestOtherIons
			remainder=solute.substr(largest.length()+multiplierWidth,solute.length()-largest.length());}
		else if(pos+largest.length()==solute.length())
			{// OtherIonslargest
			remainder=solute.substr(0,pos);}
		else{// OtherIonslargestOtherIons
			remainder=solute.substr(0,pos);
			remainder+=solute.substr(pos+largest.length()+multiplierWidth,solute.length()-pos-largest.length());}
		}
	else if(numRepeats==1)
		{// Parantheses found, but w/o a number
		if(pos==1)
			{// (largest)OtherIons
			remainder=solute.substr(largest.length()+2,solute.length()-largest.length()-2);}
		else if(pos+largest.length()==solute.length()-1)
			{// OtherIons(largest)
			remainder=solute.substr(0,pos-1);}
		else
			{// OtherIons(largest)OtherIons
			remainder=solute.substr(0,pos-1);
			remainder+=solute.substr(pos+largest.length()+1,solute.length()-pos-largest.length()-1);}
		}
	else if(numRepeats>1)
		{// Parantheses found, but w/ a number that is Counter long
		if(pos==0)
			{// (largest)?OtherIons, ?=numRepeats
			remainder=solute.substr(largest.length()+2+Counter,solute.length()-largest.length()-Counter-2);}
		else if(pos+largest.length()==solute.length()-1-Counter)
			{// OtherIons(largest)?
			remainder=solute.substr(0,pos-1);}
		else
			{// OtherIons(largest)?OtherIons
			remainder=solute.substr(0,pos-1);
			remainder+=solute.substr(pos+largest.length()+1+Counter,solute.length()-pos-largest.length()-Counter-1);}
		}
	else{cerr<<"Error in Solution::interpretSolute!\nBad solute ("<<solute<<").\n";exit(EXIT_FAILURE);}
	// Now Assess Remainder
	//cout<<remainder<<endl;
	bool ASSESSING=true;
	string remainder2="";
	while(ASSESSING)
		{ionList="";chargeList="";radiusList="";delete [] iList;delete [] cList;delete [] rList;
		// Generate List of Possible Ion(s) contained in solute
		Counter=0;
		for(int i=0;i<numIons;i++)
			{pos=remainder.find(modeledIons[i],0);
			if(pos!=string::npos)
				{// Add to list of culprits
				ionList+=modeledIons[i]+delimiter;Counter++;
				chargeList+=modeledIonsCharge[i]+delimiter;
				radiusList+=modeledIonsRadius[i]+delimiter;}}
		if(ionList.length()==0){cerr<<"Error in Solution::interpretSolute!\nNo ions identified in solute ("<<remainder<<").\n";exit(EXIT_FAILURE);}
		iList=fill_string_array(ionList,Counter,delimiter);
		cList=fill_string_array(chargeList,Counter,delimiter);
		rList=fill_string_array(radiusList,Counter,delimiter);
		// Get Largest Culprit, which is guaranteed to be an ion w/o mistake
		Max=iList[0].length(),maxPos=0;
		for(int i=1;i<Counter;i++){if(iList[i].length()>Max){maxPos=i;Max=iList[i].length();}}
		largest=iList[maxPos];
		// Update Output charge valence of first ion
		Output+=cList[maxPos]+",";
		// Check for any number(s) following largest ion
		pos=remainder.find(largest.c_str(),0);
		SEARCHING=true;Counter=0;tmp="";
		while(SEARCHING)
			{letter=remainder[pos+largest.length()+Counter];
			if(isdigit(letter)){tmp+=letter;Counter++;}
			else{SEARCHING=false;break;}}
		if(tmp.compare("")==0){multiplier=1;}
		else{multiplier=atoi(tmp.c_str());}
		// Check for preceding and/or finishing parantheses
		numRepeats=0;
		if(pos-1>=0)
			{letter=remainder[pos-1];
			if(strcmp(&letter,"(")==0)
				{// Preceding Parantheses found, find end
				pos2=remainder.find(")",pos-1);
				if(pos2==string::npos){cerr<<"Error in Solution::interpretSolute!\nUnbalanced parantheses in solute ("<<remainder<<").\n";exit(EXIT_FAILURE);}
				// Acquire Number of Repeats of Unit
				SEARCHING=true;Counter=0;tmp="";
				while(SEARCHING)
					{letter=remainder[pos2+1+Counter];
					if(isdigit(letter)){tmp+=letter;Counter++;}
					else{SEARCHING=false;break;}}
				if(tmp.compare("")==0){numRepeats=1;}
				else{numRepeats=atoi(tmp.c_str());}}
			}
		// Update Output Concentration Multiplier
		if(numRepeats==0){Output+=cnvrtNumToStrng(multiplier,0)+",";}
		else if(numRepeats>0){Output+=cnvrtNumToStrng(multiplier*numRepeats,0)+",";}
		else{cerr<<"Error in Solution::interpretSolute!\nBad parantheses in solute ("<<remainder<<").\n";exit(EXIT_FAILURE);}
		// Update Output: Ionic Radius (A)
		Output+=rList[maxPos]+",";
		//cout<<Output<<endl;
		// Assess Parantheses to define remaining ion(s)
		remainder2="";
		if(multiplier!=1){multiplierWidth=floor(log10(multiplier))+1;}
		else if(multiplier==1){multiplierWidth=0;}
		if(numRepeats==0)
			{// No Parantheses found, so simply remove largest ion
			if(pos==0)
				{// largestOtherIons
				remainder2=remainder.substr(largest.length()+multiplierWidth,remainder.length()-largest.length());}
			else if(pos+largest.length()==remainder.length())
				{// OtherIonslargest
				remainder2=remainder.substr(0,pos);}
			else{// OtherIonslargestOtherIons
				remainder2=remainder.substr(0,pos);
				remainder2+=remainder.substr(pos+largest.length()+multiplierWidth,remainder.length()-pos-largest.length());}
			}
		else if(numRepeats==1)
			{// Parantheses found, but w/o a number
			if(pos==1)
				{// (largest)OtherIons
				remainder2=remainder.substr(largest.length()+2,remainder.length()-largest.length()-2);}
			else if(pos+largest.length()==remainder.length()-1)
				{// OtherIons(largest)
				remainder2=remainder.substr(0,pos-1);}
			else
				{// OtherIons(largest)OtherIons
				remainder2=remainder.substr(0,pos-1);
				remainder2+=remainder.substr(pos+largest.length()+1,remainder.length()-pos-largest.length()-1);}
			}
		else if(numRepeats>1)
			{// Parantheses found, but w/ a number that is Counter long
			if(pos==0)
				{// (largest)?OtherIons, ?=numRepeats
				remainder2=remainder.substr(largest.length()+2+Counter,remainder.length()-largest.length()-Counter-2);}
			else if(pos+largest.length()==remainder.length()-1-Counter)
				{// OtherIons(largest)?
				remainder2=remainder.substr(0,pos-1);}
			else
				{// OtherIons(largest)?OtherIons
				remainder2=remainder.substr(0,pos-1);
				remainder2+=remainder.substr(pos+largest.length()+1+Counter,remainder.length()-pos-largest.length()-Counter-1);}
			}
		else{cerr<<"Error in Solution::interpretSolute!\nBad solute ("<<remainder<<").\n";exit(EXIT_FAILURE);}
		// Exit Loop when no more ions remain
		if(remainder2.compare("")==0){ASSESSING=false;break;}
		else{remainder=remainder2;}
		}
	return Output;}

excessVolume Solution::getFunctionGroupInformation(string group)
	{if(group.compare("CH3")==0){return CH3;}
	else if(group.compare("CH2")==0){return CH2;}
	else if(group.compare("CH")==0){return CH;}
	else if(group.compare("*1CH(CH)4*1CH")==0){return benzene;}
	else if(group.compare("*1CH2(CH2)4*1CH2")==0){return cyclohexane;}
	else if(group.compare("*1CH(CH)4*1C[CH3]")==0){return toluene;}
	else if(group.compare("OH")==0){return OH;}
	else if(group.compare("CH3C[OH]HCH3")==0){return CH3CHOHCH3;}
	else if(group.compare("H2O")==0){return H2O;}
	else if(group.compare("*1CH(CH)4*1C[OH]")==0){return phenol;}
	else if(group.compare("CH3C::O")==0){return CH3CO;}
	else if(group.compare("C[::O]H")==0){return CHO;}
	else if(group.compare("C[::O]OH")==0){return COOH;}
	else if(group.compare("CH3C[::O]OH")==0){return CH3COOH;}
	else if(group.compare("CH2O")==0){return CH2O;}
	else if(group.compare("CHCl2")==0){return CHCl2;}
	else if(group.compare("CCl3")==0){return CCl3;}
	else if(group.compare("*1CH(CH)4*1C[NH2]")==0){return aniline;}
	else if(group.compare("CNH2")==0){return CNH2;}
	else
		{cerr<<"Error in Solution::getFunctionGroupInformation()!\nUnrecognized functional group ("<<group<<")\n"; exit(EXIT_FAILURE);}
	}

double Solution::getFunctionGroupInteractionParameter(string group,string group2)
	{excessVolume eV=getFunctionGroupInformation(group);
	if(group2.compare("CH3")==0){return eV.a_CH3;}
	else if(group2.compare("CH2")==0){return eV.a_CH2;}
	else if(group2.compare("CH")==0){return eV.a_CH;}
	else if(group2.compare("*1CH(CH)4*1CH")==0){return eV.a_benzene;}
	else if(group2.compare("*1CH2(CH2)4*1CH2")==0){return eV.a_cyclohexane;}
	else if(group2.compare("*1CH(CH)4*1C[CH3]")==0){return eV.a_toluene;}
	else if(group2.compare("OH")==0){return eV.a_OH;}
	else if(group2.compare("CH3C[OH]HCH3")==0){return eV.a_CH3CHOHCH3;}
	else if(group2.compare("H2O")==0){return eV.a_H2O;}
	else if(group2.compare("*1CH(CH)4*1C[OH]")==0){return eV.a_phenol;}
	else if(group2.compare("CH3C::O")==0){return eV.a_CH3CO;}
	else if(group2.compare("C[::O]H")==0||group2.compare("CHO")==0){return eV.a_CHO;}
	else if(group2.compare("C[::O]OH")==0group2.compare("COOH")==0){return eV.a_COOH;}
	else if(group2.compare("CH3C[::O]OH")==0){return eV.a_CH3COOH;}
	else if(group2.compare("CH2O")==0){return eV.a_CH2O;}
	else if(group2.compare("CHCl2")==0){return eV.a_CHCl2;}
	else if(group2.compare("CCl3")==0){return eV.a_CCl3;}
	else if(group2.compare("*1CH(CH)4*1C[NH2]")==0){return eV.a_aniline;}
	else if(group2.compare("CNH2")==0){return eV.a_CNH2;}
	else
		{cerr<<"Error in Solution::getFunctionGroupInteractionParameter()!\nUnrecognized functional group ("<<group2<<") interacting with group ("<<group<<")\n"; exit(EXIT_FAILURE);}
	}

string Solution::getIonData(string soluteConc,string soluteName,string delimiter)
	{int numSolutes=count_delimiter(soluteName,delimiter);
	int ErrChck=count_delimiter(soluteConc,delimiter);
	if(numSolutes!=ErrChck)
		{cerr<<"Error in Solution::getIonData\nNumber of solutes ("<<soluteName<<") does not match number of salt concentrations ("<<soluteConc<<")"<<endl;exit(EXIT_FAILURE);}
	// IonData Format: charge_valence,molar_concentration,radius_(A),
	int pos,multiplier,nIons;
	string Salt="",Conc="",ionString="",Output="",tmp="";
	double iC,ionSzCompare=0,test;
	for(int i=0;i<numSolutes;i++)
		{// Get Salt Component
		if(i==0)
			{// Acquire First Salt			
			Salt=soluteName.substr(0,soluteName.find(delimiter,0));
			// Acquire First Salt Concentration
			tmp=soluteConc.substr(0,soluteConc.find(delimiter,0));
			Conc=formatNumberString(tmp);}
		else{// Acquire Other Salts
			pos=0;
			for(int j=0;j<i;j++){pos=soluteName.find(delimiter,pos+1);}
			Salt=soluteName.substr(pos+1,soluteName.find(delimiter,pos+1)-pos-1);
			// Acquire Other Salt Concentration
			pos=0;
			for(int j=0;j<i;j++){pos=soluteConc.find(delimiter,pos+1);}
			tmp=soluteConc.substr(pos+1,soluteConc.find(delimiter,pos+1)-pos-1);
			Conc=formatNumberString(tmp);}

		// Generate Ion String for Salt
		ionString=interpretSolute(Salt);		
		nIons=countDelimiter(ionString,",");
		nIons/=3;pos=0;
		// Multiply Concentration Multiplier by Molar Concentration
		for(int j=0;j<nIons;j++)
			{if(j==0){pos=ionString.find(",",pos+1);}
			else{for(int k=0;k<3;k++){pos=ionString.find(",",pos+1);}}
			tmp=ionString.substr(pos+1,ionString.find(",",pos+1)-pos-1);
			multiplier=atoi(tmp.c_str());
			iC=strtod(Conc.c_str(),NULL);
			iC*=multiplier;
			ionString.replace(pos+1,tmp.length(),formatNumberString(cnvrtNumToStrng(iC,BUFFER_SIZE)));}
		Output+=ionString;}
	//largestIon=cnvrtNumToStrng(ionSzCompare,DISTANCE_SIG_FIGS);	// Keep Ion Radius with 2 significant digits
	return Output;}

int Solution::getSolventFunctionGroups(string theSolvent,string &theFunctionGroups)
	{string delimiter=GLOBAL_DELIMITER;
	theFunctionGroups="";
	theSolvent="|"+theSolvent+"|";
	if(theSolvent.compare("|CH3C::[O]OH|")==0||theSolvent.compare("|acetic acid|")==0||theSolvent.compare("|acetic_acid|")==0)		// Acetic Acid
		{theFunctionGroups+="1CH3"+delimiter;
		theFunctionGroups+="1CH3C[::O]OH"+delimiter;
		theFunctionGroups+="1C[::O]OH"+delimiter;
		return 3;}
	else if(theSolvent.compare("|CH3C[CH3]::O|")==0||theSolvent.compare("|acetone|")==0)	// Acetone
		{theFunctionGroups+="2CH3"+delimiter;
		theFunctionGroups+="1C[::O]H"+delimiter;
		return 2;}
	else if(theSolvent.compare("|CH3C:::N|")==0||theSolvent.compare("|acetonitrile|")==0)		// Acetonitrile
		{theFunctionGroups+="1CH3"+delimiter;
		theFunctionGroups+="1CNH2"+delimiter;
		return 2;}
	else if(theSolvent.compare("|CH3(CH2)3OH|")==0||theSolvent.compare("|butanol|")==0)	// Butanol
		{theFunctionGroups+="1CH3"+delimiter;
		theFunctionGroups+="3CH2"+delimiter;
		theFunctionGroups+="1OH"+delimiter;
		return 3;}
	else if(theSolvent.compare("|CH3(CH2)3C[::O]CH3|")==0||theSolvent.compare("|butyl methyl ketone|")==0||theSolvent.compare("|butyl_methyl_ketone|")==0)	// Butyl Methyl Ketone
		{theFunctionGroups+="2CH3"+delimiter;
		theFunctionGroups+="3CH2"+delimiter;
		theFunctionGroups+="1C[::O]H"+delimiter;
		return 3;}
	else if(theSolvent.compare("|*1CHCHC[Cl]CHCH*1CH|")==0||theSolvent.compare("|chlorobenzene|")==0)	// Chlorobenzene
		{theFunctionGroups+="6CH"+delimiter;
		theFunctionGroups+="1CHCl2"+delimiter;
		return 2;}
	else if(theSolvent.compare("|CHCl3|")==0||theSolvent.compare("|chloroform|")==0)			// Chloroform
		{theFunctionGroups+="1CCl3"+delimiter;
		return 1;}
	else if(theSolvent.compare("|*1CH2(CH2)4*1CH2|")==0||theSolvent.compare("|cyclohexane|")==0)	// Cyclohexane
		{theFunctionGroups+="6CH2"+delimiter;
		return 1;}
	else if(theSolvent.compare("|CH3(CH2)3O(CH2)3CH3|")==0||theSolvent.compare("|dibutyl ether|")==0||theSolvent.compare("|dibutyl_ether|")==0)	// Dibutyl Ether
		{theFunctionGroups+="2CH3"+delimiter;
		theFunctionGroups+="6CH2"+delimiter;
		theFunctionGroups+="1CH2O"+delimiter;
		return 3;}
	else if(theSolvent.compare("|C[Cl]H2C[Cl]H2|")==0||theSolvent.compare("|1,2-dichloroethane|")==0)	// 1,2-dichloroethane
		{theFunctionGroups+="2CHCl2"+delimiter;
		theFunctionGroups+="1CH2"+delimiter;
		return 2;}
	else if(theSolvent.compare("|ClC[Cl]H2|")==0||theSolvent.compare("|dichloromethane|")==0)			// dichloromethane
		{theFunctionGroups+="2CHCl2"+delimiter;
		theFunctionGroups+="1CH2"+delimiter;
		return 2;}
	else if(theSolvent.compare("|CH3OCH2CH2OCH3|")==0||theSolvent.compare("|1,2-dimethoxyethane|")==0)			// 1,2-dimethoxyethane
		{theFunctionGroups+="2CH3"+delimiter;
		theFunctionGroups+="2CH2"+delimiter;
		theFunctionGroups+="2CHO"+delimiter;
		return 3;}
	else if(theSolvent.compare("|CH3C[::O]N[CH3]CH3|")==0||theSolvent.compare("|n,n-dimethylacetamide|")==0)			// n,n-dimethylacetamide
		{theFunctionGroups+="3CH3"+delimiter;
		theFunctionGroups+="1CHO"+delimiter;
		theFunctionGroups+="1CNH2"+delimiter;
		return 3;}
	else if(theSolvent.compare("|CH3N[CH3]C[::O]H|")==0||theSolvent.compare("|n,n-dimethylformamide|")==0)			// n,n-dimethylformamide
		{theFunctionGroups+="2CH3"+delimiter;
		theFunctionGroups+="1CHO"+delimiter;
		theFunctionGroups+="1CNH2"+delimiter;
		return 3;}
	else if(theSolvent.compare("|*1CH2OCH2CH2O*1CH2|")==0||theSolvent.compare("|1,4-dioxane|")==0)			// 1,4-dioxane
		{theFunctionGroups+="4CH2"+delimiter;
		theFunctionGroups+="2CHO"+delimiter;
		return 2;}
	else if(theSolvent.compare("|CH3S[::O]CH3|")==0||theSolvent.compare("|dimethyl sulfoxide|")==0||theSolvent.compare("|dimethyl_sulfoxide|")==0)			// dimethyl sulfoxide
		{theFunctionGroups+="2CH3"+delimiter;
		theFunctionGroups+="1CHO"+delimiter;
		return 2;}
	else if(theSolvent.compare("|CH3CH2OH|")==0||theSolvent.compare("|ethanol|")==0)		// Ethanol
		{theFunctionGroups+="1CH3"+delimiter;
		theFunctionGroups+="1CH2"+delimiter;
		theFunctionGroups+="1OH"+delimiter;
		return 3;}
	else if(theSolvent.compare("|CH3C[::O]OCH2CH3|")==0||theSolvent.compare("|ethyl acetate|")==0||theSolvent.compare("|ethyl_acetate|")==0)		// ethyl acetate
		{theFunctionGroups+="2CH3"+delimiter;
		theFunctionGroups+="1CH2"+delimiter;
		theFunctionGroups+="1COOH"+delimiter;
		return 3;}
	else if(theSolvent.compare("|CH3CH2OCH2CH3|")==0||theSolvent.compare("|diethyl ether|")==0||theSolvent.compare("|diethyl_ether|")==0)		// diethyl ether
		{theFunctionGroups+="2CH3"+delimiter;
		theFunctionGroups+="2CH2"+delimiter;
		theFunctionGroups+="1CHO"+delimiter;
		return 3;}
	else if(theSolvent.compare("|HOCH2CH2OH|")==0||theSolvent.compare("|ethylene glycol|")==0||theSolvent.compare("|ethylene_glycol|")==0)		// ethylene_glycol
		{theFunctionGroups+="2CH2"+delimiter;
		theFunctionGroups+="2OH"+delimiter;
		return 2;}
	else if(theSolvent.compare("|CH3(CH2)3C[CH2CH3]HCH2OH|")==0||theSolvent.compare("|2-ethylhexanol|")==0)		// 2-ethylhexanol
		{theFunctionGroups+="2CH3"+delimiter;
		theFunctionGroups+="5CH2"+delimiter;
		theFunctionGroups+="1CH"+delimiter;
		theFunctionGroups+="1OH"+delimiter;
		return 4;}
	else if(theSolvent.compare("|NH2C[::O]H|")==0||theSolvent.compare("|formamide|")==0)		// formamide
		{theFunctionGroups+="1CNH2"+delimiter;
		theFunctionGroups+="1CHO"+delimiter;
		return 2;}
	else if(theSolvent.compare("|HOCH2C[OH]HCH2OH|")==0||theSolvent.compare("|glycerol|")==0)		// glycerol
		{theFunctionGroups+="3OH"+delimiter;
		theFunctionGroups+="2CH2"+delimiter;
		theFunctionGroups+="1CHO"+delimiter;
		return 3;}
	else if(theSolvent.compare("|CH3(CH2)4CH3|")==0||theSolvent.compare("|hexane|")==0)		// hexane
		{theFunctionGroups+="2CH3"+delimiter;
		theFunctionGroups+="4CH2"+delimiter;
		return 2;}
	else if(theSolvent.compare("|CH3C[CH3]HCH2OH|")==0||theSolvent.compare("|isobutanol|")==0)		// isobutanol
		{theFunctionGroups+="2CH3"+delimiter;
		theFunctionGroups+="1CH2"+delimiter;
		theFunctionGroups+="1CH"+delimiter;
		theFunctionGroups+="1OH"+delimiter;
		return 4;}
	else if(theSolvent.compare("|CH3OH|")==0||theSolvent.compare("|methanol|")==0)		// methanol
		{theFunctionGroups+="1CH3"+delimiter;
		theFunctionGroups+="1OH"+delimiter;
		return 2;}
	else if(theSolvent.compare("|CH3O(CH2)2OH|")==0||theSolvent.compare("|2-methoxyethanol|")==0)		// 2-methoxyethanol
		{theFunctionGroups+="1CH3"+delimiter;
		theFunctionGroups+="2CH2"+delimiter;
		theFunctionGroups+="1CHO"+delimiter;
		theFunctionGroups+="1OH"+delimiter;
		return 4;}
	else if(theSolvent.compare("|CH3CH2C[::O]CH3|")==0||theSolvent.compare("|methyl ethyl ketone|")==0||theSolvent.compare("|methyl_ethyl_ketone|")==0)		// methyl ethyl ketone
		{theFunctionGroups+="2CH3"+delimiter;
		theFunctionGroups+="1CH2"+delimiter;
		theFunctionGroups+="1CHO"+delimiter;
		return 3;}
	else if(theSolvent.compare("|CH3C[CH3]HC[::O]CH3|")==0||theSolvent.compare("|methyl isopropyl ketone|")==0||theSolvent.compare("|methyl_isopropyl_ketone|")==0)		// methyl isopropyl ketone
		{theFunctionGroups+="3CH3"+delimiter;
		theFunctionGroups+="1CH"+delimiter;
		theFunctionGroups+="1CHO"+delimiter;
		return 3;}
	else if(theSolvent.compare("|CH3N[::O]O|")==0||theSolvent.compare("|nitromethane|")==0)		// nitromethane
		{theFunctionGroups+="1CH3"+delimiter;
		theFunctionGroups+="1CNH2"+delimiter;
		return 2;}
	else if(theSolvent.compare("|CH3C[OH]HCH2OH|")==0||theSolvent.compare("|propylene glycol|")==0||theSolvent.compare("|propylene_glycol|")==0)		// propylene_glycol
		{theFunctionGroups+="1CH3"+delimiter;
		theFunctionGroups+="2OH"+delimiter;
		theFunctionGroups+="1CH2"+delimiter;
		theFunctionGroups+="1CH"+delimiter;
		return 4;}
	else if(theSolvent.compare("|CH3(CH2)2OH|")==0||theSolvent.compare("|propanol|")==0)		// propanol
		{theFunctionGroups+="1CH3"+delimiter;
		theFunctionGroups+="2CH2"+delimiter;
		theFunctionGroups+="1OH"+delimiter;
		return 3;}
	else if(theSolvent.compare("|CH3C[OH]HCH3|")==0||theSolvent.compare("|isopropanol|")==0)		// isopropanol
		{theFunctionGroups+="2CH3"+delimiter;
		theFunctionGroups+="1CH"+delimiter;
		theFunctionGroups+="1OH"+delimiter;
		return 3;}
	else if(theSolvent.compare("|*1CH2(CH2)3*1O|")==0||theSolvent.compare("|tetrahydrofuran|")==0)		// tetrahydrofuran
		{theFunctionGroups+="3CH2"+delimiter;
		theFunctionGroups+="1CHO"+delimiter;
		return 2;}
	else if(theSolvent.compare("|*1CH(CH)2C[CH2(CH2)2*2CH2]*2C*1CH|")==0||theSolvent.compare("|tetralin|")==0)		// tetralin
		{theFunctionGroups+="6CH"+delimiter;
		theFunctionGroups+="4CH2"+delimiter;
		return 2;}
	else if(theSolvent.compare("|CH3C[CH3]HCH2C[::O]CH3|")==0||theSolvent.compare("|methyl isobutyl ketone|")==0||theSolvent.compare("|methyl_isobutyl_ketone|")==0)		// methyl_isobutyl_ketone
		{theFunctionGroups+="3CH3"+delimiter;
		theFunctionGroups+="1CH"+delimiter;
		theFunctionGroups+="1CH2"+delimiter;
		theFunctionGroups+="1CHO"+delimiter;
		return 4;}
	else if(theSolvent.compare("|*1CH(CH)4*1C[CH3]|")==0||theSolvent.compare("|toluene|")==0)		// toluene
		{theFunctionGroups+="1*1CH(CH)4*1C[CH3]"+delimiter;
		return 1;}
	else if(theSolvent.compare("|ClC[Cl]::C[Cl]H|")==0||theSolvent.compare("|trichloroethene|")==0)		// trichloroethene
		{theFunctionGroups+="2CHCl2"+delimiter;
		return 1;}
	else if(theSolvent.compare("|H2O|")==0||theSolvent.compare("|water|")==0)			// Water
		{theFunctionGroups+="1H2O"+delimiter;
		return 1;}
	else
		{cerr<<"Error in Solution::getSolventFunctionGroups()!\n"<<theSolvent<<" is not recognized as so I am exiting...\n"; exit(EXIT_FAILURE);}
	}

double Solution::calcDebyeLength(double T)
	{// Count Number of Entries
	int nIon=count_delimiter(ionData,",");
	if(nIon%3!=0){cerr<<"ERROR in calcDebyeLength!!!\nIncorrect ion data input!"<<endl;exit(EXIT_FAILURE);}
	// Extract Data
	string *All=fill_string_array(ionData,nIon,",");
	// Number of Ions
	nIon/=3;
	int *ion_val=new int[nIon];
	double *ion_conc=new double[nIon];
	// Extract Information from Ion Data
	int Counter=0;
	for(int i=0;i<3*nIon;i++)
		{if(i%3==0 && i!=0){Counter++;}
		switch(i%3)
			{case 0:
				ion_val[Counter]=atoi(All[i].c_str());
				break;
			case 1:
				ion_conc[Counter]=strtod(All[i].c_str(),NULL);
				break;
			case 2:
				//ion_rad[Counter]=strtod(All[i].c_str(),NULL);
				break;
			default:
				cerr<<"How is this even possible!"<<endl;}}
	// Calculate Solution Properties
	for(int i=0;i<nIon;i++){ion_conc[i]*=1e3;} // Convert from [mol/L] to [mol/meter^3] (or equivalently millimolar)
	double ionic_strength=0;
	for(int i=0;i<nIon;i++){ionic_strength+=ion_conc[i]*ion_val[i]*ion_val[i];}
	ionic_strength/=2;
	double trm1=dielectric*eo*kb*T;
	double trm2=2*ec*ec*Na*ionic_strength;
	double debyeLength=sqrt(trm1/trm2);
	// Clear Dynamically Allocated Memory
	delete [] All; delete [] ion_val; delete [] ion_conc;
	// Convert to Angstroms
	debyeLength*=1e10;
	return debyeLength;}
