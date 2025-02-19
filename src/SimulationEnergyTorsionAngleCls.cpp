// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


void RNASimulation::getoutenetorsioncls(double val, double mean, double lower, double upper, int rangetype, double &diff, int &bdiff){
	double tmp;
	
	if(rangetype == 1){
		if((lower <= val & val <= 360) | (0 <= val & val <= upper)){
			bdiff=0;
			diff=0;
		}else{
			bdiff=1;
			tmp = abs(val - mean);
			tmp = min(tmp, 360 - tmp);
			diff=log(tmp);
		}
		//cout<<to_string(1) + "-" + to_string(val) + "-" + to_string(mean) + "-" + to_string(lower) + "-" + to_string(upper) + "-" + to_string(bdiff) + "-" + to_string(diff)<<endl;
	}else{
		if(lower <= val & val <= upper){
			bdiff=0;
			diff=0;
		}else{
			bdiff=1;
			tmp = abs(val - mean);
			tmp = min(tmp, 360 - tmp);
			diff=log(tmp);
		}
		//cout<<to_string(0) + "-" + to_string(val) + "-" + to_string(mean) + "-" + to_string(lower) + "-" + to_string(upper) + "-" + to_string(bdiff) + "-" + to_string(diff)<<endl;
	}
}

void RNASimulation::fenergyclsdihedral(vector<Residueatoms> decstr,int numseq,double &nene,double &fene, Energyterms &lp,bool flagacc)
{
	int i,j,k,l,n;
	int ti,tj;
	double tdist;
	int istart;
	Point3d tp;
	char aat;
	int ss;
	double totalnn;
	int tmpval1;
	double tmpval2;
	int bdiff, bdiffmin, bdifftmp;
	double diff, diffmin, difftmp;
	bool flag;
	//int estartidx=BLENGTH_ENERGY_NUM+ILENGTH_ENERGY_NUM+INNERANGLE_ENERGY_NUM;

	if(!flagacc)
	{
		for(j=0;j<TORSION_CLS_ENERGY_NUM;j++){
			lp.fn[TORSION_CLS_ENERGY_ST+j]=0;
			lp.nn[TORSION_CLS_ENERGY_ST+j]=0;
			
			for(i=0;i<numseq;i++)
			{
				aat=decstr[i].aaa;

				if(j == 0 & i != 0){
					bdiffmin=10000;
					diffmin=10000.0;
					for(k=0;k<TORSIONLIST[j].size();k++){
						bdifftmp=0; difftmp=0;
						getoutenetorsioncls(decstr[i-1].dih[TorsionClsList[j][0]], TORSIONLIST[j][k][0].mean, TORSIONLIST[j][k][0].lower, TORSIONLIST[j][k][0].upper, TORSIONLIST[j][k][0].rangetype, diff, bdiff);
						bdifftmp+=bdiff;
						difftmp+=diff;

						for(l=1;l<TORSIONLIST[j][k].size();l++)
						{
							getoutenetorsioncls(decstr[i].dih[TorsionClsList[j][l]], TORSIONLIST[j][k][l].mean, TORSIONLIST[j][k][l].lower, TORSIONLIST[j][k][l].upper, TORSIONLIST[j][k][l].rangetype, diff, bdiff);
							bdifftmp+=bdiff;
							difftmp+=diff;
						}

						if(bdiffmin > bdifftmp){
							if(diffmin > difftmp){
								bdiffmin=bdifftmp;
								diffmin=difftmp;
							}
						}
					}
				}else{
					bdiffmin=10000;
					diffmin=10000.0;
					for(k=0;k<TORSIONLIST[j].size();k++){
						bdifftmp=0; difftmp=0;

						for(l=0;l<TORSIONLIST[j][k].size();l++)
						{
							getoutenetorsioncls(decstr[i].dih[TorsionClsList[j][l]], TORSIONLIST[j][k][l].mean, TORSIONLIST[j][k][l].lower, TORSIONLIST[j][k][l].upper, TORSIONLIST[j][k][l].rangetype, diff, bdiff);
							bdifftmp+=bdiff;
							difftmp+=diff;
						}

						if(bdiffmin > bdifftmp){
							if(diffmin > difftmp){
								bdiffmin=bdifftmp;
								diffmin=difftmp;
							}
						}
					}
				}
				
				lp.indf[TORSION_CLS_ENERGY_ST+j][i]=diffmin;
				lp.indn[TORSION_CLS_ENERGY_ST+j][i]=bdiffmin;
				
				lp.fn[TORSION_CLS_ENERGY_ST+j]+=lp.indf[TORSION_CLS_ENERGY_ST+j][i];
				lp.nn[TORSION_CLS_ENERGY_ST+j]+=lp.indn[TORSION_CLS_ENERGY_ST+j][i];
			}
		}
	}
	else
	{
		int iold[TORSION_CLS_ENERGY_NUM],inew[TORSION_CLS_ENERGY_NUM];
		double fold[TORSION_CLS_ENERGY_NUM],fnew[TORSION_CLS_ENERGY_NUM];
		for(i=0;i<TORSION_CLS_ENERGY_NUM;i++){
			iold[i]=0; inew[i]=0; fold[i]=0; fnew[i]=0;
		}
		
		for(j=0;j<TORSION_CLS_ENERGY_NUM;j++){

			for(i=0;i<numseq;i++) 
			{
				flag=false;
				if(lp.flap[TORSIONANGLE_ENERGY_ST+TorsionClsList[j][0]][i-1]) flag=true;
				for(l=1;l<TorsionClsList[j].size();l++)
				{
					if(lp.flap[TORSIONANGLE_ENERGY_ST+TorsionClsList[j][l]][i]) flag=true;
				}

				if(!flag) continue;

				iold[j]+=lp.indn[TORSION_CLS_ENERGY_ST+j][i];
				fold[j]+=lp.indf[TORSION_CLS_ENERGY_ST+j][i];

				aat=decstr[i].aaa;

				if(j == 0 & i != 0){
					bdiffmin=10000;
					diffmin=10000.0;
					for(k=0;k<TORSIONLIST[j].size();k++){
						bdifftmp=0;
						difftmp=0;

						getoutenetorsioncls(decstr[i-1].dih[TorsionClsList[j][0]], TORSIONLIST[j][k][0].mean, TORSIONLIST[j][k][0].lower, TORSIONLIST[j][k][0].upper, TORSIONLIST[j][k][0].rangetype, diff, bdiff);
						bdifftmp+=bdiff;
						difftmp+=diff;

						for(l=1;l<TORSIONLIST[j][k].size();l++)
						{
							getoutenetorsioncls(decstr[i].dih[TorsionClsList[j][l]], TORSIONLIST[j][k][l].mean, TORSIONLIST[j][k][l].lower, TORSIONLIST[j][k][l].upper, TORSIONLIST[j][k][l].rangetype, diff, bdiff);
							bdifftmp+=bdiff;
							difftmp+=diff;
						}

						if(bdiffmin > bdifftmp){
							if(diffmin > difftmp){
								bdiffmin=bdifftmp;
								diffmin=difftmp;
							}
						}
					}
				}else{
					bdiffmin=10000;
					diffmin=10000.0;
					for(k=0;k<TORSIONLIST[j].size();k++){
						bdifftmp=0; difftmp=0;

						for(l=0;l<TORSIONLIST[j][k].size();l++)
						{
							getoutenetorsioncls(decstr[i].dih[TorsionClsList[j][l]], TORSIONLIST[j][k][l].mean, TORSIONLIST[j][k][l].lower, TORSIONLIST[j][k][l].upper, TORSIONLIST[j][k][l].rangetype, diff, bdiff);
							bdifftmp+=bdiff;
							difftmp+=diff;
						}

						if(bdiffmin > bdifftmp){
							if(diffmin > difftmp){
								bdiffmin=bdifftmp;
								diffmin=difftmp;
							}
						}
					}
				}

				lp.indf[TORSION_CLS_ENERGY_ST+j][i]=diffmin;
				lp.indn[TORSION_CLS_ENERGY_ST+j][i]=bdiffmin;

				inew[j]+=lp.indn[TORSION_CLS_ENERGY_ST+j][i];
				fnew[j]+=lp.indf[TORSION_CLS_ENERGY_ST+j][i];
			}
			
			lp.nn[TORSION_CLS_ENERGY_ST+j]+=(inew[j]-iold[j]);
			lp.fn[TORSION_CLS_ENERGY_ST+j]+=(fnew[j]-fold[j]);
		}
	}

	fene=0; nene=0;
	for(j=0;j<TORSION_CLS_ENERGY_NUM;j++){
		fene+=lp.fn[TORSION_CLS_ENERGY_ST+j];
	}
	for(j=0;j<TORSION_CLS_ENERGY_NUM;j++){
		nene+=lp.nn[TORSION_CLS_ENERGY_ST+j];
	}
}





















