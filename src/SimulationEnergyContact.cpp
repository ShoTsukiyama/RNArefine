// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"
 
void get3gcontact8A(double dist, double u, int &enei, double &enef){
	if(dist < 8){
		enef=-u;
		enei=0;
	}
	else if(8 <= dist & dist < 20){
		enef=-0.5*u*(1-sin(M_PI*(dist/12-1.1666666666666667)));
		enei=1;
	}
	else if(20 <= dist & dist < 80){
		enef=0.5*u*(1+sin(M_PI*(dist/60-0.8333333333333334)));
		enei=1;
	}
	else{
		enef=u;
		enei=1;
	}
}

void get3gcontact20A(double dist, double u, int &enei, double &enef){
	if(dist < 20){
		enef=-u;
		enei=0;
	}
	else if(20 <= dist & dist < 40){
		enef=-0.5*u*(1-sin(M_PI*(dist/20-1.5)));
		enei=1;
	}
	else if(40 <= dist & dist < 120){
		enef=0.5*u*(1+sin(M_PI*(dist/80-1)));
		enei=1;
	}
	else{
		enef=u;
		enei=1;
	}
}

void get3gnoncontact8A(double dist, double u, int &enei, double &enef){
	if(dist < 2){
		enef=u;
		enei=1;
	}
	else if(2 <= dist & dist < 6){
		enef=0.5*u*(1-sin(M_PI*(dist/4-1)));
		enei=1;
	}
	else if(6 <= dist & dist < 8){
		enef=-0.5*u*(1+sin(M_PI*(dist/2-3.5)));
		enei=1;
	}
	else{
		enef=-u;
		enei=0;
	}
}

void get3gnoncontact20A(double dist, double u, int &enei, double &enef){
	if(dist < 6){
		enef=u;
		enei=0;
	}
	else if(6 <= dist & dist < 16){
		enef=0.5*u*(1-sin(M_PI*(dist/10-1.1)));
		enei=0;
	}
	else if(16 <= dist & dist < 20){
		enef=-0.5*u*(1+sin(M_PI*(dist/4-4.5)));
		enei=0;
	}
	else{
		enef=-u;
		enei=1;
	}
}

void RNASimulation::contactenergysingle(vector<Residueatoms> decstr, int numseq, double &fene1, double &fene2, Energyterms &lp, vector<vector<vector<double>>> distmaprpp, vector<Phycontact> &contactenerpp, vector<vector<bool>> mapbbflagrpp,bool flagacc)
{
	double totene=0;
	int i,j,k,l,n;
	Point3d tp;
	double tval;
	int aidx;
	vector<int> indn(CONTACT_ENERGY_NUM, 0);
	vector<double> indf(CONTACT_ENERGY_NUM, 0.0);
	int enei;
	double enef;
	double dist;

	if(!flagacc)
	{
		for(l=0;l<CONTACT_ENERGY_NUM;l++){
			lp.fn[CONTACT_ENERGY_ST+l]=0;
			lp.nn[CONTACT_ENERGY_ST+l]=0;
		}

		for(i=0;i<numseq;i++)
		{
			for(j=i+1;j<numseq;j++)
			{
				k=j*numseq+i;

				for(l=0;l<CONTACT_ENERGY_NUM;l++){
					contactenerpp[k].pval[l]=0;
					contactenerpp[k].pcw[l]=0;
				}

				if(CONTACTBINARY[0][i][j] == 1){
					get3gcontact8A(distmaprpp[i][j][CTARGETATOM], CONTACTSCORE[0][i][j], enei, enef);
					contactenerpp[k].pval[0]+=(enef+1);
					contactenerpp[k].pcw[0]+=enei;
				}
				//}else{
					//get3gnoncontact8A(distmaprpp[i][j][CTARGETATOM], 1-CONTACTSCORE[0][i][j], enei, enef);
					//contactenerpp[k].pval[0]+=(enef+1);
					//contactenerpp[k].pcw[0]+=enei;
				//}
				
				if(CONTACTBINARY[1][i][j] == 1){
					get3gcontact20A(distmaprpp[i][j][CTARGETATOM], CONTACTSCORE[1][i][j], enei, enef);
					contactenerpp[k].pval[1]+=(enef+1);
					contactenerpp[k].pcw[1]+=enei;
				}
				//}else{
					//get3gnoncontact20A(distmaprpp[i][j][CTARGETATOM], 1-CONTACTSCORE[1][i][j], enei, enef);
					//contactenerpp[k].pval[1]+=(enef+1);
					//contactenerpp[k].pcw[1]+=enei;
				//}
				for(l=0;l<CONTACT_ENERGY_NUM;l++){
					indf[l]+=contactenerpp[k].pval[l];
					indn[l]+=contactenerpp[k].pcw[l];
				}
			}
		}
		for(l=0;l<CONTACT_ENERGY_NUM;l++){
			lp.fn[CONTACT_ENERGY_ST+l]=indf[l];
			lp.nn[CONTACT_ENERGY_ST+l]=indn[l];
		}
	}
	else{
		vector<double> eold(CONTACT_ENERGY_NUM, 0.0);
		vector<double> enew(CONTACT_ENERGY_NUM, 0.0);
		vector<int> iold(CONTACT_ENERGY_NUM, 0);
		vector<int> inew(CONTACT_ENERGY_NUM, 0);

		for(i=0;i<numseq;i++)
		{
			for(j=i+1;j<numseq;j++)
			{
				k=j*numseq+i;

				if(mapbbflagrpp[i][j])
				{
					for(l=0;l<CONTACT_ENERGY_NUM;l++){
						eold[l]+=contactenerpp[k].pval[l];
						iold[l]+=contactenerpp[k].pcw[l];
					}

					contactenerpp[k].pval[0]=0;
					contactenerpp[k].pcw[0]=0;
					if(CONTACTBINARY[0][i][j] == 1){
						get3gcontact8A(distmaprpp[i][j][CTARGETATOM], CONTACTSCORE[0][i][j], enei, enef);
						contactenerpp[k].pval[0]+=(enef+1);
						contactenerpp[k].pcw[0]+=enei;
					}
					//}else{
						//get3gnoncontact8A(distmaprpp[i][j][CTARGETATOM], 1-CONTACTSCORE[0][i][j], enei, enef);
						//contactenerpp[k].pval[0]+=(enef+1);
						//contactenerpp[k].pcw[0]+=enei;
					//}

					if(CONTACTBINARY[1][i][j] == 1){
						get3gcontact20A(distmaprpp[i][j][CTARGETATOM], CONTACTSCORE[1][i][j], enei, enef);
						contactenerpp[k].pval[1]+=(enef+1);
						contactenerpp[k].pcw[1]+=enei;
					}
					//}else{
						//get3gnoncontact20A(distmaprpp[i][j][CTARGETATOM], 1-CONTACTSCORE[1][i][j], enei, enef);
						//contactenerpp[k].pval[1]+=(enef+1);
						//contactenerpp[k].pcw[1]+=enei;
					//}

					for(l=0;l<CONTACT_ENERGY_NUM;l++){
						enew[l]+=contactenerpp[k].pval[l];
						inew[l]+=contactenerpp[k].pcw[l];
					}
				}
			}
		}
		for(l=0;l<CONTACT_ENERGY_NUM;l++){
			lp.fn[CONTACT_ENERGY_ST+l]+=(enew[l]-eold[l]);
			lp.nn[CONTACT_ENERGY_ST+l]+=(inew[l]-iold[l]);
		}
	}
	
	fene1=lp.fn[CONTACT_ENERGY_ST];
	fene2=lp.fn[CONTACT_ENERGY_ST+1];

	//for(i=0;i<numseq;i++)
	//{
		//for(j=i+1;j<numseq;j++)
		//{
			//k=j*numseq+i;
			//lp.indn[CONTACT_ENERGY_ST][i]+=contactenerpp[k].pcw[0];
			//lp.indn[CONTACT_ENERGY_ST][j]+=contactenerpp[k].pcw[0];
			//lp.indf[CONTACT_ENERGY_ST][i]+=contactenerpp[k].pval[0];
			//lp.indf[CONTACT_ENERGY_ST][j]+=contactenerpp[k].pval[0];
		//}
	//}
}









