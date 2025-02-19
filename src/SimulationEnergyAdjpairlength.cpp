// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


void RNASimulation::fenergyadjpairlengths(vector<Residueatoms> decstr, int numseq, double &nene, double &fene, Energyterms &lp, vector<vector<vector<double>>> adjpairdistrpp, vector<vector<bool>> mapflagrpp, bool flagacc)
{
	int i,j,k,l,n;
	int ti,tj;
	double tdist;
	int istart;
	Point3d tp;
	int ss;
	double totalnn;
	int tmpval1;
	double tmpval2;
	int bdiff, bdifftmp;
	double diff, difftmp;
	bool flag;
	double tmptmp;
	char aat1, aat2;
	
	if(!flagacc)
	{
		lp.fn[ADJPAIRLENG_ENERGY_ST]=0;
		lp.nn[ADJPAIRLENG_ENERGY_ST]=0;

		for(i=1;i<numseq;i++)
		{
			aat1=decstr[i-1].aaa;
			aat2=decstr[i].aaa;
			
			bdifftmp=0; difftmp=0;
			for(j=0;j<PAIRLENGTHTARGETNUM;j++)
			{
				for(k=0;k<PAIRLENGTHTARGETNUM;k++)
				{
					getoutenelength(adjpairdistrpp[i][j][k], PAIRLENGADJPARAM[aat1][aat2][0][j][k].mean, PAIRLENGADJPARAM[aat1][aat2][0][j][k].lower, PAIRLENGADJPARAM[aat1][aat2][0][j][k].upper, PAIRLENGADJPARAM[aat1][aat2][0][j][k].slower, PAIRLENGADJPARAM[aat1][aat2][0][j][k].supper, 1, diff, bdiff);
					bdifftmp+=bdiff; difftmp+=diff;

					//if(bdiff == 1) 
					//cout<<to_string(i) + "-" + to_string(j) + "-" + to_string(k) + "-" + to_string(adjpairdistrpp[i][j][k]) + "-" + to_string(PAIRLENGADJPARAM[aat1][aat2][j][k].mean) + "-" + to_string(PAIRLENGADJPARAM[aat1][aat2][j][k].lower) + "-" + to_string(PAIRLENGADJPARAM[aat1][aat2][j][k].upper) + "-" + to_string(diff) + "-" + to_string(bdiff)<<endl;
				}
			}

			lp.indf[ADJPAIRLENG_ENERGY_ST][i]=difftmp;
			lp.indn[ADJPAIRLENG_ENERGY_ST][i]=bdifftmp;

			lp.fn[ADJPAIRLENG_ENERGY_ST]+=lp.indf[ADJPAIRLENG_ENERGY_ST][i];
			lp.nn[ADJPAIRLENG_ENERGY_ST]+=lp.indn[ADJPAIRLENG_ENERGY_ST][i];

		}
	}else{
		int iold,inew;
		double fold,fnew;
		iold=0; inew=0; fold=0; fnew=0;

		for(i=1;i<numseq;i++)
		{
			if(!mapflagrpp[i-1][i]) continue;

			iold+=lp.indn[ADJPAIRLENG_ENERGY_ST][i];
			fold+=lp.indf[ADJPAIRLENG_ENERGY_ST][i];

			aat1=decstr[i-1].aaa;
			aat2=decstr[i].aaa;

			bdifftmp=0; difftmp=0;
			for(j=0;j<PAIRLENGTHTARGETNUM;j++)
			{
				for(k=0;k<PAIRLENGTHTARGETNUM;k++)
				{
					getoutenelength(adjpairdistrpp[i][j][k], PAIRLENGADJPARAM[aat1][aat2][0][j][k].mean, PAIRLENGADJPARAM[aat1][aat2][0][j][k].lower, PAIRLENGADJPARAM[aat1][aat2][0][j][k].upper, PAIRLENGADJPARAM[aat1][aat2][0][j][k].slower, PAIRLENGADJPARAM[aat1][aat2][0][j][k].supper, 1, diff, bdiff);
					bdifftmp+=bdiff; difftmp+=diff;
				}
			}
			
			lp.indf[ADJPAIRLENG_ENERGY_ST][i]=difftmp;
			lp.indn[ADJPAIRLENG_ENERGY_ST][i]=bdifftmp;

			inew+=lp.indn[ADJPAIRLENG_ENERGY_ST][i];
			fnew+=lp.indf[ADJPAIRLENG_ENERGY_ST][i];
		}

		lp.nn[ADJPAIRLENG_ENERGY_ST]+=(inew-iold);
		lp.fn[ADJPAIRLENG_ENERGY_ST]+=(fnew-fold);
	}
	fene=lp.fn[ADJPAIRLENG_ENERGY_ST]; 
	nene=lp.nn[ADJPAIRLENG_ENERGY_ST];

}





















