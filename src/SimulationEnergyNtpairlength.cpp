// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


void RNASimulation::fenergyntpairlengths(vector<Residueatoms> decstr, int numseq, double &nene, double &fene, Energyterms &lp, vector<vector<vector<double>>> ntpairdistrpp, vector<vector<bool>> mapflagrpp, bool flagacc)
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
	char aat;
	
	if(!flagacc)
	{
		lp.fn[NTPAIRLENG_ENERGY_ST]=0;
		lp.nn[NTPAIRLENG_ENERGY_ST]=0;

		for(i=0;i<numseq;i++)
		{
			aat=decstr[i].aaa;
			
			bdifftmp=0; difftmp=0;
			for(j=0;j<PAIRLENGTHTARGETNUM;j++)
			{
				for(k=j+1;k<PAIRLENGTHTARGETNUM;k++)
				{
					getoutenelength(ntpairdistrpp[i][j][k], PAIRLENGNTPARAM[aat][0][j][k].mean, PAIRLENGNTPARAM[aat][0][j][k].lower, PAIRLENGNTPARAM[aat][0][j][k].upper, PAIRLENGNTPARAM[aat][0][j][k].slower, PAIRLENGNTPARAM[aat][0][j][k].supper, 1, diff, bdiff);
					bdifftmp+=bdiff; difftmp+=diff;

					//if(bdiff == 1) cout<<to_string(i) + "-" + to_string(j) + "-" + to_string(k) + "-" + to_string(ntpairdistrpp[i][j][k]) + "-" + to_string(PAIRLENGNTPARAM[aat][j][k].mean) + "-" + to_string(PAIRLENGNTPARAM[aat][j][k].lower) + "-" + to_string(PAIRLENGNTPARAM[aat][j][k].upper) + "-" + to_string(diff) + "-" + to_string(bdiff)<<endl;
				}
			}

			lp.indf[NTPAIRLENG_ENERGY_ST][i]=difftmp;
			lp.indn[NTPAIRLENG_ENERGY_ST][i]=bdifftmp;

			lp.fn[NTPAIRLENG_ENERGY_ST]+=lp.indf[NTPAIRLENG_ENERGY_ST][i];
			lp.nn[NTPAIRLENG_ENERGY_ST]+=lp.indn[NTPAIRLENG_ENERGY_ST][i];

		}
	}else{
		int iold,inew;
		double fold,fnew;
		iold=0; inew=0; fold=0; fnew=0;

		for(i=0;i<numseq;i++)
		{
			if(!mapflagrpp[i][i]) continue;

			iold+=lp.indn[NTPAIRLENG_ENERGY_ST][i];
			fold+=lp.indf[NTPAIRLENG_ENERGY_ST][i];

			aat=decstr[i].aaa;

			bdifftmp=0; difftmp=0;
			for(j=0;j<PAIRLENGTHTARGETNUM;j++)
			{
				for(k=j+1;k<PAIRLENGTHTARGETNUM;k++)
				{
					getoutenelength(ntpairdistrpp[i][j][k], PAIRLENGNTPARAM[aat][0][j][k].mean, PAIRLENGNTPARAM[aat][0][j][k].lower, PAIRLENGNTPARAM[aat][0][j][k].upper, PAIRLENGNTPARAM[aat][0][j][k].slower, PAIRLENGNTPARAM[aat][0][j][k].supper, 1, diff, bdiff);
					bdifftmp+=bdiff; difftmp+=diff;

					//cout<<to_string(i) + "-" + to_string(j) + "-" + to_string(k) + "-" + to_string(ntpairdistrpp[i][j][k]) + "-" + to_string(PAIRLENGNTPARAM[aat][j][k].mean) + "-" + to_string(PAIRLENGNTPARAM[aat][j][k].lower) + "-" + to_string(PAIRLENGNTPARAM[aat][j][k].upper) + "-" + to_string(diff) + "-" + to_string(bdiff)<<endl;
				}
			}
			
			lp.indf[NTPAIRLENG_ENERGY_ST][i]=difftmp;
			lp.indn[NTPAIRLENG_ENERGY_ST][i]=bdifftmp;

			fnew+=lp.indf[NTPAIRLENG_ENERGY_ST][i];
			inew+=lp.indn[NTPAIRLENG_ENERGY_ST][i];
		}

		lp.nn[NTPAIRLENG_ENERGY_ST]+=(inew-iold);
		lp.fn[NTPAIRLENG_ENERGY_ST]+=(fnew-fold);
	}
	fene=lp.fn[NTPAIRLENG_ENERGY_ST]; 
	nene=lp.nn[NTPAIRLENG_ENERGY_ST];

}





















