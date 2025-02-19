// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


void RNASimulation::fenergymultidihedral(vector<Residueatoms> decstr,int numseq,double &nene,double &fene, Energyterms &lp,bool flagacc)
{
	int i,j;
	char aat;
	int toridx1, toridx2;
	double dih1, dih2;
	int binidx1, binidx2;
	int binnum1, binnum2;
	double indftmp;
	int indntmp;
	//int estartidx=BLENGTH_ENERGY_NUM+ILENGTH_ENERGY_NUM+INNERANGLE_ENERGY_NUM;

	if(!flagacc)
	{
		for(j=0;j<MULTITORSION_ENERGY_NUM;j++)
		{
			lp.fn[MULTITORSION_ENERGY_ST+j]=0;
			lp.nn[MULTITORSION_ENERGY_ST+j]=0;

			for(i=1;i<numseq;i++)
			{
				aat=decstr[i].aaa;
				binnum2=MULTITORSIONPARAM[aat][j].binnum2;
				toridx1 = CombTorParamIdxList[aat][j][0];
				toridx2 = CombTorParamIdxList[aat][j][1];

				dih1 = decstr[i].dih[toridx1];
				dih2 = decstr[i].dih[toridx2];

				binidx1=getbin(dih1, MULTITORSIONPARAM[aat][j].boundary1);
				binidx2=getbin(dih2, MULTITORSIONPARAM[aat][j].boundary2);

				//cout<<to_string(binidx1) + "-" + to_string(binidx2)<<endl;
				//cout<<binnum2<<endl;

				lp.indn[MULTITORSION_ENERGY_ST+j][i]=MULTITORSIONPARAM[aat][j].binary[binidx1*binnum2+binidx2];
				lp.indf[MULTITORSION_ENERGY_ST+j][i]=MULTITORSIONPARAM[aat][j].energy[binidx1*binnum2+binidx2];

				lp.fn[MULTITORSION_ENERGY_ST+j]+=lp.indf[MULTITORSION_ENERGY_ST+j][i];
				lp.nn[MULTITORSION_ENERGY_ST+j]+=lp.indn[MULTITORSION_ENERGY_ST+j][i];
			}
		}
	}
	else
	{
		int iold[MULTITORSION_ENERGY_NUM],inew[MULTITORSION_ENERGY_NUM];
		double fold[MULTITORSION_ENERGY_NUM],fnew[MULTITORSION_ENERGY_NUM];
		for(i=0;i<MULTITORSION_ENERGY_NUM;i++){
			iold[i]=0; inew[i]=0; fold[i]=0; fnew[i]=0;
		}

		for(j=0;j<MULTITORSION_ENERGY_NUM;j++)
		{
			for(i=1;i<numseq;i++)
			{
				aat=decstr[i].aaa;
				binnum2=MULTITORSIONPARAM[aat][j].binnum2;
				toridx1 = CombTorParamIdxList[aat][j][0];
				toridx2 = CombTorParamIdxList[aat][j][1];

				if((!lp.flap[TORSIONANGLE_ENERGY_ST+toridx1][i]) & (!lp.flap[TORSIONANGLE_ENERGY_ST+toridx2][i])) continue;

				iold[j]+=lp.indn[MULTITORSION_ENERGY_ST+j][i];
				fold[j]+=lp.indf[MULTITORSION_ENERGY_ST+j][i];

				lp.indn[MULTITORSION_ENERGY_ST+j][i]=0;
				lp.indf[MULTITORSION_ENERGY_ST+j][i]=0;

				dih1 = decstr[i].dih[toridx1];
				dih2 = decstr[i].dih[toridx2];

				binidx1=getbin(dih1, MULTITORSIONPARAM[aat][j].boundary1);
				binidx2=getbin(dih2, MULTITORSIONPARAM[aat][j].boundary2);

				lp.indn[MULTITORSION_ENERGY_ST+j][i]+=MULTITORSIONPARAM[aat][j].binary[binidx1*binnum2+binidx2];
				lp.indf[MULTITORSION_ENERGY_ST+j][i]+=MULTITORSIONPARAM[aat][j].energy[binidx1*binnum2+binidx2];

				inew[j]+=lp.indn[MULTITORSION_ENERGY_ST+j][i];
				fnew[j]+=lp.indf[MULTITORSION_ENERGY_ST+j][i];
				
			}

			lp.nn[MULTITORSION_ENERGY_ST+j]+=(inew[j]-iold[j]);
			lp.fn[MULTITORSION_ENERGY_ST+j]+=(fnew[j]-fold[j]);
		}
	}

	fene=0; nene=0;
	for(j=0;j<MULTITORSION_ENERGY_NUM;j++){
		fene+=lp.fn[MULTITORSION_ENERGY_ST+j];
	}
	for(j=0;j<MULTITORSION_ENERGY_NUM;j++){
		nene+=lp.nn[MULTITORSION_ENERGY_ST+j];
	}
}





















