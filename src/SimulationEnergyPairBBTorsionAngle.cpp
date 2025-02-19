// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


void RNASimulation::fenergypairbbdihedral(vector<Residueatoms> decstr,int numseq,double &fene, Energyterms &lp,bool flagacc)
{
	int i,j,k;
	int ti,tj;
	double tdist;
	int istart;
	Point3d tp;
	char aat;
	int ss;
	double totalnn;
	int tmpval1;
	double tmpval2;
	double dist1, dist2;
	double minenef, enef, eneftmp;
	int minenei, enei, eneitmp;
	int binidx1, binidx2, binnum2;
	//int estartidx=BLENGTH_ENERGY_NUM+ILENGTH_ENERGY_NUM+INNERANGLE_ENERGY_NUM;

	if(!flagacc)
	{
		for(j=0;j<PairBBTorParam.size();j++){
			lp.fn[PAIRBBTORSION_ENERGY_ST+j]=0;
		}

		for(i=0;i<numseq;i++)
		{
			aat=decstr[i].aaa;

			for(j = 0;j<PairBBTorParam.size();j++){
				binidx1 = getbin(decstr[i].dih[PairBBTorParam[j][0]], TORSIONBBPAIRSPARAM[j].boundary1);
				binidx2 = getbin(decstr[i].dih[PairBBTorParam[j][1]], TORSIONBBPAIRSPARAM[j].boundary2);
				binnum2=TORSIONBBPAIRSPARAM[j].binnum2;

				lp.indf[PAIRBBTORSION_ENERGY_ST+j][i]=TORSIONBBPAIRSPARAM[j].energy[binidx1*binnum2+binidx2];
				lp.fn[PAIRBBTORSION_ENERGY_ST+j]+=lp.indf[PAIRBBTORSION_ENERGY_ST+j][i];
			}
					//cout<<to_string(i) + "-" + to_string(j) + "-" + to_string(binidx1) + "-" + to_string(binidx2) + "-" + to_string(lp.indf[PAIRBBTORSION_ENERGY_ST+j][i])<<endl;
		}
	}
	else
	{
		int iold[PairBBTorParam.size()],inew[PairBBTorParam.size()];
		double fold[PairBBTorParam.size()],fnew[PairBBTorParam.size()];
		for(i=0;i<PairBBTorParam.size();i++){
			iold[i]=0; inew[i]=0; fold[i]=0; fnew[i]=0;
		}
		bool upflag;

		for(i=0;i<numseq;i++)
		{
			aat=decstr[i].aaa;
			upflag=false;
			for(j=0;j<PairBBTorParam.size();j++)
			{
				if(lp.flap[TORSIONANGLE_ENERGY_ST+PairBBTorParam[j][0]][i]){
					upflag=true;
					break;
				}
				if(lp.flap[TORSIONANGLE_ENERGY_ST+PairBBTorParam[j][1]][i]){
					upflag=true;
					break;
				}
			}

			if(upflag)
			{
				for(j=0;j<PairBBTorParam.size();j++)
				{
					fold[j]+=lp.indf[PAIRBBTORSION_ENERGY_ST+j][i];
					lp.indf[PAIRBBTORSION_ENERGY_ST+j][i]=0;
				}

				
				for(j = 0;j<PairBBTorParam.size();j++){
					binidx1 = getbin(decstr[i].dih[PairBBTorParam[j][0]], TORSIONBBPAIRSPARAM[j].boundary1);
					binidx2 = getbin(decstr[i].dih[PairBBTorParam[j][1]], TORSIONBBPAIRSPARAM[j].boundary2);
					binnum2=TORSIONBBPAIRSPARAM[j].binnum2;

					lp.indf[PAIRBBTORSION_ENERGY_ST+j][i]+=TORSIONBBPAIRSPARAM[j].energy[binidx1*binnum2+binidx2];
					//lp.fn[PAIRBBTORSION_ENERGY_ST+j]+=lp.indf[PAIRBBTORSION_ENERGY_ST+j][i];
				}
				

				for(j = 0;j<PairBBTorParam.size();j++)
				{
					fnew[j]+=lp.indf[PAIRBBTORSION_ENERGY_ST+j][i];
				}
			}
		}

		for(j = 0;j<PairBBTorParam.size();j++)
		{
			lp.fn[PAIRBBTORSION_ENERGY_ST+j]+=(fnew[j]-fold[j]);
		}
	}

	fene=0;
	for(j=0;j<PairBBTorParam.size();j++)
	{
		fene+=lp.fn[PAIRBBTORSION_ENERGY_ST+j];
	}
}

































