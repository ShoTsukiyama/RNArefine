// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


void RNASimulation::fenergypairdihedral(vector<Residueatoms> decstr,int numseq,double &fene, Energyterms &lp,bool flagacc)
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
		for(j=0;j<PairTorParam['A'].size();j++){
			lp.fn[PAIRTORSION_ENERGY_ST+j]=0;
		}

		for(i=0;i<numseq;i++)
		{
			aat=decstr[i].aaa;

			dist1=getdisttor(decstr[i].dih[ClsTorParams[aat][0]], TORCLSPARAM[aat][0][0].mean);
			dist2=getdisttor(decstr[i].dih[ClsTorParams[aat][0]], TORCLSPARAM[aat][1][0].mean);

			if(dist1 < dist2){
				for(j = 0;j<PairTorParam[aat].size();j++){
					binidx1 = getbin(decstr[i].dih[PairTorParam[aat][j][0]], TORSIONPAIRSPARAM[aat][0][j].boundary1);
					binidx2 = getbin(decstr[i].dih[PairTorParam[aat][j][1]], TORSIONPAIRSPARAM[aat][0][j].boundary2);
					binnum2=TORSIONPAIRSPARAM[aat][0][j].binnum2;

					lp.indf[PAIRTORSION_ENERGY_ST+j][i]=TORSIONPAIRSPARAM[aat][0][j].energy[binidx1*binnum2+binidx2];
					lp.fn[PAIRTORSION_ENERGY_ST+j]+=lp.indf[PAIRTORSION_ENERGY_ST+j][i];

					//cout<<to_string(i) + "-" + to_string(j) + "-" + to_string(binidx1) + "-" + to_string(binidx2) + "-" + to_string(lp.indf[PAIRTORSION_ENERGY_ST+j][i])<<endl;
				}
			}else{
				for(j = 0;j<PairTorParam[aat].size();j++){
					binidx1 = getbin(decstr[i].dih[PairTorParam[aat][j][0]], TORSIONPAIRSPARAM[aat][1][j].boundary1);
					binidx2 = getbin(decstr[i].dih[PairTorParam[aat][j][1]], TORSIONPAIRSPARAM[aat][1][j].boundary2);
					binnum2=TORSIONPAIRSPARAM[aat][1][j].binnum2;

					lp.indf[PAIRTORSION_ENERGY_ST+j][i]=TORSIONPAIRSPARAM[aat][0][j].energy[binidx1*binnum2+binidx2];
					lp.fn[PAIRTORSION_ENERGY_ST+j]+=lp.indf[PAIRTORSION_ENERGY_ST+j][i];
				}
			}
		}
	}
	else
	{
		int iold[PairTorParam['A'].size()],inew[PairTorParam['A'].size()];
		double fold[PairTorParam['A'].size()],fnew[PairTorParam['A'].size()];
		for(i=0;i<PairTorParam['A'].size();i++){
			iold[i]=0; inew[i]=0; fold[i]=0; fnew[i]=0;
		}
		bool upflag;

		for(i=0;i<numseq;i++)
		{
			aat=decstr[i].aaa;
			upflag=false;
			for(j=0;j<PairTorParam[aat].size();j++)
			{
				if(lp.flap[TORSIONANGLE_ENERGY_ST+PairTorParam[aat][j][0]][i]){
					upflag=true;
					break;
				}
				if(lp.flap[TORSIONANGLE_ENERGY_ST+PairTorParam[aat][j][1]][i]){
					upflag=true;
					break;
				}
			}

			if(upflag)
			{
				for(j=0;j<PairTorParam[aat].size();j++)
				{
					fold[j]+=lp.indf[PAIRTORSION_ENERGY_ST+j][i];
					lp.indf[PAIRTORSION_ENERGY_ST+j][i]=0;
				}

				dist1=getdisttor(decstr[i].dih[ClsTorParams[aat][0]], TORCLSPARAM[aat][0][0].mean);
				dist2=getdisttor(decstr[i].dih[ClsTorParams[aat][0]], TORCLSPARAM[aat][1][0].mean);

				if(dist1 < dist2){
					for(j = 0;j<PairTorParam[aat].size();j++){
						binidx1 = getbin(decstr[i].dih[PairTorParam[aat][j][0]], TORSIONPAIRSPARAM[aat][0][j].boundary1);
						binidx2 = getbin(decstr[i].dih[PairTorParam[aat][j][1]], TORSIONPAIRSPARAM[aat][0][j].boundary2);
						binnum2=TORSIONPAIRSPARAM[aat][0][j].binnum2;

						lp.indf[PAIRTORSION_ENERGY_ST+j][i]+=TORSIONPAIRSPARAM[aat][0][j].energy[binidx1*binnum2+binidx2];
						//lp.fn[PAIRTORSION_ENERGY_ST+j]+=lp.indf[PAIRTORSION_ENERGY_ST+j][i];
					}
				}else{
					for(j = 0;j<PairTorParam[aat].size();j++){
						binidx1 = getbin(decstr[i].dih[PairTorParam[aat][j][0]], TORSIONPAIRSPARAM[aat][1][j].boundary1);
						binidx2 = getbin(decstr[i].dih[PairTorParam[aat][j][1]], TORSIONPAIRSPARAM[aat][1][j].boundary2);
						binnum2=TORSIONPAIRSPARAM[aat][1][j].binnum2;

						lp.indf[PAIRTORSION_ENERGY_ST+j][i]+=TORSIONPAIRSPARAM[aat][0][j].energy[binidx1*binnum2+binidx2];
						//lp.fn[PAIRTORSION_ENERGY_ST+j]+=lp.indf[PAIRTORSION_ENERGY_ST+j][i];
					}
				}

				for(j = 0;j<PairTorParam[aat].size();j++)
				{
					fnew[j]+=lp.indf[PAIRTORSION_ENERGY_ST+j][i];
				}
			}
		}

		for(j = 0;j<PairTorParam[aat].size();j++)
		{
			lp.fn[PAIRTORSION_ENERGY_ST+j]+=(fnew[j]-fold[j]);
		}
	}

	fene=0;
	for(j=0;j<PairTorParam[aat].size();j++)
	{
		fene+=lp.fn[PAIRTORSION_ENERGY_ST+j];
	}
}

































