// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

/*
void RNASimulation::fenergydihedral(vector<Residueatoms> decstr,int numseq,double &nene,double &fene, Energyterms &lp,bool flagacc)
{
	int i,j;
	int ti,tj;
	double tdist;
	int istart;
	Point3d tp;
	char aat;
	int ss;
	double totalnn;
	int tmpval1;
	double tmpval2;
	//int estartidx=BLENGTH_ENERGY_NUM+ILENGTH_ENERGY_NUM+INNERANGLE_ENERGY_NUM;

	if(!flagacc)
	{
		for(j=0;j<TORSIONANGLETARGET.size();j++){
			lp.fn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]]=0;
			lp.nn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]]=0;

			aat=decstr[0].aaa;
			if(TORSIONANGLETARGET[j] == 10){
				getoutenetorsion(decstr[0].dih[TORSIONANGLETARGET[j]], OPTORSIONPARAM[aat][0].mode, OPTORSIONPARAM[aat][0].boundary, OPTORSIONPARAM[aat][0].energy, OPTORSIONPARAM[aat][0].binary,0.01, lp.indf[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0], lp.indn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0]);
				//getoutlier(decstr[0].dih[TORSIONANGLETARGET[j]], OPTORSIONPARAM[aat][0].boundary, OPTORSIONPARAM[aat][0].energy, OPTORSIONPARAM[aat][0].binary, lp.indf[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0], lp.indn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0]);
				//getoutenetorsion(decstr[0].dih[TORSIONANGLETARGET[j]], OPTORSIONPARAM[aat][0].mode, OPTORSIONPARAM[aat][0].boundary, OPTORSIONPARAM[aat][0].energy, OPTORSIONPARAM[aat][0].binary, lp.indf[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0], lp.indn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0]);
			}else if(TORSIONANGLETARGET[j] == 11){
				getoutenetorsion(decstr[0].dih[TORSIONANGLETARGET[j]], OPTORSIONPARAM[aat][1].mode, OPTORSIONPARAM[aat][1].boundary, OPTORSIONPARAM[aat][1].energy, OPTORSIONPARAM[aat][1].binary,0.01, lp.indf[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0], lp.indn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0]);
				//getoutlier(decstr[0].dih[TORSIONANGLETARGET[j]], OPTORSIONPARAM[aat][1].boundary, OPTORSIONPARAM[aat][1].energy, OPTORSIONPARAM[aat][1].binary, lp.indf[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0], lp.indn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0]);
				//getoutenetorsion(decstr[0].dih[TORSIONANGLETARGET[j]], OPTORSIONPARAM[aat][1].mode, OPTORSIONPARAM[aat][1].boundary, OPTORSIONPARAM[aat][1].energy, OPTORSIONPARAM[aat][1].binary, lp.indf[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0], lp.indn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0]);
			}else{
				getoutenetorsion(decstr[0].dih[TORSIONANGLETARGET[j]], TORSIONPARAM[aat][TORSIONANGLETARGET[j]].mode, TORSIONPARAM[aat][TORSIONANGLETARGET[j]].boundary, TORSIONPARAM[aat][TORSIONANGLETARGET[j]].energy, TORSIONPARAM[aat][TORSIONANGLETARGET[j]].binary,0.01, lp.indf[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0], lp.indn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0]);
				//getoutlier(decstr[0].dih[TORSIONANGLETARGET[j]], TORSIONPARAM[aat][TORSIONANGLETARGET[j]].boundary, TORSIONPARAM[aat][TORSIONANGLETARGET[j]].energy, TORSIONPARAM[aat][TORSIONANGLETARGET[j]].binary, lp.indf[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0], lp.indn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0]);
				//getoutenetorsion(decstr[0].dih[TORSIONANGLETARGET[j]], TORSIONPARAM[aat][TORSIONANGLETARGET[j]].mode, TORSIONPARAM[aat][TORSIONANGLETARGET[j]].boundary, TORSIONPARAM[aat][TORSIONANGLETARGET[j]].energy, TORSIONPARAM[aat][TORSIONANGLETARGET[j]].binary, lp.indf[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0], lp.indn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0]);
			}
			lp.fn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]]+=lp.indf[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0];
			lp.nn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]]+=lp.indn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0];
			
			for(i=1;i<numseq;i++)
			{
				aat=decstr[i].aaa;
				getoutenetorsion(decstr[i].dih[TORSIONANGLETARGET[j]], TORSIONPARAM[aat][TORSIONANGLETARGET[j]].mode, TORSIONPARAM[aat][TORSIONANGLETARGET[j]].boundary, TORSIONPARAM[aat][TORSIONANGLETARGET[j]].energy, TORSIONPARAM[aat][TORSIONANGLETARGET[j]].binary,0.01, lp.indf[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][i], lp.indn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][i]);
				
				//if(lp.indn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][i] != 0){
					//cout<<to_string(1) + "-" + to_string(i) + "-" + to_string(j) + "-" + to_string(decstr[i].dih[TORSIONANGLETARGET[j]]) + "-" + to_string(lp.indf[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][i]) + "-" + to_string(lp.indn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][i])<<endl;
				//}
				//else{
					//cout<<to_string(0) + "-" + to_string(i) + "-" + to_string(j) + "-" + to_string(decstr[i].dih[TORSIONANGLETARGET[j]]) + "-" + to_string(lp.indf[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][i]) + "-" + to_string(lp.indn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][i])<<endl;
				//}
				
				lp.fn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]]+=lp.indf[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][i];
				lp.nn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]]+=lp.indn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][i];
			}
		}
	}
	else
	{
		int iold[TORSIONANGLETARGET.size()],inew[TORSIONANGLETARGET.size()];
		double fold[TORSIONANGLETARGET.size()],fnew[TORSIONANGLETARGET.size()];
		for(i=0;i<TORSIONANGLETARGET.size();i++){
			iold[i]=0; inew[i]=0; fold[i]=0; fnew[i]=0;
		}
		
		for(j=0;j<TORSIONANGLETARGET.size();j++){

			if(lp.flap[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0]){
				iold[j]+=lp.indn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0];
				fold[j]+=lp.indf[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0];

				aat=decstr[0].aaa;
				if(TORSIONANGLETARGET[j] == 10){
					getoutenetorsion(decstr[0].dih[TORSIONANGLETARGET[j]], OPTORSIONPARAM[aat][0].mode, OPTORSIONPARAM[aat][0].boundary, OPTORSIONPARAM[aat][0].energy, OPTORSIONPARAM[aat][0].binary,0.01, lp.indf[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0], lp.indn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0]);
					//getoutlier(decstr[0].dih[TORSIONANGLETARGET[j]], OPTORSIONPARAM[aat][0].boundary, OPTORSIONPARAM[aat][0].energy, OPTORSIONPARAM[aat][0].binary, lp.indf[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0], lp.indn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0]);
					//getoutenetorsion(decstr[0].dih[TORSIONANGLETARGET[j]], OPTORSIONPARAM[aat][0].mode, OPTORSIONPARAM[aat][0].boundary, OPTORSIONPARAM[aat][0].energy, OPTORSIONPARAM[aat][0].binary, lp.indf[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0], lp.indn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0]);
				}else if(TORSIONANGLETARGET[j] == 11){
					getoutenetorsion(decstr[0].dih[TORSIONANGLETARGET[j]], OPTORSIONPARAM[aat][1].mode, OPTORSIONPARAM[aat][1].boundary, OPTORSIONPARAM[aat][1].energy, OPTORSIONPARAM[aat][1].binary,0.01, lp.indf[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0], lp.indn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0]);
					//getoutlier(decstr[0].dih[TORSIONANGLETARGET[j]], OPTORSIONPARAM[aat][1].boundary, OPTORSIONPARAM[aat][1].energy, OPTORSIONPARAM[aat][1].binary, lp.indf[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0], lp.indn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0]);
					//getoutenetorsion(decstr[0].dih[TORSIONANGLETARGET[j]], OPTORSIONPARAM[aat][1].mode, OPTORSIONPARAM[aat][1].boundary, OPTORSIONPARAM[aat][1].energy, OPTORSIONPARAM[aat][1].binary, lp.indf[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0], lp.indn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0]);
				}else{
					getoutenetorsion(decstr[0].dih[TORSIONANGLETARGET[j]], TORSIONPARAM[aat][TORSIONANGLETARGET[j]].mode, TORSIONPARAM[aat][TORSIONANGLETARGET[j]].boundary, TORSIONPARAM[aat][TORSIONANGLETARGET[j]].energy, TORSIONPARAM[aat][TORSIONANGLETARGET[j]].binary,0.01, lp.indf[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0], lp.indn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0]);
					//getoutlier(decstr[0].dih[TORSIONANGLETARGET[j]], TORSIONPARAM[aat][TORSIONANGLETARGET[j]].boundary, TORSIONPARAM[aat][TORSIONANGLETARGET[j]].energy, TORSIONPARAM[aat][TORSIONANGLETARGET[j]].binary, lp.indf[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0], lp.indn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0]);
					//getoutenetorsion(decstr[0].dih[TORSIONANGLETARGET[j]], TORSIONPARAM[aat][TORSIONANGLETARGET[j]].mode, TORSIONPARAM[aat][TORSIONANGLETARGET[j]].boundary, TORSIONPARAM[aat][TORSIONANGLETARGET[j]].energy, TORSIONPARAM[aat][TORSIONANGLETARGET[j]].binary, lp.indf[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0], lp.indn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0]);
				}
				inew[j]+=lp.indn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0];
				fnew[j]+=lp.indf[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0];
				lp.flap[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][0]=false;
			}

			for(i=1;i<numseq;i++) if(lp.flap[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][i])
			{
				iold[j]+=lp.indn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][i];
				fold[j]+=lp.indf[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][i];

				aat=decstr[i].aaa;
				getoutenetorsion(decstr[i].dih[TORSIONANGLETARGET[j]], TORSIONPARAM[aat][TORSIONANGLETARGET[j]].mode, TORSIONPARAM[aat][TORSIONANGLETARGET[j]].boundary, TORSIONPARAM[aat][TORSIONANGLETARGET[j]].energy, TORSIONPARAM[aat][TORSIONANGLETARGET[j]].binary,0.01, lp.indf[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][i], lp.indn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][i]);

				inew[j]+=lp.indn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][i];
				fnew[j]+=lp.indf[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][i];
				lp.flap[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][i]=false;
			}
			
			lp.nn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]]+=(inew[j]-iold[j]);
			lp.fn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]]+=(fnew[j]-fold[j]);
		}
	}

	fene=0; nene=0;
	for(j=0;j<TORSIONANGLETARGET.size();j++){
		fene+=lp.fn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]];
	}
	for(j=0;j<TORSIONANGLETARGET.size();j++){
		nene+=lp.nn[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]];
	}
}
*/


double RNASimulation::getdisttor(double val, double mean){
	double tmp1, tmp2;

	tmp1=abs(val-mean);
	tmp2=360-tmp1;

	if(tmp1 > tmp2) return tmp2;
	else return tmp1;
}



void RNASimulation::fenergydihedral(vector<Residueatoms> decstr,int numseq,double &nene,double &fene, Energyterms &lp,bool flagacc)
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
	vector<double> minenevecftmp;
	vector<int>minenevecitmp;
	vector<double> minenevecf(6);
	vector<int> mineneveci(6);
	//int estartidx=BLENGTH_ENERGY_NUM+ILENGTH_ENERGY_NUM+INNERANGLE_ENERGY_NUM;

	if(!flagacc)
	{
		for(i=0;i<numseq;i++)
		{
			aat=decstr[i].aaa;

			dist1=getdisttor(decstr[i].dih[ClsTorParams[aat][0]], TORCLSPARAM[aat][0][0].mean);
			dist2=getdisttor(decstr[i].dih[ClsTorParams[aat][0]], TORCLSPARAM[aat][1][0].mean);

			if(dist1 < dist2){
				for(j=0;j<6;j++){
					getoutenetorangle(decstr[i].dih[ClsTorParams[aat][j]], TORCLSPARAM[aat][0][j].mean, TORCLSPARAM[aat][0][j].lower, TORCLSPARAM[aat][0][j].upper, 1, eneftmp, eneitmp);
					lp.indf[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]][i]+=eneftmp;
					lp.indn[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]][i]+=eneitmp;
					lp.fn[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]]+=lp.indf[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]][i];
					lp.nn[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]]+=lp.indn[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]][i];
				}
			}else{
				for(j=0;j<6;j++){
					getoutenetorangle(decstr[i].dih[ClsTorParams[aat][j]], TORCLSPARAM[aat][1][j].mean, TORCLSPARAM[aat][1][j].lower, TORCLSPARAM[aat][1][j].upper, 1, eneftmp, eneitmp);
					lp.indf[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]][i]+=eneftmp;
					lp.indn[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]][i]+=eneitmp;
					lp.fn[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]]+=lp.indf[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]][i];
					lp.nn[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]]+=lp.indn[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]][i];
				}
			}
		}
	}
	else
	{
		int iold[8],inew[8];
		double fold[8],fnew[8];
		for(i=0;i<8;i++){
			iold[i]=0; inew[i]=0; fold[i]=0; fnew[i]=0;
		}
		bool upflag;

		for(i=0;i<numseq;i++)
		{
			aat=decstr[i].aaa;
			upflag=false;
			for(j=0;j<6;j++)
			{
				if(lp.flap[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]][i]){
					upflag=true;
					break;
				}
			}

			if(upflag)
			{
				for(j=0;j<6;j++)
				{
					iold[j]+=lp.indn[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]][i];
					fold[j]+=lp.indf[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]][i];
					lp.indf[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]][i]=0;
					lp.indn[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]][i]=0;
				}

				dist1=getdisttor(decstr[i].dih[ClsTorParams[aat][0]], TORCLSPARAM[aat][0][0].mean);
				dist2=getdisttor(decstr[i].dih[ClsTorParams[aat][0]], TORCLSPARAM[aat][1][0].mean);

				if(dist1 < dist2){
				for(j=0;j<6;j++){
					getoutenetorangle(decstr[i].dih[ClsTorParams[aat][j]], TORCLSPARAM[aat][0][j].mean, TORCLSPARAM[aat][0][j].lower, TORCLSPARAM[aat][0][j].upper, 1, eneftmp, eneitmp);
					lp.indf[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]][i]+=eneftmp;
					lp.indn[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]][i]+=eneitmp;
					//lp.fn[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]]+=lp.indf[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]][i];
					//lp.nn[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]]+=lp.indn[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]][i];
					}
				}else{
					for(j=0;j<6;j++){
						getoutenetorangle(decstr[i].dih[ClsTorParams[aat][j]], TORCLSPARAM[aat][1][j].mean, TORCLSPARAM[aat][1][j].lower, TORCLSPARAM[aat][1][j].upper, 1, eneftmp, eneitmp);
						lp.indf[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]][i]+=eneftmp;
						lp.indn[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]][i]+=eneitmp;
						//lp.fn[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]]+=lp.indf[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]][i];
						//lp.nn[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]]+=lp.indn[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]][i];
					}
				}

				for(j=0;j<6;j++)
				{
					//lp.indn[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]][i]=mineneveci[j];
					//lp.indf[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]][i]=minenevecf[j];

					inew[j]+=lp.indn[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]][i];
					fnew[j]+=lp.indf[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]][i];

					//lp.fn[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]]+=lp.indf[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]][i];
					//lp.nn[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]]+=lp.indn[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]][i];
				}

				//getoutenetorangle(decstr[i].dih[0], 160, 300, 1, eneftmp, eneitmp);
				//lp.indf[TORSIONANGLE_ENERGY_ST][i]=eneftmp;
				//lp.indn[TORSIONANGLE_ENERGY_ST][i]=eneitmp;

				//fnew[6]+=lp.indf[TORSIONANGLE_ENERGY_ST][i];
				//inew[6]+=lp.indn[TORSIONANGLE_ENERGY_ST][i];

				//getoutenetorangle(decstr[i].dih[3], 60, 280, 1, eneftmp, eneitmp);
				//lp.indf[TORSIONANGLE_ENERGY_ST+3][i]=eneftmp;
				//lp.indn[TORSIONANGLE_ENERGY_ST+3][i]=eneitmp;

				//fnew[7]+=lp.indf[TORSIONANGLE_ENERGY_ST+3][i];
				//inew[7]+=lp.indn[TORSIONANGLE_ENERGY_ST+3][i];
			}
		}

		for(j=0;j<6;j++)
		{
			lp.nn[TORSIONANGLE_ENERGY_ST+ClsTorParams['A'][j]]+=(inew[j]-iold[j]);
			lp.fn[TORSIONANGLE_ENERGY_ST+ClsTorParams['A'][j]]+=(fnew[j]-fold[j]);
		}
		//lp.nn[TORSIONANGLE_ENERGY_ST]+=(inew[6]-iold[6]);
		//lp.fn[TORSIONANGLE_ENERGY_ST]+=(fnew[6]-fold[6]);
		//lp.nn[TORSIONANGLE_ENERGY_ST+3]+=(inew[7]-iold[7]);
		//lp.fn[TORSIONANGLE_ENERGY_ST+3]+=(fnew[7]-fold[7]);
	}

	fene=0; nene=0;
	for(j=0;j<6;j++)
	{
		fene+=lp.fn[TORSIONANGLE_ENERGY_ST+ClsTorParams['A'][j]];
	}
	//fene+=lp.fn[TORSIONANGLE_ENERGY_ST];
	//fene+=lp.fn[TORSIONANGLE_ENERGY_ST+3];
	for(j=0;j<6;j++)
	{
		nene+=lp.nn[TORSIONANGLE_ENERGY_ST+ClsTorParams['A'][j]];
	}
	//nene+=lp.nn[TORSIONANGLE_ENERGY_ST];
	//nene+=lp.nn[TORSIONANGLE_ENERGY_ST+3];
}

































