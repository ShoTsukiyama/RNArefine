// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

void RNASimulation::getoutenetorpucker(double val, double mean, double lower, double upper, double k, double &diff, int &bdiff)
{
	double tmp;

	if(lower <= val && val <= upper){
		bdiff=0;
		diff=0;
	}else{
		bdiff=1;
		tmp=abs(val-mean);
		diff=k*min(tmp, 360-tmp);
	}
}

double distanceFromPointToLine(const Point3s& A, const Point3s& B, const Point3s& P) {
    // ベクトル AB
    double ABx = B.x - A.x;
    double ABy = B.y - A.y;
    double ABz = B.z - A.z;

    // ベクトル AP
    double APx = P.x - A.x;
    double APy = P.y - A.y;
    double APz = P.z - A.z;

    // ABの大きさ
    double AB_length = std::sqrt(ABx * ABx + ABy * ABy + ABz * ABz);

    // 外積 AB × AP
    double crossX = ABy * APz - ABz * APy;
    double crossY = ABz * APx - ABx * APz;
    double crossZ = ABx * APy - ABy * APx;

    // 外積の大きさ
    double cross_length = std::sqrt(crossX * crossX + crossY * crossY + crossZ * crossZ);

    // 距離を計算
    return cross_length / AB_length;
}


void RNASimulation::fenergypucker(vector<Residueatoms> decstr,int numseq,double &nene,double &fene, Energyterms &lp, vector<vector<bool>> mapflagrpp, bool flagacc)
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

	int atomn, atomc, atomp;
	Point3s tmpn,tmpc,tmpp;
	double dist, deltator;
	double enefdist1, enefdist2, eneftor1, eneftor2, enef1, enef2;
	int enendist1, enendist2, enentor1, enentor2, enen1, enen2;
	//int estartidx=BLENGTH_ENERGY_NUM+ILENGTH_ENERGY_NUM+INNERANGLE_ENERGY_NUM;

	if(!flagacc)
	{
		lp.fn[PUCKER_ENERGY_ST]=0;
		lp.nn[PUCKER_ENERGY_ST]=0;

		for(i=0;i<numseq-1;i++)
		{
			aat=decstr[i].aaa;
			atomn=PHOPHGLYCOSIDEATOMS[aat][1];
			atomc=PHOPHGLYCOSIDEATOMS[aat][2];
			atomp=PHOPHGLYCOSIDEATOMS[aat][0];

			tmpn=decstr[i].ptv[atomn];
			tmpc=decstr[i].ptv[atomc];
			tmpp=decstr[i+1].ptv[atomp];

			dist=distanceFromPointToLine(tmpn, tmpc, tmpp);
			deltator=decstr[i].dih[5];

			getoutenelength(dist, PUCKERPARAM[0][0].mean, PUCKERPARAM[0][0].lower, PUCKERPARAM[0][0].upper, 1, enefdist1, enendist1);
			getoutenetorpucker(deltator, PUCKERPARAM[0][1].mean, PUCKERPARAM[0][1].lower, PUCKERPARAM[0][1].upper, 0.05, eneftor1, enentor1);
			enen1=enendist1+enentor1;
			enef1=enefdist1+eneftor1;

			getoutenelength(dist, PUCKERPARAM[1][0].mean, PUCKERPARAM[1][0].lower, PUCKERPARAM[1][0].upper, 1, enefdist2, enendist2);
			getoutenetorpucker(deltator, PUCKERPARAM[1][1].mean, PUCKERPARAM[1][1].lower, PUCKERPARAM[1][1].upper, 0.05, eneftor2, enentor2);
			enen2=enendist2+enentor2;
			enef2=enefdist2+eneftor2;

			if(enen1 == enen2){
				if(enef1 < enef2){
					lp.indf[PUCKER_ENERGY_ST][i]=enef1;
					lp.indn[PUCKER_ENERGY_ST][i]=enen1;
				}else{
					lp.indf[PUCKER_ENERGY_ST][i]=enef2;
					lp.indn[PUCKER_ENERGY_ST][i]=enen2;
				}
			}else if(enen1 < enen2){
				lp.indf[PUCKER_ENERGY_ST][i]=enef1;
				lp.indn[PUCKER_ENERGY_ST][i]=enen1;
			}else{
				lp.indf[PUCKER_ENERGY_ST][i]=enef2;
				lp.indn[PUCKER_ENERGY_ST][i]=enen2;
			}

			lp.fn[PUCKER_ENERGY_ST]+=lp.indf[PUCKER_ENERGY_ST][i];
			lp.nn[PUCKER_ENERGY_ST]+=lp.indn[PUCKER_ENERGY_ST][i];
		}
	}
	else
	{
		int iold[PUCKER_ENERGY_NUM],inew[PUCKER_ENERGY_NUM];
		double fold[PUCKER_ENERGY_NUM],fnew[PUCKER_ENERGY_NUM];
		for(i=0;i<PUCKER_ENERGY_NUM;i++){
			iold[i]=0; inew[i]=0; fold[i]=0; fnew[i]=0;
		}

		for(i=0;i<numseq-1;i++)
		{
			aat=decstr[i].aaa;
			
			if(lp.flap[TORSIONANGLE_ENERGY_ST+PHOPHGLYCOSIDEATOMS[aat][3]][i] | mapflagrpp[i][i+1]){

				iold[0]+=lp.indn[PUCKER_ENERGY_ST][i];
				fold[0]+=lp.indf[PUCKER_ENERGY_ST][i];

				atomn=PHOPHGLYCOSIDEATOMS[aat][1];
				atomc=PHOPHGLYCOSIDEATOMS[aat][2];
				atomp=PHOPHGLYCOSIDEATOMS[aat][0];

				tmpn=decstr[i].ptv[atomn];
				tmpc=decstr[i].ptv[atomc];
				tmpp=decstr[i+1].ptv[atomp];

				dist=distanceFromPointToLine(tmpn, tmpc, tmpp);
				deltator=decstr[i].dih[PHOPHGLYCOSIDEATOMS[aat][3]];

				getoutenelength(dist, PUCKERPARAM[0][0].mean, PUCKERPARAM[0][0].lower, PUCKERPARAM[0][0].upper, 1, enefdist1, enendist1);
				getoutenetorpucker(deltator, PUCKERPARAM[0][1].mean, PUCKERPARAM[0][1].lower, PUCKERPARAM[0][1].upper, 0.05, eneftor1, enentor1);
				enen1=enendist1+enentor1;
				enef1=enefdist1+eneftor1;

				getoutenelength(dist, PUCKERPARAM[1][0].mean, PUCKERPARAM[1][0].lower, PUCKERPARAM[1][0].upper, 1, enefdist2, enendist2);
				getoutenetorpucker(deltator, PUCKERPARAM[1][1].mean, PUCKERPARAM[1][1].lower, PUCKERPARAM[1][1].upper, 0.05, eneftor2, enentor2);
				enen2=enendist2+enentor2;
				enef2=enefdist2+eneftor2;

				if(enen1 == enen2){
					if(enef1 < enef2){
						lp.indf[PUCKER_ENERGY_ST][i]=enef1;
						lp.indn[PUCKER_ENERGY_ST][i]=enen1;
					}else{
						lp.indf[PUCKER_ENERGY_ST][i]=enef2;
						lp.indn[PUCKER_ENERGY_ST][i]=enen2;
					}
			    }else if(enen1 < enen2){
					lp.indf[PUCKER_ENERGY_ST][i]=enef1;
					lp.indn[PUCKER_ENERGY_ST][i]=enen1;
				}else{
					lp.indf[PUCKER_ENERGY_ST][i]=enef2;
					lp.indn[PUCKER_ENERGY_ST][i]=enen2;
				}

				inew[0]+=lp.indn[PUCKER_ENERGY_ST][i];
				fnew[0]+=lp.indf[PUCKER_ENERGY_ST][i];
			}
		}

		lp.fn[PUCKER_ENERGY_ST]+=(inew[0]-iold[0]);
		lp.nn[PUCKER_ENERGY_ST]+=(fnew[0]-fold[0]);
	}

	fene=0; nene=0;
	fene+=lp.fn[PUCKER_ENERGY_ST];
	nene+=lp.nn[PUCKER_ENERGY_ST];
}





















