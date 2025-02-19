// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

double RNASimulation::calcsingletempene(Residueatoms restr){
	char aat=restr.aaa;
	int ss=restr.ss;
	int i, j;
	int atom;
	vector<vector<double>> RotMatix;
    vector<double> TranVect;
	vector<double> xyztmp(3,0);
	Point3s vec;
	double minrmsd=100000000000;
	double rmsd;

	vector<double> tmp(3,0);
	vector<vector<double> > xyz_list1(AccAtomParamNum[restr.aaa][1],tmp);
	vector<vector<double> > xyz_list2(AccAtomParamNum[restr.aaa][1],tmp);

	for(i=0;i<SINGLENTTEMPLATE[ss][aat].size();i++){
		for(j=0; j<AccAtomParamNum[restr.aaa][1];j++){
			atom=SNTTEMPLATEATOMS[aat][j];
			xyz_list1[j][0]=restr.ptv[atom].x;
			xyz_list1[j][1]=restr.ptv[atom].y;
			xyz_list1[j][2]=restr.ptv[atom].z;

			xyz_list2[j][0]=SINGLENTTEMPLATE[ss][aat][i][0].ptv[atom].x;
			xyz_list2[j][1]=SINGLENTTEMPLATE[ss][aat][i][0].ptv[atom].y;
			xyz_list2[j][2]=SINGLENTTEMPLATE[ss][aat][i][0].ptv[atom].z;
		}

		bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

		//for(int j=0; j<SNTTEMPLATEATOMS[aat].size();j++){
		rmsd=0;
		for(j=0; j<AccAtomParamNum[restr.aaa][1];j++){
			atom=SNTTEMPLATEATOMS[aat][j];
			bf.Transform(xyz_list2[j],RotMatix,TranVect,xyztmp);

			vec.x=restr.ptv[atom].x-xyztmp[0];
			vec.y=restr.ptv[atom].y-xyztmp[1];
			vec.z=restr.ptv[atom].z-xyztmp[2];
			rmsd+=bf.norm(vec);
		}

		if(minrmsd > rmsd){
			minrmsd=rmsd;
		}
	}

	return minrmsd;
}


void RNASimulation::fenergysinglenttemplate(vector<Residueatoms> decstr,int numseq,double &fene, Energyterms &lp, vector<vector<bool>> mapflagrpp,bool flagacc)
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
	int tmp;
	//int estartidx=BLENGTH_ENERGY_NUM+ILENGTH_ENERGY_NUM+INNERANGLE_ENERGY_NUM;

	if(!flagacc)
	{
		lp.fn[SINGLENTTEMP_ENERGY_ST]=0;

		for(i=0;i<numseq;i++)
		{
			lp.indf[SINGLENTTEMP_ENERGY_ST][i]=calcsingletempene(decstr[i]);
			lp.fn[SINGLENTTEMP_ENERGY_ST]+=lp.indf[SINGLENTTEMP_ENERGY_ST][i];
		}
	}
	else
	{
		int iold[SINGLENTTEMP_ENERGY_NUM],inew[SINGLENTTEMP_ENERGY_NUM];
		double fold[SINGLENTTEMP_ENERGY_NUM],fnew[SINGLENTTEMP_ENERGY_NUM];
		//for(i=0;i<SINGLENTTEMP_ENERGY_NUM;i++){
		iold[0]=0; inew[0]=0; fold[0]=0; fnew[0]=0;
		//}
		bool upflag;

		for(i=0;i<numseq;i++)
		{
			aat=decstr[i].aaa;
			upflag=false;
			for(j=0;j<numseq;j++)
			{
				if(mapflagrpp[i][j]){
					upflag=true;
					break;
				}
			}

			if(upflag)
			{
				fold[0]+=lp.indf[SINGLENTTEMP_ENERGY_ST][i];
				lp.indf[SINGLENTTEMP_ENERGY_ST][i]=0;
				lp.indf[SINGLENTTEMP_ENERGY_ST][i]=calcsingletempene(decstr[i]);
				fnew[0]+=lp.indf[SINGLENTTEMP_ENERGY_ST][i];
			}
		}

		lp.fn[SINGLENTTEMP_ENERGY_ST]+=(fnew[0]-fold[0]);
	}

	fene=lp.fn[SINGLENTTEMP_ENERGY_ST];
}

































