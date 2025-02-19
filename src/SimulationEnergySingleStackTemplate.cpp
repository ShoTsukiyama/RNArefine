// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

double RNASimulation::calcstacktempene(Residueatoms restr1, Residueatoms restr2, int sttype){
	char aat1=restr1.aaa;
	char aat2=restr2.aaa;

	int i, j;
	int atom;
	vector<vector<double>> RotMatix;
    vector<double> TranVect;
	vector<double> xyztmp(3,0);
	Point3s vec;
	double minrmsd=100000000000;
	double rmsd;

	vector<double> tmp(3,0);
	vector<vector<double> > xyz_list1(AccAtomParamNum[restr1.aaa][2]+AccAtomParamNum[restr2.aaa][2],tmp);
    vector<vector<double> > xyz_list2(AccAtomParamNum[restr1.aaa][2]+AccAtomParamNum[restr2.aaa][2],tmp);
	vector<vector<double> > xyz_listtmp(15,tmp);

	//cout<<"bp template"<<endl;
	//for(int j=0; j<SNTTEMPLATEATOMS[aat1].size();j++){
		//atom=SNTTEMPLATEATOMS[aat1][j];
		//cout<< "0" << "-" << j << "-" << restr1.ptv[atom].x << "-" << restr1.ptv[atom].y << "-" << restr1.ptv[atom].z<<endl;
	//}
	//for(int j=0; j<SNTTEMPLATEATOMS[aat2].size();j++){
		//atom=SNTTEMPLATEATOMS[aat2][j];
		//cout<< "1" << "-" << j << "-" << restr2.ptv[atom].x << "-" << restr2.ptv[atom].y << "-" << restr2.ptv[atom].z<<endl;
	//}
	
	for(i=0;i<STACKTEMPLATE[sttype][aat1][aat2].size();i++)
	{
		//for(int j=0; j<SNTTEMPLATEATOMS[aat].size();j++){
		int count = 0;
		for(j=0; j<AccAtomParamNum[restr1.aaa][2];j++){
			atom=SNTTEMPLATEATOMS[aat1][j];
			xyz_list1[count][0]=restr1.ptv[atom].x;
			xyz_list1[count][1]=restr1.ptv[atom].y;
			xyz_list1[count][2]=restr1.ptv[atom].z;

			xyz_list2[count][0]=STACKTEMPLATE[sttype][aat1][aat2][i][0].ptv[atom].x;
			xyz_list2[count][1]=STACKTEMPLATE[sttype][aat1][aat2][i][0].ptv[atom].y;
			xyz_list2[count][2]=STACKTEMPLATE[sttype][aat1][aat2][i][0].ptv[atom].z;
			count+=1;
		}

		for(j=0; j<AccAtomParamNum[restr2.aaa][2];j++){
			atom=SNTTEMPLATEATOMS[aat2][j];
			xyz_list1[count][0]=restr2.ptv[atom].x;
			xyz_list1[count][1]=restr2.ptv[atom].y;
			xyz_list1[count][2]=restr2.ptv[atom].z;

			xyz_list2[count][0]=STACKTEMPLATE[sttype][aat1][aat2][i][1].ptv[atom].x;
			xyz_list2[count][1]=STACKTEMPLATE[sttype][aat1][aat2][i][1].ptv[atom].y;
			xyz_list2[count][2]=STACKTEMPLATE[sttype][aat1][aat2][i][1].ptv[atom].z;
			count+=1;
		}

		bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

		rmsd=0;
		count=0;
		//for(int j=0; j<SNTTEMPLATEATOMS[aat].size();j++){
		for(j=0; j<AccAtomParamNum[restr1.aaa][2];j++){
			atom=SNTTEMPLATEATOMS[aat1][j];
			bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);
			vec.x=restr1.ptv[atom].x-xyztmp[0];
			vec.y=restr1.ptv[atom].y-xyztmp[1];
			vec.z=restr1.ptv[atom].z-xyztmp[2];
			rmsd+=bf.norm(vec);
			count+=1;
		}

		for(j=0; j<AccAtomParamNum[restr2.aaa][2];j++){
			atom=SNTTEMPLATEATOMS[aat2][j];
			bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);

			vec.x=restr2.ptv[atom].x-xyztmp[0];
			vec.y=restr2.ptv[atom].y-xyztmp[1];
			vec.z=restr2.ptv[atom].z-xyztmp[2];
			rmsd+=bf.norm(vec);
			count+=1;
		}

		if(minrmsd > rmsd){
			minrmsd=rmsd;
		}
	}

	if(sttype == 0 && aat1 == aat2)
	{
		for(i=0;i<STACKTEMPLATE[sttype][aat1][aat2].size();i++)
		{
			int count = 0;
			for(j=0; j<AccAtomParamNum[restr2.aaa][2];j++){
				atom=SNTTEMPLATEATOMS[aat1][j];
				xyz_list1[count][0]=restr2.ptv[atom].x;
				xyz_list1[count][1]=restr2.ptv[atom].y;
				xyz_list1[count][2]=restr2.ptv[atom].z;

				xyz_list2[count][0]=STACKTEMPLATE[sttype][aat1][aat2][i][0].ptv[atom].x;
				xyz_list2[count][1]=STACKTEMPLATE[sttype][aat1][aat2][i][0].ptv[atom].y;
				xyz_list2[count][2]=STACKTEMPLATE[sttype][aat1][aat2][i][0].ptv[atom].z;
				count+=1;
			}

			for(j=0; j<AccAtomParamNum[restr1.aaa][2];j++){
				atom=SNTTEMPLATEATOMS[aat2][j];
				xyz_list1[count][0]=restr1.ptv[atom].x;
				xyz_list1[count][1]=restr1.ptv[atom].y;
				xyz_list1[count][2]=restr1.ptv[atom].z;

				xyz_list2[count][0]=STACKTEMPLATE[sttype][aat1][aat2][i][1].ptv[atom].x;
				xyz_list2[count][1]=STACKTEMPLATE[sttype][aat1][aat2][i][1].ptv[atom].y;
				xyz_list2[count][2]=STACKTEMPLATE[sttype][aat1][aat2][i][1].ptv[atom].z;
				count+=1;
			}

			bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

			rmsd=0;
			count=0;
			//for(int j=0; j<SNTTEMPLATEATOMS[aat].size();j++){
			for(j=0; j<AccAtomParamNum[restr2.aaa][2];j++){
				atom=SNTTEMPLATEATOMS[aat1][j];
				bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);
				vec.x=restr2.ptv[atom].x-xyztmp[0];
				vec.y=restr2.ptv[atom].y-xyztmp[1];
				vec.z=restr2.ptv[atom].z-xyztmp[2];
				rmsd+=bf.norm(vec);
				count+=1;
			}

			for(j=0; j<AccAtomParamNum[restr1.aaa][2];j++){
				atom=SNTTEMPLATEATOMS[aat2][j];
				bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);

				vec.x=restr1.ptv[atom].x-xyztmp[0];
				vec.y=restr1.ptv[atom].y-xyztmp[1];
				vec.z=restr1.ptv[atom].z-xyztmp[2];
				rmsd+=bf.norm(vec);
				count+=1;
			}

			if(minrmsd > rmsd){
				minrmsd=rmsd;
			}
		}
	}
	return minrmsd;
}



void RNASimulation::fenergystacktemplate(vector<Residueatoms> decstr,int numseq,double &fene,Energyterms &lp, vector<Stackcontact> &stackrpp, vector<vector<bool>> mapflagrpp,bool flagacc)
{
	int i,j;
	int basei,basej;
	int pairwiseidx;
	double dih, len;
	char aat1, aat2;
	int atom1, atom2, atom3, atom4;
	int atom11,atom12,atom13,atom21,atom22,atom23;
	Point3s tmp1,tmp2,tmp3,tmp4;
	Point3s res1N91, res1C42, res1C86, res2N91, res2C42, res2C86;
	double r, alpha, beta, gamma;
	int binidx1, binidx2, binidx3;
	int binnum1, binnum2, binnum3;
	Point3s pt;
	int stacktype;
	double indftmp;
	int indntmp;
	vector<double> indf(STACKINGTEMP_ENERGY_NUM,0);
	vector<int> indn(STACKINGTEMP_ENERGY_NUM,0);
	Point3s basecenter1, basecenter2;
	double baseangle, baseangle1, baseangle2;
	Point3s c2cvec, z1, z2;
	double tmp;

	if(!flagacc)
	{
		for(j=0;j<STACKINGTEMP_ENERGY_NUM;j++)
		{
			lp.fn[STACKINGTEMP_ENERGY_ST+j]=0;
			lp.nn[STACKINGTEMP_ENERGY_ST+j]=0;


			for(i=0;i<numseq;i++){
				lp.indf[STACKINGTEMP_ENERGY_ST+j][i]=0;
				lp.indf[STACKINGTEMP_ENERGY_ST+j][i]=0;
				lp.indn[STACKINGTEMP_ENERGY_ST+j][i]=0;
				lp.indn[STACKINGTEMP_ENERGY_ST+j][i]=0;
			}
		}

		for(i=0;i<STACKINGPAIRS.size();i++)
		{
			basei=STACKINGPAIRS[i].basepidx.first;//base idx 1
			basej=STACKINGPAIRS[i].basepidx.second;//base idx 2
			stacktype=STACKINGPAIRS[i].sttype;

			pairwiseidx=bsmapid2idx[basei*numseq+basej];


			if(basei < basej)
			{
				aat1=decstr[basei].aaa;
				aat2=decstr[basej].aaa;
                                

				indf[0]=0;
				indf[0] += calcstacktempene(decstr[basei], decstr[basej], stacktype);
				

				for(j=0;j<STACKINGTEMP_ENERGY_NUM;j++)
				{
					stackrpp[pairwiseidx].tval[j]=indf[j];
					lp.indf[STACKINGTEMP_ENERGY_ST+j][basei]+=indf[j];
					lp.indf[STACKINGTEMP_ENERGY_ST+j][basej]+=indf[j];
					lp.fn[STACKINGTEMP_ENERGY_ST+j]+=indf[j];
				}
			}
		}
	}
	else
	{
		//int idx = 0;
		double fold[STACKINGTEMP_ENERGY_NUM],fnew[STACKINGTEMP_ENERGY_NUM];
		for(i=0;i<STACKINGTEMP_ENERGY_NUM;i++){
			fold[i]=0;fnew[i]=0;
		}


		for(i=0;i<STACKINGPAIRS.size();i++)
		{
			basei=STACKINGPAIRS[i].basepidx.first;//base idx 1
			basej=STACKINGPAIRS[i].basepidx.second;//base idx 2
			stacktype=STACKINGPAIRS[i].sttype;
			
			pairwiseidx=bsmapid2idx[basei*numseq+basej];

			if(basei<basej)
			{
				if(!mapflagrpp[basei][basej]) continue;
			
				for(j=0;j<STACKINGTEMP_ENERGY_NUM;j++)
				{
					fold[j]+=stackrpp[pairwiseidx].tval[j];

					lp.indf[STACKINGTEMP_ENERGY_ST+j][basei]=lp.indf[STACKINGTEMP_ENERGY_ST+j][basei]-stackrpp[pairwiseidx].pval[j];
					lp.indf[STACKINGTEMP_ENERGY_ST+j][basej]=lp.indf[STACKINGTEMP_ENERGY_ST+j][basej]-stackrpp[pairwiseidx].pval[j];
					lp.indn[STACKINGTEMP_ENERGY_ST+j][basei]=lp.indn[STACKINGTEMP_ENERGY_ST+j][basei]-stackrpp[pairwiseidx].pcw[j];
					lp.indn[STACKINGTEMP_ENERGY_ST+j][basej]=lp.indn[STACKINGTEMP_ENERGY_ST+j][basej]-stackrpp[pairwiseidx].pcw[j];

					stackrpp[pairwiseidx].tval[j]=0;
				}

				aat1=decstr[basei].aaa;
				aat2=decstr[basej].aaa;

				indf[0]=0;
				indf[0] += calcstacktempene(decstr[basei], decstr[basej], stacktype);

				for(j=0;j<STACKINGTEMP_ENERGY_NUM;j++)
				{
					stackrpp[pairwiseidx].tval[j]+=indf[j];

					lp.indf[STACKINGTEMP_ENERGY_ST+j][basei]+=indf[j];
					lp.indf[STACKINGTEMP_ENERGY_ST+j][basej]+=indf[j];

					fnew[j]+=stackrpp[pairwiseidx].tval[j];
				}
			}
		}

		for(i=0;i<STACKINGTEMP_ENERGY_NUM;i++)
		{
			lp.fn[STACKINGTEMP_ENERGY_ST+i]+=(fnew[i]-fold[i]);
		}
	}

	fene=0;
	for(i=0;i<STACKINGTEMP_ENERGY_NUM;i++)
	{
		fene+=lp.fn[STACKINGTEMP_ENERGY_ST+i];
	}
}






















