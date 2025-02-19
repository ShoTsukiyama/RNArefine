// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


double RNASimulation::calcbptempene(Residueatoms restr1, Residueatoms restr2){
	char aat1=restr1.aaa;
	char aat2=restr2.aaa;

	int i, j;
	int atom;
	vector<vector<double>> RotMatix;
    vector<double> TranVect;
	vector<double> xyztmp(3,0.0);
	Point3s vec;
	double minrmsd=100000000;
	double rmsd;


        //vector<double> tmp(3,0.0);
	int paramnum=AccAtomParamNum[aat1][2]+AccAtomParamNum[aat2][2];
	
	//vector<vector<double>> xyz_list1(paramnum,tmp);
	vector<vector<double>> xyz_list1;
        xyz_list1.resize(paramnum, vector<double>(3, 0.0));

        //vector<vector<double>> xyz_list2(paramnum,tmp);
	vector<vector<double>> xyz_list2;
        xyz_list2.resize(paramnum, vector<double>(3, 0.0));

	//cout<<"bp template"<<endl;
	//for(int j=0; j<SNTTEMPLATEATOMS[aat1].size();j++){
		//atom=SNTTEMPLATEATOMS[aat1][j];
		//cout<< "0" << "-" << j << "-" << restr1.ptv[atom].x << "-" << restr1.ptv[atom].y << "-" << restr1.ptv[atom].z<<endl;
	//}
	//for(int j=0; j<SNTTEMPLATEATOMS[aat2].size();j++){
		//atom=SNTTEMPLATEATOMS[aat2][j];
		//cout<< "1" << "-" << j << "-" << restr2.ptv[atom].x << "-" << restr2.ptv[atom].y << "-" << restr2.ptv[atom].z<<endl;
	//}
	
	for(i=0;i<BPLENTTEMPLATE[aat1][aat2].size();i++)
	{
		//for(int j=0; j<SNTTEMPLATEATOMS[aat].size();j++){
		int count = 0;
		for(j=0; j<AccAtomParamNum[restr1.aaa][2];j++){
			atom=SNTTEMPLATEATOMS[aat1][j];
			xyz_list1[count][0]=restr1.ptv[atom].x;
			xyz_list1[count][1]=restr1.ptv[atom].y;
			xyz_list1[count][2]=restr1.ptv[atom].z;

			xyz_list2[count][0]=BPLENTTEMPLATE[aat1][aat2][i][0].ptv[atom].x;
			xyz_list2[count][1]=BPLENTTEMPLATE[aat1][aat2][i][0].ptv[atom].y;
			xyz_list2[count][2]=BPLENTTEMPLATE[aat1][aat2][i][0].ptv[atom].z;
			count+=1;
		}

		for(j=0; j<AccAtomParamNum[restr2.aaa][2];j++){
			atom=SNTTEMPLATEATOMS[aat2][j];
			xyz_list1[count][0]=restr2.ptv[atom].x;
			xyz_list1[count][1]=restr2.ptv[atom].y;
			xyz_list1[count][2]=restr2.ptv[atom].z;

			xyz_list2[count][0]=BPLENTTEMPLATE[aat1][aat2][i][1].ptv[atom].x;
			xyz_list2[count][1]=BPLENTTEMPLATE[aat1][aat2][i][1].ptv[atom].y;
			xyz_list2[count][2]=BPLENTTEMPLATE[aat1][aat2][i][1].ptv[atom].z;
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

	if(aat1 == aat2)
	{
		for(i=0;i<BPLENTTEMPLATE[aat1][aat2].size();i++)
		{
			int count = 0;
			for(j=0; j<AccAtomParamNum[restr2.aaa][2];j++){
				atom=SNTTEMPLATEATOMS[aat1][j];
				xyz_list1[count][0]=restr2.ptv[atom].x;
				xyz_list1[count][1]=restr2.ptv[atom].y;
				xyz_list1[count][2]=restr2.ptv[atom].z;

				xyz_list2[count][0]=BPLENTTEMPLATE[aat1][aat2][i][0].ptv[atom].x;
				xyz_list2[count][1]=BPLENTTEMPLATE[aat1][aat2][i][0].ptv[atom].y;
				xyz_list2[count][2]=BPLENTTEMPLATE[aat1][aat2][i][0].ptv[atom].z;
				count+=1;
			}

			for(j=0; j<AccAtomParamNum[restr1.aaa][2];j++){
				atom=SNTTEMPLATEATOMS[aat2][j];
				xyz_list1[count][0]=restr1.ptv[atom].x;
				xyz_list1[count][1]=restr1.ptv[atom].y;
				xyz_list1[count][2]=restr1.ptv[atom].z;

				xyz_list2[count][0]=BPLENTTEMPLATE[aat1][aat2][i][1].ptv[atom].x;
				xyz_list2[count][1]=BPLENTTEMPLATE[aat1][aat2][i][1].ptv[atom].y;
				xyz_list2[count][2]=BPLENTTEMPLATE[aat1][aat2][i][1].ptv[atom].z;
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


void RNASimulation::fenergynoncbasepairtemp(vector<Residueatoms> decstr,int numseq,double &fene,Energyterms &lp, vector<Basecontact> &baserpp, vector<vector<bool>> mapflagrpp,bool flagacc){
	int i,j;
	int basei,basej;
	char aat1, aat2;
	int pairwiseidx;
	vector<double> indf(BASEPAIRTEMP_ENERGY_ST,0);
	double indftmp;


	if(!flagacc)
	{
		for(j=0;j<BASEPAIRTEMP_ENERGY_NUM;j++)
		{
			lp.fn[BASEPAIRTEMP_ENERGY_ST+j]=0;
		}
		
		for(i=0;i<BASEPAIRS.size();i++)
		{
			basei=BASEPAIRS[i].basepidx.first;//base idx 1
			basej=BASEPAIRS[i].basepidx.second;//base idx 2

			aat1=decstr[basei].aaa;
			aat2=decstr[basej].aaa;

			//if((aat1 == 'A' && aat2 == 'U') | (aat2 == 'A' && aat1 == 'U') | (aat1 == 'G' && aat2 == 'C') | (aat2 == 'G' && aat1 == 'C')) continue;

			for(j=0;j<BASEPAIRTEMP_ENERGY_NUM;j++)
			{
				lp.indf[BASEPAIRTEMP_ENERGY_ST+j][basei]=0;
				lp.indf[BASEPAIRTEMP_ENERGY_ST+j][basej]=0;
			}

			pairwiseidx=bpmapid2idx[basei*numseq+basej];

			if(basei < basej)
			{
				for(j=0;j<BASEPAIRTEMP_ENERGY_ST;j++)
				{
					indf[j]=0;
				}

				indftmp=calcbptempene(decstr[basei], decstr[basej]);
				indf[0]+=indftmp;

				for(j=0;j<BASEPAIRTEMP_ENERGY_NUM;j++)
				{
					//cout<<to_string(basei) + "-" + to_string(basej) + "-" + to_string(j) + "-" + to_string(indf[j]) + "-" + to_string(indn[j])<<endl;

					baserpp[pairwiseidx].tval[j]=indf[j];

					lp.indf[BASEPAIRTEMP_ENERGY_ST+j][basei]+=indf[j];
					lp.indf[BASEPAIRTEMP_ENERGY_ST+j][basej]+=indf[j];
					lp.fn[BASEPAIRTEMP_ENERGY_ST+j]+=indf[j];
				}
			}
		}
	}
	else
	{
		//int idx = 0;
		double fold[BASEPAIRTEMP_ENERGY_NUM],fnew[BASEPAIRTEMP_ENERGY_NUM];
		for(i=0;i<BASEPAIRTEMP_ENERGY_NUM;i++){
			fold[i]=0;fnew[i]=0;
		}

		for(i=0;i<BASEPAIRS.size();i++)
		{
			basei=BASEPAIRS[i].basepidx.first;//base idx 1
			basej=BASEPAIRS[i].basepidx.second;//base idx 2

			aat1=decstr[basei].aaa;
			aat2=decstr[basej].aaa;

			//if((aat1 == 'A' && aat2 == 'U') | (aat2 == 'A' && aat1 == 'U') | (aat1 == 'G' && aat2 == 'C') | (aat2 == 'G' && aat1 == 'C')) continue;

			pairwiseidx=bpmapid2idx[basei*numseq+basej];

			if(basei<basej)
			{
				//if(!mapflagrpp[basei][basej]) continue;
			
				for(j=0;j<BASEPAIRTEMP_ENERGY_NUM;j++)
				{
					fold[j]+=baserpp[pairwiseidx].tval[j];

					lp.indf[BASEPAIRTEMP_ENERGY_ST+j][basei]=lp.indf[BASEPAIRTEMP_ENERGY_ST+j][basei]-baserpp[pairwiseidx].tval[j];
					lp.indf[BASEPAIRTEMP_ENERGY_ST+j][basej]=lp.indf[BASEPAIRTEMP_ENERGY_ST+j][basej]-baserpp[pairwiseidx].tval[j];
					
					baserpp[pairwiseidx].tval[j]=0;
				}

				for(j=0;j<BASEPAIRTEMP_ENERGY_NUM;j++)
				{
					indf[j]=0;
				}

				indftmp=calcbptempene(decstr[basei], decstr[basej]);
				indf[0]+=indftmp;

				for(j=0;j<BASEPAIRTEMP_ENERGY_NUM;j++)
				{
					//cout<<to_string(basei) + "-" + to_string(basej) + "-" + to_string(j) + "-" + to_string(indf[j]) + "-" + to_string(indn[j])<<endl;
					baserpp[pairwiseidx].tval[j]=indf[j];

					lp.indf[BASEPAIRTEMP_ENERGY_ST+j][basei]+=baserpp[pairwiseidx].tval[j];
					lp.indf[BASEPAIRTEMP_ENERGY_ST+j][basej]+=baserpp[pairwiseidx].tval[j];
					
					fnew[j]+=baserpp[pairwiseidx].tval[j];
				}
			}
		}

		for(i=0;i<BASEPAIRTEMP_ENERGY_NUM;i++)
		{
			lp.fn[BASEPAIRTEMP_ENERGY_ST+i]+=(fnew[i]-fold[i]);
		}
	}

	fene=0;
	for(i=0;i<BASEPAIRTEMP_ENERGY_NUM;i++)
	{
		fene+=lp.fn[BASEPAIRTEMP_ENERGY_ST+i];//base
	}
}






















