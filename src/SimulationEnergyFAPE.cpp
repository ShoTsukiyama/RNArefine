// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

double RNASimulation::calcFAPEthresh(double dist){
	return 3+(150/(5+exp(0.5*(dist-10))));
	//return 1+(30/(1+exp(-0.5*(dist-20))));
}

void RNASimulation::fapesingleref(vector<Atom> proseq,vector<Boneinfo> bb, int numbb,vector<Boneinfo> ibb, int inumbb,int numseq)
{
	//map<char,map<int,vector<double>>> basecoordindall = {
	//{'A',{{0,{3.063,8.025,-4.135}},{3,{6.719,6.889,-3.258}},{4,{6.776,5.642,-4.140}}}},
	//{'C',{{0,{6.913,5.099,-6.683}},{3,{9.376,2.167,-5.806}},{4,{8.750,1.087,-6.688}}}},
	//{'G',{{0,{-1.758,-8.408,1.587}},{3,{1.932,-9.427,0.710}},{4,{2.653,-8.408,1.592}}}},
	//{'U',{{0,{3.063,-8.025,4.135}},{3,{6.719,-6.889,3.258}},{4,{6.776,-5.642,4.140}}}}
	//};

	map<char,map<int,vector<double>>> basecoordindall = {
		{'A',{{0,{2.3240912,-0.80317456,-1.5209152}},{3,{-0.8031747,1.3708649,-0.56768906}},{10,{-1.520915,-0.5676891,2.088603}}}},
		{'C',{{0,{2.35021,-0.7880346,-1.5621773}},{3,{-0.78803474,1.3541472,-0.5661155}},{11,{-1.5621774,-0.56611556,2.1282897}}}},
		{'G',{{0,{2.3250363,-0.80280226,-1.5222337}},{3,{-0.8028024,1.3702725,-0.5674678}},{10,{-1.5222336,-0.5674678,2.0897026}}}},
		{'U',{{0,{2.3493478,-0.7881313,-1.5612125}},{3,{-0.7881318,1.3542684,-0.5661326}},{11,{-1.5612122,-0.5661328,2.1273444}}}}
	};

	//map<char,map<int,vector<double>>> basecoordindall = {
	//{'A',{{0,{3.063,8.025,-4.135}},{3,{6.719,6.889,-3.258}},{10,{5.666,4.158,-1.647}}}},
	//{'C',{{0,{6.913,5.099,-6.683}},{3,{9.376,2.167,-5.806}},{11,{7.014,0.438,-4.195}}}},
	//{'G',{{0,{-1.758,-8.408,1.587}},{3,{1.932,-9.427,0.710}},{10,{2.522,-6.560,-0.901}}}},
	//{'U',{{0,{3.063,-8.025,4.135}},{3,{6.719,-6.889,3.258}},{11,{5.666,-4.158,1.647}}}}
	//};

	initnpaa3=0;
	int i,j,k;

	int aidx;
	vector<double> bcoor;
	vector<double> tmp(FAPE_BASEATOM_NUM,0);
    vector<vector<double> > xyz_list1(FAPE_BASEATOM_NUM,tmp);
    vector<vector<double> > xyz_list2(FAPE_BASEATOM_NUM,tmp);
	vector<vector<double>> xyz_list3(FAPE_BASEATOM_NUM,tmp);
	vector<double> xyz_tmp(tmp);
	vector<vector<double> > RotMatix;
    vector<double> TranVect;
	map<int,vector<double>> basec;
	Point3d pcd;
	//lentmp.resize(refidx+1);
	//ntflag.resize(refidx+1);
	//ntstr.resize(refidx+1);
	
	//resind in simu [0,numseq-1]
	//ibb           resind don't need from 1
	//bb        similar range as ibb
	int bbmin=1000000;
	int bbmax=-10000000;

	//残基idの最大値(bbmax)と最小値(bbmin)を求める
	for(i=0;i<inumbb;i++)
	{
		if(bbmin>ibb[i].resind) bbmin=ibb[i].resind;
		if(bbmax<ibb[i].resind) bbmax=ibb[i].resind;
	}

	//残基数+5000の配列を用意し、0~残基数+1の要素にはiを格納
	//int *bbind=new int[5000+bbmax-bbmin+1];
	vector<int> bbind(5000+bbmax-bbmin+1);
	for(i=0;i<bbmax-bbmin+1;i++) bbind[i]=-1; 
	for(i=0;i<inumbb;i++) bbind[ibb[i].resind-bbmin]=i;

	initpaa3.resize(initpaa3.size()+numseq*numseq);
	
	k=initnpaa3;
	for(i=0;i<numbb;i++)
	{
		if(bb[i].indatom[0]==-1 || bb[i].indatom[3]==-1 || bb[i].indatom[4]==-1) continue;
		//if(bb[i].indatom[0]==-1 || bb[i].indatom[3]==-1) continue;
		//if((bb[i].resid=='A' || bb[i].resid=='G') && bb[i].indatom[10]==-1) continue;
		//if((bb[i].resid=='C' || bb[i].resid=='U') && bb[i].indatom[11]==-1) continue;
		if(bb[i].resind<bbmin || bb[i].resind>bbmax) continue;
		if(bbind[bb[i].resind-bbmin]==-1) continue;

		basec=basecoordindall[bb[i].resid];

		int count=0;
		for(auto it = basec.begin(); it != basec.end(); ++it){
        	aidx = it->first;
        	bcoor = it->second;

			//cout<<to_string(i) + "|" + to_string(aidx)<<endl;
			//cout<<to_string(proseq[bb[i].indatom[aidx]].x) + "|" + to_string(proseq[bb[i].indatom[aidx]].y) + "|" + to_string(proseq[bb[i].indatom[aidx]].z)<<endl;
			xyz_list1[count][0]=proseq[bb[i].indatom[aidx]].x;
        	xyz_list1[count][1]=proseq[bb[i].indatom[aidx]].y;
        	xyz_list1[count][2]=proseq[bb[i].indatom[aidx]].z;

			//cout<<to_string(i) +"||  "+ to_string(proseq[bb[i].indatom[aidx]].x)<<endl;

			xyz_list2[count][0]=bcoor[0];
        	xyz_list2[count][1]=bcoor[1];
        	xyz_list2[count][2]=bcoor[2];
			count+=1;
		}

		bf.RotateCoor(xyz_list2,xyz_list1,RotMatix,TranVect);

		for(j=0;j<numbb;j++)
		{
			if(bb[j].indatom[0]==-1 || bb[j].indatom[3]==-1 || bb[j].indatom[4]==-1) continue;
			//if(bb[j].indatom[0]==-1 || bb[j].indatom[3]==-1) continue;
			//if((bb[j].resid=='A' || bb[j].resid=='G') && bb[j].indatom[10]==-1) continue;
			//if((bb[j].resid=='C' || bb[j].resid=='U') && bb[j].indatom[11]==-1) continue;
			if(bb[j].resind<bbmin || bb[j].resind>bbmax) continue;
			if(bbind[bb[j].resind-bbmin]==-1) continue;
			if(fabs(bbind[bb[j].resind-bbmin]-bbind[bb[i].resind-bbmin])<1) continue;

			basec=basecoordindall[bb[j].resid];

			initpaa3[k].cood.resize(basec.size());
			count=0;
			for(auto it = basec.begin(); it != basec.end(); ++it)
			{
				aidx = it->first;
				xyz_tmp[0]=proseq[bb[j].indatom[aidx]].x;
				xyz_tmp[1]=proseq[bb[j].indatom[aidx]].y;
				xyz_tmp[2]=proseq[bb[j].indatom[aidx]].z;
				bf.Transform(xyz_tmp,RotMatix,TranVect,initpaa3[k].cood[count]);
				count+=1;
			}

		//	if(bb[j].resind-bb[i].resind<2) continue;
			
			initpaa3[k].ispair=true;
			initpaa3[k].pbin=bbind[bb[i].resind-bbmin];
			initpaa3[k].tnum=bbind[bb[j].resind-bbmin];

			pcd.x=proseq[bb[j].indatom[3]].x-proseq[bb[i].indatom[3]].x; 
			pcd.y=proseq[bb[j].indatom[3]].y-proseq[bb[i].indatom[3]].y; 
			pcd.z=proseq[bb[j].indatom[3]].z-proseq[bb[i].indatom[3]].z;
			initpaa3[k].threshold=calcFAPEthresh(bf.norm(pcd));
			k++;
		}

	}
	initnpaa3=k;
	k=initnpaa3;
	if(k<numseq*numseq) k=numseq*numseq;
	initpaa3.resize(k);
	//strpval.resize(k);
	//paa3=(pairaa *)realloc(paa3,k*sizeof(pairaa));

	//refinestructureとreference structureで共通したの残基数をlentmpに格納
	initlentmp=0; 
	for(i=0;i<numbb;i++)
	{
		if(bb[i].resind<bbmin || bb[i].resind>bbmax) continue;
		if(bbind[bb[i].resind-bbmin]==-1) continue;
		initlentmp++;
	}
	//refinestructureとreference structureで共通したの残基であるかどうか(共通:true、共通でない:false)をntflagに格納
	//refinestructureとreference structureで共通したの残基のC4'の座標をntstrに格納
	initntflag.resize(inumbb);
	initntstr.resize(inumbb*FAPE_BASEATOM_NUM);

	int maxbbind=0;
	j=0;
	for(i=0;i<numbb;i++)
	{

		if(bb[i].resind<bbmin || bb[i].resind>bbmax) continue;
		if(bbind[bb[i].resind-bbmin]==-1) continue;

	//	ntflag[bb[i].resind-1]=true;
		if(maxbbind<bbind[bb[i].resind-bbmin]) maxbbind=bbind[bb[i].resind-bbmin];
		initntflag[bbind[bb[i].resind-bbmin]]=true;
		initntstr[j]=proseq[bb[i].indatom[3]].x;
		initntstr[initlentmp+j]=proseq[bb[i].indatom[3]].y;
		initntstr[2*initlentmp+j]=proseq[bb[i].indatom[3]].z;
		j++;
	}
	initntstr.resize(FAPE_BASEATOM_NUM*initlentmp);
	initntflag.resize(maxbbind);

	bbind.clear();
	bbind.shrink_to_fit();

}

//getinitcontact:: reference structureとrefined structurenの共通の残基を特定し、共通した残基についてreference structureの残機間の距離を算出
//refrence配列のatomの情報、bb:reference配列のbackbornの情報、numbb:reference配列のbackbornの残基数、refined配列のbackborn atomの情報、refined配列のbackbornの残基数
void RNASimulation::fapemultipleref(vector<Atom> proseq,vector<Boneinfo> bb, int numbb,vector<Boneinfo> ibb, int inumbb,int numseq, int refidx)
{
	//map<char,map<int,vector<double>>> basecoordindall = {
	//{'A',{{0,{2.32409,-0.803175,-1.52092}},{3,{-0.803175,1.37086,-0.567689}},{10,{-1.52092,-0.567689,2.0886}}}},
	//{'C',{{0,{2.35021,-0.788035,-1.56218}},{3,{-0.788035,1.35415,-0.566115}},{11,{-1.56218,-0.566116,2.12829}}}},
	//{'G',{{0,{2.32504,-0.802802,-1.52223}},{3,{-0.802802,1.37027,-0.567468}},{10,{-1.52223,-0.567468,2.0897}}}},
	//{'U',{{0,{2.34935,-0.788131,-1.56121}},{3,{-0.788132,1.35427,-0.566133}},{11,{-1.56121,-0.566133,2.12734}}}}
	//};

	//map<char,map<int,vector<double>>> basecoordindall = {
	//{'A',{{0,{3.063,8.025,-4.135}},{3,{6.719,6.889,-3.258}},{4,{6.776,5.642,-4.140}}}},
	//{'C',{{0,{6.913,5.099,-6.683}},{3,{9.376,2.167,-5.806}},{4,{8.750,1.087,-6.688}}}},
	//{'G',{{0,{-1.758,-8.408,1.587}},{3,{1.932,-9.427,0.710}},{4,{2.653,-8.408,1.592}}}},
	//{'U',{{0,{3.063,-8.025,4.135}},{3,{6.719,-6.889,3.258}},{4,{6.776,-5.642,4.140}}}}
	//};

	map<char,map<int,vector<double>>> basecoordindall = {
		{'A',{{0,{2.3240912,-0.80317456,-1.5209152}},{3,{-0.8031747,1.3708649,-0.56768906}},{10,{-1.520915,-0.5676891,2.088603}}}},
		{'C',{{0,{2.35021,-0.7880346,-1.5621773}},{3,{-0.78803474,1.3541472,-0.5661155}},{11,{-1.5621774,-0.56611556,2.1282897}}}},
		{'G',{{0,{2.3250363,-0.80280226,-1.5222337}},{3,{-0.8028024,1.3702725,-0.5674678}},{10,{-1.5222336,-0.5674678,2.0897026}}}},
		{'U',{{0,{2.3493478,-0.7881313,-1.5612125}},{3,{-0.7881318,1.3542684,-0.5661326}},{11,{-1.5612122,-0.5661328,2.1273444}}}}
	};
	

	//map<char,map<int,vector<double>>> basecoordindall = {
	//{'A',{{0,{3.063,8.025,-4.135}},{3,{6.719,6.889,-3.258}},{10,{5.666,4.158,-1.647}}}},
	//{'C',{{0,{6.913,5.099,-6.683}},{3,{9.376,2.167,-5.806}},{11,{7.014,0.438,-4.195}}}},
	//{'G',{{0,{-1.758,-8.408,1.587}},{3,{1.932,-9.427,0.710}},{10,{2.522,-6.560,-0.901}}}},
	//{'U',{{0,{3.063,-8.025,4.135}},{3,{6.719,-6.889,3.258}},{11,{5.666,-4.158,1.647}}}}
	//};

	npaa3.resize(refidx+1);
	lentmp.resize(refidx+1);
	ntflag.resize(refidx+1);
	ntstr.resize(refidx+1);
	npaa3[refidx]=0;
	int i,j,k;

	int aidx;
	vector<double> bcoor;
	vector<double> tmp(FAPE_BASEATOM_NUM,0);
    vector<vector<double> > xyz_list1(FAPE_BASEATOM_NUM,tmp);
    vector<vector<double> > xyz_list2(FAPE_BASEATOM_NUM,tmp);
	vector<vector<double>> xyz_list3(FAPE_BASEATOM_NUM,tmp);
	vector<double> xyz_tmp(tmp);
	vector<vector<double> > RotMatix;
    vector<double> TranVect;
	map<int,vector<double>> basec;
	Point3d pcd;
	
	int bbmin=1000000;
	int bbmax=-10000000;

	//残基idの最大値(bbmax)と最小値(bbmin)を求める
	for(i=0;i<inumbb;i++)
	{
		if(bbmin>ibb[i].resind) bbmin=ibb[i].resind;
		if(bbmax<ibb[i].resind) bbmax=ibb[i].resind;
	}

	//残基数+5000の配列を用意し、0~残基数+1の要素にはiを格納
	//int *bbind=new int[5000+bbmax-bbmin+1];
	vector<int> bbind(5000+bbmax-bbmin+1);
	for(i=0;i<bbmax-bbmin+1;i++) bbind[i]=-1; 
	for(i=0;i<inumbb;i++) bbind[ibb[i].resind-bbmin]=i;

	if(paa3.size() < refidx + 1) paa3.resize(refidx + 1);
	paa3[refidx].resize(npaa3[refidx]+numseq*numseq);
	
	BasicFunc bf;
	k=npaa3[refidx];
	for(i=0;i<numbb;i++)
	{
		basec=basecoordindall[bb[i].resid];
		if(bb[i].indatom[0]==-1 || bb[i].indatom[3]==-1 || bb[i].indatom[4]==-1) continue;
		//if(bb[i].indatom[0]==-1 || bb[i].indatom[3]==-1) continue;
		//if((bb[i].resid=='A' || bb[i].resid=='G') && bb[i].indatom[10]==-1) continue;
		//if((bb[i].resid=='C' || bb[i].resid=='U') && bb[i].indatom[11]==-1) continue;
		if(bb[i].resind<bbmin || bb[i].resind>bbmax) continue;
		if(bbind[bb[i].resind-bbmin]==-1) continue;

		int count=0;
		for(auto it = basec.begin(); it != basec.end(); ++it){
        	aidx = it->first;
        	bcoor = it->second;

			xyz_list1[count][0]=proseq[bb[i].indatom[aidx]].x;
        	xyz_list1[count][1]=proseq[bb[i].indatom[aidx]].y;
        	xyz_list1[count][2]=proseq[bb[i].indatom[aidx]].z;

			xyz_list2[count][0]=bcoor[0];
        	xyz_list2[count][1]=bcoor[1];
        	xyz_list2[count][2]=bcoor[2];
			count+=1;
		}

		bf.RotateCoor(xyz_list2,xyz_list1,RotMatix,TranVect);

		for(j=0;j<numbb;j++)
		{
			if(bb[j].indatom[0]==-1 || bb[j].indatom[3]==-1 || bb[j].indatom[4]==-1) continue;
			//if(bb[j].indatom[0]==-1 || bb[j].indatom[3]==-1) continue;
			//if((bb[j].resid=='A' || bb[j].resid=='G') && bb[j].indatom[10]==-1) continue;
			//if((bb[j].resid=='C' || bb[j].resid=='U') && bb[j].indatom[11]==-1) continue;
			if(bb[j].resind<bbmin || bb[j].resind>bbmax) continue;
			if(bbind[bb[j].resind-bbmin]==-1) continue;
			if(fabs(bbind[bb[j].resind-bbmin]-bbind[bb[i].resind-bbmin])<1) continue;

			basec=basecoordindall[bb[j].resid];

			paa3[refidx][k].cood.resize(basec.size());
			count=0;
			for(auto it = basec.begin(); it != basec.end(); ++it)
			{
				aidx = it->first;
				xyz_tmp[0]=proseq[bb[j].indatom[aidx]].x;
				xyz_tmp[1]=proseq[bb[j].indatom[aidx]].y;
				xyz_tmp[2]=proseq[bb[j].indatom[aidx]].z;
				bf.Transform(xyz_tmp,RotMatix,TranVect,paa3[refidx][k].cood[count]);
				count+=1;
			}
			
			paa3[refidx][k].ispair=true;
			paa3[refidx][k].pbin=bbind[bb[i].resind-bbmin];
			paa3[refidx][k].tnum=bbind[bb[j].resind-bbmin];

			pcd.x=proseq[bb[j].indatom[3]].x-proseq[bb[i].indatom[3]].x; 
			pcd.y=proseq[bb[j].indatom[3]].y-proseq[bb[i].indatom[3]].y; 
			pcd.z=proseq[bb[j].indatom[3]].z-proseq[bb[i].indatom[3]].z;
			paa3[refidx][k].threshold=calcFAPEthresh(bf.norm(pcd));

			k++;
		}

	}
	npaa3[refidx]=k;
	k=npaa3[refidx];
	if(k<numseq*numseq) k=numseq*numseq;
	paa3[refidx].resize(k);

	//refinestructureとreference structureで共通したの残基数をlentmpに格納
	lentmp[refidx]=0; 
	for(i=0;i<numbb;i++)
	{
		if(bb[i].resind<bbmin || bb[i].resind>bbmax) continue;
		if(bbind[bb[i].resind-bbmin]==-1) continue;
		lentmp[refidx]++;
	}
	//refinestructureとreference structureで共通したの残基であるかどうか(共通:true、共通でない:false)をntflagに格納
	//refinestructureとreference structureで共通したの残基のC4'の座標をntstrに格納
	ntflag[refidx].resize(inumbb);
	ntstr[refidx].resize(inumbb*FAPE_BASEATOM_NUM);

	int maxbbind=0;
	j=0;
	for(i=0;i<numbb;i++)
	{

		if(bb[i].resind<bbmin || bb[i].resind>bbmax) continue;
		if(bbind[bb[i].resind-bbmin]==-1) continue;

		if(maxbbind<bbind[bb[i].resind-bbmin]) maxbbind=bbind[bb[i].resind-bbmin];
		ntflag[refidx][bbind[bb[i].resind-bbmin]]=true;
		ntstr[refidx][j]=proseq[bb[i].indatom[3]].x;
		ntstr[refidx][lentmp[refidx]+j]=proseq[bb[i].indatom[3]].y;
		ntstr[refidx][2*lentmp[refidx]+j]=proseq[bb[i].indatom[3]].z;
		j++;
	}
	ntstr[refidx].resize(FAPE_BASEATOM_NUM*lentmp[refidx]);
	ntflag[refidx].resize(maxbbind);

	bbind.clear();
	bbind.shrink_to_fit();

}

void RNASimulation::fapeenergysingle(vector<Residueatoms> decstr,int numseq,double &fene,Energyterms &lp,FAPEobj &faperpp,bool flagacc)
{
	double totene=0;
	int i,j,k,l;
	Point3d tp;
	double tval;
	int aidx;
	vector<double> bcoor;
	vector<double> tmp(FAPE_BASEATOM_NUM,0);
    vector<vector<double> > xyz_list1(FAPE_BASEATOM_NUM,tmp);
    vector<vector<double> > xyz_list2(FAPE_BASEATOM_NUM,tmp);
	vector<vector<double>> xyz_list3(FAPE_BASEATOM_NUM,tmp);
	vector<double> xyz_tmp(tmp);
	vector<vector<double> > RotMatix;
    vector<double> TranVect;
	map<int,vector<double>> basec;
	double indf=0;
	double ttmp;

	if(!flagacc)
	{
		lp.fn[FAPE_ENERGY_ST]=0;

		for(i=0;i<initnpaa3;i++)
		{
			k=initpaa3[i].tnum*numseq+initpaa3[i].pbin;

			faperpp.fapeene1[k]=0;
			for(j=0;j<FAPE_BASEATOM_NUM;j++)
			{
				tp.x=faperpp.strpval[k][j].x-initpaa3[i].cood[j].x;
				tp.y=faperpp.strpval[k][j].y-initpaa3[i].cood[j].y;
				tp.z=faperpp.strpval[k][j].z-initpaa3[i].cood[j].z;
				tval=bf.norm(tp);

				if(tval > FAPE_THRESH) tval=FAPE_THRESH;

				faperpp.fapeene1[k]+=tval;
			}
			indf+=faperpp.fapeene1[k];
		}
		lp.fn[FAPE_ENERGY_ST]=indf;
	}
	else{
		double eold=0; double enew=0;

		for(i=0;i<initnpaa3;i++)
		{
			k=initpaa3[i].tnum*numseq+initpaa3[i].pbin;

			if(faperpp.flap1[k])
			{
				eold+=faperpp.fapeene1[k];

				faperpp.fapeene1[k]=0;
				for(j=0;j<FAPE_BASEATOM_NUM;j++){
					tp.x=faperpp.strpval[k][j].x-initpaa3[i].cood[j].x;
					tp.y=faperpp.strpval[k][j].y-initpaa3[i].cood[j].y;
					tp.z=faperpp.strpval[k][j].z-initpaa3[i].cood[j].z;
					tval=bf.norm(tp);

					if(tval > FAPE_THRESH) tval=FAPE_THRESH;
					faperpp.fapeene1[k]+=tval;
				}

				enew+=faperpp.fapeene1[k];
				faperpp.flap1[k]=false;
			}
		}
		
		lp.fn[FAPE_ENERGY_ST]+=(enew-eold);
	}
	fene=lp.fn[FAPE_ENERGY_ST];
}

void RNASimulation::fapeenergymultiple(vector<Residueatoms> decstr,int numseq,double &fene,Energyterms &lp,FAPEobj &faperpp,bool flagacc)
{
	double totene=0;
	int i,j,k,l;
	Point3d tp;
	double tval;
	int aidx;
	vector<double> bcoor;
	vector<double> tmp(FAPE_BASEATOM_NUM,0);
    vector<vector<double> > xyz_list1(FAPE_BASEATOM_NUM,tmp);
    vector<vector<double> > xyz_list2(FAPE_BASEATOM_NUM,tmp);
	vector<vector<double>> xyz_list3(FAPE_BASEATOM_NUM,tmp);
	vector<double> xyz_tmp(tmp);
	vector<vector<double> > RotMatix;
    vector<double> TranVect;
	map<int,vector<double>> basec;
	double threshold;
	double indf=0;
	
	if(!flagacc)
	{
		lp.fn[FAPE_ENERGY_ST+1]=0;
		for(l=0;l<npaa3.size();l++)
		{
			for(i=0;i<npaa3[l];i++)
			{
				k=paa3[l][i].tnum*numseq+paa3[l][i].pbin;

				faperpp.fapeene2[l][k]=0;
				for(j=0;j<FAPE_BASEATOM_NUM;j++)
				{
					tp.x=faperpp.strpval[k][j].x-paa3[l][i].cood[j].x;
					tp.y=faperpp.strpval[k][j].y-paa3[l][i].cood[j].y;
					tp.z=faperpp.strpval[k][j].z-paa3[l][i].cood[j].z;
					tval=bf.norm(tp);
					if(tval > FAPE_THRESH) tval=FAPE_THRESH;
					faperpp.fapeene2[l][k]+=tval;
				}
				indf+=faperpp.fapeene2[l][k];
			}
		}
		lp.fn[FAPE_ENERGY_ST+1]=indf;
	}else{
		double eold=0; 
		double enew=0;

		for(l=0;l<npaa3.size();l++)
		{
			for(i=0;i<npaa3[l];i++)
			{
				k=paa3[l][i].tnum*numseq+paa3[l][i].pbin;

				if(faperpp.flap2[k])
				{
					eold+=faperpp.fapeene2[l][k];

					faperpp.fapeene2[l][k]=0;
					for(j=0;j<FAPE_BASEATOM_NUM;j++)
					{
						tp.x=faperpp.strpval[k][j].x-paa3[l][i].cood[j].x;
						tp.y=faperpp.strpval[k][j].y-paa3[l][i].cood[j].y;
						tp.z=faperpp.strpval[k][j].z-paa3[l][i].cood[j].z;
						tval=bf.norm(tp);
						if(tval > FAPE_THRESH) tval=FAPE_THRESH;
						faperpp.fapeene2[l][k]+=tval;
					}
					enew+=faperpp.fapeene2[l][k];
				}
			}
		}

		lp.fn[FAPE_ENERGY_ST+1]+=(enew-eold);
		for(i=0;i<numseq;i++){
			for(j=0;j<numseq;j++){
				k=i*numseq+j;
				faperpp.flap2[k]=false;
			}
		}
	}
	fene=lp.fn[FAPE_ENERGY_ST+1];
}

void RNASimulation::getfape(vector<Residueatoms> decstr,int numseq,FAPEobj &faperpp,bool flagacc){
	int i,j,k;
	int aidx;
	int count;
	vector<double> bcoor;
	vector<double> tmp(FAPE_BASEATOM_NUM,0);
    vector<vector<double> > xyz_list1(FAPE_BASEATOM_NUM,tmp);
    vector<vector<double> > xyz_list2(FAPE_BASEATOM_NUM,tmp);
	vector<double> xyz_tmp(tmp);
	vector<vector<double> > RotMatix;
    vector<double> TranVect;
	map<int,vector<double>> basec;

	//map<char,map<int,vector<double>>> basecoord = {
	//{'A',{{0,{3.063,8.025,-4.135}},{5,{6.719,6.889,-3.258}},{7,{6.776,5.642,-4.140}}}},
	//{'C',{{0,{6.913,5.099,-6.683}},{5,{9.376,2.167,-5.806}},{7,{8.750,1.087,-6.688}}}},
	//{'G',{{0,{-1.758,-8.408,1.587}},{5,{1.932,-9.427,0.710}},{7,{2.653,-8.408,1.592}}}},
	//{'U',{{0,{3.063,-8.025,4.135}},{5,{6.719,-6.889,3.258}},{7,{6.776,-5.642,4.140}}}}
	//};

	map<char,map<int,vector<double>>> basecoord = {
		{'A',{{0,{2.3240912,-0.80317456,-1.5209152}},{5,{-0.8031747,1.3708649,-0.56768906}},{12,{-1.520915,-0.5676891,2.088603}}}},
		{'C',{{0,{2.35021,-0.7880346,-1.5621773}},{5,{-0.78803474,1.3541472,-0.5661155}},{12,{-1.5621774,-0.56611556,2.1282897}}}},
		{'G',{{0,{2.3250363,-0.80280226,-1.5222337}},{5,{-0.8028024,1.3702725,-0.5674678}},{12,{-1.5222336,-0.5674678,2.0897026}}}},
		{'U',{{0,{2.3493478,-0.7881313,-1.5612125}},{5,{-0.7881318,1.3542684,-0.5661326}},{12,{-1.5612122,-0.5661328,2.1273444}}}}
	};
	
	//map<char,map<int,vector<double>>> basecoord = {
	//{'A',{{0,{3.063,8.025,-4.135}},{5,{6.719,6.889,-3.258}},{12,{5.666,4.158,-1.647}}}},
	//{'C',{{0,{6.913,5.099,-6.683}},{5,{9.376,2.167,-5.806}},{12,{7.014,0.438,-4.195}}}},
	//{'G',{{0,{-1.758,-8.408,1.587}},{5,{1.932,-9.427,0.710}},{12,{2.522,-6.560,-0.901}}}},
	//{'U',{{0,{3.063,-8.025,4.135}},{5,{6.719,-6.889,3.258}},{12,{5.666,-4.158,1.647}}}}
	//};

	if(!flagacc)
	{
		for(i=0;i<numseq;i++)
		{
			basec=basecoord[decstr[i].aaa];
			count=0;
			for(auto it = basec.begin(); it != basec.end(); ++it){
				aidx = it->first;
				bcoor = it->second;

				xyz_list1[count][0]=decstr[i].ptv[aidx].x;
				xyz_list1[count][1]=decstr[i].ptv[aidx].y;
				xyz_list1[count][2]=decstr[i].ptv[aidx].z;
				
				xyz_list2[count][0]=bcoor[0];
				xyz_list2[count][1]=bcoor[1];
				xyz_list2[count][2]=bcoor[2];
				count+=1;
			}
			bf.RotateCoor(xyz_list2,xyz_list1,RotMatix,TranVect);
			for(j=0;j<numseq;j++)
			{
				k=j*numseq+i;
				basec=basecoord[decstr[i].aaa];
				faperpp.strpval[k].resize(basec.size());
				
				count=0;
				for(auto it = basec.begin(); it != basec.end(); ++it)
				{
					aidx = it->first;
					xyz_tmp[0]=decstr[j].ptv[aidx].x;
					xyz_tmp[1]=decstr[j].ptv[aidx].y;
					xyz_tmp[2]=decstr[j].ptv[aidx].z;
					//cout<<to_string(xyz_tmp[0])+"|"+to_string(xyz_tmp[1])+"|"+to_string(xyz_tmp[2])<<endl;
					bf.Transform(xyz_tmp,RotMatix,TranVect,faperpp.strpval[k][count]);
					count+=1;
					//cout<<to_string(strpval[k][count].x)+"|"+to_string(strpval[k][count].y)+"|"+to_string(strpval[k][count].z)<<endl;
				}
			}
		}
	}else{
		for(i=0;i<numseq;i++)
		{
			basec=basecoord[decstr[i].aaa];
			count=0;
			for(auto it = basec.begin(); it != basec.end(); ++it)
			{
				aidx = it->first;
				bcoor = it->second;

				xyz_list1[count][0]=decstr[i].ptv[aidx].x;
				xyz_list1[count][1]=decstr[i].ptv[aidx].y;
				xyz_list1[count][2]=decstr[i].ptv[aidx].z;
				
				xyz_list2[count][0]=bcoor[0];
				xyz_list2[count][1]=bcoor[1];
				xyz_list2[count][2]=bcoor[2];
				count+=1;
			}
			bf.RotateCoor(xyz_list2,xyz_list1,RotMatix,TranVect);
			for(j=0;j<numseq;j++)
			{
				k=j*numseq+i;
				if(!faperpp.flap[k]) continue;

				basec=basecoord[decstr[i].aaa];
				faperpp.strpval[k].resize(basec.size());
				
				count=0;
				for(auto it = basec.begin(); it != basec.end(); ++it)
				{
					aidx = it->first;
					xyz_tmp[0]=decstr[j].ptv[aidx].x;
					xyz_tmp[1]=decstr[j].ptv[aidx].y;
					xyz_tmp[2]=decstr[j].ptv[aidx].z;
					//cout<<to_string(xyz_tmp[0])+"|"+to_string(xyz_tmp[1])+"|"+to_string(xyz_tmp[2])<<endl;
					bf.Transform(xyz_tmp,RotMatix,TranVect,faperpp.strpval[k][count]);
					count+=1;
					//cout<<to_string(strpval[k][count].x)+"|"+to_string(strpval[k][count].y)+"|"+to_string(strpval[k][count].z)<<endl;
				}
				faperpp.flap[k]=false;
			}
		}
	}
}

void RNASimulation::distsingleref(vector<Atom> proseq,vector<Boneinfo> bb, int numbb,vector<Boneinfo> ibb, int inumbb,int numseq)
{

	initdistnpaa=0;
	int i,j,k;

	int aidx;
	vector<double> bcoor;
	vector<double> tmp(FAPE_BASEATOM_NUM,0);
    vector<vector<double> > xyz_list1(FAPE_BASEATOM_NUM,tmp);
    vector<vector<double> > xyz_list2(FAPE_BASEATOM_NUM,tmp);
	vector<vector<double>> xyz_list3(FAPE_BASEATOM_NUM,tmp);
	vector<double> xyz_tmp(tmp);
	vector<vector<double> > RotMatix;
    vector<double> TranVect;
	map<int,vector<double>> basec;
	Point3d pcd;
	//lentmp.resize(refidx+1);
	//ntflag.resize(refidx+1);
	//ntstr.resize(refidx+1);
	
	//resind in simu [0,numseq-1]
	//ibb           resind don't need from 1
	//bb        similar range as ibb
	int bbmin=1000000;
	int bbmax=-10000000;

	//残基idの最大値(bbmax)と最小値(bbmin)を求める
	for(i=0;i<inumbb;i++)
	{
		if(bbmin>ibb[i].resind) bbmin=ibb[i].resind;
		if(bbmax<ibb[i].resind) bbmax=ibb[i].resind;
	}

	//残基数+5000の配列を用意し、0~残基数+1の要素にはiを格納
	//int *bbind=new int[5000+bbmax-bbmin+1];
	vector<int> bbind(5000+bbmax-bbmin+1);
	for(i=0;i<bbmax-bbmin+1;i++) bbind[i]=-1; 
	for(i=0;i<inumbb;i++) bbind[ibb[i].resind-bbmin]=i;

	initdist.resize(initdist.size()+numseq*numseq);
	
	k=initdistnpaa;
	for(i=0;i<numbb;i++)
	{
		if(bb[i].indatom[3]==-1) continue;
		if(bb[i].resind<bbmin || bb[i].resind>bbmax) continue;
		if(bbind[bb[i].resind-bbmin]==-1) continue;

		for(j=i+2;j<numbb;j++)
		{
			if(bb[i].indatom[3]==-1) continue;
			if(bb[j].resind<bbmin || bb[j].resind>bbmax) continue;
			if(bbind[bb[j].resind-bbmin]==-1) continue;
			if(fabs(bbind[bb[j].resind-bbmin]-bbind[bb[i].resind-bbmin])<1) continue;
			
			initdist[k].ispair=true;
			initdist[k].pbin=bbind[bb[i].resind-bbmin];
			initdist[k].tnum=bbind[bb[j].resind-bbmin];

			pcd.x=proseq[bb[j].indatom[3]].x-proseq[bb[i].indatom[3]].x; 
			pcd.y=proseq[bb[j].indatom[3]].y-proseq[bb[i].indatom[3]].y; 
			pcd.z=proseq[bb[j].indatom[3]].z-proseq[bb[i].indatom[3]].z;
			initdist[k].dstd=initdist[k].dist;
			if(initdist[k].dstd<3.9) initdist[k].dstd=3.9;
			initdist[k].dist=bf.norm(pcd);
			k++;
		}

	}
	initdistnpaa=k;
	k=initdistnpaa;
	if(k<numseq*numseq) k=numseq*numseq;
	initdist.resize(k);

	//refinestructureとreference structureで共通したの残基数をlentmpに格納
	initlentmp=0; 
	for(i=0;i<numbb;i++)
	{
		if(bb[i].resind<bbmin || bb[i].resind>bbmax) continue;
		if(bbind[bb[i].resind-bbmin]==-1) continue;
		initlentmp++;
	}
	//refinestructureとreference structureで共通したの残基であるかどうか(共通:true、共通でない:false)をntflagに格納
	//refinestructureとreference structureで共通したの残基のC4'の座標をntstrに格納
	initntflag.resize(inumbb);
	initntstr.resize(inumbb*FAPE_BASEATOM_NUM);

	int maxbbind=0;
	j=0;
	for(i=0;i<numbb;i++)
	{

		if(bb[i].resind<bbmin || bb[i].resind>bbmax) continue;
		if(bbind[bb[i].resind-bbmin]==-1) continue;

	//	ntflag[bb[i].resind-1]=true;
		if(maxbbind<bbind[bb[i].resind-bbmin]) maxbbind=bbind[bb[i].resind-bbmin];
		initntflag[bbind[bb[i].resind-bbmin]]=true;
		initntstr[j]=proseq[bb[i].indatom[3]].x;
		initntstr[initlentmp+j]=proseq[bb[i].indatom[3]].y;
		initntstr[2*initlentmp+j]=proseq[bb[i].indatom[3]].z;
		j++;
	}
	initntstr.resize(FAPE_BASEATOM_NUM*initlentmp);
	initntflag.resize(maxbbind);

	bbind.clear();
	bbind.shrink_to_fit();
}

double RNASimulation::fenergyinitcontact(vector<Residueatoms> decstr,int numseq)
{
	if(initdistnpaa==0) return 0.0;
	double totene=0;
	int i,j,k;
	Point3d tp;
	double tval;
	for(i=0;i<numseq;i++)
	{
		for(j=i+2;j<numseq;j++)
		{
			k=j*numseq+i;
			tp.x=decstr[j].ptv[5].x-decstr[i].ptv[5].x;
			tp.y=decstr[j].ptv[5].y-decstr[i].ptv[5].y;
			tp.z=decstr[j].ptv[5].z-decstr[i].ptv[5].z;
			initdist[k].pval=bf.norm(tp);
		}
	}
	for(i=0;i<initdistnpaa;i++)
	{
		k=initdist[i].tnum*numseq+initdist[i].pbin;
		tval=fabs(initdist[k].pval-initdist[i].dist)/initdist[i].dstd;//deltad
		totene+=tval;
	}
	return totene/double(initdistnpaa);
}








