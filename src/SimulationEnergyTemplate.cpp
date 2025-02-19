// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

double calcFAPEthresh(double dist){
	return 5+(130/(5+exp(0.5*(dist-30))));
}

void RNASimulation::distsingleref(vector<Atom> proseq,vector<Boneinfo> bb, int numbb,vector<Boneinfo> ibb, int inumbb,int numseq)
{
	initnpaa3=0;
	int i,j,k,n;
	int aidx;
	Point3d pcd;
	int bbmin=1000000;
	int bbmax=-10000000;
	int count;

	//残基idの最大値(bbmax)と最小値(bbmin)を求める
	for(i=0;i<inumbb;i++)
	{
		if(bbmin>ibb[i].resind) bbmin=ibb[i].resind;
		if(bbmax<ibb[i].resind) bbmax=ibb[i].resind;
	}

	//残基数+5000の配列を用意し、0~残基数+1の要素にはiを格納
	vector<int> bbind(5000+bbmax-bbmin+1);
	for(i=0;i<bbmax-bbmin+1;i++) bbind[i]=-1; 
	for(i=0;i<inumbb;i++) bbind[ibb[i].resind-bbmin]=i;

	initpaa3.resize(initpaa3.size()+numseq*numseq);
	
	k=initnpaa3;
	for(i=0;i<numbb;i++)
	{
		if(bb[i].indatom[0]==-1 || bb[i].indatom[3]==-1) continue;
		if(bb[i].resind<bbmin || bb[i].resind>bbmax) continue;
		if(bbind[bb[i].resind-bbmin]==-1) continue;

		for(j=i+1;j<numbb;j++)
		{
			if(bb[j].indatom[0]==-1 || bb[j].indatom[3]==-1) continue;
			if(bb[j].resind<bbmin || bb[j].resind>bbmax) continue;
			if(bbind[bb[j].resind-bbmin]==-1) continue;
			if(fabs(bbind[bb[j].resind-bbmin]-bbind[bb[i].resind-bbmin])<1) continue;

			initpaa3[k].dist.resize(DISTTARGETATOMS.size());
			initpaa3[k].threshold.resize(DISTTARGETATOMS.size());
			count=0;
			for(n=0;n<DISTTARGETATOMS.size();n++)
			{
				aidx = DISTTARGETATOMS[n];
				pcd.x=proseq[bb[i].indatom[aidx]].x-proseq[bb[j].indatom[aidx]].x; 
				pcd.y=proseq[bb[i].indatom[aidx]].y-proseq[bb[j].indatom[aidx]].y; 
				pcd.z=proseq[bb[i].indatom[aidx]].z-proseq[bb[j].indatom[aidx]].z;
				initpaa3[k].dist[n] = bf.norm(pcd);
				initpaa3[k].threshold[n]=calcFAPEthresh(initpaa3[k].dist[n]);
				count+=1;
			}

		//	if(bb[j].resind-bb[i].resind<2) continue;
			
			initpaa3[k].ispair=true;
			initpaa3[k].pbin=bbind[bb[i].resind-bbmin];
			initpaa3[k].tnum=bbind[bb[j].resind-bbmin];
			k++;
		}
	}
	initnpaa3=k;
	if(k<numseq*numseq) k=numseq*numseq;
	initpaa3.resize(k);

	bbind.clear();
	bbind.shrink_to_fit();

}

//getinitdist:: reference structureとrefined structurenの共通の残基を特定し、共通した残基についてreference structureの残機間の距離を算出
//refrence配列のatomの情報、bb:reference配列のbackbornの情報、numbb:reference配列のbackbornの残基数、refined配列のbackborn atomの情報、refined配列のbackbornの残基数
void RNASimulation::distmultipleref(vector<Atom> proseq,vector<Boneinfo> bb, int numbb,vector<Boneinfo> ibb, int inumbb,int numseq, double weight, int refidx)
{
	npaa3.resize(refidx+1);
	npaa3[refidx]=0;
	int i,j,k,n;
	int aidx;
	Point3d pcd;
	int bbmin=1000000;
	int bbmax=-10000000;
	int count;

	//残基idの最大値(bbmax)と最小値(bbmin)を求める
	for(i=0;i<inumbb;i++)
	{
		if(bbmin>ibb[i].resind) bbmin=ibb[i].resind;
		if(bbmax<ibb[i].resind) bbmax=ibb[i].resind;
	}

	//残基数+5000の配列を用意し、0~残基数+1の要素にはiを格納
	vector<int> bbind(5000+bbmax-bbmin+1);
	for(i=0;i<bbmax-bbmin+1;i++) bbind[i]=-1; 
	for(i=0;i<inumbb;i++) bbind[ibb[i].resind-bbmin]=i;

	if(paa3.size() < refidx + 1) paa3.resize(refidx + 1);
	paa3[refidx].resize(npaa3[refidx]+numseq*numseq);
	
	k=npaa3[refidx];
	for(i=0;i<numbb;i++)
	{
		if(bb[i].indatom[0]==-1 || bb[i].indatom[3]==-1) continue;
		if(bb[i].resind<bbmin || bb[i].resind>bbmax) continue;
		if(bbind[bb[i].resind-bbmin]==-1) continue;

		for(j=i+1;j<numbb;j++)
		{
			if(bb[j].indatom[0]==-1 || bb[j].indatom[3]==-1) continue;
			if(bb[j].resind<bbmin || bb[j].resind>bbmax) continue;
			if(bbind[bb[j].resind-bbmin]==-1) continue;
			if(fabs(bbind[bb[j].resind-bbmin]-bbind[bb[i].resind-bbmin])<1) continue;

			paa3[refidx][k].dist.resize(DISTTARGETATOMS.size());
			paa3[refidx][k].threshold.resize(DISTTARGETATOMS.size());
			count=0;
			for(n=0;n<DISTTARGETATOMS.size();n++)
			{
				aidx = DISTTARGETATOMS[n];
				pcd.x=proseq[bb[i].indatom[aidx]].x-proseq[bb[j].indatom[aidx]].x; 
				pcd.y=proseq[bb[i].indatom[aidx]].y-proseq[bb[j].indatom[aidx]].y; 
				pcd.z=proseq[bb[i].indatom[aidx]].z-proseq[bb[j].indatom[aidx]].z;
				paa3[refidx][k].dist[n] = bf.norm(pcd);
				paa3[refidx][k].threshold[n] = calcFAPEthresh(paa3[refidx][k].dist[n]);
				count+=1;
			}
			
			paa3[refidx][k].ispair=true;
			paa3[refidx][k].pbin=bbind[bb[i].resind-bbmin];
			paa3[refidx][k].tnum=bbind[bb[j].resind-bbmin];
			paa3[refidx][k].weight=weight;
			k++;
		}

	}
	npaa3[refidx]=k;
	if(k<numseq*numseq) k=numseq*numseq;
	paa3[refidx].resize(k);

	bbind.clear();
	bbind.shrink_to_fit();

}

void RNASimulation::distenergysingle(vector<Residueatoms> decstr,int numseq,double &fene,Energyterms &lp, vector<vector<vector<double>>> distmaprpp, vector<double> &distenerpp, vector<vector<bool>> mapbbflagrpp,bool flagacc)
{
	double totene=0;
	int i,j,k,l,n;
	Point3d tp;
	double tval;
	int aidx;
	double indf=0;
	double dist;

	if(!flagacc)
	{
		lp.fn[DIST_ENERGY_ST]=0;

		for(i=0;i<initnpaa3;i++)
		{
			k=initpaa3[i].tnum*numseq+initpaa3[i].pbin;
			distenerpp[k]=0;
			for(n=0;n<DISTTARGETATOMS.size();n++)
			{	
				dist=initpaa3[i].dist[n]-distmaprpp[initpaa3[i].pbin][initpaa3[i].tnum][n];
				tval=sqrt(dist*dist);

				//if(tval > DIST_THRESH) tval=DIST_THRESH;
				if(tval > initpaa3[i].threshold[n]) tval=initpaa3[i].threshold[n];
				distenerpp[k]+=tval;
			}
			indf+=distenerpp[k];
		}
		lp.fn[DIST_ENERGY_ST]=indf;
	}
	else{
		double eold=0; double enew=0;

		for(i=0;i<initnpaa3;i++)
		{
			k=initpaa3[i].tnum*numseq+initpaa3[i].pbin;

			if(mapbbflagrpp[initpaa3[i].pbin][initpaa3[i].tnum])
			{
				eold+=distenerpp[k];
				distenerpp[k]=0;
				for(n=0;n<DISTTARGETATOMS.size();n++)
				{	
					dist=initpaa3[i].dist[n]-distmaprpp[initpaa3[i].pbin][initpaa3[i].tnum][n];
					tval=sqrt(dist*dist);

					//if(tval > DIST_THRESH) tval=DIST_THRESH;
					if(tval > initpaa3[i].threshold[n]) tval=initpaa3[i].threshold[n];
					distenerpp[k]+=tval;
				}
				enew+=distenerpp[k];
			}
		}
		lp.fn[DIST_ENERGY_ST]+=(enew-eold);
	}
	fene=lp.fn[DIST_ENERGY_ST];
}


void RNASimulation::distenergymultiple(vector<Residueatoms> decstr,int numseq,double &fene,Energyterms &lp, vector<vector<vector<double>>> distmaprpp, vector<vector<double>> &multidistenerpp, vector<vector<bool>> mapbbflagrpp, bool flagacc)
{
	double totene=0;
	int i,j,k,l,n;
	Point3d tp;
	double tval;
	int aidx;
	double indf=0;
	double dist;

	if(!flagacc)
	{
		lp.fn[DIST_ENERGY_ST+1]=0;
		for(l=0;l<npaa3.size();l++)
		{
			for(i=0;i<npaa3[l];i++)
			{
				k=paa3[l][i].tnum*numseq+paa3[l][i].pbin;
				multidistenerpp[l][k]=0;
				for(n=0;n<DISTTARGETATOMS.size();n++)
				{	
					dist=paa3[l][i].dist[n]-distmaprpp[paa3[l][i].pbin][paa3[l][i].tnum][n];
					tval=sqrt(dist*dist);

					//if(tval > DIST_THRESH) tval=DIST_THRESH;
					if(tval > paa3[l][i].threshold[n]) tval=paa3[l][i].threshold[n];
					multidistenerpp[l][k]+=paa3[l][k].weight*tval;
				}
				indf+=multidistenerpp[l][k];
				//cout<<to_string(l) + "-" + to_string(i) + "-" + to_string(paa3[l][i].dist[0]) + "-" + to_string(paa3[l][i].dist[1]) + "-" + to_string(distmaprpp[paa3[l][i].pbin][paa3[l][i].tnum][0]) + "-" + to_string(distmaprpp[paa3[l][i].pbin][paa3[l][i].tnum][1]) + "-" + to_string(multidistenerpp[l][k])<<endl;
			}
		}
		lp.fn[DIST_ENERGY_ST+1]=indf;
	}else{
		double eold=0; double enew=0;

		for(l=0;l<npaa3.size();l++)
		{
			for(i=0;i<npaa3[l];i++)
			{
				k=paa3[l][i].tnum*numseq+paa3[l][i].pbin;

				if(mapbbflagrpp[initpaa3[i].pbin][initpaa3[i].tnum])
				{
					eold+=multidistenerpp[l][k];
					multidistenerpp[l][k]=0;
					for(n=0;n<DISTTARGETATOMS.size();n++)
					{
						dist=paa3[l][i].dist[n]-distmaprpp[paa3[l][i].pbin][paa3[l][i].tnum][n];
						tval=sqrt(dist*dist);

						//if(tval > DIST_THRESH) tval=DIST_THRESH;
						if(tval > paa3[l][i].threshold[n]) tval=paa3[l][i].threshold[n];
						multidistenerpp[l][k]+=paa3[l][k].weight*tval;
					}
					enew+=multidistenerpp[l][k];
					//cout<<to_string(l) + "-" + to_string(i) + "-" + to_string(paa3[l][i].dist[0]) + "-" + to_string(paa3[l][i].dist[1]) + "-" + to_string(distmaprpp[paa3[l][i].pbin][paa3[l][i].tnum][0]) + "-" + to_string(distmaprpp[paa3[l][i].pbin][paa3[l][i].tnum][1]) + "-" + to_string(multidistenerpp[l][k])<<endl;
				}//else{
					//cout<<to_string(l) + "-" + to_string(i) + "-" + to_string(paa3[l][i].dist[0]) + "-" + to_string(paa3[l][i].dist[1]) + "-" + to_string(distmaprpp[paa3[l][i].pbin][paa3[l][i].tnum][0]) + "-" + to_string(distmaprpp[paa3[l][i].pbin][paa3[l][i].tnum][1]) + "-" + to_string(multidistenerpp[l][k])<<endl;
				//}
			}
		}

		lp.fn[DIST_ENERGY_ST+1]+=(enew-eold);
	}
	fene=lp.fn[DIST_ENERGY_ST+1];
}








