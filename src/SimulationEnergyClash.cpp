// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

double RNASimulation::fenergyoxhydb(vector<Residueatoms> decstr, int basei, int basej, int ii, int jj, double tdist){
	double ene=0;
	Point3s tmp2, tmp3, tmp4;
	float ang;
	char aat;
	int binidx1, binidx2, binnum2;
	//if(tdist < upper & lower < tdist) ene = -exp(-0.5 * std::pow((tdist - mean) / stddev, 2));
	//else ene = 0;
	if(ii == 9){
		aat = decstr[basej].aaa;

		if(tdist > OXHYDBPARAM[aat][AtomOxhydbId2pId[aat][jj]].minval2 & tdist < OXHYDBPARAM[aat][AtomOxhydbId2pId[aat][jj]].maxval2){
			tmp2.x = decstr[basei].ptv[6].x; tmp2.y = decstr[basei].ptv[6].y; tmp2.z = decstr[basei].ptv[6].z;
			tmp3.x = decstr[basei].ptv[9].x; tmp3.y = decstr[basei].ptv[9].y; tmp3.z = decstr[basei].ptv[9].z;
			tmp4.x = decstr[basej].ptv[jj].x; tmp4.y = decstr[basej].ptv[jj].y; tmp4.z = decstr[basej].ptv[jj].z;
			coord2ang(tmp2, tmp3, tmp4, ang);

			binidx1 = getbin(ang, OXHYDBPARAM[aat][AtomOxhydbId2pId[aat][jj]].boundary1);
			binidx2 = getbin(tdist, OXHYDBPARAM[aat][AtomOxhydbId2pId[aat][jj]].boundary2);
			binnum2=OXHYDBPARAM[aat][AtomOxhydbId2pId[aat][jj]].binnum2;
			ene-=OXHYDBPARAM[aat][AtomOxhydbId2pId[aat][jj]].energy[binidx1*binnum2+binidx2];

			cout<< ii << "-" << jj << "-" << ang << "-" << tdist << "-" << OXHYDBPARAM[aat][AtomOxhydbId2pId[aat][jj]].energy[binidx1*binnum2+binidx2]<<endl;
		}
	}
	if(jj == 9){
		aat = decstr[basei].aaa;
		if(tdist > OXHYDBPARAM[aat][AtomOxhydbId2pId[aat][ii]].minval2 & tdist < OXHYDBPARAM[aat][AtomOxhydbId2pId[aat][ii]].maxval2){
			tmp2.x = decstr[basej].ptv[6].x; tmp2.y = decstr[basej].ptv[6].y; tmp2.z = decstr[basej].ptv[6].z;
			tmp3.x = decstr[basej].ptv[9].x; tmp3.y = decstr[basej].ptv[9].y; tmp3.z = decstr[basej].ptv[9].z;
			tmp4.x = decstr[basei].ptv[ii].x; tmp4.y = decstr[basei].ptv[ii].y; tmp4.z = decstr[basei].ptv[ii].z;
			coord2ang(tmp2, tmp3, tmp4, ang);

			binidx1 = getbin(ang, OXHYDBPARAM[aat][AtomOxhydbId2pId[aat][ii]].boundary1);
			binidx2 = getbin(tdist, OXHYDBPARAM[aat][AtomOxhydbId2pId[aat][ii]].boundary2);
			binnum2=OXHYDBPARAM[aat][AtomOxhydbId2pId[aat][ii]].binnum2;
			ene-=OXHYDBPARAM[aat][AtomOxhydbId2pId[aat][ii]].energy[binidx1*binnum2+binidx2];

			cout<< ii << "-" << jj << "-" << ang << "-" << tdist << "-" << OXHYDBPARAM[aat][AtomOxhydbId2pId[aat][ii]].energy[binidx1*binnum2+binidx2]<<endl;
		}
	}

	return ene;
}

void RNASimulation::fenergyclash(vector<Residueatoms> decstr,int numseq,double &nene,double &fene, double &hfene, Energyterms &lp,vector<vector<Rescontact>> &rpp, vector<vector<int>> &nmaprpp, vector<vector<vector<double>>> &distmaprpp, vector<vector<vector<double>>> &ntpairdistrpp, vector<vector<vector<double>>> &adjpairdistrpp, vector<vector<bool>> mapflagrpp, vector<vector<bool>> mapbbflagrpp,bool flagacc){
	int i,j,k;
	int ii,jj;
	int atomidx1, atomidx2;
	int inda,indb;
	char aat1, aat2, aatk;
	bool flagmi,flagmj,flagres;
	int atom1, atom2, atom3, atom4;
	double cdihe, cleng, leng[4];
	int ijtype;
	int basetype;
	//vector<vector<int>> baseparams;
	Point3d tp;
	Point3d ve;
	Point3s tmp1,tmp2,tmp3,tmp4;
	BasicFunc bf;
	double tdist,dspair1,dspair2,penal;
	double dscut=1.0;
	double minicut=0.01;
	int bptype;
	double fn;
	int nn;
	double tmpval;
	bool flagatom;
	double hbondene;

	if(!flagacc)
	{
		for(i=0;i<CLASH_ENERGY_NUM;i++)
		{
			lp.nn[CLASH_ENERGY_ST+i]=0;
			lp.fn[CLASH_ENERGY_ST+i]=0;
		}

		for(i=0;i<numseq;i++)
		{
			for(j=i;j<numseq;j++)
			//for(j=i+1;j<numseq;j++)
			{
				for(k=0;k<CLASH_ENERGY_NUM;k++)
				{
					rpp[i][j].pval[k]=0;
					rpp[i][j].pcw[k]=0;
				}
				rpp[i][j].pval[1]=0;
				nmaprpp[i][j]=0; 
				nmaprpp[j][i]=0;
			}
		}

		for(i=0;i<numseq;i++)
		{
			inda=decstr[i].resind;
			aat1=decstr[i].aaa;

			for(ii=0;ii<NNheavy[inda][0];ii++)
			{
				atomidx1=NNheavy[inda][ii+1];
				for(j=i;j<numseq;j++)
				{
					indb=decstr[j].resind;
					aat2=decstr[j].aaa;

					for(jj=0;jj<NNheavy[indb][0];jj++)
					{
						if(i == j & ii >= jj) continue;

						atomidx2=NNheavy[indb][jj+1];
						if(atomidx1 == 32 | atomidx1 == 34 | atomidx1 == 35 | atomidx1 == 36 | atomidx1 == 37 | atomidx1 == 38 | atomidx1 == 39 | atomidx1 == 40 | atomidx1 == 41 | atomidx1 == 42  | atomidx1 == 43  | atomidx1 == 44  | atomidx1 == 45) continue;
						if(atomidx2 == 32 | atomidx2 == 34 | atomidx2 == 35 | atomidx2 == 36 | atomidx2 == 37 | atomidx2 == 38 | atomidx2 == 39 | atomidx2 == 40 | atomidx2 == 41 | atomidx2 == 42  | atomidx2 == 43  | atomidx2 == 44  | atomidx2 == 45) continue;

						auto itcontactatom = std::find(DISTTARGETATOMS.begin(), DISTTARGETATOMS.end(), atomidx1);
						if((atomidx1 == atomidx2) && (itcontactatom != DISTTARGETATOMS.end())){
							tp.x=decstr[i].ptv[ii].x-decstr[j].ptv[jj].x;
							tp.y=decstr[i].ptv[ii].y-decstr[j].ptv[jj].y;
							tp.z=decstr[i].ptv[ii].z-decstr[j].ptv[jj].z;
							tdist=sqrt(tp.x*tp.x+tp.y*tp.y+tp.z*tp.z);

							int index = std::distance(DISTTARGETATOMS.begin(), itcontactatom);
							
							distmaprpp[i][j][index]=tdist;
							distmaprpp[j][i][index]=tdist;

							if(i==j){
								flagatom = EXPCLASHATOMSNT[aat1][ii][jj];
								if (flagatom) continue;
							}

							if(i+1==j){
								flagatom = EXPCLASHATOMSADJ[aat1][aat2][ii][jj];
								if (flagatom) continue;
							}

							if(fabs(tp.x)>15) continue;
							if(fabs(tp.y)>15) continue;
							if(fabs(tp.z)>15) continue;
						}
						/*
						else if((ii == 9 & Atomoxhydbtarget[aat2][jj] == 1) | (jj == 9 & Atomoxhydbtarget[aat1][ii] == 1)){
							tp.x=decstr[i].ptv[ii].x-decstr[j].ptv[jj].x;
							tp.y=decstr[i].ptv[ii].y-decstr[j].ptv[jj].y;
							tp.z=decstr[i].ptv[ii].z-decstr[j].ptv[jj].z;
							tdist=sqrt(tp.x*tp.x+tp.y*tp.y+tp.z*tp.z);

							hbondene=fenergyoxhydb(decstr, i, j, ii, jj, tdist);

							lp.fn[OXHBOND_ENERGY_ST]+=hbondene;
							rpp[i][j].pval[1]+=hbondene;

							if(i==j){
								flagatom = EXPCLASHATOMSNT[aat1][ii][jj];
								if (flagatom) continue;
							}

							if(i+1==j){
								flagatom = EXPCLASHATOMSADJ[aat1][aat2][ii][jj];
								if (flagatom) continue;
							}

							if(fabs(tp.x)>15) continue;
							if(fabs(tp.y)>15) continue;
							if(fabs(tp.z)>15) continue;
						}
						*/
						else{
							if(i==j){
								flagatom = EXPCLASHATOMSNT[aat1][ii][jj];
								if (flagatom) continue;
							}

							if(i+1==j){
								flagatom = EXPCLASHATOMSADJ[aat1][aat2][ii][jj];
								if (flagatom) continue;
							}

							tp.x=decstr[i].ptv[ii].x-decstr[j].ptv[jj].x;
							if(fabs(tp.x)>15) continue;
							tp.y=decstr[i].ptv[ii].y-decstr[j].ptv[jj].y;
							if(fabs(tp.y)>15) continue;
							tp.z=decstr[i].ptv[ii].z-decstr[j].ptv[jj].z;
							if(fabs(tp.z)>15) continue;

							tdist=sqrt(tp.x*tp.x+tp.y*tp.y+tp.z*tp.z);
						}

						dspair2=VWPARAM[atomidx1]+VWPARAM[atomidx2]-0.4;

						if(tdist<dspair2){
							penal=dspair2-tdist;
							penal=sqrt(penal*penal);

							lp.nn[CLASH_ENERGY_ST]+=1;
							lp.fn[CLASH_ENERGY_ST]+=penal;
							rpp[i][j].pcw[0]+=1;
							rpp[i][j].pval[0]+=penal;
							nmaprpp[i][j]=1; nmaprpp[j][i]=1;

						}
					}
				}
			}
		}
	}
	else
	{
		int iold[CLASH_ENERGY_NUM+OXHBOND_ENERGY_NUM],inew[CLASH_ENERGY_NUM+OXHBOND_ENERGY_NUM];
		double fold[CLASH_ENERGY_NUM+OXHBOND_ENERGY_NUM],fnew[CLASH_ENERGY_NUM+OXHBOND_ENERGY_NUM];
		for(i=0;i<CLASH_ENERGY_NUM+OXHBOND_ENERGY_NUM;i++){
			iold[i]=0;inew[i]=0;fold[i]=0;fnew[i]=0;
		}

		for(i=0;i<numseq;i++)
		{
			for(j=i;j<numseq;j++)
			{
				if(mapflagrpp[i][j]){
					iold[0]+=rpp[i][j].pcw[0];//cla
					fold[0]+=rpp[i][j].pval[0];//cla
					rpp[i][j].pcw[0]=0;
					rpp[i][j].pval[0]=0;

					fold[1]+=rpp[i][j].pval[1];//cla
					rpp[i][j].pval[1]=0;
				}
				nmaprpp[i][j]=0; 
				nmaprpp[j][i]=0;
			}
		}

		for(i=0;i<numseq;i++)
		{
			inda=decstr[i].resind;
			aat1=decstr[i].aaa;

			for(ii=0;ii<NNheavy[inda][0];ii++)
			{
				atomidx1=NNheavy[inda][ii+1];
				
				for(j=i;j<numseq;j++)
				{
					indb=decstr[j].resind;
					aat2=decstr[j].aaa;

					for(jj=0;jj<NNheavy[indb][0];jj++)
					{
						if(i == j & ii >= jj) continue;

						atomidx2=NNheavy[indb][jj+1];
						if(atomidx1 == 32 | atomidx1 == 34 | atomidx1 == 35 | atomidx1 == 36 | atomidx1 == 37 | atomidx1 == 38 | atomidx1 == 39 | atomidx1 == 40 | atomidx1 == 41 | atomidx1 == 42  | atomidx1 == 43  | atomidx1 == 44  | atomidx1 == 45) continue;
						if(atomidx2 == 32 | atomidx2 == 34 | atomidx2 == 35 | atomidx2 == 36 | atomidx2 == 37 | atomidx2 == 38 | atomidx2 == 39 | atomidx2 == 40 | atomidx2 == 41 | atomidx2 == 42  | atomidx2 == 43  | atomidx2 == 44  | atomidx2 == 45) continue;
						if(!mapflagrpp[i][j]) continue;

						auto itcontactatom = std::find(DISTTARGETATOMS.begin(), DISTTARGETATOMS.end(), atomidx1);
						if((atomidx1 == atomidx2) && (itcontactatom != DISTTARGETATOMS.end())){
							tp.x=decstr[i].ptv[ii].x-decstr[j].ptv[jj].x;
							tp.y=decstr[i].ptv[ii].y-decstr[j].ptv[jj].y;
							tp.z=decstr[i].ptv[ii].z-decstr[j].ptv[jj].z;
							tdist=sqrt(tp.x*tp.x+tp.y*tp.y+tp.z*tp.z);

							int index = std::distance(DISTTARGETATOMS.begin(), itcontactatom);

							distmaprpp[i][j][index]=tdist;
							distmaprpp[j][i][index]=tdist;

							if(i==j){
								flagatom = EXPCLASHATOMSNT[aat1][ii][jj];
								if (flagatom) continue;
							}

							if(i+1==j){
								flagatom = EXPCLASHATOMSADJ[aat1][aat2][ii][jj];
								if (flagatom) continue;
							}
							if(fabs(tp.x)>15) continue;
							if(fabs(tp.y)>15) continue;
							if(fabs(tp.z)>15) continue;
						}
						/*
						else if((ii == 9 & Atomoxhydbtarget[aat2][jj] == 1) | (jj == 9 & Atomoxhydbtarget[aat1][ii] == 1)){
							tp.x=decstr[i].ptv[ii].x-decstr[j].ptv[jj].x;
							tp.y=decstr[i].ptv[ii].y-decstr[j].ptv[jj].y;
							tp.z=decstr[i].ptv[ii].z-decstr[j].ptv[jj].z;
							tdist=sqrt(tp.x*tp.x+tp.y*tp.y+tp.z*tp.z);

							hbondene=fenergyoxhydb(decstr, i, j, ii, jj, tdist);
							rpp[i][j].pval[1]+=hbondene;

							if(i==j){
								flagatom = EXPCLASHATOMSNT[aat1][ii][jj];
								if (flagatom) continue;
							}

							if(i+1==j){
								flagatom = EXPCLASHATOMSADJ[aat1][aat2][ii][jj];
								if (flagatom) continue;
							}

							if(fabs(tp.x)>15) continue;
							if(fabs(tp.y)>15) continue;
							if(fabs(tp.z)>15) continue;
						}
						*/
						else{
							if(i==j){
								flagatom = EXPCLASHATOMSNT[aat1][ii][jj];
								if (flagatom) continue;
							}

							if(i+1==j){
								flagatom = EXPCLASHATOMSADJ[aat1][aat2][ii][jj];
								if (flagatom) continue;
							}
							tp.x=decstr[i].ptv[ii].x-decstr[j].ptv[jj].x;
							if(fabs(tp.x)>15) continue;
							tp.y=decstr[i].ptv[ii].y-decstr[j].ptv[jj].y;
							if(fabs(tp.y)>15) continue;
							tp.z=decstr[i].ptv[ii].z-decstr[j].ptv[jj].z;
							if(fabs(tp.z)>15) continue;

							tdist=sqrt(tp.x*tp.x+tp.y*tp.y+tp.z*tp.z);
						}
						dspair2=VWPARAM[atomidx1]+VWPARAM[atomidx2]-0.4;

						if(tdist<dspair2){
							penal=dspair2-tdist;
							penal=sqrt(penal*penal);

							rpp[i][j].pval[0]+=penal;
							rpp[i][j].pcw[0]+=1;
							nmaprpp[i][j]=1; nmaprpp[j][i]=1;
						}
					}//jj
				}//j
			}//ii
		}//i


		for(i=0;i<numseq;i++)
		{	
			for(j=i;j<numseq;j++)
			{	
				if(!mapflagrpp[i][j]) continue;
				inew[0]+=rpp[i][j].pcw[0];//cla
				fnew[0]+=rpp[i][j].pval[0];//cla

				fnew[1]+=rpp[i][j].pval[1];//cla
			}//j
		}//i

		for(i=0;i<CLASH_ENERGY_NUM;i++)
		{
			lp.nn[CLASH_ENERGY_ST+i]+=(inew[i]-iold[i]);
			lp.fn[CLASH_ENERGY_ST+i]+=(fnew[i]-fold[i]);
		}
		lp.fn[OXHBOND_ENERGY_ST]+=(fnew[1]-fold[1]);
	}

	nene=0; fene=0;
	for(i=0;i<CLASH_ENERGY_NUM;i++)
	{
		nene+=lp.nn[CLASH_ENERGY_ST+i];
		fene+=lp.fn[CLASH_ENERGY_ST+i];
	}

	hfene=lp.fn[OXHBOND_ENERGY_ST];

	for(i=0;i<numseq;i++){
		lp.indn[CLASH_ENERGY_ST][i]=0;
		lp.indf[CLASH_ENERGY_ST][i]=0;
	}
	for(i=0;i<numseq;i++)
	{	
		for(j=i;j<numseq;j++)
		{	
			lp.indn[CLASH_ENERGY_ST][i]+=rpp[i][j].pcw[0];
			lp.indn[CLASH_ENERGY_ST][j]+=rpp[i][j].pcw[0];
			lp.indf[CLASH_ENERGY_ST][i]+=rpp[i][j].pval[0];
			lp.indf[CLASH_ENERGY_ST][i]+=rpp[i][j].pval[0];
		}
	}
}






















