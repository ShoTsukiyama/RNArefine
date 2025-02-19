// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

void RNASimulation::getoutenebasepairhydlength(double val, double mean, double lower, double upper, double k, double &diff, int &bdiff){
	double tmp;

	if((lower <= val) && (val <= upper)){
		bdiff=0;
		diff=0;
	}else{
		bdiff=1;
		tmp=(val-mean);
		diff=k*sqrt(tmp*tmp);
	}
}

void RNASimulation::getoutenebasepairhydangle(double val, double mean, double lower, double upper, double k, double &diff, int &bdiff){
	double tmp;

	if((lower <= val) && (val <= upper)){
		bdiff=0;
		diff=0;
	}else{
		bdiff=1;
		tmp=cos((val-mean)*RADDEG);
		diff=k*exp(-tmp);
	}
}

void RNASimulation::calchbondene(vector<Residueatoms> decstr, char daat, char aaat, int dbase, int abase, int pairtype, string hbtype, int &nn, double &fn, vector<double> &hbatoms)
{
	int atom1, atom2, atom3, atom4;
	Point3s tmpd,tmpa,tmph;
	Point3s ptha, ptda, ptdh;
	double hlen1, hlen2, hang;
	int hbondidx;
	double fene, fenetmp, ene;
	int bene, benetmp;
	vector<double> hbatomstmp;

	vector<vector<int>> atomslist=BasePairHBondInfo[daat][aaat][pairtype][hbtype];
	hbondidx=HbondParamIdx[daat][aaat][pairtype][hbtype];

	fn=100000;
	for(int i=0;i<atomslist.size();i++){
		hbatomstmp.clear();
		atom1=atomslist[i][0]; atom2=atomslist[i][1]; atom3=atomslist[i][2];
		tmpd=decstr[dbase].ptv[atom1];
		tmph=decstr[dbase].ptv[atom2];
		tmpa=decstr[abase].ptv[atom3];

		hbatomstmp.assign({static_cast<double>(dbase),static_cast<double>(dbase),static_cast<double>(abase),static_cast<double>(atom1),static_cast<double>(atom2),static_cast<double>(atom3),BASEPAIRHYDPARAM[daat][aaat][pairtype][hbtype][1][0].mean,BASEPAIRHYDPARAM[daat][aaat][pairtype][hbtype][0][1].mean,BASEPAIRHYDPARAM[daat][aaat][pairtype][hbtype][0][0].mean});

		ptha.x=double(tmpa.x-tmph.x); ptha.y=double(tmpa.y-tmph.y); ptha.z=double(tmpa.z-tmph.z); //acceptor-hydrogen
		hlen1=bf.norm(ptha);

		ptda.x=double(tmpd.x-tmpa.x); ptda.y=double(tmpd.y-tmpa.y); ptda.z=double(tmpd.z-tmpa.z); //donor-acceptor
		hlen2=bf.norm(ptda);

		fene=0; bene=0;
		ptdh.x=double(tmpd.x-tmph.x); ptdh.y=double(tmpd.y-tmph.y); ptdh.z=double(tmpd.z-tmph.z); //donor-hydrogen
		
		//donor-hydrogen-acceptor
		hang=bf.angv(ptdh,ptha); 
		getoutenebasepairhydangle(hang, BASEPAIRHYDPARAM[daat][aaat][pairtype][hbtype][1][0].mean, BASEPAIRHYDPARAM[daat][aaat][pairtype][hbtype][1][0].lower, BASEPAIRHYDPARAM[daat][aaat][pairtype][hbtype][1][0].upper, 1.0, fenetmp, benetmp);
		fene+=fenetmp;
		bene+=benetmp;

		//donor-acceptor
		getoutenebasepairhydlength(hlen2, BASEPAIRHYDPARAM[daat][aaat][pairtype][hbtype][0][1].mean, BASEPAIRHYDPARAM[daat][aaat][pairtype][hbtype][0][1].lower, BASEPAIRHYDPARAM[daat][aaat][pairtype][hbtype][0][1].upper, 1.0, fenetmp, benetmp);
		fene+=fenetmp;
		bene+=benetmp;

		//hydrogen-acceptor
		getoutenebasepairhydlength(hlen1, BASEPAIRHYDPARAM[daat][aaat][pairtype][hbtype][0][0].mean, BASEPAIRHYDPARAM[daat][aaat][pairtype][hbtype][0][0].lower, BASEPAIRHYDPARAM[daat][aaat][pairtype][hbtype][0][1].upper, 1.0, fenetmp, benetmp);
		fene+=fenetmp;
		bene+=benetmp;

		//ene=getlengschematicene(hlen1, BASEPARAM[daat][aaat][pairtype].hbondcutoff[hbondidx][0], 1.0);
		//if(ene>0) bene=1;
		//fene+=ene;
		//ene=getlengschematicene(hlen2, BASEPARAM[daat][aaat][pairtype].hbondcutoff[hbondidx][1], 1.0);
		//if(ene>0) bene=1;
		//fene+=ene;
		//ene=gethbondangschematicene(hang, BASEPARAM[daat][aaat][pairtype].hbondcutoff[hbondidx][2], 1.0);
		//if(ene>0) bene=1;
		//fene+=ene;
		if(fn > fene){
			nn=bene;
			fn=fene;
			hbatoms=hbatomstmp;
		}
	}
}


void RNASimulation::fenergybasepairhyd(vector<Residueatoms> decstr,int numseq,double &nene,double &fene,Energyterms &lp, vector<Basecontact> &baserpp, vector<vector<bool>> mapflagrpp,bool flagacc){
	int i,j,k;
	int basei,basej;
	char aat1, aat2;
	double hbondfntotal1, hbondfntotal2;
	int hbondnn21, hbondnn22;
	double totalfn1, mintotalfn1;
	double totalfn2, mintotalfn2;
	int minnn2, nn2;
	int hbondnn2;
	double hbondfntotal21, hbondfntotal22;
	double fntmp;
	int nntmp;
	int pairwiseidx;
	vector<double> hatomtypetmp;
	vector<vector<double>> hatomtype1, hatomtype2, hatomtype;
	vector<double> indf(BASEHYD_ENERGY_NUM,0);
	vector<int> indn(BASEHYD_ENERGY_NUM,0);
	

	if(!flagacc)
	{
		for(j=0;j<BASEHYD_ENERGY_NUM;j++)
		{
			lp.fn[BASEHYD_ENERGY_ST+j]=0;
			lp.nn[BASEHYD_ENERGY_ST+j]=0;
		}
		
		for(i=0;i<BASEPAIRS.size();i++)
		{
			basei=BASEPAIRS[i].basepidx.first;//base idx 1
			basej=BASEPAIRS[i].basepidx.second;//base idx 2

			pairwiseidx=bpmapid2idx[basei*numseq+basej];

			if(basei < basej)
			{
				aat1=decstr[basei].aaa;
				aat2=decstr[basej].aaa;

				for(j=0;j<BASEHYD_ENERGY_NUM;j++)
				{
					indf[j]=0; indn[j]=0;
				}

				mintotalfn2=1000000;
				for(const auto& pair : BasePairADorder[aat1][aat2])
				{
					totalfn1=0; totalfn2=0; nn2=0;
					if(aat1 == aat2){
						hbondfntotal21=0, hbondfntotal22=0;
						hbondnn21=0; hbondnn22=0;
						hatomtype1.clear(); hatomtype2.clear();
						for(const auto& hb : pair.second)
						{
							calchbondene(decstr, aat1, aat2, basei, basej, pair.first, hb.first, nntmp, fntmp, hatomtypetmp);
							hbondfntotal21+=fntmp;
							hbondnn21=nntmp;
							hatomtype1.push_back(hatomtypetmp);

							calchbondene(decstr, aat2, aat1, basej, basei, pair.first, hb.first, nntmp, fntmp, hatomtypetmp);
							hbondfntotal22+=fntmp;
							hbondnn22=nntmp;
							hatomtype2.push_back(hatomtypetmp);
						}

						hbondfntotal21=hbondfntotal21;
						hbondfntotal22=hbondfntotal22;

						if(hbondfntotal21<hbondfntotal22){
							totalfn2+=hbondfntotal21;
							hatomtype=hatomtype1;
							//nn2=hbondnn21;
							if(hbondnn21 >= 1){
								nn2=1;
							}else{
								nn2=0;
							}
						}
						else{
							totalfn2+=hbondfntotal22;
							hatomtype=hatomtype2;
							//nn2=hbondnn22;
							if(hbondnn22 >= 1){
								nn2=1;
							}else{
								nn2=0;
							}
						}
					}else{
						hbondfntotal2=0;
						hbondnn2=0;

						hatomtype.clear();
						for(const auto& hb : pair.second){
							if(BasePairADorder[aat1][aat2][pair.first][hb.first][0] == 0){
								calchbondene(decstr, aat1, aat2, basei, basej, pair.first, hb.first, nntmp, fntmp, hatomtypetmp);
							}
							else{
								calchbondene(decstr, aat2, aat1, basej, basei, pair.first, hb.first, nntmp, fntmp, hatomtypetmp);
							}

							hatomtype.push_back(hatomtypetmp);
							hbondfntotal2+=fntmp;
							hbondnn2=nntmp;
						}

						hbondfntotal2=hbondfntotal2;
						totalfn2+=hbondfntotal2;

						if(hbondnn2 >= 1){
							nn2=1;
						}else{
							nn2=0;
						}
						//nn2=hbondnn2;
					}

					if(totalfn2 < mintotalfn2){
						mintotalfn2=totalfn2;
						minnn2=nn2;
						BASEPAIRS[i].hydtype=pair.first;
						BASEPAIRS[i].hatomtype=hatomtype;
					}
				}

				indf[0]+=mintotalfn2;
				indn[0]=minnn2;

				for(j=0;j<BASEHYD_ENERGY_NUM;j++)
				{
					baserpp[pairwiseidx].pval[BASE_ENERGY_NUM+j]=indf[j];
					baserpp[pairwiseidx].pval[BASE_ENERGY_NUM+j]=indf[j];
					baserpp[pairwiseidx].pcw[BASE_ENERGY_NUM+j]=indn[j];
					baserpp[pairwiseidx].pcw[BASE_ENERGY_NUM+j]=indn[j];

					lp.indf[BASEHYD_ENERGY_ST+j][basei]+=indf[j];
					lp.indf[BASEHYD_ENERGY_ST+j][basej]+=indf[j];
					lp.indn[BASEHYD_ENERGY_ST+j][basei]+=indn[j];
					lp.indn[BASEHYD_ENERGY_ST+j][basej]+=indn[j];
					lp.fn[BASEHYD_ENERGY_ST+j]+=indf[j];
					lp.nn[BASEHYD_ENERGY_ST+j]+=indn[j];
				}
			}
		}
	}
	else
	{
		//int idx = 0;
		int iold[BASEHYD_ENERGY_NUM],inew[BASEHYD_ENERGY_NUM];
		double fold[BASEHYD_ENERGY_NUM],fnew[BASEHYD_ENERGY_NUM];
		for(i=0;i<BASEHYD_ENERGY_NUM;i++){
			iold[i]=0;inew[i]=0;
			fold[i]=0;fnew[i]=0;
		}

		for(i=0;i<BASEPAIRS.size();i++)
		{
			basei=BASEPAIRS[i].basepidx.first;//base idx 1
			basej=BASEPAIRS[i].basepidx.second;//base idx 2

			pairwiseidx=bpmapid2idx[basei*numseq+basej];

			if(basei<basej)
			{

				if(!mapflagrpp[basei][basej]) continue;
			
				for(j=0;j<BASEHYD_ENERGY_NUM;j++)
				{
					fold[j]+=baserpp[pairwiseidx].pval[BASE_ENERGY_NUM+j];
					iold[j]+=baserpp[pairwiseidx].pcw[BASE_ENERGY_NUM+j];

					lp.indf[BASEHYD_ENERGY_ST+j][basei]=lp.indf[BASEHYD_ENERGY_ST+j][basei]-baserpp[pairwiseidx].pval[BASE_ENERGY_NUM+j];
					lp.indf[BASEHYD_ENERGY_ST+j][basej]=lp.indf[BASEHYD_ENERGY_ST+j][basej]-baserpp[pairwiseidx].pval[BASE_ENERGY_NUM+j];
					lp.indn[BASEHYD_ENERGY_ST+j][basei]=lp.indn[BASEHYD_ENERGY_ST+j][basei]-baserpp[pairwiseidx].pcw[BASE_ENERGY_NUM+j];
					lp.indn[BASEHYD_ENERGY_ST+j][basej]=lp.indn[BASEHYD_ENERGY_ST+j][basej]-baserpp[pairwiseidx].pcw[BASE_ENERGY_NUM+j];

					baserpp[pairwiseidx].pval[BASE_ENERGY_NUM+j]=0;
					baserpp[pairwiseidx].pcw[BASE_ENERGY_NUM+j]=0;
				}

				aat1=decstr[basei].aaa;
				aat2=decstr[basej].aaa;

				for(j=0;j<BASEHYD_ENERGY_NUM;j++)
				{
					indf[j]=0; indn[j]=0;
				}

				mintotalfn2=1000000;
				for(const auto& pair : BasePairADorder[aat1][aat2])
				{
					totalfn1=0; totalfn2=0; nn2=0;
					if(aat1 == aat2){
						hbondfntotal21=0, hbondfntotal22=0;
						hbondnn21=0; hbondnn22=0;
						hatomtype1.clear(); hatomtype2.clear();
						for(const auto& hb : pair.second)
						{
							calchbondene(decstr, aat1, aat2, basei, basej, pair.first, hb.first, nntmp, fntmp, hatomtypetmp);
							hbondfntotal21+=fntmp;
							hbondnn21=nntmp;
							hatomtype1.push_back(hatomtypetmp);

							calchbondene(decstr, aat2, aat1, basej, basei, pair.first, hb.first, nntmp, fntmp, hatomtypetmp);
							hbondfntotal22+=fntmp;
							hbondnn22=nntmp;
							hatomtype2.push_back(hatomtypetmp);
						}

						hbondfntotal21=hbondfntotal21;
						hbondfntotal22=hbondfntotal22;

						if(hbondfntotal21<hbondfntotal22){
							totalfn2+=hbondfntotal21;
							hatomtype=hatomtype1;
							//nn2=hbondnn21;
							if(hbondnn21 >= 1){
								nn2=1;
							}else{
								nn2=0;
							}
						}
						else{
							totalfn2+=hbondfntotal22;
							hatomtype=hatomtype2;
							//nn2=hbondnn22;
							if(hbondnn22 >= 1){
								nn2=1;
							}else{
								nn2=0;
							}
						}
					}else{
						hbondfntotal2=0;
						hbondnn2=0;

						hatomtype.clear();
						for(const auto& hb : pair.second){
							if(BasePairADorder[aat1][aat2][pair.first][hb.first][0] == 0){
								calchbondene(decstr, aat1, aat2, basei, basej, pair.first, hb.first, nntmp, fntmp, hatomtypetmp);
							}
							else{
								calchbondene(decstr, aat2, aat1, basej, basei, pair.first, hb.first, nntmp, fntmp, hatomtypetmp);
							}

							hatomtype.push_back(hatomtypetmp);
							hbondfntotal2+=fntmp;
							hbondnn2=nntmp;
						}

						hbondfntotal2=hbondfntotal2;
						totalfn2+=hbondfntotal2;

						if(hbondnn2 >= 1){
							nn2=1;
						}else{
							nn2=0;
						}
						//nn2=hbondnn2;
					}

					if(totalfn2 < mintotalfn2){
						mintotalfn2=totalfn2;
						minnn2=nn2;
						BASEPAIRS[i].hydtype=pair.first;
						BASEPAIRS[i].hatomtype=hatomtype;
					}
				}

				indf[0]+=mintotalfn2;
				indn[0]=minnn2;
				
				for(j=0;j<BASEHYD_ENERGY_NUM;j++)
				{
					baserpp[pairwiseidx].pval[BASE_ENERGY_NUM+j]=indf[j];
					baserpp[pairwiseidx].pcw[BASE_ENERGY_NUM+j]=indn[j];

					lp.indf[BASEHYD_ENERGY_ST+j][basei]+=baserpp[pairwiseidx].pval[BASE_ENERGY_NUM+j];
					lp.indf[BASEHYD_ENERGY_ST+j][basej]+=baserpp[pairwiseidx].pval[BASE_ENERGY_NUM+j];
					lp.indn[BASEHYD_ENERGY_ST+j][basei]+=baserpp[pairwiseidx].pcw[BASE_ENERGY_NUM+j];
					lp.indn[BASEHYD_ENERGY_ST+j][basej]+=baserpp[pairwiseidx].pcw[BASE_ENERGY_NUM+j];
					
					fnew[j]+=baserpp[pairwiseidx].pval[BASE_ENERGY_NUM+j];
					inew[j]+=baserpp[pairwiseidx].pcw[BASE_ENERGY_NUM+j];
				}
			}
		}

		for(i=0;i<BASEHYD_ENERGY_NUM;i++)
		{
			lp.fn[BASEHYD_ENERGY_ST+i]+=(fnew[i]-fold[i]);
			lp.nn[BASEHYD_ENERGY_ST+i]+=(inew[i]-iold[i]);
		}
	}

	for(i=0;i<BASEHYD_ENERGY_NUM;i++)
	{
		fene+=lp.fn[BASEHYD_ENERGY_ST+i];//base
		nene+=lp.nn[BASEHYD_ENERGY_ST+i];//base
	}
}






















