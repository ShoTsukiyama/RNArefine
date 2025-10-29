
#include "MonteCarloSimulation.hpp"

void MonteCarloSimulation::setLMPflap(double trandpos1, double trandpos2)
{
	int i,j;
	int atom1, atom2, atom3;

	for(i=trandpos1;i<trandpos2;i++){
		char base=refinestr_temp[i].aaa;

        for(j=0;j<MCLENGTHPARAM[base].size();j++){
            lengthflag[i][j] = 1;
        }

        for(j=0;j<MCINNERANGPARAM[base].size();j++){
            angleflag[i][j] = 1;
        }

        for(j=0;j<MCTORSIONPARAM[base].size();j++){
            torsionflag[i][j] = 1;
        }

        for(j=0;j<MCHYDLENGTHPARAM[base].size();j++){
            hydlengthflag[i][j] = 1;
        }

        for(j=0;j<MCHYDINNERANGPARAM[base].size();j++){
            hydangleflag[i][j] = 1;
        }

        for(j=0;j<MCHYDTORSIONPARAM[base].size();j++){
            hydtorsionflag[i][j] = 1;
        }

        for(j=0;j<AddLengthParamsList[base].size();j++){
            addlengthflag[i][j] = 1;
        }

        for(j=0;j<AddAngleParamsList[base].size();j++){
            addangleflag[i][j] = 1;
        }

		for(j=0;j<AddHydLengthParamsList[base].size();j++)
		{
			hydaddlengthflag[i][j]=1;
		}

		for(j=0;j<AddHydAngleParamsList[base].size();j++)
		{
			hydaddangleflag[i][j]=1;
		}
	}

	for(i=0;i<numbb;i++)
	{
		ntflag[i] = 1;
		for(j=i;j<numbb;j++)
		{
			if(j<trandpos1 || i>trandpos2){}
			else
			{
				mapbbflag[i][j]=1;
				mapbbflag[j][i]=1;
				mapflag[i][j]=1;
				mapflag[j][i]=1;
			}
		}
	}
}

void MonteCarloSimulation::singlemoveLMP(int k,int flagdir)
{
	BasicFunc bf;
	Point3d rp;
	double rnorm;
	char base1, base2;
	int tmp;

	int rtype = int(Random()*2.0);

	if(rtype == 0 & flagdir==0){
		base1=refinestr_temp[k-1].aaa;
		base2=refinestr_temp[k].aaa;

		int mtype1 = int(Random()*AccAtomParamNum[base1][2]);
		int mtype2 = int(Random()*AccAtomParamNum[base2][2]);

		rp=bf.setv(refinestr_temp[k].ptv[mtype2].x-refinestr_temp[k-1].ptv[mtype1].x,refinestr_temp[k].ptv[mtype2].y-refinestr_temp[k-1].ptv[mtype1].y,refinestr_temp[k].ptv[mtype2].z-refinestr_temp[k-1].ptv[mtype1].z);
		rnorm=bf.norm(rp);
		
		if(rnorm<EPS)//two in the same place
		{
			refinestr_temp[k].ptv[mtype2].x+=0.0001f; 
			refinestr_temp[k].ptv[mtype2].y+=0.0001f; 
			refinestr_temp[k].ptv[mtype2].z+=0.0001f;
			rp=bf.setv(refinestr_temp[k].ptv[mtype2].x-refinestr_temp[k-1].ptv[mtype1].x,refinestr_temp[k].ptv[mtype2].y-refinestr_temp[k-1].ptv[mtype1].y,refinestr_temp[k].ptv[mtype2].z-refinestr_temp[k-1].ptv[mtype1].z);
			rnorm=bf.norm(rp);
		}

		if(rnorm<PAIRLENGADJPARAM[base1][base2][mtype1][mtype2].lower || rnorm>PAIRLENGADJPARAM[base1][base2][mtype1][mtype2].upper)
		{
			rp=bf.scal(rp,PAIRLENGADJPARAM[base1][base2][mtype1][mtype2].mean/rnorm);
			refinestr_temp[k].ptv[mtype2].x=refinestr_temp[k-1].ptv[mtype1].x+rp.x;
			refinestr_temp[k].ptv[mtype2].y=refinestr_temp[k-1].ptv[mtype1].y+rp.y;
			refinestr_temp[k].ptv[mtype2].z=refinestr_temp[k-1].ptv[mtype1].z+rp.z;
		}
	}else if(rtype == 1 & flagdir==0){
		int mtype1 = int(Random()*6.0);
		int mtype2 = int(Random()*5.0);
		if(mtype2 >= mtype1) mtype2+=1;

		if(mtype1 > mtype2)
		{
			tmp=mtype1;
			mtype1=mtype2;
			mtype2=tmp;
		}

		base1=refinestr_temp[k].aaa;
		base2=refinestr_temp[k].aaa;

		rp=bf.setv(refinestr_temp[k].ptv[mtype2].x-refinestr_temp[k].ptv[mtype1].x,refinestr_temp[k].ptv[mtype2].y-refinestr_temp[k].ptv[mtype1].y,refinestr_temp[k].ptv[mtype2].z-refinestr_temp[k].ptv[mtype1].z);
		rnorm=bf.norm(rp);
		if(rnorm<EPS)//two in the same place
		{
			refinestr_temp[k].ptv[mtype2].x+=0.0001f; refinestr_temp[k].ptv[mtype2].y+=0.0001f; refinestr_temp[k].ptv[mtype2].z+=0.0001f;
			rp=bf.setv(refinestr_temp[k].ptv[mtype2].x-refinestr_temp[k].ptv[mtype1].x,refinestr_temp[k].ptv[mtype2].y-refinestr_temp[k].ptv[mtype1].y,refinestr_temp[k].ptv[mtype2].z-refinestr_temp[k].ptv[mtype1].z);
			rnorm=bf.norm(rp);
		}

		if(rnorm<ALLPAIRLENGNTPARAM[base1][mtype1][mtype2].lower || rnorm>ALLPAIRLENGNTPARAM[base1][mtype1][mtype2].upper)
		{
			rp=bf.scal(rp,ALLPAIRLENGNTPARAM[base1][mtype1][mtype2].mean/rnorm);
			refinestr_temp[k].ptv[mtype2].x=refinestr_temp[k].ptv[mtype1].x+rp.x;
			refinestr_temp[k].ptv[mtype2].y=refinestr_temp[k].ptv[mtype1].y+rp.y;
			refinestr_temp[k].ptv[mtype2].z=refinestr_temp[k].ptv[mtype1].z+rp.z;
		}
	}else if(rtype == 0 & flagdir==1){
		int mtype1 = int(Random()*6.0);
		int mtype2 = int(Random()*6.0);

		base1=refinestr_temp[k].aaa;
		base2=refinestr_temp[k-1].aaa;

		rp=bf.setv(refinestr_temp[k-1].ptv[mtype1].x-refinestr_temp[k].ptv[mtype2].x, refinestr_temp[k-1].ptv[mtype1].y-refinestr_temp[k].ptv[mtype2].y, refinestr_temp[k-1].ptv[mtype1].z-refinestr_temp[k].ptv[mtype2].z);
		rnorm=bf.norm(rp);
		
		if(rnorm<EPS)//two in the same place
		{
			refinestr_temp[k-1].ptv[mtype1].x+=0.0001f; refinestr_temp[k-1].ptv[mtype1].y+=0.0001f; refinestr_temp[k-1].ptv[mtype1].z+=0.0001f;
			rp=bf.setv(refinestr_temp[k-1].ptv[mtype1].x-refinestr_temp[k].ptv[mtype2].x, refinestr_temp[k-1].ptv[mtype1].y-refinestr_temp[k].ptv[mtype2].y, refinestr_temp[k-1].ptv[mtype1].z-refinestr_temp[k].ptv[mtype2].z);
			rnorm=bf.norm(rp);
		}

		if(rnorm<PAIRLENGADJPARAM[base2][base1][mtype1][mtype2].lower || rnorm>PAIRLENGADJPARAM[base2][base1][mtype1][mtype2].upper)
		{
			rp=bf.scal(rp,PAIRLENGADJPARAM[base2][base1][mtype1][mtype2].mean/rnorm);
			refinestr_temp[k-1].ptv[mtype1].x=refinestr_temp[k].ptv[mtype2].x+rp.x;
			refinestr_temp[k-1].ptv[mtype1].y=refinestr_temp[k].ptv[mtype2].y+rp.y;
			refinestr_temp[k-1].ptv[mtype1].z=refinestr_temp[k].ptv[mtype2].z+rp.z;
		}
	}else if(rtype == 1 & flagdir==1){
		int mtype1 = int(Random()*6.0);
		int mtype2 = int(Random()*5.0);
		if(mtype2 >= mtype1) mtype2+=1;

		if(mtype1 < mtype2){
			tmp=mtype1;
			mtype1=mtype2;
			mtype2=tmp;
		}

		base1=refinestr_temp[k-1].aaa;
		base2=refinestr_temp[k-1].aaa;

		rp=bf.setv(refinestr_temp[k-1].ptv[mtype2].x-refinestr_temp[k-1].ptv[mtype1].x,refinestr_temp[k-1].ptv[mtype2].y-refinestr_temp[k-1].ptv[mtype1].y,refinestr_temp[k-1].ptv[mtype2].z-refinestr_temp[k-1].ptv[mtype1].z);
		rnorm=bf.norm(rp);
		if(rnorm<EPS)//two in the same place
		{
			refinestr_temp[k-1].ptv[mtype2].x+=0.0001f; refinestr_temp[k-1].ptv[mtype2].y+=0.0001f; refinestr_temp[k-1].ptv[mtype2].z+=0.0001f;
			rp=bf.setv(refinestr_temp[k-1].ptv[mtype2].x-refinestr_temp[k-1].ptv[mtype1].x,refinestr_temp[k-1].ptv[mtype2].y-refinestr_temp[k-1].ptv[mtype1].y,refinestr_temp[k-1].ptv[mtype2].z-refinestr_temp[k-1].ptv[mtype1].z);
			rnorm=bf.norm(rp);
		}

		if(rnorm<ALLPAIRLENGNTPARAM[base1][mtype1][mtype2].lower || rnorm>ALLPAIRLENGNTPARAM[base1][mtype1][mtype2].upper)
		{
			rp=bf.scal(rp,ALLPAIRLENGNTPARAM[base1][mtype1][mtype2].mean/rnorm);
			refinestr_temp[k-1].ptv[mtype2].x=refinestr_temp[k-1].ptv[mtype1].x+rp.x;
			refinestr_temp[k-1].ptv[mtype2].y=refinestr_temp[k-1].ptv[mtype1].y+rp.y;
			refinestr_temp[k-1].ptv[mtype2].z=refinestr_temp[k-1].ptv[mtype1].z+rp.z;
		}
	}
}


int MonteCarloSimulation::selecttargetposLMP(int &trandpos, double r)
{
	double randval=Random();

	if(r > randval){
		int i,j;
		char base;
		vector<double> prob;
		prob.resize(numbb, 0.0);
		double totalprob=0.0;

		for(i=0;i<numbb;i++){
			base=refinestr[i].aaa;
			for(j=0;j<AccAtomParamNum[base][4];j++){
				const double e = std::max(0.0, pairwiseenergyatoms[i][j]);
				prob[i]+=e;
				totalprob+=e;
			}
		}

		if (!(totalprob > 0.0)) return 0;

		double wpos=Random()*totalprob;

		double currentenen=0;
		double nextene;
		double tenergy=0;
		for(i=0;i<numbb;i++){
			base=refinestr[i].aaa;
			nextene=currentenen+prob[i];
			if(currentenen <= wpos && wpos < nextene){
				trandpos=i;
				return 1;
			}
			currentenen=nextene;
		}
	}else{
		trandpos=Random()*numbb;
	}

	return 0;
}


Point3s ranv(double fac)
{
	Point3s tp;
	double dtheta,dphi,dr;
	dr=0.00001+fac*Random();
	dtheta=(1.0*Random())*PI;
	dphi=(2.0*Random()-1.0)*PI;
	tp.x=dr*sin(dtheta)*cos(dphi);
	tp.y=dr*sin(dtheta)*sin(dphi);
	tp.z=dr*cos(dtheta);
	return(tp);
}

bool MonteCarloSimulation::fmovementLMP(double oldenergy, double &newenergy)
{
	int i, j, k;
	bool flagtor;
	int trandpos, trandposst, trandposed,span;
	double mean, std, upper, lower;
	int idx;
	double rmax, rnorm;
	Point3s rp;
	int mtype1, mtype2;
	char base, base1, base2;
	bool flagdone;
	bool exitOuterLoop = false;

	int selflag = selecttargetposLMP(trandpos, 0.8);

	if(selflag == 0){
		newenergy=oldenergy+1000000000;
		return false;
	}

	base=refinestr[trandpos].aaa;

	//cout<<"=============="<<endl;
	//cout<<trandpos<<endl;

	StrCopy(refinestr, refinestr_temp);
	resetfrag(0);

	span=2+int(6*Random());
	trandposst=trandpos-int(span/2);
	trandposed=trandposst+span;
	if(trandposst<1) trandposst=1;
	if(trandposed>numbb-1) trandposed=numbb-1;
	span=trandposed-trandposst;
	if(span < 2){
		newenergy=oldenergy+1000000000;
		return false;
	}

	int threshiter=span*200;
	if(threshiter>3000) threshiter=3000;

	for(i=trandposst;i<trandposed;i++)
	{
		for(j=0;j<AccAtomParamNum[base][2];j++){
			//rmax=0.25+0.25*Random();
			rmax=0.5+0.5*Random();
			rp=ranv(rmax);
			refinestr_temp[i].ptv[j].x+=rp.x;
			refinestr_temp[i].ptv[j].y+=rp.y;
			refinestr_temp[i].ptv[j].z+=rp.z;
		}
	}

	int numiter=0;
	do{	
		for(i=trandposst;i<trandposed;i++)
		{
			for(j=0;j<100;j++){
				singlemoveLMP(i,0);
			}
		}

		for(i=trandposed;i>trandposst;i--)
		{
			for(j=0;j<100;j++){
				singlemoveLMP(i,1);
			}
		}

		flagdone=true;
		exitOuterLoop = false;
		for(k=trandposst;k<=trandposed;k++)
		{
			base1=refinestr_temp[k-1].aaa;
			base2=refinestr_temp[k].aaa;
			for(mtype1=0;mtype1<AccAtomParamNum[base1][2];mtype1++)
			{
				for(mtype2=0;mtype2<AccAtomParamNum[base2][2];mtype2++)
				{
					rp=bf.setvs(refinestr_temp[k].ptv[mtype2].x-refinestr_temp[k-1].ptv[mtype1].x,refinestr_temp[k].ptv[mtype2].y-refinestr_temp[k-1].ptv[mtype1].y,refinestr_temp[k].ptv[mtype2].z-refinestr_temp[k-1].ptv[mtype1].z);
					rnorm=bf.norm(rp);

					if(rnorm<PAIRLENGADJPARAM[base1][base2][mtype1][mtype2].lower || rnorm>PAIRLENGADJPARAM[base1][base2][mtype1][mtype2].upper)
					{
						flagdone=false;
						exitOuterLoop=true;
						break;
					}
				}
				if(exitOuterLoop) break;
			}
			if(exitOuterLoop) break;
		}

		if(flagdone){
			for(k=trandposst;k<trandposed;k++)
			{
				base1=refinestr_temp[k].aaa;
				base2=refinestr_temp[k].aaa;

				for(mtype1=0;mtype1<AccAtomParamNum[base1][2];mtype1++)
				{
					for(mtype2=mtype1+1;mtype2<AccAtomParamNum[base2][2];mtype2++)
					{
						rp=bf.setvs(refinestr_temp[k].ptv[mtype2].x-refinestr_temp[k].ptv[mtype1].x,refinestr_temp[k].ptv[mtype2].y-refinestr_temp[k].ptv[mtype1].y,refinestr_temp[k].ptv[mtype2].z-refinestr_temp[k].ptv[mtype1].z);
						rnorm=bf.norm(rp);
						if(rnorm<ALLPAIRLENGNTPARAM[base1][mtype1][mtype2].lower || rnorm>ALLPAIRLENGNTPARAM[base1][mtype1][mtype2].upper)
						{
							flagdone=false;
							exitOuterLoop=true;
							break;
						}
					}
					if(exitOuterLoop) break;
				}
				if(exitOuterLoop) break;
			}
		}
		numiter++;
	}while(!flagdone && numiter<threshiter);

	int istart=trandposst-1;
	int iend=trandposed+2;
	if(iend>=numbb) iend=numbb;

	coord2parambb(refinestr_temp,trandposst,iend);
	coord2paramsugar(refinestr_temp,trandposst,iend);
	param2coordbase(refinestr_temp,istart,iend);
	param2coordhydrogen(refinestr_temp,trandposst,iend);
	coord2paramadd(refinestr_temp,istart,iend);
	coord2paramhydrogen(refinestr_temp,istart,iend);
	coord2addparamhydrogen(refinestr_temp,istart,iend);

	setLMPflap(istart, iend);
	newenergy = calcenergymain();

	bool upflag=metropolisUpdate(oldenergy, newenergy);
	if(upflag){
		StrCopy(refinestr_temp, refinestr);
		initializeatomenergy();
		updateenergies();
	}
	
	/*
	cout<<"===================="<<endl;
	cout<<"QqqqqqqqqqqqqqqqqQ"<<endl;
	checkupdate();
	cout<<"PppppppppppppppppP"<<endl;
	resetfrag(1);
	double newenergy_tmp = calcenergymain();

	if(abs(newenergy - newenergy_tmp)/max(abs(newenergy), abs(newenergy_tmp)) > 0.002){
		cout<<"Energy Error"<<endl;
		cout<<newenergy << "-" << newenergy_tmp<<endl;
		exit(1);
	}else{
		cout<<"OK"<<endl;
	}
	*/

	return upflag;
}

