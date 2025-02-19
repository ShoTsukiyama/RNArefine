// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

void RNASimulation::singlemoveLMPfull(vector<Residueatoms> &restr,int k,int flagdir,int paramidx)
{
	BasicFunc bf;
	Point3d rp;
	double rnorm;
	char aat1, aat2;
	int tmp;

	int rtype = int(Random()*2.0);

	if(rtype == 0 & flagdir==0){
		int mtype1 = int(Random()*PAIRLENGTHTARGETMOVNUM);
		int mtype2 = int(Random()*PAIRLENGTHTARGETMOVNUM);

		aat1=restr[k-1].aaa;
		aat2=restr[k].aaa;

		rp=bf.setv(restr[k].ptv[mtype2].x-restr[k-1].ptv[mtype1].x,restr[k].ptv[mtype2].y-restr[k-1].ptv[mtype1].y,restr[k].ptv[mtype2].z-restr[k-1].ptv[mtype1].z);
		rnorm=bf.norm(rp);
		
		if(rnorm<EPS)//two in the same place
		{
			restr[k].ptv[mtype2].x+=0.0001f; restr[k].ptv[mtype2].y+=0.0001f; restr[k].ptv[mtype2].z+=0.0001f;
			rp=bf.setv(restr[k].ptv[mtype2].x-restr[k-1].ptv[mtype1].x,restr[k].ptv[mtype2].y-restr[k-1].ptv[mtype1].y,restr[k].ptv[mtype2].z-restr[k-1].ptv[mtype1].z);
			rnorm=bf.norm(rp);
		}

		if(rnorm<PAIRLENGADJPARAM[aat1][aat2][0][mtype1][mtype2].lower || rnorm>PAIRLENGADJPARAM[aat1][aat2][0][mtype1][mtype2].upper)
		{
			rp=bf.scal(rp,PAIRLENGADJPARAM[aat1][aat2][0][mtype1][mtype2].mean/rnorm);
			restr[k].ptv[mtype2].x=restr[k-1].ptv[mtype1].x+rp.x;
			restr[k].ptv[mtype2].y=restr[k-1].ptv[mtype1].y+rp.y;
			restr[k].ptv[mtype2].z=restr[k-1].ptv[mtype1].z+rp.z;
		}
	}else if(rtype == 1 & flagdir==0){
		int mtype1 = int(Random()*PAIRLENGTHTARGETMOVNUM);
		int mtype2 = int(Random()*(PAIRLENGTHTARGETMOVNUM-1));
		if(mtype2 >= mtype1) mtype2+=1;

		if(mtype1 > mtype2)
		{
			tmp=mtype1;
			mtype1=mtype2;
			mtype2=tmp;
		}
		
		aat1=restr[k].aaa;
		aat2=restr[k].aaa;

		rp=bf.setv(restr[k].ptv[mtype2].x-restr[k].ptv[mtype1].x,restr[k].ptv[mtype2].y-restr[k].ptv[mtype1].y,restr[k].ptv[mtype2].z-restr[k].ptv[mtype1].z);
		rnorm=bf.norm(rp);
		if(rnorm<EPS)//two in the same place
		{
			restr[k].ptv[mtype2].x+=0.0001f; restr[k].ptv[mtype2].y+=0.0001f; restr[k].ptv[mtype2].z+=0.0001f;
			rp=bf.setv(restr[k].ptv[mtype2].x-restr[k].ptv[mtype1].x,restr[k].ptv[mtype2].y-restr[k].ptv[mtype1].y,restr[k].ptv[mtype2].z-restr[k].ptv[mtype1].z);
			rnorm=bf.norm(rp);
		}
		if(rnorm<PAIRLENGNTPARAM[aat1][paramidx][mtype1][mtype2].lower || rnorm>PAIRLENGNTPARAM[aat1][paramidx][mtype1][mtype2].upper)
		{
			//cout<<rnorm<<endl;
			//cout<<to_string(rp.x) + "-" + to_string(rp.y) + "-" + to_string(rp.z)<<endl;
			rp=bf.scal(rp,PAIRLENGNTPARAM[aat1][paramidx][mtype1][mtype2].mean/rnorm);
			//cout<<to_string(rp.x) + "-" + to_string(rp.y) + "-" + to_string(rp.z)<<endl;
			restr[k].ptv[mtype2].x=restr[k].ptv[mtype1].x+rp.x;
			restr[k].ptv[mtype2].y=restr[k].ptv[mtype1].y+rp.y;
			restr[k].ptv[mtype2].z=restr[k].ptv[mtype1].z+rp.z;
		}
	}else if(rtype == 0 & flagdir==1){
		int mtype1 = int(Random()*PAIRLENGTHTARGETMOVNUM);
		int mtype2 = int(Random()*PAIRLENGTHTARGETMOVNUM);

		aat1=restr[k].aaa;
		aat2=restr[k-1].aaa;

		rp=bf.setv(restr[k-1].ptv[mtype1].x-restr[k].ptv[mtype2].x, restr[k-1].ptv[mtype1].y-restr[k].ptv[mtype2].y, restr[k-1].ptv[mtype1].z-restr[k].ptv[mtype2].z);
		rnorm=bf.norm(rp);
		
		if(rnorm<EPS)//two in the same place
		{
			restr[k-1].ptv[mtype1].x+=0.0001f; restr[k-1].ptv[mtype1].y+=0.0001f; restr[k-1].ptv[mtype1].z+=0.0001f;
			rp=bf.setv(restr[k-1].ptv[mtype1].x-restr[k].ptv[mtype2].x, restr[k-1].ptv[mtype1].y-restr[k].ptv[mtype2].y, restr[k-1].ptv[mtype1].z-restr[k].ptv[mtype2].z);
			rnorm=bf.norm(rp);
		}
		if(rnorm<PAIRLENGADJPARAM[aat2][aat1][0][mtype1][mtype2].lower || rnorm>PAIRLENGADJPARAM[aat2][aat1][0][mtype1][mtype2].upper)
		{
			//cout<<rnorm<<endl;
			//cout<<to_string(rp.x) + "-" + to_string(rp.y) + "-" + to_string(rp.z)<<endl;
			rp=bf.scal(rp,PAIRLENGADJPARAM[aat2][aat1][0][mtype1][mtype2].mean/rnorm);
			//cout<<to_string(rp.x) + "-" + to_string(rp.y) + "-" + to_string(rp.z)<<endl;
			restr[k-1].ptv[mtype1].x=restr[k].ptv[mtype2].x+rp.x;
			restr[k-1].ptv[mtype1].y=restr[k].ptv[mtype2].y+rp.y;
			restr[k-1].ptv[mtype1].z=restr[k].ptv[mtype2].z+rp.z;
		}
	}else if(rtype == 1 & flagdir==1){
		int mtype1 = int(Random()*PAIRLENGTHTARGETMOVNUM);
		int mtype2 = int(Random()*(PAIRLENGTHTARGETMOVNUM-1));
		if(mtype2 >= mtype1) mtype2+=1;

		if(mtype1 < mtype2){
			tmp=mtype1;
			mtype1=mtype2;
			mtype2=tmp;
		}
		
		aat1=restr[k-1].aaa;
		aat2=restr[k-1].aaa;

		rp=bf.setv(restr[k-1].ptv[mtype2].x-restr[k-1].ptv[mtype1].x,restr[k-1].ptv[mtype2].y-restr[k-1].ptv[mtype1].y,restr[k-1].ptv[mtype2].z-restr[k-1].ptv[mtype1].z);
		rnorm=bf.norm(rp);
		if(rnorm<EPS)//two in the same place
		{
			restr[k-1].ptv[mtype2].x+=0.0001f; restr[k-1].ptv[mtype2].y+=0.0001f; restr[k-1].ptv[mtype2].z+=0.0001f;
			rp=bf.setv(restr[k-1].ptv[mtype2].x-restr[k-1].ptv[mtype1].x,restr[k-1].ptv[mtype2].y-restr[k-1].ptv[mtype1].y,restr[k-1].ptv[mtype2].z-restr[k-1].ptv[mtype1].z);
			rnorm=bf.norm(rp);
		}
		if(rnorm<PAIRLENGNTPARAM[aat1][paramidx][mtype1][mtype2].lower || rnorm>PAIRLENGNTPARAM[aat1][paramidx][mtype1][mtype2].upper)
		{
			//cout<<rnorm<<endl;
			//cout<<to_string(rp.x) + "-" + to_string(rp.y) + "-" + to_string(rp.z)<<endl;
			rp=bf.scal(rp,PAIRLENGNTPARAM[aat1][paramidx][mtype1][mtype2].mean/rnorm);
			//cout<<to_string(rp.x) + "-" + to_string(rp.y) + "-" + to_string(rp.z)<<endl;
			restr[k-1].ptv[mtype2].x=restr[k-1].ptv[mtype1].x+rp.x;
			restr[k-1].ptv[mtype2].y=restr[k-1].ptv[mtype1].y+rp.y;
			restr[k-1].ptv[mtype2].z=restr[k-1].ptv[mtype1].z+rp.z;
		}
	}
}

void RNASimulation::settmplmpfullflap(double numseq, double trandpos1, double trandpos2)
{
	int i,j;
	int st, ed;

	for(i=trandpos1;i<trandpos2;i++){
		for(j=0;j<BLENGTHTARGET.size();j++){
			lp2.flap[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][i]=true;
		}

		for(j=0;j<INNERANGLETARGET.size();j++){
			lp2.flap[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][i]=true;
		}

		for(j=0;j<TORSIONANGLETARGET.size();j++){
			lp2.flap[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][i]=true;
		}

		for(j=0;j<ADDBLENGTHTARGET.size();j++)
		{
			lp2.flap[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]][i]=true;
		}
		
		for(j=0;j<ADDINNERANGLETARGET.size();j++)
		{
			lp2.flap[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]][i]=true;
		}
	}
	
	for(i=0;i<numseq;i++)
	{
		for(j=i;j<numseq;j++)
		//for(j=i+1;j<numseq;j++)
		{
			if(j<trandpos1 || i>trandpos2){}
			else
			{
				mapbbflagbk[i][j]=true;
				mapbbflagbk[j][i]=true;
				mapflagbk[i][j]=true;
				mapflagbk[j][i]=true;
			}
		}
	}
}

void RNASimulation::setlmpparam(vector<Residueatoms> &restr, int trandposst, int trandposed, vector<int> &ntparamidx, vector<int> &adjparamidx)
{
	int i, j;
	char aat1, aat2;
	int ls, minls;
	int idx;
	double randval;
	int mtype1, mtype2;
	double rmax, rnorm;
	Point3d rp;
	
	for(j=trandposst;j<trandposed;j++)
	{
		if(restr[i].dih[5] > 110){
			if(restr[i].dih[4] > 120 & restr[i].dih[4]<=250){
				idx=0;
			}else{
				idx=2;
			}
		}else{
			if(restr[i].dih[4] > 120 & restr[i].dih[4]<=250){
				idx=1;
			}else{
				idx=3;
			}
		}
		
		randval=Random();
		//if(randval > 0.98){
			//randval=Random();
			//if(randval < 0.25) ntparamidx.push_back(0);
			//else if(randval >= 0.25 & randval < 0.5) ntparamidx.push_back(1);
			//else if(randval >= 0.5 & randval < 0.75) ntparamidx.push_back(2);
			//else if(randval >= 0.75 & randval < 1.0) ntparamidx.push_back(3);
		//}else{
		ntparamidx.push_back(idx);
		//}
	}

	for(j=trandposst;j<=trandposed;j++)
	{
		adjparamidx.push_back(0);
	}
}


bool RNASimulation::fmovementlmpfull(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc)
{
	int i,j,k,inext;
	int istart,iend;
	int trandpos, trandposed, trandposst, span;
	double trandang;
	Point3d rp;
	double mati[9];
	int wpos;
	char aat1, aat2;
	bool idxtmp1, idxtmp2;
	vector<int> scpos(numseq,0);
	vector<int> accsc(numseq,0);
	vector<int> accscr(numseq,0);
	int transbase;
	int total;
	int numiter;
	double rmax, rnorm;
	int mtype1, mtype2;
	vector<int> ntparamidx;
	vector<int> adjparamidx;
	bool flagdone;
	double randval=Random();
	bool exitOuterLoop = false;

	int totwrong=lp1.nn[CLASH_ENERGY_ST]*2;

	if(totwrong>0 && randval<0.7){
		wpos=int(Random()*totwrong);

		i=0;
		for(trandpos=0;trandpos<numseq;trandpos++)
		{
			inext=i+lp1.indn[CLASH_ENERGY_ST][trandpos];
			if(i<=wpos && wpos<inext){
				break;
			}
			i=inext;
		}
	}else{
		trandpos=int((numseq)*Random());
	}

	span=2+int(6*Random());
	trandposst=trandpos-span/2;
	trandposed=trandposst+span;
	if(trandposst<1) trandposst=1;
	if(trandposed>numseq-1) trandposed=numseq-1;
	span=trandposed-trandposst;
	if(span < 2){
		newenergy=oldenergy;
		return false;
	}
	
	copystr(decstr, numseq);

	setlmpparam(restr, trandposst, trandposed, ntparamidx, adjparamidx);

	int threshiter=span*200;
	if(threshiter>3000) threshiter=3000;

	for(i=trandposst;i<trandposed;i++)
	{
		for(j=0;j<PAIRLENGTHTARGETMOVNUM;j++){
			//rmax=0.25+0.25*Random();
			rmax=0.5+0.5*Random();
			rp=bf.ranv(rmax);
			restr[i].ptv[j].x+=rp.x;
			restr[i].ptv[j].y+=rp.y;
			restr[i].ptv[j].z+=rp.z;
		}
	}

	numiter=0;
	do{	
		for(i=trandposst;i<trandposed;i++)
		{
			for(j=0;j<100;j++){
				singlemoveLMPfull(restr,i,0,ntparamidx[i-trandposst]);
			}
		}

		for(i=trandposed;i>trandposst;i--)
		{
			for(j=0;j<100;j++){
				singlemoveLMPfull(restr,i,1,ntparamidx[i-trandposst-1]);
			}
		}
		//cout<<"-----------"<<endl;

		flagdone=true;
		exitOuterLoop = false;
		for(k=trandposst;k<=trandposed;k++)
		{
			for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
			{
				for(mtype2=0;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
				{
					aat1=restr[k-1].aaa;
					aat2=restr[k].aaa;

					rp=bf.setv(restr[k].ptv[mtype2].x-restr[k-1].ptv[mtype1].x,restr[k].ptv[mtype2].y-restr[k-1].ptv[mtype1].y,restr[k].ptv[mtype2].z-restr[k-1].ptv[mtype1].z);
					rnorm=bf.norm(rp);

					//cout<<to_string(k) + "-" + to_string(mtype1) + "-" + to_string(mtype2) + "-" + to_string(rnorm) + "-" + to_string(PAIRLENGADJPARAM[aat1][aat2][mtype1][mtype2].lower) + "-" + to_string(PAIRLENGADJPARAM[aat1][aat2][mtype1][mtype2].upper)<<endl;

					if(rnorm<PAIRLENGADJPARAM[aat1][aat2][0][mtype1][mtype2].lower || rnorm>PAIRLENGADJPARAM[aat1][aat2][0][mtype1][mtype2].upper)
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
				for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
				{
					for(mtype2=mtype1+1;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
					{
						aat1=restr[k].aaa;
						aat2=restr[k].aaa;

						rp=bf.setv(restr[k].ptv[mtype2].x-restr[k].ptv[mtype1].x,restr[k].ptv[mtype2].y-restr[k].ptv[mtype1].y,restr[k].ptv[mtype2].z-restr[k].ptv[mtype1].z);
						rnorm=bf.norm(rp);

						//cout<<to_string(k) + "-" + to_string(mtype1) + "-" + to_string(mtype2) + "-" + to_string(rnorm) + "-" + to_string([aat1][mtype1][mtype2].lower) + "-" + to_string(ALLPAIRLENGNTPARAM[aat1][mtype1][mtype2].upper)<<endl;

						if(rnorm<PAIRLENGNTPARAM[aat1][ntparamidx[k-trandposst]][mtype1][mtype2].lower || rnorm>PAIRLENGNTPARAM[aat1][ntparamidx[k-trandposst]][mtype1][mtype2].upper)
						//if(rnorm<ALLPAIRLENGNTPARAM[aat1][mtype1][mtype2].lower || rnorm>ALLPAIRLENGNTPARAM[aat1][mtype1][mtype2].upper)
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
			
	istart=trandposst-1;
	iend=trandposed+2;
	if(iend>=numseq) iend=numseq;

	coord2parambb(restr,trandposst,iend,1);
	coord2paramsugar(restr,trandposst,iend,1);
	param2coordbase(restr,trandposst,iend);
	//param2coordsugar(restr,istart,numseq);
	param2coordhydrogen(restr,trandposst,numseq);

	coord2paramadd(restr,istart,iend);

	settmplmpflap(numseq, istart, iend);
	newenergy=fcalcrmsdenergy(restr, numseq, lp2, rppmatbk, bpenergybk, bsenergybk, nmapbk, distmapbk, ntpairdistbk, adjpairdistbk, distenebk, multidistenebk, contactenebk, mapflagbk, mapbbflagbk, flagacc);
	/*
	cout<<"QqqqqqqqqqqqqqqqqQ"<<endl;
	checkupdate(restr, numseq);
	cout<<"PppppppppppppppppP"<<endl;
	fsetinitial(lp2,rppmatbk, bpenergybk, bsenergybk, nmapbk, distmapbk, ntpairdistbk, adjpairdistbk, distenebk, multidistenebk, contactenebk, mapflagbk, mapbbflagbk, numseq);
	double newenergy1=fcalcrmsdenergy(restr, numseq, lp2, rppmatbk, bpenergybk, bsenergybk, nmapbk, distmapbk, ntpairdistbk, adjpairdistbk, distenebk, multidistenebk, contactenebk, mapflagbk, mapbbflagbk, false);
	if(fabs(newenergy1-newenergy) > 10){
		cout<<"Error"<<endl;
		cout<<"2:: " + to_string(newenergy1)<<endl;
		cout<<"3:: " + to_string(newenergy)<<endl;
		exit(1);
	}
	cout<<"OK"<<endl;
	*/
	
	bool upidx=metropolisUpdate(oldenergy, newenergy, Bt);

	if(upidx)
	{
		rcopystr(decstr, restr, numseq);
		return true;
	}
	else
	{
		return false;
	}
}












