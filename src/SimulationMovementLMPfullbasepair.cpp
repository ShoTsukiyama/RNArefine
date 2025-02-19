// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

void RNASimulation::singbplemoveLMPfull(vector<Residueatoms> &restr,int k,int rtype,int flagdir,int paramidx)
{
	//BasicFunc bf;
	Point3d rp;
	double rnorm;
	char aat1, aat2;
	int tmp;

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

		if(rnorm<BPATOMPAIRLENGADJPARAM[aat1][aat2][paramidx][mtype1][mtype2].lower || rnorm>BPATOMPAIRLENGADJPARAM[aat1][aat2][paramidx][mtype1][mtype2].upper)
		{
			rp=bf.scal(rp,BPATOMPAIRLENGADJPARAM[aat1][aat2][paramidx][mtype1][mtype2].mean/rnorm);
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

		if(rnorm<BPATOMPAIRLENGNTPARAM[aat1][paramidx][mtype1][mtype2].lower || rnorm>BPATOMPAIRLENGNTPARAM[aat1][paramidx][mtype1][mtype2].upper)
		{
			//cout<<rnorm<<endl;
			//cout<<to_string(rp.x) + "-" + to_string(rp.y) + "-" + to_string(rp.z)<<endl;
			rp=bf.scal(rp,BPATOMPAIRLENGNTPARAM[aat1][paramidx][mtype1][mtype2].mean/rnorm);
			//cout<<to_string(rp.x) + "-" + to_string(rp.y) + "-" + to_string(rp.z)<<endl;
			restr[k].ptv[mtype2].x=restr[k].ptv[mtype1].x+rp.x;
			restr[k].ptv[mtype2].y=restr[k].ptv[mtype1].y+rp.y;
			restr[k].ptv[mtype2].z=restr[k].ptv[mtype1].z+rp.z;
		}
	}else if(rtype == 0 & flagdir==1){
		int mtype1 = int(Random()*PAIRLENGTHTARGETMOVNUM);
		int mtype2 = int(Random()*PAIRLENGTHTARGETMOVNUM);

		aat1=restr[k].aaa;
		aat2=restr[k+1].aaa;

		rp=bf.setv(restr[k].ptv[mtype1].x-restr[k+1].ptv[mtype2].x, restr[k].ptv[mtype1].y-restr[k+1].ptv[mtype2].y, restr[k].ptv[mtype1].z-restr[k+1].ptv[mtype2].z);
		rnorm=bf.norm(rp);
		
		if(rnorm<EPS)//two in the same place
		{
			restr[k].ptv[mtype1].x+=0.0001f; restr[k].ptv[mtype1].y+=0.0001f; restr[k].ptv[mtype1].z+=0.0001f;
			rp=bf.setv(restr[k].ptv[mtype1].x-restr[k+1].ptv[mtype2].x, restr[k].ptv[mtype1].y-restr[k+1].ptv[mtype2].y, restr[k].ptv[mtype1].z-restr[k+1].ptv[mtype2].z);
			rnorm=bf.norm(rp);
		}

		if(rnorm<BPATOMPAIRLENGADJPARAM[aat2][aat1][paramidx][mtype1][mtype2].lower || rnorm>BPATOMPAIRLENGADJPARAM[aat2][aat1][paramidx][mtype1][mtype2].upper)
		{
			//cout<<rnorm<<endl;
			//cout<<to_string(rp.x) + "-" + to_string(rp.y) + "-" + to_string(rp.z)<<endl;
			rp=bf.scal(rp,BPATOMPAIRLENGADJPARAM[aat2][aat1][paramidx][mtype1][mtype2].mean/rnorm);
			//cout<<to_string(rp.x) + "-" + to_string(rp.y) + "-" + to_string(rp.z)<<endl;
			restr[k].ptv[mtype1].x=restr[k+1].ptv[mtype2].x+rp.x;
			restr[k].ptv[mtype1].y=restr[k+1].ptv[mtype2].y+rp.y;
			restr[k].ptv[mtype1].z=restr[k+1].ptv[mtype2].z+rp.z;
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

		if(rnorm<BPATOMPAIRLENGNTPARAM[aat1][paramidx][mtype1][mtype2].lower || rnorm>BPATOMPAIRLENGNTPARAM[aat1][paramidx][mtype1][mtype2].upper)
		{
			//cout<<rnorm<<endl;
			//cout<<to_string(rp.x) + "-" + to_string(rp.y) + "-" + to_string(rp.z)<<endl;
			rp=bf.scal(rp,BPATOMPAIRLENGNTPARAM[aat1][paramidx][mtype1][mtype2].mean/rnorm);
			//cout<<to_string(rp.x) + "-" + to_string(rp.y) + "-" + to_string(rp.z)<<endl;
			restr[k].ptv[mtype2].x=restr[k].ptv[mtype1].x+rp.x;
			restr[k].ptv[mtype2].y=restr[k].ptv[mtype1].y+rp.y;
			restr[k].ptv[mtype2].z=restr[k].ptv[mtype1].z+rp.z;
		}
	}
}

void RNASimulation::singbplemoveLMPfullbp(vector<Residueatoms> &restr,int i, int j,int rtype)
{
	Point3d rp;
	double rnorm;
	char aat1, aat2;
	int tmp;

	int mtype1 = int(Random()*PAIRLENGTHTARGETMOVNUM);
	int mtype2 = int(Random()*PAIRLENGTHTARGETMOVNUM);

	int basei, basej;

	if(rtype == 0){
		basei=i;
		basej=j;
	}else{
		basej=i;
		basei=j;
	}

	aat1=restr[basej].aaa;
	aat2=restr[basei].aaa;

	rp=bf.setv(restr[basei].ptv[mtype2].x-restr[basej].ptv[mtype1].x,restr[basei].ptv[mtype2].y-restr[basej].ptv[mtype1].y,restr[basei].ptv[mtype2].z-restr[basej].ptv[mtype1].z);
	rnorm=bf.norm(rp);

	if(rnorm<EPS)//two in the same place
	{
		restr[basei].ptv[mtype2].x+=0.0001f; restr[basei].ptv[mtype2].y+=0.0001f; restr[basei].ptv[mtype2].z+=0.0001f;
		rp=bf.setv(restr[basei].ptv[mtype2].x-restr[basej].ptv[mtype1].x,restr[basei].ptv[mtype2].y-restr[basej].ptv[mtype1].y,restr[basei].ptv[mtype2].z-restr[basej].ptv[mtype1].z);
		rnorm=bf.norm(rp);
	}

	//exit(1);

	if(rnorm<BPATOMPAIRLENGPARAM[aat1][aat2][mtype1][mtype2].lower || rnorm>BPATOMPAIRLENGPARAM[aat1][aat2][mtype1][mtype2].upper)
	{
		rp=bf.scal(rp,BPATOMPAIRLENGPARAM[aat1][aat2][mtype1][mtype2].mean/rnorm);
		restr[basei].ptv[mtype2].x=restr[basej].ptv[mtype1].x+rp.x;
		restr[basei].ptv[mtype2].y=restr[basej].ptv[mtype1].y+rp.y;
		restr[basei].ptv[mtype2].z=restr[basej].ptv[mtype1].z+rp.z;
	}
	
}

void RNASimulation::settmpbplmpfullflap(double numseq, double trandpos1, double trandpos2)
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

bool RNASimulation::checkadjlength(vector<Residueatoms> restr, int basei,int paramidx)
{
	int mtype1, mtype2;
	char aat1, aat2;
	bool flagdone=true;
	bool exitOuterLoop = false;
	Point3d rp;
	double rnorm;

	for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
	{
		for(mtype2=0;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
		{
			aat1=restr[basei-1].aaa;
			aat2=restr[basei].aaa;

			rp=bf.setv(restr[basei].ptv[mtype2].x-restr[basei-1].ptv[mtype1].x,restr[basei].ptv[mtype2].y-restr[basei-1].ptv[mtype1].y,restr[basei].ptv[mtype2].z-restr[basei-1].ptv[mtype1].z);
			rnorm=bf.norm(rp);

			if(rnorm<BPATOMPAIRLENGADJPARAM[aat1][aat2][paramidx][mtype1][mtype2].lower || rnorm>BPATOMPAIRLENGADJPARAM[aat1][aat2][paramidx][mtype1][mtype2].upper)
			{
				flagdone=false;
				exitOuterLoop=true;
				break;
			}
		}if(exitOuterLoop) break;
	}

	return flagdone;
}

bool RNASimulation::checkntlength(vector<Residueatoms> restr, int basei,int paramidx)
{
	int mtype1, mtype2;
	char aat1, aat2;
	bool flagdone=true;
	bool exitOuterLoop = false;
	Point3d rp;
	double rnorm;

	for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
	{
		for(mtype2=mtype1+1;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
		{
			aat1=restr[basei].aaa;
			aat2=restr[basei].aaa;

			rp=bf.setv(restr[basei].ptv[mtype2].x-restr[basei].ptv[mtype1].x,restr[basei].ptv[mtype2].y-restr[basei].ptv[mtype1].y,restr[basei].ptv[mtype2].z-restr[basei].ptv[mtype1].z);
			rnorm=bf.norm(rp);

			if(rnorm<BPATOMPAIRLENGNTPARAM[aat1][paramidx][mtype1][mtype2].lower || rnorm>BPATOMPAIRLENGNTPARAM[aat1][paramidx][mtype1][mtype2].upper)
			{
				flagdone=false;
				exitOuterLoop=true;
				break;
			}
		}if(exitOuterLoop) break;
	}

	return flagdone;
}

bool RNASimulation::checkbplength(vector<Residueatoms> restr, int basei, int basej)
{
	int mtype1, mtype2;
	char aat1, aat2;
	bool flagdone=true;
	bool exitOuterLoop = false;
	Point3d rp;
	double rnorm;

	for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
	{
		for(mtype2=0;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
		{
			aat1=restr[basej].aaa;
			aat2=restr[basei].aaa;

			rp=bf.setv(restr[basei].ptv[mtype2].x-restr[basej].ptv[mtype1].x,restr[basei].ptv[mtype2].y-restr[basej].ptv[mtype1].y,restr[basei].ptv[mtype2].z-restr[basej].ptv[mtype1].z);
			rnorm=bf.norm(rp);

			if(rnorm<BPATOMPAIRLENGPARAM[aat1][aat2][mtype1][mtype2].lower || rnorm>BPATOMPAIRLENGPARAM[aat1][aat2][mtype1][mtype2].upper)
			{
				flagdone=false;
				exitOuterLoop=true;
				break;
			}
		}if(exitOuterLoop) break;
	}

	return flagdone;
}



bool RNASimulation::fmovementlmpfullbasepair(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc)
{
	int i,j,k,inext;
	int istart,iend;
	int trandpos, trandposed, trandposst, span;
	int basei, basej;
	double trandang;
	Point3d rp;
	double mati[9];
	int wpos;
	char aat1, aat2;
	bool idxtmp1, idxtmp2;
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
	int mtype;
	int rtype;
	int pairwiseidx;

	int totwrong1=lp1.nn[CLASH_ENERGY_ST];
	int totwrong2=lp1.nn[BASE_ENERGY_ST+0] + lp1.nn[BASE_ENERGY_ST+1];

	if(totwrong1>0 && randval<0.4){
		wpos=int(Random()*totwrong1);

		i=0;
		for(trandpos=0;trandpos<numseq;trandpos++)
		{
			inext=i+lp1.indn[CLASH_ENERGY_ST][trandpos];
			if(i<=wpos && wpos<inext){
				break;
			}
			i=inext;
		}
	}else if(totwrong2>0 && randval>0.4 && randval<0.8){
		wpos=int(Random()*totwrong2);
		mtype=1;

		i=0;
		for(trandpos=0;trandpos<BASEPAIRS.size();trandpos++)
		{
			basei=BASEPAIRS[trandpos].basepidx.first;
			basej=BASEPAIRS[trandpos].basepidx.second;

			pairwiseidx=bpmapid2idx[basei*numseq+basej];

			inext=i+bpenergy[pairwiseidx].pcw[0]+bpenergy[pairwiseidx].pcw[1];
			if(i<=wpos && wpos<inext){
				break;
			}
			i=inext;
		}

		randval=Random();
		if(randval > 0.5) trandpos=basei;
		else trandpos=basej;
	}
	
	else{
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

	//cout<<to_string(trandposst) + "-" + to_string(trandposed)<<endl;

	vector<pair<int,int>> pairsf; 
	for(i=trandposst;i<trandposed;i++){
		if(bppairednt[i] == -1){
			pairsf.push_back({i,-1});
		}else if(bppairednt[i] < 1){
			break;
		}else if(bppairednt[i] > numseq-1){
			break;
		}else{
			pairsf.push_back({i,bppairednt[i]});
		}
	}

	vector<pair<int,int>> pairsb;
	for(i=trandposed-1;i>=trandposst;i--){
		if(bppairednt[i] == -1){
			pairsb.push_back({i,-1});
		}else if(bppairednt[i] < 1){
			break;
		}else if(bppairednt[i] > numseq-1){
			break;
		}else{
			pairsb.push_back({i,bppairednt[i]});
		}
	}

	vector<pair<int,int>> pairs;
	if(pairsf.size() >= pairsb.size()){
		pairs=pairsf;
	}else{
		for (i=pairsb.size()-1; i>=0; i--) {
			pairs.push_back(pairsb[i]);
		}
	}

	copystr(decstr, numseq);

	vector<vector<int>> lmpparamidx;
	setlmpbpparam(restr, pairs, lmpparamidx);

	int threshiter=span*200;
	if(threshiter>3000) threshiter=3000;

	for(i=0;i<pairs.size();i++)
	{
		basei = pairs[i].first;
		basej = pairs[i].second;

		if(basej == -1){
			for(j=0;j<PAIRLENGTHTARGETMOVNUM;j++){
				//cout<<to_string(restr[basei].ptv[j].x) + "-" + to_string(restr[basei].ptv[j].y) + "-" + to_string(restr[basei].ptv[j].z)<<endl;
				rmax=0.25+0.25*Random();
				//rmax=0.5+0.5*Random();
				rp=bf.ranv(rmax);
				restr[basei].ptv[j].x+=rp.x;
				restr[basei].ptv[j].y+=rp.y;
				restr[basei].ptv[j].z+=rp.z;
				//cout<<to_string(restr[basei].ptv[j].x) + "-" + to_string(restr[basei].ptv[j].y) + "-" + to_string(restr[basei].ptv[j].z)<<endl;
			}
		}else{
			for(j=0;j<PAIRLENGTHTARGETMOVNUM;j++){
				//cout<<to_string(restr[basei].ptv[j].x) + "-" + to_string(restr[basei].ptv[j].y) + "-" + to_string(restr[basei].ptv[j].z)<<endl;
				rmax=0.25+0.25*Random();
				//rmax=0.5+0.5*Random();
				rp=bf.ranv(rmax);
				restr[basei].ptv[j].x+=rp.x;
				restr[basei].ptv[j].y+=rp.y;
				restr[basei].ptv[j].z+=rp.z;
				//cout<<to_string(restr[basei].ptv[j].x) + "-" + to_string(restr[basei].ptv[j].y) + "-" + to_string(restr[basei].ptv[j].z)<<endl;
			}
			for(j=0;j<PAIRLENGTHTARGETMOVNUM;j++){
				//cout<<to_string(restr[basej].ptv[j].x) + "-" + to_string(restr[basej].ptv[j].y) + "-" + to_string(restr[basej].ptv[j].z)<<endl;
				
				//rmax=0.25+0.25*Random();
				rmax=0.5+0.5*Random();
				rp=bf.ranv(rmax);
				restr[basej].ptv[j].x+=rp.x;
				restr[basej].ptv[j].y+=rp.y;
				restr[basej].ptv[j].z+=rp.z;
				//cout<<to_string(restr[basej].ptv[j].x) + "-" + to_string(restr[basej].ptv[j].y) + "-" + to_string(restr[basej].ptv[j].z)<<endl;
				
			}
		}
	}

	vector<bool> checkidxtmp(numseq, false);

	numiter=0;
	do{
		for(i=0;i<pairs.size();i++)
		{
			basei = pairs[i].first;
			basej = pairs[i].second;

			//cout<<to_string(basei) + "-" + to_string(basej)<<endl;
			
			if(basej == -1){
				for(j=0;j<100;j++){
					rtype = int(Random()*2.0);
					if(rtype == 0){
						singbplemoveLMPfull(restr,basei,rtype,0,lmpparamidx[i][1]);
					}else{
						singbplemoveLMPfull(restr,basei,rtype,0,lmpparamidx[i][0]);
					}
				}
			}else{
				for(j=0;j<100;j++){
					rtype = int(Random()*2.0);
					singbplemoveLMPfullbp(restr,basei,basej,rtype);
				}

				for(j=0;j<100;j++){
					randval=Random();
					if(randval < 0.5){
						rtype = int(Random()*2.0);
						if(rtype == 0){
							singbplemoveLMPfull(restr,basei,rtype,0,lmpparamidx[i][1]);
						}else{
							singbplemoveLMPfull(restr,basei,rtype,0,lmpparamidx[i][0]);
						}
					}else{
						rtype = int(Random()*2.0);
						if(rtype == 0){
							singbplemoveLMPfull(restr,basej,rtype,0,lmpparamidx[i][4]);
						}else if(rtype == 1){
							singbplemoveLMPfull(restr,basej,rtype,0,lmpparamidx[i][3]);
						}
					}
				}
			}
		}

		for(i=pairs.size()-1;i>=0;i--)
		{
			basei = pairs[i].first;
			basej = pairs[i].second;
			if(basej == -1){
				for(j=0;j<100;j++){
					rtype = int(Random()*2.0);
					if(rtype == 0){
						singbplemoveLMPfull(restr,basei,rtype,1,lmpparamidx[i][2]);
					}else{
						singbplemoveLMPfull(restr,basei,rtype,1,lmpparamidx[i][0]);
					}
				}
			}else{
				for(j=0;j<100;j++){
					randval=Random();
					if(randval < 0.5){
						rtype = int(Random()*2.0);
						if(rtype == 0){
							singbplemoveLMPfull(restr,basei,rtype,1,lmpparamidx[i][2]);
						}else if(rtype == 1){
							singbplemoveLMPfull(restr,basei,rtype,1,lmpparamidx[i][0]);
						}
					}else{
						rtype = int(Random()*2.0);
						if(rtype == 0){
							singbplemoveLMPfull(restr,basej,rtype,1,lmpparamidx[i][5]);
						}else if(rtype == 1){
							singbplemoveLMPfull(restr,basej,rtype,1,lmpparamidx[i][3]);
						}
					}
				}

				for(j=0;j<100;j++){
					rtype = int(Random()*2.0);
					singbplemoveLMPfullbp(restr,basei,basej,rtype);
				}
			}
		}

		flagdone=true;
		exitOuterLoop = false;
		vector<bool> checkidx=checkidxtmp;

		for(i=0;i<pairs.size();i++)
		{
			basei = pairs[i].first;
			basej = pairs[i].second;

			if(basej == -1){
				if(!checkidx[basei]){
					flagdone=checkadjlength(restr, basei,lmpparamidx[i][1]);
					checkidx[basei] = true;
					if(!flagdone) break;
				}

				flagdone=checkntlength(restr, basei, lmpparamidx[i][0]);
				if(!flagdone) break;

				if(!checkidx[basei+1]){
					flagdone=checkadjlength(restr, basei+1,lmpparamidx[i][2]);
					checkidx[basei+1] = true;
					if(!flagdone) break;
				}
			}else{
				flagdone=checkbplength(restr, basei, basej);
				if(!flagdone) break;

				if(!checkidx[basei]){
					flagdone=checkadjlength(restr, basei,lmpparamidx[i][1]);
					checkidx[basei] = true;
					if(!flagdone) break;
				}

				flagdone=checkntlength(restr, basei, lmpparamidx[i][0]);
				if(!flagdone) break;

				if(!checkidx[basei+1]){
					flagdone=checkadjlength(restr, basei+1,lmpparamidx[i][2]);
					checkidx[basei+1] = true;
					if(!flagdone) break;
				}

				if(!checkidx[basej]){
					flagdone=checkadjlength(restr, basej,lmpparamidx[i][4]);
					checkidx[basej] = true;
					if(!flagdone) break;
				}

				flagdone=checkntlength(restr, basej, lmpparamidx[i][3]);
				if(!flagdone) break;

				if(!checkidx[basej+1]){
					flagdone=checkadjlength(restr, basej+1,lmpparamidx[i][5]);
					checkidx[basej+1] = true;
					if(!flagdone) break;
				}
			}
		}
		
		numiter++;
	}while(!flagdone && numiter<threshiter);
	
	vector<bool> targetidxtmp(numseq, false);
	for(i=0;i<pairs.size();i++)
	{
		basei = pairs[i].first;
		basej = pairs[i].second;

		targetidxtmp[basei] = 1;
		if(basej != -1){
			targetidxtmp[basej] = 1;
		}
	}
	vector<pair<int,int>> continuousidx;
	int start = -1;

	for(i=0;i<numseq;i++){
		if (targetidxtmp[i] == 1) {
			if (start == -1) {
				start = i;
			}
		}else {
            if (start != -1) {
                continuousidx.push_back({start, i - 1});
                start = -1;
            }
        }
	}
	if (start != -1) {
        continuousidx.push_back({start, numseq - 1});
    }

	for(i = 0; i<continuousidx.size();i++){
		istart=continuousidx[i].first-1;
		iend=continuousidx[i].second+2;
		if(iend>=numseq) iend=numseq;
		//cout<<to_string(istart) + "-" + to_string(iend)<<endl;
		coord2parambb(restr,continuousidx[i].first,iend,1);
		coord2paramsugar(restr,continuousidx[i].first,iend,1);
		param2coordbase(restr,continuousidx[i].first,iend);
		param2coordhydrogen(restr,continuousidx[i].first,iend);

		coord2paramadd(restr,istart,iend);
		settmplmpflap(numseq, istart, iend);
	}

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












