// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


void RNASimulation::setlenflap(double numseq, int trandpos, double aat, int lentype)
{
	int i,j;

	lp2.flap[BLENGTH_ENERGY_ST+lentype][trandpos]=true;

	for(j=0;j<ADDBLENGTHTARGET.size();j++)
	{
		lp2.flap[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]][trandpos]=true;
	}

	if(lentype == 5 && trandpos < numseq-1){
		lp2.flap[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[0]][trandpos+1]=true;
	}

	for(j=0;j<ADDINNERANGLETARGET.size();j++)
	{
		lp2.flap[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]][trandpos]=true;
	}

	for(i=0;i<numseq;i++)
	{
		for(j=i;j<numseq;j++)
		//for(j=i+1;j<numseq;j++)
		{
			if(j<trandpos || i>trandpos+1){}
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


void RNASimulation::setlenntflap(double numseq, int trandpos, double aat, int lentype)
{
	int i,j;

	lp2.flap[BLENGTH_ENERGY_ST+lentype][trandpos]=true;

	for(j=0;j<ADDBLENGTHTARGET.size();j++)
	{
		lp2.flap[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]][trandpos]=true;
	}

	if(lentype == 5 && trandpos < numseq-1){
		lp2.flap[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[0]][trandpos+1]=true;
	}

	for(j=0;j<ADDINNERANGLETARGET.size();j++)
	{
		lp2.flap[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]][trandpos]=true;
	}

	for(i=0;i<numseq;i++)
	{
		mapflagbk[i][trandpos]=true;
		mapflagbk[trandpos][i]=true;
	}

	/*
	for(i=0;i<numseq;i++)
	{
		for(j=i;j<numseq;j++)
		//for(j=i+1;j<numseq;j++)
		{
			if(j<trandpos || i>trandpos+1){}
			else
			{
				mapflagbk[i][j]=true;
				mapflagbk[j][i]=true;
			}
		}
	}
	*/
}

double RNASimulation::getmovlenval(double val, double mean, double std, double k){
	double newval;
	double deg;
	if(val <= mean){
		deg=mean-val;
		if(deg > std) deg = std;
		newval = val + (deg)*k*Random();
	}else{
		deg=val-mean;
		if(deg > std) deg = std;
		newval = val - (deg)*k*Random();
	}

	return newval;
}

bool RNASimulation::fmovementlen(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc)//length 
{

	int i,j;
	double r,l;
	double nr,nl;
	int trandpos;
	bool flagtor;
	int wpos;
	char aat;
	double pval;
	bool upidx;

	int lentype=BLENGTHTARGET[int(Random()*BLENGTHTARGET.size())];
	int totwrong=lp1.nn[BLENGTH_ENERGY_ST+lentype];
	int wrong=0;
	double randval=Random();
	int uptype;
	//double Btang;

	//cout<<to_string(0) + "-" + to_string(basepairs[0].basepidx.first) + "-" + to_string(basepairs[0].basepidx.second)<<endl;


	if(totwrong>0 && randval<0.4+0.3*(totwrong/numseq)){
		uptype=0;
		wpos=int(Random()*totwrong);

		i=-1;
		for(trandpos=0;trandpos<numseq;trandpos++)
		{
			if(lp1.indn[BLENGTH_ENERGY_ST+lentype][trandpos]==1) i++;
			if(i==wpos) break;
		}
	}
	else{
		uptype=1;
		trandpos=int((numseq)*Random());
	}

	copystr(decstr, numseq);
	aat=restr[trandpos].aaa;
	
	if(trandpos==0 && lentype == 10){
		if(uptype == 0){
			restr[trandpos].len[lentype]=getmovlenval(restr[trandpos].len[lentype], OPLENGTHPARAM[aat][0].mean, OPLENGTHPARAM[aat][0].std, 0.1);
			//Btang=(1.4+0.02*numseq);
		}else{
			restr[trandpos].len[lentype]=getrandomnewparam(OPLENGTHPARAM[aat][0].lower, OPLENGTHPARAM[aat][0].upper);
			//Btang=-1;
		}
	}else if(trandpos==0 && lentype == 11){
		if(uptype == 0){
			restr[trandpos].len[lentype]=getmovlenval(restr[trandpos].len[lentype], OPLENGTHPARAM[aat][1].mean, OPLENGTHPARAM[aat][1].std, 0.1);
			//Btang=(1.4+0.02*numseq);
		}else{
			restr[trandpos].len[lentype]=getrandomnewparam(OPLENGTHPARAM[aat][1].lower, OPLENGTHPARAM[aat][1].upper);
			//Btang=-1;
		}
	}else{
		if(uptype == 0){
			restr[trandpos].len[lentype]=getmovlenval(restr[trandpos].len[lentype], LENGTHPARAM[aat][lentype].mean, LENGTHPARAM[aat][lentype].std, 0.1);
			//Btang=(1.4+0.02*numseq);
		}else{
			restr[trandpos].len[lentype]=getrandomnewparam(LENGTHPARAM[aat][lentype].lower, LENGTHPARAM[aat][lentype].upper);
			//Btang=-1;
		}
	}
	
	if(lentype < 6){
		flagtor=param2coordbb(restr,trandpos,numseq);
		param2coordsugar(restr,trandpos,numseq);
		param2coordbase(restr,trandpos,numseq);
		param2coordhydrogen(restr,trandpos,numseq);
		if(!flagtor) printf("tor2str wrong in len %d\n",trandpos);

		if(lentype == 5){
			int end = min(trandpos+2, numseq);
			coord2paramadd(restr,trandpos,end);
		}else{
			coord2paramadd(restr,trandpos,trandpos+1);
		}

		setlenflap(numseq, trandpos, aat, lentype);
	}else{
		param2coordsugar(restr,trandpos,trandpos+1);
		param2coordbase(restr,trandpos,trandpos+1);
		param2coordhydrogen(restr,trandpos,trandpos+1);
		coord2paramadd(restr,trandpos,trandpos+1);
		setlenntflap(numseq, trandpos, aat, lentype);
	}

	newenergy=fcalcrmsdenergy(restr, numseq, lp2, rppmatbk, bpenergybk, bsenergybk, nmapbk, distmapbk, ntpairdistbk, adjpairdistbk, distenebk, multidistenebk, contactenebk, mapflagbk, mapbbflagbk, flagacc);
	upidx=metropolisUpdate(oldenergy, newenergy, Bt);
	//upidx=metropolisUpdate(oldenergy, newenergy, Btang);

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

	if(upidx){
		rcopystr(decstr, restr, numseq);
		return true;
	}
	else
	{
		return false;
	}
}












