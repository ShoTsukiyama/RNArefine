// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////



void RNASimulation::setangflap(double numseq, double trandpos, double aat, int angtype)
{
	int i,j;

	lp2.flap[INNERANGLE_ENERGY_ST+angtype][trandpos]=true;

	for(j=0;j<ADDBLENGTHTARGET.size();j++)
	{
		lp2.flap[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]][trandpos]=true;
	}

	if(angtype == 5 && trandpos < numseq-1){
		lp2.flap[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[0]][trandpos+1]=true;
	}

	for(j=0;j<ADDINNERANGLETARGET.size();j++)
	{
		lp2.flap[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]][trandpos]=true;
	}

	for(i=0;i<numseq;i++)
	{
		for(j=i;j<numseq;j++)
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

void RNASimulation::setangntflap(double numseq, double trandpos, double aat, int angtype)
{
	int i,j;

	lp2.flap[INNERANGLE_ENERGY_ST+angtype][trandpos]=true;

	for(j=0;j<ADDBLENGTHTARGET.size();j++)
	{
		lp2.flap[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]][trandpos]=true;
	}

	if(angtype == 5 && trandpos < numseq-1){
		lp2.flap[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[0]][trandpos+1]=true;
	}

	for(j=0;j<ADDINNERANGLETARGET.size();j++)
	{
		lp2.flap[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]][trandpos]=true;
	}

	//for(i=0;i<numseq;i++)
	//{
		//mapflagbk[i][trandpos]=true;
		//mapflagbk[trandpos][i]=true;
	//}

	
	for(i=0;i<numseq;i++)
	{
		for(j=i;j<numseq;j++)
		{
			if(j<trandpos || i>trandpos+1){}
			else
			{
				mapflagbk[i][j]=true;
				mapflagbk[j][i]=true;
			}
		}
	}
	
}

double RNASimulation::getmovangval(double val, double mean, double std, double k){
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

bool RNASimulation::fmovementang(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc)//angle 
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

	int angtype=INNERANGLETARGET[int(Random()*INNERANGLETARGET.size())];
	int totwrong=lp1.nn[INNERANGLE_ENERGY_ST+angtype];
	int wrong=0;
	double randval=Random();
	int uptype;
	//double Btang;

	if(totwrong>0 && randval<0.4+0.3*(totwrong/numseq)){
		uptype=0;
		wpos=int(Random()*totwrong);

		i=-1;
		for(trandpos=0;trandpos<numseq;trandpos++)
		{
			if(lp1.indn[INNERANGLE_ENERGY_ST+angtype][trandpos]>0) i++;
			if(i==wpos) break;
		}
	}
	else{
		uptype=1;
		trandpos=int((numseq)*Random());
	}

	copystr(decstr, numseq);
	aat=restr[trandpos].aaa;

	if(trandpos==0 && angtype == 10){
		if(uptype == 0){
			restr[trandpos].ang[angtype]=getmovangval(restr[trandpos].ang[angtype], OPINNERANGPARAM[aat][0].mean, OPINNERANGPARAM[aat][0].std, 0.1);
			//Btang=(1.4+0.02*numseq);
		}else{
			restr[trandpos].ang[angtype]=getrandomnewparam(OPINNERANGPARAM[aat][0].lower, OPINNERANGPARAM[aat][0].upper);
			//Btang=-1;
		}
	}else if(trandpos==0 && angtype == 11){
		if(uptype == 0){
			restr[trandpos].ang[angtype]=getmovangval(restr[trandpos].ang[angtype], OPINNERANGPARAM[aat][1].mean, OPINNERANGPARAM[aat][1].std, 0.1);
			//Btang=(1.4+0.02*numseq);
		}else{
			restr[trandpos].ang[angtype]=getrandomnewparam(OPINNERANGPARAM[aat][1].lower, OPINNERANGPARAM[aat][1].upper);
			//Btang=-1;
		}
	}else{
		if(uptype == 0){
			restr[trandpos].ang[angtype]=getmovangval(restr[trandpos].ang[angtype], INNERANGPARAM[aat][angtype].mean, INNERANGPARAM[aat][angtype].std, 0.1);
			//Btang=(1.4+0.02*numseq);
		}else{
			restr[trandpos].ang[angtype]=getrandomnewparam(INNERANGPARAM[aat][angtype].lower, INNERANGPARAM[aat][angtype].upper);
			//Btang=-1;
		}
	}

	if(angtype < 6){
		flagtor=param2coordbb(restr,trandpos,numseq);
		param2coordsugar(restr,trandpos,numseq);
		param2coordbase(restr,trandpos,numseq);
		param2coordhydrogen(restr,trandpos,numseq);
		if(!flagtor) printf("tor2str wrong in ang %d\n",trandpos);

		if(angtype == 5){
			int end = min(trandpos+2, numseq);
			coord2paramadd(restr,trandpos,end);
		}else{
			coord2paramadd(restr,trandpos,trandpos+1);
		}

		setangflap(numseq, trandpos, aat, angtype);
	}else{
		param2coordsugar(restr,trandpos,trandpos+1);
		param2coordbase(restr,trandpos,trandpos+1);
		param2coordhydrogen(restr,trandpos,trandpos+1);
		coord2paramadd(restr,trandpos,trandpos+1);
		setangntflap(numseq, trandpos, aat, angtype);
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














