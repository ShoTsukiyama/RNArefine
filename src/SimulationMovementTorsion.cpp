// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////



void RNASimulation::settorflap(double numseq, double trandpos, double aat, int dihtype)
{
	int i,j;

	lp2.flap[TORSIONANGLE_ENERGY_ST+dihtype][trandpos]=true;

	for(j=0;j<ADDBLENGTHTARGET.size();j++)
	{
		lp2.flap[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]][trandpos]=true;
	}

	if(dihtype == 5 && trandpos < numseq-1){
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

void RNASimulation::settorntflap(double numseq, double trandpos, double aat, int dihtype)
{
	int i,j;

	lp2.flap[TORSIONANGLE_ENERGY_ST+dihtype][trandpos]=true;

	for(j=0;j<ADDBLENGTHTARGET.size();j++)
	{
		lp2.flap[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]][trandpos]=true;
	}

	if(dihtype == 5 && trandpos < numseq-1){
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

double RNASimulation::getmovtorval(double val, double std, double k){
	double newval;
	double randval=Random();
	if(randval > 0.5){
		newval = val + std*k*Random();
	}else{
		newval = val - std*k*Random();
	}

	if(newval > 360) newval = newval - 360;
	if(newval < 0) newval = newval + 360;

	return newval;
}

bool RNASimulation::fmovementtor(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt, bool flagacc)//dihedral
{
	int i,j,inext;
	double r,l;
	double nr,nl;
	int trandpos, trandendpos;
	bool flagtor;
	int wpos;
	char aat;
	double dist1, dist2;
	double pval;
	bool upidx;
	//int totwrong=numseq+lp1.nn[TORSIONANGLE_ENERGY_ST+dihtype];

	//int tortype=int(Random()*TORSIONANGLETARGET.size());
	int dihtype=TORSIONANGLETARGET[int(Random()*TORSIONANGLETARGET.size())];
	//int totwrong=lp1.nn[TORSIONANGLE_ENERGY_ST+dihtype];
	int wrong=0;
	//if(totwrong>0 && randval<0.8){
		//wpos=int(Random()*totwrong);
		//i=-1;
		//for(trandpos=0;trandpos<numseq;trandpos++)
		//{
			//if(lp1.indn[TORSIONANGLE_ENERGY_ST+dihtype][trandpos]==1) i++;
			//if(i==wpos) break;
		//}
	//}
	//else{
	//trandpos=int((numseq)*Random());
	//}

	double randval=Random();
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

	copystr(decstr, numseq);
	aat=restr[trandpos].aaa;

	//restr[trandpos].dih[dihtype]=getnewparambinary(TORSIONPARAM[aat][dihtype].binary, TORSIONPARAM[aat][dihtype].total, TORSIONPARAM[aat][dihtype].center, TORSIONPARAM[aat][dihtype].span);

	//restr[trandpos].dih[dihtype]=getnewparam(TORSIONPARAM[aat][dihtype].distr, TORSIONPARAM[aat][dihtype].totaldist, TORSIONPARAM[aat][dihtype].center, TORSIONPARAM[aat][dihtype].span);
	if(dihtype < 6 | dihtype == 10 | dihtype == 11 | dihtype == 13 | dihtype == 14){
		restr[trandpos].dih[dihtype]=getnewparam(TORSIONPARAM[aat][dihtype].logdistr, TORSIONPARAM[aat][dihtype].totallogdist, TORSIONPARAM[aat][dihtype].center, TORSIONPARAM[aat][dihtype].span);
	}else{

		dist1=getdisttor(restr[trandpos].dih[ClsTorParams[aat][0]], TORCLSPARAM[aat][0][0].mean);
		dist2=getdisttor(restr[trandpos].dih[ClsTorParams[aat][0]], TORCLSPARAM[aat][1][0].mean);


		if(dist1 < dist2){
			if(dihtype == 6) restr[trandpos].dih[dihtype]=getrandomnewparam(TORCLSPARAM[aat][0][1].lower, TORCLSPARAM[aat][0][1].upper);
			if(dihtype == 7) restr[trandpos].dih[dihtype]=getrandomnewparam(TORCLSPARAM[aat][0][2].lower, TORCLSPARAM[aat][0][2].upper);
			if(dihtype == 8) restr[trandpos].dih[dihtype]=getrandomnewparam(TORCLSPARAM[aat][0][3].lower, TORCLSPARAM[aat][0][3].upper);
			if(dihtype == 9) restr[trandpos].dih[dihtype]=getrandomnewparam(TORCLSPARAM[aat][0][4].lower, TORCLSPARAM[aat][0][4].upper);
			if(dihtype == 12) restr[trandpos].dih[dihtype]=getrandomnewparam(TORCLSPARAM[aat][0][5].lower, TORCLSPARAM[aat][0][5].upper);
		}else{
			if(dihtype == 6) restr[trandpos].dih[dihtype]=getrandomnewparam(TORCLSPARAM[aat][1][1].lower, TORCLSPARAM[aat][1][1].upper);
			if(dihtype == 7) restr[trandpos].dih[dihtype]=getrandomnewparam(TORCLSPARAM[aat][1][2].lower, TORCLSPARAM[aat][1][2].upper);
			if(dihtype == 8) restr[trandpos].dih[dihtype]=getrandomnewparam(TORCLSPARAM[aat][1][3].lower, TORCLSPARAM[aat][1][3].upper);
			if(dihtype == 9) restr[trandpos].dih[dihtype]=getrandomnewparam(TORCLSPARAM[aat][1][4].lower, TORCLSPARAM[aat][1][4].upper);
			if(dihtype == 12) restr[trandpos].dih[dihtype]=getrandomnewparam(TORCLSPARAM[aat][1][5].lower, TORCLSPARAM[aat][1][5].upper);
		}
	}
	
	//double randval=Random();
	//if(randval > 0.5){
	//restr[trandpos].dih[dihtype]=getmovtorval(restr[trandpos].dih[dihtype], TORSIONPARAM[aat][dihtype].std, 0.1);
	//}else{
		//restr[trandpos].dih[dihtype]=getnewparam(TORSIONPARAM[aat][dihtype].logdistr, TORSIONPARAM[aat][dihtype].totallogdist, TORSIONPARAM[aat][dihtype].center, TORSIONPARAM[aat][dihtype].span);
	//}

	if(dihtype < 6){
		flagtor=param2coordbb(restr,trandpos,numseq);
		param2coordsugar(restr,trandpos,numseq);
		param2coordbase(restr,trandpos,numseq);
		param2coordhydrogen(restr,trandpos,numseq);
		if(!flagtor) printf("tor2str wrong in dih %d %d\n",trandpos, dihtype);

		if(dihtype == 5){
			int end = min(trandpos+2, numseq);
			coord2paramadd(restr,trandpos,end);
		}else{
			coord2paramadd(restr,trandpos,trandpos+1);
		}

		settorflap(numseq, trandpos, aat, dihtype);
	}else{
		param2coordsugar(restr,trandpos,trandpos+1);
		param2coordbase(restr,trandpos,trandpos+1);
		param2coordhydrogen(restr,trandpos,trandpos+1);
		coord2paramadd(restr,trandpos,trandpos+1);
		settorntflap(numseq, trandpos, aat, dihtype);
	}

	newenergy=fcalcrmsdenergy(restr, numseq, lp2, rppmatbk, bpenergybk, bsenergybk, nmapbk, distmapbk, ntpairdistbk, adjpairdistbk, distenebk, multidistenebk, contactenebk, mapflagbk, mapbbflagbk, flagacc);
	upidx=metropolisUpdate(oldenergy, newenergy, Bt);
	
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
	
	//bool upidx=metropolisUpdate(oldenergy, newenergy, Bt);
	//if(newenergy<oldenergy){
	if(upidx){
		rcopystr(decstr, restr, numseq);
		return true;
		
	}
	else
	{
		return false;
	}
}





