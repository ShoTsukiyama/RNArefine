// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////



void RNASimulation::setmultitorflap(double numseq, double trandpos, double aat)
{
	int i,j;

	for(j=0;j<ADDBLENGTHTARGET.size();j++)
	{
		lp2.flap[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]][trandpos]=true;
	}

	if(trandpos < numseq-1){
		lp2.flap[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[0]][trandpos+1]=true;
	}

	if(trandpos > 0){
		//lp2.flap[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[5]][trandpos-1]=true;
		//lp2.flap[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[4]][trandpos-1]=true;
		lp2.flap[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[3]][trandpos-1]=true;
		lp2.flap[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[4]][trandpos-1]=true;
		lp2.flap[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[5]][trandpos-1]=true;
		lp2.flap[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[6]][trandpos-1]=true;
	}

	for(j=0;j<ADDINNERANGLETARGET.size();j++)
	{
		lp2.flap[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]][trandpos]=true;
	}
	//if(trandpos > 0){
		//for(j=0;j<ADDINNERANGLETARGET.size();j++)
		//{
			//lp2.flap[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]][trandpos-1]=true;
		//}
	//}

	for(i=0;i<numseq;i++)
	{
		for(j=i;j<numseq;j++)
		//for(j=i+1;j<numseq;j++)
		{
			if(j<trandpos-1 || i>trandpos+1){}
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

double RNASimulation::getmovmultitorval(double val, double mean, double k){
	double diff, diff_inv;
	double newval;
	if(val < mean){
		diff=mean-val;
		diff_inv=360-diff;
		if(diff_inv < diff){
			newval = val - diff_inv*k*Random();
		}else{
			newval = val + diff*k*Random();
		}
	}else{
		diff=val-mean;
		diff_inv=360-diff;
		if(diff_inv < diff){
			newval = val + diff_inv*k*Random();
		}else{
			newval = val - diff*k*Random();
		}
	}

	if(newval > 360) newval = newval - 360;
	if(newval < 0) newval = newval + 360;
	return newval;
}

bool RNASimulation::fmovementmultitor(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt, bool flagacc)//dihedral
{
	int i,j,inext;
	double r,l;
	double nr,nl;
	int trandpos, dihtype, dihcls;
	bool flagtor;
	int wpos;
	char aat;
	double pval;
	bool upidx;
	//int totwrong=numseq+lp1.nn[TORSIONANGLE_ENERGY_ST+dihtype];

	//int tortype=int(Random()*TORSIONANGLETARGET.size());
	//int dihtype=TORSIONANGLETARGET[int(Random()*TORSIONANGLETARGET.size())];
	//int totwrong=lp1.nn[TORSION_CLS_ENERGY_ST];

	int wrong=0;
	double randval=Random();
	bool exitOuterLoop = false;

	/*
	if(totwrong>0 && randval<0.8){
		wpos=int(Random()*totwrong);
		i=-1;
		for(trandpos=1;trandpos<numseq;trandpos++)
		{
			for(dihtype=0;dihtype<TorsionClsList[0].size();dihtype++){
				if(toroutflag[trandpos][dihtype]) i++;
				
				if(i==wpos){
					exitOuterLoop = true;
					break;
				}
			}
			if(exitOuterLoop){
				break;
			}
		}

		dihcls=torclsidx[trandpos];
	}
	*/

	int totwrong1=0;
	for(i=0;i<BLENGTHTARGET.size();i++){
		totwrong1+=lp1.nn[BLENGTH_ENERGY_ST+BLENGTHTARGET[i]];
		totwrong1-=lp1.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[i]][0];
	}
	for(i=0;i<INNERANGLETARGET.size();i++){
		totwrong1+=lp1.nn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[i]];
		totwrong1-=lp1.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[i]][0];
	}

	int totwrong2=lp1.nn[CLASH_ENERGY_ST]*2-lp1.indn[CLASH_ENERGY_ST][0];

	//int totene=BLENGTHTARGET.size()+INNERANGLETARGET.size();
	//totene=totene*numseq-totene;
	if(totwrong1>0 && randval<0.4){
		wpos=int(Random()*totwrong1);

		i=0;
		for(trandpos=1;trandpos<numseq;trandpos++)
		{
			for(j=0;j<BLENGTHTARGET.size();j++){
				if(i<=wpos & wpos<i+lp1.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][trandpos]){
					exitOuterLoop = true;
				}
				i+=lp1.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][trandpos];
			}

			for(j=0;j<INNERANGLETARGET.size();j++){
				if(i<=wpos & wpos<i+lp1.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][trandpos]){
					exitOuterLoop = true;
				}
				i+=lp1.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][trandpos];
			}
			if(exitOuterLoop) break;
		}
		//dihcls=torclsidx[trandpos];
	}else if(totwrong2>0 && 0.4<=randval && randval<0.8){
		wpos=int(Random()*totwrong2);

		i=0;
		for(trandpos=1;trandpos<numseq;trandpos++)
		{
			inext=i+lp1.indn[CLASH_ENERGY_ST][trandpos];
			if(i<=wpos && wpos<inext){
				break;
			}
			i=inext;
		}
	}else{
		trandpos=int((numseq-1)*Random())+1;
	}

	if(trandpos >= numseq){
		trandpos=int((numseq-1)*Random())+1;
	}

	dihcls=int(TORSIONLIST[0].size()*Random());
	copystr(decstr, numseq);
	aat=restr[trandpos].aaa;

	randval=Random();
	//if(randval > 0.5){
		//restr[trandpos-1].dih[TorsionClsList[0][0]]=getmovmultitorval(restr[trandpos-1].dih[TorsionClsList[0][0]], TORSIONLIST[0][torclsidx[trandpos]][0].mean, 0.1);
		//lp2.flap[TORSIONANGLE_ENERGY_ST+TorsionClsList[0][0]][trandpos-1]=true;
		//for(i=1;i<toroutflag[trandpos].size();i++){
			//restr[trandpos].dih[TorsionClsList[0][j]]=getmovmultitorval(restr[trandpos].dih[TorsionClsList[0][j]], TORSIONLIST[0][torclsidx[trandpos]][j].mean, 0.1);
			//lp2.flap[TORSIONANGLE_ENERGY_ST+TorsionClsList[0][j]][trandpos]=true;
		//}
	//}else{

	restr[trandpos-1].dih[TorsionClsList[0][0]]=getrandomnewtor(TORSIONLIST[0][dihcls][0].lower, TORSIONLIST[0][dihcls][0].upper, TORSIONLIST[0][dihcls][0].rangetype);
	restr[trandpos-1].ang[TorsionClsList[0][0]]=getrandomnewparam(INNERANGPARAM[aat][TorsionClsList[0][0]].lower, INNERANGPARAM[aat][TorsionClsList[0][0]].upper);
	restr[trandpos-1].len[TorsionClsList[0][0]]=getrandomnewparam(LENGTHPARAM[aat][TorsionClsList[0][0]].lower, LENGTHPARAM[aat][TorsionClsList[0][0]].upper);
	lp2.flap[TORSIONANGLE_ENERGY_ST+TorsionClsList[0][0]][trandpos-1]=true;
	lp2.flap[INNERANGLE_ENERGY_ST+TorsionClsList[0][0]][trandpos-1]=true;
	lp2.flap[BLENGTH_ENERGY_ST+TorsionClsList[0][0]][trandpos-1]=true;

	for(i=1;i<TorsionClsList[0].size();i++){
		restr[trandpos].dih[TorsionClsList[0][i]]=getrandomnewtor(TORSIONLIST[0][dihcls][i].lower, TORSIONLIST[0][dihcls][i].upper, TORSIONLIST[0][dihcls][i].rangetype);
		restr[trandpos].ang[TorsionClsList[0][i]]=getrandomnewparam(INNERANGPARAM[aat][TorsionClsList[0][i]].lower, INNERANGPARAM[aat][TorsionClsList[0][i]].upper);
		restr[trandpos].len[TorsionClsList[0][i]]=getrandomnewparam(LENGTHPARAM[aat][TorsionClsList[0][i]].lower, LENGTHPARAM[aat][TorsionClsList[0][i]].upper);
		lp2.flap[TORSIONANGLE_ENERGY_ST+TorsionClsList[0][i]][trandpos]=true;
		lp2.flap[INNERANGLE_ENERGY_ST+TorsionClsList[0][i]][trandpos]=true;
		lp2.flap[BLENGTH_ENERGY_ST+TorsionClsList[0][i]][trandpos]=true;
	}
	//}

	flagtor=param2coordbb(restr,trandpos-1,numseq);
	param2coordsugar(restr,trandpos-1,numseq);
	param2coordbase(restr,trandpos-1,numseq);
	param2coordhydrogen(restr,trandpos-1,numseq);
	if(!flagtor) printf("tor2str wrong in len %d\n",trandpos);

	int end = min(trandpos+2, numseq);
	coord2paramadd(restr,trandpos-1,end);

	setmultitorflap(numseq, trandpos, aat);
	
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





