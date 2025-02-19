// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////



void RNASimulation::setmultitorsugarflap(double numseq, double trandpos, double aat)
{
	int i,j;

	for(j=0;j<ADDBLENGTHTARGET.size();j++)
	{
		lp2.flap[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]][trandpos]=true;
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

bool RNASimulation::fmovementmultitorsugar(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt, bool flagacc)//dihedral
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

	int wrong=0;
	double randval=Random();
	bool exitOuterLoop = false;

	int totwrong1=0;
	for(i=0;i<BLENGTHTARGET.size();i++){
		totwrong1+=lp1.nn[BLENGTH_ENERGY_ST+BLENGTHTARGET[i]];
	}
	for(i=0;i<INNERANGLETARGET.size();i++){
		totwrong1+=lp1.nn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[i]];
	}

	int totwrong2=lp1.nn[CLASH_ENERGY_ST]*2;

	//int totene=BLENGTHTARGET.size()+INNERANGLETARGET.size();
	//totene=totene*numseq-totene;
	if(totwrong1>0 && randval<0.4){
		wpos=int(Random()*totwrong1);

		i=0;
		for(trandpos=0;trandpos<numseq;trandpos++)
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
		for(trandpos=0;trandpos<numseq;trandpos++)
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

	dihcls=int(TORSIONLIST[1].size()*Random());
	copystr(decstr, numseq);
	aat=restr[trandpos].aaa;

	randval=Random();

	for(i=0;i<TorsionClsList[1].size();i++)
	{
		restr[trandpos].dih[TorsionClsList[1][i]]=getrandomnewtor(TORSIONLIST[1][dihcls][i].lower, TORSIONLIST[1][dihcls][i].upper, TORSIONLIST[1][dihcls][i].rangetype);
		restr[trandpos].ang[TorsionClsList[1][i]]=getrandomnewparam(INNERANGPARAM[aat][TorsionClsList[1][i]].lower, INNERANGPARAM[aat][TorsionClsList[1][i]].upper);
		restr[trandpos].len[TorsionClsList[1][i]]=getrandomnewparam(LENGTHPARAM[aat][TorsionClsList[1][i]].lower, LENGTHPARAM[aat][TorsionClsList[1][i]].upper);
		
		//cout<<to_string(i) + "-" + to_string(TorsionClsList[1][i]) + "-" + to_string(restr[trandpos].dih[TorsionClsList[1][i]])<<endl;
		
		lp2.flap[TORSIONANGLE_ENERGY_ST+TorsionClsList[1][i]][trandpos]=true;
		lp2.flap[INNERANGLE_ENERGY_ST+TorsionClsList[1][i]][trandpos]=true;
		lp2.flap[BLENGTH_ENERGY_ST+TorsionClsList[1][i]][trandpos]=true;
	}
	//}

	flagtor=param2coordbb(restr,trandpos,numseq);
	param2coordsugar(restr,trandpos,numseq);
	param2coordbase(restr,trandpos,numseq);
	param2coordhydrogen(restr,trandpos,numseq);

	if(!flagtor) printf("tor2str wrong in len %d\n",trandpos);

	int end = min(trandpos+2, numseq);

	coord2paramadd(restr,trandpos,end);

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





