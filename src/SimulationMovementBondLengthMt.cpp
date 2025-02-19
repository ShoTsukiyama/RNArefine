// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

void RNASimulation::setlenntmtflap(double numseq, int trandpos, double aat, vector<int> targetlengths)
{
	int i,j;

	for (i=0; i < targetlengths.size(); ++i) {
		lp2.flap[BLENGTH_ENERGY_ST+targetlengths[i]][trandpos]=true;

		if(targetlengths[i] == 5 && trandpos < numseq-1){
			lp2.flap[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[0]][trandpos+1]=true;
		}
	}

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

bool RNASimulation::fmovementlenmt(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc)//length 
{

	int i,j,nexti;
	double r,l;
	double nr,nl;
	int trandpos;
	bool flagtor;
	int wpos;
	char aat;
	double pval;
	bool upidx;
	bool exitOuterLoop;
	int lentype;
	int wrong=0;
	double randval=Random();
	int uptype;
	
	int totwrong=0;
	for(i=0;i<BLENGTHTARGET.size();i++){
		totwrong+=lp1.nn[BLENGTH_ENERGY_ST+BLENGTHTARGET[i]];
	}

	if(totwrong>0){
		wpos=int(Random()*totwrong);

		i=0;
		for(trandpos=0;trandpos<numseq;trandpos++)
		{
			for(j=0;j<BLENGTHTARGET.size();j++){
				nexti=i+lp1.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][trandpos];
				if(i<=wpos & wpos<nexti){
					exitOuterLoop = true;
				}
				i=nexti;
			}if(exitOuterLoop) break;
		}
	}else{
		newenergy=oldenergy;
		return false;
	}

	vector<int> targetlengths;
	for(j=0;j<BLENGTHTARGET.size();j++){
		if(lp1.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][trandpos] == 1) targetlengths.push_back(BLENGTHTARGET[j]);
	}

	copystr(decstr, numseq);
	aat=restr[trandpos].aaa;

	if(trandpos==0){
		for(j=0;j<targetlengths.size();j++){
			if(targetlengths[j] == 10){
				restr[trandpos].len[targetlengths[j]]=getrandomnewparam(OPLENGTHPARAM[aat][0].lower, OPLENGTHPARAM[aat][0].upper);
			}else if(targetlengths[j] == 11){
				restr[trandpos].len[targetlengths[j]]=getrandomnewparam(OPLENGTHPARAM[aat][1].lower, OPLENGTHPARAM[aat][1].upper);
			}else{
				restr[trandpos].len[targetlengths[j]]=getrandomnewparam(LENGTHPARAM[aat][targetlengths[j]].lower, LENGTHPARAM[aat][targetlengths[j]].upper);
			}
		}
	}else{
		for(j=0;j<targetlengths.size();j++){
			restr[trandpos].len[targetlengths[j]]=getrandomnewparam(LENGTHPARAM[aat][targetlengths[j]].lower, LENGTHPARAM[aat][targetlengths[j]].upper);
		}
	}
	
	flagtor=param2coordbb(restr,trandpos,numseq);
	param2coordsugar(restr,trandpos,numseq);
	param2coordbase(restr,trandpos,numseq);
	param2coordhydrogen(restr,trandpos,numseq);
	if(!flagtor) printf("tor2str wrong in len %d\n",trandpos);

	if(find(targetlengths.begin(), targetlengths.end(), 5) != targetlengths.end()){
		int end = min(trandpos+2, numseq);
		coord2paramadd(restr,trandpos,end);
	}else{
		coord2paramadd(restr,trandpos,trandpos+1);
	}
	setlenntmtflap(numseq, trandpos, aat, targetlengths);

	newenergy=fcalcrmsdenergy(restr, numseq, lp2, rppmatbk, bpenergybk, bsenergybk, nmapbk, distmapbk, ntpairdistbk, adjpairdistbk, distenebk, multidistenebk, contactenebk, mapflagbk, mapbbflagbk, flagacc);
	upidx=metropolisUpdate(oldenergy, newenergy, 30);
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












