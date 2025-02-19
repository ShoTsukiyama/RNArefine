// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

void RNASimulation::setmultianglenflap(double numseq, double trandpos, double aat, int span)
{
	int i,j;

	for(i=0;i<numseq;i++)
	{
		for(j=i;j<numseq;j++)
		//for(j=i+1;j<numseq;j++)
		{
			if(j<trandpos-1 || i>trandpos+span+1){}
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

bool RNASimulation::fmovementlenang(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt, bool flagacc)//dihedral
{
	int i,j;
	double r,l;
	double nr,nl;
	int trandpos, trandposed, trandposst, span, dihtype, dihcls;
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

	int totene=BLENGTHTARGET.size()*numseq+INNERANGLETARGET.size()*numseq;
	int totwrong=0;
	for(i=0;i<BLENGTHTARGET.size();i++){
		totwrong+=lp1.nn[BLENGTH_ENERGY_ST+BLENGTHTARGET[i]];
	}
	for(i=0;i<INNERANGLETARGET.size();i++){
		totwrong+=lp1.nn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[i]];
	}

	if(totwrong>0 && randval<0.5+0.5*(totwrong/totene)){
		wpos=int(Random()*totwrong);

		i=-1;
		for(trandposst=0;trandposst<numseq;trandposst++)
		{
			for(j=0;j<BLENGTHTARGET.size();j++)
			{

				if(lp1.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][trandposst] > 0) i++;
				if(i == wpos){
					exitOuterLoop = true;
					break;
				}
			}

			for(j=0;j<INNERANGLETARGET.size();j++){

				if(lp1.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][trandposst] > 0) i++;
				if(i == wpos){
					exitOuterLoop = true;
					break;
				}
			}
			if(exitOuterLoop) break;
		}
	}else{
		trandposst=int((numseq)*Random());
	}

	if(trandposst >= numseq){
		trandposst=int((numseq)*Random());
	}

	span=int(5*Random())+1;
	trandposed=min(trandposst+span, numseq);

	copystr(decstr, numseq);

	for(trandpos=trandposst;trandpos<trandposed;trandpos++){
		aat=restr[trandpos].aaa;

		//randval=Random();

		for(j=0;j<BLENGTHTARGET.size();j++){
			//randval=Random();
			if(lp1.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][trandpos] == 1){
				restr[trandpos].len[BLENGTHTARGET[j]]=getmovlenval(restr[trandpos].len[BLENGTHTARGET[j]], LENGTHPARAM[aat][BLENGTHTARGET[j]].mean, LENGTHPARAM[aat][BLENGTHTARGET[j]].std, 0.1);

				//restr[trandpos].len[BLENGTHTARGET[j]]=getrandomnewparam(LENGTHPARAM[aat][BLENGTHTARGET[j]].lower, LENGTHPARAM[aat][BLENGTHTARGET[j]].upper);

				lp2.flap[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][trandpos]=true;
			}
		}

		for(j=0;j<INNERANGLETARGET.size();j++){
			//randval=Random();
			if(lp1.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][trandpos] == 1){
				restr[trandpos].ang[INNERANGLETARGET[j]]=getmovangval(restr[trandpos].ang[INNERANGLETARGET[j]], INNERANGPARAM[aat][INNERANGLETARGET[j]].mean, INNERANGPARAM[aat][INNERANGLETARGET[j]].std, 0.1);

				//restr[trandpos].ang[INNERANGLETARGET[j]]=getrandomnewparam(INNERANGPARAM[aat][INNERANGLETARGET[j]].lower, INNERANGPARAM[aat][INNERANGLETARGET[j]].upper);

				lp2.flap[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][trandpos]=true;
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
	}

	if(trandposed < numseq){
		lp2.flap[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[0]][trandposed]=true;
	}

	flagtor=param2coordbb(restr,trandposst,numseq);
	param2coordsugar(restr,trandposst,numseq);
	param2coordbase(restr,trandposst,numseq);
	param2coordhydrogen(restr,trandposst,numseq);
	//setmultitorflap(numseq, trandpos, aat, dihtype);

	coord2paramadd(restr,trandposst,trandposed);

	setmultianglenflap(numseq, trandposst, aat, span);

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





