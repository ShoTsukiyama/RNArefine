// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////



void RNASimulation::setsugartorsflap(double numseq, int trandpos1, int trandpos2)
{
	int i,j;

	for(i=trandpos1;i<trandpos2;i++){
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

bool RNASimulation::fmovementsugartors(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt, bool flagacc)//dihedral
{
	int i,j,nexti;
	double r,l;
	double nr,nl;
	bool flagtor;
	int wpos;
	char aat;
	double pval;
	bool upidx;
	double dist1, dist2;
	int trandpos, trandposed, trandposst, span;
	bool flagdone;
	double randval=Random();
	bool exitOuterLoop = false;
	int totwrong=0;
	for(j=0;j<6;j++)
	{
		totwrong+=lp1.nn[TORSIONANGLE_ENERGY_ST+ClsTorParams['A'][j]];
	}
	
	if(totwrong>0 && randval<0.5){
		wpos=int(Random()*totwrong);

		i=0;
		for(trandpos=0;trandpos<numseq;trandpos++)
		{
			aat=decstr[i].aaa;

			for(j=0;j<6;j++)
			{
				nexti=i+lp1.indn[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]][trandpos];
				if(i<=wpos && wpos<nexti){
					exitOuterLoop = true;
					break;
				}
				i=nexti;
			}if(exitOuterLoop) break;
		}
	}else{
		trandpos=int((numseq)*Random());
	}

	copystr(decstr, numseq);
	
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

	for(i=trandposst;i<trandposed;i++)
	{
		aat=restr[i].aaa;

		for(auto& idx : SUGARRINGPARAMS){
			restr[i].len[BLENGTHTARGET[idx]]=getrandomnewparam(LENGTHPARAM[aat][BLENGTHTARGET[idx]].lower, LENGTHPARAM[aat][BLENGTHTARGET[idx]].upper);
			lp2.flap[BLENGTH_ENERGY_ST+BLENGTHTARGET[idx]][i]=true;

			restr[i].ang[INNERANGLETARGET[idx]]=getrandomnewparam(INNERANGPARAM[aat][INNERANGLETARGET[idx]].lower, INNERANGPARAM[aat][INNERANGLETARGET[idx]].upper);
			lp2.flap[INNERANGLE_ENERGY_ST+INNERANGLETARGET[idx]][i]=true;
		}

		dist1=getdisttor(restr[i].dih[ClsTorParams[aat][0]], TORCLSPARAM[aat][0][0].mean);
		dist2=getdisttor(restr[i].dih[ClsTorParams[aat][0]], TORCLSPARAM[aat][1][0].mean);

		if(dist1 < dist2){
			for(j=1;j<ClsTorParams[aat].size();j++){
				restr[i].dih[ClsTorParams[aat][j]]=getrandomnewparam(TORCLSPARAM[aat][0][j].lower, TORCLSPARAM[aat][0][j].upper);
				lp2.flap[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]][i]=true;
			}
		}else{
			for(j=1;j<ClsTorParams[aat].size();j++){
				restr[i].dih[ClsTorParams[aat][j]]=getrandomnewparam(TORCLSPARAM[aat][1][j].lower, TORCLSPARAM[aat][1][j].upper);
				lp2.flap[TORSIONANGLE_ENERGY_ST+ClsTorParams[aat][j]][i]=true;
			}
		}
	}

	param2coordsugar(restr,trandposst,trandposed);
	param2coordbase(restr,trandposst,trandposed);
	param2coordhydrogen(restr,trandposst,trandposed);
	coord2paramadd(restr,trandposst,trandposed);

	setsugartorsflap(numseq,trandposst,trandposed);

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





