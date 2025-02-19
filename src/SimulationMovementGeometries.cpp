// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


void RNASimulation::setgeoflap(double numseq, int trandpos, int paramtype, double aat)
{
	int i,j;

	lp2.flap[BLENGTH_ENERGY_ST+paramtype][trandpos]=true;
	lp2.flap[INNERANGLE_ENERGY_ST+paramtype][trandpos]=true;
	lp2.flap[TORSIONANGLE_ENERGY_ST+paramtype][trandpos]=true;

	for(j=0;j<ADDBLENGTHTARGET.size();j++)
	{
		lp2.flap[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]][trandpos]=true;
	}

	if(paramtype == 5 && trandpos < numseq-1){
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

void RNASimulation::setgeontflap(double numseq, int trandpos, int paramtype, double aat)
{
	int i,j;

	lp2.flap[BLENGTH_ENERGY_ST+paramtype][trandpos]=true;
	lp2.flap[INNERANGLE_ENERGY_ST+paramtype][trandpos]=true;
	lp2.flap[TORSIONANGLE_ENERGY_ST+paramtype][trandpos]=true;

	for(j=0;j<ADDBLENGTHTARGET.size();j++)
	{
		lp2.flap[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]][trandpos]=true;
	}

	if(paramtype == 5 && trandpos < numseq-1){
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
		{
			if(j<trandpos-1 || i>trandpos+1){}
			else
			{
				mapflagbk[i][j]=true;
				mapflagbk[j][i]=true;
			}
		}
	}
	*/
}

bool RNASimulation::fmovementgeo(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc)//length 
{
	int i,j;
	double r,l;
	double nr,nl;
	int trandpos;
	bool flagtor;
	int wpos;
	int ss;
	char aat;
	double lengval, angval, torval;
	int atom1, atom2, atom3, atom4;
	Point3s pn;
	int paramtype=MODTARGET[int(Random()*MODTARGET.size())];
	//int totwrong=lp1.nn[BLENGTH_ENERGY_ST+paramtype]+lp1.nn[INNERANGLE_ENERGY_ST+paramtype]+lp1.nn[TORSIONANGLE_ENERGY_ST+paramtype];
	int totwrong=lp1.nn[BLENGTH_ENERGY_ST+paramtype]+lp1.nn[INNERANGLE_ENERGY_ST+paramtype];
	int wrong=0;
	double randval=Random();
	int st;
	bool flagts;
	BasicFunc bf;

	if(totwrong>0 && randval<0.6){
		wpos=int(Random()*totwrong);

		i=-1;
		for(trandpos=0;trandpos<numseq;trandpos++)
		{
			if(lp1.indn[BLENGTH_ENERGY_ST+paramtype][trandpos]==1) i++;
			if(i==wpos) break;
			if(lp1.indn[INNERANGLE_ENERGY_ST+paramtype][trandpos]==1) i++;
			if(i==wpos) break;
			//if(lp1.indn[TORSIONANGLE_ENERGY_ST+paramtype][trandpos]==1) i++;
			//if(i==wpos) break;
		}
	}
	else{
		trandpos=int((numseq)*Random());
	}

	if(trandpos >= numseq){
		trandpos=int((numseq)*Random());
	}

	copystr(decstr, numseq);

	aat=restr[trandpos].aaa;
	restr[trandpos].len[paramtype]=getrandomnewparam(LENGTHPARAM[aat][paramtype].lower, LENGTHPARAM[aat][paramtype].upper);
	restr[trandpos].ang[paramtype]=getrandomnewparam(INNERANGPARAM[aat][paramtype].lower, INNERANGPARAM[aat][paramtype].upper);
	//restr[trandpos].dih[paramtype]=getnewparambinary(TORSIONPARAM[aat][paramtype].binary, TORSIONPARAM[aat][paramtype].total, TORSIONPARAM[aat][paramtype].center, TORSIONPARAM[aat][paramtype].span);
	restr[trandpos].dih[paramtype]=getnewparam(TORSIONPARAM[aat][paramtype].logdistr, TORSIONPARAM[aat][paramtype].totallogdist, TORSIONPARAM[aat][paramtype].center, TORSIONPARAM[aat][paramtype].span);
	

	//flagtor=param2coordbb(restr,trandpos,numseq);
	//cout<<to_string(trandpos) + "-" + to_string(paramtype)<<endl;
	if(paramtype < 6){
		flagtor=param2coordbb(restr,trandpos,numseq);
		param2coordsugar(restr,trandpos,numseq);
		param2coordbase(restr,trandpos,numseq);
		param2coordhydrogen(restr,trandpos,numseq);

		if(paramtype == 5){
			int end = min(trandpos+2, numseq);
			coord2paramadd(restr,trandpos,end);
		}else{
			coord2paramadd(restr,trandpos,trandpos+1);
		}

		setgeoflap(numseq, trandpos, paramtype, aat);
	}else{
		param2coordsugar(restr,trandpos,trandpos+1);
		param2coordbase(restr,trandpos,trandpos+1);
		param2coordhydrogen(restr,trandpos,trandpos+1);
		coord2paramadd(restr,trandpos,trandpos+1);
		setgeontflap(numseq, trandpos, paramtype, aat);
	}

	newenergy=fcalcrmsdenergy(restr, numseq, lp2, rppmatbk, bpenergybk, bsenergybk, nmapbk, distmapbk, ntpairdistbk, adjpairdistbk, distenebk, multidistenebk, contactenebk, mapflagbk, mapbbflagbk, flagacc);
	bool upidx=metropolisUpdate(oldenergy, newenergy, Bt);
	
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












