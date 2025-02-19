// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

void RNASimulation::setangleflaplocal(double numseq, int trandpos, double aat, int angtype)
{
	int i,j;

	lp2.flap[INNERANGLE_ENERGY_ST+angtype][trandpos]=true;

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

void RNASimulation::changelocalang(vector<Residueatoms> &decstr, int trandpos, int angtype, int numseq)
{
	char aat=decstr[trandpos].aaa;
	Point3s pn;
	bool flagts;
	int atom1, atom2, atom3, atom4;
	int moveatom;

	if(trandpos == 0){
		if(angtype == 10 || angtype == 11){
			atom1=FIRSTOPParamsList[aat][angtype][0];atom2=FIRSTOPParamsList[aat][angtype][1];atom3=FIRSTOPParamsList[aat][angtype][2];moveatom=FIRSTOPParamsList[aat][angtype][3];
		}else{
			atom1=ParamsList[aat][angtype][0];atom2=ParamsList[aat][angtype][1];atom3=ParamsList[aat][angtype][2];moveatom=ParamsList[aat][angtype][3];
		}

		flagts=bf.tor2pos22(decstr[trandpos].ptv[atom1].x,decstr[trandpos].ptv[atom1].y,decstr[trandpos].ptv[atom1].z,
							decstr[trandpos].ptv[atom2].x,decstr[trandpos].ptv[atom2].y,decstr[trandpos].ptv[atom2].z,
							decstr[trandpos].ptv[atom3].x,decstr[trandpos].ptv[atom3].y,decstr[trandpos].ptv[atom3].z,
							decstr[trandpos].dih[angtype]*RADDEG,decstr[trandpos].len[angtype],decstr[trandpos].ang[angtype]*RADDEG,
							&pn.x,&pn.y,&pn.z);
	}else{
		atom1=ParamsList[aat][angtype][0];atom2=ParamsList[aat][angtype][1];atom3=ParamsList[aat][angtype][2];moveatom=ParamsList[aat][angtype][3];
		
		flagts=bf.tor2pos22(decstr[trandpos].ptv[atom1].x,decstr[trandpos].ptv[atom1].y,decstr[trandpos].ptv[atom1].z,
							decstr[trandpos].ptv[atom2].x,decstr[trandpos].ptv[atom2].y,decstr[trandpos].ptv[atom2].z,
							decstr[trandpos].ptv[atom3].x,decstr[trandpos].ptv[atom3].y,decstr[trandpos].ptv[atom3].z,
							decstr[trandpos].dih[angtype]*RADDEG,decstr[trandpos].len[angtype],decstr[trandpos].ang[angtype]*RADDEG,
							&pn.x,&pn.y,&pn.z);
	}

	decstr[trandpos].ptv[moveatom].x=pn.x;
	decstr[trandpos].ptv[moveatom].y=pn.y;
	decstr[trandpos].ptv[moveatom].z=pn.z;

	int iend=trandpos+2;
	if(iend > numseq) iend=numseq;
	coord2paramlocal(decstr, trandpos, moveatom, trandpos, iend);
	coord2paramsugarlocal(decstr, trandpos, moveatom, trandpos, iend);
	coord2paramaddlocal(decstr, trandpos, moveatom, trandpos, iend);
}


bool RNASimulation::fmovementanglocal(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc)//angle 
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

	if(totwrong>0){
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
		newenergy=oldenergy;
        return false;
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

	changelocalang(restr, trandpos, angtype, numseq);

	param2coordbase(restr,trandpos,numseq);
	param2coordhydrogen(restr,trandpos,numseq);

	setangleflaplocal(numseq, trandpos, aat, angtype);

	newenergy=fcalcrmsdenergy(restr, numseq, lp2, rppmatbk, bpenergybk, bsenergybk, nmapbk, distmapbk, ntpairdistbk, adjpairdistbk, distenebk, multidistenebk, contactenebk, mapflagbk, mapbbflagbk, flagacc);
	upidx=metropolisUpdate(oldenergy, newenergy, 20);
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











