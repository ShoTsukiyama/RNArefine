// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


void RNASimulation::setgeohydflap(double numseq, int trandpos, int paramtype, double aat)
{
	int i,j;

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

bool RNASimulation::fmovementhydgeo(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc)//length 
{
	int i,j;
	double r,l;
	double nr,nl;
	int trandpos;
	bool flagtor;
	int wpos;
	int ss;
	char aat;
	int paramtype;
	double lengval, angval, torval;
	int atom1, atom2, atom3, atom4;
	Point3s pn;
	//int totwrong=lp1.nn[BLENGTH_ENERGY_ST+paramtype]+lp1.nn[INNERANGLE_ENERGY_ST+paramtype]+lp1.nn[TORSIONANGLE_ENERGY_ST+paramtype];
	//int totwrong=lp1.nn[BLENGTH_ENERGY_ST+paramtype]+lp1.nn[INNERANGLE_ENERGY_ST+paramtype];
	int wrong=0;
	double randval=Random();
	int st;
	bool flagts;
	BasicFunc bf;

	trandpos=int((numseq)*Random());

	copystr(decstr, numseq);

	aat=restr[trandpos].aaa;
	paramtype=int(Random()*HydParamsList[aat].size());

	//cout<<to_string(trandpos) + "-" + to_string(paramtype)<<endl;

	//restr[trandpos].hydlen[paramtype]=getrandomnewparam(LENGTHPARAM[aat][paramtype].lower, LENGTHPARAM[aat][paramtype].upper);
	restr[trandpos].hydlen[paramtype]=getnewparambinary(HYDLENGTHPARAM[aat][paramtype].binary, HYDLENGTHPARAM[aat][paramtype].total, HYDLENGTHPARAM[aat][paramtype].center, HYDLENGTHPARAM[aat][paramtype].span);
	restr[trandpos].hydang[paramtype]=getrandomnewparam(HYDINNERANGPARAM[aat][paramtype].lower, HYDINNERANGPARAM[aat][paramtype].upper);
	restr[trandpos].hyddih[paramtype]=getnewparambinary(HYDTORSIONPARAM[aat][paramtype].binary, HYDTORSIONPARAM[aat][paramtype].total, HYDTORSIONPARAM[aat][paramtype].center, HYDTORSIONPARAM[aat][paramtype].span);
	
	//cout<<to_string(trandpos) + "-" + to_string(paramtype) + "-" + to_string(restr[trandpos].hydlen[paramtype]) + "-" + to_string(restr[trandpos].hydang[paramtype]) + "-" + to_string(restr[trandpos].hyddih[paramtype])<<endl;


	param2coordhydrogen(restr,trandpos,trandpos+1);

	setgeohydflap(numseq, trandpos, paramtype, aat);
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












