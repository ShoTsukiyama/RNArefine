// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

double RNASimulation::getrandomnewparam(double l, double r){
	return l+(r-l)*Random();
}

double RNASimulation::getrandomnewtor(double l, double r, int rangetype){
	if(rangetype == 1){
		double r_ = r+360;
		double randoval = l+(r_-l)*Random();

		if(randoval >= 360){
			return randoval-360;
		}else{
			return randoval;
		}
	}else{
		return l+(r-l)*Random();
	}
}

double RNASimulation::getrandomgradparam(double crrp, double l, double r, vector<int> distr, vector<double> boundary, vector<double> energy, double span){
	double mtype=Random();
	
	if(mtype < 0.5){
		int binidx=getbin(crrp, boundary);
		if(binidx !=0 & binidx < energy.size()-1){
			if(energy[binidx] > energy[binidx+1] & energy[binidx-1] > energy[binidx+1]){
				return crrp+span*Random();
			}else if(energy[binidx] > energy[binidx-1] & energy[binidx+1] > energy[binidx-1] ){
				return crrp-span*Random();
			}else{
				return l+(r-l)*Random();
			}
		}else if(binidx ==0){
			return crrp+span*Random();
		}else{
			return crrp-span*Random();
		}
	}
	else{
		return l+(r-l)*Random();
	}
}

double RNASimulation::getnewparambinary(vector<int> binary, int total, vector<double> center, double span){
	int i;
	double pval;
	//double total;
	double sval;
	double rval;
	int trandpos;

    //total=0;
	//for(int i=0; i<binary.size();i++){
		//total+=int(1-binary[i]);
	//}

	rval=int(total*(Random()-EPS));

	sval=-1;
	for(trandpos=0;trandpos<binary.size();trandpos++)
	{
		if(binary[trandpos]==0) sval++;
		if(sval==rval) break;
	}
	if(trandpos == binary.size()){
		cout<<"Error444"<<endl;
		exit(1);
	}
	pval=center[trandpos] + (span/2)*(2*Random()-1);
	return pval;
}

double RNASimulation::getnewparam(vector<int> distr, int total, vector<double> center, double span){
	int i;
	double pval;
	double sval;
	double rval;
	int binidx=0;

	rval=total*(Random()-EPS);
	
	sval=0;
	for(int i=0; i<distr.size();i++){
		if(sval<rval & rval<=sval+distr[i]){
			binidx=i;
			break;
		}
		sval+=distr[i];
	}
	if(binidx >= distr.size()){
		cout<<binidx<<endl;
		cout<<distr.size()<<endl;
		cout<<"Error444"<<endl;
		exit(1);
	}
	pval=center[binidx] + (span/2)*(2*Random()-1);
	return pval;
}

double RNASimulation::getnewparam(vector<double> logdistr, double logtotal, vector<double> center, double span){
	int i;
	double pval;
	double sval;
	double rval;
	int binidx=0;

	rval=logtotal*(Random()-EPS);
	
	sval=0;
	for(int i=0; i<logdistr.size();i++){
		if(sval<rval & rval<=sval+logdistr[i]){
			binidx=i;
			break;
		}
		sval+=logdistr[i];
	}
	if(binidx >= logdistr.size()){
		cout<<binidx<<endl;
		cout<<logdistr.size()<<endl;
		cout<<"Error444"<<endl;
		exit(1);
	}
	pval=center[binidx] + (span/2)*(2*Random()-1);
	return pval;
}

void RNASimulation::selectrange(double l, double r, double &nl, double &nr)
{
	double valrange=r-l;
	double randval=Random();
	if(randval<0.1) randval=0.1;
	double nvalrange=randval*(r-l);
	double movedist=Random()*(valrange-nvalrange);

	nl=l+movedist;
	nr=nl+nvalrange;
}

double RNASimulation::gettemplature(double pos, double numseq, double scale){
	return (0.5*scale+0.5)-(0.5*(scale-1))*sin(((1/numseq)*(pos)-0.5)*M_PI);
}

bool RNASimulation::metropolisUpdate(double oldenergy, double newenergy, double Bt){
	double v, r;

	if(Bt < 0){
		if(newenergy<oldenergy){
			return true;
		}
		else{
			return false;
		}
	}
	else{
		if(newenergy<oldenergy){
			return true;
		}
		else{
			v=exp(-(newenergy-oldenergy)/(Bt));
			//cout<<to_string(newenergy-oldenergy) + "-" + to_string(v)<<endl;
			r=Random();
			if(r<v) return true;
			else return false;
		}
	}
}

int RNASimulation::getrandval(vector<double> tmov,double trandnum)
{
	int i;
	for(i=0;i<tmov.size()-1;i++)
	{
		if(trandnum>=tmov[i] && trandnum<tmov[i+1])
		{
			return i;
		}
	}
	return 0;
}

bool RNASimulation::fmcfragsweepcom(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,int &movtype,double Bt)
{	
	bool flagacc=true;
	bool flagres=false;
	vector<double> tmov;
	int ntype;
	double movetype=Random();

	//flagres=fmovementgradientlength2(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
	//flagres=fmovementlgm2(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
	//flagres=fmovementlenlocal(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
	//flagres=fmovementanglocal(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
	
	//movtype=ntype;
	//if(flagres) for(int i=0;i<30;i++) enelistbk[i]=enelist[i];
	//return flagres;

	if(enetype==0){
		if(BASEPAIRS.size() != 0 & STACKINGPAIRS.size() != 0) tmov={0,0.1,0.2,0.3,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.9,0.95,1.0};
		else if(BASEPAIRS.size() == 0 & STACKINGPAIRS.size() != 0) tmov={0,0.1,0.2,0.3,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.9,0.95,1.0};
		else if(BASEPAIRS.size() != 0 & STACKINGPAIRS.size() == 0) tmov={0,0.1,0.2,0.3,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.9,0.95,1.0};
		else tmov={0,0.1,0.2,0.3,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.9,0.95,1.0};
		//tmov={0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
		ntype=getrandval(tmov,movetype);

		//cout<<"============="<<endl;
		//cout<<ntype<<endl;
		//cout<<"============="<<endl;
	        
		switch(ntype)
		{
			case 0://0-0.1
				flagres=fmovementang(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				break;
			case 1://0.1-0.2
				flagres=fmovementlen(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				break;
			case 2://0.2-0.3
				flagres=fmovementtor(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				break;
			case 3://0.3-0.4
				flagres=fmovementlenang(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				break;
			case 4: //0.4-0.45
				flagres=fmovementrot(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				break;
			//case 4: //0.98 0.98
				//flagres=fmovementgeo(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				//break;
			case 5: //0.45-0.5
				flagres=fmovementmultitor(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				break;
			case 6: //0.5-0.55
				flagres=fmovementmultitorsugar(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				break;
			case 7: //0.55-0.6
				flagres=fmovementsugartors(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				break;
			case 8://0.6-0.65
				flagres=fmovementlmp(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				break;
			case 9: //0.65-0.7
				flagres=fmovementlmpfull(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				break;
			//case 10://0.7-0.8
				//flagres=fmovementlmpfullbasepair(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				//break;
			//case 10://0.7-0.8
				//flagres=fmovementgradientlength(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				//break;
			case 10://0.8-0.9
				flagres=fmovementgradientlength2(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				break;
			case 11: //0.9-0.95
				flagres=fmovementhydgeo(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				break;
			case 12://0.95-1.0
				flagres=fmovementgradhydgeo(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				break;
		}
	}else{
		if(BASEPAIRS.size() != 0 & STACKINGPAIRS.size() != 0) tmov={0,0.0625,0.125,0.1875,0.25,0.35,0.45,0.5,0.55,0.6,0.8,0.9,0.95,1.0};
		else if(BASEPAIRS.size() == 0 & STACKINGPAIRS.size() != 0) tmov={0,0.0625,0.125,0.1875,0.25,0.35,0.45,0.5,0.55,0.6,0.8,0.9,0.95,1.0};
		else if(BASEPAIRS.size() != 0 & STACKINGPAIRS.size() == 0) tmov={0,0.0625,0.125,0.1875,0.25,0.35,0.45,0.5,0.55,0.6,0.8,0.9,0.95,1.0};
		else tmov={0,0.0625,0.125,0.1875,0.25,0.35,0.45,0.5,0.55,0.6,0.8,0.9,0.95,1.0};
		//tmov={0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.45,0.55,0.6,0.65,0.7,0.8,0.9,1.0};
		ntype=getrandval(tmov,movetype);

		//cout<<"============="<<endl;
		//cout<<ntype<<endl;
		//cout<<"============="<<endl;

		//flagres=fmovementtor(decstr,numseq,oldenergy,newenergy,12,Bt,flagacc);
		//movtype=ntype;
		//if(flagres) for(int i=0;i<30;i++) enelistbk[i]=enelist[i];
		//return flagres;

		switch(ntype)
		{
			case 0://ang 0.1 0.1
				flagres=fmovementang(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				break;
			case 1:// 0.125
				flagres=fmovementangmt(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				//flagres=fmovementanglocal(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				break;
			case 2://len 0.225 0.2
				flagres=fmovementlen(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				break;
			case 3://0.25
				flagres=fmovementlenmt(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				//flagres=fmovementlenlocal(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				break;
			case 4://bet 0.35 0.3
				flagres=fmovementtor(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				break;
			case 5://bet 0.45 0.4
				flagres=fmovementlenang(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				break;
			case 6: // 0.5 0.45
				flagres=fmovementrot(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				break;
			//case 4:
				//flagres=fmovementgeo(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				//break;
			//case 7: //0.55
				//flagres=fmovementmultitor(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				//break;
			case 7: //0.55 0.5
				flagres=fmovementmultitorsugar(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				break;
			case 8: // 0.6 0.55
				flagres=fmovementsugartors(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				break;
			//case 9: //0.65 0.6
				//flagres=fmovementlmp(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				//break;
			//case 11:   //0.75
				//flagres=fmovementlmpfull(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				//break;
			//case 10://0.95
				//flagres=fmovementlmpfullbasepair(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				//break;
			//case 10:
				//flagres=fmovementgradientlength(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				//break;
			case 9://0.8
				flagres=fmovementgradientlength2(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				break;
			//case 12://0.89
				//flagres=fmovementlgm(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				//flagres=fmovementlgm(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				//break;
			case 10://0.9
				//flagres=fmovementlgm(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				flagres=fmovementlgm2(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				break;
			case 11: //0.95
				flagres=fmovementhydgeo(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				break;
			case 12: //1.0
				flagres=fmovementgradhydgeo(decstr,numseq,oldenergy,newenergy,Bt,flagacc);
				break;
		}
	}
	movtype=ntype;
	if(flagres) for(int i=0;i<30;i++) enelistbk[i]=enelist[i];
	return flagres;
}






