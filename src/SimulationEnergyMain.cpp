// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

bool RNASimulation::checkoutlier(double val, vector<double> boundary, vector<int> binary){
	int binidx=getbin(val, boundary);
	int bdiff=binary[binidx];
	if(bdiff == 1) return false;
	else return true;
}

void RNASimulation::getoutenelength(double val, double mean, double lower, double upper,  double slower, double supper, double k, double &diff, int &bdiff){
	double tmp;

	if(lower <= val && val <= upper){
		bdiff=0;
		diff=0;
	}else if(lower > val){
		bdiff=1;
		tmp=abs(mean-val);
		diff=k*tmp;
	}else if(val > upper){
		bdiff=1;
		tmp=abs(val-mean);
		diff=k*tmp;
	}

	/*
	else if(lower > val){
		bdiff=1;
		tmp=abs(lower-val);
		if(tmp < 1) diff=k*tmp;
		else diff=k*sqrt(tmp*tmp);
	}else if(val > upper){
		bdiff=1;
		tmp=abs(val-upper);
		if(tmp < 1) diff=k*tmp;
		else diff=k*sqrt(tmp*tmp);
	}
	*/
}

void RNASimulation::getouteneangle(double val, double mean, double lower, double upper,  double slower, double supper, double k, double &diff, int &bdiff){
	double tmp;

	if(lower <= val && val <= upper){
		bdiff=0;
		diff=0;
	}else if(lower > val){
		bdiff=1;
		tmp=abs(mean-val);
		diff=k*tmp;
	}else if(val > upper){
		bdiff=1;
		tmp=abs(val-mean);
		diff=k*tmp;
	}

	/*
	else if(lower > val){
		bdiff=1;
		tmp=abs(lower-val);
		if(tmp < 1) diff=k*tmp;
		else diff=k*sqrt(tmp*tmp);
	}else if(val > upper){
		bdiff=1;
		tmp=abs(val-upper);
		if(tmp < 1) diff=k*tmp;
		else diff=k*sqrt(tmp*tmp);
	}
	*/
}

void RNASimulation::getoutenetorangle(double val, double mean, double lower, double upper, double k, double &diff, int &bdiff){
	double tmp1, tmp2;

	if((lower < upper) && (lower <= val) && (val <= upper)){
		bdiff=0;
		diff=0;
	}
	else if((upper < lower) & ((lower <= val && val <= 360) | (0 <= val && val <= upper))){
		bdiff=0;
		diff=0;
	}
	else if(lower > val){
		bdiff=1;
		tmp1=lower-val;
		tmp2=val+360-upper;
		if(tmp1 < tmp2) diff=k*abs(tmp1);
		else diff=k*abs(tmp2);
	}else if(val > upper){
		bdiff=1;
		tmp1=val-upper;
		tmp2=lower+360-val;
		if(tmp1 < tmp2) diff=k*abs(tmp1);
		diff=k*abs(tmp2);
	}
}

void RNASimulation::getoutenetorsion(double val, double mean, vector<double> boundary, vector<double> energy, vector<int> binary, double k, double &diff, int &bdiff){
	double difftmp;
	int binidx=getbin(val, boundary);
	bdiff=binary[binidx];
	diff=energy[binidx];
	//if(bdiff == 1){
		//difftmp=abs(val-mean);
		//if(difftmp > 180) diff=360-difftmp;
		//else diff=difftmp;
		//diff=k*sqrt(diff*diff);
	//}else{
		//diff=0;
	//}
}


void RNASimulation::getoutene(double val, double mean, vector<double> boundary, vector<double> energy, vector<int> binary, double &diff, int &bdiff){
	int binidx=getbin(val, boundary);
	bdiff=binary[binidx];
	diff=bdiff*abs(val-mean);
}

void RNASimulation::getoutlier(double val, vector<double> boundary, vector<double> energy, vector<int> binary, double &diff, int &bdiff){
	int binidx=getbin(val, boundary);
	diff=energy[binidx];
	bdiff=binary[binidx];
}

int RNASimulation::getbin(double val, vector<double> boundary){
	for(int i=0;i<boundary.size();i++){
		if(val < boundary[i]) return i;
	}
	return boundary.size()-1;
}

void RNASimulation::resetflag(vector<vector<bool>> &mapflagrpp, vector<vector<bool>> &mapbbflagrpp){
	for(int i=0;i<mapflagrpp.size();i++){
		for(int j=0;j<mapflagrpp[i].size();j++){
			mapflagrpp[i][j]=false;
			mapflagrpp[j][i]=false;

			mapbbflagrpp[i][j]=false;
			mapbbflagrpp[j][i]=false;
		}
	}
}

double RNASimulation::fcalcrmsdenergy(vector<Residueatoms> decstr,int numseq, Energyterms &lp,vector<vector<Rescontact>> &rpp, vector<Basecontact> &bprpp, vector<Stackcontact> &bsrpp, vector<vector<int>> &nmaprpp, vector<vector<vector<double>>> &distmaprpp, vector<vector<vector<double>>> &ntpairdistrpp, vector<vector<vector<double>>> &adjpairdistrpp, vector<double> &distenerpp, vector<vector<double>> &multidistenerpp, vector<Phycontact> &contactenerpp, vector<vector<bool>> &mapflagrpp, vector<vector<bool>> &mapbbflagrpp, bool flaginit)
{
	double trmsd=0;	
	bool disp = false;
	map<int, map<string, double>> enecoeff={
		//{0, {{"All", 0.1}, {"initdist", 5.0}, {"multidist", coefftemplate}, {"contact8A", 0}, {"contact20A", 0}, {"clash", 5}, {"adjpairlength", 0}, {"adjpairlength bin", 0}, {"ntpairlength", 0}, {"ntpairlength bin", 0}, {"paired torsion", 10}, {"paired bb torsion", 10}, {"multi torsion", 0}, {"multi torsion bin", 0}, {"clustered torsion", 0}, {"clustered torsion bin", 0}, {"bond length", 100}, {"bond length bin", 100}, {"bond angle", 50}, {"bond angle bin", 50}, {"bond torsion", 0}, {"bond torsion bin", 0}, {"bp", 10}, {"bpout", 10}, {"bphyd", 10.0}, {"bphydout", 10.0}, {"stacking", 10}, {"stackingout", 10}, {"pucker", 0}, {"pucker bin", 0}, {"oxhbond", 0}}}, //1.0001
		//{1, {{"All", 0.1}, {"initdist", 2.0}, {"multidist", coefftemplate}, {"contact8A", 0}, {"contact20A", 0}, {"clash", 150}, {"adjpairlength", 0}, {"adjpairlength bin", 0}, {"ntpairlength", 0}, {"ntpairlength bin", 0}, {"paired torsion", 30}, {"paired bb torsion", 30}, {"multi torsion", 0}, {"multi torsion bin", 0}, {"clustered torsion", 0}, {"clustered torsion bin", 0}, {"bond length", 100}, {"bond length bin", 100}, {"bond angle", 5}, {"bond angle bin", 5}, {"bond torsion", 1}, {"bond torsion bin", 1}, {"bp", 40}, {"bpout", 40}, {"bphyd", 200.0}, {"bphydout", 200.0}, {"stacking", 500}, {"stackingout", 500}, {"pucker", 0}, {"pucker bin", 0}, {"oxhbond", 0}}}
		
		//{0, {{"All", 0.1}, {"initdist", 5.0}, {"multidist", coefftemplate}, {"singlenttemp", 50}, {"contact8A", 0}, {"contact20A", 0}, {"clash", 5}, {"adjpairlength", 0}, {"adjpairlength bin", 0}, {"ntpairlength", 0}, {"ntpairlength bin", 0}, {"paired torsion", 10}, {"paired bb torsion", 10}, {"multi torsion", 0}, {"multi torsion bin", 0}, {"clustered torsion", 0}, {"clustered torsion bin", 0}, {"bond length", 100}, {"bond length bin", 100}, {"bond angle", 50}, {"bond angle bin", 50}, {"bond torsion", 0}, {"bond torsion bin", 0}, {"bp", 10}, {"bpout", 10}, {"bphyd", 10.0}, {"bphydout", 10.0}, {"stacking", 10}, {"stackingout", 10}, {"pucker", 0}, {"pucker bin", 0}, {"oxhbond", 0}, {"bptmp", 50}, {"stacktmp", 50}}}, //1.0001
		//{1, {{"All", 0.1}, {"initdist", 2.0}, {"multidist", coefftemplate}, {"singlenttemp", 0}, {"contact8A", 0}, {"contact20A", 0}, {"clash", 200}, {"adjpairlength", 0}, {"adjpairlength bin", 0}, {"ntpairlength", 0}, {"ntpairlength bin", 0}, {"paired torsion", 0}, {"paired bb torsion", 30}, {"multi torsion", 0}, {"multi torsion bin", 0}, {"clustered torsion", 0}, {"clustered torsion bin", 0}, {"bond length", 200}, {"bond length bin", 200}, {"bond angle", 5}, {"bond angle bin", 5}, {"bond torsion", 0}, {"bond torsion bin", 0}, {"bp", 40}, {"bpout", 40}, {"bphyd", 200.0}, {"bphydout", 200.0}, {"stacking", 0}, {"stackingout", 0}, {"pucker", 0}, {"pucker bin", 0}, {"oxhbond", 0}, {"bptmp", 50}, {"stacktmp", 80}}}

		{0, {{"All", 0.1}, {"initdist", 5.0}, {"multidist", coefftemplate}, {"singlenttemp", 50}, {"contact8A", 0}, {"contact20A", 0}, {"clash", 5}, {"adjpairlength", 0}, {"adjpairlength bin", 0}, {"ntpairlength", 0}, {"ntpairlength bin", 0}, {"paired torsion", 10}, {"paired bb torsion", 10}, {"multi torsion", 0}, {"multi torsion bin", 0}, {"clustered torsion", 0}, {"clustered torsion bin", 0}, {"bond length", 100}, {"bond length bin", 100}, {"bond angle", 50}, {"bond angle bin", 50}, {"bond torsion", 0}, {"bond torsion bin", 0}, {"bp", 10}, {"bpout", 10}, {"bphyd", 10.0}, {"bphydout", 10.0}, {"stacking", 10}, {"stackingout", 10}, {"pucker", 0}, {"pucker bin", 0}, {"oxhbond", 0}, {"bptmp", 50}, {"stacktmp", 50}}}, //1.0001
		{1, {{"All", 0.1}, {"initdist", 2.0}, {"multidist", coefftemplate}, {"singlenttemp", 0}, {"contact8A", 0}, {"contact20A", 0}, {"clash", 200}, {"adjpairlength", 0}, {"adjpairlength bin", 0}, {"ntpairlength", 0}, {"ntpairlength bin", 0}, {"paired torsion", 0}, {"paired bb torsion", 30}, {"multi torsion", 0}, {"multi torsion bin", 0}, {"clustered torsion", 0}, {"clustered torsion bin", 0}, {"bond length", 150}, {"bond length bin", 150}, {"bond angle", 5}, {"bond angle bin", 5}, {"bond torsion", 0}, {"bond torsion bin", 0}, {"bp", 50}, {"bpout", 50}, {"bphyd", 200.0}, {"bphydout", 200.0}, {"stacking", 0}, {"stackingout", 0}, {"pucker", 0}, {"pucker bin", 0}, {"oxhbond", 0}, {"bptmp", 50}, {"stacktmp", 50}}}

	};

	//for(i=0;i<20;i++) eneval[i]=0;
	for(int i=0;i<40;i++) enelist[i]=0;
	//ps.tor2dih(decstr,numseq);

	//if(enecoeff[enetype]["initFAPE"] > 0 || enecoeff[enetype]["multiFAPE"] > 0 || enecoeff[enetype]["tempFAPE"] > 0) getfape(decstr,numseq,faperpp,flaginit);

	//必ずdistにdistenergysingleとdistenergymultipleより前に計算すること
	if(enecoeff[enetype]["clash"] > 0 | enecoeff[enetype]["oxhbond"] > 0 | enecoeff[enetype]["initdist"] > 0 | enecoeff[enetype]["multidist"] > 0){
		//calculation of clash loss
		fenergyclash(decstr,numseq,enelist[2],enelist[12],enelist[34],lp,rpp,nmaprpp,distmaprpp,ntpairdistrpp,adjpairdistrpp,mapflagrpp,mapbbflagrpp,flaginit);
		enelist[2]=enecoeff[enetype]["All"]*enecoeff[enetype]["clash"]*enelist[2];
		enelist[12]=enecoeff[enetype]["All"]*enecoeff[enetype]["clash"]*enelist[12];
		enelist[34]=enecoeff[enetype]["All"]*enecoeff[enetype]["oxhbond"]*enelist[34];
		trmsd+=enelist[2];
		trmsd+=enelist[12];
		trmsd+=enelist[34];
		
		if(disp){
			cout<<"clash"<<endl;
			cout<<enelist[2]<<endl;
			cout<<enelist[12]<<endl;

			cout<<"oxhbond"<<endl;
			cout<<enelist[34]<<endl;
		}
	}

	if(enecoeff[enetype]["initdist"] > 0){
		distenergysingle(decstr,numseq,enelist[0],lp,distmaprpp,distenerpp,mapbbflagrpp,flaginit);
		enelist[0]=enecoeff[enetype]["All"]*enecoeff[enetype]["initdist"]*enelist[0];
		trmsd+=enelist[0];

		if(disp){
			cout<<"Distance(init)"<<endl;
			cout<<enelist[0]<<endl;
		}
	}

	if(enecoeff[enetype]["multidist"] > 0){
		distenergymultiple(decstr,numseq,enelist[1],lp,distmaprpp,multidistenerpp,mapbbflagrpp,flaginit);
		enelist[1]=enecoeff[enetype]["All"]*enecoeff[enetype]["multidist"]*enelist[1];
		trmsd+=enelist[1];

		if(disp){
			cout<<"Distance(template)"<<endl;
			cout<<enelist[1]<<endl;
		}
	}

	if(enecoeff[enetype]["singlenttemp"] > 0){
		fenergysinglenttemplate(decstr,numseq,enelist[35],lp,mapflagrpp,flaginit);
		enelist[35]=enecoeff[enetype]["All"]*enecoeff[enetype]["singlenttemp"]*enelist[35];
		trmsd+=enelist[35];

		if(disp){
			cout<<"Single template"<<endl;
			cout<<enelist[35]<<endl;
		}
	}

	if(enecoeff[enetype]["contact8A"] > 0 | enecoeff[enetype]["contact20A"] > 0){
		contactenergysingle(decstr,numseq,enelist[8],enelist[18],lp,distmaprpp,contactenerpp,mapbbflagrpp,flaginit);
		enelist[8]=enecoeff[enetype]["All"]*enecoeff[enetype]["contact8A"]*enelist[8];
		enelist[18]=enecoeff[enetype]["All"]*enecoeff[enetype]["contact20A"]*enelist[18];

		trmsd+=enelist[8];
		trmsd+=enelist[18];

		if(disp){
			cout<<"Contact"<<endl;
			cout<<enelist[8]<<endl;
			cout<<enelist[18]<<endl;
		}
	}

	//必ずfenergydihedralより上に定義
	if(enecoeff[enetype]["multi torsion"] > 0){
		fenergymultidihedral(decstr,numseq,enelist[9],enelist[19],lp,flaginit);
		enelist[9]=enecoeff[enetype]["All"]*enecoeff[enetype]["multi torsion bin"]*enelist[9];
		enelist[19]=enecoeff[enetype]["All"]*enecoeff[enetype]["multi torsion"]*enelist[19];

		trmsd+=enelist[19];

		if(disp){
			cout<<"Multi torsion"<<endl;
			cout<<enelist[9]<<endl;
			cout<<enelist[19]<<endl;
		}
	}

	//必ずfenergydihedralより上に定義
	if(enecoeff[enetype]["paired bb torsion"] > 0){
		fenergypairbbdihedral(decstr,numseq,enelist[31],lp,flaginit);
		enelist[31]=enecoeff[enetype]["All"]*enecoeff[enetype]["paired bb torsion"]*enelist[31];

		trmsd+=enelist[31];

		if(disp){
			cout<<"Paired bb torsion"<<endl;
			cout<<enelist[31]<<endl;
		}
	}

	//必ずfenergydihedralより上に定義
	if(enecoeff[enetype]["paired torsion"] > 0){
		fenergypairdihedral(decstr,numseq,enelist[30],lp,flaginit);
		enelist[30]=enecoeff[enetype]["All"]*enecoeff[enetype]["paired torsion"]*enelist[30];

		trmsd+=enelist[30];

		if(disp){
			cout<<"Paired torsion"<<endl;
			cout<<enelist[30]<<endl;
		}
	}

	//必ずfenergydihedralより上に定義
	if(enecoeff[enetype]["clustered torsion"] > 0){
		fenergyclsdihedral(decstr,numseq,enelist[22],enelist[23],lp,flaginit);
		enelist[22]=enecoeff[enetype]["All"]*enecoeff[enetype]["clustered torsion bin"]*enelist[22];
		enelist[23]=enecoeff[enetype]["All"]*enecoeff[enetype]["clustered torsion"]*enelist[23];

		trmsd+=enelist[23];
		if(disp){
			cout<<"Clustered torsion"<<endl;
			cout<<enelist[22]<<endl;
			cout<<enelist[23]<<endl;
		}
	}

	//必ずfenergyclashより上に定義
	/*
	if(enecoeff[enetype]["adjpairlength"] > 0){
		fenergyadjpairlengths(decstr,numseq,enelist[26],enelist[27],lp,adjpairdistrpp,mapflagrpp,flaginit);
		enelist[26]=enecoeff[enetype]["All"]*enecoeff[enetype]["adjpairlength bin"]*enelist[26];
		enelist[27]=enecoeff[enetype]["All"]*enecoeff[enetype]["adjpairlength"]*enelist[27];

		trmsd+=enelist[26];
		trmsd+=enelist[27];

		if(disp){
			cout<<"Adjpairlength"<<endl;
			cout<<enelist[26]<<endl;
			cout<<enelist[27]<<endl;
		}
	}

	//必ずfenergyclashより上に定義
	if(enecoeff[enetype]["ntpairlength"] > 0){
		fenergyntpairlengths(decstr,numseq,enelist[28],enelist[29],lp,ntpairdistrpp,mapflagrpp,flaginit);
		enelist[28]=enecoeff[enetype]["All"]*enecoeff[enetype]["ntpairlength bin"]*enelist[28];
		enelist[29]=enecoeff[enetype]["All"]*enecoeff[enetype]["ntpairlength"]*enelist[29];

		trmsd+=enelist[28];
		trmsd+=enelist[29];

		if(disp){
			cout<<"Nt pair length"<<endl;
			cout<<enelist[28]<<endl;
			cout<<enelist[29]<<endl;
		}
	}
	

	//必ずfenergydihedralより上に定義
	if(enecoeff[enetype]["pucker"] > 0){
		fenergypucker(decstr, numseq, enelist[24], enelist[25],lp,mapflagrpp,flaginit);
		enelist[24]=enecoeff[enetype]["All"]*enecoeff[enetype]["pucker bin"]*enelist[24];
		enelist[25]=enecoeff[enetype]["All"]*enecoeff[enetype]["pucker"]*enelist[25];
		//trmsd+=enelist[24];
		trmsd+=enelist[25];

		if(disp){
			cout<<"Pucker"<<endl;
			cout<<enelist[24]<<endl;
			cout<<enelist[25]<<endl;
		}
	}
	*/

	if(enecoeff[enetype]["bond length"] > 0){
		fenergybondlength(decstr,numseq,enelist[3],enelist[13],lp,flaginit);
		enelist[3]=enecoeff[enetype]["All"]*enecoeff[enetype]["bond length bin"]*enelist[3];
		enelist[13]=enecoeff[enetype]["All"]*enecoeff[enetype]["bond length"]*enelist[13];
		trmsd+=enelist[3];
		trmsd+=enelist[13];
		
		if(disp){
			cout<<"Bond length"<<endl;
			cout<<enelist[3]<<endl;
			cout<<enelist[13]<<endl;
		}
	}

	if(enecoeff[enetype]["bond angle"] > 0){
		fenergybondangle(decstr,numseq,enelist[4],enelist[14],lp,flaginit);
		enelist[4]=enecoeff[enetype]["All"]*enecoeff[enetype]["bond angle bin"]*enelist[4];
		enelist[14]=enecoeff[enetype]["All"]*enecoeff[enetype]["bond angle"]*enelist[14];
		trmsd+=enelist[4];
		trmsd+=enelist[14];
		
		if(disp){
			cout<<"Bond angle"<<endl;
			cout<<enelist[4]<<endl;
			cout<<enelist[14]<<endl;
		}
	}

	if(enecoeff[enetype]["bond torsion"] > 0){
		fenergydihedral(decstr,numseq,enelist[5],enelist[15],lp,flaginit);
		enelist[5]=enecoeff[enetype]["All"]*enecoeff[enetype]["bond torsion bin"]*enelist[5];
		enelist[15]=enecoeff[enetype]["All"]*enecoeff[enetype]["bond torsion"]*enelist[15];
		//trmsd+=enelist[5];
		trmsd+=enelist[15];
		
		if(disp){
			cout<<"Bond torsion"<<endl;
			cout<<enelist[5]<<endl;
			cout<<enelist[15]<<endl;
		}
	}

	if(enecoeff[enetype]["bp"] > 0){
		fenergybasepair(decstr, numseq, enelist[6], enelist[16], lp, bprpp, mapflagrpp, flaginit);
		enelist[6]=enecoeff[enetype]["All"]*enecoeff[enetype]["bpout"]*enelist[6];
		enelist[16]=enecoeff[enetype]["All"]*enecoeff[enetype]["bp"]*enelist[16];
		trmsd+=enelist[6];
		trmsd+=enelist[16];

		if(disp){
			cout<<"Base pair"<<endl;
			cout<<enelist[6]<<endl;
			cout<<enelist[16]<<endl;
		}
	}

	if(enecoeff[enetype]["bphyd"] > 0){
		fenergybasepairhyd(decstr, numseq, enelist[32], enelist[33], lp, bprpp, mapflagrpp, flaginit);
		enelist[32]=enecoeff[enetype]["All"]*enecoeff[enetype]["bphydout"]*enelist[32];
		enelist[33]=enecoeff[enetype]["All"]*enecoeff[enetype]["bphyd"]*enelist[33];
		trmsd+=enelist[32];
		trmsd+=enelist[33];

		if(disp){
			cout<<"base pair (hydrogen)"<<endl;
			cout<<enelist[32]<<endl;
			cout<<enelist[33]<<endl;
		}
	}

	if(enecoeff[enetype]["stacking"] > 0){
		fenergystacking(decstr, numseq, enelist[7], enelist[17], lp, bsrpp, mapflagrpp, flaginit);
		enelist[7]=enecoeff[enetype]["All"]*enecoeff[enetype]["stackingout"]*enelist[7];
		enelist[17]=enecoeff[enetype]["All"]*enecoeff[enetype]["stacking"]*enelist[17];
		trmsd+=enelist[7];
		trmsd+=enelist[17];

		if(disp){
			cout<<"Base stacking"<<endl;
			cout<<enelist[7]<<endl;
			cout<<enelist[17]<<endl;
		}
	}
        
	if(enecoeff[enetype]["bptmp"] > 0){
		fenergynoncbasepairtemp(decstr,numseq,enelist[37],lp,bprpp,mapflagrpp,flaginit);
		enelist[37]=enecoeff[enetype]["All"]*enecoeff[enetype]["bptmp"]*enelist[37];
		trmsd+=enelist[37];

		if(disp){
			cout<<"Base pair template"<<endl;
			cout<<enelist[37]<<endl;
		}
	}

	if(enecoeff[enetype]["stacktmp"] > 0){
		fenergystacktemplate(decstr, numseq, enelist[36], lp, bsrpp, mapflagrpp, flaginit);
		enelist[36]=enecoeff[enetype]["All"]*enecoeff[enetype]["stacktmp"]*enelist[36];
		trmsd+=enelist[36];

		if(disp){
			cout<<"Stacking template"<<endl;
			cout<<enelist[36]<<endl;
		}
	}

	//if(enecoeff[enetype]["hbond"] > 0)
	//{
		//fenergyhbond(decstr,numseq,enelist[8],enelist[18], lp, hbondrpp, flaginit);
		//enelist[8]=enecoeff[enetype]["All"]*enecoeff[enetype]["hbondout"]*baseenecoeff*enelist[8];
		//enelist[18]=enecoeff[enetype]["All"]*enecoeff[enetype]["hbond"]*enelist[18]/numseq;
		//trmsd+=enelist[8];
		//trmsd+=enelist[18];

		//cout<<"hbond"<<endl;
		//cout<<enelist[8]<<endl;
		//cout<<enelist[18]<<endl;
	//}

	resetflag(mapflagrpp,mapbbflagrpp);

	return trmsd;
}









