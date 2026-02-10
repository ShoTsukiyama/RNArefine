
#include "MonteCarloSimulation.hpp"

void MonteCarloSimulation::setlengflap(double trandpos, double base, int lengtype)
{
	int i,j;
	int atom1, atom2, atom3;

	lengthflag[trandpos][lengtype]=1;

	for(i=trandpos;i<numbb;i++) ntflag[i]=1;
	
	if(lengtype == 5){
		if(trandpos < numbb-1){
			ntflag[trandpos+1] = 1;
			addlengthflag[trandpos+1][0]=1;
		}
	}

	for(i=0;i<AddLengthParamsList[base].size();i++){
		addlengthflag[trandpos][i]=1;
	}
	
	for(i=0;i<AddAngleParamsList[base].size();i++)
	{
		addangleflag[trandpos][i]=1;
	}

	for(i=0;i<AddHydLengthParamsList[base].size();i++)
	{
		hydaddlengthflag[trandpos][i]=1;
	}

	for(i=0;i<AddHydAngleParamsList[base].size();i++)
	{
		hydaddangleflag[trandpos][i]=1;
	}

	for(i=0;i<numbb;i++)
	{
		for(j=i;j<numbb;j++)
		{
			if(j<trandpos || i>trandpos+1){}
			else
			{
				mapbbflag[i][j]=1;
				mapbbflag[j][i]=1;
				mapflag[i][j]=1;
				mapflag[j][i]=1;
			}
		}
	}
}

void MonteCarloSimulation::setlengntflap(double trandpos, double base, int lengtype)
{
	int i,j;
	int atom1, atom2, atom3;

	lengthflag[trandpos][lengtype]=1;

	for(i=trandpos;i<numbb;i++) ntflag[i]=1;
	
	if(lengtype == 5){
		if(trandpos < numbb-1){
			ntflag[trandpos+1] = 1;
			addlengthflag[trandpos+1][0]=1;
		}
	}

	for(i=0;i<AddLengthParamsList[base].size();i++){
		addlengthflag[trandpos][i]=1;
	}
	
	for(i=0;i<AddAngleParamsList[base].size();i++)
	{
		addangleflag[trandpos][i]=1;
	}

	for(i=0;i<AddHydLengthParamsList[base].size();i++)
	{
		hydaddlengthflag[trandpos][i]=1;
	}

	for(i=0;i<AddHydAngleParamsList[base].size();i++)
	{
		hydaddangleflag[trandpos][i]=1;
	}
	
	for(i=0;i<numbb;i++)
	{
		for(j=i;j<numbb;j++)
		{
			if(j<trandpos || i>trandpos+1){}
			else
			{
				mapflag[i][j]=1;
				mapflag[j][i]=1;
			}
		}
	}
}

double MonteCarloSimulation::getmovlengval(double val, double mean, double std, double k){
	double newval;
	double deg;
	if(val <= mean){
		deg=mean-val;
		if(deg > std) deg = std;
		newval = val + (deg)*k*Random();
	}else{
		deg=val-mean;
		if(deg > std) deg = std;
		newval = val - (deg)*k*Random();
	}

	return newval;
}

bool MonteCarloSimulation::selecttargetleng(int &trandpos, int &lengtype, double r)
{
	double randval=Random();

	if(r > randval){
		int i, j, idx;
		int res2, res3;
		int atom2, atom3;
		char base;
		vector<vector<double>> prob;
		prob.resize(numbb);

		base=refinestr[0].aaa;
		prob[0].resize(AccAtomParamNum[base][3], 0);
		res2=0; res3=0;
		atom2=ParamsList[base][2][2];atom3=ParamsList[base][2][3];
		prob[0][2] += lengthenergyatoms[res2][atom2];
		prob[0][2] += lengthenergyatoms[res3][atom3];

		for(i=3;i<AccAtomParamNum[base][3];i++){
			if(i == 10 || i == 11){
				atom2=FIRSTOPParamsList[base][i][2]; atom3=FIRSTOPParamsList[base][i][3];
			}else{
				if(0 < i && i < AccAtomParamNum[base][2]){
					atom2=ParamsList[base][i][2];atom3=ParamsList[base][i][3];
				}else{
					idx=i-AccAtomParamNum[base][2];
					atom2=BaseParamsList[base][idx][2]; atom3=BaseParamsList[base][idx][3];
				}
			}
			prob[0][i] += lengthenergyatoms[res2][atom2];
			prob[0][i] += lengthenergyatoms[res3][atom3];
		}

		for(i=1;i<numbb;i++){
			base=refinestr[i].aaa;
			prob[i].resize(AccAtomParamNum[base][3], 0);
			
			res2=i-1; res3=i;
			atom2=ParamsList[base][0][2];atom3=ParamsList[base][0][3];
			prob[i][0] += lengthenergyatoms[res2][atom2];
			prob[i][0] += lengthenergyatoms[res3][atom3];

			res2=i; res3=i;
			atom2=ParamsList[base][1][2];atom3=ParamsList[base][1][3];
			prob[i][1] += lengthenergyatoms[res2][atom2];
			prob[i][1] += lengthenergyatoms[res3][atom3];

			for(j=2;j<AccAtomParamNum[base][3];j++){
				res2=i; res3=i;

				if(0 < j && j < AccAtomParamNum[base][2]){
					atom2=ParamsList[base][j][2];atom3=ParamsList[base][j][3];
				}else{
					idx=j-AccAtomParamNum[base][2];
					atom2=BaseParamsList[base][idx][2]; atom3=BaseParamsList[base][idx][3];
				}
				prob[i][j] += lengthenergyatoms[res2][atom2];
				prob[i][j] += lengthenergyatoms[res3][atom3];
			}
		}

		double tenergy=0;
		for(i=0;i<numbb;i++){
			base=refinestr[i].aaa;
			for(j=0;j<AccAtomParamNum[base][3];j++){
				tenergy+=prob[i][j];
			}
		}

		if(tenergy < 0.000001){
			return false;
		}
		
		double wpos=Random()*tenergy;
		double currentenen=0;
		double nextene;
		tenergy=0;
		for(i=0;i<numbb;i++){
			base=refinestr[i].aaa;
			for(j=0;j<AccAtomParamNum[base][3];j++){
				nextene=currentenen+prob[i][j];
				if(currentenen <= wpos && wpos < nextene){
					trandpos=i;
					lengtype=j;
					return true;
				}
				currentenen=nextene;
			}
		}
	}else{
		trandpos=Random()*numbb;
		char base=refinestr[trandpos].aaa;
		lengtype=int(AccAtomParamNum[base][3]*Random());
		return true;
	}
	return false;
}


bool MonteCarloSimulation::fmovementleng(double oldenergy, double &newenergy)
{
	bool flagtor;
	int trandpos, lengtype;
	bool selflag = selecttargetleng(trandpos, lengtype, 0.8);

	if(selflag == 0){
		newenergy=oldenergy+1000000000;
		return false;
	}

	char base=refinestr[trandpos].aaa;

	//cout<<"=============="<<endl;
	//cout<<trandpos << "-" << lengtype<<endl;

	double mean, std, upper, lower;
	int idx;

	int uptype=int(Random()*2);

	StrCopy(refinestr, refinestr_temp);
	resetfrag(0);

	if(trandpos==0 && lengtype == 10){
		if(uptype == 0){
			refinestr_temp[trandpos].len[lengtype]=getmovlengval(refinestr_temp[trandpos].len[lengtype], MCOPLENGTHPARAM[base][0][0], MCOPLENGTHPARAM[base][0][1], 0.1);
		}else{
			refinestr_temp[trandpos].len[lengtype]=getrandomnewparam(MCOPLENGTHPARAM[base][0][2], MCOPLENGTHPARAM[base][0][3]);
		}
	}else if(trandpos==0 && lengtype == 11){
		if(uptype == 0){
			refinestr_temp[trandpos].len[lengtype]=getmovlengval(refinestr_temp[trandpos].len[lengtype], MCOPLENGTHPARAM[base][1][0], MCOPLENGTHPARAM[base][1][1], 0.1);
		}else{
			refinestr_temp[trandpos].len[lengtype]=getrandomnewparam(MCOPLENGTHPARAM[base][1][2], MCOPLENGTHPARAM[base][1][3]);
		}
	}else{
		if(0 <= lengtype && lengtype < AccAtomParamNum[base][2]){
			mean = MCLENGTHPARAM[base][lengtype][0];
			std = MCLENGTHPARAM[base][lengtype][1];
			lower = MCLENGTHPARAM[base][lengtype][2];
			upper = MCLENGTHPARAM[base][lengtype][3];
	   	}else if(AccAtomParamNum[base][2] <= lengtype && lengtype < AccAtomParamNum[base][3]){
			idx = lengtype - AccAtomParamNum[base][2];
			mean = MCBASELENGTHPARAM[base][idx][0];
			std = MCBASELENGTHPARAM[base][idx][1];
			lower = MCBASELENGTHPARAM[base][idx][2];
			upper = MCBASELENGTHPARAM[base][idx][3];
	   	}

		if(uptype == 0){
			refinestr_temp[trandpos].len[lengtype]=getmovlengval(refinestr_temp[trandpos].len[lengtype], mean, std, 0.1);
		}else{
			refinestr_temp[trandpos].len[lengtype]=getrandomnewparam(lower, upper);
		}
	}

	if(lengtype < 6){
		flagtor=param2coordbb(refinestr_temp,trandpos,numbb);
		param2coordsugar(refinestr_temp,trandpos,numbb);
		param2coordbase(refinestr_temp,trandpos,numbb);
		param2coordhydrogen(refinestr_temp,trandpos,numbb);
		if(!flagtor) printf("tor2str wrong in leng %d\n",trandpos);

		if(lengtype == 5 || lengtype == 4){
			int end = min(trandpos+2, numbb);
			coord2paramadd(refinestr_temp,trandpos,end);
			coord2addparamhydrogen(refinestr_temp,trandpos,trandpos+1);
		}else{
			int end = min(trandpos+1, numbb);
			coord2paramadd(refinestr_temp,trandpos,end);
			coord2addparamhydrogen(refinestr_temp,trandpos,end);
		}
		setlengflap(trandpos, base, lengtype);
	}else{
		param2coordsugar(refinestr_temp,trandpos,trandpos+1);
		param2coordbase(refinestr_temp,trandpos,trandpos+1);
		param2coordhydrogen(refinestr_temp,trandpos,trandpos+1);
		coord2paramadd(refinestr_temp,trandpos,trandpos+1);
		coord2addparamhydrogen(refinestr_temp,trandpos,trandpos+1);
		setlengntflap(trandpos, base, lengtype);
	}

	newenergy = calcenergymain();
	bool upflag=metropolisUpdate(oldenergy, newenergy);
	if(upflag){
		StrCopy(refinestr_temp, refinestr);
		initializeatomenergy();
		updateenergies();
	}

	/*
	cout<<"===================="<<endl;
	cout<<"QqqqqqqqqqqqqqqqqQ"<<endl;
	checkupdate();
	cout<<"PppppppppppppppppP"<<endl;
	resetfrag(1);
	double newenergy_tmp = calcenergymain();

	if(abs(newenergy - newenergy_tmp)/max(abs(newenergy), abs(newenergy_tmp)) > 0.002){
		cout<<"Energy Error"<<endl;
		cout<<newenergy << "-" << newenergy_tmp<<endl;
		exit(1);
	}else{
		cout<<"OK"<<endl;
	}
	*/

	return upflag;
}

