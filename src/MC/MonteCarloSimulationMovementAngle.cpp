
#include "MonteCarloSimulation.hpp"

void MonteCarloSimulation::setangflap(double trandpos, double base, int angtype)
{
	int i,j;
	int atom1, atom2, atom3;

	angleflag[trandpos][angtype]=1;
	
	for(i=trandpos;i<numbb;i++) ntflag[i]=1;
	//ntflag[trandpos+1]=1;
	
	if(angtype == 5){
		if(trandpos < numbb-1){
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

	//for(int k=0;k<numbb;k++){
	for(i=0;i<AddHydLengthParamsList[base].size();i++)
	{
		hydaddlengthflag[trandpos][i]=1;
	}

	for(i=0;i<AddHydAngleParamsList[base].size();i++)
	{
		hydaddangleflag[trandpos][i]=1;
	}
	//}

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

void MonteCarloSimulation::setangntflap(double trandpos, double base, int angtype)
{
	int i,j;
	int atom1, atom2, atom3;

	angleflag[trandpos][angtype]=1;
	
	for(i=trandpos;i<numbb;i++) ntflag[i]=1;
	//ntflag[trandpos+1]=1;

	if(angtype == 5){
		if(trandpos < numbb-1){
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

	//for(int k=0;k<numbb;k++){
	for(i=0;i<AddHydLengthParamsList[base].size();i++)
	{
		hydaddlengthflag[trandpos][i]=1;
	}

	for(i=0;i<AddHydAngleParamsList[base].size();i++)
	{
		hydaddangleflag[trandpos][i]=1;
	}
	//}
	
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

double MonteCarloSimulation::getmovangval(double val, double mean, double std, double k){
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

bool MonteCarloSimulation::selecttargetang(int &trandpos, int &angtype, double r)
{
	double randval=Random();

	if(r > randval){
		int i, j, idx;
		int res1, res2, res3;
		int atom1, atom2, atom3;
		char base;
		vector<vector<double>> prob;
		prob.resize(numbb);

		base=refinestr[0].aaa;
		prob[0].resize(AccAtomParamNum[base][3], 0);
		res1=0; res2=0; res3=0;
		atom1=ParamsList[base][2][1];atom2=ParamsList[base][2][2];atom3=ParamsList[base][2][3];
		prob[0][2] += angleenergyatoms[res1][atom1];
		prob[0][2] += angleenergyatoms[res2][atom2];
		prob[0][2] += angleenergyatoms[res3][atom3];

		for(i=3;i<AccAtomParamNum[base][3];i++){
			if(i == 10 || i == 11){
				atom1=FIRSTOPParamsList[base][i][1]; atom2=FIRSTOPParamsList[base][i][2]; atom3=FIRSTOPParamsList[base][i][3];
			}else{
				if(0 < i && i < AccAtomParamNum[base][2]){
					atom1=ParamsList[base][i][1];atom2=ParamsList[base][i][2];atom3=ParamsList[base][i][3];
				}else{
					idx=i-AccAtomParamNum[base][2];
					atom1=BaseParamsList[base][idx][1]; atom2=BaseParamsList[base][idx][2]; atom3=BaseParamsList[base][idx][3];
				}
			}
			prob[0][i] += angleenergyatoms[res1][atom1];
			prob[0][i] += angleenergyatoms[res2][atom2];
			prob[0][i] += angleenergyatoms[res3][atom3];
		}

		for(i=1;i<numbb;i++){
			base=refinestr[i].aaa;
			prob[i].resize(AccAtomParamNum[base][3], 0);
			
			res1=i-1; res2=i-1; res3=i;
			atom1=ParamsList[base][0][1];atom2=ParamsList[base][0][2];atom3=ParamsList[base][0][3];
			prob[i][0] += angleenergyatoms[res1][atom1];
			prob[i][0] += angleenergyatoms[res2][atom2];
			prob[i][0] += angleenergyatoms[res3][atom3];

			res1=i-1; res2=i; res3=i;
			atom1=ParamsList[base][1][1];atom2=ParamsList[base][1][2];atom3=ParamsList[base][1][3];
			prob[i][1] += angleenergyatoms[res1][atom1];
			prob[i][1] += angleenergyatoms[res2][atom2];
			prob[i][1] += angleenergyatoms[res3][atom3];

			res1=i-1; res2=i; res3=i;
			atom1=ParamsList[base][11][1];atom2=ParamsList[base][11][2];atom3=ParamsList[base][11][3];
			prob[i][11] += angleenergyatoms[res1][atom1];
			prob[i][11] += angleenergyatoms[res2][atom2];
			prob[i][11] += angleenergyatoms[res3][atom3];

			for(j=2;j<11;j++){
				res1=i; res2=i; res3=i;
				atom1=ParamsList[base][j][1];atom2=ParamsList[base][j][2];atom3=ParamsList[base][j][3];
				prob[i][j] += angleenergyatoms[res1][atom1];
				prob[i][j] += angleenergyatoms[res2][atom2];
				prob[i][j] += angleenergyatoms[res3][atom3];
			}

			for(j=12;j<AccAtomParamNum[base][3];j++){
				res1=i; res2=i; res3=i;

				if(0 < j && j < AccAtomParamNum[base][2]){
					atom1=ParamsList[base][j][1];atom2=ParamsList[base][j][2];atom3=ParamsList[base][j][3];
				}else{
					idx=j-AccAtomParamNum[base][2];
					atom1=BaseParamsList[base][idx][1]; atom2=BaseParamsList[base][idx][2]; atom3=BaseParamsList[base][idx][3];
				}
				//atom2=ParamsList[base][j][1];atom3=ParamsList[base][j][2];atom4=ParamsList[base][j][3];
				prob[i][j] += angleenergyatoms[res1][atom1];
				prob[i][j] += angleenergyatoms[res2][atom2];
				prob[i][j] += angleenergyatoms[res3][atom3];
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
					angtype=j;
					return true;
				}
				currentenen=nextene;
			}
		}
	}else{
		trandpos=Random()*numbb;
		char base=refinestr[trandpos].aaa;
		angtype=int(AccAtomParamNum[base][3]*Random());
		return true;
	}

	return false;
}

bool MonteCarloSimulation::fmovementang(double oldenergy, double &newenergy)
{
	bool flagtor;
	//int trandpos=Random()*numbb;
	//char base=refinestr[trandpos].aaa;
	//int angtype=int(AccAtomParamNum[base][3]*Random());

	int trandpos, angtype;
	bool selflag = selecttargetang(trandpos, angtype, 0.8);

	if(selflag == 0){
		newenergy=oldenergy+1000000000;
		return false;
	}


	char base=refinestr[trandpos].aaa;

	//cout<<"=============="<<endl;
	//cout<<trandpos << "-" << angtype<<endl;

	double mean, std, upper, lower;
	int idx;

	int uptype=int(Random()*2);

	StrCopy(refinestr, refinestr_temp);
	resetfrag(0);

	if(trandpos==0 && angtype == 10){
		if(uptype == 0){
			refinestr_temp[trandpos].ang[angtype]=getmovangval(refinestr_temp[trandpos].ang[angtype], MCOPINNERANGPARAM[base][0][0], MCOPINNERANGPARAM[base][0][1], 0.1);
		}else{
			refinestr_temp[trandpos].ang[angtype]=getrandomnewparam(MCOPINNERANGPARAM[base][0][2], MCOPINNERANGPARAM[base][0][3]);
		}
	}else if(trandpos==0 && angtype == 11){
		if(uptype == 0){
			refinestr_temp[trandpos].ang[angtype]=getmovangval(refinestr_temp[trandpos].ang[angtype], MCOPINNERANGPARAM[base][1][0], MCOPINNERANGPARAM[base][1][1], 0.1);
		}else{
			refinestr_temp[trandpos].ang[angtype]=getrandomnewparam(MCOPINNERANGPARAM[base][1][2], MCOPINNERANGPARAM[base][1][3]);
		}
	}else{
		if(0 <= angtype && angtype < AccAtomParamNum[base][2]){
			mean = MCINNERANGPARAM[base][angtype][0];
			std = MCINNERANGPARAM[base][angtype][1];
			lower = MCINNERANGPARAM[base][angtype][2];
			upper = MCINNERANGPARAM[base][angtype][3];
	   	}else if(AccAtomParamNum[base][2] <= angtype && angtype < AccAtomParamNum[base][3]){
			idx = angtype - AccAtomParamNum[base][2];
			mean = MCBASEINNERANGPARAM[base][idx][0];
			std = MCBASEINNERANGPARAM[base][idx][1];
			lower = MCBASEINNERANGPARAM[base][idx][2];
			upper = MCBASEINNERANGPARAM[base][idx][3];
	   	}
		
		if(uptype == 0){
			refinestr_temp[trandpos].ang[angtype]=getmovangval(refinestr_temp[trandpos].ang[angtype], mean, std, 0.1);
		}else{
			refinestr_temp[trandpos].ang[angtype]=getrandomnewparam(lower, upper);
		}
	}

	if(angtype < 6){
		flagtor=param2coordbb(refinestr_temp,trandpos,numbb);
		param2coordsugar(refinestr_temp,trandpos,numbb);
		param2coordbase(refinestr_temp,trandpos,numbb);
		param2coordhydrogen(refinestr_temp,trandpos,numbb);
		if(!flagtor) printf("tor2str wrong in ang %d\n",trandpos);

		if(angtype == 5 || angtype == 4){
			int end = min(trandpos+2, numbb);
			coord2paramadd(refinestr_temp,trandpos,end);
			coord2addparamhydrogen(refinestr_temp,trandpos,end);
		}else{
			int end = min(trandpos+1, numbb);
			coord2paramadd(refinestr_temp,trandpos,end);
			coord2addparamhydrogen(refinestr_temp,trandpos,end);
		}
		setangflap(trandpos, base, angtype);
	}else{
		param2coordsugar(refinestr_temp,trandpos,trandpos+1);
		param2coordbase(refinestr_temp,trandpos,trandpos+1);
		param2coordhydrogen(refinestr_temp,trandpos,trandpos+1);
		coord2paramadd(refinestr_temp,trandpos,trandpos+1);
		coord2addparamhydrogen(refinestr_temp,trandpos,trandpos+1);
		setangntflap(trandpos, base, angtype);
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

