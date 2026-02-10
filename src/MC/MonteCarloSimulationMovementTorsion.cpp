
#include "MonteCarloSimulation.hpp"

void MonteCarloSimulation::settorflap(double trandpos, double base, int tortype)
{
	int i,j;
	int atom1, atom2, atom3;

	torsionflag[trandpos][tortype]=1;

	for(i=trandpos;i<numbb;i++) ntflag[i]=1;
	
	if(tortype == 5 && trandpos+1 < numbb){
		addlengthflag[trandpos+1][0]=1;
		ntflag[trandpos+1] = 1;
	}

	for(i=0;i<AddLengthParamsList[base].size();i++){
		addlengthflag[trandpos][i]=1;
	}
	
	for(i=0;i<AddAngleParamsList[base].size();i++)
	{
		addangleflag[trandpos][i]=1;
	}

	for(i=0;i<AddHydLengthParamsList[base].size();i++){
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

void MonteCarloSimulation::settorntflap(double trandpos, double base, int tortype)
{
	int i,j;
	int atom1, atom2, atom3;

	torsionflag[trandpos][tortype]=1;

	ntflag[trandpos] = 1;

	if(tortype == 5 && trandpos+1 < numbb){
		addlengthflag[trandpos+1][0]=1;
		ntflag[trandpos+1] = 1;
	}
	
	for(i=0;i<AddLengthParamsList[base].size();i++){
		addlengthflag[trandpos][i]=1;
	}
	
	for(i=0;i<AddAngleParamsList[base].size();i++)
	{
		addangleflag[trandpos][i]=1;
	}

	for(i=0;i<AddHydLengthParamsList[base].size();i++){
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

bool MonteCarloSimulation::selecttargettor(int &trandpos, int &tortype, double r)
{
	double randval=Random();

	if(r > randval){
		int i, j, idx;
		int res1, res2, res3, res4;
		int atom1, atom2, atom3, atom4;
		char base;
		vector<vector<double>> prob;
		prob.resize(numbb);

		base=refinestr[0].aaa;
		prob[0].resize(AccAtomParamNum[base][3], 0);
		res1=0; res2=0; res3=0; res4=0;
		atom1=ParamsList[base][2][0];atom2=ParamsList[base][2][1];atom3=ParamsList[base][2][2];atom4=ParamsList[base][2][3];
		prob[0][2] += torsionenergyatoms[res1][atom1];
		prob[0][2] += torsionenergyatoms[res2][atom2];
		prob[0][2] += torsionenergyatoms[res3][atom3];
		prob[0][2] += torsionenergyatoms[res4][atom4];

		for(i=3;i<AccAtomParamNum[base][3];i++){
			if(i == 10 || i == 11){
				atom1=FIRSTOPParamsList[base][i][0]; atom2=FIRSTOPParamsList[base][i][1]; atom3=FIRSTOPParamsList[base][i][2]; atom4=FIRSTOPParamsList[base][i][3];
			}else{
				if(0 < i && i < AccAtomParamNum[base][2]){
					atom1=ParamsList[base][i][0];atom2=ParamsList[base][i][1];atom3=ParamsList[base][i][2];atom4=ParamsList[base][i][3];
				}else{
					idx=i-AccAtomParamNum[base][2];
					atom1=BaseParamsList[base][idx][0]; atom2=BaseParamsList[base][idx][1]; atom3=BaseParamsList[base][idx][2]; atom4=BaseParamsList[base][idx][3];
				}
			}
			prob[0][i] += torsionenergyatoms[res1][atom1];
			prob[0][i] += torsionenergyatoms[res2][atom2];
			prob[0][i] += torsionenergyatoms[res3][atom3];
			prob[0][i] += torsionenergyatoms[res4][atom4];
		}

		for(i=1;i<numbb;i++){
			base=refinestr[i].aaa;
			prob[i].resize(AccAtomParamNum[base][3], 0);
			
			res1=i-1; res2=i-1; res3=i-1; res4=i;
			atom1=ParamsList[base][0][0];atom2=ParamsList[base][0][1];atom3=ParamsList[base][0][2];atom4=ParamsList[base][0][3];
			prob[i][0] += torsionenergyatoms[res1][atom1];
			prob[i][0] += torsionenergyatoms[res2][atom2];
			prob[i][0] += torsionenergyatoms[res3][atom3];
			prob[i][0] += torsionenergyatoms[res4][atom4];

			res1=i-1; res2=i-1; res3=i; res4=i;
			atom1=ParamsList[base][1][0];atom2=ParamsList[base][1][1];atom3=ParamsList[base][1][2];atom4=ParamsList[base][1][3];
			prob[i][1] += torsionenergyatoms[res1][atom1];
			prob[i][1] += torsionenergyatoms[res2][atom2];
			prob[i][1] += torsionenergyatoms[res3][atom3];
			prob[i][1] += torsionenergyatoms[res4][atom4];

			res1=i-1; res2=i; res3=i; res4=i;
			atom1=ParamsList[base][2][0];atom2=ParamsList[base][2][1];atom3=ParamsList[base][2][2];atom4=ParamsList[base][2][3];
			prob[i][1] += torsionenergyatoms[res1][atom1];
			prob[i][1] += torsionenergyatoms[res2][atom2];
			prob[i][1] += torsionenergyatoms[res3][atom3];
			prob[i][1] += torsionenergyatoms[res4][atom4];

			res1=i-1; res2=i; res3=i; res4=i;
			atom1=ParamsList[base][10][0];atom2=ParamsList[base][10][1];atom3=ParamsList[base][10][2];atom4=ParamsList[base][10][3];
			prob[i][10] += torsionenergyatoms[res1][atom1];
			prob[i][10] += torsionenergyatoms[res2][atom2];
			prob[i][10] += torsionenergyatoms[res3][atom3];
			prob[i][10] += torsionenergyatoms[res4][atom4];

			res1=i; res2=i-1; res3=i; res4=i;
			atom1=ParamsList[base][11][0];atom2=ParamsList[base][11][1];atom3=ParamsList[base][11][2];atom4=ParamsList[base][11][3];
			prob[i][11] += torsionenergyatoms[res1][atom1];
			prob[i][11] += torsionenergyatoms[res2][atom2];
			prob[i][11] += torsionenergyatoms[res3][atom3];
			prob[i][11] += torsionenergyatoms[res4][atom4];

			for(j=2;j<10;j++){
				res1=i; res2=i; res3=i; res4=i;
				atom1=ParamsList[base][j][0];atom2=ParamsList[base][j][1];atom3=ParamsList[base][j][2];atom3=ParamsList[base][j][3];
				prob[i][j] += torsionenergyatoms[res1][atom1];
				prob[i][j] += torsionenergyatoms[res2][atom2];
				prob[i][j] += torsionenergyatoms[res3][atom3];
				prob[i][j] += torsionenergyatoms[res4][atom4];
			}

			for(j=10;j<AccAtomParamNum[base][3];j++){
				res1=i; res2=i; res3=i; res4=i;

				if(0 < j && j < AccAtomParamNum[base][2]){
					atom1=ParamsList[base][j][0];atom2=ParamsList[base][j][1];atom3=ParamsList[base][j][2];atom4=ParamsList[base][j][3];
				}else{
					idx=j-AccAtomParamNum[base][2];
					atom1=BaseParamsList[base][idx][0]; atom2=BaseParamsList[base][idx][1]; atom3=BaseParamsList[base][idx][2];atom3=BaseParamsList[base][idx][3];
				}
				//atom2=ParamsList[base][j][1];atom3=ParamsList[base][j][2];atom4=ParamsList[base][j][3];
				prob[i][j] += torsionenergyatoms[res1][atom1];
				prob[i][j] += torsionenergyatoms[res2][atom2];
				prob[i][j] += torsionenergyatoms[res3][atom3];
				prob[i][j] += torsionenergyatoms[res4][atom4];
			}
		}

		double tenergy=0;
		for(i=0;i<numbb;i++){
			base=refinestr[i].aaa;
			for(j=0;j<AccAtomParamNum[base][3];j++){
				tenergy+=prob[i][j];
			}
		}

		if(tenergy > -0.000001){
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
				if(currentenen >= wpos && wpos > nextene){
					trandpos=i;
					tortype=j;
					return true;
				}
				currentenen=nextene;
			}
		}
	}else{
		trandpos=Random()*numbb;
		char base=refinestr[trandpos].aaa;
		tortype=int(AccAtomParamNum[base][3]*Random());
		return true;
	}
	return false;
}

double MonteCarloSimulation::sampletorsion(vector<double> &logdistr, double logtotal, vector<double> &center, double span)
{
	int i;
	double pval, sval, rval;
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
	pval=center[binidx] + (span/2)*(2*Random()-1);
	return pval;
}

bool MonteCarloSimulation::fmovementtor(double oldenergy, double &newenergy)
{
	bool flagtor;
	int trandpos, tortype;
	bool selflag = selecttargettor(trandpos, tortype, 0.8);

	if(selflag == 0){
		newenergy=oldenergy+1000000000;
		return false;
	}

	char base=refinestr[trandpos].aaa;

	//cout<<"=============="<<endl;
	//cout<<trandpos << "-" << tortype<<endl;

	double mean, std, upper, lower;
	int idx;

	StrCopy(refinestr, refinestr_temp);
	resetfrag(0);

	if(trandpos==0 && tortype == 10){
		refinestr_temp[trandpos].dih[tortype]=sampletorsion(MCOPTORSIONDISTRIBUTION[base][0].energy, MCOPTORSIONDISTRIBUTION[base][0].totalenergy, MCOPTORSIONDISTRIBUTION[base][0].center, MCOPTORSIONDISTRIBUTION[base][0].span);
	}else if(trandpos==0 && tortype == 11){
		refinestr_temp[trandpos].dih[tortype]=sampletorsion(MCOPTORSIONDISTRIBUTION[base][1].energy, MCOPTORSIONDISTRIBUTION[base][1].totalenergy, MCOPTORSIONDISTRIBUTION[base][1].center, MCOPTORSIONDISTRIBUTION[base][1].span);
	}else{
		if(0 <= tortype && tortype < AccAtomParamNum[base][2]){
			refinestr_temp[trandpos].dih[tortype]=sampletorsion(MCTORSIONDISTRIBUTION[base][tortype].energy, MCTORSIONDISTRIBUTION[base][tortype].totalenergy, MCTORSIONDISTRIBUTION[base][tortype].center, MCTORSIONDISTRIBUTION[base][tortype].span);
	   	}else if(AccAtomParamNum[base][2] <= tortype && tortype < AccAtomParamNum[base][3]){
			idx = tortype - AccAtomParamNum[base][2];
			refinestr_temp[trandpos].dih[tortype]=sampletorsion(MCBASETORSIONDISTRIBUTION[base][idx].energy, MCBASETORSIONDISTRIBUTION[base][idx].totalenergy, MCBASETORSIONDISTRIBUTION[base][idx].center, MCBASETORSIONDISTRIBUTION[base][idx].span);
	   	}
	}

	if(tortype < 6){
		flagtor=param2coordbb(refinestr_temp,trandpos,numbb);
		param2coordsugar(refinestr_temp,trandpos,numbb);
		param2coordbase(refinestr_temp,trandpos,numbb);
		param2coordhydrogen(refinestr_temp,trandpos,numbb);
		if(!flagtor) printf("tor2str wrong in tor %d\n",trandpos);

		if(tortype == 5 || tortype == 4){
			int end = min(trandpos+2, numbb);
			coord2paramadd(refinestr_temp,trandpos,end);
			coord2addparamhydrogen(refinestr_temp,trandpos,trandpos+1);
		}else{
			int end = min(trandpos+1, numbb);
			coord2paramadd(refinestr_temp,trandpos,end);
			coord2addparamhydrogen(refinestr_temp,trandpos,end);
		}
		settorflap(trandpos, base, tortype);
	}else{
		int end = min(trandpos+1, numbb);
		param2coordsugar(refinestr_temp,trandpos,end);
		param2coordbase(refinestr_temp,trandpos,end);
		param2coordhydrogen(refinestr_temp,trandpos,end);
		coord2addparamhydrogen(refinestr_temp,trandpos,end);
		coord2paramadd(refinestr_temp,trandpos,end);
		settorntflap(trandpos, base, tortype);
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

