
#include "MonteCarloSimulation.hpp"

void MonteCarloSimulation::setcombtorflap(double trandpos, double base)
{
	int i,j;
	int atom1, atom2, atom3;

	for(i=trandpos;i<numbb;i++) ntflag[i] = 1;

	for(j=0;j<torcombmovidx.size();j++){
		torsionflag[trandpos][torcombmovidx[j]]=1;
	}
	
	if(trandpos < numbb-1){
		ntflag[trandpos+1] = 1;
		addlengthflag[trandpos+1][0]=1;
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


void MonteCarloSimulation::selecttargetcombtor(int &trandpos, int &tortype, double r)
{
	double randval=Random();

	if(r > randval){
		int i, j, idx;
		int res1, res2, res3, res4;
		int atom1, atom2, atom3, atom4;
		char base;
		vector<vector<double>> prob;
		prob.resize(numbb);

		for(i=0;i<numbb;i++){
			base=refinestr[i].aaa;

			prob[i].resize(torcombmovidx.size(), 0);
			for(j=0;j<torcombmovidx.size();j++){
				res1=i; res2=i; res3=i; res4=i;
				atom1=ParamsList[base][torcombmovidx[j]][0];atom2=ParamsList[base][torcombmovidx[j]][1];atom3=ParamsList[base][torcombmovidx[j]][2];atom3=ParamsList[base][torcombmovidx[j]][3];
				prob[i][j] += torsionenergyatoms[res1][atom1];
				prob[i][j] += torsionenergyatoms[res2][atom2];
				prob[i][j] += torsionenergyatoms[res3][atom3];
				prob[i][j] += torsionenergyatoms[res4][atom4];
			}
		}

		double tenergy=0;
		for(i=0;i<numbb;i++){
			base=refinestr[i].aaa;
			for(j=0;j<torcombmovidx.size();j++){
				tenergy+=prob[i][j];
			}
		}
		
		double wpos=Random()*tenergy;
		double currentenen=0;
		double nextene;
		tenergy=0;
		for(i=0;i<numbb;i++){
			base=refinestr[i].aaa;
			for(j=0;j<torcombmovidx.size();j++){
				nextene=currentenen+prob[i][j];
				if(currentenen >= wpos && wpos > nextene){
					trandpos=i;
					tortype=torcombmovidx[j];
					return;
				}
				currentenen=nextene;
			}
		}
	}else{
		trandpos=Random()*numbb;
		char base=refinestr[trandpos].aaa;
		tortype=torcombmovidx[int(torcombmovidx.size()*Random())];
	}
}

/*
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
	cout<<pval<<endl;
	return pval;
}
*/

double MonteCarloSimulation::samplesecondtor2bin(vector<double> &prob, double totalprob)
{
	double rval=totalprob*(Random()-EPS);
	double sval=0;
	int binidx;

	for(int i=0;i<prob.size();i++){
		if(sval<rval & rval<=sval+prob[i]){
			binidx=i;
			break;
		}
		sval+=prob[i];
	}
	if(binidx == prob.size()) binidx=prob.size()-1;
	return binidx;
}

double MonteCarloSimulation::samplesecondtor2(vector<double> &center, double span, int binidx)
{	
	return center[binidx] + (span/2)*(2*Random()-1);
}

bool MonteCarloSimulation::fmovementcombtor(double oldenergy, double &newenergy)
{
	bool flagtor;
	int trandpos, tor1, tor2, paramidx, paramorder, binnum1, binnum2;
	selecttargetcombtor(trandpos, tor1, 0.8);
	char base=refinestr[trandpos].aaa;

	double mean, std, upper, lower;
	int idx;

	StrCopy(refinestr, refinestr_temp);
	resetfrag(0);

	refinestr_temp[trandpos].dih[tor1]=sampletorsion(MCTORSIONDISTRIBUTION[base][tor1].energy, MCTORSIONDISTRIBUTION[base][tor1].totalenergy, MCTORSIONDISTRIBUTION[base][tor1].center, MCTORSIONDISTRIBUTION[base][tor1].span);
	
	vector<vector<int>> torsion_combs = torcombmov[base][tor1];
	for(int j=0;j<torsion_combs.size();j++){
		tor2 = torsion_combs[j][0];
		paramidx = torsion_combs[j][1];
		paramorder = torsion_combs[j][2];

		if(paramorder == 0){
			binnum1 = TORSIONPAIRSPARAM[base][paramidx].binnum1;
			binnum2 = TORSIONPAIRSPARAM[base][paramidx].binnum2;

			int binidx=getbin(refinestr_temp[trandpos].dih[tor1], TORSIONPAIRSPARAM[base][paramidx].boundary1);
			
			vector<double> tor2_probs(binnum2);
			double totalprob=0;
			for(int i=0;i<binnum2;i++){
				tor2_probs[i] = TORSIONPAIRSPARAM[base][paramidx].energy[binidx*binnum2+i];
				totalprob+=tor2_probs[i];
			}
			int abinidx=samplesecondtor2bin(tor2_probs, totalprob);

			refinestr_temp[trandpos].dih[tor2]=samplesecondtor2(TORSIONPAIRSPARAM[base][paramidx].center2, TORSIONPAIRSPARAM[base][paramidx].span2, abinidx);
		}else{
			binnum1 = TORSIONPAIRSPARAM[base][paramidx].binnum1;
			binnum2 = TORSIONPAIRSPARAM[base][paramidx].binnum2;

			int binidx=getbin(refinestr_temp[trandpos].dih[tor1], TORSIONPAIRSPARAM[base][paramidx].boundary2);

			vector<double> tor2_probs(binnum1);
			double totalprob=0;
			for(int i=0;i<binnum1;i++){
				tor2_probs[i] = TORSIONPAIRSPARAM[base][paramidx].energy[i*binnum2+binidx];
				totalprob+=tor2_probs[i];
			}
			int abinidx=samplesecondtor2bin(tor2_probs, totalprob);
			refinestr_temp[trandpos].dih[tor2]=samplesecondtor2(TORSIONPAIRSPARAM[base][paramidx].center1, TORSIONPAIRSPARAM[base][paramidx].span1, abinidx);
		}
	}

	flagtor=param2coordbb(refinestr_temp,trandpos,numbb);
	param2coordsugar(refinestr_temp,trandpos,numbb);
	param2coordbase(refinestr_temp,trandpos,numbb);
	param2coordhydrogen(refinestr_temp,trandpos,numbb);
	if(!flagtor) printf("tor2str wrong in tor %d\n",trandpos);

	int end = min(trandpos+2, numbb);
	coord2paramadd(refinestr_temp,trandpos,end);
	coord2addparamhydrogen(refinestr_temp,trandpos,end);
	
	setcombtorflap(trandpos, base);
	
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

