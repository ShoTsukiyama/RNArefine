// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

bool RNASimulation::fmovementgradientangle(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc)
{
	int i,j,k,inext;
	int istart,iend;
	int trandpos, trandposed, trandposst, span;
	int basei, basej;
	double trandang;
	Point3d rp;
	double mati[9];
	int wpos;
	char aat1, aat2, aat;
	bool idxtmp1, idxtmp2;
	int transbase;
	int total;
	int numiter;
	double rmax, rnorm;
	int mtype1, mtype2;
	vector<int> ntparamidx;
	vector<int> adjparamidx;
	bool flagdone;
	double randval=Random();
	bool exitOuterLoop = false;
	int mtype;
	int rtype;
	int pairwiseidx;
	vector<bool> movenodeinfo;
	vector<pair<int, int>> edgeinfo;
	vector<int> edgetype;
	vector<int> paraminfo;

	double learning_rate;
	int epochs;
	double tolerance;

	//double anglecoeff=0.01;
	//double initcoeff=1.0;
	//double paircoeff=10.0;

	double anglecoeff=0.1;
	double initcoeff=1.0;

	//learning_rate=0.0001;
	//epochs=10000;
	//tolerance=0.000001;

	learning_rate=0.00001;
	epochs=5000;
	tolerance=0.00001;

	randval=Random();

	int totwrong1=lp1.nn[CLASH_ENERGY_ST]*2;
	int totwrong2=lp1.nn[BASE_ENERGY_ST+0] + lp1.nn[BASE_ENERGY_ST+1];

	if(totwrong1>0 && randval<0.4){
		wpos=int(Random()*totwrong1);

		i=0;
		for(trandpos=0;trandpos<numseq;trandpos++)
		{
			inext=i+lp1.indn[CLASH_ENERGY_ST][trandpos];
			if(i<=wpos && wpos<inext){
				break;
			}
			i=inext;
		}
	}else if(totwrong2>0 && randval>0.4 && randval<0.8){
		wpos=int(Random()*totwrong2);
		mtype=1;

		i=0;
		for(trandpos=0;trandpos<BASEPAIRS.size();trandpos++)
		{
			basei=BASEPAIRS[trandpos].basepidx.first;
			basej=BASEPAIRS[trandpos].basepidx.second;

			pairwiseidx=bpmapid2idx[basei*numseq+basej];

			inext=i+bpenergy[pairwiseidx].pcw[0]+bpenergy[pairwiseidx].pcw[1];
			if(i<=wpos && wpos<inext){
				break;
			}
			i=inext;
		}

		randval=Random();
		if(randval > 0.5) trandpos=basei;
		else trandpos=basej;
	}
	
	else{
		trandpos=int((numseq)*Random());
	}

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

	//cout<<to_string(trandposst) + "-" + to_string(trandposed)<<endl;

	vector<pair<int,int>> pairsf; 
	for(i=trandposst;i<trandposed;i++){
		if(bppairednt[i] == -1){
			pairsf.push_back({i,-1});
		}else if(bppairednt[i] < 1){
			break;
		}else if(bppairednt[i] > numseq-1){
			break;
		}else{
			pairsf.push_back({i,bppairednt[i]});
		}
	}

	vector<pair<int,int>> pairsb;
	for(i=trandposed-1;i>=trandposst;i--){
		if(bppairednt[i] == -1){
			pairsb.push_back({i,-1});
		}else if(bppairednt[i] < 1){
			break;
		}else if(bppairednt[i] > numseq-1){
			break;
		}else{
			pairsb.push_back({i,bppairednt[i]});
		}
	}

	vector<pair<int,int>> pairs;
	if(pairsf.size() >= pairsb.size()){
		pairs=pairsf;
	}else{
		for (i=pairsb.size()-1; i>=0; i--) {
			pairs.push_back(pairsb[i]);
		}
	}

	copystr(decstr, numseq);

	vector<vector<int>> lmpparamidx;
	setlmpbpparam(restr, pairs, lmpparamidx);

	int threshiter=span*200;
	if(threshiter>3000) threshiter=3000;

	/*
	for(i=0;i<pairs.size();i++)
	{
		basei = pairs[i].first;
		basej = pairs[i].second;

		//cout<<basei<<endl;
		//cout<<basej<<endl;

		
		if(basej == -1){
			for(j=0;j<PAIRLENGTHTARGETMOVNUM;j++){
				//cout<<to_string(restr[basei].ptv[j].x) + "-" + to_string(restr[basei].ptv[j].y) + "-" + to_string(restr[basei].ptv[j].z)<<endl;
				//rmax=0.01+0.01*Random();
				rmax=0.25+0.25*Random();
				rp=bf.ranv(rmax);
				restr[basei].ptv[j].x+=rp.x;
				restr[basei].ptv[j].y+=rp.y;
				restr[basei].ptv[j].z+=rp.z;
				//cout<<to_string(restr[basei].ptv[j].x) + "-" + to_string(restr[basei].ptv[j].y) + "-" + to_string(restr[basei].ptv[j].z)<<endl;
			}
		}else{
			for(j=0;j<PAIRLENGTHTARGETMOVNUM;j++){
				//cout<<to_string(restr[basei].ptv[j].x) + "-" + to_string(restr[basei].ptv[j].y) + "-" + to_string(restr[basei].ptv[j].z)<<endl;
				//rmax=0.01+0.01*Random();
				rmax=0.25+0.25*Random();
				rp=bf.ranv(rmax);
				restr[basei].ptv[j].x+=rp.x;
				restr[basei].ptv[j].y+=rp.y;
				restr[basei].ptv[j].z+=rp.z;
				//cout<<to_string(restr[basei].ptv[j].x) + "-" + to_string(restr[basei].ptv[j].y) + "-" + to_string(restr[basei].ptv[j].z)<<endl;
			}
			for(j=0;j<PAIRLENGTHTARGETMOVNUM;j++){
				//cout<<to_string(restr[basej].ptv[j].x) + "-" + to_string(restr[basej].ptv[j].y) + "-" + to_string(restr[basej].ptv[j].z)<<endl;
				
				//rmax=0.01+0.01*Random();
				rmax=0.25+0.25*Random();
				rp=bf.ranv(rmax);
				restr[basej].ptv[j].x+=rp.x;
				restr[basej].ptv[j].y+=rp.y;
				restr[basej].ptv[j].z+=rp.z;
				//cout<<to_string(restr[basej].ptv[j].x) + "-" + to_string(restr[basej].ptv[j].y) + "-" + to_string(restr[basej].ptv[j].z)<<endl;
				
			}
		}
	}
	*/
	

	getgraphinfoLMP(pairs, numseq, lmpparamidx, movenodeinfo, edgeinfo, edgetype, paraminfo);

	vector<int> movementresidues;
	map<int, vector<Point3s>> gradptv;
	map<int, vector<double>> initlengths;
	map<int, vector<Point3s>> initdists;

	map<int, vector<double>> angles;
	map<int, vector<double>> anglediffs;
	map<int, vector<bool>> angleouterflags;
	map<int, vector<double>> theta;
	map<int, vector<double>> normPQs;
	map<int, vector<double>> normRQs; 
	map<int, vector<double>> dotps;
	map<int, vector<Point3s>> vecPQs; 
	map<int, vector<Point3s>> vecRQs;
	int paramnum=AccAtomParamNum['A'][2]+AddAngleParamsList['A'].size();
	

	int initcount, anglecount;
	double mse_prev=0;
    double mse=0;
	double mseinit_prev=0;
	double mseangle_prev=0;
	double mseinit=0;
	double mseangle=0;

	initcount=0;
	for(i=0;i<pairs.size();i++)
	{
		basei = pairs[i].first; basej = pairs[i].second;
		if(basej == -1){
			movementresidues.push_back(basei);
		}else{
			movementresidues.push_back(basei);
			movementresidues.push_back(basej);
		}
	}

	initcount=0; anglecount=0;
	for(i=0;i<movementresidues.size();i++)
	{
		basei = movementresidues[i];
		initlengths[basei].resize(PAIRLENGTHTARGETMOVNUM, 0.0);
		initdists[basei].resize(PAIRLENGTHTARGETMOVNUM);
		calcmseinitlength(restr[basei].ptv, decstr[basei].ptv, initlengths[basei], initdists[basei], mseinit_prev, initcount, initcoeff);

		angles[basei].resize(paramnum, 0.0);
		anglediffs[basei].resize(paramnum, 0.0);
		angleouterflags[basei].resize(paramnum, false);
		theta[basei].resize(paramnum, 0.0);
		normPQs[basei].resize(paramnum, 0.0);
		normRQs[basei].resize(paramnum, 0.0);
		dotps[basei].resize(paramnum, 0.0);
		vecPQs[basei].resize(paramnum);
		vecRQs[basei].resize(paramnum);

		calcmseangles(restr, basei, angles[basei], anglediffs[basei], angleouterflags[basei], theta[basei], normPQs[basei], normRQs[basei], dotps[basei], vecPQs[basei], vecRQs[basei], mseangle_prev, anglecount, anglecoeff);
	}
	mseinit_prev=mseinit_prev/initcount;
	mseangle_prev=mseangle_prev/anglecount;
	mse_prev=mseinit_prev+mseangle_prev;

	//cout << "Initial" << " mse: " << mse_prev << " mse(pair): " << msepair_prev << " mse(init): " << mseinit_prev  << " mse(angle): " << mseangle_prev <<endl;

	/*
	for(i=0;i<movenodeinfo.size();i++){
		if(movenodeinfo[i]){
			//cout<<i<<endl;
			for(mtype=0;mtype<PAIRLENGTHTARGETMOVNUM;mtype++)
			{
				cout<< restr[i].ptv[mtype].x << "-" << restr[i].ptv[mtype].y << "-" << restr[i].ptv[mtype].z << endl;
			}
		}
	}
	*/

	for (int epoch = 0; epoch < epochs; ++epoch) {
		initializegrad(gradptv, movenodeinfo, mse);
		mseinit=0;
		mseangle=0;

		for(i=0;i<movementresidues.size();i++)
		{
			basei = movementresidues[i];
			calcgradientinitlength(restr[basei].ptv, decstr[basei].ptv, basei, gradptv, initlengths[basei], initdists[basei], initcount, initcoeff);
			calcgradientangles(restr, basei, gradptv, movenodeinfo, anglediffs[basei], angleouterflags[basei], theta[basei], normPQs[basei], normRQs[basei], dotps[basei], vecPQs[basei], vecRQs[basei], anglecount, anglecoeff);
		}

		for(i=0;i<movenodeinfo.size();i++){
			if(movenodeinfo[i]){
				for(mtype=0;mtype<PAIRLENGTHTARGETMOVNUM;mtype++)
				{
					restr[i].ptv[mtype].x -= learning_rate*gradptv[i][mtype].x;
					restr[i].ptv[mtype].y -= learning_rate*gradptv[i][mtype].y;
					restr[i].ptv[mtype].z -= learning_rate*gradptv[i][mtype].z;
				}
			}
		}

		initcount=0; anglecount=0;
		for(i=0;i<movementresidues.size();i++)
		{
			basei = movementresidues[i];
			initlengths[basei].resize(PAIRLENGTHTARGETMOVNUM, 0.0);
			initdists[basei].resize(PAIRLENGTHTARGETMOVNUM);
			calcmseinitlength(restr[basei].ptv, decstr[basei].ptv, initlengths[basei], initdists[basei], mseinit, initcount, initcoeff);


			angles[basei].resize(paramnum, 0.0);
			anglediffs[basei].resize(paramnum, 0.0);
			angleouterflags[basei].resize(paramnum, false);
			theta[basei].resize(paramnum, 0.0);
			normPQs[basei].resize(paramnum, 0.0);
			normRQs[basei].resize(paramnum, 0.0);
			dotps[basei].resize(paramnum, 0.0);
			vecPQs[basei].resize(paramnum);
			vecRQs[basei].resize(paramnum);

			calcmseangles(restr, basei, angles[basei], anglediffs[basei], angleouterflags[basei], theta[basei], normPQs[basei], normRQs[basei], dotps[basei], vecPQs[basei], vecRQs[basei], mseangle, anglecount, anglecoeff);
		}
		mseinit=mseinit/initcount;
		mseangle=mseangle/anglecount;
		mse=mseinit+mseangle;

        //cout << "epoch:" << epoch << " mse: " << mse << " mse(init): " << mseinit  << " mse(angle): " << mseangle << " prev mse: " << mse_prev << " prev mse(init): " << mseinit_prev << " prev mse(angle): " << mseangle_prev <<endl;
		//cout << "epoch:" << epoch << " mse: " << mse << " mse(init): " << mseinit  << " mse(angle): " << mseangle <<endl;


		 if(abs(mse_prev - mse) < tolerance){
			//cout << "Convergence reached at epoch " << epoch << ", MSE: " << mse << endl;
			/*
			for(i=0;i<movenodeinfo.size();i++){
				if(movenodeinfo[i]){
					//cout<<i<<endl;
					for(mtype=0;mtype<PAIRLENGTHTARGETMOVNUM;mtype++)
					{
						cout<< restr[i].ptv[mtype].x << "-" << restr[i].ptv[mtype].y << "-" << restr[i].ptv[mtype].z << endl;
					}
				}
			}
			*/
			break;
		 }

		mseangle_prev = mseangle;
		mseinit_prev = mseinit;
		mse_prev = mse;   
	}

	vector<bool> targetidxtmp(numseq, false);
	for(i=0;i<pairs.size();i++)
	{
		basei = pairs[i].first;
		basej = pairs[i].second;

		targetidxtmp[basei] = 1;
		if(basej != -1){
			targetidxtmp[basej] = 1;
		}
	}
	vector<pair<int,int>> continuousidx;
	int start = -1;

	for(i=0;i<numseq;i++){
		if (targetidxtmp[i] == 1) {
			if (start == -1) {
				start = i;
			}
		}else {
            if (start != -1) {
                continuousidx.push_back({start, i - 1});
                start = -1;
            }
        }
	}
	if (start != -1) {
        continuousidx.push_back({start, numseq - 1});
    }

	for(i = 0; i<continuousidx.size();i++){
		istart=continuousidx[i].first-1;
		iend=continuousidx[i].second+2;
		if(iend>=numseq) iend=numseq;
		//cout<<to_string(istart) + "-" + to_string(iend)<<endl;
		coord2parambb(restr,continuousidx[i].first,iend,1);
		coord2paramsugar(restr,continuousidx[i].first,iend,1);
		param2coordbase(restr,continuousidx[i].first,iend);
		param2coordhydrogen(restr,continuousidx[i].first,iend);

		coord2paramadd(restr,istart,iend);
		settmplmpflap(numseq, istart, iend);
	}

	newenergy=fcalcrmsdenergy(restr, numseq, lp2, rppmatbk, bpenergybk, bsenergybk, nmapbk, distmapbk, ntpairdistbk, adjpairdistbk, distenebk, multidistenebk, contactenebk, mapflagbk, mapbbflagbk, flagacc);

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

	//exit(1);
	
	bool upidx=metropolisUpdate(oldenergy, newenergy, Bt);

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












