// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////




void RNASimulation::settmpbpgradientlmpfullflap(double numseq, double trandpos1, double trandpos2)
{
	int i,j;
	int st, ed;

	for(i=trandpos1;i<trandpos2;i++){
		for(j=0;j<BLENGTHTARGET.size();j++){
			lp2.flap[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][i]=true;
		}

		for(j=0;j<INNERANGLETARGET.size();j++){
			lp2.flap[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][i]=true;
		}

		for(j=0;j<TORSIONANGLETARGET.size();j++){
			lp2.flap[TORSIONANGLE_ENERGY_ST+TORSIONANGLETARGET[j]][i]=true;
		}

		for(j=0;j<ADDBLENGTHTARGET.size();j++)
		{
			lp2.flap[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]][i]=true;
		}
		
		for(j=0;j<ADDINNERANGLETARGET.size();j++)
		{
			lp2.flap[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]][i]=true;
		}
	}
	
	for(i=0;i<numseq;i++)
	{
		for(j=i;j<numseq;j++)
		//for(j=i+1;j<numseq;j++)
		{
			if(j<trandpos1 || i>trandpos2){}
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

bool RNASimulation::fmovementgradientlength(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc)
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
	map<int, vector<Point3s>> snttemplatestr;

	double learning_rate;
	int epochs;
	double tolerance;
	int tmp;

	//double anglecoeff=0.01;
	//double initcoeff=1.0;
	//double paircoeff=10.0;

	//double initcoeff=50.0;
	//double paircoeffadj=5000.0;
	//double paircoeffnt=5000.0;
	//double paircoeffbp=0.3;
	//double anglecoeff=0;

	//double initcoeff=50.0;
	//double paircoeffadj=5000.0;
	//double paircoeffnt=5000.0;
	//double paircoeffbp=0.3;
	//double anglecoeff=0.1;

	double initcoeff;
	double tmpcoeff;
	double anglecoeff;
	double paircoeffadj;
	double paircoeffnt;
	double paircoeffbp;

	randval=Random();
	if(randval < 0.5){
		initcoeff=0;
		tmpcoeff=5000.0;
		anglecoeff=500;
		paircoeffadj=5000;
		paircoeffnt=5000;
		paircoeffbp=5.0;
	}else{
		initcoeff=0;
		tmpcoeff=1000.0;
		anglecoeff=0;
		paircoeffadj=5000;
		paircoeffnt=5000;
		paircoeffbp=5.0;
	}
	

	/*
	if(randval < 0.333){
		learning_rate=0.01;
		epochs=10000;
		tolerance=0.01;
	}else if(randval >= 0.333 & randval < 0.667){
		learning_rate=0.001;
		epochs=10000;
		tolerance=0.01;
	}else{
		learning_rate=0.0005;
		epochs=10000;
		tolerance=0.01;
	}
	*/

	//learning_rate=0.0001;
	//epochs=10000;
	//tolerance=0.000001;

	//learning_rate=0.0001;
	//epochs=1000;
	//tolerance=0.000001;

	//learning_rate=0.0001;
	//epochs=5000;
	//tolerance=0.000001;
	//learning_rate=0.00001;
	//epochs=5000;
	//tolerance=0.000001;

	learning_rate=0.000001;
	epochs=100+900*Random();
	tolerance=0.0000001;

	randval=Random();

	int totwrong1=lp1.nn[CLASH_ENERGY_ST]*2;
	int totwrong2=lp1.nn[BASE_ENERGY_ST+0] + lp1.nn[BASE_ENERGY_ST+1];
	int totwrong3=0;
	for(i=0;i<BLENGTHTARGET.size();i++){
		totwrong3+=lp1.nn[BLENGTH_ENERGY_ST+BLENGTHTARGET[i]];
		//totwrong3-=lp1.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[i]][0];
	}
	for(i=0;i<INNERANGLETARGET.size();i++){
		totwrong3+=lp1.nn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[i]];
		//totwrong3-=lp1.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[i]][0];
	}

	if(totwrong3>0 && randval<0.5){
		wpos=int(Random()*totwrong3-EPS);

		i=0;
		for(trandpos=0;trandpos<numseq;trandpos++)
		{
			for(j=0;j<BLENGTHTARGET.size();j++){
				inext=i+lp1.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][trandpos];
				if(i<=wpos<inext){
					exitOuterLoop = true;
					break;
				}
				i=inext;
			}

			for(j=0;j<INNERANGLETARGET.size();j++){
				inext=i+lp1.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][trandpos];
				if(i<=wpos<inext){
					exitOuterLoop = true;
					break;
				}
				i=inext;
			}
			if(exitOuterLoop) break;
		}
	}else if(totwrong1>0 && randval>=0.5 && randval<0.7){
		wpos=int(Random()*totwrong1-EPS);

		i=0;
		for(trandpos=0;trandpos<numseq;trandpos++)
		{
			inext=i+lp1.indn[CLASH_ENERGY_ST][trandpos];
			if(i<=wpos && wpos<inext){
				break;
			}
			i=inext;
		}
	}else if(totwrong2>0 && randval>=0.7 && randval<0.9){
		wpos=int(Random()*totwrong2-EPS);

		i=0;
		for(trandpos=0;trandpos<BASEPAIRS.size();trandpos++)
		{
			basei=BASEPAIRS[trandpos].basepidx.first;
			basej=BASEPAIRS[trandpos].basepidx.second;

			if(basei > basej){
				tmp=basei;
				basei=basej;
				basej=tmp;
			}
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
	}else{
		trandpos=int((numseq-2)*Random())+1;
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

	vector<pair<int,int>> pairsf; 
	for(i=trandposst;i<trandposed-1;i++){
		if(bppairednt[i] == -1){
			pairsf.push_back({i,-1});
		}else if(bppairednt[i] < 1){
			break;
		}else if(bppairednt[i] >= numseq-1){
			break;
		}else{
			pairsf.push_back({i,bppairednt[i]});
		}
	}

	vector<pair<int,int>> pairsb;
	for(i=trandposed-2;i>=trandposst;i--){
		if(bppairednt[i] == -1){
			pairsb.push_back({i,-1});
		}else if(bppairednt[i] < 1){
			break;
		}else if(bppairednt[i] >= numseq-1){
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

	//int threshiter=span*200;
	//if(threshiter>3000) threshiter=3000;

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
				rmax+=0.01+0.01*Random();
				//rmax=0.25+0.25*Random();
				rp=bf.ranv(rmax);
				//cout<< restr[basei].ptv[j].x << "-" << restr[basei].ptv[j].y << "-" << restr[basei].ptv[j].z << "-"<<endl;
				restr[basei].ptv[j].x+=rp.x;
				restr[basei].ptv[j].y+=rp.y;
				restr[basei].ptv[j].z+=rp.z;
				//cout<< restr[basei].ptv[j].x << "-" << restr[basei].ptv[j].y << "-" << restr[basei].ptv[j].z << "-"<<endl;
				//cout<<to_string(restr[basei].ptv[j].x) + "-" + to_string(restr[basei].ptv[j].y) + "-" + to_string(restr[basei].ptv[j].z)<<endl;
			}
		}else{
			for(j=0;j<PAIRLENGTHTARGETMOVNUM;j++){
				//cout<<to_string(restr[basei].ptv[j].x) + "-" + to_string(restr[basei].ptv[j].y) + "-" + to_string(restr[basei].ptv[j].z)<<endl;
				rmax+=0.01+0.01*Random();
				//rmax=0.25+0.25*Random();
				rp=bf.ranv(rmax);
				//cout<< restr[basei].ptv[j].x << "-" << restr[basei].ptv[j].y << "-" << restr[basei].ptv[j].z << "-"<<endl;
				restr[basei].ptv[j].x+=rp.x;
				restr[basei].ptv[j].y+=rp.y;
				restr[basei].ptv[j].z+=rp.z;
				//cout<< restr[basei].ptv[j].x << "-" << restr[basei].ptv[j].y << "-" << restr[basei].ptv[j].z << "-"<<endl;
				//cout<<to_string(restr[basei].ptv[j].x) + "-" + to_string(restr[basei].ptv[j].y) + "-" + to_string(restr[basei].ptv[j].z)<<endl;
			}
			for(j=0;j<PAIRLENGTHTARGETMOVNUM;j++){
				//cout<<to_string(restr[basej].ptv[j].x) + "-" + to_string(restr[basej].ptv[j].y) + "-" + to_string(restr[basej].ptv[j].z)<<endl;
				
				rmax+=0.01+0.01*Random();
				//rmax=0.25+0.25*Random();
				rp=bf.ranv(rmax);
				//cout<< restr[basej].ptv[j].x << "-" << restr[basej].ptv[j].y << "-" << restr[basej].ptv[j].z << "-"<<endl;
				restr[basej].ptv[j].x+=rp.x;
				restr[basej].ptv[j].y+=rp.y;
				restr[basej].ptv[j].z+=rp.z;
				//cout<< restr[basej].ptv[j].x << "-" << restr[basej].ptv[j].y << "-" << restr[basej].ptv[j].z << "-"<<endl;
				//cout<<to_string(restr[basej].ptv[j].x) + "-" + to_string(restr[basej].ptv[j].y) + "-" + to_string(restr[basej].ptv[j].z)<<endl;
				
			}
		}
	}
	*/

	getgraphinfoLMP(pairs, numseq, lmpparamidx, movenodeinfo, edgeinfo, edgetype, paraminfo);

	map<int, Residueatoms> snttemplates;
	for(i=0;i<movenodeinfo.size();i++){
		Residueatoms snttemplate;
		if(movenodeinfo[i]){
			//cout<<"============="<<endl;
			//cout<< i << "-" << restr[i].ptv[0].x << "-" << restr[i].ptv[0].y << "-" << restr[i].ptv[0].z<<endl;
			selectsnttemplate(restr[i], snttemplate);
			snttemplates[i]=snttemplate;
			//cout<< i << "-" << restr[i].ptv[0].x << "-" << restr[i].ptv[0].y << "-" << restr[i].ptv[0].z<<endl;
			//cout<< i << "-" << snttemplates[i].ptv[0].x << "-" << snttemplates[i].ptv[0].y << "-" << snttemplates[i].ptv[0].z<<endl;
		}
	}

	vector<int> movementresidues;
	vector<double> coeffs;
	vector<vector<vector<double>>> lengths;
	vector<vector<vector<double>>> diffs;
	vector<vector<vector<Point3s>>> dists;
	vector<vector<vector<bool>>> outlineridx;
	map<int, vector<Point3s>> gradptv;
	map<int, vector<double>> initlengths, tmplengths;
	map<int, vector<Point3s>> initdists, tmpdists;

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

	int countadj, countnt, countbp, initcount, tmpcount, anglecount;
	double mse_prev=0;
    double mse=0;
	double mseinit_prev=0;
	double msetmp_prev=0;
	double mseangle_prev=0;
	double mselengthadj_prev=0;
	double mselengthnt_prev=0;
	double mselengthbp_prev=0;
	double mseinit=0;
	double msetmp=0;
	double mseangle=0;
	double mselengthadj=0;
	double mselengthnt=0;
	double mselengthbp=0;

	lengths.resize(edgeinfo.size(), vector<vector<double>>(PAIRLENGTHTARGETMOVNUM, vector<double>(PAIRLENGTHTARGETMOVNUM, 0.0)));
	diffs.resize(edgeinfo.size(), vector<vector<double>>(PAIRLENGTHTARGETMOVNUM, vector<double>(PAIRLENGTHTARGETMOVNUM, 0.0)));
	dists.resize(edgeinfo.size(), vector<vector<Point3s>>(PAIRLENGTHTARGETMOVNUM, vector<Point3s>(PAIRLENGTHTARGETMOVNUM)));
	outlineridx.resize(edgeinfo.size(), vector<vector<bool>>(PAIRLENGTHTARGETMOVNUM, vector<bool>(PAIRLENGTHTARGETMOVNUM, false)));

	countadj=0; countnt=0; countbp=0;
	for(i=0;i<edgeinfo.size();i++){
		basei=edgeinfo[i].first; basej=edgeinfo[i].second;
		aat1 = restr[basei].aaa; aat2 = restr[basej].aaa;

		//cout<< basei << "-" << basej<<endl;

		if(edgetype[i] == 1){
			calcmseadjlength(restr[basei].ptv, restr[basej].ptv, basei, basej, aat1, aat2, restr[basei].resind, restr[basej].resind, paraminfo[i], lengths[i], diffs[i], dists[i], mselengthadj_prev, countadj, outlineridx[i], paircoeffadj, false);
			coeffs.push_back(paircoeffadj);
		}else if(edgetype[i] == 0){
			calcmsentlength(restr[basei].ptv, restr[basej].ptv, aat1, paraminfo[i], lengths[i], diffs[i], dists[i], mselengthnt_prev, countnt, outlineridx[i], paircoeffnt, false);
			coeffs.push_back(paircoeffnt);
		}else if(edgetype[i] == 2){
			calcmsebplength(restr[basei].ptv, restr[basej].ptv, aat1, aat2, lengths[i], diffs[i], dists[i], mselengthbp_prev, countbp, outlineridx[i], paircoeffbp, false);
			coeffs.push_back(paircoeffbp);
		}
	}

	if(countadj > 0) mselengthadj_prev=mselengthadj_prev/countadj;
	if(countnt > 0) mselengthnt_prev=mselengthnt_prev/countnt;
	if(countbp > 0) mselengthbp_prev=mselengthbp_prev/countbp;

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

	initcount=0; tmpcount=0; anglecount=0;
	for(i=0;i<movementresidues.size();i++)
	{
		basei = movementresidues[i];

		initlengths[basei].resize(PAIRLENGTHTARGETMOVNUM, 0.0);
		initdists[basei].resize(PAIRLENGTHTARGETMOVNUM);
		calcmseinitlength(restr[basei].ptv, decstr[basei].ptv, initlengths[basei], initdists[basei], mseinit_prev, initcount, initcoeff);

		tmplengths[basei].resize(PAIRLENGTHTARGETMOVNUM, 0.0);
		tmpdists[basei].resize(PAIRLENGTHTARGETMOVNUM);
		calcmsetmplength(restr[basei].ptv, snttemplates[basei].ptv, tmplengths[basei], tmpdists[basei], msetmp_prev, tmpcount, tmpcoeff);

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
	msetmp_prev=msetmp_prev/tmpcount;
	mseangle_prev=mseangle_prev/anglecount;
	mse_prev=mselengthadj_prev+mselengthnt_prev+mselengthbp_prev+mseinit_prev+msetmp_prev+mseangle_prev;

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

	//cout << "Initial:" << " mse: " << mse_prev << " mse(length adj): " << mselengthadj_prev << " mse(length nt): " << mselengthnt_prev << " mse (length bp): " << mselengthbp_prev << " mse(init): " << mseinit_prev<<endl;


	for (int epoch = 0; epoch < epochs; ++epoch) {
		initializegrad(gradptv, movenodeinfo, mse);
		mseinit=0; msetmp=0;
		mseangle=0;
		mselengthadj=0; mselengthnt=0; mselengthbp=0;

		for(i=0;i<edgeinfo.size();i++){
			basei=edgeinfo[i].first; basej=edgeinfo[i].second;
			if(edgetype[i] == 1){
				calcgradientlength(restr[basei].ptv, restr[basej].ptv, basei, basej, gradptv, movenodeinfo, lengths[i], diffs[i], dists[i], countadj, outlineridx[i], coeffs[i]);
			}else if(edgetype[i] == 0){
				calcgradientlength(restr[basei].ptv, restr[basej].ptv, basei, basej, gradptv, movenodeinfo, lengths[i], diffs[i], dists[i], countnt, outlineridx[i], coeffs[i]);
			}else if(edgetype[i] == 2){
				calcgradientlength(restr[basei].ptv, restr[basej].ptv, basei, basej, gradptv, movenodeinfo, lengths[i], diffs[i], dists[i], countbp, outlineridx[i], coeffs[i]);
			}
		}
		for(i=0;i<movementresidues.size();i++)
		{
			basei = movementresidues[i];
			calcgradientinitlength(restr[basei].ptv, decstr[basei].ptv, basei, gradptv, initlengths[basei], initdists[basei], initcount, initcoeff);
		
			calcgradienttmplength(restr[basei].ptv, snttemplates[basei].ptv, basei, gradptv, tmplengths[basei], tmpdists[basei], tmpcount, tmpcoeff);
		
			calcgradientangles(restr, basei, gradptv, movenodeinfo, anglediffs[basei], angleouterflags[basei], theta[basei], normPQs[basei], normRQs[basei], dotps[basei], vecPQs[basei], vecRQs[basei], anglecount, anglecoeff);
		}

		for(i=0;i<movenodeinfo.size();i++){
			if(movenodeinfo[i]){
				for(mtype=0;mtype<PAIRLENGTHTARGETMOVNUM;mtype++)
				{
					//cout<< restr[i].ptv[mtype].x << "-" << restr[i].ptv[mtype].y << "-" << restr[i].ptv[mtype].z << "-"<<endl;
					restr[i].ptv[mtype].x -= learning_rate*gradptv[i][mtype].x;
					restr[i].ptv[mtype].y -= learning_rate*gradptv[i][mtype].y;
					restr[i].ptv[mtype].z -= learning_rate*gradptv[i][mtype].z;
					//cout<< restr[i].ptv[mtype].x << "-" << restr[i].ptv[mtype].y << "-" << restr[i].ptv[mtype].z << "-"<<endl;
				}
			}
		}

		lengths.resize(edgeinfo.size(), vector<vector<double>>(PAIRLENGTHTARGETMOVNUM, vector<double>(PAIRLENGTHTARGETMOVNUM, 0.0)));
		diffs.resize(edgeinfo.size(), vector<vector<double>>(PAIRLENGTHTARGETMOVNUM, vector<double>(PAIRLENGTHTARGETMOVNUM, 0.0)));
		dists.resize(edgeinfo.size(), vector<vector<Point3s>>(PAIRLENGTHTARGETMOVNUM, vector<Point3s>(PAIRLENGTHTARGETMOVNUM)));
		outlineridx.resize(edgeinfo.size(), vector<vector<bool>>(PAIRLENGTHTARGETMOVNUM, vector<bool>(PAIRLENGTHTARGETMOVNUM, false)));
		countadj=0; countnt=0; countbp=0;
		for(i=0;i<edgeinfo.size();i++){
			basei=edgeinfo[i].first; basej=edgeinfo[i].second;
			aat1 = restr[basei].aaa; aat2 = restr[basej].aaa;

			if(edgetype[i] == 1){
				calcmseadjlength(restr[basei].ptv, restr[basej].ptv, basei, basej, aat1, aat2, restr[basei].resind, restr[basej].resind, paraminfo[i], lengths[i], diffs[i], dists[i], mselengthadj, countadj, outlineridx[i], paircoeffadj,false);
			}else if(edgetype[i] == 0){
				calcmsentlength(restr[basei].ptv, restr[basej].ptv, aat1, paraminfo[i], lengths[i], diffs[i], dists[i], mselengthnt, countnt, outlineridx[i], paircoeffnt,false);
			}else if(edgetype[i] == 2){
				calcmsebplength(restr[basei].ptv, restr[basej].ptv, aat1, aat2, lengths[i], diffs[i], dists[i], mselengthbp, countbp, outlineridx[i], paircoeffbp,false);
			}
		}

		if(countadj > 0) mselengthadj=mselengthadj/countadj;
		if(countnt > 0) mselengthnt=mselengthnt/countnt;
		if(countbp > 0) mselengthbp=mselengthbp/countbp;

		initcount=0; tmpcount=0; anglecount=0;
		for(i=0;i<movementresidues.size();i++)
		{
			basei = movementresidues[i];
			initlengths[basei].resize(PAIRLENGTHTARGETMOVNUM, 0.0);
			initdists[basei].resize(PAIRLENGTHTARGETMOVNUM);
			calcmseinitlength(restr[basei].ptv, decstr[basei].ptv, initlengths[basei], initdists[basei], mseinit, initcount, initcoeff);
		
			tmplengths[basei].resize(PAIRLENGTHTARGETMOVNUM, 0.0);
			tmpdists[basei].resize(PAIRLENGTHTARGETMOVNUM);
			calcmsetmplength(restr[basei].ptv, snttemplates[basei].ptv, tmplengths[basei], tmpdists[basei], msetmp, tmpcount, tmpcoeff);
		
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
		msetmp=msetmp/tmpcount;
		if(anglecount > 0) mseangle=mseangle/anglecount;
		mse=mselengthadj+mselengthnt+mselengthbp+mseinit+msetmp+mseangle;

        //cout << "epoch:" << epoch << " mse: " << mse << " mse(length adj): " << mselengthadj << " mse(length nt): " << mselengthnt << " mse (length bp): " << mselengthbp << " mse(init): " << mseinit << " prev mse: " << mse_prev << " prev mse(length adj): " << mselengthadj_prev << " prev mse(length nt): " << mselengthnt_prev << " prev mse (length bp): " << mselengthbp_prev << " prev mse(init): " << mseinit_prev<<endl;
		//cout << "epoch:" << epoch << " mse: " << mse << " mse(length adj): " << mselengthadj << " mse(length nt): " << mselengthnt << " mse (length bp): " << mselengthbp << " mse(init): " << mseinit << " mse(tmp): " << msetmp << " mse(angle): " << mseangle<<endl;


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

		mselengthadj_prev=mselengthadj;
		mselengthnt_prev=mselengthnt;
		mselengthbp_prev=mselengthbp;
		mseinit_prev = mseinit;
		msetmp_prev = msetmp;
		mseangle_prev = mseangle;
		mse_prev = mse;   
	}
	//cout << "Final:" << " mse: " << mse << " mse(length adj): " << mselengthadj << " mse(length nt): " << mselengthnt << " mse (length bp): " << mselengthbp << " mse(init): " << mseinit<<endl;


	//exit(1);



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












