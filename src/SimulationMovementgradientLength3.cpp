// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

void RNASimulation::selectbptemplate(Residueatoms restr1, Residueatoms restr2, vector<Residueatoms> &dnttemplate){
	char aat1=restr1.aaa;
	char aat2=restr2.aaa;

	dnttemplate.resize(2);

	int i, j;
	Residueatoms templatetemp;
	int atom;
	vector<vector<double>> RotMatix;
    vector<double> TranVect;
	vector<vector<double>> RotMatixtmp;
    vector<double> TranVecttmp;
	vector<double> xyztmp(3,0);
	Point3s vec;
	double minrmsd=100000000000;
	int minidx=-1;
	int count;

	vector<double> tmp(3,0);
	vector<vector<double> > xyz_list1(3,tmp);
    vector<vector<double> > xyz_list2(3,tmp);
	vector<vector<double> > xyz_listtmp(15,tmp);

	for(i=0;i<BPLENTTEMPLATE[aat1][aat2].size();i++){
		count=0;
		for(j=AccAtomParamNum[restr2.aaa][1]; j<AccAtomParamNum[restr2.aaa][2];j++){
			atom=SNTTEMPLATEATOMS[aat2][j];
			xyz_list1[count][0]=restr2.ptv[atom].x;
			xyz_list1[count][1]=restr2.ptv[atom].y;
			xyz_list1[count][2]=restr2.ptv[atom].z;

			xyz_list2[count][0]=BPLENTTEMPLATE[aat1][aat2][i][1].ptv[atom].x;
			xyz_list2[count][1]=BPLENTTEMPLATE[aat1][aat2][i][1].ptv[atom].y;
			xyz_list2[count][2]=BPLENTTEMPLATE[aat1][aat2][i][1].ptv[atom].z;
			count+=1;
		}

		bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

		count=0;
		for(j=AccAtomParamNum[restr1.aaa][1]; j<AccAtomParamNum[restr1.aaa][2];j++){
			atom=SNTTEMPLATEATOMS[aat1][j];
			xyz_list2[count][0]=BPLENTTEMPLATE[aat1][aat2][i][0].ptv[atom].x;
			xyz_list2[count][1]=BPLENTTEMPLATE[aat1][aat2][i][0].ptv[atom].y;
			xyz_list2[count][2]=BPLENTTEMPLATE[aat1][aat2][i][0].ptv[atom].z;

			bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);
			templatetemp.ptv[atom].x=xyztmp[0];
			templatetemp.ptv[atom].y=xyztmp[1];
			templatetemp.ptv[atom].z=xyztmp[2];

			count+=1;
		}

		double rmsd=0;
		for(j=AccAtomParamNum[restr1.aaa][2]-3; j<AccAtomParamNum[restr1.aaa][2];j++){
			atom=SNTTEMPLATEATOMS[aat1][j];
			vec.x=restr1.ptv[atom].x-templatetemp.ptv[atom].x;
			vec.y=restr1.ptv[atom].y-templatetemp.ptv[atom].y;
			vec.z=restr1.ptv[atom].z-templatetemp.ptv[atom].z;
			rmsd+=bf.norm(vec);
		}

		if(minrmsd > rmsd){
			RotMatixtmp=RotMatix;
			TranVecttmp=TranVect;
			minidx=i;
			minrmsd=rmsd;
		}
	}

	for(int j=0; j<SNTTEMPLATEATOMS[aat1].size();j++){
		atom=SNTTEMPLATEATOMS[aat1][j];
		xyz_listtmp[j][0]=BPLENTTEMPLATE[aat1][aat2][minidx][0].ptv[atom].x;
		xyz_listtmp[j][1]=BPLENTTEMPLATE[aat1][aat2][minidx][0].ptv[atom].y;
		xyz_listtmp[j][2]=BPLENTTEMPLATE[aat1][aat2][minidx][0].ptv[atom].z;

		bf.Transform(xyz_listtmp[j],RotMatixtmp,TranVecttmp,xyztmp);
		dnttemplate[0].ptv[atom].x=xyztmp[0];
		dnttemplate[0].ptv[atom].y=xyztmp[1];
		dnttemplate[0].ptv[atom].z=xyztmp[2];
	}

	for(int j=0; j<SNTTEMPLATEATOMS[aat2].size();j++){
		atom=SNTTEMPLATEATOMS[aat2][j];
		xyz_listtmp[j][0]=BPLENTTEMPLATE[aat1][aat2][minidx][1].ptv[atom].x;
		xyz_listtmp[j][1]=BPLENTTEMPLATE[aat1][aat2][minidx][1].ptv[atom].y;
		xyz_listtmp[j][2]=BPLENTTEMPLATE[aat1][aat2][minidx][1].ptv[atom].z;

		bf.Transform(xyz_listtmp[j],RotMatixtmp,TranVecttmp,xyztmp);
		dnttemplate[1].ptv[atom].x=xyztmp[0];
		dnttemplate[1].ptv[atom].y=xyztmp[1];
		dnttemplate[1].ptv[atom].z=xyztmp[2];
	}
}

bool RNASimulation::checkoutlinerlgm(double val, double lower, double upper){

	if(lower <= val && val <= upper){
		return false;
	}else if(lower > val){
		return true;
	}else if(val > upper){
		return true;
	}

	/*
	if(lower <= val && val <= upper){
		bdiff=0;
		diff=0;
	}else if(lower > val){
		bdiff=1;
		tmp=abs(mean-val);
		if(tmp < 1) diff=k*tmp;
		else diff=k*sqrt(tmp*tmp);
	}else if(val > upper){
		bdiff=1;
		tmp=abs(val-mean);
		if(tmp < 1) diff=k*tmp;
		else diff=k*sqrt(tmp*tmp);
	}
	*/
	return false;
}

void RNASimulation::initgeo(vector<Residueatoms> &restr, int st, int ed, int type){
	int i,j;
	char aat;
	bool flagtor;
	double randoval;

	if(st == 0){
		aat=restr[0].aaa;
		for(j=0;j<INNERANGLETARGET.size();j++){
			if(lp1.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0] > 0){
				//randoval=Random();
				if(j == 10){
					if(checkoutlinerlgm(restr[0].ang[j], OPINNERANGPARAM[aat][0].lower, OPINNERANGPARAM[aat][0].upper)){
					//if(randoval < 0.1){
						//cout<<"aaa"<<endl;
						//restr[0].ang[j]=getrandomnewparam(OPINNERANGPARAM[aat][0].lower, OPINNERANGPARAM[aat][0].upper);
						restr[0].ang[j]=getmovangval(restr[0].ang[j], OPINNERANGPARAM[aat][0].mean, OPINNERANGPARAM[aat][0].std, 0.5);
					//}
					}
				}else if(j == 11){
					if(checkoutlinerlgm(restr[0].ang[j], OPINNERANGPARAM[aat][1].lower, OPINNERANGPARAM[aat][1].upper)){
					//if(randoval < 0.1){
						//cout<<"aaa"<<endl;
						//restr[0].ang[j]=getrandomnewparam(OPINNERANGPARAM[aat][1].lower, OPINNERANGPARAM[aat][1].upper);
						restr[0].ang[j]=getmovangval(restr[0].ang[j], OPINNERANGPARAM[aat][1].mean, OPINNERANGPARAM[aat][1].std, 0.5);
					//}
					}
				}else{
					if(checkoutlinerlgm(restr[0].ang[j], INNERANGPARAM[aat][j].lower, INNERANGPARAM[aat][j].upper)){
					//if(randoval < 0.1){
						//cout<<"aaa"<<endl;
						//restr[0].ang[j]=getrandomnewparam(INNERANGPARAM[aat][j].lower, INNERANGPARAM[aat][j].upper);
						restr[0].ang[j]=getmovangval(restr[0].ang[j], INNERANGPARAM[aat][j].mean, INNERANGPARAM[aat][j].std, 0.5);
					//}
					}
				}
			}
		}

		for(j=0;j<BLENGTHTARGET.size();j++){
			if(lp1.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0] > 0){
				//randoval=Random();
				if(j == 10){
					if(checkoutlinerlgm(restr[0].len[j], OPLENGTHPARAM[aat][0].lower, OPLENGTHPARAM[aat][0].upper)){
					//if(randoval < 0.1){
						//cout<<"aaa"<<endl;
						//restr[0].len[j]=getrandomnewparam(OPLENGTHPARAM[aat][0].lower, OPLENGTHPARAM[aat][0].upper);
						restr[0].len[j]=getmovlenval(restr[0].len[j], OPLENGTHPARAM[aat][0].mean, OPLENGTHPARAM[aat][0].std, 1);
					//}
					}
				}else if(j == 11){
					if(checkoutlinerlgm(restr[0].len[j], OPLENGTHPARAM[aat][1].lower, OPLENGTHPARAM[aat][1].upper)){
					//if(randoval < 0.1){
						//cout<<"aaa"<<endl;
						//restr[0].len[j]=getrandomnewparam(OPLENGTHPARAM[aat][1].lower, OPLENGTHPARAM[aat][1].upper);
						restr[0].len[j]=getmovlenval(restr[0].len[j], OPLENGTHPARAM[aat][1].mean, OPLENGTHPARAM[aat][1].std, 1);
					//}
					}
				}else{
					if(checkoutlinerlgm(restr[0].len[j], LENGTHPARAM[aat][j].lower, LENGTHPARAM[aat][j].upper)){
					//if(randoval < 0.1){
						//cout<<"aaa"<<endl;
						//restr[0].len[j]=getrandomnewparam(LENGTHPARAM[aat][j].lower, LENGTHPARAM[aat][j].upper);
						restr[0].len[j]=getmovlenval(restr[0].len[j], LENGTHPARAM[aat][j].mean, LENGTHPARAM[aat][j].std, 1);
					//}
					}
				}
			}
		}

		//flagtor=param2coordbb(restr,0,1);
		//if(type > 0){
			//param2coordsugar(restr,0,1);
			//param2coordbase(restr,0,1);
			//param2coordhydrogen(restr,0,1);
		//}
	}

	for(int i=st;i<ed;i++){
		aat=restr[i].aaa;
		for(j=0;j<INNERANGLETARGET.size();j++){
			randoval=Random();
			if(lp1.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][i] > 0){
				if(checkoutlinerlgm(restr[i].ang[j], INNERANGPARAM[aat][j].lower, INNERANGPARAM[aat][j].upper)){
				//if(randoval < 0.1){
					//cout<<"aaa"<<endl;
					//restr[i].ang[j]=getrandomnewparam(INNERANGPARAM[aat][j].lower, INNERANGPARAM[aat][j].upper);
					restr[i].ang[j]=getmovangval(restr[i].ang[j], INNERANGPARAM[aat][j].mean, INNERANGPARAM[aat][j].std, 0.5);
				//}
				}
			}
		}

		for(j=0;j<BLENGTHTARGET.size();j++){
			randoval=Random();
			if(lp1.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][i] > 0){
				if(checkoutlinerlgm(restr[i].len[j], LENGTHPARAM[aat][j].lower, LENGTHPARAM[aat][j].upper)){
				//if(randoval < 0.1){
					//cout<<"aaa"<<endl;
					//restr[i].len[j]=getrandomnewparam(LENGTHPARAM[aat][j].lower, LENGTHPARAM[aat][j].upper);
					restr[i].len[j]=getmovlenval(restr[i].len[j], LENGTHPARAM[aat][j].mean, LENGTHPARAM[aat][j].std, 1);
				//}
				}
			}
		}
	}

	flagtor=param2coordbb(restr,st,ed);
	if(type > 0){
		param2coordsugar(restr,st,ed);
		param2coordbase(restr,st,ed);
		param2coordhydrogen(restr,st,ed);
	}
}

double RNASimulation::distfromene(vector<Residueatoms> restr, vector<Residueatoms> decstr, int st, int ed, int type){
	Point3s vec;
	int i,j;
	double rnorm;
	int mtype1, mtype2;
	double rmsd1=0;
	double rmsd2=0;
	int resnumst, resnum;
	int basei, basej;
	double diff;
	
	if(type==0){
		resnumst=0;
		resnum=6;
	}

	
	for(i=st;i<ed;i++){
		if(type==1){
			resnumst=AccAtomParamNum[restr[i].aaa][1];
			resnum=AccAtomParamNum[restr[i].aaa][3];
		}

		for(j=resnumst;j<resnum;j++){
			vec.x = restr[i].ptv[j].x - decstr[i].ptv[j].x; 
			vec.y = restr[i].ptv[j].y - decstr[i].ptv[j].y; 
			vec.z = restr[i].ptv[j].z - decstr[i].ptv[j].z;

			rmsd1+=bf.norm(vec);
		}
	}


	/*
	for(i=0;i<bpairs.size();i++){
		basei = bpairs[i].first; basej = bpairs[i].second;
		if(basei >= 0 & basej >= 0){
			for(mtype1=PAIRLENGTHTARGETMOVNUM-3;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
			{
					
					vec.x = restr[basei].ptv[mtype1].x - bptemplates[basei][0].ptv[mtype1].x;
					vec.y = restr[basei].ptv[mtype1].y - bptemplates[basei][0].ptv[mtype1].y;
					vec.z = restr[basei].ptv[mtype1].z - bptemplates[basei][0].ptv[mtype1].z;

					rnorm=bf.norm(vec);

					//cout<<basei << "-" << basej << "-" << mtype1 << "-" << mtype2 << "-" << rnorm << "-" << bptemplateslength[basei].interdistmap[mtype1][mtype2] << "-" << diff<<endl;
					rmsd2 += sqrt(rnorm*rnorm);
			}
		}
	}
	*/

	return rmsd1;
}

double RNASimulation::goldenSearchcalcene(vector<Residueatoms> restr, vector<Residueatoms> decstr, double torval, int dihtype, int trandpos, int trandposst, int trandposed){
	vector<Residueatoms> restrtmp;
	bool flagtor;
	restrtmp.resize(restr.size());
	restrtmp=restr;

	restrtmp[trandpos].dih[dihtype]=torval;
	flagtor=param2coordbb(restrtmp,trandpos,trandposed+1);
	param2coordsugar(restrtmp,trandpos,trandposed+1);
	param2coordbase(restrtmp,trandpos,trandposed+1);

	double ene=distfromene(restrtmp, decstr, trandposst, trandposed+1, 1);

	return ene;
}

const double PHI = (1 + std::sqrt(5)) / 2;
double RNASimulation::goldenSearchlgm(vector<Residueatoms> restr, vector<Residueatoms> decstr, int trandpos, int dihtype, int trandposst, int trandposed, double a, double b, double tolerance){
	double c = b - (b - a) / PHI;
    double d = a + (b - a) / PHI;
	char aat=restr[trandpos].aaa;
	double cene, dene;

	while (std::abs(b - a) > tolerance) {
		cene=goldenSearchcalcene(restr, decstr, c, dihtype, trandpos, trandposst, trandposed);
		dene=goldenSearchcalcene(restr, decstr, d, dihtype, trandpos, trandposst, trandposed);

		if (cene < dene) {
            b = d;
        } else {
            a = c;
        }

		c = b - (b - a) / PHI;
        d = a + (b - a) / PHI;
	}

	return (a + b) / 2;
}


bool RNASimulation::fmovementlgm(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc)
{
	int i,j,k,inext;
	int trandpos, trandposed, trandposst;
	int wpos;
	char aat;
	double randval=Random();
	bool exitOuterLoop = false;
	map<int, vector<Point3s>> snttemplatestr;
	int dihtype;
	bool flagtor;
	int dihcls;

	randval=Random();

	int totwrong3=0;
	for(i=0;i<BLENGTHTARGET.size();i++){
		totwrong3+=lp1.nn[BLENGTH_ENERGY_ST+BLENGTHTARGET[i]];
		//totwrong3-=lp1.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[i]][0];
	}
	for(i=0;i<INNERANGLETARGET.size();i++){
		totwrong3+=lp1.nn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[i]];
		//totwrong3-=lp1.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[i]][0];
	}

	if(totwrong3>0){
		wpos=int(Random()*totwrong3-EPS);

		i=0;
		for(trandpos=0;trandpos<numseq;trandpos++)
		{
			for(j=0;j<BLENGTHTARGET.size();j++){
				inext=i+lp1.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][trandpos];
				if(i<=wpos & wpos<inext){
					exitOuterLoop = true;
					break;
				}
				i=inext;
			}

			for(j=0;j<INNERANGLETARGET.size();j++){
				inext=i+lp1.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][trandpos];
				if(i<=wpos & wpos<inext){
					exitOuterLoop = true;
					break;
				}
				i=inext;
			}
			if(exitOuterLoop) break;
		}
	}else{
		newenergy=oldenergy;
		return false;
	}

	i=0;
	for(trandposed=trandpos;trandposed<numseq;trandposed++)
	{
		exitOuterLoop = false;
		for(j=0;j<BLENGTHTARGET.size();j++){
			if(lp1.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][trandposed] == 1){
				exitOuterLoop = true;
				break;
			}
		}

		for(j=0;j<INNERANGLETARGET.size();j++){
			if(lp1.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][trandposed] == 1){
				exitOuterLoop = true;
				break;
			}
		}

		if(!exitOuterLoop){
			break;
		}
	}

	if(trandposed >= numseq) trandposed=numseq-1;

	i=0;
	for(trandposst=trandpos;trandposst>=0;trandposst--)
	{
		exitOuterLoop = false;
		for(j=0;j<BLENGTHTARGET.size();j++){
			if(lp1.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][trandposst] == 1){
				exitOuterLoop = true;
				break;
			}
		}

		for(j=0;j<INNERANGLETARGET.size();j++){
			if(lp1.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][trandposst] == 1){
				exitOuterLoop = true;
				break;
			}
		}

		if(!exitOuterLoop){
			break;
		}
	}

	if(trandposst < 0) trandposst=0;

	int trandposedtmp=trandposed;
	int trandpossttmp=trandposst;

	int residuenum;
	double maxspan=2+int(6*Random());
	if(trandposed-trandposst+1 > maxspan){
		residuenum=trandposed-trandposst-maxspan+2;
		wpos=int(Random()*residuenum-EPS);
		trandposst=wpos+trandposst;
		trandposed=trandposst+maxspan-1;
	}

	//cout<<trandposst << "-" << trandpos << "-" << trandposed<<endl;

	if(trandposedtmp<trandposed | trandpossttmp>trandposst){
		cout<<"Error"<<endl;
		exit(1);
	}

	copystr(decstr, numseq);

	/*
	vector<pair<int,int>> bpairs; 
	map<int, vector<Residueatoms>> bptemplates;
	map<int, TemplateDistMap> bptemplateslength;
	int basei, basej;
	char aat1, aat2;
	for(i=trandposst;i<trandposed+1;i++){
		if(bppairednt[i] >= 0){

			basei=i; 
			basej=bppairednt[i];
			aat1 = restr[basei].aaa; 
			aat2 = restr[basej].aaa;

			bpairs.push_back({i,bppairednt[i]});
			selectbptemplate(restr[basei], restr[basej], bptemplates[i]);
		}else{
			bpairs.push_back({i,-1});
		}
	}
	*/

	//initgeo(restr, trandposst, trandposed+1, 0);
	residuenum=trandposed-trandposst+1;

	//for(i=trandposst;i<trandposed+1;i++){
		//for(int j=0;j<AccAtomParamNum[decstr[i].aaa][3];j++){
			//cout<<decstr[i].ptv[j].x << "-" << decstr[i].ptv[j].y << "-" << decstr[i].ptv[j].z<<endl;
		//}
	//}

	//cout<<"==========="<<endl;

	//for(i=trandposst;i<trandposed+1;i++){
		//for(int j=0;j<AccAtomParamNum[restr[i].aaa][3];j++){
			//cout<<restr[i].ptv[j].x << "-" << restr[i].ptv[j].y << "-" << restr[i].ptv[j].z<<endl;
		//}
	//}


	vector<Residueatoms> restrtmp;
	restrtmp.resize(restr.size());
	restrtmp=restr;
	double enetmp,ene;
	//double ene=distfromene(restr, decstr, trandposst, trandposed+1, 0);

	
	//backbone refinement
	int maxiter=100000+1000*residuenum;
	int iter=0;

	/*
	for(iter=0;iter<maxiter;iter++){
		dihtype=TORSIONANGLETARGET[int(Random()*6)];
		wpos=int(Random()*residuenum-EPS)+trandposst;
		aat=restrtmp[wpos].aaa;
		restrtmp=restr;

		restrtmp[wpos].dih[dihtype]=getrandomnewparam(TORSIONPARAM[aat][dihtype].lower, TORSIONPARAM[aat][dihtype].upper);
		flagtor=param2coordbb(restrtmp,wpos,trandposed+1);

		enetmp=distfromene(restrtmp, decstr, bpairs, bptemplates, bptemplateslength, trandposst, trandposed+1, 0);

		if(ene>enetmp){
			ene=enetmp;
			restr=restrtmp;
			cout<<ene<<endl;
		}
	}
	*/

	
	initgeo(restr, trandposst, trandposed+1, 1);
	
	//maxiter=1500+150*residuenum;
	//maxiter=15000+1500*residuenum;
	//maxiter=1000+100*residuenum;
	maxiter=100+10*residuenum;
	//maxiter=500+50*residuenum;
	iter=0;
	ene=distfromene(restr, decstr, trandposst, trandposed+1, 1);
	for(iter=0;iter<maxiter;iter++){
		dihtype=TORSIONANGLETARGET[int(Random()*TORSIONANGLETARGET.size())];
		wpos=int(Random()*residuenum-EPS)+trandposst;
		aat=restrtmp[wpos].aaa;
		restrtmp=restr;

		//restrtmp[wpos].dih[dihtype]=getrandomnewparam(TORSIONPARAM[aat][dihtype].lower, TORSIONPARAM[aat][dihtype].upper);
		
		/*
		if(iter%1000 == 0 & iter != 0){
			initgeo(restrtmp, trandposst, trandposed+1, 1);

			flagtor=param2coordbb(restrtmp,wpos,trandposed+1);
			param2coordsugar(restrtmp,wpos,trandposed+1);
			param2coordbase(restrtmp,wpos,trandposed+1);

			enetmp=distfromene(restrtmp, decstr, trandposst, trandposed+1, 1);
			ene=enetmp;
			restr=restrtmp;
		}
		*/
		
		restrtmp[wpos].dih[dihtype]=goldenSearchlgm(restrtmp, decstr, wpos, dihtype, trandposst, trandposed, TORSIONPARAM[aat][dihtype].lower, TORSIONPARAM[aat][dihtype].upper, 1);
		//restrtmp[wpos].dih[dihtype]=getrandomnewparam(TORSIONPARAM[aat][dihtype].lower, TORSIONPARAM[aat][dihtype].upper);
		//cout<<iter<<endl;

		//
		/*
		randval=Random();
		if(randval < 0.334){
			dihcls=int(TORSIONLIST[1].size()*Random());
			for(i=0;i<TorsionClsList[1].size();i++)
			{
				restr[wpos].dih[TorsionClsList[1][i]]=getrandomnewtor(TORSIONLIST[1][dihcls][i].lower, TORSIONLIST[1][dihcls][i].upper, TORSIONLIST[1][dihcls][i].rangetype);
			}
		}else if(0.334<=randval & randval <0.667){
			dihcls=int(TORSIONLIST[0].size()*Random());
			for(i=1;i<TorsionClsList[0].size();i++){
				restr[wpos].dih[TorsionClsList[0][i]]=getrandomnewtor(TORSIONLIST[0][dihcls][i].lower, TORSIONLIST[0][dihcls][i].upper, TORSIONLIST[0][dihcls][i].rangetype);
			}
		}else{
			restrtmp[wpos].dih[dihtype]=getrandomnewparam(TORSIONPARAM[aat][dihtype].lower, TORSIONPARAM[aat][dihtype].upper);
		}
		*/

		flagtor=param2coordbb(restrtmp,wpos,trandposed+1);
		param2coordsugar(restrtmp,wpos,trandposed+1);
		param2coordbase(restrtmp,wpos,trandposed+1);

		enetmp=distfromene(restrtmp, decstr, trandposst, trandposed+1, 1);

		if(ene>enetmp){
			ene=enetmp;
			restr=restrtmp;
			//cout<<ene<<endl;
		}
	}
	

	flagtor=param2coordbb(restr,trandposst,numseq);
	param2coordsugar(restr,trandposst,numseq);
	param2coordbase(restr,trandposst,numseq);
	param2coordhydrogen(restr,trandposst, numseq);

	if(trandposst == 0){
		coord2parambb(restr,trandposst, trandposst+1,1);
		coord2paramsugar(restr,trandposst, trandposst+1,1);
		param2coordbase(restr,trandposst, trandposst+1);
		param2coordhydrogen(restr,trandposst, trandposst+1);
	}
	coord2parambb(restr,trandposed, trandposed+1,1);
	coord2paramsugar(restr,trandposed, trandposed+1,1);
	param2coordbase(restr,trandposed, trandposed+1);
	
	coord2paramadd(restr,trandposst, trandposed+1);
	settmplmpflap(numseq, trandposst, trandposed+1);

	//cout<<"==========="<<endl;

	//for(i=trandposst;i<trandposed+1;i++){
		//for(int j=0;j<AccAtomParamNum[restr[i].aaa][3];j++){
			//cout<<restr[i].ptv[j].x << "-" << restr[i].ptv[j].y << "-" << restr[i].ptv[j].z<<endl;
		//}
	//}

	//exit(1);

	newenergy=fcalcrmsdenergy(restr, numseq, lp2, rppmatbk, bpenergybk, bsenergybk, nmapbk, distmapbk, ntpairdistbk, adjpairdistbk, distenebk, multidistenebk, contactenebk, mapflagbk, mapbbflagbk, flagacc);
	//exit(1);
	
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
	//bool upidx=metropolisUpdate(oldenergy, newenergy, 30*Random());

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












