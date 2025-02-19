// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

double RNASimulation::getbptorsionene(double val, vector<double> boundary, vector<double> energy){
	int binidx=getbin(val, boundary);
	return energy[binidx];
}

void RNASimulation::getoutenebaselength(double val, double mean, double lower, double upper, double k, double &diff, int &bdiff){
	double tmp;

	if(lower <= val && val <= upper){
		bdiff=0;
		diff=0;
	}else{
		bdiff=1;
		tmp=val-mean;
		diff=k*sqrt(tmp*tmp);
	}
}

void RNASimulation::getplanebp(Point3s pointA, Point3s pointB, Point3s pointC, double &a, double &b, double &c, double &d){
	Point3s AB;
    Point3s AC;

	AB.x=pointB.x - pointA.x; AB.y=pointB.y - pointA.y; AB.z=pointB.z - pointA.z;
	AC.x=pointC.x - pointA.x; AC.y=pointC.y - pointA.y; AC.z=pointC.z - pointA.z;

	Point3s normal = bf.basecross(AB,AC);

	a = normal.x;
    b = normal.y;
    c = normal.z;

	d = -(a * pointA.x + b * pointA.y + c * pointA.z);
}

double RNASimulation::getdistfromplane(double a, double b, double c, double d, Point3s pointP){
	double numerator = abs(a * pointP.x + b * pointP.y + c * pointP.z + d);
    double denominator = sqrt(a * a + b * b + c * c);
    
    return numerator / denominator;
}

double RNASimulation::bpplaneene(vector<Residueatoms> decstr, int basei, int basej){
	double a, b, c, d;
	Point3s pointA, pointB, pointC;
	char aat1, aat2;
	double ene=0.0;
	aat1 = decstr[basei].aaa;
	aat2 = decstr[basej].aaa;

	vector<int> baseatoms1=BaseHeavyAtom[aat1];
	vector<int> baseatoms2=BaseHeavyAtom[aat2];

	pointA=decstr[basei].ptv[baseatoms1[0]]; pointB=decstr[basei].ptv[baseatoms1[1]]; pointC=decstr[basei].ptv[baseatoms1[2]];
	getplanebp(pointA, pointB, pointC, a, b, c, d);

	for(int i=0; i<BaseHeavyAtom[aat2].size();i++){
		ene+=getdistfromplane(a, b, c, d, decstr[basej].ptv[baseatoms2[i]]);
	}

	pointA=decstr[basej].ptv[baseatoms2[0]]; pointB=decstr[basej].ptv[baseatoms2[1]]; pointC=decstr[basej].ptv[baseatoms2[2]];
	getplanebp(pointA, pointB, pointC, a, b, c, d);
	
	for(int i=0; i<BaseHeavyAtom[aat1].size();i++){
		ene+=getdistfromplane(a, b, c, d, decstr[basei].ptv[baseatoms1[i]]);
	}

	return ene;
}


void RNASimulation::fenergybasepair(vector<Residueatoms> decstr,int numseq,double &nene,double &fene,Energyterms &lp, vector<Basecontact> &baserpp, vector<vector<bool>> mapflagrpp,bool flagacc){
	int i,j,k;
	int basei,basej;
	char aat1, aat2;
	int atom1, atom2, atom3, atom4;
	int atom11,atom12,atom13,atom21,atom22,atom23;
	Point3s tmp1,tmp2,tmp3,tmp4;
	Point3s res1N91, res1C42, res1C86, res2N91, res2C42, res2C86;
	double r, alpha, beta, gamma;
	int binidx1, binidx2, binidx3;
	int binnum1, binnum2, binnum3;
	Point3s pt;
	double dih, len;
	int pairwiseidx;
	vector<double> indf(BASE_ENERGY_NUM,0);
	vector<int> indn(BASE_ENERGY_NUM,0);
	double indftmp;
	int indntmp;

	//binnum1=BASEPAIRROT[0].binnum1;
	binnum2=BASEPAIRROT[0].binnum2;
	binnum3=BASEPAIRROT[0].binnum3;

	if(!flagacc)
	{
		for(j=0;j<BASE_ENERGY_NUM;j++)
		{
			lp.fn[BASE_ENERGY_ST+j]=0;
			lp.nn[BASE_ENERGY_ST+j]=0;
		}
		
		for(i=0;i<BASEPAIRS.size();i++)
		{
			basei=BASEPAIRS[i].basepidx.first;//base idx 1
			basej=BASEPAIRS[i].basepidx.second;//base idx 2

			for(j=0;j<BASE_ENERGY_NUM;j++)
			{
				lp.indf[BASE_ENERGY_ST+j][basei]=0;
				lp.indf[BASE_ENERGY_ST+j][basej]=0;
				lp.indn[BASE_ENERGY_ST+j][basei]=0;
				lp.indn[BASE_ENERGY_ST+j][basej]=0;
			}

			pairwiseidx=bpmapid2idx[basei*numseq+basej];

			if(basei < basej)
			{
				aat1=decstr[basei].aaa;
				aat2=decstr[basej].aaa;

				for(j=0;j<BASE_ENERGY_NUM;j++)
				{
					indf[j]=0; indn[j]=0;
				}

				//atom11=BASEAXISBP[aat1][aat2][0][0];
				//atom12=BASEAXISBP[aat1][aat2][0][1];
				//atom13=BASEAXISBP[aat1][aat2][0][2];
				//atom21=BASEAXISBP[aat1][aat2][0][3];
				//atom22=BASEAXISBP[aat1][aat2][0][4];
				//atom23=BASEAXISBP[aat1][aat2][0][5];

				//res1N91=decstr[basei].ptv[atom11]; res1C42=decstr[basei].ptv[atom12]; res1C86=decstr[basei].ptv[atom13];
				//res2N91=decstr[basej].ptv[atom21]; res2C42=decstr[basej].ptv[atom22]; res2C86=decstr[basej].ptv[atom23];

				//getbaseangles(res1N91, res1C42, res1C86, res2N91, res2C42, res2C86, r, alpha, beta, gamma);

				//getoutenebaselength(r, BASEPAIRPARAM[0][7].mean, BASEPAIRPARAM[0][7].lower, BASEPAIRPARAM[0][7].upper, 1, indftmp, indntmp);
				//indf[0]+=indftmp; indn[0]+=indntmp;

				//if(indn[0] == 1){
					//cout<<"a1-" + to_string(r) + "-" + to_string(BASEPAIRPARAM[0][7].lower) + "-" + to_string(BASEPAIRPARAM[0][7].upper) + "-" + to_string(indf[0]) + "-" + to_string(indn[0])<<endl;
				//}else{
					//cout<<"b1-" + to_string(r) + "-" + to_string(BASEPAIRPARAM[0][7].lower) + "-" + to_string(BASEPAIRPARAM[0][7].upper) + "-" + to_string(indf[0]) + "-" + to_string(indn[0])<<endl;
				//}

				//binidx1=getbin(alpha, BASEPAIRROT[0].boundary1);
				//binidx2=getbin(beta, BASEPAIRROT[0].boundary2);
				//binidx3=getbin(gamma, BASEPAIRROT[0].boundary3);

				//cout<<to_string(binidx1) + "-" + to_string(binidx2) + "-" + to_string(binidx3)<<endl;
				//cout<<to_string(binnum2) + "-" + to_string(binnum3)<<endl;

				//indftmp=BASEPAIRROT[0].energy[binidx1*binnum2*binnum3+binidx2*binnum3+binidx3];
				//indntmp=BASEPAIRROT[0].binary[binidx1*binnum2*binnum3+binidx2*binnum3+binidx3];
				//indf[1]+=indftmp; indn[1]+=indntmp;

				//if(indn[1] == 1){
					//cout<<"a2-" + to_string(r) + "-" + to_string(alpha) + "-" + to_string(beta) + "-" + to_string(gamma) + "-" + to_string(indf[1]) + "-" + to_string(indn[1])<<endl;
				//}else{
					//cout<<"b2-" + to_string(r) + "-" + to_string(alpha) + "-" + to_string(beta) + "-" + to_string(gamma) + "-" + to_string(indf[1]) + "-" + to_string(indn[1])<<endl;
				//}


				//torsion C1'-N9-N9-C1'
				atom1=BasePairTorList[aat1][aat2][0][0];
				atom2=BasePairTorList[aat1][aat2][0][1];
				atom3=BasePairTorList[aat1][aat2][0][2];
				atom4=BasePairTorList[aat1][aat2][0][3];
				tmp1=decstr[basei].ptv[atom1];tmp2=decstr[basei].ptv[atom2];tmp3=decstr[basej].ptv[atom3];tmp4=decstr[basej].ptv[atom4];
				dih=bf.phi(tmp1.x,tmp1.y,tmp1.z,tmp2.x,tmp2.y,tmp2.z,tmp3.x,tmp3.y,tmp3.z,tmp4.x,tmp4.y,tmp4.z);
				if(dih>360) dih=dih-360;
				else if(dih<0) dih=dih+360;
				getoutenetorangle(dih, BASEPAIRPARAM[1][0].mean, BASEPAIRPARAM[1][0].lower, BASEPAIRPARAM[1][0].upper, 1, indftmp, indntmp);
				indf[3]+=indftmp; indn[3]+=indntmp;
				//indftmp=getbptorsionene(dih, BASEPAIRPARAM[1][2].boundary, BASEPAIRPARAM[1][2].energy);
				//indf[0]+=indftmp;

				//cout<<to_string(basei) + "-" + to_string(basej) + "-" + to_string(dih) + "-" + to_string(indftmp)<<endl;

				//torsion O3'-C3'-C3'-O3'
				//atom1=BasePairTorList[aat1][aat2][1][0];
				//atom2=BasePairTorList[aat1][aat2][1][1];
				//atom3=BasePairTorList[aat1][aat2][1][2];
				//atom4=BasePairTorList[aat1][aat2][1][3];
				//tmp1=decstr[basei].ptv[atom1];tmp2=decstr[basei].ptv[atom2];tmp3=decstr[basej].ptv[atom3];tmp4=decstr[basej].ptv[atom4];
				//dih=bf.phi(tmp1.x,tmp1.y,tmp1.z,tmp2.x,tmp2.y,tmp2.z,tmp3.x,tmp3.y,tmp3.z,tmp4.x,tmp4.y,tmp4.z);
				//if(dih>360) dih=dih-360;
				//else if(dih<0) dih=dih+360;
				//indftmp=getbptorsionene(dih, BASEPAIRPARAM[1][3].boundary, BASEPAIRPARAM[1][3].energy);
				//indf[0]+=indftmp;

				//cout<<to_string(basei) + "-" + to_string(basej) + "-" + to_string(dih) + "-" + to_string(indftmp)<<endl;

				//C1'-C1'
				//atom1=BasePairLenList[aat1][aat2][6][0];
				//atom2=BasePairLenList[aat1][aat2][6][1];
				//tmp1=decstr[basei].ptv[atom1];
				//tmp2=decstr[basej].ptv[atom2];
				//pt.x=double(tmp1.x-tmp2.x); pt.y=double(tmp1.y-tmp2.y); pt.z=double(tmp1.z-tmp2.z);
				//len=bf.norm(pt);
				//getoutenebaselength(len, BASEPAIRPARAM[0][6].mean, BASEPAIRPARAM[0][6].lower, BASEPAIRPARAM[0][6].upper, 1, indftmp, indntmp);
				//indf[2]+=indftmp; indn[2]+=indntmp;
				
				if((aat1 == 'A' & aat2 == 'U') | (aat1 == 'U' & aat2 == 'A') | (aat1 == 'G' & aat2 == 'C') | (aat1 == 'C' & aat2 == 'G')){
					indf[2]+=bpplaneene(decstr, basei, basej);
				}
				

				//cout<<to_string(BASEPAIRPARAM[0][3].lower) + "-" + to_string(BASEPAIRPARAM[0][3].upper)<<endl;
				//cout<<to_string(basei) + "-" + to_string(basej) + "-" + to_string(len) + "-" + to_string(indftmp)<<endl;

				//C3'-C3'
				//atom1=BasePairLenList[aat1][aat2][1][0];
				//atom2=BasePairLenList[aat1][aat2][1][1];
				//tmp1=decstr[basei].ptv[atom1];
				//tmp2=decstr[basej].ptv[atom2];
				//pt.x=double(tmp1.x-tmp2.x); pt.y=double(tmp1.y-tmp2.y); pt.z=double(tmp1.z-tmp2.z);
				//len=bf.norm(pt);
				//getoutenebaselength(len, BASEPAIRPARAM[0][4].mean, BASEPAIRPARAM[0][4].lower, BASEPAIRPARAM[0][4].upper , 1, indftmp, indntmp);
				//indf[0]+=indftmp; indn[0]+=indntmp;

				//cout<<to_string(BASEPAIRPARAM[0][4].lower) + "-" + to_string(BASEPAIRPARAM[0][4].upper)<<endl;
				//cout<<to_string(basei) + "-" + to_string(basej) + "-" + to_string(len) + "-" + to_string(indftmp)<<endl;

				for(j=0;j<BASE_ENERGY_NUM;j++)
				{
					//cout<<to_string(basei) + "-" + to_string(basej) + "-" + to_string(j) + "-" + to_string(indf[j]) + "-" + to_string(indn[j])<<endl;

					baserpp[pairwiseidx].pval[j]=indf[j];
					//baserpp[pairwiseidx].pval[j]=indf[j];
					baserpp[pairwiseidx].pcw[j]=indn[j];
					//baserpp[pairwiseidx].pcw[j]=indn[j];

					lp.indf[BASE_ENERGY_ST+j][basei]+=indf[j];
					lp.indf[BASE_ENERGY_ST+j][basej]+=indf[j];
					lp.indn[BASE_ENERGY_ST+j][basei]+=indn[j];
					lp.indn[BASE_ENERGY_ST+j][basej]+=indn[j];
					lp.fn[BASE_ENERGY_ST+j]+=indf[j];
					lp.nn[BASE_ENERGY_ST+j]+=indn[j];
				}
			}
		}
	}
	else
	{
		//int idx = 0;
		int iold[BASE_ENERGY_NUM],inew[BASE_ENERGY_NUM];
		double fold[BASE_ENERGY_NUM],fnew[BASE_ENERGY_NUM];
		for(i=0;i<BASE_ENERGY_NUM;i++){
			iold[i]=0;inew[i]=0;
			fold[i]=0;fnew[i]=0;
		}

		for(i=0;i<BASEPAIRS.size();i++)
		{
			basei=BASEPAIRS[i].basepidx.first;//base idx 1
			basej=BASEPAIRS[i].basepidx.second;//base idx 2

			pairwiseidx=bpmapid2idx[basei*numseq+basej];

			if(basei<basej)
			{

				//if(!mapflagrpp[basei][basej]) continue;
			
				for(j=0;j<BASE_ENERGY_NUM;j++)
				{
					fold[j]+=baserpp[pairwiseidx].pval[j];
					iold[j]+=baserpp[pairwiseidx].pcw[j];

					lp.indf[BASE_ENERGY_ST+j][basei]=lp.indf[BASE_ENERGY_ST+j][basei]-baserpp[pairwiseidx].pval[j];
					lp.indf[BASE_ENERGY_ST+j][basej]=lp.indf[BASE_ENERGY_ST+j][basej]-baserpp[pairwiseidx].pval[j];
					lp.indn[BASE_ENERGY_ST+j][basei]=lp.indn[BASE_ENERGY_ST+j][basei]-baserpp[pairwiseidx].pcw[j];
					lp.indn[BASE_ENERGY_ST+j][basej]=lp.indn[BASE_ENERGY_ST+j][basej]-baserpp[pairwiseidx].pcw[j];

					baserpp[pairwiseidx].pval[j]=0;
					baserpp[pairwiseidx].pcw[j]=0;
				}

				aat1=decstr[basei].aaa;
				aat2=decstr[basej].aaa;

				for(j=0;j<BASE_ENERGY_NUM;j++)
				{
					indf[j]=0; indn[j]=0;
				}

				//atom11=BASEAXISBP[aat1][aat2][0][0];
				//atom12=BASEAXISBP[aat1][aat2][0][1];
				//atom13=BASEAXISBP[aat1][aat2][0][2];
				//atom21=BASEAXISBP[aat1][aat2][0][3];
				//atom22=BASEAXISBP[aat1][aat2][0][4];
				//atom23=BASEAXISBP[aat1][aat2][0][5];

				//res1N91=decstr[basei].ptv[atom11]; res1C42=decstr[basei].ptv[atom12]; res1C86=decstr[basei].ptv[atom13];
				//res2N91=decstr[basej].ptv[atom21]; res2C42=decstr[basej].ptv[atom22]; res2C86=decstr[basej].ptv[atom23];

				//getbaseangles(res1N91, res1C42, res1C86, res2N91, res2C42, res2C86, r, alpha, beta, gamma);

				//getoutenebaselength(r, BASEPAIRPARAM[0][7].mean, BASEPAIRPARAM[0][7].lower, BASEPAIRPARAM[0][7].upper, 1, indftmp, indntmp);
				//indf[0]+=indftmp; indn[0]+=indntmp;

				//if(indn[0] == 1) cout<< 0 <<"-"<< indf[0]<<endl;

				//if(indn[0] == 1){
					//cout<<"a1-" + to_string(r) + "-" + to_string(BASEPAIRPARAM[0][7].lower) + "-" + to_string(BASEPAIRPARAM[0][7].upper) + "-" + to_string(indf[0]) + "-" + to_string(indn[0])<<endl;
				//}else{
					//cout<<"b1-" + to_string(r) + "-" + to_string(BASEPAIRPARAM[0][7].lower) + "-" + to_string(BASEPAIRPARAM[0][7].upper) + "-" + to_string(indf[0]) + "-" + to_string(indn[0])<<endl;
				//}

				//binidx1=getbin(alpha, BASEPAIRROT[0].boundary1);
				//binidx2=getbin(beta, BASEPAIRROT[0].boundary2);
				//binidx3=getbin(gamma, BASEPAIRROT[0].boundary3);

				//cout<<to_string(binidx1) + "-" + to_string(binidx2) + "-" + to_string(binidx3)<<endl;
				//cout<<to_string(binnum2) + "-" + to_string(binnum3)<<endl;

				//indftmp=BASEPAIRROT[0].energy[binidx1*binnum2*binnum3+binidx2*binnum3+binidx3];
				//indntmp=BASEPAIRROT[0].binary[binidx1*binnum2*binnum3+binidx2*binnum3+binidx3];
				//indf[1]+=indftmp; indn[1]+=indntmp;

				//if(indn[1] == 1) cout<< 1 <<"-"<< indf[1]<<endl;

				//torsion C1'-N9-N9-C1'
				atom1=BasePairTorList[aat1][aat2][0][0];
				atom2=BasePairTorList[aat1][aat2][0][1];
				atom3=BasePairTorList[aat1][aat2][0][2];
				atom4=BasePairTorList[aat1][aat2][0][3];
				tmp1=decstr[basei].ptv[atom1];tmp2=decstr[basei].ptv[atom2];tmp3=decstr[basej].ptv[atom3];tmp4=decstr[basej].ptv[atom4];
				dih=bf.phi(tmp1.x,tmp1.y,tmp1.z,tmp2.x,tmp2.y,tmp2.z,tmp3.x,tmp3.y,tmp3.z,tmp4.x,tmp4.y,tmp4.z);
				if(dih>360) dih=dih-360;
				else if(dih<0) dih=dih+360;
				getoutenetorangle(dih, BASEPAIRPARAM[1][0].mean, BASEPAIRPARAM[1][0].lower, BASEPAIRPARAM[1][0].upper, 1, indftmp, indntmp);
				indf[3]+=indftmp; indn[3]+=indntmp;

				//if(indn[3] == 1) cout<< 3 <<"-"<< indf[3] << "-" << dih << "-" << BASEPAIRPARAM[1][0].lower << "-" << BASEPAIRPARAM[1][0].upper <<endl;

				//indf[2]+=getbptorsionene(dih, BASEPAIRPARAM[1][0].boundary, BASEPAIRPARAM[1][0].energy);

				//torsion O3'-C3'-C3'-O3'
				//atom1=BasePairTorList[aat1][aat2][1][0];
				//atom2=BasePairTorList[aat1][aat2][1][1];
				//atom3=BasePairTorList[aat1][aat2][1][2];
				//atom4=BasePairTorList[aat1][aat2][1][3];
				//tmp1=decstr[basei].ptv[atom1];tmp2=decstr[basei].ptv[atom2];tmp3=decstr[basej].ptv[atom3];tmp4=decstr[basej].ptv[atom4];
				//dih=bf.phi(tmp1.x,tmp1.y,tmp1.z,tmp2.x,tmp2.y,tmp2.z,tmp3.x,tmp3.y,tmp3.z,tmp4.x,tmp4.y,tmp4.z);
				//if(dih>360) dih=dih-360;
				//else if(dih<0) dih=dih+360;
				//indf[0]+=getbptorsionene(dih, BASEPAIRPARAM[1][3].boundary, BASEPAIRPARAM[1][3].energy);

				//C1'-C1'
				//atom1=BasePairLenList[aat1][aat2][6][0];
				//atom2=BasePairLenList[aat1][aat2][6][1];
				//tmp1=decstr[basei].ptv[atom1];
				//tmp2=decstr[basej].ptv[atom2];
				//pt.x=double(tmp1.x-tmp2.x); pt.y=double(tmp1.y-tmp2.y); pt.z=double(tmp1.z-tmp2.z);
				//len=bf.norm(pt);
				//getoutenelength(len, BASEPAIRPARAM[0][6].mean, BASEPAIRPARAM[0][6].lower, BASEPAIRPARAM[0][6].upper, 1, indftmp, indntmp);
				//indf[2]+=indftmp; indn[2]+=indntmp;

				if((aat1 == 'A' & aat2 == 'U') | (aat1 == 'U' & aat2 == 'A') | (aat1 == 'G' & aat2 == 'C') | (aat1 == 'C' & aat2 == 'G')){
					indf[2]+=bpplaneene(decstr, basei, basej);
				}

				//cout<< 2 <<"-"<< indf[2]<<endl;

				//if(indn[2] == 1) cout<< 2 <<"-"<< indf[2]<<endl;

				//C3'-C3'
				//atom1=BasePairLenList[aat1][aat2][1][0];
				//atom2=BasePairLenList[aat1][aat2][1][1];
				//tmp1=decstr[basei].ptv[atom1];
				//tmp2=decstr[basej].ptv[atom2];
				//pt.x=double(tmp1.x-tmp2.x); pt.y=double(tmp1.y-tmp2.y); pt.z=double(tmp1.z-tmp2.z);
				//len=bf.norm(pt);
				//getoutenelength(len, BASEPAIRPARAM[0][4].mean, BASEPAIRPARAM[0][4].lower, BASEPAIRPARAM[0][4].upper, 1, indftmp, indntmp);
				//indf[0]+=indftmp; indn[0]+=indntmp;

				for(j=0;j<BASE_ENERGY_NUM;j++)
				{
					//cout<<to_string(basei) + "-" + to_string(basej) + "-" + to_string(j) + "-" + to_string(indf[j]) + "-" + to_string(indn[j])<<endl;
					baserpp[pairwiseidx].pval[j]=indf[j];
					baserpp[pairwiseidx].pcw[j]=indn[j];

					lp.indf[BASE_ENERGY_ST+j][basei]+=baserpp[pairwiseidx].pval[j];
					lp.indf[BASE_ENERGY_ST+j][basej]+=baserpp[pairwiseidx].pval[j];
					lp.indn[BASE_ENERGY_ST+j][basei]+=baserpp[pairwiseidx].pcw[j];
					lp.indn[BASE_ENERGY_ST+j][basej]+=baserpp[pairwiseidx].pcw[j];
					
					fnew[j]+=baserpp[pairwiseidx].pval[j];
					inew[j]+=baserpp[pairwiseidx].pcw[j];
				}
			}
		}

		for(i=0;i<BASE_ENERGY_NUM;i++)
		{
			lp.fn[BASE_ENERGY_ST+i]+=(fnew[i]-fold[i]);
			lp.nn[BASE_ENERGY_ST+i]+=(inew[i]-iold[i]);
		}
	}

	fene=0; nene=0;
	for(i=0;i<BASE_ENERGY_NUM;i++)
	{
		fene+=lp.fn[BASE_ENERGY_ST+i];//base
		nene+=lp.nn[BASE_ENERGY_ST+i];//base
	}
}






















