// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

double RNASimulation::getsttorsionene(double val, vector<double> boundary, vector<double> energy){
	int binidx=getbin(val, boundary);
	return energy[binidx];
}

void RNASimulation::getoutenestacklength(double val, double mean, double lower, double upper, double k, double &diff, int &bdiff){
	double tmp;

	if(lower <= val && val <= upper){
		bdiff=0;
		diff=0;
	}else{
		bdiff=1;
		tmp=val-mean;
		diff=k*(tmp*tmp);
	}
}

void RNASimulation::fenergystacking(vector<Residueatoms> decstr,int numseq,double &nene,double &fene,Energyterms &lp, vector<Stackcontact> &stackrpp, vector<vector<bool>> mapflagrpp,bool flagacc){
	int i,j;
	int basei,basej;
	int pairwiseidx;
	double dih, len;
	char aat1, aat2;
	int atom1, atom2, atom3, atom4;
	int atom11,atom12,atom13,atom21,atom22,atom23;
	Point3s tmp1,tmp2,tmp3,tmp4;
	Point3s res1N91, res1C42, res1C86, res2N91, res2C42, res2C86;
	double r, alpha, beta, gamma;
	int binidx1, binidx2, binidx3;
	int binnum1, binnum2, binnum3;
	Point3s pt;
	int stacktype;
	double indftmp;
	int indntmp;
	vector<double> indf(STACKING_ENERGY_NUM,0);
	vector<int> indn(STACKING_ENERGY_NUM,0);
	Point3s basecenter1, basecenter2;
	double baseangle, baseangle1, baseangle2;
	Point3s c2cvec, z1, z2;
	double tmp;

	if(!flagacc)
	{
		for(j=0;j<STACKING_ENERGY_NUM;j++)
		{
			lp.fn[STACKING_ENERGY_ST+j]=0;
			lp.nn[STACKING_ENERGY_ST+j]=0;


			for(i=0;i<numseq;i++){
				lp.indf[STACKING_ENERGY_ST+j][i]=0;
				lp.indf[STACKING_ENERGY_ST+j][i]=0;
				lp.indn[STACKING_ENERGY_ST+j][i]=0;
				lp.indn[STACKING_ENERGY_ST+j][i]=0;
			}
		}
		
		for(i=0;i<STACKINGPAIRS.size();i++)
		{
			basei=STACKINGPAIRS[i].basepidx.first;//base idx 1
			basej=STACKINGPAIRS[i].basepidx.second;//base idx 2
			stacktype=STACKINGPAIRS[i].sttype;

			binnum2=STACKINGROT[stacktype][0].binnum2;
			binnum3=STACKINGROT[stacktype][0].binnum3;

			pairwiseidx=bsmapid2idx[basei*numseq+basej];

			if(basei < basej)
			{
				aat1=decstr[basei].aaa;
				aat2=decstr[basej].aaa;

				for(j=0;j<STACKING_ENERGY_NUM;j++)
				{
					indf[j]=0; indn[j]=0;
				}

				atom11=BASEAXISST[aat1][aat2][0][0];
				atom12=BASEAXISST[aat1][aat2][0][1];
				atom13=BASEAXISST[aat1][aat2][0][2];
				atom21=BASEAXISST[aat1][aat2][0][3];
				atom22=BASEAXISST[aat1][aat2][0][4];
				atom23=BASEAXISST[aat1][aat2][0][5];

				res1N91=decstr[basei].ptv[atom11]; res1C42=decstr[basei].ptv[atom12]; res1C86=decstr[basei].ptv[atom13];
				res2N91=decstr[basej].ptv[atom21]; res2C42=decstr[basej].ptv[atom22]; res2C86=decstr[basej].ptv[atom23];

				//getbaseangles(res1N91, res1C42, res1C86, res2N91, res2C42, res2C86, r, alpha, beta, gamma);
				//getbaselength(res1N91, res2N91,r);

				//getoutenestacklength(r, STACKINGPARAM[stacktype][0][1].mean, STACKINGPARAM[stacktype][0][1].lower, STACKINGPARAM[stacktype][0][1].upper, 1, indftmp, indntmp);
				//indf[0]+=indftmp; indn[0]=indntmp;

				//if(r > 5.5){
					//tmp=r-5.5;
					//indf[0]+=sqrt(tmp*tmp);
					//indn[0]+=1;
				//}

				getnormalvec(res1N91, res1C42, res1C86, z1);
				getnormalvec(res2N91, res2C42, res2C86, z2);
				baseangle = angleBetweenVectors(z1, z2);
				if(baseangle > 90) baseangle=180-baseangle;
				if(baseangle > 30){
					tmp=baseangle-30;
					indf[1]+=0.5*sqrt(tmp*tmp);
					indn[1]+=1;
				}

				getbasecenterpoint(decstr[basei].ptv, aat1, basecenter1);
				getbasecenterpoint(decstr[basej].ptv, aat2, basecenter2);
				getcenter2centervec(basecenter1, basecenter2, c2cvec);
				baseangle1 = angleBetweenVectors(z1, c2cvec);
				baseangle2 = angleBetweenVectors(z2, c2cvec);

				if(baseangle1 > 90) baseangle1 = 180-baseangle1;
				if(baseangle1 > 40){
					tmp=baseangle-40;
					indf[2]+=0.5*sqrt(tmp*tmp);
					indn[2]=1;
				}

				if(baseangle2 > 90) baseangle2 = 180-baseangle2;
				if(baseangle2 > 40){
					tmp=baseangle-40;
					indf[2]+=0.5*sqrt(tmp*tmp);
					indn[2]=1;
				}

				r=bf.norm(c2cvec);
				if(r > 5.5){
					tmp=r-5.5;
					indf[0]+=sqrt(tmp*tmp);
					indn[0]+=1;
				}

				//if(baseangle > 90) baseangle=180-baseangle;
				//getouteneangle(baseangle, STACKINGPARAM[stacktype][2][0].mean, STACKINGPARAM[stacktype][2][0].lower, STACKINGPARAM[stacktype][2][0].upper, 1, indftmp, indntmp);
				//indf[1]=indftmp; indn[1]=indntmp;

				//if(indn[0] == 1){
					//cout<<"a1-" + to_string(r) + "-" + to_string(STACKINGPARAM[stacktype][0][1].lower) + "-" + to_string(STACKINGPARAM[stacktype][0][1].upper) + "-" + to_string(indf[0]) + "-" + to_string(indn[0])<<endl;
				//}else{
					//cout<<"b1-" + to_string(r) + "-" + to_string(STACKINGPARAM[stacktype][0][1].lower) + "-" + to_string(STACKINGPARAM[stacktype][0][1].upper) + "-" + to_string(indf[0]) + "-" + to_string(indn[0])<<endl;
				//}

				//binidx1=getbin(alpha, STACKINGROT[stacktype][0].boundary1);
				//binidx2=getbin(beta, STACKINGROT[stacktype][0].boundary2);
				//binidx3=getbin(gamma, STACKINGROT[stacktype][0].boundary3);

				//indftmp=STACKINGROT[stacktype][0].energy[binidx1*binnum2*binnum3+binidx2*binnum3+binidx3];
				//indntmp=STACKINGROT[stacktype][0].binary[binidx1*binnum2*binnum3+binidx2*binnum3+binidx3];
				//indf[1]+=indftmp; indn[1]+=indntmp;

				//if(indn[1] == 1){
					//cout<<"a2-" + to_string(r) + "-" + to_string(alpha) + "-" + to_string(beta) + "-" + to_string(gamma) + "-" + to_string(indf[1]) + "-" + to_string(indn[1])<<endl;
				//}else{
					//cout<<"b2-" + to_string(r) + "-" + to_string(alpha) + "-" + to_string(beta) + "-" + to_string(gamma) + "-" + to_string(indf[1]) + "-" + to_string(indn[1])<<endl;
				//}

				//torsion C1'-N9-N9-C1'
				//atom1=BaseStackingTorList[aat1][aat2][0][0];
				//atom2=BaseStackingTorList[aat1][aat2][0][1];
				//atom3=BaseStackingTorList[aat1][aat2][0][2];
				//atom4=BaseStackingTorList[aat1][aat2][0][3];
				//tmp1=decstr[basei].ptv[atom1];tmp2=decstr[basei].ptv[atom2];tmp3=decstr[basej].ptv[atom3];tmp4=decstr[basej].ptv[atom4];
				//dih=bf.phi(tmp1.x,tmp1.y,tmp1.z,tmp2.x,tmp2.y,tmp2.z,tmp3.x,tmp3.y,tmp3.z,tmp4.x,tmp4.y,tmp4.z);
				//if(dih>360) dih=dih-360;
				//else if(dih<0) dih=dih+360;
				//getoutenetorangle(dih, STACKINGPARAM[stacktype][1][0].mean, STACKINGPARAM[stacktype][1][0].lower, STACKINGPARAM[stacktype][1][0].upper, 1, indftmp, indntmp);
				//indftmp=getsttorsionene(dih, STACKINGPARAM[stacktype][1][0].boundary, STACKINGPARAM[stacktype][1][0].energy);
				//indf[0]+=indftmp;

				//cout<<to_string(basei) + "-" + to_string(basej) + "-" + to_string(stacktype) + "-" + to_string(dih) + "-" + to_string(indftmp)<<endl;

				//torsion O3'-C3'-C3'-O3'
				//atom1=BaseStackingTorList[aat1][aat2][1][0];
				//atom2=BaseStackingTorList[aat1][aat2][1][1];
				//atom3=BaseStackingTorList[aat1][aat2][1][2];
				//atom4=BaseStackingTorList[aat1][aat2][1][3];
				//tmp1=decstr[basei].ptv[atom1];tmp2=decstr[basei].ptv[atom2];tmp3=decstr[basej].ptv[atom3];tmp4=decstr[basej].ptv[atom4];
				//dih=bf.phi(tmp1.x,tmp1.y,tmp1.z,tmp2.x,tmp2.y,tmp2.z,tmp3.x,tmp3.y,tmp3.z,tmp4.x,tmp4.y,tmp4.z);
				//if(dih>360) dih=dih-360;
				//else if(dih<0) dih=dih+360;
				//indftmp=getsttorsionene(dih, STACKINGPARAM[stacktype][1][3].boundary, STACKINGPARAM[stacktype][1][3].energy);
				//indf[0]+=indftmp;

				//N9-N9
				//atom1=BaseStackingLenList[aat1][aat2][0][0];
				//atom2=BaseStackingLenList[aat1][aat2][0][1];
				//tmp1=decstr[basei].ptv[atom1];
				//tmp2=decstr[basej].ptv[atom2];
				//pt.x=double(tmp1.x-tmp2.x); pt.y=double(tmp1.y-tmp2.y); pt.z=double(tmp1.z-tmp2.z);
				//len=bf.norm(pt);
				//getoutenestacklength(len, STACKINGPARAM[stacktype][0][0].mean, STACKINGPARAM[stacktype][0][0].lower, STACKINGPARAM[stacktype][0][0].upper, 1, indftmp, indntmp);
				//getoutenelength(len, STACKINGPARAM[stacktype][0][7].mean, STACKINGPARAM[stacktype][0][7].lower, STACKINGPARAM[stacktype][0][7].upper, 1, indftmp, indntmp);
				//indf[0]+=indftmp; indn[0]+=indntmp;

				//C3'-C3'
				//atom1=BaseStackingLenList[aat1][aat2][1][0];
				//atom2=BaseStackingLenList[aat1][aat2][1][1];
				//tmp1=decstr[basei].ptv[atom1];
				//tmp2=decstr[basej].ptv[atom2];
				//pt.x=double(tmp1.x-tmp2.x); pt.y=double(tmp1.y-tmp2.y); pt.z=double(tmp1.z-tmp2.z);
				//len=bf.norm(pt);
				//getoutenestacklength(len, STACKINGPARAM[stacktype][0][3].mean, STACKINGPARAM[stacktype][0][3].lower, STACKINGPARAM[stacktype][0][3].upper, 1, indftmp, indntmp);
				//indf[0]+=indftmp; indn[0]+=indntmp;
				
				for(j=0;j<STACKING_ENERGY_NUM;j++)
				{
					stackrpp[pairwiseidx].pval[j]=indf[j];
					stackrpp[pairwiseidx].pval[j]=indf[j];
					stackrpp[pairwiseidx].pcw[j]=indn[j];
					stackrpp[pairwiseidx].pcw[j]=indn[j];

					lp.indf[STACKING_ENERGY_ST+j][basei]+=indf[j];
					lp.indf[STACKING_ENERGY_ST+j][basej]+=indf[j];
					lp.indn[STACKING_ENERGY_ST+j][basei]+=indn[j];
					lp.indn[STACKING_ENERGY_ST+j][basej]+=indn[j];
					lp.fn[STACKING_ENERGY_ST+j]+=indf[j];
					lp.nn[STACKING_ENERGY_ST+j]+=indn[j];
				}
			}
		}
	}
	else
	{
		//int idx = 0;
		int iold[STACKING_ENERGY_NUM],inew[STACKING_ENERGY_NUM];
		double fold[STACKING_ENERGY_NUM],fnew[STACKING_ENERGY_NUM];
		for(i=0;i<STACKING_ENERGY_NUM;i++){
			iold[i]=0;inew[i]=0;
			fold[i]=0;fnew[i]=0;
		}


		for(i=0;i<STACKINGPAIRS.size();i++)
		{
			basei=STACKINGPAIRS[i].basepidx.first;//base idx 1
			basej=STACKINGPAIRS[i].basepidx.second;//base idx 2
			stacktype=STACKINGPAIRS[i].sttype;

			binnum2=STACKINGROT[stacktype][0].binnum2;
			binnum3=STACKINGROT[stacktype][0].binnum3;
			
			pairwiseidx=bsmapid2idx[basei*numseq+basej];

			if(basei<basej)
			{
				if(!mapflagrpp[basei][basej]) continue;
			
				for(j=0;j<STACKING_ENERGY_NUM;j++)
				{
					fold[j]+=stackrpp[pairwiseidx].pval[j];
					iold[j]+=stackrpp[pairwiseidx].pcw[j];

					lp.indf[STACKING_ENERGY_ST+j][basei]=lp.indf[STACKING_ENERGY_ST+j][basei]-stackrpp[pairwiseidx].pval[j];
					lp.indf[STACKING_ENERGY_ST+j][basej]=lp.indf[STACKING_ENERGY_ST+j][basej]-stackrpp[pairwiseidx].pval[j];
					lp.indn[STACKING_ENERGY_ST+j][basei]=lp.indn[STACKING_ENERGY_ST+j][basei]-stackrpp[pairwiseidx].pcw[j];
					lp.indn[STACKING_ENERGY_ST+j][basej]=lp.indn[STACKING_ENERGY_ST+j][basej]-stackrpp[pairwiseidx].pcw[j];

					stackrpp[pairwiseidx].pval[j]=0;
					stackrpp[pairwiseidx].pcw[j]=0;
				}

				aat1=decstr[basei].aaa;
				aat2=decstr[basej].aaa;

				for(j=0;j<STACKING_ENERGY_NUM;j++)
				{
					indf[j]=0; indn[j]=0;
				}

				atom11=BASEAXISST[aat1][aat2][0][0];
				atom12=BASEAXISST[aat1][aat2][0][1];
				atom13=BASEAXISST[aat1][aat2][0][2];
				atom21=BASEAXISST[aat1][aat2][0][3];
				atom22=BASEAXISST[aat1][aat2][0][4];
				atom23=BASEAXISST[aat1][aat2][0][5];

				res1N91=decstr[basei].ptv[atom11]; res1C42=decstr[basei].ptv[atom12]; res1C86=decstr[basei].ptv[atom13];
				res2N91=decstr[basej].ptv[atom21]; res2C42=decstr[basej].ptv[atom22]; res2C86=decstr[basej].ptv[atom23];

				//getbaseangles(res1N91, res1C42, res1C86, res2N91, res2C42, res2C86, r, alpha, beta, gamma);
				//getbaselength(res1N91, res2N91, r);

				//getoutenestacklength(r, STACKINGPARAM[stacktype][0][1].mean, STACKINGPARAM[stacktype][0][1].lower, STACKINGPARAM[stacktype][0][1].upper, 1, indftmp, indntmp);
				//indf[0]=indftmp; indn[0]=indntmp;

				//cout<< 0 << "-" << indf[0]<<endl;
				//if(r > 5.5){
					//tmp=r-5.5;
					//indf[0]+=sqrt(tmp*tmp);
					//indn[0]+=1;
				//}

				getnormalvec(res1N91, res1C42, res1C86, z1);
				getnormalvec(res2N91, res2C42, res2C86, z2);
				baseangle = angleBetweenVectors(z1, z2);
				if(baseangle > 90) baseangle=180-baseangle;
				if(baseangle > 30){
					tmp=baseangle-30;
					indf[1]+=0.5*sqrt(tmp*tmp);
					indn[1]+=1;
				}

				getbasecenterpoint(decstr[basei].ptv, aat1, basecenter1);
				getbasecenterpoint(decstr[basej].ptv, aat2, basecenter2);
				getcenter2centervec(basecenter1, basecenter2, c2cvec);
				baseangle1 = angleBetweenVectors(z1, c2cvec);
				baseangle2 = angleBetweenVectors(z2, c2cvec);

				if(baseangle1 > 90) baseangle1 = 180-baseangle1;
				if(baseangle1 > 40){
					tmp=baseangle-40;
					indf[2]+=0.5*sqrt(tmp*tmp);
					indn[2]=1;
				}

				if(baseangle2 > 90) baseangle2 = 180-baseangle2;
				if(baseangle2 > 40){
					tmp=baseangle-40;
					indf[2]+=0.5*sqrt(tmp*tmp);
					indn[2]=1;
				}

				r=bf.norm(c2cvec);
				if(r > 5.5){
					tmp=r-5.5;
					indf[0]+=sqrt(tmp*tmp);
					indn[0]+=1;
				}

				//getbasecenterpoint(decstr[basei].ptv, aat1, basecenter1);
				//getbasecenterpoint(decstr[basej].ptv, aat2, basecenter2);

				//getnormalvecangle(basecenter1, res1C42, res1C86, basecenter2, res2C42, res2C86, baseangle);

				//if(baseangle > 90) baseangle=180-baseangle;
				//getouteneangle(baseangle, STACKINGPARAM[stacktype][2][0].mean, STACKINGPARAM[stacktype][2][0].lower, STACKINGPARAM[stacktype][2][0].upper, 1, indftmp, indntmp);
				//indf[1]=indftmp; indn[1]=indntmp;

				//if(indntmp == 1){
					//cout<<"aaa"<<endl;
					//cout<< 1 << "-" <<baseangle << "-" << STACKINGPARAM[stacktype][2][0].lower << "-" << STACKINGPARAM[stacktype][2][0].upper << "-" << indf[1]<<endl;
				//}
				//if(indn[0] == 1){
					//cout<<"a1-" + to_string(r) + "-" + to_string(STACKINGPARAM[stacktype][0][1].lower) + "-" + to_string(STACKINGPARAM[stacktype][0][1].upper) + "-" + to_string(indf[0]) + "-" + to_string(indn[0])<<endl;
				//}else{
					//cout<<"b1-" + to_string(r) + "-" + to_string(STACKINGPARAM[stacktype][0][1].lower) + "-" + to_string(STACKINGPARAM[stacktype][0][1].upper) + "-" + to_string(indf[0]) + "-" + to_string(indn[0])<<endl;
				//}

				//binidx1=getbin(alpha, STACKINGROT[stacktype][0].boundary1);
				//binidx2=getbin(beta, STACKINGROT[stacktype][0].boundary2);
				//binidx3=getbin(gamma, STACKINGROT[stacktype][0].boundary3);

				//indftmp=STACKINGROT[stacktype][0].energy[binidx1*binnum2*binnum3+binidx2*binnum3+binidx3];
				//indntmp=STACKINGROT[stacktype][0].binary[binidx1*binnum2*binnum3+binidx2*binnum3+binidx3];
				//indf[1]+=indftmp; indn[1]+=indntmp;

				//cout<< 1 << "-" << indf[1]<<endl;

				//torsion C1'-N9-N9-C1'
				//atom1=BaseStackingTorList[aat1][aat2][0][0];
				//atom2=BaseStackingTorList[aat1][aat2][0][1];
				//atom3=BaseStackingTorList[aat1][aat2][0][2];
				//atom4=BaseStackingTorList[aat1][aat2][0][3];
				//tmp1=decstr[basei].ptv[atom1];tmp2=decstr[basei].ptv[atom2];tmp3=decstr[basej].ptv[atom3];tmp4=decstr[basej].ptv[atom4];
				//dih=bf.phi(tmp1.x,tmp1.y,tmp1.z,tmp2.x,tmp2.y,tmp2.z,tmp3.x,tmp3.y,tmp3.z,tmp4.x,tmp4.y,tmp4.z);
				//if(dih>360) dih=dih-360;
				//else if(dih<0) dih=dih+360;
				//indftmp=getsttorsionene(dih, STACKINGPARAM[stacktype][1][0].boundary, STACKINGPARAM[stacktype][1][0].energy);
				//indf[0]+=indftmp;

				//torsion O3'-C3'-C3'-O3'
				//atom1=BaseStackingTorList[aat1][aat2][1][0];
				//atom2=BaseStackingTorList[aat1][aat2][1][1];
				//atom3=BaseStackingTorList[aat1][aat2][1][2];
				//atom4=BaseStackingTorList[aat1][aat2][1][3];
				//tmp1=decstr[basei].ptv[atom1];tmp2=decstr[basei].ptv[atom2];tmp3=decstr[basej].ptv[atom3];tmp4=decstr[basej].ptv[atom4];
				//dih=bf.phi(tmp1.x,tmp1.y,tmp1.z,tmp2.x,tmp2.y,tmp2.z,tmp3.x,tmp3.y,tmp3.z,tmp4.x,tmp4.y,tmp4.z);
				//if(dih>360) dih=dih-360;
				//else if(dih<0) dih=dih+360;
				//indftmp=getsttorsionene(dih, STACKINGPARAM[stacktype][1][3].boundary, STACKINGPARAM[stacktype][1][3].energy);
				//indf[0]+=indftmp;

				//N9-N9
				//atom1=BaseStackingLenList[aat1][aat2][0][0];
				//atom2=BaseStackingLenList[aat1][aat2][0][1];
				//tmp1=decstr[basei].ptv[atom1];
				//tmp2=decstr[basej].ptv[atom2];
				//pt.x=double(tmp1.x-tmp2.x); pt.y=double(tmp1.y-tmp2.y); pt.z=double(tmp1.z-tmp2.z);
				//len=bf.norm(pt);
				//getoutenelength(len, STACKINGPARAM[stacktype][0][7].mean, STACKINGPARAM[stacktype][0][7].lower, STACKINGPARAM[stacktype][0][7].upper, 1, indftmp, indntmp);
				//getoutenestacklength(len, STACKINGPARAM[stacktype][0][0].mean, STACKINGPARAM[stacktype][0][0].lower, STACKINGPARAM[stacktype][0][0].upper, 1, indftmp, indntmp);
				//indf[0]+=indftmp; indn[0]+=indntmp;

				//C3'-C3'
				//atom1=BaseStackingLenList[aat1][aat2][1][0];
				//atom2=BaseStackingLenList[aat1][aat2][1][1];
				//tmp1=decstr[basei].ptv[atom1];
				//tmp2=decstr[basej].ptv[atom2];
				//pt.x=double(tmp1.x-tmp2.x); pt.y=double(tmp1.y-tmp2.y); pt.z=double(tmp1.z-tmp2.z);
				//len=bf.norm(pt);
				//getoutenelength(len, STACKINGPARAM[stacktype][0][3].mean, STACKINGPARAM[stacktype][0][3].lower, STACKINGPARAM[stacktype][0][3].upper, 1, indftmp, indntmp);
				//indf[0]+=indftmp; indn[0]+=indntmp;

				for(j=0;j<STACKING_ENERGY_NUM;j++)
				{
					stackrpp[pairwiseidx].pval[j]+=indf[j];
					stackrpp[pairwiseidx].pcw[j]+=indn[j];

					lp.indf[STACKING_ENERGY_ST+j][basei]+=indf[j];
					lp.indf[STACKING_ENERGY_ST+j][basej]+=indf[j];
					lp.indn[STACKING_ENERGY_ST+j][basei]+=indn[j];
					lp.indn[STACKING_ENERGY_ST+j][basej]+=indn[j];

					fnew[j]+=stackrpp[pairwiseidx].pval[j];
					inew[j]+=stackrpp[pairwiseidx].pcw[j];
				}
			}
		}

		for(i=0;i<STACKING_ENERGY_NUM;i++)
		{
			lp.fn[STACKING_ENERGY_ST+i]+=(fnew[i]-fold[i]);
			lp.nn[STACKING_ENERGY_ST+i]+=(inew[i]-iold[i]);
		}
	}

	fene=0;
	nene=0;
	for(i=0;i<STACKING_ENERGY_NUM;i++)
	{
		fene+=lp.fn[STACKING_ENERGY_ST+i];//stacking
		nene+=lp.nn[STACKING_ENERGY_ST+i];//stacking
	}
}






















