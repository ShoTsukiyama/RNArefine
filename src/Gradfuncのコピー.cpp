// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

void RNASimulation::initgrad(map<int, vector<Point3s>> &gradptv, vector<bool> movenodeinfo)
{
	int mtype1;
    for(int i=0;i<movenodeinfo.size();i++){
		gradptv[i].resize(PAIRLENGTHTARGETMOVNUM);
		//if(movenodeinfo[i]){
		for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
		{
			gradptv[i][mtype1].x = 0;
			gradptv[i][mtype1].y = 0;
			gradptv[i][mtype1].z = 0;
		}
		//}
	}
}

void RNASimulation::initializegrad(map<int, vector<Point3s>> &gradptv, vector<bool> movenodeinfo, double &mse)
{
	int mtype1;
    for(int i=0;i<movenodeinfo.size();i++){
		gradptv[i].resize(PAIRLENGTHTARGETMOVNUM);
		//if(movenodeinfo[i]){
			for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
			{
				gradptv[i][mtype1].x = 0;
				gradptv[i][mtype1].y = 0;
				gradptv[i][mtype1].z = 0;
			}
		//}
	}
    mse=0;
}


void RNASimulation::calcmseangle(Point3s P, Point3s Q, Point3s R, double &angle, double angleideal, double anglerangelower, double anglerangeupper, bool &outflag, double &mse, double &theta, double &diff, double &normPQ, double &normRQ, double &dotp, Point3s &vecPQ, Point3s &vecRQ, double coeff){
	int i;
	Point3s PQ;
    Point3s RQ;

	PQ.x = Q.x - P.x; PQ.y = Q.y - P.y; PQ.z = Q.z - P.z;
    RQ.x = Q.x - R.x; RQ.y = Q.y - R.y; RQ.z = Q.z - R.z;
        
	vecPQ.x=PQ.x; vecPQ.y=PQ.y; vecPQ.z=PQ.z;
	vecRQ.x=RQ.x; vecRQ.y=RQ.y; vecRQ.z=RQ.z;

	normPQ = bf.norm(PQ) + EPS;
    normRQ = bf.norm(RQ) + EPS;

    dotp=bf.dotv(PQ,RQ);
    double t3=dotp/normPQ/normRQ;
	if(t3<-1.0)
	{
		t3=-1.0;
	}
	else if(t3>1.0)
	{
		t3=1.0;
	}

    theta=acos(t3);
        
    angle = theta*RAD2DEG;

	//if(angle < anglerangelower | angle > anglerangeupper){
		//diff=angle-angleideal;
    	//mse += coeff*pow(diff, 2);
		//outflag=true;


	
	if(angle < anglerangelower | angle > anglerangeupper){
		diff=angle-angleideal;
		mse += coeff*pow(diff, 2);
		outflag=true;
	}else{
		outflag=false;
	}


	//diff=angle-angleideal;
	//mse += coeff*pow(diff, 2);
	//outflag=true;
}

void RNASimulation::calcmseangles(vector<Residueatoms> restr, int targetnt, vector<double> &angle, vector<double> &diff, vector<bool> &outerflags, vector<double> &theta, vector<double> &normPQs, vector<double> &normRQs, vector<double> &dotps, vector<Point3s> &vecPQs, vector<Point3s> &vecRQs, double &mse, int &count, double coeff){
	int i;
	bool flag;
	Point3s tmp1,tmp2,tmp3;
	char aat=restr[targetnt].aaa;
	int atom1, atom2, atom3;
	int repidx;

	//cout<< targetnt << endl;

	if(targetnt != 0){
		for(i=0; i<AccAtomParamNum[aat][2]; i++){
		//for(i=0; i<2; i++){
			atom1=ParamsList[aat][i][1];atom2=ParamsList[aat][i][2];atom3=ParamsList[aat][i][3];
			if(i == 0){
				tmp1=restr[targetnt-1].ptv[atom1];tmp2=restr[targetnt-1].ptv[atom2];tmp3=restr[targetnt].ptv[atom3];
			}else if(i == 1 | i == 11){
				tmp1=restr[targetnt-1].ptv[atom1];tmp2=restr[targetnt].ptv[atom2];tmp3=restr[targetnt].ptv[atom3];
			}else{
				tmp1=restr[targetnt].ptv[atom1];tmp2=restr[targetnt].ptv[atom2];tmp3=restr[targetnt].ptv[atom3];
			}

			calcmseangle(tmp1, tmp2, tmp3, angle[i], INNERANGPARAM[aat][INNERANGLETARGET[i]].mean, INNERANGPARAM[aat][INNERANGLETARGET[i]].lower, INNERANGPARAM[aat][INNERANGLETARGET[i]].upper, flag, mse, theta[i], diff[i], normPQs[i], normRQs[i], dotps[i], vecPQs[i], vecRQs[i], coeff);
			
			//if(flag){
				//cout<< "aaa: " << i << "-" << angle[i] << "-" << INNERANGPARAM[aat][INNERANGLETARGET[i]].mean<<endl;
				//cout<< atom1 << "-" << atom2 << "-" << atom3<<endl;
			//}
			//cout<< targetnt << "-" << i << "-" << mse<<endl;
			
			outerflags[i]=flag;
			count+=1;
		}

		/*
		for(i=0;i<AddAngleParamsList[aat].size();i++)
		{
			atom1=AddAngleParamsList[aat][i][0]; atom2=AddAngleParamsList[aat][i][1]; atom3=AddAngleParamsList[aat][i][2];
			if(i==0){
				tmp1=restr[targetnt-1].ptv[atom1]; tmp2=restr[targetnt].ptv[atom2]; tmp3=restr[targetnt].ptv[atom3];
			}else{
				tmp1=restr[targetnt].ptv[atom1]; tmp2=restr[targetnt].ptv[atom2]; tmp3=restr[targetnt].ptv[atom3];
			}
			repidx=AccAtomParamNum[aat][2]+i;
			calcmseangle(tmp1, tmp2, tmp3, angle[repidx], ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[i]].mean, ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[i]].lower, ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[i]].upper, flag, mse, theta[repidx], diff[repidx], normPQs[repidx], normRQs[repidx], dotps[repidx], vecPQs[repidx], vecRQs[repidx], coeff);
			
			//cout<< "bbb: " << repidx << "-" << angle[i] << "-" << ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[i]].mean<<endl;
			
			outerflags[repidx]=flag;

			//cout<< targetnt << "-" << i << "-" << mse<<endl;

			count+=1;
		}
		*/
	}else{
		for(int i=0; i<AccAtomParamNum[aat][2]; i++){
		//for(int i=0; i<2; i++){
			outerflags[i]=false;
		}
		//for(i=0;i<AddAngleParamsList[aat].size();i++)
		//{
			//outerflags[AccAtomParamNum[aat][2]+i]=false;
		//}
	}

}

void RNASimulation::calcmsefullangles(vector<Residueatoms> restr, int targetnt, vector<double> &angle, vector<double> &diff, vector<bool> &outerflags, vector<double> &theta, vector<double> &normPQs, vector<double> &normRQs, vector<double> &dotps, vector<Point3s> &vecPQs, vector<Point3s> &vecRQs, double &mse, int &count, double coeff){
	int i;
	bool flag;
	Point3s tmp1,tmp2,tmp3;
	char aat=restr[targetnt].aaa;
	int atom1, atom2, atom3;
	int repidx;

	//cout<< targetnt << endl;

	if(targetnt != 0){
		for(i=0; i<AccAtomParamNum[aat][2]; i++){
			atom1=ParamsList[aat][i][1];atom2=ParamsList[aat][i][2];atom3=ParamsList[aat][i][3];
			if(i == 0){
				tmp1=restr[targetnt-1].ptv[atom1];tmp2=restr[targetnt-1].ptv[atom2];tmp3=restr[targetnt].ptv[atom3];
			}else if(i == 1 | i == 11){
				tmp1=restr[targetnt-1].ptv[atom1];tmp2=restr[targetnt].ptv[atom2];tmp3=restr[targetnt].ptv[atom3];
			}else{
				tmp1=restr[targetnt].ptv[atom1];tmp2=restr[targetnt].ptv[atom2];tmp3=restr[targetnt].ptv[atom3];
			}

			calcmseangle(tmp1, tmp2, tmp3, angle[i], INNERANGPARAM[aat][INNERANGLETARGET[i]].mean, INNERANGPARAM[aat][INNERANGLETARGET[i]].lower, INNERANGPARAM[aat][INNERANGLETARGET[i]].upper, flag, mse, theta[i], diff[i], normPQs[i], normRQs[i], dotps[i], vecPQs[i], vecRQs[i], coeff);
			
			//cout<< i << "-" << targetnt << "-" << angle[i] << "-" << INNERANGPARAM[aat][INNERANGLETARGET[i]].mean << "-" << mse<<endl;
			//if(flag){
				//cout<< "aaa: " << i << "-" << angle[i] << "-" << INNERANGPARAM[aat][INNERANGLETARGET[i]].mean<<endl;
				//cout<< atom1 << "-" << atom2 << "-" << atom3<<endl;
			//}
			//cout<< targetnt << "-" << i << "-" << mse<<endl;
			
			outerflags[i]=flag;
			count+=1;
		}

		
		for(i=0;i<AddAngleParamsList[aat].size();i++)
		{
			atom1=AddAngleParamsList[aat][i][0]; atom2=AddAngleParamsList[aat][i][1]; atom3=AddAngleParamsList[aat][i][2];
			if(i==0){
				tmp1=restr[targetnt-1].ptv[atom1]; tmp2=restr[targetnt].ptv[atom2]; tmp3=restr[targetnt].ptv[atom3];
			}else{
				tmp1=restr[targetnt].ptv[atom1]; tmp2=restr[targetnt].ptv[atom2]; tmp3=restr[targetnt].ptv[atom3];
			}
			repidx=AccAtomParamNum[aat][2]+i;
			calcmseangle(tmp1, tmp2, tmp3, angle[repidx], ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[i]].mean, ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[i]].lower, ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[i]].upper, flag, mse, theta[repidx], diff[repidx], normPQs[repidx], normRQs[repidx], dotps[repidx], vecPQs[repidx], vecRQs[repidx], coeff);
			
			//cout<< "bbb: " << repidx << "-" << angle[i] << "-" << ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[i]].mean<<endl;
			
			outerflags[repidx]=flag;

			//cout<< targetnt << "-" << i << "-" << mse<<endl;

			count+=1;
		}
	}else{
		//for(int i=0; i<AccAtomParamNum[aat][2]; i++){
		for(int i=0; i<AccAtomParamNum[aat][2]+AddAngleParamsList[aat].size(); i++){
			outerflags[i]=false;
		}
		//for(i=0;i<AddAngleParamsList[aat].size();i++)
		//{
			//outerflags[AccAtomParamNum[aat][2]+i]=false;
		//}
	}

}

void RNASimulation::calcgradientangle(Point3s P, Point3s Q, Point3s R, map<int, vector<Point3s>> &gradptv, int basep, int baseq, int baser, int mttypep, int mttypeq, int mttyper, bool &flagP, bool &flagQ, bool &flagR, double theta, double diff, double normPQ, double normRQ, double dotp, Point3s vecPQ, Point3s vecRQ, double count, double coeff){
	double ddotp, dnorm;
	double dmse = 2 * (diff);
	double normPQRQ2=pow(normPQ * normRQ,2);
    double dcos=(-1/sin(theta));

	double normRQ_PQ=(normRQ/normPQ);
    double normPQ_RQ=(normPQ/normRQ);

	//double k = coeff*dmse*RAD2DEG*dcos/(normPQRQ2*count);
	double k = coeff*dmse*RAD2DEG*dcos/(normPQRQ2);

	if(flagP){
		ddotp=-vecRQ.x;
        dnorm=normRQ_PQ * (-vecPQ.x);
        gradptv[basep][mttypep].x += k * (ddotp*(normPQ * normRQ) - dotp*dnorm);
		ddotp=-vecRQ.y;
        dnorm=normRQ_PQ * (-vecPQ.y);
        gradptv[basep][mttypep].y += k * (ddotp*(normPQ * normRQ) - dotp*dnorm);
		ddotp=-vecRQ.z;
        dnorm=normRQ_PQ * (-vecPQ.z);
        gradptv[basep][mttypep].z += k * (ddotp*(normPQ * normRQ) - dotp*dnorm);
	}else if(flagQ){
		ddotp=vecPQ.x+vecRQ.x;
        dnorm=normRQ_PQ * (vecPQ.x) +  normPQ_RQ * (vecRQ.x);
        gradptv[baseq][mttypeq].x += k * (ddotp*(normPQ * normRQ) - dotp*dnorm);
		ddotp=vecPQ.y+vecRQ.y;
        dnorm=normRQ_PQ * (vecPQ.y) +  normPQ_RQ * (vecRQ.y);
        gradptv[baseq][mttypeq].y += k * (ddotp*(normPQ * normRQ) - dotp*dnorm);
		ddotp=vecPQ.z+vecRQ.z;
        dnorm=normRQ_PQ * (vecPQ.z) +  normPQ_RQ * (vecRQ.z);
        gradptv[baseq][mttypeq].z += k * (ddotp*(normPQ * normRQ) - dotp*dnorm);
	}else if(flagR){
		ddotp=-vecPQ.x;
        dnorm=normPQ_RQ * (-vecRQ.x);
        gradptv[baser][mttyper].x += k * (ddotp*(normPQ * normRQ) - dotp*dnorm);
		ddotp=-vecPQ.y;
        dnorm=normPQ_RQ * (-vecRQ.y);
        gradptv[baser][mttyper].y += k * (ddotp*(normPQ * normRQ) - dotp*dnorm);
		ddotp=-vecPQ.z;
        dnorm=normPQ_RQ * (-vecRQ.z);
        gradptv[baser][mttyper].z += k * (ddotp*(normPQ * normRQ) - dotp*dnorm);
	}
}

void RNASimulation::calcgradientangles(vector<Residueatoms> restr, int targetnt, map<int, vector<Point3s>> &gradptv, vector<bool> movenodeinfo, vector<double> diff, vector<bool> outerflags, vector<double> theta, vector<double> normPQs, vector<double> normRQs, vector<double> dotps, vector<Point3s> vecPQs, vector<Point3s> vecRQs, int count, double coeff){
	int basep, baseq, baser;
	bool flagp, flagq, flagr;
	int i;
	Point3s tmp1,tmp2,tmp3;
	char aat=restr[targetnt].aaa;
	int atom1, atom2, atom3;
	int repidx;

	if(targetnt != 0){
		for(i=0; i<AccAtomParamNum[aat][2]; i++){
		//for(i=0; i<2; i++){
			if(outerflags[i]){
				atom1=ParamsList[aat][i][1];atom2=ParamsList[aat][i][2];atom3=ParamsList[aat][i][3];
				if(i == 0){
					tmp1=restr[targetnt-1].ptv[atom1];tmp2=restr[targetnt-1].ptv[atom2];tmp3=restr[targetnt].ptv[atom3];
					basep=targetnt-1; baseq=targetnt-1; baser=targetnt;
				}else if(i == 1 | i == 11){
					tmp1=restr[targetnt-1].ptv[atom1];tmp2=restr[targetnt].ptv[atom2];tmp3=restr[targetnt].ptv[atom3];
					basep=targetnt-1; baseq=targetnt; baser=targetnt;
				}else{
					tmp1=restr[targetnt].ptv[atom1];tmp2=restr[targetnt].ptv[atom2];tmp3=restr[targetnt].ptv[atom3];
					basep=targetnt; baseq=targetnt; baser=targetnt;
				}	
				flagp=movenodeinfo[basep];
				flagq=movenodeinfo[baseq];
				flagr=movenodeinfo[baser];

				//if(i==0){
					//cout<<"bbb " << flagp << "-" << flagq << "-" << flagr<<endl;
					//cout<<"bbb " << basep << "-" << baseq << "-" << baser<<endl;
				//}

				calcgradientangle(tmp1, tmp2, tmp3, gradptv, basep, baseq, baser, atom1, atom2, atom3, flagp, flagq, flagr, theta[i], diff[i], normPQs[i], normRQs[i], dotps[i], vecPQs[i], vecRQs[i], count, coeff);
			}
		}

		/*
		for(i=0;i<AddAngleParamsList[aat].size();i++)
		{
			repidx=AccAtomParamNum[aat][2]+i;
			if(outerflags[repidx]){
				atom1=AddAngleParamsList[aat][i][0]; atom2=AddAngleParamsList[aat][i][1]; atom3=AddAngleParamsList[aat][i][2];
				if(i==0){
					tmp1=restr[targetnt-1].ptv[atom1]; tmp2=restr[targetnt].ptv[atom2]; tmp3=restr[targetnt].ptv[atom3];
					basep=targetnt-1; baseq=targetnt; baser=targetnt;
				}else{
					tmp1=restr[targetnt].ptv[atom1]; tmp2=restr[targetnt].ptv[atom2]; tmp3=restr[targetnt].ptv[atom3];
					basep=targetnt; baseq=targetnt; baser=targetnt;
				}

				flagp=movenodeinfo[basep];
				flagq=movenodeinfo[baseq];
				flagr=movenodeinfo[baser];
				calcgradientangle(tmp1, tmp2, tmp3, gradptv, basep, baseq, baser, atom1, atom2, atom3, flagp, flagq, flagr, theta[repidx], diff[repidx], normPQs[repidx], normRQs[repidx], dotps[repidx], vecPQs[repidx], vecRQs[repidx], count, coeff);
			}
		}
		*/
	}
}

void RNASimulation::calcgradientfullangles(vector<Residueatoms> restr, int targetnt, map<int, vector<Point3s>> &gradptv, vector<bool> movenodeinfo, vector<double> diff, vector<bool> outerflags, vector<double> theta, vector<double> normPQs, vector<double> normRQs, vector<double> dotps, vector<Point3s> vecPQs, vector<Point3s> vecRQs, int count, double coeff){
	int basep, baseq, baser;
	bool flagp, flagq, flagr;
	int i;
	Point3s tmp1,tmp2,tmp3;
	char aat=restr[targetnt].aaa;
	int atom1, atom2, atom3;
	int repidx;

	if(targetnt != 0){
		for(i=0; i<AccAtomParamNum[aat][2]; i++){
		//for(i=0; i<2; i++){
			if(outerflags[i]){
				atom1=ParamsList[aat][i][1];atom2=ParamsList[aat][i][2];atom3=ParamsList[aat][i][3];
				if(i == 0){
					tmp1=restr[targetnt-1].ptv[atom1];tmp2=restr[targetnt-1].ptv[atom2];tmp3=restr[targetnt].ptv[atom3];
					basep=targetnt-1; baseq=targetnt-1; baser=targetnt;
				}else if(i == 1 | i == 11){
					tmp1=restr[targetnt-1].ptv[atom1];tmp2=restr[targetnt].ptv[atom2];tmp3=restr[targetnt].ptv[atom3];
					basep=targetnt-1; baseq=targetnt; baser=targetnt;
				}else{
					tmp1=restr[targetnt].ptv[atom1];tmp2=restr[targetnt].ptv[atom2];tmp3=restr[targetnt].ptv[atom3];
					basep=targetnt; baseq=targetnt; baser=targetnt;
				}	
				flagp=movenodeinfo[basep];
				flagq=movenodeinfo[baseq];
				flagr=movenodeinfo[baser];

				//if(i==0){
					//cout<<"bbb " << flagp << "-" << flagq << "-" << flagr<<endl;
					//cout<<"bbb " << basep << "-" << baseq << "-" << baser<<endl;
				//}

				calcgradientangle(tmp1, tmp2, tmp3, gradptv, basep, baseq, baser, atom1, atom2, atom3, flagp, flagq, flagr, theta[i], diff[i], normPQs[i], normRQs[i], dotps[i], vecPQs[i], vecRQs[i], count, coeff);
			}
		}

		
		for(i=0;i<AddAngleParamsList[aat].size();i++)
		{
			repidx=AccAtomParamNum[aat][2]+i;
			if(outerflags[repidx]){
				atom1=AddAngleParamsList[aat][i][0]; atom2=AddAngleParamsList[aat][i][1]; atom3=AddAngleParamsList[aat][i][2];
				if(i==0){
					tmp1=restr[targetnt-1].ptv[atom1]; tmp2=restr[targetnt].ptv[atom2]; tmp3=restr[targetnt].ptv[atom3];
					basep=targetnt-1; baseq=targetnt; baser=targetnt;
				}else{
					tmp1=restr[targetnt].ptv[atom1]; tmp2=restr[targetnt].ptv[atom2]; tmp3=restr[targetnt].ptv[atom3];
					basep=targetnt; baseq=targetnt; baser=targetnt;
				}

				flagp=movenodeinfo[basep];
				flagq=movenodeinfo[baseq];
				flagr=movenodeinfo[baser];
				calcgradientangle(tmp1, tmp2, tmp3, gradptv, basep, baseq, baser, atom1, atom2, atom3, flagp, flagq, flagr, theta[repidx], diff[repidx], normPQs[repidx], normRQs[repidx], dotps[repidx], vecPQs[repidx], vecRQs[repidx], count, coeff);
			}
		}
	}
}

void RNASimulation::calcmsentlength(Point3s *P, Point3s *Q, char aat,int paramidx, vector<vector<double>> &length, vector<vector<double>> &diff, vector<vector<Point3s>> &dists, double &mse, int &count, vector<vector<bool>> &outlineridx, double coeff, bool showres)
{
	int i;
	int mtype1, mtype2;
	double rnorm;
	Point3s PQ;

	for(i=1;i<AccAtomParamNum[aat][2];i++){
		mtype1=ParamsList[aat][i][2];mtype2=ParamsList[aat][i][3];
		PQ.x = P[mtype1].x - Q[mtype2].x; 
		PQ.y = P[mtype1].y - Q[mtype2].y; 
		PQ.z = P[mtype1].z - Q[mtype2].z;

		dists[mtype1][mtype2].x=PQ.x;
		dists[mtype1][mtype2].y=PQ.y;
		dists[mtype1][mtype2].z=PQ.z;

		rnorm=bf.norm(PQ);

		length[mtype1][mtype2] = rnorm;

		if(length[mtype1][mtype2] < LENGTHPARAM[aat][BLENGTHTARGET[i]].lower | length[mtype1][mtype2] > LENGTHPARAM[aat][BLENGTHTARGET[i]].upper){
			diff[mtype1][mtype2] = length[mtype1][mtype2]-LENGTHPARAM[aat][BLENGTHTARGET[i]].mean;
			
			//cout<< length[mtype1][mtype2] << "-" << LENGTHPARAM[aat][BLENGTHTARGET[i]].mean<<endl;
			mse += coeff*diff[mtype1][mtype2]*diff[mtype1][mtype2];
			if(showres){
				cout<<"----------"<<endl;
				cout<<rnorm<<endl;
				cout<<LENGTHPARAM[aat][BLENGTHTARGET[i]].mean<<endl;
				cout<<diff[mtype1][mtype2]*diff[mtype1][mtype2]<<endl;
			}
			outlineridx[mtype1][mtype2]=true;
		}
		count+=1;
	}

	for(i=0;i<AddLengthParamsList[aat].size();i++)
	{
		mtype1=AddLengthParamsList[aat][i][0]; mtype2=AddLengthParamsList[aat][i][1];
		PQ.x = P[mtype1].x - Q[mtype2].x; 
		PQ.y = P[mtype1].y - Q[mtype2].y; 
		PQ.z = P[mtype1].z - Q[mtype2].z;

		dists[mtype1][mtype2].x=PQ.x;
		dists[mtype1][mtype2].y=PQ.y;
		dists[mtype1][mtype2].z=PQ.z;

		rnorm=bf.norm(PQ);

		length[mtype1][mtype2] = rnorm;

		//if(length[mtype1][mtype2] < ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[i]].lower | length[mtype1][mtype2] > ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[i]].upper){
			diff[mtype1][mtype2] = length[mtype1][mtype2]-ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[i]].mean;
			
			//cout<< length[mtype1][mtype2] << "-" << ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[i]].mean<<endl;
			//BPATOMPAIRLENGNTPARAM[aat][paramidx][mtype1][mtype2].lower;
			mse += coeff*diff[mtype1][mtype2]*diff[mtype1][mtype2];
			if(showres){
				cout<<"----------"<<endl;
				cout<<rnorm<<endl;
				cout<<ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[i]].mean<<endl;
				cout<<diff[mtype1][mtype2]*diff[mtype1][mtype2]<<endl;
			}
			outlineridx[mtype1][mtype2]=true;
		//}
		count+=1;
	}

	/*
	for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
	{
		for(mtype2=mtype1+1;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
		{

			PQ.x = P[mtype1].x - Q[mtype2].x; 
			PQ.y = P[mtype1].y - Q[mtype2].y; 
			PQ.z = P[mtype1].z - Q[mtype2].z;

			dists[mtype1][mtype2].x=PQ.x;
			dists[mtype1][mtype2].y=PQ.y;
			dists[mtype1][mtype2].z=PQ.z;

			rnorm=bf.norm(PQ);
	*/
			/*
			if(rnorm<BPATOMPAIRLENGNTPARAM[aat][paramidx][mtype1][mtype2].lower){
				length[mtype1][mtype2] = rnorm;
				diff[mtype1][mtype2] = length[mtype1][mtype2]-BPATOMPAIRLENGNTPARAM[aat][paramidx][mtype1][mtype2].lower;
				mse += coeff*diff[mtype1][mtype2]*diff[mtype1][mtype2];
				if(showres){
					cout<<"----------"<<endl;
					cout<<rnorm<<endl;
					cout<<BPATOMPAIRLENGNTPARAM[aat][paramidx][mtype1][mtype2].lower<<endl;
					cout<<diff[mtype1][mtype2]*diff[mtype1][mtype2]<<endl;
				}
				outlineridx[mtype1][mtype2]=true;
			}else if(rnorm>BPATOMPAIRLENGNTPARAM[aat][paramidx][mtype1][mtype2].upper){
				length[mtype1][mtype2] = rnorm;
				diff[mtype1][mtype2] = length[mtype1][mtype2]-BPATOMPAIRLENGNTPARAM[aat][paramidx][mtype1][mtype2].upper;
				mse += coeff*diff[mtype1][mtype2]*diff[mtype1][mtype2];
				if(showres){
					cout<<"----------"<<endl;
					cout<<rnorm<<endl;
					cout<<BPATOMPAIRLENGNTPARAM[aat][paramidx][mtype1][mtype2].lower<<endl;
					cout<<diff[mtype1][mtype2]*diff[mtype1][mtype2]<<endl;
				}
				outlineridx[mtype1][mtype2]=true;
			}
			*/

	/*
			length[mtype1][mtype2] = rnorm;
			diff[mtype1][mtype2] = length[mtype1][mtype2]-BPATOMPAIRLENGNTPARAM[aat][paramidx][mtype1][mtype2].lower;
			mse += coeff*diff[mtype1][mtype2]*diff[mtype1][mtype2];
			if(showres){
				cout<<"----------"<<endl;
				cout<<rnorm<<endl;
				cout<<BPATOMPAIRLENGNTPARAM[aat][paramidx][mtype1][mtype2].lower<<endl;
				cout<<diff[mtype1][mtype2]*diff[mtype1][mtype2]<<endl;
			}
			outlineridx[mtype1][mtype2]=true;

			count+=1;
		}
		
	}
	*/
}

void RNASimulation::calcmseadjlength(Point3s *P, Point3s *Q, int basei, int basej, char aat1, char aat2, int inda, int indb,int paramidx, vector<vector<double>> &length, vector<vector<double>> &diff, vector<vector<Point3s>> &dists, double &mse, int &count, vector<vector<bool>> &outlineridx, double coeff, bool showres)
{
	int mtype1, mtype2;
	double rnorm;
	Point3s PQ;
	double dspair2;
	int atomidx1, atomidx2;

	for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
	{
		for(mtype2=0;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
		{
			
			
			if(basei < basej & mtype1 == 5 & mtype2 == 0){
				PQ.x = P[mtype1].x - Q[mtype2].x; 
				PQ.y = P[mtype1].y - Q[mtype2].y; 
				PQ.z = P[mtype1].z - Q[mtype2].z;

				dists[mtype1][mtype2].x=PQ.x;
				dists[mtype1][mtype2].y=PQ.y;
				dists[mtype1][mtype2].z=PQ.z;

				rnorm=bf.norm(PQ);
				length[mtype1][mtype2] = rnorm;
				if(length[mtype1][mtype2] < LENGTHPARAM[aat2][0].lower | length[mtype1][mtype2] > LENGTHPARAM[aat2][0].upper){
					diff[mtype1][mtype2] = length[mtype1][mtype2]-LENGTHPARAM[aat2][0].mean;
					
					mse += coeff*diff[mtype1][mtype2]*diff[mtype1][mtype2];
					if(showres){
						cout<<"----------"<<endl;
						cout<<rnorm<<endl;
						cout<<LENGTHPARAM[aat2][0].mean<<endl;
						cout<<diff[mtype1][mtype2]*diff[mtype1][mtype2]<<endl;
					}
					outlineridx[mtype1][mtype2]=true;
				}
				count+=1;
			}else if(basei > basej & mtype1 == 0 & mtype2 == 5){
				PQ.x = P[mtype1].x - Q[mtype2].x; 
				PQ.y = P[mtype1].y - Q[mtype2].y; 
				PQ.z = P[mtype1].z - Q[mtype2].z;

				dists[mtype1][mtype2].x=PQ.x;
				dists[mtype1][mtype2].y=PQ.y;
				dists[mtype1][mtype2].z=PQ.z;

				rnorm=bf.norm(PQ);
				length[mtype1][mtype2] = rnorm;
				if(length[mtype1][mtype2] < LENGTHPARAM[aat1][0].lower | length[mtype1][mtype2] > LENGTHPARAM[aat1][0].upper){
					diff[mtype1][mtype2] = length[mtype1][mtype2]-LENGTHPARAM[aat1][0].mean;
					
					mse += coeff*diff[mtype1][mtype2]*diff[mtype1][mtype2];
					if(showres){
						cout<<"----------"<<endl;
						cout<<rnorm<<endl;
						cout<<LENGTHPARAM[aat1][0].mean<<endl;
						cout<<diff[mtype1][mtype2]*diff[mtype1][mtype2]<<endl;
					}
					outlineridx[mtype1][mtype2]=true;
				}
				count+=1;
			}
			
			
			
			/*
			if(EXPCLASHATOMSADJ[aat1][aat2][mtype1][mtype2]) continue;

			PQ.x = P[mtype1].x - Q[mtype2].x; 
			PQ.y = P[mtype1].y - Q[mtype2].y; 
			PQ.z = P[mtype1].z - Q[mtype2].z;

			dists[mtype1][mtype2].x=PQ.x;
			dists[mtype1][mtype2].y=PQ.y;
			dists[mtype1][mtype2].z=PQ.z;

			rnorm=bf.norm(PQ);
			atomidx1=NNheavy[inda][mtype1+1];
			atomidx2=NNheavy[indb][mtype2+1];
			dspair2=VWPARAM[atomidx1]+VWPARAM[atomidx2]-0.4;

			if(rnorm<dspair2){
				length[mtype1][mtype2] = rnorm;
				diff[mtype1][mtype2] = length[mtype1][mtype2]-dspair2;
				mse += coeff*diff[mtype1][mtype2]*diff[mtype1][mtype2];
				if(showres){
					cout<<"----------"<<endl;
					cout<<rnorm<<endl;
					cout<<dspair2<<endl;
					cout<<diff[mtype1][mtype2]*diff[mtype1][mtype2]<<endl;
				}
				outlineridx[mtype1][mtype2]=true;
			}
			count+=1;
			*/
			


			
			
			/*
			if(rnorm<BPATOMPAIRLENGADJPARAM[aat1][aat2][paramidx][mtype1][mtype2].lower | rnorm>BPATOMPAIRLENGADJPARAM[aat1][aat2][paramidx][mtype1][mtype2].upper){
				length[mtype1][mtype2] = rnorm;
				diff[mtype1][mtype2] = length[mtype1][mtype2]-BPATOMPAIRLENGADJPARAM[aat1][aat2][paramidx][mtype1][mtype2].mean;
				mse += coeff*diff[mtype1][mtype2]*diff[mtype1][mtype2];
				if(showres){
					cout<<"----------"<<endl;
					cout<<rnorm<<endl;
					cout<<BPATOMPAIRLENGADJPARAM[aat1][aat2][paramidx][mtype1][mtype2].lower<<endl;
					cout<<diff[mtype1][mtype2]*diff[mtype1][mtype2]<<endl;
				}
				outlineridx[mtype1][mtype2]=true;
			}
			
			
			else if(rnorm>BPATOMPAIRLENGADJPARAM[aat1][aat2][paramidx][mtype1][mtype2].upper){
				length[mtype1][mtype2] = rnorm;
				diff[mtype1][mtype2] = length[mtype1][mtype2]-BPATOMPAIRLENGADJPARAM[aat1][aat2][paramidx][mtype1][mtype2].upper;
				mse += coeff*diff[mtype1][mtype2]*diff[mtype1][mtype2];
				if(showres){
					cout<<"----------"<<endl;
					cout<<rnorm<<endl;
					cout<<BPATOMPAIRLENGADJPARAM[aat1][aat2][paramidx][mtype1][mtype2].upper<<endl;
					cout<<diff[mtype1][mtype2]*diff[mtype1][mtype2]<<endl;
				}
				outlineridx[mtype1][mtype2]=true;
			}
			*/
			

			/*
			length[mtype1][mtype2] = rnorm;
			diff[mtype1][mtype2] = length[mtype1][mtype2]-BPATOMPAIRLENGADJPARAM[aat1][aat2][paramidx][mtype1][mtype2].mean;
			mse += coeff*diff[mtype1][mtype2]*diff[mtype1][mtype2];
			if(showres){
				cout<<"----------"<<endl;
				cout<<rnorm<<endl;
				cout<<BPATOMPAIRLENGADJPARAM[aat1][aat2][paramidx][mtype1][mtype2].upper<<endl;
				cout<<diff[mtype1][mtype2]*diff[mtype1][mtype2]<<endl;
			}
			outlineridx[mtype1][mtype2]=true;
			*/
		}
	}
}

void RNASimulation::calcmsebplength(Point3s *P, Point3s *Q, char aat1, char aat2, vector<vector<double>> &length, vector<vector<double>> &diff, vector<vector<Point3s>> &dists, double &mse, int &count, vector<vector<bool>> &outlineridx, double coeff, bool showres)
{
	int mtype1, mtype2;
	double rnorm;
	Point3s PQ;

	for(mtype1=PAIRLENGTHTARGETMOVNUM-3;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
	{
		for(mtype2=PAIRLENGTHTARGETMOVNUM-3;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
		{
			PQ.x = P[mtype1].x - Q[mtype2].x; 
			PQ.y = P[mtype1].y - Q[mtype2].y; 
			PQ.z = P[mtype1].z - Q[mtype2].z;

			dists[mtype1][mtype2].x=PQ.x;
			dists[mtype1][mtype2].y=PQ.y;
			dists[mtype1][mtype2].z=PQ.z;

			rnorm=bf.norm(PQ);

			/*
			if(rnorm<BPATOMPAIRLENGPARAM[aat1][aat2][mtype1][mtype2].lower){
				length[mtype1][mtype2] = rnorm;
				diff[mtype1][mtype2] = length[mtype1][mtype2]-BPATOMPAIRLENGPARAM[aat1][aat2][mtype1][mtype2].lower;
				mse += coeff*diff[mtype1][mtype2]*diff[mtype1][mtype2];
				if(showres){
					cout<<"----------"<<endl;
					cout<<rnorm<<endl;
					cout<<BPATOMPAIRLENGPARAM[aat1][aat2][mtype1][mtype2].lower<<endl;
					cout<<diff[mtype1][mtype2]*diff[mtype1][mtype2]<<endl;
				}
				outlineridx[mtype1][mtype2]=true;
			}else if(rnorm>BPATOMPAIRLENGPARAM[aat1][aat2][mtype1][mtype2].upper){
				length[mtype1][mtype2] = rnorm;
				diff[mtype1][mtype2] = length[mtype1][mtype2]-BPATOMPAIRLENGPARAM[aat1][aat2][mtype1][mtype2].upper;
				mse += coeff*diff[mtype1][mtype2]*diff[mtype1][mtype2];
				if(showres){
					cout<<"----------"<<endl;
					cout<<rnorm<<endl;
					cout<<BPATOMPAIRLENGPARAM[aat1][aat2][mtype1][mtype2].upper<<endl;
					cout<<diff[mtype1][mtype2]*diff[mtype1][mtype2]<<endl;
				}
				outlineridx[mtype1][mtype2]=true;
			}
			*/

			length[mtype1][mtype2] = rnorm;

			if(length[mtype1][mtype2] < BPATOMPAIRLENGPARAM[aat1][aat2][mtype1][mtype2].lower | length[mtype1][mtype2] > BPATOMPAIRLENGPARAM[aat1][aat2][mtype1][mtype2].upper){
				diff[mtype1][mtype2] = length[mtype1][mtype2]-BPATOMPAIRLENGPARAM[aat1][aat2][mtype1][mtype2].mean;
				mse += coeff*diff[mtype1][mtype2]*diff[mtype1][mtype2];
				if(showres){
					cout<<"----------"<<endl;
					cout<<rnorm<<endl;
					cout<<BPATOMPAIRLENGPARAM[aat1][aat2][mtype1][mtype2].lower<<endl;
					cout<<diff[mtype1][mtype2]*diff[mtype1][mtype2]<<endl;
				}
				outlineridx[mtype1][mtype2]=true;
			}
			count+=1;
		}
	}
}




void RNASimulation::calcmseinitlength(Point3s *P, Point3s *Q, vector<double> &length, vector<Point3s> &dists, double &mse, int &count, double coeff)
{
	int mtype1;
	double rnorm;
	Point3s PQ;

	for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
	//for(mtype1=0;mtype1<6;mtype1++)
	{
		PQ.x = P[mtype1].x - Q[mtype1].x; 
		PQ.y = P[mtype1].y - Q[mtype1].y; 
		PQ.z = P[mtype1].z - Q[mtype1].z;

		dists[mtype1].x=PQ.x;
		dists[mtype1].y=PQ.y;
		dists[mtype1].z=PQ.z;

		rnorm=bf.norm(PQ);

		length[mtype1] = rnorm;
		mse += coeff*rnorm;
		count+=1;
	}
}

void RNASimulation::calcmsetmplength(Point3s *P, Point3s *Q, vector<double> &length, vector<Point3s> &dists, double &mse, int &count, double coeff)
{
	int mtype1;
	double rnorm;
	Point3s PQ;

	for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
	//for(mtype1=0;mtype1<6;mtype1++)
	{
		PQ.x = P[mtype1].x - Q[mtype1].x; 
		PQ.y = P[mtype1].y - Q[mtype1].y; 
		PQ.z = P[mtype1].z - Q[mtype1].z;

		dists[mtype1].x=PQ.x;
		dists[mtype1].y=PQ.y;
		dists[mtype1].z=PQ.z;

		rnorm=bf.norm(PQ);

		length[mtype1] = rnorm;
		mse += coeff*rnorm;
		count+=1;
	}
}

void RNASimulation::calcgradientlength(Point3s *P, Point3s *Q, int basei, int basej, map<int, vector<Point3s>> &gradptv, vector<bool> movenodeinfo, vector<vector<double>> lengths, vector<vector<double>> diff, vector<vector<Point3s>> dists, int count, vector<vector<bool>> &outlineridx, double coeff)
{
	int mtype1, mtype2;
	double rnorm;
	Point3s PQ;

	for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
	{
		for(mtype2=0;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
		{
			if(basei == basej & mtype1 == mtype2) continue;

			if(movenodeinfo[basei] && outlineridx[mtype1][mtype2] && lengths[mtype1][mtype2] > 0){
				//gradptv[basei][mtype1].x += coeff * diff[mtype1][mtype2] * 2 * (dists[mtype1][mtype2].x)/(count*lengths[mtype1][mtype2]);
            	//gradptv[basei][mtype1].y += coeff * diff[mtype1][mtype2] * 2 * (dists[mtype1][mtype2].y)/(count*lengths[mtype1][mtype2]);
            	//gradptv[basei][mtype1].z += coeff * diff[mtype1][mtype2] * 2 * (dists[mtype1][mtype2].z)/(count*lengths[mtype1][mtype2]);
				
				gradptv[basei][mtype1].x += coeff * diff[mtype1][mtype2] * 2 * (dists[mtype1][mtype2].x)/(lengths[mtype1][mtype2]);
            	gradptv[basei][mtype1].y += coeff * diff[mtype1][mtype2] * 2 * (dists[mtype1][mtype2].y)/(lengths[mtype1][mtype2]);
            	gradptv[basei][mtype1].z += coeff * diff[mtype1][mtype2] * 2 * (dists[mtype1][mtype2].z)/(lengths[mtype1][mtype2]);
				

				//gradptv[basei][mtype1].x += coeff * diff[mtype1][mtype2] * 2 * (P[mtype1].x - Q[mtype2].x)/(count*lengths[mtype1][mtype2]);
            	//gradptv[basei][mtype1].y += coeff * diff[mtype1][mtype2] * 2 * (P[mtype1].y - Q[mtype2].y)/(count*lengths[mtype1][mtype2]);
            	//gradptv[basei][mtype1].z += coeff * diff[mtype1][mtype2] * 2 * (P[mtype1].z - Q[mtype2].z)/(count*lengths[mtype1][mtype2]);
			}

			if(movenodeinfo[basej] && outlineridx[mtype1][mtype2] && lengths[mtype1][mtype2] > 0){
				//cout<< coeff << "-" << diff[mtype1][mtype2] << "-" << P[mtype1].x << "-" << Q[mtype2].x << "-" << P[mtype1].x << "-" << count << "-" << lengths[mtype1][mtype2] << endl;
				//if(lengths[mtype1][mtype2] == 0){
					//cout<<"Error"<<endl;
					//cout<< P[mtype1].x << "-" << P[mtype1].y << "-" << P[mtype1].z<<endl;
					//cout<< Q[mtype2].x << "-" << Q[mtype2].y << "-" << Q[mtype2].z<<endl;
					//exit(1);
				//}
				//gradptv[basej][mtype2].x += coeff * diff[mtype1][mtype2] * 2 * (-dists[mtype1][mtype2].x)/(count*lengths[mtype1][mtype2]);
            	//gradptv[basej][mtype2].y += coeff * diff[mtype1][mtype2] * 2 * (-dists[mtype1][mtype2].y)/(count*lengths[mtype1][mtype2]);
            	//gradptv[basej][mtype2].z += coeff * diff[mtype1][mtype2] * 2 * (-dists[mtype1][mtype2].z)/(count*lengths[mtype1][mtype2]);

				gradptv[basej][mtype2].x += coeff * diff[mtype1][mtype2] * 2 * (-dists[mtype1][mtype2].x)/(lengths[mtype1][mtype2]);
            	gradptv[basej][mtype2].y += coeff * diff[mtype1][mtype2] * 2 * (-dists[mtype1][mtype2].y)/(lengths[mtype1][mtype2]);
            	gradptv[basej][mtype2].z += coeff * diff[mtype1][mtype2] * 2 * (-dists[mtype1][mtype2].z)/(lengths[mtype1][mtype2]);
				
				//cout<< gradptv[basej][mtype2].x << "-" << gradptv[basej][mtype2].y << "-" << gradptv[basej][mtype2].z<<endl;
			}
		}
	}
}

void RNASimulation::calcgradientinitlength(Point3s *P, Point3s *Q, int basei, map<int, vector<Point3s>> &gradptv, vector<double> lengths, vector<Point3s> dists, int count, double coeff)
{
	int mtype1;
	Point3s PQ;

	for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
	//for(mtype1=0;mtype1<6;mtype1++)
	{
		//gradptv[basei][mtype1].x += coeff * (dists[mtype1].x)/(count*lengths[mtype1]);
        //gradptv[basei][mtype1].y += coeff * (dists[mtype1].y)/(count*lengths[mtype1]);
        //gradptv[basei][mtype1].z += coeff * (dists[mtype1].z)/(count*lengths[mtype1]);

		gradptv[basei][mtype1].x += coeff * (dists[mtype1].x)/(lengths[mtype1]);
        gradptv[basei][mtype1].y += coeff * (dists[mtype1].y)/(lengths[mtype1]);
        gradptv[basei][mtype1].z += coeff * (dists[mtype1].z)/(lengths[mtype1]);

		//gradptv[basei][mtype1].x += coeff * (P[mtype1].x - Q[mtype1].x)/(count*lengths[mtype1]);
        //gradptv[basei][mtype1].y += coeff * (P[mtype1].y - Q[mtype1].y)/(count*lengths[mtype1]);
        //gradptv[basei][mtype1].z += coeff * (P[mtype1].z - Q[mtype1].z)/(count*lengths[mtype1]);
	}
}

void RNASimulation::calcgradienttmplength(Point3s *P, Point3s *Q, int basei, map<int, vector<Point3s>> &gradptv, vector<double> lengths, vector<Point3s> dists, int count, double coeff)
{
	int mtype1;
	Point3s PQ;

	for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
	//for(mtype1=0;mtype1<6;mtype1++)
	{
		gradptv[basei][mtype1].x += coeff * (dists[mtype1].x)/(count*lengths[mtype1]);
        gradptv[basei][mtype1].y += coeff * (dists[mtype1].y)/(count*lengths[mtype1]);
        gradptv[basei][mtype1].z += coeff * (dists[mtype1].z)/(count*lengths[mtype1]);

		//gradptv[basei][mtype1].x += coeff * (P[mtype1].x - Q[mtype1].x)/(count*lengths[mtype1]);
        //gradptv[basei][mtype1].y += coeff * (P[mtype1].y - Q[mtype1].y)/(count*lengths[mtype1]);
        //gradptv[basei][mtype1].z += coeff * (P[mtype1].z - Q[mtype1].z)/(count*lengths[mtype1]);
	}
}

/*
void RNASimulation::calcmsetempadj(Point3s *restr1, Point3s *restr2, char aat1, char aat2, vector<Residueatoms> strtemplate, vector<double> &length, vector<Point3s> &dists, double &mse, int &count, double coeff){
	int mtype1, mtype2;
	double rnorm;
	Point3s PQ;

	int i, j, k;
	vector<vector<double>> RotMatix;
    vector<double> TranVect;
	vector<double> xyztmp(3,0);
	double minrmsd=100000000000;
	double rmsd;

	int parammun = AccAtomParamNum[aat1][2] + AccAtomParamNum[aat2][2];
	vector<double> tmp(3,0);
	vector<vector<double> > xyz_list1(parammun,tmp);
    vector<vector<double> > xyz_list2(parammun,tmp);
	vector<Point3s> diststmp(parammun);
	vector<double> lengthtmp(parammun, 0.0);

	//cout<<"bbb1"<<endl;

	for(i=0;i<DILENTTEMPLATE[aat1][aat2].size();i++){
		k = 0;
		for(j=0; j<AccAtomParamNum[aat1][2];j++){
			xyz_list1[k][0]=restr1[j].x;
			xyz_list1[k][1]=restr1[j].y;
			xyz_list1[k][2]=restr1[j].z;

			xyz_list2[k][0]=DILENTTEMPLATE[aat1][aat2][i][0].ptv[j].x;
			xyz_list2[k][1]=DILENTTEMPLATE[aat1][aat2][i][0].ptv[j].y;
			xyz_list2[k][2]=DILENTTEMPLATE[aat1][aat2][i][0].ptv[j].z;
			k+=1;
		}

		for(j=0; j<AccAtomParamNum[aat2][2];j++){
			xyz_list1[k][0]=restr2[j].x;
			xyz_list1[k][1]=restr2[j].y;
			xyz_list1[k][2]=restr2[j].z;

			xyz_list2[k][0]=DILENTTEMPLATE[aat1][aat2][i][1].ptv[j].x;
			xyz_list2[k][1]=DILENTTEMPLATE[aat1][aat2][i][1].ptv[j].y;
			xyz_list2[k][2]=DILENTTEMPLATE[aat1][aat2][i][1].ptv[j].z;
			k+=1;
		}

		//cout<<"bbb2"<<endl;

		bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

		//cout<<"bbb3"<<endl;

		rmsd=0;
		k=0;
		for(j=0; j<AccAtomParamNum[aat1][2];j++){
			bf.Transform(xyz_list2[k],RotMatix,TranVect,xyztmp);
			diststmp[k].x=restr1[j].x-xyztmp[0];
			diststmp[k].y=restr1[j].y-xyztmp[1];
			diststmp[k].z=restr1[j].z-xyztmp[2];
			lengthtmp[k]=bf.norm(diststmp[k]);
			rmsd+=lengthtmp[k];
			k+=1;
		}
		for(j=0; j<AccAtomParamNum[aat2][2];j++){
			bf.Transform(xyz_list2[k],RotMatix,TranVect,xyztmp);
			diststmp[k].x=restr2[j].x-xyztmp[0];
			diststmp[k].y=restr2[j].y-xyztmp[1];
			diststmp[k].z=restr2[j].z-xyztmp[2];
			lengthtmp[k]=bf.norm(diststmp[k]);
			rmsd+=lengthtmp[k];
			k+=1;
		}

		//cout<<"bbb4"<<endl;

		if(minrmsd > rmsd){
			minrmsd=rmsd;

			k=0;
			for(j=0;j<AccAtomParamNum[aat1][2];j++){
				dists[j].x=diststmp[k].x;
				dists[j].y=diststmp[k].y;
				dists[j].z=diststmp[k].z;
				length[j]=lengthtmp[k];
				k+=1;
			}
			for(j=0;j<AccAtomParamNum[aat2][2];j++){
				dists[AccAtomParamNum[aat1][2]+j].x=diststmp[k].x;
				dists[AccAtomParamNum[aat1][2]+j].y=diststmp[k].y;
				dists[AccAtomParamNum[aat1][2]+j].z=diststmp[k].z;
				length[j]=lengthtmp[k];
				k+=1;
			}
		}

		//cout<<"bbb5"<<endl;
	}
	count=parammun;
	mse=minrmsd;
}
*/


void RNASimulation::calcmsetempadj(Point3s *restr1, Point3s *restr2, char aat1, char aat2, vector<Residueatoms> strtemplate, vector<double> &length, vector<Point3s> &dists, double &mse, int &count, double coeff){
	int mtype1, mtype2;
	double rnorm;
	Point3s PQ;

	int i, j, k;
	vector<vector<double>> RotMatix;
    vector<double> TranVect;
	vector<double> xyztmp(3,0);
	//double minrmsd=100000000000;
	double rmsd;

	int parammun = AccAtomParamNum[aat1][2] + AccAtomParamNum[aat2][2];
	vector<double> tmp(3,0);
	vector<vector<double> > xyz_list1(parammun,tmp);
    vector<vector<double> > xyz_list2(parammun,tmp);
	//vector<Point3s> diststmp(parammun);
	//vector<double> lengthtmp(parammun, 0.0);

	//cout<<"bbb1"<<endl;
	
	//for(i=0;i<DILENTTEMPLATE[aat1][aat2].size();i++){
	k = 0;
	for(j=0; j<AccAtomParamNum[aat1][2];j++){
		xyz_list1[k][0]=restr1[j].x;
		xyz_list1[k][1]=restr1[j].y;
		xyz_list1[k][2]=restr1[j].z;

		xyz_list2[k][0]=strtemplate[0].ptv[j].x;
		xyz_list2[k][1]=strtemplate[0].ptv[j].y;
		xyz_list2[k][2]=strtemplate[0].ptv[j].z;

		//cout<<"==========="<<endl;
		//cout<<xyz_list1[k][0]<< "-" << xyz_list1[k][1] << "-" << xyz_list1[k][2]<<endl;
		//cout<<xyz_list2[k][0]<< "-" << xyz_list2[k][1] << "-" << xyz_list2[k][2]<<endl;
		k+=1;
	}

	//cout<<"--------------"<<endl;

	for(j=0; j<AccAtomParamNum[aat2][2];j++){
		xyz_list1[k][0]=restr2[j].x;
		xyz_list1[k][1]=restr2[j].y;
		xyz_list1[k][2]=restr2[j].z;

		xyz_list2[k][0]=strtemplate[1].ptv[j].x;
		xyz_list2[k][1]=strtemplate[1].ptv[j].y;
		xyz_list2[k][2]=strtemplate[1].ptv[j].z;

		//cout<<"==========="<<endl;
		//cout<<xyz_list1[k][0]<< "-" << xyz_list1[k][1] << "-" << xyz_list1[k][2]<<endl;
		//cout<<xyz_list2[k][0]<< "-" << xyz_list2[k][1] << "-" << xyz_list2[k][2]<<endl;
		k+=1;
	}

	//cout<<"bbb2"<<endl;

	bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

	//cout<<"bbb3"<<endl;

	rmsd=0;
	k=0;
	for(j=0; j<AccAtomParamNum[aat1][2];j++){
		bf.Transform(xyz_list2[k],RotMatix,TranVect,xyztmp);
		dists[j].x=restr1[j].x-xyztmp[0];
		dists[j].y=restr1[j].y-xyztmp[1];
		dists[j].z=restr1[j].z-xyztmp[2];
		length[j]=bf.norm(dists[j]);
		rmsd+=length[j];
		k+=1;
	}
	for(j=0; j<AccAtomParamNum[aat2][2];j++){
		bf.Transform(xyz_list2[k],RotMatix,TranVect,xyztmp);
		dists[AccAtomParamNum[aat1][2]+j].x=restr2[j].x-xyztmp[0];
		dists[AccAtomParamNum[aat1][2]+j].y=restr2[j].y-xyztmp[1];
		dists[AccAtomParamNum[aat1][2]+j].z=restr2[j].z-xyztmp[2];
		length[AccAtomParamNum[aat1][2]+j]=bf.norm(dists[AccAtomParamNum[aat1][2]+j]);
		rmsd+=length[AccAtomParamNum[aat1][2]+j];
		k+=1;
	}
	
	/*
	rmsd=0;
	k=0;
	for(j=0; j<AccAtomParamNum[aat1][2];j++){
		//bf.Transform(xyz_list2[k],RotMatix,TranVect,xyztmp);
		dists[j].x=restr1[j].x-strtemplate[0].ptv[j].x;
		dists[j].y=restr1[j].y-strtemplate[0].ptv[j].y;
		dists[j].z=restr1[j].z-strtemplate[0].ptv[j].z;
		length[j]=bf.norm(dists[j]);
		rmsd+=length[j];
		k+=1;
	}
	for(j=0; j<AccAtomParamNum[aat2][2];j++){
		//bf.Transform(xyz_list2[k],RotMatix,TranVect,xyztmp);
		dists[AccAtomParamNum[aat1][2]+j].x=restr2[j].x-strtemplate[1].ptv[j].x;
		dists[AccAtomParamNum[aat1][2]+j].y=restr2[j].y-strtemplate[1].ptv[j].y;
		dists[AccAtomParamNum[aat1][2]+j].z=restr2[j].z-strtemplate[1].ptv[j].z;
		length[AccAtomParamNum[aat1][2]+j]=bf.norm(dists[AccAtomParamNum[aat1][2]+j]);
		rmsd+=length[AccAtomParamNum[aat1][2]+j];
		k+=1;
	}
	*/

	count=parammun;
	mse=coeff*rmsd;
}

void RNASimulation::calcmsetempbp(Point3s *restr1, Point3s *restr2, char aat1, char aat2, vector<Residueatoms> strtemplate, vector<double> &length, vector<Point3s> &dists, double &mse, int &count, double coeff){
	int mtype1, mtype2;
	double rnorm;
	Point3s PQ;

	int i, j, k;
	vector<vector<double>> RotMatix;
    vector<double> TranVect;
	vector<double> xyztmp(3,0);
	//double minrmsd=100000000000;
	double rmsd;
	int atom;

	
	map<char, vector<int>> SuppTargetAtoms={
		{'A',{6,7,8,9,12,13,14}},
		{'C',{6,7,8,9,12,13,14}},
		{'G',{6,7,8,9,12,13,14}},
		{'U',{6,7,8,9,12,13,14}}
	};

	//int parammun = AccAtomParamNum[aat1][2] + AccAtomParamNum[aat2][2] - AccAtomParamNum[aat1][1] - AccAtomParamNum[aat2][1];
	int parammun = SuppTargetAtoms[aat1].size() + SuppTargetAtoms[aat2].size();
	vector<double> tmp(3,0);
	vector<vector<double> > xyz_list1(parammun,tmp);
    vector<vector<double> > xyz_list2(parammun,tmp);
	//vector<Point3s> diststmp(parammun);
	//vector<double> lengthtmp(parammun, 0.0);

	//cout<<"bbb1"<<endl;

	//for(i=0;i<DILENTTEMPLATE[aat1][aat2].size();i++){
	
	k = 0;
	//for(j=AccAtomParamNum[aat1][1]; j<AccAtomParamNum[aat1][2];j++){
	for(j=0; j<SuppTargetAtoms[aat1].size();j++){
		atom=SuppTargetAtoms[aat1][j];
		xyz_list1[k][0]=restr1[atom].x;
		xyz_list1[k][1]=restr1[atom].y;
		xyz_list1[k][2]=restr1[atom].z;

		xyz_list2[k][0]=strtemplate[0].ptv[atom].x;
		xyz_list2[k][1]=strtemplate[0].ptv[atom].y;
		xyz_list2[k][2]=strtemplate[0].ptv[atom].z;
		k+=1;
	}

	//for(j=AccAtomParamNum[aat1][1]; j<AccAtomParamNum[aat2][2];j++){
	for(j=0; j<SuppTargetAtoms[aat2].size();j++){
		atom=SuppTargetAtoms[aat2][j];
		xyz_list1[k][0]=restr2[atom].x;
		xyz_list1[k][1]=restr2[atom].y;
		xyz_list1[k][2]=restr2[atom].z;

		xyz_list2[k][0]=strtemplate[1].ptv[atom].x;
		xyz_list2[k][1]=strtemplate[1].ptv[atom].y;
		xyz_list2[k][2]=strtemplate[1].ptv[atom].z;
		k+=1;
	}

	//cout<<"bbb2"<<endl;

	bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

	//cout<<"bbb3"<<endl;

	rmsd=0;
	k=0;
	//for(j=AccAtomParamNum[aat1][1]; j<AccAtomParamNum[aat1][2];j++){
	for(j=0; j<SuppTargetAtoms[aat1].size();j++){
		atom=SuppTargetAtoms[aat1][j];
		bf.Transform(xyz_list2[k],RotMatix,TranVect,xyztmp);
		dists[j].x=restr1[atom].x-xyztmp[0];
		dists[j].y=restr1[atom].y-xyztmp[1];
		dists[j].z=restr1[atom].z-xyztmp[2];
		length[j]=bf.norm(dists[j]);
		rmsd+=length[j];
		k+=1;
	}
	//for(j=AccAtomParamNum[aat1][1]; j<AccAtomParamNum[aat2][2];j++){
	for(j=0; j<SuppTargetAtoms[aat2].size();j++){
		atom=SuppTargetAtoms[aat2][j];
		bf.Transform(xyz_list2[k],RotMatix,TranVect,xyztmp);
		dists[SuppTargetAtoms[aat1].size()+j].x=restr2[atom].x-xyztmp[0];
		dists[SuppTargetAtoms[aat1].size()+j].y=restr2[atom].y-xyztmp[1];
		dists[SuppTargetAtoms[aat1].size()+j].z=restr2[atom].z-xyztmp[2];
		length[SuppTargetAtoms[aat1].size()+j]=bf.norm(dists[SuppTargetAtoms[aat1].size()+j]);
		rmsd+=length[SuppTargetAtoms[aat1].size()+j];
		k+=1;
	}
	
	/*
	rmsd=0;
	k=0;
	//for(j=AccAtomParamNum[aat1][1]; j<AccAtomParamNum[aat1][2];j++){
	for(j=0; j<SuppTargetAtoms[aat1].size();j++){
		atom=SuppTargetAtoms[aat1][j];
		//bf.Transform(xyz_list2[k],RotMatix,TranVect,xyztmp);
		dists[j].x=restr1[atom].x-strtemplate[0].ptv[atom].x;
		dists[j].y=restr1[atom].y-strtemplate[0].ptv[atom].y;
		dists[j].z=restr1[atom].z-strtemplate[0].ptv[atom].z;
		length[j]=bf.norm(dists[j]);
		rmsd+=length[j];
		k+=1;
	}
	//for(j=AccAtomParamNum[aat1][1]; j<AccAtomParamNum[aat2][2];j++){
	for(j=0; j<SuppTargetAtoms[aat2].size();j++){
		atom=SuppTargetAtoms[aat2][j];
		//bf.Transform(xyz_list2[k],RotMatix,TranVect,xyztmp);
		dists[SuppTargetAtoms[aat1].size()+j].x=restr2[atom].x-strtemplate[1].ptv[atom].x;
		dists[SuppTargetAtoms[aat1].size()+j].y=restr2[atom].y-strtemplate[1].ptv[atom].y;
		dists[SuppTargetAtoms[aat1].size()+j].z=restr2[atom].z-strtemplate[1].ptv[atom].z;
		length[SuppTargetAtoms[aat1].size()+j]=bf.norm(dists[SuppTargetAtoms[aat1].size()+j]);
		rmsd+=length[SuppTargetAtoms[aat1].size()+j];
		k+=1;
	}
	*/
		//cout<<"bbb5"<<endl;
	//}
	count=parammun;
	mse=coeff*rmsd;
}


//void RNASimulation::calcmseditmplength(Point3s *P1, Point3s *Q1, Point3s *P2, Point3s *Q2, vector<vector<double>> &length, vector<vector<Point3s>> &dists, double &mse, int &count, double coeff)
void RNASimulation::calcmseditmplength(Point3s *P1, Point3s *Q1, Point3s *P2, Point3s *Q2, TemplateDistMap templateslength, vector<vector<double>> &length, vector<vector<double>> &diff, vector<vector<Point3s>> &dists, double &mse, int &count, double coeff)
{
	int mtype1, mtype2;
	double rnorm;
	Point3s PQ;

	for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
	{
		for(mtype2=0;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
		{
			PQ.x = P1[mtype1].x - P2[mtype2].x; 
			PQ.y = P1[mtype1].y - P2[mtype2].y; 
			PQ.z = P1[mtype1].z - P2[mtype2].z;

			dists[mtype1][mtype2].x=PQ.x;
			dists[mtype1][mtype2].y=PQ.y;
			dists[mtype1][mtype2].z=PQ.z;

			rnorm=bf.norm(PQ);

			length[mtype1][mtype2] = rnorm;
			
			diff[mtype1][mtype2] = length[mtype1][mtype2]-templateslength.interdistmap[mtype1][mtype2];

			mse += coeff*diff[mtype1][mtype2]*diff[mtype1][mtype2];
			count+=1;
		}
	}

	/*
	for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
	//for(mtype1=0;mtype1<6;mtype1++)
	{
		PQ.x = P1[mtype1].x - Q1[mtype1].x; 
		PQ.y = P1[mtype1].y - Q1[mtype1].y; 
		PQ.z = P1[mtype1].z - Q1[mtype1].z;

		dists[mtype1].x=PQ.x;
		dists[mtype1].y=PQ.y;
		dists[mtype1].z=PQ.z;

		rnorm=bf.norm(PQ);

		length[mtype1] = rnorm;
		mse += coeff*rnorm;
		count+=1;
	}

	for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
	//for(mtype1=0;mtype1<6;mtype1++)
	{
		PQ.x = P2[mtype1].x - Q2[mtype1].x; 
		PQ.y = P2[mtype1].y - Q2[mtype1].y; 
		PQ.z = P2[mtype1].z - Q2[mtype1].z;

		dists[PAIRLENGTHTARGETMOVNUM+mtype1].x=PQ.x;
		dists[PAIRLENGTHTARGETMOVNUM+mtype1].y=PQ.y;
		dists[PAIRLENGTHTARGETMOVNUM+mtype1].z=PQ.z;

		rnorm=bf.norm(PQ);

		length[PAIRLENGTHTARGETMOVNUM+mtype1] = rnorm;
		mse += coeff*rnorm;
		count+=1;
	}
	*/
}

void RNASimulation::calcgradienttempadj(int basei, int basej, char aat1, char aat2, map<int, vector<Point3s>> &gradptv, vector<double> lengths, vector<Point3s> dists, int count, double coeff)
{
	int i;
	int atom;
	for(i=0;i<AccAtomParamNum[aat1][2];i++)
	{
		if(lengths[i] > 0){
			//gradptv[basei][i].x += coeff * (dists[i].x)/(count*lengths[i]);
			//gradptv[basei][i].y += coeff * (dists[i].y)/(count*lengths[i]);
			//gradptv[basei][i].z += coeff * (dists[i].z)/(count*lengths[i]);
			gradptv[basei][i].x += coeff * (dists[i].x)/(lengths[i]);
			gradptv[basei][i].y += coeff * (dists[i].y)/(lengths[i]);
			gradptv[basei][i].z += coeff * (dists[i].z)/(lengths[i]);
		}
	}

	for(i=0;i<AccAtomParamNum[aat2][2];i++)
	{
		if(lengths[AccAtomParamNum[aat1][2]+i] > 0){
			//gradptv[basej][i].x += coeff * (dists[AccAtomParamNum[aat1][2]+i].x)/(count*lengths[AccAtomParamNum[aat1][2]+i]);
			//gradptv[basej][i].y += coeff * (dists[AccAtomParamNum[aat1][2]+i].y)/(count*lengths[AccAtomParamNum[aat1][2]+i]);
			//gradptv[basej][i].z += coeff * (dists[AccAtomParamNum[aat1][2]+i].z)/(count*lengths[AccAtomParamNum[aat1][2]+i]);
			gradptv[basej][i].x += coeff * (dists[AccAtomParamNum[aat1][2]+i].x)/(lengths[AccAtomParamNum[aat1][2]+i]);
			gradptv[basej][i].y += coeff * (dists[AccAtomParamNum[aat1][2]+i].y)/(lengths[AccAtomParamNum[aat1][2]+i]);
			gradptv[basej][i].z += coeff * (dists[AccAtomParamNum[aat1][2]+i].z)/(lengths[AccAtomParamNum[aat1][2]+i]);
		}
	}
}

void RNASimulation::calcgradienttempbp(int basei, int basej, char aat1, char aat2, map<int, vector<Point3s>> &gradptv, vector<double> lengths, vector<Point3s> dists, int count, double coeff)
{
	int i;
	int atom;

	map<char, vector<int>> SuppTargetAtoms={
		{'A',{6,7,8,9,12,13,14}},
		{'C',{6,7,8,9,12,13,14}},
		{'G',{6,7,8,9,12,13,14}},
		{'U',{6,7,8,9,12,13,14}}
	};

	//for(i=AccAtomParamNum[aat1][1];i<AccAtomParamNum[aat1][2];i++)
	for(i=0;i<SuppTargetAtoms[aat1].size();i++)
	{
		atom=SuppTargetAtoms[aat1][i];
		if(lengths[i] > 0){
			//gradptv[basei][i].x += coeff * (dists[i].x)/(count*lengths[i]);
			//gradptv[basei][i].y += coeff * (dists[i].y)/(count*lengths[i]);
			//gradptv[basei][i].z += coeff * (dists[i].z)/(count*lengths[i]);
			gradptv[basei][atom].x += coeff * (dists[i].x)/(lengths[i]);
			gradptv[basei][atom].y += coeff * (dists[i].y)/(lengths[i]);
			gradptv[basei][atom].z += coeff * (dists[i].z)/(lengths[i]);
		}
	}

	//for(i=AccAtomParamNum[aat2][1];i<AccAtomParamNum[aat2][2];i++)
	for(i=0;i<SuppTargetAtoms[aat2].size();i++)
	{
		atom=SuppTargetAtoms[aat2][i];
		if(lengths[SuppTargetAtoms[aat1].size()+i] > 0){
			//gradptv[basej][i].x += coeff * (dists[AccAtomParamNum[aat1][2]+i].x)/(count*lengths[AccAtomParamNum[aat1][2]+i]);
			//gradptv[basej][i].y += coeff * (dists[AccAtomParamNum[aat1][2]+i].y)/(count*lengths[AccAtomParamNum[aat1][2]+i]);
			//gradptv[basej][i].z += coeff * (dists[AccAtomParamNum[aat1][2]+i].z)/(count*lengths[AccAtomParamNum[aat1][2]+i]);
			gradptv[basej][atom].x += coeff * (dists[SuppTargetAtoms[aat1].size()+i].x)/(lengths[SuppTargetAtoms[aat1].size()+i]);
			gradptv[basej][atom].y += coeff * (dists[SuppTargetAtoms[aat1].size()+i].y)/(lengths[SuppTargetAtoms[aat1].size()+i]);
			gradptv[basej][atom].z += coeff * (dists[SuppTargetAtoms[aat1].size()+i].z)/(lengths[SuppTargetAtoms[aat1].size()+i]);
		}
	}
}

//void RNASimulation::calcgradientditmplength(int basei, int basej, map<int, vector<Point3s>> &gradptv, vector<double> lengths, vector<Point3s> dists, int count, double coeff)
void RNASimulation::calcgradientditmplength(int basei, int basej, map<int, vector<Point3s>> &gradptv, vector<bool> movenodeinfo, vector<vector<double>> lengths, vector<vector<double>> diff, vector<vector<Point3s>> dists, int count, double coeff)
{
	int mtype1, mtype2;
	Point3s PQ;

	
	for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
	{
		for(mtype2=0;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
		{
			//if(basei == basej & mtype1 == mtype2) continue;

			if(lengths[mtype1][mtype2] > 0){
				//gradptv[basei][mtype1].x += coeff * diff[mtype1][mtype2] * 2 * (dists[mtype1][mtype2].x)/(count*lengths[mtype1][mtype2]);
				//gradptv[basei][mtype1].y += coeff * diff[mtype1][mtype2] * 2 * (dists[mtype1][mtype2].y)/(count*lengths[mtype1][mtype2]);
				//gradptv[basei][mtype1].z += coeff * diff[mtype1][mtype2] * 2 * (dists[mtype1][mtype2].z)/(count*lengths[mtype1][mtype2]);

				gradptv[basei][mtype1].x += coeff * diff[mtype1][mtype2] * 2 * (dists[mtype1][mtype2].x)/(lengths[mtype1][mtype2]);
				gradptv[basei][mtype1].y += coeff * diff[mtype1][mtype2] * 2 * (dists[mtype1][mtype2].y)/(lengths[mtype1][mtype2]);
				gradptv[basei][mtype1].z += coeff * diff[mtype1][mtype2] * 2 * (dists[mtype1][mtype2].z)/(lengths[mtype1][mtype2]);

				//cout<<"g-" << mtype1 << "-" << mtype2 <<  "-" << gradptv[basei][mtype1].x <<  "-" << diff[mtype1][mtype2]  <<  "-" << dists[mtype1][mtype2].x  <<  "-" << lengths[mtype1][mtype2] <<  "-" << count<<endl;
				//gradptv[basei][mtype1].x += coeff * diff[mtype1][mtype2] * 2 * (P[mtype1].x - Q[mtype2].x)/(count*lengths[mtype1][mtype2]);
				//gradptv[basei][mtype1].y += coeff * diff[mtype1][mtype2] * 2 * (P[mtype1].y - Q[mtype2].y)/(count*lengths[mtype1][mtype2]);
				//gradptv[basei][mtype1].z += coeff * diff[mtype1][mtype2] * 2 * (P[mtype1].z - Q[mtype2].z)/(count*lengths[mtype1][mtype2]);
			}

			if(lengths[mtype1][mtype2] > 0){
				//cout<< coeff << "-" << diff[mtype1][mtype2] << "-" << P[mtype1].x << "-" << Q[mtype2].x << "-" << P[mtype1].x << "-" << count << "-" << lengths[mtype1][mtype2] << endl;
				//if(lengths[mtype1][mtype2] == 0){
					//cout<<"Error"<<endl;
					//cout<< P[mtype1].x << "-" << P[mtype1].y << "-" << P[mtype1].z<<endl;
					//cout<< Q[mtype2].x << "-" << Q[mtype2].y << "-" << Q[mtype2].z<<endl;
					//exit(1);
				//}
				//gradptv[basej][mtype2].x += coeff * diff[mtype1][mtype2] * 2 * (-dists[mtype1][mtype2].x)/(count*lengths[mtype1][mtype2]);
				//gradptv[basej][mtype2].y += coeff * diff[mtype1][mtype2] * 2 * (-dists[mtype1][mtype2].y)/(count*lengths[mtype1][mtype2]);
				//gradptv[basej][mtype2].z += coeff * diff[mtype1][mtype2] * 2 * (-dists[mtype1][mtype2].z)/(count*lengths[mtype1][mtype2]);
				
				gradptv[basej][mtype2].x += coeff * diff[mtype1][mtype2] * 2 * (-dists[mtype1][mtype2].x)/(lengths[mtype1][mtype2]);
				gradptv[basej][mtype2].y += coeff * diff[mtype1][mtype2] * 2 * (-dists[mtype1][mtype2].y)/(lengths[mtype1][mtype2]);
				gradptv[basej][mtype2].z += coeff * diff[mtype1][mtype2] * 2 * (-dists[mtype1][mtype2].z)/(lengths[mtype1][mtype2]);
				
				//cout<< gradptv[basej][mtype2].x << "-" << gradptv[basej][mtype2].y << "-" << gradptv[basej][mtype2].z<<endl;
			}
		}
	}

	/*
	for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
	//for(mtype1=0;mtype1<6;mtype1++)
	{
		gradptv[basei][mtype1].x += coeff * (dists[mtype1].x)/(count*lengths[mtype1]);
        gradptv[basei][mtype1].y += coeff * (dists[mtype1].y)/(count*lengths[mtype1]);
        gradptv[basei][mtype1].z += coeff * (dists[mtype1].z)/(count*lengths[mtype1]);

		//cout<< gradptv[basei][mtype1].x << "-" << gradptv[basei][mtype1].y << "-" << gradptv[basei][mtype1].z<<endl;
	}

	for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
	//for(mtype1=0;mtype1<6;mtype1++)
	{
		gradptv[basej][mtype1].x += coeff * (dists[PAIRLENGTHTARGETMOVNUM+mtype1].x)/(count*lengths[PAIRLENGTHTARGETMOVNUM+mtype1]);
        gradptv[basej][mtype1].y += coeff * (dists[PAIRLENGTHTARGETMOVNUM+mtype1].y)/(count*lengths[PAIRLENGTHTARGETMOVNUM+mtype1]);
        gradptv[basej][mtype1].z += coeff * (dists[PAIRLENGTHTARGETMOVNUM+mtype1].z)/(count*lengths[PAIRLENGTHTARGETMOVNUM+mtype1]);
	
		//cout<< gradptv[basej][mtype1].x << "-" << gradptv[basej][mtype1].y << "-" << gradptv[basej][mtype1].z<<endl;
	}
	*/
}

void RNASimulation::calcmsecnnditmpintralength(vector<Residueatoms> restr, int base, vector<vector<double>> intralength, vector<vector<double>> interlength, vector<vector<double>> &length, vector<vector<double>> &diff, vector<vector<Point3s>> &dists, double &mse, int &count, double coeff)
{
	int i;
	int mtype1, mtype2;
	double rnorm;
	int basei, basej;
	char aat1, aat2, aat;
	Point3s PQ;
	double upperl, lowerl;

	aat=restr[base].aaa;

	for(i=0;i<AccAtomParamNum[aat][2];i++){
		if(i==0 | i==1 | i==11){
			if(base == 0) continue;
			basei=base-1; basej=base;

			aat1=restr[basei].aaa;
			aat2=restr[basej].aaa;

			mtype1=ParamsList[aat][i][1];
			mtype2=ParamsList[aat][i][3];

			PQ.x = restr[basei].ptv[mtype1].x - restr[basej].ptv[mtype2].x;
			PQ.y = restr[basei].ptv[mtype1].y - restr[basej].ptv[mtype2].y;
			PQ.z = restr[basei].ptv[mtype1].z - restr[basej].ptv[mtype2].z;

			dists[mtype1][mtype2].x=PQ.x;
			dists[mtype1][mtype2].y=PQ.y;
			dists[mtype1][mtype2].z=PQ.z;

			rnorm=bf.norm(PQ);

			length[mtype1][mtype2] = rnorm;

			//cout<< i << "-" << mtype1 << "-" << mtype2 << "-" << length[mtype1][mtype2] << "-" << interlength[mtype2][mtype1]<<endl;

			//if(length[mtype1][mtype2] < LENGTHPARAM[aat][BLENGTHTARGET[i]].lower | length[mtype1][mtype2] > LENGTHPARAM[aat][BLENGTHTARGET[i]].upper){
			diff[mtype1][mtype2] = length[mtype1][mtype2]-interlength[mtype2][mtype1];
			mse += coeff*diff[mtype1][mtype2]*diff[mtype1][mtype2];
			//outlineridx[mtype1][mtype2]=true;
			//}
			count+=1;

		}else{
			basei=base; basej=base;

			mtype1=ParamsList[aat][i][1];
			mtype2=ParamsList[aat][i][3];

			PQ.x = restr[basei].ptv[mtype1].x - restr[basej].ptv[mtype2].x;
			PQ.y = restr[basei].ptv[mtype1].y - restr[basej].ptv[mtype2].y;
			PQ.z = restr[basei].ptv[mtype1].z - restr[basej].ptv[mtype2].z;

			dists[mtype1][mtype2].x=PQ.x;
			dists[mtype1][mtype2].y=PQ.y;
			dists[mtype1][mtype2].z=PQ.z;

			rnorm=bf.norm(PQ);
			length[mtype1][mtype2] = rnorm;

			//cout<< i << "-" << mtype1 << "-" << mtype2 << "-" << length[mtype1][mtype2] << "-" << intralength[mtype1][mtype2]<<endl;

			diff[mtype1][mtype2] = length[mtype1][mtype2]-intralength[mtype1][mtype2];
			mse += coeff*diff[mtype1][mtype2]*diff[mtype1][mtype2];

			//if(length[mtype1][mtype2] < LENGTHPARAM[aat][BLENGTHTARGET[i]].lower | length[mtype1][mtype2] > LENGTHPARAM[aat][BLENGTHTARGET[i]].upper){
				//diff[mtype1][mtype2] = length[mtype1][mtype2]-LENGTHPARAM[aat][BLENGTHTARGET[i]].mean;
				//mse += coeff*diff[mtype1][mtype2]*diff[mtype1][mtype2];
				//outlineridx[mtype1][mtype2]=true;
			//}
			count+=1;
		}
	}


	/*
	for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
	{
		for(mtype2=mtype1+1;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
		{
			PQ.x = P[mtype1].x - P[mtype2].x; 
			PQ.y = P[mtype1].y - P[mtype2].y; 
			PQ.z = P[mtype1].z - P[mtype2].z;

			dists[mtype1][mtype2].x=PQ.x;
			dists[mtype1][mtype2].y=PQ.y;
			dists[mtype1][mtype2].z=PQ.z;

			rnorm=bf.norm(PQ);

			length[mtype1][mtype2] = rnorm;
			
			diff[mtype1][mtype2] = length[mtype1][mtype2]-templatesintralength[mtype1][mtype2];

			mse += coeff*diff[mtype1][mtype2]*diff[mtype1][mtype2];
			count+=1;
		}
	}
	*/
}

void RNASimulation::calcmseditmpintralength(Point3s *P, vector<vector<double>> templatesintralength, vector<vector<double>> &length, vector<vector<double>> &diff, vector<vector<Point3s>> &dists, double &mse, int &count, double coeff)
{
	int mtype1, mtype2;
	double rnorm;
	Point3s PQ;

	for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
	{
		for(mtype2=mtype1+1;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
		{
			PQ.x = P[mtype1].x - P[mtype2].x; 
			PQ.y = P[mtype1].y - P[mtype2].y; 
			PQ.z = P[mtype1].z - P[mtype2].z;

			dists[mtype1][mtype2].x=PQ.x;
			dists[mtype1][mtype2].y=PQ.y;
			dists[mtype1][mtype2].z=PQ.z;

			rnorm=bf.norm(PQ);

			length[mtype1][mtype2] = rnorm;
			
			diff[mtype1][mtype2] = length[mtype1][mtype2]-templatesintralength[mtype1][mtype2];

			mse += coeff*diff[mtype1][mtype2]*diff[mtype1][mtype2];
			count+=1;
		}
	}
}


void RNASimulation::calcgracnnentditmpintra(int base, map<int, vector<Point3s>> &gradptv, vector<bool> movenodeinfo, vector<vector<double>> lengths, vector<vector<double>> diff, vector<vector<Point3s>> dists, int count, double coeff)
{
	int i;
	int mtype1, mtype2;
	Point3s PQ;
	int basei, basej;
	char aat1, aat2, aat;

	aat=restr[base].aaa;

	for(i=0;i<AccAtomParamNum[aat][2];i++){

		if(i==0 | i==1 | i==11){
			if(base == 0) continue;
			basei=base-1; basej=base;

			aat1=restr[basei].aaa;
			aat2=restr[basej].aaa;

			mtype1=ParamsList[aat][i][1];
			mtype2=ParamsList[aat][i][3];
		}else{
			basei=base; basej=base;

			mtype1=ParamsList[aat][i][1];
			mtype2=ParamsList[aat][i][3];
		}
	
		if(lengths[mtype1][mtype2] > 0){
			gradptv[basei][mtype1].x += coeff * diff[mtype1][mtype2] * 2 * (dists[mtype1][mtype2].x)/(lengths[mtype1][mtype2]);
			gradptv[basei][mtype1].y += coeff * diff[mtype1][mtype2] * 2 * (dists[mtype1][mtype2].y)/(lengths[mtype1][mtype2]);
			gradptv[basei][mtype1].z += coeff * diff[mtype1][mtype2] * 2 * (dists[mtype1][mtype2].z)/(lengths[mtype1][mtype2]);
		}
		if(lengths[mtype1][mtype2] > 0){
			gradptv[basej][mtype2].x += coeff * diff[mtype1][mtype2] * 2 * (-dists[mtype1][mtype2].x)/(lengths[mtype1][mtype2]);
			gradptv[basej][mtype2].y += coeff * diff[mtype1][mtype2] * 2 * (-dists[mtype1][mtype2].y)/(lengths[mtype1][mtype2]);
			gradptv[basej][mtype2].z += coeff * diff[mtype1][mtype2] * 2 * (-dists[mtype1][mtype2].z)/(lengths[mtype1][mtype2]);
		}
	}
		
}

void RNASimulation::calcgradientditmplength(int basei, map<int, vector<Point3s>> &gradptv, vector<bool> movenodeinfo, vector<vector<double>> lengths, vector<vector<double>> diff, vector<vector<Point3s>> dists, int count, double coeff)
{
	int mtype1, mtype2;
	Point3s PQ;

	
	for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
	{
		for(mtype2=mtype1+1;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
		{
			if(lengths[mtype1][mtype2] > 0){
				gradptv[basei][mtype1].x += coeff * diff[mtype1][mtype2] * 2 * (dists[mtype1][mtype2].x)/(lengths[mtype1][mtype2]);
				gradptv[basei][mtype1].y += coeff * diff[mtype1][mtype2] * 2 * (dists[mtype1][mtype2].y)/(lengths[mtype1][mtype2]);
				gradptv[basei][mtype1].z += coeff * diff[mtype1][mtype2] * 2 * (dists[mtype1][mtype2].z)/(lengths[mtype1][mtype2]);
			}

			if(lengths[mtype1][mtype2] > 0){
				gradptv[basei][mtype2].x += coeff * diff[mtype1][mtype2] * 2 * (-dists[mtype1][mtype2].x)/(lengths[mtype1][mtype2]);
				gradptv[basei][mtype2].y += coeff * diff[mtype1][mtype2] * 2 * (-dists[mtype1][mtype2].y)/(lengths[mtype1][mtype2]);
				gradptv[basei][mtype2].z += coeff * diff[mtype1][mtype2] * 2 * (-dists[mtype1][mtype2].z)/(lengths[mtype1][mtype2]);
			}
		}
	}
}

/*
void RNASimulation::selectdnttemplate(Residueatoms restr1, Residueatoms restr2, vector<Residueatoms> &dnttemplate, TemplateDistMap &dnttemplatelength){
	char aat1=restr1.aaa;
	char aat2=restr2.aaa;

	dnttemplate.resize(2);

	int i, j, l, m;
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
	double rmsd;

	vector<double> tmp(3,0);
	//vector<vector<double> > xyz_list1(SNTTEMPLATEATOMS[aat].size(),tmp);
    //vector<vector<double> > xyz_list2(SNTTEMPLATEATOMS[aat].size(),tmp);
	//vector<vector<double> > xyz_listtmp(SNTTEMPLATEATOMS[aat].size(),tmp);
	vector<vector<double> > xyz_list1(30,tmp);
    vector<vector<double> > xyz_list2(30,tmp);
	vector<vector<double> > xyz_listtmp(15,tmp);

	//cout<<"dnt template"<<endl;
	//for(int j=0; j<SNTTEMPLATEATOMS[aat1].size();j++){
		//atom=SNTTEMPLATEATOMS[aat1][j];
		//cout<< "0" << "-" << j << "-"  << restr1.ptv[atom].x << "-" << restr1.ptv[atom].y << "-" << restr1.ptv[atom].z<<endl;
	//}
	//for(int j=0; j<SNTTEMPLATEATOMS[aat2].size();j++){
		//atom=SNTTEMPLATEATOMS[aat1][j];
		//cout<< "1" << "-" << j << "-"  << restr2.ptv[atom].x << "-" << restr2.ptv[atom].y << "-" << restr2.ptv[atom].z<<endl;
	//}

	for(i=0;i<DILENTTEMPLATE[aat1][aat2].size();i++){
		//for(int j=0; j<SNTTEMPLATEATOMS[aat].size();j++){
		count=0;
		for(j=0; j<15;j++){
			xyz_list1[count][0]=restr1.ptv[j].x;
			xyz_list1[count][1]=restr1.ptv[j].y;
			xyz_list1[count][2]=restr1.ptv[j].z;

			xyz_list2[count][0]=DILENTTEMPLATE[aat1][aat2][i][0].ptv[j].x;
			xyz_list2[count][1]=DILENTTEMPLATE[aat1][aat2][i][0].ptv[j].y;
			xyz_list2[count][2]=DILENTTEMPLATE[aat1][aat2][i][0].ptv[j].z;

			count+=1;
		}

		for(j=0; j<15;j++){
			xyz_list1[count][0]=restr2.ptv[j].x;
			xyz_list1[count][1]=restr2.ptv[j].y;
			xyz_list1[count][2]=restr2.ptv[j].z;

			xyz_list2[count][0]=DILENTTEMPLATE[aat1][aat2][i][1].ptv[j].x;
			xyz_list2[count][1]=DILENTTEMPLATE[aat1][aat2][i][1].ptv[j].y;
			xyz_list2[count][2]=DILENTTEMPLATE[aat1][aat2][i][1].ptv[j].z;

			count+=1;
		}

		TranVect.clear();
		TranVect.resize(3, 0.0);
		RotMatix.clear();
		RotMatix.resize(3, vector<double>(3, 0.0));

		bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

		count=0;
		rmsd=0;
		//for(int j=0; j<SNTTEMPLATEATOMS[aat].size();j++){
		for(j=0;j<15;j++){
			xyztmp[0]=0.0; xyztmp[1]=0.0; xyztmp[2]=0.0;
			bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);
			vec.x=restr1.ptv[j].x-xyztmp[0];
			vec.y=restr1.ptv[j].y-xyztmp[1];
			vec.z=restr1.ptv[j].z-xyztmp[2];
			rmsd+=bf.norm(vec);
			count+=1;
		}

		for(j=0;j<15;j++){
			xyztmp[0]=0.0; xyztmp[1]=0.0; xyztmp[2]=0.0;
			bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);
			vec.x=restr2.ptv[j].x-xyztmp[0];
			vec.y=restr2.ptv[j].y-xyztmp[1];
			vec.z=restr2.ptv[j].z-xyztmp[2];
			rmsd+=bf.norm(vec);
			count+=1;
		}

		if(minrmsd > rmsd){
			RotMatixtmp=RotMatix;
			TranVecttmp=TranVect;
			minidx=i;
			minrmsd=rmsd;
		}
	}

	for(int j=0; j<15;j++){
		xyz_listtmp[j][0]=DILENTTEMPLATE[aat1][aat2][minidx][0].ptv[j].x;
		xyz_listtmp[j][1]=DILENTTEMPLATE[aat1][aat2][minidx][0].ptv[j].y;
		xyz_listtmp[j][2]=DILENTTEMPLATE[aat1][aat2][minidx][0].ptv[j].z;

		bf.Transform(xyz_listtmp[j],RotMatixtmp,TranVecttmp,xyztmp);
		dnttemplate[0].ptv[j].x=xyztmp[0];
		dnttemplate[0].ptv[j].y=xyztmp[1];
		dnttemplate[0].ptv[j].z=xyztmp[2];
	}

	dnttemplate[0].len[1]=DILENTTEMPLATE[aat1][aat2][minidx][0].len[1];
	dnttemplate[0].len[2]=DILENTTEMPLATE[aat1][aat2][minidx][0].len[2];
	dnttemplate[0].ang[2]=DILENTTEMPLATE[aat1][aat2][minidx][0].ang[2];
	for(j=3;j<AccAtomParamNum[aat1][0];j++){
		dnttemplate[0].len[j]=DILENTTEMPLATE[aat1][aat2][minidx][0].len[j];
		dnttemplate[0].ang[j]=DILENTTEMPLATE[aat1][aat2][minidx][0].ang[j];
		dnttemplate[0].dih[j]=DILENTTEMPLATE[aat1][aat2][minidx][0].dih[j];
	}

	for(int j=0; j<15;j++){
		xyz_listtmp[j][0]=DILENTTEMPLATE[aat1][aat2][minidx][1].ptv[j].x;
		xyz_listtmp[j][1]=DILENTTEMPLATE[aat1][aat2][minidx][1].ptv[j].y;
		xyz_listtmp[j][2]=DILENTTEMPLATE[aat1][aat2][minidx][1].ptv[j].z;

		bf.Transform(xyz_listtmp[j],RotMatixtmp,TranVecttmp,xyztmp);
		dnttemplate[1].ptv[j].x=xyztmp[0];
		dnttemplate[1].ptv[j].y=xyztmp[1];
		dnttemplate[1].ptv[j].z=xyztmp[2];
	}

	for(j=0;j<AccAtomParamNum[aat2][0];j++){
		dnttemplate[1].len[j]=DILENTTEMPLATE[aat1][aat2][minidx][1].len[j];
		dnttemplate[1].ang[j]=DILENTTEMPLATE[aat1][aat2][minidx][1].ang[j];
		dnttemplate[1].dih[j]=DILENTTEMPLATE[aat1][aat2][minidx][1].dih[j];
	}

	dnttemplatelength=DILENTTEMPLATELENGTH[aat1][aat2][minidx];
}
*/

void RNASimulation::selectdnttemplate(Residueatoms restr1, Residueatoms restr2, int basei, int basej, vector<Residueatoms> &dnttemplate, TemplateDistMap &dnttemplatelength){
	char aat1=restr1.aaa;
	char aat2=restr2.aaa;

	dnttemplate.resize(2);
	int mtype1, mtype2;
	int i, j, l, m;
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
	double rmsd;

	vector<double> tmp(3,0);
	//vector<vector<double> > xyz_list1(SNTTEMPLATEATOMS[aat].size(),tmp);
    //vector<vector<double> > xyz_list2(SNTTEMPLATEATOMS[aat].size(),tmp);
	//vector<vector<double> > xyz_listtmp(SNTTEMPLATEATOMS[aat].size(),tmp);
	vector<vector<double> > xyz_list1(30,tmp);
    vector<vector<double> > xyz_list2(30,tmp);
	vector<vector<double> > xyz_listtmp(15,tmp);
	vector<vector<Residueatoms>> templatestructures;

	if(stpairednt[basei] != basej){
		if(basei < basej){
			templatestructures = DILENTTEMPLATE[aat1][aat2];
		}else{
			//cout<<basei << "-" << basej<<endl;
			templatestructures.resize(DILENTTEMPLATE[aat2][aat1].size());
			for(i=0;i<DILENTTEMPLATE[aat2][aat1].size();i++){
				templatestructures[i].resize(2);
				templatestructures[i][0]=DILENTTEMPLATE[aat2][aat1][i][1];
				templatestructures[i][1]=DILENTTEMPLATE[aat2][aat1][i][0];
			}
		}
		//templatestructures = STACKTEMPLATE[2][aat1][aat2];
	}else{
		if(basei+1 == basej){
			templatestructures = STACKTEMPLATE[1][aat1][aat2];
		}else if(basej+1 == basei){
			templatestructures.resize(STACKTEMPLATE[1][aat2][aat1].size());
			for(i=0;i<STACKTEMPLATE[1][aat2][aat1].size();i++){
				templatestructures[i].resize(2);
				templatestructures[i][0]=STACKTEMPLATE[1][aat2][aat1][i][1];
				templatestructures[i][1]=STACKTEMPLATE[1][aat2][aat1][i][0];
			}
		}else{
			templatestructures = STACKTEMPLATE[0][aat1][aat2];
		}
	}

	//cout<<"dnt template"<<endl;
	//for(int j=0; j<SNTTEMPLATEATOMS[aat1].size();j++){
		//atom=SNTTEMPLATEATOMS[aat1][j];
		//cout<< "0" << "-" << j << "-"  << restr1.ptv[atom].x << "-" << restr1.ptv[atom].y << "-" << restr1.ptv[atom].z<<endl;
	//}
	//for(int j=0; j<SNTTEMPLATEATOMS[aat2].size();j++){
		//atom=SNTTEMPLATEATOMS[aat1][j];
		//cout<< "1" << "-" << j << "-"  << restr2.ptv[atom].x << "-" << restr2.ptv[atom].y << "-" << restr2.ptv[atom].z<<endl;
	//}

	for(i=0;i<templatestructures.size();i++){
		//for(int j=0; j<SNTTEMPLATEATOMS[aat].size();j++){
		count=0;
		for(j=0; j<15;j++){
			xyz_list1[count][0]=restr1.ptv[j].x;
			xyz_list1[count][1]=restr1.ptv[j].y;
			xyz_list1[count][2]=restr1.ptv[j].z;

			xyz_list2[count][0]=templatestructures[i][0].ptv[j].x;
			xyz_list2[count][1]=templatestructures[i][0].ptv[j].y;
			xyz_list2[count][2]=templatestructures[i][0].ptv[j].z;

			//cout<<xyz_list2[count][0] << "-" << xyz_list2[count][1] << "-" << xyz_list2[count][2]<<endl;

			count+=1;
		}

		for(j=0; j<15;j++){
			xyz_list1[count][0]=restr2.ptv[j].x;
			xyz_list1[count][1]=restr2.ptv[j].y;
			xyz_list1[count][2]=restr2.ptv[j].z;

			xyz_list2[count][0]=templatestructures[i][1].ptv[j].x;
			xyz_list2[count][1]=templatestructures[i][1].ptv[j].y;
			xyz_list2[count][2]=templatestructures[i][1].ptv[j].z;

			//cout<<xyz_list2[count][0] << "-" << xyz_list2[count][1] << "-" << xyz_list2[count][2]<<endl;

			count+=1;
		}
		
		//if(stpairednt[basei] == basej && basej+1 == basei) exit(1);

		TranVect.clear();
		TranVect.resize(3, 0.0);
		RotMatix.clear();
		RotMatix.resize(3, vector<double>(3, 0.0));
		bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

		count=0;
		rmsd=0;
		//for(int j=0; j<SNTTEMPLATEATOMS[aat].size();j++){
		for(j=0;j<15;j++){
			xyztmp[0]=0.0; xyztmp[1]=0.0; xyztmp[2]=0.0;
			bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);
			vec.x=restr1.ptv[j].x-xyztmp[0];
			vec.y=restr1.ptv[j].y-xyztmp[1];
			vec.z=restr1.ptv[j].z-xyztmp[2];
			rmsd+=bf.norm(vec);
			count+=1;
		}

		for(j=0;j<15;j++){
			xyztmp[0]=0.0; xyztmp[1]=0.0; xyztmp[2]=0.0;
			bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);
			vec.x=restr2.ptv[j].x-xyztmp[0];
			vec.y=restr2.ptv[j].y-xyztmp[1];
			vec.z=restr2.ptv[j].z-xyztmp[2];
			rmsd+=bf.norm(vec);
			count+=1;
		}

		if(minrmsd > rmsd){
			RotMatixtmp=RotMatix;
			TranVecttmp=TranVect;
			minidx=i;
			minrmsd=rmsd;
		}
	}

	for(int j=0; j<15;j++){
		xyz_listtmp[j][0]=templatestructures[minidx][0].ptv[j].x;
		xyz_listtmp[j][1]=templatestructures[minidx][0].ptv[j].y;
		xyz_listtmp[j][2]=templatestructures[minidx][0].ptv[j].z;

		bf.Transform(xyz_listtmp[j],RotMatixtmp,TranVecttmp,xyztmp);
		dnttemplate[0].ptv[j].x=xyztmp[0];
		dnttemplate[0].ptv[j].y=xyztmp[1];
		dnttemplate[0].ptv[j].z=xyztmp[2];
	}

	dnttemplate[0].len[1]=templatestructures[minidx][0].len[1];
	dnttemplate[0].len[2]=templatestructures[minidx][0].len[2];
	dnttemplate[0].ang[2]=templatestructures[minidx][0].ang[2];
	for(j=3;j<AccAtomParamNum[aat1][0];j++){
		dnttemplate[0].len[j]=templatestructures[minidx][0].len[j];
		dnttemplate[0].ang[j]=templatestructures[minidx][0].ang[j];
		dnttemplate[0].dih[j]=templatestructures[minidx][0].dih[j];
	}

	for(int j=0; j<15;j++){
		xyz_listtmp[j][0]=templatestructures[minidx][1].ptv[j].x;
		xyz_listtmp[j][1]=templatestructures[minidx][1].ptv[j].y;
		xyz_listtmp[j][2]=templatestructures[minidx][1].ptv[j].z;

		bf.Transform(xyz_listtmp[j],RotMatixtmp,TranVecttmp,xyztmp);
		dnttemplate[1].ptv[j].x=xyztmp[0];
		dnttemplate[1].ptv[j].y=xyztmp[1];
		dnttemplate[1].ptv[j].z=xyztmp[2];
	}

	for(j=0;j<AccAtomParamNum[aat2][0];j++){
		dnttemplate[1].len[j]=templatestructures[minidx][1].len[j];
		dnttemplate[1].ang[j]=templatestructures[minidx][1].ang[j];
		dnttemplate[1].dih[j]=templatestructures[minidx][1].dih[j];
	}

	if(stpairednt[basei] != basej){
		if(basei < basej){
			dnttemplatelength=DILENTTEMPLATELENGTH[aat1][aat2][minidx];
		}else{
			for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
			{
				for(mtype2=0;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
				{
					dnttemplatelength.interdistmap[mtype1][mtype2]=DILENTTEMPLATELENGTH[aat2][aat1][minidx].interdistmap[mtype2][mtype1];
				}
			}
			dnttemplatelength.intradistmap[0]=DILENTTEMPLATELENGTH[aat2][aat1][minidx].intradistmap[1];
			dnttemplatelength.intradistmap[1]=DILENTTEMPLATELENGTH[aat2][aat1][minidx].intradistmap[0];
		}
	}else{
		if(basei+1 == basej){
			dnttemplatelength=STLENTTEMPLATELENGTH[1][aat1][aat2][minidx];
		}else if(basej+1 == basei){
			for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
			{
				for(mtype2=0;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
				{
					dnttemplatelength.interdistmap[mtype1][mtype2]=STLENTTEMPLATELENGTH[1][aat2][aat1][minidx].interdistmap[mtype2][mtype1];
				}
			}
			dnttemplatelength.intradistmap[0]=STLENTTEMPLATELENGTH[1][aat2][aat1][minidx].intradistmap[1];
			dnttemplatelength.intradistmap[1]=STLENTTEMPLATELENGTH[1][aat2][aat1][minidx].intradistmap[0];
		}else{
			dnttemplatelength=STLENTTEMPLATELENGTH[0][aat1][aat2][minidx];
		}
	}
}

void RNASimulation::selectbptemplate(Residueatoms restr1, Residueatoms restr2, vector<Residueatoms> &dnttemplate, TemplateDistMap &bptemplatelength){
	char aat1=restr1.aaa;
	char aat2=restr2.aaa;

	dnttemplate.resize(2);

	int i, j;
	Residueatoms templatetemp1, templatetemp2;
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
	double rmsd;
	vector<double> tmp(3,0);
	vector<vector<double> > xyz_list1(15,tmp);
    vector<vector<double> > xyz_list2(15,tmp);
	vector<vector<double> > xyz_listtmp(15,tmp);

	//cout<<"bp template"<<endl;
	//for(int j=0; j<SNTTEMPLATEATOMS[aat1].size();j++){
		//atom=SNTTEMPLATEATOMS[aat1][j];
		//cout<< "0" << "-" << j << "-" << restr1.ptv[atom].x << "-" << restr1.ptv[atom].y << "-" << restr1.ptv[atom].z<<endl;
	//}
	//for(int j=0; j<SNTTEMPLATEATOMS[aat2].size();j++){
		//atom=SNTTEMPLATEATOMS[aat2][j];
		//cout<< "1" << "-" << j << "-" << restr2.ptv[atom].x << "-" << restr2.ptv[atom].y << "-" << restr2.ptv[atom].z<<endl;
	//}
	

	for(i=0;i<BPLENTTEMPLATE[aat1][aat2].size();i++){
		//for(int j=0; j<SNTTEMPLATEATOMS[aat].size();j++){
		count=0;
		for(j=12; j<15;j++){
			atom=SNTTEMPLATEATOMS[aat1][j];
			xyz_list1[count][0]=restr1.ptv[atom].x;
			xyz_list1[count][1]=restr1.ptv[atom].y;
			xyz_list1[count][2]=restr1.ptv[atom].z;

			xyz_list2[count][0]=BPLENTTEMPLATE[aat1][aat2][i][0].ptv[atom].x;
			xyz_list2[count][1]=BPLENTTEMPLATE[aat1][aat2][i][0].ptv[atom].y;
			xyz_list2[count][2]=BPLENTTEMPLATE[aat1][aat2][i][0].ptv[atom].z;
			count+=1;
		}

		for(j=12; j<15;j++){
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

		rmsd=0;
		count=0;
		for(j=12;j<15;j++){
			atom=SNTTEMPLATEATOMS[aat1][j];
			bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);
			vec.x=restr1.ptv[atom].x-xyztmp[0];
			vec.y=restr1.ptv[atom].y-xyztmp[1];
			vec.z=restr1.ptv[atom].z-xyztmp[2];
			rmsd+=bf.norm(vec);
		}
		count=0;
		for(j=12; j<15;j++){
			atom=SNTTEMPLATEATOMS[aat2][j];
			bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);
			vec.x=restr2.ptv[atom].x-xyztmp[0];
			vec.y=restr2.ptv[atom].y-xyztmp[1];
			vec.z=restr2.ptv[atom].z-xyztmp[2];
			rmsd+=bf.norm(vec);
		}

		if(minrmsd > rmsd){
			RotMatixtmp=RotMatix;
			TranVecttmp=TranVect;
			minidx=i;
			minrmsd=rmsd;
		}
	}

	//cout<<"============"<<endl;

	for(int j=0; j<15;j++){
		atom=SNTTEMPLATEATOMS[aat1][j];
		xyz_listtmp[j][0]=BPLENTTEMPLATE[aat1][aat2][minidx][0].ptv[atom].x;
		xyz_listtmp[j][1]=BPLENTTEMPLATE[aat1][aat2][minidx][0].ptv[atom].y;
		xyz_listtmp[j][2]=BPLENTTEMPLATE[aat1][aat2][minidx][0].ptv[atom].z;

		bf.Transform(xyz_listtmp[j],RotMatixtmp,TranVecttmp,xyztmp);
		dnttemplate[0].ptv[atom].x=xyztmp[0];
		dnttemplate[0].ptv[atom].y=xyztmp[1];
		dnttemplate[0].ptv[atom].z=xyztmp[2];
		//cout<< "0" << "-" << j << "-" << dnttemplate[0].ptv[atom].x << "-" << dnttemplate[0].ptv[atom].y << "-" << dnttemplate[0].ptv[atom].z<<endl;
	}

	dnttemplate[0].len[1]=BPLENTTEMPLATE[aat1][aat2][minidx][0].len[1];
	dnttemplate[0].len[2]=BPLENTTEMPLATE[aat1][aat2][minidx][0].len[2];
	dnttemplate[0].ang[2]=BPLENTTEMPLATE[aat1][aat2][minidx][0].ang[2];
	for(j=3;j<AccAtomParamNum[aat1][0];j++){
		dnttemplate[0].len[j]=BPLENTTEMPLATE[aat1][aat2][minidx][0].len[j];
		dnttemplate[0].ang[j]=BPLENTTEMPLATE[aat1][aat2][minidx][0].ang[j];
		dnttemplate[0].dih[j]=BPLENTTEMPLATE[aat1][aat2][minidx][0].dih[j];
	}

	for(int j=0; j<15;j++){
	//for(j=0; j<6;j++){
		atom=SNTTEMPLATEATOMS[aat2][j];
		xyz_listtmp[j][0]=BPLENTTEMPLATE[aat1][aat2][minidx][1].ptv[atom].x;
		xyz_listtmp[j][1]=BPLENTTEMPLATE[aat1][aat2][minidx][1].ptv[atom].y;
		xyz_listtmp[j][2]=BPLENTTEMPLATE[aat1][aat2][minidx][1].ptv[atom].z;

		bf.Transform(xyz_listtmp[j],RotMatixtmp,TranVecttmp,xyztmp);
		dnttemplate[1].ptv[atom].x=xyztmp[0];
		dnttemplate[1].ptv[atom].y=xyztmp[1];
		dnttemplate[1].ptv[atom].z=xyztmp[2];
		//cout<< "1" << "-" << j << "-" << dnttemplate[1].ptv[atom].x << "-" << dnttemplate[1].ptv[atom].y << "-" << dnttemplate[1].ptv[atom].z<<endl;
	}

	for(j=0;j<AccAtomParamNum[aat2][0];j++){
		dnttemplate[1].len[j]=BPLENTTEMPLATE[aat1][aat2][minidx][1].len[j];
		dnttemplate[1].ang[j]=BPLENTTEMPLATE[aat1][aat2][minidx][1].ang[j];
		dnttemplate[1].dih[j]=BPLENTTEMPLATE[aat1][aat2][minidx][1].dih[j];
	}

	bptemplatelength=BPLENTTEMPLATELENGTH[aat1][aat2][minidx];

	//exit(1);
}

/*
void RNASimulation::selectbptemplate(Residueatoms restr1, Residueatoms restr2, vector<Residueatoms> &dnttemplate, TemplateDistMap &bptemplatelength){
	char aat1=restr1.aaa;
	char aat2=restr2.aaa;

	dnttemplate.resize(2);

	int i, j;
	Residueatoms templatetemp1, templatetemp2;
	int atom;
	vector<vector<double>> RotMatix;
    vector<double> TranVect;
	vector<vector<double>> RotMatixtmp;
    vector<double> TranVecttmp;
	vector<double> xyztmp(3,0);
	Point3s vec;
	double minrmsd=100000000000;
	int minidx=-1;

	vector<double> tmp(3,0);
	vector<vector<double> > xyz_list1(15,tmp);
    vector<vector<double> > xyz_list2(15,tmp);
	vector<vector<double> > xyz_listtmp(15,tmp);

	//cout<<"bp template"<<endl;
	//for(int j=0; j<SNTTEMPLATEATOMS[aat1].size();j++){
		//atom=SNTTEMPLATEATOMS[aat1][j];
		//cout<< "0" << "-" << j << "-" << restr1.ptv[atom].x << "-" << restr1.ptv[atom].y << "-" << restr1.ptv[atom].z<<endl;
	//}
	//for(int j=0; j<SNTTEMPLATEATOMS[aat2].size();j++){
		//atom=SNTTEMPLATEATOMS[aat2][j];
		//cout<< "1" << "-" << j << "-" << restr2.ptv[atom].x << "-" << restr2.ptv[atom].y << "-" << restr2.ptv[atom].z<<endl;
	//}
	

	for(i=0;i<BPLENTTEMPLATE[aat1][aat2].size();i++){
		//for(int j=0; j<SNTTEMPLATEATOMS[aat].size();j++){
		for(j=0; j<6;j++){
			atom=SNTTEMPLATEATOMS[aat1][j];
			xyz_list1[j][0]=restr1.ptv[atom].x;
			xyz_list1[j][1]=restr1.ptv[atom].y;
			xyz_list1[j][2]=restr1.ptv[atom].z;

			xyz_list2[j][0]=BPLENTTEMPLATE[aat1][aat2][i][0].ptv[atom].x;
			xyz_list2[j][1]=BPLENTTEMPLATE[aat1][aat2][i][0].ptv[atom].y;
			xyz_list2[j][2]=BPLENTTEMPLATE[aat1][aat2][i][0].ptv[atom].z;
		}

		for(j=0; j<6;j++){
			atom=SNTTEMPLATEATOMS[aat2][j];
			xyz_list1[6+j][0]=restr2.ptv[atom].x;
			xyz_list1[6+j][1]=restr2.ptv[atom].y;
			xyz_list1[6+j][2]=restr2.ptv[atom].z;

			xyz_list2[6+j][0]=BPLENTTEMPLATE[aat1][aat2][i][1].ptv[atom].x;
			xyz_list2[6+j][1]=BPLENTTEMPLATE[aat1][aat2][i][1].ptv[atom].y;
			xyz_list2[6+j][2]=BPLENTTEMPLATE[aat1][aat2][i][1].ptv[atom].z;
		}

		bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);


		//for(int j=0; j<SNTTEMPLATEATOMS[aat].size();j++){
		for(j=0; j<6;j++){
			atom=SNTTEMPLATEATOMS[aat1][j];
			bf.Transform(xyz_list2[j],RotMatix,TranVect,xyztmp);
			templatetemp1.ptv[atom].x=xyztmp[0];
			templatetemp1.ptv[atom].y=xyztmp[1];
			templatetemp1.ptv[atom].z=xyztmp[2];
		}

		for(j=0; j<6;j++){
			atom=SNTTEMPLATEATOMS[aat2][j];
			bf.Transform(xyz_list2[6+j],RotMatix,TranVect,xyztmp);
			templatetemp2.ptv[atom].x=xyztmp[0];
			templatetemp2.ptv[atom].y=xyztmp[1];
			templatetemp2.ptv[atom].z=xyztmp[2];
		}

		double rmsd=0;
		//for(int j=0; j<SNTTEMPLATEATOMS[aat].size();j++){
		for(j=0; j<6;j++){
			atom=SNTTEMPLATEATOMS[aat1][j];
			vec.x=restr1.ptv[atom].x-templatetemp1.ptv[atom].x;
			vec.y=restr1.ptv[atom].y-templatetemp1.ptv[atom].y;
			vec.z=restr1.ptv[atom].z-templatetemp1.ptv[atom].z;
			rmsd+=bf.norm(vec);
		}

		for(j=0; j<6;j++){
			atom=SNTTEMPLATEATOMS[aat2][j];
			vec.x=restr2.ptv[atom].x-templatetemp2.ptv[atom].x;
			vec.y=restr2.ptv[atom].y-templatetemp2.ptv[atom].y;
			vec.z=restr2.ptv[atom].z-templatetemp2.ptv[atom].z;
			rmsd+=bf.norm(vec);
		}

		if(minrmsd > rmsd){
			RotMatixtmp=RotMatix;
			TranVecttmp=TranVect;
			minidx=i;
			minrmsd=rmsd;
		}
	}

	//cout<<"============"<<endl;

	for(int j=0; j<SNTTEMPLATEATOMS[aat1].size();j++){
	//for(j=0; j<6;j++){
		atom=SNTTEMPLATEATOMS[aat1][j];
		xyz_listtmp[j][0]=BPLENTTEMPLATE[aat1][aat2][minidx][0].ptv[atom].x;
		xyz_listtmp[j][1]=BPLENTTEMPLATE[aat1][aat2][minidx][0].ptv[atom].y;
		xyz_listtmp[j][2]=BPLENTTEMPLATE[aat1][aat2][minidx][0].ptv[atom].z;

		bf.Transform(xyz_listtmp[j],RotMatixtmp,TranVecttmp,xyztmp);
		dnttemplate[0].ptv[atom].x=xyztmp[0];
		dnttemplate[0].ptv[atom].y=xyztmp[1];
		dnttemplate[0].ptv[atom].z=xyztmp[2];
		//cout<< "0" << "-" << j << "-" << dnttemplate[0].ptv[atom].x << "-" << dnttemplate[0].ptv[atom].y << "-" << dnttemplate[0].ptv[atom].z<<endl;
	}

	dnttemplate[0].len[1]=BPLENTTEMPLATE[aat1][aat2][minidx][0].len[1];
	dnttemplate[0].len[2]=BPLENTTEMPLATE[aat1][aat2][minidx][0].len[2];
	dnttemplate[0].ang[2]=BPLENTTEMPLATE[aat1][aat2][minidx][0].ang[2];
	for(j=3;j<AccAtomParamNum[aat1][0];j++){
		dnttemplate[0].len[j]=BPLENTTEMPLATE[aat1][aat2][minidx][0].len[j];
		dnttemplate[0].ang[j]=BPLENTTEMPLATE[aat1][aat2][minidx][0].ang[j];
		dnttemplate[0].dih[j]=BPLENTTEMPLATE[aat1][aat2][minidx][0].dih[j];
	}

	for(int j=0; j<SNTTEMPLATEATOMS[aat2].size();j++){
	//for(j=0; j<6;j++){
		atom=SNTTEMPLATEATOMS[aat2][j];
		xyz_listtmp[j][0]=BPLENTTEMPLATE[aat1][aat2][minidx][1].ptv[atom].x;
		xyz_listtmp[j][1]=BPLENTTEMPLATE[aat1][aat2][minidx][1].ptv[atom].y;
		xyz_listtmp[j][2]=BPLENTTEMPLATE[aat1][aat2][minidx][1].ptv[atom].z;

		bf.Transform(xyz_listtmp[j],RotMatixtmp,TranVecttmp,xyztmp);
		dnttemplate[1].ptv[atom].x=xyztmp[0];
		dnttemplate[1].ptv[atom].y=xyztmp[1];
		dnttemplate[1].ptv[atom].z=xyztmp[2];
		//cout<< "1" << "-" << j << "-" << dnttemplate[1].ptv[atom].x << "-" << dnttemplate[1].ptv[atom].y << "-" << dnttemplate[1].ptv[atom].z<<endl;
	}

	for(j=0;j<AccAtomParamNum[aat2][0];j++){
		dnttemplate[1].len[j]=BPLENTTEMPLATE[aat1][aat2][minidx][1].len[j];
		dnttemplate[1].ang[j]=BPLENTTEMPLATE[aat1][aat2][minidx][1].ang[j];
		dnttemplate[1].dih[j]=BPLENTTEMPLATE[aat1][aat2][minidx][1].dih[j];
	}

	bptemplatelength=BPLENTTEMPLATELENGTH[aat1][aat2][minidx];

	//exit(1);
}

*/

void RNASimulation::setgradientlmpbpparam(vector<Residueatoms> &restr, vector<pair<int,int>> pairs, vector<vector<int>> &lmpparamidx)
{
	int i, j;
	char aat1, aat2;
	int ls, minls;
	int idx;
	double randval;
	int mtype1, mtype2;
	double rmax, rnorm;
	Point3d rp;
	int basei, basej;
	int ntidx1, adjprevidx1, adjnextidx1, ntidx2, adjprevidx2, adjnextidx2;

	for(i=0;i<pairs.size();i++){
		basei = pairs[i].first;
		basej = pairs[i].second;
		//cout<<to_string(basei) + "-" + to_string(basej)<<endl;
		if(basej == -1){
			ntidx1=0;
			ntidx2=-1;
			adjprevidx2=-1;
			adjnextidx2=-1;

			if(bppairednt[basei-1] == -1){
				adjprevidx1=3;
			}else{
				adjprevidx1=2;
			}

			if(bppairednt[basei+1] == -1){
				adjnextidx1=3;
			}else{
				adjnextidx1=2;
			}

		}else{
			ntidx1=1;
			ntidx2=1;

			if(bppairednt[basei-1] == -1 | bppairednt[basei-1] != bppairednt[basej+1]){
				adjprevidx1=1;
				adjnextidx2=1;
			}else{
				adjprevidx1=0;
				adjnextidx2=0;
			}

			if(bppairednt[basei+1] == -1 | bppairednt[basei+1] != bppairednt[basej-1]){
				adjnextidx1=1;
				adjprevidx2=1;
			}else{
				adjnextidx1=0;
				adjprevidx2=0;
			}
		}
		lmpparamidx.push_back({ntidx1, adjprevidx1, adjnextidx1, ntidx2, adjprevidx2, adjnextidx2});
	}
}


void RNASimulation::setlmpbpparam(vector<Residueatoms> &restr, vector<pair<int,int>> pairs, vector<vector<int>> &lmpparamidx)
{
	int i, j;
	char aat1, aat2;
	int ls, minls;
	int idx;
	double randval;
	int mtype1, mtype2;
	double rmax, rnorm;
	Point3d rp;
	int basei, basej;
	int ntidx1, adjprevidx1, adjnextidx1, ntidx2, adjprevidx2, adjnextidx2;

	for(i=0;i<pairs.size();i++){
		basei = pairs[i].first;
		basej = pairs[i].second;
		//cout<<to_string(basei) + "-" + to_string(basej)<<endl;
		if(basej == -1){
			ntidx1=0;
			ntidx2=-1;
			adjprevidx2=-1;
			adjnextidx2=-1;

			if(bppairednt[basei-1] == -1){
				adjprevidx1=3;
			}else{
				adjprevidx1=2;
			}

			if(bppairednt[basei+1] == -1){
				adjnextidx1=3;
			}else{
				adjnextidx1=2;
			}

		}else{
			ntidx1=1;
			ntidx2=1;

			if(bppairednt[basei-1] == -1 | bppairednt[basei-1] != bppairednt[basej+1]){
				adjprevidx1=1;
				adjnextidx2=1;
			}else{
				adjprevidx1=0;
				adjnextidx2=0;
			}

			if(bppairednt[basei+1] == -1 | bppairednt[basei+1] != bppairednt[basej-1]){
				adjnextidx1=1;
				adjprevidx2=1;
			}else{
				adjnextidx1=0;
				adjprevidx2=0;
			}
		}
		lmpparamidx.push_back({ntidx1, adjprevidx1, adjnextidx1, ntidx2, adjprevidx2, adjnextidx2});
	}
}

void RNASimulation::selectsnttemplate(Residueatoms restr, Residueatoms &snttemplate){
	char aat=restr.aaa;
	Residueatoms templateorigin, templatetemp;
	int atom;
	vector<vector<double>> RotMatix;
    vector<double> TranVect;
	vector<vector<double>> RotMatixtmp;
    vector<double> TranVecttmp;
	vector<double> xyztmp(3,0);
	Point3s vec;
	double minrmsd=100000000000;
	int minidx=-1;

	vector<double> tmp(3,0);
	vector<vector<double> > xyz_list1(15,tmp);
    vector<vector<double> > xyz_list2(15,tmp);
	vector<vector<double> > xyz_listtmp(15,tmp);
	//vector<vector<double> > xyz_list1(6,tmp);
    //vector<vector<double> > xyz_list2(6,tmp);
	//vector<vector<double> > xyz_listtmp(6,tmp);

	for(int i=0;i<SINGLENTTEMPLATE[restr.ss][aat].size();i++){
		templateorigin=SINGLENTTEMPLATE[restr.ss][aat][i][0];

		for(int j=0; j<15;j++){
		//for(int j=0; j<6;j++){
			atom=SNTTEMPLATEATOMS[aat][j];
			xyz_list1[j][0]=restr.ptv[atom].x;
			xyz_list1[j][1]=restr.ptv[atom].y;
			xyz_list1[j][2]=restr.ptv[atom].z;

			xyz_list2[j][0]=SINGLENTTEMPLATE[restr.ss][aat][i][0].ptv[atom].x;
			xyz_list2[j][1]=SINGLENTTEMPLATE[restr.ss][aat][i][0].ptv[atom].y;
			xyz_list2[j][2]=SINGLENTTEMPLATE[restr.ss][aat][i][0].ptv[atom].z;
		}

		bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);


		for(int j=0; j<15;j++){
		//for(int j=0; j<6;j++){
			atom=SNTTEMPLATEATOMS[aat][j];
			bf.Transform(xyz_list2[j],RotMatix,TranVect,xyztmp);
			templatetemp.ptv[atom].x=xyztmp[0];
			templatetemp.ptv[atom].y=xyztmp[1];
			templatetemp.ptv[atom].z=xyztmp[2];
		}

		double rmsd=0;
		for(int j=0; j<15;j++){
		//for(int j=0; j<6;j++){
			atom=SNTTEMPLATEATOMS[aat][j];
			vec.x=restr.ptv[atom].x-templatetemp.ptv[atom].x;
			vec.y=restr.ptv[atom].y-templatetemp.ptv[atom].y;
			vec.z=restr.ptv[atom].z-templatetemp.ptv[atom].z;
			rmsd+=bf.norm(vec);
		}

		if(minrmsd > rmsd){
			RotMatixtmp=RotMatix;
			TranVecttmp=TranVect;
			minidx=i;
			minrmsd=rmsd;
		}
	}

	for(int j=0; j<15;j++){
	//for(int j=0; j<6;j++){
		atom=SNTTEMPLATEATOMS[aat][j];
		xyz_listtmp[j][0]=SINGLENTTEMPLATE[restr.ss][aat][minidx][0].ptv[atom].x;
		xyz_listtmp[j][1]=SINGLENTTEMPLATE[restr.ss][aat][minidx][0].ptv[atom].y;
		xyz_listtmp[j][2]=SINGLENTTEMPLATE[restr.ss][aat][minidx][0].ptv[atom].z;

		bf.Transform(xyz_listtmp[j],RotMatixtmp,TranVecttmp,xyztmp);
		snttemplate.ptv[atom].x=xyztmp[0];
		snttemplate.ptv[atom].y=xyztmp[1];
		snttemplate.ptv[atom].z=xyztmp[2];
	}
}

void RNASimulation::getgraphinfoLMP(vector<pair<int,int>> pairs, int numseq, vector<vector<int>> lmpparamidx, vector<bool> &movenodeinfo, vector<pair<int, int>> &edgeinfo, vector<int> &edgetype, vector<int> &paraminfo){

	vector<bool> checkidx(numseq, false);
	int basei, basej;
	int i;

	movenodeinfo.clear();
	movenodeinfo.resize(numseq, false);

	for(i=0;i<pairs.size();i++)
	{
		basei = pairs[i].first;
		basej = pairs[i].second;

		if(basej == -1){
			movenodeinfo[basei] = true;

			if(!checkidx[basei]){
				edgeinfo.push_back({basei, basei-1});
				paraminfo.push_back(lmpparamidx[i][1]);
				edgetype.push_back(1);
				checkidx[basei] = true;
			}

			edgeinfo.push_back({basei, basei});
			paraminfo.push_back(lmpparamidx[i][0]);
			edgetype.push_back(0);

			if(!checkidx[basei+1]){
				edgeinfo.push_back({basei+1, basei});
				paraminfo.push_back(lmpparamidx[i][2]);
				edgetype.push_back(1);
				checkidx[basei+1] = true;
			}
		}else{
			movenodeinfo[basei] = true;
			movenodeinfo[basej] = true;

			edgeinfo.push_back({basei, basej});
			paraminfo.push_back(-1);
			edgetype.push_back(2);

			if(!checkidx[basei]){
				edgeinfo.push_back({basei, basei-1});
				paraminfo.push_back(lmpparamidx[i][1]);
				edgetype.push_back(1);
				checkidx[basei] = true;
			}

			edgeinfo.push_back({basei, basei});
			paraminfo.push_back(lmpparamidx[i][0]);
			edgetype.push_back(0);

			if(!checkidx[basei+1]){
				edgeinfo.push_back({basei+1, basei});
				paraminfo.push_back(lmpparamidx[i][2]);
				edgetype.push_back(1);
				checkidx[basei+1] = true;
			}
			
			if(!checkidx[basej]){
				edgeinfo.push_back({basej, basej-1});
				paraminfo.push_back(lmpparamidx[i][4]);
				edgetype.push_back(1);
				checkidx[basej] = true;
			}

			edgeinfo.push_back({basej, basej});
			paraminfo.push_back(lmpparamidx[i][3]);
			edgetype.push_back(0);

			if(!checkidx[basej+1]){
				edgeinfo.push_back({basej+1, basej});
				paraminfo.push_back(lmpparamidx[i][5]);
				edgetype.push_back(1);
				checkidx[basej+1] = true;
			}
		}
	}
}

void RNASimulation::getgraphinfoLMP2(vector<pair<int,int>> pairs, int numseq, vector<vector<int>> lmpparamidx, vector<bool> &movenodeinfo, vector<pair<int, int>> &edgeinfo, vector<int> &edgetype, vector<int> &paraminfo){

	vector<bool> checkidx(numseq, false);
	int basei, basej;
	int i;

	movenodeinfo.clear();
	movenodeinfo.resize(numseq, false);

	for(i=0;i<pairs.size();i++)
	{
		basei = pairs[i].first;
		basej = pairs[i].second;

		if(basej == -1){
			movenodeinfo[basei] = true;

			if(!checkidx[basei]){
				edgeinfo.push_back({basei, basei-1});
				paraminfo.push_back(lmpparamidx[i][1]);
				edgetype.push_back(1);
				checkidx[basei] = true;
			}

			edgeinfo.push_back({basei, basei});
			paraminfo.push_back(lmpparamidx[i][0]);
			edgetype.push_back(0);

			if(!checkidx[basei+1]){
				edgeinfo.push_back({basei+1, basei});
				paraminfo.push_back(lmpparamidx[i][2]);
				edgetype.push_back(1);
				checkidx[basei+1] = true;
			}
		}else{
			movenodeinfo[basei] = true;

			edgeinfo.push_back({basei, basej});
			paraminfo.push_back(-1);
			edgetype.push_back(2);

			if(!checkidx[basei]){
				edgeinfo.push_back({basei, basei-1});
				paraminfo.push_back(lmpparamidx[i][1]);
				edgetype.push_back(1);
				checkidx[basei] = true;
			}

			edgeinfo.push_back({basei, basei});
			paraminfo.push_back(lmpparamidx[i][0]);
			edgetype.push_back(0);

			if(!checkidx[basei+1]){
				edgeinfo.push_back({basei+1, basei});
				paraminfo.push_back(lmpparamidx[i][2]);
				edgetype.push_back(1);
				checkidx[basei+1] = true;
			}
		}
	}
}

void RNASimulation::initgradgeo(vector<Residueatoms> &recrrstr, vector<bool> movenodeinfo, int numseq){
	int i,j;
	char aat;
	bool flagtor;
	double randoval;

	if(movenodeinfo[0] == 0){
		aat=recrrstr[0].aaa;
		for(j=0;j<INNERANGLETARGET.size();j++){
			if(lp1.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0] > 0){
				//randoval=Random();
				if(j == 10){
					if(checkoutlinerlgm(recrrstr[0].ang[j], OPINNERANGPARAM[aat][0].lower, OPINNERANGPARAM[aat][0].upper)){
					//if(randoval < 0.1){
						//cout<<"aaa"<<endl;
						//recrrstr[0].ang[j]=getrandomnewparam(OPINNERANGPARAM[aat][0].lower, OPINNERANGPARAM[aat][0].upper);
						recrrstr[0].ang[j]=getmovangval(recrrstr[0].ang[j], OPINNERANGPARAM[aat][0].mean, OPINNERANGPARAM[aat][0].std, 1);
					//}
					}
				}else if(j == 11){
					if(checkoutlinerlgm(recrrstr[0].ang[j], OPINNERANGPARAM[aat][1].lower, OPINNERANGPARAM[aat][1].upper)){
					//if(randoval < 0.1){
						//cout<<"aaa"<<endl;
						//recrrstr[0].ang[j]=getrandomnewparam(OPINNERANGPARAM[aat][1].lower, OPINNERANGPARAM[aat][1].upper);
						recrrstr[0].ang[j]=getmovangval(recrrstr[0].ang[j], OPINNERANGPARAM[aat][1].mean, OPINNERANGPARAM[aat][1].std, 1);
					//}
					}
				}else{
					if(checkoutlinerlgm(recrrstr[0].ang[j], INNERANGPARAM[aat][j].lower, INNERANGPARAM[aat][j].upper)){
					//if(randoval < 0.1){
						//cout<<"aaa"<<endl;
						//recrrstr[0].ang[j]=getrandomnewparam(INNERANGPARAM[aat][j].lower, INNERANGPARAM[aat][j].upper);
						recrrstr[0].ang[j]=getmovangval(recrrstr[0].ang[j], INNERANGPARAM[aat][j].mean, INNERANGPARAM[aat][j].std, 1);
					//}
					}
				}
			}
		}

		for(j=0;j<BLENGTHTARGET.size();j++){
			//if(lp1.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0] > 0){
				//randoval=Random();
				if(j == 10){
					if(checkoutlinerlgm(recrrstr[0].len[j], OPLENGTHPARAM[aat][0].lower, OPLENGTHPARAM[aat][0].upper)){
					//if(randoval < 0.1){
						//cout<<"aaa"<<endl;
						//recrrstr[0].len[j]=getrandomnewparam(OPLENGTHPARAM[aat][0].lower, OPLENGTHPARAM[aat][0].upper);
						recrrstr[0].len[j]=getmovlenval(recrrstr[0].len[j], OPLENGTHPARAM[aat][0].mean, OPLENGTHPARAM[aat][0].std, 1);
					//}
					}
				}else if(j == 11){
					if(checkoutlinerlgm(recrrstr[0].len[j], OPLENGTHPARAM[aat][1].lower, OPLENGTHPARAM[aat][1].upper)){
					//if(randoval < 0.1){
						//cout<<"aaa"<<endl;
						//recrrstr[0].len[j]=getrandomnewparam(OPLENGTHPARAM[aat][1].lower, OPLENGTHPARAM[aat][1].upper);
						recrrstr[0].len[j]=getmovlenval(recrrstr[0].len[j], OPLENGTHPARAM[aat][1].mean, OPLENGTHPARAM[aat][1].std, 1);
					//}
					}
				}else{
					if(checkoutlinerlgm(recrrstr[0].len[j], LENGTHPARAM[aat][j].lower, LENGTHPARAM[aat][j].upper)){
					//if(randoval < 0.1){
						//cout<<"aaa"<<endl;
						//recrrstr[0].len[j]=getrandomnewparam(LENGTHPARAM[aat][j].lower, LENGTHPARAM[aat][j].upper);
						recrrstr[0].len[j]=getmovlenval(recrrstr[0].len[j], LENGTHPARAM[aat][j].mean, LENGTHPARAM[aat][j].std, 1);
					//}
					}
				}
			//}
		}

		flagtor=param2coordbb(recrrstr,0,1);
		param2coordsugar(recrrstr,0,1);

		//flagtor=param2coordbb(recrrstr,0,1);
		//if(type > 0){
			//param2coordsugar(recrrstr,0,1);
			//param2coordbase(recrrstr,0,1);
			//param2coordhydrogen(recrrstr,0,1);
		//}
	}

	for(int i=1;i<numseq;i++){
		if(movenodeinfo[i]){
			aat=recrrstr[i].aaa;
			for(j=0;j<INNERANGLETARGET.size();j++){
				randoval=Random();
				//if(lp1.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][i] > 0){
					if(checkoutlinerlgm(recrrstr[i].ang[j], INNERANGPARAM[aat][j].lower, INNERANGPARAM[aat][j].upper)){
					//if(randoval < 0.1){
						//cout<<"aaa"<<endl;
						//recrrstr[i].ang[j]=getrandomnewparam(INNERANGPARAM[aat][j].lower, INNERANGPARAM[aat][j].upper);
						recrrstr[i].ang[j]=getmovangval(recrrstr[i].ang[j], INNERANGPARAM[aat][j].mean, INNERANGPARAM[aat][j].std, 1);
					//}
					}
				//}
			}

			for(j=0;j<BLENGTHTARGET.size();j++){
				randoval=Random();
				//if(lp1.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][i] > 0){
					if(checkoutlinerlgm(recrrstr[i].len[j], LENGTHPARAM[aat][j].lower, LENGTHPARAM[aat][j].upper)){
					//if(randoval < 0.1){
						//cout<<"aaa"<<endl;
						//recrrstr[i].len[j]=getrandomnewparam(LENGTHPARAM[aat][j].lower, LENGTHPARAM[aat][j].upper);
						recrrstr[i].len[j]=getmovlenval(recrrstr[i].len[j], LENGTHPARAM[aat][j].mean, LENGTHPARAM[aat][j].std, 1);
					//}
					}
				//}
			}
		}

		flagtor=param2coordbb(recrrstr,i,i+1);
		param2coordsugar(recrrstr,i,i+1);
	}

	//flagtor=param2coordbb(recrrstr,st,ed);
	//param2coordsugar(recrrstr,st,ed);
	//param2coordbase(recrrstr,st,ed);
	//param2coordhydrogen(recrrstr,st,ed);
}

void RNASimulation::initrandommove(vector<Residueatoms> &recrrstr, vector<bool> movenodeinfo, int numseq){
	double rmax, rnorm;
	Point3d rp;
	double randval=Random();

	//if(randval > 0.5){
	for(int i=0;i<numseq;i++)
	{
		if(movenodeinfo[i]){
			for(int j=0;j<PAIRLENGTHTARGETMOVNUM;j++){
				//rmax=0.5+0.5*Random();
				//rp=bf.ranv(rmax);
				//recrrstr[i].ptv[j].x+=rp.x;
				//recrrstr[i].ptv[j].y+=rp.y;
				//recrrstr[i].ptv[j].z+=rp.z;

				rmax=0.5+0.5*Random();
				rp=bf.ranv(rmax);
				recrrstr[i].ptv[j].x+=rp.x;
				recrrstr[i].ptv[j].y+=rp.y;
				recrrstr[i].ptv[j].z+=rp.z;
			}
		}
	}
	//}
}

















