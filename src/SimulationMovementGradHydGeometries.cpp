// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


void RNASimulation::gethbhydanglemse(Point3s P, Point3s &Q, Point3s R, double &angles, double angleideal, double &mse, double &theta, double &diff, double &normPQ, double &normRQ, double &dotp, Point3s &vecPQ, Point3s &vecRQ, double coeff){
    int i;
    Point3s PQ;
    Point3s RQ;
    //double normPQ, normRQ;
    //double dotp:
    //double theta, degree, diff;

    
    PQ.x = Q.x - P.x; 
	PQ.y = Q.y - P.y; 
	PQ.z = Q.z - P.z;
    RQ.x = Q.x - R.x; 
	RQ.y = Q.y - R.y; 
	RQ.z = Q.z - R.z;

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
        
    angles = theta*RAD2DEG;

    diff=angles-angleideal;
    mse += coeff*pow(diff, 2);
}

void RNASimulation::gethbhydlengthmse(Point3s P, Point3s &Q, double &lengths, double lengthideal, double &mse, double &diff, Point3s vecPQ, double coeff){
    int i;
    Point3s PQ;

    //PQ.x = P.x - Q.x; PQ.y = P.y - Q.y; PQ.z = P.z - Q.z;
	PQ.x=-vecPQ.x; PQ.y=-vecPQ.y; PQ.z=-vecPQ.z;
    lengths = bf.norm(PQ);
    diff = lengths-lengthideal;
    mse += coeff*diff*diff;
}

void RNASimulation::initializegrad(Point3s &gradQ, double &mse){
    gradQ.x=0;
    gradQ.y=0;
    gradQ.z=0;
    mse=0;
}

void RNASimulation::calchbhydanglemsegradient(Point3s P, Point3s &Q, Point3s R, double &angles, double angleideal, Point3s& gradQ, double theta, double diff, double normPQ, double normRQ, double dotp, Point3s &vecPQ, Point3s &vecRQ, double coeff)
{
    int i;
    double ddotp, dnorm;
    double normRQ_PQ, normPQ_RQ, normPQRQ2, dcos, dmse, k;

    normRQ_PQ=(normRQ/normPQ);
	normPQ_RQ=(normPQ/normRQ);
	normPQRQ2=pow(normPQ * normRQ,2);
	dcos=(-1/sin(theta));
	dmse = 2 * (diff);

	k = dmse*RAD2DEG*dcos/normPQRQ2;

	ddotp=2*vecPQ.x-R.x;
	//dnorm=normRQ_PQ * (Q.x - P.x) +  normPQ_RQ * (Q.x - R.x);
	dnorm=normRQ_PQ * (vecPQ.x) +  normPQ_RQ * (vecRQ.x);
	gradQ.x += coeff * k * (ddotp*(normPQ * normRQ) - dotp*dnorm);

	ddotp=2*vecPQ.y-R.y;
	//dnorm=normRQ_PQ * (Q.y - P.y) + normPQ_RQ * (Q.y - R.y);
	dnorm=normRQ_PQ * (vecPQ.y) + normPQ_RQ * (vecRQ.y);
	gradQ.y += coeff * k * (ddotp*(normPQ * normRQ) - dotp*dnorm);

	ddotp=2*vecPQ.z-R.z;
	//dnorm=normRQ_PQ * (Q.z - P.z) + normPQ_RQ * (Q.z - R.z);
	dnorm=normRQ_PQ * (vecPQ.z) + normPQ_RQ * (vecRQ.z);
	gradQ.z += coeff * k * (ddotp*(normPQ * normRQ) - dotp*dnorm);
}

void RNASimulation::calchbhydlengthmsegradient(Point3s P, Point3s &Q, double &lengths, double lengthideal, Point3s &gradQ, double &diff, Point3s vecPQ, double coeff){
	//gradQ.x += coeff * diff * 2 * (Q.x - P.x)/lengths;
	//gradQ.y += coeff * diff * 2 * (Q.y - P.y)/lengths;
	//gradQ.z += coeff * diff * 2 * (Q.z - P.z)/lengths;
	gradQ.x += coeff * diff * 2 * (vecPQ.x)/lengths;
	gradQ.y += coeff * diff * 2 * (vecPQ.y)/lengths;
	gradQ.z += coeff * diff * 2 * (vecPQ.z)/lengths;
}


void RNASimulation::applyhbhydSGD(Point3s P, Point3s &Q, Point3s R, double &angles, double &lengths1, double &lengths2, double angleideal, double lengthideal1, double lengthideal2, double learning_rate, int epochs, double tolerance){
	int i;
	int epoch;
    Point3s gradQ;
    double mse_prev=0;
    double mse=0;
    double theta, diffangle, difflength1, difflength2, normPQ, normRQ, dotp;
	Point3s vecPQ, vecRQ;

	//cout<< "initial coordinate" <<endl;
    //cout<< "Px:" << P.x << " Py" << P.y << " Pz" << P.z<<endl;
    //cout<< "Qx:" << Q.x << " Qy" << Q.y << " Qz" << Q.z<<endl;
    //cout<< "Rx:" << R.x << " Ry" << R.y << " Rz" << R.z<<endl;

	gethbhydanglemse(P, Q, R, angles, angleideal, mse_prev, theta, diffangle, normPQ, normRQ, dotp, vecPQ, vecRQ, 0.1);
    gethbhydlengthmse(R, Q, lengths1, lengthideal1, mse_prev, difflength1, vecRQ, 1);
	gethbhydlengthmse(P, Q, lengths2, lengthideal2, mse_prev, difflength2, vecPQ, 1);

	for (epoch = 0; epoch < epochs; ++epoch) {
		initializegrad(gradQ, mse);
		calchbhydanglemsegradient(P, Q, R, angles, angleideal, gradQ, theta, diffangle, normPQ, normRQ, dotp, vecPQ, vecRQ, 0.1);
        calchbhydlengthmsegradient(R, Q, lengths1, lengthideal1, gradQ, difflength1, vecRQ, 1);
		calchbhydlengthmsegradient(P, Q, lengths2, lengthideal2, gradQ, difflength2, vecPQ, 1);

		Q.x -= learning_rate*gradQ.x;
    	Q.y -= learning_rate*gradQ.y;
        Q.z -= learning_rate*gradQ.z;

		gethbhydanglemse(P, Q, R, angles, angleideal, mse, theta, diffangle, normPQ, normRQ, dotp, vecPQ, vecRQ, 0.1);
        gethbhydlengthmse(R, Q, lengths1, lengthideal1, mse, difflength1, vecRQ, 1);
		gethbhydlengthmse(P, Q, lengths2, lengthideal2, mse, difflength2, vecPQ, 1);

		//cout << "mse: " << mse << " ang: " << angles << " length1: " << lengths1 << " length2: " << lengths2 << " gradx: " << gradQ.x << " grady: " << gradQ.y << " gradz: " << gradQ.z << endl;
		//cout << "mse: " << mse << " ang: " << angles << " length1: " << lengths1 << " length2: " << lengths2 << endl;

		if(abs(mse_prev - mse) < tolerance){
            //cout << "Convergence reached at epoch " << epoch << ", MSE: " << mse << endl;
            //cout<< "final coordinate" <<endl;
            //cout<< "Px:" << P.x << " Py" << P.y << " Pz" << P.z<<endl;
            //cout<< "Qx:" << Q.x << " Qy" << Q.y << " Qz" << Q.z<<endl;
            //cout<< "Rx:" << R.x << " Ry" << R.y << " Rz" << R.z<<endl;
			break;
        }

        mse_prev = mse;
	}

	//cout << "mse: " << mse << " idealang: " << angleideal << " ideallength1: " << lengthideal1 << " ideallength2: " << lengthideal2 << endl;

	//cout << "Convergence reached at epoch " << epoch << ", MSE: " << mse << endl;
	//cout<< "Px:" << P.x << " Py" << P.y << " Pz" << P.z<<endl;
	//cout<< "Qx:" << Q.x << " Qy" << Q.y << " Qz" << Q.z<<endl;
	//cout<< "Rx:" << R.x << " Ry" << R.y << " Rz" << R.z<<endl;

}

bool RNASimulation::fmovementgradhydgeo(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc)//length 
{
	int i,j;
	double r,l;
	double nr,nl;
	int trandpos;
	bool flagtor;
	int wpos;
	int ss;
	char aat;
	double lengval, angval, torval;
	Point3s pn;
	//int totwrong=lp1.nn[BLENGTH_ENERGY_ST+paramtype]+lp1.nn[INNERANGLE_ENERGY_ST+paramtype]+lp1.nn[TORSIONANGLE_ENERGY_ST+paramtype];
	//int totwrong=lp1.nn[BLENGTH_ENERGY_ST+paramtype]+lp1.nn[INNERANGLE_ENERGY_ST+paramtype];
	int wrong=0;
	double randval;
	int st;
	bool flagts;
	BasicFunc bf;
	vector<vector<double>> hatomtype;
	int hydtype;
	char aat1, aat2;
	int basei, basej;
	Point3s tmpd,tmpa,tmph;
	int atom1, atom2, atom3;
	int dbase, hbase, abase;
	double learning_rate=0.0000001;
	int epochs=5000;
	double tolerance=0.0001;
	double actangle, idealangle, actlength1, actlength2, ideallength1, ideallength2;
	int inext;

	int totwrong=lp1.nn[BASEHYD_ENERGY_ST]*4;
	randval=Random();

	if(BASEPAIRS.size() == 0){
		newenergy=oldenergy;
		return false;
	}

	if(totwrong>0 && randval<0.9){
		wpos=int(Random()*totwrong);

		i=0;
		for(trandpos=0;trandpos<BASEPAIRS.size();trandpos++)
		{
			basei=BASEPAIRS[trandpos].basepidx.first;//base idx 1
			basej=BASEPAIRS[trandpos].basepidx.second;//base idx 2

			inext=i+lp1.indn[BASEHYD_ENERGY_ST][basei];
			if(i<=wpos && wpos<inext){
				break;
			}
			i=inext;

			inext=i+lp1.indn[BASEHYD_ENERGY_ST][basej];
			if(i<=wpos && wpos<inext){
				break;
			}
			i=inext;
		}
	}
	else{
		trandpos=int(BASEPAIRS.size()*Random()-EPS);
	}

	if(trandpos >= BASEPAIRS.size()){
		trandpos=int(BASEPAIRS.size()*Random()-EPS);
	}

	copystr(decstr, numseq);
	
	basei=BASEPAIRS[trandpos].basepidx.first;//base idx 1
	basej=BASEPAIRS[trandpos].basepidx.second;//base idx 2

	aat1=decstr[basei].aaa;
	aat2=decstr[basej].aaa;
	hatomtype=BASEPAIRS[trandpos].hatomtype;

	for(j=0;j<hatomtype.size();j++){
		//cout<<to_string(i) + "-" + to_string(hatomtype[j][0]) + "-" + to_string(hatomtype[j][1]) + "-" + to_string(hatomtype[j][2]) + "-" + to_string(hatomtype[j][3]) + "-" + to_string(hatomtype[j][4]) + "-" + to_string(hatomtype[j][5]) + "-" + to_string(hatomtype[j][6]) + "-" + to_string(hatomtype[j][7]) + "-" + to_string(hatomtype[j][8])<<endl;

		dbase=hatomtype[j][0];
		hbase=hatomtype[j][1];
		abase=hatomtype[j][2];

		atom1=hatomtype[j][3];
		atom2=hatomtype[j][4];
		atom3=hatomtype[j][5];

		tmpd=restr[dbase].ptv[atom1];
		tmph=restr[hbase].ptv[atom2];
		tmpa=restr[abase].ptv[atom3];

		idealangle=hatomtype[j][6];
		ideallength1=hatomtype[j][7];
		ideallength2=HYDLENGTHPARAM[decstr[hbase].aaa][HYDATOMID2HYDPARAMID[decstr[hbase].aaa][atom2]].mode;

		applyhbhydSGD(tmpd, tmph, tmpa, actangle, actlength1, actlength2, idealangle, ideallength1, ideallength2, learning_rate, epochs, tolerance);

		restr[hbase].ptv[atom2].x=tmph.x;
		restr[hbase].ptv[atom2].y=tmph.y;
		restr[hbase].ptv[atom2].z=tmph.z;
	}

	coord2paramhydrogenBase(restr,basei,basei+1,0);
	coord2paramhydrogenBase(restr,basej,basej+1,0);

	for(j=0;j<numseq;j++)
	{
		mapflagbk[j][basei]=true;
		mapflagbk[basei][j]=true;
		mapflagbk[j][basej]=true;
		mapflagbk[basej][j]=true;
	}

	newenergy=fcalcrmsdenergy(restr, numseq, lp2, rppmatbk, bpenergybk, bsenergybk, nmapbk, distmapbk, ntpairdistbk, adjpairdistbk, distenebk, multidistenebk, contactenebk, mapflagbk, mapbbflagbk, flagacc);
	bool upidx=metropolisUpdate(oldenergy, newenergy, Bt);
	
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
	
	if(upidx){
		rcopystr(decstr, restr, numseq);
		return true;
	}
	else
	{
		return false;
	}
}












