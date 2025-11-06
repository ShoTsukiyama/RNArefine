#include "MonteCarloSimulation.hpp"

void MonteCarloSimulation::getmolConnect(int res1, int res2, vector<vector<int>> &molConnect, map<string, vector<vector<bool>>> &intraAtomConnect, const ExtParams &extparameters)
{
    int i, j, k;
    char base1, base2, base;
	string base1_str, base2_str, base_str;
    int atomnum1, atomnum2, atomnum;
	int idx1, idx2;
    vector<string> nts={"A","C","G","U"};
    map<char, int> atomnumnts={{'A', 33}, {'C', 31}, {'G', 34}, {'U', 30}};

    for(i=0;i<nts.size();i++){
        base=nts[i][0];
		base_str=string(1, base);
        atomnum = atomnumnts[base];
        intraAtomConnect[base_str].resize(atomnum, vector<bool>(atomnum, false));
        for(k=0;k<extparameters.molConnectone.at(base_str).size();k++){
			idx1 = atomid2mcidx[base][extparameters.molConnectone.at(base_str)[k][0]];
			idx2 = atomid2mcidx[base][extparameters.molConnectone.at(base_str)[k][1]];
            intraAtomConnect[base_str][idx1][idx2] = true;
            intraAtomConnect[base_str][idx2][idx1] = true;
        }
    }

	atomnum1 = sizeof(refinestr[res1].ptv) / sizeof(refinestr[res1].ptv[0]);
	atomnum2 = sizeof(refinestr[res2].ptv) / sizeof(refinestr[res2].ptv[0]);
	molConnect.clear();
    molConnect.resize(atomnum1, vector<int>(atomnum2, 4));

	if(res1 == res2){
		base1 = refinestr[res1].aaa;
		base2 = refinestr[res2].aaa;
		base1_str=string(1, base1);
		base2_str=string(1, base2);

		for(k=0;k<extparameters.molConnectthree.at(base1_str).size();k++){
			idx1 = atomid2mcidx[base1][extparameters.molConnectthree.at(base1_str)[k][0]];
			idx2 = atomid2mcidx[base2][extparameters.molConnectthree.at(base1_str)[k][3]];
			molConnect[idx1][idx2] = 3;
			molConnect[idx2][idx1] = 3;
		}
		for(k=0;k<extparameters.molConnecttwo.at(base1_str).size();k++){
			idx1 = atomid2mcidx[base1][extparameters.molConnecttwo.at(base1_str)[k][0]];
			idx2 = atomid2mcidx[base2][extparameters.molConnecttwo.at(base1_str)[k][2]];
			molConnect[idx1][idx2] = 2;
			molConnect[idx2][idx1] = 2;
		}
		for(k=0;k<extparameters.molConnectone.at(base1_str).size();k++){
			idx1 = atomid2mcidx[base1][extparameters.molConnectone.at(base1_str)[k][0]];
			idx2 = atomid2mcidx[base2][extparameters.molConnectone.at(base1_str)[k][1]];
			molConnect[idx1][idx2] = 1;
			molConnect[idx2][idx1] = 1;
		}
		for(k=0;k<atomnum1;k++){
			molConnect[k][k] = 0;
		}
	}
	if(res2 == res1+1){
		base1 = refinestr[res1].aaa;
		base2 = refinestr[res2].aaa;
		base1_str=string(1, base1);
		base2_str=string(1, base2);

		for(k=0;k<extparameters.molAdjConnectthree.at(base2_str).size();k++){
			idx1 = atomid2mcidx[base1][extparameters.molAdjConnectthree.at(base2_str)[k][0]];
			idx2 = atomid2mcidx[base2][extparameters.molAdjConnectthree.at(base2_str)[k][3]];
			molConnect[idx1][idx2] = 3;
		}
		for(k=0;k<extparameters.molAdjConnecttwo.at(base2_str).size();k++){
			idx1 = atomid2mcidx[base1][extparameters.molAdjConnecttwo.at(base2_str)[k][0]];
			idx2 = atomid2mcidx[base2][extparameters.molAdjConnecttwo.at(base2_str)[k][2]];
			molConnect[idx1][idx2] = 2;
		}
		for(k=0;k<extparameters.molAdjConnectone.at(base2_str).size();k++){
			idx1 = atomid2mcidx[base1][extparameters.molAdjConnectone.at(base2_str)[k][0]];
			idx2 = atomid2mcidx[base2][extparameters.molAdjConnectone.at(base2_str)[k][1]];
			molConnect[idx1][idx2] = 1;
		}
	}
}

double ClashMCeneobj::calcenergyclash() const
{
	if(*flag){
		Point3s vec12;
		double newene;
		double deltaene;

		vec12.x = teAtomB->x - teAtomA->x;
		vec12.y = teAtomB->y - teAtomA->y;
		vec12.z = teAtomB->z - teAtomA->z;

		double length12=sqrt(vec12.x*vec12.x+vec12.y*vec12.y+vec12.z*vec12.z);

		if(length12 < vdW){
			newene = vdW-length12;
			newene = sqrt(newene * newene);
			return newene;
		}else{
			return 0;
		}
	}else{
		return prevene;
	}
}

void ClashMCeneobj::setprevene(const double &_prevene) const
{
	prevene=_prevene;
	*atomenergy1+=_prevene;
	*atomenergy2+=_prevene;
}

double InitMCdistobj::calcenergyinit() const
{
	if(*flag){
		Point3s vec;
		vec.x = (*teAtomB).x - (*teAtomA).x;
		vec.y = (*teAtomB).y - (*teAtomA).y;
		vec.z = (*teAtomB).z - (*teAtomA).z;

		double dist = sqrt(vec.x*vec.x+vec.y*vec.y+vec.z*vec.z);
		double diff = dist - initdist;
		diff = sqrt(diff*diff);

		if(diff > cutoff){
			return K*cutoff;
		}else{
			return K*diff;
		}
	}else{
		return prevene;
	}
}

void InitMCdistobj::setprevene(const double &_prevene) const
{
	prevene=_prevene;
	*atomenergy1+=_prevene;
	*atomenergy2+=_prevene;
}

double InitMCPosobj::calcenergyinit() const
{
	if(*flag){
		Point3s vec;
		vec.x = X0 - (*teAtom).x;
		vec.y = Y0 - (*teAtom).y;
		vec.z = Z0 - (*teAtom).z;

		double diff = sqrt(vec.x*vec.x+vec.y*vec.y+vec.z*vec.z+1e-8);

		if(diff > cutoff){
			return K*cutoff;
		}else{
			return K*diff;
		}
	}else{
		return prevene;
	}
}

void InitMCPosobj::setprevene(const double &_prevene) const
{
	prevene=_prevene;
	*atomenergy+=_prevene;
}


double ElectrosMCeneobj::calcenergyelectros() const
{
	if(*flag){
		Point3s vec12;
		double newene;
		double deltaene;

		vec12.x = teAtomB->x - teAtomA->x;
		vec12.y = teAtomB->y - teAtomA->y;
		vec12.z = teAtomB->z - teAtomA->z;

		double length12=sqrt(vec12.x*vec12.x+vec12.y*vec12.y+vec12.z*vec12.z);

		return dQQep/length12;
	}else{
		return prevene;
	}
}

void ElectrosMCeneobj::setprevene(const double &_prevene) const
{
	prevene=_prevene;
	*atomenergy1+=min(exp(_prevene/10), 1000.0);
	*atomenergy2+=min(exp(_prevene/10), 1000.0);
}

double LJMCeneobj::calcenergylj() const
{
	if(*flag){
		Point3s vec12;
		double newene;
		double deltaene;

		vec12.x = teAtomB->x - teAtomA->x;
		vec12.y = teAtomB->y - teAtomA->y;
		vec12.z = teAtomB->z - teAtomA->z;

		double length12SQR=vec12.x*vec12.x+vec12.y*vec12.y+vec12.z*vec12.z;

		if(length12SQR > CUTOFF_SQR) return 0.0;
		return (LJ6/(length12SQR*length12SQR*length12SQR))+(LJ12/(length12SQR*length12SQR*length12SQR*length12SQR*length12SQR*length12SQR));
	}else{
		return prevene;
	}
}

void LJMCeneobj::setprevene(const double &_prevene) const
{
	prevene=_prevene;
	*atomenergy1+=min(exp(_prevene/10), 1000.0);
	*atomenergy2+=min(exp(_prevene/10), 1000.0);
}

double HBondMCeneobj::calcenergyhydbond() const
{
	//cout<<teAtomDon->x << "-" << teAtomDon->y << "-" << teAtomDon->z<<endl;
	//cout<<teAtomHyd->x << "-" << teAtomHyd->y << "-" << teAtomHyd->z<<endl;
	//cout<<teAtomAcc->x << "-" << teAtomAcc->y << "-" << teAtomAcc->z<<endl;

	//cout<<idealahleng<< "-" << idealadleng<<endl;

	if(*flag){
		Point3s vec21, vec23, vec13, grad1, grad2, grad3, gradsum;

		vec21.x = teAtomDon->x - teAtomHyd->x;
		vec21.y = teAtomDon->y - teAtomHyd->y;
		vec21.z = teAtomDon->z - teAtomHyd->z;

		vec23.x = teAtomAcc->x - teAtomHyd->x;
		vec23.y = teAtomAcc->y - teAtomHyd->y;
		vec23.z = teAtomAcc->z - teAtomHyd->z;

		vec13.x = teAtomDon->x - teAtomAcc->x;
		vec13.y = teAtomDon->y - teAtomAcc->y;
		vec13.z = teAtomDon->z - teAtomAcc->z;

		double length21=sqrt(vec21.x*vec21.x+vec21.y*vec21.y+vec21.z*vec21.z+1e-8);
		double length23=sqrt(vec23.x*vec23.x+vec23.y*vec23.y+vec23.z*vec23.z+1e-8);
		double length13=sqrt(vec13.x*vec13.x+vec13.y*vec13.y+vec13.z*vec13.z+1e-8);

		vec21.x /= length21; vec21.y /= length21; vec21.z /= length21;
		vec23.x /= length23; vec23.y /= length23; vec23.z /= length23;
		
		double innerprot = vec21.x*vec23.x+vec21.y*vec23.y+vec21.z*vec23.z;
		double dcosangle = innerprot;
		if     (dcosangle >  1.) dcosangle =  1.;
		else if(dcosangle < -1.) dcosangle = -1.;

		if(dcosangle < idealahdang) dcosangle = -1;

		double diffahleng = length23-idealahleng;
		double diffadleng = length13-idealadleng;

		//cout<<-dK*exp(-1*(diffahleng*diffahleng + diffadleng*diffadleng) + dXiTheta*dcosangle)<<endl;
		//cout<<dK<<endl;
		return -dK*exp(-1*(diffahleng*diffahleng + diffadleng*diffadleng) - dXiTheta*dcosangle);
	}else{
		return prevene;
	}
}

void HBondMCeneobj::setprevene(const double &_prevene) const
{
	prevene=_prevene;
	*atomenergy2+=min(exp(_prevene),1000.0);
	*atomenergy3+=min(exp(_prevene),1000.0);
}

double HBond2MCeneobj::calcenergyhydbond() const
{
	//cout<<teAtomDon->x << "-" << teAtomDon->y << "-" << teAtomDon->z<<endl;
	//cout<<teAtomHyd->x << "-" << teAtomHyd->y << "-" << teAtomHyd->z<<endl;
	//cout<<teAtomAcc->x << "-" << teAtomAcc->y << "-" << teAtomAcc->z<<endl;

	//cout<<idealahleng<< "-" << idealadleng<<endl;

	if(*flag){
		Point3s vec23_1, vec23_2, vec23, vec13, vec21, grad1, grad2, grad3, gradsum;

		vec23_1.x = teAtomAcc->x - teAtomHyd1->x;
		vec23_1.y = teAtomAcc->y - teAtomHyd1->y;
		vec23_1.z = teAtomAcc->z - teAtomHyd1->z;

		vec23_2.x = teAtomAcc->x - teAtomHyd2->x;
		vec23_2.y = teAtomAcc->y - teAtomHyd2->y;
		vec23_2.z = teAtomAcc->z - teAtomHyd2->z;

		double length23_1=sqrt(vec23_1.x*vec23_1.x+vec23_1.y*vec23_1.y+vec23_1.z*vec23_1.z+1e-8);
		double length23_2=sqrt(vec23_2.x*vec23_2.x+vec23_2.y*vec23_2.y+vec23_2.z*vec23_2.z+1e-8);

		int hydatom;
		if(length23_1 < length23_2){
			vec21.x = teAtomDon->x - teAtomHyd1->x;
			vec21.y = teAtomDon->y - teAtomHyd1->y;
			vec21.z = teAtomDon->z - teAtomHyd1->z;

			vec23.x = vec23_1.x;
			vec23.y = vec23_1.y;
			vec23.z = vec23_1.z;

			hydatom = -1;
		}else{
			vec21.x = teAtomDon->x - teAtomHyd2->x;
			vec21.y = teAtomDon->y - teAtomHyd2->y;
			vec21.z = teAtomDon->z - teAtomHyd2->z;

			vec23.x = vec23_2.x;
			vec23.y = vec23_2.y;
			vec23.z = vec23_2.z;

			hydatom = 1;
		}

		vec13.x = teAtomDon->x - teAtomAcc->x;
		vec13.y = teAtomDon->y - teAtomAcc->y;
		vec13.z = teAtomDon->z - teAtomAcc->z;

		double length21=sqrt(vec21.x*vec21.x+vec21.y*vec21.y+vec21.z*vec21.z+1e-8);
		double length23=sqrt(vec23.x*vec23.x+vec23.y*vec23.y+vec23.z*vec23.z+1e-8);
		double length13=sqrt(vec13.x*vec13.x+vec13.y*vec13.y+vec13.z*vec13.z+1e-8);

		vec21.x /= length21; vec21.y /= length21; vec21.z /= length21;
		vec23.x /= length23; vec23.y /= length23; vec23.z /= length23;
		vec13.x /= length13; vec13.y /= length13; vec13.z /= length13;

		double innerprot = vec21.x*vec23.x+vec21.y*vec23.y+vec21.z*vec23.z;
		double dcosangle = innerprot;
		if     (dcosangle >  1.) dcosangle =  1.;
		else if(dcosangle < -1.) dcosangle = -1.;

		bool angle_active = (dcosangle >= idealahdang);
		if (!angle_active) dcosangle = -1.0;

		double idealahdiff = length23-idealahleng;
		double idealaddiff = length13-idealadleng;

		double denergydepth = -dK*exp(-(dXiR * (idealahdiff*idealahdiff + idealaddiff*idealaddiff)) - dXiTheta*dcosangle);

		//cout<<dK<<endl;
		return denergydepth;

	}else{
		return prevene;
	}
}

void HBond2MCeneobj::setprevene(const double &_prevene) const
{
	prevene=_prevene;
	*atomenergy2+=min(exp(_prevene),1000.0);
	*atomenergy3+=min(exp(_prevene),1000.0);
}

/*
void BasePairMCobj::CalcCenterOfMass() const 
{
	// It's not precisely the center of mass (because atoms have no masses) but rather geometrical center.
	Xm = 0., Ym = 0., Zm = 0.;
	for(unsigned int i = 0; i < n; ++i)
	{
		double Xk, Yk, Zk;
		Xk = AtomPtrs[i]->x;
		Yk = AtomPtrs[i]->y;
		Zk = AtomPtrs[i]->z;

		Xm += Xk;
		Ym += Yk;
		Zm += Zk;
	}
	Xm /= double(n);	// Mean coordinates
	Ym /= double(n);
	Zm /= double(n);  
}

double BasePairMCobj::TotalSQRDistance(double _a, double _b) const
{
	double SumDelta = 0.;
	for(unsigned int i = 0; i < n; ++i)
	{
		double xk, yk, zk;
		xk = AtomPtrs[i]->x;
		yk = AtomPtrs[i]->y;
		zk = AtomPtrs[i]->z;
		xk -= Xm;		// Reduced coords
		yk -= Ym;
		zk -= Zm;
		SumDelta += (_a*xk + _b*yk + zk)*(_a*xk + _b*yk + zk);
	}
	return SumDelta/(_a*_a + _b*_b + 1.);
}

//---------------------------------------------------------------------------------------------------
double BasePairMCobj::TotalSQRDistance() const
{
	double SumDelta = 0.;
	for(unsigned int i = 0; i < n; ++i)
	{
		double xk, yk, zk;
		xk = AtomPtrs[i]->x;
		yk = AtomPtrs[i]->y;
		zk = AtomPtrs[i]->z;
		xk -= Xm;		// Reduced coords
		yk -= Ym;
		zk -= Zm;
		SumDelta += (a*xk + b*yk + zk)*(a*xk + b*yk + zk);
	}
	return SumDelta/(a*a + b*b + 1.);
}


double BasePairMCobj::Fit3DPlane() const
{
	CalcCenterOfMass();
	
	double Sxx = 0., Syy = 0., Szz = 0., Sxy = 0., Sxz = 0., Syz = 0.;
	const unsigned int n = AtomPtrs.size();  
	for(unsigned int i = 0; i < n; ++i)
	{
		double xk, yk, zk;
		xk = AtomPtrs[i]->x;
		yk = AtomPtrs[i]->y;
		zk = AtomPtrs[i]->z;

		xk -= Xm;		// Reduced coords
		yk -= Ym;
		zk -= Zm;
		
		Sxx += xk*xk;
		Syy += yk*yk;
		Szz += zk*zk;
		Sxy += xk*yk;
		Sxz += xk*zk;
		Syz += yk*zk;
	}
	
	const double c0 = Syz*(Sxy*Sxy - Sxz*Sxz) + Sxy*Sxz*(Szz - Syy),
			c1 = Sxy*Sxy*Sxy + Sxy*(Sxz*Sxz - 2.*Syz*Syz - Szz*Szz) + Sxy*(Sxx*Szz + Syy*Szz - Sxx*Syy) + 
				Sxz*Syz*(Syy + Szz - 2.*Sxx),
			c2 = Syz*Syz*Syz + Syz*(Sxz*Sxz - 2.*Sxy*Sxy - Sxx*Sxx) + Syz*(Sxx*Szz + Sxx*Syy - Syy*Szz) +
				Sxy*Sxz*(Sxx + Syy - 2.*Szz),
			c3 = Sxy*(Syz*Syz - Sxz*Sxz) + Sxz*Syz*(Sxx - Syy);

	const double r = c2/c3,
			s = c1/c3,
			t = c0/c3;
	
	const double p = s - r*r/3.,
			q = 2.*r*r*r/27. - r*s/3. + t,
			R = q*q/4. + p*p*p/27.;

	double SumDelta;	// Sum over all AtomPtrs of squared distances to the plane
	
	if(R >= 0.)
	{
		a = -r/3. + pow(-q/2. + sqrt(R), 0.333333333) + pow(-q/2. - sqrt(R), 0.333333333);
		b = (Sxy*Syz*a*a + (Syz*Syz - Sxy*Sxy)*a - Sxy*Syz)/
		((Syz*(Sxx - Syy) - Sxy*Sxz)*a + Sxy*(Syy - Szz) + Sxz*Syz);
		SumDelta = TotalSQRDistance(a, b);
	}
	else
	{
		const double rho = sqrt(-p*p*p/27.),
			phi = acos(-q/(2.*rho));
		
		const double a1 = -r/3. + 2.*pow(rho, 0.333333333)*cos(phi/3.),
			a2 = -r/3. + 2.*pow(rho, 0.333333333)*cos((phi + 2.*PI_MACRO)/3.),
			a3 = -r/3. + 2.*pow(rho, 0.333333333)*cos((phi + 4.*PI_MACRO)/3.);
			
		const double b1 = (Sxy*Syz*a1*a1 + (Syz*Syz - Sxy*Sxy)*a1 - Sxy*Syz)/
				((Syz*(Sxx - Syy) - Sxy*Sxz)*a1 + Sxy*(Syy - Szz) + Sxz*Syz),
			b2 = (Sxy*Syz*a2*a2 + (Syz*Syz - Sxy*Sxy)*a2 - Sxy*Syz)/
				((Syz*(Sxx - Syy) - Sxy*Sxz)*a2 + Sxy*(Syy - Szz) + Sxz*Syz),
			b3 = (Sxy*Syz*a3*a3 + (Syz*Syz - Sxy*Sxy)*a3 - Sxy*Syz)/
				((Syz*(Sxx - Syy) - Sxy*Sxz)*a3 + Sxy*(Syy - Szz) + Sxz*Syz);
			
		const double SumDelta1 = TotalSQRDistance(a1, b1),
			SumDelta2 = TotalSQRDistance(a2, b2),
			SumDelta3 = TotalSQRDistance(a3, b3);
			
		a = a1, b = b1; SumDelta = SumDelta1;
		if(SumDelta2 < SumDelta1) {a = a2; b = b2; SumDelta = SumDelta2;}
		if(SumDelta3 < SumDelta1 && SumDelta3 < SumDelta2) {a = a3; b = b3; SumDelta = SumDelta3;}
	} 
	
	return SumDelta;  
}

BasePairMCobj::BasePairMCobj(const vector<Point3s*>& _AtomPtrs, vector<double*> &_atomenergies, int* _flag)
{
	if(!_AtomPtrs.empty())
	{
		AtomPtrs = _AtomPtrs;
		n = AtomPtrs.size();
		Fit3DPlane();  
	}  
	atomenergies=_atomenergies;
	flag=_flag;
}

double BasePairMCobj::calcenergybasepair() const
{
	if(*flag){
		Fit3DPlane();
  		double dEnergy = TotalSQRDistance();
  		return (isnan(dEnergy)) ? 0. : dEnergy;
	}else{
		return prevene;
	}
}

void BasePairMCobj::setprevene(const double &_prevene) const
{
	prevene=_prevene;

	for(int i=0;i<atomenergies.size();i++){
		*atomenergies[i] +=_prevene;
	}
}
*/

double BasePairMCobj::getdistfromplane(double a, double b, double c, double d, Point3s pointP) const
{
	double numerator = abs(a * pointP.x + b * pointP.y + c * pointP.z + d);
    double denominator = sqrt(a * a + b * b + c * c);
    
    return numerator / denominator;
}

void BasePairMCobj::getplanebp(Point3s pointA, Point3s pointB, Point3s pointC, double &a, double &b, double &c, double &d) const
{
	Point3s AB;
    Point3s AC;

	AB.x=pointB.x - pointA.x; AB.y=pointB.y - pointA.y; AB.z=pointB.z - pointA.z;
	AC.x=pointC.x - pointA.x; AC.y=pointC.y - pointA.y; AC.z=pointC.z - pointA.z;

	Point3s normal;
	normal.x = AB.y * AC.z - AB.z * AC.y;
	normal.y = AB.z * AC.x - AB.x * AC.z;
	normal.z = AB.x * AC.y - AB.y * AC.x;

	a = normal.x;
    b = normal.y;
    c = normal.z;

	d = -(a * pointA.x + b * pointA.y + c * pointA.z);
}

double BasePairMCobj::calcenergybasepair() const
{
	if(*flag){
		double a, b, c, d;
		double ene=0;
		getplanebp(*baseA1, *baseA2, *baseA3, a, b, c, d);
		for(int i=0;i<AtomPtrsB.size();i++) ene+=getdistfromplane(a, b, c, d, *AtomPtrsB[i]);

		getplanebp(*baseB1, *baseB2, *baseB3, a, b, c, d);
		for(int i=0;i<AtomPtrsA.size();i++) ene+=getdistfromplane(a, b, c, d, *AtomPtrsA[i]);

		return ene;
	}else{
		return prevene;
	}
}

void BasePairMCobj::setprevene(const double &_prevene) const
{
	prevene=_prevene;

	for(int i=0;i<atomenergiesA.size();i++){
		*atomenergiesA[i] +=_prevene;
	}

	for(int i=0;i<atomenergiesB.size();i++){
		*atomenergiesB[i] +=_prevene;
	}
}

BasePairMCobj2::BasePairMCobj2(const vector<Point3s*> _AtomPtrs, const vector<double*> _atomenergies, double _dK, int* _flag)
{
  if(!_AtomPtrs.empty())
  {
    AtomPtrs = _AtomPtrs;
	atomenergies = _atomenergies;
    n = AtomPtrs.size();
	dK = _dK;
	flag = _flag;
    Fit3DPlane();  
  }  
}

void BasePairMCobj2::CalcCenterOfMass()
{// It's not precisely the center of mass (because atoms have no masses) but rather geometrical center.
  Xm = 0., Ym = 0., Zm = 0.;
  for(unsigned int i = 0; i < n; ++i)
  {
    double Xk, Yk, Zk;
	Xk = AtomPtrs[i]->x;
	Yk = AtomPtrs[i]->y;
	Zk = AtomPtrs[i]->z;

    Xm += Xk;
    Ym += Yk;
    Zm += Zk;
  }
  Xm /= double(n);	// Mean coordinates
  Ym /= double(n);
  Zm /= double(n);  
}

double BasePairMCobj2::TotalSQRDistance(double _a, double _b) const
{
  double SumDelta = 0.;
  for(unsigned int i = 0; i < n; ++i)
  {
    double xk, yk, zk;
	xk = AtomPtrs[i]->x;
	yk = AtomPtrs[i]->y;
	zk = AtomPtrs[i]->z;
    xk -= Xm;		// Reduced coords
    yk -= Ym;
    zk -= Zm;
    SumDelta += (_a*xk + _b*yk + zk)*(_a*xk + _b*yk + zk);
  }
  return SumDelta/(_a*_a + _b*_b + 1.);
}

double BasePairMCobj2::TotalSQRDistance() const
{
  double SumDelta = 0.;
  for(unsigned int i = 0; i < n; ++i)
  {
    double xk, yk, zk;
	xk = AtomPtrs[i]->x;
	yk = AtomPtrs[i]->y;
	zk = AtomPtrs[i]->z;
    xk -= Xm;		// Reduced coords
    yk -= Ym;
    zk -= Zm;
    SumDelta += (a*xk + b*yk + zk)*(a*xk + b*yk + zk);
  }
  return SumDelta/(a*a + b*b + 1.);
}

double BasePairMCobj2::Fit3DPlane()
{/******************************************************************************************************
 * Returns the same as TotalSQRDistance(), but for optimal {a, b} (also calculates the optimal {a, b}
 * and stores them as this->a and this->b). It's a solution to 3rd degree polynomial.
 * Save your time and don't stare at this for too long.
 ******************************************************************************************************/
  CalcCenterOfMass();
  
  double Sxx = 0., Syy = 0., Szz = 0., Sxy = 0., Sxz = 0., Syz = 0.;
  const unsigned int n = AtomPtrs.size();  
  for(unsigned int i = 0; i < n; ++i)
  {
    double xk, yk, zk;
	xk = AtomPtrs[i]->x;
	yk = AtomPtrs[i]->y;
	zk = AtomPtrs[i]->z;

    xk -= Xm;		// Reduced coords
    yk -= Ym;
    zk -= Zm;
    
    Sxx += xk*xk;
    Syy += yk*yk;
    Szz += zk*zk;
    Sxy += xk*yk;
    Sxz += xk*zk;
    Syz += yk*zk;
  }
  
  const double c0 = Syz*(Sxy*Sxy - Sxz*Sxz) + Sxy*Sxz*(Szz - Syy),
	       c1 = Sxy*Sxy*Sxy + Sxy*(Sxz*Sxz - 2.*Syz*Syz - Szz*Szz) + Sxy*(Sxx*Szz + Syy*Szz - Sxx*Syy) + 
		    Sxz*Syz*(Syy + Szz - 2.*Sxx),
	       c2 = Syz*Syz*Syz + Syz*(Sxz*Sxz - 2.*Sxy*Sxy - Sxx*Sxx) + Syz*(Sxx*Szz + Sxx*Syy - Syy*Szz) +
		    Sxy*Sxz*(Sxx + Syy - 2.*Szz),
	       c3 = Sxy*(Syz*Syz - Sxz*Sxz) + Sxz*Syz*(Sxx - Syy);

  const double r = c2/c3,
	       s = c1/c3,
	       t = c0/c3;
  
  const double p = s - r*r/3.,
	       q = 2.*r*r*r/27. - r*s/3. + t,
	       R = q*q/4. + p*p*p/27.;

  double SumDelta;	// Sum over all AtomPtrs of squared distances to the plane
  
  if(R >= 0.)
  {
    a = -r/3. + pow(-q/2. + sqrt(R+1e-8), 0.333333333) + pow(-q/2. - sqrt(R+1e-8), 0.333333333);
    b = (Sxy*Syz*a*a + (Syz*Syz - Sxy*Sxy)*a - Sxy*Syz)/
	((Syz*(Sxx - Syy) - Sxy*Sxz)*a + Sxy*(Syy - Szz) + Sxz*Syz);
    SumDelta = TotalSQRDistance(a, b);
  }
  else
  {
    const double rho = sqrt((-p*p*p/27.)+1e-8),
		 phi = acos(-q/(2.*rho));
    
    const double a1 = -r/3. + 2.*pow(rho, 0.333333333)*cos(phi/3.),
		 a2 = -r/3. + 2.*pow(rho, 0.333333333)*cos((phi + 2.*PI_MACRO)/3.),
		 a3 = -r/3. + 2.*pow(rho, 0.333333333)*cos((phi + 4.*PI_MACRO)/3.);
		 
    const double b1 = (Sxy*Syz*a1*a1 + (Syz*Syz - Sxy*Sxy)*a1 - Sxy*Syz)/
		      ((Syz*(Sxx - Syy) - Sxy*Sxz)*a1 + Sxy*(Syy - Szz) + Sxz*Syz),
		 b2 = (Sxy*Syz*a2*a2 + (Syz*Syz - Sxy*Sxy)*a2 - Sxy*Syz)/
		      ((Syz*(Sxx - Syy) - Sxy*Sxz)*a2 + Sxy*(Syy - Szz) + Sxz*Syz),
		 b3 = (Sxy*Syz*a3*a3 + (Syz*Syz - Sxy*Sxy)*a3 - Sxy*Syz)/
		      ((Syz*(Sxx - Syy) - Sxy*Sxz)*a3 + Sxy*(Syy - Szz) + Sxz*Syz);
		 
    const double SumDelta1 = TotalSQRDistance(a1, b1),
		 SumDelta2 = TotalSQRDistance(a2, b2),
		 SumDelta3 = TotalSQRDistance(a3, b3);
		 
    a = a1, b = b1; SumDelta = SumDelta1;
    if(SumDelta2 < SumDelta1) {a = a2; b = b2; SumDelta = SumDelta2;}
    if(SumDelta3 < SumDelta1 && SumDelta3 < SumDelta2) {a = a3; b = b3; SumDelta = SumDelta3;}
  } 

/*// If you need the explicit equation for the plane: A*x + B*y + C*z + 1 = 0
  double lambda = 0.;
  for(unsigned int i = 0; i < n; ++i)
  {
    double Xk, Yk, Zk;
    AtomPtrs[i]->GetCoords(Xk, Yk, Zk);
    lambda += a*Xk + b*Yk + Zk;
  }
  lambda /= double(n);
  
  double A = -a/lambda,
	 B = -b/lambda,
	 C = -1./lambda;*/
  
  return SumDelta;  
}


double BasePairMCobj2::calcenergybasepair()
{
	if(*flag){
		Fit3DPlane();
		double ene = dK*TotalSQRDistance();

		return (isnan(ene)) ? 0. : ene;
	}else{
		return prevene;
	}
}

void BasePairMCobj2::setprevene(const double &_prevene) const
{
	prevene=_prevene;

	for(int i=0;i<atomenergies.size();i++){
		*atomenergies[i] +=_prevene;
	}
}

double BasePairMCobj3::calcenergy() const
{
	if(*flag){
		double newene;
		double deltaene;

		if(*bleng < *lower || *upper < *bleng){
			double tmp = *bleng - *ideal;
			
			//if(sqrt(tmp*tmp) > 1.0) cout<< remark << "-" << sqrt(tmp*tmp)<<endl;

			return sqrt(tmp*tmp);
		}else{
			return 0;
		}
	}else{
		return prevene;
	}
}

void BasePairMCobj3::setprevene(const double &_prevene) const
{
	prevene=_prevene;
	*atomenergy1+=_prevene;
	*atomenergy2+=_prevene;
}

double getstackenergy(const vector<Point3s*> AtomPtrsA, const vector<Point3s*> AtomPtrsB, const char &baseA, const char &baseB, const int &sttype){
	int i,j;
	int atom, atom1, atom2;
	int count;
	vector<vector<double>> RotMatix;
    vector<double> TranVect;
	int minidx, minorder;
	double minrmsd=1.0/0.0;
	double rmsd;
	vector<double> xyztmp(3,0);
	Point3s vec;
	BasicFunc bf;

	
    map<char, vector<int>> targetatoms={
		{'A',{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21}},
        {'C',{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19}},
        {'G',{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}},
        {'U',{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19}}
	};
	
	/*
	map<char, vector<int>> targetatoms={
		{'A',{12,13,14,15,16,17,18,19,20,21}},
        {'C',{12,13,14,15,16,17,18,19}},
        {'G',{12,13,14,15,16,17,18,19,20,21,22}},
        {'U',{12,13,14,15,16,17,18,19}}
	};
	*/
	
	map<char, vector<int>> targetatoms_dist={
	    {'A',{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21}},
        {'C',{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19}},
        {'G',{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}},
        {'U',{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19}}
    };
	
	/*
	static map<char, vector<int>> targetatoms={
            {'A',{12, 13, 14, 15, 16, 17, 18, 19, 21}},
            {'C',{12, 13, 14, 15, 16, 17}},
            {'G',{12, 13, 14, 15, 16, 17, 18, 19, 20}},
            {'U',{12, 13, 14, 15, 16, 17}}
    };
	*/
	
	/*
	static map<char, vector<int>> targetatoms={
		{'A',{12,13,14,15,16,17,18,19,20,21}},
        {'C',{12,13,14,15,16,17,18,19}},
        {'G',{12,13,14,15,16,17,18,19,20,21,22}},
        {'U',{12,13,14,15,16,17,18,19}}
	};
	*/

	/*
	vector<vector<double> > xyz_list1A(targetatoms[baseA].size(), vector<double>(3,0));
    vector<vector<double> > xyz_list2A(targetatoms[baseA].size(), vector<double>(3,0));
	vector<vector<double> > xyz_list3A(targetatoms[baseB].size(), vector<double>(3,0));

	vector<vector<double> > xyz_list1B(targetatoms[baseB].size(), vector<double>(3,0));
    vector<vector<double> > xyz_list2B(targetatoms[baseB].size(), vector<double>(3,0));
	vector<vector<double> > xyz_list3B(targetatoms[baseA].size(), vector<double>(3,0));

	for(i=0;i<STACKTEMPLATE[sttype][baseA][baseB].size();i++)
	{
		count = 0;
		for(j=0; j<targetatoms[baseA].size();j++){
			atom=targetatoms[baseA][j];
			xyz_list1A[count][0]=AtomPtrsA[atom]->x;
			xyz_list1A[count][1]=AtomPtrsA[atom]->y;
			xyz_list1A[count][2]=AtomPtrsA[atom]->z;

			xyz_list2A[count][0]=STACKTEMPLATE[sttype][baseA][baseB][i][0].ptv[atom].x;
			xyz_list2A[count][1]=STACKTEMPLATE[sttype][baseA][baseB][i][0].ptv[atom].y;
			xyz_list2A[count][2]=STACKTEMPLATE[sttype][baseA][baseB][i][0].ptv[atom].z;
			count+=1;
		}

		bf.RotateCoor(xyz_list1A,xyz_list2A,RotMatix,TranVect);

		rmsd=0;
		count=0;
		for(j=0; j<targetatoms[baseB].size();j++){
			atom=targetatoms[baseB][j];
			xyz_list3A[count][0]=STACKTEMPLATE[sttype][baseA][baseB][i][1].ptv[atom].x;
			xyz_list3A[count][1]=STACKTEMPLATE[sttype][baseA][baseB][i][1].ptv[atom].y;
			xyz_list3A[count][2]=STACKTEMPLATE[sttype][baseA][baseB][i][1].ptv[atom].z;

			bf.Transform(xyz_list3A[count],RotMatix,TranVect,xyztmp);
			vec.x=AtomPtrsB[atom]->x-xyztmp[0];
			vec.y=AtomPtrsB[atom]->y-xyztmp[1];
			vec.z=AtomPtrsB[atom]->z-xyztmp[2];

			rmsd+=bf.norm(vec);
			count+=1;
		}

		count = 0;
		for(j=0; j<targetatoms[baseB].size();j++){
			atom=targetatoms[baseB][j];
			xyz_list1B[count][0]=AtomPtrsB[atom]->x;
			xyz_list1B[count][1]=AtomPtrsB[atom]->y;
			xyz_list1B[count][2]=AtomPtrsB[atom]->z;

			xyz_list2B[count][0]=STACKTEMPLATE[sttype][baseA][baseB][i][1].ptv[atom].x;
			xyz_list2B[count][1]=STACKTEMPLATE[sttype][baseA][baseB][i][1].ptv[atom].y;
			xyz_list2B[count][2]=STACKTEMPLATE[sttype][baseA][baseB][i][1].ptv[atom].z;
			count+=1;
		}

		bf.RotateCoor(xyz_list1B,xyz_list2B,RotMatix,TranVect);

		count = 0;
		for(j=0; j<targetatoms[baseA].size();j++){
			atom=targetatoms[baseA][j];
			xyz_list3B[count][0]=STACKTEMPLATE[sttype][baseA][baseB][i][0].ptv[atom].x;
			xyz_list3B[count][1]=STACKTEMPLATE[sttype][baseA][baseB][i][0].ptv[atom].y;
			xyz_list3B[count][2]=STACKTEMPLATE[sttype][baseA][baseB][i][0].ptv[atom].z;

			bf.Transform(xyz_list3B[count],RotMatix,TranVect,xyztmp);
			vec.x=AtomPtrsA[atom]->x-xyztmp[0];
			vec.y=AtomPtrsA[atom]->y-xyztmp[1];
			vec.z=AtomPtrsA[atom]->z-xyztmp[2];

			rmsd+=bf.norm(vec);
			count+=1;
		}

		if(minrmsd > rmsd) minrmsd=rmsd;
	}

	if(sttype == 0 && baseA == baseB)
	{
		for(i=0;i<STACKTEMPLATE[sttype][baseA][baseB].size();i++)
		{
			count = 0;
			for(j=0; j<targetatoms[baseA].size();j++){
				atom=targetatoms[baseA][j];
				xyz_list1A[count][0]=AtomPtrsA[atom]->x;
				xyz_list1A[count][1]=AtomPtrsA[atom]->y;
				xyz_list1A[count][2]=AtomPtrsA[atom]->z;

				xyz_list2A[count][0]=STACKTEMPLATE[sttype][baseA][baseB][i][1].ptv[atom].x;
				xyz_list2A[count][1]=STACKTEMPLATE[sttype][baseA][baseB][i][1].ptv[atom].y;
				xyz_list2A[count][2]=STACKTEMPLATE[sttype][baseA][baseB][i][1].ptv[atom].z;
				count+=1;
			}

			bf.RotateCoor(xyz_list1A,xyz_list2A,RotMatix,TranVect);

			rmsd=0;
			count=0;
			for(j=0; j<targetatoms[baseB].size();j++){
				atom=targetatoms[baseB][j];
				xyz_list3A[count][0]=STACKTEMPLATE[sttype][baseA][baseB][i][0].ptv[atom].x;
				xyz_list3A[count][1]=STACKTEMPLATE[sttype][baseA][baseB][i][0].ptv[atom].y;
				xyz_list3A[count][2]=STACKTEMPLATE[sttype][baseA][baseB][i][0].ptv[atom].z;

				bf.Transform(xyz_list3A[count],RotMatix,TranVect,xyztmp);
				vec.x=AtomPtrsB[atom]->x-xyztmp[0];
				vec.y=AtomPtrsB[atom]->y-xyztmp[1];
				vec.z=AtomPtrsB[atom]->z-xyztmp[2];

				rmsd+=bf.norm(vec);
				count+=1;
			}

			count = 0;
			for(j=0; j<targetatoms[baseB].size();j++){
				atom=targetatoms[baseB][j];
				xyz_list1B[count][0]=AtomPtrsB[atom]->x;
				xyz_list1B[count][1]=AtomPtrsB[atom]->y;
				xyz_list1B[count][2]=AtomPtrsB[atom]->z;

				xyz_list2B[count][0]=STACKTEMPLATE[sttype][baseA][baseB][i][0].ptv[atom].x;
				xyz_list2B[count][1]=STACKTEMPLATE[sttype][baseA][baseB][i][0].ptv[atom].y;
				xyz_list2B[count][2]=STACKTEMPLATE[sttype][baseA][baseB][i][0].ptv[atom].z;
				count+=1;
			}

			bf.RotateCoor(xyz_list1B,xyz_list2B,RotMatix,TranVect);

			count = 0;
			for(j=0; j<targetatoms[baseA].size();j++){
				atom=targetatoms[baseA][j];
				xyz_list3B[count][0]=STACKTEMPLATE[sttype][baseA][baseB][i][1].ptv[atom].x;
				xyz_list3B[count][1]=STACKTEMPLATE[sttype][baseA][baseB][i][1].ptv[atom].y;
				xyz_list3B[count][2]=STACKTEMPLATE[sttype][baseA][baseB][i][1].ptv[atom].z;

				bf.Transform(xyz_list3B[count],RotMatix,TranVect,xyztmp);
				vec.x=AtomPtrsA[atom]->x-xyztmp[0];
				vec.y=AtomPtrsA[atom]->y-xyztmp[1];
				vec.z=AtomPtrsA[atom]->z-xyztmp[2];

				rmsd+=bf.norm(vec);
				count+=1;
			}

			
			if(minrmsd > rmsd){
				minrmsd=rmsd;
			}
		}
	}

	return minrmsd;
	*/



	vector<vector<double> > xyz_list1(targetatoms[baseA].size()+targetatoms[baseB].size(), vector<double>(3,0));
    vector<vector<double> > xyz_list2(targetatoms[baseA].size()+targetatoms[baseB].size(), vector<double>(3,0));

	for(i=0;i<STACKTEMPLATE[sttype][baseA][baseB].size();i++)
	{
		count = 0;
		for(j=0; j<targetatoms[baseA].size();j++){
			atom=targetatoms[baseA][j];
			xyz_list1[count][0]=AtomPtrsA[atom]->x;
			xyz_list1[count][1]=AtomPtrsA[atom]->y;
			xyz_list1[count][2]=AtomPtrsA[atom]->z;

			xyz_list2[count][0]=STACKTEMPLATE[sttype][baseA][baseB][i][0].ptv[atom].x;
			xyz_list2[count][1]=STACKTEMPLATE[sttype][baseA][baseB][i][0].ptv[atom].y;
			xyz_list2[count][2]=STACKTEMPLATE[sttype][baseA][baseB][i][0].ptv[atom].z;
			count+=1;
		}

		for(j=0; j<targetatoms[baseB].size();j++){
			atom=targetatoms[baseB][j];
			xyz_list1[count][0]=AtomPtrsB[atom]->x;
			xyz_list1[count][1]=AtomPtrsB[atom]->y;
			xyz_list1[count][2]=AtomPtrsB[atom]->z;

			xyz_list2[count][0]=STACKTEMPLATE[sttype][baseA][baseB][i][1].ptv[atom].x;
			xyz_list2[count][1]=STACKTEMPLATE[sttype][baseA][baseB][i][1].ptv[atom].y;
			xyz_list2[count][2]=STACKTEMPLATE[sttype][baseA][baseB][i][1].ptv[atom].z;
			count+=1;
		}

		bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

		rmsd=0;
		count=0;
		for(j=0; j<targetatoms[baseA].size();j++){
			atom=targetatoms[baseA][j];
			bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);
			vec.x=AtomPtrsA[atom]->x-xyztmp[0];
			vec.y=AtomPtrsA[atom]->y-xyztmp[1];
			vec.z=AtomPtrsA[atom]->z-xyztmp[2];

			rmsd+=bf.norm(vec);
			count+=1;
		}

		for(j=0; j<targetatoms[baseB].size();j++){
			atom=targetatoms[baseB][j];
			bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);

			vec.x=AtomPtrsB[atom]->x-xyztmp[0];
			vec.y=AtomPtrsB[atom]->y-xyztmp[1];
			vec.z=AtomPtrsB[atom]->z-xyztmp[2];

			rmsd+=bf.norm(vec);
			count+=1;
		}

		if(minrmsd > rmsd){
			minrmsd=rmsd;
			minidx=i;
			minorder=0;
		}
	}

	if(sttype == 0 && baseA == baseB)
	{
		for(i=0;i<STACKTEMPLATE[sttype][baseA][baseB].size();i++)
		{
			int count = 0;
			for(j=0; j<targetatoms[baseB].size();j++){
				atom=targetatoms[baseB][j];
				xyz_list1[count][0]=AtomPtrsB[atom]->x;
				xyz_list1[count][1]=AtomPtrsB[atom]->y;
				xyz_list1[count][2]=AtomPtrsB[atom]->z;

				xyz_list2[count][0]=STACKTEMPLATE[sttype][baseA][baseB][i][0].ptv[atom].x;
				xyz_list2[count][1]=STACKTEMPLATE[sttype][baseA][baseB][i][0].ptv[atom].y;
				xyz_list2[count][2]=STACKTEMPLATE[sttype][baseA][baseB][i][0].ptv[atom].z;
				count+=1;
			}

			for(j=0; j<targetatoms[baseA].size();j++){
				atom=targetatoms[baseA][j];
				xyz_list1[count][0]=AtomPtrsA[atom]->x;
				xyz_list1[count][1]=AtomPtrsA[atom]->y;
				xyz_list1[count][2]=AtomPtrsA[atom]->z;

				xyz_list2[count][0]=STACKTEMPLATE[sttype][baseA][baseB][i][1].ptv[atom].x;
				xyz_list2[count][1]=STACKTEMPLATE[sttype][baseA][baseB][i][1].ptv[atom].y;
				xyz_list2[count][2]=STACKTEMPLATE[sttype][baseA][baseB][i][1].ptv[atom].z;
				count+=1;
			}

			bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

			rmsd=0;
			count=0;
			//for(int j=0; j<SNTTEMPLATEATOMS[aat].size();j++){
			for(j=0; j<targetatoms[baseB].size();j++){
				atom=targetatoms[baseB][j];
				bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);
				vec.x=AtomPtrsB[atom]->x-xyztmp[0];
				vec.y=AtomPtrsB[atom]->y-xyztmp[1];
				vec.z=AtomPtrsB[atom]->z-xyztmp[2];
				rmsd+=bf.norm(vec);
				count+=1;
			}

			for(j=0; j<targetatoms[baseA].size();j++){
				atom=targetatoms[baseA][j];
				bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);

				vec.x=AtomPtrsA[atom]->x-xyztmp[0];
				vec.y=AtomPtrsA[atom]->y-xyztmp[1];
				vec.z=AtomPtrsA[atom]->z-xyztmp[2];
				rmsd+=bf.norm(vec);
				count+=1;
			}

			if(minrmsd > rmsd){
				minrmsd=rmsd;
				minidx=i;
				minorder=1;
			}
		}
	}

	//cout<<minrmsd<<endl;

	double energy=0;
	double leng, ileng;
	if(minorder == 0){
        for(int i=0;i<targetatoms_dist[baseA].size();i++){
            atom1=targetatoms_dist[baseA][i];
            for(int j=0;j<targetatoms_dist[baseB].size();j++){
                atom2=targetatoms_dist[baseB][j];

				double dx = AtomPtrsA[atom1]->x - AtomPtrsB[atom2]->x;
                double dy = AtomPtrsA[atom1]->y - AtomPtrsB[atom2]->y;
                double dz = AtomPtrsA[atom1]->z - AtomPtrsB[atom2]->z;
                leng = sqrt(dx*dx + dy*dy+ dz*dz+1e-12);
                
				ileng=STACKDISTMAP[sttype][baseA][baseB][minidx][atom1][atom2];

				double diff = leng - ileng;
				//cout<<atom1 << "-" << atom2 << "-" << sqrt(diff*diff) << "-" << leng << "-" << ileng<<endl;
				energy+=sqrt(diff*diff);
            }
        }
    }else{
        for(int i=0;i<targetatoms_dist[baseA].size();i++){
            atom1=targetatoms_dist[baseA][i];
            for(int j=0;j<targetatoms_dist[baseB].size();j++){
                atom2=targetatoms_dist[baseB][j];

				double dx = AtomPtrsA[atom1]->x - AtomPtrsB[atom2]->x;
                double dy = AtomPtrsA[atom1]->y - AtomPtrsB[atom2]->y;
                double dz = AtomPtrsA[atom1]->z - AtomPtrsB[atom2]->z;
                leng = sqrt(dx*dx + dy*dy+ dz*dz+1e-12);
                
				ileng=STACKDISTMAP[sttype][baseA][baseB][minidx][atom2][atom1];

				double diff = leng - ileng;
				//cout<<atom1 << "-" << atom2 << "-" << sqrt(diff*diff) << "-" << leng << "-" << ileng<<endl;
				energy+=sqrt(diff*diff);
            }
        }
    }

	//double energy=minrmsd;
	//if(energy > 80) energy=80;

	return energy;
}

double StackPairMCobj::calcenergystacking() const
{
	if(*flag){
		double a, b, c, d;
		//cout<<baseA << "-" << baseB << "-" << K << endl;
		double ene=K*getstackenergy(AtomPtrsA, AtomPtrsB, baseA, baseB, sttype);
	
		return ene;
	}else{
		return prevene;
	}
}

void StackPairMCobj::setprevene(const double &_prevene) const
{
	prevene=_prevene;

	for(int i=0;i<atomenergiesA.size();i++){
		*atomenergiesA[i] +=_prevene;
	}

	for(int i=0;i<atomenergiesB.size();i++){
		*atomenergiesB[i] +=_prevene;
	}
}

double getstackgeoenergy(const vector<Point3s*> AtomPtrsA, const vector<Point3s*> AtomPtrsB, const char &baseA, const char &baseB, const int &sttype){
	int i,j;
	int atom, atom1, atom2;
	int count;
	vector<vector<double>> RotMatix;
    vector<double> TranVect;
	int minidx, minorder;
	double minrmsd=1.0/0.0;
	double rmsd;
	vector<double> xyztmp(3,0);
	Point3s vec;
	BasicFunc bf;

	map<char, vector<int>> targetatoms={
		{'A',{12,13,14,15,16,17,18,19,20,21}},
        {'C',{12,13,14,15,16,17,18,19}},
        {'G',{12,13,14,15,16,17,18,19,20,21,22}},
        {'U',{12,13,14,15,16,17,18,19}}
	};

	double length_thresh = STACKPARAM[sttype][baseA][baseB][0][0]+3*STACKPARAM[sttype][baseA][baseB][0][1];
	//double length_lower = STACKPARAM[sttype][baseA][baseB][0][2];
	//double length_upper = STACKPARAM[sttype][baseA][baseB][0][3];
	double angle_thresh = cos((STACKPARAM[sttype][baseA][baseB][1][0]+3*STACKPARAM[sttype][baseA][baseB][1][1])*RADDEG);
	//double angle_lower = STACKPARAM[sttype][baseA][baseB][1][2];
	//double angle_upper = STACKPARAM[sttype][baseA][baseB][1][3];

	Point3s center_dA, center_dB;
	int nA, nB;

	center_dA.x=0.0; center_dA.y=0.0; center_dA.z=0.0; nA=0;
	for(int i=0;i<targetatoms[baseA].size();i++){
		atom=targetatoms[baseA][i];
		center_dA.x+=AtomPtrsA[atom]->x;
		center_dA.y+=AtomPtrsA[atom]->y;
		center_dA.z+=AtomPtrsA[atom]->z;
		nA+=1;
	}
	center_dA.x/=nA; center_dA.y/=nA; center_dA.z/=nA;

	center_dB.x=0.0; center_dB.y=0.0; center_dB.z=0.0; nB=0;
	for(int i=0;i<targetatoms[baseB].size();i++){
		atom=targetatoms[baseB][i];
		center_dB.x+=AtomPtrsB[atom]->x;
		center_dB.y+=AtomPtrsB[atom]->y;
		center_dB.z+=AtomPtrsB[atom]->z;
		nB+=1;
	}
	center_dB.x/=nB; center_dB.y/=nB; center_dB.z/=nB;

	double dx, dy, dz, leng;
	dx = center_dA.x - center_dB.x;
    dy = center_dA.y - center_dB.y;
    dz = center_dA.z - center_dB.z;
	leng = sqrt(dx*dx + dy*dy+ dz*dz+1e-12);

	//double energy1 = 0;
	//if(leng < length_lower | length_upper < leng){
		//double diff = leng - length_mean;
		//energy1 = diff*diff;
	//}

	double diff_leng;
	if (leng < length_thresh) diff_leng = 0.0;
	else diff_leng = leng - length_thresh;

	Point3s res1C42, res1C86, z1;
	Point3s res2C42, res2C86, z2;
	res1C42.x = AtomPtrsA[13]->x; res1C42.y = AtomPtrsA[13]->y; res1C42.z = AtomPtrsA[13]->z;
	res1C86.x = AtomPtrsA[14]->x; res1C86.y = AtomPtrsA[14]->y; res1C86.z = AtomPtrsA[14]->z;
	res2C42.x = AtomPtrsB[13]->x; res2C42.y = AtomPtrsB[13]->y; res2C42.z = AtomPtrsB[13]->z;
	res2C86.x = AtomPtrsB[14]->x; res2C86.y = AtomPtrsB[14]->y; res2C86.z = AtomPtrsB[14]->z;

	bf.getnormalvec(center_dA, res1C42, res1C86, z1);
	bf.getnormalvec(center_dB, res2C42, res2C86, z2);
	double costheta = bf.angleBetweenVectors(z1, z2);

	if(     costheta > 1. ) costheta = 1.;
  	else if(costheta < -1.) costheta = -1.;

	if(costheta < 0) costheta *= -1;

	double diff_ang;

	if (costheta > angle_thresh) diff_ang = 0.0;
	else diff_ang = costheta - angle_thresh;

	//cout<<diff_ang*diff_ang << "-" << diff_leng*diff_leng << "-" << -exp(-(diff_ang*diff_ang + diff_leng*diff_leng))<<endl;
 
	//if(baseangle > 90) baseangle = 180-baseangle;

	//double energy2 = 0;
	//if(baseangle < angle_lower | angle_upper < baseangle){
		//double diff = (baseangle - angle_mean) * (M_PI/180.0);
		//energy2 = diff*diff;
	//}

	//cout<<energy1<<" "<<energy2<<endl;

	//double energy=energy1+energy2;

	//return diff_ang*diff_ang + diff_leng*diff_leng;
	return -exp(-(diff_ang*diff_ang + diff_leng*diff_leng));
}


double StackPairGeoMCobj::calcenergystacking() const
{
	if(*flag){
		double a, b, c, d;
		double ene=K*getstackgeoenergy(AtomPtrsA, AtomPtrsB, baseA, baseB, sttype);
	
		return ene;
	}else{
		return prevene;
	}
}

void StackPairGeoMCobj::setprevene(const double &_prevene) const
{
	prevene=_prevene;

	for(int i=0;i<atomenergiesA.size();i++){
		*atomenergiesA[i] +=_prevene;
	}

	for(int i=0;i<atomenergiesB.size();i++){
		*atomenergiesB[i] +=_prevene;
	}
}

double calcenergyhydrogenbond(const Point3s &teAtomHyd, const Point3s &teAtomDon, const Point3s &teAtomAcc)
{
	double dREq=1.55;
	double dK=1.0;
	double dXiR=2.0;
	double dXiTheta=2.0;
	double dThreshold=-1.0;
	Point3s vec21, vec23, vec13, grad1, grad2, grad3, gradsum;

	vec21.x = teAtomDon.x - teAtomHyd.x;
	vec21.y = teAtomDon.y - teAtomHyd.y;
	vec21.z = teAtomDon.z - teAtomHyd.z;

	vec23.x = teAtomAcc.x - teAtomHyd.x;
	vec23.y = teAtomAcc.y - teAtomHyd.y;
	vec23.z = teAtomAcc.z - teAtomHyd.z;

	double length21=sqrt(vec21.x*vec21.x+vec21.y*vec21.y+vec21.z*vec21.z);
	double length23=sqrt(vec23.x*vec23.x+vec23.y*vec23.y+vec23.z*vec23.z);

	vec21.x /= length21; vec21.y /= length21; vec21.z /= length21;
	vec23.x /= length23; vec23.y /= length23; vec23.z /= length23;
	
	double innerprot = vec21.x*vec23.x+vec21.y*vec23.y+vec21.z*vec23.z;

	double dcosangle = innerprot;
	if     (dcosangle >  1.) dcosangle =  1.;
  	else if(dcosangle < -1.) dcosangle = -1.;

	return -dK*exp(-dXiR*(length23-dREq)*(length23-dREq) - dXiTheta*dcosangle);
}

bool MonteCarloSimulation::checkhbond2(Point3s* teAtomHyd1, Point3s* teAtomHyd2, Point3s* teAtomDon, Point3s* teAtomAcc)
{
	Point3s vec23_1, vec23_2, vec23, vec13, vec21;

	vec23_1.x = teAtomAcc->x - teAtomHyd1->x;
	vec23_1.y = teAtomAcc->y - teAtomHyd1->y;
	vec23_1.z = teAtomAcc->z - teAtomHyd1->z;

	vec23_2.x = teAtomAcc->x - teAtomHyd2->x;
	vec23_2.y = teAtomAcc->y - teAtomHyd2->y;
	vec23_2.z = teAtomAcc->z - teAtomHyd2->z;

	double length23_1=sqrt(vec23_1.x*vec23_1.x+vec23_1.y*vec23_1.y+vec23_1.z*vec23_1.z+1e-8);
	double length23_2=sqrt(vec23_2.x*vec23_2.x+vec23_2.y*vec23_2.y+vec23_2.z*vec23_2.z+1e-8);

	int hydatom;
	if(length23_1 < length23_2){
		vec21.x = teAtomDon->x - teAtomHyd1->x;
		vec21.y = teAtomDon->y - teAtomHyd1->y;
		vec21.z = teAtomDon->z - teAtomHyd1->z;

		vec23.x = vec23_1.x; vec23.y = vec23_1.y; vec23.z = vec23_1.z;

		hydatom = -1;
	}else{
		vec21.x = teAtomDon->x - teAtomHyd2->x;
		vec21.y = teAtomDon->y - teAtomHyd2->y;
		vec21.z = teAtomDon->z - teAtomHyd2->z;

		vec23.x = vec23_2.x; vec23.y = vec23_2.y; vec23.z = vec23_2.z;

		hydatom = 1;
	}

	vec13.x = teAtomDon->x - teAtomAcc->x;
	vec13.y = teAtomDon->y - teAtomAcc->y;
	vec13.z = teAtomDon->z - teAtomAcc->z;

	double length21=sqrt(vec21.x*vec21.x+vec21.y*vec21.y+vec21.z*vec21.z+1e-8);
	double length23=sqrt(vec23.x*vec23.x+vec23.y*vec23.y+vec23.z*vec23.z+1e-8);
	double length13=sqrt(vec13.x*vec13.x+vec13.y*vec13.y+vec13.z*vec13.z+1e-8);

	vec21.x /= length21; vec21.y /= length21; vec21.z /= length21;
	vec23.x /= length23; vec23.y /= length23; vec23.z /= length23;
	vec13.x /= length13; vec13.y /= length13; vec13.z /= length13;

	double innerprot = vec21.x*vec23.x+vec21.y*vec23.y+vec21.z*vec23.z;
	double dcosangle = innerprot;
	if     (dcosangle >  1.) dcosangle =  1.;
	else if(dcosangle < -1.) dcosangle = -1.;

	double angle = acos(dcosangle)*DEGRAD;

	//cout<<angle << "-" << length13 << "-" << (length13 <= 4.5)<<endl;

	//return ((angle >= 95) && (length13 <= 4.5));
	return (length13 <= 4.0);
}

bool MonteCarloSimulation::checkhbond(Point3s* teAtomHyd, Point3s* teAtomDon, Point3s* teAtomAcc)
{
	Point3s vec21, vec23, vec13;

	vec21.x = teAtomDon->x - teAtomHyd->x;
	vec21.y = teAtomDon->y - teAtomHyd->y;
	vec21.z = teAtomDon->z - teAtomHyd->z;

	vec23.x = teAtomAcc->x - teAtomHyd->x;
	vec23.y = teAtomAcc->y - teAtomHyd->y;
	vec23.z = teAtomAcc->z - teAtomHyd->z;

	vec13.x = teAtomDon->x - teAtomAcc->x;
	vec13.y = teAtomDon->y - teAtomAcc->y;
	vec13.z = teAtomDon->z - teAtomAcc->z;

	double length21=sqrt(vec21.x*vec21.x+vec21.y*vec21.y+vec21.z*vec21.z+1e-8);
	double length23=sqrt(vec23.x*vec23.x+vec23.y*vec23.y+vec23.z*vec23.z+1e-8);
	double length13=sqrt(vec13.x*vec13.x+vec13.y*vec13.y+vec13.z*vec13.z+1e-8);

	vec21.x /= length21; vec21.y /= length21; vec21.z /= length21;
	vec23.x /= length23; vec23.y /= length23; vec23.z /= length23;
	
	double innerprot = vec21.x*vec23.x+vec21.y*vec23.y+vec21.z*vec23.z;
	double dcosangle = innerprot;
	if     (dcosangle >  1.) dcosangle =  1.;
	else if(dcosangle < -1.) dcosangle = -1.;

	double angle = acos(dcosangle)*DEGRAD;

	//cout<<angle << "-" << length13 << "-" << (length13 <= 4.5)<<endl;

	//return ((angle >= 95) && (length13 <= 4.5));
	return (length13 <= 4.0);
}

void MonteCarloSimulation::setpairwiseeneobj(const ExtParams &extparameters, const FragmentStructure &inputtedstrall, const vector<vector<bool>> &inputtedstrallflag, const vector<FragmentStructure> &templatestrall, const vector<vector<vector<bool>>> &templatestrallflag){
	int inda,indb;
	char base1, base2;
	int idx1, idx2;
	string base1_str, base2_str;
	double vdW, dQQ, dQQep;
	int k,l,m;
	int count=0;
	int count_st=0;
	char AccBase, DonBase;
	Point3d p12;
	BasicFunc bf;
	int sttype;
	double adlength;
	bool flag;

	vector<vector<int>> molConnect;
	map<string, vector<vector<bool>>> intraAtomConnect;

	for(int i=0;i<numbb;i++) looplist.push_back(i);

	bpscoresum=0;
	int bpcount=0;
	for(int i=0;i<numbb;i++) for(int j=i+1;j<numbb;j++) if(cbpmap[i][j] > 0) bpcount+=1;
	for(int i=0;i<numbb;i++) for(int j=i+1;j<numbb;j++) if(noncbpmap[i][j] > 0) bpcount+=1;

	if(bpcount != 0){
		BasePtrsA.resize(bpcount);
		BasePtrsB.resize(bpcount);
		BasePtrs.resize(2*bpcount);
		BaseAtomenergiesA.resize(bpcount);
		BaseAtomenergiesB.resize(bpcount);
		BaseAtomenergies.resize(2*bpcount);
	}

	int stackcount=0;
	for(int i=0;i<numbb;i++) for(int j=i+1;j<numbb;j++) if(stmap[i][j] > 0) stackcount+=1;
	if(stackcount != 0){
		StackPtrsA.resize(stackcount);
		StackPtrsB.resize(stackcount);
		StackAtomenergiesA.resize(stackcount);
		StackAtomenergiesB.resize(stackcount);
	}

	clashmcene.clear();
	electrosene.clear();
	for(int i=0;i<numbb;i++)
	{
		base1=refinestr_temp[i].aaa;
		base1_str=string(1, base1);
		for(int j=i;j<numbb;j++)
		{
			base2=refinestr_temp[j].aaa;
			base2_str=string(1, base2);
			getmolConnect(i, j, molConnect, intraAtomConnect, extparameters);

			if(i != j){
				if(stmap[i][j] > 0){
					stackscoresum+=stmap[i][j];
					for(int ii=0;ii<AccAtomParamNum[base1][3];ii++){
						StackPtrsA[count_st].push_back(&refinestr_temp[i].ptv[ii]);
					}
					for(int jj=0;jj<AccAtomParamNum[base2][3];jj++){
						StackPtrsB[count_st].push_back(&refinestr_temp[j].ptv[jj]);
					}
					
					for(int ii=0;ii<AccAtomParamNum[base1][3];ii++){
						StackAtomenergiesA[count_st].push_back(&pairwiseenergyatoms[i][ii]);
					}
					for(int jj=0;jj<AccAtomParamNum[base2][3];jj++){
						StackAtomenergiesB[count_st].push_back(&pairwiseenergyatoms[j][jj]);
					}

					if(i+1==j) sttype=1;
					else sttype=0;

					StackPairMCobj stackmc = StackPairMCobj(StackPtrsA[count_st], StackPtrsB[count_st], base1, base2, sttype, StackAtomenergiesA[count_st], StackAtomenergiesB[count_st], &mapflag[i][j], stmap[i][j]);
					stackene.push_back(stackmc);

					StackPairGeoMCobj stackgeomc = StackPairGeoMCobj(StackPtrsA[count_st], StackPtrsB[count_st], base1, base2, sttype, StackAtomenergiesA[count_st], StackAtomenergiesB[count_st], &mapflag[i][j], stmap[i][j]);
					stackgeoene.push_back(stackgeomc);



					count_st+=1;
					
					stacklist.push_back({i,j,sttype});
				}

				if(cbpmap[i][j] > 0)
				{
					bpscoresum+=cbpmap[i][j];
				//if(bpmap[i][j] == 2){
					/*
					for(int ii=AccAtomParamNum[base1][1];ii<AccAtomParamNum[base1][3];ii++){
						BasePtrsA[count].push_back(&refinestr_temp[i].ptv[ii]);
					}
					for(int jj=AccAtomParamNum[base2][1];jj<AccAtomParamNum[base2][3];jj++){
						BasePtrsB[count].push_back(&refinestr_temp[j].ptv[jj]);
					}
					
					for(int ii=AccAtomParamNum[base1][1];ii<AccAtomParamNum[base1][3];ii++){
						BaseAtomenergiesA[count].push_back(&pairwiseenergyatoms[i][ii]);
					}
					for(int jj=AccAtomParamNum[base2][1];jj<AccAtomParamNum[base2][3];jj++){
						BaseAtomenergiesB[count].push_back(&pairwiseenergyatoms[j][jj]);
					}
					
					BasePairMCobj bpmc = BasePairMCobj(BasePtrsA[count], BasePtrsB[count], &refinestr_temp[i].ptv[12], &refinestr_temp[i].ptv[13], &refinestr_temp[i].ptv[14], &refinestr_temp[j].ptv[12], &refinestr_temp[j].ptv[13], &refinestr_temp[j].ptv[14], BaseAtomenergiesA[count], BaseAtomenergiesB[count], &mapflag[i][j]);
					bpene.push_back(bpmc);
					*/

                    //if((base1_str == "G" && base2_str == "C") || (base1_str == "C" && base2_str == "G") || (base1_str == "A" && base2_str == "U") || (base1_str == "U" && base2_str == "A")){
						//cout<< i << "__" << j<<endl;
					for(int ii=AccAtomParamNum[base1][1];ii<AccAtomParamNum[base1][3];ii++){
						//cout<<ii<<endl;
						BasePtrs[count].push_back(&refinestr_temp[i].ptv[ii]);
					}
					for(int jj=AccAtomParamNum[base2][1];jj<AccAtomParamNum[base2][3];jj++){
						//cout<<jj<<endl;
						BasePtrs[count].push_back(&refinestr_temp[j].ptv[jj]);
					}

					for(int ii=AccAtomParamNum[base1][1];ii<AccAtomParamNum[base1][3];ii++){
						BaseAtomenergies[count].push_back(&pairwiseenergyatoms[i][ii]);
					}
					for(int jj=AccAtomParamNum[base2][1];jj<AccAtomParamNum[base2][3];jj++){
						BaseAtomenergies[count].push_back(&pairwiseenergyatoms[j][jj]);
					}

					//cout<<"aa-" << i << "-" << j << "-" << cbpmap[i][j]<<endl;
					BasePairMCobj2 bpmc = BasePairMCobj2(BasePtrs[count], BaseAtomenergies[count], cbpmap[i][j], &mapflag[i][j]);
					bpene2.push_back(bpmc);
					
					count+=1;
					//}

					//if((base1_str == "G" && base2_str == "C") || (base1_str == "C" && base2_str == "G") || (base1_str == "A" && base2_str == "U") || (base1_str == "U" && base2_str == "A") || (base1_str == "G" && base2_str == "U") || (base1_str == "U" && base2_str == "G"))
					//{
					vector<vector<vector<int>>> hbondinfo = BpHbondInfo[base1_str+base2_str];
					for(k=0;k<hbondinfo.size();k++)
					{
						if(hbondinfo[k].size() > 1)
						{
							if(hbondinfo[k][0][3] == 0){
								const Point3s &ptH1 = refinestr_temp[j].ptv[hbondinfo[k][0][1]];
								const Point3s &ptH2 = refinestr_temp[j].ptv[hbondinfo[k][1][1]];
								const Point3s &ptAcc = refinestr_temp[i].ptv[hbondinfo[k][0][0]];
								const Point3s &ptDon = refinestr_temp[j].ptv[hbondinfo[k][0][2]];

								AccBase = base1;
								DonBase = base2;
								
								HBond2MCeneobj hbondmc = HBond2MCeneobj::create(&refinestr_temp[j].ptv[hbondinfo[k][0][1]], &refinestr_temp[j].ptv[hbondinfo[k][1][1]], &refinestr_temp[j].ptv[hbondinfo[k][0][2]], &refinestr_temp[i].ptv[hbondinfo[k][0][0]], &pairwiseenergyatoms[j][hbondinfo[k][0][2]], &pairwiseenergyatoms[i][hbondinfo[k][0][0]], BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean, cbpmap[i][j], &mapflag[i][j]);
								hbond2ene.push_back(hbondmc);
							}else{
								const Point3s &ptH1 = refinestr_temp[i].ptv[hbondinfo[k][0][1]];
								const Point3s &ptH2 = refinestr_temp[i].ptv[hbondinfo[k][1][1]];
								const Point3s &ptAcc = refinestr_temp[j].ptv[hbondinfo[k][0][0]];
								const Point3s &ptDon = refinestr_temp[i].ptv[hbondinfo[k][0][2]];

								AccBase = base2;
								DonBase = base1;

								HBond2MCeneobj hbondmc = HBond2MCeneobj::create(&refinestr_temp[i].ptv[hbondinfo[k][0][1]], &refinestr_temp[i].ptv[hbondinfo[k][1][1]], &refinestr_temp[i].ptv[hbondinfo[k][0][2]], &refinestr_temp[j].ptv[hbondinfo[k][0][0]], &pairwiseenergyatoms[i][hbondinfo[k][0][2]], &pairwiseenergyatoms[j][hbondinfo[k][0][0]], BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean, cbpmap[i][j], &mapflag[i][j]);
								hbond2ene.push_back(hbondmc);
							}
						}else{

							if(hbondinfo[k][0][3] == 0){
								const Point3s &ptH = refinestr_temp[j].ptv[hbondinfo[k][0][1]];
								const Point3s &ptAcc = refinestr_temp[i].ptv[hbondinfo[k][0][0]];
								const Point3s &ptDon = refinestr_temp[j].ptv[hbondinfo[k][0][2]];

								AccBase = base1;
								DonBase = base2;
								HBondMCeneobj hbondmc = HBondMCeneobj::create(&refinestr_temp[j].ptv[hbondinfo[k][0][1]], &refinestr_temp[j].ptv[hbondinfo[k][0][2]], &refinestr_temp[i].ptv[hbondinfo[k][0][0]], &pairwiseenergyatoms[j][hbondinfo[k][0][2]], &pairwiseenergyatoms[i][hbondinfo[k][0][0]], BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean, cbpmap[i][j], &mapflag[i][j]);
								hbondene.push_back(hbondmc);
							}else{
								const Point3s &ptH = refinestr_temp[i].ptv[hbondinfo[k][0][1]];
								const Point3s &ptAcc = refinestr_temp[j].ptv[hbondinfo[k][0][0]];
								const Point3s &ptDon = refinestr_temp[i].ptv[hbondinfo[k][0][2]];

								AccBase = base2;
								DonBase = base1;
								HBondMCeneobj hbondmc = HBondMCeneobj::create(&refinestr_temp[i].ptv[hbondinfo[k][0][1]], &refinestr_temp[i].ptv[hbondinfo[k][0][2]], &refinestr_temp[j].ptv[hbondinfo[k][0][0]], &pairwiseenergyatoms[i][hbondinfo[k][0][2]], &pairwiseenergyatoms[j][hbondinfo[k][0][0]], BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean, cbpmap[i][j], &mapflag[i][j]);
								hbondene.push_back(hbondmc);
							}

						}
					}

					//if((base1_str == "G" && base2_str == "C") || (base1_str == "C" && base2_str == "G") || (base1_str == "A" && base2_str == "U") || (base1_str == "U" && base2_str == "A") || (base1_str == "G" && base2_str == "U") || (base1_str == "U" && base2_str == "G"))
					//{
					bplist.push_back({i,j});
					//}

					looplist.erase(std::remove(looplist.begin(), looplist.end(), i), looplist.end());
					looplist.erase(std::remove(looplist.begin(), looplist.end(), j), looplist.end());

				}

				if(noncbpmap[i][j] > 0)
				{
					for(int ii=AccAtomParamNum[base1][1];ii<AccAtomParamNum[base1][3];ii++){
						BasePtrs[count].push_back(&refinestr_temp[i].ptv[ii]);
					}
					for(int jj=AccAtomParamNum[base2][1];jj<AccAtomParamNum[base2][3];jj++){
						BasePtrs[count].push_back(&refinestr_temp[j].ptv[jj]);
					}

					for(int ii=AccAtomParamNum[base1][1];ii<AccAtomParamNum[base1][3];ii++){
						BaseAtomenergies[count].push_back(&pairwiseenergyatoms[i][ii]);
					}
					for(int jj=AccAtomParamNum[base2][1];jj<AccAtomParamNum[base2][3];jj++){
						BaseAtomenergies[count].push_back(&pairwiseenergyatoms[j][jj]);
					}

					//cout<<"bb-" << i << "-" << j << "-" << noncbpmap[i][j]<<endl;
					BasePairMCobj2 bpmc = BasePairMCobj2(BasePtrs[count], BaseAtomenergies[count], noncbpmap[i][j], &mapflag[i][j]);
					bpene2.push_back(bpmc);
					
					count+=1;

					vector<vector<vector<int>>> hbondinfo = BpHbondInfo[base1_str+base2_str];
					for(k=0;k<hbondinfo.size();k++)
					{
						if(hbondinfo[k].size() > 1)
						{
							if(hbondinfo[k][0][3] == 0){
								AccBase = base1;
								DonBase = base2;

								flag = checkhbond2(&refinestr_temp[j].ptv[hbondinfo[k][0][1]], &refinestr_temp[j].ptv[hbondinfo[k][1][1]], &refinestr_temp[j].ptv[hbondinfo[k][0][2]], &refinestr_temp[i].ptv[hbondinfo[k][0][0]]);
								
								if(!flag) continue;

								HBond2MCeneobj hbondmc = HBond2MCeneobj::create(&refinestr_temp[j].ptv[hbondinfo[k][0][1]], &refinestr_temp[j].ptv[hbondinfo[k][1][1]], &refinestr_temp[j].ptv[hbondinfo[k][0][2]], &refinestr_temp[i].ptv[hbondinfo[k][0][0]], &pairwiseenergyatoms[j][hbondinfo[k][0][2]], &pairwiseenergyatoms[i][hbondinfo[k][0][0]], BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean, noncbpmap[i][j], &mapflag[i][j]);
								hbond2ene.push_back(hbondmc);
							}else{
								AccBase = base2;
								DonBase = base1;

								flag = checkhbond2(&refinestr_temp[i].ptv[hbondinfo[k][0][1]], &refinestr_temp[i].ptv[hbondinfo[k][1][1]], &refinestr_temp[i].ptv[hbondinfo[k][0][2]], &refinestr_temp[j].ptv[hbondinfo[k][0][0]]);

								if(!flag) continue;

								HBond2MCeneobj hbondmc = HBond2MCeneobj::create(&refinestr_temp[i].ptv[hbondinfo[k][0][1]], &refinestr_temp[i].ptv[hbondinfo[k][1][1]], &refinestr_temp[i].ptv[hbondinfo[k][0][2]], &refinestr_temp[j].ptv[hbondinfo[k][0][0]], &pairwiseenergyatoms[i][hbondinfo[k][0][2]], &pairwiseenergyatoms[j][hbondinfo[k][0][0]], BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean, noncbpmap[i][j], &mapflag[i][j]);
								hbond2ene.push_back(hbondmc);
							}
						}else{

							if(hbondinfo[k][0][3] == 0){
								AccBase = base1;
								DonBase = base2;

								flag = checkhbond(&refinestr_temp[j].ptv[hbondinfo[k][0][1]], &refinestr_temp[j].ptv[hbondinfo[k][0][2]], &refinestr_temp[i].ptv[hbondinfo[k][0][0]]);

								if(!flag) continue;

								HBondMCeneobj hbondmc = HBondMCeneobj::create(&refinestr_temp[j].ptv[hbondinfo[k][0][1]], &refinestr_temp[j].ptv[hbondinfo[k][0][2]], &refinestr_temp[i].ptv[hbondinfo[k][0][0]], &pairwiseenergyatoms[j][hbondinfo[k][0][2]], &pairwiseenergyatoms[i][hbondinfo[k][0][0]], BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean, noncbpmap[i][j], &mapflag[i][j]);
								hbondene.push_back(hbondmc);
							}else{
								AccBase = base2;
								DonBase = base1;

								flag = checkhbond(&refinestr_temp[i].ptv[hbondinfo[k][0][1]], &refinestr_temp[i].ptv[hbondinfo[k][0][2]], &refinestr_temp[j].ptv[hbondinfo[k][0][0]]);

								if(!flag) continue;

								HBondMCeneobj hbondmc = HBondMCeneobj::create(&refinestr_temp[i].ptv[hbondinfo[k][0][1]], &refinestr_temp[i].ptv[hbondinfo[k][0][2]], &refinestr_temp[j].ptv[hbondinfo[k][0][0]], &pairwiseenergyatoms[i][hbondinfo[k][0][2]], &pairwiseenergyatoms[j][hbondinfo[k][0][0]], BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean, noncbpmap[i][j], &mapflag[i][j]);
								hbondene.push_back(hbondmc);
							}
						}
					}

					bplist.push_back({i,j});
					looplist.erase(std::remove(looplist.begin(), looplist.end(), i), looplist.end());
					looplist.erase(std::remove(looplist.begin(), looplist.end(), j), looplist.end());
				}


				for(int ii=0;ii<atomid2mcidx[base1].size();ii++)
				{
					for(int jj=0;jj<atomid2mcidx[base2].size();jj++)
					{
						if(molConnect[ii][jj] >= 2){
							k = mcidx2atomid[base1][ii];
							l = mcidx2atomid[base2][jj];
							
							//double vdW2 = extparameters.atomfeatures.at(base1_str).at(k)[1] + extparameters.atomfeatures.at(base2_str).at(l)[1];
							//double vdW3 = MCVWPARAM[base1][ii]+MCVWPARAM[base2][jj]-0.4;
							vdW = MCVWPARAM[base1][ii]+MCVWPARAM[base2][jj]-0.4;

							//cout<<k << "-" << l << "-" << vdW << "-" << vdW2 << "-" << vdW3<<endl;
							
							ClashMCeneobj clmc = ClashMCeneobj::create(&refinestr_temp[i].ptv[ii], &refinestr_temp[j].ptv[jj], &pairwiseenergyatoms[i][ii], &pairwiseenergyatoms[j][jj], &mapflag[i][j], vdW);
							clashmcene.push_back(clmc);
						}

						if(molConnect[ii][jj] >= 3){
							if(molConnect[ii][jj] == 3){
								dQQ = (COULOMB*extparameters.atomfeatures.at(base1_str).at(k)[0]*extparameters.atomfeatures.at(base2_str).at(l)[0])/1.2;
							}else{
								dQQ = COULOMB*extparameters.atomfeatures.at(base1_str).at(k)[0]*extparameters.atomfeatures.at(base2_str).at(l)[0];
							}
							dQQep = dQQ/EpsilonRNA;

							ElectrosMCeneobj elmc = ElectrosMCeneobj::create(&refinestr_temp[i].ptv[ii], &refinestr_temp[j].ptv[jj], &pairwiseenergyatoms[i][ii], &pairwiseenergyatoms[j][jj], &mapflag[i][j], dQQep);
							electrosene.push_back(elmc);
						}

						const Point3s &pt1 = refinestr_temp[i].ptv[ii];
						const Point3s &pt2 = refinestr_temp[j].ptv[jj];

						Point3s vec12;
						vec12.x = pt2.x - pt1.x;
						vec12.y = pt2.y - pt1.y;
						vec12.z = pt2.z - pt1.z;

						double length12=vec12.x*vec12.x+vec12.y*vec12.y+vec12.z*vec12.z;

						if(length12 < CUTOFF_SQR+10){
							double vdWInfodR01 = extparameters.atomfeatures.at(base1_str).at(k)[1];
							double vdWInfodR02 = extparameters.atomfeatures.at(base2_str).at(l)[1];
							double vdWInfodEpsilon1 = extparameters.atomfeatures.at(base1_str).at(k)[2];
							double vdWInfodEpsilon2 = extparameters.atomfeatures.at(base2_str).at(l)[2];
							double star = (vdWInfodR01 + vdWInfodR02);
							double epsilon =  sqrt(vdWInfodEpsilon1*vdWInfodEpsilon2);
							double LJ6 = -2*epsilon*pow(star,6);
							double LJ12 = epsilon*pow(star,12);

							//if(molConnect[ii][jj] == 3){
								//LJ6 /=2.0;
								//LJ12 /=2.0;
							//}

							LJMCeneobj ljmc = LJMCeneobj::create(&refinestr_temp[i].ptv[ii], &refinestr_temp[j].ptv[jj], &pairwiseenergyatoms[i][ii], &pairwiseenergyatoms[j][jj], &mapflag[i][j], LJ6, LJ12);
							ljene.push_back(ljmc);
						}
					
						/*
						if((HBondHydrogenAtoms[base1_str][k] == 1) && (HBondAcceptorAtoms[base2_str][l] == 1) && (molConnect[ii][jj] == 4)){
							for(m=0;m<atomid2mcidx[base1].size();m++){
								if(intraAtomConnect[base1_str][ii][m]){
									const Point3s &ptH = refinestr_temp[i].ptv[ii];
									const Point3s &ptAcc = refinestr_temp[j].ptv[jj];
									const Point3s &ptDon = refinestr_temp[i].ptv[m];

									p12=bf.setv(ptAcc.x-ptDon.x,ptAcc.y-ptDon.y,ptAcc.z-ptDon.z);
									double hleng=bf.norm(p12);

									//double energy=calcenergyhydrogenbond(ptH, ptDon, ptAcc);
									//if(energy < HBONDThreshold){
									if(hleng < HBONDADLengThreshold){
										HBondMCeneobj hbondmc = HBondMCeneobj::create(&refinestr_temp[i].ptv[ii], &refinestr_temp[i].ptv[m], &refinestr_temp[j].ptv[jj], &pairwiseenergyatoms[i][m], &pairwiseenergyatoms[j][jj], &mapflag[i][j]);
										hbondene.push_back(hbondmc);
									}
								}
							}
						}

						
						if((HBondHydrogenAtoms[base2_str][l] == 1) && (HBondAcceptorAtoms[base1_str][k] == 1) && (molConnect[ii][jj] == 4)){
							for(m=0;m<atomid2mcidx[base2].size();m++){
								if(intraAtomConnect[base2_str][jj][m]){
									const Point3s &ptH = refinestr_temp[j].ptv[jj];
									const Point3s &ptAcc = refinestr_temp[i].ptv[ii];
									const Point3s &ptDon = refinestr_temp[j].ptv[m];

									//double energy=calcenergyhydrogenbond(ptH, ptDon, ptAcc);
									p12=bf.setv(ptAcc.x-ptDon.x,ptAcc.y-ptDon.y,ptAcc.z-ptDon.z);
									double hleng=bf.norm(p12);

									//if(energy < HBONDThreshold){
									if(hleng < HBONDADLengThreshold){
										HBondMCeneobj hbondmc = HBondMCeneobj::create(&refinestr_temp[j].ptv[jj], &refinestr_temp[j].ptv[m], &refinestr_temp[i].ptv[ii], &pairwiseenergyatoms[j][m], &pairwiseenergyatoms[i][ii], &mapflag[i][j]);
										hbondene.push_back(hbondmc);
									}
								}
							}
						}
						*/

						/*
						if(ii == 0 && jj == 0){
							if(k < inputtedstrallflag[i].size() && l < inputtedstrallflag[j].size()){
								if(inputtedstrallflag[i][k] && inputtedstrallflag[j][l]){
									Point3s initatomA = inputtedstrall.residues[i].atoms[k].coord;
									Point3s initatomB = inputtedstrall.residues[j].atoms[l].coord;

									Point3s initvec;
									initvec.x = initatomA.x - initatomB.x;
									initvec.y = initatomA.y - initatomB.y;
									initvec.z = initatomA.z - initatomB.z;

									double initdist = sqrt(initvec.x*initvec.x+initvec.y*initvec.y+initvec.z*initvec.z);
									
									InitMCdistobj initmc = InitMCdistobj::create(&refinestr_temp[i].ptv[ii], &refinestr_temp[j].ptv[jj], initdist, 20, &pairwiseenergyatoms[i][ii], &pairwiseenergyatoms[j][jj], 1.0/(1.0+templatestrall.size()), &mapflag[i][j]);
									initene.push_back(initmc);
								}
							}
						}
						*/
					}
				}

				for(int m=0;m<template_target_atoms.size();m++){
					int k = template_target_atoms_idx[m];
					int l = template_target_atoms_idx[m];
					int ii = atomid2mcidx[base1][k];
					int jj = atomid2mcidx[base2][l];
					if(m < (*MCinitdistflag)[i][j].size()){
						if((*MCinitdistflag)[i][j][m]){

							Point3s initvec;
							initvec.x = refinestr_temp[i].ptv[ii].x - refinestr_temp[j].ptv[jj].x;
							initvec.y = refinestr_temp[i].ptv[ii].y - refinestr_temp[j].ptv[jj].y;
							initvec.z = refinestr_temp[i].ptv[ii].z - refinestr_temp[j].ptv[jj].z;
							double initdist = sqrt(initvec.x*initvec.x+initvec.y*initvec.y+initvec.z*initvec.z);

							//InitMCdistobj initmc = InitMCdistobj::create(&refinestr_temp[i].ptv[ii], &refinestr_temp[j].ptv[jj], (*MCinitdistmap)[i][j][m], 20, &pairwiseenergyatoms[i][ii], &pairwiseenergyatoms[j][jj], (*MCinitdistscore)[i][j]*(1.0/(1.0+templatestrall.size())), &mapflag[i][j]);
							
							InitMCdistobj initmc = InitMCdistobj::create(&refinestr_temp[i].ptv[ii], &refinestr_temp[j].ptv[jj], (*MCinitdistmap)[i][j][m], 20, &pairwiseenergyatoms[i][ii], &pairwiseenergyatoms[j][jj], (1.0/(1.0+templatestrall.size())), &mapflag[i][j]);
							initene.push_back(initmc);
						}
					}
				}

				for(int n=0;n<templatestrall.size();n++){
					for(int m=0;m<template_target_atoms.size();m++){
						int k = template_target_atoms_idx[m];
						int l = template_target_atoms_idx[m];
						int ii = atomid2mcidx[base1][k];
						int jj = atomid2mcidx[base2][l];

						if(m < (*MCtemplatedistflag)[n][i][j].size()){
							if((*MCtemplatedistflag)[n][i][j][m]){
								InitMCdistobj initmc = InitMCdistobj::create(&refinestr_temp[i].ptv[ii], &refinestr_temp[j].ptv[jj], (*MCtemplatedistmap)[n][i][j][m], 20, &pairwiseenergyatoms[i][ii], &pairwiseenergyatoms[j][jj], 1.0/(1.0+templatestrall.size()), &mapflag[i][j]);
								templateene.push_back(initmc);

							}
						}
					}
				}

			}else{
				for(int ii=0;ii<atomid2mcidx[base1].size();ii++)
				{
					for(int jj=ii+1;jj<atomid2mcidx[base2].size();jj++)
					{
						k = mcidx2atomid[base1][ii];
						l = mcidx2atomid[base2][jj];

						if(molConnect[ii][jj] >= 3){
							if(molConnect[ii][jj] == 3){
								dQQ = (COULOMB*extparameters.atomfeatures.at(base1_str).at(k)[0]*extparameters.atomfeatures.at(base2_str).at(l)[0])/1.2;
							}else{
								dQQ = COULOMB*extparameters.atomfeatures.at(base1_str).at(k)[0]*extparameters.atomfeatures.at(base2_str).at(l)[0];
							}
							dQQep = dQQ/EpsilonRNA;

							ElectrosMCeneobj elmc = ElectrosMCeneobj::create(&refinestr_temp[i].ptv[ii], &refinestr_temp[j].ptv[jj], &pairwiseenergyatoms[i][ii], &pairwiseenergyatoms[j][jj], &mapflag[i][j], dQQep);
							electrosene.push_back(elmc);
						}

						const Point3s &pt1 = refinestr_temp[i].ptv[ii];
						const Point3s &pt2 = refinestr_temp[j].ptv[jj];

						Point3s vec12;
						vec12.x = pt2.x - pt1.x;
						vec12.y = pt2.y - pt1.y;
						vec12.z = pt2.z - pt1.z;

						double length12=vec12.x*vec12.x+vec12.y*vec12.y+vec12.z*vec12.z;

						if(length12 < CUTOFF_SQR+5){
							double vdWInfodR01 = extparameters.atomfeatures.at(base1_str).at(k)[1];
							double vdWInfodR02 = extparameters.atomfeatures.at(base2_str).at(l)[1];
							double vdWInfodEpsilon1 = extparameters.atomfeatures.at(base1_str).at(k)[2];
							double vdWInfodEpsilon2 = extparameters.atomfeatures.at(base2_str).at(l)[2];
							double star = (vdWInfodR01 + vdWInfodR02);
							double epsilon =  sqrt(vdWInfodEpsilon1*vdWInfodEpsilon2);
							double LJ6 = -2*epsilon*pow(star,6);
							double LJ12 = epsilon*pow(star,12);

							if(molConnect[ii][jj] == 3){
								LJ6 /=2.0;
								LJ12 /=2.0;
							}

							LJMCeneobj ljmc = LJMCeneobj::create(&refinestr_temp[i].ptv[ii], &refinestr_temp[j].ptv[jj], &pairwiseenergyatoms[i][ii], &pairwiseenergyatoms[j][jj], &mapflag[i][j], LJ6, LJ12);
							ljene.push_back(ljmc);
						}

						/*
						if(k < inputtedstrallflag[i].size() && l < inputtedstrallflag[j].size()){
							if(inputtedstrallflag[i][k] && inputtedstrallflag[j][l]){
								Point3s initatomA = inputtedstrall.residues[i].atoms[k].coord;
								Point3s initatomB = inputtedstrall.residues[j].atoms[l].coord;

								Point3s initvec;
								initvec.x = initatomA.x - initatomB.x;
								initvec.y = initatomA.y - initatomB.y;
								initvec.z = initatomA.z - initatomB.z;

								double initdist = sqrt(initvec.x*initvec.x+initvec.y*initvec.y+initvec.z*initvec.z);

								InitMCdistobj initmc = InitMCdistobj::create(&refinestr_temp[i].ptv[ii], &refinestr_temp[j].ptv[jj], initdist, 20, &pairwiseenergyatoms[i][ii], &pairwiseenergyatoms[j][jj], &mapflag[i][j]);
								initene.push_back(initmc);
							}
						}
						*/
					}
				}
			}
		}
	}
}

double MonteCarloSimulation::calceneclash(){
	double totalene = 0;
	double ene;
	clashenergies.clear();

	for (std::list<ClashMCeneobj>::const_iterator it = clashmcene.begin(); it != clashmcene.end(); ++it) {
		ene=it->calcenergyclash();
		clashenergies.push_back(ene);
		totalene+=ene;
	}

	return totalene;
}

double MonteCarloSimulation::calceneelectros(){
	double totalene = 0;
	double ene;
	electrosenergies.clear();
	
	for (std::list<ElectrosMCeneobj>::const_iterator it = electrosene.begin(); it != electrosene.end(); ++it) {
		ene=it->calcenergyelectros();
		electrosenergies.push_back(ene);
		totalene+=ene;
	}
	return totalene;
}

double MonteCarloSimulation::calcenelj(){
	double totalene = 0;
	double ene;
	ljenergies.clear();
	
	for (std::list<LJMCeneobj>::const_iterator it = ljene.begin(); it != ljene.end(); ++it) {
		ene=it->calcenergylj();
		ljenergies.push_back(ene);
		totalene+=ene;
	}
	return totalene;
}

double MonteCarloSimulation::calcenehbond(){
	double totalene = 0;
	double ene;
	hbondenergies.clear();
	
	for (std::list<HBondMCeneobj>::const_iterator it = hbondene.begin(); it != hbondene.end(); ++it) {
		ene=it->calcenergyhydbond();
		hbondenergies.push_back(ene);
		totalene+=ene;
	}

	hbond2energies.clear();
	for (std::list<HBond2MCeneobj>::const_iterator it = hbond2ene.begin(); it != hbond2ene.end(); ++it) {
		ene=it->calcenergyhydbond();
		hbond2energies.push_back(ene);
		totalene+=ene;
	}

	return totalene;
}

void MonteCarloSimulation::updateeneclash(){
	int count=0;

	for (std::list<ClashMCeneobj>::const_iterator it = clashmcene.begin(); it != clashmcene.end(); ++it) {
		it->setprevene(clashenergies[count]);
		count += 1;
	}
}

void MonteCarloSimulation::updateeneelectros(){
	int count=0;
	
	for (std::list<ElectrosMCeneobj>::const_iterator it = electrosene.begin(); it != electrosene.end(); ++it) {
		it->setprevene(electrosenergies[count]);
		count += 1;
	}
}

void MonteCarloSimulation::updateenelj(){
	int count=0;
	
	for (std::list<LJMCeneobj>::const_iterator it = ljene.begin(); it != ljene.end(); ++it) {
		it->setprevene(ljenergies[count]);
		count += 1;
	}
}

void MonteCarloSimulation::updateenehbond(){
	int count=0;
	for (std::list<HBondMCeneobj>::const_iterator it = hbondene.begin(); it != hbondene.end(); ++it) {
		it->setprevene(hbondenergies[count]);
		count += 1;
	}

	count=0;
	for (std::list<HBond2MCeneobj>::const_iterator it = hbond2ene.begin(); it != hbond2ene.end(); ++it) {
		it->setprevene(hbond2energies[count]);
		count += 1;
	}
}

double MonteCarloSimulation::calceneinitpos(){
	double totalene = 0;
	double ene;
	initposenergies.clear();
	
	for (std::list<InitMCPosobj>::const_iterator it = initposene.begin(); it != initposene.end(); ++it) {
		ene=it->calcenergyinit();
		initposenergies.push_back(ene);
		totalene+=ene;
	}
	return totalene;
}

void MonteCarloSimulation::updateeneinitpos(){
	int count=0;
	
	for (std::list<InitMCPosobj>::const_iterator it = initposene.begin(); it != initposene.end(); ++it) {
		it->setprevene(initposenergies[count]);
		count += 1;
	}
}

double MonteCarloSimulation::calceneinit(){
	double totalene = 0;
	double ene;
	initenergies.clear();
	
	for (std::list<InitMCdistobj>::const_iterator it = initene.begin(); it != initene.end(); ++it) {
		ene=it->calcenergyinit();
		initenergies.push_back(ene);
		totalene+=ene;
	}
	return totalene;
}

void MonteCarloSimulation::updateeneinit(){
	int count=0;
	
	for (std::list<InitMCdistobj>::const_iterator it = initene.begin(); it != initene.end(); ++it) {
		it->setprevene(initenergies[count]);
		count += 1;
	}
}

double MonteCarloSimulation::calcenetemplate(){
	double totalene = 0;
	double ene;
	templateenergies.clear();
	
	for (std::list<InitMCdistobj>::const_iterator it = templateene.begin(); it != templateene.end(); ++it) {
		ene=it->calcenergyinit();
		templateenergies.push_back(ene);
		totalene+=ene;
	}
	return totalene;
}

void MonteCarloSimulation::updateenetemplate(){
	int count=0;
	
	for (std::list<InitMCdistobj>::const_iterator it = templateene.begin(); it != templateene.end(); ++it) {
		it->setprevene(templateenergies[count]);
		count += 1;
	}
}

double MonteCarloSimulation::calcenebp(){
	double totalene = 0;
	double ene;
	bpenergies.clear();
	int count=0;
	
	/*
	for (std::list<BasePairMCobj>::const_iterator it = bpene.begin(); it != bpene.end(); ++it) {
		ene=it->calcenergybasepair();
		bpenergies.push_back(ene);
		totalene+=ene;
		count+=1;
	}
	*/

	for (std::list<BasePairMCobj2>::iterator it = bpene2.begin(); it != bpene2.end(); ++it) {
		ene=it->calcenergybasepair();
		bpenergies.push_back(ene);
		totalene+=ene;
		count+=1;
	}
	return totalene;
}

void MonteCarloSimulation::updateenebp(){
	int count=0;
	
	/*
	for (std::list<BasePairMCobj>::const_iterator it = bpene.begin(); it != bpene.end(); ++it) {
		it->setprevene(bpenergies[count]);
		count += 1;
	}
	*/

	for (std::list<BasePairMCobj2>::const_iterator it = bpene2.begin(); it != bpene2.end(); ++it) {
		it->setprevene(bpenergies[count]);
		count += 1;
	}
}

double MonteCarloSimulation::calcenebaseplane(){
	double totalene = 0;
	double ene;
	basepenergies.clear();
	int count=0;

	for (std::list<BasePairMCobj2>::iterator it = baseplane.begin(); it != baseplane.end(); ++it) {
		ene=it->calcenergybasepair();
		basepenergies.push_back(ene);
		totalene+=ene;
		count+=1;
	}
	return totalene;
}

void MonteCarloSimulation::updateenebaseplane(){
	int count=0;

	for (std::list<BasePairMCobj2>::const_iterator it = baseplane.begin(); it != baseplane.end(); ++it) {
		it->setprevene(basepenergies[count]);
		count += 1;
	}
}

double MonteCarloSimulation::calcenestack(){
	double totalene = 0;
	double ene;
	stackenergies.clear();
	int count=0;
	
	for (std::list<StackPairMCobj>::const_iterator it = stackene.begin(); it != stackene.end(); ++it) {
		ene=it->calcenergystacking();
		stackenergies.push_back(ene);
		totalene+=ene;
		count+=1;
	}
	return totalene;
}

void MonteCarloSimulation::updateenestack(){
	int count=0;
	
	for (std::list<StackPairMCobj>::const_iterator it = stackene.begin(); it != stackene.end(); ++it) {
		it->setprevene(stackenergies[count]);
		count += 1;
	}
}

double MonteCarloSimulation::calcenestackgeo(){
	double totalene = 0;
	double ene;
	stackgeoenergies.clear();
	int count=0;
	
	for (std::list<StackPairGeoMCobj>::const_iterator it = stackgeoene.begin(); it != stackgeoene.end(); ++it) {
		ene=it->calcenergystacking();
		stackgeoenergies.push_back(ene);
		totalene+=ene;
		count+=1;
	}
	return totalene;
}

void MonteCarloSimulation::updateenestackgeo(){
	int count=0;
	
	for (std::list<StackPairGeoMCobj>::const_iterator it = stackgeoene.begin(); it != stackgeoene.end(); ++it) {
		it->setprevene(stackgeoenergies[count]);
		count += 1;
	}
}


















































