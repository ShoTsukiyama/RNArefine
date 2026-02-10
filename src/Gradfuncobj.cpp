// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include <cmath>
#include "Gradfuncobj.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

double InitPairobj::calcenergyinitpair() const
{
	Point3s vec12;

	vec12.x = teAtomB->x - teAtomA->x;
	vec12.y = teAtomB->y - teAtomA->y;
	vec12.z = teAtomB->z - teAtomA->z;
	
	double length12=sqrt(vec12.x*vec12.x+vec12.y*vec12.y+vec12.z*vec12.z+1e-8);

	double dist = length12 - ideal;
	double diff = dist*dist;

	if(diff < 400) return KI*diff;
	else return KI*400;

	//return KI*dist*dist;
}

void InitPairobj::calcgradientinitpair() const
{
	Point3s vec12;

	vec12.x = teAtomB->x - teAtomA->x;
	vec12.y = teAtomB->y - teAtomA->y;
	vec12.z = teAtomB->z - teAtomA->z;

	double length12=sqrt(vec12.x*vec12.x+vec12.y*vec12.y+vec12.z*vec12.z+1e-8)+1e-8;

	vec12.x /= length12;
	vec12.y /= length12;
	vec12.z /= length12;

	double diff = length12 - ideal;

	if ((diff * diff) >= 400.0) return; 
	
	double coeff = 2.0 * KI * diff;

	teAtomAgrad->x -= coeff*vec12.x;
	teAtomAgrad->y -= coeff*vec12.y;
	teAtomAgrad->z -= coeff*vec12.z;
	teAtomBgrad->x += coeff*vec12.x;
	teAtomBgrad->y += coeff*vec12.y;
	teAtomBgrad->z += coeff*vec12.z;
}

double InitPosobj::calcenergyinit() const
{
	Point3s vec;
	vec.x = X0 - teAtomA->x;
	vec.y = Y0 - teAtomA->y;
	vec.z = Z0 - teAtomA->z;

	//cout<<"=================="<<endl;
	//cout<<X0 << "-" << Y0 << "-" << Z0<<endl;
	//cout<<teAtomA->x << "-" << teAtomA->y << "-" << teAtomA->z<<endl;

	double dist = sqrt(vec.x*vec.x+vec.y*vec.y+vec.z*vec.z+1e-8);

	//cout<<dist<<endl;
	//if(std::isinf(Kr*dist*dist)){
		//cout<<"Error energy init"<<endl;
		//cout<<X0 << "-" << Y0 << "-" << Z0<<endl;
		//cout<<teAtomA->x << "-" << teAtomA->y << "-" << teAtomA->z<<endl;
		//cout<<dist<<endl;
		//std::exit(EXIT_FAILURE);
	//}

	return Kr*dist*dist;
}

void InitPosobj::calcgradientinit() const
{
	teAtomAgrad->x -= 2.0*Kr*(X0 - teAtomA->x);
	teAtomAgrad->y -= 2.0*Kr*(Y0 - teAtomA->y);
	teAtomAgrad->z -= 2.0*Kr*(Z0 - teAtomA->z);

	//if (std::isinf(teAtomAgrad->x) || std::isinf(teAtomAgrad->y) || std::isinf(teAtomAgrad->z)) {
		//cout<<"Error grad init"<<endl;
		//cout<<X0 << "-" << Y0 << "-" << Z0<<endl;
		//cout<<teAtomA->x << "-" << teAtomA->y << "-" << teAtomA->z<<endl;
		//std::exit(EXIT_FAILURE);
	//}
}

double TElectrobj::calcenergyelectros() const
{
	Point3s vec12;

	vec12.x = teAtomB->x - teAtomA->x;
	vec12.y = teAtomB->y - teAtomA->y;
	vec12.z = teAtomB->z - teAtomA->z;

	double length12=vec12.x*vec12.x+vec12.y*vec12.y+vec12.z*vec12.z;

	//if(std::isinf(dQQep/(sqrt(length12)+1e-8))){
		//cout<<"Error energy electros"<<endl;
		//std::exit(EXIT_FAILURE);
	//}

	return dQQep/(sqrt(length12+1e-8)+1e-8);
}

void TElectrobj::calcgradientelectros() const
{
	Point3s vec1, vec2, vec12;

	vec12.x = teAtomB->x - teAtomA->x;
	vec12.y = teAtomB->y - teAtomA->y;
	vec12.z = teAtomB->z - teAtomA->z;

	double length12SQR=vec12.x*vec12.x+vec12.y*vec12.y+vec12.z*vec12.z;
	double length12 = sqrt(length12SQR+1e-8);
	double factor = dQQep/(length12SQR*length12+1e-8);

	vec12.x *= factor;
	vec12.y *= factor;
	vec12.z *= factor;

	teAtomAgrad->x += vec12.x;
	teAtomAgrad->y += vec12.y;
	teAtomAgrad->z += vec12.z;
	teAtomBgrad->x -= vec12.x;
	teAtomBgrad->y -= vec12.y;
	teAtomBgrad->z -= vec12.z;

	//if (std::isinf(teAtomAgrad->x) || std::isinf(teAtomAgrad->y) || std::isinf(teAtomAgrad->z) || std::isinf(teAtomBgrad->x) || std::isinf(teAtomBgrad->y) || std::isinf(teAtomBgrad->z)) {
		//cout<<"Error electros"<<endl;
		//std::exit(EXIT_FAILURE);
	//}
}

double LennardsJobj::calcenergylj() const
{
	Point3s vec12;

	vec12.x = teAtomB->x - teAtomA->x;
	vec12.y = teAtomB->y - teAtomA->y;
	vec12.z = teAtomB->z - teAtomA->z;

	double length12SQR=vec12.x*vec12.x+vec12.y*vec12.y+vec12.z*vec12.z;

	if(length12SQR > CUTOFF_SQR) return 0.0;

	//if(std::isinf((LJ6/(length12SQR*length12SQR*length12SQR+1e-8))+(LJ12/(length12SQR*length12SQR*length12SQR*length12SQR*length12SQR*length12SQR+1e-8)))){
		//cout<<"Error energy lj"<<endl;
		//std::exit(EXIT_FAILURE);
	//}

	//if(((LJ6/(length12SQR*length12SQR*length12SQR+1e-12))+(LJ12/(length12SQR*length12SQR*length12SQR*length12SQR*length12SQR*length12SQR+1e-12))) > 10){
		//cout<<(LJ6/(length12SQR*length12SQR*length12SQR+1e-12)) << "-" << (LJ12/(length12SQR*length12SQR*length12SQR*length12SQR*length12SQR*length12SQR+1e-12))<<endl;
	//}

	return K*((LJ6/(length12SQR*length12SQR*length12SQR+1e-12))+(LJ12/(length12SQR*length12SQR*length12SQR*length12SQR*length12SQR*length12SQR+1e-12)));
}

void LennardsJobj::calcgradientlj() const
{
	Point3s vec12;

	vec12.x = teAtomB->x - teAtomA->x;
	vec12.y = teAtomB->y - teAtomA->y;
	vec12.z = teAtomB->z - teAtomA->z;

	double length12SQR=vec12.x*vec12.x+vec12.y*vec12.y+vec12.z*vec12.z;

	if(length12SQR > CUTOFF_SQR) return;

	double factor = K*((6.0*LJ6/(length12SQR*length12SQR*length12SQR*length12SQR+1e-8))+(12.0*LJ12/(length12SQR*length12SQR*length12SQR*length12SQR*length12SQR*length12SQR*length12SQR+1e-8)));

	vec12.x *= factor;
	vec12.y *= factor;
	vec12.z *= factor;

	teAtomAgrad->x += vec12.x;
	teAtomAgrad->y += vec12.y;
	teAtomAgrad->z += vec12.z;
	teAtomBgrad->x -= vec12.x;
	teAtomBgrad->y -= vec12.y;
	teAtomBgrad->z -= vec12.z;

	//if (std::isinf(teAtomAgrad->x) || std::isinf(teAtomAgrad->y) || std::isinf(teAtomAgrad->z) || std::isinf(teAtomBgrad->x) || std::isinf(teAtomBgrad->y) || std::isinf(teAtomBgrad->z)) {
		//cout<<"Error LJ"<<endl;
		//std::exit(EXIT_FAILURE);
	//}
}

double BondLengthobj::calcenergylength() const
{
	Point3s vec12;

	vec12.x = teAtomB->x - teAtomA->x;
	vec12.y = teAtomB->y - teAtomA->y;
	vec12.z = teAtomB->z - teAtomA->z;
	
	double length12=sqrt(vec12.x*vec12.x+vec12.y*vec12.y+vec12.z*vec12.z+1e-8);

	double dist = length12 - ideal;

	//if(fabs(dist) < 0.5*stdval) return 0.0;
	//cout<<fabs(dist) << "-" << ideal << "-" << stdval<<endl;

	return KI*dist*dist;
}

void BondLengthobj::calcgradientlength() const
{
	Point3s vec12;

	vec12.x = teAtomB->x - teAtomA->x;
	vec12.y = teAtomB->y - teAtomA->y;
	vec12.z = teAtomB->z - teAtomA->z;

	double length12=sqrt(vec12.x*vec12.x+vec12.y*vec12.y+vec12.z*vec12.z+1e-8)+1e-8;

	vec12.x /= length12;
	vec12.y /= length12;
	vec12.z /= length12;

	double dist = length12 - ideal;
	//if(fabs(dist) < 0.5*stdval) return;

	double coeff=2.*KI*dist;

	teAtomAgrad->x -= coeff*vec12.x;
	teAtomAgrad->y -= coeff*vec12.y;
	teAtomAgrad->z -= coeff*vec12.z;
	teAtomBgrad->x += coeff*vec12.x;
	teAtomBgrad->y += coeff*vec12.y;
	teAtomBgrad->z += coeff*vec12.z;

	//if (std::isinf(teAtomAgrad->x) || std::isinf(teAtomAgrad->y) || std::isinf(teAtomAgrad->z) || std::isinf(teAtomBgrad->x) || std::isinf(teAtomBgrad->y) || std::isinf(teAtomBgrad->z)) {
		//cout<<"Error length"<<endl;
		//std::exit(EXIT_FAILURE);
	//}

}

Point3s BondAngleobj::basecross(Point3s p1, Point3s p2)  const
{
	Point3s temp;
	temp.x=p1.y * p2.z - p1.z * p2.y;
	temp.y=p1.z * p2.x - p1.x * p2.z;
	temp.z=p1.x * p2.y - p1.y * p2.x;
	return temp;
}

Point3s basecross(Point3s p1, Point3s p2)
{
	Point3s temp;
	temp.x=p1.y * p2.z - p1.z * p2.y;
	temp.y=p1.z * p2.x - p1.x * p2.z;
	temp.z=p1.x * p2.y - p1.y * p2.x;
	return temp;
}

double BondAngleobj::calcenergyangle() const
{
	Point3s vec21, vec23;

	vec21.x = teAtomA->x - teAtomB->x;
	vec21.y = teAtomA->y - teAtomB->y;
	vec21.z = teAtomA->z - teAtomB->z;

	vec23.x = teAtomC->x - teAtomB->x;
	vec23.y = teAtomC->y - teAtomB->y;
	vec23.z = teAtomC->z - teAtomB->z;

	double length21=sqrt(vec21.x*vec21.x+vec21.y*vec21.y+vec21.z*vec21.z+1e-8);
	double length23=sqrt(vec23.x*vec23.x+vec23.y*vec23.y+vec23.z*vec23.z+1e-8);

	vec21.x /= length21; vec21.y /= length21; vec21.z /= length21;
	vec23.x /= length23; vec23.y /= length23; vec23.z /= length23;

	double  costheta = vec21.x*vec23.x+vec21.y*vec23.y+vec21.z*vec23.z;
	if(     costheta > 1. ) costheta = 1.;
  	else if(costheta < -1.) costheta = -1.;

	double dist = acos(costheta) - ideal;

	//if(fabs(dist) < 0.5*stdval) return 0.0;
	//cout<<fabs(dist) << "-" << ideal << "-" << stdval<<endl;

	return KI*dist*dist;
}

void BondAngleobj::calcgradientangle() const
{
	Point3s vec21, vec23, outerprot, grad1, grad2, grad3;

	vec21.x = teAtomA->x - teAtomB->x;
	vec21.y = teAtomA->y - teAtomB->y;
	vec21.z = teAtomA->z - teAtomB->z;

	vec23.x = teAtomC->x - teAtomB->x;
	vec23.y = teAtomC->y - teAtomB->y;
	vec23.z = teAtomC->z - teAtomB->z;

	double length21=sqrt(vec21.x*vec21.x+vec21.y*vec21.y+vec21.z*vec21.z+1e-8);
	double length23=sqrt(vec23.x*vec23.x+vec23.y*vec23.y+vec23.z*vec23.z+1e-8);

	vec21.x /= length21; vec21.y /= length21; vec21.z /= length21;
	vec23.x /= length23; vec23.y /= length23; vec23.z /= length23;

	double  costheta = vec21.x*vec23.x+vec21.y*vec23.y+vec21.z*vec23.z;
	if(     costheta > 1. ) costheta = 1.;
  	else if(costheta < -1.) costheta = -1.;

	double dist = acos(costheta) - ideal;
	//if(fabs(dist) < 0.5*stdval) return;
	double energydepth = 2.*KI*(dist);

	outerprot = basecross(vec21, vec23);
	double sintheta = sqrt(outerprot.x*outerprot.x + outerprot.y*outerprot.y + outerprot.z*outerprot.z+1e-8);
	
	if(fabs(sintheta) < 0.00001 || isnan(sintheta)) return;	// Protection against numerical instabilities.

	double inerprot = vec21.x*vec23.x + vec21.y*vec23.y + vec21.z*vec23.z;
	
	grad1.x = (inerprot*vec21.x - vec23.x)/length21;
	grad1.y = (inerprot*vec21.y - vec23.y)/length21;
	grad1.z = (inerprot*vec21.z - vec23.z)/length21;

	grad3.x = (inerprot*vec23.x - vec21.x)/length23;
	grad3.y = (inerprot*vec23.y - vec21.y)/length23;
	grad3.z = (inerprot*vec23.z - vec21.z)/length23;

	double coef = energydepth/sintheta;
	grad1.x *= coef; grad1.y *= coef; grad1.z *= coef;
	grad3.x *= coef; grad3.y *= coef; grad3.z *= coef;
	grad2.x = -(grad1.x+grad3.x); grad2.y = -(grad1.y+grad3.y); grad2.z = -(grad1.z+grad3.z); 

	teAtomAgrad->x += grad1.x;
	teAtomAgrad->y += grad1.y;
	teAtomAgrad->z += grad1.z;
	teAtomCgrad->x += grad3.x;
	teAtomCgrad->y += grad3.y;
	teAtomCgrad->z += grad3.z;
	teAtomBgrad->x += grad2.x;
	teAtomBgrad->y += grad2.y;
	teAtomBgrad->z += grad2.z;

	//if (std::isinf(teAtomAgrad->x) || std::isinf(teAtomAgrad->y) || std::isinf(teAtomAgrad->z) || std::isinf(teAtomBgrad->x) || std::isinf(teAtomBgrad->y) || std::isinf(teAtomBgrad->z) || std::isinf(teAtomCgrad->x) || std::isinf(teAtomCgrad->y) || std::isinf(teAtomCgrad->z)) {
		//cout<<"Error angle"<<endl;
		//std::exit(EXIT_FAILURE);
	//}
}

double BondDihedralobj::calcenergydihedral()  const
{
	Point3s vec12, vec23, vec34, vecnormal123, vecnormal234, outerprot, a, b, grad1, grad2, grad3, grad4;
	Point3s vec123, vec234, outerprot23a, outerprota123, outerprot34b, outerprotb234, outerprot12a, outerprot23b;

	vec12.x = teAtomB->x - teAtomA->x;
	vec12.y = teAtomB->y - teAtomA->y;
	vec12.z = teAtomB->z - teAtomA->z;
	
	vec23.x = teAtomC->x - teAtomB->x;
	vec23.y = teAtomC->y - teAtomB->y;
	vec23.z = teAtomC->z - teAtomB->z;

	vec34.x = teAtomD->x - teAtomC->x;
	vec34.y = teAtomD->y - teAtomC->y;
	vec34.z = teAtomD->z - teAtomC->z;

	vecnormal123 = basecross(vec12,vec23);
	vecnormal234 = basecross(vec23,vec34);

	double length123=sqrt(max(double(vecnormal123.x*vecnormal123.x+vecnormal123.y*vecnormal123.y+vecnormal123.z*vecnormal123.z),1e-8));
	double length234=sqrt(max(double(vecnormal234.x*vecnormal234.x+vecnormal234.y*vecnormal234.y+vecnormal234.z*vecnormal234.z),1e-8));

	vecnormal123.x /= length123; vecnormal123.y /= length123; vecnormal123.z /= length123;
	vecnormal234.x /= length234; vecnormal234.y /= length234; vecnormal234.z /= length234;

	double cosphi=vecnormal123.x*vecnormal234.x+vecnormal123.y*vecnormal234.y+vecnormal123.z*vecnormal234.z;
	
	if(cosphi>1.) cosphi=1.;
	else if(cosphi<-1.) cosphi=-1.;
	
	outerprot=basecross(vecnormal123,vecnormal234);
	double dval = acos(cosphi) * ((vec23.x*outerprot.x+vec23.y*outerprot.y+vec23.z*outerprot.z) >= 0. ? 1.:-1.);
	
	if (dval < 0.0) dval += 2.0 * M_PI;

	double diff = 1 + cos(rotn*dval - ideal);

	return K * diff;
}


void BondDihedralobj::calcgradientdihedral() const
{
	Point3s vec12, vec23, vec34, vecnormal123, vecnormal234, outerprot, a, b; 
    Point3s outerprot23a, outerprot12a, outerprot34d, outerprot23b;

    // 1. ベクトルの計算 (変更なし)
    vec12.x = teAtomB->x - teAtomA->x;
    vec12.y = teAtomB->y - teAtomA->y;
    vec12.z = teAtomB->z - teAtomA->z;
    
    vec23.x = teAtomC->x - teAtomB->x;
    vec23.y = teAtomC->y - teAtomB->y;
    vec23.z = teAtomC->z - teAtomB->z;

    vec34.x = teAtomD->x - teAtomC->x;
    vec34.y = teAtomD->y - teAtomC->y;
    vec34.z = teAtomD->z - teAtomC->z;

	//cout<<teAtomA->x << "," << teAtomA->y << "," << teAtomA->z<<endl;
	//cout<<teAtomB->x << "," << teAtomB->y << "," << teAtomB->z<<endl;
	//cout<<teAtomC->x << "," << teAtomC->y << "," << teAtomC->z<<endl;
	//cout<<teAtomD->x << "," << teAtomD->y << "," << teAtomD->z<<endl;

    vecnormal123 = basecross(vec12,vec23);
	vecnormal234 = basecross(vec23,vec34);

	double length123=sqrt(max(double(vecnormal123.x*vecnormal123.x+vecnormal123.y*vecnormal123.y+vecnormal123.z*vecnormal123.z),1e-8));
	double length234=sqrt(max(double(vecnormal234.x*vecnormal234.x+vecnormal234.y*vecnormal234.y+vecnormal234.z*vecnormal234.z),1e-8));

	vecnormal123.x /= length123; vecnormal123.y /= length123; vecnormal123.z /= length123;
	vecnormal234.x /= length234; vecnormal234.y /= length234; vecnormal234.z /= length234;

	double innerprot=vecnormal123.x*vecnormal234.x+vecnormal123.y*vecnormal234.y+vecnormal123.z*vecnormal234.z;

	if(innerprot>1.) innerprot=1.; 
	else if(innerprot<-1.) innerprot=-1.; 

	outerprot=basecross(vecnormal123,vecnormal234);

	double dval = acos(innerprot) * ((vec23.x*outerprot.x+vec23.y*outerprot.y+vec23.z*outerprot.z) >= 0. ? 1.:-1.);

	double sinphi = sin(dval);
    
	if(fabs(sinphi) < 1e-8 || isnan(sinphi)) return;

	//cout<<tor+1 << "-" << minideals.size()<<endl;
	double diff = rotn * dval - ideal;

	//cout<<1 + cos(rotn*dval - ideal)<<endl;

	//cout<<dval*DEGRAD << "-" << +cos(diff) << "-" << cos(diff)<<endl;

	//cout<<dval*DEGRAD << "-" << ideal*DEGRAD+180 << "-" << cos(rotn*dval - ideal)<<endl;

	//double diff = dval - (*minideals[tor])*RADDEG;
	//double s  = sin(diff);
	//double c  = cos(diff);
	//double delta = atan2(s, c);

	//double ene = K * (1.0 - cos(delta));
	double gradient_E_to_phi = K * rotn * sin(diff);

	//cout<<"---------"<<endl;
	//cout<< tor+1 << "-" << 1 + cos(diff) << "-" << dval*DEGRAD << "-" << phi0 * DEGRAD<<endl;
	//cout<<"---------"<<endl;

	//double gradient_scale = -(K * sin(delta)) / sinphi;
	double gradient_scale = gradient_E_to_phi / sinphi;

	double inv_len123_sq = 1.0 / sqrt(length123 * length123);
	double inv_len234_sq = 1.0 / sqrt(length234 * length234);

	a.x = gradient_scale * (vecnormal234.x - vecnormal123.x * innerprot)*inv_len123_sq;
	a.y = gradient_scale * (vecnormal234.y - vecnormal123.y * innerprot)*inv_len123_sq;
	a.z = gradient_scale * (vecnormal234.z - vecnormal123.z * innerprot)*inv_len123_sq;

	b.x = gradient_scale * (vecnormal123.x - vecnormal234.x * innerprot)*inv_len234_sq;
	b.y = gradient_scale * (vecnormal123.y - vecnormal234.y * innerprot)*inv_len234_sq;
	b.z = gradient_scale * (vecnormal123.z - vecnormal234.z * innerprot)*inv_len234_sq;

	outerprot12a = basecross(a,vec12);     // vec(A->B) x a
	outerprot23a = basecross(vec23,a);     // vec(B->C) x a
	outerprot23b = basecross(b,vec23);     // vec(B->C) x b
	outerprot34d = basecross(vec34,b);     // vec(C->D) x b

	teAtomAgrad->x += (-outerprot23a.x); 
	teAtomAgrad->y += (-outerprot23a.y); 
	teAtomAgrad->z += (-outerprot23a.z);

	teAtomBgrad->x += (-outerprot34d.x-outerprot12a.x+outerprot23a.x); 
	teAtomBgrad->y += (-outerprot34d.y-outerprot12a.y+outerprot23a.y);
	teAtomBgrad->z += (-outerprot34d.z-outerprot12a.z+outerprot23a.z); 

	teAtomCgrad->x += (outerprot12a.x+outerprot34d.x-outerprot23b.x);
	teAtomCgrad->y += (outerprot12a.y+outerprot34d.y-outerprot23b.y);
	teAtomCgrad->z += (outerprot12a.z+outerprot34d.z-outerprot23b.z); 

	teAtomDgrad->x += (outerprot23b.x); 
	teAtomDgrad->y += (outerprot23b.y); 
	teAtomDgrad->z += (outerprot23b.z);

	//cout<<teAtomAgrad->x << "-" << teAtomAgrad->y << "-" << teAtomAgrad->z<<endl;
	//cout<<teAtomBgrad->x << "-" << teAtomBgrad->y << "-" << teAtomBgrad->z<<endl;
	//cout<<teAtomCgrad->x << "-" << teAtomCgrad->y << "-" << teAtomCgrad->z<<endl;
	//cout<<teAtomDgrad->x << "-" << teAtomDgrad->y << "-" << teAtomDgrad->z<<endl;
}

MultiBondDihedralobj::MultiBondDihedralobj(const vector<Point3s*>& _AtomPtrs, const vector<Point3s*>& _Atomgrad, const vector<vector<vector<double*>>> &_ideals, const double &_K)
    : teAtoms(_AtomPtrs), teAtomgrad(_Atomgrad), ideals(_ideals), K(_K)
{
    param_cls_num = ideals.size();
}

double distanceFromPointToLine(const Point3s& A, const Point3s& B, const Point3s& P){
    // ベクトル AB
    double ABx = B.x - A.x;
    double ABy = B.y - A.y;
    double ABz = B.z - A.z;

    // ベクトル AP
    double APx = P.x - A.x;
    double APy = P.y - A.y;
    double APz = P.z - A.z;

    // ABの大きさ
    double AB_length = max(double(sqrt(ABx * ABx + ABy * ABy + ABz * ABz)), 1e-8);

    // 外積 AB × AP
    double crossX = ABy * APz - ABz * APy;
    double crossY = ABz * APx - ABx * APz;
    double crossZ = ABx * APy - ABy * APx;

    // 外積の大きさ
    double cross_length = max(double(sqrt(crossX * crossX + crossY * crossY + crossZ * crossZ)), 1e-8);

    // 距離を計算
    return cross_length / AB_length;
}

void distanceFromPointToLineGrad(const Point3s& A, const Point3s& B, const Point3s& P, Point3s& gradA, Point3s& gradB, Point3s& gradP, const double &ideal, const double K){
    double ABx = B.x - A.x;
    double ABy = B.y - A.y;
    double ABz = B.z - A.z;

	// ベクトル AP
    double APx = P.x - A.x;
    double APy = P.y - A.y;
    double APz = P.z - A.z;

	double AB_length_sq = ABx * ABx + ABy * ABy + ABz * ABz;

	if(AB_length_sq < 1e-16) return;

	double AB_length = max(double(sqrt(AB_length_sq)), 1e-8);

	if (AB_length_sq < 1e-16) return; 

	double AB_dot_AP = ABx * APx + ABy * APy + ABz * APz;

	double t = AB_dot_AP / AB_length_sq;

	Point3s R_perp;
    R_perp.x = APx - t * ABx;
    R_perp.y = APy - t * ABy;
    R_perp.z = APz - t * ABz;

	double d_sq = R_perp.x * R_perp.x + R_perp.y * R_perp.y + R_perp.z * R_perp.z;

	if (d_sq < 1e-8) return;

	double d = sqrt(d_sq);

	//cout<<"0 " << d << "-" << ideal<<endl;

	double diff_dist = d - ideal;
	double common_factor = 2.0 * K * diff_dist;

	double inv_d = 1.0 / d;
    Point3s u_d;
    u_d.x = R_perp.x * inv_d;
    u_d.y = R_perp.y * inv_d;
    u_d.z = R_perp.z * inv_d;

	double u_dot_v = u_d.x * ABx + u_d.y * ABy + u_d.z * ABz;
    double t_prime = u_dot_v / AB_length_sq;

	Point3s Grad_d_p;
    Grad_d_p.x = u_d.x - t_prime * ABx;
    Grad_d_p.y = u_d.y - t_prime * ABy;
    Grad_d_p.z = u_d.z - t_prime * ABz;

	Point3s Grad_d_n;
    Grad_d_n.x = (t - 1.0) * u_d.x + t_prime * ABx;
    Grad_d_n.y = (t - 1.0) * u_d.y + t_prime * ABy;
    Grad_d_n.z = (t - 1.0) * u_d.z + t_prime * ABz;

	Point3s Grad_d_c;
    Grad_d_c.x = -t * u_d.x + t_prime * ABx;
    Grad_d_c.y = -t * u_d.y + t_prime * ABy;
    Grad_d_c.z = -t * u_d.z + t_prime * ABz;

	// P_p (teAtoms[10]) の勾配
    gradP.x += common_factor * Grad_d_p.x;
    gradP.y += common_factor * Grad_d_p.y;
    gradP.z += common_factor * Grad_d_p.z;

    // P_n (teAtoms[8]) の勾配
    gradA.x += common_factor * Grad_d_n.x;
    gradA.y += common_factor * Grad_d_n.y;
    gradA.z += common_factor * Grad_d_n.z;

    // P_c (teAtoms[5]) の勾配
    gradB.x += common_factor * Grad_d_c.x;
    gradB.y += common_factor * Grad_d_c.y;
    gradB.z += common_factor * Grad_d_c.z;

	//cout<<common_factor * Grad_d_n.x << "," << common_factor * Grad_d_n.y << "," << common_factor * Grad_d_n.z <<endl;
	//cout<<common_factor * Grad_d_c.x << "," << common_factor * Grad_d_c.y << "," << common_factor * Grad_d_c.z <<endl;
	//cout<<common_factor * Grad_d_p.x << "," << common_factor * Grad_d_p.y << "," << common_factor * Grad_d_p.z <<endl;
}

double MultiBondDihedralobj::calcenergymultidihedral()  const
{
	int atom1, atom2, atom3, atom4;
	Point3s teAtomA, teAtomB, teAtomC, teAtomD;
	Point3s vec12, vec23, vec34, vecnormal123, vecnormal234, outerprot, a, b, grad1, grad2, grad3, grad4;
	Point3s vec123, vec234, outerprot23a, outerprota123, outerprot34b, outerprotb234, outerprot12a, outerprot23b;

	vector<double> torsions(torsion_atom_idx.size(),0.0);
	double baseplengths;

	for(int tor=0;tor<torsion_atom_idx.size();tor++){
		if(tor == 1) continue;
		atom1 = torsion_atom_idx[tor][0];
		atom2 = torsion_atom_idx[tor][1];
		atom3 = torsion_atom_idx[tor][2];
		atom4 = torsion_atom_idx[tor][3];

		teAtomA.x = teAtoms[atom1]->x; teAtomA.y = teAtoms[atom1]->y; teAtomA.z = teAtoms[atom1]->z;
		teAtomB.x = teAtoms[atom2]->x; teAtomB.y = teAtoms[atom2]->y; teAtomB.z = teAtoms[atom2]->z;
		teAtomC.x = teAtoms[atom3]->x; teAtomC.y = teAtoms[atom3]->y; teAtomC.z = teAtoms[atom3]->z;
		teAtomD.x = teAtoms[atom4]->x; teAtomD.y = teAtoms[atom4]->y; teAtomD.z = teAtoms[atom4]->z;

		vec12.x = teAtomB.x - teAtomA.x;
		vec12.y = teAtomB.y - teAtomA.y;
		vec12.z = teAtomB.z - teAtomA.z;
		
		vec23.x = teAtomC.x - teAtomB.x;
		vec23.y = teAtomC.y - teAtomB.y;
		vec23.z = teAtomC.z - teAtomB.z;

		vec34.x = teAtomD.x - teAtomC.x;
		vec34.y = teAtomD.y - teAtomC.y;
		vec34.z = teAtomD.z - teAtomC.z;

		vecnormal123 = basecross(vec12,vec23);
		vecnormal234 = basecross(vec23,vec34);

		double length123=sqrt(max(double(vecnormal123.x*vecnormal123.x+vecnormal123.y*vecnormal123.y+vecnormal123.z*vecnormal123.z),1e-8));
		double length234=sqrt(max(double(vecnormal234.x*vecnormal234.x+vecnormal234.y*vecnormal234.y+vecnormal234.z*vecnormal234.z),1e-8));

		vecnormal123.x /= length123; vecnormal123.y /= length123; vecnormal123.z /= length123;
		vecnormal234.x /= length234; vecnormal234.y /= length234; vecnormal234.z /= length234;

		double cosphi=vecnormal123.x*vecnormal234.x+vecnormal123.y*vecnormal234.y+vecnormal123.z*vecnormal234.z;

		if(cosphi>1.) cosphi=1.;
		else if(cosphi<-1.) cosphi=-1.;
		
		outerprot=basecross(vecnormal123,vecnormal234);
		double dval = acos(cosphi) * ((vec23.x*outerprot.x+vec23.y*outerprot.y+vec23.z*outerprot.z) >= 0. ? 1.:-1.);

		if (dval < 0.0) dval += 2.0 * M_PI;

		torsions[tor] = dval;
	}

	Point3s ptn, ptc, ptp;
	ptn.x = teAtoms[8]->x;  ptn.y = teAtoms[8]->y;  ptn.z = teAtoms[8]->z;
	ptc.x = teAtoms[5]->x;  ptc.y = teAtoms[5]->y;  ptc.z = teAtoms[5]->z;
	ptp.x = teAtoms[10]->x; ptp.y = teAtoms[10]->y; ptp.z = teAtoms[10]->z;
		
	baseplengths = distanceFromPointToLine(ptn, ptc, ptp);

	double minenergy = 1.0/0.0;
	double diff;
	double idealval;
	for(int i=0;i<ideals.size();i++){
		double energy = 0;

		diff = baseplengths - (*ideals[i][0][0]);

		//cout<<baseplengths << "-" << (*ideals[i][0][0]) << "-" << diff<<endl;

		energy += (K/100) * diff * diff;

		for(int j=1;j<ideals[i].size();j++){
			if(j == 2) continue;

			diff = 1 + cos((*ideals[i][j][0])*torsions[j-1] - (*ideals[i][j][1])*RADDEG);

			//double s  = sin(diff);
			//double c  = cos(diff);
			//double delta = atan2(s, c); 
			//energy += K * (1.0 - cos(delta));
			energy += K * diff;
		}

		if(minenergy > energy) minenergy = energy;
	}

	//cout<<minenergy<<endl;
	
	return minenergy;	
}

void MultiBondDihedralobj::calcgradientmultidihedral()
{
	int atom1, atom2, atom3, atom4;
	Point3s teAtomA, teAtomB, teAtomC, teAtomD;
	Point3s vec12, vec23, vec34;
	Point3s vecnormal123, vecnormal234, outerprot, a, b; 
	Point3s outerprot23a, outerprot12a, outerprot34d, outerprot23b;
	double baseplengths;

	vector<double> torsions(torsion_atom_idx.size(),0.0);

	for(int tor=0;tor<torsion_atom_idx.size();tor++)
	{
		if(tor == 1) continue;
		atom1 = torsion_atom_idx[tor][0];
		atom2 = torsion_atom_idx[tor][1];
		atom3 = torsion_atom_idx[tor][2];
		atom4 = torsion_atom_idx[tor][3];

		teAtomA.x = teAtoms[atom1]->x; teAtomA.y = teAtoms[atom1]->y; teAtomA.z = teAtoms[atom1]->z;
		teAtomB.x = teAtoms[atom2]->x; teAtomB.y = teAtoms[atom2]->y; teAtomB.z = teAtoms[atom2]->z;
		teAtomC.x = teAtoms[atom3]->x; teAtomC.y = teAtoms[atom3]->y; teAtomC.z = teAtoms[atom3]->z;
		teAtomD.x = teAtoms[atom4]->x; teAtomD.y = teAtoms[atom4]->y; teAtomD.z = teAtoms[atom4]->z;

		vec12.x = teAtomB.x - teAtomA.x;
		vec12.y = teAtomB.y - teAtomA.y;
		vec12.z = teAtomB.z - teAtomA.z;
		
		vec23.x = teAtomC.x - teAtomB.x;
		vec23.y = teAtomC.y - teAtomB.y;
		vec23.z = teAtomC.z - teAtomB.z;

		vec34.x = teAtomD.x - teAtomC.x;
		vec34.y = teAtomD.y - teAtomC.y;
		vec34.z = teAtomD.z - teAtomC.z;

		vecnormal123 = basecross(vec12,vec23);
		vecnormal234 = basecross(vec23,vec34);

		double length123=sqrt(max(double(vecnormal123.x*vecnormal123.x+vecnormal123.y*vecnormal123.y+vecnormal123.z*vecnormal123.z),1e-8));
		double length234=sqrt(max(double(vecnormal234.x*vecnormal234.x+vecnormal234.y*vecnormal234.y+vecnormal234.z*vecnormal234.z),1e-8));

		vecnormal123.x /= length123; vecnormal123.y /= length123; vecnormal123.z /= length123;
		vecnormal234.x /= length234; vecnormal234.y /= length234; vecnormal234.z /= length234;

		double innerprot=vecnormal123.x*vecnormal234.x+vecnormal123.y*vecnormal234.y+vecnormal123.z*vecnormal234.z;

		if(innerprot>1.) innerprot=1.; 
		else if(innerprot<-1.) innerprot=-1.; 

		outerprot=basecross(vecnormal123,vecnormal234);

		double dval = acos(innerprot) * ((vec23.x*outerprot.x+vec23.y*outerprot.y+vec23.z*outerprot.z) >= 0. ? 1.:-1.);

		if (dval < 0.0) dval += 2.0 * M_PI;

		torsions[tor] = dval;
	}

	Point3s ptn, ptc, ptp;
	ptn.x = teAtoms[8]->x;  ptn.y = teAtoms[8]->y;  ptn.z = teAtoms[8]->z;
	ptc.x = teAtoms[5]->x;  ptc.y = teAtoms[5]->y;  ptc.z = teAtoms[5]->z;
	ptp.x = teAtoms[10]->x; ptp.y = teAtoms[10]->y; ptp.z = teAtoms[10]->z;

	//cout<<"================="<<endl;
	//cout<<ptn.x << "," << ptn.y << "," << ptn.z <<endl;
	//cout<<ptc.x << "," << ptc.y << "," << ptc.z <<endl;
	//cout<<ptp.x << "," << ptp.y << "," << ptp.z <<endl;
		
	baseplengths = distanceFromPointToLine(ptn, ptc, ptp);

	//cout<<baseplengths<<endl;

	double minenergy = 1.0/0.0;
	int minidx = -1;
	double diff;

	for(int i=0;i<ideals.size();i++)
	{
		double energy = 0;

		diff = baseplengths - (*ideals[i][0][0]);

		//cout<<"aaa-" << (K/100) * diff * diff<<endl;

		energy += (K/100) * diff * diff;

		for(int j=1;j<ideals[i].size();j++)
		{
			if(j == 2) continue;
			//cout<<(*ideals[i][j][1])*RADDEG<<endl;
			diff = 1 + cos((*ideals[i][j][0])*torsions[j-1] - (*ideals[i][j][1])*RADDEG);

			//cout<<"bbb-" << diff<<endl;

			energy += K * diff;
		}

		if(minenergy > energy){
			minenergy = energy;
			minidx = i;
		}
	}
	//exit(1);

	//cout<<minenergy<<endl;

	if (minidx == -1) return;
	
	vector<vector<double*>> minideals = ideals[minidx];

	//cout<<torsion_atom_idx.size() << "================"<<endl;
	for(int tor=0;tor<torsion_atom_idx.size();tor++){
		if(tor == 1) continue;
		atom1 = torsion_atom_idx[tor][0];
		atom2 = torsion_atom_idx[tor][1];
		atom3 = torsion_atom_idx[tor][2];
		atom4 = torsion_atom_idx[tor][3];

		teAtomA.x = teAtoms[atom1]->x; teAtomA.y = teAtoms[atom1]->y; teAtomA.z = teAtoms[atom1]->z;
		teAtomB.x = teAtoms[atom2]->x; teAtomB.y = teAtoms[atom2]->y; teAtomB.z = teAtoms[atom2]->z;
		teAtomC.x = teAtoms[atom3]->x; teAtomC.y = teAtoms[atom3]->y; teAtomC.z = teAtoms[atom3]->z;
		teAtomD.x = teAtoms[atom4]->x; teAtomD.y = teAtoms[atom4]->y; teAtomD.z = teAtoms[atom4]->z;

		vec12.x = teAtomB.x - teAtomA.x;
		vec12.y = teAtomB.y - teAtomA.y;
		vec12.z = teAtomB.z - teAtomA.z;
		
		vec23.x = teAtomC.x - teAtomB.x;
		vec23.y = teAtomC.y - teAtomB.y;
		vec23.z = teAtomC.z - teAtomB.z;

		vec34.x = teAtomD.x - teAtomC.x;
		vec34.y = teAtomD.y - teAtomC.y;
		vec34.z = teAtomD.z - teAtomC.z;

		vecnormal123 = basecross(vec12,vec23);
		vecnormal234 = basecross(vec23,vec34);

		double length123=sqrt(max(double(vecnormal123.x*vecnormal123.x+vecnormal123.y*vecnormal123.y+vecnormal123.z*vecnormal123.z),1e-8));
		double length234=sqrt(max(double(vecnormal234.x*vecnormal234.x+vecnormal234.y*vecnormal234.y+vecnormal234.z*vecnormal234.z),1e-8));

		vecnormal123.x /= length123; vecnormal123.y /= length123; vecnormal123.z /= length123;
		vecnormal234.x /= length234; vecnormal234.y /= length234; vecnormal234.z /= length234;

		double innerprot=vecnormal123.x*vecnormal234.x+vecnormal123.y*vecnormal234.y+vecnormal123.z*vecnormal234.z;

		if(innerprot>1.) innerprot=1.; 
		else if(innerprot<-1.) innerprot=-1.; 

		outerprot=basecross(vecnormal123,vecnormal234);

		double dval = acos(innerprot) * ((vec23.x*outerprot.x+vec23.y*outerprot.y+vec23.z*outerprot.z) >= 0. ? 1.:-1.);

		double sinphi = sin(dval);
    
		if(fabs(sinphi) < 1e-8 || isnan(sinphi)) continue;

		//cout<<tor+1 << "-" << minideals.size()<<endl;
		double n = (*minideals[tor+1][0]);
        double phi0 = (*minideals[tor+1][1])*RADDEG;
        double diff = n * dval - phi0;

		//double diff = dval - (*minideals[tor])*RADDEG;
		//double s  = sin(diff);
		//double c  = cos(diff);
		//double delta = atan2(s, c);

		//double ene = K * (1.0 - cos(delta));
		double gradient_E_to_phi = K * n * sin(diff);

		//cout<<"---------"<<endl;
		//cout<< tor+1 << "-" << 1 + cos(diff) << "-" << dval*DEGRAD << "-" << phi0 * DEGRAD<<endl;
		//cout<<"---------"<<endl;

		//double gradient_scale = -(K * sin(delta)) / sinphi;
		double gradient_scale = gradient_E_to_phi / sinphi;

		double inv_len123_sq = 1.0 / sqrt(length123 * length123);
		double inv_len234_sq = 1.0 / sqrt(length234 * length234);
		//double inv_len123_sq = 1.0 / length123 * length123;
		//double inv_len234_sq = 1.0 / length234 * length234;

		a.x = gradient_scale * (vecnormal234.x - vecnormal123.x * innerprot)*inv_len123_sq;
		a.y = gradient_scale * (vecnormal234.y - vecnormal123.y * innerprot)*inv_len123_sq;
		a.z = gradient_scale * (vecnormal234.z - vecnormal123.z * innerprot)*inv_len123_sq;

		b.x = gradient_scale * (vecnormal123.x - vecnormal234.x * innerprot)*inv_len234_sq;
		b.y = gradient_scale * (vecnormal123.y - vecnormal234.y * innerprot)*inv_len234_sq;
		b.z = gradient_scale * (vecnormal123.z - vecnormal234.z * innerprot)*inv_len234_sq;

		outerprot12a = basecross(a,vec12);     // vec(A->B) x a
		outerprot23a = basecross(vec23,a);     // vec(B->C) x a
		outerprot23b = basecross(b,vec23);     // vec(B->C) x b
		outerprot34d = basecross(vec34,b);     // vec(C->D) x b

		teAtomgrad[atom1]->x += (-outerprot23a.x); 
		teAtomgrad[atom1]->y += (-outerprot23a.y); 
		teAtomgrad[atom1]->z += (-outerprot23a.z);

		teAtomgrad[atom2]->x += (-outerprot34d.x-outerprot12a.x+outerprot23a.x); 
		teAtomgrad[atom2]->y += (-outerprot34d.y-outerprot12a.y+outerprot23a.y);
		teAtomgrad[atom2]->z += (-outerprot34d.z-outerprot12a.z+outerprot23a.z); 

		teAtomgrad[atom3]->x += (outerprot12a.x+outerprot34d.x-outerprot23b.x);
		teAtomgrad[atom3]->y += (outerprot12a.y+outerprot34d.y-outerprot23b.y);
		teAtomgrad[atom3]->z += (outerprot12a.z+outerprot34d.z-outerprot23b.z); 

		teAtomgrad[atom4]->x += (outerprot23b.x); 
		teAtomgrad[atom4]->y += (outerprot23b.y); 
		teAtomgrad[atom4]->z += (outerprot23b.z);
	}

	//cout<<teAtomgrad[8]->x << "," << teAtomgrad[8]->y << "," << teAtomgrad[8]->z <<endl;
	//cout<<teAtomgrad[5]->x << "," << teAtomgrad[5]->y << "," << teAtomgrad[5]->z <<endl;
	//cout<<teAtomgrad[10]->x << "," << teAtomgrad[10]->y << "," << teAtomgrad[10]->z <<endl;
	distanceFromPointToLineGrad(ptn, ptc, ptp, *teAtomgrad[8], *teAtomgrad[5], *teAtomgrad[10], (*minideals[0][0]), (K/100));

	//cout<<"=================="<<endl;

	//cout<<teAtomgrad[8]->x << "," << teAtomgrad[8]->y << "," << teAtomgrad[8]->z <<endl;
	//cout<<teAtomgrad[5]->x << "," << teAtomgrad[5]->y << "," << teAtomgrad[5]->z <<endl;
	//cout<<teAtomgrad[10]->x << "," << teAtomgrad[10]->y << "," << teAtomgrad[10]->z <<endl;

	//exit(1);
}

Point3s BondImpropersobj::basecross(Point3s p1, Point3s p2)  const
{
	Point3s temp;
	temp.x=p1.y * p2.z - p1.z * p2.y;
	temp.y=p1.z * p2.x - p1.x * p2.z;
	temp.z=p1.x * p2.y - p1.y * p2.x;
	return temp;
}

double BondImpropersobj::calcenergyimpropers() const
{
	Point3s vec12, vec23, vec34, vecnormal123, vecnormal234, outerprot, a, b, grad1, grad2, grad3, grad4;
	Point3s vec123, vec234, outerprot23a, outerprota123, outerprot34b, outerprotb234, outerprot12a, outerprot23b;

	vec12.x = teAtomA->x - teAtomC->x;
	vec12.y = teAtomA->y - teAtomC->y;
	vec12.z = teAtomA->z - teAtomC->z;
	
	vec23.x = teAtomB->x - teAtomA->x;
	vec23.y = teAtomB->y - teAtomA->y;
	vec23.z = teAtomB->z - teAtomA->z;

	vec34.x = teAtomD->x - teAtomB->x;
	vec34.y = teAtomD->y - teAtomB->y;
	vec34.z = teAtomD->z - teAtomB->z;

	vecnormal123 = basecross(vec12,vec23);
	vecnormal234 = basecross(vec23,vec34);

	double length123=sqrt(vecnormal123.x*vecnormal123.x+vecnormal123.y*vecnormal123.y+vecnormal123.z*vecnormal123.z+1e-8);
	double length234=sqrt(vecnormal234.x*vecnormal234.x+vecnormal234.y*vecnormal234.y+vecnormal234.z*vecnormal234.z+1e-8);

	vecnormal123.x /= length123; vecnormal123.y /= length123; vecnormal123.z /= length123;
	vecnormal234.x /= length234; vecnormal234.y /= length234; vecnormal234.z /= length234;

	double cosphi=vecnormal123.x*vecnormal234.x+vecnormal123.y*vecnormal234.y+vecnormal123.z*vecnormal234.z;
	
	if(cosphi>1.) cosphi=1.;
  	else if(cosphi<-1.) cosphi=-1.;
	
	outerprot=basecross(vecnormal123,vecnormal234);
	double dval = acos(cosphi) * ((vec23.x*outerprot.x+vec23.y*outerprot.y+vec23.z*outerprot.z) >= 0. ? 1.:-1.);
	
	return KI*(1+cos(PN*dval-Phase*RADDEG));
}

void BondImpropersobj::calcgradientimpropers() const
{
	Point3s vec12, vec23, vec34, vecnormal123, vecnormal234, outerprot, a, b, grad1, grad2, grad3, grad4;
	Point3s vec123, vec234, outerprot23a, outerprota123, outerprot34b, outerprotb234, outerprot12a, outerprot23b;

	vec12.x = teAtomA->x - teAtomC->x;
	vec12.y = teAtomA->y - teAtomC->y;
	vec12.z = teAtomA->z - teAtomC->z;
	
	vec23.x = teAtomB->x - teAtomA->x;
	vec23.y = teAtomB->y - teAtomA->y;
	vec23.z = teAtomB->z - teAtomA->z;

	vec34.x = teAtomD->x - teAtomB->x;
	vec34.y = teAtomD->y - teAtomB->y;
	vec34.z = teAtomD->z - teAtomB->z;

	vecnormal123 = basecross(vec12,vec23);
	vecnormal234 = basecross(vec23,vec34);

	double length123=sqrt(vecnormal123.x*vecnormal123.x+vecnormal123.y*vecnormal123.y+vecnormal123.z*vecnormal123.z+1e-8);
	double length234=sqrt(vecnormal234.x*vecnormal234.x+vecnormal234.y*vecnormal234.y+vecnormal234.z*vecnormal234.z+1e-8);

	vecnormal123.x /= length123; vecnormal123.y /= length123; vecnormal123.z /= length123;
	vecnormal234.x /= length234; vecnormal234.y /= length234; vecnormal234.z /= length234;

	double innerprot=vecnormal123.x*vecnormal234.x+vecnormal123.y*vecnormal234.y+vecnormal123.z*vecnormal234.z;
	outerprot=basecross(vecnormal123,vecnormal234);

	double phi = acos(innerprot) * ((vec23.x*outerprot.x+vec23.y*outerprot.y+vec23.z*outerprot.z) >= 0. ? 1.:-1.);
	double sinphi = sin(phi);

	if(fabs(sinphi) < 0.00001 || isnan(sinphi)) return;

	double energydepth = KI*PN*sin(phi*PN-Phase*RADDEG)/sinphi;
	double cosdivleng123 = cos(phi)/length123;
	double cosdivleng234 = cos(phi)/length234;

	a.x = -vecnormal234.x+vecnormal123.x*cosdivleng123;
	a.y = -vecnormal234.y+vecnormal123.y*cosdivleng123;
	a.z = -vecnormal234.z+vecnormal123.z*cosdivleng123;

	b.x = -vecnormal123.x+vecnormal234.x*cosdivleng234;
	b.y = -vecnormal123.y+vecnormal234.y*cosdivleng234;
	b.z = -vecnormal123.z+vecnormal234.z*cosdivleng234;

	vec123.x = vec12.x+vec23.x; vec123.y = vec12.y+vec23.y; vec123.z = vec12.z+vec23.z;
	vec234.x = vec23.x+vec34.x; vec234.y = vec23.y+vec34.y; vec234.z = vec23.z+vec34.z;

	outerprot23a=basecross(vec23,a);
	outerprota123=basecross(a,vec123);
	outerprot34b=basecross(vec34,b);
	outerprotb234=basecross(b,vec234);
	outerprot12a=basecross(vec12,a);
	outerprot23b=basecross(vec23,b);

	teAtomAgrad->x += energydepth*outerprot23a.x; 
	teAtomAgrad->y += energydepth*outerprot23a.y; 
	teAtomAgrad->z += energydepth*outerprot23a.z;

	teAtomBgrad->x += energydepth*(outerprota123.x+outerprot34b.x);
	teAtomBgrad->y += energydepth*(outerprota123.y+outerprot34b.y);
	teAtomBgrad->z += energydepth*(outerprota123.z+outerprot34b.z);

	teAtomCgrad->x += energydepth*(outerprotb234.x+outerprot12a.x);
	teAtomCgrad->y += energydepth*(outerprotb234.y+outerprot12a.y);
	teAtomCgrad->z += energydepth*(outerprotb234.z+outerprot12a.z);
	
	teAtomDgrad->x += energydepth*outerprot23b.x;
	teAtomDgrad->y += energydepth*outerprot23b.y;
	teAtomDgrad->z += energydepth*outerprot23b.z;
}

void Hbondobj::calcgradienthbond() const
{
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
	vec13.x /= length13; vec13.y /= length13; vec13.z /= length13;

	double innerprot = vec21.x*vec23.x+vec21.y*vec23.y+vec21.z*vec23.z;
	double dcosangle = innerprot;
	if     (dcosangle >  1.) dcosangle =  1.;
  	else if(dcosangle < -1.) dcosangle = -1.;

	double dcosangle_diff = dcosangle - idealahdang;
	double dcosangle_diff_beta = beta*dcosangle_diff;
	double dcosangle_softplus = (1/beta)*log1p(exp(dcosangle_diff_beta));

	double idealahdiff = length23-idealahleng;
	double idealaddiff = length13-idealadleng;

	double denergydepth = -dK*exp(-dXiR * (idealahdiff*idealahdiff + idealaddiff*idealaddiff) - dXiTheta*dcosangle_softplus*dcosangle_softplus + dREq);

	//cout<<"=============="<<endl;
	//cout<<-dXiR * (idealahdiff*idealahdiff + idealaddiff*idealaddiff) - dXiTheta*dcosangle_softplus*dcosangle_softplus + dREq<<endl;
	//cout<<length13 << "-" << idealadleng<<endl;
	//cout<<length23 << "-" << idealahleng<<endl;
	//cout<<dcosangle << "-" << idealahdang<<endl;

	const double angle_term_deriv = -dXiTheta * 2.0 *  dcosangle_softplus * (1.0 / (1.0 + exp(-dcosangle_diff_beta)));

	double coeff1 = angle_term_deriv / (length21 + 1e-8);
	double coeff2 = angle_term_deriv / (length23 + 1e-8);
	double coeff3 = innerprot-1.0;
	double coeff4 = -2.*dXiR;

	grad1.x = coeff1*(-dcosangle*vec21.x + vec23.x) + coeff4*idealaddiff*vec13.x;
	grad1.y = coeff1*(-dcosangle*vec21.y + vec23.y) + coeff4*idealaddiff*vec13.y;
	grad1.z = coeff1*(-dcosangle*vec21.z + vec23.z) + coeff4*idealaddiff*vec13.z;

	grad3.x = coeff4*idealahdiff*vec23.x + coeff2*(-dcosangle*vec23.x + vec21.x) - coeff4*idealaddiff*vec13.x;
	grad3.y = coeff4*idealahdiff*vec23.y + coeff2*(-dcosangle*vec23.y + vec21.y) - coeff4*idealaddiff*vec13.y;
	grad3.z = coeff4*idealahdiff*vec23.z + coeff2*(-dcosangle*vec23.z + vec21.z) - coeff4*idealaddiff*vec13.z;

	grad2.x = -coeff1*(-dcosangle*vec21.x + vec23.x)-coeff2*(-dcosangle*vec23.x + vec21.x)-coeff4*idealahdiff*vec23.x;
	grad2.y = -coeff1*(-dcosangle*vec21.y + vec23.y)-coeff2*(-dcosangle*vec23.y + vec21.y)-coeff4*idealahdiff*vec23.y;
	grad2.z = -coeff1*(-dcosangle*vec21.z + vec23.z)-coeff2*(-dcosangle*vec23.z + vec21.z)-coeff4*idealahdiff*vec23.z;

	gradsum.x = grad1.x + grad2.x + grad3.x;
	gradsum.y = grad1.y + grad2.y + grad3.y;
	gradsum.z = grad1.z + grad2.z + grad3.z;

	double gradsumleng = sqrt(gradsum.x*gradsum.x+gradsum.y*gradsum.y+gradsum.z*gradsum.z+1e-8);

	if(isnan(gradsumleng)) return;

	teAtomDongrad->x += denergydepth*grad1.x;
	teAtomDongrad->y += denergydepth*grad1.y;
	teAtomDongrad->z += denergydepth*grad1.z;
	teAtomHydgrad->x += denergydepth*grad2.x;
	teAtomHydgrad->y += denergydepth*grad2.y;
	teAtomHydgrad->z += denergydepth*grad2.z;
	teAtomAccgrad->x += denergydepth*grad3.x;
	teAtomAccgrad->y += denergydepth*grad3.y;
	teAtomAccgrad->z += denergydepth*grad3.z;
}

double Hbondobj::calcenergyhbond() const
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

	double length21=max(double(sqrt(vec21.x*vec21.x+vec21.y*vec21.y+vec21.z*vec21.z)),1e-8);
	double length23=max(double(sqrt(vec23.x*vec23.x+vec23.y*vec23.y+vec23.z*vec23.z)),1e-8);
	double length13=max(double(sqrt(vec13.x*vec13.x+vec13.y*vec13.y+vec13.z*vec13.z)),1e-8);

	vec21.x /= length21; vec21.y /= length21; vec21.z /= length21;
	vec23.x /= length23; vec23.y /= length23; vec23.z /= length23;
	
	double innerprot = vec21.x*vec23.x+vec21.y*vec23.y+vec21.z*vec23.z;
	double dcosangle = innerprot;
	if     (dcosangle >  1.) dcosangle =  1.;
  	else if(dcosangle < -1.) dcosangle = -1.;

	double dcosangle_diff =  dcosangle-idealahdang;
	double dcosangle_softplus = (1/beta)*log1p(exp(beta*dcosangle_diff));
    
	double idealahdiff = length23-idealahleng;
	double idealaddiff = length13-idealadleng;

	return -dK*exp(-dXiR * (idealahdiff*idealahdiff + idealaddiff*idealaddiff) - dXiTheta*dcosangle_softplus*dcosangle_softplus + dREq);
}

void Hbond2hydobj::calcgradienthbond() const
{
	Point3s vec23_1, vec23_2, vec23, vec13, vec21, grad1, grad2, grad3, gradsum;

	vec23_1.x = teAtomAcc->x - teAtomHyd1->x;
	vec23_1.y = teAtomAcc->y - teAtomHyd1->y;
	vec23_1.z = teAtomAcc->z - teAtomHyd1->z;

	vec23_2.x = teAtomAcc->x - teAtomHyd2->x;
	vec23_2.y = teAtomAcc->y - teAtomHyd2->y;
	vec23_2.z = teAtomAcc->z - teAtomHyd2->z;

	double length23_1=max(double(sqrt(vec23_1.x*vec23_1.x+vec23_1.y*vec23_1.y+vec23_1.z*vec23_1.z)),1e-8);
	double length23_2=max(double(sqrt(vec23_2.x*vec23_2.x+vec23_2.y*vec23_2.y+vec23_2.z*vec23_2.z)),1e-8);

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

	//bool angle_active = (dcosangle >= idealahdang);
	//if (!angle_active) dcosangle = -1.0;

	double dcosangle_diff = dcosangle - idealahdang;
	double dcosangle_diff_beta = beta*dcosangle_diff;
	double dcosangle_softplus = (1/beta)*log1p(exp(dcosangle_diff_beta));
    //double dcosangle_diff_sq = dcosangle_diff * dcosangle_diff;

	//cout<<dK << "-" << dcosangle << "-" <<  idealahdang << "-" << dcosangle_softplus << "-" << length23 << "-" << idealahleng << "-" << length13 << "-" << idealadleng<<endl;

	double idealahdiff = length23-idealahleng;
	double idealaddiff = length13-idealadleng;

	double denergydepth = -dK*exp(-dXiR * (idealahdiff*idealahdiff + idealaddiff*idealaddiff) - dXiTheta*dcosangle_softplus*dcosangle_softplus + dREq);

	const double angle_term_deriv = -dXiTheta * 2.0 *  dcosangle_softplus * (1.0 / (1.0 + exp(-dcosangle_diff_beta)));

	//double coeff1 = angle_active ? ( -dXiTheta/(length21 + 1e-8) ) : 0.0;
	double coeff1 = angle_term_deriv / (length21 + 1e-8);
	//double coeff2 = angle_active ? ( -dXiTheta/(length23 + 1e-8) ) : 0.0;
	double coeff2 = angle_term_deriv / (length23 + 1e-8);
	double coeff3 = innerprot-1.0;
	double coeff4 = -2.*dXiR;

	grad1.x = coeff1*(-dcosangle*vec21.x + vec23.x) + coeff4*idealaddiff*vec13.x;
	grad1.y = coeff1*(-dcosangle*vec21.y + vec23.y) + coeff4*idealaddiff*vec13.y;
	grad1.z = coeff1*(-dcosangle*vec21.z + vec23.z) + coeff4*idealaddiff*vec13.z;

	grad3.x = coeff4*idealahdiff*vec23.x + coeff2*(-dcosangle*vec23.x + vec21.x) - coeff4*idealaddiff*vec13.x;
	grad3.y = coeff4*idealahdiff*vec23.y + coeff2*(-dcosangle*vec23.y + vec21.y) - coeff4*idealaddiff*vec13.y;
	grad3.z = coeff4*idealahdiff*vec23.z + coeff2*(-dcosangle*vec23.z + vec21.z) - coeff4*idealaddiff*vec13.z;

	grad2.x = -coeff1*(-dcosangle*vec21.x + vec23.x)-coeff2*(-dcosangle*vec23.x + vec21.x)-coeff4*idealahdiff*vec23.x;
	grad2.y = -coeff1*(-dcosangle*vec21.y + vec23.y)-coeff2*(-dcosangle*vec23.y + vec21.y)-coeff4*idealahdiff*vec23.y;
	grad2.z = -coeff1*(-dcosangle*vec21.z + vec23.z)-coeff2*(-dcosangle*vec23.z + vec21.z)-coeff4*idealahdiff*vec23.z;

	gradsum.x = grad1.x + grad2.x + grad3.x;
	gradsum.y = grad1.y + grad2.y + grad3.y;
	gradsum.z = grad1.z + grad2.z + grad3.z;

	double gradsumleng = sqrt(gradsum.x*gradsum.x+gradsum.y*gradsum.y+gradsum.z*gradsum.z+1e-8);

	if(isnan(gradsumleng)) return;

	teAtomDongrad->x += denergydepth*grad1.x;
	teAtomDongrad->y += denergydepth*grad1.y;
	teAtomDongrad->z += denergydepth*grad1.z;
	teAtomAccgrad->x += denergydepth*grad3.x;
	teAtomAccgrad->y += denergydepth*grad3.y;
	teAtomAccgrad->z += denergydepth*grad3.z;

	if(hydatom < 0){
		teAtomHydgrad1->x += denergydepth*grad2.x;
		teAtomHydgrad1->y += denergydepth*grad2.y;
		teAtomHydgrad1->z += denergydepth*grad2.z;
	}else{
		teAtomHydgrad2->x += denergydepth*grad2.x;
		teAtomHydgrad2->y += denergydepth*grad2.y;
		teAtomHydgrad2->z += denergydepth*grad2.z;
	}

}

double Hbond2hydobj::calcenergyhbond() const
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

	//bool angle_active = (dcosangle >= idealahdang);
	//if (!angle_active) dcosangle = -1.0;

	double dcosangle_diff = dcosangle - idealahdang;
	double dcosangle_softplus = (1/beta)*log1p(exp(beta*dcosangle_diff));
    //double dcosangle_diff_sq = dcosangle_diff * dcosangle_diff;

	double idealahdiff = length23-idealahleng;
	double idealaddiff = length13-idealadleng;

	//cout<<-dK*exp((-dXiR * (idealahdiff*idealahdiff + idealaddiff*idealaddiff) + dXiTheta*dcosangle))<<endl;

	return -dK*exp((-dXiR * (idealahdiff*idealahdiff + idealaddiff*idealaddiff) - dXiTheta*dcosangle_softplus*dcosangle_softplus + dREq));
}

/*
bool BondHbondobj::bIsNan(const double _dVal){
    volatile double _dVolatVal = _dVal; 
    return _dVolatVal != _dVolatVal;
}
*/

void BasePairobj::CalcCenterOfMass()
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

double BasePairobj::TotalSQRDistance(double _a, double _b) const
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
double BasePairobj::TotalSQRDistance() const
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


double BasePairobj::Fit3DPlane()
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

BasePairobj::BasePairobj(const vector<Point3s*>& _AtomPtrs, const vector<Point3s*>& _Atomgrad, const double &_K)
{
  if(!_AtomPtrs.empty())
  {
    AtomPtrs = _AtomPtrs;
	Atomgrad = _Atomgrad;
	K = _K;
    n = AtomPtrs.size();
    Fit3DPlane();  
  }  

	/*
    for(unsigned int i = 0; i < n; ++i)
      for(unsigned int j = i; j < n; ++j)
      {
	if(AtomPtrs[i]->sGetResIdent() != AtomPtrs[j]->sGetResIdent() &&
	    Textutils::bIsWCHydrogenBond(AtomPtrs[i]->sGetName(), AtomPtrs[i]->sGetResName(), 
					 AtomPtrs[j]->sGetName(), AtomPtrs[j]->sGetResName()) )
	{
	  BaseSprings.push_back(TSpring(AtomPtrs[i], AtomPtrs[j], BondParam(10., 1.65)));
	}
      }
  }
  else std::cerr << "Internal Warning: Attempt to build an empty pair of bases." << std::endl;
  */
}

double BasePairobj::calcenergybasepair() const
{
  double dEnergy = K*TotalSQRDistance();

  //cout<<dEnergy<<endl;

  //for(std::list<TSpring>::const_iterator I = BaseSprings.begin(), end = BaseSprings.end(); I != end; ++I)
    //dEnergy += I->dCalcEnergy(); 
  
  return (isnan(dEnergy)) ? 0. : dEnergy;
}

void BasePairobj::calcgradientbasepair()
{
  Fit3DPlane();
  
  if(!isnan(a) && !isnan(b))
    for(unsigned int i = 0; i < AtomPtrs.size(); ++i)
    {
      double xk, yk, zk;
	  xk = AtomPtrs[i]->x;
	  yk = AtomPtrs[i]->y;
	  zk = AtomPtrs[i]->z;
      xk -= Xm;		// Reduced coords
      yk -= Ym;
      zk -= Zm;
	  
	  Atomgrad[i]->x += (2.*K*(a*xk + b*yk + zk)/(a*a + b*b + 1.))*a;
	  Atomgrad[i]->y += (2.*K*(a*xk + b*yk + zk)/(a*a + b*b + 1.))*b;
	  Atomgrad[i]->z += (2.*K*(a*xk + b*yk + zk)/(a*a + b*b + 1.));
    }
  
  //for(std::list<TSpring>::const_iterator I = BaseSprings.begin(), end = BaseSprings.end(); I != end; ++I)
    //I->CalcGradients();
}

int TemplateSpuobj::selecttemplate()
{
	vector<vector<double> > xyz_list1(AtomPtrs.size(), vector<double>(3,0));
    vector<vector<double> > xyz_list2(AtomPtrs.size(), vector<double>(3,0));
	vector<vector<double>> RotMatix;
    vector<double> TranVect;
	vector<double> xyztmp(3,0);
	Point3s vec;
	double rmsd;
	double minrmsd=100000000000;
	int minfragidx;

	for(int i=0;i<TempAtomPtrs.size();i++)
	{	
		for(int j=0;j<AtomPtrs.size();j++){
			xyz_list1[j][0] = AtomPtrs[j]->x;
			xyz_list1[j][1] = AtomPtrs[j]->y;
			xyz_list1[j][2] = AtomPtrs[j]->z;

			xyz_list2[j][0] = TempAtomPtrs[i][j]->x;
			xyz_list2[j][1] = TempAtomPtrs[i][j]->y;
			xyz_list2[j][2] = TempAtomPtrs[i][j]->z;
		}

		bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

		rmsd=0;
		for(int j=0;j<AtomPtrs.size();j++){
			bf.Transform(xyz_list2[j],RotMatix,TranVect,xyztmp);

			vec.x=xyz_list1[j][0]-xyztmp[0];
			vec.y=xyz_list1[j][1]-xyztmp[1];
			vec.z=xyz_list1[j][2]-xyztmp[2];

			rmsd+=bf.norm(vec);
		}

		if(minrmsd > rmsd){
            minrmsd=rmsd;
            minfragidx=i;
        }
	}

	return minfragidx;
}

TemplateSpuobj::TemplateSpuobj(const vector<Point3s*>& _AtomPtrs, const vector<Point3s*>& _Atomgrad, const vector<vector<Point3s*>> &_TempAtomPtrs, const double &_dK)
{
  if(!_AtomPtrs.empty())
  {
    AtomPtrs = _AtomPtrs;
	Atomgrad = _Atomgrad;
	TempAtomPtrs = _TempAtomPtrs;
	dK = _dK;
	templateidx = selecttemplate();
  }  
}

double Stackingobj::calcenergystack() const
{
	Point3s vec12;

	vec12.x = teAtomB->x - teAtomA->x;
	vec12.y = teAtomB->y - teAtomA->y;
	vec12.z = teAtomB->z - teAtomA->z;
	
	double length12=sqrt(vec12.x*vec12.x+vec12.y*vec12.y+vec12.z*vec12.z+1e-8);

	double dist = length12 - ideal;

	return KI*dist*dist;
}

void Stackingobj::calcgradientstack() const
{
	Point3s vec12;

	vec12.x = teAtomB->x - teAtomA->x;
	vec12.y = teAtomB->y - teAtomA->y;
	vec12.z = teAtomB->z - teAtomA->z;

	double length12=sqrt(vec12.x*vec12.x+vec12.y*vec12.y+vec12.z*vec12.z+1e-8)+1e-8;

	vec12.x /= length12;
	vec12.y /= length12;
	vec12.z /= length12;

	double coeff=2.*KI*(length12 - ideal);

	teAtomAgrad->x -= coeff*vec12.x;
	teAtomAgrad->y -= coeff*vec12.y;
	teAtomAgrad->z -= coeff*vec12.z;
	teAtomBgrad->x += coeff*vec12.x;
	teAtomBgrad->y += coeff*vec12.y;
	teAtomBgrad->z += coeff*vec12.z;
}

double BPdistobj::calcenergybp() const
{
	Point3s vec12;

	vec12.x = teAtomB->x - teAtomA->x;
	vec12.y = teAtomB->y - teAtomA->y;
	vec12.z = teAtomB->z - teAtomA->z;
	
	double length12=sqrt(vec12.x*vec12.x+vec12.y*vec12.y+vec12.z*vec12.z+1e-8);

	double dist = length12 - ideal;

	return KI*dist*dist;
}

void BPdistobj::calcgradientbp() const
{
	Point3s vec12;

	vec12.x = teAtomB->x - teAtomA->x;
	vec12.y = teAtomB->y - teAtomA->y;
	vec12.z = teAtomB->z - teAtomA->z;

	double length12=sqrt(vec12.x*vec12.x+vec12.y*vec12.y+vec12.z*vec12.z+1e-8)+1e-8;

	vec12.x /= length12;
	vec12.y /= length12;
	vec12.z /= length12;

	double coeff=2.*KI*(length12 - ideal);

	teAtomAgrad->x -= coeff*vec12.x;
	teAtomAgrad->y -= coeff*vec12.y;
	teAtomAgrad->z -= coeff*vec12.z;
	teAtomBgrad->x += coeff*vec12.x;
	teAtomBgrad->y += coeff*vec12.y;
	teAtomBgrad->z += coeff*vec12.z;
}

StackTemplateobj::StackTemplateobj(FragmentResidue* _ResA, FragmentResidue* _ResB, vector<Point3s>* _ResAgrad, vector<Point3s>* _ResBgrad, const vector<ModelUnit>* _ideal_rna, const vector<int>* _indices, const vector<int>* _orders, const int &_templatenum, double _K)
{
	ResA = _ResA;
	ResB = _ResB;
	ResAgrad = _ResAgrad;
	ResBgrad = _ResBgrad;
	ideal_rna = _ideal_rna;
	indices =  _indices;
	orders = _orders;
	templatenum = _templatenum;
	K = _K;

	supResA.resize(ResA->atoms.size(), vector<double>(3,0));
	supResB.resize(ResB->atoms.size(), vector<double>(3,0));

	xyz_list1.resize(ResA->atoms.size()+ResB->atoms.size(), vector<double>(3,0));
    xyz_list2.resize(ResA->atoms.size()+ResB->atoms.size(), vector<double>(3,0));

	base1 = ResA->restype;
	base2 = ResB->restype;
}

double StackTemplateobj::selecttemplate()
{
	int i,j;
	vector<vector<double>> RotMatix, RotMatixmin;
    vector<double> TranVect, TranVectmin;
	vector<double> xyztmp(3,0);
	vector<double> xyz(3,0);
	Point3s vec;
	double rmsd;
	double minrmsd=100000000000;
	int count, idx, ord, minfragidx, minfragord;
	int atom;
	int base1_temp, base2_temp;

	for(i=0; i<templatenum; i++){
		idx = (*indices)[i];
		ord = (*orders)[i];

		if(ord == 0){
			base1_temp = 0;
			base2_temp = 1;
		}else if(ord == 1){
			base1_temp = 1;
			base2_temp = 0;
		}

		count = 0;
    	for(j=0; j<targetatoms[base1].size();j++){
			atom=targetatoms[base1][j];
			xyz_list1[count][0]=ResA->atoms[atom].coord.x;
			xyz_list1[count][1]=ResA->atoms[atom].coord.y;
			xyz_list1[count][2]=ResA->atoms[atom].coord.z;

			xyz_list2[count][0]=(*ideal_rna)[idx].chains[0].residues[base1_temp].atoms[atom].xyz[0];
			xyz_list2[count][1]=(*ideal_rna)[idx].chains[0].residues[base1_temp].atoms[atom].xyz[1];
			xyz_list2[count][2]=(*ideal_rna)[idx].chains[0].residues[base1_temp].atoms[atom].xyz[2];

			count+=1;
    	}

		for(j=0; j<targetatoms[base2].size();j++){
			atom=targetatoms[base2][j];
			xyz_list1[count][0]=ResB->atoms[atom].coord.x;
			xyz_list1[count][1]=ResB->atoms[atom].coord.y;
			xyz_list1[count][2]=ResB->atoms[atom].coord.z;

			xyz_list2[count][0]=(*ideal_rna)[idx].chains[0].residues[base2_temp].atoms[atom].xyz[0];
			xyz_list2[count][1]=(*ideal_rna)[idx].chains[0].residues[base2_temp].atoms[atom].xyz[1];
			xyz_list2[count][2]=(*ideal_rna)[idx].chains[0].residues[base2_temp].atoms[atom].xyz[2];

			count+=1;
		}

		bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

		rmsd=0;
		count=0;
		for(j=0; j<targetatoms[base1].size();j++){
			atom=targetatoms[base1][j];
			bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);

			vec.x=ResA->atoms[atom].coord.x-xyztmp[0];
			vec.y=ResA->atoms[atom].coord.y-xyztmp[1];
			vec.z=ResA->atoms[atom].coord.z-xyztmp[2];

			rmsd+=bf.norm(vec);
			count+=1;
		}

		for(j=0; j<targetatoms[base2].size();j++){
			atom=targetatoms[base2][j];
			bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);

			vec.x=ResB->atoms[atom].coord.x-xyztmp[0];
			vec.y=ResB->atoms[atom].coord.y-xyztmp[1];
			vec.z=ResB->atoms[atom].coord.z-xyztmp[2];

			rmsd+=bf.norm(vec);
			count+=1;
		}

		if(minrmsd > rmsd){
            minrmsd=rmsd;
            minfragidx=idx;
			minfragord=ord;
			RotMatixmin=RotMatix;
			TranVectmin=TranVect;
        }
	}

	if(minfragord == 0){
		base1_temp = 0;
		base2_temp = 1;
	}else if(minfragord == 1){
		base1_temp = 1;
		base2_temp = 0;
	}

	for(j=0; j<targetatoms[base1].size();j++){
		atom=targetatoms[base1][j];
		xyz[0]=(*ideal_rna)[minfragidx].chains[0].residues[base1_temp].atoms[atom].xyz[0];
		xyz[1]=(*ideal_rna)[minfragidx].chains[0].residues[base1_temp].atoms[atom].xyz[1];
		xyz[2]=(*ideal_rna)[minfragidx].chains[0].residues[base1_temp].atoms[atom].xyz[2];

		bf.Transform(xyz,RotMatixmin,TranVectmin,xyztmp);
		supResA[atom][0] = xyztmp[0];
		supResA[atom][1] = xyztmp[1];
		supResA[atom][2] = xyztmp[2];
	}

	for(j=0; j<targetatoms[base2].size();j++){
		atom=targetatoms[base2][j];
		xyz[0]=(*ideal_rna)[minfragidx].chains[0].residues[base2_temp].atoms[atom].xyz[0];
		xyz[1]=(*ideal_rna)[minfragidx].chains[0].residues[base2_temp].atoms[atom].xyz[1];
		xyz[2]=(*ideal_rna)[minfragidx].chains[0].residues[base2_temp].atoms[atom].xyz[2];

		bf.Transform(xyz,RotMatixmin,TranVectmin,xyztmp);
		supResB[atom][0] = xyztmp[0];
		supResB[atom][1] = xyztmp[1];
		supResB[atom][2] = xyztmp[2];
	}

	return minrmsd;
}

double StackTemplateobj::calcenergystacktemp()
{
	double dist = selecttemplate();
	return K*dist;
}

void StackTemplateobj::calcgradientstacktemp() const
{
	Point3s vec12;
	int atom;

	for(int j=0; j<targetatoms.at(base1).size();j++){
		atom=targetatoms.at(base1)[j];

		(*ResAgrad)[atom].x -= 2.0*K*(supResA[atom][0] - ResA->atoms[atom].coord.x);
		(*ResAgrad)[atom].y -= 2.0*K*(supResA[atom][1] - ResA->atoms[atom].coord.y);
		(*ResAgrad)[atom].z -= 2.0*K*(supResA[atom][2] - ResA->atoms[atom].coord.z);

		//cout<<"=========="<<endl;
		//cout<<supResA[atom][0] << "-" << supResA[atom][1] << "-" << supResA[atom][2]<<endl;
		//cout<<ResA->atoms[atom].coord.x << "-" << ResA->atoms[atom].coord.y << "-" << ResA->atoms[atom].coord.z<<endl;
	}

	for(int j=0; j<targetatoms.at(base2).size();j++){
		atom=targetatoms.at(base2)[j];

		(*ResBgrad)[atom].x -= 2.0*K*(supResB[atom][0] - ResB->atoms[atom].coord.x);
		(*ResBgrad)[atom].y -= 2.0*K*(supResB[atom][1] - ResB->atoms[atom].coord.y);
		(*ResBgrad)[atom].z -= 2.0*K*(supResB[atom][2] - ResB->atoms[atom].coord.z);

		//cout<<"=========="<<endl;
		//cout<<supResB[atom][0] << "-" << supResB[atom][1] << "-" << supResB[atom][2]<<endl;
		//cout<<ResB->atoms[atom].coord.x << "-" << ResB->atoms[atom].coord.y << "-" << ResB->atoms[atom].coord.z<<endl;
	}
}

StackTemplateobj2::StackTemplateobj2(FragmentResidue* _ResA, FragmentResidue* _ResB, vector<Point3s>* _ResAgrad, vector<Point3s>* _ResBgrad, const vector<ModelUnit>* _ideal_rna, const vector<int>* _indices, const vector<int>* _orders, const int &_templatenum, double _K)
{
	ResA = _ResA;
	ResB = _ResB;
	ResAgrad = _ResAgrad;
	ResBgrad = _ResBgrad;
	ideal_rna = _ideal_rna;
	indices =  _indices;
	orders = _orders;
	templatenum = _templatenum;
	K = _K;

	supResA.resize(ResA->atoms.size(), vector<double>(3,0));
	supResB.resize(ResB->atoms.size(), vector<double>(3,0));

	base1 = ResA->restype;
	base2 = ResB->restype;

	xyz_listA1.resize(targetatoms[base1].size(), vector<double>(3,0));
    xyz_listA2.resize(targetatoms[base1].size(), vector<double>(3,0));
	xyz_listA3.resize(targetatoms[base2].size(), vector<double>(3,0));

	xyz_listB1.resize(targetatoms[base2].size(), vector<double>(3,0));
    xyz_listB2.resize(targetatoms[base2].size(), vector<double>(3,0));
	xyz_listB3.resize(targetatoms[base1].size(), vector<double>(3,0));
}

double StackTemplateobj2::selecttemplate()
{
	int i,j;
	vector<vector<double>> RotMatixA, RotMatixB, RotMatixminA, RotMatixminB;
    vector<double> TranVectA, TranVectB, TranVectminA, TranVectminB;
	vector<double> xyztmp(3,0);
	vector<double> xyz(3,0);
	Point3s vec;
	double rmsd;
	double minrmsd=100000000000;
	int count, idx, ord, minfragidx, minfragord;
	int atom;
	int base1_temp, base2_temp;

	for(i=0; i<templatenum; i++){
		idx = (*indices)[i];
		ord = (*orders)[i];

		if(ord == 0){
			base1_temp = 0;
			base2_temp = 1;
		}else if(ord == 1){
			base1_temp = 1;
			base2_temp = 0;
		}

		count = 0;
    	for(j=0; j<targetatoms[base1].size();j++){
			atom=targetatoms[base1][j];

			xyz_listA1[count][0]=ResA->atoms[atom].coord.x;
			xyz_listA1[count][1]=ResA->atoms[atom].coord.y;
			xyz_listA1[count][2]=ResA->atoms[atom].coord.z;

			xyz_listA2[count][0]=(*ideal_rna)[idx].chains[0].residues[base1_temp].atoms[atom].xyz[0];
			xyz_listA2[count][1]=(*ideal_rna)[idx].chains[0].residues[base1_temp].atoms[atom].xyz[1];
			xyz_listA2[count][2]=(*ideal_rna)[idx].chains[0].residues[base1_temp].atoms[atom].xyz[2];

			count+=1;
    	}

		bf.RotateCoor(xyz_listA1,xyz_listA2,RotMatixA,TranVectA);

		rmsd=0;
		count=0;
		for(j=0; j<targetatoms[base2].size();j++)
		{
			atom=targetatoms[base2][j];
			xyz_listA3[count][0]=(*ideal_rna)[idx].chains[0].residues[base2_temp].atoms[atom].xyz[0];
			xyz_listA3[count][1]=(*ideal_rna)[idx].chains[0].residues[base2_temp].atoms[atom].xyz[1];
			xyz_listA3[count][2]=(*ideal_rna)[idx].chains[0].residues[base2_temp].atoms[atom].xyz[2];

			bf.Transform(xyz_listA3[count],RotMatixA,TranVectA,xyztmp);

			vec.x=ResB->atoms[atom].coord.x-xyztmp[0];
			vec.y=ResB->atoms[atom].coord.y-xyztmp[1];
			vec.z=ResB->atoms[atom].coord.z-xyztmp[2];

			rmsd+=bf.norm(vec);
			count+=1;
		}

		count = 0;
    	for(j=0; j<targetatoms[base2].size();j++)
		{
			atom=targetatoms[base2][j];
			xyz_listB1[count][0]=ResB->atoms[atom].coord.x;
			xyz_listB1[count][1]=ResB->atoms[atom].coord.y;
			xyz_listB1[count][2]=ResB->atoms[atom].coord.z;

			xyz_listB2[count][0]=(*ideal_rna)[idx].chains[0].residues[base2_temp].atoms[atom].xyz[0];
			xyz_listB2[count][1]=(*ideal_rna)[idx].chains[0].residues[base2_temp].atoms[atom].xyz[1];
			xyz_listB2[count][2]=(*ideal_rna)[idx].chains[0].residues[base2_temp].atoms[atom].xyz[2];

			count+=1;
    	}

		bf.RotateCoor(xyz_listB1,xyz_listB2,RotMatixB,TranVectB);

		count=0;
		for(j=0; j<targetatoms[base1].size();j++)
		{
			atom=targetatoms[base1][j];
			xyz_listB3[count][0]=(*ideal_rna)[idx].chains[0].residues[base1_temp].atoms[atom].xyz[0];
			xyz_listB3[count][1]=(*ideal_rna)[idx].chains[0].residues[base1_temp].atoms[atom].xyz[1];
			xyz_listB3[count][2]=(*ideal_rna)[idx].chains[0].residues[base1_temp].atoms[atom].xyz[2];

			bf.Transform(xyz_listB3[count],RotMatixB,TranVectB,xyztmp);

			vec.x=ResA->atoms[atom].coord.x-xyztmp[0];
			vec.y=ResA->atoms[atom].coord.y-xyztmp[1];
			vec.z=ResA->atoms[atom].coord.z-xyztmp[2];

			rmsd+=bf.norm(vec);
			count+=1;
		}

		if(minrmsd > rmsd)
		{
            minrmsd=rmsd;
            minfragidx=idx;
			minfragord=ord;
			RotMatixminA=RotMatixA;
			TranVectminA=TranVectA;
			RotMatixminB=RotMatixB;
			TranVectminB=TranVectB;
        }
	}

	if(minfragord == 0){
		base1_temp = 0;
		base2_temp = 1;
	}else if(minfragord == 1){
		base1_temp = 1;
		base2_temp = 0;
	}

	for(j=0; j<targetatoms[base2].size();j++){
		atom=targetatoms[base2][j];
		xyz[0]=(*ideal_rna)[minfragidx].chains[0].residues[base2_temp].atoms[atom].xyz[0];
		xyz[1]=(*ideal_rna)[minfragidx].chains[0].residues[base2_temp].atoms[atom].xyz[1];
		xyz[2]=(*ideal_rna)[minfragidx].chains[0].residues[base2_temp].atoms[atom].xyz[2];

		bf.Transform(xyz,RotMatixminA,TranVectminA,xyztmp);
		supResB[atom][0] = xyztmp[0];
		supResB[atom][1] = xyztmp[1];
		supResB[atom][2] = xyztmp[2];
	}

	for(j=0; j<targetatoms[base1].size();j++){
		atom=targetatoms[base1][j];
		xyz[0]=(*ideal_rna)[minfragidx].chains[0].residues[base1_temp].atoms[atom].xyz[0];
		xyz[1]=(*ideal_rna)[minfragidx].chains[0].residues[base1_temp].atoms[atom].xyz[1];
		xyz[2]=(*ideal_rna)[minfragidx].chains[0].residues[base1_temp].atoms[atom].xyz[2];

		bf.Transform(xyz,RotMatixminB,TranVectminB,xyztmp);
		supResA[atom][0] = xyztmp[0];
		supResA[atom][1] = xyztmp[1];
		supResA[atom][2] = xyztmp[2];
	}

	return minrmsd;
}

double StackTemplateobj2::calcenergystacktemp()
{
	double dist = selecttemplate();
	return K*dist;
}

void StackTemplateobj2::calcgradientstacktemp() const
{
	Point3s vec12;
	int atom;

	for(int j=0; j<targetatoms.at(base1).size();j++){
		atom=targetatoms.at(base1)[j];
		(*ResAgrad)[atom].x += 2.0*K*(ResA->atoms[atom].coord.x - supResA[atom][0]);
		(*ResAgrad)[atom].y += 2.0*K*(ResA->atoms[atom].coord.y - supResA[atom][1]);
		(*ResAgrad)[atom].z += 2.0*K*(ResA->atoms[atom].coord.z - supResA[atom][2]);
	}

	for(int j=0; j<targetatoms.at(base2).size();j++){
		atom=targetatoms.at(base2)[j];
		(*ResBgrad)[atom].x += 2.0*K*(ResB->atoms[atom].coord.x - supResB[atom][0]);
		(*ResBgrad)[atom].y += 2.0*K*(ResB->atoms[atom].coord.y - supResB[atom][1]);
		(*ResBgrad)[atom].z += 2.0*K*(ResB->atoms[atom].coord.z - supResB[atom][2]);
	}
}

void getmsegrad(const double &length, double &lengthgrad, const double &ideal)
{
	double dist = length - ideal;
	lengthgrad = 2.*dist;
}

void getdistgrad(const double &x1, const double &y1, const double &z1, const double &x2, const double &y2, const double &z2, double &dx1, double &dy1, double &dz1, double &dx2, double &dy2, double &dz2)
{
	Point3s vec12;
	vec12.x = x2 - x1; vec12.y = y2 - y1; vec12.z = z2 - z1;
	double length12=1.0/sqrt(vec12.x*vec12.x+vec12.y*vec12.y+vec12.z*vec12.z + 1e-12);
	vec12.x *= length12;
	vec12.y *= length12;
	vec12.z *= length12;

	dx1 = -vec12.x; dy1 = -vec12.y; dz1 = -vec12.z;
	dx2 =  vec12.x; dy2 =  vec12.y; dz2 =  vec12.z;
}

void getcentergrad(double &dx, const int &n)
{
	dx = 1.0/n;
}

void getoutprodgrad(const Point3s &gradprev, const Point3s &vecA, const Point3s &vecB, Point3s &gradA, Point3s &gradB)
{
	gradA.x = vecB.y * gradprev.z - vecB.z * gradprev.y;
    gradA.y = vecB.z * gradprev.x - vecB.x * gradprev.z;
	gradA.z = vecB.x * gradprev.y - vecB.y * gradprev.x;

	gradB.x = gradprev.y * vecA.z - gradprev.z * vecA.y;
    gradB.y = gradprev.z * vecA.x - gradprev.x * vecA.z;
	gradB.z = gradprev.x * vecA.y - gradprev.y * vecA.x;
}


StackGeoobj::StackGeoobj(Point3s* _ResA, Point3s* _ResB, vector<Point3s>* _ResAgrad, vector<Point3s>* _ResBgrad, char &_baseA, char &_baseB, const double &_ideallength, const double &_idealangle, double _K1, double _K2)
{
	ResA = _ResA;
	ResB = _ResB;
	ResAgrad = _ResAgrad;
	ResBgrad = _ResBgrad;
	baseA = _baseA;
	baseB = _baseB;
	K1 = _K1;
	K2 = _K2;
	ideallength = _ideallength;
	idealangle = _idealangle*RADDEG;

	atomnumA=targetatoms[baseA].size();
	atomnumB=targetatoms[baseB].size();
}

StackGeoobjMC::StackGeoobjMC(vector<Point3s*>& _ResA, vector<Point3s*>& _ResB, vector<Point3s*>& _ResAgrad, vector<Point3s*>& _ResBgrad, char &_baseA, char &_baseB, const double &_ideallength, const double &_idealangle, double _K1, double _K2)
{
	ResA = _ResA;
	ResB = _ResB;
	ResAgrad = _ResAgrad;
	ResBgrad = _ResBgrad;
	baseA = _baseA;
	baseB = _baseB;
	K1 = _K1;
	K2 = _K2;
	ideallength = _ideallength;
	//idealangle = _idealangle*RADDEG;
	idealangle = cos(_idealangle*RADDEG);

	atomnumA=targetatoms[baseA].size();
	atomnumB=targetatoms[baseB].size();

	gradA.resize(ResAgrad.size());
	gradB.resize(ResBgrad.size());
}

double StackGeoobjMC::calcenergy() const
{
	Point3s center_dA, center_dB;
	int atom;
	double dx, dy, dz, leng;

	const auto& taA = targetatoms.at(baseA);
	const auto& taB = targetatoms.at(baseB);

	center_dA.x=0.0; center_dA.y=0.0; center_dA.z=0.0;
	//for(int i=0;i<targetatoms[baseA].size();i++){
	for(int i=0;i<(int)taA.size();i++)
	{
		//atom=targetatoms[baseA][i];
		atom=taA[i];
		center_dA.x+=ResA[atom]->x;
		center_dA.y+=ResA[atom]->y;
		center_dA.z+=ResA[atom]->z;
	}
	center_dA.x/=atomnumA; center_dA.y/=atomnumA; center_dA.z/=atomnumA;

	center_dB.x=0.0; center_dB.y=0.0; center_dB.z=0.0;
	//for(int i=0;i<targetatoms[baseB].size();i++){
	for(int i=0;i<(int)taB.size();i++)
	{
		//atom=targetatoms[baseB][i];
		atom=taB[i];
		center_dB.x+=ResB[atom]->x;
		center_dB.y+=ResB[atom]->y;
		center_dB.z+=ResB[atom]->z;
	}
	center_dB.x/=atomnumB; center_dB.y/=atomnumB; center_dB.z/=atomnumB;

	dx = center_dA.x - center_dB.x;
    dy = center_dA.y - center_dB.y;
    dz = center_dA.z - center_dB.z;
	leng = sqrt(dx*dx + dy*dy+ dz*dz + 1e-12);

	//double energy1 = 0;
	//if(leng > 6.0){

	//cout<<leng << "-" << ideallength<<endl;

	double diff_leng;
	if (leng < ideallength) diff_leng = 0.0;
	else diff_leng = leng - ideallength;
	//double energy1 = K1*diff*diff;
	//}

	Point3s vecAC42, vecAC86, z1;
	Point3s vecBC42, vecBC86, z2;
	Point3s crossvecA, crossvecB;
	vecAC42.x = ResA[1]->x - center_dA.x; vecAC42.y = ResA[1]->y - center_dA.y; vecAC42.z = ResA[1]->z - center_dA.z;
	vecAC86.x = ResA[2]->x - center_dA.x; vecAC86.y = ResA[2]->y - center_dA.y; vecAC86.z = ResA[2]->z - center_dA.z;
	vecBC42.x = ResB[1]->x - center_dB.x; vecBC42.y = ResB[1]->y - center_dB.y; vecBC42.z = ResB[1]->z - center_dB.z;
	vecBC86.x = ResB[2]->x - center_dB.x; vecBC86.y = ResB[2]->y - center_dB.y; vecBC86.z = ResB[2]->z - center_dB.z;

	crossvecA.x = vecAC42.y * vecAC86.z - vecAC42.z * vecAC86.y;
    crossvecA.y = vecAC42.z * vecAC86.x - vecAC42.x * vecAC86.z;
	crossvecA.z = vecAC42.x * vecAC86.y - vecAC42.y * vecAC86.x;

	crossvecB.x = vecBC42.y * vecBC86.z - vecBC42.z * vecBC86.y;
    crossvecB.y = vecBC42.z * vecBC86.x - vecBC42.x * vecBC86.z;
	crossvecB.z = vecBC42.x * vecBC86.y - vecBC42.y * vecBC86.x;

	double  b;
	b = 1/sqrt((double)(crossvecA.x*crossvecA.x + crossvecA.y*crossvecA.y + crossvecA.z*crossvecA.z) + 1e-12);
  	crossvecA.x *= b; crossvecA.y *= b; crossvecA.z *= b;

	b = 1/sqrt((double)(crossvecB.x*crossvecB.x + crossvecB.y*crossvecB.y + crossvecB.z*crossvecB.z) + 1e-12);
  	crossvecB.x *= b; crossvecB.y *= b; crossvecB.z *= b;

	//double length21=sqrt(crossvecA.x*crossvecA.x+crossvecA.y*crossvecA.y+crossvecA.z*crossvecA.z+1e-8)+1e-8;
	//double length23=sqrt(crossvecB.x*crossvecB.x+crossvecB.y*crossvecB.y+crossvecB.z*crossvecB.z+1e-8)+1e-8;

	//cout<<length21 << "-" << length23<<endl;

	//crossvecA.x /= length21; crossvecA.y /= length21; crossvecA.z /= length21;
	//crossvecB.x /= length23; crossvecB.y /= length23; crossvecB.z /= length23;

	double  costheta = crossvecA.x*crossvecB.x+crossvecA.y*crossvecB.y+crossvecA.z*crossvecB.z;
	if(     costheta > 1. ) costheta = 1.;
  	else if(costheta < -1.) costheta = -1.;

	if(costheta < 0) costheta *= -1;

	//double angle_act = acos(costheta);
	//double angle = (angle_act <= M_PI_2) ? angle_act : (M_PI - angle_act);

	//cout<<idealangle << "-" << angle<<endl;

	//double energy2;
	//if(angle > 30){
	//if (angle > idealangle) diff = 0.0;
	//else diff = angle - idealangle;

	//cout<<costheta << "-" << idealangle << "-" << acos(costheta)*(180/M_PI) << "-" << acos(idealangle)*(180/M_PI)<<endl;

	double diff_ang;
	if (costheta > idealangle) diff_ang = 0.0;
	else diff_ang = costheta - idealangle;
	//double energy2 = K2*diff*diff;
	//}else{
		//energy2 = 0;
	//}

	double energy = -K1*exp(-(diff_ang*diff_ang + diff_leng*diff_leng));

	return energy;

}

void StackGeoobjMC::calcgradient()
{
	for(int i=0;i<ResAgrad.size();i++){
		gradA[i].x = 0; gradA[i].y = 0; gradA[i].z = 0;
	}
	for(int i=0;i<ResBgrad.size();i++){
		gradB[i].x = 0; gradB[i].y = 0; gradB[i].z = 0;
	}

	Point3s center_dA, center_dB;
	int atom;
	double dx, dy, dz, length, dmseleng, dmseangle, dangle;
	double ddistAx, ddistAy, ddistAz, ddistBx, ddistBy, ddistBz;
	double dcenterAx, dcenterAy, dcenterAz, dcenterBx, dcenterBy, dcenterBz;

	const auto& taA = targetatoms.at(baseA);
	const auto& taB = targetatoms.at(baseB);

	center_dA.x=0.0; center_dA.y=0.0; center_dA.z=0.0;
	//for(int i=0;i<targetatoms[baseA].size();i++){
	for(int i=0;i<(int)taA.size();i++){
		//atom=targetatoms[baseA][i];
		atom=taA[i];
		center_dA.x+=ResA[atom]->x;
		center_dA.y+=ResA[atom]->y;
		center_dA.z+=ResA[atom]->z;
	}
	center_dA.x/=atomnumA; center_dA.y/=atomnumA; center_dA.z/=atomnumA;

	center_dB.x=0.0; center_dB.y=0.0; center_dB.z=0.0;
	//for(int i=0;i<targetatoms[baseB].size();i++){
	for(int i=0;i<(int)taB.size();i++){
		//atom=targetatoms[baseB][i];
		atom=taB[i];
		center_dB.x+=ResB[atom]->x;
		center_dB.y+=ResB[atom]->y;
		center_dB.z+=ResB[atom]->z;
	}
	center_dB.x/=atomnumB; center_dB.y/=atomnumB; center_dB.z/=atomnumB;

	dx = center_dA.x - center_dB.x;
    dy = center_dA.y - center_dB.y;
    dz = center_dA.z - center_dB.z;
	length = sqrt(dx*dx + dy*dy+ dz*dz + 1e-12);

	getcentergrad(dcenterAx, atomnumA);
	getcentergrad(dcenterAy, atomnumA);
	getcentergrad(dcenterAz, atomnumA);
	getcentergrad(dcenterBx, atomnumB);
	getcentergrad(dcenterBy, atomnumB);
	getcentergrad(dcenterBz, atomnumB);

	double diff_leng;
	if (length < ideallength) diff_leng = 0.0;
	else diff_leng = length - ideallength;

	if(length >= ideallength){
		getmsegrad(length, dmseleng, ideallength);
		getdistgrad(center_dA.x, center_dA.y, center_dA.z, center_dB.x, center_dB.y, center_dB.z, ddistAx, ddistAy, ddistAz, ddistBx, ddistBy, ddistBz);

		//for(int i=0;i<targetatoms[baseA].size();i++){
		for(int i=0;i<(int)taA.size();i++){
			//atom=targetatoms[baseA][i];
			atom=taA[i];
			gradA[atom].x -= dmseleng*ddistAx*dcenterAx;
			gradA[atom].y -= dmseleng*ddistAy*dcenterAy;
			gradA[atom].z -= dmseleng*ddistAz*dcenterAz;
		}

		//for(int i=0;i<targetatoms[baseB].size();i++){
		for(int i=0;i<(int)taB.size();i++){
			//atom=targetatoms[baseB][i];
			atom=taB[i];
			gradB[atom].x -= dmseleng*ddistBx*dcenterBx;
			gradB[atom].y -= dmseleng*ddistBy*dcenterBy;
			gradB[atom].z -= dmseleng*ddistBz*dcenterBz;
		}
	}

	Point3s vecAC42, vecAC86, z1;
	Point3s vecBC42, vecBC86, z2;
	Point3s crossvecA, crossvecB, crossvecA_norm, crossvecB_norm;
	vecAC42.x = ResA[1]->x - center_dA.x; vecAC42.y = ResA[1]->y - center_dA.y; vecAC42.z = ResA[1]->z - center_dA.z;
	vecAC86.x = ResA[2]->x - center_dA.x; vecAC86.y = ResA[2]->y - center_dA.y; vecAC86.z = ResA[2]->z - center_dA.z;
	vecBC42.x = ResB[1]->x - center_dB.x; vecBC42.y = ResB[1]->y - center_dB.y; vecBC42.z = ResB[1]->z - center_dB.z;
	vecBC86.x = ResB[2]->x - center_dB.x; vecBC86.y = ResB[2]->y - center_dB.y; vecBC86.z = ResB[2]->z - center_dB.z;

	crossvecA.x = vecAC42.y * vecAC86.z - vecAC42.z * vecAC86.y;
    crossvecA.y = vecAC42.z * vecAC86.x - vecAC42.x * vecAC86.z;
	crossvecA.z = vecAC42.x * vecAC86.y - vecAC42.y * vecAC86.x;

	crossvecB.x = vecBC42.y * vecBC86.z - vecBC42.z * vecBC86.y;
    crossvecB.y = vecBC42.z * vecBC86.x - vecBC42.x * vecBC86.z;
	crossvecB.z = vecBC42.x * vecBC86.y - vecBC42.y * vecBC86.x;

	double  bA, bB;
	bA = 1/sqrt((double)(crossvecA.x*crossvecA.x + crossvecA.y*crossvecA.y + crossvecA.z*crossvecA.z) + 1e-12);
  	crossvecA_norm.x = crossvecA.x * bA; crossvecA_norm.y = crossvecA.y * bA; crossvecA_norm.z = crossvecA.z * bA;

	bB = 1/sqrt((double)(crossvecB.x*crossvecB.x + crossvecB.y*crossvecB.y + crossvecB.z*crossvecB.z) + 1e-12);
	crossvecB_norm.x = crossvecB.x * bB; crossvecB_norm.y = crossvecB.y * bB; crossvecB_norm.z = crossvecB.z * bB;

	//double length21=sqrt(crossvecA.x*crossvecA.x+crossvecA.y*crossvecA.y+crossvecA.z*crossvecA.z+1e-8)+1e-8;
	//double length23=sqrt(crossvecB.x*crossvecB.x+crossvecB.y*crossvecB.y+crossvecB.z*crossvecB.z+1e-8)+1e-8;

	//crossvecA.x /= length21; crossvecA.y /= length21; crossvecA.z /= length21;
	//crossvecB.x /= length23; crossvecB.y /= length23; crossvecB.z /= length23;

	double  costheta = crossvecA_norm.x*crossvecB_norm.x+crossvecA_norm.y*crossvecB_norm.y+crossvecA_norm.z*crossvecB_norm.z;
	if(     costheta > 1. ) costheta = 1.;
  	else if(costheta < -1.) costheta = -1.;

	//double angle_act = acos(costheta);
	//double angle = (angle_act <= M_PI_2) ? angle_act : (M_PI - angle_act);

	double costheta_temp = costheta;
	if(costheta < 0) costheta_temp *= -1;

	double diff_ang;
	if (costheta_temp > idealangle) diff_ang = 0.0;
	else diff_ang = costheta_temp - idealangle;

	//cout<<costheta << "-" << costheta_temp<<endl;

	Point3s crossgradA, crossgradB;
	Point3s gradAC42, gradAC86;
	Point3s gradBC42, gradBC86;
	Point3s g_centerA, g_centerB;
	
	//if(angle >= idealangle){
	if(costheta_temp <= idealangle){
		getmsegrad(costheta_temp, dmseangle, idealangle);
		//dangle = (costheta >= 0) ? 1.0 : -1.0;
		double dabs = (std::fabs(costheta) < 1e-8) ? 0.0 : (costheta >= 0.0 ? 1.0 : -1.0);

		crossgradA.x = (crossvecB_norm.x - costheta*crossvecA_norm.x)*bA;
		crossgradA.y = (crossvecB_norm.y - costheta*crossvecA_norm.y)*bA;
		crossgradA.z = (crossvecB_norm.z - costheta*crossvecA_norm.z)*bA;

		crossgradB.x = (crossvecA_norm.x - costheta*crossvecB_norm.x)*bB;
		crossgradB.y = (crossvecA_norm.y - costheta*crossvecB_norm.y)*bB;
		crossgradB.z = (crossvecA_norm.z - costheta*crossvecB_norm.z)*bB;
		
		getoutprodgrad(crossgradA, vecAC42, vecAC86, gradAC42, gradAC86);
		getoutprodgrad(crossgradB, vecBC42, vecBC86, gradBC42, gradBC86);

		double coeff=dmseangle*dabs;
		//double coeff=K2*dmseangle;

		gradA[1].x -= coeff*gradAC42.x;  gradA[1].y -= coeff*gradAC42.y;  gradA[1].z -= coeff*gradAC42.z;
		gradA[2].x -= coeff*gradAC86.x;  gradA[2].y -= coeff*gradAC86.y;  gradA[2].z -= coeff*gradAC86.z;
		gradB[1].x -= coeff*gradBC42.x;  gradB[1].y -= coeff*gradBC42.y;  gradB[1].z -= coeff*gradBC42.z;
		gradB[2].x -= coeff*gradBC86.x;  gradB[2].y -= coeff*gradBC86.y;  gradB[2].z -= coeff*gradBC86.z;

		g_centerA.x = 0; g_centerA.y = 0; g_centerA.z = 0;
		g_centerA.x -= gradAC42.x; g_centerA.y -= gradAC42.y; g_centerA.z -= gradAC42.z;
		g_centerA.x -= gradAC86.x; g_centerA.y -= gradAC86.y; g_centerA.z -= gradAC86.z;

		g_centerB.x = 0; g_centerB.y = 0; g_centerB.z = 0;
		g_centerB.x -= gradBC42.x; g_centerB.y -= gradBC42.y; g_centerB.z -= gradBC42.z;
		g_centerB.x -= gradBC86.x; g_centerB.y -= gradBC86.y; g_centerB.z -= gradBC86.z;

		//for(int i=0;i<targetatoms[baseA].size();i++){
		for(int i=0;i<(int)taA.size();i++){
			//atom=targetatoms[baseA][i];
			atom=taA[i];
			gradA[atom].x += coeff*dcenterAx*g_centerA.x;
			gradA[atom].y += coeff*dcenterAy*g_centerA.y;
			gradA[atom].z += coeff*dcenterAz*g_centerA.z;
		}

		//for(int i=0;i<targetatoms[baseB].size();i++){
		for(int i=0;i<(int)taB.size();i++){
			//atom=targetatoms[baseB][i];
			atom=taB[i];
			gradB[atom].x += coeff*dcenterBx*g_centerB.x;
			gradB[atom].y += coeff*dcenterBy*g_centerB.y;
			gradB[atom].z += coeff*dcenterBz*g_centerB.z;
		}
	}

	double denergydepth = -K1*exp(-(diff_ang*diff_ang + diff_leng*diff_leng));

	for(int i=0;i<(int)taA.size();i++){
		atom=taA[i];
		ResAgrad[atom]->x += denergydepth*gradA[atom].x;
		ResAgrad[atom]->y += denergydepth*gradA[atom].y;
		ResAgrad[atom]->z += denergydepth*gradA[atom].z;
	}

	for(int i=0;i<(int)taB.size();i++){
		atom=taB[i];
		ResBgrad[atom]->x += denergydepth*gradB[atom].x;
		ResBgrad[atom]->y += denergydepth*gradB[atom].y;
		ResBgrad[atom]->z += denergydepth*gradB[atom].z;
	}

}


SingleNtTemplateobj::SingleNtTemplateobj(FragmentResidue* _Res, vector<Point3s>* _Resgrad, const vector<ModelUnit>* _ideal_rna, const vector<int>* _indices, const int &_templatenum, double _K)
{
	Res = _Res;
	Resgrad = _Resgrad;
	ideal_rna = _ideal_rna;
	indices =  _indices;
	templatenum = _templatenum;
	K = _K;

	supRes.resize(Res->atoms.size(), vector<double>(3,0));

	base = Res->restype;

	xyz_list1.resize(targetatoms[base].size(), vector<double>(3,0));
    xyz_list2.resize(targetatoms[base].size(), vector<double>(3,0));
	xyz_list3.resize(targetatoms[base].size(), vector<double>(3,0));
}

double SingleNtTemplateobj::selecttemplate()
{
	int i,j;
	vector<vector<double>> RotMatix, RotMatixmin;
    vector<double> TranVect, TranVectmin;
	vector<double> xyztmp(3,0);
	vector<double> xyz(3,0);
	Point3s vec;
	double rmsd;
	double minrmsd=1.0/0.0;
	int count, idx, ord, minfragidx, minfragord;
	int atom;

	for(i=0; i<templatenum; i++){
		idx = (*indices)[i];

		count = 0;
    	for(j=0; j<targetatoms[base].size();j++){
			atom=targetatoms[base][j];

			xyz_list1[count][0]=Res->atoms[atom].coord.x;
			xyz_list1[count][1]=Res->atoms[atom].coord.y;
			xyz_list1[count][2]=Res->atoms[atom].coord.z;

			xyz_list2[count][0]=(*ideal_rna)[idx].chains[0].residues[0].atoms[atom].xyz[0];
			xyz_list2[count][1]=(*ideal_rna)[idx].chains[0].residues[0].atoms[atom].xyz[1];
			xyz_list2[count][2]=(*ideal_rna)[idx].chains[0].residues[0].atoms[atom].xyz[2];

			count+=1;
    	}

		bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

		rmsd=0;
		count=0;
		for(j=0; j<targetatoms[base].size();j++)
		{
			atom=targetatoms[base][j];

			bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);

			vec.x=xyz_list1[count][0]-xyztmp[0];
			vec.y=xyz_list1[count][1]-xyztmp[1];
			vec.z=xyz_list1[count][2]-xyztmp[2];

			rmsd+=bf.norm(vec);
			count+=1;
		}

		if(minrmsd > rmsd)
		{
            minrmsd=rmsd;
            minfragidx=idx;
			minfragord=ord;
			RotMatixmin=RotMatix;
			TranVectmin=TranVect;
        }
	}

	for(j=0; j<targetatoms[base].size();j++){
		atom=targetatoms[base][j];
		xyz[0]=(*ideal_rna)[minfragidx].chains[0].residues[0].atoms[atom].xyz[0];
		xyz[1]=(*ideal_rna)[minfragidx].chains[0].residues[0].atoms[atom].xyz[1];
		xyz[2]=(*ideal_rna)[minfragidx].chains[0].residues[0].atoms[atom].xyz[2];

		bf.Transform(xyz,RotMatixmin,TranVectmin,xyztmp);
		supRes[atom][0] = xyztmp[0];
		supRes[atom][1] = xyztmp[1];
		supRes[atom][2] = xyztmp[2];
	}

	return minrmsd;
}

double SingleNtTemplateobj::calcenergystacktemp()
{
	double dist = selecttemplate();
	return K*dist;
}

void SingleNtTemplateobj::calcgradientstacktemp() const
{
	Point3s vec12;
	int atom;

	for(int j=0; j<targetatoms.at(base).size();j++){
		atom=targetatoms.at(base)[j];
		(*Resgrad)[atom].x += 2.0*K*(Res->atoms[atom].coord.x - supRes[atom][0]);
		(*Resgrad)[atom].y += 2.0*K*(Res->atoms[atom].coord.y - supRes[atom][1]);
		(*Resgrad)[atom].z += 2.0*K*(Res->atoms[atom].coord.z - supRes[atom][2]);
	}
}

SingleNtTemplateobjMC::SingleNtTemplateobjMC(Point3s* _Res, vector<Point3s>* _Resgrad, char &_base, const vector<vector<Residueatoms>>* _ideal_rna, const vector<int>* _indices, const int &_templatenum, double _K)
{
	Res = _Res;
	Resgrad = _Resgrad;
	base = _base;
	ideal_rna = _ideal_rna;
	indices =  _indices;
	templatenum = _templatenum;
	K = _K;

	atomnum=targetatomnum[base];

	supRes.resize(atomnum, vector<double>(3,0));

	xyz_list1.resize(suptargetatoms[base].size(), vector<double>(3,0));
    xyz_list2.resize(suptargetatoms[base].size(), vector<double>(3,0));
}


double SingleNtTemplateobjMC::selecttemplate()
{
	int i,j;
	vector<vector<double>> RotMatix, RotMatixmin;
    vector<double> TranVect, TranVectmin;
	vector<double> xyztmp(3,0);
	vector<double> xyz(3,0);
	Point3s vec;
	double rmsd;
	double minrmsd=1.0/0.0;
	int count, idx, ord, minfragidx, minfragord;
	int atom;

	for(i=0; i<templatenum; i++){
		idx = (*indices)[i];

		count = 0;
    	for(j=0; j<suptargetatoms[base].size();j++)
		{
			atom=suptargetatoms[base][j];

			xyz_list1[count][0]=Res[atom].x;
			xyz_list1[count][1]=Res[atom].y;
			xyz_list1[count][2]=Res[atom].z;

			xyz_list2[count][0]=(*ideal_rna)[idx][0].ptv[atom].x;
			xyz_list2[count][1]=(*ideal_rna)[idx][0].ptv[atom].y;
			xyz_list2[count][2]=(*ideal_rna)[idx][0].ptv[atom].z;

			count+=1;
    	}

		bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

		rmsd=0;
		count=0;
		for(j=0; j<suptargetatoms[base].size();j++)
		{
			atom=suptargetatoms[base][j];
			bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);

			vec.x=xyz_list1[count][0]-xyztmp[0];
			vec.y=xyz_list1[count][1]-xyztmp[1];
			vec.z=xyz_list1[count][2]-xyztmp[2];

			rmsd+=bf.norm(vec);
			count+=1;
		}

		if(minrmsd > rmsd)
		{
            minrmsd=rmsd;
            minfragidx=idx;
			minfragord=ord;
			RotMatixmin=RotMatix;
			TranVectmin=TranVect;
        }
	}

	for(j=0; j<reptargetatoms[base].size();j++){
		atom=reptargetatoms[base][j];
		xyz[0]=(*ideal_rna)[minfragidx][0].ptv[atom].x;
		xyz[1]=(*ideal_rna)[minfragidx][0].ptv[atom].y;
		xyz[2]=(*ideal_rna)[minfragidx][0].ptv[atom].z;

		bf.Transform(xyz,RotMatixmin,TranVectmin,xyztmp);

		supRes[atom][0] = xyztmp[0];
		supRes[atom][1] = xyztmp[1];
		supRes[atom][2] = xyztmp[2];
	}

	return minrmsd;
}

double SingleNtTemplateobjMC::calcenergystacktemp()
{
	double dist = selecttemplate();
	return K*dist;
}

void SingleNtTemplateobjMC::calcgradientstacktemp() const
{
	Point3s vec12;
	int atom;

	for(int j=0; j<reptargetatoms.at(base).size();j++){
		atom=reptargetatoms.at(base)[j];
		(*Resgrad)[atom].x += 2.0*K*(Res[atom].x - supRes[atom][0]);
		(*Resgrad)[atom].y += 2.0*K*(Res[atom].y - supRes[atom][1]);
		(*Resgrad)[atom].z += 2.0*K*(Res[atom].z - supRes[atom][2]);
	}
}

StackTemplateobj2MC::StackTemplateobj2MC(Point3s* _ResA, Point3s* _ResB, vector<Point3s>* _ResAgrad, vector<Point3s>* _ResBgrad, char &_base1, char &_base2, const vector<vector<Residueatoms>>* _ideal_rna, const vector<int>* _indices, const vector<int>* _orders, const int &_templatenum, double _K)
{
	ResA = _ResA;
	ResB = _ResB;
	ResAgrad = _ResAgrad;
	ResBgrad = _ResBgrad;
	base1 = _base1;
	base2 = _base2;
	ideal_rna = _ideal_rna;
	indices =  _indices;
	orders = _orders;
	templatenum = _templatenum;
	K = _K;

	atomnumA=targetatomnum[base1];
	atomnumB=targetatomnum[base2];

	supResA.resize(atomnumA, vector<double>(3,0));
	supResB.resize(atomnumB, vector<double>(3,0));

	xyz_listA1.resize(targetatoms[base1].size(), vector<double>(3,0));
    xyz_listA2.resize(targetatoms[base1].size(), vector<double>(3,0));
	xyz_listA3.resize(targetatoms[base2].size(), vector<double>(3,0));

	xyz_listB1.resize(targetatoms[base2].size(), vector<double>(3,0));
    xyz_listB2.resize(targetatoms[base2].size(), vector<double>(3,0));
	xyz_listB3.resize(targetatoms[base1].size(), vector<double>(3,0));
}

double StackTemplateobj2MC::selecttemplate()
{
	int i,j;
	vector<vector<double>> RotMatixA, RotMatixB, RotMatixminA, RotMatixminB;
    vector<double> TranVectA, TranVectB, TranVectminA, TranVectminB;
	vector<double> xyztmp(3,0);
	vector<double> xyz(3,0);
	Point3s vec;
	double rmsd;
	double minrmsd=100000000000;
	int count, idx, ord, minfragidx, minfragord;
	int atom;
	int base1_temp, base2_temp;

	for(i=0; i<templatenum; i++){
		idx = (*indices)[i];
		ord = (*orders)[i];

		if(ord == 0){
			base1_temp = 0;
			base2_temp = 1;
		}else if(ord == 1){
			base1_temp = 1;
			base2_temp = 0;
		}

		count = 0;
    	for(j=0; j<targetatoms[base1].size();j++){
			atom=targetatoms[base1][j];

			xyz_listA1[count][0]=ResA[atom].x;
			xyz_listA1[count][1]=ResA[atom].y;
			xyz_listA1[count][2]=ResA[atom].z;

			xyz_listA2[count][0]=(*ideal_rna)[idx][base1_temp].ptv[atom].x;
			xyz_listA2[count][1]=(*ideal_rna)[idx][base1_temp].ptv[atom].y;
			xyz_listA2[count][2]=(*ideal_rna)[idx][base1_temp].ptv[atom].z;

			count+=1;
    	}

		bf.RotateCoor(xyz_listA1,xyz_listA2,RotMatixA,TranVectA);

		rmsd=0;
		count=0;
		for(j=0; j<targetatoms[base2].size();j++)
		{
			atom=targetatoms[base2][j];
			xyz_listA3[count][0]=(*ideal_rna)[idx][base2_temp].ptv[atom].x;
			xyz_listA3[count][1]=(*ideal_rna)[idx][base2_temp].ptv[atom].y;
			xyz_listA3[count][2]=(*ideal_rna)[idx][base2_temp].ptv[atom].z;

			bf.Transform(xyz_listA3[count],RotMatixA,TranVectA,xyztmp);

			vec.x=ResB[atom].x-xyztmp[0];
			vec.y=ResB[atom].y-xyztmp[1];
			vec.z=ResB[atom].z-xyztmp[2];

			rmsd+=bf.norm(vec);
			count+=1;
		}

		count = 0;
    	for(j=0; j<targetatoms[base2].size();j++)
		{
			atom=targetatoms[base2][j];
			xyz_listB1[count][0]=ResB[atom].x;
			xyz_listB1[count][1]=ResB[atom].y;
			xyz_listB1[count][2]=ResB[atom].z;

			xyz_listB2[count][0]=(*ideal_rna)[idx][base2_temp].ptv[atom].x;
			xyz_listB2[count][1]=(*ideal_rna)[idx][base2_temp].ptv[atom].y;
			xyz_listB2[count][2]=(*ideal_rna)[idx][base2_temp].ptv[atom].z;

			count+=1;
    	}

		bf.RotateCoor(xyz_listB1,xyz_listB2,RotMatixB,TranVectB);

		count=0;
		for(j=0; j<targetatoms[base1].size();j++)
		{
			atom=targetatoms[base1][j];
			xyz_listB3[count][0]=(*ideal_rna)[idx][base1_temp].ptv[atom].x;
			xyz_listB3[count][1]=(*ideal_rna)[idx][base1_temp].ptv[atom].y;
			xyz_listB3[count][2]=(*ideal_rna)[idx][base1_temp].ptv[atom].z;

			bf.Transform(xyz_listB3[count],RotMatixB,TranVectB,xyztmp);

			vec.x=ResA[atom].x-xyztmp[0];
			vec.y=ResA[atom].y-xyztmp[1];
			vec.z=ResA[atom].z-xyztmp[2];

			rmsd+=bf.norm(vec);
			count+=1;
		}

		if(minrmsd > rmsd)
		{
            minrmsd=rmsd;
            minfragidx=idx;
			minfragord=ord;
			RotMatixminA=RotMatixA;
			TranVectminA=TranVectA;
			RotMatixminB=RotMatixB;
			TranVectminB=TranVectB;
        }
	}

	if(minfragord == 0){
		base1_temp = 0;
		base2_temp = 1;
	}else if(minfragord == 1){
		base1_temp = 1;
		base2_temp = 0;
	}

	for(j=0; j<targetatoms[base2].size();j++){
		atom=targetatoms[base2][j];
		xyz[0]=(*ideal_rna)[minfragidx][base2_temp].ptv[atom].x;
		xyz[1]=(*ideal_rna)[minfragidx][base2_temp].ptv[atom].y;
		xyz[2]=(*ideal_rna)[minfragidx][base2_temp].ptv[atom].z;
		bf.Transform(xyz,RotMatixminA,TranVectminA,xyztmp);
		supResB[atom][0] = xyztmp[0];
		supResB[atom][1] = xyztmp[1];
		supResB[atom][2] = xyztmp[2];
	}

	for(j=0; j<targetatoms[base1].size();j++){
		atom=targetatoms[base1][j];
		xyz[0]=(*ideal_rna)[minfragidx][base1_temp].ptv[atom].x;
		xyz[1]=(*ideal_rna)[minfragidx][base1_temp].ptv[atom].y;
		xyz[2]=(*ideal_rna)[minfragidx][base1_temp].ptv[atom].z;

		bf.Transform(xyz,RotMatixminB,TranVectminB,xyztmp);
		supResA[atom][0] = xyztmp[0];
		supResA[atom][1] = xyztmp[1];
		supResA[atom][2] = xyztmp[2];
	}

	return minrmsd;
}

double StackTemplateobj2MC::calcenergystacktemp()
{
	double dist = selecttemplate();
	return K*dist;
}

void StackTemplateobj2MC::calcgradientstacktemp() const
{
	Point3s vec12;
	int atom;

	for(int j=0; j<targetatoms.at(base1).size();j++){
		atom=targetatoms.at(base1)[j];
		(*ResAgrad)[atom].x += 2.0*K*(ResA[atom].x - supResA[atom][0]);
		(*ResAgrad)[atom].y += 2.0*K*(ResA[atom].y - supResA[atom][1]);
		(*ResAgrad)[atom].z += 2.0*K*(ResA[atom].z - supResA[atom][2]);
	}

	for(int j=0; j<targetatoms.at(base2).size();j++){
		atom=targetatoms.at(base2)[j];
		(*ResBgrad)[atom].x += 2.0*K*(ResB[atom].x - supResB[atom][0]);
		(*ResBgrad)[atom].y += 2.0*K*(ResB[atom].y - supResB[atom][1]);
		(*ResBgrad)[atom].z += 2.0*K*(ResB[atom].z - supResB[atom][2]);
	}
}


double Clashobj::calcenergyclash() const
{
	Point3s vec12;

	vec12.x = teAtomB->x - teAtomA->x;
	vec12.y = teAtomB->y - teAtomA->y;
	vec12.z = teAtomB->z - teAtomA->z;
	
	double length12=sqrt(vec12.x*vec12.x+vec12.y*vec12.y+vec12.z*vec12.z+1e-8);

	if(length12 > vwd) return 0;

	double dist = length12-vwd;

	return K*dist*dist;
}

void Clashobj::calcgradientclash() const
{
	Point3s vec12;

	vec12.x = teAtomB->x - teAtomA->x;
	vec12.y = teAtomB->y - teAtomA->y;
	vec12.z = teAtomB->z - teAtomA->z;

	double length12=sqrt(vec12.x*vec12.x+vec12.y*vec12.y+vec12.z*vec12.z+1e-8)+1e-8;

	if(length12 > vwd) return;

	vec12.x /= length12;
	vec12.y /= length12;
	vec12.z /= length12;

	double coeff=2.*K*(length12 - vwd);

	teAtomAgrad->x -= coeff*vec12.x;
	teAtomAgrad->y -= coeff*vec12.y;
	teAtomAgrad->z -= coeff*vec12.z;
	teAtomBgrad->x += coeff*vec12.x;
	teAtomBgrad->y += coeff*vec12.y;
	teAtomBgrad->z += coeff*vec12.z;
}















