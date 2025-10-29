// BasicFunc.h: interface for the BasicFunc class.
//
//////////////////////////////////////////////////////////////////////
#pragma once

//#if !defined(AFX_BASICFUNC_H__131A9948_54B2_4D0B_87A8_BBC321EAFBA3__INCLUDED_)
//#define AFX_BASICFUNC_H__131A9948_54B2_4D0B_87A8_BBC321EAFBA3__INCLUDED_

//#if _MSC_VER > 1000
#pragma once
//#endif // _MSC_VER > 1000
#include "CommonParams.hpp"

class BasicFunc  
{
public:
	BasicFunc();  
	
	double norm(Point3d p1);
	double norm(Point3s p1);

	double dotv(Point3d p1,Point3d p2);
	double dotv(Point3s p1,Point3s p2);
	
	double phi(double xi,double yi,double zi,
	           double xj,double yj,double zj,
			   double xk,double yk,double zk,
			   double xl,double yl,double zl);
	double angv(Point3d p1,Point3d p2);
	double angv(Point3s p1,Point3s p2);

	void basenormalize(Point3s &p1);
	Point3s basesetv(double tx,double ty,double tz);
	Point3s basecross(Point3s p1, Point3s p2);
	void getcoordsystem(Point3s res1N, Point3s res1C1, Point3s res1C2, Point3s &x, Point3s &y, Point3s &z);

	Point3d setv(double tx,double ty,double tz);
	
	Point3s setvs(double tx,double ty,double tz);
	Point3s calcsubtract(Point3s p1, Point3s p2);
	Point3s calcadd(Point3s p1, Point3s p2);
	void normalize(Point3s &p1);
	void normalize(double a[3]);
	bool tor2pos22(float xi,float yi,float zi,float xj,float yj,float zj,float xk,
		    float yk,float zk,float tang,float tleng,float tinner,float *xl,float *yl,float *zl);
	//void transposeRotationMatrix(vector<vector<double> > &RotMatix);
	void Transform(vector<double> xyz, vector<vector<double> > rotationMatrix, vector<double> TranVect, vector<double> &xyzres);
	void Transform(vector<double> xyz_tmp, vector<vector<double> > rotationMatrix, vector<double> TranVect,Point3d &tp);
	void transposeRotationMatrix(vector<vector<double> > &RotMatix);
	void RotateCoor(const vector<vector<double> > target_stru,const vector<vector<double> > refrence_stru, 
			vector<vector<double> > &RotMatix, vector<double> &TranVect);
	void VectorTimesMatrix(double *v1, const vector<vector<double> > &rotation_matrix,  double *v2);
	void ChangeCoor(vector<double> &v1, const vector<vector<double> > &rotation_matrix, 
			const vector<double> &translate_vec,Point3d &v2);
	Point3d unit(Point3d p1);
	Point3s units(Point3s p1);
	Point3d prod(Point3d p1,Point3d p2);
	Point3d scal(Point3d p1,double f);
	Point3d mmat(double mat[9],Point3d p1);
	Point3d ranv(double fac);
	Point3d crossang(Point3d p1,Point3d p2);
	Point3s cross(Point3s p1, Point3s p2);
	void cross(double a[3], double b[3], double c[3]);

	int maxnormal(Point3d p1);
	int minnormal(Point3d p1);

	void  q2rot(double *q,double rmax[]);

	void getnormalvec(Point3s center1, Point3s res11, Point3s res12, Point3s &z);
	double angleBetweenVectors(Point3s a, Point3s b);

	virtual ~BasicFunc();

};

//#endif // !defined(AFX_BASICFUNC_H__131A9948_54B2_4D0B_87A8_BBC321EAFBA3__INCLUDED_)
