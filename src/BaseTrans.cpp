
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#include "SimulationBase.hpp"

using namespace std;
using namespace Eigen;

/*
void basenormalize(Point3s &p1)
{
  double  b = 1/sqrt((double)(p1.x*p1.x + p1.y*p1.y + p1.z*p1.z));
  p1.x *= b;
  p1.y *= b;
  p1.z *= b;
}

Point3s basesetv(double tx,double ty,double tz)
{
	Point3s temp;
	temp.x=tx;
	temp.y=ty;
	temp.z=tz;
	return temp;
}

Point3s calcsubtract(Point3s p1, Point3s p2)
{
  Point3s p3;
  p3=basesetv(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
  return p3;
}

Point3s calcadd(Point3s p1, Point3s p2)
{
  Point3s p3;
  p3=basesetv(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z);
  return p3;
}

Point3s basecross(Point3s p1, Point3s p2)
{
  return basesetv(
            p1.y * p2.z - p1.z * p2.y,
            p1.z * p2.x - p1.x * p2.z,
            p1.x * p2.y - p1.y * p2.x
            );
}

void RNASimulation::getcoordsystem(Point3s res1N, Point3s res1C1, Point3s res1C2, Point3s &x, Point3s &y, Point3s &z)
{
  Point3s u,v;

  u=calcsubtract(res1C1, res1N);
  v=calcsubtract(res1C2, res1N);

  basenormalize(u);
  basenormalize(v);  

  y=calcadd(u,v);
  basenormalize(y);

  z=basecross(u,v);
  basenormalize(z);

  x=basecross(y,z);
  basenormalize(x);
}
*/

double RNASimulation::angleBetweenVectors(Point3s a, Point3s b){
    BasicFunc bf;
    double dot = bf.dotv(a, b);
    double normA = bf.norm(a);
    double normB = bf.norm(b);
    double cosTheta = dot / (normA * normB);
    double angle = acos(cosTheta);
    double angleInDegrees = angle * (180.0 / M_PI);
    if(angleInDegrees > 180) angleInDegrees=angleInDegrees-180;
    if(angleInDegrees < 0) angleInDegrees=angleInDegrees+180;
    return angleInDegrees;
}

//void gettransfparams(Point3s res1N91, Point3s res1C42, Point3s res1C86, Point3s res2N91, Point3s res2C42, Point3s res2C86, double &r1, double &phi1, double &theta1, double &r2, double &phi2, double &theta2, double &alpha1, double &alpha2, double &beta1, double &beta2, double &gamma1, double &gamma2)
void RNASimulation::getbaseangles(Point3s res1N91, Point3s res1C42, Point3s res1C86, Point3s res2N91, Point3s res2C42, Point3s res2C86, double &r, double &alpha, double &beta, double &gamma)
{
  Point3s d;
  Point3s x1, y1, z1;
  Point3s x2, y2, z2;

  bf.getcoordsystem(res1N91, res1C42, res1C86, x1, y1, z1);
  bf.getcoordsystem(res2N91, res2C42, res2C86, x2, y2, z2);

  d.x=res2N91.x-res1N91.x ; d.y=res2N91.y-res1N91.y; d.z=res2N91.z-res1N91.z;
  r = sqrt(d.x * d.x + d.y * d.y + d.z * d.z);

  alpha=angleBetweenVectors(x1, x2);
  beta=angleBetweenVectors(y1, y2);
  gamma=angleBetweenVectors(z1, z2);
}

void RNASimulation::getbaselength(Point3s res1N91, Point3s res2N91, double &r)
{
  Point3s d;
  d.x=res2N91.x-res1N91.x ; d.y=res2N91.y-res1N91.y; d.z=res2N91.z-res1N91.z;
  r = sqrt(d.x * d.x + d.y * d.y + d.z * d.z);
}

void RNASimulation::getbasecenterpoint(Point3s *ptv, char aat, Point3s &center){
  double count;

  center.x=0; center.y=0; center.z=0; count=0;
  for(int i=0;i<BASEFULLATOM[aat].size();i++){
    center.x+=ptv[BASEFULLATOM[aat][i]].x;
    center.y+=ptv[BASEFULLATOM[aat][i]].y;
    center.z+=ptv[BASEFULLATOM[aat][i]].z;
    count+=1;
  }

  center.x=center.x/count;
  center.y=center.y/count;
  center.z=center.z/count;
}

void RNASimulation::getcenter2centervec(Point3s center1, Point3s center2, Point3s &vec){
  vec.x = center2.x-center1.x;
  vec.y = center2.y-center1.y;
  vec.z = center2.z-center1.z;
}

void RNASimulation::getnormalvec(Point3s center1, Point3s res11, Point3s res12, Point3s &z){
  Point3s u,v;

  u=bf.calcsubtract(res11, center1);
  v=bf.calcsubtract(res12, center1);

  z=bf.basecross(u,v);
  bf.basenormalize(z);

  //baseangle = angleBetweenVectors(z1, z2);
}





