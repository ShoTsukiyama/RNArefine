/*************************************************************************
Geometry Functions for Structural Bioinformatics
Author: Cao Yang
contact: cao@scu.edu.cn
Date:2008.7.11.
*************************************************************************/

#ifndef Geometry_hpp
#define Geometry_hpp

#pragma once

#include <vector>
#include <cmath>
#include <iostream>
#include <cstdlib> 

using namespace std;

typedef vector<vector<double> > matrix;

void ShowMyvector(const vector<double> &cc);
void crossproduct(const vector<double> &c1, const vector<double> &c2, vector<double> &cc);
double innerproduct(const vector<double> &c1, const vector<double> &c2);
double innerproduct_univ(const vector<double> &c1, const vector<double> &c2);
void vectorsum(const vector<double> &c1, const vector<double> &c2, vector<double> &cc);
void vectorsum_univ(const vector<double> &c1, const vector<double> &c2, vector<double> &cc);
void subtract(const vector<double> &c1, const vector<double> &c2, vector<double> &cc);
void subtract_univ(const vector<double> &c1, const vector<double> &c2, vector<double> &cc);
bool norm(const vector<double> &c, vector<double> &cc);
bool norm_univ(const vector<double> &c, vector<double> &cc);
bool norm_warnless(vector<double> &c);
bool norm_univ_warnless(vector<double> &c);
double vectorlength(const vector<double> &c);
double vectorlength_univ(const vector<double> &c);
double vectorlength2(const vector<double> &c);
double vectorlength2_univ(const vector<double> &c);
void multi(double coefficient, vector<double> &c, vector<double> &cc);
void multi_univ(double coefficient, vector<double> &c, vector<double> &cc);
void multi(double coefficient, vector<double> &cc);
void multi_univ(double coefficient, vector<double> &cc);
double deg2rad(double deg);
double rad2deg(double rad);
double Points2Distance2(const vector<double> &c1, const vector<double> &c2);
double Points2Distance(const vector<double> &c1, const vector<double> &c2);
double Points2Angle(const vector<double> &c1, const vector<double> &c2, const vector<double> &c3);
double Points4Angle(const vector<double> &c1, const vector<double> &c2, const vector<double> &c3, const vector<double> &c4);
double Points2Dihedral(const vector<double> &c1, const vector<double> &c2, const vector<double> &c3, const vector<double> &c4);
void SetMatrix(matrix &sm, int m, int n);
void MatrixTimesMatrix(const matrix &a, const matrix &v, matrix &x, int m, int n, int l);
bool TransVectorTimesVector(const vector<double> &trans, const vector<double> &vctor, matrix &mtx);
bool MatrixTimesTransVector(const matrix &mtx, const vector<double> &tvt, vector<double> &vct);
void RealTimesMatrix(double at, const matrix &mx, matrix &mc);
bool MatrixAddMatrix(const matrix &ma, const matrix &mb, matrix &mc);
bool norm2(const vector<double> &c, vector<double> &cc);
bool norm2_warnless(vector<double> &c);
bool RotationMatrixA(const vector<double> &axis, double angle, matrix &romtx);
bool RotationMatrixB(const vector<double> &axis, double angle, matrix &romtx);
bool CoordinateRotation(vector<double> &pointA, const vector<double> &axis, double angle, vector<double> &pointB);
bool CoordinateRotation(const vector<double> &pointA, const vector<double> &axisA, const vector<double> &axisB, double angle, vector<double> &pointB);
bool GroupRotation(const vector<double> &axisA, const vector<double> &axisB, double angle, vector<vector<double> > &pointB);
bool GroupRotation(const vector<double> &axisA, const vector<double> &axisB, double angle, vector<vector<double> > &pointB, const vector<short>& index);
bool GroupTranslation(const vector<double> &trans, vector<vector<double> > &pointB);

#endif


