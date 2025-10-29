#ifndef Superpose_HPP
#define Superpose_HPP 1
/*********************************************************
Compute the best rotation matrix and translation vector that optimize 
overlap between two sets of coordinates 
Date: 20100118
*********************************************************/
#include <vector>
#include <cmath>
#include <cstdlib>
#include "BasicFunc.hpp"

/* vector functions using c arrays */

#define ROTATE(a,i,j,k,l) { g = a[i][j]; \
                            h = a[k][l]; \
                            a[i][j] = g-s*(h+g*tau); \
                            a[k][l] = h+s*(g-h*tau); }

using namespace std;

int jacobi3(double a[3][3], double d[3], double v[3][3], int* n_rot);
void normalize(double a[3]);
double dot_prod(double a[3], double b[3]);
void cross(double a[3], double b[3], double c[3]);
int diagonalize_symmetric(double matrix[3][3], 
  double eigen_vec[3][3], 
  double eigenval[3]);
void VectorTimesMatrix(double *v1, const vector<vector<double> > &rotation_matrix,  double *v2);
void ChangeCoor(vector<double> &v1, const vector<vector<double> > &rotation_matrix, const vector<double> &translate_vec,vector<double> &v2);
void RotateCoor(const vector<vector<double> > &target_stru,const vector<vector<double> > &refrence_stru, vector<vector<double> > &RotMatix, vector<double> &TranVect);

#endif
