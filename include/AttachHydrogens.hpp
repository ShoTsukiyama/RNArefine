#ifndef ATTACHHYDROGENS_HPP
#define ATTACHHYDROGENS_HPP 1

#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <list>
#include "MissingRNAatom.hpp"
#include "PDBParser.hpp"
#include "Gradfuncobj.hpp"

#define EPS 1e-7

//extern map<string, vector<string>> full_atoms;
//extern map<string, vector<string>> bb_atoms;

using namespace std;

double norm(Point3s p1);
int maxnormal(Point3s p1);
Point3s unit(Point3s p1);
Point3s prod(Point3s p1,Point3s p2);
Point3s scal(Point3s p1,double f);
void q2rot(double *q,double rmax[]);
Point3s mmat(double mat[9],Point3s p1);
bool tor2pos22(float xi,float yi,float zi,float xj,float yj,float zj,float xk,float yk,float zk,float tang,float tleng,float tinner,float *xl,float *yl,float *zl);
void AttachHydrogen(FragmentStructure &stru, const map<string, vector<vector<double>>> &hydparam, const map<string, vector<vector<int>>> &hydparamidx);
void AttachHydrogen(FragmentStructure &stru, const vector<vector<vector<double>>> &hydrogencoords, const map<string, vector<vector<double>>> &hydparam, const map<string, vector<vector<int>>> &hydparamidx);
void UpdateHydrogen(FragmentStructure &stru, const map<string, vector<vector<double>>> &hydparam, const map<string, vector<vector<int>>> &hydparamidx);
void RemoveHydrogens(FragmentStructure &stru);

#endif