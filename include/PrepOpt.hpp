#ifndef PREPOPT_HPP
#define PREPOPT_HPP 1

#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <list>
#include <deque>
#include "MissingRNAatom.hpp"
#include "PDBParser.hpp"
#include "Gradfuncobj.hpp"
#include "CommonParams.hpp"
#include "Logger.hpp"

using namespace std;

extern map<string, vector<string>> full_atoms;
extern map<string, vector<string>> bb_atoms;

void ChangeCoor(const vector<double> &v1, const vector<vector<double> > &rotation_matrix, const vector<double> &translate_vec,Point3s &v2);
void prepopt(const int c1, const int r1,
const int c2, const int r2, const ModelUnit &pdb_entry,
map<string, vector<ModelUnit>> &ideal_rna, bool stack_prev, bool stack_next, 
vector<vector<double>> minRotMatix, vector<double> minTranVect, int mintempidx,
FragmentStructure &superimposedstr, FragmentStructure &inputtedstr, vector<vector<bool>> &inputtedstrflag);
void prepoptloop(const int c, const int r, const ModelUnit &pdb_entry,
    map<string, vector<ModelUnit>> &ideal_rna, bool stack_prev, bool stack_next, 
    vector<vector<double>> minRotMatix, vector<double> minTranVect, int mintempidx,
    FragmentStructure &superimposedstr, FragmentStructure &inputtedstr, vector<vector<bool>> &inputtedstrflag);
void prepgrad(const FragmentStructure &stru, vector<vector<Point3s>> &gradptv);
void initgrad(const FragmentStructure &stru, vector<vector<Point3s>> &gradptv);
void setinitobj(FragmentStructure &stru, const FragmentStructure &targetstr, const vector<vector<bool>> &targetflag, vector<vector<Point3s>> &gradient, list<InitPosobj> &binitobj, double K, vector<bool> &moveresidues);
void setinitpairobj(FragmentStructure &stru, const FragmentStructure &targetstr, const vector<vector<bool>> &targetflag, vector<vector<Point3s>> &gradient, list<InitPairobj> &initpairobj, double K, vector<bool> &moveresidues);
void setbasepairobj(FragmentStructure &stru, vector<vector<Point3s>> &gradient, vector<vector<int>> &bindmap, const vector<vector<double>> &cbpmap, const vector<vector<double>> &noncbpmap, list<BasePairobj> &bpobj);

//void selectstackfrag(FragmentStructure &stru, vector<vector<Point3s>> &gradient, int res1, int res2, const vector<ModelUnit> &ideal_rna, const vector<vector<vector<double>>> &distmap, const int &sttype, int &minfragidx, double &minrmsd, int &baseorder, double K, list<Stackingobj> &stackobj);
void selectstackfrag(FragmentStructure &stru, vector<vector<Point3s>> &gradient, int res1, int res2, const vector<ModelUnit> &ideal_rna, vector<int> &indices, vector<int> &orders, int ord, const int &sttype, double K, list<StackTemplateobj> &stackobj);
void setstackingobj(FragmentStructure &stru, vector<vector<Point3s>> &gradient, const map<string, map<string, vector<ModelUnit>>> &ideal_rna, const vector<vector<pair<int, double>>> &stackmap, const map<string, map<string, vector<vector<vector<double>>>>> &stackdistmap, double K, list<Stackingobj> &stackobj);
void setstackinggeoobj(FragmentStructure &stru, vector<vector<Point3s>> &gradient, const map<int, map<string, map<string, vector<vector<double>>>>> &stackparam, const vector<vector<pair<int, double>>> &stackmap, vector<vector<Point3s*>> &AtomPtrslistA, vector<vector<Point3s*>> &AtomgradlistA, vector<vector<Point3s*>> &AtomPtrslistB, vector<vector<Point3s*>> &AtomgradlistB, double K, list<StackGeoobjMC> &stackgeoobj);

void stack_superimpose(FragmentStructure &stru, const ModelUnit &ideal_rna_one, int res1, int res2, int order, int ord, double &rmsd);
void selectstackfrag(FragmentStructure &stru, vector<vector<Point3s>> &gradient, int res1, int res2, const vector<ModelUnit> &ideal_rna, int stidx, vector<vector<int>> &indices, vector<vector<int>> &orders, int ord, const int &sttype, double K, list<StackTemplateobj2> &stackobj);
void setstackingobj(FragmentStructure &stru, vector<vector<Point3s>> &gradient, const map<string, map<string, vector<ModelUnit>>> &ideal_rna, const vector<vector<pair<int, double>>> &stackmap, const map<string, map<string, vector<vector<vector<double>>>>> &stackdistmap, vector<vector<int>> &indices, vector<vector<int>> &orders, double K, list<StackTemplateobj2> &stackobj);

void setbpobj(FragmentStructure &stru, vector<vector<Point3s>> &gradient, const map<string, map<string, vector<ModelUnit>>> &ideal_rna, vector<vector<int>> &bpmap, const map<string, vector<vector<vector<double>>>> &bpdistmap, double K, list<BPdistobj> &bpdistobj);
void selectbpfrag(FragmentStructure &stru, vector<vector<Point3s>> &gradient, int res1, int res2, const vector<ModelUnit> &ideal_rna, const vector<vector<vector<double>>> &distmap, int &minfragidx, double &minrmsd, const double &K, list<BPdistobj> &bpdistobj);

//void setlengthobj(FragmentStructure &stru, vector<vector<Point3s>> &gradient, int trandposst, int trandposed, vector<vector<int>> &bindmap, const map<string, vector<vector<double>>> &blengthparam, const map<string, vector<vector<double>>> &adjblengthparam, const map<string, vector<vector<double>>> &addlengthparam, const map<string, vector<vector<double>>> &adjblengthparam, double K, list<BondLengthobj> &blengthobj, vector<bool> &moveresidues);
void setlengthobj(FragmentStructure &stru, vector<vector<Point3s>> &gradient, int trandposst, int trandposed, vector<vector<int>> &bindmap, const map<string, vector<vector<double>>> &blengthparam, const map<string, vector<vector<double>>> &adjblengthparam, const map<string, vector<vector<double>>> &hydparam, const map<string, vector<vector<int>>> &hydparamidx, const map<string, vector<vector<double>>> &addhydparam, const map<string, vector<vector<int>>> &addhydparamidx, double K, list<BondLengthobj> &blengthobj, vector<bool> &moveresidues);
//void setangleobj(FragmentStructure &stru, vector<vector<Point3s>> &gradient, int trandposst, int trandposed, vector<vector<int>> &bindmap, const map<string, vector<vector<double>>> &ntangleparam, const map<string, vector<vector<int>>> &ntangleparamidx, const map<string, vector<vector<double>>> &adjangleparam, const map<string, vector<vector<int>>> &adjangleparamidx, const map<string, vector<vector<double>>> &adjangleparam, const map<string, vector<vector<int>>> &adjangleparamidx, double K, list<BondAngleobj> &bangleobj, vector<bool> &moveresidues);
void setangleobj(FragmentStructure &stru, vector<vector<Point3s>> &gradient, int trandposst, int trandposed, vector<vector<int>> &bindmap, const map<string, vector<vector<double>>> &ntangleparam, const map<string, vector<vector<int>>> &ntangleparamidx, const map<string, vector<vector<double>>> &adjangleparam, const map<string, vector<vector<int>>> &adjangleparamidx, const map<string, vector<vector<double>>> &hydangleparam, const map<string, vector<vector<int>>> &hydangleparamidx, double K, list<BondAngleobj> &bangleobj, vector<bool> &moveresidues);
//void setdihedralobj(FragmentStructure &stru, vector<vector<Point3s>> &gradient, int trandposst, int trandposed, const map<string, vector<pair<vector<int>, vector<double>>>> &hydtorsionparam, list<BondDihedralobj> &dihedralobj, vector<bool> &moveresidues);
void setdihedralobj(FragmentStructure &stru, vector<vector<Point3s>> &gradient, int trandposst, int trandposed, const map<string, vector<pair<vector<int>, vector<double>>>> &hydtorsionparam, map<string, vector<pair<vector<int>, vector<double>>>> torsionparam, double K, list<BondDihedralobj> &dihedralobj, vector<bool> &moveresidues);

void setimpropersobj(FragmentStructure &stru, vector<vector<Point3s>> &gradient, int trandposst, int trandposed, list<BondAngleobj> &bangleobj, list<BondImpropersobj> &improperobj);
void getgradientinit(list<InitPosobj> binitobj);
double getenergyinit(list<InitPosobj> binitobj);
void getgradientinitpair(list<InitPairobj> binitpairobj);
double getenergyinitpair(list<InitPairobj> binitpairobj);
void getgradienttemplatepair(list<InitPairobj> btemplatepairobj);
double getenergytemplatepair(list<InitPairobj> btemplatepairobj);
void getgradientbasepair(list<BasePairobj> &basepairobj);
double getenergybasepair(list<BasePairobj> &basepairobj);

void getgradientbaseplane(list<BasePairobj> &baseplaneobj);
double getenergybaseplane(list<BasePairobj> &baseplaneobj);

void getgradientbasepairdist(list<BPdistobj> &bpdistobj);
double getenergybasepairdist(list<BPdistobj> &bpdistobj);

void getgradientstacking(list<Stackingobj> &stackobj);
double getenergystacking(list<Stackingobj> &stackobj);

void getgradienthbond(list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj);
double getenergyhbond(list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj);
void getgradientblength(list<BondLengthobj> &blengthobj);
double getenergyblength(const list<BondLengthobj> &blengthobj);
void getgradientbangle(list<BondAngleobj> &bangleobj);
void getgradienthydtorsion(list<BondDihedralobj> &dihedralobj);
double getenergybangle(list<BondAngleobj> &bangleobj);
double getenergyhydtor(list<BondDihedralobj> &hydtorobj);
void getgradientbimpropers(list<BondImpropersobj> &improperobj);
double getenergybimpropers(list<BondImpropersobj> &improperobj);
void getgradientbdihedral(list<BondDihedralobj> &dihedralobj);
double getenergybdihedral(list<BondDihedralobj> &dihedralobj);
void getgradientelectros(list<TElectrobj> &telecteobj);
double getenergyelectros(list<TElectrobj> &telecteobj);
void getgradientlj(list<LennardsJobj> &ljobj);
double getenergylj(list<LennardsJobj> &ljobj);
void getgradient(list<BondLengthobj> &blengthobj, list<BondAngleobj> &bangleobj, list<BondDihedralobj> &dihedralobj, list<BasePairobj> &baseplaneobj, list<BasePairobj> &basepairobj, list<BPdistobj> &bpdistobj, list<Stackingobj> &stackobj, list<StackGeoobjMC> &stackgeoobj, list<InitPosobj> &binitobj, list<InitPairobj> &binitpairobj, list<InitPairobj> &btemplatepairobj, list<TElectrobj> &telecteobj, list<LennardsJobj> &ljobj, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, list<Clashobj> &clashobj, list<SingleNtTemplateobj> &nttempobj);
void getgradientloop(list<BondLengthobj> &blengthobj, list<BondAngleobj> &bangleobj, list<BondDihedralobj> &dihedralobj, list<InitPosobj> &binitobj, list<InitPairobj> &binitpairobj, list<InitPairobj> &btemplatepairobj, list<TElectrobj> &telecteobj, list<LennardsJobj> &ljobj, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, list<Clashobj> &clashobj, list<SingleNtTemplateobj> &nttempobj);
double getcovalentenergystem(list<BondLengthobj> &blengthobj, list<BondAngleobj> &bangleobj, list<BondDihedralobj> &dihedralobj, list<BasePairobj> &baseplaneobj, list<BasePairobj> &basepairobj, list<BPdistobj> &bpdistobj, list<Stackingobj> &stackobj, list<StackGeoobjMC> &stackgeoobj, list<InitPosobj> &binitobj, list<InitPairobj> &binitpairobj, list<InitPairobj> &btemplatepairobj, list<TElectrobj> &telecteobj, list<LennardsJobj> &ljobj, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, list<Clashobj> &clashobj, list<SingleNtTemplateobj> &nttempobj);
double getcovalentenergystemloop(list<BondLengthobj> &blengthobj, list<BondAngleobj> &bangleobj, list<BondDihedralobj> &dihedralobj, list<BasePairobj> &baseplaneobj, list<BasePairobj> &basepairobj, list<Stackingobj> &stackobj, list<StackGeoobjMC> &stackgeoobj, list<InitPosobj> &binitobj, list<InitPairobj> &binitpairobj, list<InitPairobj> &btemplatepairobj, list<TElectrobj> &telecteobj, list<LennardsJobj> &ljobj, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, list<Clashobj> &clashobj, list<SingleNtTemplateobj> &nttempobj);
double calcaveragenorm(vector<vector<Point3s>> &gradient, double clip_val = -1.0);
void updatestructure(FragmentStructure &stru, const vector<vector<Point3s>> &gradient, double step_size, vector<bool> &moveresidues);
double goldenelectionsearch(FragmentStructure &stru, const vector<vector<Point3s>> &gradient, list<BondLengthobj> &blengthobj, list<BondAngleobj> &bangleobj, list<BondDihedralobj> &dihedralobj, list<BasePairobj> &baseplaneobj, list<BasePairobj> &basepairobj, list<BPdistobj> &bpdistobj, list<Stackingobj> &stackobj, list<StackGeoobjMC> &stackgeoobj, list<InitPosobj> &binitobj, list<InitPairobj> &binitpairobj, list<InitPairobj> &btemplatepairobj, list<TElectrobj> &telecteobj, list<LennardsJobj> &ljobj, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, list<Clashobj> &clashobj, list<SingleNtTemplateobj> &nttempobj, vector<bool> &moveresidues);
double goldenelectionsearchloop(FragmentStructure &stru, const vector<vector<Point3s>> &gradient, list<BondLengthobj> &blengthobj, list<BondAngleobj> &bangleobj, list<BondDihedralobj> &dihedralobj, list<BasePairobj> &baseplaneobj, list<BasePairobj> &bpobj, list<Stackingobj> &stackobj, list<StackGeoobjMC> &stackgeoobj, list<InitPosobj> &binitobj, list<InitPairobj> &binitpairobj, list<InitPairobj> &btemplatepairobj, list<TElectrobj> &telecteobj, list<LennardsJobj> &ljobj, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, list<Clashobj> &clashobj, list<SingleNtTemplateobj> &nttempobj, vector<bool> &moveresidues);
double linesearch(FragmentStructure &stru, const vector<vector<Point3s>> &gradient, list<BondLengthobj> &blengthobj, list<BondAngleobj> &bangleobj, list<BondDihedralobj> &dihedralobj, list<BasePairobj> &baseplaneobj, list<BasePairobj> &basepairobj, list<BPdistobj> &bpdistobj, list<Stackingobj> &stackobj, list<StackGeoobjMC> &stackgeoobj, list<InitPosobj> &binitobj, list<InitPairobj> &binitpairobj, list<InitPairobj> &btemplatepairobj, list<TElectrobj> &telecteobj, list<LennardsJobj> &ljobj, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, list<Clashobj> &clashobj, list<SingleNtTemplateobj> &nttempobj, vector<bool> &moveresidues, double alpha0 = 1.0, double rho = 0.5, double beta = 0.1);
void getljpairs(FragmentStructure &superimposedstr, vector<vector<Point3s>> &gradient, const map<string, map<int, vector<double>>> &atomfeatures, int res1, int res2, vector<vector<int>> molConnect, list<LennardsJobj> &ljobj);
void getclashpairs(FragmentStructure &superimposedstr, vector<vector<Point3s>> &gradient, const map<string, map<int, vector<double>>> &atomfeatures, int res1, int res2, vector<vector<int>> molConnect, double K, list<Clashobj> &clashobj);

void getclashpairs(FragmentStructure &superimposedstr, vector<vector<Point3s>> &gradient, const map<string, vector<double>> &clashparam, int res1, int res2, vector<vector<int>> molConnect, double K, list<Clashobj> &clashobj);


void getelectrospairs(FragmentStructure &superimposedstr, vector<vector<Point3s>> &gradient, const map<string, map<int, vector<double>>> &atomfeatures, int res1, int res2, vector<vector<int>> molConnect, list<TElectrobj> &telecteobj);
double calcenergyhbond(const Point3s &teAtomHyd, const Point3s &teAtomDon, const Point3s &teAtomAcc);
void gethbondpairs(FragmentStructure &superimposedstr, vector<vector<Point3s>> &gradient, int res1, int res2, vector<vector<int>> molConnect, vector<vector<bool>> intraatomConnect1, vector<vector<bool>> intraatomConnect2, list<Hbondobj> &hbondobj);

/*
void getmolConnect(FragmentStructure &superimposedstr, vector<vector<Point3s>> &gradient, const vector<vector<int>> &stackmap, const map<string, map<int, vector<double>>> &atomfeatures,
const map<string, vector<vector<int>>> &molConnectone, const map<string, vector<vector<int>>> &molConnecttwo, const map<string, vector<vector<int>>> &molConnectthree, 
const map<string, vector<vector<int>>> &molAdjConnectone, const map<string, vector<vector<int>>> &molAdjConnecttwo, const map<string, vector<vector<int>>> &molAdjConnectthree, 
const map<string, vector<vector<int>>> &molAdjConnectoneidx, const map<string, vector<vector<int>>> &molAdjConnecttwoidx, const map<string, vector<vector<int>>> &molAdjConnectthreeidx, vector<vector<int>> &bindmap, vector<bool> &moveresidues,
list<TElectrobj> &telecteobj, list<LennardsJobj> &ljobj, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, list<Clashobj> &clashobj, map<string, double> enecoeff);

void getmolConnect(FragmentStructure &superimposedstr, vector<vector<Point3s>> &gradient, const vector<vector<int>> &stackmap, const map<string, map<int, vector<double>>> &atomfeatures, const map<string, vector<double>> &clashparam,
const map<string, vector<vector<int>>> &molConnectone, const map<string, vector<vector<int>>> &molConnecttwo, const map<string, vector<vector<int>>> &molConnectthree, 
const map<string, vector<vector<int>>> &molAdjConnectone, const map<string, vector<vector<int>>> &molAdjConnecttwo, const map<string, vector<vector<int>>> &molAdjConnectthree, 
const map<string, vector<vector<int>>> &molAdjConnectoneidx, const map<string, vector<vector<int>>> &molAdjConnecttwoidx, const map<string, vector<vector<int>>> &molAdjConnectthreeidx, vector<vector<int>> &bindmap, vector<bool> &moveresidues,
list<TElectrobj> &telecteobj, list<LennardsJobj> &ljobj, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, list<Clashobj> &clashobj, map<string, double> enecoeff);
*/

void getmolConnect(FragmentStructure &superimposedstr, vector<vector<Point3s>> &gradient, const vector<vector<pair<int, double>>> &stackmap, const map<string, map<int, vector<double>>> &atomfeatures, const map<string, vector<double>> &clashparam, const map<char, map<char, map<int, map<int, HydroBondParameters>>>> &hbondparam,
    const map<string, vector<vector<int>>> &molConnectone, const map<string, vector<vector<int>>> &molConnecttwo, const map<string, vector<vector<int>>> &molConnectthree, 
    const map<string, vector<vector<int>>> &molAdjConnectone, const map<string, vector<vector<int>>> &molAdjConnecttwo, const map<string, vector<vector<int>>> &molAdjConnectthree, 
    const map<string, vector<vector<int>>> &molAdjConnectoneidx, const map<string, vector<vector<int>>> &molAdjConnecttwoidx, const map<string, vector<vector<int>>> &molAdjConnectthreeidx, const vector<vector<int>> &bindmap, const vector<vector<double>> &cbpmap, const vector<vector<double>> &noncbpmap, vector<bool> &moveresidues,
    list<TElectrobj> &telecteobj, list<LennardsJobj> &ljobj, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, list<Clashobj> &clashobj, map<string, double> enecoeff);

/*
void getmolConnectbp(FragmentStructure &superimposedstr, vector<vector<Point3s>> &gradient, int bp_len, const vector<vector<int>> &stackmap, const map<string, map<int, vector<double>>> &atomfeatures,
const map<string, vector<vector<int>>> &molConnectone, const map<string, vector<vector<int>>> &molConnecttwo, const map<string, vector<vector<int>>> &molConnectthree, 
const map<string, vector<vector<int>>> &molAdjConnectone, const map<string, vector<vector<int>>> &molAdjConnecttwo, const map<string, vector<vector<int>>> &molAdjConnectthree, 
const map<string, vector<vector<int>>> &molAdjConnectoneidx, const map<string, vector<vector<int>>> &molAdjConnecttwoidx, const map<string, vector<vector<int>>> &molAdjConnectthreeidx, vector<vector<int>> &bindmap, vector<bool> &moveresidues,
list<TElectrobj> &telecteobj, list<LennardsJobj> &ljobj, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, list<Clashobj> &clashobj, map<string, double> enecoeff);

void getmolConnectbp(FragmentStructure &superimposedstr, vector<vector<Point3s>> &gradient, int bp_len, const vector<vector<int>> &stackmap, const map<string, map<int, vector<double>>> &atomfeatures, const map<string, vector<double>> &clashparam,
const map<string, vector<vector<int>>> &molConnectone, const map<string, vector<vector<int>>> &molConnecttwo, const map<string, vector<vector<int>>> &molConnectthree, 
const map<string, vector<vector<int>>> &molAdjConnectone, const map<string, vector<vector<int>>> &molAdjConnecttwo, const map<string, vector<vector<int>>> &molAdjConnectthree, 
const map<string, vector<vector<int>>> &molAdjConnectoneidx, const map<string, vector<vector<int>>> &molAdjConnecttwoidx, const map<string, vector<vector<int>>> &molAdjConnectthreeidx, vector<vector<int>> &bindmap, vector<bool> &moveresidues,
list<TElectrobj> &telecteobj, list<LennardsJobj> &ljobj, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, list<Clashobj> &clashobj, map<string, double> enecoeff);
*/

void getmolConnectbp(FragmentStructure &superimposedstr, vector<vector<Point3s>> &gradient, int bp_len, const vector<vector<int>> &stackmap, const map<string, map<int, vector<double>>> &atomfeatures, const map<string, vector<double>> &clashparam, const map<char, map<char, map<int, map<int, HydroBondParameters>>>> &hbondparam,
    const map<string, vector<vector<int>>> &molConnectone, const map<string, vector<vector<int>>> &molConnecttwo, const map<string, vector<vector<int>>> &molConnectthree, 
    const map<string, vector<vector<int>>> &molAdjConnectone, const map<string, vector<vector<int>>> &molAdjConnecttwo, const map<string, vector<vector<int>>> &molAdjConnectthree, 
    const map<string, vector<vector<int>>> &molAdjConnectoneidx, const map<string, vector<vector<int>>> &molAdjConnecttwoidx, const map<string, vector<vector<int>>> &molAdjConnectthreeidx, vector<vector<int>> &bindmap, const vector<vector<double>> &cbpmap, const vector<vector<double>> &noncbpmap, vector<bool> &moveresidues,
    list<TElectrobj> &telecteobj, list<LennardsJobj> &ljobj, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, list<Clashobj> &clashobj, map<string, double> enecoeff);

void lbfgsupdate(const vector<vector<Point3s>> &grad, const deque<vector<vector<Point3s>>> &s_list, const deque<vector<vector<Point3s>>> &y_list, const deque<double> &rho_list, const vector<bool> &movenodeinfo, vector<vector<Point3s>> &z);
void updatelbfgsfeatures(const FragmentStructure &superimposedstr, const FragmentStructure &superimposedstr_prev, const vector<vector<Point3s>> &grad, const vector<vector<Point3s>> &grad_prev, deque<vector<vector<Point3s>>> &s_list, deque<vector<vector<Point3s>>> &y_list, deque<double> &rho_list, const vector<bool> &movenodeinfo, int m);
//void fragmentoptimizationbp(FragmentStructure &superimposedstr, const FragmentStructure &inputtedstr, const vector<vector<bool>> &inputtedstrflag, const vector<FragmentStructure> &templatestrall, const vector<vector<vector<bool>>> &templatestrallflag, const map<string, map<string, vector<ModelUnit>>> &ideal_rna, int bplen, 
//double &startene, double &endene, const ExtParams &extparameters, vector<vector<int>> &bindmap, vector<vector<int>> &stackmap, const map<string, map<string, vector<vector<vector<double>>>>> &stackdistmap, vector<bool> &moveresidues, const int &epochs);
//void fragmentoptimizationloop(FragmentStructure &superimposedstr, const FragmentStructure &inputtedstr, const vector<vector<bool>> &inputtedstrflag, const vector<FragmentStructure> &templatestrall, const vector<vector<vector<bool>>> &templatestrallflag, const map<string, map<string, vector<ModelUnit>>> &ideal_rna, int looplen, 
//double &startene, double &endene, const ExtParams &extparameters, vector<vector<int>> &bindmap, vector<vector<int>> &stackmap, const map<string, map<string, vector<vector<vector<double>>>>> &stackdistmap, vector<bool> &moveresidues, const int &epochs);
//void fragmentoptimizationbpopt(
    //FragmentStructure &superimposedstr, const FragmentStructure &inputtedstr, const vector<vector<bool>> &inputtedstrflag, const vector<FragmentStructure> &templatestrall, const vector<vector<vector<bool>>> &templatestrallflag, 
    //double &startene, double &endene, const ExtParams &extparameters, vector<vector<int>> &bindmap, vector<vector<int>> &stackmap, vector<bool> &moveresidues, const int &epochs);
//void fragmentoptimizationallopt(
    //FragmentStructure &superimposedstr, const FragmentStructure &inputtedstr, const vector<vector<bool>> &inputtedstrflag, const vector<FragmentStructure> &templatestrall, const vector<vector<vector<bool>>> &templatestrallflag, 
    //double &startene, double &endene, const ExtParams &extparameters, vector<vector<int>> &bindmap, vector<vector<int>> &stackmap, vector<bool> &moveresidues, const int &epochs);
void fragmentoptimization(
    FragmentStructure &superimposedstr, const FragmentStructure &inputtedstr, const vector<vector<bool>> &inputtedstrflag, const vector<FragmentStructure> &templatestrall, const vector<vector<vector<bool>>> &templatestrallflag, const map<string, map<string, vector<ModelUnit>>> &ideal_rna, 
    double &startene, double &endene, const ExtParams &extparameters, vector<vector<int>> &bindmap, const vector<vector<pair<int, double>>> &stackmap, const vector<vector<double>> &cbpmap, const vector<vector<double>> &noncbpmap, const map<string, vector<vector<vector<double>>>> &bpdistmap, const map<string, map<string, vector<vector<vector<double>>>>> &stackdistmap, vector<bool> &moveresidues, const double maxtime, const int &epochs, Logger &logger);

    

#endif