#ifndef MONTECARLOSIMULATION_HPP
#define MONTECARLOSIMULATION_HPP 1

#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iterator>
#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "MissingRNAatom.hpp"
#include "AttachHydrogens.hpp"
#include "CommonParams.hpp"
#include "StructureOptimization.hpp"
#include "Gradfuncobj.hpp"
#include "BasicFunc.hpp"
#include "Logger.hpp"

typedef struct Binnedparamtype{
    double minval, maxval, span, mean, std, mode, upper, lower, totalenergy;
    int binnum, rangetype, total;
    vector<int> distr;
    vector<double> energy;
    vector<double> prob;
    vector<double> boundary;
	vector<double> center;
}Binnedparamtype;

typedef struct D2Binnedparamtype{
	double minval1, minval2;
	double maxval1, maxval2;
	double span1, span2;
	int binnum1, binnum2;
	vector<int> distr;
	vector<double> energy;
	vector<int> binary;
	vector<double> boundary1, boundary2;
	vector<double> center1, center2;
}D2Binnedparamtype;

typedef struct PairLengInfo
{
	double mean, std, upper, lower;
}PairLengInfo;

extern map<char, vector<vector<double>>> MCLENGTHPARAM;
extern map<char, vector<vector<double>>> MCINNERANGPARAM;
extern map<char, vector<vector<double>>> MCTORSIONPARAM;

extern map<char, vector<vector<double>>> MCBASELENGTHPARAM;
extern map<char, vector<vector<double>>> MCBASEINNERANGPARAM;
extern map<char, vector<vector<double>>> MCBASETORSIONPARAM;

extern map<char, vector<vector<double>>> MCADDINNERANGPARAM;
extern map<char, vector<vector<double>>> MCADDLENGTHPARAM;

extern map<char, vector<vector<double>>> MCOPTORSIONPARAM;
extern map<char, vector<vector<double>>> MCOPINNERANGPARAM;
extern map<char, vector<vector<double>>> MCOPLENGTHPARAM;

extern map<char, vector<vector<double>>> MCHYDTORSIONPARAM;
extern map<char, vector<vector<double>>> MCHYDINNERANGPARAM;
extern map<char, vector<vector<double>>> MCHYDLENGTHPARAM;

extern map<char, vector<vector<double>>> MCHYDADDLENGTHPARAM;
extern map<char, vector<vector<double>>> MCHYDADDANGLEPARAM;

extern map<char, vector<double>> MCVWPARAM;

extern map<char, vector<Binnedparamtype>> MCTORSIONDISTRIBUTION;
extern map<char, vector<Binnedparamtype>> MCBASETORSIONDISTRIBUTION;
extern map<char, vector<Binnedparamtype>> MCOPTORSIONDISTRIBUTION;
extern map<char, vector<Binnedparamtype>> MCHYDTORSIONDISTRIBUTION;

extern map<char, vector<D2Binnedparamtype>> TORSIONPAIRSPARAM;

extern map<char, vector<vector<PairLengInfo>>> ALLPAIRLENGNTPARAM;
extern map<char, map<char, vector<vector<PairLengInfo>>>> PAIRLENGADJPARAM;

extern map<int, map<char, map<char, vector<vector<Residueatoms>>>>> STACKTEMPLATE;
extern map<int, map<char, map<char, vector<vector<double>>>>> STACKPARAM;
extern map<int, map<char, map<char, vector<vector<vector<double>>>>>> STACKDISTMAP;
extern map<char, map<char, vector<vector<vector<double>>>>> BPDISTMAP;
extern map<char, map<char, vector<vector<Residueatoms>>>> BPTEMPLATE;
extern map<char, map<char, vector<vector<double>>>> BPPARAM;
extern map<char, vector<vector<Residueatoms>>> LOOPNTTEMPLATE;
extern map<char, vector<vector<Residueatoms>>> NTTEMPLATE;
extern map<char, map<char, map<int, map<int, HydroBondParameters>>>> BPHBOND;

extern map<string, vector<vector<int>>> molConnectone;
extern map<string, vector<vector<int>>> molConnecttwo; 
extern map<string, vector<vector<int>>> molConnectthree;
extern map<string, vector<vector<int>>> molAdjConnectone;
extern map<string, vector<vector<int>>> molAdjConnecttwo;
extern map<string, vector<vector<int>>> molAdjConnectthree;
extern map<string, vector<vector<int>>> molAdjConnectoneidx;
extern map<string, vector<vector<int>>> molAdjConnecttwoidx;
extern map<string, vector<vector<int>>> molAdjConnectthreeidx;
extern map<string, map<int, vector<double>>> AtomFeatures;

extern map<string, vector<vector<bool>>> intraAtomConnect;
extern map<string, map<string, vector<vector<int>>>> intramolConnect;
extern map<string, vector<vector<int>>> ntmolConnect;

static map<char, map<int, vector<int>>>ParamsList={
	{'A',{{0,{3,4,5,0}},{1,{4,5,0,1}},{2,{5,0,1,2}},{3,{0,1,2,3}},{4,{1,2,3,4}},{5,{2,3,4,5}},{6,{5,3,4,6}},{7,{4,6,3,7}},{8,{3,7,6,8}},{9,{7,8,6,9}},{10,{5,1,0,10}},{11,{10,5,0,11}},{12,{3,7,8,12}},{13,{7,8,12,13}},{14,{13,8,12,14}}}},
	{'C',{{0,{3,4,5,0}},{1,{4,5,0,1}},{2,{5,0,1,2}},{3,{0,1,2,3}},{4,{1,2,3,4}},{5,{2,3,4,5}},{6,{5,3,4,6}},{7,{4,6,3,7}},{8,{3,7,6,8}},{9,{7,8,6,9}},{10,{5,1,0,10}},{11,{10,5,0,11}},{12,{3,7,8,12}},{13,{7,8,12,13}},{14,{13,8,12,14}}}},
	{'G',{{0,{3,4,5,0}},{1,{4,5,0,1}},{2,{5,0,1,2}},{3,{0,1,2,3}},{4,{1,2,3,4}},{5,{2,3,4,5}},{6,{5,3,4,6}},{7,{4,6,3,7}},{8,{3,7,6,8}},{9,{7,8,6,9}},{10,{5,1,0,10}},{11,{10,5,0,11}},{12,{3,7,8,12}},{13,{7,8,12,13}},{14,{13,8,12,14}}}},
	{'U',{{0,{3,4,5,0}},{1,{4,5,0,1}},{2,{5,0,1,2}},{3,{0,1,2,3}},{4,{1,2,3,4}},{5,{2,3,4,5}},{6,{5,3,4,6}},{7,{4,6,3,7}},{8,{3,7,6,8}},{9,{7,8,6,9}},{10,{5,1,0,10}},{11,{10,5,0,11}},{12,{3,7,8,12}},{13,{7,8,12,13}},{14,{13,8,12,14}}}}
};

static map<char, map<int, vector<int>>>AddAngleParamsList={
	{'A',{{0,{5,0,10}},{1,{10,0,11}},{2,{1,0,11}},{3,{4,6,9}},{4,{2,3,7}},{5,{5,4,6}},{6,{4,6,8}},{7,{4,3,7}},{8,{8,7,3}},{9,{6,8,7}},{10,{6,8,12}},{11,{13,12,14}},{12,{19,18,14}},{13,{17,19,18}},{14,{13,19,18}},{15,{15,17,20}},{16,{16,15,17}},{17,{21,16,15}},{18,{21,13,12}}}},
	{'C',{{0,{5,0,10}},{1,{10,0,11}},{2,{1,0,11}},{3,{4,6,9}},{4,{2,3,7}},{5,{5,4,6}},{6,{4,6,8}},{7,{4,3,7}},{8,{8,7,3}},{9,{6,8,7}},{10,{6,8,12}},{11,{14,12,13}},{12,{17,13,18}},{13,{15,17,13}},{14,{16,15,17}},{15,{16,15,19}}}},
	{'G',{{0,{5,0,10}},{1,{10,0,11}},{2,{1,0,11}},{3,{4,6,9}},{4,{2,3,7}},{5,{5,4,6}},{6,{4,6,8}},{7,{4,3,7}},{8,{8,7,3}},{9,{6,8,7}},{10,{6,8,12}},{11,{13,12,14}},{12,{19,18,14}},{13,{17,19,18}},{14,{13,19,18}},{15,{15,17,21}},{16,{16,15,17}},{17,{20,16,15}},{18,{20,16,22}},{19,{20,13,12}}}},
	{'U',{{0,{5,0,10}},{1,{10,0,11}},{2,{1,0,11}},{3,{4,6,9}},{4,{2,3,7}},{5,{5,4,6}},{6,{4,6,8}},{7,{4,3,7}},{8,{8,7,3}},{9,{6,8,7}},{10,{6,8,12}},{11,{13,12,14}},{12,{18,13,17}},{13,{15,17,13}},{14,{16,15,17}},{15,{19,15,16}}}}
};

static map<char, map<int, vector<int>>>HydParamsList={
	{'A',{{0,{3,1,2,22}},{1,{3,1,2,23}},{2,{2,4,3,24}},{3,{5,3,4,25}},{4,{8,4,6,26}},{5,{4,6,9,27}},{6,{6,7,8,28}},{7,{8,12,14,29}},{8,{19,17,20,30}},{9,{19,17,20,31}},{10,{13,21,16,32}}}},
	{'C',{{0,{3,1,2,20}},{1,{3,1,2,21}},{2,{2,4,3,22}},{3,{5,3,4,23}},{4,{8,4,6,24}},{5,{4,6,9,25}},{6,{6,7,8,26}},{7,{13,12,14,27}},{8,{12,14,16,28}},{9,{16,15,19,29}},{10,{16,15,19,30}}}},
	{'G',{{0,{3,1,2,23}},{1,{3,1,2,24}},{2,{2,4,3,25}},{3,{5,3,4,26}},{4,{8,4,6,27}},{5,{4,6,9,28}},{6,{6,7,8,29}},{7,{8,12,14,30}},{8,{20,16,22,31}},{9,{20,16,22,32}},{10,{20,16,15,33}}}},
	{'U',{{0,{3,1,2,20}},{1,{3,1,2,21}},{2,{2,4,3,22}},{3,{5,3,4,23}},{4,{8,4,6,24}},{5,{4,6,9,25}},{6,{6,7,8,26}},{7,{8,12,14,27}},{8,{12,14,16,28}},{9,{12,13,17,29}}}}
};

static map<char, map<int, vector<int>>>BaseParamsList={
	{'A',{{0,{8,12,14,18}},{1,{8,12,13,19}},{2,{12,13,19,17}},{3,{18,19,13,21}},{4,{13,19,17,15}},{5,{19,13,21,16}},{6,{15,19,17,20}}}},
	{'C',{{0,{8,12,13,17}},{1,{8,12,14,16}},{2,{12,14,16,15}},{3,{8,12,13,18}},{4,{13,17,15,19}}}},
	{'G',{{0,{8,12,14,18}},{1,{8,12,13,19}},{2,{12,13,19,17}},{3,{18,19,13,20}},{4,{13,19,17,15}},{5,{19,13,20,16}},{6,{15,19,17,21}},{7,{20,15,16,22}}}},
	{'U',{{0,{8,12,13,17}},{1,{8,12,14,16}},{2,{12,14,16,15}},{3,{8,12,13,18}},{4,{13,17,15,19}}}}
};

static map<char, map<int, vector<int>>>FIRSTOPParamsList={
	{'A',{{10,{2,1,0,10}},{11,{10,1,0,11}}}},
	{'C',{{10,{2,1,0,10}},{11,{10,1,0,11}}}},
	{'G',{{10,{2,1,0,10}},{11,{10,1,0,11}}}},
	{'U',{{10,{2,1,0,10}},{11,{10,1,0,11}}}}
};

static map<char, vector<int>> BaseHydParamIdx={
	{'A',{7,8,9,10}},
	{'C',{7,8,9,10}},
	{'G',{7,8,9,10}},
	{'U',{7,8,9}}
};

static map<char, vector<int>> AccAtomParamNum={
	{'A',{6,12,15,22,33}},
	{'C',{6,12,15,20,31}},
	{'G',{6,12,15,23,34}},
	{'U',{6,12,15,20,30}}
};

static map<char, map<int, vector<int>>>AddLengthParamsList={
	{'A',{{0,{4,0}}, {1,{5,10}}, {2,{5,11}}, {3,{5,1}}, {4,{0,2}}, {5,{1,10}}, {6,{1,11}}, {7,{1,3}}, {8,{2,7}}, {9,{2,4}}, {10,{3,5}}, {11,{3,6}}, {12,{3,8}}, {13,{4,7}}, {14,{4,9}}, {15,{4,8}}, {16,{5,6}}, {17,{6,7}}, {18,{8,9}}, {19,{10,11}}, {20,{7,8}}, {21,{6,12}}, {22,{7,12}}, {23,{8,14}}, {24,{8,13}}, {25,{14,13}}, {26,{15,16}}, {27,{19,18}}, {28,{12,18}}, {29,{12,19}}, {30,{12,21}}, {31,{14,19}}, {32,{18,17}}, {33,{18,13}}, {34,{19,20}}, {35,{19,15}}, {36,{19,21}}, {37,{17,16}}, {38,{17,13}}, {39,{15,20}}, {40,{15,21}}, {41,{16,13}}}},
	{'C',{{0,{4,0}}, {1,{5,10}}, {2,{5,11}}, {3,{5,1}}, {4,{0,2}}, {5,{1,10}}, {6,{1,11}}, {7,{1,3}}, {8,{2,7}}, {9,{2,4}}, {10,{3,5}}, {11,{3,6}}, {12,{3,8}}, {13,{4,7}}, {14,{4,9}}, {15,{4,8}}, {16,{5,6}}, {17,{6,7}}, {18,{8,9}}, {19,{10,11}}, {20,{7,8}}, {21,{6,12}}, {22,{7,12}}, {23,{8,14}}, {24,{8,13}}, {25,{14,13}}, {26,{15,17}}, {27,{12,18}}, {28,{12,17}}, {29,{12,16}}, {30,{13,15}}, {31,{18,17}}, {32,{17,19}}, {33,{17,16}}, {34,{15,14}}, {35,{19,16}}}},
	{'G',{{0,{4,0}}, {1,{5,10}}, {2,{5,11}}, {3,{5,1}}, {4,{0,2}}, {5,{1,10}}, {6,{1,11}}, {7,{1,3}}, {8,{2,7}}, {9,{2,4}}, {10,{3,5}}, {11,{3,6}}, {12,{3,8}}, {13,{4,7}}, {14,{4,9}}, {15,{4,8}}, {16,{5,6}}, {17,{6,7}}, {18,{8,9}}, {19,{10,11}}, {20,{7,8}}, {21,{6,12}}, {22,{7,12}}, {23,{8,14}}, {24,{8,13}}, {25,{14,13}}, {26,{18,19}}, {27,{15,16}}, {28,{12,18}}, {29,{12,19}}, {30,{12,20}}, {31,{14,19}}, {32,{18,17}}, {33,{18,13}}, {34,{19,21}}, {35,{19,15}}, {36,{19,20}}, {37,{17,16}}, {38,{17,13}}, {39,{21,15}}, {40,{15,22}}, {41,{15,20}}, {42,{16,13}}, {43,{22,20}}}},
	{'U',{{0,{4,0}}, {1,{5,10}}, {2,{5,11}}, {3,{5,1}}, {4,{0,2}}, {5,{1,10}}, {6,{1,11}}, {7,{1,3}}, {8,{2,7}}, {9,{2,4}}, {10,{3,5}}, {11,{3,6}}, {12,{3,8}}, {13,{4,7}}, {14,{4,9}}, {15,{4,8}}, {16,{5,6}}, {17,{6,7}}, {18,{8,9}}, {19,{10,11}}, {20,{7,8}}, {21,{6,12}}, {22,{7,12}}, {23,{8,14}}, {24,{8,13}}, {25,{14,13}}, {26,{15,17}}, {27,{12,18}}, {28,{12,17}}, {29,{12,16}}, {30,{13,15}}, {31,{18,17}}, {32,{17,19}}, {33,{17,16}}, {34,{15,14}}, {35,{19,16}}}}
};

static map<char, map<int, vector<int>>>AddHydLengthParamsList={
	{'A',{{0,{22,1}}, {1,{22,3}}, {2,{23,1}}, {3,{23,3}}, {4,{24,2}}, {5,{24,4}}, {6,{25,3}}, {7,{25,5}}, {8,{26,4}}, {9,{26,8}}, {10,{27,6}}, {11,{28,6}}, {12,{28,7}}, {13,{29,18}}, {14,{29,12}}, {15,{30,17}}, {16,{31,17}}, {17,{32,15}}, {18,{32,21}}, {19,{28,12}}}},
	{'C',{{0,{20,1}}, {1,{20,3}}, {2,{21,1}}, {3,{21,3}}, {4,{22,2}}, {5,{22,4}}, {6,{23,3}}, {7,{23,5}}, {8,{24,4}}, {9,{24,8}}, {10,{25,6}}, {11,{26,6}}, {12,{26,7}}, {13,{27,12}}, {14,{27,16}}, {15,{28,15}}, {16,{28,14}}, {17,{29,15}}, {18,{30,15}}, {19,{26,12}}}},
	{'G',{{0,{23,1}}, {1,{23,3}}, {2,{24,1}}, {3,{24,3}}, {4,{25,2}}, {5,{25,4}}, {6,{26,3}}, {7,{26,5}}, {8,{27,4}}, {9,{27,8}}, {10,{28,6}}, {11,{29,6}}, {12,{29,7}}, {13,{30,18}}, {14,{30,12}}, {15,{31,16}}, {16,{32,16}}, {17,{33,16}}, {18,{33,17}}, {19,{29,12}}}},
	{'U',{{0,{20,1}}, {1,{20,3}}, {2,{21,1}}, {3,{21,3}}, {4,{22,2}}, {5,{22,4}}, {6,{23,3}}, {7,{23,5}}, {8,{24,4}}, {9,{24,8}}, {10,{25,6}}, {11,{26,6}}, {12,{26,7}}, {13,{27,12}}, {14,{27,16}}, {15,{28,15}}, {16,{28,14}}, {17,{29,13}}, {18,{29,15}}, {19,{26,12}}}}
};

static map<char, map<int, vector<int>>>AddHydAngleParamsList={
	{'A',{{0,{3,2,22}}, {1,{23,2,22}}, {2,{3,2,23}}, {3,{2,3,24}}, {4,{7,3,24}}, {5,{6,4,25}}, {6,{5,4,25}}, {7,{8,6,26}}, {8,{9,6,26}}, {9,{6,8,28}}, {10,{12,8,28}}, {11,{18,14,29}}, {12,{31,20,30}}, {13,{15,16,32}}}},
	{'C',{{0,{3,2,20}}, {1,{21,2,20}}, {2,{3,2,21}}, {3,{2,3,22}}, {4,{7,3,22}}, {5,{6,4,23}}, {6,{5,4,23}}, {7,{8,6,24}}, {8,{9,6,24}}, {9,{6,8,26}}, {10,{12,8,26}}, {11,{16,14,27}}, {12,{15,16,28}}, {13,{30,19,29}}}},
	{'G',{{0,{3,2,23}}, {1,{24,2,23}}, {2,{3,2,24}}, {3,{2,3,25}}, {4,{7,3,25}}, {5,{6,4,26}}, {6,{5,4,26}}, {7,{8,6,27}}, {8,{9,6,27}}, {9,{6,8,29}}, {10,{12,8,29}}, {11,{18,14,30}}, {12,{32,22,31}}, {13,{17,15,33}}}},
	{'U',{{0,{3,2,20}}, {1,{21,2,20}}, {2,{3,2,21}}, {3,{2,3,22}}, {4,{7,3,22}}, {5,{6,4,23}}, {6,{5,4,23}}, {7,{8,6,24}}, {8,{9,6,24}}, {9,{6,8,26}}, {10,{12,8,26}}, {11,{16,14,27}}, {12,{15,16,28}}, {13,{15,17,29}}}}
};

static map<char, map<int, int>> atomid2mcidx = {
    {'A', {{0,0}, {1,10}, {2,11}, {3,1}, {4,2}, {5,3}, {6,7}, {7,4}, {8,5}, {9,6}, {10,9}, {11,8}, {12,12}, {13,14}, {14,18}, {15,19}, {16,17}, {17,20}, {18,15}, {19,16}, {20,21}, {21,13}, {22,22}, {23,23}, {24,24}, {25,25}, {26,26}, {27,27}, {28,28}, {29,29}, {30,30}, {31,31}, {32,32}}}, 
    {'C', {{0,0}, {1,10}, {2,11}, {3,1}, {4,2}, {5,3}, {6,7}, {7,4}, {8,5}, {9,6}, {10,9}, {11,8}, {12,12}, {13,13}, {14,18}, {15,17}, {16,15}, {17,19}, {18,16}, {19,14}, {20,20}, {21,21}, {22,22}, {23,23}, {24,24}, {25,25}, {26,26}, {27,27}, {28,28}, {29,29}, {30,30}}},
    {'G', {{0,0}, {1,10}, {2,11}, {3,1}, {4,2}, {5,3}, {6,7}, {7,4}, {8,5}, {9,6}, {10,9}, {11,8}, {12,12}, {13,14}, {14,18}, {15,19}, {16,17}, {17,21}, {18,15}, {19,16}, {20,22}, {21,20}, {22,13}, {23,23}, {24,24}, {25,25}, {26,26}, {27,27}, {28,28}, {29,29}, {30,30}, {31,31}, {32,32}, {33,33}}}, 
    {'U', {{0,0}, {1,10}, {2,11}, {3,1}, {4,2}, {5,3}, {6,7}, {7,4}, {8,5}, {9,6}, {10,9}, {11,8}, {12,12}, {13,13}, {14,18}, {15,17}, {16,15}, {17,19}, {18,16}, {19,14}, {20,20}, {21,21}, {22,22}, {23,23}, {24,24}, {25,25}, {26,26}, {27,27}, {28,28}, {29,29}}}
};

static map<char, map<int, int>> mcidx2atomid = {
    {'A', {{0,0}, {10,1}, {11,2}, {1,3}, {2,4}, {3,5}, {7,6}, {4,7}, {5,8}, {6,9}, {9,10}, {8,11}, {12,12}, {14,13}, {18,14}, {19,15}, {17,16}, {20,17}, {15,18}, {16,19}, {21,20}, {13,21}, {22,22}, {23,23}, {24,24}, {25,25}, {26,26}, {27,27}, {28,28}, {29,29}, {30,30}, {31,31}, {32,32}}},
    {'C', {{0,0}, {10,1}, {11,2}, {1,3}, {2,4}, {3,5}, {7,6}, {4,7}, {5,8}, {6,9}, {9,10}, {8,11}, {12,12}, {13,13}, {18,14}, {17,15}, {15,16}, {19,17}, {16,18}, {14,19}, {20,20}, {21,21}, {22,22}, {23,23}, {24,24}, {25,25}, {26,26}, {27,27}, {28,28}, {29,29}, {30,30}}},
    {'G', {{0,0}, {10,1}, {11,2}, {1,3}, {2,4}, {3,5}, {7,6}, {4,7}, {5,8}, {6,9}, {9,10}, {8,11}, {12,12}, {14,13}, {18,14}, {19,15}, {17,16}, {21,17}, {15,18}, {16,19}, {22,20}, {20,21}, {13,22}, {23,23}, {24,24}, {25,25}, {26,26}, {27,27}, {28,28}, {29,29}, {30,30}, {31,31}, {32,32}, {33,33}}},
    {'U', {{0,0}, {10,1}, {11,2}, {1,3}, {2,4}, {3,5}, {7,6}, {4,7}, {5,8}, {6,9}, {9,10}, {8,11}, {12,12}, {13,13}, {18,14}, {17,15}, {15,16}, {19,17}, {16,18}, {14,19}, {20,20}, {21,21}, {22,22}, {23,23}, {24,24}, {25,25}, {26,26}, {27,27}, {28,28}, {29,29}}}
};

static map<char, vector<vector<int>>> combtorsion = {
    {'A', {{5,7},{5,8},{5,9},{5,12},{7,8},{7,9},{7,12},{8,9},{8,12},{9,12}}},
    {'C', {{5,7},{5,8},{5,9},{5,12},{7,8},{7,9},{7,12},{8,9},{8,12},{9,12}}},
    {'G', {{5,7},{5,8},{5,9},{5,12},{7,8},{7,9},{7,12},{8,9},{8,12},{9,12}}},
    {'U', {{5,7},{5,8},{5,9},{5,12},{7,8},{7,9},{7,12},{8,9},{8,12},{9,12}}}
};

static map<char, map<int, vector<vector<int>>>> torcombmov = {
    {'A', {{5, {{7,0,0}, {8,1,0}, {9,2,0}, {12,3,0}}},{7, {{5,0,1}, {8,4,0}, {9,5,0}, {12,6,0}}},{8, {{5,1,1}, {7,4,1}, {9,7,0}, {12,8,0}}},{9, {{5,2,1}, {7,5,1}, {8,7,1}, {12,9,0}}},{12,{{5,3,1}, {7,6,1}, {8,8,1}, {9,9,1}}}}},
    {'C', {{5, {{7,0,0}, {8,1,0}, {9,2,0}, {12,3,0}}},{7, {{5,0,1}, {8,4,0}, {9,5,0}, {12,6,0}}},{8, {{5,1,1}, {7,4,1}, {9,7,0}, {12,8,0}}},{9, {{5,2,1}, {7,5,1}, {8,7,1}, {12,9,0}}},{12,{{5,3,1}, {7,6,1}, {8,8,1}, {9,9,1}}}}},
    {'G', {{5, {{7,0,0}, {8,1,0}, {9,2,0}, {12,3,0}}},{7, {{5,0,1}, {8,4,0}, {9,5,0}, {12,6,0}}},{8, {{5,1,1}, {7,4,1}, {9,7,0}, {12,8,0}}},{9, {{5,2,1}, {7,5,1}, {8,7,1}, {12,9,0}}},{12,{{5,3,1}, {7,6,1}, {8,8,1}, {9,9,1}}}}},
    {'U', {{5, {{7,0,0}, {8,1,0}, {9,2,0}, {12,3,0}}},{7, {{5,0,1}, {8,4,0}, {9,5,0}, {12,6,0}}},{8, {{5,1,1}, {7,4,1}, {9,7,0}, {12,8,0}}},{9, {{5,2,1}, {7,5,1}, {8,7,1}, {12,9,0}}},{12,{{5,3,1}, {7,6,1}, {8,8,1}, {9,9,1}}}}}
};

static vector<int> torcombmovidx={5,7,8,9,12};

/*
static map<string, vector<int>> HBondHydrogenAtoms = {
    {"A",{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,0}},
    {"C",{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1}},
    {"G",{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,1}},
    {"U",{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1}}
};
*/

static map<string, vector<int>> HBondHydrogenAtoms = {
    {"A",{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0}},
    {"C",{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0}},
    {"G",{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0}},
    {"U",{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0}}
};

static map<string, vector<int>> HBondAcceptorAtoms = {
    {"A",{0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0}},
    {"C",{0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}},
    {"G",{0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0}},
    {"U",{0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0}}
};

static map<string, vector<int>> HBondDonerAtoms = {
    {"A",{0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}},
    {"C",{0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}},
    {"G",{0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}},
    {"U",{0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}}
};

//Acc idx, Hyd idx, Don idx, order
static map<string, vector<vector<vector<int>>>> BpHbondInfo = {
    {"AA",{{{21,30,20,0}, {21,31,20,0}}, {{15,30,20,0}, {15,31,20,0}}, {{18,30,20,0}, {18,31,20,0}}}},
    {"AC",{{{17,30,20,1}, {17,31,20,1}}, {{18,29,19,0}, {18,30,19,0}}, {{18,30,20,1}, {18,31,20,1}}, {{21,29,19,0}, {21,30,19,0}}, {{15,29,19,0}, {15,30,19,0}}}},
    {"CA",{{{17,30,20,0}, {17,31,20,0}}, {{18,29,19,1}, {18,30,19,1}}, {{18,30,20,0}, {18,31,20,0}}, {{21,29,19,1}, {21,30,19,1}}, {{15,29,19,1}, {15,30,19,1}}}},
    {"AG",{{{18,30,20,1}, {18,31,20,1}}, {{15,31,22,0}, {15,32,22,0}}, {{18,31,22,0}, {18,32,22,0}}, {{21,31,22,0}, {21,32,22,0}}, {{15,33,15,0}}, {{21,30,20,1}, {21,31,20,1}}, {{18,33,15,0}}, {{21,33,15,0}}}},
    {"GA",{{{18,30,20,0}, {18,31,20,0}}, {{15,31,22,1}, {15,32,22,1}}, {{18,31,22,1}, {18,32,22,1}}, {{21,31,22,1}, {21,32,22,1}}, {{15,33,15,1}}, {{21,30,20,0}, {21,31,20,0}}, {{18,33,15,1}}, {{21,33,15,1}}}},
    {"AU",{{{19,30,20,1}, {19,31,20,1}}, {{15,29,17,0}}}},
    {"UA",{{{19,30,20,0}, {19,31,20,0}}, {{15,29,17,1}}}},
    {"GG",{{{18,31,22,0}, {18,32,22,0}}, {{21,33,15,0}}, {{21,31,22,0},{21,32,22,0}}, {{18,33,15,0}}}},
    {"GC",{{{21,29,19,0}, {21,30,19,0}}, {{17,33,15,1}}, {{18,31,22,1},{18,32,22,1}}}},
    {"CG",{{{21,29,19,1}, {21,30,19,1}}, {{17,33,15,0}}, {{18,31,22,0},{18,32,22,0}}}},
    {"GU",{{{19,31,22,1}, {19,32,22,1}}, {{21,29,17,0}}}},
    {"UG",{{{19,31,22,0}, {19,32,22,0}}, {{21,29,17,1}}}},
    {"CC",{{{18,29,19,0}, {18,30,19,0}}, {{17,29,19,0},{17,30,19,0}}}},
    {"UC",{{{18,29,19,0}, {18,30,19,0}}, {{19,29,19,0},{19,30,19,0}}, {{17,29,17,1}}, {{18,29,17,1}}}},
    {"UC",{{{18,29,19,1}, {18,30,19,1}}, {{19,29,19,1},{19,30,19,1}}, {{17,29,17,0}}, {{18,29,17,0}}}},
    {"UU",{{{18,29,17,0}}, {{19,29,17,0}}}}
};

class ClashMCeneobj
{
        private:
        Point3s* teAtomA;
        Point3s* teAtomB;
        double *atomenergy1, *atomenergy2;
        double vdW;
        int *flag;
        mutable double prevene;

        public:
        ClashMCeneobj(Point3s* _teAtomA, Point3s* _teAtomB, double* _atomenergy1, double* _atomenergy2, int* _flag, double _vdW =-1.0)
            : teAtomA(_teAtomA), teAtomB(_teAtomB), atomenergy1(_atomenergy1), atomenergy2(_atomenergy2), flag(_flag), vdW(_vdW), prevene(0.0) {}
    
    double calcenergyclash() const;

    static ClashMCeneobj create(Point3s* _teAtomA, Point3s* _teAtomB, double* _atomenergy1, double* _atomenergy2, int* _flag, double _vdW = -1.0) {
        return ClashMCeneobj(_teAtomA, _teAtomB, _atomenergy1, _atomenergy2, _flag, _vdW);
    }

    void setprevene(const double &_prevene) const;
};

class ElectrosMCeneobj
{
        private:
        Point3s* teAtomA;
        Point3s* teAtomB;
        double *atomenergy1, *atomenergy2;
        double dQQep;
        int* flag;
        mutable double prevene;

        public:
        ElectrosMCeneobj(Point3s* _teAtomA, Point3s* _teAtomB, double* _atomenergy1, double* _atomenergy2, int* _flag, double _dQQep =-1.0)
            : teAtomA(_teAtomA), teAtomB(_teAtomB), flag(_flag), atomenergy1(_atomenergy1), atomenergy2(_atomenergy2), dQQep(_dQQep), prevene(0.0) {}
    
    double calcenergyelectros() const;

    static ElectrosMCeneobj create(Point3s* _teAtomA, Point3s* _teAtomB, double* _atomenergy1, double* _atomenergy2, int* _flag, double _dQQep = -1.0) {
        return ElectrosMCeneobj(_teAtomA, _teAtomB, _atomenergy1, _atomenergy2, _flag, _dQQep);
    }

    void setprevene(const double &_prevene) const;
};

class LJMCeneobj
{
        private:
        Point3s* teAtomA;
        Point3s* teAtomB;
        double* atomenergy1, *atomenergy2;
        double LJ6, LJ12;
        int* flag;
        mutable double prevene;

        public:
        LJMCeneobj(Point3s* _teAtomA, Point3s* _teAtomB, double* _atomenergy1, double* _atomenergy2, int* _flag, double _LJ6 =-1.0, double _LJ12 =-1.0)
            : teAtomA(_teAtomA), teAtomB(_teAtomB), atomenergy1(_atomenergy1), atomenergy2(_atomenergy2), flag(_flag), LJ6(_LJ6), LJ12(_LJ12), prevene(0.0) {}
    
    double calcenergylj() const;

    static LJMCeneobj create(Point3s* _teAtomA, Point3s* _teAtomB, double* _atomenergy1, double* _atomenergy2, int* _flag, double _LJ6 = -1.0, double _LJ12 = -1.0) {
        return LJMCeneobj(_teAtomA, _teAtomB, _atomenergy1, _atomenergy2, _flag, _LJ6, _LJ12);
    }

    void setprevene(const double &_prevene) const;
};

class HBondMCeneobj
{
        private:
        Point3s* teAtomHyd;
	    Point3s* teAtomDon;
	    Point3s* teAtomAcc;
        double *atomenergy2, *atomenergy3;
        double idealadleng, idealahleng, idealahdang;
        double dK;
        int* flag;
        mutable double prevene;
        double dREq=1.55;
        double dXiR=2.0;
        double dXiTheta=2.0;
        double dThreshold=-1.0;

        public:
        HBondMCeneobj(Point3s* _teAtomHyd, Point3s* _teAtomDon, Point3s* _teAtomAcc, double* _atomenergy2, double* _atomenergy3, const double &_idealadleng, const double &_idealahleng, const double &_idealahdang, double _dK, int* _flag)
            : teAtomHyd(_teAtomHyd), teAtomDon(_teAtomDon), teAtomAcc(_teAtomAcc), atomenergy2(_atomenergy2), atomenergy3(_atomenergy3), idealadleng(_idealadleng), idealahleng(_idealahleng), idealahdang(cos(_idealahdang*RADDEG)), dK(_dK), flag(_flag), prevene(0.0), dREq(1.55), dXiR(2.0), dXiTheta(2.0), dThreshold(-1.0) {}
    
            double calcenergyhydbond() const;

    static HBondMCeneobj create(Point3s* _teAtomHyd, Point3s* _teAtomDon, Point3s* _teAtomAcc, double* _atomenergy2, double* _atomenergy3, const double &_idealadleng, const double &_idealahleng, const double &_idealahdang, double _dK, int* _flag) {
        return HBondMCeneobj(_teAtomHyd, _teAtomDon, _teAtomAcc, _atomenergy2, _atomenergy3, _idealadleng, _idealahleng, _idealahdang, _dK, _flag);
    }

    void setprevene(const double &_prevene) const;
};

class HBond2MCeneobj
{
        private:
        Point3s* teAtomHyd1;
        Point3s* teAtomHyd2;
	    Point3s* teAtomDon;
	    Point3s* teAtomAcc;
        double *atomenergy2, *atomenergy3;
        double idealadleng, idealahleng, idealahdang;
        double dK;
        int* flag;
        mutable double prevene;
        double dREq=1.55;
        double dXiR=2.0;
        double dXiTheta=2.0;
        double dThreshold=-1.0;

        public:
        HBond2MCeneobj(Point3s* _teAtomHyd1, Point3s* _teAtomHyd2, Point3s* _teAtomDon, Point3s* _teAtomAcc, double* _atomenergy2, double* _atomenergy3, const double &_idealadleng, const double &_idealahleng, const double &_idealahdang, double _dK, int* _flag)
            : teAtomHyd1(_teAtomHyd1), teAtomHyd2(_teAtomHyd2), teAtomDon(_teAtomDon), teAtomAcc(_teAtomAcc), atomenergy2(_atomenergy2), atomenergy3(_atomenergy3), idealadleng(_idealadleng), idealahleng(_idealahleng), idealahdang(cos(_idealahdang*RADDEG)), dK(_dK), flag(_flag), prevene(0.0), dREq(1.55), dXiR(2.0), dXiTheta(2.0), dThreshold(-1.0) {}
    
            double calcenergyhydbond() const;

    static HBond2MCeneobj create(Point3s* _teAtomHyd1, Point3s* _teAtomHyd2, Point3s* _teAtomDon, Point3s* _teAtomAcc, double* _atomenergy2, double* _atomenergy3, const double &_idealadleng, const double &_idealahleng, const double &_idealahdang, double _dK, int* _flag) {
        return HBond2MCeneobj(_teAtomHyd1, _teAtomHyd2, _teAtomDon, _teAtomAcc, _atomenergy2, _atomenergy3, _idealadleng, _idealahleng, _idealahdang, _dK, _flag);
    }

    void setprevene(const double &_prevene) const;
};

class InitMCdistobj
{
  private:
  	Point3s *teAtomA;
    Point3s *teAtomB;
    double* atomenergy1, *atomenergy2;
    int* flag;
	const double initdist;
    const double cutoff;
    mutable double prevene;
    const double K;

  public:
  InitMCdistobj(Point3s* _teAtomA, Point3s* _teAtomB, const double _initdist, const double _cutoff, double* _atomenergy1, double* _atomenergy2, const double _K, int* _flag)
  		: teAtomA(_teAtomA), teAtomB(_teAtomB), initdist(_initdist), cutoff(_cutoff), atomenergy1(_atomenergy1), atomenergy2(_atomenergy2), K(_K), flag(_flag), prevene(0.0) {}
    
    double calcenergyinit() const;
    void setprevene(const double &_prevene) const;

	static InitMCdistobj create(Point3s* _teAtomA, Point3s* _teAtomB, const double _initdist, const double _cutoff, double* _atomenergy1, double* _atomenergy2, const double _K, int* _flag) {
        return InitMCdistobj(_teAtomA, _teAtomB, _initdist, _cutoff, _atomenergy1, _atomenergy2, _K, _flag);
    }
};

class InitMCPosobj
{
  private:
  	Point3s* teAtom;
	const double X0, Y0, Z0, K;
    double* atomenergy;
    int* flag;
    const double cutoff;
    mutable double prevene;

  public:
  InitMCPosobj(Point3s* _teAtom, const double _cutoff, double* _atomenergy, double _X0, double _Y0, double _Z0, double _K, int* _flag)
  		: teAtom(_teAtom), cutoff(_cutoff), atomenergy(_atomenergy), X0(_X0), Y0(_Y0), Z0(_Z0), K(_K), flag(_flag), prevene(0.0) {}
    
    double calcenergyinit() const;
    void setprevene(const double &_prevene) const;

	static InitMCPosobj create(Point3s* _teAtom, const double _cutoff, double* _atomenergy, double X0, double Y0, double Z0, double K, int* _flag) {
        return InitMCPosobj(_teAtom, _cutoff, _atomenergy, X0, Y0, Z0, K, _flag);
    }
};

/*
class BasePairMCobj
{
	private:
		vector<Point3s*> AtomPtrs;
		mutable double Xm,Ym,Zm;
		mutable double a,b;
		unsigned int n;
        int *flag;
        vector<double*> atomenergies;
        mutable double prevene;

		double Fit3DPlane() const;
		double TotalSQRDistance(double _a, double _b) const;
    	double TotalSQRDistance() const;
    	void CalcCenterOfMass() const;

	public:
        BasePairMCobj(const vector<Point3s*>& _AtomPtrs, vector<double*> &_atomenergies, int* _flag);

		double calcenergybasepair() const;
        void setprevene(const double &_prevene) const;
};
*/

class BasePairMCobj
{
	private:
		vector<Point3s*> AtomPtrsA;
        vector<Point3s*> AtomPtrsB;
        int *flag;
        vector<double*> atomenergiesA;
        vector<double*> atomenergiesB;
        Point3s *baseA1, *baseA2, *baseA3;
        Point3s *baseB1, *baseB2, *baseB3;
        mutable double prevene;

	public:
        BasePairMCobj(const vector<Point3s*>& _AtomPtrsA, const vector<Point3s*>& _AtomPtrsB, Point3s*_baseA1, Point3s*_baseA2, Point3s*_baseA3, Point3s*_baseB1, Point3s*_baseB2, Point3s*_baseB3, vector<double*> &_atomenergiesA, vector<double*> &_atomenergiesB, int* _flag)
            : AtomPtrsA(_AtomPtrsA), AtomPtrsB(_AtomPtrsB), baseA1(_baseA1), baseA2(_baseA2), baseA3(_baseA3), baseB1(_baseB1), baseB2(_baseB2), baseB3(_baseB3), atomenergiesA(_atomenergiesA), atomenergiesB(_atomenergiesB), flag(_flag), prevene(0.0) {}

		double calcenergybasepair() const;
        void setprevene(const double &_prevene) const;
        double getdistfromplane(double a, double b, double c, double d, Point3s pointP) const;
        void getplanebp(Point3s pointA, Point3s pointB, Point3s pointC, double &a, double &b, double &c, double &d) const;

        static BasePairMCobj create(const vector<Point3s*>& _AtomPtrsA, const vector<Point3s*>& _AtomPtrsB, Point3s*_baseA1, Point3s*_baseA2, Point3s*_baseA3, Point3s*_baseB1, Point3s*_baseB2, Point3s*_baseB3, vector<double*> &_atomenergiesA, vector<double*> &_atomenergiesB, int* _flag) {
            return BasePairMCobj(_AtomPtrsA, _AtomPtrsB, _baseA1, _baseA2, _baseA3, _baseB1, _baseB2, _baseB3, _atomenergiesA, _atomenergiesB, _flag);
        }
};

class BasePairMCobj2
{
	private:
		vector<Point3s*> AtomPtrs;
        int *flag;
        vector<double*> atomenergies;
        mutable double prevene;
        double dK;
        double Xm,Ym,Zm;
		double a,b;
		unsigned int n;

	public:
        BasePairMCobj2(const vector<Point3s*> _AtomPtrs, const vector<double*> _atomenergies, double dK, int* _flag);

		double Fit3DPlane();
        void CalcCenterOfMass();
        double TotalSQRDistance(double _a, double _b) const;
        double TotalSQRDistance() const;
        double calcenergybasepair();
        void setprevene(const double &_prevene) const;
};

class BasePairMCobj3
{
    private:
  	    float* bleng;
        double* ideal;
        double* upper;
        double* lower;
        double *atomenergy1, *atomenergy2;
        int* flag;
        mutable double prevene;
        string remark;

    public:
    BasePairMCobj3(float* _bleng, double* _atomenergy1, double* _atomenergy2, int* _flag, double* _ideal, double* _upper, double* _lower)
  		: bleng(_bleng), atomenergy1(_atomenergy1), atomenergy2(_atomenergy2), flag(_flag), ideal(_ideal), upper(_upper), lower(_lower), prevene(0.0) {}
    
    double calcenergy() const;
    void setprevene(const double &_prevene) const;

	static BasePairMCobj3 create(float* _bleng, double* _atomenergy1, double* _atomenergy2, int* _flag, double* _ideal, double* _upper, double* _lower) {
        return BasePairMCobj3(_bleng, _atomenergy1, _atomenergy2, _flag, _ideal, _upper, _lower);
    }
};

class StackPairMCobj
{
	private:
		vector<Point3s*> AtomPtrsA;
        vector<Point3s*> AtomPtrsB;
        int *flag;
        vector<double*> atomenergiesA;
        vector<double*> atomenergiesB;
        int sttype;
        char baseA, baseB;
        mutable double prevene;
        double K;

	public:
    StackPairMCobj(const vector<Point3s*>& _AtomPtrsA, const vector<Point3s*>& _AtomPtrsB, const char &_baseA, const char &_baseB, const int _sttype, vector<double*> &_atomenergiesA, vector<double*> &_atomenergiesB, int* _flag, double _K)
            : AtomPtrsA(_AtomPtrsA), AtomPtrsB(_AtomPtrsB), baseA(_baseA), baseB(_baseB), sttype(_sttype), atomenergiesA(_atomenergiesA), atomenergiesB(_atomenergiesB), flag(_flag), prevene(0.0), K(_K) {}

		double calcenergystacking() const;
        void setprevene(const double &_prevene) const;
        double getdistfromplane(double a, double b, double c, double d, Point3s pointP) const;
        void getplanebp(Point3s pointA, Point3s pointB, Point3s pointC, double &a, double &b, double &c, double &d) const;

        static StackPairMCobj create(const vector<Point3s*>& _AtomPtrsA, const vector<Point3s*>& _AtomPtrsB, const char &_baseA, const char &_baseB, const int _sttype, vector<double*> &_atomenergiesA, vector<double*> &_atomenergiesB, int* _flag, double _K) {
            return StackPairMCobj(_AtomPtrsA, _AtomPtrsB, _baseA, _baseB, _sttype, _atomenergiesA, _atomenergiesB, _flag, _K);
        }
};

class StackPairGeoMCobj
{
	private:
		vector<Point3s*> AtomPtrsA;
        vector<Point3s*> AtomPtrsB;
        int *flag;
        vector<double*> atomenergiesA;
        vector<double*> atomenergiesB;
        int sttype;
        char baseA, baseB;
        mutable double prevene;
        double K;

	public:
    StackPairGeoMCobj(const vector<Point3s*>& _AtomPtrsA, const vector<Point3s*>& _AtomPtrsB, const char &_baseA, const char &_baseB, const int _sttype, vector<double*> &_atomenergiesA, vector<double*> &_atomenergiesB, int* _flag, double _K)
            : AtomPtrsA(_AtomPtrsA), AtomPtrsB(_AtomPtrsB), baseA(_baseA), baseB(_baseB), sttype(_sttype), atomenergiesA(_atomenergiesA), atomenergiesB(_atomenergiesB), flag(_flag), prevene(0.0), K(_K) {}

		double calcenergystacking() const;
        void setprevene(const double &_prevene) const;
        double getdistfromplane(double a, double b, double c, double d, Point3s pointP) const;
        void getplanebp(Point3s pointA, Point3s pointB, Point3s pointC, double &a, double &b, double &c, double &d) const;

        static StackPairGeoMCobj create(const vector<Point3s*>& _AtomPtrsA, const vector<Point3s*>& _AtomPtrsB, const char &_baseA, const char &_baseB, const int _sttype, vector<double*> &_atomenergiesA, vector<double*> &_atomenergiesB, int* _flag, double _K) {
            return StackPairGeoMCobj(_AtomPtrsA, _AtomPtrsB, _baseA, _baseB, _sttype, _atomenergiesA, _atomenergiesB, _flag, _K);
        }
};

class BondLengthMCeneobj
{
    private:
  	    float* bleng;
        double* ideal;
        double* upper;
        double* lower;
        double *atomenergy1, *atomenergy2;
        int* flag;
        mutable double prevene;
        string remark;

    public:
    BondLengthMCeneobj(float* _bleng, double* _atomenergy1, double* _atomenergy2, int* _flag, double* _ideal, double* _upper, double* _lower, string _remark)
  		: bleng(_bleng), atomenergy1(_atomenergy1), atomenergy2(_atomenergy2), flag(_flag), ideal(_ideal), upper(_upper), lower(_lower), prevene(0.0), remark(_remark) {}
    
    double calcenergyblength() const;
    void setprevene(const double &_prevene) const;

	static BondLengthMCeneobj create(float* _bleng, double* _atomenergy1, double* _atomenergy2, int* _flag, double* _ideal, double* _upper, double* _lower, string _remark) {
        return BondLengthMCeneobj(_bleng, _atomenergy1, _atomenergy2, _flag, _ideal, _upper, _lower, _remark);
    }
};

class BondAngleMCeneobj
{
    private:
        float* bangle;
        double* ideal;
        double* upper;
        double* lower;
        double *atomenergy1, *atomenergy2, *atomenergy3;
        int* flag;
        mutable double prevene;
        string remark;

    public:
    BondAngleMCeneobj(float* _bangle, double* _atomenergy1, double* _atomenergy2, double* _atomenergy3, int* _flag, double* _ideal, double* _upper, double* _lower, string _remark)
  		: bangle(_bangle), atomenergy1(_atomenergy1), atomenergy2(_atomenergy2), atomenergy3(_atomenergy3), flag(_flag), ideal(_ideal), upper(_upper), lower(_lower), prevene(0.0), remark(_remark) {}
    
	double calcenergybangle() const;
    void setprevene(const double &_prevene) const;

	static BondAngleMCeneobj create(float* _bangle, double* _atomenergy1, double* _atomenergy2, double* _atomenergy3, int* _flag, double* _ideal, double* _upper, double* _lower, string _remark) {
        return BondAngleMCeneobj(_bangle, _atomenergy1, _atomenergy2, _atomenergy3, _flag, _ideal, _upper, _lower, _remark);
    }
};

class BondTorsioneneobj
{
    private:
        float* btorsion;
        double *atomenergy1, *atomenergy2, *atomenergy3, *atomenergy4;
        int* flag;
        mutable double prevene;
        vector<double> *boundary;
        vector<double> *energy;

    public:
    BondTorsioneneobj(float* _btorsion, double* _atomenergy1, double* _atomenergy2, double* _atomenergy3, double* _atomenergy4, int* _flag, vector<double> *_boundary, vector<double> *_energy)
    : btorsion(_btorsion), atomenergy1(_atomenergy1), atomenergy2(_atomenergy2), atomenergy3(_atomenergy3), atomenergy4(_atomenergy4), flag(_flag), boundary(_boundary), energy(_energy), prevene(0.0) {}

    double calcenergybtorsion() const;
    void setprevene(const double &_prevene) const;

    static BondTorsioneneobj create(float* _btorsion, double* _atomenergy1, double* _atomenergy2, double* _atomenergy3, double* _atomenergy4, int* _flag, vector<double> *_boundary, vector<double> *_energy){
        return BondTorsioneneobj(_btorsion, _atomenergy1, _atomenergy2, _atomenergy3, _atomenergy4, _flag, _boundary, _energy);
    }
};

class CombTorsionseneobj
{
    private:
        float *btorsion1, *btorsion2;
        double *atomenergy11, *atomenergy12, *atomenergy13, *atomenergy14;
        double *atomenergy21, *atomenergy22, *atomenergy23, *atomenergy24;
        int *flag1, *flag2;
        mutable double prevene;
        vector<double> *boundary1, *boundary2;
        vector<double> *energy;
        int *binnum1, *binnum2;

    public:
    CombTorsionseneobj(float* _btorsion1, float* _btorsion2, double* _atomenergy11, double* _atomenergy12, double* _atomenergy13, double* _atomenergy14, double* _atomenergy21, double* _atomenergy22, double* _atomenergy23, double* _atomenergy24, int* _flag1, int* _flag2, vector<double> *_boundary1, vector<double> *_boundary2, vector<double> *_energy, int *_binnum1, int *_binnum2)
    : btorsion1(_btorsion1), btorsion2(_btorsion2), atomenergy11(_atomenergy11), atomenergy12(_atomenergy12), atomenergy13(_atomenergy13), atomenergy14(_atomenergy14), atomenergy21(_atomenergy21), atomenergy22(_atomenergy22), atomenergy23(_atomenergy23), atomenergy24(_atomenergy24), flag1(_flag1), flag2(_flag2), boundary1(_boundary1), boundary2(_boundary2), energy(_energy), binnum1(_binnum1), binnum2(_binnum2), prevene(0.0) {}

    double calcenergycombtorsion() const;
    void setprevene(const double &_prevene) const;

    static CombTorsionseneobj create(float *_btorsion1, float *_btorsion2, double* _atomenergy11, double* _atomenergy12, double* _atomenergy13, double* _atomenergy14, double* _atomenergy21, double* _atomenergy22, double* _atomenergy23, double* _atomenergy24, int* _flag1, int* _flag2, vector<double> *_boundary1, vector<double> *_boundary2, vector<double> *_energy, int *_binnum1, int *_binnum2){
        return CombTorsionseneobj(_btorsion1, _btorsion2, _atomenergy11, _atomenergy12, _atomenergy13, _atomenergy14, _atomenergy21, _atomenergy22, _atomenergy23, _atomenergy24, _flag1, _flag2, _boundary1, _boundary2, _energy, _binnum1, _binnum2);
    }
};

class BondDihedralMCeneobj
{
    private:
    float* bdangle;
        double KI, PN, Phase, iDiv;
        double* ideal;
        double* upper;
        double* lower;
        double *atomenergy1, *atomenergy2, *atomenergy3, *atomenergy4;
        int* flag;
        mutable double prevene;

    public:
    int nt1, nt2, nt3, nt4;
  	BondDihedralMCeneobj(float* _bdangle, double* _atomenergy1, double* _atomenergy2, double* _atomenergy3, double* _atomenergy4, int _nt4, int* _flag, double* _ideal, double* _upper, double* _lower, double _KI, double _PN, double _Phase, double _iDiv)
  		: bdangle(_bdangle), atomenergy1(_atomenergy1), atomenergy2(_atomenergy2), atomenergy3(_atomenergy3), atomenergy4(_atomenergy4), nt4(_nt4), flag(_flag), ideal(_ideal), upper(_upper), lower(_lower), KI(_KI), PN(_PN), Phase(_Phase), iDiv(_iDiv) {}
    
	double calcenergybdihedral() const;
    void setprevene(const double &_prevene) const;

	static BondDihedralMCeneobj create(float* _bdangle, double* _atomenergy1, double* _atomenergy2, double* _atomenergy3, double* _atomenergy4, int _nt4, int* _flag, double* _ideal, double* _upper, double * _lower, double _KI, double _PN, double _Phase, double _iDiv) {
        return BondDihedralMCeneobj(_bdangle, _atomenergy1, _atomenergy2, _atomenergy3, _atomenergy4, _nt4, _flag, _ideal, _upper, _lower, _KI, _PN, _Phase, _iDiv);
    }
};

class MonteCarloSimulation
{
private:
    int numbb;
    double temperature;
    int lsweepnum;
    //Logger* logger = nullptr;
    Logger logger = Logger();
    vector<Residueatoms> refinestr, refinestr_temp, initstr;
    vector<vector<Point3s*>> BasePtrsA, BasePtrsB, BasePtrs;
    vector<vector<Point3s*>> StackPtrsA, StackPtrsB;
    vector<vector<double*>> BaseAtomenergiesA, BaseAtomenergiesB, BaseAtomenergies;
    vector<vector<double*>> StackAtomenergiesA, StackAtomenergiesB;
    vector<vector<vector<int>>> stembps;
    map<int, vector<pair<int, int>>> loopnts;

    const FragmentStructure* MCinputtedstrallPtr;
    const vector<vector<bool>>* MCinputtedstrallflagPtr;
    const vector<FragmentStructure>* MCtemplatestrallPtr;
    const vector<vector<vector<bool>>>* MCtemplatestrallflagPtr;

    const vector<vector<vector<double>>>* MCinitdistmap;
    const vector<vector<vector<bool>>>* MCinitdistflag;
    const vector<vector<double>>* MCinitdistscore;
    //const vector<double>* MCinitntscore;
    const vector<vector<vector<vector<double>>>>* MCtemplatedistmap;
    const vector<vector<vector<vector<bool>>>>* MCtemplatedistflag;

    double bpscoresum;
    double stackscoresum;
    vector<vector<int>> bpmap;
    vector<vector<double>> cbpmap;
    vector<vector<double>> noncbpmap;
    vector<vector<double>> stmap;
    vector<vector<int>> bplist;
    vector<vector<int>> stacklist;
    vector<int> looplist;

    vector<vector<double>> originpoint;

    vector<int> ntflag;
    vector<vector<int>> mapflag;
    vector<vector<int>> mapbbflag;
    vector<vector<int>> lengthflag;
    vector<vector<int>> angleflag;
    vector<vector<int>> torsionflag;
    vector<vector<int>> hydlengthflag;
    vector<vector<int>> hydangleflag;
    vector<vector<int>> hydtorsionflag;
    vector<vector<int>> addlengthflag;
    vector<vector<int>> addangleflag;
    vector<vector<int>> hydaddlengthflag;
    vector<vector<int>> hydaddangleflag;

    list<ClashMCeneobj> clashmcene;
    list<ElectrosMCeneobj> electrosene;
    list<LJMCeneobj> ljene;
    list<HBondMCeneobj> hbondene;
    list<HBond2MCeneobj> hbond2ene;
    list<InitMCdistobj> initene;
    list<InitMCdistobj> templateene;
    list<InitMCPosobj> initposene;
    list<BasePairMCobj> bpene;
    list<BasePairMCobj2> bpene2;
    list<BasePairMCobj2> baseplane;
    list<BasePairMCobj3> bpene3;
    list<StackPairMCobj> stackene;
    list<StackPairGeoMCobj> stackgeoene;

    vector<double> clashenergies;
    vector<double> electrosenergies;
    vector<double> ljenergies;
    vector<double> hbondenergies;
    vector<double> hbond2energies;
    vector<double> initenergies;
    vector<double> initposenergies;
    vector<double> templateenergies;
    vector<double> bpenergies;
    vector<double> basepenergies;
    vector<double> bpenergies2;
    vector<double> stackenergies;
    vector<double> stackgeoenergies;

    list<BondLengthMCeneobj> bondlengthmcene;
    list<BondAngleMCeneobj> bondanglemcene;
    list<BondTorsioneneobj> bondtorsionmcene;
    list<BondLengthMCeneobj> addbondlengthmcene;
    list<BondLengthMCeneobj> addhydbondlengthmcene;
    list<BondAngleMCeneobj> addhydbondanglemcene;
    list<BondAngleMCeneobj> addbondanglemcene;
    list<CombTorsionseneobj> combtorsionmcene;

    vector<double> lengthenergies;
    vector<double> angleenergies;
    vector<double> torsionenergies;
    vector<double> combtorsionenergies;

    vector<vector<double>> lengthenergyatoms;
    vector<vector<double>> angleenergyatoms;
    vector<vector<double>> torsionenergyatoms;
    vector<vector<double>> pairwiseenergyatoms;

    //map<string, double> weightenergy= {{"Init", 1}, {"Length energy", 80}, {"Angle energy", 1.0}, {"Torsion energy", 1.0}, {"Comb torsion energy", 1.0}, {"Clash energy", 10.0}, {"Electros energy", 0.00001}, {"LJ energy", 0.000001}, {"Hbond energy", 1.0}, {"BP energy", 10}, {"Stack energy", 50}};
	//map<string, double> weightenergy= {{"Init", 0.1}, {"Template", 0.1}, {"Length energy", 70.0}, {"Angle energy", 1.0}, {"Torsion energy", 1.0}, {"Comb torsion energy", 1.0}, {"Clash energy", 50.0}, {"Electros energy", 0}, {"LJ energy", 0}, {"Hbond energy", 1.0}, {"BP energy", 15.0}, {"Stack energy", 15.0}};
	//map<string, double> weightenergy= {{"Init", 1.5}, {"Template", 1.5}, {"Length energy", 70.0}, {"Angle energy", 10.0}, {"Torsion energy", 1.0}, {"Comb torsion energy", 1.0}, {"Clash energy", 100.0}, {"Electros energy", 0}, {"LJ energy", 0}, {"Hbond energy", 1.0}, {"BP energy", 15.0}, {"Stack energy", 30.0}};
	//map<string, double> weightenergy= {{"InitPos", 0}, {"Init", 3.0}, {"Template", 3.0}, {"Length energy", 60.0}, {"Angle energy", 3.0}, {"Torsion energy", 1.0}, {"Comb torsion energy", 1.0}, {"Clash energy", 120.0}, {"Electros energy", 0}, {"LJ energy", 0}, {"Hbond energy", 1.0}, {"BP energy", 20.0}, {"Stack energy", 40.0}};
    //map<string, double> weightenergy= {{"InitPos", 0}, {"Init", 3.0}, {"Template", 3.0}, {"Length energy", 60.0}, {"Angle energy", 3.0}, {"Torsion energy", 0.5}, {"Comb torsion energy", 0.5}, {"Clash energy", 120.0}, {"Electros energy", 0}, {"LJ energy", 0}, {"Hbond energy", 3.0}, {"BP energy", 25.0}, {"Stack energy", 40.0}};
    
    //map<string, double> weightenergy= {{"InitPos", 0}, {"Init", 3.0}, {"Template", 3.0}, {"Length energy", 60.0}, {"Angle energy", 3.0}, {"Torsion energy", 0.1}, {"Comb torsion energy", 0.1}, {"Clash energy", 120.0}, {"Electros energy", 0}, {"LJ energy", 0}, {"Hbond energy", 20.0}, {"BP energy", 20.0}, {"Stack energy", 40.0}, {"Stack geometrical energy", 40.0}, {"Baseplane energy", 0.0}};
    
    map<string, double> weightenergy= {{"InitPos", 0}, {"Init", 3.0}, {"Template", 3.0}, {"Length energy", 60.0}, {"Angle energy", 1.0}, {"Torsion energy", 1.0}, {"Comb torsion energy", 0}, {"Clash energy", 120.0}, {"Electros energy", 0}, {"LJ energy", 0}, {"Hbond energy", 20.0}, {"BP energy", 0.0}, {"Stack energy", 30.0}, {"Stack geometrical energy", 0.0}, {"Baseplane energy", 10.0}};

    //map<string, double> weightenergy= {{"Init", 3.0}, {"Template", 3.0}, {"Length energy", 70.0}, {"Angle energy", 1.0}, {"Torsion energy", 1.0}, {"Comb torsion energy", 1.0}, {"Clash energy", 150.0}, {"Electros energy", 0}, {"LJ energy", 0}, {"Hbond energy", 1.0}, {"BP energy", 15.0}, {"Stack energy", 50.0}};

    //map<string, double> weightenergy= {{"Init", 1}, {"Length energy", 1}, {"Angle energy", 1}, {"Torsion energy", 1}, {"Comb torsion energy", 1}, {"Clash energy", 1}, {"Electros energy", 1}, {"LJ energy", 1}, {"Hbond energy", 1}, {"BP energy", 1}};

public:
    BasicFunc bf;
    void setInitStr(const FragmentStructure &pdb_entry, int c);
    void printStrParams();
    //void LoadParams(string pramdir);
    void setStr(const FragmentStructure &pdb_entry, int c);
    void Init(const ExtParams &extparameters, const FragmentStructure &inputtedstrall, const vector<vector<bool>> &inputtedstrallflag, const vector<FragmentStructure> &templatestrall, const vector<vector<vector<bool>>> &templatestrallflag, const vector<vector<vector<double>>> &initdistmap, const vector<vector<vector<bool>>> &initdistflag, const vector<vector<double>> &initdistscores, const vector<vector<vector<vector<double>>>> &templatedistmap, const vector<vector<vector<vector<bool>>>> &templatedistflag, const vector<vector<int>> &bindmap, const vector<vector<pair<int, double>>> &stmap, const vector<vector<double>> &cbpmap_input, const vector<vector<double>> &noncbpmap_input);
    void returnStr(FragmentStructure &pdb_entry, int c);
    void setOrigin();
    void resetfrag(int index=0);
    void StrCopy(vector<Residueatoms> strA, vector<Residueatoms> &strB);
    void getmolConnect(int res1, int res2, vector<vector<int>> &molConnect, map<string, vector<vector<bool>>> &intraAtomConnect, const ExtParams &extparameters);
    void checkupdate();
    
    //Simulation main
    double Simulation(double maxtime, int startepochidx);

    //Energy main
    double calcenergymain();
    void initializeatomenergy();

    //pair-wise energies
    bool checkhbond2(Point3s* teAtomHyd1, Point3s* teAtomHyd2, Point3s* teAtomDon, Point3s* teAtomAcc);
    bool checkhbond(Point3s* teAtomHyd, Point3s* teAtomDon, Point3s* teAtomAcc);
    void setpairwiseeneobj(const ExtParams &extparameters, const FragmentStructure &inputtedstrall, const vector<vector<bool>> &inputtedstrallflag, const vector<FragmentStructure> &templatestrall, const vector<vector<vector<bool>>> &templatestrallflag);
    //electros
    double calceneelectros();
    void updateeneelectros();
    //clash
    double calceneclash();
    void updateeneclash();
    //LJ potential
    double calcenelj();
    void updateenelj();
    //hbond
    double calcenehbond();
    void updateenehbond();
    //Init
    double calceneinit();
    void updateeneinit();
    //Init pos
    double calceneinitpos();
    void updateeneinitpos();
    //Template
    double calcenetemplate();
    void updateenetemplate();
    //Base pair
    double calcenebp();
    void updateenebp();
    //Base plane
    double calcenebaseplane();
    void updateenebaseplane();
    //Stacking template
    double calcenestack();
    void updateenestack();
    //Stacking geometries
    double calcenestackgeo();
    void updateenestackgeo();

    //geometrical energies
    void setgeometricaleneobj();
    //length
    double calcenelength();
    void updateenelength();
    //angle
    double calceneangle();
    void updateeneangle();
    //torsion
    double calcenetorsion();
    void updateenetorsion();

    //combination torsion
    double calcenecombtorsion();
    void updateenecombtorsion();

    //movement 
    bool calcmovementmain(double oldenergy, double &newenergy, int &moveidx);
    void updateenergies();
    void setsimulationparams(double _temperature, int _lsweepnum);
    double gettemperature();
    bool metropolisUpdate(double oldenergy, double newenergy);
    double getrandomnewparam(double l, double r);
    int getrandval(vector<double> tmov,double trandnum);
    //angle
    bool fmovementang(double oldenergy, double &newenergy);
    void selecttargetang(int &trandpos, int &angtype, double r);
    double getmovangval(double val, double mean, double std, double k);
    void setangntflap(double trandpos, double base, int angtype);
    void setangflap( double trandpos, double base, int angtype);
    //length
    bool fmovementleng(double oldenergy, double &newenergy);
    void selecttargetleng(int &trandpos, int &lengtype, double r);
    double getmovlengval(double val, double mean, double std, double k);
    void setlengflap(double trandpos, double base, int lengtype);
    void setlengntflap(double trandpos, double base, int lengtype);
    //torsion
    bool fmovementtor(double oldenergy, double &newenergy);
    void selecttargettor(int &trandpos, int &tortype, double r);
    double sampletorsion(vector<double> &logdistr, double logtotal, vector<double> &center, double span);
    void settorflap(double trandpos, double base, int tortype);
    void settorntflap(double trandpos, double base, int tortype);

    //combination torsion
    bool fmovementcombtor(double oldenergy, double &newenergy);
    void selecttargetcombtor(int &trandpos, int &tortype, double r);
    double samplesecondtor2bin(vector<double> &prob, double totalprob);
    double samplesecondtor2(vector<double> &center, double span, int binidx);
    void setcombtorflap(double trandpos, double base);

    //LMP
    bool fmovementLMP(double oldenergy, double &newenergy);
    int selecttargetposLMP(int &trandpos, double r);
    void setLMPflap(double trandpos, double base);
    void singlemoveLMP(int k,int flagdir);

    //Rotate
    bool fmovementrot(double oldenergy, double &newenergy);
    void setrotateflap(double trandpos1, double trandpos2);

    //fragment optimization
    void MCprepgrad(vector<vector<Point3s>> &gradptv);
    void MCfragmentoptimization(
        const map<int, vector<Point3s>> &initstructure, const vector<vector<int>> &bindingmap,
        const vector<bool> &moveresidues, const vector<vector<int>> &contactmap, const vector<vector<double>> &stackingmap,
        const vector<vector<int>> &basepairlist, const vector<vector<int>> &noncbasepairlist,
        const vector<double> &basepairscorelist, const vector<double> &noncbasepairscorelist);
    void MCsetinitobj(const map<int, vector<Point3s>> &initstructure, vector<vector<Point3s>> &gradient, list<InitPosobj> &binitobj, double K);
    void MCsetinitpairobj(const map<int, vector<Point3s>> &initstructure, const vector<bool> &moveresidues, vector<vector<Point3s>> &gradient, list<InitPairobj> &binitdistobj, double K);
    void MCsettemplatepairobj(const vector<bool> &moveresidues, vector<vector<Point3s>> &gradient, list<InitPairobj> &binitdistobj, double K);
    void MCsetlengthobj(const vector<vector<int>> &bindingmap, vector<vector<Point3s>> &gradient, double K, list<BondLengthobj> &blengthobj);
    void MCsetangleobj(const vector<vector<int>> &bindingmap, vector<vector<Point3s>> &gradient, double K, list<BondAngleobj> &bangleobj);
    
    //void MCsetstackingobj(vector<vector<Point3s>> &gradient, const vector<vector<int>> &stackmap, double K, list<Stackingobj> &stackobj);
    double getstackenergy(const int &stack_base1, const int &stack_base2, const int &sttype, vector<double> &templateenergies);
    double getstackenergy(const int &stack_base1, const int &stack_base2, const int &sttype, vector<double> &templateenergies, int &mintempidx, int &minorder);
    
    void selectstacktemplateidx(const int &stack_base1, const int &stack_base2, const int &sttype, const vector<double> &templateenergies, int &templateidx, int &temporder);
    void MCselectstackfrag(vector<vector<Point3s>> &gradient, int res1, int res2, const int &sttype, int &minfragidx, double &minrmsd, int &baseorder, double K, list<Stackingobj> &stackobj);

    void MCselectstackfrag(vector<vector<Point3s>> &gradient, int res1, int res2, const int &sttype, int stidx, vector<vector<int>> &indices, vector<vector<int>> &orders, double K, list<Stackingobj> &stackobj);
    void MCsetstackingobj(vector<vector<Point3s>> &gradient, const vector<vector<double>> &stackmap, vector<vector<int>> &indices, vector<vector<int>> &orders, double K, list<Stackingobj> &stackobj);
    double MClinesearch(const vector<vector<Point3s>> &gradient, list<BondLengthobj> &blengthobj, list<BondAngleobj> &bangleobj, list<BondDihedralobj> &dihedralobj, list<BasePairobj> &bplaneobj, list<BasePairobj> &basepairobj, list<BPdistobj> &bpdistobj, list<Stackingobj> &stackobj, list<StackGeoobjMC> &stackgeoobj, list<InitPosobj> &binitobj, list<InitPairobj> &binitdistobj, list<InitPairobj> &btemplatedistobj, list<TElectrobj> &telecteobj, list<LennardsJobj> &ljobj, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, list<Clashobj> &clashobj, list<SingleNtTemplateobjMC> &nttempobj, const vector<bool> &moveresidues, double alpha0, double rho, double beta);

    double MCgoldenelectionsearch(const vector<vector<Point3s>> &gradient, list<BondLengthobj> &blengthobj, list<BondAngleobj> &bangleobj, list<BondDihedralobj> &dihedralobj, list<BasePairobj> &bplaneobj, list<BasePairobj> &basepairobj, list<BPdistobj> &bpdistobj, list<Stackingobj> &stackobj, list<StackGeoobjMC> &stackgeoobj, list<InitPosobj> &binitobj, list<InitPairobj> &binitdistobj, list<InitPairobj> &btemplatedistobj, list<TElectrobj> &telecteobj, list<LennardsJobj> &ljobj, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, list<Clashobj> &clashobj, list<SingleNtTemplateobjMC> &nttempobj, const vector<bool> &moveresidues);

    double getbaseenergy(const int &nt_base, vector<double> &templateenergies, int &templateidx);
    void selectnttemplateidx(const int &nt_base, const vector<double> &templateenergies, int &templateidx);

    void selectbptemplateidx(const int &bp_base1, const int &bp_base2, const vector<double> &templateenergies, int &templateidx);
    double getbasepairenergy(const int &bp_base1, const int &bp_base2, vector<double> &templateenergies, int &templateidx);
    void MCsetbaseplaneobj(vector<vector<Point3s>> &gradient, vector<vector<Point3s*>> &BaseAtomPtrslist, vector<vector<Point3s*>> &BaseAtomgradlist, const vector<bool> &moveresidues, const double &K, list<BasePairobj> &bplaneobj);
    void MCsetbasepairobj(vector<vector<Point3s>> &gradient, const vector<vector<int>> &basepairlist, const vector<vector<int>> &noncbasepairlist, const vector<double> &basepairscorelist, const vector<double> &noncbasepairscorelist, vector<vector<Point3s*>> &AtomPtrslist, vector<vector<Point3s*>> &Atomgradlist, const double &K, list<BasePairobj> &bpobj);
    void MCsetbpdistobj(vector<vector<Point3s>> &gradient, const vector<vector<int>> &basepairlist, vector<vector<int>> &indices, double K, list<BPdistobj> &stackobj);
    void MCsetNonbondedobj(vector<vector<Point3s>> &gradient, const vector<vector<int>> &contactmap, const vector<vector<int>> &bindingmap, const vector<vector<int>> &cbasepairlist, const vector<vector<int>> &noncbasepairlist, const vector<double> &basepairscorelist, const vector<double> &noncbasepairscorelist,list<TElectrobj> &telecteobj, list<LennardsJobj> &ljobj, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, list<Clashobj> &clashobj, map<string, double> &enecoeff);
    void MCupdatelbfgsfeatures(const vector<vector<Point3s>> &refinestr_temp_prev, const vector<vector<Point3s>> &grad, const vector<vector<Point3s>> &grad_prev, deque<vector<vector<Point3s>>> &s_list, deque<vector<vector<Point3s>>> &y_list, deque<double> &rho_list, const vector<bool> &movenodeinfo, int m);
    void MCupdatestructure(const vector<vector<Point3s>> &gradient, double step_size, const vector<bool> &moveresidues);
    double MClinesearch(const vector<vector<Point3s>> &gradient, list<BondLengthobj> &blengthobj, list<BondAngleobj> &bangleobj, list<BondDihedralobj> &dihedralobj, list<BasePairobj> &bplaneobj, list<BasePairobj> &basepairobj, list<BPdistobj> &bpdistobj, list<StackTemplateobj2MC> &stackobj, list<StackGeoobjMC> &stackgeoobj, list<InitPosobj> &binitobj, list<InitPairobj> &binitdistobj, list<InitPairobj> &btemplatedistobj, list<TElectrobj> &telecteobj, list<LennardsJobj> &ljobj, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, list<Clashobj> &clashobj, list<SingleNtTemplateobjMC> &nttempobj, const vector<bool> &moveresidues, double alpha0, double rho, double beta);
    double MCgoldenelectionsearch(const vector<vector<Point3s>> &gradient, list<BondLengthobj> &blengthobj, list<BondAngleobj> &bangleobj, list<BondDihedralobj> &dihedralobj, list<BasePairobj> &bplaneobj, list<BasePairobj> &basepairobj, list<BPdistobj> &bpdistobj, list<StackTemplateobj2MC> &stackobj, list<StackGeoobjMC> &stackgeoobj, list<InitPosobj> &binitobj, list<InitPairobj> &binitdistobj, list<InitPairobj> &btemplatedistobj, list<TElectrobj> &telecteobj, list<LennardsJobj> &ljobj, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, list<Clashobj> &clashobj, list<SingleNtTemplateobjMC> &nttempobj, const vector<bool> &moveresidues);
    
    bool fmovementfrag(double oldenergy, double &newenergy);
    int selecttargetfrag(const int &fraglen);
    void setnonadjfragflap(int istart1, int iend1, int istart2, int iend2);
    void setfragflap(int istart, int iend);
    void setadjfragflap(int istart, int iend);

    void MCstack_superimpose(int res1, int res2, char base1, char base2, const int &sttype, double &rmsd, const int &order, const int &tmpidx);
    void MCselectstackfrag(vector<vector<Point3s>> &gradient, int res1, int res2, const int &sttype, int stidx, vector<vector<int>> &indices, vector<vector<int>> &orders, double K, list<StackTemplateobj2MC> &stackobj);
    void MCsetstackingobj(vector<vector<Point3s>> &gradient, const vector<vector<double>> &stackmap, vector<vector<int>> &indices, vector<vector<int>> &orders, double K, list<StackTemplateobj2MC> &stackobj);
    void MCsetstackinggeoobj(vector<vector<Point3s>> &gradient, const vector<vector<double>> &stackmap, vector<vector<Point3s*>> &basestackptrsA, vector<vector<Point3s*>> &basestackptrsB, vector<vector<Point3s*>> basestackgradA, vector<vector<Point3s*>> basestackgradB, double K, list<StackGeoobjMC> &stackobj);

    void MCbp_superimpose(int res1, int res2, char base1, char base2, double &rmsd, const int &tmpidx);
    void MCselectbpfrag(vector<vector<Point3s>> &gradient, int res1, int res2, int bpidx, vector<vector<int>> &indices, double K, list<BPdistobj> &bpdistobj);

    void MCsetnttemplateobj(vector<vector<Point3s>> &gradient, vector<vector<int>> &indices, double K, list<SingleNtTemplateobjMC> &nttempobj);
    void MCselectntfrag(vector<vector<Point3s>> &gradient, int res, int ntidx, vector<vector<int>> &indices, double K, list<SingleNtTemplateobjMC> &nttempobj);
    void MCnt_superimpose(int res, char base, double &rmsd, const int &tmpidx);
    //void MCselectntfrag(vector<vector<Point3s>> &gradient, int res, int ntidx, vector<vector<int>> &indices, double K, list<SingleNtTemplateobjMC> &nttempobj);    

    //Movement for stack
    int selecttargetstack(int &basei, int &basej, double &bpscore);
    bool replacestackfrag(int &stack_base1, int &stack_base2, int &sttype, int &templateidx, int &temporder);
    bool fmovementstackfrag(double oldenergy, double &newenergy);
    void setadjstackfragflap(int istart, int iend);
    void setnonadjstackfragflap(int istart1, int iend1, int istart2, int iend2);

    //Movement for base pair
    //int selecttargetbp();
    int selecttargetbp(int &basei, int &basej, double &bpscore);
    bool fmovementbpfrag(double oldenergy, double &newenergy);
    bool replacebpfrag(int &bp_base1, int &bp_base2, int &selected_fragidx);
    void setbpfragflap(int istart1, int iend1, int istart2, int iend2);

    //movement for loop
    int selecttargetloopnt();
    bool fmovementloopfrag(double oldenergy, double &newenergy);
    bool replacentfrag(int &base, int &selected_fragidx);
    void setloopntfragflap(int istart, int iend);

    //coord2param
    void coord2params(const Point3s &tmp1, const Point3s &tmp2, const Point3s &tmp3, const Point3s &tmp4, float &len, float &ang, float &dih);
    void coord2dangle(const Point3s &tmp1, const Point3s &tmp2, const Point3s &tmp3, const Point3s &tmp4, float &dih);
    void coord2anglen(const Point3s &tmp2, const Point3s &tmp3, const Point3s &tmp4, float &len, float &ang);
    void coord2ang(const Point3s &tmp2, const Point3s &tmp3, const Point3s &tmp4, float &ang);
    void coord2len(const Point3s &tmp3, const Point3s &tmp4, float &len);
    void coord2parambb(vector<Residueatoms> &decstr,int start, int end);
    void coord2paramsugar(vector<Residueatoms> &decstr,int start,int end);
    void coord2paramadd(vector<Residueatoms> &decstr,int start,int end);
    void coord2paramhydrogen(vector<Residueatoms> &decstr,int start,int end);
    void coord2addparamhydrogen(vector<Residueatoms> &decstr,int start,int end);
    void coord2paramhydrogenBase(vector<Residueatoms> &decstr,int start,int end);
    void coord2parambase(vector<Residueatoms> &decstr,int start,int end);

    //param2coord
    bool param2coordbb(vector<Residueatoms> &decstr,int istart,int iend);
    void param2coordsugar(vector<Residueatoms> &decstr,int lstart,int lend);
    bool param2coordbase(vector<Residueatoms> &decstr, int istart, int iend);
    void param2coordhydrogen(vector<Residueatoms> &decstr,int lstart,int lend);

    //logger
    void setLogger(bool logflag, string logfile="");

    MonteCarloSimulation(const FragmentStructure &pdb_entry, const int &c, const ExtParams &extparameters, const FragmentStructure &inputtedstrall, const vector<vector<bool>> &inputtedstrallflag, const vector<FragmentStructure> &templatestrall, const vector<vector<vector<bool>>> &templatestrallflag, const vector<vector<vector<double>>> &initdistmap, const vector<vector<vector<bool>>> &initdistflag, const vector<vector<double>> &initdistscores, const vector<vector<vector<vector<double>>>> &templatedistmap, const vector<vector<vector<vector<bool>>>> &templatedistflag, const vector<vector<int>> &bindmap, const vector<vector<pair<int, double>>> &stmap, const vector<vector<double>> &cbpmap_input, const vector<vector<double>> &noncbpmap_input);
    
    MonteCarloSimulation(const MonteCarloSimulation& other) {
        refinestr = other.refinestr;
        refinestr_temp = other.refinestr_temp;
        initstr = other.initstr;
        stembps = other.stembps;
        loopnts = other.loopnts;
        clashmcene = other.clashmcene;
    }
    //virtual ~MonteCarloSimulation();
};

int getbin(double val, const vector<double> &boundary);
void LoadParams(string pramdir);
double MonteCarloSimulationProcess(ModelUnit &pdb_entry, const size_t c, const string &pramdir,
    const vector<pair<double,vector<size_t>>> &bp_vec, const map<int, vector<pair<int, int>>> &loopnts, const vector<pair<double,vector<size_t>>> &noncbp_vec,
    const vector<vector<int>> &stemconnectnts, const map<int, vector<vector<vector<int>>>> &connectnts, const vector<vector<pair<int, double>>> &stmap,
    const FragmentStructure &inputtedstrall, const vector<vector<bool>> &inputtedstrallflag, const vector<double> &score_vec,
    const vector<FragmentStructure> &templatestrall, const vector<vector<vector<bool>>> &templatestrallflag,
    const ExtParams &extparameters, vector<vector<vector<double>>> &hydrogencoords, double maxtime, const string &outdir, const bool &logflag, const int option=0);

double MonteCarloSimulationProcess(ModelUnit &pdb_entry, const string &pramdir, 
    const vector<pair<double,vector<size_t>>> &bp_vec, const map<int, vector<pair<int, int>>> &loopnts, const vector<pair<double,vector<size_t>>> &noncbp_vec,
    const vector<vector<int>> &stemconnectnts, const map<int, vector<vector<vector<int>>>> &connectnts, const vector<vector<pair<int, double>>> &stmap,
    const FragmentStructure &inputtedstrall, const vector<vector<bool>> &inputtedstrallflag, const vector<double> &score_vec,
    const vector<FragmentStructure> &templatestrall, const vector<vector<vector<bool>>> &templatestrallflag,
    const ExtParams &extparameters, vector<vector<vector<double>>> &hydrogencoords, double maxtime, const string &outdir, const bool &logflag, const int option=0);

#endif
