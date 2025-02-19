#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
//#include <string.h>
#include <math.h>
#include <sys/stat.h>
//#include <vector>
using namespace std;

#pragma once

#define currentsystem 1//0 windows 1 linux
#if currentsystem==0
#include <windows.h>
#else
#include <sys/time.h>
#endif

#define MAX_CHAIN_NUM 128
#define PARAM_NUM 15
#define ADDLENGTH_NUM 1
#define ADDINNERANGLE_NUM 8
#define ATOM_RES_NUM 34
#define HYD_PARAM_NUM 11
#define ATOM_TYPE_NUM 18
//#define ADD_ATOM_NUM 12
//#define BS_ATOM_NUM 14
#define EPS 1e-7
#define PI (4.0 * atan(1.0))
#define RADDEG PI/180.0
#define RAD2DEG 180.0/M_PI
#define DIST_THRESH 900
#define CONTACT_THRESH 0.5
#define MRANGETYPE 0

#define PAIRLENGTHTARGETNUM 15
#define PAIRLENGTHTARGETMOVNUM 15
#define CLSTORPARAMNUM 6

#define STREAMS    256 
#define DEFAULT    123456789  /* initial seed, use 0 < DEFAULT < MODULUS  */
#define NN 624
#define MODULUS    2147483647 /* DON'T CHANGE THIS VALUE 				  */
#define MULTIPLIER 48271      /* DON'T CHANGE THIS VALUE                  */
#define A256       22925      /* jump multiplier, DON'T CHANGE THIS VALUE */

#define HBONDCENTERATOM 5
//#define ENERGY_TYPE_NUM 200
#define TOTAL_ENERGY_NUM 200
#define BLENGTH_ENERGY_NUM 15
#define ADDBLENGTH_ENERGY_NUM 1
#define INNERANGLE_ENERGY_NUM 15
#define ADDINNERANGLE_ENERGY_NUM 8
#define TORSIONANGLE_ENERGY_NUM 15
#define CLASH_ENERGY_NUM 1
#define BASE_ENERGY_NUM 4
#define BASEHYD_ENERGY_NUM 2
#define STACKING_ENERGY_NUM 3
#define CONTACT_ENERGY_NUM 2
#define DIST_ENERGY_NUM 2
#define TORSION_CLS_ENERGY_NUM 2
#define PUCKER_ENERGY_NUM 1
#define ADJPAIRLENG_ENERGY_NUM 1
#define NTPAIRLENG_ENERGY_NUM 1
#define MULTITORSION_ENERGY_NUM 8
#define PAIRTORSION_ENERGY_NUM 6
#define PAIRBBTORSION_ENERGY_NUM 5
#define OXHBOND_ENERGY_NUM 1
#define SINGLENTTEMP_ENERGY_NUM 1
#define STACKINGTEMP_ENERGY_NUM 1
#define BASEPAIRTEMP_ENERGY_NUM 1

#define BLENGTH_ENERGY_ST 0
#define ADDBLENGTH_ENERGY_ST 15
#define INNERANGLE_ENERGY_ST 16
#define ADDINNERANGLE_ENERGY_ST 31
#define TORSIONANGLE_ENERGY_ST 39
#define CLASH_ENERGY_ST 54
#define BASE_ENERGY_ST 55
#define BASEHYD_ENERGY_ST 59
#define STACKING_ENERGY_ST 61
#define CONTACT_ENERGY_ST 64
#define DIST_ENERGY_ST 66
#define TORSION_CLS_ENERGY_ST 68
#define PUCKER_ENERGY_ST 70
#define ADJPAIRLENG_ENERGY_ST 71
#define NTPAIRLENG_ENERGY_ST 72
#define MULTITORSION_ENERGY_ST 73
#define PAIRTORSION_ENERGY_ST 81
#define PAIRBBTORSION_ENERGY_ST 87
#define OXHBOND_ENERGY_ST 92
#define SINGLENTTEMP_ENERGY_ST 93
#define STACKINGTEMP_ENERGY_ST 94
#define BASEPAIRTEMP_ENERGY_ST 95

static unsigned long mt[NN]; /* the array for the state vector  */
static int mti=NN+1; /* mti==N+1 means mt[N] is not initialized */
static long seed[STREAMS] = {DEFAULT};  /* current state of each stream   */
static int  stream        = 0;          /* stream index, 0 is the default */
static int  initialized   = 0;          /* test for stream initialization */

//static map< string, double> outlinerstd={
	//{"length", 8.0}, {"inner", 8.0}, {"torsion", 4.0}
//};

//static map< string, double> outlinerstdsugar={
	//{"length", 5.0}, {"inner", 5.0}, {"torsion", 4.0}
//};

static map<string, map<string, vector<double>>>outlinethread={
	{
		"geo",{
			{"length",{5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0}}, 
			{"inner",{5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0}}, 
			{"torsion",{5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0}},
			{"oppleng",{5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0}}
			}
	},
	{
		"sgeo",{
			{"length",{12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0}}, 
			{"inner",{12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0}}, 
			{"torsion",{12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0}},
			{"oppleng",{12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0}}
			}
	},
	{
		"geo_add",{
			//{"length", {4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0}}, 
			//{"inner", {4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0}},
			//{"torsion", {4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0}}
			{"length", {5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0}}, 
			{"inner", {5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0}},
			{"torsion", {5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0}}
			}
	},
	{
		"sgeo_add",{
			//{"length", {4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0}}, 
			//{"inner", {4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0}},
			//{"torsion", {4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0}}
			{"length", {12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0}}, 
			{"inner", {12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0}},
			{"torsion", {12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0}}
			}
	},
	{
		"geo_first",{
			{"length", {5.0, 5.0}}, 
			{"inner", {5.0, 5.0}}, 
			{"torsion", {5.0, 5.0}}
			}
	},
	{
		"sgeo_first",{
			{"length", {12.0, 12.0}}, 
			{"inner", {12.0, 12.0}}, 
			{"torsion", {12.0, 12.0}}
			}
	},
	{
		"geo_hyd",{
			{"length", {5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0}}, 
			{"inner", {5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0}}, 
			{"torsion", {5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0}}
			}
	},
	{
		"sgeo_hyd",{
			{"length", {12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0}}, 
			{"inner", {12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0}}, 
			{"torsion", {12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0}}
			}
	},
	{
		"cls_tor",{
			{"length", {3.0, 3.0, 3.0, 3.0, 3.0, 3.0}}, 
			{"inner", {3.0, 3.0, 3.0, 3.0, 3.0, 3.0}}, 
			{"torsion", {3.0, 3.0, 3.0, 3.0, 3.0, 3.0}}
			}
	}
};

static map< string, double> outlinerbpstd={
	{"length", 4}, {"inner", 4}, {"torsion", 1.0}
};

static map< string, double> outlinerbphydstd={
	{"length", 4}, {"angle", 4}
};

static map< string, double> outlinerbsstd={
	{"length", 4}, {"inner", 4}, {"torsion", 4}, {"baseangle", 2.5}
};

static map< string, double> outlinerpucker{
	{"length", 4}, {"torsion", 4}
};

static map< string, double> outlineroshyd{
	{"length", 1}, {"angle", 1}
};

static map< string, double> outpairlengstd={{"nt", 4.0}, {"adj", 4.0}};

static map< int, double> outlinertorcls={{0, 3.0}, {1, 3.0}};

static map< string, int> Atom2id={
	{" P  ", 0},{" O5'", 1},{" C5'", 2},{" C4'", 3},{" C3'", 4},{" O3'", 5},{" C2'", 6},{" O4'", 7},{" C1'", 8},{" O2'",9},{" OP1",10},{" OP2",11},{" N9 ", 12},{" N1 ",13},
	{" C4 ",14},{" C2 ",15},{" C8 ",16},{" C6 ",17},{" N7 ",18},{" C5 ",19},{" N6 ",20},{" N3 ",21},{" O2 ",22},{" N4 ",23},{" O6 ",24},{" N2 ",25},{" O4 ",26},
	{" H5'",27},{"H5''",28},{" H4'",29},{" H3'",30},{" H2'",31},{"HO2'",32},{" H1'",33},{" H8 ",34},{" H61",35},{" H62",36},{" H2 ",37},{" H6 ",38},{" H5 ",39},{" H3 ",40},{" H21",41},
	{" H22",42},{" H1 ",43},{" H41",44},{" H42",45},{"  H8",34},{"  H2",37},{"  H6",38},{"  H5",39},{"  H3",40},{"  H1",43}
};

static map<char, map<string, int>> Atom2cid={
	{
		'A',{{" P  ", 0},{" O5'", 1},{" C5'", 2},{" C4'", 3},{" C3'", 4},{" O3'", 5},{" C2'", 6},{" O4'", 7},{" C1'", 8},{" O2'",9},{" OP1",10},{" OP2",11},{" N9 ",12},{" C4 ",13},{" C8 ",14}}
	},
	{
		'C',{{" P  ", 0},{" O5'", 1},{" C5'", 2},{" C4'", 3},{" C3'", 4},{" O3'", 5},{" C2'", 6},{" O4'", 7},{" C1'", 8},{" O2'",9},{" OP1",10},{" OP2",11},{" N1 ",12},{" C2 ",13},{" C6 ",14}}
	},
	{
		'G',{{" P  ", 0},{" O5'", 1},{" C5'", 2},{" C4'", 3},{" C3'", 4},{" O3'", 5},{" C2'", 6},{" O4'", 7},{" C1'", 8},{" O2'",9},{" OP1",10},{" OP2",11},{" N9 ",12},{" C4 ",13},{" C8 ",14}}
	},
	{
		'U',{{" P  ", 0},{" O5'", 1},{" C5'", 2},{" C4'", 3},{" C3'", 4},{" O3'", 5},{" C2'", 6},{" O4'", 7},{" C1'", 8},{" O2'",9},{" OP1",10},{" OP2",11},{" N1 ",12},{" C2 ",13},{" C6 ",14}}
	}
};

static map<int, string> Id2atom={
	{0," P  "},{1," O5'"},{2," C5'"},{3," C4'"},{4," C3'"},{5," O3'"},{6," C2'"},{7," O4'"},{8," C1'"},{9," O2'"},{10," OP1"},{11," OP2"},{12," N9 "},{13," N1 "},{14," C4 "},{15," C2 "},
	{16," C8 "},{17," C6 "},{18," N7 "},{19," C5 ",},{20," N6 "},{21," N3 "},{22," O2 "},{23," N4 "},{24," O6 "},{25," N2 "},{26," O4 "},
	{27," H5'"},{28,"H5''"},{29," H4'"},{30," H3'"},{31," H2'"},{32,"HO2'"},{33," H1'"},{34," H8 "},{35," H61"},{36," H62"},{37," H2 "},{38," H6 "},{39," H5 "},{40," H3 "},{41," H21"},
	{42," H22"},{43," H1 "},{44," H41"},{45," H42"}
};

//static map< string, int> HeavyAtom2id={
	//{" P  ", 0},{" O5'", 1},{" C5'", 2},{" C4'", 3},{" C3'", 4},{" O3'", 5}
//};

//static map< int, string> id2HeavyAtom={
	//{0, " P  "},{1, " O5'"},{2, " C5'"},{3, " C4'"},{4, " C3'"},{5, " O3'"}
//};

static map< string, char> PdbRes2ResIdchar={
	{"  A",'0'},{"  C",'1'},{"  G",'2'},{"  U",'3'},{" DA",'0'},{" DC",'1'},{" DG",'2'},{" DT",'3'},{"RA3",'0'},{"RC3",'1'},{"RG3",'2'},{"RU3",'3'}
};

static map< string, string> nonstPdbRes2stPdbRes={
	{" DA","  A"},{" DC","  C"},{" DG","  G"},{" DT","  U"},{"RA3","  A"},{"RC3","  C"},{"RG3","  G"},{"RU3","  U"}
};

static map< char, char> ResIdchar2Reschar={
	{'0','A'},{'1','C'},{'2','G'},{'3','U'},{'4','X'}
};

static map< int, char> Idint2Reschar={
	{0,'A'},{1,'C'},{2,'G'},{3,'U'},{4,'X'}
};

static map< char, int> Reschar2Idint={
	{'A',0},{'C',1},{'G',2},{'U',3},{'X',4}
};

static map< char, const char*> Ntd3={
	{0, "  A"},{1, "  C"},{2, "  G"},{3, "  U"},{4,"  X"}
};

static vector<int> DISTTARGETATOMS={0,3};
static int CTARGETATOM = 1;

static vector<vector<int>> NNheavy={
	{33,0,1,2,3,4,5,6,7,8,9,10,11,12,14,16,13,15,17,18,19,20,21,27,28,29,30,31,32,33,34,35,36,37},
	{31,0,1,2,3,4,5,6,7,8,9,10,11,13,15,17,14,19,21,22,23,27,28,29,30,31,32,33,38,39,44,45},
	{34,0,1,2,3,4,5,6,7,8,9,10,11,12,14,16,13,15,17,18,19,21,24,25,27,28,29,30,31,32,33,34,41,42,43},
	{30,0,1,2,3,4,5,6,7,8,9,10,11,13,15,17,14,19,21,22,26,27,28,29,30,31,32,33,38,39,40}
};

static vector<vector<int>> NNheavyIdx={
	{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29},
	{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27},
	{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30},
	{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27}
};

static vector<vector<int>> SaveNN={
	{22,0,1,2,3,4,5,6,7,8,9,10,11,12,14,16,13,15,17,18,19,20,21},
	{20,0,1,2,3,4,5,6,7,8,9,10,11,13,15,17,14,19,21,22,23},
	{23,0,1,2,3,4,5,6,7,8,9,10,11,12,14,16,13,15,17,18,19,21,24,25},
	{20,0,1,2,3,4,5,6,7,8,9,10,11,13,15,17,14,19,21,22,26}
};

static vector<vector<int>> NNAll={
	{33,0,1,2,3,4,5,6,7,8,9,10,11,12,14,16,13,15,17,18,19,20,21,27,28,29,30,31,32,33,34,35,36,37},
	{31,0,1,2,3,4,5,6,7,8,9,10,11,13,15,17,14,19,21,22,23,27,28,29,30,31,32,33,38,39,44,45},
	{34,0,1,2,3,4,5,6,7,8,9,10,11,12,14,16,13,15,17,18,19,21,24,25,27,28,29,30,31,32,33,34,41,42,43},
	{30,0,1,2,3,4,5,6,7,8,9,10,11,13,15,17,14,19,21,22,26,27,28,29,30,31,32,33,38,39,40}
};

static map<char, vector<int>>Atomoxhydbtarget={
	{'A',{0,1,0,0,0,1,0,1,0,1,1,1,1,0,0,1,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0}},
	{'C',{0,1,0,0,0,1,0,1,0,1,1,1,1,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0}},
	{'G',{0,1,0,0,0,1,0,1,0,1,1,1,1,0,0,1,0,0,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0}},
	{'U',{0,1,0,0,0,1,0,1,0,1,1,1,1,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0}}
};

static map<char, vector<int>> BaseHeavyAtom={
	{'A',{12,13,14,15,16,17,18,19,20,21}},
	{'C',{12,13,14,15,16,17,18,19}},
	{'G',{12,13,14,15,16,17,18,19,20,21,22}},
	{'U',{12,13,14,15,16,17,18,19}}
};

static map<char, map<int, int>> AtomOxhydbId2pId={
	{'A',{{1,0}, {5,1}, {7,2}, {9,3}, {10,4}, {11,4}, {15,5}, {21,6}, {18,7}, {20,8}}},
	{'C',{{1,0}, {5,1}, {7,2}, {9,3}, {10,4}, {11,4}, {18,5}, {17,6}, {19,7}}},
	{'G',{{1,0}, {5,1}, {7,2}, {9,3}, {10,4}, {11,4}, {20,5}, {21,6}, {18,7}, {22,8}, {15,9}}},
	{'U',{{1,0}, {5,1}, {7,2}, {9,3}, {10,4}, {11,4}, {18,5}, {19,6}, {17,7}}}
};

static map<char, map<int, vector<int>>>FIRSTOPParamsList={
	{
		'A',{{10,{2,1,0,10}},{11,{10,1,0,11}}}
	},
	{
		'C',{{10,{2,1,0,10}},{11,{10,1,0,11}}}
	},
	{
		'G',{{10,{2,1,0,10}},{11,{10,1,0,11}}}
	},
	{
		'U',{{10,{2,1,0,10}},{11,{10,1,0,11}}}
	},
	{
		'X',{{10,{2,1,0,10}},{11,{10,1,0,11}}}
	}
};

static map<char, map<int, int>>ClsTorParams={
	{'A',{{0,5},{1,6},{2,7},{3,8},{4,9},{5,12}}},
	{'C',{{0,5},{1,6},{2,7},{3,8},{4,9},{5,12}}},
	{'G',{{0,5},{1,6},{2,7},{3,8},{4,9},{5,12}}},
	{'U',{{0,5},{1,6},{2,7},{3,8},{4,9},{5,12}}},
	{'X',{{0,5},{1,6},{2,7},{3,8},{4,9},{5,12}}}
};

static map<char, map<int, vector<int>>> PairTorParam={
	{'A',{{0,{5,7}},{1,{8,9}},{2,{8,12}}}},
	{'C',{{0,{5,7}},{1,{8,9}},{2,{8,12}}}},
	{'G',{{0,{5,7}},{1,{8,9}},{2,{8,12}}}},
	{'U',{{0,{5,7}},{1,{8,9}},{2,{8,12}}}}
};

static map<char, vector<int>> PairTorParamSugar={
	{'A',{5,6,7,8,9,12}},
	{'C',{5,6,7,8,9,12}},
	{'G',{5,6,7,8,9,12}},
	{'U',{5,6,7,8,9,12}}
};

static map<int, vector<int>> PairBBTorParam={{0,{0,1}},{1,{1,2}},{2,{1,4}},{3,{2,4}},{4,{3,4}}};

static vector<vector<int>> TorsionClsList={
	{5,0,1,2,3,4,5}, {5,6,7,9,8,12}
};

static map<char, vector<int>>SNTTEMPLATEATOMS={
	{'A',{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21}},
	{'C',{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19}},
	{'G',{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}},
	{'U',{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19}}
};

static map<char, vector<int>>BPTEMPLATEATOMS={
	{'A',{8,12,13,14}},
	{'C',{8,12,13,14}},
	{'G',{8,12,13,14}},
	{'U',{8,12,13,14}}
};

static map<char, map<int, vector<int>>>ParamsList={
	{
		'A',{
			{0,{3,4,5,0}},
			{1,{4,5,0,1}},
			{2,{5,0,1,2}},
			{3,{0,1,2,3}},
			{4,{1,2,3,4}},
			{5,{2,3,4,5}},
			{6,{5,3,4,6}},
			{7,{4,6,3,7}},
			{8,{3,7,6,8}},
			{9,{7,8,6,9}},
			{10,{5,1,0,10}},
			{11,{10,5,0,11}},
			{12,{3,7,8,12}},
			{13,{7,8,12,13}},
			{14,{13,8,12,14}}
		}
	},
	{
		'C',{
			{0,{3,4,5,0}},
			{1,{4,5,0,1}},
			{2,{5,0,1,2}},
			{3,{0,1,2,3}},
			{4,{1,2,3,4}},
			{5,{2,3,4,5}},
			{6,{5,3,4,6}},
			{7,{4,6,3,7}},
			{8,{3,7,6,8}},
			{9,{7,8,6,9}},
			{10,{5,1,0,10}},
			{11,{10,5,0,11}},
			{12,{3,7,8,12}},
			{13,{7,8,12,13}},
			{14,{13,8,12,14}}
		}
	},
	{
		'G',{
			{0,{3,4,5,0}},
			{1,{4,5,0,1}},
			{2,{5,0,1,2}},
			{3,{0,1,2,3}},
			{4,{1,2,3,4}},
			{5,{2,3,4,5}},
			{6,{5,3,4,6}},
			{7,{4,6,3,7}},
			{8,{3,7,6,8}},
			{9,{7,8,6,9}},
			{10,{5,1,0,10}},
			{11,{10,5,0,11}},
			{12,{3,7,8,12}},
			{13,{7,8,12,13}},
			{14,{13,8,12,14}}
		}
	},
	{
		'U',{
			{0,{3,4,5,0}},
			{1,{4,5,0,1}},
			{2,{5,0,1,2}},
			{3,{0,1,2,3}},
			{4,{1,2,3,4}},
			{5,{2,3,4,5}},
			{6,{5,3,4,6}},
			{7,{4,6,3,7}},
			{8,{3,7,6,8}},
			{9,{7,8,6,9}},
			{10,{5,1,0,10}},
			{11,{10,5,0,11}},
			{12,{3,7,8,12}},
			{13,{7,8,12,13}},
			{14,{13,8,12,14}}
		}
	}
};

static map<char, map<int, vector<int>>>AddLengthParamsList={
	{'A',{{0,{7,8}}}},
	{'C',{{0,{7,8}}}},
	{'G',{{0,{7,8}}}},
	{'U',{{0,{7,8}}}}
};

static map<char, map<int, vector<int>>>AddAngleParamsList={
	{
		'A',{
			{0,{5,0,10}},
			{1,{10,0,11}},
			{2,{1,0,11}},
			{3,{4,6,9}},
			{4,{2,3,7}},
			{5,{5,4,6}},
			{6,{4,6,8}},
			{7,{6,8,12}}
		}
	},
	{
		'C',{
			{0,{5,0,10}},
			{1,{10,0,11}},
			{2,{1,0,11}},
			{3,{4,6,9}},
			{4,{2,3,7}},
			{5,{5,4,6}},
			{6,{4,6,8}},
			{7,{6,8,12}}
		}
	},
	{
		'G',{
			{0,{5,0,10}},
			{1,{10,0,11}},
			{2,{1,0,11}},
			{3,{4,6,9}},
			{4,{2,3,7}},
			{5,{5,4,6}},
			{6,{4,6,8}},
			{7,{6,8,12}}
		}
	},
	{
		'U',{
			{0,{5,0,10}},
			{1,{10,0,11}},
			{2,{1,0,11}},
			{3,{4,6,9}},
			{4,{2,3,7}},
			{5,{5,4,6}},
			{6,{4,6,8}},
			{7,{6,8,12}}
		}
	}
};

static map<char, map<int, vector<int>>>HydParamsList={
	{
		'A',{
			{0,{3,1,2,22}},
			{1,{3,1,2,23}},
			{2,{2,4,3,24}},
			{3,{5,3,4,25}},
			{4,{8,4,6,26}},
			{5,{4,6,9,27}},
			{6,{6,7,8,28}},
			{7,{8,12,14,29}},
			{8,{19,17,20,30}},
			{9,{19,17,20,31}},
			{10,{13,21,16,32}}
		}
	},
	{
		'C',{
			{0,{3,1,2,20}},
			{1,{3,1,2,21}},
			{2,{2,4,3,22}},
			{3,{5,3,4,23}},
			{4,{8,4,6,24}},
			{5,{4,6,9,25}},
			{6,{6,7,8,26}},
			{7,{13,12,14,27}},
			{8,{12,14,16,28}},
			{9,{16,15,19,29}},
			{10,{16,15,19,30}}
		}
	},
	{
		'G',{
			{0,{3,1,2,23}},
			{1,{3,1,2,24}},
			{2,{2,4,3,25}},
			{3,{5,3,4,26}},
			{4,{8,4,6,27}},
			{5,{4,6,9,28}},
			{6,{6,7,8,29}},
			{7,{8,12,14,30}},
			{8,{20,16,22,31}},
			{9,{20,16,22,32}},
			{10,{20,16,15,33}}
		}
	},
	{
		'U',{
			{0,{3,1,2,20}},
			{1,{3,1,2,21}},
			{2,{2,4,3,22}},
			{3,{5,3,4,23}},
			{4,{8,4,6,24}},
			{5,{4,6,9,25}},
			{6,{6,7,8,26}},
			{7,{8,12,14,27}},
			{8,{12,14,16,28}},
			{9,{12,13,17,29}}
		}
	}
};

static map<char, map<int, int>> HYDATOMID2HYDPARAMID={
	{
		'A',{{22,0},{23,1},{24,2},{25,3},{26,4},{27,5},{28,6},{29,7},{30,8},{31,9},{32,10}}
	},
	{
		'C',{{20,0},{21,1},{22,2},{23,3},{24,4},{25,5},{26,6},{27,7},{28,8},{29,9},{30,10}}
	},
	{
		'G',{{23,0},{24,1},{25,2},{26,3},{27,4},{28,5},{29,6},{30,7},{31,8},{32,9},{33,10}}
	},
	{
		'U',{{20,0},{21,1},{22,2},{23,3},{24,4},{25,5},{26,6},{27,7},{28,8},{29,9}}
	}
};

static vector<int> SUGARRINGPARAMS={6,7,8,9,12};

static map<char, vector<int>> BASESATOMS={
	{'A',{13,12,14}},
	{'C',{13,12,14}},
	{'G',{13,12,14}},
	{'U',{13,12,14}}
};

static map<char, vector<int>> BASEATOM={
	{'A',{15,16,17,18,19,20,21}},
	{'C',{15,16,17,18,19}},
	{'G',{15,16,17,18,19,20,21,22}},
	{'U',{15,16,17,18,19}}
};

static map<char, vector<int>> BASEFULLATOM={
	{'A',{12,13,14,15,16,17,18,19,20,21}},
	{'C',{12,13,14,15,16,17,18,19}},
	{'G',{12,13,14,15,16,17,18,19,20,21,22}},
	{'U',{12,13,14,15,16,17,18,19}}
};

static map<char, vector<int>> PHOPHGLYCOSIDEATOMS={
	//id of phosphoric acid, id of nitrogen atom, id of C1', id of torsion angle of delta
	{'A',{0,12,7,5}},
	{'C',{0,12,7,5}},
	{'G',{0,12,7,5}},
	{'U',{0,12,7,5}}
};

static map<char, vector<int>> AccAtomParamNum={
	{'A',{6,12,15,22,33}},
	{'C',{6,12,15,20,31}},
	{'G',{6,12,15,23,34}},
	{'U',{6,12,15,20,30}}
};

static map<char, vector<int>> BaseHydParamIdx={
	{'A',{7,8,9,10}},
	{'C',{7,8,9,10}},
	{'G',{7,8,9,10}},
	{'U',{7,8,9}}
};


static map<char, vector<vector<int>>> CombTorParamIdxList={
	{'A',{{0,1},{1,2},{2,4},{5,7},{5,6},{7,8},{8,9},{8,12}}},
	{'C',{{0,1},{1,2},{2,4},{5,7},{5,6},{7,8},{8,9},{8,12}}},
	{'G',{{0,1},{1,2},{2,4},{5,7},{5,6},{7,8},{8,9},{8,12}}},
	{'U',{{0,1},{1,2},{2,4},{5,7},{5,6},{7,8},{8,9},{8,12}}}
};


static map<char, map<char, vector<vector<int>>>>BASEAXISBP={
	{
		'A',{
			{'A',{{12,13,14,12,13,14}}},
			{'C',{{12,13,14,12,13,14}}},
			{'G',{{12,13,14,12,13,14}}},
			{'U',{{12,13,14,12,13,14}}}
		}
	},
	{
		'C',{
			{'A',{{12,13,14,12,13,14}}},
			{'C',{{12,13,14,12,13,14}}},
			{'G',{{12,13,14,12,13,14}}},
			{'U',{{12,13,14,12,13,14}}}
		}
	},
	{
		'G',{
			{'A',{{12,13,14,12,13,14}}},
			{'C',{{12,13,14,12,13,14}}},
			{'G',{{12,13,14,12,13,14}}},
			{'U',{{12,13,14,12,13,14}}}
		}
	},
	{
		'U',{
			{'A',{{12,13,14,12,13,14}}},
			{'C',{{12,13,14,12,13,14}}},
			{'G',{{12,13,14,12,13,14}}},
			{'U',{{12,13,14,12,13,14}}}
		}
	}
};

static map<char, map<char, vector<vector<int>>>>BASEAXISST={
	{
		'A',{
			{'A',{{12,13,14,12,13,14}}},
			{'C',{{12,13,14,12,13,14}}},
			{'G',{{12,13,14,12,13,14}}},
			{'U',{{12,13,14,12,13,14}}}
		}
	},
	{
		'C',{
			{'A',{{12,13,14,12,13,14}}},
			{'C',{{12,13,14,12,13,14}}},
			{'G',{{12,13,14,12,13,14}}},
			{'U',{{12,13,14,12,13,14}}}
		}
	},
	{
		'G',{
			{'A',{{12,13,14,12,13,14}}},
			{'C',{{12,13,14,12,13,14}}},
			{'G',{{12,13,14,12,13,14}}},
			{'U',{{12,13,14,12,13,14}}}
		}
	},
	{
		'U',{
			{'A',{{12,13,14,12,13,14}}},
			{'C',{{12,13,14,12,13,14}}},
			{'G',{{12,13,14,12,13,14}}},
			{'U',{{12,13,14,12,13,14}}}
		}
	}
};

static map<char, map<char, vector<vector<int>>>>BasePairTorList={
	{
		'A',{
			{'A',{{8,12,12,8}}},
			{'C',{{8,12,12,8}}},
			{'G',{{8,12,12,8}}},
			{'U',{{8,12,12,8}}}
		}
	},
	{
		'C',{
			{'A',{{8,12,12,8}}},
			{'C',{{8,12,12,8}}},
			{'G',{{8,12,12,8}}},
			{'U',{{8,12,12,8}}}
		}
	},
	{
		'G',{
			{'A',{{8,12,12,8}}},
			{'C',{{8,12,12,8}}},
			{'G',{{8,12,12,8}}},
			{'U',{{8,12,12,8}}}
		}
	},
	{
		'U',{
			{'A',{{8,12,12,8}}},
			{'C',{{8,12,12,8}}},
			{'G',{{8,12,12,8}}},
			{'U',{{8,12,12,8}}}
		}
	}
};

static map<char, map<char, vector<vector<int>>>>BasePairLenList={
	{
		'A',{
			{'A',{{0,0},{1,1},{2,2},{3,3},{4,4},{5,5},{8,8},{12,12}}},
			{'C',{{0,0},{1,1},{2,2},{3,3},{4,4},{5,5},{8,8},{12,12}}},
			{'G',{{0,0},{1,1},{2,2},{3,3},{4,4},{5,5},{8,8},{12,12}}},
			{'U',{{0,0},{1,1},{2,2},{3,3},{4,4},{5,5},{8,8},{12,12}}}
		}
	},
	{
		'C',{
			{'A',{{0,0},{1,1},{2,2},{3,3},{4,4},{5,5},{8,8},{12,12}}},
			{'C',{{0,0},{1,1},{2,2},{3,3},{4,4},{5,5},{8,8},{12,12}}},
			{'G',{{0,0},{1,1},{2,2},{3,3},{4,4},{5,5},{8,8},{12,12}}},
			{'U',{{0,0},{1,1},{2,2},{3,3},{4,4},{5,5},{8,8},{12,12}}}
		}
	},
	{
		'G',{
			{'A',{{0,0},{1,1},{2,2},{3,3},{4,4},{5,5},{8,8},{12,12}}},
			{'C',{{0,0},{1,1},{2,2},{3,3},{4,4},{5,5},{8,8},{12,12}}},
			{'G',{{0,0},{1,1},{2,2},{3,3},{4,4},{5,5},{8,8},{12,12}}},
			{'U',{{0,0},{1,1},{2,2},{3,3},{4,4},{5,5},{8,8},{12,12}}}
		}
	},
	{
		'U',{
			{'A',{{0,0},{1,1},{2,2},{3,3},{4,4},{5,5},{8,8},{12,12}}},
			{'C',{{0,0},{1,1},{2,2},{3,3},{4,4},{5,5},{8,8},{12,12}}},
			{'G',{{0,0},{1,1},{2,2},{3,3},{4,4},{5,5},{8,8},{12,12}}},
			{'U',{{0,0},{1,1},{2,2},{3,3},{4,4},{5,5},{8,8},{12,12}}}
		}
	}
};

static map<char, map<char, vector<vector<int>>>>BaseStackingTorList={
	{
		'A',{
			{'A',{{8,12,12,8}}},
			{'C',{{8,12,12,8}}},
			{'G',{{8,12,12,8}}},
			{'U',{{8,12,12,8}}}
		}
	},
	{
		'C',{
			{'A',{{8,12,12,8}}},
			{'C',{{8,12,12,8}}},
			{'G',{{8,12,12,8}}},
			{'U',{{8,12,12,8}}}
		}
	},
	{
		'G',{
			{'A',{{8,12,12,8}}},
			{'C',{{8,12,12,8}}},
			{'G',{{8,12,12,8}}},
			{'U',{{8,12,12,8}}}
		}
	},
	{
		'U',{
			{'A',{{8,12,12,8}}},
			{'C',{{8,12,12,8}}},
			{'G',{{8,12,12,8}}},
			{'U',{{8,12,12,8}}}
		}
	}
};

static map<char, map<char, vector<vector<int>>>>BaseStackingLenList={
	{
		'A',{
			{'A',{{12,12}}},
			{'C',{{12,12}}},
			{'G',{{12,12}}},
			{'U',{{12,12}}}
		}
	},
	{
		'C',{
			{'A',{{12,12}}},
			{'C',{{12,12}}},
			{'G',{{12,12}}},
			{'U',{{12,12}}}
		}
	},
	{
		'G',{
			{'A',{{12,12}}},
			{'C',{{12,12}}},
			{'G',{{12,12}}},
			{'U',{{12,12}}}
		}
	},
	{
		'U',{
			{'A',{{12,12}}},
			{'C',{{12,12}}},
			{'G',{{12,12}}},
			{'U',{{12,12}}}
		}
	}
};

static map<char, map<char, map<int, map<string, int>>>> HbondParamInfo={
	{'A',{
		{'A',{
			{0,{{"N6-N3",0}}},
			{4,{{"N6-N7",0}}},
			{1,{{"N6-N1",0}}},
			{2,{{"N6-N1",0}}},
			{3,{{"N6-N1",0}, {"N6-N7",1}}},
			{5,{{"N6-N7",0}}}
			}
		},
		{'C',{
			{2,{{"N6-O2",0}}},
			{3,{{"N6-N3",0},{"N4-N7",1}}},
			{1,{{"N6-N3",0}}},
			{0,{{"N4-N1",0}}},
			{5,{{"N4-N3",0}}}
			}
		},
		{'G',{
			{3,{{"N2-N7",0},{"N6-N3",1}}},
			{0,{{"N2-N7",0}}},
			{2,{{"N6-O6",0},{"N1-N1",1}}}
			}
		},
		{'U',{
			{0,{{"N3-N1",0},{"N6-O4",1}}}
			}
		}
		}
	},
	{'C',{
		{'A',{
			{2,{{"N6-O2",0}}},
			{3,{{"N6-N3",0},{"N4-N7",1}}},
			{1,{{"N6-N3",0}}},
			{0,{{"N4-N1",0}}},
			{5,{{"N4-N3",0}}}
			}
		},
		{'C',{
			{0,{{"N4-O2",0}}},
			{1,{{"N4-N3",0}}}
			}
		},
		{'G',{
			{0,{{"N1-N3",0},{"N2-O2",1},{"N4-O6",2}}}
			}
		},
		{'U',{
			{3,{{"N4-O2",0}}},
			{1,{{"N4-O4",0}}},
			{0,{{"N3-N3",0},{"N4-O4",1}}},
			{2,{{"N3-O2",0}}},
			}
		}
		}
	},
	{'G',{
		{'A',{
			{3,{{"N2-N7",0},{"N6-N3",1}}},
			{0,{{"N2-N7",0}}},
			{2,{{"N6-O6",0},{"N1-N1",1}}}
			}
		},
		{'C',{
			{0,{{"N1-N3",0},{"N2-O2",1},{"N4-O6",2}}}
			}
		},
		{'G',{
			{0,{{"N2-O6",0}}},
			{2,{{"N1-O6",0}}},
			{1,{{"N2-N7",0},{"N1-O6",1}}},
			{4,{{"N2-N7",0}}},
			{3,{{"N2-N3",0}}},
			{7,{{"N1-O6",0}}}
			}
		},
		{'U',{
			{0,{{"N1-O2",0}, {"N3-O6",1}}}
			}
		}
		}
	},
	{'U',{
		{'A',{
			{0,{{"N3-N1",0},{"N6-O4",1}}}
			}
		},
		{'C',{
			{3,{{"N4-O2",0}}},
			{1,{{"N4-O4",0}}},
			{0,{{"N3-N3",0},{"N4-O4",1}}},
			{2,{{"N3-O2",0}}},
			}
		},
		{'G',{
			{0,{{"N1-O2",0}, {"N3-O6",1}}}
			}
		},
		{'U',{
			{2,{{"N3-O4",0}, {"N3-O2",1}}},
			{0,{{"N3-O4",0}}},
			}
		}
		}
	}
};

static map<char, map<char, map<int, map<string, vector<int>>>>> BasePairADorder={
	{'A',{
		{'A',{
			{0,{{"N6-N3",{-1,-1}}}},
			{4,{{"N6-N7",{-1,-1}}, {"N6-N7",{-1,-1}}}},
			{1,{{"N6-N1",{-1,-1}}}},
			{2,{{"N6-N1",{-1,-1}}, {"N6-N1",{-1,-1}}}},
			{3,{{"N6-N1",{-1,-1}}, {"N6-N7",{-1,-1}}}},
			{5,{{"N6-N7",{-1,-1}}}},
			}
		},
		{'C',{
			{2,{{"N6-O2",{0,1}}}},
			{3,{{"N6-N3",{0,1}},{"N4-N7",{1,0}}}},
			{1,{{"N6-N3",{0,1}}}},
			{0,{{"N4-N1",{1,0}}}},
			{5,{{"N4-N3",{1,0}}}}
			}
		},
		{'G',{
			{3,{{"N2-N7",{1,0}},{"N6-N3",{0,1}}}},
			{0,{{"N2-N7",{1,0}}}},
			{2,{{"N6-O6",{0,1}},{"N1-N1",{1,0}}}}
			}
		},
		{'U',{
			{0,{{"N3-N1",{1,0}},{"N6-O4",{0,1}}}}
			}
		}
		}
	},
	{'C',{
		{'A',{
			{2,{{"N6-O2",{1,0}}}},
			{3,{{"N6-N3",{1,0}},{"N4-N7",{0,1}}}},
			{1,{{"N6-N3",{1,0}}}},
			{0,{{"N4-N1",{0,1}}}},
			{5,{{"N4-N3",{0,1}}}}
			}
		},
		{'C',{
			{0,{{"N4-O2",{-1,-1}}}},
			{1,{{"N4-N3",{-1,-1}}}}
			}
		},
		{'G',{
			{0,{{"N1-N3",{1,0}},{"N2-O2",{1,0}},{"N4-O6",{0,1}}}}
			}
		},
		{'U',{
			{3,{{"N4-O2",{0,1}}}},
			{1,{{"N4-O4",{0,1}}}},
			{0,{{"N3-N3",{1,0}},{"N4-O4",{0,1}}}},
			{2,{{"N3-O2",{1,0}}}}
			}
		}
		}
	},
	{'G',{
		{'A',{
			{3,{{"N2-N7",{0,1}},{"N6-N3",{1,0}}}},
			{0,{{"N2-N7",{0,1}}}},
			{2,{{"N6-O6",{1,0}},{"N1-N1",{0,1}}}}
			}
		},
		{'C',{
			{0,{{"N1-N3",{0,1}},{"N2-O2",{0,1}},{"N4-O6",{1,0}}}}
			}
		},
		{'G',{
			{0,{{"N2-O6",{-1,-1}}}},
			{2,{{"N1-O6",{-1,-1}}}},
			{1,{{"N2-N7",{-1,-1}},{"N1-O6",{-1,-1}}}},
			{4,{{"N2-N7",{-1,-1}}}},
			{3,{{"N2-N3",{-1,-1}},{"N2-N3",{-1,-1}}}},
			{7,{{"N1-O6",{-1,-1}},{"N1-O6",{-1,-1}}}}
			}
		},
		{'U',{
			{0,{{"N1-O2",{0,1}}, {"N3-O6",{1,0}}}}
			}
		}
		}
	},
	{'U',{
		{'A',{
			{0,{{"N3-N1",{0,1}},{"N6-O4",{1,0}}}}
			}
		},
		{'C',{
			{3,{{"N4-O2",{1,0}}}},
			{1,{{"N4-O4",{1,0}}}},
			{0,{{"N3-N3",{0,1}},{"N4-O4",{1,0}}}},
			{2,{{"N3-O2",{0,1}}}}
			}
		},
		{'G',{
			{0,{{"N1-O2",{1,0}}, {"N3-O6",{0,1}}}}
			}
		},
		{'U',{
			{2,{{"N3-O4",{-1,-1}}, {"N3-O2",{-1,-1}}}},
			{0,{{"N3-O4",{-1,-1}}}}
			}
		}
		}
	}
};

static map<char, map<char, map<int, map<string, vector<vector<int>>>>>> BasePairHBondInfo={
	{'A',{
		{'A',{
			{0,{{"N6-N3",{{20,30,21},{20,31,21}}}}},
			{4,{{"N6-N7",{{20,30,18},{20,31,18}}}, {"N6-N7",{{20,30,18},{20,31,18}}}}},
			{1,{{"N6-N1",{{20,30,15},{20,31,15}}}}},
			{2,{{"N6-N1",{{20,30,15},{20,31,15}}}, {"N6-N1",{{20,30,15},{20,31,15}}}}},
			{3,{{"N6-N1",{{20,30,15},{20,31,15}}}, {"N6-N7",{{20,30,18},{20,31,18}}}}},
			{5,{{"N6-N7",{{20,30,18},{20,31,18}}}}}
			}
		},
		{'C',{
			{2,{{"N6-O2",{{20,30,18},{20,31,18}}}}},
			{3,{{"N6-N3",{{20,30,17},{20,31,17}}},{"N4-N7",{{19,29,18},{19,30,18}}}}},
			{1,{{"N6-N3",{{20,30,17},{20,31,17}}}}},
			{0,{{"N4-N1",{{19,29,15},{19,30,15}}}}},
			{5,{{"N4-N3",{{19,29,21},{19,30,21}}}}}
			}
		},
		{'G',{
			{3,{{"N2-N7",{{22,31,18},{22,32,18}}},{"N6-N3",{{20,30,20},{20,31,20}}}}},
			{0,{{"N2-N7",{{22,31,18},{22,32,18}}}}},
			{2,{{"N6-O6",{{20,30,21},{20,31,21}}},{"N1-N1",{{15,33,15}}}}}
			}
		},
		{'U',{
			{0,{{"N3-N1",{{17,29,15}}},{"N6-O4",{{20,30,19},{20,31,19}}}}}
			}
		}
		}
	},
	{'C',{
		{'A',{
			{2,{{"N6-O2",{{20,30,18},{20,31,18}}}}},
			{3,{{"N6-N3",{{20,30,17},{20,31,17}}},{"N4-N7",{{19,29,18},{19,30,18}}}}},
			{1,{{"N6-N3",{{20,30,17},{20,31,17}}}}},
			{0,{{"N4-N1",{{19,29,15},{19,30,15}}}}},
			{5,{{"N4-N3",{{19,29,21},{19,30,21}}}}}
			}
		},
		{'C',{
			{0,{{"N4-O2",{{19,29,18},{19,30,18}}}}},
			{1,{{"N4-N3",{{19,29,17},{19,30,17}}}}}
			}
		},
		{'G',{
			{0,{{"N1-N3",{{15,33,17}}},{"N2-O2",{{22,31,18},{22,32,18}}},{"N4-O6",{{19,29,21},{19,30,21}}}}}
			}
		},
		{'U',{
			{3,{{"N4-O2",{{19,29,18},{19,30,18}}}}},
			{1,{{"N4-O4",{{19,29,19},{19,30,19}}}}},
			{0,{{"N3-N3",{{17,29,17}}},{"N4-O4",{{19,29,19},{19,30,19}}}}},
			{2,{{"N3-O2",{{17,29,18}}}}}
			}
		}
		}
	},
	{'G',{
		{'A',{
			{3,{{"N2-N7",{{22,31,18},{22,32,18}}},{"N6-N3",{{20,30,20},{20,31,20}}}}},
			{0,{{"N2-N7",{{22,31,18},{22,32,18}}}}},
			{2,{{"N6-O6",{{20,30,21},{20,31,21}}},{"N1-N1",{{15,33,15}}}}}
			}
		},
		{'C',{
			{0,{{"N1-N3",{{15,33,17}}},{"N2-O2",{{22,31,18},{22,32,18}}},{"N4-O6",{{19,29,21},{19,30,21}}}}}
			}
		},
		{'G',{
			{0,{{"N2-O6",{{22,31,21},{22,32,21}}}}},
			{2,{{"N1-O6",{{15,33,21}}}}},
			{1,{{"N2-N7",{{22,31,18},{22,32,18}}},{"N1-O6",{{15,33,21}}}}},
			{4,{{"N2-N7",{{22,31,18},{22,32,18}}}}},
			{3,{{"N2-N3",{{22,31,20},{22,32,20}}},{"N2-N3",{{22,31,20},{22,32,20}}}}},
			{7,{{"N1-O6",{{15,33,21}}},{"N1-O6",{{15,33,21}}}}}
			}
		},
		{'U',{
			{0,{{"N1-O2",{{15,33,18}}}, {"N3-O6",{{17,29,21}}}}}
			}
		}
		}
	},
	{'U',{
		{'A',{
			{0,{{"N3-N1",{{17,29,15}}},{"N6-O4",{{20,30,19},{20,31,19}}}}}
			}
		},
		{'C',{
			{3,{{"N4-O2",{{19,29,18},{19,30,18}}}}},
			{1,{{"N4-O4",{{19,29,19},{19,30,19}}}}},
			{0,{{"N3-N3",{{17,29,17}}},{"N4-O4",{{19,29,19},{19,30,19}}}}},
			{2,{{"N3-O2",{{17,29,18}}}}}
			}
		},
		{'G',{
			{0,{{"N1-O2",{{15,33,18}}}, {"N3-O6",{{17,29,21}}}}}
			}
		},
		{'U',{
			{2,{{"N3-O4",{{17,29,19}}}, {"N3-O2",{{17,29,18}}}}},
			{0,{{"N3-O4",{{17,29,19}}}}}
			}
		}
		}
	}
};

static map<char, map<char, map<int, map<string, int>>>> HbondParamIdx={
	{'A',{
		{'A',{
			{0,{{"N6-N3",0}}},
			{4,{{"N6-N7",0}, {"N6-N7",1}}},
			{1,{{"N6-N1",0}}},
			{2,{{"N6-N1",0}, {"N6-N1",1}}},
			{3,{{"N6-N1",0}, {"N6-N7",1}}},
			{5,{{"N6-N7",0}}}
			}
		},
		{'C',{
			{2,{{"N6-O2",0}}},
			{3,{{"N6-N3",0},{"N4-N7",1}}},
			{1,{{"N6-N3",0}}},
			{0,{{"N4-N1",0}}},
			{5,{{"N4-N3",0}}}
			}
		},
		{'G',{
			{3,{{"N2-N7",0},{"N6-N3",1}}},
			{0,{{"N2-N7",0}}},
			{2,{{"N6-O6",0},{"N1-N1",1}}}
			}
		},
		{'U',{
			{0,{{"N3-N1",0},{"N6-O4",1}}}
			}
		}
		}
	},
	{'C',{
		{'A',{
			{2,{{"N6-O2",0}}},
			{3,{{"N6-N3",0},{"N4-N7",1}}},
			{1,{{"N6-N3",0}}},
			{0,{{"N4-N1",0}}},
			{5,{{"N4-N3",0}}}
			}
		},
		{'C',{
			{0,{{"N4-O2",0}}},
			{1,{{"N4-N3",0}}}
			}
		},
		{'G',{
			{0,{{"N1-N3",0},{"N2-O2",1},{"N4-O6",2}}}
			}
		},
		{'U',{
			{3,{{"N4-O2",0}}},
			{1,{{"N4-O4",0}}},
			{0,{{"N3-N3",0},{"N4-O4",1}}},
			{2,{{"N3-O2",0}}},
			}
		}
		}
	},
	{'G',{
		{'A',{
			{3,{{"N2-N7",0},{"N6-N3",1}}},
			{0,{{"N2-N7",0}}},
			{2,{{"N6-O6",0},{"N1-N1",1}}}
			}
		},
		{'C',{
			{0,{{"N1-N3",0},{"N2-O2",1},{"N4-O6",2}}}
			}
		},
		{'G',{
			{0,{{"N2-O6",0}}},
			{2,{{"N1-O6",0}}},
			{1,{{"N2-N7",0},{"N1-O6",1}}},
			{4,{{"N2-N7",0}}},
			{3,{{"N2-N3",0},{"N2-N3",1}}},
			{7,{{"N1-O6",0},{"N1-O6",1}}}
			}
		},
		{'U',{
			{0,{{"N1-O2",0}, {"N3-O6",1}}}
			}
		}
		}
	},
	{'U',{
		{'A',{
			{0,{{"N3-N1",0},{"N6-O4",1}}}
			}
		},
		{'C',{
			{3,{{"N4-O2",0}}},
			{1,{{"N4-O4",0}}},
			{0,{{"N3-N3",0},{"N4-O4",0}}},
			{2,{{"N3-O2",0}}},
			}
		},
		{'G',{
			{0,{{"N1-O2",0}, {"N3-O6",1}}}
			}
		},
		{'U',{
			{2,{{"N3-O4",0}, {"N3-O2",1}}},
			{0,{{"N3-O4",0}}},
			}
		}
		}
	}
};

static map<char, map<char,int>> BaseTypeList={
	{'A',{{'A',0},{'C',1},{'G',2},{'U',3}}},
	{'C',{{'A',1},{'C',4},{'G',5},{'U',6}}},
	{'G',{{'A',2},{'C',5},{'G',7},{'U',8}}},
	{'U',{{'A',3},{'C',6},{'G',8},{'U',9}}}
};

static vector<int> MODTARGET={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
static vector<int> BLENGTHTARGET={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
static vector<int> SUBBLENGTHTARGET={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
static vector<int> ADDBLENGTHTARGET={0};
static vector<int> INNERANGLETARGET={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
static vector<int> ADDINNERANGLETARGET={0,1,2,3,4,5,6,7};
static vector<int> TORSIONANGLETARGET={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
static vector<int> BASECOORDTORSIONANGLETARGET={0,1,2,3,4,5,6,7,8,9,10,11,12,13};
static vector<int> BBATOMIDX={0,1,2,3,4,5,6,7,8,9,10,11,12,13};

typedef struct Point3s
{
	float x,y,z;
    Point3s() : x(0.0f), y(0.0f), z(0.0f) {}
    Point3s(const Point3s& other) {
        x = other.x;
        y = other.y;
        z = other.z;
    }
}Point3s;

typedef struct Point3d
{
	double x,y,z;
    Point3d() : x(0.0), y(0.0), z(0.0) {}
    Point3d(const Point3d& other) {
        x = other.x;
        y = other.y;
        z = other.z;
    }
}Point3d;

typedef struct Atom	 // Structure For An Object
{
	int        simpletype;
	int        modelid;					
	string     seqno;           //原子の通しid
	string     detailtype;  // 原子名
	int        detail;         // 原子の種類のid
	string     residue;     // 残基名
	char       residueid; // 残基の種類のid
	char       chainid;        //chainのid
	int        resno; 		//残基の通しid
	float	   x,y,z;		//atomの座標
	char       alt,ins;		//ins:残基の挿入　コード
} Atom;	

typedef struct Segment
{
	int init;//modelの最初のatomのindex
	int term;//modelの最後のatomのindex
}Segment;

typedef struct Prochain
{
	Segment chainseg; //該当chainのatomの開始位置と終了位置の情報
	char chainid;	//該当chainのid
}Prochain;

typedef struct Promodel
{
	Prochain procha[MAX_CHAIN_NUM];//各chainの情報
	int  nchain;	//chain数
}Promodel;

typedef struct Boneinfo
{
	map<int, int> indatom;
	int resind;//残基の通しid
	char resid;//残基名
	int istart,iend;//各残基の最初と最後のatomにおける原子の通しid
}Boneinfo;

typedef struct Cutoffparamtype{
	double lin,lmiddle,lout,rin,rmiddle,rout;
	int ptype;
}Cutoffparamtype;

typedef struct Binnedparamtype{
	double minval;
	double maxval;
	double span;
	int binnum;
	double mean;
	double std;
	double mode;
	double upper, supper;
	double lower, slower;
	int rangetype, srangetype;
	int total;
	int totaldist;
	double totallogdist;
	vector<int> distr;
	vector<double> logdistr;
	vector<double> energy;
	vector<int> binary;
	vector<double> boundary;
	vector<double> center;
	//vector<vector<double>> tmpdist;
}Binnedparamtype;

typedef struct Torsionclsparamtype{
	double mean;
	double std;
	double upper;
	double lower;
	int rangetype;
}Torsionclsparamtype;

typedef struct Hydoparamtype{
	double minval1, minval2;
	double maxval1, maxval2;
	double span1, span2;
	int binnum1, binnum2;
	vector<int> distr;
	vector<double> energy;
	vector<int> binary;
	vector<double> boundary1, boundary2;
	vector<double> center1, center2;
}Hydoparamtype;

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

typedef struct D3Binnedparamtype{
	double minval1, minval2, minval3;
	double maxval1, maxval2, maxval3;
	double span1, span2, span3;
	int binnum1, binnum2, binnum3;
	vector<int> distr;
	vector<double> energy;
	vector<int> binary;
	vector<double> boundary1, boundary2, boundary3;
	vector<double> center1, center2, center3;
}D3Binnedparamtype;

typedef struct Baseparamstype{
	Binnedparamtype tor;
	Binnedparamtype leng;
	map<int, Cutoffparamtype> lengcutoff;
	map<int, Cutoffparamtype> torcutoff;
	vector<map<int, Cutoffparamtype>> hbondcutoff;
}Baseparamstype;

typedef struct Phycontact
{
	double pval[CONTACT_ENERGY_NUM+OXHBOND_ENERGY_NUM];
	int pcw[CONTACT_ENERGY_NUM+OXHBOND_ENERGY_NUM];

	Phycontact() {
		std::fill(pval, pval + CONTACT_ENERGY_NUM+OXHBOND_ENERGY_NUM, 0);
		std::fill(pcw, pcw + CONTACT_ENERGY_NUM+OXHBOND_ENERGY_NUM, 0);
	}

	Phycontact(const Phycontact& other) {
		std::copy(other.pval, other.pval + CONTACT_ENERGY_NUM+OXHBOND_ENERGY_NUM, pval);
        std::copy(other.pcw, other.pcw + CONTACT_ENERGY_NUM+OXHBOND_ENERGY_NUM, pcw);
	}
}Phycontact;

typedef struct Rescontact
{
	double pval[CLASH_ENERGY_NUM+OXHBOND_ENERGY_NUM];//atomのclashのペナルティを格納する変数(0:resAの主鎖とresBの主鎖の間のatomのclash数、1:resAの側鎖とresBの主鎖の間のatomのclash数、2:resAの主鎖とresBの側鎖の間のatomのclash数、3:resAの側鎖とresBの側鎖の間のatomのclash数)
	int pcw[CLASH_ENERGY_NUM+OXHBOND_ENERGY_NUM];//atomのclashの数を格納する変数

	Rescontact() {
        std::fill(pval, pval + CLASH_ENERGY_NUM+OXHBOND_ENERGY_NUM, 0.0);
        std::fill(pcw, pcw + CLASH_ENERGY_NUM+OXHBOND_ENERGY_NUM, 0);
    }

	Rescontact(const Rescontact& other) {
        // 配列のコピー
        std::copy(other.pval, other.pval + CLASH_ENERGY_NUM+OXHBOND_ENERGY_NUM, pval);
        std::copy(other.pcw, other.pcw + CLASH_ENERGY_NUM+OXHBOND_ENERGY_NUM, pcw);
    }
}Rescontact;

typedef struct Basecontact
{
	double pval[BASE_ENERGY_NUM+BASEHYD_ENERGY_NUM];
	int pcw[BASE_ENERGY_NUM+BASEHYD_ENERGY_NUM];
	double tval[BASEPAIRTEMP_ENERGY_NUM];

	Basecontact() {
        std::fill(pval, pval + BASE_ENERGY_NUM+BASEHYD_ENERGY_NUM, 0.0);
        std::fill(pcw, pcw + BASE_ENERGY_NUM+BASEHYD_ENERGY_NUM, 0);
		std::fill(tval, tval + BASEPAIRTEMP_ENERGY_NUM, 0);
    }

	Basecontact(const Basecontact& other) {
        // 配列のコピー
        std::copy(other.pval, other.pval + BASE_ENERGY_NUM+BASEHYD_ENERGY_NUM, pval);
        std::copy(other.pcw, other.pcw + BASE_ENERGY_NUM+BASEHYD_ENERGY_NUM, pcw);
		std::copy(other.tval, other.tval + BASEPAIRTEMP_ENERGY_NUM, tval);
    }
}Basecontact;

typedef struct Stackcontact
{
	double pval[STACKING_ENERGY_NUM];
	int pcw[STACKING_ENERGY_NUM];
	double tval[STACKINGTEMP_ENERGY_NUM];

	Stackcontact() {
        std::fill(pval, pval + STACKING_ENERGY_NUM, 0.0);
        std::fill(pcw, pcw + STACKING_ENERGY_NUM, 0);
		std::fill(tval, tval + STACKINGTEMP_ENERGY_NUM, 0);
    }

	Stackcontact(const Stackcontact& other) {
        // 配列のコピー
        std::copy(other.pval, other.pval + STACKING_ENERGY_NUM, pval);
        std::copy(other.pcw, other.pcw + STACKING_ENERGY_NUM, pcw);
		std::copy(other.tval, other.tval + STACKINGTEMP_ENERGY_NUM, tval);
    }
}Stackcontact;

typedef struct Baseinfo
{
	int baseidx;
	int basetype;
	double dihe;
	double clength;
	double baselength1;
	double baselength2;
	double baselength3;
}Baseinfo;

typedef struct PairLengInfo
{
	double mean, std, upper, lower;
}PairLengInfo;

typedef struct TorClsInfo
{
	double mean, std, upper, lower;
}TorClsInfo;

typedef struct Residueatoms
{
	int residueid;//sequence
	float dih[PARAM_NUM],len[PARAM_NUM],ang[PARAM_NUM],oplen[PARAM_NUM];
	float addlen[ADDLENGTH_NUM], addang[ADDINNERANGLE_NUM];
	float hyddih[HYD_PARAM_NUM],hydlen[HYD_PARAM_NUM],hydang[HYD_PARAM_NUM];
	Point3s ptv[ATOM_RES_NUM];//heavily atomの座標の情報
	
	char aaa;//残基名
	int resind; //残基の種類のid
	int residx;
	int ss;
	Baseinfo base;

	Residueatoms() 
        : residueid(0), 
          aaa('\0'), 
          resind(0), 
          residx(0), 
          ss(0), 
          base() // Baseinfoのデフォルトコンストラクタを呼び出す
    {
        // 配列の初期化
        std::fill(dih, dih + PARAM_NUM, 0.0f);
        std::fill(len, len + PARAM_NUM, 0.0f);
        std::fill(ang, ang + PARAM_NUM, 0.0f);
		std::fill(oplen, oplen + PARAM_NUM, 0.0f);
		std::fill(addlen, addlen + ADDLENGTH_NUM, 0.0f);
		std::fill(addang, addang + ADDINNERANGLE_NUM, 0.0f);
		std::fill(hyddih, hyddih + HYD_PARAM_NUM, 0.0f);
        std::fill(hydlen, hydlen + HYD_PARAM_NUM, 0.0f);
        std::fill(hydang, hydang + HYD_PARAM_NUM, 0.0f);
        std::fill(ptv, ptv + ATOM_RES_NUM, Point3s()); // Point3sのデフォルトコンストラクタを呼び出す
    }

	Residueatoms(const Residueatoms& other) 
        : residueid(other.residueid), 
          aaa(other.aaa), 
          resind(other.resind), 
          residx(other.residx), 
          ss(other.ss), 
          base(other.base)
    {
        std::copy(other.dih, other.dih + PARAM_NUM, dih);
        std::copy(other.len, other.len + PARAM_NUM, len);
        std::copy(other.ang, other.ang + PARAM_NUM, ang);
		std::copy(other.oplen, other.oplen + PARAM_NUM, oplen);
		std::copy(other.addlen, other.addlen + ADDLENGTH_NUM, addlen);
        std::copy(other.addang, other.addang + ADDINNERANGLE_NUM, addang);
		std::copy(other.hyddih, other.hyddih + HYD_PARAM_NUM, hyddih);
        std::copy(other.hydlen, other.hydlen + HYD_PARAM_NUM, hydlen);
        std::copy(other.hydang, other.hydang + HYD_PARAM_NUM, hydang);
        std::copy(other.ptv, other.ptv + ATOM_RES_NUM, ptv);
    }
}Residueatoms;

typedef struct TemplateDistMap
{
	vector<vector<vector<double>>> intradistmap;
	vector<vector<double>> interdistmap;

	// デフォルトコンストラクタ
    TemplateDistMap() {
        // 必要ならば、メンバ変数をデフォルトで初期化する処理を追加
        // 例: intradistmapのサイズを決める場合
        intradistmap.resize(2, std::vector<std::vector<double>>(PAIRLENGTHTARGETMOVNUM, std::vector<double>(PAIRLENGTHTARGETMOVNUM, 0.0)));
        interdistmap.resize(PAIRLENGTHTARGETMOVNUM, std::vector<double>(PAIRLENGTHTARGETMOVNUM, 0.0));
    }

    // コピーコンストラクタ
    TemplateDistMap(const TemplateDistMap& other)
        : intradistmap(other.intradistmap), interdistmap(other.interdistmap) {
        // メンバ変数を他のインスタンスからコピー
    }
}TemplateDistMap;

typedef struct Pairaa
{
	bool ispair;//2つの残基が隣接した残基でない2つの残基である場合はtrue、それ以外はfalse
	vector<double> dist;//2つの残基間の距離
	double dstd;//2つの残基間の距離 (3.9より小さい距離は3.9に変換)
	Point3d pval;
	vector<Point3d> cood;
	int pbin,tnum;//2つの残基の通しid (tnumのidの方がpbinのidより大きい)
	vector<double> threshold;
	double weight;
	double ratio[20];

	Pairaa() 
        : ispair(false), 
          dstd(0.0), 
          pbin(0), 
          tnum(0),
		  weight(0)
    {
        std::fill(ratio, ratio + 20, 0.0); // ratio配列を初期化
    }

    // コピーコンストラクタ
    Pairaa(const Pairaa& other) 
        : ispair(other.ispair), 
          dist(other.dist),
          dstd(other.dstd), 
          pval(other.pval), 
          cood(other.cood),
          pbin(other.pbin), 
          tnum(other.tnum), 
          threshold(other.threshold),
		  weight(other.weight)
    {
        std::copy(other.ratio, other.ratio + 20, ratio); // ratio配列をコピー
    }
}Pairaa;

typedef struct BPBasicInfo
{
	pair<int, int> basepidx;
	int bptype;
	int sstype;
	double sc;
	int hydtype;
	vector<vector<double>> hatomtype;
}BPBasicInfo;

typedef struct STBasicInfo
{
	pair<int, int> basepidx;
	int sttype;
	double sc;
}STBasicInfo;

typedef struct CTBasicInfo
{
	pair<int, int> basepidx;
	double sc;
}CTBasicInfo;

typedef struct Energyterms
{
	int nn[TOTAL_ENERGY_NUM];//wrong pos
	double fn[TOTAL_ENERGY_NUM];//totval
	vector<vector<bool>> flap;//pos concern
	vector<vector<int>> indn;//ene int
	vector<vector<double>> indf;//ene dou

    Energyterms() {
        std::fill(std::begin(nn), std::end(nn), 0);
        std::fill(std::begin(fn), std::end(fn), 0.0);
        flap = vector<vector<bool>>(TOTAL_ENERGY_NUM);
        indn = vector<vector<int>>(TOTAL_ENERGY_NUM);
        indf = vector<vector<double>>(TOTAL_ENERGY_NUM);
    }

    Energyterms(const Energyterms& other) {
        std::copy(std::begin(other.nn), std::end(other.nn), std::begin(nn));
        std::copy(std::begin(other.fn), std::end(other.fn), std::begin(fn));
        flap = other.flap;
        indn = other.indn;
        indf = other.indf;
    }
}Energyterms;

typedef struct Gradlength
{
	double mseadj, msent, msebp; 
	double mseadjprev, msentprev, msebpprev;
	int countadj, countnt, countbp;

	vector<vector<vector<double>>> lengths;
	vector<vector<vector<double>>> diffs;
	vector<vector<vector<Point3s>>> dists;
	vector<vector<vector<bool>>> outlineridx;

	Gradlength() 
        : mseadj(0.0), msent(0.0), msebp(0.0),
          mseadjprev(0.0), msentprev(0.0), msebpprev(0.0),
          countadj(0), countnt(0), countbp(0) {}

	Gradlength(const Gradlength& other) 
        : mseadj(other.mseadj), msent(other.msent), msebp(other.msebp),
          mseadjprev(other.mseadjprev), msentprev(other.msentprev), msebpprev(other.msebpprev),
          countadj(other.countadj), countnt(other.countnt), countbp(other.countbp),
          lengths(other.lengths), diffs(other.diffs), dists(other.dists), outlineridx(other.outlineridx) {}
}Gradlength;

typedef struct Gradtemp
{
	double mse, mseprev;
	int count;
	map<int, vector<vector<double>>> diffs; //initlengths, tmplengths, bptmplengths
	map<int, vector<vector<Point3s>>> dists; //initdists, tmpdists, bptmpdists
	map<int, vector<vector<double>>> lengths;

	map<int, vector<vector<double>>> intradiffs1; 
	map<int, vector<vector<Point3s>>> intradists1; 
	map<int, vector<vector<double>>> intralengths1;
	map<int, vector<vector<bool>>> outlineridx1;

	map<int, vector<vector<double>>> intradiffs2; 
	map<int, vector<vector<Point3s>>> intradists2; 
	map<int, vector<vector<double>>> intralengths2;
	map<int, vector<vector<bool>>> outlineridx2;

	Gradtemp() : mse(0.0), mseprev(0.0) {}

    Gradtemp(const Gradtemp& other) 
        : mse(other.mse), mseprev(other.mseprev), diffs(other.diffs), dists(other.dists), lengths(other.lengths), 
		  intradiffs1(other.intradiffs1), intradists1(other.intradists1), intralengths1(other.intralengths1), outlineridx1(other.outlineridx1), 
		  intradiffs2(other.intradiffs2), intradists2(other.intradists2), intralengths2(other.intralengths2), outlineridx2(other.outlineridx2) {}

}Gradtemp;

typedef struct Gradtemplate
{
	double mse, mseprev;
	int count;
	map<int, vector<Point3s>> dists; //initdists, tmpdists, bptmpdists
	map<int, vector<double>> lengths;

	Gradtemplate() : mse(0.0), mseprev(0.0) {}

    Gradtemplate(const Gradtemplate& other) 
        : mse(other.mse), mseprev(other.mseprev), dists(other.dists), lengths(other.lengths) {}

}Gradtemplate;

typedef struct GradtempInit
{
	double mse, mseprev;
	int count;
	map<int, vector<double>> tmpdiff; //initlengths, tmplengths, bptmplengths
	map<int, vector<Point3s>> tmpdists; //initdists, tmpdists, bptmpdists

	GradtempInit() : mse(0.0), mseprev(0.0) {}

    GradtempInit(const GradtempInit& other) 
        : mse(other.mse), mseprev(other.mseprev), tmpdiff(other.tmpdiff), tmpdists(other.tmpdists) {}

}GradtempInit;

typedef struct Gradangle
{
	double mse, mseprev;
	int count;
	map<int, vector<double>> angles;
	map<int, vector<double>> anglediffs;
	map<int, vector<bool>> angleouterflags;
	map<int, vector<double>> theta;
	map<int, vector<double>> normPQs;
	map<int, vector<double>> normRQs; 
	map<int, vector<double>> dotps;
	map<int, vector<Point3s>> vecPQs; 
	map<int, vector<Point3s>> vecRQs;

	Gradangle() : mse(0.0), mseprev(0.0) {}

    Gradangle(const Gradangle& other) 
        : mse(other.mse), mseprev(other.mseprev), angles(other.angles),
          anglediffs(other.anglediffs), angleouterflags(other.angleouterflags), theta(other.theta),
          normPQs(other.normPQs), normRQs(other.normRQs), dotps(other.dotps),
          vecPQs(other.vecPQs), vecRQs(other.vecRQs) {}
}Gradangle;

typedef struct Gradmain2
{
	double mse, mseprev;
	Gradlength lengthobj;
	Gradtemplate adjtempobj;
	Gradtemp adjtemplengthobj;
	GradtempInit inittempobj;
	Gradtemplate bptempobj;
	Gradangle angleobj;

	Gradmain2() : mse(0.0), mseprev(0.0) {}

    Gradmain2(const Gradmain2& other) 
        : mse(other.mse), mseprev(other.mseprev),
          lengthobj(other.lengthobj), adjtempobj(other.adjtempobj),
          inittempobj(other.inittempobj), bptempobj(other.bptempobj), angleobj(other.angleobj), adjtemplengthobj(other.adjtemplengthobj) {}
}Gradmain2;

double Random(void);
void sgenrand(unsigned long seed);
void SelectStream(int index);
void PlantSeeds(long x);
void PutSeed(long x);






















