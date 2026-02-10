/* Purpose: Refine bond lengths and angles
Ideal bond and pseudobond lengths are taken from the Amber forcefield parameters.
*/

#ifndef COMMONPARAMS_HPP
#define COMMONPARAMS_HPP 1

#include <math.h>
#include <stdbool.h>
#include <vector>
#include <map>
#include <utility>
#include <cstdlib>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

#define PARAM_NUM 23
#define ADDLENGTH_NUM 44
#define ADDINNERANGLE_NUM 20
#define ATOM_RES_NUM 34
#define HYD_PARAM_NUM 11
#define HYD_ADDPARAM_NUM 20
#define HYD_ADDANGLE_NUM 14
#define ATOM_TYPE_NUM 18

#define COULOMB 332.063711
#define EpsilonRNA 20
#define CUTOFF 15
#define CUTOFF_SQR CUTOFF*CUTOFF
#define GLOBALCUTOFF 30
#define GLOBALCUTOFF_SQR GLOBALCUTOFF*GLOBALCUTOFF
#define HBONDThreshold -1
#define HBONDADLengThreshold 10
#define PI_MACRO (4.0 * atan(1.0))
#define RADDEG (PI_MACRO/180.0)
#define DEGRAD (180.0 / PI_MACRO)
#define HBOND_CONST 4.0

#define EPS 1e-7
#define PI (4.0 * atan(1.0))
#define MODULUS    2147483647
#define MULTIPLIER 48271
#define STREAMS    256 
#define DEFAULT    123456789
static long seed[STREAMS] = {DEFAULT};
static int  stream        = 0;

using namespace std;

static map<string, int> heavily_atom_num ={{"A", 22}, {"C", 20}, {"G", 23}, {"U", 20}};

typedef struct HydroBondParameters
{
    double adlength_mean, adlength_std;
	double dhlength_mean, dhlength_std;
	double ahlength_mean, ahlength_std;
	double adhangle_mean, adhangle_std;
	double dhaangle_mean, dhaangle_std;
	double hadangle_mean, hadangle_std;
}HydroBondParameters;

typedef struct BBHydroBondParameters
{
    double adlength_mean, adlength_std;
    double ahlength_mean, ahlength_std;
	double dhaangle_mean, dhaangle_std;
    double haaaangle_mean, haaaangle_std;
    double daaaangle_mean, daaaangle_std;
}BBHydroBondParameters;

static map<string, vector<string>> full_atoms_pdb={
	{"  A", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'", " N9 ", " C8 ", " N7 ", " C5 ", " C6 ", " N6 ", " N1 ", " C2 ", " N3 ", " C4 ", "H5'1", "H5'2", " H4'", " H3'", " H2'", "HO'2", " H1'", "  H8", " H61", " H62", "  H2"}},
	{"  C", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'", " N1 ", " C2 ", " O2 ", " N3 ", " C4 ", " N4 ", " C5 ", " C6 ", "H5'1", "H5'2", " H4'", " H3'", " H2'", "HO'2", " H1'", "  H6", "  H5", " H41", " H42"}},
	{"  G", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'", " N9 ", " C8 ", " N7 ", " C5 ", " C6 ", " O6 ", " N1 ", " C2 ", " N2 ", " N3 ", " C4 ", "H5'1", "H5'2", " H4'", " H3'", " H2'", "HO'2", " H1'", "  H8", " H21", " H22", "  H1"}},
	{"  U", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'", " N1 ", " C2 ", " O2 ", " N3 ", " C4 ", " O4 ", " C5 ", " C6 ", "H5'1", "H5'2", " H4'", " H3'", " H2'", "HO'2", " H1'", "  H6", "  H5", "  H3"}},
};

static map<string, vector<string>> coarse_grained_atoms_pdb={
	{"  A", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'", " N9 ", " C8 ", " N7 ", " C5 ", " C6 ", " N6 ", " N1 ", " C2 ", " N3 ", " C4 "}},
	{"  C", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'", " N1 ", " C2 ", " O2 ", " N3 ", " C4 ", " N4 ", " C5 ", " C6 "}},
	{"  G", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'", " N9 ", " C8 ", " N7 ", " C5 ", " C6 ", " O6 ", " N1 ", " C2 ", " N2 ", " N3 ", " C4 "}},
	{"  U", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'", " N1 ", " C2 ", " O2 ", " N3 ", " C4 ", " O4 ", " C5 ", " C6 "}},
};

static map<string, vector<string>> full_atoms={
    {"A", {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", "N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4", "H5'1", "H5'2", "H4'", "H3'", "H2'", "HO'2", "H1'", "H8", "H61", "H62", "H2"}},
    {"C", {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", "N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6", "H5'1", "H5'2", "H4'", "H3'", "H2'", "HO'2", "H1'", "H6", "H5", "H41", "H42"}},
    {"G", {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", "N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4", "H5'1", "H5'2", "H4'", "H3'", "H2'", "HO'2", "H1'", "H8", "H21", "H22", "H1"}},
    {"U", {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", "N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6", "H5'1", "H5'2", "H4'", "H3'", "H2'", "HO'2", "H1'", "H6", "H5", "H3"}},
};

static map<string, vector<string>> coarse_grained_atoms={
    {"A", {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", "N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"}},
    {"C", {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", "N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"}},
    {"G", {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", "N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4"}},
    {"U", {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", "N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"}},
};

static map<string, vector<string>> hyd_full_atoms_name={
    {"A", {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", "N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4", "H5'1", "H5'2", "H4'", "H3'", "H2'", "HO'2", "H1'", "H8", "H61", "H62", "H2"}},
    {"C", {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", "N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6", "H5'1", "H5'2", "H4'", "H3'", "H2'", "HO'2", "H1'", "H6", "H5", "H41", "H42"}},
    {"G", {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", "N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4", "H5'1", "H5'2", "H4'", "H3'", "H2'", "HO'2", "H1'", "H8", "H21", "H22", "H1"}},
    {"U", {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", "N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6", "H5'1", "H5'2", "H4'", "H3'", "H2'", "HO'2", "H1'", "H6", "H5", "H3"}},
};

static map<string, vector<string>> full_atoms_name={
    {"A", {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", "N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"}},
    {"C", {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", "N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"}},
    {"G", {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", "N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4"}},
    {"U", {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", "N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"}},
};

static map<string, vector<string>> bb_atoms_name={
    {"A", {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"}},
    {"C", {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"}},
    {"G", {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"}},
    {"U", {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"}},
};

static map<string, vector<string>> bb_atoms={
    {"A", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'"}},
    {"C", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'"}},
    {"G", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'"}},
    {"U", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'"}},
};

static map<string, map<string, int>> hydfullatomname2idx={
    {"A", {{"P", 0}, {"OP1",1}, {"OP2",2}, {"O5'",3}, {"C5'",4}, {"C4'",5}, {"O4'",6}, {"C3'",7}, {"O3'",8}, {"C2'",9}, {"O2'",10}, {"C1'",11}, {"N9", 12}, {"C8", 13}, {"N7", 14}, {"C5", 15}, {"C6", 16}, {"N6", 17}, {"N1", 18}, {"C2", 19}, {"N3", 20}, {"C4",21}, {"H5'1",22}, {"H5'2",23}, {"H4'", 24}, {"H3'", 25}, {"H2'", 26}, {"HO'2", 27}, {"H1'", 28}, {"H8", 29}, {"H61", 30}, {"H62", 31}, {"H2", 32}}},
    {"C", {{"P", 0}, {"OP1",1}, {"OP2",2}, {"O5'",3}, {"C5'",4}, {"C4'",5}, {"O4'",6}, {"C3'",7}, {"O3'",8}, {"C2'",9}, {"O2'",10}, {"C1'",11}, {"N1", 12}, {"C2", 13}, {"O2", 14}, {"N3", 15}, {"C4", 16}, {"N4", 17}, {"C5", 18}, {"C6", 19}, {"H5'1", 20}, {"H5'2", 21}, {"H4'", 22}, {"H3'", 23}, {"H2'", 24}, {"HO'2", 25}, {"H1'", 26}, {"H6", 27}, {"H5", 28}, {"H41", 29}, {"H42",30}}},
    {"G", {{"P", 0}, {"OP1",1}, {"OP2",2}, {"O5'",3}, {"C5'",4}, {"C4'",5}, {"O4'",6}, {"C3'",7}, {"O3'",8}, {"C2'",9}, {"O2'",10}, {"C1'",11}, {"N9", 12}, {"C8", 13}, {"N7", 14}, {"C5", 15}, {"C6", 16}, {"O6", 17}, {"N1", 18}, {"C2", 19}, {"N2", 20}, {"N3", 21}, {"C4", 22}, {"H5'1", 23}, {"H5'2", 24}, {"H4'", 25}, {"H3'", 26}, {"H2'", 27}, {"HO'2", 28}, {"H1'", 29}, {"H8", 30}, {"H21", 31}, {"H22", 32}, {"H1", 33}}},
    {"U", {{"P", 0}, {"OP1",1}, {"OP2",2}, {"O5'",3}, {"C5'",4}, {"C4'",5}, {"O4'",6}, {"C3'",7}, {"O3'",8}, {"C2'",9}, {"O2'",10}, {"C1'",11}, {"N1", 12}, {"C2", 13}, {"O2", 14}, {"N3", 15}, {"C4", 16}, {"O4", 17}, {"C5", 18}, {"C6", 19}, {"H5'1", 20}, {"H5'2", 21}, {"H4'", 22}, {"H3'", 23}, {"H2'", 24}, {"HO'2", 25}, {"H1'", 26}, {"H6", 27}, {"H5", 28}, {"H3", 29}}},
};

static vector<string> template_target_atoms = {" P  ", " C4'"}; 
static vector<int> template_target_atoms_idx = {0, 5}; 

typedef struct ExtParams
{
	map<string, vector<vector<double>>> blengthparam;
    map<string, vector<vector<double>>> blengthstdval;
	map<string, vector<vector<double>>> adjblengthparam;
    map<string, vector<vector<double>>> adjblengthstdval;

	map<string, vector<vector<double>>> ntangleparam;
	map<string, vector<vector<int>>> ntangleparamidx; 
	map<string, vector<vector<double>>> adjangleparam; 
	map<string, vector<vector<int>>> adjangleparamidx;
    
    map<string, vector<vector<vector<double>>>> sugartorsion;
    map<string, vector<vector<double>>> bbtorsion;
    map<string, vector<pair<vector<int>, vector<double>>>> btorsionparam;

	map<string, vector<vector<double>>> hydparam;
	map<string, vector<vector<int>>> hydparamidx;
    map<string, vector<vector<double>>> hydaddangle;
	map<string, vector<vector<int>>> hydaddangleidx;
    map<string, vector<vector<double>>> hydaddlength;
	map<string, vector<vector<int>>> hydaddlengthidx;
    map<string, vector<pair<vector<int>, vector<double>>>> hydtorsionparam;

	map<string, vector<vector<int>>> molConnectone;
	map<string, vector<vector<int>>> molConnecttwo;
	map<string, vector<vector<int>>> molConnectthree;

	map<string, vector<vector<int>>> molAdjConnectone;
	map<string, vector<vector<int>>> molAdjConnecttwo;
	map<string, vector<vector<int>>> molAdjConnectthree;
	map<string, vector<vector<int>>> molAdjConnectoneidx;
	map<string, vector<vector<int>>> molAdjConnecttwoidx;
	map<string, vector<vector<int>>> molAdjConnectthreeidx;
	map<string, map<int, vector<double>>> atomfeatures;
    map<string, vector<double>> clashparam;

    map<int, map<string, map<string, vector<vector<double>>>>> stackparam;
    map<string, map<string, vector<vector<double>>>> bpparam;
    map<char, map<char, map<int, map<int, HydroBondParameters>>>> bphbond;
    map<char, map<char, map<int, map<int, BBHydroBondParameters>>>> bbhbond;
}ExtParams;

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

typedef struct FragmentAtoms
{
    Point3s coord;
	FragmentAtoms() = default;
    FragmentAtoms(const FragmentAtoms& other) : coord(other.coord) {}
}FragmentAtoms;

typedef struct FragmentResidue
{
    vector<FragmentAtoms> atoms;
    string restype;
	FragmentResidue() = default;
    FragmentResidue(const FragmentResidue& other)
        : atoms(other.atoms), restype(other.restype) {}
}FragmentResidue;

typedef struct FragmentStructure
{   
    vector<FragmentResidue> residues;

	FragmentStructure() = default;
    FragmentStructure(const FragmentStructure& other)
        : residues(other.residues) {}
}FragmentStructure;

double Random(void);

#endif
