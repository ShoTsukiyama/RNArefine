/* calculate backbone torsion angles (omega, psi, phi, kappa, alpha) */
#ifndef cssr_HPP
#define cssr_HPP 1
#include <cstring>
#include <sstream>
#include <vector>
#include <string>
#include <set>
#include <map>
#include "PDBParser.hpp"
#include "GeometryTools.hpp"

using namespace std;

void filter_bp(vector<pair<double,vector<size_t> > >&bp_vec);
bool bp_tor_score(
    const vector<double> &c1, const vector<double> &c2,
    const vector<double> &c3, const vector<double> &c4,
    const double mu, const double sd, const double tol, const double weight,
    double &nominator, double &denominator);
bool bp_ang_score(
    const vector<double> &c1, const vector<double> &c2,
    const vector<double> &c3, const vector<double> &c4,
    const double mu, const double sd, const double tol, const double weight,
    double &nominator, double &denominator);
bool bp_len_score(const vector<double>&c1, const vector<double>&c2, const double mu, const double sd, const double tol, const double weight, double &nominator, double &denominator);
bool bp_nn_score(const bool previnextj, const bool nextiprevj,
    const bool has_prev_ci,       const bool has_next_ci,
    const bool has_prev_cj,       const bool has_next_cj,
    const vector<double> &prev_ci, const vector<double> &next_ci,
    const vector<double> &prev_cj, const vector<double> &next_cj,
    const double mu, const double sd, const double tol, const double weight,
    double &nominator, double &denominator);
double mean_vec(const vector<double>&vec);
void cssr(const ModelUnit &pdb_entry, vector<vector<size_t> >&res_str_vec,
    vector<pair<double,vector<size_t> > >&bp_vec, 
    const bool interchain=false);
void getsselement(const ModelUnit &pdb_entry, vector<pair<double,vector<size_t> > > bp_vec, vector<vector<vector<int>>> &stembps, map<int, vector<pair<int, int>>> &loopnts, vector<vector<int>> &stemconnectnts, map<int, vector<vector<vector<int>>>> &connectnts);
string getseq(ModelUnit pdb_entry);
bool loadDistmap(string distpath, string seq, vector<vector<double>> &distmap, vector<vector<bool>> &distmapflag, int &distseqst, int &distseqend);
bool loadContactmap(string contpath, string seq, vector<vector<double>> &contactmap, vector<vector<bool>> &contactmapflag, int &contactseqst, int &contactseqend);
bool loadBP(string bppath, string seq, map<string, int> &chainstr2chainidx, map<int, map<int, int>> &pdbidx2seqidx, vector<pair<double,vector<size_t>>> &bp_vec, vector<pair<double,vector<size_t>>> &noncbp_vec, int &bpseqst, int &bpseqend);
bool loadST(string stpath, string seq, map<string, int> &chainstr2chainidx, map<int, map<int, int>> &pdbidx2seqidx, vector<pair<double,vector<size_t> > > &stack_vec, int &stackseqst, int &stackseqend);
void getstmap(const vector<pair<double,vector<size_t>>> &stack_vec, vector<vector<pair<int, double>>> &stackmap, const int &length);
void getmultidistmap(const ModelUnit &pdb_entry, const string &distdir, const string &seq, vector<vector<vector<double>>> &distmaps, vector<vector<vector<bool>>> &distmapflags, vector<vector<int>> &distmapranges);
void get_edge_loop(vector<pair<double,vector<size_t>>> &bp_vec, vector<pair<double,vector<size_t>>> &noncbp_vec, int numbb, vector<double> &nt_weight);

#endif
