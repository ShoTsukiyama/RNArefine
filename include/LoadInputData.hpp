
#ifndef LoadInputData_HPP
#define LoadInputData_HPP 1

#include <iostream>
#include <filesystem>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include <sstream>
#include "PDBParser.hpp"
#include "GetInteraction.hpp"

using namespace std;
namespace fs = std::filesystem;

void get_multitemplate_info(const ModelUnit &pdb_entry, const string &tempstrdir, vector<FragmentStructure> &templatestrall, vector<vector<vector<bool>>> &templatestrallflag, int &atomic_detail, int &allowX, vector<string> neccatoms = {"  P ", " C4'", " C3'"});
void gettargetfilesfiles(string tdir, vector<string> &targetfiles, string ext=".pdb");
bool check_targetpdb(ModelUnit pdb_entry, vector<string> neccatoms = {"  P ", " C4'", " C3'"});

void get_multistack_info(const ModelUnit &pdb_entry, const string &stackdir, const string &seq, map<string, int> &chainstr2chainidx, map<int, map<int, int>> &pdbidx2seqidx, vector<vector<pair<double,vector<size_t>>>> &stack_list, vector<vector<vector<double>>> &stack_scoremat, vector<vector<vector<bool>>> &stack_flagmat, vector<vector<int>> &stack_ranges);
void get_multibp_info(const ModelUnit &pdb_entry, const string &bpdir, const string &seq, map<string, int> &chainstr2chainidx, map<int, map<int, int>> &pdbidx2seqidx, vector<vector<pair<double,vector<size_t>>>> &bps_list, vector<vector<pair<double,vector<size_t>>>> &noncbps_list, vector<vector<vector<double>>> &bps_scoremat, vector<vector<vector<bool>>> &bps_flagmat, vector<vector<int>> &bps_ranges);

//void get_multibp_info(const ModelUnit &pdb_entry, const string &bpdir, const string &seq, vector<pair<double,vector<size_t>>> &bps_list, vector<vector<double>> &bps_scoremat, double threshold=0.25);
//void get_multistack_info(const ModelUnit &pdb_entry, const string &stackdir, const string &seq, vector<pair<double,vector<size_t>>> &stack_list, vector<vector<double>> &stack_scoremat, double threshold=0.25);

bool loadContactmap(string contpath, string seq, vector<vector<double>> &contactmap, vector<vector<bool>> &contactmapflag, int &contactseqst, int &contactseqend);
void getmulticontactmap(const ModelUnit &pdb_entry, const string &contactdir, const string &seq, vector<vector<vector<double>>> &contactmaps, vector<vector<vector<bool>>> &contactmapflags, vector<vector<int>> &contactmapranges);

#endif
