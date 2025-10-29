/* IdealRNA.hpp - Idealized RNA nucleotide conformations */

#ifndef IdealRNA_HPP
#define IdealRNA_HPP 1

#include <string>
#include <map>
#include <vector>
#include <filesystem>
#include "PDBParser.hpp"

namespace fs = std::filesystem;
using namespace std;

void standardize_atom_order(ModelUnit &pep);
vector<string> getfilelist(string path);
map<string, vector<ModelUnit>> loadtemplatestrs(string basepath, vector<string> subdirs);
map<string, map<string, vector<ModelUnit>>>parse_ideal_rna(string paramdir);
bool load_pairwiseinfo(string filepath, string res1, string res2, vector<vector<double>> &pmat);
map<string, vector<DNTTEMPLATEUnit>> loadpwstrs(string basepath, vector<string> subdirs);
map<string, map<string, vector<DNTTEMPLATEUnit>>>parse_ideal_rna_distmap(string paramdir);
void get_distmap(const map<string, vector<ModelUnit>> &rna_str, map<string, vector<vector<vector<double>>>> &stackdistmap);
void calc_interaction_distmap(const ModelUnit &rna_str, vector<vector<double>> &pmat);
void get_score_map(string seq, const vector<pair<double,vector<size_t>>> &bp_vec, vector<double> &score_vec);

#endif
