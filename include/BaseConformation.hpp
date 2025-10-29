#ifndef BaseConformation_HPP
#define BaseConformation_HPP 1

#include <iostream>
#include <string>
#include "PDBParser.hpp"
#include "Superpose.hpp"
#include "GetInteraction.hpp"

#include  "PDBParser.hpp"
#include "MissingRNAatom.hpp"

using namespace std;

bool SelectBPTemplate(ModelUnit &pdb_entry, map<string, vector<ModelUnit>> &ideal_rna, const vector<pair<double,vector<size_t>>> bp_vec, map<int, int> &bp_template_idx);
int fixBaseConformation(ResidueUnit &residue, ModelUnit&ideal_rna);
int fixBasePairConformation(ResidueUnit &residue1, ResidueUnit &residue2, ModelUnit&ideal_rna);
int fixBaseConformation(ModelUnit &pdb_entry, map<string, vector<ModelUnit>> &ideal_rna, const vector<pair<double,vector<size_t> > >&bp_vec, map<int, int> &bp_template_idx);

#endif