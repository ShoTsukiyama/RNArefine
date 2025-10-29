#ifndef FragmentRefine_HPP
#define FragmentRefine_HPP 1

#include <iostream>
#include "PDBParser.hpp"
#include "Superpose.hpp"
#include "MissingRNAatom.hpp"

bool fillFragmentInBasePair(const int c1, const int r1,
    const int c2, const int r2, const ModelUnit &pdb_entry, ChainUnit &fill_chain,
    map<string, vector<ModelUnit>> &ideal_rna, int &templateidx, bool fullatom);
bool fillFragmentLocal(const int c, const int r,
    const ModelUnit &pdb_entry, ChainUnit &fill_chain, map<string, vector<ModelUnit>> &ideal_rna, int &templateidx, bool fullatom);
bool FragRef(ModelUnit &pdb_entry, const size_t c,
    map<string, map<string, vector<ModelUnit>>> &ideal_rna,
    const vector<pair<double,vector<size_t>>>&bp_vec, map<int, int> &bp_template_idx, map<int, int> &nonbp_template_idx,
    const int option=0);   
bool FragRef(ModelUnit &pdb_entry, 
    map<string, map<string, vector<ModelUnit>>> &ideal_rna,
    const vector<pair<double,vector<size_t> > >&bp_vec, map<int, int> &bp_template_idx, map<int, int> &nonbp_template_idx,
    const int option=0);

#endif
