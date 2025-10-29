#ifndef RefineRNA_HPP
#define RefineRNA_HPP 1

#include <iostream>
#include "PDBParser.hpp"
#include "BondLengths.hpp"
#include "BaseConformation.hpp"
#include "GetInteraction.hpp"
#include "AtomicClashes.hpp"
#include "AdjBondLengths.hpp"

void RefineRNA(ModelUnit &pdb_entry, map<string, map<string, vector<ModelUnit>>> &ideal_rna, const vector<pair<double,vector<size_t> > >&bp_vec, map<int, int> &bp_template_idx, map<int, int> &nonbp_template_idx, double alpha, double beta, double tolerance, const string outdir, const bool &logflag);

#endif
