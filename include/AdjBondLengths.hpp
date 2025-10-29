/* Purpose: Refine bond lengths and angles
Ideal bond and pseudobond lengths are taken from the Amber forcefield parameters.
*/

#ifndef AdjBondLengths_HPP
#define AdjBondLengths_HPP 1

#include <math.h>
#include <stdbool.h>
#include <vector>
#include <cstdlib>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <map>

#include "GeometryTools.hpp"
#include "PDBParser.hpp"
#include "Superpose.hpp"
#include "BondLengths.hpp"

void get_template_idx(const ModelUnit &pdb_entry, const int c1, const int c2, const int r1, const int r2, map<string, map<string, vector<DNTTEMPLATEUnit>>> ideal_rna, int &templateidx);
int fix_adj_bond_lengths(ModelUnit &pep, map<size_t, vector<int>> loop_vec, map<string, map<string, vector<DNTTEMPLATEUnit>>> ideal_rna, double tolerance);

#endif
