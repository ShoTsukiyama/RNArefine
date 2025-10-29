#ifndef PreEnergyRefinement_HPP
#define PreEnergyRefinement_HPP 1

#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include "PDBParser.hpp"
#include "Superpose.hpp"
#include "IdealRNA.hpp"
#include "PrepOpt.hpp"
#include "AttachHydrogens.hpp"
#include "CommonParams.hpp"

using namespace std;

void UpdateNtStructure(const int c1, const int r1, ModelUnit &pdb_entry, FragmentResidue &spresidue);

void ExtractNtPre(const int c, const int r, const ModelUnit &pdb_entry, FragmentResidue &spresidue);

bool PreEnergyRefine(ModelUnit &pdb_entry, const vector<FragmentStructure> &templatestrall, const vector<vector<vector<bool>>> &templatestrallflag, const map<string, map<string, vector<ModelUnit>>> &ideal_rna, const size_t c,
    vector<vector<int>> &stackmap, const map<string, map<string, vector<vector<vector<double>>>>> &stackdistmap, const vector<vector<vector<int>>> &stembps, const map<int, vector<pair<int, int>>> &loopnts, const map<int, vector<vector<vector<int>>>> &connectnts,
    ExtParams &extparameters);

bool PreEnergyRefine(ModelUnit &pdb_entry, const vector<FragmentStructure> &templatestrall, const vector<vector<vector<bool>>> &templatestrallflag, const map<string, map<string, vector<ModelUnit>>> &ideal_rna, vector<vector<int>> &stackmap, const map<string, map<string, vector<vector<vector<double>>>>> &stackdistmap,
    const vector<vector<vector<int>>> &stembps, const map<int, vector<pair<int, int>>> &loopnts, const map<int, vector<vector<vector<int>>>> &connectnts,
    ExtParams &extparameters);

#endif
