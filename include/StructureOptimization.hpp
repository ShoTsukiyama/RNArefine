#ifndef STRUCTUREOPTIMIZATION_HPP
#define STRUCTUREOPTIMIZATION_HPP 1

#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include "MissingRNAatom.hpp"
#include "AttachHydrogens.hpp"
#include "CommonParams.hpp"

using namespace std;

void UpdateStructure(const int c, const int r, ModelUnit &pdb_entry, FragmentResidue &spresidue);
void GetHydrogenInfo(const int c, const int r, ModelUnit pdb_entry, FragmentResidue &spresidue, vector<vector<double>> &hydrogencoords);
bool StructureOptimization(ModelUnit &pdb_entry, const size_t c, const map<string, map<string, vector<ModelUnit>>> &ideal_rna,
    const FragmentStructure &inputtedstrall, const vector<vector<bool>> &inputtedstrallflag,
    const vector<FragmentStructure> &templatestrall, const vector<vector<vector<bool>>> &templatestrallflag,
    const vector<pair<double,vector<size_t>>> &bp_vec, const map<int, vector<pair<int, int>>> &loopnts, const vector<pair<double,vector<size_t>>> &noncbp_vec, 
    const vector<vector<int>> &stemconnectnts, const map<int, vector<vector<vector<int>>>> &connectnts,
    const vector<vector<pair<int, double>>> &stackmap, const map<string, vector<vector<vector<double>>>> &bpdistmap, const map<string, map<string, vector<vector<vector<double>>>>> &stackdistmap, const ExtParams &extparameters, const vector<double> &nt_weight, 
    const vector<vector<vector<double>>> &hydrogencoords, const double maxtime, Logger &logger, const int option=0);
bool StructureOptimization(ModelUnit &pdb_entry, const map<string, map<string, vector<ModelUnit>>> &ideal_rna, 
    const FragmentStructure &inputtedstrall, const vector<vector<bool>> &inputtedstrallflag,
    const vector<FragmentStructure> &templatestrall, const vector<vector<vector<bool>>> &templatestrallflag,
    const vector<pair<double,vector<size_t>>> &bp_vec, const map<int, vector<pair<int, int>>> &loopnts, const vector<pair<double,vector<size_t>>> &noncbp_vec, 
    const vector<vector<int>> &stemconnectnts, const map<int, vector<vector<vector<int>>>> &connectnts,
    const vector<vector<pair<int, double>>> &stackmap, const map<string, vector<vector<vector<double>>>> &bpdistmap, const map<string, map<string, vector<vector<vector<double>>>>> &stackdistmap, const ExtParams &extparameters, const vector<double> &nt_weight, 
    const vector<vector<vector<double>>> &hydrogencoords, const double maxtime, const string &outdir, const bool &logflag, const int option=0);

#endif
