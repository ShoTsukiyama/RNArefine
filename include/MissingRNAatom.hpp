#ifndef MiSSINGRNAATOM_HPP
#define MiSSINGRNAATOM_HPP 1

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
#include "Logger.hpp"

using namespace std;

bool checkMissingRNAatom(ResidueUnit &residue, 
    vector<string> &missing_atom_vec,
    vector<size_t> &nonstd_atom_vec, Logger &logger, const int option=0);
size_t countUniqAtom(vector<vector<double> >&xyz_list2);
bool AssembleBasePairFrag(const int c1, const int r1, const int c2, const int r2, 
    const ModelUnit &pdb_entry, map<string, vector<ModelUnit>> &ideal_rna,
    FragmentResidue &spresidue1, FragmentResidue &spresidue2, 
    int &templateidx);
void UpdateStemStructure(const int c1, const int r1, const int c2, const int r2, 
    const ModelUnit &pdb_entry, ChainUnit &fill_chain,
    vector<string>missing_atom_vec_r1, vector<string>missing_atom_vec_r2,
    FragmentResidue &spresidue1, FragmentResidue &spresidue2, 
    bool fullatom);
void ExtractNt(const int c, const int r, const ModelUnit &pdb_entry, FragmentResidue &spresidue);
void GetInitStr(int c, const ModelUnit &pdb_entry, FragmentStructure &initresidues, vector<vector<bool>> &inputtedstrflags);
bool AssembleNonBpNtFrag(const int c, const int r,
    const ModelUnit &pdb_entry, map<string, vector<ModelUnit>> &ideal_rna, FragmentResidue &spresidue, int templateidx);
void UpdateLoopStructure(const int c, const int r,
    const ModelUnit &pdb_entry, ChainUnit &fill_chain,
    vector<string>missing_atom_vec, FragmentResidue &spresidue,
    bool fullatom);
void AddMissingAtoms(ModelUnit &pdb_entry, int c, const ChainUnit &fill_chain, const map<size_t,vector<string> > &fill_residue_map, const vector<bool> &updateflag, int option);
bool MissingRNAatom(ModelUnit &pdb_entry, const size_t c,
    map<string, map<string, vector<ModelUnit>>> &ideal_rna,
    const vector<vector<vector<int>>> &stembps, const map<int, vector<pair<int, int>>> &loopnts, const map<int, vector<vector<vector<int>>>> &connectnts,
    map<int, int> &bp_template_idx, map<int, int> &nonbp_template_idx,
    FragmentStructure &inputtedstrall, vector<vector<bool>> &inputtedstrallflag,
    ExtParams &extparameters, Logger &logger, const int option=0);
bool MissingRNAatom(ModelUnit &pdb_entry, 
    map<string, map<string, vector<ModelUnit>>> &ideal_rna,
    const vector<vector<vector<int>>> &stembps, const map<int, vector<pair<int, int>>> &loopnts, const map<int, vector<vector<vector<int>>>> &connectnts,
    map<int, int> &bp_template_idx, map<int, int> &nonbp_template_idx,
    FragmentStructure &inputtedstrall, vector<vector<bool>> &inputtedstrallflag,
    ExtParams &extparameters, const string &outdir, const bool &logflag, const int option=0);


#endif
