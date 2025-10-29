/* Purpose: parse PDB file into data structure similar to Bio.PDB in biopython
 * (model - chain - residue - atom). */
#ifndef PDBParser_HPP
#define PDBParser_HPP 1

#pragma once

#include <vector>
#include <cstdlib>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <filesystem>

#include "pstream.h"
#include "CommonParams.hpp"

using namespace std;

struct AtomUnit    // struct for each atom entry
{
    string name;       // atom name
    vector<double> xyz; // coordinate
    bool movable;       // 0 if atom is in input, 1 if atom is added by the program
    bool input;
};

struct ResidueUnit // struct for each residue
{
    bool het;               // true - HETATM, false - ATOM
    int resi;               // residue sequence number
    char icode;             // insertion code
    string resn;            // residue name
    vector<AtomUnit> atoms; // list of atoms
};

struct ChainUnit  // struct for each chain
{
    string chainID_full;          // chain ID, might be more than 1 char
    char chainID;                 // short chain ID, must be 1 char
    string sequence;              // sequence converted from CA coordinate
    string sarst;                 // SARST (Structural similarity search
                                  // Aided by Ramachandran Sequential 
                                  // Transformation) code
    string ss;                    // secondary structure
    vector<ResidueUnit> residues; // list of residues
};

struct ModelUnit  // struct for each model in mult-model PDB
{
    vector<ChainUnit> chains; // list of chains
};

struct DNTTEMPLATEUnit
{
    ModelUnit model;
    vector<vector<double>> pwdist_mean;
    vector<vector<double>> pwdist_std;
};

int parse_pdb_line(const string line,ModelUnit &pep, ChainUnit &chain,
    ResidueUnit &residue, AtomUnit &atom, map<char,string> &chainIDmap,
    const int atomic_detail=2,const int allowX=1);
ModelUnit read_pdb_structure(const char *filename, map<string, int> &chainstr2chainidx, map<int, map<int, int>> &pdbidx2seqidx, const int atomic_detail=2,const int allowX=1);
string write_pdb_structure(ChainUnit &chain,int &i);
void write_pdb_structure(const char *filename,ChainUnit &chain, string &initlabel);
string write_pdb_structure(ModelUnit &pep, string &initlabel);
void write_pdb_structure(const char *filename,ModelUnit &pep, string &initlabel);
void standardize_pdb_order(ModelUnit &pep);
void reindex_pdb(const int startindex,ChainUnit& chain);
void reindex_pdb(const int startindex,ModelUnit& pep);
char aa3to1(const string resn,const int convertX=2);
string pdb2fasta(ChainUnit& chain);
string aa1to3(const char aa);
string pdb2fasta(ModelUnit& pep,const string PDBid="",const int ShowSeqLen=0);
int has_atom_name(ResidueUnit residue,string name=" CA ");
void remove_sidechain(ResidueUnit& residue,int atomic_detail=1);
void remove_sidechain(ChainUnit& chain,int atomic_detail=1);
void remove_sidechain(ModelUnit& pep,int atomic_detail=1);

#endif
