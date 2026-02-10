/* IdealRNA.hpp - Idealized RNA nucleotide conformations */

#ifndef LOADGEOPARAMS_HPP
#define LOADGEOPARAMS_HPP 1

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <filesystem>
#include "PDBParser.hpp"
#include "CommonParams.hpp"

using namespace std;
namespace fs = std::filesystem;

bool loadlength(string basepath, map<string, vector<vector<double>>> &ntlengthparam, map<string, vector<vector<double>>> &ntlengthstdval, map<string, vector<vector<double>>> &adjlengthparam, map<string, vector<vector<double>>> &adjlengthstdval);
bool loadangle(string basepath, map<string, vector<vector<double>>> &ntangleparam, map<string, vector<vector<double>>> &ntanglestdval, map<string, vector<vector<int>>> &ntangleparamidx, map<string, vector<vector<double>>> &adjangleparam, map<string, vector<vector<double>>> &adjanglestdval, map<string, vector<vector<int>>> &adjangleparamidx);
bool loadhydtorsion(string basepath, map<string, vector<pair<vector<int>, vector<double>>>> &hydtorsionparam);
bool loadhydrogenparam(string basepath, map<string, vector<vector<double>>> &hydparam, map<string, vector<vector<int>>> &hydparamidx);
bool getConnect(string basepath, map<string, vector<vector<int>>> &molConnectone, map<string, vector<vector<int>>> &molConnecttwo, map<string, vector<vector<int>>> &molConnectthree);
bool getConnectAdj(string basepath, map<string, vector<vector<int>>> &molAdjConnectone, map<string, vector<vector<int>>> &molAdjConnecttwo, map<string, vector<vector<int>>> &molAdjConnectthree, map<string, vector<vector<int>>> &molAdjConnectoneidx, map<string, vector<vector<int>>> &molAdjConnecttwoidx, map<string, vector<vector<int>>> &molAdjConnectthreeidx);
bool loadfeatures(string basepath, map<string, map<int, vector<double>>> &atomfeatures);
void loadparameters(string paramdir, ExtParams &extparameters);
bool loadlength(string basepath, map<string, vector<vector<double>>> &ntlengthparam, map<string, vector<vector<double>>> &adjlengthparam);
bool loadangle(string basepath, map<string, vector<vector<double>>> &ntangleparam, map<string, vector<vector<int>>> &ntangleparamidx, map<string, vector<vector<double>>> &adjangleparam, map<string, vector<vector<int>>> &adjangleparamidx);
bool loadhydrogenparam(string basepath, map<string, vector<vector<double>>> &hydparam, map<string, vector<vector<int>>> &hydparamidx);
bool loadhydrogenaddangle(string basepath, map<string, vector<vector<double>>> &hydparam, map<string, vector<vector<int>>> &hydparamidx);
bool loadhydrogenaddlength(string basepath, map<string, vector<vector<double>>> &hydparam, map<string, vector<vector<int>>> &hydparamidx);
bool getConnect(string basepath, map<string, vector<vector<int>>> &molConnectone, map<string, vector<vector<int>>> &molConnecttwo, map<string, vector<vector<int>>> &molConnectthree);
bool getConnectAdj(string basepath, map<string, vector<vector<int>>> &molAdjConnectone, map<string, vector<vector<int>>> &molAdjConnecttwo, map<string, vector<vector<int>>> &molAdjConnectthree, map<string, vector<vector<int>>> &molAdjConnectoneidx, map<string, vector<vector<int>>> &molAdjConnecttwoidx, map<string, vector<vector<int>>> &molAdjConnectthreeidx);
bool loadfeatures(string basepath, map<string, map<int, vector<double>>> &atomfeatures);
void loadparameters(string paramdir, ExtParams &extparameters);
bool loadtorsion(string basepath, map<string, vector<pair<vector<int>, vector<double>>>> &hydtorsion);
bool loadstackparameters(string basepath, map<int, map<string, map<string, vector<vector<double>>>>> &stparameter);
bool loadbpparameters(string basepath, map<string, map<string, vector<vector<double>>>> &stparameter);
bool loadbphbondparameters(string basepath, map<char, map<char, map<int, map<int, HydroBondParameters>>>> &stparameter);
bool loadbbhbondparameters(string basepath, map<char, map<int, map<int, BBHydroBondParameters>>> &stparameter);
#endif












