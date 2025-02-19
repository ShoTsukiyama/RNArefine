
#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <map>
#include <fstream>
#include <sstream>
#include <filesystem>
#include "CommonPara.hpp"
#include "RNArefine.hpp"
#include "BasicFunc.hpp"

using namespace std;

class LoadParams
{
public:
	LoadParams();

	BasicFunc bf;

	bool loadvwd(string datafile);
	bool loadclashexpatoms(string basepath, map<char, vector<vector<bool>>> &stparameter1, map<char, map<char, vector<vector<bool>>>> &stparameter2, map<char, map<char, vector<vector<bool>>>> &stparameter3);

	void paramloadmain(string basePath, int numbb, string seq);

	bool load2ddistribution(string binfilepath, string energyfilepath, string binaryfilepath, string distributionfilepath, vector<D2Binnedparamtype> &stparameter);
	bool load3ddistribution(string binfilepath, string energyfilepath, string binaryfilepath, string distributionfilepath, vector<D3Binnedparamtype> &stparameter);

	bool loadstackingrot(string basepath, map<int, vector<D3Binnedparamtype>> &stparameter);
	bool loadbasepairrot(string basepath, vector<D3Binnedparamtype> &stparameter);

	bool loadstackingparameters(string basepath, map<int, map<int, vector<Binnedparamtype>>> &stparameter);
	bool loadstandardparameters(string basepath, map< char, vector<Binnedparamtype>> &stparameter, string paramname, string paramsetname);

	bool loadbaepairhydrogenparameters(string basepath, map<char, map<char, map<int, map<string, vector<vector<Binnedparamtype>>>>>> &stparameter);
	
	bool loadbasepairparameters(string basepath, map<int, vector<Binnedparamtype>> &stparameter);
	bool loadbasebininfo(string bininfofile, Binnedparamtype &param);

	void loadparamvecd(string filepath, vector<double> &vec);
	void loadparamveci(string filepath, vector<int> &vec);
	void loadparamvec2dd(string filepath, size_t dim1_size, size_t dim2_size, vector<vector<double>> &vec_2d);
	void loadparamvec2di(string filepath, size_t dim1_size, size_t dim2_size, vector<vector<int>> &vec_2d);

	bool loadbaseparamtype1cutoff(string file, map<int, Cutoffparamtype> &cutoffparam);
	bool loadbaseparamtype2cutoff(string file, map<int, Cutoffparamtype> &cutoffparam);
	bool loadbasehbondcutoff(string file, map<int, Cutoffparamtype> &cutoffparam);

	bool loadcombtorsionparameter(string basepath, map<char, vector<D2Binnedparamtype>> &stparameter);
	bool loadbasetemplate(string basepath, map<char, vector<vector<double>>> &stparameter);
	
	bool loadtorsionlist(string basepath, vector<vector<vector<Torsionclsparamtype>>> &stparameter);

	bool loadpuckerparameters(string basepath, vector<vector<Binnedparamtype>> &stparameter);

	//bool loadpairlengnt(string basepath, map<char, vector<vector<PairLengInfo>>> &stparameter);
	//bool loadpairlengadj(string basepath, map<char, map<char, vector<vector<PairLengInfo>>>> &stparameter);
	bool loadpairlengnt(string basepath, map<char, vector<vector<vector<PairLengInfo>>>> &stparameter);
	bool loadpairlengadj(string basepath, map<char, map<char, vector<vector<vector<PairLengInfo>>>>> &stparameter);

	bool loadsinglepairlengadj(string basepath, map<char, map<char, vector<vector<PairLengInfo>>>> &stparameter);
	bool loadsinglepairlengnt(string basepath, map<char, vector<vector<PairLengInfo>>> &stparameter);

	bool loadtorclsparam(string basepath, map<char, vector<vector<TorClsInfo>>> &stparameter);

	bool loadramachandranplot(string binfilepath, string ramachandranpath, vector<D2Binnedparamtype> &stparameter);
	bool loadtorcls2dparam(string basepath, map<char, vector<vector<D2Binnedparamtype>>> &stparameter);
	bool loadtorclsbb2dparam(string basepath, vector<D2Binnedparamtype> &stparameter);

	bool loadbppairlengnt(string basepath, map<char, vector<vector<vector<PairLengInfo>>>> &stparameter);
	bool loadbppairlengbp(string basepath, map<char, map<char, vector<vector<PairLengInfo>>>> &stparameter);
	bool loadbppairlengadj(string basepath, map<char, map<char, vector<vector<vector<PairLengInfo>>>>> &stparameter);

	bool loadsinglenttemplates(string basepath, map<int, map<char, vector<vector<Residueatoms>>>> &stparameter);
	void gettemplatedistmap(vector<Residueatoms> stparameter, TemplateDistMap &distparameter);
	bool loaddinttemplates(string basepath, map<char, map<char, vector<vector<Residueatoms>>>> &stparameter, map<char, map<char, vector<TemplateDistMap>>> &distparameter);
	bool loadnptemplates(string basepath, map<char, map<char, vector<vector<Residueatoms>>>> &stparameter, map<char, map<char, vector<TemplateDistMap>>> &distparameter);

	bool loadstacktemplates(string basepath, map<int, map<char, map<char, vector<vector<Residueatoms>>>>> &stparameter, map<char, map<char, map<char, vector<TemplateDistMap>>>> &distparameter);

	void coord2anglen(Point3s tmp2, Point3s tmp3, Point3s tmp4, float &len, float &ang);
	void coord2len(Point3s tmp3, Point3s tmp4, float &len);
	void coord2params(Point3s tmp1, Point3s tmp2, Point3s tmp3, Point3s tmp4, float &len, float &ang, float &dih);
	bool gettemplategeoparams(vector<Residueatoms> &stparameter);

	bool loadoxhydbparams(string basepath, map<char, vector<D2Binnedparamtype>> &stparameter);

	virtual ~LoadParams();
};

