
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

using namespace std;

class LoadDeepParams
{
public:
	LoadDeepParams();

	bool loadcontactmapfile(string filepath, vector<vector<double>> &contactscore, vector<vector<int>> &contactbinary);
	bool loadcontactmap(string basepath1, string basepath2, int seqnum, vector<vector<vector<double>>> &contactscore, vector<vector<vector<int>>> &contactbinary);
	bool loadBP(string bppath, string seq, vector<BPBasicInfo> &basepairs);
	bool loadST(string stpath, vector<STBasicInfo> &stackingpairs);

	//void paramloadmain(string ct1path, string ct2path, int numbb, string seq);
	//void paramloadmain(string bppath, string stpath, string ct1path, string ct2path, int numbb, string seq);
	void paramloadmain(string bppath, string stpath, int numbb, string seq);

	virtual ~LoadDeepParams();
};

