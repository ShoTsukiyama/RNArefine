// TestingFunc.h: interface for the TestingFunc class.
//
//////////////////////////////////////////////////////////////////////

//#if !defined(AFX_TESTINGFUNC_H__5ADC12DB_C425_45A5_A3B8_7DBBC34CAE4A__INCLUDED_)
//#define AFX_TESTINGFUNC_H__5ADC12DB_C425_45A5_A3B8_7DBBC34CAE4A__INCLUDED_
#pragma once
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <map>
#include <iterator>
using namespace std;

#include "CommonPara.hpp"

//#if _MSC_VER > 1000
//#endif // _MSC_VER > 1000

class LoadPDB
{
public:
    LoadPDB();

    string pdbname;
	//char *aminoinchain[maxnumchains];// aminos in each chain

	int numpromod;
	int allocpromod;//許容されるPDBファイル内の最大MODEL数を格納する変数
	int allocproseq;//許容されるatomのレコード数を格納する変数
	bool flagfirstmodel;
	bool flagbeginmodel;
	int numproseq;//currentのatomのレコードを格納するindexを格納する変数
	bool flagter;//TARレコードであるかどうかのindexを格納する変数

	//vector<atom> proseq;//atomの情報を格納する変数
	//vector<promodel> promod;//あるMODELの情報を格納する変数

	int numbb, effbb;
	//vector<boneinfo> bb;
	//vector<dihedral> diheang;

    int dealline(string st,vector<Atom> &proseq, vector<Promodel> &promod);
    bool loadpdb(string pdbfile, vector<Atom> &proseq, vector<Promodel> &promod);
    void extractbb(vector<Atom> &proseq, vector<Promodel> &promod, vector<Boneinfo> &bb, int model, int chain);
    virtual ~LoadPDB();

};

//#endif // !defined(AFX_TESTINGFUNC_H__5ADC12DB_C425_45A5_A3B8_7DBBC34CAE4A__INCLUDED_)


