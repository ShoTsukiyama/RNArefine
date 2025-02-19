// TestingFunc.h: interface for the TestingFunc class.
//
//////////////////////////////////////////////////////////////////////

//#if !defined(AFX_TESTINGFUNC_H__5ADC12DB_C425_45A5_A3B8_7DBBC34CAE4A__INCLUDED_)
//#define AFX_TESTINGFUNC_H__5ADC12DB_C425_45A5_A3B8_7DBBC34CAE4A__INCLUDED_
#pragma once
#include <vector>
#include <string>
#include <map>
using namespace std;

#include "LoadPDB.hpp"
#include "LoadParams.hpp"
#include "LoadDeepParams.hpp"
//#include "LoadBP.hpp"
#include "ParseSeq.hpp"
#include "SimulationBase.hpp"
#include "StrVis.hpp"
//#include "LoadST.hpp"

//extern vector<vector<Pairaa>> paa3;
//extern vector<Pairaa> initpaa3;
//extern vector<int> npaa3;
//extern int initnpaa3;
extern map< char, vector<Binnedparamtype>> TORSIONPARAM;
extern map< char, vector<Binnedparamtype>> INNERANGPARAM;
extern map< char, vector<Binnedparamtype>> LENGTHPARAM;
extern map< char, vector<Binnedparamtype>> SUBLENGTHPARAM;

extern map< char, vector<Binnedparamtype>> ADDINNERANGPARAM;
extern map< char, vector<Binnedparamtype>> ADDLENGTHPARAM;

extern map< char, vector<Binnedparamtype>> OPTORSIONPARAM;
extern map< char, vector<Binnedparamtype>> OPINNERANGPARAM;
extern map< char, vector<Binnedparamtype>> OPLENGTHPARAM;

extern map< char, vector<Binnedparamtype>> HYDTORSIONPARAM;
extern map< char, vector<Binnedparamtype>> HYDINNERANGPARAM;
extern map< char, vector<Binnedparamtype>> HYDLENGTHPARAM;

extern map<int, vector<Binnedparamtype>> BASEPAIRPARAM;
extern map<int, map<int, vector<Binnedparamtype>>> STACKINGPARAM;

extern vector<D3Binnedparamtype> BASEPAIRROT;
extern map<int, vector<D3Binnedparamtype>> STACKINGROT;

extern map<char, map<char, map<int, map<string, vector<vector<Binnedparamtype>>>>>> BASEPAIRHYDPARAM;

extern vector<double>VWPARAM;

extern map< char, vector<D2Binnedparamtype>> MULTITORSIONPARAM;
extern vector<vector<vector<Torsionclsparamtype>>> TORSIONLIST;
extern map< char, vector<vector<D2Binnedparamtype>>> TORSIONPAIRSPARAM;
extern vector<D2Binnedparamtype> TORSIONBBPAIRSPARAM;

extern map<char, vector<vector<bool>>> EXPCLASHATOMSNT;
extern map<char, map<char, vector<vector<bool>>>> EXPCLASHATOMSADJ;
extern map<char, map<char, vector<vector<bool>>>> EXPCLASHATOMSANY;

extern map<char, vector<vector<double>>> BASETEMPLATE; 

extern vector<vector<Binnedparamtype>> PUCKERPARAM;

extern vector<vector<vector<double>>> CONTACTSCORE;
extern vector<vector<vector<int>>> CONTACTBINARY;
extern vector<BPBasicInfo> BASEPAIRS;
extern vector<STBasicInfo> STACKINGPAIRS;

extern map<char, vector<vector<vector<PairLengInfo>>>> PAIRLENGNTPARAM;
extern map<char, map<char, vector<vector<vector<PairLengInfo>>>>> PAIRLENGADJPARAM;

extern map<char, vector<vector<vector<PairLengInfo>>>> BPATOMPAIRLENGNTPARAM;
extern map<char, map<char, vector<vector<vector<PairLengInfo>>>>> BPATOMPAIRLENGADJPARAM;
extern map<char, map<char, vector<vector<PairLengInfo>>>> BPATOMPAIRLENGPARAM;

extern map<char, vector<vector<PairLengInfo>>> ALLPAIRLENGNTPARAM;
extern map<char, map<char, vector<vector<PairLengInfo>>>> ALLPAIRLENGADJPARAM;

extern map<char, vector<vector<TorClsInfo>>> TORCLSPARAM;

extern map<int, map<char, vector<vector<Residueatoms>>>> SINGLENTTEMPLATE;
extern map<char, map<char, vector<vector<Residueatoms>>>> DILENTTEMPLATE;
extern map<char, map<char, vector<vector<Residueatoms>>>> BPLENTTEMPLATE;

extern map<char, map<char, map<char, vector<TemplateDistMap>>>> STLENTTEMPLATELENGTH;
extern map<char, map<char, vector<TemplateDistMap>>> DILENTTEMPLATELENGTH;
extern map<char, map<char, vector<TemplateDistMap>>> BPLENTTEMPLATELENGTH;

extern map<char, vector<D2Binnedparamtype>> OXHYDBPARAM;

extern map<int, map<char, map<char, vector<vector<Residueatoms>>>>> STACKTEMPLATE;

//#if _MSC_VER > 1000
//#endif // _MSC_VER > 1000

class RNArefineFunc
{
public:
    RNArefineFunc();
    //void rnarefine(string exepath, string targetstrfile, vector<string>templatefiles, string ssdir,  string outfile, string paramdir, int reftype, double cuttime1, double cuttime2);//single traj
    void rnarefine(string exepath, string targetstrfile, vector<string>templatefiles, double templatecoeff, vector<double>templateweights, string bppath, string stpath, string outfile, string paramdir, int reftype, double cuttime1, double cuttime2);
    //void rnarefine(string exepath, string targetstrfile, vector<string>templatefiles, string ssdir,  string outfile, string paramdir,double cuttime);//single traj
    bool judgetermination(map<int, int> tempidx2repidx, vector<double> energies);
    void swapmovement(map<int, int> &repidx2tempidx, map<int, int> &tempidx2repidx, vector<double> Temperatures, vector<double> energies);
    //void threadprocess(int start, int end, vector<RNASimulation> &rsims, vector<vector<Residueatoms>> &decstrs, vector<ofstream> &outputFiles, vector<int> &ms, vector<int> &fnows, vector<int> &dupds, vector<double> &oldenergies, vector<double> &newenergies, string outfile, vector<double> Temperatures, map<int, int> repidx2tempidx, int cycleidx, int numbb, vector<int> bbind, clock_t remcstart, clock_t remcfinish);
    virtual ~RNArefineFunc();

};

//#endif // !defined(AFX_TESTINGFUNC_H__5ADC12DB_C425_45A5_A3B8_7DBBC34CAE4A__INCLUDED_)


