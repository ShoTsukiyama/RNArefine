// EMRefinement.h: interface for the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

#pragma once
#include <stdio.h>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <deque>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "ParseSeq.hpp"
#include "LoadPDB.hpp"
#include "CommonPara.hpp"
#include "BasicFunc.hpp"
#include "RNArefine.hpp"

using namespace std;
using namespace Eigen;

class RNASimulation 
{
public:
	RNASimulation();

	double TK;
	double enelist[40],enelistbk[40];
	double contwt;
	int curcycle;
	double baseenecoeff=1;
	int enetype=0;
	//int baseflag=0;
	double initfape;
	double fapecoff;

	vector<double> dihevec;
	vector<int> lpedpos;

	//vector<BPBasicInfo> basepairs;
	//vector<STBasicInfo> stackingpairs;
	map<int, map<int, int>> bpmap;

	double sseneratio[2];

	vector<vector<Pairaa>> paa3;
	vector<Pairaa> initpaa3;
	vector<int> npaa3;
	int initnpaa3;

	//double initfape;
	//double fapecoff;
	//FAPEobj fapemat, fapematbk;
	vector<vector<Rescontact>> rppmat, rppmatbk; //clashの情報を格納するための変数
	vector<vector<int>> nmap, nmapbk;
	vector<Residueatoms> restr;
	vector<vector<vector<double>>> distmap, distmapbk;
	vector<vector<vector<double>>> ntpairdist, ntpairdistbk;
	vector<vector<vector<double>>> adjpairdist, adjpairdistbk;
	vector<double> distene, distenebk;
	vector<Phycontact> contactene, contactenebk;
	vector<vector<double>> multidistene, multidistenebk;
	vector<vector<bool>> mapflag, mapflagbk;
	vector<vector<bool>> mapbbflag, mapbbflagbk;
	Energyterms lp1,lp2;

	map<int, int> bpmapid2idx;
	map<int, int> bsmapid2idx;
	map<int, int> bppairednt;
	map<int, int> stpairednt;
	vector<Basecontact> bpenergy, bpenergybk;
	vector<Stackcontact> bsenergy, bsenergybk;
	double coefftemplate;

	//vector<vector<bool>> toroutflag;
	//vector<int> torclsidx;

	//vector<vector<double>> ntstr,tpstr;//ntstr:refined strutureと共通の残基におけるreference structureのC4'のx、y、z座標(x1,x2,x2...y1,y2,y3...z1,z2,z3)
	//vector<vector<bool>> ntflag;//reference structureの各残基について、refined structureと共通した残基であるかどうか（共通:true、共通でない:false）
	//vector<double> initntstr,inittpstr;//ntstr:refined strutureと共通の残基におけるreference structureのC4'のx、y、z座標(x1,x2,x2...y1,y2,y3...z1,z2,z3)
	//vector<bool> initntflag;//reference structureの各残基について、refined structureと共通した残基であるかどうか（共通:true、共通でない:false）

	//vector<vector<Pairaa>> paa3;
	//vector<Pairaa> initpaa3;
	//vector<int> npaa3;
	//int initnpaa3;
	
	int tempnum;
	//vector<vector<Pairaa>> tmppaa3;
	//vector<int> tmpnpaa3;
	//vector<Stemtemplates> stemstr;
	//map<int, vector<Looptemplates>> loopstr;

	//ParseSeq ps;
	BasicFunc bf;
	//LoadPDB pp;

	int numdim;
	vector<int> lentmp;
	int initlentmp;

	//bool loadST(string stpath);
	//bool loadBP(string bppath, string seq);

	void setStr(vector<Residueatoms> &decstr, vector<Atom> proseq, vector<Boneinfo> bb, string seq);
	void fsetinitialarray(int numseq);
	double gettemplature(double pos, double numseq, double scale);
	bool checkoutlier(double val, vector<double> boundary, vector<int> binary);
	void checkupdate(vector<Residueatoms> decstr, int numseq);

	//Param2Coord
	bool param2coordbb(vector<Residueatoms> &decstr,int istart,int iend);
	void param2coordsugar(vector<Residueatoms> &decstr,int lstart,int lend);
	bool param2coordbase(vector<Residueatoms> &decstr, int istart, int iend);
	void param2coordhydrogen(vector<Residueatoms> &decstr,int lstart,int lend);

	//Coord2Param
	void coord2params(Point3s tmp1, Point3s tmp2, Point3s tmp3, Point3s tmp4, float &len, float &ang, float &dih);
	void coord2dangle(Point3s tmp1, Point3s tmp2, Point3s tmp3, Point3s tmp4, float &dih);
	void coord2ang(Point3s tmp2, Point3s tmp3, Point3s tmp4, float &ang);
	void coord2anglen(Point3s tmp2, Point3s tmp3, Point3s tmp4, float &len, float &ang);
	void coord2parambblocal(vector<Residueatoms> &decstr,int pos);
	void coord2len(Point3s tmp3, Point3s tmp4, float &len);
	void coord2parambb(vector<Residueatoms> &decstr,int start, int end,int simmode);
	void coord2paramsugar(vector<Residueatoms> &decstr,int start,int end,int simmode);
	void coord2paramadd(vector<Residueatoms> &decstr,int start,int end);
	void coord2paramhydrogen(vector<Residueatoms> &decstr,int start,int end,int simmode);
	void coord2paramhydrogenBase(vector<Residueatoms> &decstr,int start,int end,int simmode);

	//Save model
	bool writepdb(vector<Residueatoms> decstr,int numseq,vector<int> bbres, string outname,int ftype);

	//main function of energy calculation
	void getoutenelength(double val, double mean, double lower, double upper,  double slower, double supper, double k, double &diff, int &bdiff);
	void getouteneangle(double val, double mean, double lower, double upper,  double slower, double supper, double k, double &diff, int &bdiff);
	void getoutenetorangle(double val, double mean, double lower, double upper, double k, double &diff, int &bdiff);
	
	//other functions for energy calculations
	void resetflag(vector<vector<bool>> &mapflagrpp, vector<vector<bool>> &mapbbflagrpp);
	int getbin(double val, vector<double> boundary);
	void getoutene(double val, double mean, vector<double> boundary, vector<double> energy, vector<int> binary, double &diff, int &bdiff);
	void getoutenetorsion(double val, double mode, vector<double> boundary, vector<double> energy, vector<int> binary, double k, double &diff, int &bdiff);
	void getoutlier(double val, vector<double> boundary, vector<double> energy, vector<int> binary, double &diff, int &bdiff);
	double fcalcrmsdenergy(vector<Residueatoms> decstr,int numseq, Energyterms &lp,vector<vector<Rescontact>> &rpp, vector<Basecontact> &bprpp, vector<Stackcontact> &bsrpp, vector<vector<int>> &nmaprpp, vector<vector<vector<double>>> &distmaprpp, vector<vector<vector<double>>> &ntpairdistrpp, vector<vector<vector<double>>> &adjpairdistrpp, vector<double> &distenerpp, vector<vector<double>> &multidistenerpp, vector<Phycontact> &contactenerpp, vector<vector<bool>> &mapflagrpp, vector<vector<bool>> &mapbbflagrpp, bool flaginit);

	//FAPE
	//double calcFAPEthresh(double dist);
	//void fapesingleref(vector<Atom> proseq,vector<Boneinfo> bb, int numbb,vector<Boneinfo> ibb, int inumbb,int numseq);
	//void fapemultipleref(vector<Atom> proseq,vector<Boneinfo> bb, int numbb,vector<Boneinfo> ibb, int inumbb,int numseq, int refidx);
	//void getfape(vector<Residueatoms> decstr,int numseq,FAPEobj &faperpp,bool flagacc);
	//void fapeenergysingle(vector<Residueatoms> decstr,int numseq,double &fene,Energyterms &lp,FAPEobj &faperpp,bool flagacc);
	//void fapeenergymultiple(vector<Residueatoms> decstr,int numseq,double &fene,Energyterms &lp,FAPEobj &faperpp,bool flagacc);

	//energy for dist
	void distsingleref(vector<Atom> proseq,vector<Boneinfo> bb, int numbb,vector<Boneinfo> ibb, int inumbb,int numseq);
	void distmultipleref(vector<Atom> proseq,vector<Boneinfo> bb, int numbb,vector<Boneinfo> ibb, int inumbb,int numseq, double weight, int refidx);
	void distenergysingle(vector<Residueatoms> decstr,int numseq,double &fene,Energyterms &lp, vector<vector<vector<double>>> distmaprpp, vector<double> &distenerpp, vector<vector<bool>> mapbbflagrpp,bool flagacc);
	void distenergymultiple(vector<Residueatoms> decstr,int numseq,double &fene,Energyterms &lp, vector<vector<vector<double>>> distmaprpp, vector<vector<double>> &multidistenerpp, vector<vector<bool>> mapbbflagrpp, bool flagacc);

	//energy for contact
	//void contactenergysingle(vector<Residueatoms> decstr,int numseq,double &fene,Energyterms &lp, vector<vector<vector<double>>> distmaprpp, vector<double> &contactenerpp, vector<vector<bool>> mapbbflagrpp,bool flagacc);
	void contactenergysingle(vector<Residueatoms> decstr, int numseq, double &fene1, double &fene2, Energyterms &lp, vector<vector<vector<double>>> distmaprpp, vector<Phycontact> &contactenerpp, vector<vector<bool>> mapbbflagrpp,bool flagacc);

	//energy for clash and hydrogen bonding of oxigen
	void fenergyclash(vector<Residueatoms> decstr,int numseq,double &nene,double &fene, double &hfene ,Energyterms &lp,vector<vector<Rescontact>> &rpp, vector<vector<int>> &nmaprpp, vector<vector<vector<double>>> &distrpp, vector<vector<vector<double>>> &ntpairdistrpp, vector<vector<vector<double>>> &adjpairdistrpp, vector<vector<bool>> mapflagrpp, vector<vector<bool>> mapbbflagrpp,bool flagacc);
	double fenergyoxhydb(vector<Residueatoms> decstr, int basei, int basej, int ii, int jj, double tdist);

	//energy for multiple torsions
	void fenergymultidihedral(vector<Residueatoms> decstr,int numseq,double &nene,double &fene, Energyterms &lp,bool flagacc);

	//energy for bond length
	void fenergybondlength(vector<Residueatoms> decstr,int numseq,double &nene1,double &fene1,Energyterms &lp,bool flagacc);
	
	//energy for inner angles
	void fenergybondangle(vector<Residueatoms> decstr, int numseq,double &nene,double &fene,Energyterms &lp, bool flagacc);
	
	//energy for torsion angles
	double getdisttor(double val, double mean);
	void fenergydihedral(vector<Residueatoms> decstr,int numseq,double &nene,double &fene,Energyterms &lp,bool flagacc);

	//energy for paired torsion angles
	void fenergypairdihedral(vector<Residueatoms> decstr,int numseq,double &fene, Energyterms &lp,bool flagacc);

	//energy for paired backbone torsion angles
	void fenergypairbbdihedral(vector<Residueatoms> decstr,int numseq,double &fene, Energyterms &lp,bool flagacc);

	//energy for adjacent pair length
	//void fenergyadjpairlengths(vector<Residueatoms> decstr,int numseq,double &nene,double &fene, Energyterms &lp, vector<vector<vector<double>>> adjpairdistrpp, vector<vector<bool>> mapflagrpp, bool flagacc);

	//energy for nucleotide pair length
	//void fenergyntpairlengths(vector<Residueatoms> decstr,int numseq,double &nene,double &fene, Energyterms &lp, vector<vector<vector<double>>> ntpairdistrpp, vector<vector<bool>> mapflagrpp, bool flagacc);

	//energy for clustered torsion angles
	void getoutenetorsioncls(double val, double mean, double lower, double upper, int rangetype, double &diff, int &bdiff);
	void fenergyclsdihedral(vector<Residueatoms> decstr,int numseq,double &nene,double &fene, Energyterms &lp,bool flagacc);

	//energy for base-paring
	void getbaseangles(Point3s res1N91, Point3s res1C42, Point3s res1C86, Point3s res2N91, Point3s res2C42, Point3s res2C86, double &r, double &alpha, double &beta, double &gamma);
	void getbaselength(Point3s res1N91, Point3s res2N91, double &r);
	//void getnormalvecangle(Point3s center1, Point3s res11, Point3s res12, Point3s center2, Point3s res21, Point3s res22, double &baseangle);
	void getnormalvec(Point3s center1, Point3s res11, Point3s res12, Point3s &z);
	void getbasecenterpoint(Point3s *ptv, char aat, Point3s &center);
	void getcenter2centervec(Point3s center1, Point3s center2, Point3s &vec);
	void getoutenebaselength(double val, double mean, double lower, double upper, double k, double &diff, int &bdiff);
	int selectmovebase(int &mnt1, int &mnt2, int numseq);
	void fenergybasepair(vector<Residueatoms> decstr,int numseq,double &nene,double &fene,Energyterms &lp, vector<Basecontact> &baserpp, vector<vector<bool>> mapflagrpp,bool flagacc);
	void getplanebp(Point3s pointA, Point3s pointB, Point3s pointC, double &a, double &b, double &c, double &d);
	double getdistfromplane(double a, double b, double c, double d, Point3s pointP);
	double getbptorsionene(double val, vector<double> boundary, vector<double> energy);
	double bpplaneene(vector<Residueatoms> decstr, int basei, int basej);
	double angleBetweenVectors(Point3s a, Point3s b);

	//energy for base-staking
	void getoutenestacklength(double val, double mean, double lower, double upper, double k, double &diff, int &bdiff);
	void fenergystacking(vector<Residueatoms> decstr,int numseq,double &nene,double &fene,Energyterms &lp, vector<Stackcontact> &stackrpp, vector<vector<bool>> mapflagrpp, bool flagacc);
	double getsttorsionene(double val, vector<double> boundary, vector<double> energy);

	//energy by single template structure
	double calcsingletempene(Residueatoms restr);
	void fenergysinglenttemplate(vector<Residueatoms> decstr,int numseq,double &fene, Energyterms &lp, vector<vector<bool>> mapflagrpp, bool flagacc);
	
	//energy by stacking template structure
	double calcstacktempene(Residueatoms restr1, Residueatoms restr2, int sttype);
	void fenergystacktemplate(vector<Residueatoms> decstr,int numseq,double &fene,Energyterms &lp, vector<Stackcontact> &stackrpp, vector<vector<bool>> mapflagrpp,bool flagacc);

	//energy by base pair template structure
	double calcbptempene(Residueatoms restr1, Residueatoms restr2);
	void fenergynoncbasepairtemp(vector<Residueatoms> decstr,int numseq,double &fene,Energyterms &lp, vector<Basecontact> &baserpp, vector<vector<bool>> mapflagrpp,bool flagacc);

	//energy for pucker
	//void getoutenetorpucker(double val, double mean, double lower, double upper, double k, double &diff, int &bdiff);
	//void fenergypucker(vector<Residueatoms> decstr,int numseq,double &nene,double &fene, Energyterms &lp, vector<vector<bool>> mapflagrpp, bool flagacc);

	//energy for hydrogen bonding of base-paring
	void getoutenebasepairhydlength(double val, double mean, double lower, double upper, double k, double &diff, int &bdiff);
	void getoutenebasepairhydangle(double val, double mean, double lower, double upper, double k, double &diff, int &bdiff);
	void calchbondene(vector<Residueatoms> decstr, char daat, char aaat, int dbase, int abase, int pairtype, string hbtype, int &nn, double &fn, vector<double> &hbatoms);
	void fenergybasepairhyd(vector<Residueatoms> decstr,int numseq,double &nene,double &fene,Energyterms &lp, vector<Basecontact> &baserpp, vector<vector<bool>> mapflagrpp,bool flagacc);

	//main functions of movements
	bool fmcfragsweepcom(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,int &movtype,double Bt);
	void selectrange(double l, double r, double &nl, double &nr);
	double getnewparam(vector<int> distr, int total, vector<double> center, double span);
	double getnewparam(vector<double> logdistr, double logtotal, vector<double> center, double span);
	double getnewparambinary(vector<int> binary, int total, vector<double> center, double span);
	double getrandomnewparam(double l, double r);
	double getrandomnewtor(double l, double r, int rangetype);

	//Movement of torsion angles
	double getmovtorval(double val, double std, double k);
	bool fmovementtor(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt, bool flagacc);
	void settorflap(double numseq, double trandpos, double aat, int dihtype);
	void settorntflap(double numseq, double trandpos, double aat, int dihtype);

	//Movement of angles
	double getmovangval(double val, double mean, double std, double k);
	bool fmovementang(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc);
	void setangflap(double numseq, double trandpos, double aat, int angtype);
	void setangntflap(double numseq, double trandpos, double aat, int angtype);

	//Movement of angles in local structure
	bool fmovementanglocal(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc);
	void changelocalang(vector<Residueatoms> &decstr, int trandpos, int angtype, int numseq);
	void setangleflaplocal(double numseq, int trandpos, double aat, int angtype);

	//Movement of angles with Metropolis method
	void setangntmtflap(double numseq, double trandpos, double aat, vector<int> targetangles);
	bool fmovementangmt(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc);

	//Movement of bond length
	double getmovlenval(double val, double mean, double std, double k);
	bool fmovementlen(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc);
	void setlenflap(double numseq, int trandpos, double aat, int lentype);
	void setlenntflap(double numseq, int trandpos, double aat, int lentype);

	//Movement of length in local structure
	void setlenflaplocal(double numseq, int trandpos, double aat, int lentype);
	void coord2paramlocal(vector<Residueatoms> &decstr, int trandpos, int moveatom, int start, int end);
	void coord2paramsugarlocal(vector<Residueatoms> &decstr, int trandpos, int moveatom, int start, int end);
	void changelocallen(vector<Residueatoms> &decstr, int trandpos, int lentype, int numseq);
	void coord2paramaddlocal(vector<Residueatoms> &decstr, int trandpos, int moveatom,int start,int end);
	bool fmovementlenlocal(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc);

	//Movement of bond length with Metropolis method
	bool fmovementlenmt(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc);
	void setlenntmtflap(double numseq, int trandpos, double aat, vector<int> targetlengths);

	//Movement of angles and bond length
	void setmultianglenflap(double numseq, double trandpos, double aat, int span);
	bool fmovementlenang(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt, bool flagacc);

	//Movement of multiple torsion angles
	double getmovmultitorval(double val, double mean, double k);
	bool fmovementmultitor(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt, bool flagacc);
	void setmultitorflap(double numseq, double trandpos, double aat);

	//Movement of multiple sugar torsion angles
	void setmultitorsugarflap(double numseq, double trandpos, double aat);
	bool fmovementmultitorsugar(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt, bool flagacc);

	//Movement of geometries
	void setgeoflap(double numseq, int trandpos, int paramtype, double aat);
	void setgeontflap(double numseq, int trandpos, int paramtype, double aat);
	bool fmovementgeo(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc);

	//Movement of hydrogen geometries
	void setgeohydflap(double numseq, int trandpos, int paramtype, double aat);
	bool fmovementhydgeo(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc);

	//Movement of hydrogen geometries based on gradients
	void gethbhydanglemse(Point3s P, Point3s &Q, Point3s R, double &angles, double angleideal, double &mse, double &theta, double &diff, double &normPQ, double &normRQ, double &dotp, Point3s &vecPQ, Point3s &vecRQ, double coeff);
	void gethbhydlengthmse(Point3s P, Point3s &Q, double &lengths, double lengthideal, double &mse, double &diff, Point3s vecPQ, double coeff);
	void initializegrad(Point3s &gradQ, double &mse);
	void calchbhydanglemsegradient(Point3s P, Point3s &Q, Point3s R, double &angles, double angleideal, Point3s& gradQ, double theta, double diff, double normPQ, double normRQ, double dotp, Point3s &vecPQ, Point3s &vecRQ, double coeff);
	void calchbhydlengthmsegradient(Point3s P, Point3s &Q, double &lengths, double lengthideal, Point3s &gradQ, double &diff, Point3s vecPQ, double coeff);
	void applyhbhydSGD(Point3s P, Point3s &Q, Point3s R, double &angles, double &lengths1, double &lengths2, double angleideal, double lengthideal1, double lengthideal2, double learning_rate, int epochs, double tolerance);
	bool fmovementgradhydgeo(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc);
	void calcmseangles(vector<Residueatoms> restr, int targetnt, vector<double> &angle, vector<double> &diff, vector<bool> &outerflags, vector<double> &theta, vector<double> &normPQs, vector<double> &normRQs, vector<double> &dotps, vector<Point3s> &vecPQs, vector<Point3s> &vecRQs, double &mse, int &count, double coeff);
	void calcmseangle(Point3s P, Point3s Q, Point3s R, double &angle, double angleideal, double anglerangelower, double anglerangeupper, bool &outflag, double &mse, double &theta, double &diff, double &normPQ, double &normRQ, double &dotp, Point3s &vecPQ, Point3s &vecRQ, double coeff);
	void calcmsefullangles(vector<Residueatoms> restr, int targetnt, vector<double> &angle, vector<double> &diff, vector<bool> &outerflags, vector<double> &theta, vector<double> &normPQs, vector<double> &normRQs, vector<double> &dotps, vector<Point3s> &vecPQs, vector<Point3s> &vecRQs, double &mse, int &count, double coeff);
	void calcgradientfullangles(vector<Residueatoms> restr, int targetnt, map<int, vector<Point3s>> &gradptv, vector<bool> movenodeinfo, vector<double> diff, vector<bool> outerflags, vector<double> theta, vector<double> normPQs, vector<double> normRQs, vector<double> dotps, vector<Point3s> vecPQs, vector<Point3s> vecRQs, int count, double coeff);

	//Movement rotate
	bool fmovementrot(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc);
	void settmprotateflap(double numseq, double trandpos1, double trandpos2, char aat1, char aat2);

	//Movement of torsions in sugar rings
	bool fmovementsugartors(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt, bool flagacc);
	void setsugartorsflap(double numseq, int trandpos1, int trandpos2);

	//Movement LMP
	//void singlemoveLMP(vector<Residueatoms> &restr,int k,int flagdir, int ntidx, int adjidx);
	//void setlmpparam(vector<Residueatoms> &restr, int trandposst, int trandposed, vector<int> &ntparamidx, vector<int> &adjparamidx);
	//void settmplmpflap(double numseq, double trandpos1, double trandpos2);
	//bool fmovementlmp(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc);
	void singlemoveLMP(vector<Residueatoms> &restr,int k,int flagdir);
	void settmplmpflap(double numseq, double trandpos1, double trandpos2);
	bool fmovementlmp(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc);

	void singlemoveLMPfull(vector<Residueatoms> &restr,int k,int flagdir,int paramidx);
	void settmplmpfullflap(double numseq, double trandpos1, double trandpos2);
	bool fmovementlmpfull(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc);
	void setlmpparam(vector<Residueatoms> &restr, int trandposst, int trandposed, vector<int> &ntparamidx, vector<int> &adjparamidx);

	//Movement LMP for bp
	void singbplemoveLMPfull(vector<Residueatoms> &restr,int k,int rtype,int flagdir,int paramidx);
	void singbplemoveLMPfullbp(vector<Residueatoms> &restr,int i, int j,int rtype);
	void settmpbplmpfullflap(double numseq, double trandpos1, double trandpos2);
	void setlmpbpparam(vector<Residueatoms> &restr, vector<pair<int,int>> pairs, vector<vector<int>> &lmpparamidx);
	bool checkadjlength(vector<Residueatoms> restr, int basei,int paramidx);
	bool checkntlength(vector<Residueatoms> restr, int basei,int paramidx);
	bool checkbplength(vector<Residueatoms> restr, int basei, int basej);
	bool fmovementlmpfullbasepair(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc);

	//Movement gradient-based LMP for bp
	void getgraphinfoLMP(vector<pair<int,int>> pairs, int numseq, vector<vector<int>> lmpparamidx, vector<bool> &movenodeinfo, vector<pair<int, int>> &edgeinfo, vector<int> &edgetype, vector<int> &paraminfo);
	void settmpbpgradientlmpfullflap(double numseq, double trandpos1, double trandpos2);
	void setgradientlmpbpparam(vector<Residueatoms> &restr, vector<pair<int,int>> pairs, vector<vector<int>> &lmpparamidx);
	bool fmovementgradientlength(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc);
	bool fmovementgradientangle(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc);
	void calcmsentlength(Point3s *P, Point3s *Q, char aat,int paramidx, vector<vector<double>> &length, vector<vector<double>> &diff, vector<vector<Point3s>> &dists, double &mse, int &count, vector<vector<bool>> &outlineridx, double coeff, bool showres);
	void calcmseadjlength(Point3s *P, Point3s *Q, int basei, int basej, char aat1, char aat2, int inda, int indb, int paramidx, vector<vector<double>> &length, vector<vector<double>> &diff, vector<vector<Point3s>> &dists, double &mse, int &count, vector<vector<bool>> &outlineridx, double coeff, bool showres);
	void calcmsebplength(Point3s *P, Point3s *Q, char aat1, char aat2, vector<vector<double>> &length, vector<vector<double>> &diff, vector<vector<Point3s>> &dists, double &mse, int &count, vector<vector<bool>> &outlineridx, double coeff, bool showres);
	void initializegrad(map<int, vector<Point3s>> &gradptv, vector<bool> movenodeinfo, double &mse);
	void calcgradientlength(Point3s *P, Point3s *Q, int basei, int basej, map<int, vector<Point3s>> &gradptv, vector<bool> movenodeinfo, vector<vector<double>> lengths, vector<vector<double>> diff, vector<vector<Point3s>> dists, int count, vector<vector<bool>> &outlineridx, double coeff);
	void calcmseinitlength(Point3s *P, Point3s *Q, vector<double> &length, vector<Point3s> &dists, double &mse, int &count, double coeff);
	void calcgradientinitlength(Point3s *P, Point3s *Q, int basei, map<int, vector<Point3s>> &gradptv, vector<double> lengths, vector<Point3s> dists, int count, double coeff);
	void calcgradientangle(Point3s P, Point3s Q, Point3s R, map<int, vector<Point3s>> &gradptv, int basep, int baseq, int baser, int mttypep, int mttypeq, int mttyper, bool &flagP, bool &flagQ, bool &flagR, double theta, double diff, double normPQ, double normRQ, double dotp, Point3s vecPQ, Point3s vecRQ, double count, double coeff);
	void calcgradientangles(vector<Residueatoms> restr, int targetnt, map<int, vector<Point3s>> &gradptv, vector<bool> movenodeinfo, vector<double> diff, vector<bool> outerflags, vector<double> theta, vector<double> normPQs, vector<double> normRQs, vector<double> dotps, vector<Point3s> vecPQs, vector<Point3s> vecRQs, int count, double coeff);
	void selectsnttemplate(Residueatoms restr, Residueatoms &snttemplate);
	void calcmsetmplength(Point3s *P, Point3s *Q, vector<double> &length, vector<Point3s> &dists, double &mse, int &count, double coeff);
	void calcgradienttmplength(Point3s *P, Point3s *Q, int basei, map<int, vector<Point3s>> &gradptv, vector<double> lengths, vector<Point3s> dists, int count, double coeff);

	void calcmsetempadj(Point3s *restr1, Point3s *restr2, char aat1, char aat2, vector<Residueatoms> strtemplate, vector<double> &length, vector<Point3s> &dists, double &mse, int &count, double coeff);
	void calcgradienttempadj(int basei, int basej, char aat1, char aat2, map<int, vector<Point3s>> &gradptv, vector<double> lengths, vector<Point3s> dists, int count, double coeff);

	void calcmsetempbp(Point3s *restr1, Point3s *restr2, char aat1, char aat2, vector<Residueatoms> strtemplate, vector<double> &length, vector<Point3s> &dists, double &mse, int &count, double coeff);
	void calcgradienttempbp(int basei, int basej, char aat1, char aat2, map<int, vector<Point3s>> &gradptv, vector<double> lengths, vector<Point3s> dists, int count, double coeff);

	double getgraddotproduct(map<int, vector<Point3s>> grad, map<int, vector<Point3s>> z, vector<bool> movenodeinfo);
	void strupdate(vector<Residueatoms> &restr, double step_size, map<int, vector<Point3s>> z, vector<bool> movenodeinfo);

	void selectdnttemplate(Residueatoms restr1, Residueatoms restr2, int basei, int basej, vector<Residueatoms> &dnttemplate, TemplateDistMap &dnttemplatelength);
	void selectbptemplate(Residueatoms restr1, Residueatoms restr2, vector<Residueatoms> &dnttemplate, TemplateDistMap &bptemplatelength);
	bool fmovementgradientlength2(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc);
	//bool fmovementgradientlength3(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc);
	//void calcmseditmplength(Point3s *P1, Point3s *Q1, Point3s *P2, Point3s *Q2, vector<double> &length, vector<vector<double>> &diff, vector<Point3s> &dists, double &mse, int &count, double coeff);
	void calcmseditmplength(Point3s *P1, Point3s *Q1, Point3s *P2, Point3s *Q2, TemplateDistMap templateslength, vector<vector<double>> &length, vector<vector<double>> &diff, vector<vector<Point3s>> &dists, double &mse, int &count, double coeff);
	//void calcgradientditmplength(int basei, int basej, map<int, vector<Point3s>> &gradptv, vector<double> lengths, vector<Point3s> dists, int count, double coeff);
	void calcgradientditmplength(int basei, int basej, map<int, vector<Point3s>> &gradptv, vector<bool> movenodeinfo, vector<vector<double>> lengths, vector<vector<double>> diff, vector<vector<Point3s>> dists, int count, double coeff);
	void getgraphinfoLMP2(vector<pair<int,int>> pairs, int numseq, vector<vector<int>> lmpparamidx, vector<bool> &movenodeinfo, vector<pair<int, int>> &edgeinfo, vector<int> &edgetype, vector<int> &paraminfo);

	void calcene2(vector<Residueatoms> restr, vector<Residueatoms> decstr, Gradmain2 &gradobj, vector<pair<int, int>> edgeinfo, vector<int> edgetype, vector<int> paraminfo, vector<int> movementresidues, map<int, vector<Residueatoms>> dnttemplates, map<int, vector<Residueatoms>> bptemplates, map<int, TemplateDistMap> dnttemplateslength, map<int, TemplateDistMap> bptemplateslength, map<string, double> enecoeff);
	//void calcene3(vector<Residueatoms> restr, vector<Residueatoms> decstr, Gradmain2 &gradobj, vector<pair<int, int>> edgeinfo, vector<int> edgetype, vector<int> paraminfo, vector<int> movementresidues, map<int, vector<Residueatoms>> dnttemplates, map<int, vector<Residueatoms>> bptemplates, map<int, TemplateDistMap> dnttemplateslength, map<int, TemplateDistMap> bptemplateslength, map<string, double> enecoeff);
	void calcgrad2(vector<Residueatoms> restr, vector<Residueatoms> decstr, map<int, vector<Point3s>> &gradptv, Gradmain2 gradobj, vector<pair<int, int>> edgeinfo, vector<int> edgetype, vector<bool> movenodeinfo, vector<int> movementresidues, map<string, double> enecoeff);
	//void calcgrad3(vector<Residueatoms> restr, vector<Residueatoms> decstr, map<int, vector<Point3s>> &gradptv, Gradmain2 gradobj , vector<pair<int, int>> edgeinfo, vector<int> edgetype, vector<bool> movenodeinfo, vector<int> movementresidues, map<string, double> enecoeff);
	void initgrad2obj(Gradmain2 &gradobj, int index);
	double lbfgsgradnorm(map<int, vector<Point3s>> grad, vector<bool> movenodeinfo);
	void lbfgsupdate(map<int, vector<Point3s>> grad, deque<map<int, vector<Point3s>>> s_list, deque<map<int, vector<Point3s>>> y_list, deque<double> rho_list, vector<bool> movenodeinfo, map<int, vector<Point3s>> &z);
	void initgrad(map<int, vector<Point3s>> &gradptv, vector<bool> movenodeinfo);
	double linesearch(vector<Residueatoms>& restr, vector<Residueatoms> decstr, map<int, vector<Point3s>>& grad, map<int, vector<Point3s>> z, Gradmain2 &gradobj, vector<pair<int, int>> edgeinfo, vector<int> edgetype, vector<int> paraminfo, vector<bool> movenodeinfo, vector<int> movementresidues, map<int, vector<Residueatoms>> dnttemplates, map<int, vector<Residueatoms>> bptemplates, map<int, TemplateDistMap> dnttemplateslength, map<int, TemplateDistMap> bptemplateslength, map<string, double> enecoeff, double alpha, double beta, double gamma);
	//double linesearch3(vector<Residueatoms>& restr, vector<Residueatoms> decstr, map<int, vector<Point3s>>& grad, map<int, vector<Point3s>> z, Gradmain2 &gradobj, vector<pair<int, int>> edgeinfo, vector<int> edgetype, vector<int> paraminfo, vector<bool> movenodeinfo, vector<int> movementresidues, map<int, vector<Residueatoms>> dnttemplates, map<int, vector<Residueatoms>> bptemplates, map<int, TemplateDistMap> dnttemplateslength, map<int, TemplateDistMap> bptemplateslength, map<string, double> enecoeff, double alpha, double beta, double gamma);
	void calcmseditmpintralength(Point3s *P, vector<vector<double>> templatesintralength, vector<vector<double>> &length, vector<vector<double>> &diff, vector<vector<Point3s>> &dists, double &mse, int &count, double coeff);
	void calcgradientditmplength(int basei, map<int, vector<Point3s>> &gradptv, vector<bool> movenodeinfo, vector<vector<double>> lengths, vector<vector<double>> diff, vector<vector<Point3s>> dists, int count, double coeff);

	void initgradgeo(vector<Residueatoms> &recrrstr, vector<bool> movenodeinfo, int numseq);
	void initrandommove(vector<Residueatoms> &recrrstr, vector<bool> movenodeinfo, int numseq);

	//lgm
	bool fmovementlgm(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc);
	double distfromene(vector<Residueatoms> restr, vector<Residueatoms> decstr, int st, int ed, int type);
	void initgeo(vector<Residueatoms> &restr, int st, int ed, int type);
	void selectbptemplate(Residueatoms restr1, Residueatoms restr2, vector<Residueatoms> &dnttemplate);
	double goldenSearchlgm(vector<Residueatoms> restr, vector<Residueatoms> decstr, int trandpos, int dihtype, int trandposst, int trandposed, double a, double b, double tolerance = 1e-6);
	double goldenSearchcalcene(vector<Residueatoms> restr, vector<Residueatoms> decstr, double torval, int dihtype, int trandpos, int trandposst, int trandposed);
	bool checkoutlinerlgm(double val, double lower, double upper);

	//lgm2
	bool fmovementlgm2(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc);
	double goldenSearchlgm2(vector<Residueatoms> restr, vector<Residueatoms> decstr, map<int, Point3s> initpoints1, map<int, Point3s> initpoints2, vector<pair<int,int>> bpairs, map<int, TemplateDistMap> bptemplateslength, map<int, TemplateDistMap> dnttemplateslength, int trandpos, int dihtype, int trandposst, int trandposed, double a, double b, double tolerance, int numseq);
	double goldenSearchcalcene2(vector<Residueatoms> restr, vector<Residueatoms> decstr, map<int, Point3s> initpoints1, map<int, Point3s> initpoints2, vector<pair<int,int>> bpairs, map<int, TemplateDistMap> bptemplateslength, map<int, TemplateDistMap> dnttemplateslength, double torval, int dihtype, int trandpos, int trandposst, int trandposed, int numseq);
	double distfromene2(vector<Residueatoms> recrrstr, vector<Residueatoms> decstr, map<int, Point3s> initpoints1, map<int, Point3s> initpoints2, vector<pair<int,int>> bpairs, map<int, TemplateDistMap> bptemplateslength, map<int, TemplateDistMap> dnttemplateslength, int st, int ed, int numseq);
	void getcenterpoints(vector<Residueatoms> decstr, map<int, Point3s> &centerpoints, int st, int ed, vector<int> atomlist);
	void initgeo2(vector<Residueatoms> &restr, int st, int ed, int type);
	double eneinitpoints(vector<Residueatoms> recrrstr, map<int, Point3s> initpoints1, map<int, Point3s> initpoints2, int st, int ed);
	void eneedgegeoparams(vector<Residueatoms> recrrstr, int ed, int numseq, double &feneangle, double &fenelength);
	void eneaddgeoparams(vector<Residueatoms> recrrstr, int st, int ed, double &feneangle, double &fenelength);
	void selectbptemplate2(Residueatoms restr1, Residueatoms restr2, vector<Residueatoms> &dnttemplate, TemplateDistMap &bptemplatelength);
	void calbplength(Residueatoms res1, Residueatoms res2, vector<vector<double>> &interdistmap);
	double compbplengthmap(vector<vector<double>> map1, vector<vector<double>> map2);
	double enebptemplatelength(vector<Residueatoms> recrrstr, vector<pair<int,int>> bpairs, map<int, TemplateDistMap> bptemplateslength);
	double enebpplane(vector<Residueatoms> decstr, int basei, int basej);
	double enebptor(vector<Residueatoms> decstr, int basei, int basej);
	void selectdnttemplate2(Residueatoms restr1, Residueatoms restr2, int basei, int basej, vector<Residueatoms> &dnttemplate, TemplateDistMap &dnttemplatelength);
	double eneadjtemplatelength(vector<Residueatoms> recrrstr, int st, int ed, int numseq, map<int, TemplateDistMap> dnttemplateslength);

	void initrandommovegss(Residueatoms &recrrstr);
	void resetflaggss(vector<bool> &residueflag, vector<vector<bool>> &atomflag, int numseq);
	double calcinitstrgss(vector<Residueatoms> restr, vector<Residueatoms> decstr, vector<int> nodes, vector<double> &energy, vector<bool> residueflag, int uptype);
	void initenergygss(map<string, vector<double>> &edgeenergy, map<string, vector<double>> &nodeenergy, vector<vector<bool>> &lengthoutflag, vector<vector<bool>> &angleoutflag, vector<pair<int, int>> edges, vector<int> nodes, int numseq);
	double calcenegss(vector<Residueatoms> restr, vector<Residueatoms> decstr, int numseq, vector<pair<int, int>> edgeinfo, vector<int> edgetype, vector<int> paraminfo, vector<int> nodes, map<int, vector<Residueatoms>> dnttemplates, map<int, vector<Residueatoms>> bptemplates, map<int, TemplateDistMap> dnttemplateslength, map<int, TemplateDistMap> bptemplateslength, map<string, double> enecoeff, vector<bool> residueflag, vector<vector<bool>> atomflag, map<string, vector<double>> &edgeenergy, map<string, vector<double>> &nodeenergy, vector<vector<bool>> &lengthoutflag, vector<vector<bool>> &angleoutflag, int uptype);
	void getflaggss(vector<bool> &residueflag, vector<vector<bool>> &atomflag, int trandpos, char aat, int paramtype);
	void changeparamgss(vector<Residueatoms> &recrrstr, int paramtype, int trandpos);
	void updateparamgss(vector<Residueatoms> &recrrstr, vector<vector<bool>> atomflag, int trandpos, int numseq);
	double goldenSearchcalcene(vector<Residueatoms> recrrstr, vector<Residueatoms> decstr, double torval, int paramtype, int paramidx, int trandpos, int numseq, vector<pair<int, int>> edgeinfo, vector<int> edgetype, vector<int> paraminfo, vector<int> nodes, map<int, vector<Residueatoms>> dnttemplates, map<int, vector<Residueatoms>> bptemplates, map<int, TemplateDistMap> dnttemplateslength, map<int, TemplateDistMap> bptemplateslength, map<string, double> enecoeff, vector<bool> residueflag, vector<vector<bool>> atomflag, map<string, vector<double>> edgeenergy, map<string, vector<double>> nodeenergy, vector<vector<bool>> lengthoutflag, vector<vector<bool>> angleoutflag);
	double goldenSearchdih(vector<Residueatoms> &recrrstr, vector<Residueatoms> decstr, int trandpos, int dihtype, double tolerance, int numseq, vector<pair<int, int>> edgeinfo, vector<int> edgetype, vector<int> paraminfo, vector<int> nodes, map<int, vector<Residueatoms>> dnttemplates, map<int, vector<Residueatoms>> bptemplates, map<int, TemplateDistMap> dnttemplateslength, map<int, TemplateDistMap> bptemplateslength, map<string, double> enecoeff, map<string, vector<double>> &edgeenergy, map<string, vector<double>> &nodeenergy, vector<vector<bool>> &lengthoutflag, vector<vector<bool>> &angleoutflag);
	double goldenSearchang(vector<Residueatoms> &recrrstr, vector<Residueatoms> decstr, int trandpos, int angtype, double tolerance, int numseq, vector<pair<int, int>> edgeinfo, vector<int> edgetype, vector<int> paraminfo, vector<int> nodes, map<int, vector<Residueatoms>> dnttemplates, map<int, vector<Residueatoms>> bptemplates, map<int, TemplateDistMap> dnttemplateslength, map<int, TemplateDistMap> bptemplateslength, map<string, double> enecoeff, map<string, vector<double>> &edgeenergy, map<string, vector<double>> &nodeenergy, vector<vector<bool>> &lengthoutflag, vector<vector<bool>> &angleoutflag);
	double goldenSearchlen(vector<Residueatoms> &recrrstr, vector<Residueatoms> decstr, int trandpos, int lentype, double tolerance, int numseq, vector<pair<int, int>> edgeinfo, vector<int> edgetype, vector<int> paraminfo, vector<int> nodes, map<int, vector<Residueatoms>> dnttemplates, map<int, vector<Residueatoms>> bptemplates, map<int, TemplateDistMap> dnttemplateslength, map<int, TemplateDistMap> bptemplateslength, map<string, double> enecoeff, map<string, vector<double>> &edgeenergy, map<string, vector<double>> &nodeenergy, vector<vector<bool>> &lengthoutflag, vector<vector<bool>> &angleoutflag);
	double calcmseditmplengthgss(vector<Residueatoms> restr, vector<pair<int, int>> edgeinfo, vector<int> edgetype, map<int, TemplateDistMap> templateslength, vector<double> &energy, vector<bool> residueflag, int uptype);
	double getouteneanglegss(double val, double mean, double lower, double upper, double slower, double supper, bool &outflag);
	double getoutenelength(double val, double mean, double lower, double upper, double slower, double supper, bool &outflag);
	double calcenebondanglegss(vector<Residueatoms> restr, vector<int> nodes, vector<double> &energy, vector<vector<bool>> &angleoutflag, vector<bool> residueflag, int numseq, int uptype);
	double calcenebondlengthgss(vector<Residueatoms> restr, vector<int> nodes, vector<double> &energy, vector<vector<bool>> &lengthoutflag, vector<bool> residueflag, int numseq, int uptype);
	double calcmsebptmplengthgss(vector<Residueatoms> restr, vector<pair<int, int>> edgeinfo, vector<int> edgetype, map<int, TemplateDistMap> templateslength, vector<double> &energy, vector<bool> residueflag, int uptype);
	double calcmsesuptempgss(Point3s *restr1, Point3s *restr2, char aat1, char aat2, vector<Residueatoms> strtemplate);
	double calcmseditmpgss(vector<Residueatoms> restr, vector<pair<int, int>> edgeinfo, vector<int> edgetype, map<int, vector<Residueatoms>> dnttemplates, vector<double> &energy, vector<bool> residueflag, int uptype);
	bool getmovetargetlengthgss(vector<vector<bool>> lengthoutflag, vector<int> nodes, vector<int> mnodes, int &targetpos, int &targetparam);
	bool getmovetargetanglegss(vector<vector<bool>> angoutflag, vector<int> nodes, vector<int> mnodes, int &targetpos, int &targetparam);

	void initopleng(vector<Residueatoms> &recrrstr, vector<int> nodes);
	//void calcmsecnnditmpintralength(vector<Residueatoms> restr, int base, vector<vector<double>> intralength, vector<vector<double>> interlength, vector<vector<double>> &length, vector<vector<double>> &diff, vector<vector<Point3s>> &dists, double &mse, int &count, double coeff);
	void calcmsecnnditmpintralength(vector<Residueatoms> restr, int base, vector<vector<double>> &length, vector<vector<double>> &diff, vector<vector<Point3s>> &dists, vector<vector<bool>> &outlineridx, double &mse, int &count, double coeff);
	//void calcgracnnentditmpintra(int base, map<int, vector<Point3s>> &gradptv, vector<bool> movenodeinfo, vector<vector<double>> lengths, vector<vector<double>> diff, vector<vector<Point3s>> dists, int count, double coeff);
	void calcgracnnentditmpintra(int base, map<int, vector<Point3s>> &gradptv, vector<bool> movenodeinfo, vector<vector<double>> lengths, vector<vector<double>> diff, vector<vector<Point3s>> dists, vector<vector<bool>> outlineridx , int count, double coeff);
	

	double calcenenbondlengthgss(vector<Residueatoms> restr, vector<int> nodes, vector<double> &energy, vector<bool> residueflag, int uptype);

	//other functions for movements
	bool metropolisUpdate(double oldenergy, double newenergy, double Bt);
	int getrandval(vector<double> tmov,double trandnum);
	void copystr(vector<Residueatoms> decstr, int numseq);
	void rcopystr(vector<Residueatoms> &decstr, vector<Residueatoms> restr, int numseq);
	void fsetinitial(Energyterms &lp, vector<vector<Rescontact>> &rpp, vector<Basecontact> &bprpp, vector<Stackcontact> &bsrpp, vector<vector<int>> &nmaprpp, vector<vector<vector<double>>> &distmaprpp, vector<vector<vector<double>>> &ntpairdistrpp, vector<vector<vector<double>>> &adjpairdistrpp, vector<double> &distenerpp, vector<vector<double>> &multidistenerpp, vector<Phycontact> &contactenerpp, vector<vector<bool>> &mapflagrpp, vector<vector<bool>> &mapbbflagrpp, int numseq);
	
	double getrandomgradparam(double crrp, double l, double r, vector<int> distr, vector<double> boundary, vector<double> energy, double span);

	RNASimulation(const RNASimulation& other) {
        // 配列のコピー
        std::copy(std::begin(other.enelist), std::end(other.enelist), enelist);
        std::copy(std::begin(other.enelistbk), std::end(other.enelistbk), enelistbk);
        
        // 基本データ型のコピー
        contwt = other.contwt;
        curcycle = other.curcycle;
        baseenecoeff = other.baseenecoeff;
        enetype = other.enetype;
		bpmap = other.bpmap;
        initfape = other.initfape;
        fapecoff = other.fapecoff;
        sseneratio[0] = other.sseneratio[0];
        sseneratio[1] = other.sseneratio[1];
        tempnum = other.tempnum;
        numdim = other.numdim;
        initlentmp = other.initlentmp;
		coefftemplate = other.coefftemplate;

        // std::vector のコピー
        dihevec = other.dihevec;
        lpedpos = other.lpedpos;
        //basepairs = other.basepairs;
        //stackingpairs = other.stackingpairs;
        rppmat = other.rppmat;
        rppmatbk = other.rppmatbk;
        nmap = other.nmap;
        nmapbk = other.nmapbk;
        restr = other.restr;
        distmap = other.distmap;
        distmapbk = other.distmapbk;
		ntpairdist = other.ntpairdist;
		ntpairdistbk = other.ntpairdistbk;
		adjpairdist = other.adjpairdist;
		adjpairdistbk = other.adjpairdistbk;
        distene = other.distene;
		contactene = other.contactene;
		contactenebk = other.contactenebk;
        distenebk = other.distenebk;
        multidistene = other.multidistene;
        multidistenebk = other.multidistenebk;
        mapflag = other.mapflag;
        mapflagbk = other.mapflagbk;
		mapbbflag = other.mapbbflag;
        mapbbflagbk = other.mapbbflagbk;
        bpmapid2idx = other.bpmapid2idx;
        bsmapid2idx = other.bsmapid2idx;
		bppairednt = other.bppairednt;
		stpairednt = other.stpairednt;
        bpenergy = other.bpenergy;
        bpenergybk = other.bpenergybk;
        bsenergy = other.bsenergy;
        bsenergybk = other.bsenergybk;

        // その他のメンバーのコピー
        lp1 = other.lp1;
        lp2 = other.lp2;
        initpaa3 = other.initpaa3;
        paa3 = other.paa3;
        npaa3 = other.npaa3;
        initnpaa3 = other.initnpaa3;
        lentmp = other.lentmp;
    }

	virtual ~RNASimulation();

};

