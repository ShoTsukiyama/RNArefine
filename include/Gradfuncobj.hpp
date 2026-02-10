
#ifndef GRADFUNCOBJ_HPP
#define GRADFUNCOBJ_HPP 1

#pragma once
#include <list>
#include <vector>
#include <string>
#include "CommonParams.hpp"
#include "BasicFunc.hpp"
#include "PDBParser.hpp"

using namespace std;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

typedef struct Residueatoms
{
	int residueid;//sequence
	float dih[PARAM_NUM],len[PARAM_NUM],ang[PARAM_NUM];
	float addlen[ADDLENGTH_NUM], addang[ADDINNERANGLE_NUM];
	float hyddih[HYD_PARAM_NUM],hydlen[HYD_PARAM_NUM],hydang[HYD_PARAM_NUM];
    float addhydlen[HYD_ADDPARAM_NUM],addhydang[HYD_ADDANGLE_NUM];
	Point3s ptv[ATOM_RES_NUM];//heavily atomの座標の情報
	
	char aaa;//残基名
	int resind; //残基の種類のid
	int residx;

	Residueatoms() 
        : residueid(0), 
          aaa('\0'), 
          resind(0), 
          residx(0)
    {
        // 配列の初期化
        fill(dih, dih + PARAM_NUM, 0.0f);
        fill(len, len + PARAM_NUM, 0.0f);
        fill(ang, ang + PARAM_NUM, 0.0f);
		fill(addlen, addlen + ADDLENGTH_NUM, 0.0f);
		fill(addang, addang + ADDINNERANGLE_NUM, 0.0f);
		fill(hyddih, hyddih + HYD_PARAM_NUM, 0.0f);
        fill(hydlen, hydlen + HYD_PARAM_NUM, 0.0f);
        fill(hydang, hydang + HYD_PARAM_NUM, 0.0f);
        fill(addhydlen, addhydlen + HYD_ADDPARAM_NUM, 0.0f);
        fill(addhydang, addhydang + HYD_ADDANGLE_NUM, 0.0f);
        fill(ptv, ptv + ATOM_RES_NUM, Point3s()); // Point3sのデフォルトコンストラクタを呼び出す
    }

	Residueatoms(const Residueatoms& other) 
        : residueid(other.residueid), 
          aaa(other.aaa), 
          resind(other.resind), 
          residx(other.residx)
    {
        copy(other.dih, other.dih + PARAM_NUM, dih);
        copy(other.len, other.len + PARAM_NUM, len);
        copy(other.ang, other.ang + PARAM_NUM, ang);
		copy(other.addlen, other.addlen + ADDLENGTH_NUM, addlen);
        copy(other.addang, other.addang + ADDINNERANGLE_NUM, addang);
		copy(other.hyddih, other.hyddih + HYD_PARAM_NUM, hyddih);
        copy(other.hydlen, other.hydlen + HYD_PARAM_NUM, hydlen);
        copy(other.hydang, other.hydang + HYD_PARAM_NUM, hydang);
        copy(other.addhydlen, other.addhydlen + HYD_ADDPARAM_NUM, addhydlen);
        copy(other.addhydang, other.addhydang + HYD_ADDANGLE_NUM, addhydang);
        copy(other.ptv, other.ptv + ATOM_RES_NUM, ptv);
    }
}Residueatoms;

class InitPairobj
{
  private:
  	Point3s* teAtomA;
	Point3s* teAtomB;
	Point3s* teAtomAgrad;
	Point3s* teAtomBgrad;
    double ideal;
	double KI;

  public:
  InitPairobj(Point3s* _teAtomA, Point3s* _teAtomB, Point3s* _teAtomAgrad, Point3s* _teAtomBgrad, double _ideal, double _KI)
  		: teAtomA(_teAtomA), teAtomB(_teAtomB), teAtomAgrad(_teAtomAgrad), teAtomBgrad(_teAtomBgrad), ideal(_ideal), KI(_KI) {}
    
    double calcenergyinitpair() const;
    void   calcgradientinitpair() const;

	static InitPairobj create(Point3s* _teAtomA, Point3s* _teAtomB, Point3s* _teAtomAgrad, Point3s* _teAtomBgrad, double _ideal, double _KI) {
        return InitPairobj(_teAtomA, _teAtomB, _teAtomAgrad, _teAtomBgrad, _ideal, _KI);
    }
};

class InitPosobj
{
  private:
  	Point3s* teAtomA;
	Point3s* teAtomAgrad;
	const double X0, Y0, Z0, Kr;

  public:
  InitPosobj(Point3s* _teAtomA, Point3s* _teAtomAgrad, double _X0, double _Y0, double _Z0, double _Kr)
  		: teAtomA(_teAtomA), teAtomAgrad(_teAtomAgrad), X0(_X0), Y0(_Y0), Z0(_Z0), Kr(_Kr) {}
    
    double calcenergyinit() const;
    void   calcgradientinit() const;

	static InitPosobj create(Point3s* _teAtomA, Point3s* _teAtomAgrad, double X0, double Y0, double Z0, double Kr) {
        return InitPosobj(_teAtomA, _teAtomAgrad, X0, Y0, Z0, Kr);
    }
};

class TElectrobj
{
  private:
  	Point3s* teAtomA;
	Point3s* teAtomB;
	Point3s* teAtomAgrad;
	Point3s* teAtomBgrad;
    double dQQep;

  public:
  	TElectrobj(Point3s* _teAtomA, Point3s* _teAtomB, Point3s* _teAtomAgrad, Point3s* _teAtomBgrad, double _dQQep =-1.0)
  		: teAtomA(_teAtomA), teAtomB(_teAtomB), teAtomAgrad(_teAtomAgrad), teAtomBgrad(_teAtomBgrad), dQQep(_dQQep) {}
    
    double calcenergyelectros() const;
    void   calcgradientelectros() const;

	static TElectrobj create(Point3s* _teAtomA, Point3s* _teAtomB, Point3s* _teAtomAgrad, Point3s* _teAtomBgrad, double _dQQep = -1.0) {
        return TElectrobj(_teAtomA, _teAtomB, _teAtomAgrad, _teAtomBgrad, _dQQep);
    }
};

class LennardsJobj
{
  private:
  	Point3s* teAtomA;
	Point3s* teAtomB;
	Point3s* teAtomAgrad;
	Point3s* teAtomBgrad;
    double LJ6, LJ12;
	double K;

  public:
  LennardsJobj(Point3s* _teAtomA, Point3s* _teAtomB, Point3s* _teAtomAgrad, Point3s* _teAtomBgrad, double _LJ6, double _LJ12, double _K)
  		: teAtomA(_teAtomA), teAtomB(_teAtomB), teAtomAgrad(_teAtomAgrad), teAtomBgrad(_teAtomBgrad), LJ6(_LJ6), LJ12(_LJ12), K(_K) {}
    
    double calcenergylj() const;
    void   calcgradientlj() const;

	static LennardsJobj create(Point3s* _teAtomA, Point3s* _teAtomB, Point3s* _teAtomAgrad, Point3s* _teAtomBgrad, double _LJ6, double _LJ12, double _K) {
        return LennardsJobj(_teAtomA, _teAtomB, _teAtomAgrad, _teAtomBgrad, _LJ6, _LJ12, _K);
    }
};

class BondLengthobj
{
  private:
  	Point3s* teAtomA;
	Point3s* teAtomB;
	Point3s* teAtomAgrad;
	Point3s* teAtomBgrad;
    double ideal, stdval;
	double KI;

  public:
  BondLengthobj(Point3s* _teAtomA, Point3s* _teAtomB, Point3s* _teAtomAgrad, Point3s* _teAtomBgrad, double _ideal, double _stdval, double _KI)
  		: teAtomA(_teAtomA), teAtomB(_teAtomB), teAtomAgrad(_teAtomAgrad), teAtomBgrad(_teAtomBgrad), ideal(_ideal), stdval(_stdval), KI(_KI) {}
    
    double calcenergylength() const;
    void   calcgradientlength() const;

	static BondLengthobj create(Point3s* _teAtomA, Point3s* _teAtomB, Point3s* _teAtomAgrad, Point3s* _teAtomBgrad, double _ideal, double _stdval, double _KI) {
        return BondLengthobj(_teAtomA, _teAtomB, _teAtomAgrad, _teAtomBgrad, _ideal, _stdval, _KI);
    }
};

class BondAngleobj
{
  private:
  	Point3s* teAtomA;
	Point3s* teAtomB;
	Point3s* teAtomC;
	Point3s* teAtomAgrad;
	Point3s* teAtomBgrad;
	Point3s* teAtomCgrad;
    double ideal, stdval, KI;

  public:
  	BondAngleobj(Point3s* _teAtomA, Point3s* _teAtomB, Point3s* _teAtomC, Point3s* _teAtomAgrad, Point3s* _teAtomBgrad, Point3s* _teAtomCgrad, double _ideal, double _stdval, double _KI)
  		: teAtomA(_teAtomA), teAtomB(_teAtomB), teAtomC(_teAtomC), teAtomAgrad(_teAtomAgrad), teAtomBgrad(_teAtomBgrad), teAtomCgrad(_teAtomCgrad), ideal(_ideal*RADDEG), stdval(_stdval*RADDEG), KI(_KI) {}
    
	Point3s basecross(Point3s p1, Point3s p2) const;
    double calcenergyangle() const;
    void   calcgradientangle() const;

	static BondAngleobj create(Point3s* _teAtomA, Point3s* _teAtomB, Point3s* _teAtomC, Point3s* _teAtomAgrad, Point3s* _teAtomBgrad, Point3s* _teAtomCgrad, double _ideal, double _stdval, double _KI) {
        return BondAngleobj(_teAtomA, _teAtomB, _teAtomC, _teAtomAgrad, _teAtomBgrad, _teAtomCgrad, _ideal, _stdval, _KI);
    }
};

class MultiBondDihedralobj
{
	private:
		vector<Point3s*> teAtoms;
		vector<Point3s*> teAtomgrad;
		vector<vector<vector<double*>>> ideals;
		int param_cls_num;
		double K;

		vector<vector<int>> torsion_atom_idx = {{0,1,2,3}, {1,2,3,10}, {1,4,5,6}, {4,5,6,2}, {5,6,2,1}, {6,2,1,4}, {2,1,4,5}, {4,5,8,9}};
	public:
		MultiBondDihedralobj(const vector<Point3s*>& _AtomPtrs, const vector<Point3s*>& _Atomgrad, const vector<vector<vector<double*>>> &_ideals, const double &_K);

		double calcenergymultidihedral() const;
		void calcgradientmultidihedral();
};

class BondDihedralobj
{
  private:
  	Point3s* teAtomA;
	Point3s* teAtomB;
	Point3s* teAtomC;
	Point3s* teAtomD;
	Point3s* teAtomAgrad;
	Point3s* teAtomBgrad;
	Point3s* teAtomCgrad;
	Point3s* teAtomDgrad;
    double K, ideal, rotn;

  public:
  	BondDihedralobj(Point3s* _teAtomA, Point3s* _teAtomB, Point3s* _teAtomC, Point3s* _teAtomD, Point3s* _teAtomAgrad, Point3s* _teAtomBgrad, Point3s* _teAtomCgrad, Point3s* _teAtomDgrad, double _K, double _ideal, double _rotn)
  		: teAtomA(_teAtomA), teAtomB(_teAtomB), teAtomC(_teAtomC), teAtomD(_teAtomD), teAtomAgrad(_teAtomAgrad), teAtomBgrad(_teAtomBgrad), teAtomCgrad(_teAtomCgrad), teAtomDgrad(_teAtomDgrad), K(_K), ideal(_ideal*RADDEG), rotn(_rotn) {}
    
    double calcenergydihedral() const;
    void   calcgradientdihedral() const;

	static BondDihedralobj create(Point3s* _teAtomA, Point3s* _teAtomB, Point3s* _teAtomC, Point3s* _teAtomD, Point3s* _teAtomAgrad, Point3s* _teAtomBgrad, Point3s* _teAtomCgrad, Point3s* _teAtomDgrad, double _K, double _ideal, double _rotn) {
        return BondDihedralobj(_teAtomA, _teAtomB, _teAtomC, _teAtomD, _teAtomAgrad, _teAtomBgrad, _teAtomCgrad, _teAtomDgrad, _K, _ideal, _rotn);
    }
};

class BondImpropersobj
{
  private:
  	Point3s* teAtomA;
	Point3s* teAtomB;
	Point3s* teAtomC;
	Point3s* teAtomD;
	Point3s* teAtomAgrad;
	Point3s* teAtomBgrad;
	Point3s* teAtomCgrad;
	Point3s* teAtomDgrad;
    double KI, PN, Phase, iDiv;

  public:
  	BondImpropersobj(Point3s* _teAtomA, Point3s* _teAtomB, Point3s* _teAtomC, Point3s* _teAtomD, Point3s* _teAtomAgrad, Point3s* _teAtomBgrad, Point3s* _teAtomCgrad, Point3s* _teAtomDgrad, double _KI, double _PN, double _Phase)
  		: teAtomA(_teAtomA), teAtomB(_teAtomB), teAtomC(_teAtomC), teAtomD(_teAtomD), teAtomAgrad(_teAtomAgrad), teAtomBgrad(_teAtomBgrad), teAtomCgrad(_teAtomCgrad), teAtomDgrad(_teAtomDgrad), KI(_KI), PN(_PN), Phase(_Phase) {}
    
	Point3s basecross(Point3s p1, Point3s p2) const;
    double calcenergyimpropers()  const;
    void   calcgradientimpropers() const;

	static BondImpropersobj create(Point3s* _teAtomA, Point3s* _teAtomB, Point3s* _teAtomC, Point3s* _teAtomD, Point3s* _teAtomAgrad, Point3s* _teAtomBgrad, Point3s* _teAtomCgrad, Point3s* _teAtomDgrad, double _KI, double _PN, double _Phase) {
        return BondImpropersobj(_teAtomA, _teAtomB, _teAtomC, _teAtomD, _teAtomAgrad, _teAtomBgrad, _teAtomCgrad, _teAtomDgrad, _KI, _PN, _Phase);
    }
};

class Hbondobj
{
  private:
  	Point3s* teAtomHyd;
	Point3s* teAtomDon;
	Point3s* teAtomAcc;
	Point3s* teAtomHydgrad;
	Point3s* teAtomDongrad;
	Point3s* teAtomAccgrad;
	double idealadleng, idealahleng, idealahdang;
	const double dK;
    double dREq;
	double dXiR;
	double dXiTheta;
	double dThreshold;
	double beta=25;

  public:
  	Hbondobj(Point3s* _teAtomHyd, Point3s* _teAtomDon, Point3s* _teAtomAcc, Point3s* _teAtomHydgrad, Point3s* _teAtomDongrad, Point3s* _teAtomAccgrad, const double &_idealadleng, const double &_idealahleng, const double &_idealahdang, const double &_dK)
  		: teAtomHyd(_teAtomHyd), teAtomDon(_teAtomDon), teAtomAcc(_teAtomAcc), teAtomHydgrad(_teAtomHydgrad), teAtomDongrad(_teAtomDongrad), teAtomAccgrad(_teAtomAccgrad), idealadleng(_idealadleng), idealahleng(_idealahleng), idealahdang(cos(_idealahdang*RADDEG)),
		  dK(_dK), dREq(0.0), dXiR(1.0), dXiTheta(1.0), dThreshold(-1.0) {}
    
    double calcenergyhbond() const;
    void   calcgradienthbond() const;

	//bool bIsNan(const double _dVal);

	static Hbondobj create(Point3s* _teAtomHyd, Point3s* _teAtomDon, Point3s* _teAtomAcc, Point3s* _teAtomHydgrad, Point3s* _teAtomDongrad, Point3s* _teAtomAccgrad, const double &_idealadleng, const double &_idealahleng, const double &_idealahdang, const double &_dK) {
        return Hbondobj(_teAtomHyd, _teAtomDon, _teAtomAcc, _teAtomHydgrad, _teAtomDongrad, _teAtomAccgrad, _idealadleng, _idealahleng, _idealahdang, _dK);
    }
};

class Hbond2hydobj
{
  private:
  	Point3s* teAtomHyd1;
	Point3s* teAtomHyd2;
	Point3s* teAtomDon;
	Point3s* teAtomAcc;
	Point3s* teAtomHydgrad1;
	Point3s* teAtomHydgrad2;
	Point3s* teAtomDongrad;
	Point3s* teAtomAccgrad;
	double idealadleng, idealahleng, idealahdang;
	const double dK;
    double dREq;
	double dXiR;
	double dXiTheta;
	double dThreshold;
	double beta=25;

  public:
  	Hbond2hydobj(Point3s* _teAtomHyd1, Point3s* _teAtomHyd2, Point3s* _teAtomDon, Point3s* _teAtomAcc, Point3s* _teAtomHydgrad1, Point3s* _teAtomHydgrad2, Point3s* _teAtomDongrad, Point3s* _teAtomAccgrad, const double &_idealadleng, const double &_idealahleng, const double &_idealahdang, const double &_dK)
  		: teAtomHyd1(_teAtomHyd1), teAtomHyd2(_teAtomHyd2), teAtomDon(_teAtomDon), teAtomAcc(_teAtomAcc), teAtomHydgrad1(_teAtomHydgrad1), teAtomHydgrad2(_teAtomHydgrad2), teAtomDongrad(_teAtomDongrad), teAtomAccgrad(_teAtomAccgrad), idealadleng(_idealadleng), idealahleng(_idealahleng), idealahdang(cos(_idealahdang*RADDEG)),
		  dK(_dK), dREq(0.0), dXiR(1.0), dXiTheta(1.0), dThreshold(-1.0) {}
    
    double calcenergyhbond() const;
    void   calcgradienthbond() const;

	//bool bIsNan(const double _dVal);

	static Hbond2hydobj create(Point3s* _teAtomHyd1, Point3s* _teAtomHyd2, Point3s* _teAtomDon, Point3s* _teAtomAcc, Point3s* _teAtomHydgrad1, Point3s* _teAtomHydgrad2, Point3s* _teAtomDongrad, Point3s* _teAtomAccgrad, const double &_idealadleng, const double &_idealahleng, const double &_idealahdang, const double &_dK) {
        return Hbond2hydobj(_teAtomHyd1, _teAtomHyd2, _teAtomDon, _teAtomAcc, _teAtomHydgrad1, _teAtomHydgrad2, _teAtomDongrad, _teAtomAccgrad, _idealadleng, _idealahleng, _idealahdang, _dK);
    }
};

class BasePairobj
{
	private:
		vector<Point3s*> AtomPtrs;
		vector<Point3s*> Atomgrad;
		double Xm,Ym,Zm;
		double a,b;
		unsigned int n;
		double K;

		double Fit3DPlane();
		double TotalSQRDistance(double _a, double _b) const;
    	double TotalSQRDistance() const;
    	void CalcCenterOfMass();

		//double dK = 5.0;

	public:
		BasePairobj(const vector<Point3s*>& _AtomPtrs, const vector<Point3s*>& _Atomgrad, const double &_K);
			//: AtomPtrs(_AtomPtrs), Atomgrad(_Atomgrad) {}

		double calcenergybasepair() const;
    	void calcgradientbasepair();
};

class TemplateSpuobj
{
	private:
		vector<Point3s*> AtomPtrs;
		vector<Point3s*> Atomgrad;
		vector<vector<Point3s*>> TempAtomPtrs;
		int templateidx;
		BasicFunc bf;
		double dK = 1.0;

	public:
		TemplateSpuobj(const vector<Point3s*>& _AtomPtrs, const vector<Point3s*>& _Atomgrad, const vector<vector<Point3s*>> &_TempAtomPtrs, const double &_dK);
		int selecttemplate();
};

class Stackingobj
{
  private:
  	Point3s* teAtomA;
	Point3s* teAtomB;
	Point3s* teAtomAgrad;
	Point3s* teAtomBgrad;
    double ideal;
	double KI;

  public:
  Stackingobj(Point3s* _teAtomA, Point3s* _teAtomB, Point3s* _teAtomAgrad, Point3s* _teAtomBgrad, double _ideal, double _KI)
  		: teAtomA(_teAtomA), teAtomB(_teAtomB), teAtomAgrad(_teAtomAgrad), teAtomBgrad(_teAtomBgrad), ideal(_ideal), KI(_KI) {}
    
    double calcenergystack() const;
    void   calcgradientstack() const;

	static Stackingobj create(Point3s* _teAtomA, Point3s* _teAtomB, Point3s* _teAtomAgrad, Point3s* _teAtomBgrad, double _ideal, double _KI) {
        return Stackingobj(_teAtomA, _teAtomB, _teAtomAgrad, _teAtomBgrad, _ideal, _KI);
    }
};

class BPdistobj
{
  private:
  	Point3s* teAtomA;
	Point3s* teAtomB;
	Point3s* teAtomAgrad;
	Point3s* teAtomBgrad;
    double ideal;
	double KI;

  public:
  BPdistobj(Point3s* _teAtomA, Point3s* _teAtomB, Point3s* _teAtomAgrad, Point3s* _teAtomBgrad, double _ideal, double _KI)
  		: teAtomA(_teAtomA), teAtomB(_teAtomB), teAtomAgrad(_teAtomAgrad), teAtomBgrad(_teAtomBgrad), ideal(_ideal), KI(_KI) {}
    
    double calcenergybp() const;
    void   calcgradientbp() const;

	static BPdistobj create(Point3s* _teAtomA, Point3s* _teAtomB, Point3s* _teAtomAgrad, Point3s* _teAtomBgrad, double _ideal, double _KI) {
        return BPdistobj(_teAtomA, _teAtomB, _teAtomAgrad, _teAtomBgrad, _ideal, _KI);
    }
};

class StackTemplateobj
{
  private:
  	FragmentResidue* ResA;
	FragmentResidue* ResB;
	vector<Point3s>* ResAgrad;
	vector<Point3s>* ResBgrad;
	const vector<ModelUnit>* ideal_rna;
	int templatenum;
	vector<vector<double>> supResA;
	vector<vector<double>> supResB;
	const vector<int>* indices;
	const vector<int>* orders;
	BasicFunc bf;

	vector<vector<double> > xyz_list1;
	vector<vector<double> > xyz_list2;

	string base1, base2;

	map<string, vector<int>> targetatoms={
        {"A",{12, 13, 14, 15, 16, 17, 18, 19, 20, 21}},
        {"C",{12, 13, 14, 15, 16, 17, 18, 19}},
        {"G",{12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22}},
        {"U",{12, 13, 14, 15, 16, 17, 18, 19}}
    };

	double K;

  	public:
  	StackTemplateobj(FragmentResidue* _ResA, FragmentResidue* _ResB, vector<Point3s>* _ResAgrad, vector<Point3s>* _ResBgrad, const vector<ModelUnit>* _ideal_rna, const vector<int>* _indices, const vector<int>* _orders, const int &_templatenum, double _K);
    
	double selecttemplate();
    double calcenergystacktemp();
    void   calcgradientstacktemp() const;
};

class StackTemplateobj2
{
  private:
  	FragmentResidue* ResA;
	FragmentResidue* ResB;
	vector<Point3s>* ResAgrad;
	vector<Point3s>* ResBgrad;
	const vector<ModelUnit>* ideal_rna;
	int templatenum;
	vector<vector<double>> supResA;
	vector<vector<double>> supResB;
	const vector<int>* indices;
	const vector<int>* orders;
	BasicFunc bf;

	vector<vector<double> > xyz_listA1;
	vector<vector<double> > xyz_listA2;
	vector<vector<double> > xyz_listA3;

	vector<vector<double> > xyz_listB1;
	vector<vector<double> > xyz_listB2;
	vector<vector<double> > xyz_listB3;

	string base1, base2;

	map<string, vector<int>> targetatoms={
        {"A",{12, 13, 14, 15, 16, 17, 18, 19, 20, 21}},
        {"C",{12, 13, 14, 15, 16, 17, 18, 19}},
        {"G",{12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22}},
        {"U",{12, 13, 14, 15, 16, 17, 18, 19}}
    };

	double K;

  	public:
  	StackTemplateobj2(FragmentResidue* _ResA, FragmentResidue* _ResB, vector<Point3s>* _ResAgrad, vector<Point3s>* _ResBgrad, const vector<ModelUnit>* _ideal_rna, const vector<int>* _indices, const vector<int>* _orders, const int &_templatenum, double _K);
    
	double selecttemplate();
    double calcenergystacktemp();
    void   calcgradientstacktemp() const;
};

class SingleNtTemplateobj
{
 private:
  	FragmentResidue* Res;
	vector<Point3s>* Resgrad;
	const vector<ModelUnit>* ideal_rna;
	int templatenum;
	vector<vector<double>> supRes;
	const vector<int>* indices;
	BasicFunc bf;

	vector<vector<double> > xyz_list1;
	vector<vector<double> > xyz_list2;
	vector<vector<double> > xyz_list3;

	string base;

	map<string, vector<int>> targetatoms={
        {"A",{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21}},
        {"C",{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}},
        {"G",{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22}},
        {"U",{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}}
    };

	double K;

  	public:
  	SingleNtTemplateobj(FragmentResidue* _Res, vector<Point3s>* _Resgrad, const vector<ModelUnit>* _ideal_rna, const vector<int>* _indices, const int &_templatenum, double _K);
    
	double selecttemplate();
    double calcenergystacktemp();
    void   calcgradientstacktemp() const;
};

class SingleNtTemplateobjMC
{
  private:
  	Point3s* Res;
	vector<Point3s>* Resgrad;
	const vector<vector<Residueatoms>>* ideal_rna;
	int templatenum;
	vector<vector<double>> supRes;
	const vector<int>* indices;
	BasicFunc bf;

	vector<vector<double> > xyz_list1;
	vector<vector<double> > xyz_list2;

	char base;
	int atomnum;

	/*
	map<char, vector<int>> reptargetatoms={
        {'A',{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21}},
        {'C',{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}},
        {'G',{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22}},
        {'U',{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}}
    };
	
	map<char, vector<int>> suptargetatoms={
        {'A',{0, 1, 2, 5, 12, 13}},
        {'C',{0, 1, 2, 5, 12, 13}},
        {'G',{0, 1, 2, 5, 12, 13}},
        {'U',{0, 1, 2, 5, 12, 13}}
    };
	*/
	map<char, vector<int>> reptargetatoms={
        {'A',{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21}},
        {'C',{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}},
        {'G',{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22}},
        {'U',{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}}
    };
	
	map<char, vector<int>> suptargetatoms={
        {'A',{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21}},
        {'C',{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}},
        {'G',{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22}},
        {'U',{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}}
    };


	map<char, int> targetatomnum={{'A',22},{'C',20},{'G',23},{'U',20}};

	double K;

  	public:
  	SingleNtTemplateobjMC(Point3s* _Res, vector<Point3s>* _Resgrad, char &_base, const vector<vector<Residueatoms>>* _ideal_rna, const vector<int>* _indices, const int &_templatenum, double _K);
    
	double selecttemplate();
    double calcenergystacktemp();
    void   calcgradientstacktemp() const;
};

class StackTemplateobj2MC
{
  private:
  	Point3s* ResA;
	Point3s* ResB;
	vector<Point3s>* ResAgrad;
	vector<Point3s>* ResBgrad;
	const vector<vector<Residueatoms>>* ideal_rna;
	int templatenum;
	vector<vector<double>> supResA;
	vector<vector<double>> supResB;
	const vector<int>* indices;
	const vector<int>* orders;
	BasicFunc bf;

	vector<vector<double> > xyz_listA1;
	vector<vector<double> > xyz_listA2;
	vector<vector<double> > xyz_listA3;

	vector<vector<double> > xyz_listB1;
	vector<vector<double> > xyz_listB2;
	vector<vector<double> > xyz_listB3;

	char base1, base2;
	int atomnumA, atomnumB;

	map<char, vector<int>> targetatoms={
        {'A',{12, 13, 14, 15, 16, 17, 18, 19, 20, 21}},
        {'C',{12, 13, 14, 15, 16, 17, 18, 19}},
        {'G',{12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22}},
        {'U',{12, 13, 14, 15, 16, 17, 18, 19}}
    };
	

	map<char, int> targetatomnum={{'A',22},{'C',20},{'G',23},{'U',20}};

	double K;

  	public:
  	StackTemplateobj2MC(Point3s* _ResA, Point3s* _ResB, vector<Point3s>* _ResAgrad, vector<Point3s>* _ResBgrad, char &_base1, char &_base2, const vector<vector<Residueatoms>>* _ideal_rna, const vector<int>* _indices, const vector<int>* _orders, const int &_templatenum, double _K);
    
	double selecttemplate();
    double calcenergystacktemp();
    void   calcgradientstacktemp() const;
};

class StackGeoobj
{
  private:
  	Point3s* ResA;
	Point3s* ResB;
	vector<Point3s>* ResAgrad;
	vector<Point3s>* ResBgrad;
	char baseA, baseB;
	int atomnumA, atomnumB;
	BasicFunc bf;
	double K1, K2;
	double ideallength, idealangle;

	map<char, vector<int>> targetatoms={
		{'A',{0,1,2,3,4,5,6,7,8,9}},
        {'C',{0,1,2,3,4,5,6,7}},
        {'G',{0,1,2,3,4,5,6,7,8,9,10}},
        {'U',{0,1,2,3,4,5,6,7}}
	};

  	public:
  	StackGeoobj(Point3s* _ResA, Point3s* _ResB, vector<Point3s>* _ResAgrad, vector<Point3s>* _ResBgrad, char &_baseA, char &_baseB, const double &_ideallength, const double &_idealangle, double _K1, double _K2);
    
	double calcenergy();
	void calcgradient();
    
};

class StackGeoobjMC
{
  private:
  	vector<Point3s*> ResA;
	vector<Point3s*> ResB;
	vector<Point3s*> ResAgrad;
	vector<Point3s*> ResBgrad;
	vector<Point3s> gradA;
	vector<Point3s> gradB;
	char baseA, baseB;
	int atomnumA, atomnumB;
	BasicFunc bf;
	double K1, K2;
	double ideallength, idealangle;

	map<char, vector<int>> targetatoms={
		{'A',{0,1,2,3,4,5,6,7,8,9}},
        {'C',{0,1,2,3,4,5,6,7}},
        {'G',{0,1,2,3,4,5,6,7,8,9,10}},
        {'U',{0,1,2,3,4,5,6,7}}
	};

  	public:
  	StackGeoobjMC(vector<Point3s*>& _ResA, vector<Point3s*>& _ResB, vector<Point3s*>& _ResAgrad, vector<Point3s*>& _ResBgrad, char &_baseA, char &_baseB, const double &_ideallength, const double &_idealangle, double _K1, double _K2);
    
	double calcenergy() const;
	void calcgradient();
    
};

class Clashobj
{
  private:
  	Point3s* teAtomA;
	Point3s* teAtomB;
	Point3s* teAtomAgrad;
	Point3s* teAtomBgrad;
    double vwd;
	double K;

  public:
  Clashobj(Point3s* _teAtomA, Point3s* _teAtomB, Point3s* _teAtomAgrad, Point3s* _teAtomBgrad, double _vwd, double _K)
  		: teAtomA(_teAtomA), teAtomB(_teAtomB), teAtomAgrad(_teAtomAgrad), teAtomBgrad(_teAtomBgrad), vwd(_vwd), K(_K) {}
    
    double calcenergyclash() const;
    void   calcgradientclash() const;

	static Clashobj create(Point3s* _teAtomA, Point3s* _teAtomB, Point3s* _teAtomAgrad, Point3s* _teAtomBgrad, double _vwd, double _K) {
        return Clashobj(_teAtomA, _teAtomB, _teAtomAgrad, _teAtomBgrad, _vwd, _K);
    }
};

#endif













