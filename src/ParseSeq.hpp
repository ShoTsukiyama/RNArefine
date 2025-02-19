// ParseSeq.h: interface for the ParseSeq class.
//
//////////////////////////////////////////////////////////////////////

//#if !defined(AFX_PARSESEQ_H__B85D878F_294B_49A8_B196_3D23A1C37A54__INCLUDED_)
//#define AFX_PARSESEQ_H__B85D878F_294B_49A8_B196_3D23A1C37A54__INCLUDED_

//#if _MSC_VER > 1000
#pragma once
#include <stdlib.h>

#include "RNArefine.hpp"
#include "CommonPara.hpp"
#include "LoadPDB.hpp"
//#endif // _MSC_VER > 1000

using namespace std;

class ParseSeq  
{
public:
	ParseSeq();
	
	string pdb2seq(vector<Boneinfo> bb);
	void writeseq(string outfile, string seq, string header);
	virtual ~ParseSeq();
};

//#endif // !defined(AFX_PARSESEQ_H__B85D878F_294B_49A8_B196_3D23A1C37A54__INCLUDED_)
