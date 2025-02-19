// ParseSeq.cpp: implementation of the ParseSeq class.
//
//////////////////////////////////////////////////////////////////////
//#include <string>
//#include <map>

#include "ParseSeq.hpp"
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ParseSeq::ParseSeq()
{
	int i,j;
}

string ParseSeq::pdb2seq(vector<Boneinfo> bb)
{
	string seqdata;

	int indres;
	for(int i=0;i<bb.size();i++)
	{
		indres=Reschar2Idint[bb[i].resid];
		if(indres>3) indres=4;
		seqdata+=Idint2Reschar[indres];	
	}

	return seqdata;
}

void ParseSeq::writeseq(string outfile, string seq, string header)
{
	ofstream outputFile(outfile.c_str());

    if (!outputFile.is_open()) {
        cerr << "Error: Unable to open file " << outfile << endl;
        return;
    }

    outputFile << ">" << header << endl;
    outputFile << seq << endl;

    outputFile.close();
}

ParseSeq::~ParseSeq()
{
}

