// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

bool RNASimulation::writepdb(vector<Residueatoms> decstr,int numseq,vector<int> bbres, string outname,int ftype)
{
//	int j,k,l;
	int i, k,l;
	int totnu=1;
	char tpname[100];
	int indaa;
	FILE *file;

	if(ftype==0)
	{
		file=fopen(outname.c_str(),"wt");
	}
	else if(ftype==1)
	{
		file=fopen(outname.c_str(),"a+");
	}
	if(!file)
	{
		printf("cannot write file %s\n",outname.c_str());
		return false;
	}

	fprintf(file,"MODEL 1\n");
	for(i=0;i<numseq;i++)//pp.numbb
	{
		indaa=decstr[i].resind;
		for(k=0;k<SaveNN[indaa][0];k++)
		//for(k=0;k<NNAll[indaa][0];k++)
		{
			tpname[0]=' ';tpname[1]=' ';tpname[2]=' ';tpname[3]=' ';
			for(l=0;l<Id2atom[SaveNN[indaa][k+1]].length();l++)
			//for(l=0;l<Id2atom[NNAll[indaa][k+1]].length();l++)
			{
				tpname[l]=Id2atom[SaveNN[indaa][k+1]].substr(l,1)[0];
				//tpname[l]=Id2atom[NNAll[indaa][k+1]].substr(l,1)[0];
			}
			fprintf(file,"ATOM  %5d %c%c%c%c %c%c%c %c%4d    %8.3f%8.3f%8.3f %5.2f %5.2f\n",
				totnu,tpname[0],tpname[1],tpname[2],tpname[3],Ntd3[indaa][0],Ntd3[indaa][1],Ntd3[indaa][2],'A',
				bbres[i],decstr[i].ptv[k].x,decstr[i].ptv[k].y,decstr[i].ptv[k].z,1.0,0.0);
			totnu++;
		}
	}
	fprintf(file,"TER\nENDMDL\n");
	fclose(file);
	return true;
}

























