// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

void RNASimulation::fenergybondlength(vector<Residueatoms> decstr,int numseq,double &nene,double &fene, Energyterms &lp, bool flagacc)
{
	int i,j;
	double tdist;
	double eval;
	int istart;
	Point3d tp;
	int binidx;
	char aat;
	int ss,ss1,ss2, eleidx;
	double lengdiff;
	double tmp1;
	int tmp2;

	if(!flagacc)
	{
		for(j=0;j<BLENGTHTARGET.size();j++)
		{
			lp.nn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]]=0;
			lp.fn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]]=0;

			aat=decstr[0].aaa;
			if(BLENGTHTARGET[j] == 10){
				getoutenelength(decstr[0].len[BLENGTHTARGET[j]], OPLENGTHPARAM[aat][0].mean, OPLENGTHPARAM[aat][0].lower, OPLENGTHPARAM[aat][0].upper, OPLENGTHPARAM[aat][0].slower, OPLENGTHPARAM[aat][0].supper, 1, lp.indf[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0], lp.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0]);
			}else if(BLENGTHTARGET[j] == 11){
				getoutenelength(decstr[0].len[BLENGTHTARGET[j]], OPLENGTHPARAM[aat][1].mean, OPLENGTHPARAM[aat][1].lower, OPLENGTHPARAM[aat][1].upper, OPLENGTHPARAM[aat][1].slower, OPLENGTHPARAM[aat][1].supper, 1, lp.indf[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0], lp.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0]);
			}else{
				getoutenelength(decstr[0].len[BLENGTHTARGET[j]], LENGTHPARAM[aat][BLENGTHTARGET[j]].mean, LENGTHPARAM[aat][BLENGTHTARGET[j]].lower, LENGTHPARAM[aat][BLENGTHTARGET[j]].upper, LENGTHPARAM[aat][BLENGTHTARGET[j]].slower, LENGTHPARAM[aat][BLENGTHTARGET[j]].supper, 1, lp.indf[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0], lp.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0]);
			}
			lp.nn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]]+=lp.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0];
			lp.fn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]]+=lp.indf[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0];

			//if(lp.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0] == 1){
				//cout<<to_string(i) + "-" + to_string(j) + "-" + to_string(decstr[0].len[BLENGTHTARGET[j]]) + "-" + to_string(LENGTHPARAM[aat][BLENGTHTARGET[j]].lower) + "-" + to_string(LENGTHPARAM[aat][BLENGTHTARGET[j]].upper)<<endl;
			//}

			for(i=1;i<numseq;i++)
			{
				aat=decstr[i].aaa;
				
				getoutenelength(decstr[i].len[BLENGTHTARGET[j]], LENGTHPARAM[aat][BLENGTHTARGET[j]].mean, LENGTHPARAM[aat][BLENGTHTARGET[j]].lower, LENGTHPARAM[aat][BLENGTHTARGET[j]].upper, LENGTHPARAM[aat][BLENGTHTARGET[j]].slower, LENGTHPARAM[aat][BLENGTHTARGET[j]].supper, 1, lp.indf[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][i], lp.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][i]);
				
				//if(lp.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][i] == 1){
					//cout<<to_string(i) + "-" + to_string(j) + "-" + to_string(decstr[i].len[BLENGTHTARGET[j]]) + "-" + to_string(LENGTHPARAM[aat][BLENGTHTARGET[j]].lower) + "-" + to_string(LENGTHPARAM[aat][BLENGTHTARGET[j]].upper)<<endl;
				//}

				lp.nn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]]+=lp.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][i];
				lp.fn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]]+=lp.indf[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][i];
			}
		}	

		for(j=0;j<ADDBLENGTHTARGET.size();j++)
		{
			lp.nn[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]]=0;
			lp.fn[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]]=0;

			for(i=0;i<numseq;i++)
			{
				aat=decstr[i].aaa;
				getoutenelength(decstr[i].addlen[ADDBLENGTHTARGET[j]], ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[j]].mean, ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[j]].lower, ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[j]].upper, ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[j]].slower, ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[j]].supper, 1, lp.indf[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]][i], lp.indn[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]][i]);
				
				lp.nn[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]]+=lp.indn[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]][i];
				lp.fn[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]]+=lp.indf[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]][i];
			}
		}
	}
	else
	{
		int iold[BLENGTHTARGET.size()],inew[BLENGTHTARGET.size()];
		double fold[BLENGTHTARGET.size()],fnew[BLENGTHTARGET.size()];
		for(i=0;i<BLENGTHTARGET.size();i++)
		{
			iold[i]=0; inew[i]=0; fold[i]=0; fnew[i]=0;
		}
		
		for(j=0;j<BLENGTHTARGET.size();j++)
		{
			if(lp.flap[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0]){
				iold[j]+=lp.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0];
				fold[j]+=lp.indf[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0];

				aat=decstr[0].aaa;
				if(BLENGTHTARGET[j] == 10){
					getoutenelength(decstr[0].len[BLENGTHTARGET[j]], OPLENGTHPARAM[aat][0].mean, OPLENGTHPARAM[aat][0].lower, OPLENGTHPARAM[aat][0].upper, OPLENGTHPARAM[aat][0].slower, OPLENGTHPARAM[aat][0].supper, 1, lp.indf[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0], lp.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0]);
					//getoutlier(decstr[0].len[BLENGTHTARGET[j]], OPLENGTHPARAM[aat][0].boundary, OPLENGTHPARAM[aat][0].energy, OPLENGTHPARAM[aat][0].binary, lp.indf[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0], lp.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0]);
					//getoutene(decstr[0].len[BLENGTHTARGET[j]], OPLENGTHPARAM[aat][0].mode, OPLENGTHPARAM[aat][0].boundary, OPLENGTHPARAM[aat][0].energy, OPLENGTHPARAM[aat][0].binary, lp.indf[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0], lp.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0]);
				}else if(BLENGTHTARGET[j] == 11){
					getoutenelength(decstr[0].len[BLENGTHTARGET[j]], OPLENGTHPARAM[aat][1].mean, OPLENGTHPARAM[aat][1].lower, OPLENGTHPARAM[aat][1].upper, OPLENGTHPARAM[aat][1].slower, OPLENGTHPARAM[aat][1].supper, 1, lp.indf[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0], lp.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0]);
					//getoutlier(decstr[0].len[BLENGTHTARGET[j]], OPLENGTHPARAM[aat][1].boundary, OPLENGTHPARAM[aat][1].energy, OPLENGTHPARAM[aat][1].binary, lp.indf[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0], lp.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0]);
					//getoutene(decstr[0].len[BLENGTHTARGET[j]], OPLENGTHPARAM[aat][1].mode, OPLENGTHPARAM[aat][1].boundary, OPLENGTHPARAM[aat][1].energy, OPLENGTHPARAM[aat][1].binary, lp.indf[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0], lp.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0]);
				}else{
					getoutenelength(decstr[0].len[BLENGTHTARGET[j]], LENGTHPARAM[aat][BLENGTHTARGET[j]].mean, LENGTHPARAM[aat][BLENGTHTARGET[j]].lower, LENGTHPARAM[aat][BLENGTHTARGET[j]].upper, LENGTHPARAM[aat][BLENGTHTARGET[j]].slower, LENGTHPARAM[aat][BLENGTHTARGET[j]].supper, 1, lp.indf[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0], lp.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0]);
					//getoutlier(decstr[0].len[BLENGTHTARGET[j]], LENGTHPARAM[aat][BLENGTHTARGET[j]].boundary, LENGTHPARAM[aat][BLENGTHTARGET[j]].energy, LENGTHPARAM[aat][BLENGTHTARGET[j]].binary, lp.indf[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0], lp.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0]);
					//getoutene(decstr[0].len[BLENGTHTARGET[j]], LENGTHPARAM[aat][BLENGTHTARGET[j]].mode, LENGTHPARAM[aat][BLENGTHTARGET[j]].boundary, LENGTHPARAM[aat][BLENGTHTARGET[j]].energy, LENGTHPARAM[aat][BLENGTHTARGET[j]].binary, lp.indf[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0], lp.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0]);
				}
				inew[j]+=lp.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0];
				fnew[j]+=lp.indf[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0];
				lp.flap[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0]=false;
			}

			for(i=1;i<numseq;i++) if(lp.flap[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][i]){
				iold[j]+=lp.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][i];
				fold[j]+=lp.indf[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][i];
				aat=decstr[i].aaa;

				getoutenelength(decstr[i].len[BLENGTHTARGET[j]], LENGTHPARAM[aat][BLENGTHTARGET[j]].mean, LENGTHPARAM[aat][BLENGTHTARGET[j]].lower, LENGTHPARAM[aat][BLENGTHTARGET[j]].upper, LENGTHPARAM[aat][BLENGTHTARGET[j]].slower, LENGTHPARAM[aat][BLENGTHTARGET[j]].supper, 1, lp.indf[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][i], lp.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][i]);

				//if(lp.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][i] != 0){
					//cout<<to_string(i) + "-" + to_string(j) + "-" + to_string(decstr[i].len[BLENGTHTARGET[j]]) + "-" + to_string(lp.indf[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][i]) + "-" + to_string(lp.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][i])<<endl;
				//}

				inew[j]+=lp.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][i];
				fnew[j]+=lp.indf[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][i];
				lp.flap[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][i]=false;
			}

			lp.nn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]]+=(inew[j]-iold[j]);
			lp.fn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]]+=(fnew[j]-fold[j]);
		}

		int addiold[ADDBLENGTHTARGET.size()],addinew[ADDBLENGTHTARGET.size()];
		double addfold[ADDBLENGTHTARGET.size()],addfnew[ADDBLENGTHTARGET.size()];
		for(i=0;i<ADDBLENGTHTARGET.size();i++)
		{
			addiold[i]=0; addinew[i]=0; addfold[i]=0; addfnew[i]=0;
		}

		for(j=0;j<ADDBLENGTHTARGET.size();j++)
		{
			for(i=0;i<numseq;i++)
			{
				if(lp.flap[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]][i]){
					addiold[j]+=lp.indn[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]][i];
					addfold[j]+=lp.indf[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]][i];

					aat=decstr[i].aaa;
					getoutenelength(decstr[i].addlen[ADDBLENGTHTARGET[j]], ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[j]].mean, ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[j]].lower, ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[j]].upper, ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[j]].slower, ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[j]].supper, 1, lp.indf[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]][i], lp.indn[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]][i]);
					
					addinew[j]+=lp.indn[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]][i];
					addfnew[j]+=lp.indf[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]][i];
					lp.flap[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]][i]=false;
				}
			}

			lp.nn[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]]+=(addinew[j]-addiold[j]);
			lp.fn[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]]+=(addfnew[j]-addfold[j]);
		}
	}

	fene=0;
	for(j=0;j<BLENGTHTARGET.size();j++){
		fene+=lp.fn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]];
	}
	for(j=0;j<ADDBLENGTHTARGET.size();j++)
	{
		fene+=lp.fn[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]];
	}

	nene=0;
	for(j=0;j<BLENGTHTARGET.size();j++){
		nene+=lp.nn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]];
	}
	for(j=0;j<ADDBLENGTHTARGET.size();j++){
		nene+=lp.nn[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]];
	}
}



















