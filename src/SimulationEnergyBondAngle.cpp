// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

void RNASimulation::fenergybondangle(vector<Residueatoms> decstr, int numseq,double &nene,double &fene, Energyterms &lp, bool flagacc)
{
	int i,j;
	int istart;
	Point3d tp1,tp2;
	char aat;
	int ss;
	double totalnn;
	double diffang;

	//int estartidx=BLENGTH_ENERGY_NUM+ILENGTH_ENERGY_NUM;

	if(!flagacc)
	{
		for(j=0;j<INNERANGLETARGET.size();j++){
			lp.fn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]]=0;
			lp.nn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]]=0;

			aat=decstr[0].aaa;
			if(INNERANGLETARGET[j] == 10){
				getouteneangle(decstr[0].ang[INNERANGLETARGET[j]], OPINNERANGPARAM[aat][0].mean, OPINNERANGPARAM[aat][0].lower, OPINNERANGPARAM[aat][0].upper, OPINNERANGPARAM[aat][0].slower, OPINNERANGPARAM[aat][0].supper, 1, lp.indf[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0], lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0]);
				//getoutlier(decstr[0].ang[INNERANGLETARGET[j]], OPINNERANGPARAM[aat][0].boundary, OPINNERANGPARAM[aat][0].energy, OPINNERANGPARAM[aat][0].binary, lp.indf[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0], lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0]);
				//getoutene(decstr[0].ang[INNERANGLETARGET[j]], OPINNERANGPARAM[aat][0].mode, OPINNERANGPARAM[aat][0].boundary, OPINNERANGPARAM[aat][0].energy, OPINNERANGPARAM[aat][0].binary, lp.indf[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0], lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0]);
			
				//if(lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0] != 0){
					//cout<<to_string(i) + "-" + to_string(j) + "-" + to_string(decstr[0].ang[INNERANGLETARGET[j]]) + "-" + to_string(OPINNERANGPARAM[aat][0].mean) + "-" + to_string(OPINNERANGPARAM[aat][0].lower) + "-" + to_string(OPINNERANGPARAM[aat][0].upper) + "-" + to_string(lp.indf[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0]) + "-" + to_string(lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0])<<endl;
				//}
			
			}else if(INNERANGLETARGET[j] == 11){
				getouteneangle(decstr[0].ang[INNERANGLETARGET[j]], OPINNERANGPARAM[aat][1].mean, OPINNERANGPARAM[aat][1].lower, OPINNERANGPARAM[aat][1].upper, OPINNERANGPARAM[aat][1].slower, OPINNERANGPARAM[aat][1].supper, 1, lp.indf[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0], lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0]);
				//getoutlier(decstr[0].ang[INNERANGLETARGET[j]], OPINNERANGPARAM[aat][1].boundary, OPINNERANGPARAM[aat][1].energy, OPINNERANGPARAM[aat][1].binary, lp.indf[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0], lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0]);
				//getoutene(decstr[0].ang[INNERANGLETARGET[j]], OPINNERANGPARAM[aat][1].mode, OPINNERANGPARAM[aat][1].boundary, OPINNERANGPARAM[aat][1].energy, OPINNERANGPARAM[aat][1].binary, lp.indf[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0], lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0]);
			
				//if(lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0] != 0){
					//cout<<to_string(i) + "-" + to_string(j) + "-" + to_string(decstr[0].ang[INNERANGLETARGET[j]]) + "-" + to_string(OPINNERANGPARAM[aat][1].mean) + "-" + to_string(OPINNERANGPARAM[aat][1].lower) + "-" + to_string(OPINNERANGPARAM[aat][1].upper) + "-" + to_string(lp.indf[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0]) + "-" + to_string(lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0])<<endl;
				//}
			
			}else{
				getouteneangle(decstr[0].ang[INNERANGLETARGET[j]], INNERANGPARAM[aat][INNERANGLETARGET[j]].mean, INNERANGPARAM[aat][INNERANGLETARGET[j]].lower, INNERANGPARAM[aat][INNERANGLETARGET[j]].upper, INNERANGPARAM[aat][INNERANGLETARGET[j]].slower, INNERANGPARAM[aat][INNERANGLETARGET[j]].supper, 1, lp.indf[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0], lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0]);
				//getoutlier(decstr[0].ang[INNERANGLETARGET[j]], INNERANGPARAM[aat][INNERANGLETARGET[j]].boundary, INNERANGPARAM[aat][INNERANGLETARGET[j]].energy, INNERANGPARAM[aat][INNERANGLETARGET[j]].binary, lp.indf[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0], lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0]);
				//getoutene(decstr[0].ang[INNERANGLETARGET[j]], INNERANGPARAM[aat][INNERANGLETARGET[j]].mode, INNERANGPARAM[aat][INNERANGLETARGET[j]].boundary, INNERANGPARAM[aat][INNERANGLETARGET[j]].energy, INNERANGPARAM[aat][INNERANGLETARGET[j]].binary, lp.indf[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0], lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0]);
			
				//if(lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0] != 0){
					//cout<<to_string(i) + "-" + to_string(j) + "-" + to_string(decstr[0].ang[INNERANGLETARGET[j]]) + "-" + to_string(INNERANGPARAM[aat][INNERANGLETARGET[j]].mean) + "-" + to_string(INNERANGPARAM[aat][INNERANGLETARGET[j]].lower) + "-" + to_string(INNERANGPARAM[aat][INNERANGLETARGET[j]].upper) + "-" + to_string(lp.indf[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0]) + "-" + to_string(lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0])<<endl;
				//}
			
			}
			lp.fn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]]+=lp.indf[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0];
			lp.nn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]]+=lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0];

			for(i=1;i<numseq;i++)
			{
				aat=decstr[i].aaa;
				getouteneangle(decstr[i].ang[INNERANGLETARGET[j]], INNERANGPARAM[aat][INNERANGLETARGET[j]].mean, INNERANGPARAM[aat][INNERANGLETARGET[j]].lower, INNERANGPARAM[aat][INNERANGLETARGET[j]].upper, INNERANGPARAM[aat][INNERANGLETARGET[j]].slower, INNERANGPARAM[aat][INNERANGLETARGET[j]].supper, 1, lp.indf[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][i], lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][i]);

				//if(lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][i] != 0){
					//cout<<"normal-" + to_string(i) + "-" + to_string(j) + "-" + to_string(decstr[i].ang[INNERANGLETARGET[j]]) + "-" + to_string(INNERANGPARAM[aat][INNERANGLETARGET[j]].mean) + "-" + to_string(INNERANGPARAM[aat][INNERANGLETARGET[j]].lower) + "-" + to_string(INNERANGPARAM[aat][INNERANGLETARGET[j]].upper) + "-" + to_string(lp.indf[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][i]) + "-" + to_string(lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][i])<<endl;
				//}

				lp.fn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]]+=lp.indf[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][i];
				lp.nn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]]+=lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][i];
			}
		}

		for(j=0;j<ADDINNERANGLETARGET.size();j++)
		{
			lp.fn[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]]=0;
			lp.nn[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]]=0;

			if(ADDINNERANGLETARGET[j] == 0) istart=1;
			else istart=0;

			for(i=istart;i<numseq;i++)
			{
				aat=decstr[i].aaa;
				getouteneangle(decstr[i].addang[ADDINNERANGLETARGET[j]], ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].mean, ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].lower, ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].upper, ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].slower, ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].supper, 1, lp.indf[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]][i], lp.indn[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]][i]);

				//if(lp.indn[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]][i] != 0){
					//cout<<"add-" + to_string(i) + "-" + to_string(j) + "-" + to_string(decstr[i].addang[ADDINNERANGLETARGET[j]]) + "-" + to_string(ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].mean) + "-" + to_string(ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].lower) + "-" + to_string(ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].upper) + "-" + to_string(lp.indf[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]][i]) + "-" + to_string(lp.indn[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]][i])<<endl;
				//}

				lp.fn[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]]+=lp.indf[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]][i];
				lp.nn[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]]+=lp.indn[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]][i];
			}
		}
	}

	else
	{
		int iold[INNERANGLETARGET.size()],inew[INNERANGLETARGET.size()];
		double fold[INNERANGLETARGET.size()],fnew[INNERANGLETARGET.size()];
		for(i=0;i<INNERANGLETARGET.size();i++){
			iold[i]=0; inew[i]=0; fold[i]=0; fnew[i]=0;
		}
		
		for(j=0;j<INNERANGLETARGET.size();j++){

			if(lp.flap[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0]){
				iold[j]+=lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0];
				fold[j]+=lp.indf[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0];

				aat=decstr[0].aaa;
				if(INNERANGLETARGET[j] == 10){
					getouteneangle(decstr[0].ang[INNERANGLETARGET[j]], OPINNERANGPARAM[aat][0].mean, OPINNERANGPARAM[aat][0].lower, OPINNERANGPARAM[aat][0].upper, OPINNERANGPARAM[aat][0].slower, OPINNERANGPARAM[aat][0].supper, 1, lp.indf[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0], lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0]);
					//getoutlier(decstr[0].ang[INNERANGLETARGET[j]], OPINNERANGPARAM[aat][0].boundary, OPINNERANGPARAM[aat][0].energy, OPINNERANGPARAM[aat][0].binary, lp.indf[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0], lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0]);
					//getoutene(decstr[0].ang[INNERANGLETARGET[j]], OPINNERANGPARAM[aat][0].mode, OPINNERANGPARAM[aat][0].boundary, OPINNERANGPARAM[aat][0].energy, OPINNERANGPARAM[aat][0].binary, lp.indf[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0], lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0]);
				
				}else if(INNERANGLETARGET[j] == 11){
					getouteneangle(decstr[0].ang[INNERANGLETARGET[j]], OPINNERANGPARAM[aat][1].mean, OPINNERANGPARAM[aat][1].lower, OPINNERANGPARAM[aat][1].upper, OPINNERANGPARAM[aat][1].slower, OPINNERANGPARAM[aat][1].supper, 1, lp.indf[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0], lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0]);
					//getoutlier(decstr[0].ang[INNERANGLETARGET[j]], OPINNERANGPARAM[aat][1].boundary, OPINNERANGPARAM[aat][1].energy, OPINNERANGPARAM[aat][1].binary, lp.indf[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0], lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0]);
					//getoutene(decstr[0].ang[INNERANGLETARGET[j]], OPINNERANGPARAM[aat][1].mode, OPINNERANGPARAM[aat][1].boundary, OPINNERANGPARAM[aat][1].energy, OPINNERANGPARAM[aat][1].binary, lp.indf[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0], lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0]);
				}else{
					getouteneangle(decstr[0].ang[INNERANGLETARGET[j]], INNERANGPARAM[aat][INNERANGLETARGET[j]].mean, INNERANGPARAM[aat][INNERANGLETARGET[j]].lower, INNERANGPARAM[aat][INNERANGLETARGET[j]].upper, INNERANGPARAM[aat][INNERANGLETARGET[j]].slower, INNERANGPARAM[aat][INNERANGLETARGET[j]].supper, 1, lp.indf[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0], lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0]);
					//getoutlier(decstr[0].ang[INNERANGLETARGET[j]], INNERANGPARAM[aat][INNERANGLETARGET[j]].boundary, INNERANGPARAM[aat][INNERANGLETARGET[j]].energy, INNERANGPARAM[aat][INNERANGLETARGET[j]].binary, lp.indf[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0], lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0]);
					//getoutene(decstr[0].ang[INNERANGLETARGET[j]], INNERANGPARAM[aat][INNERANGLETARGET[j]].mode, INNERANGPARAM[aat][INNERANGLETARGET[j]].boundary, INNERANGPARAM[aat][INNERANGLETARGET[j]].energy, INNERANGPARAM[aat][INNERANGLETARGET[j]].binary, lp.indf[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0], lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0]);
				}

				inew[j]+=lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0];
				fnew[j]+=lp.indf[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0];
				lp.flap[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0]=false;
			}

			for(i=1;i<numseq;i++) if(lp.flap[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][i])
			{
				iold[j]+=lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][i];
				fold[j]+=lp.indf[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][i];

				aat=decstr[i].aaa;
				getouteneangle(decstr[i].ang[INNERANGLETARGET[j]], INNERANGPARAM[aat][INNERANGLETARGET[j]].mean, INNERANGPARAM[aat][INNERANGLETARGET[j]].lower, INNERANGPARAM[aat][INNERANGLETARGET[j]].upper, INNERANGPARAM[aat][INNERANGLETARGET[j]].slower, INNERANGPARAM[aat][INNERANGLETARGET[j]].supper, 1, lp.indf[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][i], lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][i]);

				//if(lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][i] != 0){
					//cout<<"normal-" + to_string(i) + "-" + to_string(j) + "-" + to_string(decstr[i].ang[INNERANGLETARGET[j]]) + "-" + to_string(INNERANGPARAM[aat][INNERANGLETARGET[j]].mean) + "-" + to_string(INNERANGPARAM[aat][INNERANGLETARGET[j]].lower) + "-" + to_string(INNERANGPARAM[aat][INNERANGLETARGET[j]].upper) + "-" + to_string(lp.indf[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][i]) + "-" + to_string(lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][i])<<endl;
				//}

				inew[j]+=lp.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][i];
				fnew[j]+=lp.indf[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][i];
				lp.flap[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][i]=false;
			}
			lp.nn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]]+=(inew[j]-iold[j]);
			lp.fn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]]+=(fnew[j]-fold[j]);
		}

		int addiold[ADDINNERANGLETARGET.size()],addinew[ADDINNERANGLETARGET.size()];
		double addfold[ADDINNERANGLETARGET.size()],addfnew[ADDINNERANGLETARGET.size()];
		for(i=0;i<ADDINNERANGLETARGET.size();i++)
		{
			addiold[i]=0; addinew[i]=0; addfold[i]=0; addfnew[i]=0;
		}

		for(j=0;j<ADDINNERANGLETARGET.size();j++)
		{
			if(ADDINNERANGLETARGET[j] == 0) istart=1;
			else istart=0;

			for(i=istart;i<numseq;i++) if(lp.flap[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]][i])
			{
				addiold[j]+=lp.indn[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]][i];
				addfold[j]+=lp.indf[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]][i];
				
				aat=decstr[i].aaa;
				getouteneangle(decstr[i].addang[ADDINNERANGLETARGET[j]], ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].mean, ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].lower, ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].upper, ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].slower, ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].supper, 1, lp.indf[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]][i], lp.indn[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]][i]);

				//if(lp.indn[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]][i] != 0){
					//cout<<"add-" + to_string(i) + "-" + to_string(j) + "-" + to_string(decstr[i].addang[ADDINNERANGLETARGET[j]]) + "-" + to_string(ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].mean) + "-" + to_string(ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].lower) + "-" + to_string(ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].upper) + "-" + to_string(lp.indf[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]][i]) + "-" + to_string(lp.indn[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]][i])<<endl;
				//}

				//if(lp.indn[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]][i] == 1){
					//cout<< i << "-" << j << "-" << lp.indf[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]][i]<<endl;
				//}

				addinew[j]+=lp.indn[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]][i];
				addfnew[j]+=lp.indf[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]][i];
				lp.flap[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]][i]=false;
			}

			lp.nn[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]]+=(addinew[j]-addiold[j]);
			lp.fn[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]]+=(addfnew[j]-addfold[j]);
		}
	}
	
	fene=0;
	for(j=0;j<INNERANGLETARGET.size();j++){
		fene+=lp.fn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]];
	}
	for(j=0;j<ADDINNERANGLETARGET.size();j++){
		fene+=1.5*lp.fn[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]];
	}

	nene=0;
	for(j=0;j<INNERANGLETARGET.size();j++){
		nene+=lp.nn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]];
	}
	for(j=0;j<ADDINNERANGLETARGET.size();j++){
		nene+=lp.nn[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]];
	}
}


















