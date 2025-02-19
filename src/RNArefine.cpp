// TestingFunc.cpp: implementation of the TestingFunc class.
//
//////////////////////////////////////////////////////////////////////
#include <string>
#include <iostream>
#include <sys/stat.h>
#include <stdlib.h>
#include <filesystem>
#include <tuple>
#include <filesystem>
#include <thread>
#include <mutex>
#include <chrono>
#include <cmath>

#include "RNArefine.hpp"

using namespace std;
namespace fs = std::filesystem;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

RNArefineFunc::RNArefineFunc()
{

}

bool outputene(RNASimulation rsim, string outpath, double enesum){
	bool flagps;
	char buffer[1024];
	ofstream outputFile;
	outputFile.open(outpath); 
	if (!outputFile.is_open()) {
		cout << "Error:: simulation results cannot be written.\n";
		return false;
	}

	//snprintf(buffer, sizeof(buffer), "oene,%10.3f\nfapeinit,%10.3f\nfapedeeptmp,%10.3f\ncontact,%10.3f\nclashout,%10.3f\nclashloss,%10.3f\nlengthout,%10.3f\nlengthloss,%10.3f\nangleout,%10.3f\nangleloss,%10.3f\nmultitorsionout,%10.3f\nmultitorsionloss,%10.3f\naddlengthout,%10.3f\naddlengthloss,%10.3f\naddangleout,%10.3f\naddangleloss,%10.3f\naddtorsionout,%10.3f\naddtorsionloss,%10.3f\nbaseout,%10.3f\nbaseloss,%10.3f\nstackout,%10.3f\nstackloss,%10.3f\nhbond,%10.3f\n",
			//enesum,rsim.enelist[0], rsim.enelist[1], rsim.enelist[8], rsim.enelist[2], rsim.enelist[12], rsim.enelist[3], rsim.enelist[13], rsim.enelist[4], rsim.enelist[14], rsim.enelist[9], rsim.enelist[19], rsim.enelist[19], rsim.enelist[20], rsim.enelist[21], rsim.enelist[22], rsim.enelist[23], rsim.enelist[24], 
			//rsim.enelist[6], rsim.enelist[16], rsim.enelist[7], rsim.enelist[17], rsim.enelist[18]);

	snprintf(buffer, sizeof(buffer), "oene,%10.3f\nfapeinit,%10.3f\nfapedeeptmp,%10.3f\ncontact8A,%10.3f\ncontact20A,%10.3f\nclashout,%10.3f\nclashloss,%10.3f\nadjpairlengthout,%10.3f\nadjpairlengthloss,%10.3f\nntpairlengthout,%10.3f\nntpairlengthloss,%10.3f\nlengthout,%10.3f\nlengthloss,%10.3f\nangleout,%10.3f\nangleloss,%10.3f\ntorsionout,%10.3f\ntorsionloss,%10.3f\nmultitorsionout,%10.3f\nmultitorsionloss,%10.3f\npairedtorsionloss,%10.3f\npairedbbtorsionloss,%10.3f\nclusteredtorsionout,%10.3f\nclusteredtorsionloss,%10.3f\nbpout,%10.3f\nbploss,%10.3f\nbphydout,%10.3f\nbphydloss,%10.3f\nstackout,%10.3f\nstackloss,%10.3f\npuckerout,%10.3f\npuckerloss,%10.3f\nsinglenttemploss,%10.3f\nbptemploss,%10.3f\nstackingtemploss,%10.3f\n",
			enesum,rsim.enelist[0], rsim.enelist[1], rsim.enelist[8], rsim.enelist[18], rsim.enelist[2], rsim.enelist[12], rsim.enelist[26], rsim.enelist[27], rsim.enelist[28], rsim.enelist[29], rsim.enelist[3], rsim.enelist[13], rsim.enelist[4], rsim.enelist[14], rsim.enelist[5], rsim.enelist[15], rsim.enelist[9], rsim.enelist[19], rsim.enelist[30], rsim.enelist[31], rsim.enelist[22], rsim.enelist[23], rsim.enelist[6], rsim.enelist[16], rsim.enelist[32], rsim.enelist[33], rsim.enelist[7], rsim.enelist[17], rsim.enelist[24], rsim.enelist[25], rsim.enelist[35], rsim.enelist[37], rsim.enelist[36]);

	std::string formattedString = buffer;
	outputFile << formattedString;
	outputFile.flush();

	outputFile.close();
	return true;
}

bool RNArefineFunc::judgetermination(map<int, int> tempidx2repidx, vector<double> energies){
	int i;
	int replicanum=tempidx2repidx.size();
	int lowrepline=int(replicanum*0.75);

	double mean = 0.0;
    double variance = 0.0;
	double stdval = 0.0;

	double sum = 0.0;
    double sumSquaredDiff = 0.0;
    int count = 0;
	double value;

	for(i=0;i<replicanum;i++){
		if(i>=lowrepline){
			value=energies[tempidx2repidx[i]];
			sum += value;
			sumSquaredDiff += value * value;
			count+=1;
		}
		//eneall+=energies[tempidx2repidx[i]];
	}

	mean = sum / count;
	variance = sumSquaredDiff / count - mean * mean;
	stdval = sqrt(variance);
	
	if(stdval < mean*0.0001) return true;
	return false;
}

void RNArefineFunc::swapmovement(map<int, int> &repidx2tempidx, map<int, int> &tempidx2repidx, vector<double> Temperatures, vector<double> energies){
	int i;
	int rep1, rep2;
	int temp1, temp2;
	double dene, dbeta, p;
	double randval=Random();

	if(randval > 0.5){
		for(i=0;i<int(energies.size()/2);i++){
			temp1=2*i;
			temp2=2*i+1;
			rep1=tempidx2repidx[temp1];
			rep2=tempidx2repidx[temp2];
	
			dene=energies[rep2]-energies[rep1];
			dbeta=(1/Temperatures[temp2])-(1/Temperatures[temp1]);
			p=exp(dene*dbeta);
			randval=Random();
			if(p>randval){
				repidx2tempidx[rep1]=temp2;
				repidx2tempidx[rep2]=temp1;
				tempidx2repidx[temp1]=rep2;
				tempidx2repidx[temp2]=rep1;
			}
		}
	}else{
		for(i=1;i<int(energies.size()/2);i++){
			temp1=2*i-1;
			temp2=2*i;
			rep1=tempidx2repidx[temp1];
			rep2=tempidx2repidx[temp2];

			dene=energies[rep2]-energies[rep1];
			dbeta=(1/Temperatures[temp2])-(1/Temperatures[temp1]);
			p=exp(dene*dbeta);
			randval=Random();
			if(p>randval){
				repidx2tempidx[rep1]=temp2;
				repidx2tempidx[rep2]=temp1;
				tempidx2repidx[temp1]=rep2;
				tempidx2repidx[temp2]=rep1;
			}
		}
	}
}

/*
void RNArefineFunc::repsingleprocess(vector<Residueatoms> &decstr, RNASimulation &rsim, ofstream &outputFile, int &ms, int &fnow, int &dupd, int cycleidx, int repidx, string refinestr, int numbb, vector<int> bbind, int tempidx, double oldenergy, double newenergy, clock_t remcstart, clock_t remcfinish){
	int movetype;
	bool flagps;
	char buffer[1024];
	
	rsim.curcycle=k;
	flagps=rsim.fmcfragsweepcom(decstr,numbb,oldenergy,newenergy,movetype,Temperatures[tempidx]);
	rsim.baseenecoeff=rsim.baseenecoeff*1.0001;

	if(flagps){
		ms+=1;
		//refinestr=outfile+"/replicas/replica_"+to_string(i+1)+"/decoy.pdb";
		rsim.writepdb(decstr,numbb,bbind,refinestr,0);
		//StrVis(exepath, refinestr, outfile+"/replicas/init_decoy.pdb", outfile+"/replicas/replica_"+to_string(i+1)+"/image", ms[i]);
		oldenergy=newenergy;
		fnow=0;
	}
	else{
		fnow++;
		dupd++;
	}
	remcfinish=clock(); 
	printf("doing cycle %10d replica %10d acc %d mov %2d eff %3d fail %3d cfail %3d time %4d oene %10.3f nene %10.3f\n",
		cycleidx,repidx,flagps,movetype,ms,dupd,fnow,int((remcfinish-remcstart)/CLOCKS_PER_SEC),oldenergy,newenergy);

	snprintf(buffer, sizeof(buffer), "doing cycle %10d eff %3d acc %d fail %3d cfail %3d time %4d oene %10.3f nene %10.3f mov %2d fapeinit %10.3f fapedeeptmp %10.3f fapetmp %10.3f clashout %10.3f clashloss %10.3f lengthout %10.3f lengthloss %10.3f adjlengthout %10.3f adjlengthloss %10.3f angleout %10.3f angleloss %10.3f torsionout %10.3f torsionloss %10.3f baseout %10.3f baseloss %10.3f\n",
		cycleidx,ms,flagps,dupd,fnow,int((remcfinish-remcstart)/CLOCKS_PER_SEC),oldenergy,newenergy,movetype,
		rsim.enelist[0], rsim.enelist[1], rsim.enelist[18], rsim.enelist[2], rsim.enelist[12], rsim.enelist[3], rsim.enelist[13], rsim.enelist[4], rsim.enelist[14],
		rsim.enelist[5], rsim.enelist[15], rsim.enelist[6], rsim.enelist[16], rsim.enelist[7], rsim.enelist[17]);
	std::string formattedString = buffer;
	outputFile << formattedString;
	outputFile.flush();
}
*/

void threadprocess(int start, int end, vector<RNASimulation> &rsims, vector<vector<Residueatoms>> &decstrs, vector<ofstream> &outputFiles, vector<int> &ms, vector<int> &fnows, vector<int> &dupds, vector<double> &oldenergies, vector<double> &newenergies, string outfile, vector<double> Temperatures, map<int, int> repidx2tempidx, int cycleidx, int numbb, vector<int> bbind, std::chrono::high_resolution_clock::time_point remcstart){
//void threadprocess(int i, vector<RNASimulation> &rsims, vector<vector<Residueatoms>> &decstrs, vector<ofstream> &outputFiles, vector<int> &ms, vector<int> &fnows, vector<int> &dupds, vector<double> &oldenergies, vector<double> &newenergies, string outfile, vector<double> Temperatures, map<int, int> repidx2tempidx, int cycleidx, int numbb, vector<int> bbind, clock_t remcstart, clock_t remcfinish){
	int movetype=-1;
	bool flagps;
	char buffer[1024];
	string refinestr;
	double tmp;

	for (int i = start; i < end; ++i)
	{
		rsims[i].curcycle=cycleidx;
		flagps=rsims[i].fmcfragsweepcom(decstrs[i],numbb,oldenergies[i],newenergies[i],movetype,Temperatures[repidx2tempidx[i]]);
		tmp = oldenergies[i];

		if(flagps){
			ms[i]+=1;
			refinestr=outfile+"/replicas/replica_"+to_string(i+1)+"/structure.pdb";
			rsims[i].writepdb(decstrs[i],numbb,bbind,refinestr,0);
			//StrVis(exepath, refinestr, outfile+"/replicas/init_decoy.pdb", outfile+"/replicas/replica_"+to_string(i+1)+"/image", ms[i]);
			oldenergies[i]=newenergies[i];
			fnows[i]=0;
		}
		else{
			fnows[i]++;
			dupds[i]++;
		}
		
		auto remcfinish=std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed_seconds = remcfinish - remcstart;
		int elapsed_seconds_int = static_cast<int>(elapsed_seconds.count());

		printf("doing cycle %10d replica %10d acc %d mov %2d eff %3d fail %3d cfail %3d time %4d oene %10.3f nene %10.3f templidx %4d tmpl %10.3f\n",
			cycleidx+1,i+1,flagps,movetype,ms[i],dupds[i],fnows[i],elapsed_seconds_int,tmp,newenergies[i],repidx2tempidx[i]+1, Temperatures[repidx2tempidx[i]]);


		//snprintf(buffer, sizeof(buffer), "doing cycle %10d eff %3d acc %d fail %3d cfail %3d time %4d oene %10.3f nene %10.3f mov %2d templidx %4d tmpl %10.3f fapeinit %10.3f fapedeeptmp %10.3f contact %10.3f clashout %10.3f clashloss %10.3f lengthout %10.3f lengthloss %10.3f angleout %10.3f angleloss %10.3f torsionout %10.3f torsionloss %10.3f addlengthout %10.3f addlengthloss %10.3f addangleout %10.3f addangleloss %10.3f addtorsionout %10.3f addtorsionloss %10.3f baseout %10.3f baseloss %10.3f stackout %10.3f stackloss %10.3f hbond %10.3f\n",
			//cycleidx+1,ms[i],flagps,dupds[i],fnows[i],elapsed_seconds_int,tmp,newenergies[i],movetype, repidx2tempidx[i]+1, Temperatures[repidx2tempidx[i]],
			//rsims[i].enelist[0], rsims[i].enelist[1], rsims[i].enelist[8], rsims[i].enelist[2], rsims[i].enelist[12], rsims[i].enelist[3], rsims[i].enelist[13], rsims[i].enelist[4], rsims[i].enelist[14], rsims[i].enelist[5], rsims[i].enelist[15], rsims[i].enelist[19], rsims[i].enelist[20], rsims[i].enelist[21], rsims[i].enelist[22], rsims[i].enelist[23], rsims[i].enelist[24],
			//rsims[i].enelist[6], rsims[i].enelist[16], rsims[i].enelist[7], rsims[i].enelist[17], rsims[i].enelist[18]);

		
		if(flagps){
			snprintf(buffer, sizeof(buffer), "doing cycle %10d eff %3d acc %d fail %3d cfail %3d time %4d oene %10.3f nene %10.3f mov %2d templidx %4d tmpl %10.3f fapeinit %10.3f fapedeeptmp %10.3f contact8A %10.3f contact20A %10.3f clashout %10.3f clashloss %10.3f adjpairlengthout %10.3f adjpairlengthloss %10.3f ntpairlengthout %10.3f ntpairlengthloss %10.3f lengthout %10.3f lengthloss %10.3f angleout %10.3f angleloss %10.3f torsionout %10.3f torsionloss %10.3f multitorsionout %10.3f multitorsionloss %10.3f pairedtorsionloss %10.3f pairedbbtorsionloss %10.3f clusteredtorsionout %10.3f clusteredtorsionloss %10.3f baseout %10.3f baseloss %10.3f bphydout %10.3f bphydloss %10.3f stackout %10.3f stackloss %10.3f puckerout %10.3f puckerloss %10.3f singlenttemploss %10.3f bptemploss %10.3f stackingtemploss %10.3f\n",
				cycleidx+1,ms[i],flagps,dupds[i],fnows[i],elapsed_seconds_int,tmp,newenergies[i],movetype, repidx2tempidx[i]+1, Temperatures[repidx2tempidx[i]],
				rsims[i].enelist[0], rsims[i].enelist[1], rsims[i].enelist[8], rsims[i].enelist[18], rsims[i].enelist[2], rsims[i].enelist[12], rsims[i].enelist[26], rsims[i].enelist[27], rsims[i].enelist[28], rsims[i].enelist[29], rsims[i].enelist[3], rsims[i].enelist[13], rsims[i].enelist[4], rsims[i].enelist[14], rsims[i].enelist[5], rsims[i].enelist[15], rsims[i].enelist[9], rsims[i].enelist[19], rsims[i].enelist[30], rsims[i].enelist[31], rsims[i].enelist[22], rsims[i].enelist[23],rsims[i].enelist[6], rsims[i].enelist[16],rsims[i].enelist[32], rsims[i].enelist[33], rsims[i].enelist[7], rsims[i].enelist[17], rsims[i].enelist[24], rsims[i].enelist[25], rsims[i].enelist[35], rsims[i].enelist[37], rsims[i].enelist[36]);
			std::string formattedString = buffer;
			outputFiles[i] << formattedString;
			outputFiles[i].flush();
		}
	}
}

void enethreadprocess(int start, int end, vector<RNASimulation> &rsims, vector<vector<Residueatoms>> &decstrs, vector<ofstream> &outputFiles, vector<int> &ms, vector<int> &fnows, vector<int> &dupds, vector<double> &oldenergies, vector<double> &newenergies, vector<bool> &fidx, int &finum, int fnum, string outfile, int cycleidx, int numbb, vector<int> bbind, std::chrono::high_resolution_clock::time_point remcstart){
//void threadprocess(int i, vector<RNASimulation> &rsims, vector<vector<Residueatoms>> &decstrs, vector<ofstream> &outputFiles, vector<int> &ms, vector<int> &fnows, vector<int> &dupds, vector<double> &oldenergies, vector<double> &newenergies, string outfile, vector<double> Temperatures, map<int, int> repidx2tempidx, int cycleidx, int numbb, vector<int> bbind, clock_t remcstart, clock_t remcfinish){
	int movetype;
	bool flagps;
	char buffer[1024];
	string refinestr;
	double Bt=-1;
	double tmp;

	for (int i = start; i < end; ++i)
	{
		if(fidx[i]) continue;
		rsims[i].curcycle=cycleidx;
		flagps=rsims[i].fmcfragsweepcom(decstrs[i],numbb,oldenergies[i],newenergies[i],movetype,Bt);
		
		tmp=oldenergies[i];
		
		if(flagps){
			ms[i]+=1;
			refinestr=outfile+"/enemin/replica_"+to_string(i+1)+"/structure.pdb";
			rsims[i].writepdb(decstrs[i],numbb,bbind,refinestr,0);
			//StrVis(exepath, refinestr, outfile+"/replicas/init_decoy.pdb", outfile+"/replicas/replica_"+to_string(i+1)+"/image", ms[i]);
			oldenergies[i]=newenergies[i];
			fnows[i]=0;
		}
		else{
			fnows[i]++;
			dupds[i]++;
		}
		if(fnows[i]>=fnum){
			fidx[i]=true;
			finum+=1;
		}
		
		auto remcfinish=std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed_seconds = remcfinish - remcstart;
		int elapsed_seconds_int = static_cast<int>(elapsed_seconds.count());

		printf("doing cycle %10d replica %10d acc %d mov %2d eff %3d fail %3d cfail %3d time %4d oene %10.3f nene %10.3f\n",
			cycleidx+1,i+1,flagps,movetype,ms[i],dupds[i],fnows[i],elapsed_seconds_int,tmp,newenergies[i]);

		//snprintf(buffer, sizeof(buffer), "doing cycle %10d eff %3d acc %d fail %3d cfail %3d time %4d oene %10.3f nene %10.3f mov %2d fapeinit %10.3f fapedeeptmp %10.3f contact %10.3f clashout %10.3f clashloss %10.3f lengthout %10.3f lengthloss %10.3f angleout %10.3f angleloss %10.3f torsionout %10.3f torsionloss %10.3f addlengthout %10.3f addlengthloss %10.3f addangleout %10.3f addangleloss %10.3f addtorsionout %10.3f addtorsionloss %10.3f baseout %10.3f baseloss %10.3f stackout %10.3f stackloss %10.3f hbond %10.3f\n",
			//cycleidx+1,ms[i],flagps,dupds[i],fnows[i],elapsed_seconds_int,tmp,newenergies[i],movetype,
			//rsims[i].enelist[0], rsims[i].enelist[1], rsims[i].enelist[8], rsims[i].enelist[2], rsims[i].enelist[12], rsims[i].enelist[3], rsims[i].enelist[13], rsims[i].enelist[4], rsims[i].enelist[14], rsims[i].enelist[5], rsims[i].enelist[15], rsims[i].enelist[19], rsims[i].enelist[20], rsims[i].enelist[21], rsims[i].enelist[22], rsims[i].enelist[23], rsims[i].enelist[24],
			//rsims[i].enelist[6], rsims[i].enelist[16], rsims[i].enelist[7], rsims[i].enelist[17], rsims[i].enelist[18]);

		
		if(flagps){
			snprintf(buffer, sizeof(buffer), "doing cycle %10d eff %3d acc %d fail %3d cfail %3d time %4d oene %10.3f nene %10.3f mov %2d fapeinit %10.3f fapedeeptmp %10.3f contact8A %10.3f contact20A %10.3f clashout %10.3f clashloss %10.3f adjpairlengthout %10.3f adjpairlengthloss %10.3f ntpairlengthout %10.3f ntpairlengthloss %10.3f lengthout %10.3f lengthloss %10.3f angleout %10.3f angleloss %10.3f torsionout %10.3f torsionloss %10.3f multitorsionout %10.3f multitorsionloss %10.3f pairedtorsionloss %10.3f pairedbbtorsionloss %10.3f clusteredtorsionout %10.3f clusteredtorsionloss %10.3f baseout %10.3f baseloss %10.3f bphydout %10.3f bphydloss %10.3f stackout %10.3f stackloss %10.3f puckerout %10.3f puckerloss %10.3f singlenttemploss %10.3f bptemploss %10.3f stackingtemploss %10.3f\n",
				cycleidx+1,ms[i],flagps,dupds[i],fnows[i],elapsed_seconds_int,tmp,newenergies[i],movetype,
				rsims[i].enelist[0], rsims[i].enelist[1], rsims[i].enelist[8], rsims[i].enelist[18], rsims[i].enelist[2], rsims[i].enelist[12], rsims[i].enelist[26], rsims[i].enelist[27], rsims[i].enelist[28], rsims[i].enelist[29], rsims[i].enelist[3], rsims[i].enelist[13], rsims[i].enelist[4], rsims[i].enelist[14], rsims[i].enelist[5], rsims[i].enelist[15], rsims[i].enelist[9], rsims[i].enelist[19], rsims[i].enelist[30], rsims[i].enelist[31], rsims[i].enelist[22], rsims[i].enelist[23], rsims[i].enelist[6], rsims[i].enelist[16], rsims[i].enelist[32], rsims[i].enelist[33], rsims[i].enelist[7], rsims[i].enelist[17], rsims[i].enelist[24], rsims[i].enelist[25], rsims[i].enelist[35], rsims[i].enelist[37], rsims[i].enelist[36]);
			std::string formattedString = buffer;
			outputFiles[i] << formattedString;
			outputFiles[i].flush();
		}
	}
}

void RNArefineFunc::rnarefine(string exepath, string targetstrfile, vector<string>templatefiles, double templatecoeff, vector<double>templateweights, string bppath, string stpath, string outfile, string paramdir, int reftype, double cuttime1, double cuttime2)
{
	bool flagps;
	int i,j,k;
	double randval;

	//variables in ss refinment
	string refinestr;
	map<int, int> initsslist, bplist, commonsslist;
	//vector<map<int, int>> tmpinitsslist, tmpnewsslist;
	map<string, int> sstype;
	string oldbrk, newbrk, sselement;
	//vector<string> tmpoldbrk, tmpnewbrk, tmpoldeleres, tmpneweleres;
	map<int, int> idx2residx;
	map<int, int> residx2idx;
	string ssstr;

	//variables for target structure
	int modelidx=0;
	int chainidx=0;
	int tmodelidx=0;
	int thainidx=0;
	char targetchainid;
	string seq;
	int numbb;
	vector<Atom> proseq;
	vector<Promodel> promod;
	vector<Boneinfo> bb;
	vector<int> bbind;

	//variables for template structure
	vector<string> tseqs;
	vector<vector<Atom>> tproseq;
	vector<vector<Promodel>> tpromod;
	vector<vector<Boneinfo>> tbb;
	vector<vector<int>> tbbind;

	//objects for refinements
	LoadPDB  lpdb;
	//LoadST lst;
	//LoadBP lbp;
	LoadParams lparams;
	LoadDeepParams ldparams;
	ParseSeq tps;
	RNASimulation rsim;
	//LoadBaseTemplate lbt;
	vector<RNASimulation> rsims;

	//vector<BPBasicInfo> basepairs;
	//vector<STBasicInfo> stackingpairs;
	
	cout<<"Loading pdb file" << endl;
	flagps=lpdb.loadpdb(targetstrfile, proseq, promod);
	if(!flagps) return;
	targetchainid=promod[modelidx].procha[chainidx].chainid;
	lpdb.extractbb(proseq, promod, bb, modelidx, chainidx);

	numbb=bb.size();
	seq=tps.pdb2seq(bb);
	for(i=0;i<numbb;i++){
		residx2idx[bb[i].resind]=i;
		idx2residx[i]=bb[i].resind;
	}

	//cout<<"Loading stacking information" << endl;
	//ldparams.loadST(stpath);
	//cout<<"Loading base pair information" << endl;
	//ldparams.loadBP(bppath, seq);

	bbind.resize(numbb);
	for(i=0;i<numbb;i++) bbind[i]=bb[i].resind;

	cout<<"Loading templates" << endl;
	if(templatefiles.size() != 0){
		//tmpinitsslist.resize(templatefiles.size());
		//tmpnewsslist.resize(templatefiles.size());
		//tmpoldbrk.resize(templatefiles.size());
		//tmpnewbrk.resize(templatefiles.size());
		//tmpoldeleres.resize(templatefiles.size());
		//tmpneweleres.resize(templatefiles.size());
		tseqs.resize(templatefiles.size());
		tproseq.resize(templatefiles.size());
		tpromod.resize(templatefiles.size());
		tbb.resize(templatefiles.size());
		tbbind.resize(templatefiles.size());
	}

	for(i=0;i<templatefiles.size();i++)
	{
		flagps=lpdb.loadpdb(templatefiles[i], tproseq[i], tpromod[i]);
		if(!flagps) return;
		lpdb.extractbb(tproseq[i], tpromod[i], tbb[i], tmodelidx, thainidx);
		tseqs[i] = tps.pdb2seq(tbb[i]);
		tbbind[i].resize(tbb[i].size());
		for(j=0;j<tbb[i].size();j++) tbbind[i][j]=tbb[i][j].resind;
	}


	//rsim.fapesingleref(proseq,bb,numbb,bb,numbb,seq.length());
	rsim.distsingleref(proseq,bb,numbb,bb,numbb,seq.length());
	for(i=0;i<templatefiles.size();i++){
		//rsim.fapemultipleref(tproseq[i],tbb[i],tbb[i].size(),bb,numbb,seq.length(),i);
		rsim.distmultipleref(tproseq[i],tbb[i],tbb[i].size(),bb,numbb,seq.length(),templateweights[i],i);
	}
	rsim.coefftemplate=templatecoeff;
	
	//rsim.getdihemoveratio();

	//for(i=0;i<templatefiles.size();i++){
		//lparams.gettemplateparams(tproseq[i],tbb[i],tbb[i].size(),bb,numbb,seq.length(),i);
	//}
	

	cout<<"Loading pre-defined parameters" << endl;
	lparams.paramloadmain(paramdir,numbb,seq);


	cout<<"Loading predicted parameters" << endl;
	ldparams.paramloadmain(bppath, stpath ,numbb, seq);

	cout<<"Start refinement process" << endl;
	rsim.fsetinitialarray(numbb); //fapemultipleref以下に定義すること
	vector<Residueatoms> decstr(numbb);
	rsim.setStr(decstr, proseq, bb, seq);
        
	int simmode = 1;
	//calculating parameters includeing lengths, inner angles, and dihedral angles (sugar rings)
	rsim.coord2paramsugar(decstr,0,numbb,simmode);
	//calculating parameters includeing lengths, inner angles, and dihedral angles (backbone)
	rsim.coord2parambb(decstr,0,numbb,simmode);
	//calculating parameters includeing lengths, inner angles, and dihedral angles (additional parameters)
	rsim.coord2paramadd(decstr,0,numbb);
	//generating coordinates based on parameters (backbone)
	rsim.param2coordbb(decstr,0,numbb);
	//generating coordinates based on parameters (sugar)
	rsim.param2coordsugar(decstr,0,numbb);
	//generating coordinates based on parameters (base)
	rsim.param2coordbase(decstr,0,numbb);

	rsim.coord2paramhydrogen(decstr,0,numbb,2);
	rsim.param2coordhydrogen(decstr,0,numbb);
	
	//calculating parameters includeing lengths, inner angles, and dihedral angles (aditional geometrical parameters)
	//rsim.coord2addparam(decstr,0,numbb);
	rsim.writepdb(decstr,numbb,bbind,outfile+"/init.pdb",0);
	//exit(1);
	rsim.enetype=0;
	//rsim.getfape(decstr,numbb);
	//rsim.initfape=rsim.fapeenergysingle(decstr,numbb);
	double initenergy=rsim.fcalcrmsdenergy(decstr,numbb,rsim.lp1,rsim.rppmat,rsim.bpenergy,rsim.bsenergy,rsim.nmap,rsim.distmap,rsim.ntpairdist,rsim.adjpairdist,rsim.distene,rsim.multidistene,rsim.contactene,rsim.mapflag,rsim.mapbbflag,false);
	//outputene(rsim, outfile+"/initial_energy.csv", initenergy);
	//exit(1);
	
	auto remcstart=std::chrono::high_resolution_clock::now();
	auto remcfinish=std::chrono::high_resolution_clock::now();
	//variables for refinement
	double wtterm=20;
	if(wtterm<0) wtterm=0;
	else if(wtterm>100) wtterm=100;
	wtterm*=20000/100.0;
	double wtinit=wtterm/4.0;
	//連続でエネルギーが更新されなかった場合の終了条件を設定
	int fnum;
	fnum=500*numbb;
	//simulationの最大時間の設定（cuttimeの単位は時間）
	double timethresh1=cuttime1*3600;//時間の閾値(この時間を超えるとsimulationが修了される)
	double timethresh2=cuttime2*3600;//時間の閾値(この時間を超えるとsimulationが修了される)
	//cout<<timethresh1<<endl;
	//cout<<timethresh2<<endl;
	//timethresh=24*3600;
	//if(cuttime>0.0 && cuttime<24) timethresh=cuttime*3600;
	//simulationの最大回数の設定（配列の長さに反比例）
	//secondary structureの再予測の間隔
	int ssterm=1000;
	int outterm=100;
	char buffer[1024];
	int movetype;
	double elapsed_seconds_double;
	int lsweepnum, sweepnum;
	if(reftype == 0 || reftype == 1 || reftype == 2)
	{
		sweepnum=int(50000000/numbb);
		//sweepnum=5;
		fs::create_directories(outfile+"/replicas");
		rsim.writepdb(decstr,numbb,bbind,outfile+"/replicas/init_structure.pdb",0);

		int replicanum=40;
		//int swapinterval=int(30*std::pow(numbb,0.5));
		int swapinterval=int(10*std::pow(numbb,0.5));
		//cout<<swapinterval<<endl;
		vector<double> oldenergies(replicanum);
		vector<double> newenergies(replicanum);
		vector<vector<Residueatoms>> decstrs;
		vector<int> ms(replicanum,0);
		vector<int> fnows(replicanum,0);
		vector<int> dupds(replicanum,0);
		vector<double> Temperatures(replicanum);
		double Temperatureexp;
		map<int, int> tempidx2repidx;
		map<int, int> repidx2tempidx;
		vector<ofstream> outputFiles(replicanum);
		//double Tmax=(2.4+0.016*numbb)*10;
		//double Tmin=(0.6+0.00067*numbb)*10;
		//double Tmax=(2.4+0.016*numbb)*10;
		//double Tmin=(0.1+0.00067*numbb)*10;
		double Tmax=(1.2+0.032*numbb)*10;
		//double Tmin=(0.1+0.00067*numbb)*10;
		double Tmin=(0.05+0.00067*numbb)*10;

		for(i=0;i<replicanum;i++){
			rsims.push_back(rsim);
			decstrs.push_back(decstr);
			Temperatureexp=(double(replicanum)-double(i+1))/(double(replicanum)-1);
			Temperatures[i]=Tmin*std::pow(Tmax/Tmin, Temperatureexp);
			//fs::create_directories(outfile+"/replicas/replica_"+to_string(i+1)+"/image");
			fs::create_directories(outfile+"/replicas/replica_"+to_string(i+1));
			refinestr=outfile+"/replicas/replica_"+to_string(i+1)+"/structure.pdb";
			rsims[i].writepdb(decstr,numbb,bbind,refinestr,0);
			//StrVis(exepath, refinestr, outfile+"/replicas/init_decoy.pdb", outfile+"/replicas/replica_"+to_string(i+1)+"/image", 0);
			outputFiles[i].open(outfile+"/replicas/replica_"+to_string(i+1)+"/simulation_results.txt"); 
			if (!outputFiles[i].is_open()) {
				cout << "Error:: simulation results cannot be written.\n";
				exit(1);
			}
			oldenergies[i]=initenergy;
			tempidx2repidx[i]=i;
			repidx2tempidx[i]=i;
		}

		//for(i=0;i<replicanum;i++){
			//for(j=0;j<rsims[i].basepairs.size();j++){
				//cout<<to_string(j) + "-" + to_string(rsims[i].basepairs[j].basepidx.first) + "-" + to_string(rsims[i].basepairs[j].basepidx.second)<<endl;
			//}
		//}
		//exit(1);

		const int numThreads = 1;
		int step = int(replicanum / numThreads);
		remcstart=std::chrono::high_resolution_clock::now();
		//clock();
		for(k=0;k<sweepnum;k++)
		{
			std::vector<std::thread> threads;
			for (i = 0; i < numThreads; ++i) {
				int start = i * step;
				int end = (i == numThreads - 1) ? replicanum : (i + 1) * step;
				
				threads.emplace_back(threadprocess, start, end, std::ref(rsims), std::ref(decstrs), std::ref(outputFiles), std::ref(ms), std::ref(fnows), std::ref(dupds), std::ref(oldenergies), std::ref(newenergies), outfile, Temperatures, repidx2tempidx, k, numbb, bbind, remcstart);
				//threads.emplace_back(threadprocess, i, std::ref(rsims), std::ref(decstrs), std::ref(outputFiles), std::ref(ms), std::ref(fnows), std::ref(dupds), std::ref(oldenergies), std::ref(newenergies), outfile, Temperatures, repidx2tempidx, k, numbb, bbind, remcstart, remcfinish);
			}

			for (auto& thread : threads) thread.join();
			cout<<"\n"<<endl;

			remcfinish=std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed_seconds = remcfinish - remcstart;
			elapsed_seconds_double = elapsed_seconds.count();
			if(elapsed_seconds_double>timethresh1){
				cout<<"Terminate condition 1"<<endl;
				break;//too much time
			}
			if(judgetermination(tempidx2repidx,oldenergies) && k > 1000){
				cout<<"Terminate condition 2"<<endl;
				break;
			}

			if(k!=0 && k%swapinterval == 0){
				swapmovement(repidx2tempidx, tempidx2repidx, Temperatures, oldenergies);
			}
		}

		int minrepidx=-1;
		double minrepene=10000000000;
		for(i=0;i<replicanum;i++){
			outputFiles[i].close();
			rsims[i].writepdb(decstrs[i],numbb,bbind,outfile+"/replicas/replica_"+to_string(i+1)+"/final_structure.pdb",0);
		
			if(minrepene>oldenergies[i]){
				minrepene=oldenergies[i];
				minrepidx=i;
			}
		}
		cout<<"Final structure is the replica "+to_string(minrepidx+1)+" (energy is "+to_string(minrepene)+")"<<endl;
		rsims[minrepidx].writepdb(decstrs[minrepidx],numbb,bbind,outfile+"/replicas/final_structure_tmp.pdb",0);

		if(reftype == 0){
			vector<bool> fidx(replicanum,false);
			int finum=0;
			lsweepnum=int(50000000/numbb);
			//lsweepnum=5;

			for(i=0;i<replicanum;i++){
				fs::create_directories(outfile+"/enemin/replica_"+to_string(i+1));
				refinestr=outfile+"/enemin/replica_"+to_string(i+1)+"/structure.pdb";
				ms[i]=0;
				fnows[i]=0;
				dupds[i]=0;
				rsims[i].enetype=1;
				rsims[i].writepdb(decstrs[i],numbb,bbind,refinestr,0);
				//StrVis(exepath, refinestr, outfile+"/replicas/init_decoy.pdb", outfile+"/replicas/replica_"+to_string(i+1)+"/image", 0);
				outputFiles[i].open(outfile+"/enemin/replica_"+to_string(i+1)+"/simulation_results.txt"); 
				if (!outputFiles[i].is_open()) {
					cout << "Error:: simulation results cannot be written.\n";
					exit(1);
				}
				oldenergies[i]=rsims[i].fcalcrmsdenergy(decstrs[i],numbb,rsims[i].lp1,rsims[i].rppmat,rsims[i].bpenergy,rsims[i].bsenergy,rsims[i].nmap,rsims[i].distmap,rsims[i].ntpairdist,rsims[i].adjpairdist,rsims[i].distene,rsims[i].multidistene,rsims[i].contactene,rsims[i].mapflag,rsims[i].mapbbflag,false);
				newenergies[i]=oldenergies[i];
			}

			const int numThreads = 1;
			int step = int(replicanum / numThreads);
			remcstart=std::chrono::high_resolution_clock::now();
			for(k=0;k<lsweepnum;k++)
			{
				std::vector<std::thread> threads;
				for (i = 0; i < numThreads; ++i) {
					int start = i * step;
					int end = (i == numThreads - 1) ? replicanum : (i + 1) * step;
					
					threads.emplace_back(enethreadprocess, start, end, std::ref(rsims), std::ref(decstrs), std::ref(outputFiles), std::ref(ms), std::ref(fnows), std::ref(dupds), std::ref(oldenergies), std::ref(newenergies), std::ref(fidx), std::ref(finum), fnum, outfile, k, numbb, bbind, remcstart);
				}

				for (auto& thread : threads) thread.join();
				cout<<"\n"<<endl;

				remcfinish=std::chrono::high_resolution_clock::now();
				std::chrono::duration<double> elapsed_seconds = remcfinish - remcstart;
				elapsed_seconds_double = elapsed_seconds.count();
				if(elapsed_seconds_double>timethresh2) break;//too much time
				if(finum == replicanum) break;
			}

			int minrepidx=-1;
			double minrepene=10000000000;
			for(i=0;i<replicanum;i++){
				outputFiles[i].close();
				rsims[i].writepdb(decstrs[i],numbb,bbind,outfile+"/enemin/replica_"+to_string(i+1)+"/final_structure.pdb",0);
			
				if(minrepene>oldenergies[i]){
					minrepene=oldenergies[i];
					minrepidx=i;
				}
			}
			cout<<"Final structure is the replica "+to_string(minrepidx+1)+" (energy is "+to_string(minrepene)+")"<<endl;
			rsims[minrepidx].writepdb(decstrs[minrepidx],numbb,bbind,outfile+"/enemin/final_structure_tmp.pdb",0);
			rsims[minrepidx].writepdb(decstrs[minrepidx],numbb,bbind,outfile+"/final_structure.pdb",0);
		}

		if(reftype == 1){
			rsim=rsims[minrepidx];
			decstr.clear();
			decstr.shrink_to_fit();
			decstr.resize(decstrs[minrepidx].size());
			std::copy(decstrs[minrepidx].begin(), decstrs[minrepidx].end(), decstr.begin());
		}
	}


	if(reftype == 1 || reftype == 3){
		lsweepnum=int(50000000/numbb);
		int m=0;
		int dupd=0;
		double Bt=-1;
		double oldenergy, newenergy;
		int fnow=0;
		string formattedString;
		fs::create_directories(outfile+"/enemin");
		ofstream outputFile;
		outputFile.open(outfile+"/enemin/simulation_results.txt"); 
		if (!outputFile.is_open()) {
			cout << "Error:: simulation results cannot be written.\n";
			exit(1);
		}

		rsim.enetype=1;	
		oldenergy=rsim.fcalcrmsdenergy(decstr,numbb,rsim.lp1,rsim.rppmat,rsim.bpenergy,rsim.bsenergy,rsim.nmap,rsim.distmap,rsim.ntpairdist,rsim.adjpairdist,rsim.distene,rsim.multidistene,rsim.contactene,rsim.mapflag,rsim.mapbbflag,false);
		
		refinestr=outfile+"/enemin/init_structure.pdb";
		rsim.writepdb(decstr,numbb,bbind,refinestr,0);

		int index=0;
		double Tmax, Tmin;
		double u=0.005;
		double d=0.003;
		double Tmax0=50.0;
		double Tmin0=20.0;
		double ramda = 0.00015;
		rsim.TK=50.0;

		remcstart=std::chrono::high_resolution_clock::now();
		for(k=0;k<lsweepnum;k++)
		{
			Tmax=Tmax0*exp(-ramda*k); 
			Tmin=Tmin0*exp(-ramda*k);

			//cout<<k << "-"<< Tmax << "-" << Tmin<< "-" << index << "-" << rsim.TK <<endl;

			flagps=rsim.fmcfragsweepcom(decstr,numbb,oldenergy,newenergy,movetype,Bt);

			//if(abs(Tmax-Tmin) < 0.1) index=-1;
			if(rsim.TK < Tmin) index=1;
        	else if(rsim.TK > Tmax) index=0;

			//if(index == -1){
				//rsim.TK=-1;
			//}else 
			if(index == 1){
				rsim.TK += u;
			}else if(index == 0){
				rsim.TK -= d;
			}

			if(newenergy < oldenergy){
				m++;
				refinestr=outfile+"/enemin/refine_structure.pdb";
				rsim.writepdb(decstr,numbb,bbind,refinestr,0);
			}
			if(flagps)
			{
				oldenergy=newenergy;
				fnow=0;
			}
			else{
				fnow++;
				dupd++;
			}
			if(fnow>=fnum) break;
			remcfinish=std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed_seconds = remcfinish - remcstart;
			elapsed_seconds_double = elapsed_seconds.count();

			printf("doing cycle %10d mov %2d eff %3d fail %3d cfail %3d time %4d oene %10.3f nene %10.3f\n",
					k,movetype,m,dupd,fnow,static_cast<int>(elapsed_seconds_double),oldenergy,newenergy);
			
			//snprintf(buffer, sizeof(buffer), "doing cycle %10d eff %3d fail %3d cfail %3d time %4d oene %10.3f nene %10.3f mov %2d fapeinit %10.3f fapedeeptmp %10.3f contact %10.3f clashout %10.3f clashloss %10.3f lengthout %10.3f lengthloss %10.3f angleout %10.3f angleloss %10.3f torsionout %10.3f torsionloss %10.3f addlengthout %10.3f addlengthloss %10.3f addangleout %10.3f addangleloss %10.3f addtorsionout %10.3f addtorsionloss %10.3f baseout %10.3f baseloss %10.3f stackout %10.3f stackloss %10.3f hbond %10.3f\n",
					//k,m,dupd,fnow,static_cast<int>(elapsed_seconds_double),oldenergy,newenergy,movetype,
					//rsim.enelist[0], rsim.enelist[1], rsim.enelist[8], rsim.enelist[2], rsim.enelist[12], rsim.enelist[3], rsim.enelist[13], rsim.enelist[4], rsim.enelist[14], rsim.enelist[5], rsim.enelist[15], rsim.enelist[19], rsim.enelist[20], rsim.enelist[21], rsim.enelist[22], rsim.enelist[23], rsim.enelist[24], 
					//rsim.enelist[6], rsim.enelist[16], rsim.enelist[7], rsim.enelist[17], rsim.enelist[18]);

			
			if(flagps)
			{
				snprintf(buffer, sizeof(buffer), "doing cycle %10d eff %3d fail %3d cfail %3d time %4d oene %10.3f nene %10.3f mov %2d fapeinit %10.3f fapedeeptmp %10.3f contact8A %10.3f contact20A %10.3f clashout %10.3f clashloss %10.3f adjpairlengthout %10.3f adjpairlengthloss %10.3f ntpairlengthout %10.3f ntpairlengthloss %10.3f lengthout %10.3f lengthloss %10.3f angleout %10.3f angleloss %10.3f torsionout %10.3f torsionloss %10.3f multitorsionout %10.3f multitorsionloss %10.3f pairedtorsionloss %10.3f pairedbbtorsionloss %10.3f clusteredtorsionout %10.3f clusteredtorsionloss %10.3f baseout %10.3f baseloss %10.3f bphydout %10.3f bphydloss %10.3f stackout %10.3f stackloss %10.3f puckerout %10.3f puckerloss %10.3f singlenttemploss %10.3f bptemploss %10.3f stackingtemploss %10.3f\n",
						k+1,m,dupd,fnow,static_cast<int>(elapsed_seconds_double),oldenergy,newenergy,movetype,
						rsim.enelist[0], rsim.enelist[1], rsim.enelist[8], rsim.enelist[18], rsim.enelist[2], rsim.enelist[12], rsim.enelist[26], rsim.enelist[27], rsim.enelist[28], rsim.enelist[29], rsim.enelist[3], rsim.enelist[13], rsim.enelist[4], rsim.enelist[14], rsim.enelist[5], rsim.enelist[15], rsim.enelist[9], rsim.enelist[19], rsim.enelist[30], rsim.enelist[31], rsim.enelist[22], rsim.enelist[23],  rsim.enelist[6], rsim.enelist[16],  rsim.enelist[32], rsim.enelist[33], rsim.enelist[7], rsim.enelist[17], rsim.enelist[24], rsim.enelist[25], rsim.enelist[35], rsim.enelist[37], rsim.enelist[36]);
				formattedString = buffer;
				outputFile << formattedString;
				outputFile.flush();
			}
			if(elapsed_seconds_double>timethresh2) break;//too much time
		}

		outputFile.close();
		rsim.writepdb(decstr,numbb,bbind,outfile+"/final_structure.pdb",0);
	}
}

RNArefineFunc::~RNArefineFunc()
{

}

	
