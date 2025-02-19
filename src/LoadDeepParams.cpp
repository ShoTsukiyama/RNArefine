// CallPython.cpp: implementation of the CallPython class.
//
//////////////////////////////////////////////////////////////////////
//#include <string>
//#include <map>
#include "LoadDeepParams.hpp"

using namespace std;
namespace fs = std::filesystem;

vector<vector<vector<double>>> CONTACTSCORE;
vector<vector<vector<int>>> CONTACTBINARY;

vector<BPBasicInfo> BASEPAIRS;
vector<STBasicInfo> STACKINGPAIRS;

//#include "ParseSS.hpp"
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

LoadDeepParams::LoadDeepParams()
{
}

bool LoadDeepParams::loadcontactmapfile(string filepath, vector<vector<double>> &contactscore, vector<vector<int>> &contactbinary)
{
	int count;
	string line, nt;
	bool flag1, flag2, flag3;
	int idx1, idx2;
	double sc;

	ifstream file(filepath);
	if (!file.is_open()) {
        cerr << "base paring information cannot be opened [" << filepath << "]" << endl;
        return false;
    }

	idx1=0;
	while (getline(file, line)) {
        stringstream ss(line);
        idx2=0;
		while (getline(ss, nt, ',')){
			try{
				contactscore[idx1][idx2]=stof(nt);
				if(contactscore[idx1][idx2] > CONTACT_THRESH) contactbinary[idx1][idx2]=1;
				//cout<<to_string(idx1) + "-" + to_string(idx2) + "-" + to_string(contactscore[idx1][idx2]) + "-" + to_string(contactbinary[idx1][idx2])<<endl;
			}catch (const exception& e) {
            	cerr << "Error processing line: " << line << " - " << e.what() << endl;
            	continue; // エラーが発生したらスキップ
       		}
			idx2+=1;
		}
		idx1+=1;
    }
	file.close();
	return true;
}

bool LoadDeepParams::loadcontactmap(string basepath1, string basepath2, int seqnum, vector<vector<vector<double>>> &contactscore, vector<vector<vector<int>>> &contactbinary)
{
	bool flag1, flag2;

	contactscore.resize(2, vector<vector<double>>(seqnum, vector<double>(seqnum, 0.0)));
	contactbinary.resize(2, vector<vector<int>>(seqnum, vector<int>(seqnum, 0)));

	//cutoff 8
	flag1 = loadcontactmapfile(basepath1, contactscore[0], contactbinary[0]);
	//cutoff 20
	flag2 = loadcontactmapfile(basepath2, contactscore[1], contactbinary[1]);

	if(flag1 && flag2) return true;
	else return false;
}

bool LoadDeepParams::loadST(string stpath, vector<STBasicInfo> &stackingpairs){
	int count;
	string line, nt;
	bool flag1, flag2, flag3;
	int idx1, idx2;
	double sc;
	STBasicInfo stacking;
	ifstream file(stpath);
	if (!file.is_open()) {
        cerr << "stacking information cannot be opened [" << stpath << "]" << endl;
        return false;
    }

	getline(file, line);
	while (getline(file, line)) {
        stringstream ss(line);
        flag1=false;
		flag2=false;
		flag3=false;        
		count = 0;
		while (getline(ss, nt, ',')){
			if(count==0){
				idx1=stoi(nt);
				flag1=true;
			}
			if(count==1){
				idx2=stoi(nt);
				flag2=true;
			}
			if(count==2){
				sc=stof(nt);
				flag3=true;
			}
			count+=1;
		}

		if(flag1 & flag2 & flag3){
			stacking.basepidx.first=idx1;
			stacking.basepidx.second=idx2;
			stacking.sc=sc;

			if((idx1 == idx2 + 1) | (idx2 == idx1 + 1)) stacking.sttype=1;
			else stacking.sttype=0;

			stackingpairs.push_back(stacking);
		}
    }
    
    // ファイルを閉じる
    file.close();

	return true;
}

bool LoadDeepParams::loadBP(string bppath, string seq, vector<BPBasicInfo> &basepairs){
	int count;
	string line, nt;
	bool flag1, flag2, flag3, flag4, flag5;
	int idx1, idx2;
	double sc;
	BPBasicInfo basepair;
	ifstream file(bppath);
	if (!file.is_open()) {
        cerr << "base paring information cannot be opened [" << bppath << "]" << endl;
        return false;
    }

	getline(file, line);
	while (getline(file, line)) {
        stringstream ss(line);
        flag1=false; 
		flag2=false; 
		flag3=false;
		count = 0;
		while (getline(ss, nt, ',')){
			if(count==0){
				idx1=stoi(nt);
				flag1=true;
			}
			if(count==1){
				idx2=stoi(nt);
				flag2=true;
			}
			if(count==2){
				sc=stof(nt);
				flag3=true;
			}
			count+=1;
		}

		if(flag1 & flag2 & flag3){
			if(idx1 < idx2){
				basepair.basepidx.first=idx1;
				basepair.basepidx.second=idx2;
				basepair.bptype=BaseTypeList[seq[idx1]][seq[idx2]];
				basepair.sc=sc;
			}else{
				basepair.basepidx.first=idx2;
				basepair.basepidx.second=idx1;
				basepair.bptype=BaseTypeList[seq[idx2]][seq[idx1]];
				basepair.sc=sc;
			}
			basepairs.push_back(basepair);
		}
    }

    
    // ファイルを閉じる
    file.close();

	return true;
}

//void LoadDeepParams::paramloadmain(string bppath, string stpath, string ct1path, string ct2path, int numbb, string seq)
void LoadDeepParams::paramloadmain(string bppath, string stpath, int numbb, string seq)
{
	bool flagps;
	
	//flagps=loadcontactmap(ct1path, ct2path, numbb, CONTACTSCORE, CONTACTBINARY);
	flagps=loadST(stpath, STACKINGPAIRS);
	flagps=loadBP(bppath, seq, BASEPAIRS);
}

LoadDeepParams::~LoadDeepParams()
{

}

