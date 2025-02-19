// CallPython.cpp: implementation of the CallPython class.
//
//////////////////////////////////////////////////////////////////////
//#include <string>
//#include <map>
#include "LoadParams.hpp"

using namespace std;
namespace fs = std::filesystem;


map< char, vector<Binnedparamtype>> TORSIONPARAM;
map< char, vector<Binnedparamtype>> INNERANGPARAM;
map< char, vector<Binnedparamtype>> LENGTHPARAM;
map< char, vector<Binnedparamtype>> SUBLENGTHPARAM;

map< char, vector<Binnedparamtype>> ADDINNERANGPARAM;
map< char, vector<Binnedparamtype>> ADDLENGTHPARAM;

map< char, vector<Binnedparamtype>> OPTORSIONPARAM;
map< char, vector<Binnedparamtype>> OPINNERANGPARAM;
map< char, vector<Binnedparamtype>> OPLENGTHPARAM;

map< char, vector<Binnedparamtype>> HYDTORSIONPARAM;
map< char, vector<Binnedparamtype>> HYDINNERANGPARAM;
map< char, vector<Binnedparamtype>> HYDLENGTHPARAM;

map<int, vector<Binnedparamtype>> BASEPAIRPARAM;
map<int, map<int, vector<Binnedparamtype>>> STACKINGPARAM;

vector<D3Binnedparamtype> BASEPAIRROT;
map<int, vector<D3Binnedparamtype>> STACKINGROT;

map<char, map<char, map<int, map<string, vector<vector<Binnedparamtype>>>>>> BASEPAIRHYDPARAM;

map< char, vector<D2Binnedparamtype>> MULTITORSIONPARAM;
map< char, vector<vector<D2Binnedparamtype>>> TORSIONPAIRSPARAM;
vector<D2Binnedparamtype> TORSIONBBPAIRSPARAM;

vector<vector<vector<Torsionclsparamtype>>> TORSIONLIST;

map<char, vector<vector<bool>>> EXPCLASHATOMSNT;
map<char, map<char, vector<vector<bool>>>> EXPCLASHATOMSADJ;
map<char, map<char, vector<vector<bool>>>> EXPCLASHATOMSANY;

map<char, vector<vector<double>>> BASETEMPLATE; 

vector<vector<Binnedparamtype>> PUCKERPARAM;

map<char, vector<vector<vector<PairLengInfo>>>> PAIRLENGNTPARAM;
map<char, map<char, vector<vector<vector<PairLengInfo>>>>> PAIRLENGADJPARAM;

map<char, vector<vector<vector<PairLengInfo>>>> BPATOMPAIRLENGNTPARAM;
map<char, map<char, vector<vector<vector<PairLengInfo>>>>> BPATOMPAIRLENGADJPARAM;
map<char, map<char, vector<vector<PairLengInfo>>>> BPATOMPAIRLENGPARAM;

map<char, vector<vector<PairLengInfo>>> ALLPAIRLENGNTPARAM;

map<char, vector<vector<TorClsInfo>>> TORCLSPARAM;

vector<double>VWPARAM;

map<int, map<char, vector<vector<Residueatoms>>>> SINGLENTTEMPLATE;
map<char, map<char, vector<vector<Residueatoms>>>> DILENTTEMPLATE;
map<char, map<char, vector<vector<Residueatoms>>>> BPLENTTEMPLATE;

map<int, map<char, map<char, vector<vector<Residueatoms>>>>> STACKTEMPLATE;

map<char, map<char, map<char, vector<TemplateDistMap>>>> STLENTTEMPLATELENGTH;
map<char, map<char, vector<TemplateDistMap>>> DILENTTEMPLATELENGTH;
map<char, map<char, vector<TemplateDistMap>>> BPLENTTEMPLATELENGTH;

map<char, vector<D2Binnedparamtype>> OXHYDBPARAM;

//#include "ParseSS.hpp"
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

LoadParams::LoadParams()
{
}

bool LoadParams::loadbaseparamtype1cutoff(string file, map<int, Cutoffparamtype> &cutoffparam)
{
	string	oneline;
	int count;
	int paramnum;
	double lin,lmiddle,lout,rin,rmiddle,rout;
	string field;
	ifstream inputFile(file);

	if(!inputFile.is_open())
	{
		printf( "Parameter file can not be opened (dihedral angle in base) %s\n", file.c_str());
		return false;
	}

	while(getline(inputFile, oneline))
	{
		stringstream ssbase(oneline);
		count=0;
		while (getline(ssbase, field, ',')) {
			if(count==0) paramnum=stoi(field);
			if(count==1) lout=stof(field);
			if(count==2) lmiddle=stof(field);
			if(count==3) lin=stof(field);
			if(count==4) rin=stof(field);
			if(count==5) rmiddle=stof(field);
			if(count==6) rout=stof(field);
			count+=1;
		}
		cutoffparam[paramnum].lout=lout;
		cutoffparam[paramnum].lmiddle=lmiddle;
		cutoffparam[paramnum].lin=lin;
		cutoffparam[paramnum].rin=rin;
		cutoffparam[paramnum].rmiddle=rmiddle;
		cutoffparam[paramnum].rout=rout;
	}
	return true;
}

bool LoadParams::loadbaseparamtype2cutoff(string file, map<int, Cutoffparamtype> &cutoffparam)
{
	string	oneline;
	int count;
	int paramnum;
	double lin,lmiddle,lout,rin,rmiddle,rout;
	string field;
	ifstream inputFile(file);

	if(!inputFile.is_open())
	{
		printf( "Parameter file can not be opened (dihedral angle in base) %s\n", file.c_str());
		return false;
	}

	while(getline(inputFile, oneline))
	{
		stringstream ssbase(oneline);
		count=0;
		while (getline(ssbase, field, ',')) {
			if(count==0) paramnum=stoi(field);
			if(count==1) lin=stof(field);
			if(count==2) lmiddle=stof(field);
			if(count==3) lout=stof(field);
			if(count==4) rout=stof(field);
			if(count==5) rmiddle=stof(field);
			if(count==6) rin=stof(field);
			count+=1;
		}
		cutoffparam[paramnum].lin=lin;
		cutoffparam[paramnum].lmiddle=lmiddle;
		cutoffparam[paramnum].lout=lout;
		cutoffparam[paramnum].rout=rout;
		cutoffparam[paramnum].rmiddle=rmiddle;
		cutoffparam[paramnum].rin=rin;
	}
	return true;
}

bool LoadParams::loadbasehbondcutoff(string file, map<int, Cutoffparamtype> &cutoffparam)
{
	string	oneline;
	int count;
	int paramnum;
	double lin,lmiddle,lout,rin,rmiddle,rout;
	string field;
	ifstream inputFile(file);

	if(!inputFile.is_open())
	{
		printf( "Parameter file can not be opened (dihedral angle in base) %s\n", file.c_str());
		return false;
	}

	while(getline(inputFile, oneline))
	{
		stringstream ssbase(oneline);
		count=0;
		while (getline(ssbase, field, ',')) {
			if(count==0) paramnum=stoi(field);
			if(count==1) lout=stof(field);
			if(count==2) lmiddle=stof(field);
			if(count==3) lin=stof(field);
			if(count==4) rin=stof(field);
			if(count==5) rmiddle=stof(field);
			if(count==6) rout=stof(field);
			count+=1;
		}
		cutoffparam[paramnum].lout=lout;
		cutoffparam[paramnum].lmiddle=lmiddle;
		cutoffparam[paramnum].lin=lin;
		cutoffparam[paramnum].rin=rin;
		cutoffparam[paramnum].rmiddle=rmiddle;
		cutoffparam[paramnum].rout=rout;
	}
	return true;
}

bool LoadParams::loadvwd(string datafile)
{
	int i;
	FILE *file;
	char	oneline[600];
	char nt;
	char atomname[5];
	atomname[0]=' ';atomname[1]=' ';atomname[2]=' ';atomname[3]=' ';atomname[4]='\0';
	string atomnamest;
	double val;
	int atomid;
	bool flagps;

	file = fopen(datafile.c_str(), "rt");
	if(file==NULL)
	{
		printf( "Parameter file can not be opened %s\n", datafile.c_str());
		return false;
	}	

	VWPARAM.resize(Atom2id.size());
	for(i=0; i<Atom2id.size();i++) VWPARAM[i]=-1;	

	while(!feof(file))
	{
		flagps=fgets(oneline,600,file);
		if(!flagps) break;
		sscanf(oneline,"%c%c%c%c%lf", &atomname[0], &atomname[1], &atomname[2], &atomname[3], &val);
		//sscanf(oneline,"%c%[^\n]%c%c%c%c%[^\n]%lf", &nt, &atomname[0], &atomname[1], &atomname[2], &atomname[3], &val);

		atomnamest=atomname;
		atomid=Atom2id[atomnamest];
		VWPARAM[atomid]=val;
	}
	return true;
}

bool LoadParams::loadbasebininfo(string bininfofile, Binnedparamtype &param)
{
	string	oneline;
	int ssidx, paramtype;
	double minval, maxval, span, binnum, mean, stdev;
	int count;
	string field;
	ifstream inputFile(bininfofile);

	if(!inputFile.is_open())
	{
		printf( "Parameter file can not be opened (dihedral angle in base) %s\n", bininfofile.c_str());
		return false;
	}

	while(getline(inputFile, oneline))
	{
		stringstream ssbase(oneline);
		count=0;
		while (getline(ssbase, field, ',')) {
			if(count==0) param.minval=stof(field); //minval
			if(count==1) param.maxval=stof(field); //maxval
			if(count==2) param.span=stof(field); //span
			if(count==3) param.binnum=stoi(field); //binnum
			if(count==4) param.mean=stof(field); //mean
			if(count==5) param.std=stof(field); //stdev
			count+=1;
		}
	}
	return true;
}

void LoadParams::loadparamvecd(string filepath, vector<double> &vec){

    // バイナリファイルからデータを読み込む
    std::ifstream input(filepath, std::ios::binary);
    if (!input) {
        std::cerr << "Error opening file (" + filepath + ")!" << std::endl;
        exit(1);
    }

    // ファイルサイズを取得
    input.seekg(0, std::ios::end);
    std::streamsize size = input.tellg();
    input.seekg(0, std::ios::beg);

    // バッファーを確保してファイルから読み込む
    vec.resize(size / sizeof(double));
    if (!input.read(reinterpret_cast<char*>(vec.data()), size)) {
        std::cerr << "Error reading file!" << std::endl;
        exit(1);
    }
}

void LoadParams::loadparamveci(string filepath, vector<int> &vec){

    // バイナリファイルからデータを読み込む
    std::ifstream input(filepath, std::ios::binary);
    if (!input) {
        std::cerr << "Error opening file (" + filepath + ")!" << std::endl;
        exit(1);
    }

    // ファイルサイズを取得
    input.seekg(0, std::ios::end);
    std::streamsize size = input.tellg();
    input.seekg(0, std::ios::beg);

    // バッファーを確保してファイルから読み込む
    vec.resize(size / sizeof(int));

	vector<double> vec_tmp;
	vec_tmp.resize(size / sizeof(double));
    if (!input.read(reinterpret_cast<char*>(vec_tmp.data()), size)) {
        std::cerr << "Error reading file!" << std::endl;
        exit(1);
    }

	for(int i=0;i<vec_tmp.size();i++) vec[i]=vec_tmp[i];
}

void LoadParams::loadparamvec2dd(string filepath, size_t dim1_size, size_t dim2_size, vector<vector<double>> &vec_2d){

    // バイナリファイルからデータを読み込む
    std::ifstream input(filepath, std::ios::binary);
    if (!input) {
        std::cerr << "Error opening file (" + filepath + ")!" << std::endl;
        exit(1);
    }

    // ファイルサイズを取得
    input.seekg(0, std::ios::end);
    std::streamsize size = input.tellg();
    input.seekg(0, std::ios::beg);

    // バッファーを確保してファイルから読み込む
    vec_2d.resize(dim1_size, vector<double>(dim2_size));
    //if (!input.read(reinterpret_cast<char*>(vec_2d.data()[0].data()), size)) {
        //std::cerr << "Error reading file!" << std::endl;
        //exit(1);
    //}

	for (size_t i = 0; i < dim1_size; ++i) {
        if (!input.read(reinterpret_cast<char*>(vec_2d[i].data()), dim2_size * sizeof(double))) {
            std::cerr << "Error reading file!" << std::endl;
            exit(1);
        }
    }
}

void LoadParams::loadparamvec2di(string filepath, size_t dim1_size, size_t dim2_size, vector<vector<int>> &vec_2d){

    // バイナリファイルからデータを読み込む
    std::ifstream input(filepath, std::ios::binary);
    if (!input) {
        std::cerr << "Error opening file (" + filepath + ")!" << std::endl;
        exit(1);
    }

    // ファイルサイズを取得
    input.seekg(0, std::ios::end);
    std::streamsize size = input.tellg();
    input.seekg(0, std::ios::beg);

	vector<vector<double>> vec_2d_tmp;
    // バッファーを確保してファイルから読み込む
    vec_2d.resize(dim1_size, vector<int>(dim2_size));
	vec_2d_tmp.resize(dim1_size, vector<double>(dim2_size));
    //if (!input.read(reinterpret_cast<char*>(vec_2d.data()[0].data()), size)) {
        //std::cerr << "Error reading file!" << std::endl;
        //exit(1);
    //}

	for (size_t i = 0; i < dim1_size; ++i) {
        if (!input.read(reinterpret_cast<char*>(vec_2d_tmp[i].data()), dim2_size * sizeof(double))) {
            std::cerr << "Error reading file!" << std::endl;
            exit(1);
        }
    }

	//vec_2d.resize(dim1_size, vector<int>(dim2_size));
	for(int i=0;i<vec_2d.size();i++){
		for(int j=0;j<vec_2d[i].size();j++){
			vec_2d[i][j]=vec_2d_tmp[i][j];
		}
	}

}

bool LoadParams::load2ddistribution(string binfilepath, string energyfilepath, string binaryfilepath, string distributionfilepath, vector<D2Binnedparamtype> &stparameter){
	string	oneline;
	int paramtype, paramnum;
	double minval1, maxval1, span1, binnum1, minval2, maxval2, span2, binnum2;
	int count;
	string field;

	paramnum=0;
	ifstream inputFile(binfilepath);
		
	if(!inputFile.is_open())
	{
		cout<<"Bin info file can not be opened (" + binfilepath + ")\n"<<endl;
		return false;
	}

	while(getline(inputFile, oneline))
	{
		stringstream ssbase(oneline);
		count=0;
		while (getline(ssbase, field, ','))
		{
			if(count==0) paramtype=stoi(field); //paramtype
			if(count==1) minval1=stof(field); //minval
			if(count==2) maxval1=stof(field); //maxval
			if(count==3) span1=stof(field); //span
			if(count==4) binnum1=stoi(field); //binnum
			if(count==5) minval2=stof(field); //minval
			if(count==6) maxval2=stof(field); //maxval
			if(count==7) span2=stof(field); //span
			if(count==8) binnum2=stoi(field); //binnum
			count+=1;
		}
		
		if(stparameter.size() < paramtype+1) stparameter.resize(paramtype+1);
			stparameter[paramtype].minval1=minval1;
			stparameter[paramtype].maxval1=maxval1;
			stparameter[paramtype].span1=span1;
			stparameter[paramtype].binnum1=binnum1;
			stparameter[paramtype].minval2=minval2;
			stparameter[paramtype].maxval2=maxval2;
			stparameter[paramtype].span2=span2;
			stparameter[paramtype].binnum2=binnum2;

			stparameter[paramtype].energy.resize(binnum1*binnum2);
			stparameter[paramtype].binary.resize(binnum1*binnum2);
			stparameter[paramtype].distr.resize(binnum1*binnum2);
			stparameter[paramtype].boundary1.resize(binnum1-1);
			stparameter[paramtype].center1.resize(binnum1);
			stparameter[paramtype].boundary2.resize(binnum2-1);
			stparameter[paramtype].center2.resize(binnum2);
		paramnum+=1;
	}
	

	vector<vector<double>> vec2denergies;
	loadparamvec2dd(energyfilepath, paramnum, binnum1*binnum2, vec2denergies);
	vector<vector<double>> vec2dbinary;
	loadparamvec2dd(binaryfilepath, paramnum, binnum1*binnum2, vec2dbinary);
	vector<vector<int>> vec2ddistribution;
	loadparamvec2di(distributionfilepath, paramnum, binnum1*binnum2, vec2ddistribution);

	for(int j=0;j<paramnum;j++)
	{
		stparameter[j].energy=vec2denergies[j];
		for(int k=0; k<vec2dbinary[j].size();k++) stparameter[j].binary[k]=vec2dbinary[j][k];
		for(int k=0; k<vec2ddistribution[j].size();k++) stparameter[j].distr[k]=vec2ddistribution[j][k];
		for(int k=0; k<binnum1-1; k++) stparameter[j].boundary1[k]=stparameter[j].minval1+k*stparameter[j].span1;
		for(int k=0; k<binnum2-1; k++) stparameter[j].boundary2[k]=stparameter[j].minval2+k*stparameter[j].span2;
		for(int k=0; k<binnum1; k++) stparameter[j].center1[k]=stparameter[j].minval1+(k-0.5)*stparameter[j].span1;
		for(int k=0; k<binnum2; k++) stparameter[j].center2[k]=stparameter[j].minval2+(k-0.5)*stparameter[j].span2;
	}

	return true;
}

bool LoadParams::loadramachandranplot(string binfilepath, string ramachandranpath, vector<D2Binnedparamtype> &stparameter){
	string	oneline;
	int paramtype, paramnum;
	double minval1, maxval1, span1, binnum1, minval2, maxval2, span2, binnum2;
	int count;
	string field;

	paramnum=0;
	ifstream inputFile(binfilepath);
		
	if(!inputFile.is_open())
	{
		cout<<"Bin info file can not be opened (" + binfilepath + ")\n"<<endl;
		return false;
	}

	while(getline(inputFile, oneline))
	{
		stringstream ssbase(oneline);
		count=0;
		while (getline(ssbase, field, ','))
		{
			if(count==0) paramtype=stoi(field); //paramtype
			if(count==1) minval1=stof(field); //minval
			if(count==2) maxval1=stof(field); //maxval
			if(count==3) span1=stof(field); //span
			if(count==4) binnum1=stoi(field); //binnum
			if(count==5) minval2=stof(field); //minval
			if(count==6) maxval2=stof(field); //maxval
			if(count==7) span2=stof(field); //span
			if(count==8) binnum2=stoi(field); //binnum
			count+=1;
		}
		
		if(stparameter.size() < paramtype+1) stparameter.resize(paramtype+1);
			stparameter[paramtype].minval1=minval1;
			stparameter[paramtype].maxval1=maxval1;
			stparameter[paramtype].span1=span1;
			stparameter[paramtype].binnum1=binnum1;
			stparameter[paramtype].minval2=minval2;
			stparameter[paramtype].maxval2=maxval2;
			stparameter[paramtype].span2=span2;
			stparameter[paramtype].binnum2=binnum2;

			stparameter[paramtype].energy.resize(binnum1*binnum2);
			stparameter[paramtype].binary.resize(binnum1*binnum2);
			stparameter[paramtype].distr.resize(binnum1*binnum2);
			stparameter[paramtype].boundary1.resize(binnum1-1);
			stparameter[paramtype].center1.resize(binnum1);
			stparameter[paramtype].boundary2.resize(binnum2-1);
			stparameter[paramtype].center2.resize(binnum2);
		paramnum+=1;
	}
	

	vector<vector<double>> vec2denergies;
	loadparamvec2dd(ramachandranpath, paramnum, binnum1*binnum2, vec2denergies);
	
	for(int j=0;j<paramnum;j++)
	{
		stparameter[j].energy=vec2denergies[j];
		for(int k=0; k<binnum1-1; k++) stparameter[j].boundary1[k]=stparameter[j].minval1+k*stparameter[j].span1;
		for(int k=0; k<binnum2-1; k++) stparameter[j].boundary2[k]=stparameter[j].minval2+k*stparameter[j].span2;
		for(int k=0; k<binnum1; k++) stparameter[j].center1[k]=stparameter[j].minval1+(k-0.5)*stparameter[j].span1;
		for(int k=0; k<binnum2; k++) stparameter[j].center2[k]=stparameter[j].minval2+(k-0.5)*stparameter[j].span2;
	}

	return true;
}

bool LoadParams::load3ddistribution(string binfilepath, string energyfilepath, string binaryfilepath, string distributionfilepath, vector<D3Binnedparamtype> &stparameter){
	string	oneline;
	int paramtype, paramnum;
	double minval1, maxval1, span1, binnum1, minval2, maxval2, span2, binnum2, minval3, maxval3, span3, binnum3;
	int count;
	string field;

	paramnum=0;
	ifstream inputFile(binfilepath);
		
	if(!inputFile.is_open())
	{
		cout<<"Bin info file can not be opened (" + binfilepath + ")\n"<<endl;
		return false;
	}

	while(getline(inputFile, oneline))
	{
		stringstream ssbase(oneline);
		count=0;
		while (getline(ssbase, field, ','))
		{
			if(count==0) paramtype=stoi(field); //paramtype
			if(count==1) minval1=stof(field); //minval
			if(count==2) maxval1=stof(field); //maxval
			if(count==3) span1=stof(field); //span
			if(count==4) binnum1=stoi(field); //binnum
			if(count==5) minval2=stof(field); //minval
			if(count==6) maxval2=stof(field); //maxval
			if(count==7) span2=stof(field); //span
			if(count==8) binnum2=stoi(field); //binnum
			if(count==9) minval3=stof(field); //minval
			if(count==10) maxval3=stof(field); //maxval
			if(count==11) span3=stof(field); //span
			if(count==12) binnum3=stoi(field); //binnum
			count+=1;
		}
		
		if(stparameter.size() < paramtype+1) stparameter.resize(paramtype+1);
			stparameter[paramtype].minval1=minval1;
			stparameter[paramtype].maxval1=maxval1;
			stparameter[paramtype].span1=span1;
			stparameter[paramtype].binnum1=binnum1;
			stparameter[paramtype].minval2=minval2;
			stparameter[paramtype].maxval2=maxval2;
			stparameter[paramtype].span2=span2;
			stparameter[paramtype].binnum2=binnum2;
			stparameter[paramtype].minval3=minval3;
			stparameter[paramtype].maxval3=maxval3;
			stparameter[paramtype].span3=span3;
			stparameter[paramtype].binnum3=binnum3;

			stparameter[paramtype].energy.resize(binnum1*binnum2*binnum3);
			stparameter[paramtype].binary.resize(binnum1*binnum2*binnum3);
			stparameter[paramtype].distr.resize(binnum1*binnum2*binnum3);
			stparameter[paramtype].boundary1.resize(binnum1-1);
			stparameter[paramtype].center1.resize(binnum1);
			stparameter[paramtype].boundary2.resize(binnum2-1);
			stparameter[paramtype].center2.resize(binnum2);
			stparameter[paramtype].boundary3.resize(binnum3-1);
			stparameter[paramtype].center3.resize(binnum3);
	}
	paramnum+=1;

	vector<vector<double>> vec2denergies;
	loadparamvec2dd(energyfilepath, paramnum, binnum1*binnum2*binnum3, vec2denergies);
	vector<vector<double>> vec2dbinary;
	loadparamvec2dd(binaryfilepath, paramnum, binnum1*binnum2*binnum3, vec2dbinary);
	vector<vector<int>> vec2ddistribution;
	loadparamvec2di(distributionfilepath, paramnum, binnum1*binnum2*binnum3, vec2ddistribution);

	for(int j=0;j<paramnum;j++)
	{
		stparameter[j].energy=vec2denergies[j];
		for(int k=0; k<vec2dbinary[j].size();k++) stparameter[j].binary[k]=vec2dbinary[j][k];
		for(int k=0; k<vec2ddistribution[j].size();k++) stparameter[j].distr[k]=vec2ddistribution[j][k];

		for(int k=0; k<binnum1-1; k++) stparameter[j].boundary1[k]=stparameter[j].minval1+k*stparameter[j].span1;
		for(int k=0; k<binnum2-1; k++) stparameter[j].boundary2[k]=stparameter[j].minval2+k*stparameter[j].span2;
		for(int k=0; k<binnum3-1; k++) stparameter[j].boundary3[k]=stparameter[j].minval3+k*stparameter[j].span3;
		for(int k=0; k<binnum1; k++) stparameter[j].center1[k]=stparameter[j].minval1+(k-0.5)*stparameter[j].span1;
		for(int k=0; k<binnum2; k++) stparameter[j].center2[k]=stparameter[j].minval2+(k-0.5)*stparameter[j].span2;
		for(int k=0; k<binnum3; k++) stparameter[j].center3[k]=stparameter[j].minval3+(k-0.5)*stparameter[j].span3;
	}
	return true;
}

bool LoadParams::loadstackingrot(string basepath, map<int, vector<D3Binnedparamtype>> &stparameter)
{
	map<int, string> sstype = {{0, "nonadj"},{1, "adj"}};
	for(int i=0;i<2;i++)
	{
		load3ddistribution(basepath+"/"+sstype[i]+"/rot3d/bin_info.txt", basepath+"/"+sstype[i]+"/rot3d/bin_energies.bin", basepath+"/"+sstype[i]+"/rot3d/bin_binary.bin", basepath+"/"+sstype[i]+"/rot3d/bin_distribution.bin", stparameter[i]);
	}
	return true;
}

bool LoadParams::loadbasepairrot(string basepath, vector<D3Binnedparamtype> &stparameter)
{
	
	load3ddistribution(basepath+"/all/rot3d/bin_info.txt", basepath+"/all/rot3d/bin_energies.bin", basepath+"/all/rot3d/bin_binary.bin", basepath+"/all/rot3d/bin_distribution.bin", stparameter);

	//vector<string> nts = {"A","C","G","U"};
	//for(int i=0;i<4;i++)
	//{
		//for(int j=i;j<4;j++)
		//{
			//load3ddistribution(basepath+"/"+nts[i]+"_"+nts[j]+"/rot3d/bin_info.txt", basepath+"/"+nts[i]+"_"+nts[j]+"/rot3d/bin_energies.bin", basepath+"/"+nts[i]+"_"+nts[j]+"/rot3d/bin_binary.bin", basepath+"/"+nts[i]+"_"+nts[j]+"/rot3d/bin_distribution.bin", stparameter[i]);
		//}
	//}
	return true;
}

bool LoadParams::loadstackingparameters(string basepath, map<int, map<int, vector<Binnedparamtype>>> &stparameter)
{
	string	oneline;
	int ssidx, paramtype, paramnum, binnum;
	double minval, maxval, span, mean, stdval, mode;
	int count;
	string field;

	map<int, string> sstype = {{0, "nonadj"},{1, "adj"}};
	map<int, string> paramname = {{0, "length"}, {1, "torsion"}, {2, "baseangle"}};

	for(int i=0;i<2;i++)
	{
		for(int j=0; j<paramname.size(); j++)
		{
			ifstream inputFile(basepath+"/"+sstype[i]+"/"+paramname[j]+"/bin_info.txt");
			if(!inputFile.is_open())
			{
				cout<<"Bin info file can not be opened (" + paramname[j] + ")\n"<<endl;
				return false;
			}

			paramnum=0;
			while(getline(inputFile, oneline))
			{
				stringstream ssbase(oneline);
				count=0;
				while (getline(ssbase, field, ','))
				{
					//if(count==0) ssidx=stoi(field); //ssidx
					if(count==0) paramtype=stoi(field); //paramtype
					if(count==1) minval=stof(field); //minval
					if(count==2) maxval=stof(field); //maxval
					if(count==3) span=stof(field); //span
					if(count==4) binnum=stoi(field); //binnum
					if(count==5) mean=stof(field); //mean
					if(count==6) stdval=stof(field); //stdval
					if(count==7) mode=stof(field); //stdval
					count+=1;
				}

				if(stparameter[i][j].size() < paramtype+1) stparameter[i][j].resize(paramtype+1);
				stparameter[i][j][paramtype].minval=minval;
				stparameter[i][j][paramtype].maxval=maxval;
				stparameter[i][j][paramtype].span=span;
				stparameter[i][j][paramtype].binnum=binnum;
				stparameter[i][j][paramtype].mean=mean;
				stparameter[i][j][paramtype].std=stdval;
				stparameter[i][j][paramtype].mode=mode;
				stparameter[i][j][paramtype].upper=mean+outlinerbsstd[paramname[j]]*stdval;
				stparameter[i][j][paramtype].lower=mean-outlinerbsstd[paramname[j]]*stdval;

				if(stparameter[i][j][paramtype].lower < 0) stparameter[i][j][paramtype].lower=0;

				//cout<< i << "-" << j << "-" << mean << "-" << mean << "-" << outlinerbsstd[paramname[j]]*stdval << "-" << stdval<<endl;
				//cout<< stparameter[i][j][paramtype].upper << "-" << stparameter[i][j][paramtype].lower<<endl;

				stparameter[i][j][paramtype].rangetype=0;
				if(paramname[j] == "torsion"){
					if(stparameter[i][j][paramtype].upper > 360){
						if(stparameter[i][j][paramtype].lower < 0){
							stparameter[i][j][paramtype].lower=0;
							stparameter[i][j][paramtype].upper=360;
						}else{
							stparameter[i][j][paramtype].upper = stparameter[i][j][paramtype].upper-360;
							stparameter[i][j][paramtype].rangetype=1;
						}
					}
					else if(stparameter[i][j][paramtype].lower < 0){
						stparameter[i][j][paramtype].lower = stparameter[i][j][paramtype].lower+360;
						stparameter[i][j][paramtype].rangetype=1;
					}
				}

				//cout<<to_string(stparameter[i][j][paramtype].minval) + "-" + to_string(stparameter[i][j][paramtype].maxval) + "-" + to_string(stparameter[i][j][paramtype].span) + "-" + to_string(stparameter[i][j][paramtype].binnum) + "-" + to_string(stparameter[i][j][paramtype].mean) + "-" + to_string(stparameter[i][j][paramtype].std) + "-" + to_string(stparameter[i][j][paramtype].mode)<<endl;

				stparameter[i][j][paramtype].energy.resize(binnum);
				stparameter[i][j][paramtype].binary.resize(binnum);
				stparameter[i][j][paramtype].distr.resize(binnum);
				stparameter[i][j][paramtype].boundary.resize(binnum-1);
				stparameter[i][j][paramtype].center.resize(binnum);
				paramnum+=1;
			}

			vector<vector<double>> vec2denergies;
			loadparamvec2dd(basepath+"/"+sstype[i]+"/"+paramname[j]+"/bin_energies.bin", paramnum, binnum, vec2denergies);
			vector<vector<double>> vec2dbinary;
			loadparamvec2dd(basepath+"/"+sstype[i]+"/"+paramname[j]+"/bin_binary.bin", paramnum, binnum, vec2dbinary);
			vector<vector<int>> vec2ddistribution;
			loadparamvec2di(basepath+"/"+sstype[i]+"/"+paramname[j]+"/bin_distribution.bin", paramnum, binnum, vec2ddistribution);
		
			for(int k=0;k<paramnum;k++){
				stparameter[i][j][k].energy=vec2denergies[k];

				stparameter[i][j][k].total=0;
				for(int l=0; l<vec2dbinary[k].size();l++){
					stparameter[i][j][k].binary[l]=vec2dbinary[k][l];
					stparameter[i][j][k].total+=int(1-stparameter[i][j][k].binary[l]);
				}

				stparameter[i][j][k].totaldist=0;
				for(int l=0; l<vec2ddistribution[k].size();l++){
					stparameter[i][j][k].distr[l]=vec2ddistribution[k][l];
					stparameter[i][j][k].totaldist+=vec2ddistribution[k][l];
				}
				for(int l=0; l<binnum-1; l++) stparameter[i][j][k].boundary[l]=stparameter[i][j][k].minval+l*stparameter[i][j][k].span;
				for(int l=0; l<binnum; l++) stparameter[i][j][k].center[l]=stparameter[i][j][k].minval+(l-0.5)*stparameter[i][j][k].span;
			}
		}
	}
	//exit(1);
	return true;
}

bool LoadParams::loadbasepairparameters(string basepath, map<int, vector<Binnedparamtype>> &stparameter)
{
	string	oneline;
	int ssidx, paramtype, paramnum, binnum;
	double minval, maxval, span, mean, stdval;
	int count;
	string field;

	map<int, string> paramname = {{0, "length"}, {1, "torsion"}};
	for(int j=0; j<paramname.size(); j++)
	{
		ifstream inputFile(basepath+"/all/"+paramname[j]+"/bin_info.txt");
		if(!inputFile.is_open())
		{
			cout<<"Bin info file can not be opened (" + paramname[j] + ")\n"<<endl;
			return false;
		}

		paramnum=0;
		while(getline(inputFile, oneline))
		{
			stringstream ssbase(oneline);
			count=0;
			while (getline(ssbase, field, ','))
			{
				//if(count==0) ssidx=stoi(field); //ssidx
				if(count==0) paramtype=stoi(field); //paramtype
				if(count==1) minval=stof(field); //minval
				if(count==2) maxval=stof(field); //maxval
				if(count==3) span=stof(field); //span
				if(count==4) binnum=stoi(field); //binnum
				if(count==5) mean=stof(field); //mean
				if(count==6) stdval=stof(field); //stdval
				count+=1;
			}

			if(stparameter[j].size() < paramtype+1) stparameter[j].resize(paramtype+1);
			stparameter[j][paramtype].minval=minval;
			stparameter[j][paramtype].maxval=maxval;
			stparameter[j][paramtype].span=span;
			stparameter[j][paramtype].binnum=binnum;
			stparameter[j][paramtype].mean=mean;
			stparameter[j][paramtype].std=stdval;


			stparameter[j][paramtype].upper=mean+outlinerbpstd[paramname[j]]*stdval;
			stparameter[j][paramtype].lower=mean-outlinerbpstd[paramname[j]]*stdval;

			stparameter[j][paramtype].rangetype=0;
			if(paramname[j] == "torsion"){
				if(stparameter[j][paramtype].upper > 360){
					if(stparameter[j][paramtype].lower < 0){
						stparameter[j][paramtype].lower=0;
						stparameter[j][paramtype].upper=360;
					}else{
						stparameter[j][paramtype].upper = stparameter[j][paramtype].upper-360;
						stparameter[j][paramtype].rangetype=1;
					}
				}
				else if(stparameter[j][paramtype].lower < 0){
					stparameter[j][paramtype].lower = stparameter[j][paramtype].lower+360;
					stparameter[j][paramtype].rangetype=1;
				}
			}

			//cout<<to_string(stparameter[j][paramtype].minval) + "-" + to_string(stparameter[j][paramtype].maxval) + "-" + to_string(stparameter[j][paramtype].span) + "-" + to_string(stparameter[j][paramtype].binnum) + "-" + to_string(stparameter[j][paramtype].mean) + "-" + to_string(stparameter[j][paramtype].std)<<endl;

			stparameter[j][paramtype].energy.resize(binnum);
			stparameter[j][paramtype].binary.resize(binnum);
			stparameter[j][paramtype].distr.resize(binnum);
			stparameter[j][paramtype].boundary.resize(binnum-1);
			stparameter[j][paramtype].center.resize(binnum);
			paramnum+=1;
		}

		vector<vector<double>> vec2denergies;
		loadparamvec2dd(basepath+"/all/"+paramname[j]+"/bin_energies.bin", paramnum, binnum, vec2denergies);
		vector<vector<double>> vec2dbinary;
		loadparamvec2dd(basepath+"/all/"+paramname[j]+"/bin_binary.bin", paramnum, binnum, vec2dbinary);
		vector<vector<int>> vec2ddistribution;
		loadparamvec2di(basepath+"/all/"+paramname[j]+"/bin_distribution.bin", paramnum, binnum, vec2ddistribution);
	
		for(int k=0;k<paramnum;k++){
			stparameter[j][k].energy=vec2denergies[k];

			stparameter[j][k].total=0;
			for(int l=0; l<vec2dbinary[k].size();l++){
				stparameter[j][k].binary[l]=vec2dbinary[k][l];
				stparameter[j][k].total+=int(1-stparameter[j][k].binary[l]);
			}
			//cout<<stparameter[j][k].total<<endl;

			stparameter[j][k].totaldist=0;
			for(int l=0; l<vec2ddistribution[k].size();l++){
				stparameter[j][k].distr[l]=vec2ddistribution[k][l];
				stparameter[j][k].totaldist+=vec2ddistribution[k][l];
			}
			//cout<<stparameter[j][k].totaldist<<endl;
	
			for(int l=0; l<binnum-1; l++) stparameter[j][k].boundary[l]=stparameter[j][k].minval+l*stparameter[j][k].span;
			for(int l=0; l<binnum; l++) stparameter[j][k].center[l]=stparameter[j][k].minval+(l-0.5)*stparameter[j][k].span;
		}
	}
	return true;
}

bool LoadParams::loadstandardparameters(string basepath, map< char, vector<Binnedparamtype>> &stparameter, string paramname, string paramsetname)
{
	string	oneline;
	int ssidx, paramtype, paramnum, binnum;
	double minval, maxval, span, mean, stdval, mode;
	int count;
	string field;
	vector<string> nts = {"A","C","G","U"};
	for(int i=0;i<4;i++)
	{
		paramnum=0;
		ifstream inputFile(basepath+"/"+nts[i][0]+"/bin_info_" + paramname + ".txt");
		
		if(!inputFile.is_open())
		{
			//cout<<basepath+"/"+nts[i][0]+"/bin_info_" + paramname + ".txt"<<endl;
			cout<<"Bin info file can not be opened (" + paramname + ")\n"<<endl;
			return false;
		}

		while(getline(inputFile, oneline))
		{
			stringstream ssbase(oneline);
			count=0;
			while (getline(ssbase, field, ','))
			{
				//if(count==0) ssidx=stoi(field); //ssidx
				if(count==0) paramtype=stoi(field); //paramtype
				if(count==1) minval=stof(field); //minval
				if(count==2) maxval=stof(field); //maxval
				if(count==3) span=stof(field); //span
				if(count==4) binnum=stoi(field); //binnum
				if(count==5) mean=stof(field); //mean
				if(count==6) stdval=stof(field); //stdval
				if(count==7) mode=stof(field); //stdval
				count+=1;
			}

			if(stparameter[nts[i][0]].size() < paramtype+1) stparameter[nts[i][0]].resize(paramtype+1);
			stparameter[nts[i][0]][paramtype].minval=minval;
			stparameter[nts[i][0]][paramtype].maxval=maxval;
			stparameter[nts[i][0]][paramtype].span=span;
			stparameter[nts[i][0]][paramtype].binnum=binnum;
			stparameter[nts[i][0]][paramtype].mean=mean;
			stparameter[nts[i][0]][paramtype].std=stdval;
			stparameter[nts[i][0]][paramtype].mode=mode;

			stparameter[nts[i][0]][paramtype].upper=mean+outlinethread[paramsetname][paramname][paramtype]*stdval;
			stparameter[nts[i][0]][paramtype].lower=mean-outlinethread[paramsetname][paramname][paramtype]*stdval;
			
			stparameter[nts[i][0]][paramtype].supper=mean+outlinethread["s"+paramsetname][paramname][paramtype]*stdval;
			stparameter[nts[i][0]][paramtype].slower=mean-outlinethread["s"+paramsetname][paramname][paramtype]*stdval;
			
			stparameter[nts[i][0]][paramtype].rangetype=0;
			stparameter[nts[i][0]][paramtype].srangetype=0;
			if(paramname == "torsion"){
				if(stparameter[nts[i][0]][paramtype].upper > 360){
					if(stparameter[nts[i][0]][paramtype].lower < 0){
						stparameter[nts[i][0]][paramtype].lower=0;
						stparameter[nts[i][0]][paramtype].upper=360;
					}else{
						stparameter[nts[i][0]][paramtype].upper = stparameter[nts[i][0]][paramtype].upper-360;
						stparameter[nts[i][0]][paramtype].rangetype=1;
					}
				}
				else if(stparameter[nts[i][0]][paramtype].lower < 0){
					stparameter[nts[i][0]][paramtype].lower = stparameter[nts[i][0]][paramtype].lower+360;
					stparameter[nts[i][0]][paramtype].rangetype=1;
				}

				if(stparameter[nts[i][0]][paramtype].supper > 360){
					if(stparameter[nts[i][0]][paramtype].slower < 0){
						stparameter[nts[i][0]][paramtype].slower=0;
						stparameter[nts[i][0]][paramtype].supper=360;
					}else{
						stparameter[nts[i][0]][paramtype].supper = stparameter[nts[i][0]][paramtype].supper-360;
						stparameter[nts[i][0]][paramtype].srangetype=1;
					}
				}
				else if(stparameter[nts[i][0]][paramtype].slower < 0){
					stparameter[nts[i][0]][paramtype].slower = stparameter[nts[i][0]][paramtype].slower+360;
					stparameter[nts[i][0]][paramtype].srangetype=1;
				}
			}

			//cout<<to_string(stparameter[nts[i][0]][paramtype].minval) + "-" + to_string(stparameter[nts[i][0]][paramtype].maxval) + "-" + to_string(stparameter[nts[i][0]][paramtype].span) + "-" + to_string(stparameter[nts[i][0]][paramtype].binnum) + "-" + to_string(stparameter[nts[i][0]][paramtype].mean) + "-" + to_string(stparameter[nts[i][0]][paramtype].std) + "-" + to_string(stparameter[nts[i][0]][paramtype].mode) + "-" + to_string(stparameter[nts[i][0]][paramtype].upper) + "-" + to_string(stparameter[nts[i][0]][paramtype].lower) + "-" + to_string(stparameter[nts[i][0]][paramtype].rangetype)<<endl;

			//distrvec.resize(binnum);
			//elementwisedih[elenum][nt][paramnum].distr.resize(binnum);
			stparameter[nts[i][0]][paramtype].energy.resize(binnum);
			stparameter[nts[i][0]][paramtype].binary.resize(binnum);
			stparameter[nts[i][0]][paramtype].logdistr.resize(binnum);
			stparameter[nts[i][0]][paramtype].distr.resize(binnum);
			stparameter[nts[i][0]][paramtype].boundary.resize(binnum-1);
			stparameter[nts[i][0]][paramtype].center.resize(binnum);
			paramnum+=1;
		}

		vector<vector<double>> vec2denergies;
		//loadparamtensor3dd(basepath+"/"+nts[i]+"/bin_energies_" + paramname + ".bin", stparameter.size(), paramnum, binnum, vec3denergies);
		loadparamvec2dd(basepath+"/"+nts[i][0]+"/bin_energies_" + paramname + ".bin", paramnum, binnum, vec2denergies);
		vector<vector<double>> vec2dbinary;
		loadparamvec2dd(basepath+"/"+nts[i][0]+"/bin_binary_" + paramname + ".bin", paramnum, binnum, vec2dbinary);
		vector<vector<int>> vec2ddistribution;
		loadparamvec2di(basepath+"/"+nts[i][0]+"/bin_distribution_" + paramname + ".bin", paramnum, binnum, vec2ddistribution);
		for(int j=0;j<paramnum;j++){
			stparameter[nts[i][0]][j].energy=vec2denergies[j];
			//if(j==0){
				//for(int k=0;k<vec2denergies[j].size();k++){
					//cout<<to_string(k)+"-"+to_string(stparameter[nts[i][0]][j].energy[k])+"-"+to_string(vec2ddistribution[j][k])+"|";
				//}
				//cout<<""<<endl;
			//}
			//exit(1);
			stparameter[nts[i][0]][j].total=0;
			for(int k=0; k<vec2dbinary[j].size();k++){
				stparameter[nts[i][0]][j].binary[k]=vec2dbinary[j][k];
				stparameter[nts[i][0]][j].total+=int(1-stparameter[nts[i][0]][j].binary[k]);
			}

			//cout<<stparameter[nts[i][0]][j].total<<endl;

			stparameter[nts[i][0]][j].totaldist=0;
			stparameter[nts[i][0]][j].totallogdist=0;
			for(int k=0; k<vec2ddistribution[j].size();k++){
				stparameter[nts[i][0]][j].distr[k]=vec2ddistribution[j][k];

				if(vec2ddistribution[j][k] == 0) stparameter[nts[i][0]][j].logdistr[k]=0;
				else stparameter[nts[i][0]][j].logdistr[k]=log(static_cast<double>(vec2ddistribution[j][k]));
				//cout<<to_string(k) + "-" + to_string(stparameter[nts[i][0]][j].distr[k]) + "-" + to_string(stparameter[nts[i][0]][j].logdistr[k])<<endl;
				stparameter[nts[i][0]][j].totaldist+=vec2ddistribution[j][k];
				stparameter[nts[i][0]][j].totallogdist+=stparameter[nts[i][0]][j].logdistr[k];
			}
			//cout<<stparameter[nts[i][0]][j].totaldist<<endl;

			for(int k=0; k<binnum-1; k++) stparameter[nts[i][0]][j].boundary[k]=stparameter[nts[i][0]][j].minval+k*stparameter[nts[i][0]][j].span;
			for(int k=0; k<binnum; k++) stparameter[nts[i][0]][j].center[k]=stparameter[nts[i][0]][j].minval+(k-0.5)*stparameter[nts[i][0]][j].span;
		}
	}
	return true;
}

bool LoadParams::loadcombtorsionparameter(string basepath, map<char, vector<D2Binnedparamtype>> &stparameter){
	vector<string> nts = {"A","C","G","U"};
	
	for(int i=0;i<4;i++)
	{
		load2ddistribution(basepath+"/" + nts[i] + "/bin_info.txt", basepath+"/" + nts[i] + "/bin_energies.bin", basepath+"/" + nts[i] + "/bin_binary.bin", basepath+"/" + nts[i] + "/bin_distribution.bin", stparameter[nts[i][0]]);
	}

	return true;
}

bool LoadParams::loadclashexpatoms(string basepath, map<char, vector<vector<bool>>> &stparameter1, map<char, map<char, vector<vector<bool>>>> &stparameter2, map<char, map<char, vector<vector<bool>>>> &stparameter3){
	int i, j, k, l;
	vector<char> nts={'A','C','G','U'};
	int atomnum, atomnum1, atomnum2;

	for(i=0;i<nts.size();i++){
		atomnum=NNAll[Reschar2Idint[nts[i]]][0];

		string str(1, nts[i]);

		vector<vector<int>> vec2ddistribution;
		loadparamvec2di(basepath+"/interresidue/"+str+"/clash_exp_mat_"+str+".bin", atomnum, atomnum, vec2ddistribution);

		stparameter1[nts[i]].resize(atomnum, vector<bool>(atomnum));
		for(j=0;j<vec2ddistribution.size();j++){
			for(k=0;k<vec2ddistribution[j].size();k++){
				//cout<<to_string(j) + "-" + to_string(k) + "-" + to_string(vec2ddistribution[j][k])<<endl;
				stparameter1[nts[i]][j][k]=static_cast<bool>(vec2ddistribution[j][k]);
			}
		}
	}

	for(i=0;i<nts.size();i++){
		for(j=0;j<nts.size();j++){
			atomnum1=NNAll[Reschar2Idint[nts[i]]][0];
			atomnum2=NNAll[Reschar2Idint[nts[j]]][0];

			string str1(1, nts[i]);
			string str2(1, nts[j]);

			vector<vector<int>> vec2ddistribution;
			loadparamvec2di(basepath+"/adjacent/"+str1+"_"+str2+"/clash_exp_mat_"+str1+"_"+str2+".bin", atomnum1, atomnum2, vec2ddistribution);

			stparameter2[nts[i]][nts[j]].resize(atomnum1, vector<bool>(atomnum2));
			for(k=0;k<vec2ddistribution.size();k++){
				for(l=0;l<vec2ddistribution[k].size();l++){
					//cout<<to_string(k) + "-" + to_string(l) + "-" + to_string(vec2ddistribution[k][l])<<endl;
					stparameter2[nts[i]][nts[j]][k][l]=static_cast<bool>(vec2ddistribution[k][l]);
				}
			}
		}
	}

	for(i=0;i<nts.size();i++){
		for(j=0;j<nts.size();j++){
			atomnum1=NNAll[Reschar2Idint[nts[i]]][0];
			atomnum2=NNAll[Reschar2Idint[nts[j]]][0];

			string str1(1, nts[i]);
			string str2(1, nts[j]);

			vector<vector<int>> vec2ddistribution;
			loadparamvec2di(basepath+"/any/"+str1+"_"+str2+"/clash_exp_mat_"+str1+"_"+str2+".bin", atomnum1, atomnum2, vec2ddistribution);

			stparameter3[nts[i]][nts[j]].resize(atomnum1, vector<bool>(atomnum2));
			for(k=0;k<vec2ddistribution.size();k++){
				for(l=0;l<vec2ddistribution[k].size();l++){
					//cout<<to_string(k) + "-" + to_string(l) + "-" + to_string(vec2ddistribution[k][l])<<endl;
					stparameter3[nts[i]][nts[j]][k][l]=static_cast<bool>(vec2ddistribution[k][l]);
				}
			}
		}
	}
	//exit(1);
	return true;
}

bool LoadParams::loadbasetemplate(string basepath, map<char, vector<vector<double>>> &stparameter){
	int i, j;
	string oneline;
	string field;
	vector<string> nts = {"A","C","G","U"};
	vector<double> templateinfo;
	int count;

	for(i=0;i<nts.size();i++)
	{
		ifstream inputFile(basepath+"/"+nts[i][0]+"/base_templates_"+nts[i][0]+".txt");
		
		if(!inputFile.is_open())
		{
			cout<<"Info file can not be opened (base template)\n"<<endl;
			return false;
		}

		while(getline(inputFile, oneline))
		{
			vector<double> temp;
			count=0;
			stringstream ssbase(oneline);
			while (getline(ssbase, field, ','))
			{
				if(count == 0){
					temp.push_back(stof(field));
					count+=1;
				}
				else if(count == 1){
					temp.push_back(stof(field));
					count+=1;
				}
				else if(count == 2){
					temp.push_back(stof(field));
					stparameter[nts[i][0]].push_back(temp);
					temp.clear();
					count=0;
				}
			}
		}
	}
	return true;
}

bool LoadParams::loadbaepairhydrogenparameters(string basepath, map<char, map<char, map<int, map<string, vector<vector<Binnedparamtype>>>>>> &stparameter)
{
	int i,j,k,l,n;
	string	oneline;
	int ssidx, paramtype, paramnum, binnum;
	double minval, maxval, span, mean, stdval, mode;
	int count;
	string field;
	vector<string> nts = {"A","U","G","C"};
	map<int, string> paramname = {{0, "length"}, {1, "angle"}};

	for(int i=0;i<4;i++)
	{
		for(int j=i;j<4;j++)
		{
			for(const auto& pair : HbondParamInfo[nts[i][0]][nts[j][0]])
			{
				for(const auto& hb : HbondParamInfo[nts[i][0]][nts[j][0]][pair.first])
				{
					stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first].resize(2);
					for(n=0; n<2; n++){
						paramnum=0;
						ifstream inputFile(basepath+"/"+nts[i][0]+"_"+nts[j][0]+"/"+to_string(pair.first)+"/bin_info_"+ paramname[n] + "_" + hb.first + ".txt");
						
						if(!inputFile.is_open())
						{
							cout<<"Bin info file can not be opened (" + paramname[n] + ")\n"<<endl;
							return false;
						}

						while(getline(inputFile, oneline))
						{
							stringstream ssbase(oneline);
							count=0;
							while (getline(ssbase, field, ','))
							{
								//if(count==0) ssidx=stoi(field); //ssidx
								if(count==0) paramtype=stoi(field); //paramtype
								if(count==1) minval=stof(field); //minval
								if(count==2) maxval=stof(field); //maxval
								if(count==3) span=stof(field); //span
								if(count==4) binnum=stoi(field); //binnum
								if(count==5) mean=stof(field); //mean
								if(count==6) stdval=stof(field); //stdval
								if(count==7) mode=stof(field); //stdval
								count+=1;
							}

							if(stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n].size() < paramtype+1) stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n].resize(paramtype+1);
							stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][paramtype].minval=minval;
							stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][paramtype].maxval=maxval;
							stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][paramtype].span=span;
							stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][paramtype].binnum=binnum;
							stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][paramtype].mean=mean;
							stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][paramtype].std=stdval;
							stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][paramtype].mode=mode;
							stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][paramtype].upper=mean+outlinerbphydstd[paramname[n]]*stdval;
							stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][paramtype].lower=mean-outlinerbphydstd[paramname[n]]*stdval;
							
							stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][paramtype].energy.resize(binnum);
							stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][paramtype].binary.resize(binnum);
							stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][paramtype].logdistr.resize(binnum);
							stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][paramtype].distr.resize(binnum);
							stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][paramtype].boundary.resize(binnum-1);
							stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][paramtype].center.resize(binnum);
							paramnum+=1;
						}

						vector<vector<double>> vec2denergies;
						//loadparamtensor3dd(basepath+"/"+nts[i]+"/bin_energies_" + paramname + ".bin", stparameter.size(), paramnum, binnum, vec3denergies);
						loadparamvec2dd(basepath+"/"+nts[i][0]+"_"+nts[j][0]+"/"+to_string(pair.first)+"/bin_energies_"+ paramname[n] + "_" + hb.first + ".bin", paramnum, binnum, vec2denergies);
						vector<vector<double>> vec2dbinary;
						loadparamvec2dd(basepath+"/"+nts[i][0]+"_"+nts[j][0]+"/"+to_string(pair.first)+"/bin_binary_"+ paramname[n] + "_" + hb.first + ".bin", paramnum, binnum, vec2dbinary);
						vector<vector<int>> vec2ddistribution;
						loadparamvec2di(basepath+"/"+nts[i][0]+"_"+nts[j][0]+"/"+to_string(pair.first)+"/bin_distribution_"+ paramname[n] + "_" + hb.first + ".bin", paramnum, binnum, vec2ddistribution);
						
						for(int l=0;l<paramnum;l++){
							stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][l].energy=vec2denergies[l];
							
							stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][l].total=0;
							for(int k=0; k<vec2dbinary[l].size();k++){
								stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][l].binary[k]=vec2dbinary[l][k];
								stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][l].total+=int(1-stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][l].binary[k]);
							}

							stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][l].totaldist=0;
							stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][l].totallogdist=0;
							for(int k=0; k<vec2ddistribution[l].size();k++){
								stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][l].distr[k]=vec2ddistribution[l][k];

								if(vec2ddistribution[l][k] == 0) stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][l].logdistr[k]=0;
								else stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][l].logdistr[k]=log(static_cast<double>(vec2ddistribution[l][k]));
								stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][l].totaldist+=vec2ddistribution[l][k];
								stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][l].totallogdist+=stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][l].logdistr[k];
							}
							//cout<<stparameter[nts[i][0]][j].totaldist<<endl;

							for(int k=0; k<binnum-1; k++) stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][l].boundary[k]=stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][l].minval+k*stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][l].span;
							for(int k=0; k<binnum; k++) stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][l].center[k]=stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][l].minval+(k-0.5)*stparameter[nts[i][0]][nts[j][0]][pair.first][hb.first][n][l].span;
						}
					}
				}
				if(nts[i] != nts[j]){
					stparameter[nts[j][0]][nts[i][0]]=stparameter[nts[i][0]][nts[j][0]];
				}
			}
		}
	}

	return true;
}

bool LoadParams::loadtorsionlist(string basepath, vector<vector<vector<Torsionclsparamtype>>> &stparameter){

	int i, j, k;
	double val;
	string	oneline, field;
	ifstream inputFile;
	vector<string> paramlist = {"backbone", "sugar"};
	stparameter.resize(2);

	for(i=0;i<paramlist.size();i++)
	{
		inputFile.open(basepath+"/"+paramlist[i]+"/cls_torsion_mean.csv");
				
		if(!inputFile.is_open())
		{
			cout<<"Bin info file can not be opened (torsion list (mean))\n"<<endl;
			return false;
		}

		vector<vector<double>> vallist;
		while(getline(inputFile, oneline))
		{
			stringstream ssbase(oneline);
			vector<double> vals;
			while (getline(ssbase, field, ','))
			{
				val=stof(field);
				vals.push_back(val);
			}
			vallist.push_back(vals);
		}

		stparameter[i].resize(vallist.size());
		for(j=0;j<vallist.size();j++){
			stparameter[i][j].resize(vallist[j].size());
			for(k=0;k<vallist[j].size();k++){
				stparameter[i][j][k].mean = vallist[j][k];
			}
		}

		inputFile.close();

		inputFile.open(basepath+"/"+paramlist[i]+"/cls_torsion_std.csv");
		if(!inputFile.is_open())
		{
			cout<<"Bin info file can not be opened (torsion list (std)\n"<<endl;
			return false;
		}

		vallist.clear();
		//vals.clear();

		while(getline(inputFile, oneline))
		{
			vector<double> vals;
			stringstream ssbase(oneline);
			while (getline(ssbase, field, ','))
			{
				val=stof(field);
				vals.push_back(val);
			}
			vallist.push_back(vals);
		}

		inputFile.close();

		stparameter[i].resize(vallist.size());
		for(j=0;j<vallist.size();j++){
			stparameter[i][j].resize(vallist[j].size());
			for(k=0;k<vallist[j].size();k++){
				stparameter[i][j][k].std = vallist[j][k];
				stparameter[i][j][k].upper=stparameter[i][j][k].mean+outlinertorcls[i]*stparameter[i][j][k].std;
				stparameter[i][j][k].lower=stparameter[i][j][k].mean-outlinertorcls[i]*stparameter[i][j][k].std;
			}
		}

		for(j=0;j<stparameter[i].size();j++){
			for(k=0;k<stparameter[i][j].size();k++){
				stparameter[i][j][k].rangetype=0;
			
				if(stparameter[i][j][k].upper > 360){
					if(stparameter[i][j][k].lower < 0){
						stparameter[i][j][k].lower=0;
						stparameter[i][j][k].upper=360;
					}else{
						stparameter[i][j][k].upper = stparameter[i][j][k].upper-360;
						stparameter[i][j][k].rangetype=1;
					}
				}
				else if(stparameter[i][j][k].lower < 0){
					stparameter[i][j][k].lower = stparameter[i][j][k].lower+360;
					stparameter[i][j][k].rangetype=1;
				}
				//cout<<to_string(i) + "-" + to_string(j) + "-" + to_string(k) + "-" + to_string(stparameter[i][j][k].mean) + "-" + to_string(stparameter[i][j][k].std) + "-" + to_string(stparameter[i][j][k].upper) + "-" + to_string(stparameter[i][j][k].lower) + "-" + to_string(stparameter[i][j][k].rangetype)<<endl;
			}
		}
	}

	return true;
}

bool LoadParams::loadpuckerparameters(string basepath, vector<vector<Binnedparamtype>> &stparameter){
	int i, j, k;
	string	oneline;
	int ssidx, paramtype, classnum, paramnum, binnum;
	double minval, maxval, span, mean, stdval, mode;
	int count;
	string field;

	stparameter.resize(2);

	map<int, string> paramname = {{0, "length"}, {1, "torsion"}};
	for(i=0;i<paramname.size();i++){
		paramnum=0;
		ifstream inputFile(basepath + "/"  + paramname[i] + "/bin_info_" + paramname[i] + ".txt");
		
		if(!inputFile.is_open())
		{
			cout<<"Bin info file can not be opened (pucker; " + paramname[i] + ")\n"<<endl;
			return false;
		}

		while(getline(inputFile, oneline))
		{
			stringstream ssbase(oneline);
			count=0;
			while (getline(ssbase, field, ','))
			{
				//if(count==0) ssidx=stoi(field); //ssidx
				if(count==0) classnum=stoi(field); //paramtype
				if(count==1) minval=stof(field); //minval
				if(count==2) maxval=stof(field); //maxval
				if(count==3) span=stof(field); //span
				if(count==4) binnum=stoi(field); //binnum
				if(count==5) mean=stof(field); //mean
				if(count==6) stdval=stof(field); //stdval
				if(count==7) mode=stof(field); //stdval
				count+=1;
			}

			if(stparameter.size() < classnum+1) stparameter.resize(classnum+1);
			if(stparameter[classnum].size() < i+1) stparameter[classnum].resize(i+1);

			stparameter[classnum][i].minval=minval;
			stparameter[classnum][i].maxval=maxval;
			stparameter[classnum][i].span=span;
			stparameter[classnum][i].binnum=binnum;
			stparameter[classnum][i].mean=mean;
			stparameter[classnum][i].std=stdval;
			stparameter[classnum][i].mode=mode;
			stparameter[classnum][i].upper=mean+outlinerpucker[paramname[i]]*stdval;
			stparameter[classnum][i].lower=mean-outlinerpucker[paramname[i]]*stdval;

			//cout<<to_string(stparameter[classnum][i].minval) + "-" + to_string(stparameter[classnum][i].maxval) + "-" + to_string(stparameter[classnum][i].span) + "-" + to_string(stparameter[classnum][i].binnum) + "-" + to_string(stparameter[classnum][i].mean) + "-" + to_string(stparameter[classnum][i].std) + "-" + to_string(stparameter[classnum][i].mode) + "-" + to_string(stparameter[classnum][i].upper) + "-" + to_string(stparameter[classnum][i].lower) + "-" + to_string(stparameter[classnum][i].rangetype)<<endl;

			//distrvec.resize(binnum);
			//elementwisedih[elenum][nt][paramnum].distr.resize(binnum);
			stparameter[classnum][i].energy.resize(binnum);
			stparameter[classnum][i].binary.resize(binnum);
			stparameter[classnum][i].logdistr.resize(binnum);
			stparameter[classnum][i].distr.resize(binnum);
			stparameter[classnum][i].boundary.resize(binnum-1);
			stparameter[classnum][i].center.resize(binnum);

			paramnum+=1;
		}

		vector<vector<double>> vec2denergies;
		//loadparamtensor3dd(basepath+"/"+nts[i]+"/bin_energies_" + paramname + ".bin", stparameter.size(), paramnum, binnum, vec3denergies);
		loadparamvec2dd(basepath + "/"  + paramname[i] + "/bin_energies_" + paramname[i] + ".bin", paramnum, binnum, vec2denergies);
		vector<vector<double>> vec2dbinary;
		loadparamvec2dd(basepath + "/"  + paramname[i] + "/bin_binary_" + paramname[i] + ".bin", paramnum, binnum, vec2dbinary);
		vector<vector<int>> vec2ddistribution;
		loadparamvec2di(basepath + "/"  + paramname[i] + "/bin_distribution_" + paramname[i] + ".bin", paramnum, binnum, vec2ddistribution);
		for(int j=0;j<paramnum;j++){
			stparameter[j][i].energy=vec2denergies[j];
			//if(j==0){
				//for(int k=0;k<vec2denergies[j].size();k++){
					//cout<<to_string(k)+"-"+to_string(stparameter[j][i].energy[k])+"-"+to_string(vec2ddistribution[j][k])+"|";
				//}
				//cout<<""<<endl;
			//}
			//exit(1);
			stparameter[j][i].total=0;
			for(int k=0; k<vec2dbinary[j].size();k++){
				stparameter[j][i].binary[k]=vec2dbinary[j][k];
				stparameter[j][i].total+=int(1-stparameter[j][i].binary[k]);
			}

			//cout<<stparameter[j][i].total<<endl;

			stparameter[j][i].totaldist=0;
			stparameter[j][i].totallogdist=0;
			for(int k=0; k<vec2ddistribution[j].size();k++){
				stparameter[j][i].distr[k]=vec2ddistribution[j][k];

				if(vec2ddistribution[j][k] == 0) stparameter[j][i].logdistr[k]=0;
				else stparameter[j][i].logdistr[k]=log(static_cast<double>(vec2ddistribution[j][k]));
				//cout<<to_string(k) + "-" + to_string(stparameter[j][i].distr[k]) + "-" + to_string(stparameter[j][i].logdistr[k])<<endl;
				stparameter[j][i].totaldist+=vec2ddistribution[j][k];
				stparameter[j][i].totallogdist+=stparameter[j][i].logdistr[k];
			}
			//cout<<stparameter[j][i].totaldist<<endl;

			for(int k=0; k<binnum-1; k++) stparameter[j][i].boundary[k]=stparameter[j][i].minval+k*stparameter[j][i].span;
			for(int k=0; k<binnum; k++) stparameter[j][i].center[k]=stparameter[j][i].minval+(k-0.5)*stparameter[j][i].span;
		}
	}

	return true;

}


bool LoadParams::loadpairlengnt(string basepath, map<char, vector<vector<vector<PairLengInfo>>>> &stparameter){
	int i;
	vector<string> paramtype = {"type1", "type2", "type3", "type4"};
	vector<string> nts = {"A","U","G","C"};
	string oneline;
	int count, paramidx;
	int atom1, atom2;
	string field;
	double mean, std;
	int pidx;
	int atomidx1, atomidx2;
	double vwradi;
	double lowertmp;
	bool flagatom;

	for(i=0; i<nts.size();i++){
		pidx=0;
		stparameter[nts[i][0]].resize(paramtype.size(), vector<vector<PairLengInfo>>(PAIRLENGTHTARGETNUM, vector<PairLengInfo>(PAIRLENGTHTARGETNUM)));

		for(auto& param: paramtype){
			ifstream inputFile(basepath+"/"+param+"/nt_bbsugar_pair_length_"+nts[i]+".csv");

			if(!inputFile.is_open())
			{
				cout<<"CSV info file can not be opened (grouped adjacent pair atom lengths)\n"<<endl;
				return false;
			}

			//stparameter[nts[i][0]].resize(PAIRLENGTHTARGETNUM, vector<PairLengInfo>(PAIRLENGTHTARGETNUM));

			paramidx=0;
			while(getline(inputFile, oneline))
			{
				stringstream ssbase(oneline);
				count=0;
				while (getline(ssbase, field, ','))
				{
					//if(count==0) ssidx=stoi(field); //ssidx
					if(count==0) mean=stof(field); //paramtype
					if(count==1) std=stof(field); //minval
					count+=1;
				}

				atom1=paramidx/PAIRLENGTHTARGETNUM;
				atom2=paramidx%PAIRLENGTHTARGETNUM;

				//cout<<nts[i]<<endl;
				//cout<<to_string(paramidx) + "-" + to_string(atom1) + "-" + to_string(atom2) + "-" + to_string(mean) + "-" + to_string(std)<<endl;

				stparameter[nts[i][0]][pidx][atom1][atom2].mean=mean;
				stparameter[nts[i][0]][pidx][atom1][atom2].std=std;
				stparameter[nts[i][0]][pidx][atom1][atom2].upper=mean+outpairlengstd["nt"]*std;
				//stparameter[nts[i][0]][pidx][atom1][atom2].lower=mean-outpairlengstd["nt"]*std;

				flagatom = EXPCLASHATOMSNT[nts[i][0]][atom1][atom2];
				lowertmp=mean-outpairlengstd["adj"]*std;
				if(!flagatom){
					atomidx1=NNheavy[Reschar2Idint[nts[i][0]]][atom1+1];
					atomidx2=NNheavy[Reschar2Idint[nts[i][0]]][atom2+1];
					vwradi=VWPARAM[atomidx1]+VWPARAM[atomidx2]-0.4;
					if(vwradi > lowertmp) stparameter[nts[i][0]][pidx][atom1][atom2].lower=vwradi;
					else stparameter[nts[i][0]][pidx][atom1][atom2].lower=lowertmp;
				}else{
					stparameter[nts[i][0]][pidx][atom1][atom2].lower=lowertmp;
				}

				paramidx+=1;
			}
			pidx+=1;
		}
	}
	return true;
}

bool LoadParams::loadpairlengadj(string basepath, map<char, map<char, vector<vector<vector<PairLengInfo>>>>> &stparameter){
	int i, j;
	vector<string> paramtype = {"type1"};
	vector<string> nts = {"A","U","G","C"};
	string oneline;
	int count, paramidx;
	int atom1, atom2;
	string field;
	double mean, std;
	int pidx;
	int atomidx1, atomidx2;
	double vwradi;
	double lowertmp;
	bool flagatom;

	for(i=0; i<nts.size();i++){
		for(j=0; j<nts.size();j++){
			
			pidx=0;
			stparameter[nts[i][0]][nts[j][0]].resize(paramtype.size(), vector<vector<PairLengInfo>>(PAIRLENGTHTARGETNUM, vector<PairLengInfo>(PAIRLENGTHTARGETNUM)));

			for(auto& param: paramtype){
				ifstream inputFile(basepath+"/"+param+"/adj_bbsugar_pair_length_"+nts[i]+nts[j]+".csv");

				if(!inputFile.is_open())
				{
					cout<<"CSV info file can not be opened (grouped pair atom lengths)\n"<<endl;
					return false;
				}

				//stparameter[nts[i][0]][nts[j][0]].resize(PAIRLENGTHTARGETNUM, vector<PairLengInfo>(PAIRLENGTHTARGETNUM));

				paramidx=0;
				while(getline(inputFile, oneline))
				{
					stringstream ssbase(oneline);
					count=0;
					while (getline(ssbase, field, ','))
					{
						//if(count==0) ssidx=stoi(field); //ssidx
						if(count==0) mean=stof(field); //paramtype
						if(count==1) std=stof(field); //minval
						count+=1;
					}

					atom1=paramidx/PAIRLENGTHTARGETNUM;
					atom2=paramidx%PAIRLENGTHTARGETNUM;

					//cout<<nts[i]+"-"+nts[j]<<endl;
					//cout<<to_string(paramidx) + "-" + to_string(atom1) + "-" + to_string(atom2) + "-" + to_string(mean) + "-" + to_string(std)<<endl;

					stparameter[nts[i][0]][nts[j][0]][pidx][atom1][atom2].mean=mean;
					stparameter[nts[i][0]][nts[j][0]][pidx][atom1][atom2].std=std;
					stparameter[nts[i][0]][nts[j][0]][pidx][atom1][atom2].upper=mean+outpairlengstd["adj"]*std;
					//stparameter[nts[i][0]][nts[j][0]][pidx][atom1][atom2].lower=mean-outpairlengstd["adj"]*std;

					flagatom = EXPCLASHATOMSADJ[nts[i][0]][nts[j][0]][atom1][atom2];
					lowertmp=mean-outpairlengstd["adj"]*std;
					if(!flagatom){
						atomidx1=NNheavy[Reschar2Idint[nts[i][0]]][atom1+1];
						atomidx2=NNheavy[Reschar2Idint[nts[j][0]]][atom2+1];
						vwradi=VWPARAM[atomidx1]+VWPARAM[atomidx2]-0.4;
						if(vwradi > lowertmp) stparameter[nts[i][0]][nts[j][0]][pidx][atom1][atom2].lower=vwradi;
						else stparameter[nts[i][0]][nts[j][0]][pidx][atom1][atom2].lower=lowertmp;
					}else{
						stparameter[nts[i][0]][nts[j][0]][pidx][atom1][atom2].lower=lowertmp;
					}

					paramidx+=1;
				}
				pidx+=1;
			}
		}
	}
	return true;
}

bool LoadParams::loadbppairlengnt(string basepath, map<char, vector<vector<vector<PairLengInfo>>>> &stparameter){
	int i;
	vector<string> paramtype = {"nonbp", "bp"};
	vector<string> nts = {"A","U","G","C"};
	string oneline;
	int count, paramidx;
	int atom1, atom2;
	string field;
	double mean, std;
	int pidx;
	int atomidx1, atomidx2;
	double vwradi;
	double lowertmp;
	bool flagatom;

	for(i=0; i<nts.size();i++){
		pidx=0;
		stparameter[nts[i][0]].resize(paramtype.size(), vector<vector<PairLengInfo>>(PAIRLENGTHTARGETNUM, vector<PairLengInfo>(PAIRLENGTHTARGETNUM)));

		for(auto& param: paramtype){
			ifstream inputFile(basepath+"/"+param+"/nt_atom_pair_length_"+nts[i]+".csv");

			if(!inputFile.is_open())
			{
				cout<<"CSV info file can not be opened (grouped pair atom lengths in base pair nucleotide)\n"<<endl;
				return false;
			}

			paramidx=0;
			while(getline(inputFile, oneline))
			{
				stringstream ssbase(oneline);
				count=0;
				while (getline(ssbase, field, ','))
				{
					if(count==0) mean=stof(field); //paramtype
					if(count==1) std=stof(field); //minval
					count+=1;
				}

				atom1=paramidx/PAIRLENGTHTARGETNUM;
				atom2=paramidx%PAIRLENGTHTARGETNUM;

				stparameter[nts[i][0]][pidx][atom1][atom2].mean=mean;
				stparameter[nts[i][0]][pidx][atom1][atom2].std=std;
				stparameter[nts[i][0]][pidx][atom1][atom2].upper=mean+outpairlengstd["adj"]*std;
				//stparameter[nts[i][0]][pidx][atom1][atom2].lower=mean-outpairlengstd["adj"]*std;

				flagatom = EXPCLASHATOMSNT[nts[i][0]][atom1][atom2];
				lowertmp=mean-outpairlengstd["adj"]*std;
				if(!flagatom){
					atomidx1=NNheavy[Reschar2Idint[nts[i][0]]][atom1+1];
					atomidx2=NNheavy[Reschar2Idint[nts[i][0]]][atom2+1];
					vwradi=VWPARAM[atomidx1]+VWPARAM[atomidx2]-0.4;
					if(vwradi > lowertmp) stparameter[nts[i][0]][pidx][atom1][atom2].lower=vwradi;
					else stparameter[nts[i][0]][pidx][atom1][atom2].lower=lowertmp;
				}else{
					stparameter[nts[i][0]][pidx][atom1][atom2].lower=lowertmp;
				}

				paramidx+=1;
			}
			pidx+=1;
		}
	}
	return true;
}

bool LoadParams::loadbppairlengbp(string basepath, map<char, map<char, vector<vector<PairLengInfo>>>> &stparameter){
	int i,j;
	vector<string> nts = {"A","U","G","C"};
	string oneline;
	int count, paramidx;
	int atom1, atom2;
	string field;
	double mean, std;
	int pidx;
	int atomidx1, atomidx2;
	double vwradi;
	double lowertmp;
	bool flagatom;

	for(i=0; i<nts.size();i++){
		for(j=0; j<nts.size();j++){
			stparameter[nts[i][0]][nts[j][0]].resize(PAIRLENGTHTARGETNUM, vector<PairLengInfo>(PAIRLENGTHTARGETNUM));

			ifstream inputFile(basepath+"/bp_atom_pair_length_"+nts[i]+nts[j]+".csv");

			if(!inputFile.is_open())
			{
				cout<<"CSV info file can not be opened (grouped pair atom lengths in base pairs)\n"<<endl;
				return false;
			}

			paramidx=0;
			while(getline(inputFile, oneline))
			{
				stringstream ssbase(oneline);
				count=0;
				while (getline(ssbase, field, ','))
				{
					if(count==0) mean=stof(field); //paramtype
					if(count==1) std=stof(field); //minval
					count+=1;
				}

				atom1=paramidx/PAIRLENGTHTARGETNUM;
				atom2=paramidx%PAIRLENGTHTARGETNUM;

				stparameter[nts[i][0]][nts[j][0]][atom1][atom2].mean=mean;
				stparameter[nts[i][0]][nts[j][0]][atom1][atom2].std=std;
				stparameter[nts[i][0]][nts[j][0]][atom1][atom2].upper=mean+outpairlengstd["adj"]*std;
				//stparameter[nts[i][0]][nts[j][0]][atom1][atom2].lower=mean-outpairlengstd["adj"]*std;

				lowertmp=mean-outpairlengstd["adj"]*std;
				atomidx1=NNheavy[Reschar2Idint[nts[i][0]]][atom1+1];
				atomidx2=NNheavy[Reschar2Idint[nts[j][0]]][atom2+1];
				vwradi=VWPARAM[atomidx1]+VWPARAM[atomidx2]-0.4;
				if(vwradi > lowertmp) stparameter[nts[i][0]][nts[j][0]][atom1][atom2].lower=vwradi;
				else stparameter[nts[i][0]][nts[j][0]][atom1][atom2].lower=lowertmp;

				paramidx+=1;
			}
		}
	}
	return true;
}

bool LoadParams::loadbppairlengadj(string basepath, map<char, map<char, vector<vector<vector<PairLengInfo>>>>> &stparameter){
	int i, j;
	vector<string> paramtype = {"stem", "closing", "opening", "unpair"};
	vector<string> nts = {"A","U","G","C"};
	string oneline;
	int count, paramidx;
	int atom1, atom2;
	string field;
	double mean, std;
	int atomidx1, atomidx2;
	double vwradi;
	double lowertmp;
	bool flagatom;
	int pidx;

	for(i=0; i<nts.size();i++){
		for(j=0; j<nts.size();j++){
			
			pidx=0;
			stparameter[nts[i][0]][nts[j][0]].resize(paramtype.size(), vector<vector<PairLengInfo>>(PAIRLENGTHTARGETNUM, vector<PairLengInfo>(PAIRLENGTHTARGETNUM)));

			for(auto& param: paramtype){
				ifstream inputFile(basepath+"/"+param+"/adj_atom_pair_length_"+nts[i]+nts[j]+".csv");

				if(!inputFile.is_open())
				{
					cout<<"CSV info file can not be opened (grouped pair atom lengths)\n"<<endl;
					return false;
				}

				//stparameter[nts[i][0]][nts[j][0]].resize(PAIRLENGTHTARGETNUM, vector<PairLengInfo>(PAIRLENGTHTARGETNUM));

				paramidx=0;
				while(getline(inputFile, oneline))
				{
					stringstream ssbase(oneline);
					count=0;
					while (getline(ssbase, field, ','))
					{
						//if(count==0) ssidx=stoi(field); //ssidx
						if(count==0) mean=stof(field); //paramtype
						if(count==1) std=stof(field); //minval
						count+=1;
					}

					atom1=paramidx/PAIRLENGTHTARGETNUM;
					atom2=paramidx%PAIRLENGTHTARGETNUM;

					//cout<<nts[i]+"-"+nts[j]<<endl;
					//cout<<to_string(paramidx) + "-" + to_string(atom1) + "-" + to_string(atom2) + "-" + to_string(mean) + "-" + to_string(std)<<endl;

					stparameter[nts[i][0]][nts[j][0]][pidx][atom1][atom2].mean=mean;
					stparameter[nts[i][0]][nts[j][0]][pidx][atom1][atom2].std=std;
					stparameter[nts[i][0]][nts[j][0]][pidx][atom1][atom2].upper=mean+outpairlengstd["adj"]*std;

					flagatom = EXPCLASHATOMSADJ[nts[i][0]][nts[j][0]][atom1][atom2];
					lowertmp=mean-outpairlengstd["adj"]*std;
					if(!flagatom){
						atomidx1=NNheavy[Reschar2Idint[nts[i][0]]][atom1+1];
						atomidx2=NNheavy[Reschar2Idint[nts[j][0]]][atom2+1];
						vwradi=VWPARAM[atomidx1]+VWPARAM[atomidx2]-0.4;
						if(vwradi > lowertmp) stparameter[nts[i][0]][nts[j][0]][pidx][atom1][atom2].lower=vwradi;
						else stparameter[nts[i][0]][nts[j][0]][pidx][atom1][atom2].lower=lowertmp;
					}else{
						stparameter[nts[i][0]][nts[j][0]][pidx][atom1][atom2].lower=lowertmp;
					}
					paramidx+=1;
				}
				pidx+=1;
			}
		}
	}
	return true;
}


bool LoadParams::loadsinglepairlengnt(string basepath, map<char, vector<vector<PairLengInfo>>> &stparameter){
	int i;
	vector<string> nts = {"A","U","G","C"};
	string oneline;
	int count, paramidx;
	int atom1, atom2;
	string field;
	double mean, std;
	int atomidx1, atomidx2;
	double vwradi;
	double lowertmp;
	bool flagatom;

	for(i=0; i<nts.size();i++){
		ifstream inputFile(basepath+"/all/nt_bbsugar_pair_length_"+nts[i]+".csv");

		if(!inputFile.is_open())
		{
			cout<<"CSV info file can not be opened (single adjacent atom lengths)\n"<<endl;
			return false;
		}

		stparameter[nts[i][0]].resize(PAIRLENGTHTARGETNUM, vector<PairLengInfo>(PAIRLENGTHTARGETNUM));

		paramidx=0;
		while(getline(inputFile, oneline))
		{
			stringstream ssbase(oneline);
			count=0;
			while (getline(ssbase, field, ','))
			{
				//if(count==0) ssidx=stoi(field); //ssidx
				if(count==0) mean=stof(field); //paramtype
				if(count==1) std=stof(field); //minval
				count+=1;
			}

			atom1=paramidx/PAIRLENGTHTARGETNUM;
			atom2=paramidx%PAIRLENGTHTARGETNUM;

			//cout<<nts[i]<<endl;
			//cout<<to_string(paramidx) + "-" + to_string(atom1) + "-" + to_string(atom2) + "-" + to_string(mean) + "-" + to_string(std)<<endl;

			stparameter[nts[i][0]][atom1][atom2].mean=mean;
			stparameter[nts[i][0]][atom1][atom2].std=std;
			stparameter[nts[i][0]][atom1][atom2].upper=mean+outpairlengstd["nt"]*std;
			//stparameter[nts[i][0]][atom1][atom2].lower=mean-outpairlengstd["nt"]*std;

			flagatom = EXPCLASHATOMSNT[nts[i][0]][atom1][atom2];
			lowertmp=mean-outpairlengstd["adj"]*std;
			if(!flagatom){
				atomidx1=NNheavy[Reschar2Idint[nts[i][0]]][atom1+1];
				atomidx2=NNheavy[Reschar2Idint[nts[i][0]]][atom2+1];
				vwradi=VWPARAM[atomidx1]+VWPARAM[atomidx2]-0.4;
				if(vwradi > lowertmp) stparameter[nts[i][0]][atom1][atom2].lower=vwradi;
				else stparameter[nts[i][0]][atom1][atom2].lower=lowertmp;
			}else{
				stparameter[nts[i][0]][atom1][atom2].lower=lowertmp;
			}

			paramidx+=1;
		}
	}
	return true;
}

bool LoadParams::loadsinglepairlengadj(string basepath, map<char, map<char, vector<vector<PairLengInfo>>>> &stparameter){
	int i, j;
	vector<string> nts = {"A","U","G","C"};
	string oneline;
	int count, paramidx;
	int atom1, atom2;
	string field;
	double mean, std;
	int atomidx1, atomidx2;
	double vwradi;
	double lowertmp;
	bool flagatom;

	for(i=0; i<nts.size();i++){
		for(j=0; j<nts.size();j++){
			ifstream inputFile(basepath+"/all/adj_bbsugar_pair_length_"+nts[i]+nts[j]+".csv");

			if(!inputFile.is_open())
			{
				cout<<"CSV info file can not be opened (single pair atom lengths)\n"<<endl;
				return false;
			}

			stparameter[nts[i][0]][nts[j][0]].resize(PAIRLENGTHTARGETNUM, vector<PairLengInfo>(PAIRLENGTHTARGETNUM));

			paramidx=0;
			while(getline(inputFile, oneline))
			{
				stringstream ssbase(oneline);
				count=0;
				while (getline(ssbase, field, ','))
				{
					//if(count==0) ssidx=stoi(field); //ssidx
					if(count==0) mean=stof(field); //paramtype
					if(count==1) std=stof(field); //minval
					count+=1;
				}

				atom1=paramidx/PAIRLENGTHTARGETNUM;
				atom2=paramidx%PAIRLENGTHTARGETNUM;

				//cout<<nts[i]+"-"+nts[j]<<endl;
				//cout<<to_string(paramidx) + "-" + to_string(atom1) + "-" + to_string(atom2) + "-" + to_string(mean) + "-" + to_string(std)<<endl;

				stparameter[nts[i][0]][nts[j][0]][atom1][atom2].mean=mean;
				stparameter[nts[i][0]][nts[j][0]][atom1][atom2].std=std;
				stparameter[nts[i][0]][nts[j][0]][atom1][atom2].upper=mean+outpairlengstd["adj"]*std;
				stparameter[nts[i][0]][nts[j][0]][atom1][atom2].lower=mean-outpairlengstd["adj"]*std;

				flagatom = EXPCLASHATOMSADJ[nts[i][0]][nts[j][0]][atom1][atom2];
				lowertmp=mean-outpairlengstd["adj"]*std;
				if(!flagatom){
					atomidx1=NNheavy[Reschar2Idint[nts[i][0]]][atom1+1];
					atomidx2=NNheavy[Reschar2Idint[nts[j][0]]][atom2+1];
					vwradi=VWPARAM[atomidx1]+VWPARAM[atomidx2]-0.4;
					if(vwradi > lowertmp) stparameter[nts[i][0]][nts[j][0]][atom1][atom2].lower=vwradi;
					else stparameter[nts[i][0]][nts[j][0]][atom1][atom2].lower=lowertmp;
				}else{
					stparameter[nts[i][0]][nts[j][0]][atom1][atom2].lower=lowertmp;
				}

				paramidx+=1;
			}
		}
	}
	return true;
}

bool LoadParams::loadtorclsparam(string basepath, map<char, vector<vector<TorClsInfo>>> &stparameter){
	int i, j;
	string oneline;
	int count, paramidx, paramgroupidx;
	string field;
	double mean, std;
	vector<string> nts = {"A","U","G","C"};
	vector<string> paramtype = {"group_1", "group_2"};

	for(i=0; i<nts.size();i++){
		stparameter[nts[i][0]].resize(paramtype.size(), vector<TorClsInfo>(CLSTORPARAMNUM));

		paramgroupidx=0;
		for(auto& param: paramtype){
			ifstream inputFile(basepath+"/"+nts[i]+"/"+param+"/parameters.txt");

			if(!inputFile.is_open())
			{
				cout<<"CSV info file can not be opened (clustered torsion)\n"<<endl;
				return false;
			}

			//stparameter[nts[i][0]][nts[j][0]].resize(PAIRLENGTHTARGETNUM, vector<PairLengInfo>(PAIRLENGTHTARGETNUM));

			paramidx=0;
			while(getline(inputFile, oneline))
			{
				stringstream ssbase(oneline);
				count=0;
				while (getline(ssbase, field, ','))
				{
					if(count==0) mean=stof(field); //paramtype
					if(count==1) std=stof(field); //minval
					count+=1;
				}

				stparameter[nts[i][0]][paramgroupidx][paramidx].mean=mean;
				stparameter[nts[i][0]][paramgroupidx][paramidx].std=std;
				stparameter[nts[i][0]][paramgroupidx][paramidx].upper=mean+outlinethread["cls_tor"]["torsion"][paramidx]*std;
				stparameter[nts[i][0]][paramgroupidx][paramidx].lower=mean-outlinethread["cls_tor"]["torsion"][paramidx]*std;

				paramidx+=1;
			}
			paramgroupidx+=1;
		}
	}
	return true;
}

bool LoadParams::loadtorcls2dparam(string basepath, map<char, vector<vector<D2Binnedparamtype>>> &stparameter){
	vector<string> nts = {"A","C","G","U"};
	
	for(int i=0;i<4;i++)
	{
		stparameter[nts[i][0]].resize(2);
		//load2ddistribution(basepath+"/" + nts[i] + "/bin_info.txt", basepath+"/" + nts[i] + "/bin_energies.bin", basepath+"/" + nts[i] + "/bin_binary.bin", basepath+"/" + nts[i] + "/bin_distribution.bin", stparameter[nts[i][0]]);
		loadramachandranplot(basepath+"/" + nts[i] + "/group_0/bin_info.txt", basepath+"/" + nts[i] + "/group_0/bin_rm.bin", stparameter[nts[i][0]][0]);
		loadramachandranplot(basepath+"/" + nts[i] + "/group_0/bin_info.txt", basepath+"/" + nts[i] + "/group_1/bin_rm.bin", stparameter[nts[i][0]][1]);
	}

	return true;
}

bool LoadParams::loadtorclsbb2dparam(string basepath, vector<D2Binnedparamtype> &stparameter){

	loadramachandranplot(basepath+"/bin_info.txt", basepath+"/bin_rm.bin", stparameter);
	
	return true;
}

bool LoadParams::loadsinglenttemplates(string basepath, map<int, map<char, vector<vector<Residueatoms>>>> &stparameter){
	vector<string> nts = {"A","C","G","U"};
	string tdir,pdbpath;
	vector<string>refname;
	FILE *filein;
	string oneline_str;
	char	oneline[600];	
	string st;
	bool flagf;

	for(int ss = 0; ss<2; ss++){
		for(int i=0;i<4;i++)
		{
			tdir=basepath+"/"+to_string(ss)+"/"+nts[i];
			refname.clear();
			for (const auto& entry : fs::directory_iterator(tdir)) {
				if (fs::is_regular_file(entry) && entry.path().extension() == ".pdb") {
					refname.push_back(entry.path());
				}
			}

			stparameter[ss][nts[i][0]].resize(refname.size());
			for(int j=0; j<refname.size();j++){
				vector<Residueatoms> residues;
				Residueatoms residue;
				int seqidx=-1;
				pdbpath=refname[j];

				filein = fopen(pdbpath.c_str(), "r");
				if(filein==NULL)
				{
					printf( "PDB file can not be opened %s\n", pdbpath.c_str());
					return false;
				}

				char oldline[600];oldline[0]='\0';
				bool modelend=false;

				while(!feof(filein) && !modelend)
				{	
					flagf=fgets(oneline,600,filein);
					if(!flagf) break;
					if(strcmp(oneline,oldline)!=0 || oldline[0]=='\0')
					{
						st=oneline;

						if (st.substr(0,3) == "TER")
						{
							residues.push_back(residue);
							modelend=true;
						}else if(st.substr(0,4)=="ATOM")
						{
							int resno = atoi(st.substr(22,4).c_str());
							if(seqidx != -1 & seqidx != resno){
								residues.push_back(residue);
								residue = Residueatoms();
								residue.residx=resno-1;
								residue.aaa=st.substr(19,1)[0];
								residue.resind=Reschar2Idint[residue.aaa];
								seqidx=resno;
							}else if(seqidx == -1){
								residue = Residueatoms();
								residue.residx=resno-1;
								residue.aaa=st.substr(19,1)[0];
								residue.resind=Reschar2Idint[residue.aaa];
								seqidx=resno;
							}

							residue.ptv[Atom2cid[residue.aaa][st.substr(12,4)]].x=float(atof(st.substr(30,8).c_str())); 
							residue.ptv[Atom2cid[residue.aaa][st.substr(12,4)]].y=float(atof(st.substr(38,8).c_str())); 
							residue.ptv[Atom2cid[residue.aaa][st.substr(12,4)]].z=float(atof(st.substr(46,8).c_str())); 
						}
					}
					strcpy(oldline,oneline);
				}

				stparameter[ss][nts[i][0]][j]=residues;
				fclose(filein);
			}
		}
	}

	return true;
}

bool LoadParams::loadstacktemplates(string basepath, map<int, map<char, map<char, vector<vector<Residueatoms>>>>> &stparameter, map<char, map<char, map<char, vector<TemplateDistMap>>>> &distparameter){
	vector<string> nts = {"A","C","G","U"};
	string tdir,pdbpath;
	vector<string>refname;
	FILE *filein;
	string oneline_str;
	char	oneline[600];	
	string st;
	bool flagf;

	for(int ss=0; ss<3; ss++){
		for(int nt1=0;nt1<4;nt1++)
		{
			for(int nt2=0;nt2<4;nt2++)
			{
				if(ss == 0 && nt1 > nt2) continue;
				tdir=basepath+"/" + to_string(ss) + "/"+nts[nt1]+"_"+nts[nt2];
				refname.clear();
				for (const auto& entry : fs::directory_iterator(tdir)) {
					if (fs::is_regular_file(entry) && entry.path().extension() == ".pdb") {
						refname.push_back(entry.path());
					}
				}

				if(ss == 0 && nt1 != nt2){
					stparameter[ss][nts[nt1][0]][nts[nt2][0]].resize(refname.size());
					stparameter[ss][nts[nt2][0]][nts[nt1][0]].resize(refname.size());
					distparameter[ss][nts[nt1][0]][nts[nt2][0]].resize(refname.size());
					distparameter[ss][nts[nt2][0]][nts[nt1][0]].resize(refname.size());
				}
				else{
					stparameter[ss][nts[nt1][0]][nts[nt2][0]].resize(refname.size());
					distparameter[ss][nts[nt1][0]][nts[nt2][0]].resize(refname.size());
				}

				for(int j=0; j<refname.size();j++){
					vector<Residueatoms> residues;
					Residueatoms residue;
					int seqidx=-1;
					pdbpath=refname[j];

					filein = fopen(pdbpath.c_str(), "r");
					if(filein==NULL)
					{
						printf( "PDB file can not be opened %s\n", pdbpath.c_str());
						return false;
					}

					char oldline[600];oldline[0]='\0';
					bool modelend=false;

					while(!feof(filein) && !modelend)
					{	
						flagf=fgets(oneline,600,filein);
						if(!flagf) break;
						if(strcmp(oneline,oldline)!=0 || oldline[0]=='\0')
						{
							st=oneline;

							if (st.substr(0,3) == "TER")
							{
								residues.push_back(residue);
								modelend=true;
							}else if(st.substr(0,4)=="ATOM")
							{
								int resno = atoi(st.substr(22,4).c_str());
								if(seqidx != -1 & seqidx != resno){
									residues.push_back(residue);
									residue = Residueatoms();
									residue.residx=resno-1;
									residue.aaa=st.substr(19,1)[0];
									residue.resind=Reschar2Idint[residue.aaa];
									seqidx=resno;
								}else if(seqidx == -1){
									residue = Residueatoms();
									residue.residx=resno-1;
									residue.aaa=st.substr(19,1)[0];
									residue.resind=Reschar2Idint[residue.aaa];
									seqidx=resno;
								}

								residue.ptv[Atom2cid[residue.aaa][st.substr(12,4)]].x=float(atof(st.substr(30,8).c_str())); 
								residue.ptv[Atom2cid[residue.aaa][st.substr(12,4)]].y=float(atof(st.substr(38,8).c_str())); 
								residue.ptv[Atom2cid[residue.aaa][st.substr(12,4)]].z=float(atof(st.substr(46,8).c_str())); 
							}
						}
						strcpy(oldline,oneline);
					}

					if(ss == 0 && nt1 != nt2){
						vector<Residueatoms> invresidues;
						invresidues.resize(residues.size());
						for(int l=0;l<residues.size();l++){
							invresidues[l]=residues[residues.size()-l-1];
						}
						stparameter[ss][nts[nt1][0]][nts[nt2][0]][j]=residues;
						stparameter[ss][nts[nt2][0]][nts[nt1][0]][j]=invresidues;

						gettemplatedistmap(stparameter[ss][nts[nt1][0]][nts[nt2][0]][j], distparameter[ss][nts[nt1][0]][nts[nt2][0]][j]);
						gettemplatedistmap(stparameter[ss][nts[nt2][0]][nts[nt1][0]][j], distparameter[ss][nts[nt2][0]][nts[nt1][0]][j]);
						residues.clear();

						invresidues.clear();
					}
					else{
						stparameter[ss][nts[nt1][0]][nts[nt2][0]][j]=residues;
						gettemplatedistmap(stparameter[ss][nts[nt1][0]][nts[nt2][0]][j], distparameter[ss][nts[nt1][0]][nts[nt2][0]][j]);

						residues.clear();
					}
					
					fclose(filein);
				}
			}
		}
	}


	return true;
}

void LoadParams::coord2anglen(Point3s tmp2, Point3s tmp3, Point3s tmp4, float &len, float &ang)
{
	BasicFunc bf;
	Point3d p12,p23;
	p12=bf.setv(tmp2.x-tmp3.x,tmp2.y-tmp3.y,tmp2.z-tmp3.z);
	p23=bf.setv(tmp4.x-tmp3.x,tmp4.y-tmp3.y,tmp4.z-tmp3.z);
	len=bf.norm(p23);
	ang=bf.angv(p12,p23);
}

void LoadParams::coord2len(Point3s tmp3, Point3s tmp4, float &len)
{
	BasicFunc bf;
	Point3d p23;
	p23=bf.setv(tmp4.x-tmp3.x,tmp4.y-tmp3.y,tmp4.z-tmp3.z);
	len=bf.norm(p23);
}

void LoadParams::coord2params(Point3s tmp1, Point3s tmp2, Point3s tmp3, Point3s tmp4, float &len, float &ang, float &dih)
{
	BasicFunc bf;
	Point3d p12,p23;
	p12=bf.setv(tmp2.x-tmp3.x,tmp2.y-tmp3.y,tmp2.z-tmp3.z);
	p23=bf.setv(tmp4.x-tmp3.x,tmp4.y-tmp3.y,tmp4.z-tmp3.z);
	len=bf.norm(p23);
	ang=bf.angv(p12,p23);
	dih=bf.phi(tmp1.x,tmp1.y,tmp1.z,
			   tmp2.x,tmp2.y,tmp2.z,
			   tmp3.x,tmp3.y,tmp3.z,
			   tmp4.x,tmp4.y,tmp4.z);
}

bool LoadParams::gettemplategeoparams(vector<Residueatoms> &stparameter){
	int i, j;
	char aat;
	int atom1, atom2, atom3, atom4;
	Point3s tmp1,tmp2,tmp3,tmp4;
	BasicFunc bf;

	stparameter.resize(2);

	aat=stparameter[0].aaa;
	atom3=ParamsList[aat][1][2];atom4=ParamsList[aat][1][3];
	tmp3=stparameter[0].ptv[atom3];tmp4=stparameter[0].ptv[atom4];
	coord2len(tmp3, tmp4, stparameter[0].len[1]);
		
	atom2=ParamsList[aat][2][1];atom3=ParamsList[aat][2][2];atom4=ParamsList[aat][2][3];
	tmp2=stparameter[0].ptv[atom2];tmp3=stparameter[0].ptv[atom3];tmp4=stparameter[0].ptv[atom4];
	coord2anglen(tmp2, tmp3, tmp4, stparameter[0].len[2], stparameter[0].ang[2]);

	for(j=3;j<AccAtomParamNum[aat][0];j++){
		atom1=ParamsList[aat][j][0];atom2=ParamsList[aat][j][1];atom3=ParamsList[aat][j][2];atom4=ParamsList[aat][j][3];
		tmp1=stparameter[0].ptv[atom1];tmp2=stparameter[0].ptv[atom2];tmp3=stparameter[0].ptv[atom3];tmp4=stparameter[0].ptv[atom4];
				
		coord2params(tmp1, tmp2, tmp3, tmp4, stparameter[0].len[j], stparameter[0].ang[j], stparameter[0].dih[j]);
	}

	for(i=1;i<stparameter.size();i++)
	{
		aat=stparameter[i].aaa;
		atom1=ParamsList[aat][0][0];atom2=ParamsList[aat][0][1];atom3=ParamsList[aat][0][2];atom4=ParamsList[aat][0][3];
		tmp1=stparameter[i-1].ptv[atom1];tmp2=stparameter[i-1].ptv[atom2];tmp3=stparameter[i-1].ptv[atom3];tmp4=stparameter[i].ptv[atom4];
		coord2params(tmp1, tmp2, tmp3, tmp4, stparameter[i].len[j], stparameter[i].ang[j], stparameter[i].dih[j]);

		atom1=ParamsList[aat][1][0];atom2=ParamsList[aat][1][1];atom3=ParamsList[aat][1][2];atom4=ParamsList[aat][1][3];
		tmp1=stparameter[i-1].ptv[atom1];tmp2=stparameter[i-1].ptv[atom2];tmp3=stparameter[i].ptv[atom3];tmp4=stparameter[i].ptv[atom4];
		coord2params(tmp1, tmp2, tmp3, tmp4, stparameter[i].len[1], stparameter[i].ang[1], stparameter[i].dih[1]);
		
		atom1=ParamsList[aat][2][0];atom2=ParamsList[aat][2][1];atom3=ParamsList[aat][2][2];atom4=ParamsList[aat][2][3];
		tmp1=stparameter[i-1].ptv[atom1];tmp2=stparameter[i].ptv[atom2];tmp3=stparameter[i].ptv[atom3];tmp4=stparameter[i].ptv[atom4];
		coord2params(tmp1, tmp2, tmp3, tmp4, stparameter[i].len[2], stparameter[i].ang[2], stparameter[i].dih[2]);
		
		for(j=3;j<AccAtomParamNum[aat][0];j++){
			atom1=ParamsList[aat][j][0];atom2=ParamsList[aat][j][1];atom3=ParamsList[aat][j][2];atom4=ParamsList[aat][j][3];
			tmp1=stparameter[i].ptv[atom1];tmp2=stparameter[i].ptv[atom2];tmp3=stparameter[i].ptv[atom3];tmp4=stparameter[i].ptv[atom4];
			
			coord2params(tmp1, tmp2, tmp3, tmp4, stparameter[i].len[j], stparameter[i].ang[j], stparameter[i].dih[j]);
		}
	}
	return true;
}

void LoadParams::gettemplatedistmap(vector<Residueatoms> stparameter, TemplateDistMap &distparameter){
	int mtype1, mtype2;
	Point3s vec;
	double rnorm;

	distparameter.intradistmap.resize(2);
	distparameter.intradistmap[0].resize(PAIRLENGTHTARGETMOVNUM);
	distparameter.intradistmap[1].resize(PAIRLENGTHTARGETMOVNUM);
	for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++){
		distparameter.intradistmap[0][mtype1].resize(PAIRLENGTHTARGETMOVNUM);
		distparameter.intradistmap[1][mtype1].resize(PAIRLENGTHTARGETMOVNUM);
		for(mtype2=mtype1;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++){
			vec.x = stparameter[0].ptv[mtype1].x - stparameter[0].ptv[mtype2].x; 
			vec.y = stparameter[0].ptv[mtype1].y - stparameter[0].ptv[mtype2].y; 
			vec.z = stparameter[0].ptv[mtype1].z - stparameter[0].ptv[mtype2].z;
			distparameter.intradistmap[0][mtype1][mtype2]=bf.norm(vec);
			distparameter.intradistmap[0][mtype2][mtype1]=bf.norm(vec);

			vec.x = stparameter[1].ptv[mtype1].x - stparameter[1].ptv[mtype2].x; 
			vec.y = stparameter[1].ptv[mtype1].y - stparameter[1].ptv[mtype2].y; 
			vec.z = stparameter[1].ptv[mtype1].z - stparameter[1].ptv[mtype2].z;
			distparameter.intradistmap[1][mtype1][mtype2]=bf.norm(vec);
			distparameter.intradistmap[1][mtype2][mtype1]=bf.norm(vec);
		}
	}

	distparameter.interdistmap.resize(PAIRLENGTHTARGETMOVNUM);
	for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++){
		distparameter.interdistmap[mtype1].resize(PAIRLENGTHTARGETMOVNUM);
		for(mtype2=0;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++){
			vec.x = stparameter[0].ptv[mtype1].x - stparameter[1].ptv[mtype2].x; 
			vec.y = stparameter[0].ptv[mtype1].y - stparameter[1].ptv[mtype2].y; 
			vec.z = stparameter[0].ptv[mtype1].z - stparameter[1].ptv[mtype2].z;
			distparameter.interdistmap[mtype1][mtype2]=bf.norm(vec);
		}
	}
}

bool LoadParams::loaddinttemplates(string basepath, map<char, map<char, vector<vector<Residueatoms>>>> &stparameter, map<char, map<char, vector<TemplateDistMap>>> &distparameter){
	vector<string> nts = {"A","C","G","U"};
	string tdir,pdbpath;
	vector<string>refname;
	FILE *filein;
	string oneline_str;
	char	oneline[600];	
	string st;
	bool flagf;

	for(int i=0;i<4;i++)
	{
		for(int k=0;k<4;k++)
		{
			tdir=basepath+"/"+nts[i]+"_"+nts[k];

			refname.clear();
			for (const auto& entry : fs::directory_iterator(tdir)) {
				if (fs::is_regular_file(entry) && entry.path().extension() == ".pdb") {
					refname.push_back(entry.path());
				}
			}

			stparameter[nts[i][0]][nts[k][0]].resize(refname.size());
			distparameter[nts[i][0]][nts[k][0]].resize(refname.size());
			for(int j=0; j<refname.size();j++){
				vector<Residueatoms> residues;
				Residueatoms residue;
				int seqidx=-1;
				pdbpath=refname[j];

				filein = fopen(pdbpath.c_str(), "r");
				if(filein==NULL)
				{
					printf( "PDB file can not be opened %s\n", pdbpath.c_str());
					fclose(filein);
					return false;
				}

				char oldline[600];oldline[0]='\0';
				bool modelend=false;

				while(!feof(filein) && !modelend)
				{	
					flagf=fgets(oneline,600,filein);
					if(!flagf) break;
					if(strcmp(oneline,oldline)!=0 || oldline[0]=='\0')
					{
						st=oneline;

						if (st.substr(0,3) == "TER")
						{
							residues.push_back(residue);
							modelend=true;
						}else if(st.substr(0,4)=="ATOM")
						{
							int resno = atoi(st.substr(22,4).c_str());
							if(seqidx != -1 & seqidx != resno){
								residues.push_back(residue);
								residue = Residueatoms();
								residue.residx=resno-1;
								residue.aaa=st.substr(19,1)[0];
								residue.resind=Reschar2Idint[residue.aaa];
								seqidx=resno;
							}else if(seqidx == -1){
								residue = Residueatoms();
								residue.residx=resno-1;
								residue.aaa=st.substr(19,1)[0];
								residue.resind=Reschar2Idint[residue.aaa];
								seqidx=resno;
							}

							residue.ptv[Atom2cid[residue.aaa][st.substr(12,4)]].x=float(atof(st.substr(30,8).c_str())); 
							residue.ptv[Atom2cid[residue.aaa][st.substr(12,4)]].y=float(atof(st.substr(38,8).c_str())); 
							residue.ptv[Atom2cid[residue.aaa][st.substr(12,4)]].z=float(atof(st.substr(46,8).c_str())); 
						}
					}
					strcpy(oldline,oneline);
				}

				stparameter[nts[i][0]][nts[k][0]][j]=residues;
				gettemplategeoparams(stparameter[nts[i][0]][nts[k][0]][j]);
				gettemplatedistmap(stparameter[nts[i][0]][nts[k][0]][j], distparameter[nts[i][0]][nts[k][0]][j]);

				residues.clear();

				fclose(filein);
			}
		}
	}

	/*
	for(int i=0;i<4;i++)
	{
		cout<< stparameter[nts[i][0]].size()<<endl;
		for(int j=0;j<stparameter[nts[i][0]].size();j++){
			cout<< stparameter[nts[i][0]][j].size()<<endl;
			for(int k=0;k<stparameter[nts[i][0]][j].size();k++){
				cout<< j << "-" << k << "-" << stparameter[nts[i][0]][j][k].residx << "-"  << stparameter[nts[i][0]][j][k].aaa << "-" << stparameter[nts[i][0]][j][k].ptv[0].x << "-" << stparameter[nts[i][0]][j][k].ptv[0].y << "-" << stparameter[nts[i][0]][j][k].ptv[0].z <<endl; 
			}
		}
	}
	*/

	return true;
}

bool LoadParams::loadnptemplates(string basepath, map<char, map<char, vector<vector<Residueatoms>>>> &stparameter, map<char, map<char, vector<TemplateDistMap>>> &distparameter){
	vector<string> nts = {"A","C","G","U"};
	string tdir,pdbpath;
	vector<string>refname;
	FILE *filein;
	string oneline_str;
	char	oneline[600];	
	string st;
	bool flagf;

	for(int i=0;i<4;i++)
	{
		for(int k=i;k<4;k++)
		{
			tdir=basepath+"/"+nts[i]+"_"+nts[k];

			refname.clear();
			for (const auto& entry : fs::directory_iterator(tdir)) {
				if (fs::is_regular_file(entry) && entry.path().extension() == ".pdb") {
					refname.push_back(entry.path());
				}
			}

			stparameter[nts[i][0]][nts[k][0]].resize(refname.size());
			stparameter[nts[k][0]][nts[i][0]].resize(refname.size());
			distparameter[nts[i][0]][nts[k][0]].resize(refname.size());
			distparameter[nts[k][0]][nts[i][0]].resize(refname.size());

			for(int j=0; j<refname.size();j++){
				vector<Residueatoms> residues;
				Residueatoms residue;
				int seqidx=-1;
				pdbpath=refname[j];

				filein = fopen(pdbpath.c_str(), "r");
				if(filein==NULL)
				{
					printf( "PDB file can not be opened %s\n", pdbpath.c_str());
					return false;
				}

				char oldline[600];oldline[0]='\0';
				bool modelend=false;

				while(!feof(filein) && !modelend)
				{	
					flagf=fgets(oneline,600,filein);
					if(!flagf) break;
					if(strcmp(oneline,oldline)!=0 || oldline[0]=='\0')
					{
						st=oneline;

						if (st.substr(0,3) == "TER")
						{
							residues.push_back(residue);
							modelend=true;
						}else if(st.substr(0,4)=="ATOM")
						{
							int resno = atoi(st.substr(22,4).c_str());
							if(seqidx != -1 & seqidx != resno){
								residues.push_back(residue);
								residue = Residueatoms();
								residue.residx=resno-1;
								residue.aaa=st.substr(19,1)[0];
								residue.resind=Reschar2Idint[residue.aaa];
								seqidx=resno;
							}else if(seqidx == -1){
								residue = Residueatoms();
								residue.residx=resno-1;
								residue.aaa=st.substr(19,1)[0];
								residue.resind=Reschar2Idint[residue.aaa];
								seqidx=resno;
							}

							residue.ptv[Atom2cid[residue.aaa][st.substr(12,4)]].x=float(atof(st.substr(30,8).c_str())); 
							residue.ptv[Atom2cid[residue.aaa][st.substr(12,4)]].y=float(atof(st.substr(38,8).c_str())); 
							residue.ptv[Atom2cid[residue.aaa][st.substr(12,4)]].z=float(atof(st.substr(46,8).c_str())); 
						}
					}
					strcpy(oldline,oneline);
				}

				vector<Residueatoms> invresidues;
				invresidues.resize(residues.size());
				for(int l=0;l<residues.size();l++){
					invresidues[l]=residues[residues.size()-l-1];
				}

				stparameter[nts[i][0]][nts[k][0]][j]=residues;
				stparameter[nts[k][0]][nts[i][0]][j]=invresidues;

				gettemplategeoparams(stparameter[nts[i][0]][nts[k][0]][j]);
				gettemplategeoparams(stparameter[nts[k][0]][nts[i][0]][j]);

				gettemplatedistmap(stparameter[nts[i][0]][nts[k][0]][j], distparameter[nts[i][0]][nts[k][0]][j]);
				gettemplatedistmap(stparameter[nts[k][0]][nts[i][0]][j], distparameter[nts[k][0]][nts[i][0]][j]);

				residues.clear();
				invresidues.clear();
				fclose(filein);
			}
		}
	}

	/*
	for(int i=0;i<4;i++)
	{
		cout<< stparameter[nts[i][0]].size()<<endl;
		for(int j=0;j<stparameter[nts[i][0]].size();j++){
			cout<< stparameter[nts[i][0]][j].size()<<endl;
			for(int k=0;k<stparameter[nts[i][0]][j].size();k++){
				cout<< j << "-" << k << "-" << stparameter[nts[i][0]][j][k].residx << "-"  << stparameter[nts[i][0]][j][k].aaa << "-" << stparameter[nts[i][0]][j][k].ptv[0].x << "-" << stparameter[nts[i][0]][j][k].ptv[0].y << "-" << stparameter[nts[i][0]][j][k].ptv[0].z <<endl; 
			}
		}
	}
	*/

	return true;
}


bool LoadParams::loadoxhydbparams(string basepath, map<char, vector<D2Binnedparamtype>> &stparameter){
	vector<string> nts = {"A","C","G","U"};
	vector<int> param_num = {9, 8, 10, 8};
	
	for(int i=0;i<4;i++)
	{
		loadramachandranplot(basepath+"/" + nts[i] + "/bin_info.txt", basepath+"/" + nts[i] + "/bin_rm.bin", stparameter[nts[i][0]]);
	}

	return true;
}

void LoadParams::paramloadmain(string basePath, int numbb, string seq)
{
	bool flagps;

	flagps=loadstandardparameters(basePath+"/geo", TORSIONPARAM, "torsion", "geo");
	if(flagps==false) exit(1);
	flagps=loadstandardparameters(basePath+"/geo", INNERANGPARAM, "inner", "geo");
	if(flagps==false) exit(1);
	flagps=loadstandardparameters(basePath+"/geo", LENGTHPARAM, "length", "geo");
	if(flagps==false) exit(1);

	flagps=loadstandardparameters(basePath+"/oppleng", SUBLENGTHPARAM, "oppleng", "geo");
	if(flagps==false) exit(1);

	flagps=loadstandardparameters(basePath+"/geo_add", ADDINNERANGPARAM, "inner", "geo_add");
	if(flagps==false) exit(1);
	flagps=loadstandardparameters(basePath+"/geo_add", ADDLENGTHPARAM, "length", "geo_add");
	if(flagps==false) exit(1);

	flagps=loadcombtorsionparameter(basePath+"/torcomb", MULTITORSIONPARAM);
	if(flagps==false) exit(1);

	flagps=loadstandardparameters(basePath+"/geo_first", OPTORSIONPARAM, "torsion", "geo_first");
	if(flagps==false) exit(1);
	flagps=loadstandardparameters(basePath+"/geo_first", OPINNERANGPARAM, "inner", "geo_first");
	if(flagps==false) exit(1);
	flagps=loadstandardparameters(basePath+"/geo_first", OPLENGTHPARAM, "length", "geo_first");
	if(flagps==false) exit(1);

	flagps=loadstandardparameters(basePath+"/geo_hyd", HYDTORSIONPARAM, "torsion", "geo_hyd");
	if(flagps==false) exit(1);
	flagps=loadstandardparameters(basePath+"/geo_hyd", HYDINNERANGPARAM, "inner", "geo_hyd");
	if(flagps==false) exit(1);
	flagps=loadstandardparameters(basePath+"/geo_hyd", HYDLENGTHPARAM, "length", "geo_hyd");
	if(flagps==false) exit(1);

	flagps=loadstackingrot(basePath+"/stacking", STACKINGROT);
	if(flagps==false) exit(1);
	flagps=loadbasepairrot(basePath+"/basepair", BASEPAIRROT);
	if(flagps==false) exit(1);

	flagps=loadbaepairhydrogenparameters(basePath+"/basepair_hyd", BASEPAIRHYDPARAM);
	if(flagps==false) exit(1);

	flagps=loadstackingparameters(basePath+"/stacking", STACKINGPARAM);
	if(flagps==false) exit(1);
	flagps=loadbasepairparameters(basePath+"/basepair", BASEPAIRPARAM);
	if(flagps==false) exit(1);

	//Parameters in clash
	flagps=loadvwd(basePath+"/clash/vdw.txt");
	if(flagps==false) exit(1);

	flagps=loadclashexpatoms(basePath+"/clash", EXPCLASHATOMSNT, EXPCLASHATOMSADJ, EXPCLASHATOMSANY);
	if(flagps==false) exit(1);

	flagps=loadbasetemplate(basePath+"/base_template", BASETEMPLATE);
	if(flagps==false) exit(1);

	flagps=loadtorsionlist(basePath+"/torsion_list", TORSIONLIST);
	if(flagps==false) exit(1);

	flagps=loadpuckerparameters(basePath+"/pucker", PUCKERPARAM);
	if(flagps==false) exit(1);

	flagps=loadpairlengnt(basePath+"/pairleng/nt_bb", PAIRLENGNTPARAM);
	if(flagps==false) exit(1);

	flagps=loadpairlengadj(basePath+"/pairleng/adj_bb", PAIRLENGADJPARAM);
	if(flagps==false) exit(1);

	flagps=loadsinglepairlengnt(basePath+"/pairleng/nt_bb", ALLPAIRLENGNTPARAM);
	if(flagps==false) exit(1);

	flagps=loadtorclsparam(basePath+"/cls_tor", TORCLSPARAM);
	if(flagps==false) exit(1);

	flagps=loadtorcls2dparam(basePath+"/cls_tor_2d", TORSIONPAIRSPARAM);
	if(flagps==false) exit(1);

	flagps=loadtorclsbb2dparam(basePath+"/cls_tor_2d_bb", TORSIONBBPAIRSPARAM);
	if(flagps==false) exit(1);

	flagps=loadbppairlengnt(basePath+"/pairelengbp/nt", BPATOMPAIRLENGNTPARAM);
	if(flagps==false) exit(1);

	flagps=loadbppairlengbp(basePath+"/pairelengbp/bp", BPATOMPAIRLENGPARAM);
	if(flagps==false) exit(1);

	flagps=loadbppairlengadj(basePath+"/pairelengbp/adj", BPATOMPAIRLENGADJPARAM);
	if(flagps==false) exit(1);

	flagps=loadsinglenttemplates(basePath+"/single_nt_template", SINGLENTTEMPLATE);
	if(flagps==false) exit(1);

	flagps=loaddinttemplates(basePath+"/di_nt_templates", DILENTTEMPLATE, DILENTTEMPLATELENGTH);
	if(flagps==false) exit(1);

	flagps=loadnptemplates(basePath+"/bp_nt_templates", BPLENTTEMPLATE, BPLENTTEMPLATELENGTH);
	if(flagps==false) exit(1);

	flagps=loadoxhydbparams(basePath+"/hbond_o2'", OXHYDBPARAM);
	if(flagps==false) exit(1);

	flagps=loadstacktemplates(basePath+"/stacking_template", STACKTEMPLATE, STLENTTEMPLATELENGTH);
	if(flagps==false) exit(1);
}

LoadParams::~LoadParams()
{

}

