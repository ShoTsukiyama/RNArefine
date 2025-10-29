#include "MonteCarloSimulation.hpp"

map<char, vector<vector<double>>> MCLENGTHPARAM;
map<char, vector<vector<double>>> MCINNERANGPARAM;
map<char, vector<vector<double>>> MCTORSIONPARAM;

map<char, vector<vector<double>>> MCBASELENGTHPARAM;
map<char, vector<vector<double>>> MCBASEINNERANGPARAM;
map<char, vector<vector<double>>> MCBASETORSIONPARAM;

map<char, vector<vector<double>>> MCADDINNERANGPARAM;
map<char, vector<vector<double>>> MCADDLENGTHPARAM;

map<char, vector<vector<double>>> MCOPTORSIONPARAM;
map<char, vector<vector<double>>> MCOPINNERANGPARAM;
map<char, vector<vector<double>>> MCOPLENGTHPARAM;

map<char, vector<vector<double>>> MCHYDTORSIONPARAM;
map<char, vector<vector<double>>> MCHYDINNERANGPARAM;
map<char, vector<vector<double>>> MCHYDLENGTHPARAM;

map<char, vector<vector<double>>> MCHYDADDLENGTHPARAM;
map<char, vector<vector<double>>> MCHYDADDANGLEPARAM;

map<char, vector<double>> MCVWPARAM;

map<char, vector<Binnedparamtype>> MCTORSIONDISTRIBUTION;
map<char, vector<Binnedparamtype>> MCBASETORSIONDISTRIBUTION;
map<char, vector<Binnedparamtype>> MCOPTORSIONDISTRIBUTION;
map<char, vector<Binnedparamtype>> MCHYDTORSIONDISTRIBUTION;

map<char, vector<D2Binnedparamtype>> TORSIONPAIRSPARAM;

map<char, vector<vector<PairLengInfo>>> ALLPAIRLENGNTPARAM;
map<char, map<char, vector<vector<PairLengInfo>>>> PAIRLENGADJPARAM;

map<int, map<char, map<char, vector<vector<Residueatoms>>>>> STACKTEMPLATE;
map<int, map<char, map<char, vector<vector<double>>>>> STACKPARAM;
map<char, map<char, vector<vector<double>>>> BPPARAM;
map<int, map<char, map<char, vector<vector<vector<double>>>>>> STACKDISTMAP;
map<char, map<char, vector<vector<vector<double>>>>> BPDISTMAP;
map<char, map<char, vector<vector<Residueatoms>>>> BPTEMPLATE;
map<char, vector<vector<Residueatoms>>> LOOPNTTEMPLATE;
map<char, vector<vector<Residueatoms>>> NTTEMPLATE;
map<char, map<char, map<int, map<int, HydroBondParameters>>>> BPHBOND;

map<string, vector<vector<int>>> molConnectone;
map<string, vector<vector<int>>> molConnecttwo; 
map<string, vector<vector<int>>> molConnectthree;
map<string, vector<vector<int>>> molAdjConnectone;
map<string, vector<vector<int>>> molAdjConnecttwo;
map<string, vector<vector<int>>> molAdjConnectthree;
map<string, vector<vector<int>>> molAdjConnectoneidx;
map<string, vector<vector<int>>> molAdjConnecttwoidx;
map<string, vector<vector<int>>> molAdjConnectthreeidx;
map<string, map<int, vector<double>>> AtomFeatures;

map<string, vector<vector<bool>>> intraAtomConnect;
map<string, map<string, vector<vector<int>>>> intramolConnect;
map<string, vector<vector<int>>> ntmolConnect;

void loadparamvec2dd(string filepath, size_t dim1_size, size_t dim2_size, vector<vector<double>> &vec_2d){

    std::ifstream input(filepath, std::ios::binary);
    if (!input) {
        std::cerr << "Error opening file (" + filepath + ")!" << std::endl;
        exit(1);
    }

    input.seekg(0, std::ios::end);
    std::streamsize size = input.tellg();
    input.seekg(0, std::ios::beg);
    vec_2d.resize(dim1_size, vector<double>(dim2_size));
	for (size_t i = 0; i < dim1_size; ++i) {
        if (!input.read(reinterpret_cast<char*>(vec_2d[i].data()), dim2_size * sizeof(double))) {
            std::cerr << "Error reading file!" << std::endl;
            exit(1);
        }
    }
}

void loadparamvec2di(string filepath, size_t dim1_size, size_t dim2_size, vector<vector<int>> &vec_2d){

    std::ifstream input(filepath, std::ios::binary);
    if (!input) {
        std::cerr << "Error opening file (" + filepath + ")!" << std::endl;
        exit(1);
    }

    input.seekg(0, std::ios::end);
    std::streamsize size = input.tellg();
    input.seekg(0, std::ios::beg);

	vector<vector<double>> vec_2d_tmp;
    vec_2d.resize(dim1_size, vector<int>(dim2_size));
	vec_2d_tmp.resize(dim1_size, vector<double>(dim2_size));

	for (size_t i = 0; i < dim1_size; ++i) {
        if (!input.read(reinterpret_cast<char*>(vec_2d_tmp[i].data()), dim2_size * sizeof(double))) {
            std::cerr << "Error reading file!" << std::endl;
            exit(1);
        }
    }

	for(int i=0;i<vec_2d.size();i++){
		for(int j=0;j<vec_2d[i].size();j++){
			vec_2d[i][j]=vec_2d_tmp[i][j];
		}
	}
}

bool loaddistributionparam(string basepath, map<char, vector<Binnedparamtype>> &distributionparam, string paramname)
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
			cout<<"Bin info file can not be opened (" + paramname + ")\n"<<endl;
			return false;
		}

		while(getline(inputFile, oneline))
		{
			stringstream ssbase(oneline);
			count=0;
			while (getline(ssbase, field, ','))
			{
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

			if(distributionparam[nts[i][0]].size() < paramtype+1) distributionparam[nts[i][0]].resize(paramtype+1);
			distributionparam[nts[i][0]][paramtype].minval=minval;
			distributionparam[nts[i][0]][paramtype].maxval=maxval;
			distributionparam[nts[i][0]][paramtype].span=span;
			distributionparam[nts[i][0]][paramtype].binnum=binnum;
			distributionparam[nts[i][0]][paramtype].mean=mean;
			distributionparam[nts[i][0]][paramtype].std=stdval;
			distributionparam[nts[i][0]][paramtype].mode=mode;

			distributionparam[nts[i][0]][paramtype].upper=mean+3.0*stdval;
			distributionparam[nts[i][0]][paramtype].lower=mean-3.0*stdval;
			
			distributionparam[nts[i][0]][paramtype].rangetype=0;
			if(paramname == "torsion"){
				if(distributionparam[nts[i][0]][paramtype].upper > 360){
					if(distributionparam[nts[i][0]][paramtype].lower < 0){
						distributionparam[nts[i][0]][paramtype].lower=0;
						distributionparam[nts[i][0]][paramtype].upper=360;
					}else{
						distributionparam[nts[i][0]][paramtype].upper = distributionparam[nts[i][0]][paramtype].upper-360;
						distributionparam[nts[i][0]][paramtype].rangetype=1;
					}
				}
				else if(distributionparam[nts[i][0]][paramtype].lower < 0){
					distributionparam[nts[i][0]][paramtype].lower = distributionparam[nts[i][0]][paramtype].lower+360;
					distributionparam[nts[i][0]][paramtype].rangetype=1;
				}
			}

			distributionparam[nts[i][0]][paramtype].energy.resize(binnum);
			distributionparam[nts[i][0]][paramtype].distr.resize(binnum);
			distributionparam[nts[i][0]][paramtype].boundary.resize(binnum-1);
			distributionparam[nts[i][0]][paramtype].center.resize(binnum);
			paramnum+=1;
		}

		vector<vector<double>> vec2denergies;
		loadparamvec2dd(basepath+"/"+nts[i][0]+"/bin_energies_" + paramname + ".bin", paramnum, binnum, vec2denergies);
		vector<vector<int>> vec2ddistribution;
		loadparamvec2di(basepath+"/"+nts[i][0]+"/bin_distribution_" + paramname + ".bin", paramnum, binnum, vec2ddistribution);
		
		for(int j=0;j<paramnum;j++){
			distributionparam[nts[i][0]][j].energy=vec2denergies[j];
			distributionparam[nts[i][0]][j].total=0;
			distributionparam[nts[i][0]][j].totalenergy=0;
			for(int k=0; k<vec2ddistribution[j].size();k++){
				distributionparam[nts[i][0]][j].distr[k]=vec2ddistribution[j][k];
				distributionparam[nts[i][0]][j].total+=vec2ddistribution[j][k];
				distributionparam[nts[i][0]][j].totalenergy+=vec2denergies[j][k];
			}

			for(int k=0; k<binnum-1; k++) distributionparam[nts[i][0]][j].boundary[k]=distributionparam[nts[i][0]][j].minval+k*distributionparam[nts[i][0]][j].span;
			for(int k=0; k<binnum; k++) distributionparam[nts[i][0]][j].center[k]=distributionparam[nts[i][0]][j].minval+(k-0.5)*distributionparam[nts[i][0]][j].span;
		}
	}
	return true;

}

bool loadramachandranplot(string binfilepath, map<char, vector<D2Binnedparamtype>> &stparameter){
	string	oneline;
	int paramtype, paramnum;
	double minval1, maxval1, span1, binnum1, minval2, maxval2, span2, binnum2;
	int count;
	string field;
	vector<string> nts = {"A","C","G","U"};

	for(int i=0;i<4;i++)
	{
		ifstream inputFile(binfilepath+"/"+nts[i]+"/bin_info.txt");
		
		if(!inputFile.is_open())
		{
			cout<<"Bin info file can not be opened (" + binfilepath+"/"+nts[i]+"/bin_info.txt" + ")\n"<<endl;
			return false;
		}

		paramnum=0;
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

			if(stparameter[nts[i][0]].size() < paramtype+1) stparameter[nts[i][0]].resize(paramtype+1);
			stparameter[nts[i][0]][paramtype].minval1=minval1;
			stparameter[nts[i][0]][paramtype].maxval1=maxval1;
			stparameter[nts[i][0]][paramtype].span1=span1;
			stparameter[nts[i][0]][paramtype].binnum1=binnum1;
			stparameter[nts[i][0]][paramtype].minval2=minval2;
			stparameter[nts[i][0]][paramtype].maxval2=maxval2;
			stparameter[nts[i][0]][paramtype].span2=span2;
			stparameter[nts[i][0]][paramtype].binnum2=binnum2;

			stparameter[nts[i][0]][paramtype].energy.resize(binnum1*binnum2);
			stparameter[nts[i][0]][paramtype].binary.resize(binnum1*binnum2);
			stparameter[nts[i][0]][paramtype].distr.resize(binnum1*binnum2);
			stparameter[nts[i][0]][paramtype].boundary1.resize(binnum1-1);
			stparameter[nts[i][0]][paramtype].center1.resize(binnum1);
			stparameter[nts[i][0]][paramtype].boundary2.resize(binnum2-1);
			stparameter[nts[i][0]][paramtype].center2.resize(binnum2);
			paramnum+=1;
		}
		
		vector<vector<double>> vec2denergies;
		loadparamvec2dd(binfilepath+"/"+nts[i]+"/bin_energies.bin", paramnum, binnum1*binnum2, vec2denergies);
		
		for(int j=0;j<paramnum;j++)
		{
			stparameter[nts[i][0]][j].energy=vec2denergies[j];
			for(int k=0; k<binnum1-1; k++) stparameter[nts[i][0]][j].boundary1[k]=stparameter[nts[i][0]][j].minval1+k*stparameter[nts[i][0]][j].span1;
			for(int k=0; k<binnum2-1; k++) stparameter[nts[i][0]][j].boundary2[k]=stparameter[nts[i][0]][j].minval2+k*stparameter[nts[i][0]][j].span2;
			for(int k=0; k<binnum1; k++) stparameter[nts[i][0]][j].center1[k]=stparameter[nts[i][0]][j].minval1+(k-0.5)*stparameter[nts[i][0]][j].span1;
			for(int k=0; k<binnum2; k++) stparameter[nts[i][0]][j].center2[k]=stparameter[nts[i][0]][j].minval2+(k-0.5)*stparameter[nts[i][0]][j].span2;
		
			//for(int k=0; k<binnum1-1; k++){
				//for(int l=0; l<binnum1-1; l++){
					//cout<<stparameter[nts[i][0]][j].energy[k+binnum1*l]<<endl;
				//}
			//}
		}
	}

	return true;
}

bool loadgeometricalparameters(string basepath, map<char, vector<vector<double>>> &parameter, string paramname, double k)
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

            double lower = mean-k*stdval;
            lower = max(lower, 0.0);
            double upper = mean+k*stdval;
            parameter[nts[i][0]].push_back({mean, stdval, lower, upper, mode});
        }

    }
    return true;
}

bool loadvwd(string basepath, map<char, vector<double>> &parameter)
{
	int i, j;
	FILE *file;
	vector<string> nts = {"A","C","G","U"};
	char	oneline[600];
	char atomname[5];
	atomname[0]=' ';atomname[1]=' ';atomname[2]=' ';atomname[3]=' ';atomname[4]='\0';
	string atomnamest;
	double val;
	int atomid;
	bool flagps;

	for(i=0;i<nts.size();i++){
		string nt = nts[i];

		file = fopen((basepath+"/vdw_"+nt+".txt").c_str(), "rt");
		if(file==NULL)
		{
			printf( "Parameter file can not be opened %s\n", basepath.c_str());
			return false;
		}	
		parameter[nt[0]].resize(AccAtomParamNum[nt[0]][4]);
		for(j=0; j<AccAtomParamNum[nt[0]][4];j++) parameter[nt[0]][j]=-1;	

		int idx=0;
		while(!feof(file))
		{
			flagps=fgets(oneline,600,file);
			if(!flagps) break;
			sscanf(oneline,"%c%c%c%c%lf", &atomname[0], &atomname[1], &atomname[2], &atomname[3], &val);
			parameter[nt[0]][idx]=val;
			idx+=1;
		}
		fclose(file);
	}

	return true;
}

bool loadpairlengadj(string basepath, map<char, map<char, vector<vector<PairLengInfo>>>> &stparameter){
	int i, j;
	vector<string> nts = {"A","U","G","C"};
	string oneline;
	int count, paramidx;
	int atom1, atom2;
	string field;
	double mean, std;

	for(i=0; i<nts.size();i++){
		for(j=0; j<nts.size();j++){
			
			stparameter[nts[i][0]][nts[j][0]].resize(AccAtomParamNum[nts[i][0]][2], vector<PairLengInfo>(AccAtomParamNum[nts[i][0]][2]));

			ifstream inputFile(basepath+"/adj/adj_bbsugar_pair_length_"+nts[i]+nts[j]+".csv");

			if(!inputFile.is_open())
			{
				cout<<"CSV info file can not be opened (grouped pair atom lengths)\n"<<endl;
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

				atom1=paramidx/AccAtomParamNum[nts[i][0]][2];
				atom2=paramidx%AccAtomParamNum[nts[j][0]][2];

				stparameter[nts[i][0]][nts[j][0]][atom1][atom2].mean=mean;
				stparameter[nts[i][0]][nts[j][0]][atom1][atom2].std=std;
				stparameter[nts[i][0]][nts[j][0]][atom1][atom2].upper=mean+3.0*std;
				
				stparameter[nts[i][0]][nts[j][0]][atom1][atom2].lower=min(0.0, mean-3.0*std);

				paramidx+=1;
			}
		}
	}
	return true;
}

bool loadpairlengnt(string basepath, map<char, vector<vector<PairLengInfo>>> &stparameter){
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

	for(i=0; i<nts.size();i++)
	{
		stparameter[nts[i][0]].resize(AccAtomParamNum[nts[i][0]][2], vector<PairLengInfo>(AccAtomParamNum[nts[i][0]][2]));

		ifstream inputFile(basepath+"/nt/nt_bbsugar_pair_length_"+nts[i]+".csv");

		if(!inputFile.is_open())
		{
			cout<<"CSV info file can not be opened (grouped adjacent pair atom lengths)\n"<<endl;
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

			atom1=paramidx/AccAtomParamNum[nts[i][0]][2];
			atom2=paramidx%AccAtomParamNum[nts[i][0]][2];

			stparameter[nts[i][0]][atom1][atom2].mean=mean;
			stparameter[nts[i][0]][atom1][atom2].std=std;
			stparameter[nts[i][0]][atom1][atom2].upper=mean+3.0*std;
			stparameter[nts[i][0]][atom1][atom2].lower=min(0.0,mean-3.0*std);

			paramidx+=1;
		}
	}
	return true;
}


bool loadstacktemplates(string basepath, map<int, map<char, map<char, vector<vector<Residueatoms>>>>> &stparameter){
	vector<string> nts = {"A","C","G","U"};
	string tdir,pdbpath;
	vector<string>refname;
	FILE *filein;
	string oneline_str;
	char	oneline[600];	
	string st;
	bool flagf;

	map< char, int> Reschar2Idint={{'A',0},{'C',1},{'G',2},{'U',3},{'X',4}};

	map<char, map<string, int>> Atom2cid={
		{'A',{{" P  ", 0},{" O5'", 1},{" C5'", 2},{" C4'", 3},{" C3'", 4},{" O3'", 5},{" C2'", 6},{" O4'", 7},{" C1'", 8},{" O2'",9},{" OP1",10},{" OP2",11},{" N9 ",12},{" C4 ",13},{" C8 ",14}, {" N1 ",15}, {" C2 ",16}, {" C6 ",17}, {" N7 ",18}, {" C5 ",19}, {" N6 ",20}, {" N3 ",21}}},
		{'C',{{" P  ", 0},{" O5'", 1},{" C5'", 2},{" C4'", 3},{" C3'", 4},{" O3'", 5},{" C2'", 6},{" O4'", 7},{" C1'", 8},{" O2'",9},{" OP1",10},{" OP2",11},{" N1 ",12},{" C2 ",13},{" C6 ",14}, {" C4 ",15}, {" C5 ",16}, {" N3 ",17}, {" O2 ",18}, {" N4 ",19}}},
		{'G',{{" P  ", 0},{" O5'", 1},{" C5'", 2},{" C4'", 3},{" C3'", 4},{" O3'", 5},{" C2'", 6},{" O4'", 7},{" C1'", 8},{" O2'",9},{" OP1",10},{" OP2",11},{" N9 ",12},{" C4 ",13},{" C8 ",14}, {" N1 ",15}, {" C2 ",16}, {" C6 ",17}, {" N7 ",18}, {" C5 ",19}, {" N3 ",20}, {" O6 ",21}, {" N2 ",22}}},
		{'U',{{" P  ", 0},{" O5'", 1},{" C5'", 2},{" C4'", 3},{" C3'", 4},{" O3'", 5},{" C2'", 6},{" O4'", 7},{" C1'", 8},{" O2'",9},{" OP1",10},{" OP2",11},{" N1 ",12},{" C2 ",13},{" C6 ",14}, {" C4 ",15}, {" C5 ",16}, {" N3 ",17}, {" O2 ",18}, {" O4 ",19}}}
	};

	for(int ss=0; ss<2; ss++){
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
				}
				else{
					stparameter[ss][nts[nt1][0]][nts[nt2][0]].resize(refname.size());
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

						residues.clear();
						invresidues.clear();
					}
					else{
						stparameter[ss][nts[nt1][0]][nts[nt2][0]][j]=residues;

						residues.clear();
					}
					
					fclose(filein);
				}
			}
		}
	}


	return true;
}

bool loadstackparams(string basepath, map<int, map<char, map<char, vector<vector<double>>>>> &stparameter){
	vector<string> nts = {"A","U","G","C"};
	string tdir,pdbpath;
	vector<string>refname;
	FILE *filein;
	string oneline_str;
	string	oneline;	
	string st;
	int count;
	int paramtype, binnum;
	double minval, maxval, span, mean, stdval, mode;
	string field;
	bool flagf;

	map< char, int> Reschar2Idint={{'A',0},{'C',1},{'G',2},{'U',3}};

	for(int ss=0; ss<2; ss++){
		for(int nt1=0;nt1<4;nt1++)
		{
			for(int nt2=nt1;nt2<4;nt2++)
			{
				stparameter[ss][nts[nt1][0]][nts[nt2][0]].resize(2);
				stparameter[ss][nts[nt2][0]][nts[nt1][0]].resize(2);
				for(int i=0;i<2;i++){
					stparameter[ss][nts[nt1][0]][nts[nt2][0]][i].resize(4);
					stparameter[ss][nts[nt2][0]][nts[nt1][0]][i].resize(4);
				}
			}
		}
	}

	for(int ss=0; ss<2; ss++){
		for(int nt1=0;nt1<4;nt1++)
		{
			for(int nt2=0;nt2<4;nt2++)
			{
				if(ss == 0 && nt1 > nt2) continue;

				tdir=basepath+"/" + to_string(ss) + "/"+nts[nt1]+"_"+nts[nt2];
				
				ifstream inputFile_leng(tdir+"/length/bin_info.txt");
				if(!inputFile_leng.is_open())
				{
					cout<<"Bin info file can not be opened (Stacking parameters [length])\n"<<endl;
					return false;
				}

				while(getline(inputFile_leng, oneline))
				{
					stringstream ssbase(oneline);
					count=0;
					while (getline(ssbase, field, ','))
					{
						if(count==0) paramtype=stoi(field);
						if(count==1) minval=stof(field); //minval
						if(count==2) maxval=stof(field); //maxval
						if(count==3) span=stof(field); //span
						if(count==4) binnum=stoi(field); //binnum
						if(count==5) mean=stof(field); //mean
						if(count==6) stdval=stof(field); //stdval
						count+=1;
					}

					if(ss == 0){
						stparameter[ss][nts[nt1][0]][nts[nt2][0]][0][0] = mean;
						stparameter[ss][nts[nt1][0]][nts[nt2][0]][0][1] = stdval;
						stparameter[ss][nts[nt1][0]][nts[nt2][0]][0][2] = max(0.0, mean-4*stdval);
						stparameter[ss][nts[nt1][0]][nts[nt2][0]][0][3] = mean+2*stdval;
						stparameter[ss][nts[nt2][0]][nts[nt1][0]][0][0] = mean;
                        stparameter[ss][nts[nt2][0]][nts[nt1][0]][0][1] = stdval;
						stparameter[ss][nts[nt2][0]][nts[nt1][0]][0][2] = max(0.0, mean-4*stdval);
						stparameter[ss][nts[nt2][0]][nts[nt1][0]][0][3] = mean+2*stdval;

					}else{
						stparameter[ss][nts[nt1][0]][nts[nt2][0]][0][0] = mean;
						stparameter[ss][nts[nt1][0]][nts[nt2][0]][0][1] = stdval;
						stparameter[ss][nts[nt1][0]][nts[nt2][0]][0][2] = max(0.0, mean-4*stdval);
						stparameter[ss][nts[nt1][0]][nts[nt2][0]][0][3] = mean+2*stdval;
					}
				}

				ifstream inputFile_ang(tdir+"/angle/bin_info.txt");
				if(!inputFile_ang.is_open())
				{
					cout<<"Bin info file can not be opened (Stacking parameters [angle])\n"<<endl;
					return false;
				}

				while(getline(inputFile_ang, oneline))
				{
					stringstream ssbase(oneline);
					count=0;
					while (getline(ssbase, field, ','))
					{
						if(count==0) paramtype=stoi(field);
						if(count==1) minval=stof(field); //minval
						if(count==2) maxval=stof(field); //maxval
						if(count==3) span=stof(field); //span
						if(count==4) binnum=stoi(field); //binnum
						if(count==5) mean=stof(field); //mean
						if(count==6) stdval=stof(field); //stdval
						count+=1;
					}

					if(ss == 0){
						stparameter[ss][nts[nt1][0]][nts[nt2][0]][1][0] = mean;
						stparameter[ss][nts[nt1][0]][nts[nt2][0]][1][1] = stdval;
						stparameter[ss][nts[nt1][0]][nts[nt2][0]][1][2] = max(0.0, mean-4*stdval);
						stparameter[ss][nts[nt1][0]][nts[nt2][0]][1][3] = min(90.0, mean+4*stdval);
						stparameter[ss][nts[nt2][0]][nts[nt1][0]][1][0] = mean;
                        stparameter[ss][nts[nt2][0]][nts[nt1][0]][1][1] = stdval;
						stparameter[ss][nts[nt2][0]][nts[nt1][0]][1][2] = max(0.0, mean-4*stdval);
						stparameter[ss][nts[nt2][0]][nts[nt1][0]][1][3] = min(90.0, mean+4*stdval);
					}else{
						stparameter[ss][nts[nt1][0]][nts[nt2][0]][1][0] = mean;
						stparameter[ss][nts[nt1][0]][nts[nt2][0]][1][1] = stdval;
						stparameter[ss][nts[nt1][0]][nts[nt2][0]][1][2] = max(0.0, mean-4*stdval);
						stparameter[ss][nts[nt1][0]][nts[nt2][0]][1][3] = min(90.0, mean+4*stdval);
					}
				}
			}
		}
	}

	return true;
}

bool loadbpparams(string basepath, map<char, map<char, vector<vector<double>>>> &stparameter){
	vector<string> nts = {"A","U","G","C"};
	string tdir,pdbpath;
	vector<string>refname;
	FILE *filein;
	string oneline_str;
	string	oneline;	
	string st;
	int count;
	double mean, stdval;
	string field;
	bool flagf;

	for(int nt1=0;nt1<4;nt1++)
	{
		for(int nt2=nt1;nt2<4;nt2++)
		{
			stparameter[nts[nt1][0]][nts[nt2][0]].resize(1);
			stparameter[nts[nt2][0]][nts[nt1][0]].resize(1);
			for(int i=0;i<1;i++){
				stparameter[nts[nt1][0]][nts[nt2][0]][i].resize(2);
				stparameter[nts[nt2][0]][nts[nt1][0]][i].resize(2);
			}
		}
	}

	for(int nt1=0;nt1<4;nt1++)
	{
		for(int nt2=nt1;nt2<4;nt2++)
		{

			tdir=basepath+"/" + "/"+nts[nt1]+"_"+nts[nt2];
			
			ifstream inputFile(tdir+"/bp_length.txt");
			if(!inputFile.is_open())
			{
				cout<<"Bin info file can not be opened (Basepair parameters [length])\n"<<endl;
				return false;
			}

			while(getline(inputFile, oneline))
			{
				stringstream ssbase(oneline);
				count=0;
				while (getline(ssbase, field, ','))
				{
					if(count==0) mean=stof(field); //mean
					if(count==1) stdval=stof(field); //stdval
					count+=1;
				}

				stparameter[nts[nt1][0]][nts[nt2][0]][0][0] = mean;
				stparameter[nts[nt1][0]][nts[nt2][0]][0][1] = stdval;
				stparameter[nts[nt2][0]][nts[nt1][0]][0][0] = mean;
				stparameter[nts[nt2][0]][nts[nt1][0]][0][1] = stdval;
			}
		}
	}

	/*
	for(int nt1=0;nt1<4;nt1++)
	{
		for(int nt2=0;nt2<4;nt2++)
		{
			cout<<nts[nt1][0] << "-" << nts[nt2][0] << "-" << stparameter[nts[nt1][0]][nts[nt2][0]][0][0] << "-" << stparameter[nts[nt1][0]][nts[nt2][0]][0][1]<<endl;
		}
	}
	*/

	return true;
}

bool loadbphbondparams(string basepath, map<char, map<char, map<int, map<int, HydroBondParameters>>>> &stparameter)
{
	vector<string> nts = {"A","U","G","C"};
	string tdir,pdbpath;
	vector<string>refname;
	FILE *filein;
	string oneline_str;
	string	oneline;	
	string st;
	int count;
	double mean, stdval;
	string field;
	string acceptor_base, donor_base;
	string acceptor_atom, donor_atom;
	double adlength_mean, adlength_std;
	double dhlength_mean, dhlength_std;
	double ahlength_mean, ahlength_std;
	double adhangle_mean, adhangle_std;
	double dhaangle_mean, dhaangle_std;
	double hadangle_mean, hadangle_std;
	bool flagf;

	for(int nt1=0;nt1<4;nt1++)
	{
		for(int nt2=nt1;nt2<4;nt2++)
		{
			tdir=basepath+"/"+nts[nt1]+"_"+nts[nt2];
			
			ifstream inputFile(tdir+"/bp_hbond.txt");
			if(!inputFile.is_open())
			{
				cout<<"Bin info file can not be opened (Basepair hbond parameters)\n"<<endl;
				return false;
			}

			while(getline(inputFile, oneline))
			{
				stringstream ssbase(oneline);
				count=0;
				while (getline(ssbase, field, ','))
				{
					if(count==0) acceptor_base=field; 
					if(count==1) donor_base=field; 
					if(count==2) acceptor_atom=field; 
					if(count==3) donor_atom=field; 
					if(count==4) adlength_mean=stof(field); 
					if(count==5) adlength_std=stof(field); 
					if(count==6) dhlength_mean=stof(field);
					if(count==7) dhlength_std=stof(field);
					if(count==8) ahlength_mean=stof(field);
					if(count==9) ahlength_std=stof(field);
					if(count==10) adhangle_mean=stof(field);
					if(count==11) adhangle_std=stof(field);
					if(count==12) dhaangle_mean=stof(field);
					if(count==13) dhaangle_std=stof(field);
					if(count==14) hadangle_mean=stof(field);
					if(count==15) hadangle_std=stof(field);
					count+=1;
				}

				stparameter[acceptor_base[0]][donor_base[0]][atomid2mcidx[acceptor_base[0]][hydfullatomname2idx[acceptor_base][acceptor_atom]]][atomid2mcidx[donor_base[0]][hydfullatomname2idx[donor_base][donor_atom]]].adlength_mean = adlength_mean;
				stparameter[acceptor_base[0]][donor_base[0]][atomid2mcidx[acceptor_base[0]][hydfullatomname2idx[acceptor_base][acceptor_atom]]][atomid2mcidx[donor_base[0]][hydfullatomname2idx[donor_base][donor_atom]]].adlength_std = adlength_std;
				stparameter[acceptor_base[0]][donor_base[0]][atomid2mcidx[acceptor_base[0]][hydfullatomname2idx[acceptor_base][acceptor_atom]]][atomid2mcidx[donor_base[0]][hydfullatomname2idx[donor_base][donor_atom]]].dhlength_mean = dhlength_mean;
				stparameter[acceptor_base[0]][donor_base[0]][atomid2mcidx[acceptor_base[0]][hydfullatomname2idx[acceptor_base][acceptor_atom]]][atomid2mcidx[donor_base[0]][hydfullatomname2idx[donor_base][donor_atom]]].dhlength_std = dhlength_std;
				stparameter[acceptor_base[0]][donor_base[0]][atomid2mcidx[acceptor_base[0]][hydfullatomname2idx[acceptor_base][acceptor_atom]]][atomid2mcidx[donor_base[0]][hydfullatomname2idx[donor_base][donor_atom]]].ahlength_mean = ahlength_mean;
				stparameter[acceptor_base[0]][donor_base[0]][atomid2mcidx[acceptor_base[0]][hydfullatomname2idx[acceptor_base][acceptor_atom]]][atomid2mcidx[donor_base[0]][hydfullatomname2idx[donor_base][donor_atom]]].ahlength_std = ahlength_std;
				stparameter[acceptor_base[0]][donor_base[0]][atomid2mcidx[acceptor_base[0]][hydfullatomname2idx[acceptor_base][acceptor_atom]]][atomid2mcidx[donor_base[0]][hydfullatomname2idx[donor_base][donor_atom]]].adhangle_mean = adhangle_mean;
				stparameter[acceptor_base[0]][donor_base[0]][atomid2mcidx[acceptor_base[0]][hydfullatomname2idx[acceptor_base][acceptor_atom]]][atomid2mcidx[donor_base[0]][hydfullatomname2idx[donor_base][donor_atom]]].adhangle_std = adhangle_std;
				stparameter[acceptor_base[0]][donor_base[0]][atomid2mcidx[acceptor_base[0]][hydfullatomname2idx[acceptor_base][acceptor_atom]]][atomid2mcidx[donor_base[0]][hydfullatomname2idx[donor_base][donor_atom]]].dhaangle_mean = dhaangle_mean;
				stparameter[acceptor_base[0]][donor_base[0]][atomid2mcidx[acceptor_base[0]][hydfullatomname2idx[acceptor_base][acceptor_atom]]][atomid2mcidx[donor_base[0]][hydfullatomname2idx[donor_base][donor_atom]]].dhaangle_std = dhaangle_std;
				stparameter[acceptor_base[0]][donor_base[0]][atomid2mcidx[acceptor_base[0]][hydfullatomname2idx[acceptor_base][acceptor_atom]]][atomid2mcidx[donor_base[0]][hydfullatomname2idx[donor_base][donor_atom]]].hadangle_mean = hadangle_mean;
				stparameter[acceptor_base[0]][donor_base[0]][atomid2mcidx[acceptor_base[0]][hydfullatomname2idx[acceptor_base][acceptor_atom]]][atomid2mcidx[donor_base[0]][hydfullatomname2idx[donor_base][donor_atom]]].hadangle_std = hadangle_std;
			}
		}
	}

	/*
	for(int nt1=0;nt1<4;nt1++)
	{
		for(int nt2=0;nt2<4;nt2++)
		{
			cout<<nts[nt1][0] << "-" << nts[nt2][0] << "-" << stparameter[nts[nt1][0]][nts[nt2][0]][0][0] << "-" << stparameter[nts[nt1][0]][nts[nt2][0]][0][1]<<endl;
		}
	}
	*/

	return true;
}


bool loadbptemplates(string basepath, map<char, map<char, vector<vector<Residueatoms>>>> &stparameter){
	vector<string> nts = {"A","C","G","U"};
	string tdir,pdbpath;
	vector<string>refname;
	FILE *filein;
	string bptype;
	string oneline_str;
	char	oneline[600];	
	string st;
	bool flagf;

	map< char, int> Reschar2Idint={{'A',0},{'C',1},{'G',2},{'U',3},{'X',4}};

	static map<char, map<string, int>> atom2cid={
		{'A',{{" P  ", 0},{" O5'", 1},{" C5'", 2},{" C4'", 3},{" C3'", 4},{" O3'", 5},{" C2'", 6},{" O4'", 7},{" C1'", 8},{" O2'",9},{" OP1",10},{" OP2",11},{" N9 ",12},{" C4 ",13},{" C8 ",14},{" N1 ",15},{" C2 ",16},{" C6 ",17},{" N7 ",18},{" C5 ",19},{" N6 ",20},{" N3 ",21}}},
		{'C',{{" P  ", 0},{" O5'", 1},{" C5'", 2},{" C4'", 3},{" C3'", 4},{" O3'", 5},{" C2'", 6},{" O4'", 7},{" C1'", 8},{" O2'",9},{" OP1",10},{" OP2",11},{" N1 ",12},{" C2 ",13},{" C6 ",14},{" C4 ",15},{" C5 ",16},{" N3 ",17},{" O2 ",18},{" N4 ",19}}},
		{'G',{{" P  ", 0},{" O5'", 1},{" C5'", 2},{" C4'", 3},{" C3'", 4},{" O3'", 5},{" C2'", 6},{" O4'", 7},{" C1'", 8},{" O2'",9},{" OP1",10},{" OP2",11},{" N9 ",12},{" C4 ",13},{" C8 ",14},{" N1 ",15},{" C2 ",16},{" C6 ",17},{" N7 ",18},{" C5 ",19},{" N3 ",20},{" O6 ",21},{" N2 ",22}}},
		{'U',{{" P  ", 0},{" O5'", 1},{" C5'", 2},{" C4'", 3},{" C3'", 4},{" O3'", 5},{" C2'", 6},{" O4'", 7},{" C1'", 8},{" O2'",9},{" OP1",10},{" OP2",11},{" N1 ",12},{" C2 ",13},{" C6 ",14},{" C4 ",15},{" C5 ",16},{" N3 ",17},{" O2 ",18},{" O4 ",19}}}
	};

	for(int i=0;i<4;i++)
	{
		for(int k=0;k<4;k++)
		{
			bptype=nts[i]+"-"+nts[k];
			if(bptype != "A-U" && bptype != "G-C" && bptype != "G-U") continue;

			tdir=basepath+"/"+nts[i]+"-"+nts[k];

			refname.clear();
			for (const auto& entry : fs::directory_iterator(tdir)) {
				if (fs::is_regular_file(entry) && entry.path().extension() == ".pdb") {
					refname.push_back(entry.path());
				}
			}

			stparameter[nts[i][0]][nts[k][0]].resize(refname.size());
			stparameter[nts[k][0]][nts[i][0]].resize(refname.size());

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

							residue.ptv[atom2cid[residue.aaa][st.substr(12,4)]].x=float(atof(st.substr(30,8).c_str())); 
							residue.ptv[atom2cid[residue.aaa][st.substr(12,4)]].y=float(atof(st.substr(38,8).c_str())); 
							residue.ptv[atom2cid[residue.aaa][st.substr(12,4)]].z=float(atof(st.substr(46,8).c_str())); 
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

				residues.clear();
				invresidues.clear();
				fclose(filein);
			}
		}
	}

	return true;
}

bool loadsinglenttemplates(string basepath, map<char, vector<vector<Residueatoms>>> &stparameter){
	vector<string> nts = {"A","C","G","U"};
	string tdir,pdbpath;
	vector<string>refname;
	FILE *filein;
	string bptype;
	string oneline_str;
	char	oneline[600];	
	string st;
	bool flagf;

	map< char, int> Reschar2Idint={{'A',0},{'C',1},{'G',2},{'U',3},{'X',4}};

	static map<char, map<string, int>> atom2cid={
		{'A',{{" P  ", 0},{" O5'", 1},{" C5'", 2},{" C4'", 3},{" C3'", 4},{" O3'", 5},{" C2'", 6},{" O4'", 7},{" C1'", 8},{" O2'",9},{" OP1",10},{" OP2",11},{" N9 ",12},{" C4 ",13},{" C8 ",14},{" N1 ",15},{" C2 ",16},{" C6 ",17},{" N7 ",18},{" C5 ",19},{" N6 ",20},{" N3 ",21}}},
		{'C',{{" P  ", 0},{" O5'", 1},{" C5'", 2},{" C4'", 3},{" C3'", 4},{" O3'", 5},{" C2'", 6},{" O4'", 7},{" C1'", 8},{" O2'",9},{" OP1",10},{" OP2",11},{" N1 ",12},{" C2 ",13},{" C6 ",14},{" C4 ",15},{" C5 ",16},{" N3 ",17},{" O2 ",18},{" N4 ",19}}},
		{'G',{{" P  ", 0},{" O5'", 1},{" C5'", 2},{" C4'", 3},{" C3'", 4},{" O3'", 5},{" C2'", 6},{" O4'", 7},{" C1'", 8},{" O2'",9},{" OP1",10},{" OP2",11},{" N9 ",12},{" C4 ",13},{" C8 ",14},{" N1 ",15},{" C2 ",16},{" C6 ",17},{" N7 ",18},{" C5 ",19},{" N3 ",20},{" O6 ",21},{" N2 ",22}}},
		{'U',{{" P  ", 0},{" O5'", 1},{" C5'", 2},{" C4'", 3},{" C3'", 4},{" O3'", 5},{" C2'", 6},{" O4'", 7},{" C1'", 8},{" O2'",9},{" OP1",10},{" OP2",11},{" N1 ",12},{" C2 ",13},{" C6 ",14},{" C4 ",15},{" C5 ",16},{" N3 ",17},{" O2 ",18},{" O4 ",19}}}
	};

	for(int i=0;i<4;i++)
	{
		tdir=basepath+"/"+nts[i];

		refname.clear();
		for (const auto& entry : fs::directory_iterator(tdir)) {
			if (fs::is_regular_file(entry) && entry.path().extension() == ".pdb") {
				refname.push_back(entry.path());
			}
		}

		stparameter[nts[i][0]].resize(refname.size());

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

						residue.ptv[atom2cid[residue.aaa][st.substr(12,4)]].x=float(atof(st.substr(30,8).c_str())); 
						residue.ptv[atom2cid[residue.aaa][st.substr(12,4)]].y=float(atof(st.substr(38,8).c_str())); 
						residue.ptv[atom2cid[residue.aaa][st.substr(12,4)]].z=float(atof(st.substr(46,8).c_str())); 
					}
				}
				strcpy(oldline,oneline);
			}

			stparameter[nts[i][0]][j]=residues;

			residues.clear();
			fclose(filein);
		}

	}

	return true;
}

bool MCgetConnect(string basepath, map<string, vector<vector<int>>> &molConnectone, map<string, vector<vector<int>>> &molConnecttwo, map<string, vector<vector<int>>> &molConnectthree)
{
    vector<string> nts = {"A","C","G","U"};
    int count;
    string	oneline, field;
    string atom1, atom2, atom3, atom4;
    string base;
    int idx1, idx2, idx3, idx4;
    ifstream inputFile;

    for(int i=0;i<nts.size();i++)
	{
        base = nts[i];

        inputFile.open(basepath+"/"+base+"_one.connect");
        if(!inputFile)
		{
			cout<<"Parameter file can not be opened (Connect one-" << base << ")\n"<<endl;
			return false;
		}

        while(getline(inputFile, oneline))
		{
            stringstream ssbase(oneline);
			count=0;
			while (getline(ssbase, field, ','))
			{
                //if(count==0) ssidx=stoi(field); //ssidx
				if(count==0) atom1=field; 
				if(count==1) atom2=field;
                count+=1;
            }
            if(count == 2){
                auto it1 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom1);
                idx1 = -1;
                if(it1 != hyd_full_atoms_name[base].end()){
                    idx1=atomid2mcidx[base[0]][distance(hyd_full_atoms_name[base].begin(), it1)];
                }else{
                    cout<<atom1<<endl;
                    cout<<"Error1"<<endl;
                    exit(1);
                }

                auto it2 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom2);
                idx2 = -1;
                if(it2 != hyd_full_atoms_name[base].end()){
                    idx2=atomid2mcidx[base[0]][distance(hyd_full_atoms_name[base].begin(), it2)];
                }else{
                    cout<<atom2<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }
                if(idx1 >= 0 && idx2 >= 0) molConnectone[base].push_back({idx1, idx2});
            }
        }
        inputFile.close();
        inputFile.clear();

        inputFile.open(basepath+"/"+base+"_two.connect");
        if(!inputFile)
		{
			cout<<"Parameter file can not be opened (Connect two-" << base << ")\n"<<endl;
			return false;
		}

        while(getline(inputFile, oneline))
		{
            stringstream ssbase(oneline);
			count=0;
			while (getline(ssbase, field, ','))
			{
                //if(count==0) ssidx=stoi(field); //ssidx
				if(count==0) atom1=field; 
				if(count==1) atom2=field;
                if(count==2) atom3=field;
                count+=1;
            }
            if(count == 3){
                auto it1 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom1);
                idx1 = -1;
                if(it1 != hyd_full_atoms_name[base].end()){
                    idx1=atomid2mcidx[base[0]][distance(hyd_full_atoms_name[base].begin(), it1)];
                }else{
                    cout<<atom1<<endl;
                    cout<<"Error1"<<endl;
                    exit(1);
                }

                auto it2 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom2);
                idx2 = -1;
                if(it2 != hyd_full_atoms_name[base].end()){
                    idx2=atomid2mcidx[base[0]][distance(hyd_full_atoms_name[base].begin(), it2)];
                }else{
                    cout<<atom2<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                auto it3 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom3);
                idx3 = -1;
                if(it3 != hyd_full_atoms_name[base].end()){
                    idx3=atomid2mcidx[base[0]][distance(hyd_full_atoms_name[base].begin(), it3)];
                }else{
                    cout<<atom3<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                if(idx1 >= 0 && idx2 >= 0 && idx3 >= 0) molConnecttwo[base].push_back({idx1, idx2, idx3});
            }
        }
        inputFile.close();
        inputFile.clear();

        inputFile.open(basepath+"/"+base+"_three.connect");
        if(!inputFile)
		{
			cout<<"Parameter file can not be opened (Connect three-" << base << ")\n"<<endl;
			return false;
		}

        while(getline(inputFile, oneline))
		{
            stringstream ssbase(oneline);
			count=0;
			while (getline(ssbase, field, ','))
			{
                //if(count==0) ssidx=stoi(field); //ssidx
				if(count==0) atom1=field; 
				if(count==1) atom2=field;
                if(count==2) atom3=field;
                if(count==3) atom4=field;
                count+=1;
            }
            if(count == 4){
                auto it1 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom1);
                idx1 = -1;
                if(it1 != hyd_full_atoms_name[base].end()){
                    idx1=atomid2mcidx[base[0]][distance(hyd_full_atoms_name[base].begin(), it1)];
                }else{
                    cout<<atom1<<endl;
                    cout<<"Error1"<<endl;
                    exit(1);
                }

                auto it2 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom2);
                idx2 = -1;
                if(it2 != hyd_full_atoms_name[base].end()){
                    idx2=atomid2mcidx[base[0]][distance(hyd_full_atoms_name[base].begin(), it2)];
                }else{
                    cout<<atom2<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                auto it3 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom3);
                idx3 = -1;
                if(it3 != hyd_full_atoms_name[base].end()){
                    idx3=atomid2mcidx[base[0]][distance(hyd_full_atoms_name[base].begin(), it3)];
                }else{
                    cout<<atom3<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                auto it4 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom4);
                idx4 = -1;
                if(it4 != hyd_full_atoms_name[base].end()){
                    idx4=atomid2mcidx[base[0]][distance(hyd_full_atoms_name[base].begin(), it4)];
                }else{
                    cout<<atom4<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                if(idx1 >= 0 && idx2 >= 0 && idx3 >= 0 && idx4 >=0) molConnectthree[base].push_back({idx1, idx2, idx3, idx4});
            }
        }
        inputFile.close();
        inputFile.clear();
    }
    return true;
}

bool MCgetConnectAdj(string basepath, map<string, vector<vector<int>>> &molAdjConnectone, map<string, vector<vector<int>>> &molAdjConnecttwo, map<string, vector<vector<int>>> &molAdjConnectthree, map<string, vector<vector<int>>> &molAdjConnectoneidx, map<string, vector<vector<int>>> &molAdjConnecttwoidx, map<string, vector<vector<int>>> &molAdjConnectthreeidx)
{
    vector<string> nts = {"A","C","G","U"};
    int count;
    string	oneline, field;
    string atom1, atom2, atom3, atom4;
    string base;
    int idx1, idx2, idx3, idx4;
    int residx1, residx2, residx3, residx4;
    ifstream inputFile;

    for(int i=0;i<nts.size();i++)
	{
        base = nts[i];

        inputFile.open(basepath+"/"+base+"_one_adj.connect");
        if(!inputFile)
		{
			cout<<"Parameter file can not be opened (Adj connect one-" << base << ")\n"<<endl;
			return false;
		}

        while(getline(inputFile, oneline))
		{
            stringstream ssbase(oneline);
			count=0;
			while (getline(ssbase, field, ','))
			{
                //if(count==0) ssidx=stoi(field); //ssidx
				if(count==0) atom1=field; 
				if(count==1) atom2=field;
                if(count==2) residx1=stoi(field); 
				if(count==3) residx2=stoi(field);
                count+=1;
            }
            if(count == 4){
                auto it1 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom1);
                idx1 = -1;
                if(it1 != hyd_full_atoms_name[base].end()){
                    idx1=atomid2mcidx[base[0]][distance(hyd_full_atoms_name[base].begin(), it1)];
                }else{
                    cout<<atom1<<endl;
                    cout<<"Error1"<<endl;
                    exit(1);
                }

                auto it2 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom2);
                idx2 = -1;
                if(it2 != hyd_full_atoms_name[base].end()){
                    idx2=atomid2mcidx[base[0]][distance(hyd_full_atoms_name[base].begin(), it2)];
                }else{
                    cout<<atom2<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                molAdjConnectone[base].push_back({idx1, idx2});
                molAdjConnectoneidx[base].push_back({residx1, residx2});
            }
        }
        inputFile.close();
        inputFile.clear();

        inputFile.open(basepath+"/"+base+"_two_adj.connect");
        if(!inputFile)
		{
			cout<<"Parameter file can not be opened (Adj connect two-" << base << ")\n"<<endl;
			return false;
		}

        while(getline(inputFile, oneline))
		{
            stringstream ssbase(oneline);
			count=0;
			while (getline(ssbase, field, ','))
			{
                //if(count==0) ssidx=stoi(field); //ssidx
				if(count==0) atom1=field; 
				if(count==1) atom2=field;
                if(count==2) atom3=field;
                if(count==3) residx1=stoi(field); 
				if(count==4) residx2=stoi(field);
                if(count==5) residx3=stoi(field);
                count+=1;
            }
            if(count == 6){
                auto it1 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom1);
                idx1 = -1;
                if(it1 != hyd_full_atoms_name[base].end()){
                    idx1=atomid2mcidx[base[0]][distance(hyd_full_atoms_name[base].begin(), it1)];
                }else{
                    cout<<atom1<<endl;
                    cout<<"Error1"<<endl;
                    exit(1);
                }

                auto it2 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom2);
                idx2 = -1;
                if(it2 != hyd_full_atoms_name[base].end()){
                    idx2=atomid2mcidx[base[0]][distance(hyd_full_atoms_name[base].begin(), it2)];
                }else{
                    cout<<atom2<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                auto it3 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom3);
                idx3 = -1;
                if(it3 != hyd_full_atoms_name[base].end()){
                    idx3=atomid2mcidx[base[0]][distance(hyd_full_atoms_name[base].begin(), it3)];
                }else{
                    cout<<atom3<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                molAdjConnecttwo[base].push_back({idx1, idx2, idx3});
                molAdjConnecttwoidx[base].push_back({residx1, residx2, residx3});
            }
        }
        inputFile.close();
        inputFile.clear();

        inputFile.open(basepath+"/"+base+"_three_adj.connect");
        if(!inputFile)
		{
			cout<<"Parameter file can not be opened (Adj connect three-" << base << ")\n"<<endl;
			return false;
		}

        while(getline(inputFile, oneline))
		{
            stringstream ssbase(oneline);
			count=0;
			while (getline(ssbase, field, ','))
			{
                //if(count==0) ssidx=stoi(field); //ssidx
				if(count==0) atom1=field; 
				if(count==1) atom2=field;
                if(count==2) atom3=field;
                if(count==3) atom4=field;
                if(count==4) residx1=stoi(field); 
				if(count==5) residx2=stoi(field);
                if(count==6) residx3=stoi(field);
                if(count==7) residx4=stoi(field);
                count+=1;
            }
            if(count == 8){
                auto it1 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom1);
                idx1 = -1;
                if(it1 != hyd_full_atoms_name[base].end()){
                    idx1=atomid2mcidx[base[0]][distance(hyd_full_atoms_name[base].begin(), it1)];
                }else{
                    cout<<atom1<<endl;
                    cout<<"Error1"<<endl;
                    exit(1);
                }

                auto it2 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom2);
                idx2 = -1;
                if(it2 != hyd_full_atoms_name[base].end()){
                    idx2=atomid2mcidx[base[0]][distance(hyd_full_atoms_name[base].begin(), it2)];
                }else{
                    cout<<atom2<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                auto it3 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom3);
                idx3 = -1;
                if(it3 != hyd_full_atoms_name[base].end()){
                    idx3=atomid2mcidx[base[0]][distance(hyd_full_atoms_name[base].begin(), it3)];
                }else{
                    cout<<atom3<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                auto it4 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom4);
                idx4 = -1;
                if(it4 != hyd_full_atoms_name[base].end()){
                    idx4=atomid2mcidx[base[0]][distance(hyd_full_atoms_name[base].begin(), it4)];
                }else{
                    cout<<atom4<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                molAdjConnectthree[base].push_back({idx1, idx2, idx3, idx4});
                molAdjConnectthreeidx[base].push_back({residx1, residx2, residx3, residx4});
            }
        }
        inputFile.close();
        inputFile.clear();
    }
    return true;
}

bool MCloadfeatures(string basepath, map<string, map<int, vector<double>>> &atomfeatures)
{
    int count;
    string	oneline, field;
    string atom;
    double Charge, vdWR0, vdWEps;
    string base;
    int idx;
    double mean, stdval;
    vector<string> nts = {"A","C","G","U"};
    for(int i=0;i<nts.size();i++)
	{
        base = nts[i];
        ifstream inputFile(basepath+"/"+base+".fea");
        if(!inputFile.is_open())
		{
			cout<<"Parameter file can not be opened (feature-" << base << ")\n"<<endl;
			return false;
		}

        while(getline(inputFile, oneline))
		{
            stringstream ssbase(oneline);
			count=0;
			while (getline(ssbase, field, ','))
			{
                //if(count==0) ssidx=stoi(field); //ssidx
				if(count==0) atom=field; 
				if(count==1) Charge=stof(field); 
				if(count==2) vdWR0=stof(field); 
                if(count==3) vdWEps=stof(field); 
				count+=1;
            }

            if(count == 4){
                auto it = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom);
                idx = -1;
                if(it != hyd_full_atoms_name[base].end()){
                    idx=atomid2mcidx[base[0]][distance(hyd_full_atoms_name[base].begin(), it)];
                }else{
                    cout<<atom<<endl;
                    cout<<"Error1"<<endl;
                    exit(1);
                }
                atomfeatures[base][idx]={Charge, vdWR0, vdWEps};
            }
        }
    }

    return true;
}

void MCsetconnect(map<string, vector<vector<bool>>> &intraatomconnect, map<string, map<string, vector<vector<int>>>> &intramolconnect, map<string, vector<vector<int>>> &ntmolconnect){
	vector<string> nts={"A","C","G","U"};
	map<string, int> atomnumnts={{"A", 33}, {"C", 31}, {"G", 34}, {"U", 30}};
	string base, base1_str, base2_str;
	char base1, base2;
	int atomnum1, atomnum2, atomnum;

	for(int i=0;i<nts.size();i++){
        base=nts[i];
        atomnum = atomnumnts[base];
        intraatomconnect[base].resize(atomnum, vector<bool>(atomnum, false));
        for(int k=0;k<molConnectone.at(base).size();k++){
            intraatomconnect[base][molConnectone.at(base)[k][0]][molConnectone.at(base)[k][1]] = true;
            intraatomconnect[base][molConnectone.at(base)[k][1]][molConnectone.at(base)[k][0]] = true;
        }
    }

	for(int i=0;i<nts.size();i++){
		base1_str=nts[i];
		base1=base1_str[0];
		atomnum1=AccAtomParamNum[base1][4];
		for(int j=0;j<nts.size();j++){
			base2_str=nts[j];
			base2=base2_str[0];
			atomnum2=AccAtomParamNum[base2][4];

			intramolconnect[base1_str][base2_str].resize(atomnum1, vector<int>(atomnum2, 4));
			for(int k=0;k<molAdjConnectthree.at(base1_str).size();k++){
				intramolconnect[base1_str][base2_str][molAdjConnectthree.at(base1_str)[k][0]][molAdjConnectthree.at(base1_str)[k][3]] = 3;
			}
			for(int k=0;k<molAdjConnecttwo.at(base1_str).size();k++){
				intramolconnect[base1_str][base2_str][molAdjConnecttwo.at(base1_str)[k][0]][molAdjConnecttwo.at(base1_str)[k][2]] = 2;
			}
			for(int k=0;k<molAdjConnectone.at(base1_str).size();k++){
				intramolconnect[base1_str][base2_str][molAdjConnectone.at(base1_str)[k][0]][molAdjConnectone.at(base1_str)[k][1]] = 1;
			}
		}

		ntmolconnect[base1_str].resize(atomnum1, vector<int>(atomnum1, 4));
		for(int k=0;k<molConnectthree.at(base1_str).size();k++){
			ntmolconnect[base1_str][molConnectthree.at(base1_str)[k][0]][molConnectthree.at(base1_str)[k][3]] = 3;
			ntmolconnect[base1_str][molConnectthree.at(base1_str)[k][3]][molConnectthree.at(base1_str)[k][0]] = 3;
		}
		for(int k=0;k<molConnecttwo.at(base1_str).size();k++){
			ntmolconnect[base1_str][molConnecttwo.at(base1_str)[k][0]][molConnecttwo.at(base1_str)[k][2]] = 2;
			ntmolconnect[base1_str][molConnecttwo.at(base1_str)[k][2]][molConnecttwo.at(base1_str)[k][0]] = 2;
		}
		for(int k=0;k<molConnectone.at(base1_str).size();k++){
			ntmolconnect[base1_str][molConnectone.at(base1_str)[k][0]][molConnectone.at(base1_str)[k][1]] = 1;
			ntmolconnect[base1_str][molConnectone.at(base1_str)[k][1]][molConnectone.at(base1_str)[k][0]] = 1;
		}
		for(int k=0;k<atomnum1;k++){
			ntmolconnect[base1_str][k][k] = 0;
		}
	}

}

void MCcalc_interaction_distmap(const vector<Residueatoms> &rna_str, vector<vector<double>> &pmat){

    double dx, dy, dz, length;
	map<char, int> atomnumnts={{'A', 22}, {'C', 20}, {'G', 23}, {'U', 20}};
    pmat.resize(atomnumnts[rna_str[0].aaa]);
    for(int i=0;i<atomnumnts[rna_str[0].aaa];i++){
        pmat[i].resize(atomnumnts[rna_str[1].aaa]);
        for(int j=0;j<atomnumnts[rna_str[1].aaa];j++){
			dx = rna_str[0].ptv[i].x - rna_str[1].ptv[j].x;
			dy = rna_str[0].ptv[i].y - rna_str[1].ptv[j].y;
			dz = rna_str[0].ptv[i].z - rna_str[1].ptv[j].z;
            pmat[i][j] = sqrt(dx*dx + dy*dy + dz*dz+1e-7);
        }
    }
}

void MCget_distmap(const map<char, map<char, vector<vector<Residueatoms>>>> &rna_str, map<char, map<char, vector<vector<vector<double>>>>> &stackdistmap){
    for (const auto& pair1 : rna_str) {
		for(const auto& pair2 : pair1.second){
			stackdistmap[pair1.first][pair2.first].resize(pair2.second.size());
			for(int i=0;i<pair2.second.size();i++){
				MCcalc_interaction_distmap(pair2.second[i], stackdistmap[pair1.first][pair2.first][i]);
			}
		}        
    }
}

void LoadParams(string pramdir){
    loadgeometricalparameters(pramdir+"/geo", MCTORSIONPARAM, "torsion", 3.5);
    loadgeometricalparameters(pramdir+"/geo", MCINNERANGPARAM, "inner", 3.5);
    loadgeometricalparameters(pramdir+"/geo", MCLENGTHPARAM, "length", 3.5);
    
    loadgeometricalparameters(pramdir+"/geo_base", MCBASEINNERANGPARAM, "inner", 3.5);
    loadgeometricalparameters(pramdir+"/geo_base", MCBASELENGTHPARAM, "length", 3.5);
    loadgeometricalparameters(pramdir+"/geo_base", MCBASETORSIONPARAM, "torsion", 3.5);
    
    loadgeometricalparameters(pramdir+"/geo_add", MCADDINNERANGPARAM, "inner", 3.5);
    loadgeometricalparameters(pramdir+"/geo_add", MCADDLENGTHPARAM, "length", 3.5);

    loadgeometricalparameters(pramdir+"/geo_first", MCOPTORSIONPARAM, "torsion", 3.5);
    loadgeometricalparameters(pramdir+"/geo_first", MCOPINNERANGPARAM, "inner", 3.5);
    loadgeometricalparameters(pramdir+"/geo_first", MCOPLENGTHPARAM, "length", 3.5);
    
    loadgeometricalparameters(pramdir+"/geo_hyd", MCHYDINNERANGPARAM, "inner", 3.5);
    loadgeometricalparameters(pramdir+"/geo_hyd", MCHYDLENGTHPARAM, "length", 3.5);
    loadgeometricalparameters(pramdir+"/geo_hyd", MCHYDTORSIONPARAM, "torsion", 3.5);
    
    loadgeometricalparameters(pramdir+"/geo_hyd_add", MCHYDADDLENGTHPARAM, "length", 3.5);
	loadgeometricalparameters(pramdir+"/geo_hyd_add", MCHYDADDANGLEPARAM, "inner", 3.5);

	loaddistributionparam(pramdir+"/geo", MCTORSIONDISTRIBUTION, "torsion");
	loaddistributionparam(pramdir+"/geo_base", MCBASETORSIONDISTRIBUTION, "torsion");
	loaddistributionparam(pramdir+"/geo_first", MCOPTORSIONDISTRIBUTION, "torsion");
	loaddistributionparam(pramdir+"/geo_hyd", MCHYDTORSIONDISTRIBUTION, "torsion");

	loadramachandranplot(pramdir+"/tor_comb", TORSIONPAIRSPARAM);

	loadpairlengadj(pramdir+"/pairleng", PAIRLENGADJPARAM);
	loadpairlengnt(pramdir+"/pairleng", ALLPAIRLENGNTPARAM);

    loadvwd(pramdir+"/clash", MCVWPARAM);

	loadstackparams(pramdir+"/stacking", STACKPARAM);
	loadbpparams(pramdir+"/basepair", BPPARAM);
	loadbphbondparams(pramdir+"/basepair_hbond", BPHBOND);

	loadstacktemplates(pramdir+"/stacking_template", STACKTEMPLATE);
	loadbptemplates(pramdir+"/bp_template", BPTEMPLATE);
	loadsinglenttemplates(pramdir+"/loop_template", LOOPNTTEMPLATE);
	loadsinglenttemplates(pramdir+"/single_nt_template", NTTEMPLATE);

	MCgetConnect(pramdir+"/connect", molConnectone, molConnecttwo, molConnectthree);
	MCgetConnectAdj(pramdir+"/connect", molAdjConnectone, molAdjConnecttwo, molAdjConnectthree, molAdjConnectoneidx, molAdjConnecttwoidx, molAdjConnectthreeidx);
	MCloadfeatures(pramdir+"/features", AtomFeatures);

	MCsetconnect(intraAtomConnect, intramolConnect, ntmolConnect);

	MCget_distmap(STACKTEMPLATE[0], STACKDISTMAP[0]);
	MCget_distmap(STACKTEMPLATE[1], STACKDISTMAP[1]);
	MCget_distmap(BPTEMPLATE, BPDISTMAP);
}














