#include <iostream>
#include <vector>
#include <filesystem>
#include <cstdlib>
#include <unordered_map>
#include <fstream>
#include <sstream>

//#include "sspred.hpp"
#include "RNArefine.hpp"
//#include <experimental/filesystem>

using namespace std;
namespace fs = std::filesystem;

//vector<string> gettemplatefiles(string tdir)
void gettemplatefiles(string tdir, vector<string> &templatestrs)
{
	//vector<string>refname;
	try {
        for (const auto& entry : fs::directory_iterator(tdir)) {
            if (fs::is_regular_file(entry) && entry.path().extension() == ".pdb") {
				templatestrs.push_back(entry.path());
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Error:: template file cannot be loaded." << e.what() << std::endl;
    }
}

void gettemplateweights(vector<string> templatestrs, string weightfile, vector<double> &weightslist){
	
	int i;
	map<string, double> weights;
	ifstream inputFile(weightfile);
	weightslist.resize(templatestrs.size());
	string oneline, field;
	string filename;
	int count;
	double n;

	if(!inputFile.is_open())
	{
		cout<<"Warning: The template weight file specified by -libcwg does not exist (" + weightfile + ")\n"<<endl;
		cout << "Default value (1/<template number>) is set.\n" <<endl;

		n=1/templatestrs.size();
		for(i=0; i<templatestrs.size();i++){
			weightslist[i]=n;
		}
		return;
	}
	while(getline(inputFile, oneline))
	{
		stringstream ssbase(oneline);
		count=0;
		while (getline(ssbase, field, ','))
		{
			if(count==0) filename=field; //paramtype
			if(count==1) n=stof(field); //minval
			count+=1;
		}

		if(count < 2){
			cout<<"Warning: The format of the template weight file specified by -libcwg is wrong.\n"<<endl;
			cout << "Default value (1/<template number>) is set.\n" <<endl;
			n=1/templatestrs.size();
			for(i=0; i<templatestrs.size();i++){
				weightslist[i]=n;
			}
			return;
		}else{
			weights[filename]=n;
		}
	}

	double totalweight=0.0;
	for(i=0; i<templatestrs.size();i++){
		if (weights.find(templatestrs[i]) != weights.end()) {
			weightslist[i]=weights[templatestrs[i]];
			totalweight+=weightslist[i];
		}else{
			cout<<"Warning: The weight of "<< templatestrs[i] <<" is not found in the template weight file pecified by -libcwg.\n"<<endl;
			cout << "Default value (1/<template number>) is set.\n" <<endl;
			n=1/templatestrs.size();
			for(i=0; i<templatestrs.size();i++){
				weightslist[i]=n;
			}
			return;
		}
	}

	for(i=0; i<templatestrs.size();i++){
		weightslist[i]=weightslist[i]/totalweight;
	}
}


void parseArgs(int argc, char* argv[], unordered_map<string, string> &args) {

	for(int i = 1; i < argc; ++i) {
        string arg = argv[i];
		if(arg == "-h"){
			//cout<<"Usage: -i [path of input PDB (required)] -o [directory to output (required)] -libdir [directory to template structures] -libwg [weights of template structures] -bp [path of base pairing information (required)] -st [path of stacking information (required)] -r [refinment type (default 0)] -c1 [cutoff time 1 (default 10)] -c2 [cutoff time 2 (default 5)]\n";
			cout<< "./RNArefiner -i [path of input PDB (required)] -o [directory to output (required)] -bp [path of base pairing information (required)] -st [path of stacking information (required)] [other options]"<<endl;
			cout<<""<<endl;
			cout<<"[other options]"<<endl;
			cout<<"-libdir: Directory of additional template structures (All PDB files in the directory are loaded; default: None)"<<endl;
			cout<<"-libwg: Weights of an additional template structures-based energy  (-libdir must be specified; default: 10 if -libdir is specified, otherwise None)"<<endl;
			cout<<"-libcwg: Weights of each template structure (-libdir must be specified; default: 1/<template number> if -libdir is specified, otherwise None)"<<endl;
			cout<<"-time: Simulation time (hours; default: 10)"<<endl;
			//cout<<"-r: refinment type"<<endl;
			//cout<<"<refinment type>"<<endl;
			//cout<<"0: Replica-exchange Monte Carlo Simulation + Energy minimization (all)"<<endl;
			//cout<<"1: Replica-exchange Monte Carlo Simulation + Energy minimization (the minimum)"<<endl;
			//cout<<"2: Replica-exchange Monte Carlo Simulation"<<endl;
			//cout<<"3: Energy minimization"<<endl;
			//cout<<"<cutoff time (hours)>"<<endl;
			//cout<<"cutoff time 1: Cutoff time of replica-exchange Monte Carlo Simulation"<<endl;
			//cout<<"cutoff time 2: Cutoff time of Energy minimization"<<endl;
			exit(1);
		}
		if(arg == "-i" || arg == "-o" || arg == "-libdir" || arg == "-libwg" || arg == "-libcwg" || arg == "-bp" || arg == "-st" || arg == "-time") {
        //if(arg == "-i" || arg == "-o" || arg == "-libdir" || arg == "-libwg" || arg == "-bp" || arg == "-st" || arg == "-ct1" || arg == "-ct2" || arg == "-r" || arg == "-c1" || arg == "-c2") {
            if (i + 1 < argc) {
				if(string(argv[i + 1]) == "-i" || string(argv[i + 1]) == "-o" || string(argv[i + 1]) == "-libdir" || string(argv[i + 1]) == "-libwg" || string(argv[i + 1]) == "-libcwg" || string(argv[i + 1]) == "-bp" || string(argv[i + 1]) == "-st" || string(argv[i + 1]) == "-time"){
					cout<< "Error:" << arg << " is not specified"<<endl;
					exit(1);
				}else{
                	args[arg] = argv[i + 1];
               		++i;
				}
            } else {
                cout<< "Error:" << arg << " is not specified"<<endl;
				exit(1);
            }
        }else{
			cout<< "Warning: No argument of " << arg <<endl;
		}
    }

	if(args.find("-i") == args.end() || args.find("-bp") == args.end() || args.find("-st") == args.end() || args.find("-o") == args.end()) {
        cerr << "Error: -i, -bp, -st and -o arguments are required.\n" <<endl;
		cerr << "Check -h.\n" << endl;
        exit(1);
    }

	if(args.find("-libdir") == args.end() & args.find("-libcwg") != args.end()){
		cerr << "-libdir must be specified when -libcwg is specified.\n" <<endl;
		exit(1);
	}

	if(args.find("-libdir") == args.end() & args.find("-libwg") != args.end()){
		cerr << "-libdir must be specified when -libwg is specified.\n" <<endl;
		exit(1);
	}

	if(args.find("-r") == args.end()) args["-r"]="3";
	if(args.find("-time") == args.end()) args["-time"]="10";
	//if(args.find("-c1") == args.end()) args["-c1"]="10";
	//if(args.find("-c2") == args.end()) args["-c2"]="5";
	if(args.find("-libdir") == args.end()){
		args["-libdir"]=args["-o"]+"/template";
		fs::create_directories(args["-libdir"]);
		args["-libwg"]="";
		args["-libcwg"]="";
	}else{
		if(args.find("-libwg") == args.end()){
			args["-libwg"]="10";
		}
	}

	args["-r"]="3";
	args["-c1"]="0";
	args["-c2"]=args["-time"];
}

int main(int argc, char** argv)
{
	//TestingFunc tf;
	unordered_map<string, string> args;
	parseArgs(argc, argv, args);
	
	string infile=args["-i"];
	string bppath=args["-bp"];
	string stpath=args["-st"];
	string outdir=args["-o"];
	string templatedir=args["-libdir"];

	int reftype;
	try {
		reftype = stoi(args["-r"]);
	} catch(const invalid_argument& e) {
		cout << "Warning: Refinment type must be integer. Default value (3) is set.\n"<<endl;
		reftype = 3;
	}

	double cuttime1;
	try {
		cuttime1 = stof(args["-c1"]);
	} catch(const invalid_argument& e) {
		cout << "Warning: Simulation time (-time) must be double. Default value (0) is set.\n"<<endl;
		cuttime1 = 0;
	}

	double cuttime2;
	try {
		cuttime2 = stof(args["-c2"]);
	} catch(const invalid_argument& e) {
		cout << "Warning: Simulation time (-time) must be double. Default value (10) is set.\n" <<endl;
		cuttime2 = 10;
	}

	//fs::create_directories(templatedir);
	//fs::create_directories(outdir+"/image");

	//const char* envVarName = "DSSR_PATH";
	//string sssubdir=std::getenv(envVarName);
	//string pythonpath="python";

	string filePath=fs::canonical(fs::path(argv[0])).string();
	string basePath=fs::path(filePath).parent_path().string();
	string paramdir=fs::canonical(fs::path(basePath+"/../param")).string();
	//string mcannotate=basePath+"/../downloads/MC-Annotate_linux";
	//string stackprocpath=basePath+"/../stackingproc";

	vector<string>templatefiles;
	vector<double>templateweights;
	double templatecoeff;
	gettemplatefiles(templatedir, templatefiles);
	if(templatefiles.size() != 0){
		try {
			templatecoeff=stof(args["-libwg"]);
		}catch(const invalid_argument& e) {
			cout << "Warning: Weight for the energy of template structures (-libwg) must be double. Default value (10) is set.\n" <<endl;
			templatecoeff=10;
		}
		gettemplateweights(templatefiles, args["-libcwg"], templateweights);
	}else{
		templatecoeff=0;
	}

	fs::path currentPath = fs::current_path();
	string exepath = currentPath.string();

	RNArefineFunc rf;
	rf.rnarefine(exepath,infile,templatefiles,templatecoeff,templateweights,bppath,stpath,outdir,paramdir,reftype,cuttime1,cuttime2);
	
	return 1;
}
