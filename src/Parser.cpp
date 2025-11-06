
#include "Parser.hpp"

const string docstring=""
"[[RNAref3D]]\n"
"./RNAref3D global -i [path of input PDB (required)] -o [output directory (required)] -bp [path of base pairing information (required)] -st [path of stacking information (required)] [other options]\n"
"./RNAref3D local -i [path of input PDB (required)] -o [output directory (required)] -strlib [Directory of template and fragment structures] -distlib [Directory of distance maps] -contlib [Directory of contact maps] -bplib [Directory of base pairing information] -stacklib [Directory of stacking information]\n"
"[other options]\n"
"-m: Refinement mode (default: 3)\n"
"   1: Full atomic structure generation only\n"
"   2: MC-based refinement\n"
"   3: L-BFGS-based refinement\n"
"   4: MC-based refinement + L-BFGS-based refinement\n"
"-mtime: MC-based refinement time [seconds] (default: 7200)\n"
"-ltime: L-BFGS-based refinement time [seconds] (default: 3600)\n"
"-verbose: Enables detailed logging. Log is written to output/directory/logs/output.log.\n";

void parseArgs(int argc, char* argv[], unordered_map<string, string> &args) {

    args["-verbose"]="f";

    //string refinement_type=argv[1];

    //if(refinement_type == "local"){
        //args["-refinement_type"]=refinement_type;

    for(int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if(arg == "-h"){
            cout<<docstring<<endl;
            exit(1);
        }

        if(arg == "-verbose") args["-verbose"]="t";
        else if(arg == "-i" || arg == "-o" || arg == "-bp" || arg == "-st" || arg == "-m" || arg == "-mtime" || arg == "-ltime"){
            if (i + 1 < argc) {
                if(string(argv[i + 1]) == "-i" || string(argv[i + 1]) == "-o" || string(argv[i + 1]) == "-bp" || string(argv[i + 1]) == "-st" || string(argv[i + 1]) == "-m" || string(argv[i + 1]) == "-mtime" || string(argv[i + 1]) == "-ltime"){
                    cerr<< "Error:" << arg << " is not specified"<<endl;
                    exit(1);
                }else{
                    args[arg] = argv[i + 1];
                    ++i;

                    if(arg == "-i"){
                        if (!fs::exists(args["-i"])) {
                            cerr<<"Error: No such file: " << args["-i"] <<endl;
                            exit(1);
                        }
                    }else if(arg == "-bp"){
                        if (!fs::exists(args["-bp"])) {
                            cerr<<"Error: No such file: " << args["-bp"] <<endl;
                            exit(1);
                        }
                    }else if(arg == "-st"){
                        if (!fs::exists(args["-st"])) {
                            cerr<<"Error: No such file: " << args["-st"] <<endl;
                            exit(1);
                        }
                    }else if(arg == "-o"){
                        if (!fs::exists(args["-o"]) ||  !fs::is_directory(args["-o"])) {
                            cerr<<"Warning: No output directory: " << args["-o"] <<endl;
                            fs::create_directories(args["-o"]);
                            cerr<<"Info: Create output directory"<<endl;
                        }
                    }
                }
            } else {
                cerr<< "Error:" << arg << " is not specified"<<endl;
                exit(1);
            }
        }else{
            cerr<< "Warning: No argument of " << arg <<endl;
        }
    }

    if(args.find("-i") == args.end() || args.find("-bp") == args.end() || args.find("-st") == args.end() || args.find("-o") == args.end()) {
        cerr << "Error: -i, -bp, -st and -o arguments are required.\n" <<endl;
        cerr << "Check -h.\n" << endl;
        exit(1);
    }

    if(args.find("-m") == args.end()){
        args["-m"]="3";
    }else{
        try {
            size_t opt;
            stoi(args["-m"], &opt);
            if(opt < 0 || opt > 3){
                cerr << "Error: Specify a integer value (0-3) for the refinement mode\n" <<endl;
                exit(1);
            }
        }catch (...) {
            cerr << "Error: Specify a integer value (0-3) for the refinement mode\n" <<endl;
            exit(1);
        }
    }

    if(args.find("-mtime") == args.end())
    {
        args["-mtime"]="7200";
    }else{
        try {
            size_t opt;
            stod(args["-mtime"], &opt);
        }catch (...) {
            cerr << "MC-based refinement time must be a double or integer value\n" <<endl;
            exit(1);
        }
    }

    if(args.find("-ltime") == args.end())
    {
        args["-ltime"]="3600";
    }else{
        try {
            size_t opt;
            stod(args["-ltime"], &opt);
        }catch (...) {
            cerr << "L-BFGS-based refinement time must be a double or integer value\n" <<endl;
            exit(1);
        }
    }

    //}
    
    /*
    else if(refinement_type == "global"){
        args["-refinement_type"]=refinement_type;

        for(int i = 2; i < argc; ++i) {
            string arg = argv[i];
            if(arg == "-h"){
                cout<<docstring<<endl;
                exit(1);
            }

            if(arg == "-verbose") args["-verbose"]="t";
            else if(arg == "-i" || arg == "-o" || arg == "-bplib" || arg == "-stacklib" || arg == "-strlib" || arg == "-distlib" || arg == "-contlib") {
                if (i + 1 < argc) {
                    if(string(argv[i + 1]) == "-i" || string(argv[i + 1]) == "-o" || string(argv[i + 1]) == "-strlib" || string(argv[i + 1]) == "-distlib" || string(argv[i + 1]) == "-contlib" || string(argv[i + 1]) == "-bplib" || string(argv[i + 1]) == "-stacklib"){
                        cerr<< "Error: No option name of " << arg <<endl;
                        exit(1);
                    }else{
                        args[arg] = argv[i + 1];
                        ++i;

                        if(arg == "-i"){
                            if (!fs::exists(args["-i"])) {
                                cerr<<"Error: No such file: " << args["-i"] <<endl;
                                exit(1);
                            }
                        }else if(arg == "-bplib"){
                            if (!fs::exists(args["-bplib"]) ||  !fs::is_directory(args["-bplib"])) {
                                cerr<<"Error: No such file: " << args["-bplib"] <<endl;
                                exit(1);
                            }
                        }else if(arg == "-stacklib"){
                            if (!fs::exists(args["-stacklib"]) ||  !fs::is_directory(args["-stacklib"])) {
                                cerr<<"Error: No such file: " << args["-stacklib"] <<endl;
                                exit(1);
                            }
                        }else if(arg == "-strlib"){
                            if (!fs::exists(args["-strlib"]) ||  !fs::is_directory(args["-strlib"])) {
                                cerr<<"Error: No such file: " << args["-strlib"] <<endl;
                                exit(1);
                            }
                        }else if(arg == "-distlib"){
                            if (!fs::exists(args["-distlib"]) ||  !fs::is_directory(args["-distlib"])) {
                                cerr<<"Error: No such file: " << args["-distlib"] <<endl;
                                exit(1);
                            }
                        }else if(arg == "-contlib"){
                            if (!fs::exists(args["-contlib"]) ||  !fs::is_directory(args["-contlib"])) {
                                cerr<<"Error: No such file: " << args["-contlib"] <<endl;
                                exit(1);
                            }
                        }else if(arg == "-o"){
                            if (!fs::exists(args["-o"]) ||  !fs::is_directory(args["-o"])) {
                                cerr<<"Warning: No output directory: " << args["-o"] <<endl;
                                fs::create_directories(args["-o"]);
                                cout<<"Info: Create output directory"<<endl;
                            }
                        }
                    }
                }else {
                    cerr<< "Error:" << arg << " is not specified"<<endl;
                    exit(1);
                }
            }
        }

        if(args.find("-i") == args.end() || args.find("-o") == args.end()) {
            cerr << "Error: -i and -o arguments are required.\n" <<endl;
            cerr << "Check -h.\n" << endl;
            exit(1);
        }

        if(args.find("-strlib") == args.end() && args.find("-distlib") == args.end() && args.find("-contlib") == args.end() && args.find("-bplib") == args.end() && args.find("-stacklib") == args.end()){
            cerr << "Error: One of the option from -strlib, -distlib, -contlib, -bplib, or -stacklib must be specified.\n" <<endl;
            cerr << "Check -h.\n" << endl;
            exit(1);
        }

        if(args.find("-strlib") == args.end()){
            args["-strlib"]="";
        }

        if(args.find("-distlib") == args.end()){
            args["-distlib"]="";
        }

        if(args.find("-contlib") == args.end()){
            args["-contlib"]="";
        }

        if(args.find("-bplib") == args.end()){
            args["-bplib"]="";
        }

        if(args.find("-stacklib") == args.end()){
            args["-stacklib"]="";
        }
        
    }else{
        cerr << "Refinement type is not correct. Specify local or global.";
        exit(1);
    }
    */
}






