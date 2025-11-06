#include <iostream>
#include <filesystem>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include <sstream>
#include "MissingRNAatom.hpp"
#include "GetInteraction.hpp"
#include "PDBParser.hpp"
#include "RefineRNA.hpp"
#include "FragmentRefine.hpp"
#include "ResuperimpositionRNAatom.hpp"
#include "LoadGeoParams.hpp"
#include "CommonParams.hpp"
#include "StructureOptimization.hpp"
#include "MonteCarloSimulation.hpp"
#include "PreEnergyRefinement.hpp"
#include "Logger.hpp"
#include "Parser.hpp"
#include "LoadInputData.hpp"

using namespace std;
namespace fs = std::filesystem;

bool makedir(string path){
    fs::path dirPath = path;

    // ディレクトリを作成
    if (fs::create_directory(dirPath)) {
        return true;
    } else {
        return false;
    }

    return false;
}


int main(int argc,char **argv)
{
    unordered_map<string, string> args;
    parseArgs(argc, argv, args);

    //string refinetype=args["-refinement_type"];

    string outdir=args["-o"];
    Logger baseiclog = Logger();
    baseiclog.setPrint(true);
    bool logflag=false;
    if(args["-verbose"] == "t"){
        fs::create_directories(outdir+"/logs");
        logflag=true;
        baseiclog.setVerbose(true);
        baseiclog.setLogFile(outdir+"/logs/basic.log");
        baseiclog.log("Loggin:: On\n");
    }

    string infile=args["-i"];
    int atomic_detail=2;
    int allowX=0;
    map<int, map<int, int>> pdbidx2seqidx;
    map<string, int> chainstr2chainidx;
    
    ModelUnit pdb_entry=read_pdb_structure(infile.c_str(), chainstr2chainidx, pdbidx2seqidx,atomic_detail,allowX);

    check_targetpdb(pdb_entry);
    string seq=getseq(pdb_entry);

    string execfilePath=fs::canonical(fs::path(argv[0])).string();
	string basePath=fs::path(execfilePath).parent_path().string();
	string paramdir=fs::canonical(fs::path(basePath+"/../init_param")).string();

    /*
    if(refinetype == "global"){
        string strlib_dir=args["-strlib"];
        string distlib_dir=args["-distlib"];
        string contlib_dir=args["-contlib"];
        string bp_dir=args["-bplib"];
        string stack_dir=args["-stacklib"];

        cout<<"Loading parameters..."<<endl;
        vector<FragmentStructure> templatestrall;
        vector<vector<vector<bool>>> templatestrallflag;
        get_multitemplate_info(pdb_entry, strlib_dir, templatestrall, templatestrallflag, atomic_detail, allowX);

        //vector<pair<double,vector<size_t>>> bps_list;
        //vector<vector<double>> bps_scoremat;
        //get_multibp_info(pdb_entry, bp_dir, seq, bps_list, bps_scoremat, 0.25);

        //vector<pair<double,vector<size_t>>> stack_list;
        //vector<vector<double>> stack_scoremat;
        //get_multistack_info(pdb_entry, stack_dir, seq, stack_list, stack_scoremat, 0.25);

        vector<vector<pair<double,vector<size_t>>>> bps_list;
        vector<vector<pair<double,vector<size_t>>>> noncbps_list;
        vector<vector<vector<double>>> bps_scoremat;
        vector<vector<vector<bool>>> bps_flagmat;
        vector<vector<int>> bps_ranges;
        get_multibp_info(pdb_entry, bp_dir, seq, chainstr2chainidx, pdbidx2seqidx, bps_list, noncbps_list, bps_scoremat, bps_flagmat, bps_ranges);

        vector<vector<pair<double,vector<size_t>>>> stack_list;
        vector<vector<vector<double>>> stack_scoremat;
        vector<vector<vector<bool>>> stack_flagmat;
        vector<vector<int>> stack_ranges;
        get_multistack_info(pdb_entry, stack_dir, seq, chainstr2chainidx, pdbidx2seqidx, stack_list, stack_scoremat, stack_flagmat, stack_ranges);

        vector<vector<vector<double>>> distmaps;
        vector<vector<vector<bool>>> distmapflags;
        vector<vector<int>> distmapranges;
        getmultidistmap(pdb_entry, distlib_dir, seq, distmaps, distmapflags, distmapranges);

        vector<vector<vector<double>>> contactmaps;
        vector<vector<vector<bool>>> contactmapflags;
        vector<vector<int>> contactmapranges;
        getmulticontactmap(pdb_entry, contlib_dir, seq, contactmaps, contactmapflags, contactmapranges);


    }
    */    
    //else if(refinetype == "local"){
    string bpfile=args["-bp"];
    string stfile=args["-st"];

    double tolerance = 0.025;
    int option = stoi(args["-m"]); 
    double m_maxtime = stod(args["-mtime"]);
    double e_maxtime = stod(args["-ltime"]);

    vector<FragmentStructure> templatestrall;
    vector<vector<vector<bool>>> templatestrallflag;

    map<int, int> bp_template_idx, nonbp_template_idx;
    vector<vector<size_t> >res_str_vec;
    vector<pair<double,vector<size_t>>> bp_vec, bp_vec_unique;
    vector<pair<double,vector<size_t>>> noncbp_vec;
    vector<pair<double,vector<size_t>>> stack_vec;
    vector<vector<pair<int, double>>> stackmap;
    double energy;
    int spnum;
    size_t c, r, bp;
    int r1,r2,c1,c2;
    vector<vector<vector<int>>> stembps;
    map<int, vector<pair<int, int>>> loopnts;
    map<int, vector<vector<vector<int>>>> connectnts;
    vector<vector<int>> stemconnectnts;
    double final_energy;
    double alpha=0.4; 
    double beta=0.4;

    int stackseqst, stackseqend;
    int bpseqst, bpseqend;
    
    bool bpflag = loadBP(bpfile, seq, chainstr2chainidx, pdbidx2seqidx, bp_vec, noncbp_vec, bpseqst, bpseqend);
    if(!bpflag) exit(1);
        
    bool stflag = loadST(stfile, seq, chainstr2chainidx, pdbidx2seqidx, stack_vec, stackseqst, stackseqend);
    if(!stflag) exit(1);
        
    bp_vec_unique = bp_vec;
    filter_bp(bp_vec_unique);
    
    getsselement(pdb_entry, bp_vec_unique, stembps, loopnts, stemconnectnts, connectnts);
    getstmap(stack_vec, stackmap, pdb_entry.chains[0].residues.size());

    map<string, map<string, vector<ModelUnit>>> ideal_rna=parse_ideal_rna(paramdir);

    map<string, map<string, vector<vector<vector<double>>>>> stackdistmap;
    get_distmap(ideal_rna["stacking_adj"], stackdistmap["stacking_adj"]);
    get_distmap(ideal_rna["stacking_nonadj"], stackdistmap["stacking_nonadj"]);

    map<string, vector<vector<vector<double>>>> bpdistmap;
    get_distmap(ideal_rna["bp"], bpdistmap);

    vector<double> score_vec;
    get_score_map(seq, bp_vec, score_vec);

    ExtParams extparameters;
    FragmentStructure inputtedstrall;
    vector<vector<bool>> inputtedstrallflag;

    loadparameters(paramdir, extparameters);

    baseiclog.log("Full-atomic structure generation...\n");
    auto start_pre = std::chrono::high_resolution_clock::now();
    MissingRNAatom(pdb_entry, ideal_rna, stembps, loopnts, connectnts, bp_template_idx, nonbp_template_idx, inputtedstrall, inputtedstrallflag, extparameters, outdir, logflag, option);
    auto end_pre = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_pre = end_pre - start_pre;
    baseiclog.log("Processing time: " + to_string(elapsed_pre.count()) + " seconds\n");

    standardize_pdb_order(pdb_entry);

    string remark="";
    write_pdb_structure((outdir+"/superimpose.pdb").c_str(), pdb_entry, remark);

    //baseiclog.log("Saving to " +outdir+"/full_atomic.pdb\n\n");
    //write_pdb_structure((outdir+"/full_atomic.pdb").c_str(),pdb_entry);
    
    /*
    baseiclog.log("Pre-energy-based refinement...\n");
    auto start_pene = std::chrono::high_resolution_clock::now();
    PreEnergyRefine(pdb_entry, ideal_rna, stackmap, stackdistmap, stembps, loopnts, connectnts, extparameters);
    auto end_pene = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_pene = end_pene - start_pene;
    baseiclog.log("Processing time: " + to_string(elapsed_pene.count()) + " seconds\n");
    //baseiclog.log("Saving to " +outdir+"/pre_refine.pdb\n\n");
    //write_pdb_structure((outdir+"/pre_refine.pdb").c_str(),pdb_entry);
    */

    int iter=1;
    if(option>=1){
        size_t c,r,a;
        size_t aStart=0;
        for (c=0; c<pdb_entry.chains.size(); c++)
            for (r=0;r<pdb_entry.chains[c].residues.size();r++)
                for (a=aStart;a<pdb_entry.chains[c].residues[r].atoms.size();a++)
                    pdb_entry.chains[c].residues[r].atoms[a].movable=1;
        
        vector<vector<vector<double>>> hydrogencoords;
        if(option == 1 || option == 3){    
            baseiclog.log("REMC-based refinement...\n");
            baseiclog.log("MC refinement...\n");
            //double m_maxtime = 7200;
            auto start_remc = std::chrono::high_resolution_clock::now();
            final_energy=MonteCarloSimulationProcess(pdb_entry, paramdir+"/MC", bp_vec, loopnts, noncbp_vec, stemconnectnts, connectnts, stackmap, inputtedstrall, inputtedstrallflag, score_vec, templatestrall, templatestrallflag, extparameters, hydrogencoords, m_maxtime, outdir, logflag, option);
            auto end_remc = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed_remc = end_remc - start_remc;
            baseiclog.log("Processing time: " + to_string(elapsed_remc.count()) + " seconds, Final energy: " + to_string(final_energy) + "\n");
            //baseiclog.log("Saving to " +outdir+"/mc_"+to_string(iter)+".pdb\n");
            //write_pdb_structure((outdir+"/mc_"+to_string(iter)+".pdb").c_str(),pdb_entry);
        }

        if(option == 2 || option == 3){
            baseiclog.log("Energy-based refinement...\n");
            //double e_maxtime = 3600;
            auto start_ene = std::chrono::high_resolution_clock::now();
            StructureOptimization(pdb_entry, ideal_rna, inputtedstrall, inputtedstrallflag, templatestrall, templatestrallflag, bp_vec, loopnts, noncbp_vec, stemconnectnts, connectnts, stackmap, bpdistmap, stackdistmap, extparameters, hydrogencoords, e_maxtime, outdir, logflag, option);
            auto end_ene = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed_ene = end_ene - start_ene;
            baseiclog.log("Processing time: " + to_string(elapsed_ene.count()) + " seconds\n");
            //baseiclog.log("Saving to " +outdir+"/lbfgs_"+to_string(iter)+".pdb\n\n");
            //write_pdb_structure((outdir+"/lbfgs_"+to_string(iter)+".pdb").c_str(),pdb_entry);
            //exit(1);
        }
    }
    

    //write_pdb_structure(outfile.c_str(), pdb_entry);

    baseiclog.log("Saving to " + outdir + "\n");
    string final_energy_str=to_string(final_energy);
    write_pdb_structure((outdir+"/rnaref3d.pdb").c_str(), pdb_entry, final_energy_str);

    map<string, map<string, vector<ModelUnit>>> ().swap(ideal_rna);
    vector<ChainUnit>().swap(pdb_entry.chains);
    vector<vector<size_t> >().swap(res_str_vec);
    vector<pair<double,vector<size_t> > > ().swap(bp_vec);
    //}
    return 0;

}
