
#include "StructureOptimization.hpp"

void UpdateStructure(const int c, const int r, ModelUnit &pdb_entry, FragmentResidue &spresidue)
{
    int a, a1, a2;
    string base1, base2;
    int atom;

    base1=pdb_entry.chains[c].residues[r].resn[2];
    for (a=0;a<pdb_entry.chains[c].residues[r].atoms.size();a++)
    {
        if(pdb_entry.chains[c].residues[r].atoms[a].movable)
        {
            pdb_entry.chains[c].residues[r].atoms[a].xyz[0] = spresidue.atoms[a].coord.x;
            pdb_entry.chains[c].residues[r].atoms[a].xyz[1] = spresidue.atoms[a].coord.y;
            pdb_entry.chains[c].residues[r].atoms[a].xyz[2] = spresidue.atoms[a].coord.z;
        }
    }
}

void GetHydrogenInfo(const int c, const int r, ModelUnit pdb_entry, FragmentResidue &spresidue, vector<vector<double>> &hydrogencoords)
{
    int a, a1, a2;
    string base;
    int idx;
    int atom;

    base=pdb_entry.chains[c].residues[r].resn[2];
    hydrogencoords.resize(full_atoms[base].size()-heavily_atom_num[base], vector<double>(3, 0));
    for (a=heavily_atom_num[base];a<full_atoms[base].size();a++)
    {
        idx = a-heavily_atom_num[base];
        hydrogencoords[idx][0] = spresidue.atoms[a].coord.x;
        hydrogencoords[idx][1] = spresidue.atoms[a].coord.y;
        hydrogencoords[idx][2] = spresidue.atoms[a].coord.z;
        //cout<<idx << "-" << hydrogencoords[idx][0] << "-" << hydrogencoords[idx][1] << "-" << hydrogencoords[idx][2]<<endl;
    }
}

bool StructureOptimization(ModelUnit &pdb_entry, const size_t c, const map<string, map<string, vector<ModelUnit>>> &ideal_rna,
    const FragmentStructure &inputtedstrall, const vector<vector<bool>> &inputtedstrallflag,
    const vector<FragmentStructure> &templatestrall, const vector<vector<vector<bool>>> &templatestrallflag,
    const vector<pair<double,vector<size_t>>> &bp_vec, const map<int, vector<pair<int, int>>> &loopnts, const vector<pair<double,vector<size_t>>> &noncbp_vec, 
    const vector<vector<int>> &stemconnectnts, const map<int, vector<vector<vector<int>>>> &connectnts,
    const vector<vector<pair<int, double>>> &stackmap, const map<string, vector<vector<vector<double>>>> &bpdistmap, const map<string, map<string, vector<vector<vector<double>>>>> &stackdistmap, const ExtParams &extparameters, const vector<double> &init_weight, const vector<vector<vector<double>>> &hydrogencoords, const double maxtime, Logger &logger, const int option)
    //map<string, vector<vector<double>>> &blengthparam, map<string, vector<vector<double>>> &adjblengthparam, 
    //map<string, vector<vector<double>>> &ntangleparam, map<string, vector<vector<int>>> &ntangleparamidx, map<string, vector<vector<double>>> &adjangleparam, map<string, vector<vector<int>>> &adjangleparamidx, 
    //map<string, vector<vector<double>>> &hydparam, map<string, vector<vector<int>>> &hydparamidx, const int option=0)
{
    size_t r,i,j,k;
    int r1,r2,c1,c2;
    size_t a;
    string base;
    int atom;
    double startene, endene;
    string restype;
    int resnum=pdb_entry.chains[c].residues.size();
    FragmentStructure superimposedstr;
    FragmentStructure inputtedstr_temp;
    vector<vector<bool>> inputtedstrflag_temp;
    vector<vector<int>> bindmap(resnum, vector<int>(resnum, 0));
    vector<vector<double>> cbpmap(resnum, vector<double>(resnum, 0));
    vector<vector<double>> noncbpmap(resnum, vector<double>(resnum, 0));
    vector<bool> moveresidues(resnum, false);
    double score;

    static map<string, vector<int>> init_target_atoms_idx={
        {"A", {0, 3, 4, 5, 7, 8}},
        {"C", {0, 3, 4, 5, 7, 8}},
        {"G", {0, 3, 4, 5, 7, 8}},
        {"U", {0, 3, 4, 5, 7, 8}}
    };
        
    superimposedstr.residues.resize(resnum);
    inputtedstr_temp.residues.resize(resnum);
    inputtedstrflag_temp.resize(resnum);

    for(j=0;j<pdb_entry.chains[c].residues.size();j++)
    {
        base=pdb_entry.chains[c].residues[j].resn[2];
        inputtedstrflag_temp[j].resize(full_atoms[base].size(), false);
    }

    for(int j=0;j<pdb_entry.chains[c].residues.size()-1;j++)
    {
        ExtractNt(c, j, pdb_entry, superimposedstr.residues[j]);
        //ExtractNt(c, j, pdb_entry, inputtedstr.residues[j]);

        //inputtedstr.residues[j] = inputtedstrall.residues[j];
        //inputtedstrflag[j] = inputtedstrallflag[j];

        //base=pdb_entry.chains[c].residues[j].resn[2];
        //inputtedstrflag[j].resize(heavily_atom_num[base], false);


        base = inputtedstrall.residues[j].restype;
        inputtedstr_temp.residues[j].restype = base;
        for(k=0;k<init_target_atoms_idx[base].size();k++){
            atom = init_target_atoms_idx[base][k];

            inputtedstr_temp.residues[j].atoms.resize(full_atoms[base].size());
            inputtedstrflag_temp[j].resize(full_atoms[base].size());

            if(inputtedstrallflag[j][atom]){
                inputtedstr_temp.residues[j].atoms[atom].coord.x = inputtedstrall.residues[j].atoms[atom].coord.x;
                inputtedstr_temp.residues[j].atoms[atom].coord.y = inputtedstrall.residues[j].atoms[atom].coord.y;
                inputtedstr_temp.residues[j].atoms[atom].coord.z = inputtedstrall.residues[j].atoms[atom].coord.z;
                inputtedstrflag_temp[j][atom] = true;
            }
        }

        bindmap[j][j+1] = 1;
        moveresidues[j] = true;
    }

    
    j=pdb_entry.chains[c].residues.size()-1;
    ExtractNt(c, j, pdb_entry, superimposedstr.residues[j]);
    //ExtractNt(c, j, pdb_entry, inputtedstr.residues[j]);

    base = inputtedstrall.residues[j].restype;
    inputtedstr_temp.residues[j].restype = base;
    for(k=0;k<init_target_atoms_idx[base].size();k++){
        atom = init_target_atoms_idx[base][k];

        inputtedstr_temp.residues[j].atoms.resize(full_atoms[base].size());
        inputtedstrflag_temp[j].resize(full_atoms[base].size());

        if(inputtedstrallflag[j][atom]){
            inputtedstr_temp.residues[j].atoms[atom].coord.x = inputtedstrall.residues[j].atoms[atom].coord.x;
            inputtedstr_temp.residues[j].atoms[atom].coord.y = inputtedstrall.residues[j].atoms[atom].coord.y;
            inputtedstr_temp.residues[j].atoms[atom].coord.z = inputtedstrall.residues[j].atoms[atom].coord.z;
            inputtedstrflag_temp[j][atom] = true;
        }
    }
    
    moveresidues[j] = true;

    for(i=0;i<bp_vec.size();i++){
        score = bp_vec[i].first;
        c1 = bp_vec[i].second[0];
        r1 = bp_vec[i].second[1];
        c2 = bp_vec[i].second[2];
        r2 = bp_vec[i].second[3];

        if(c1 == c && c2 == c){
            cbpmap[r1][r2]=score;
            cbpmap[r2][r1]=score;
        }
    }

    for(i=0;i<noncbp_vec.size();i++){
        score = noncbp_vec[i].first;
        c1 = noncbp_vec[i].second[0];
        r1 = noncbp_vec[i].second[1];
        c2 = noncbp_vec[i].second[2];
        r2 = noncbp_vec[i].second[3];

        if(c1 == c && c2 == c){
            noncbpmap[r1][r2]=score;
            noncbpmap[r2][r1]=score;
        }
    }

    //if(option == 3){
        //AttachHydrogen(superimposedstr, hydrogencoords, extparameters.hydparam, extparameters.hydparamidx);
    //}else{
    AttachHydrogen(superimposedstr, extparameters.hydparam, extparameters.hydparamidx);
    //}
    //UpdateHydrogen(superimposedstr, extparameters.hydparam, extparameters.hydparamidx);
    fragmentoptimization(superimposedstr, inputtedstr_temp, inputtedstrflag_temp, templatestrall, templatestrallflag, init_weight, ideal_rna, startene, endene, extparameters, bindmap, stackmap, cbpmap, noncbpmap, bpdistmap, stackdistmap, moveresidues, maxtime, 30000, logger);
    //cout<<"Fragment optimization is finished (whole structure)" << "starting energy: " << startene << ", final energy: " << endene<<endl;

    RemoveHydrogens(superimposedstr);

    for(int j=0;j<pdb_entry.chains[c].residues.size();j++){
        UpdateStructure(c, j, pdb_entry, superimposedstr.residues[j]);
        UpdateStructure(c, j, pdb_entry, superimposedstr.residues[j]);
    }

    return true;
}

bool StructureOptimization(ModelUnit &pdb_entry, const map<string, map<string, vector<ModelUnit>>> &ideal_rna, 
    const FragmentStructure &inputtedstrall, const vector<vector<bool>> &inputtedstrallflag,
    const vector<FragmentStructure> &templatestrall, const vector<vector<vector<bool>>> &templatestrallflag,
    const vector<pair<double,vector<size_t>>> &bp_vec, const map<int, vector<pair<int, int>>> &loopnts, const vector<pair<double,vector<size_t>>> &noncbp_vec,
    const vector<vector<int>> &stemconnectnts, const map<int, vector<vector<vector<int>>>> &connectnts,
    const vector<vector<pair<int, double>>> &stackmap, const map<string, vector<vector<vector<double>>>> &bpdistmap, const map<string, map<string, vector<vector<vector<double>>>>> &stackdistmap, const ExtParams &extparameters, const vector<double> &init_weight, const vector<vector<vector<double>>> &hydrogencoords, const double maxtime, const string &outdir, const bool &logflag, const int option)
    //map<string, vector<vector<double>>> &blengthparam, map<string, vector<vector<double>>> &adjblengthparam, 
    //map<string, vector<vector<double>>> &ntangleparam, map<string, vector<vector<int>>> &ntangleparamidx, map<string, vector<vector<double>>> &adjangleparam, map<string, vector<vector<int>>> &adjangleparamidx, 
    //map<string, vector<vector<double>>> &hydparam, map<string, vector<vector<int>>> &hydparamidx, const int option=0)
{
    Logger logger = Logger();
    logger.setPrint(true);
    if(logflag){
        logger.setVerbose(true);
        logger.setLogFile(outdir+"/logs/l-bfgs_refinement.log");
    }

    size_t c;
    for (c=0;c<pdb_entry.chains.size();c++)
    StructureOptimization(pdb_entry, c, ideal_rna, inputtedstrall, inputtedstrallflag, templatestrall, templatestrallflag, bp_vec, loopnts, noncbp_vec, stemconnectnts, connectnts, stackmap, bpdistmap, stackdistmap, extparameters, init_weight, hydrogencoords, maxtime, logger, option);

    return true;
}


