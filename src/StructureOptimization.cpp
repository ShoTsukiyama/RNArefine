
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
    const vector<vector<pair<int, double>>> &stackmap, const map<string, vector<vector<vector<double>>>> &bpdistmap, const map<string, map<string, vector<vector<vector<double>>>>> &stackdistmap, const ExtParams &extparameters, const vector<vector<vector<double>>> &hydrogencoords, const double maxtime, Logger &logger, const int option)
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
    
    //inputtedstr.residues[j] = inputtedstrall.residues[j];
    //inputtedstrflag[j] = inputtedstrallflag[j];
    moveresidues[j] = true;

    /*
    for(i=0;i<stembps.size();i++){
        for (j=0;j<stembps[i].size();j++)
        {
            c1= stembps[i][j][0];
            r1= stembps[i][j][1];
            c2= stembps[i][j][2];
            r2= stembps[i][j][3];

            if(c1 == c && c2 == c){
               bindmap[r1][r2]=2;
               bindmap[r2][r1]=2;
            }
        }
    }
    */

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

    AttachHydrogen(superimposedstr, hydrogencoords, extparameters.hydparam, extparameters.hydparamidx);
    //UpdateHydrogen(superimposedstr, extparameters.hydparam, extparameters.hydparamidx);
    fragmentoptimization(superimposedstr, inputtedstr_temp, inputtedstrflag_temp, templatestrall, templatestrallflag, ideal_rna, startene, endene, extparameters, bindmap, stackmap, cbpmap, noncbpmap, bpdistmap, stackdistmap, moveresidues, maxtime, 30000, logger);
    //cout<<"Fragment optimization is finished (whole structure)" << "starting energy: " << startene << ", final energy: " << endene<<endl;

    RemoveHydrogens(superimposedstr);

    for(int j=0;j<pdb_entry.chains[c].residues.size();j++){
        UpdateStructure(c, j, pdb_entry, superimposedstr.residues[j]);
        UpdateStructure(c, j, pdb_entry, superimposedstr.residues[j]);
    }

            /*

            double startene, endene;

            for (i=0;i<stembps.size();i++)
            {

                int baseleftprev_c = stemconnectnts[i][0];
                int baseleftprev_r = stemconnectnts[i][1];
                int baseleftnext_c = stemconnectnts[i][2];
                int baseleftnext_r = stemconnectnts[i][3];
                int baserightprev_c = stemconnectnts[i][4];
                int baserightprev_r = stemconnectnts[i][5];
                int baserightnext_c = stemconnectnts[i][6];
                int baserightnext_r = stemconnectnts[i][7];

                int addres=0; 
                int addres_lp=0; 
                int addres_ln=0; 
                int addres_rp=0; 
                int addres_rn=0;

                if(baseleftprev_r >= 0) addres_lp+=1;
                if(baseleftnext_r >= 0) addres_ln+=1;
                if(baserightprev_r >= 0) addres_rp+=1;
                if(baserightnext_r >= 0) addres_rn+=1;

                addres=addres_lp+addres_ln+addres_rp+addres_rn;

                int edgeleft_r = stembps[i][0][1];
                int bp_len = stembps[i].size();
                int bpnt_len = bp_len*2+addres;
                int edgeright_r = stembps[i][bp_len-1][3];

                FragmentStructure superimposedstr;
                FragmentStructure inputtedstr;
                vector<vector<bool>> inputtedstrflag;
                vector<vector<int>> bindmap(bpnt_len, vector<int>(bpnt_len, 0));
                vector<bool> moveresidues(bpnt_len, false);

                superimposedstr.residues.resize(bpnt_len);
                inputtedstr.residues.resize(bpnt_len);
                inputtedstrflag.resize(bpnt_len);

                int rel_pos, rel_pos_1, rel_pos_2;

                if(addres_lp > 0){
                    rel_pos=0;
                    ExtractNt(baseleftprev_c, baseleftprev_r, pdb_entry, superimposedstr.residues[rel_pos]);
                    inputtedstr.residues[rel_pos] = inputtedstrall.residues[baseleftprev_r];
                    inputtedstrflag[rel_pos] = inputtedstrallflag[baseleftprev_r];
                    bindmap[rel_pos][rel_pos+1] = 1;
                }

                if(addres_ln > 0){
                    rel_pos=addres_lp+bp_len;
                    ExtractNt(baseleftnext_c, baseleftnext_r, pdb_entry, superimposedstr.residues[rel_pos]);
                    inputtedstr.residues[rel_pos] = inputtedstrall.residues[baseleftnext_r];
                    inputtedstrflag[rel_pos] = inputtedstrallflag[baseleftnext_r];
                    bindmap[rel_pos-1][rel_pos] = 1;
                }

                if(addres_rp > 0){
                    rel_pos=addres_lp+bp_len+addres_ln;
                    ExtractNt(baserightprev_c, baserightprev_r, pdb_entry, superimposedstr.residues[rel_pos]);
                    inputtedstr.residues[rel_pos] = inputtedstrall.residues[baserightprev_r];
                    inputtedstrflag[rel_pos] = inputtedstrallflag[baserightprev_r];
                    bindmap[rel_pos][rel_pos+1] = 1;
                }

                if(addres_rn > 0){
                    rel_pos=addres_lp+2*bp_len+addres_ln+addres_rn;
                    ExtractNt(baserightnext_c, baserightnext_r, pdb_entry, superimposedstr.residues[rel_pos]);
                    inputtedstr.residues[rel_pos] = inputtedstrall.residues[baserightnext_r];
                    inputtedstrflag[rel_pos] = inputtedstrallflag[baserightnext_r];
                    bindmap[rel_pos-1][rel_pos] = 1;
                }

                for (j=0;j<stembps[i].size();j++)
                {
                    c1= stembps[i][j][0];
                    r1= stembps[i][j][1];
                    c2= stembps[i][j][2];
                    r2= stembps[i][j][3];

                    cout<<r1 << "-" << r2<<endl;

                    rel_pos_1 = addres_lp+r1-edgeleft_r;
                    rel_pos_2 = addres_lp+addres_ln+addres_rp+r2-edgeright_r+bp_len;

                    ExtractNt(c1, r1, pdb_entry, superimposedstr.residues[rel_pos_1]);
                    ExtractNt(c2, r2, pdb_entry, superimposedstr.residues[rel_pos_2]);
                    inputtedstr.residues[rel_pos_1] = inputtedstrall.residues[r1];
                    inputtedstrflag[rel_pos_1] = inputtedstrallflag[r1];
                    inputtedstr.residues[rel_pos_2] = inputtedstrall.residues[r2];
                    inputtedstrflag[rel_pos_2] = inputtedstrallflag[r2];

                    moveresidues[rel_pos_1] = true;
                    moveresidues[rel_pos_2] = true;

                    bindmap[rel_pos_1][rel_pos_2]=2;
                    bindmap[rel_pos_2][rel_pos_1]=2;

                    if(j != bp_len-1 || addres_ln > 0){
                        //cout<<rel_pos_1 << "-" << rel_pos_1+1<<endl;
                        bindmap[rel_pos_1][rel_pos_1+1]=1;
                    }
                    if(j != bp_len-1 || addres_rn > 0){
                        //cout<<rel_pos_2 << "-" << rel_pos_2+1<<endl;
                        bindmap[rel_pos_2][rel_pos_2+1]=1;
                    }
                }

                for(j=0;j<bindmap.size();j++){
                    for(int k=0;k<bindmap[j].size();k++){
                        cout<<bindmap[j][k] << ",";
                    }
                    cout<<""<<endl;
                }

                AttachHydrogen(superimposedstr, extparameters.hydparam, extparameters.hydparamidx);
                
                //UpdateHydrogen(superimposedstr, extparameters.hydparam, extparameters.hydparamidx);
                
                fragmentoptimizationbpopt(superimposedstr, inputtedstr, inputtedstrflag, startene, endene, extparameters, bindmap, moveresidues, 2500);
                cout<<"Fragment optimization is finished (stem " << i+1 << "), " << "starting energy: " << startene << ", final energy: " << endene<<endl;
    
                RemoveHydrogens(superimposedstr);

                for (j=0;j<stembps[i].size();j++)
                {
                    c1= stembps[i][j][0];
                    r1= stembps[i][j][1];
                    c2= stembps[i][j][2];
                    r2= stembps[i][j][3];

                    rel_pos_1 = addres_lp+r1-edgeleft_r;
                    rel_pos_2 = addres_lp+addres_ln+addres_rp+r2-edgeright_r+bp_len;

                    UpdateStructure(c1, r1, pdb_entry, superimposedstr.residues[rel_pos_1]);
                    UpdateStructure(c2, r2, pdb_entry, superimposedstr.residues[rel_pos_2]);
                
                }
            }

            for (auto& it : loopnts)
            {
                c1=it.first;
                for(i=0;i<it.second.size();i++){
                    int edgeleft_r = it.second[i].first;
                    int edgeright_r = it.second[i].second;

                    int outedge_left_r = connectnts.at(c1)[i][0][0];
                    int outedge_left_c = connectnts.at(c1)[i][0][2];
                    int outedge_left_pair_r = connectnts.at(c1)[i][0][1];
                    int outedge_left_pair_c = connectnts.at(c1)[i][0][3];
                    int pairtype_left = connectnts.at(c1)[i][0][4];

                    int outedge_right_r = connectnts.at(c1)[i][1][0];
                    int outedge_right_c = connectnts.at(c1)[i][1][2];
                    int outedge_right_pair_r = connectnts.at(c1)[i][1][1];
                    int outedge_right_pair_c = connectnts.at(c1)[i][1][3];
                    int pairtype_right = connectnts.at(c1)[i][1][4];

                    int addres=0;
                    int addres_f=0;
                    int addres_e=0;

                    if(pairtype_left == 2){
                        addres_f+=1;
                        addres_e+=1;
                    }else if(pairtype_left == 1 && pairtype_right == 1){
                        addres_f+=2;
                        addres_e+=2;
                    }else if(pairtype_left == 1){
                        addres_f+=2;
                    }else if(pairtype_right == 1){
                        addres_e+=2;
                    }
                    addres = addres_f+addres_e;

                    FragmentStructure superimposedstr;
                    FragmentStructure inputtedstr;
                    vector<vector<bool>> inputtedstrflag;
                    int loop_len = edgeright_r-edgeleft_r+1+addres;

                    vector<vector<int>> bindmap(loop_len, vector<int>(loop_len, 0));
                    vector<bool> moveresidues(loop_len, false);
                    superimposedstr.residues.resize(loop_len);
                    inputtedstr.residues.resize(loop_len);
                    inputtedstrflag.resize(loop_len);

                    if(pairtype_left == 2)
                    {
                        ExtractNt(c1, outedge_left_r, pdb_entry, superimposedstr.residues[0]);
                        inputtedstr.residues[0] = inputtedstrall.residues[outedge_left_r];
                        inputtedstrflag[0] = inputtedstrallflag[outedge_left_r];
                        
                        ExtractNt(c1, outedge_left_pair_r, pdb_entry, superimposedstr.residues[loop_len-1]);
                        inputtedstr.residues[loop_len-1] = inputtedstrall.residues[outedge_left_pair_r];
                        inputtedstrflag[loop_len-1] = inputtedstrallflag[outedge_left_pair_r];
                        
                        bindmap[0][loop_len-1] = 2;
                        bindmap[loop_len-1][0] = 2;
                        bindmap[0][1] = 1;

                        cout<<outedge_left_r<<endl;
                        cout<<outedge_left_pair_r<<endl;
                    }else if(pairtype_left == 1)
                    {
                        ExtractNt(outedge_left_pair_c, outedge_left_pair_r, pdb_entry, superimposedstr.residues[0]);
                        inputtedstr.residues[0] = inputtedstrall.residues[outedge_left_pair_r];
                        inputtedstrflag[0] = inputtedstrallflag[outedge_left_pair_r];

                        ExtractNt(outedge_left_c, outedge_left_r, pdb_entry, superimposedstr.residues[1]);
                        inputtedstr.residues[1] = inputtedstrall.residues[outedge_left_r];
                        inputtedstrflag[1] = inputtedstrallflag[outedge_left_r];
                        
                        bindmap[0][1] = 2;
                        bindmap[1][0] = 2;
                        bindmap[1][2] = 1;

                        cout<<outedge_left_pair_r<<endl;
                        cout<<outedge_left_r<<endl;
                    }

                    for(j=edgeleft_r;j<edgeright_r+1;j++)
                    {
                        ExtractNt(c1, j, pdb_entry, superimposedstr.residues[addres_f+j-edgeleft_r]);

                        inputtedstr.residues[addres_f+j-edgeleft_r] = inputtedstrall.residues[j];
                        inputtedstrflag[addres_f+j-edgeleft_r] = inputtedstrallflag[j];
                        
                        moveresidues[addres_f+j-edgeleft_r] = true;
                        if(j != edgeright_r){
                            bindmap[addres_f+j-edgeleft_r][addres_f+j-edgeleft_r+1] = 1;
                        }

                        cout<<j<<endl;
                    }

                    if(pairtype_right == 1)
                    {
                        //AssembleBasePairFrag(outedge_right_c,outedge_right_r,outedge_right_pair_c,outedge_right_pair_r, pdb_entry, ideal_rna["bp"], superimposedstr.residues[loop_len-2], superimposedstr.residues[loop_len-1], inputtedstr.residues[loop_len-2], inputtedstr.residues[loop_len-1], inputtedstrflag[loop_len-2], inputtedstrflag[loop_len-1], templateidx);
                        ExtractNt(outedge_right_c, outedge_right_r, pdb_entry, superimposedstr.residues[loop_len-2]);
                        inputtedstr.residues[loop_len-2] = inputtedstrall.residues[outedge_right_r];
                        inputtedstrflag[loop_len-2] = inputtedstrallflag[outedge_right_r];

                        ExtractNt(outedge_right_pair_c, outedge_right_pair_r, pdb_entry, superimposedstr.residues[loop_len-1]);
                        inputtedstr.residues[loop_len-1] = inputtedstrall.residues[outedge_right_pair_r];
                        inputtedstrflag[loop_len-1] = inputtedstrallflag[outedge_right_pair_r];

                        bindmap[loop_len-2][loop_len-1] = 2;
                        bindmap[loop_len-1][loop_len-2] = 2;
                        bindmap[loop_len-3][loop_len-2] = 1;

                        cout<<outedge_right_r<<endl;
                        cout<<outedge_right_pair_r<<endl;
                    }

                    AttachHydrogen(superimposedstr, extparameters.hydparam, extparameters.hydparamidx);
                    //UpdateHydrogen(superimposedstr, extparameters.hydparam, extparameters.hydparamidx);

                    fragmentoptimizationloop(superimposedstr, inputtedstr, inputtedstrflag, loop_len, startene, endene, extparameters, bindmap, moveresidues, 2500);
                    cout<<"Fragment optimization is finished (loop " << i+1 << "), " << "starting energy: " << startene << ", final energy: " << endene<<endl;

                    RemoveHydrogens(superimposedstr);

                    for(j=edgeleft_r;j<edgeright_r+1;j++){
                        UpdateStructure(c1, j, pdb_entry, superimposedstr.residues[addres_f+j-edgeleft_r]);
                        
                    }
                }
            }
            */

    return true;
}

bool StructureOptimization(ModelUnit &pdb_entry, const map<string, map<string, vector<ModelUnit>>> &ideal_rna, 
    const FragmentStructure &inputtedstrall, const vector<vector<bool>> &inputtedstrallflag,
    const vector<FragmentStructure> &templatestrall, const vector<vector<vector<bool>>> &templatestrallflag,
    const vector<pair<double,vector<size_t>>> &bp_vec, const map<int, vector<pair<int, int>>> &loopnts, const vector<pair<double,vector<size_t>>> &noncbp_vec,
    const vector<vector<int>> &stemconnectnts, const map<int, vector<vector<vector<int>>>> &connectnts,
    const vector<vector<pair<int, double>>> &stackmap, const map<string, vector<vector<vector<double>>>> &bpdistmap, const map<string, map<string, vector<vector<vector<double>>>>> &stackdistmap, const ExtParams &extparameters, const vector<vector<vector<double>>> &hydrogencoords, const double maxtime, const string &outdir, const bool &logflag, const int option)
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
    StructureOptimization(pdb_entry, c, ideal_rna, inputtedstrall, inputtedstrallflag, templatestrall, templatestrallflag, bp_vec, loopnts, noncbp_vec, stemconnectnts, connectnts, stackmap, bpdistmap, stackdistmap, extparameters, hydrogencoords, maxtime, logger, option);

    return true;
}


