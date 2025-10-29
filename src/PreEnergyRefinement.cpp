
#include "PreEnergyRefinement.hpp"

void UpdateNtStructure(const int c1, const int r1, ModelUnit &pdb_entry, FragmentResidue &spresidue)
{
    for (int a=0;a<pdb_entry.chains[c1].residues[r1].atoms.size();a++)
    {
        pdb_entry.chains[c1].residues[r1].atoms[a].xyz[0]=spresidue.atoms[a].coord.x;
        pdb_entry.chains[c1].residues[r1].atoms[a].xyz[1]=spresidue.atoms[a].coord.y;
        pdb_entry.chains[c1].residues[r1].atoms[a].xyz[2]=spresidue.atoms[a].coord.z;
    }
}

void ExtractNtPre(const int c, const int r, const ModelUnit &pdb_entry, FragmentResidue &spresidue)
{
    int a;
    int atom;
    string base, base_pdb;

    base=pdb_entry.chains[c].residues[r].resn[2];
    base_pdb=pdb_entry.chains[c].residues[r].resn;

    FragmentAtoms coordtmp;
    coordtmp.coord.x = 0; coordtmp.coord.y = 0; coordtmp.coord.z = 0;
    spresidue.restype = base;
    spresidue.atoms.resize(coarse_grained_atoms_pdb[base_pdb].size(), coordtmp);
    
    for (a=0;a<pdb_entry.chains[c].residues[r].atoms.size();a++){
        spresidue.atoms[a].coord.x = pdb_entry.chains[c].residues[r].atoms[a].xyz[0];
        spresidue.atoms[a].coord.y = pdb_entry.chains[c].residues[r].atoms[a].xyz[1];
        spresidue.atoms[a].coord.z = pdb_entry.chains[c].residues[r].atoms[a].xyz[2];
    }
}

bool PreEnergyRefine(ModelUnit &pdb_entry, const vector<FragmentStructure> &templatestrall, const vector<vector<vector<bool>>> &templatestrallflag, const map<string, map<string, vector<ModelUnit>>> &ideal_rna, const size_t c,
    vector<vector<int>> &stackmap, const map<string, map<string, vector<vector<vector<double>>>>> &stackdistmap, const vector<vector<vector<int>>> &stembps, const map<int, vector<pair<int, int>>> &loopnts, const map<int, vector<vector<vector<int>>>> &connectnts,
    ExtParams &extparameters)
    //map<string, vector<vector<double>>> &blengthparam, map<string, vector<vector<double>>> &adjblengthparam, 
    //map<string, vector<vector<double>>> &ntangleparam, map<string, vector<vector<int>>> &ntangleparamidx, map<string, vector<vector<double>>> &adjangleparam, map<string, vector<vector<int>>> &adjangleparamidx, 
    //map<string, vector<vector<double>>> &hydparam, map<string, vector<vector<int>>> &hydparamidx, const int option=0)
{
    size_t i,j,k;
    int r1,r2,c1,c2,rr1,rr2;
    double startene, endene;
    string base;

    for (i=0;i<stembps.size();i++)
    {
        FragmentStructure superimposedstr;
        FragmentStructure inputtedstr;
        vector<vector<bool>> inputtedstrflag;
        int edgeleft_r = stembps[i][0][1];
        int bp_len = stembps[i].size();
        int bpnt_len = bp_len*2;
        vector<vector<int>> bindmap(bpnt_len, vector<int>(bpnt_len, 0));
        vector<vector<int>> stackingmap(bpnt_len, vector<int>(bpnt_len, 0));
        vector<bool> moveresidues(bpnt_len, false);

        int edgeright_r = stembps[i][bp_len-1][3];

        superimposedstr.residues.resize(bpnt_len);
        inputtedstr.residues.resize(bpnt_len);
        inputtedstrflag.resize(bpnt_len);

        vector<vector<int>> ntidx;
        for (j=0;j<stembps[i].size();j++)
        {
            c1= stembps[i][j][0];
            r1= stembps[i][j][1];
            c2= stembps[i][j][2];
            r2= stembps[i][j][3];

            ExtractNtPre(c1, r1, pdb_entry, superimposedstr.residues[r1-edgeleft_r]);
            ExtractNtPre(c2, r2, pdb_entry, superimposedstr.residues[r2-edgeright_r+bp_len]);

            ExtractNtPre(c1, r1, pdb_entry, inputtedstr.residues[r1-edgeleft_r]);
            ExtractNtPre(c2, r2, pdb_entry, inputtedstr.residues[r2-edgeright_r+bp_len]);

            base=pdb_entry.chains[c1].residues[r1].resn[2];
            inputtedstrflag[r1-edgeleft_r].resize(heavily_atom_num[base], false);
            inputtedstrflag[r1-edgeleft_r][0] = true;

            base=pdb_entry.chains[c2].residues[r2].resn[2];
            inputtedstrflag[r2-edgeright_r+bp_len].resize(heavily_atom_num[base], false);
            inputtedstrflag[r2-edgeright_r+bp_len][0] = true;

            //bp_template_idx[bp]=templateidx;

            moveresidues[r1-edgeleft_r] = true;
            moveresidues[r2-edgeright_r+bp_len] = true;

            bindmap[r1-edgeleft_r][r2-edgeright_r+bp_len]=2;
            bindmap[r2-edgeright_r+bp_len][r1-edgeleft_r]=2;

            if(j != bp_len-1)
            {
                bindmap[j][j+1]=1;
                bindmap[bp_len+j][bp_len+j+1]=1;
            }

            ntidx.push_back({r1,r1-edgeleft_r});
            ntidx.push_back({r2,r2-edgeright_r+bp_len});
        }

        for(j=0;j<ntidx.size();j++){
            for(k=j;k<ntidx.size();k++){
                r1 = ntidx[j][0];
                r2 = ntidx[k][0];
                rr1 = ntidx[j][1];
                rr2 = ntidx[k][1];

                stackingmap[rr1][rr2] = stackmap[r1][r2];
                stackingmap[rr2][rr1] = stackmap[r2][r1];
            }
        }

        AttachHydrogen(superimposedstr, extparameters.hydparam, extparameters.hydparamidx);   
        //fragmentoptimizationbp(superimposedstr, inputtedstr, inputtedstrflag, templatestrall, templatestrallflag, ideal_rna, bp_len, startene, endene, extparameters, bindmap, stackingmap, stackdistmap, moveresidues, 5000);
        cout<<"Fragment optimization is finished (stem " << i+1 << "), " << "starting energy: " << startene << ", final energy: " << endene<<endl;
        RemoveHydrogens(superimposedstr);

        for (j=0;j<stembps[i].size();j++)
        {
            c1= stembps[i][j][0];
            r1= stembps[i][j][1];
            c2= stembps[i][j][2];
            r2= stembps[i][j][3];

            UpdateNtStructure(c1, r1, pdb_entry, superimposedstr.residues[r1-edgeleft_r]);
            UpdateNtStructure(c2, r2, pdb_entry, superimposedstr.residues[r2-edgeright_r+bp_len]);
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
            vector<vector<int>> stackingmap(loop_len, vector<int>(loop_len, 0));
            vector<bool> moveresidues(loop_len, false);
            superimposedstr.residues.resize(loop_len);
            inputtedstr.residues.resize(loop_len);
            inputtedstrflag.resize(loop_len);
            for(j=0;j<loop_len;j++){
                base=pdb_entry.chains[c].residues[j].resn[2];
                inputtedstrflag[j].resize(heavily_atom_num[base], false);
            }

            vector<vector<int>> ntidx;
            if(pairtype_left == 2)
            {
                ExtractNtPre(c1, outedge_left_r, pdb_entry, superimposedstr.residues[0]);
                //inputtedstr.residues[0] = inputtedstrall.residues[outedge_left_r];
                //inputtedstrflag[0] = inputtedstrallflag[outedge_left_r];
                ExtractNtPre(c1, outedge_left_r, pdb_entry, inputtedstr.residues[0]);

                base=pdb_entry.chains[c1].residues[outedge_left_r].resn[2];
                inputtedstrflag[0].resize(heavily_atom_num[base], false);
                inputtedstrflag[0][0] = true;
            
                ExtractNtPre(c1, outedge_left_pair_r, pdb_entry, superimposedstr.residues[loop_len-1]);
                //inputtedstr.residues[loop_len-1] = inputtedstrall.residues[outedge_left_pair_r];
                //inputtedstrflag[loop_len-1] = inputtedstrallflag[outedge_left_pair_r];
                ExtractNtPre(c1, outedge_left_pair_r, pdb_entry, inputtedstr.residues[loop_len-1]);

                base=pdb_entry.chains[c1].residues[outedge_left_pair_r].resn[2];
                inputtedstrflag[loop_len-1].resize(heavily_atom_num[base], false);
                inputtedstrflag[loop_len-1][0] = true;
                
                bindmap[0][loop_len-1] = 2;
                bindmap[loop_len-1][0] = 2;
                bindmap[0][1] = 1;

                ntidx.push_back({outedge_left_r,0});
                ntidx.push_back({outedge_left_pair_r,loop_len-1});

            }else if(pairtype_left == 1)
            {
                ExtractNtPre(outedge_left_pair_c, outedge_left_pair_r, pdb_entry, superimposedstr.residues[0]);
                //inputtedstr.residues[0] = inputtedstrall.residues[outedge_left_pair_r];
                //inputtedstrflag[0] = inputtedstrallflag[outedge_left_pair_r];
                ExtractNtPre(outedge_left_pair_c, outedge_left_pair_r, pdb_entry, inputtedstr.residues[0]);

                base=pdb_entry.chains[outedge_left_pair_c].residues[outedge_left_pair_r].resn[2];
                inputtedstrflag[0].resize(heavily_atom_num[base], false);
                inputtedstrflag[0][0] = true;

                ExtractNtPre(outedge_left_c, outedge_left_r, pdb_entry, superimposedstr.residues[1]);
                //inputtedstr.residues[1] = inputtedstrall.residues[outedge_left_r];
                //inputtedstrflag[1] = inputtedstrallflag[outedge_left_r];
                ExtractNtPre(outedge_left_c, outedge_left_r, pdb_entry, inputtedstr.residues[1]);

                base=pdb_entry.chains[outedge_left_c].residues[outedge_left_r].resn[2];
                inputtedstrflag[1].resize(heavily_atom_num[base], false);
                inputtedstrflag[1][0] = true;
                
                bindmap[0][1] = 2;
                bindmap[1][0] = 2;
                bindmap[1][2] = 1;

                ntidx.push_back({outedge_left_pair_r,0});
                ntidx.push_back({outedge_left_r,1});
            }

            for(j=edgeleft_r;j<edgeright_r+1;j++)
            {
                ExtractNtPre(c1, j, pdb_entry, superimposedstr.residues[addres_f+j-edgeleft_r]);
                //inputtedstr.residues[addres_f+j-edgeleft_r] = inputtedstrall.residues[j];
                //inputtedstrflag[addres_f+j-edgeleft_r] = inputtedstrallflag[j];
                ExtractNtPre(c1, j, pdb_entry, inputtedstr.residues[addres_f+j-edgeleft_r]);

                base=pdb_entry.chains[c1].residues[j].resn[2];
                inputtedstrflag[addres_f+j-edgeleft_r].resize(heavily_atom_num[base], false);
                inputtedstrflag[addres_f+j-edgeleft_r][0] = true;
                
                moveresidues[addres_f+j-edgeleft_r] = true;
                if(j != edgeright_r){
                    bindmap[addres_f+j-edgeleft_r][addres_f+j-edgeleft_r+1] = 1;
                }

                r1 = j;
                rr1 = addres_f+j-edgeleft_r;
                ntidx.push_back({r1, rr1});
            }

            if(pairtype_right == 1)
            {
                ExtractNtPre(outedge_right_c, outedge_right_r, pdb_entry, superimposedstr.residues[loop_len-2]);
                //inputtedstr.residues[loop_len-2] = inputtedstrall.residues[outedge_right_r];
                //inputtedstrflag[loop_len-2] = inputtedstrallflag[outedge_right_r];
                ExtractNtPre(outedge_right_c, outedge_right_r, pdb_entry, inputtedstr.residues[loop_len-2]);

                base=pdb_entry.chains[outedge_right_c].residues[outedge_right_r].resn[2];
                inputtedstrflag[loop_len-2].resize(heavily_atom_num[base], false);
                inputtedstrflag[loop_len-2][0] = true;

                ExtractNtPre(outedge_right_pair_c, outedge_right_pair_r, pdb_entry, superimposedstr.residues[loop_len-1]);
                //inputtedstr.residues[loop_len-1] = inputtedstrall.residues[outedge_right_pair_r];
                //inputtedstrflag[loop_len-1] = inputtedstrallflag[outedge_right_pair_r];
                ExtractNtPre(outedge_right_pair_c, outedge_right_pair_r, pdb_entry, inputtedstr.residues[loop_len-1]);

                base=pdb_entry.chains[outedge_right_pair_c].residues[outedge_right_pair_r].resn[2];
                inputtedstrflag[loop_len-1].resize(heavily_atom_num[base], false);
                inputtedstrflag[loop_len-1][0] = true;

                bindmap[loop_len-2][loop_len-1] = 2;
                bindmap[loop_len-1][loop_len-2] = 2;
                bindmap[loop_len-3][loop_len-2] = 1;

                ntidx.push_back({outedge_right_r, loop_len-2});
                ntidx.push_back({outedge_right_pair_r, loop_len-1});
            }

            for(j=0;j<ntidx.size();j++){
                for(k=j;k<ntidx.size();k++){
                    r1 = ntidx[j][0];
                    r2 = ntidx[k][0];
                    rr1 = ntidx[j][1];
                    rr2 = ntidx[k][1];
    
                    stackingmap[rr1][rr2] = stackmap[r1][r2];
                    stackingmap[rr2][rr1] = stackmap[r2][r1];
                }
            }



            //for(j=0;j<loop_len;j++){
                //base=pdb_entry.chains[c].residues[j].resn[2];
                //for(int k=0;k<heavily_atom_num[base];k++){
                    //cout<<inputtedstrflag[j][k] << "-" << inputtedstr.residues[j].atoms[k].coord.x << "-" <<  superimposedstr.residues[j].atoms[k].coord.x<<endl;
                //}
            //}
            //exit(1);

            AttachHydrogen(superimposedstr, extparameters.hydparam, extparameters.hydparamidx);
            //fragmentoptimizationloop(superimposedstr, inputtedstr, inputtedstrflag, templatestrall, templatestrallflag, ideal_rna, loop_len, startene, endene, extparameters, bindmap, stackingmap, stackdistmap, moveresidues, 5000);
            cout<<"Fragment optimization is finished (loop " << i+1 << "), " << "starting energy: " << startene << ", final energy: " << endene<<endl;
            RemoveHydrogens(superimposedstr);

            for(j=edgeleft_r;j<edgeright_r+1;j++){
                UpdateNtStructure(c1, j, pdb_entry, superimposedstr.residues[addres_f+j-edgeleft_r]);
            }
        }
    }

    return true;
}

bool PreEnergyRefine(ModelUnit &pdb_entry, const vector<FragmentStructure> &templatestrall, const vector<vector<vector<bool>>> &templatestrallflag, const map<string, map<string, vector<ModelUnit>>> &ideal_rna, vector<vector<int>> &stackmap, const map<string, map<string, vector<vector<vector<double>>>>> &stackdistmap,
    const vector<vector<vector<int>>> &stembps, const map<int, vector<pair<int, int>>> &loopnts, const map<int, vector<vector<vector<int>>>> &connectnts,
    ExtParams &extparameters)
    //map<string, vector<vector<double>>> &blengthparam, map<string, vector<vector<double>>> &adjblengthparam, 
    //map<string, vector<vector<double>>> &ntangleparam, map<string, vector<vector<int>>> &ntangleparamidx, map<string, vector<vector<double>>> &adjangleparam, map<string, vector<vector<int>>> &adjangleparamidx, 
    //map<string, vector<vector<double>>> &hydparam, map<string, vector<vector<int>>> &hydparamidx, const int option=0)
{
    size_t c;
    for (c=0;c<pdb_entry.chains.size();c++)
        PreEnergyRefine(pdb_entry, templatestrall, templatestrallflag, ideal_rna, c, stackmap, stackdistmap, stembps, loopnts, connectnts, extparameters);


    return true;
}


