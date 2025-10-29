
#include "IdealRNA.hpp"

void standardize_atom_order(ModelUnit &pep)
{
    /* Define a map of maps for each residue:
    first key = residue name (string)
    second key = atom name (string)
    value = PDB standard order number (int) */

    // Define maps
    map< string, map< string, int> > NTPs;
    map< string, int> A;

    // ATP
    A[" P  "]= 0; A[" OP1"]= 1; A[" OP2"]= 2; A[" O5'"]= 3; A[" C5'"]= 4; A[" C4'"]= 5; A[" O4'"]= 6; A[" C3'"]= 7; A[" O3'"]= 8; A[" C2'"]= 9; 
    A[" O2'"]= 10; A[" C1'"]= 11; A[" N9 "]= 12; A[" C8 "]= 13; A[" N7 "]= 14; A[" C5 "]= 15; A[" C6 "]= 16; A[" N6 "]= 17; A[" N1 "]= 18; A[" C2 "]= 19;
    A[" N3 "]= 20; A[" C4 "]= 21;
    // add residue to map of NTPs
    NTPs["  A"]= A;
    
    // CTP
    map< string, int> C;
    C[" P  "]= 0; C[" OP1"]= 1; C[" OP2"]= 2; C[" O5'"]= 3; C[" C5'"]= 4; C[" C4'"]= 5; C[" O4'"]= 6; C[" C3'"]= 7; C[" O3'"]= 8; C[" C2'"]= 9;
    C[" O2'"]= 10; C[" C1'"]= 11; C[" N1 "]= 12; C[" C2 "]= 13; C[" O2 "]= 14; C[" N3 "]= 15; C[" C4 "]= 16; C[" N4 "]= 17; C[" C5 "]= 18; C[" C6 "]= 19;
    // add residue to map of NTPs
    NTPs["  C"]= C;

    // GTP
    map< string, int> G;
    G[" P  "]= 0; G[" OP1"]= 1; G[" OP2"]= 2; G[" O5'"]= 3; G[" C5'"]= 4; G[" C4'"]= 5; G[" O4'"]= 6; G[" C3'"]= 7; G[" O3'"]= 8; G[" C2'"]= 9; 
    G[" O2'"]= 10; G[" C1'"]= 11; G[" N9 "]= 12; G[" C8 "]= 13; G[" N7 "]= 14; G[" C5 "]= 15; G[" C6 "]= 16; G[" O6 "]= 17; G[" N1 "]= 18; G[" C2 "]= 19;
    G[" N2 "]= 20; G[" N3 "]= 21; G[" C4 "]= 22;
    // add residue to map of NTPs
    NTPs["  G"]= G;

    // UTP
    map< string, int> U;
    U[" P  "]= 0; U[" OP1"]= 1; U[" OP2"]= 2; U[" O5'"]= 3; U[" C5'"]= 4; U[" C4'"]= 5; U[" O4'"]= 6; U[" C3'"]= 7; U[" O3'"]=8; U[" C2'"]= 9;
    U[" O2'"]= 10; U[" C1'"]= 11; U[" N1 "]= 12; U[" C2 "]= 13; U[" O2 "]= 14; U[" N3 "]= 15; U[" C4 "]= 16; U[" O4 "]=17; U[" C5 "]= 18; U[" C6 "]= 19;
    // add residue to map of NTPs
    NTPs["  U"]= U;
    
    // Loop through chains
    for (int c=0;c<pep.chains.size();c++)
    {
        // Loop through residues in chain
        for (int r=0;r<pep.chains[c].residues.size();r++)
        {
            // Set up a temporary atom and residue
            AtomUnit temp_atom; 
            temp_atom.xyz.assign(3,0);
            ResidueUnit temp_residue;
            temp_residue.atoms.assign(23,temp_atom); // Longest residue length is 23 atoms, and each atom has x,y,z coordinates

            // Loop through atoms in each residue
            for (int a=0;a<pep.chains[c].residues[r].atoms.size();a++)
            {
                if(NTPs[pep.chains[c].residues[r].resn].count(pep.chains[c].residues[r].atoms[a].name) == 0) continue;
                // Find correct PDB order for residue-atom pair in the map of maps
                int order = NTPs[pep.chains[c].residues[r].resn][pep.chains[c].residues[r].atoms[a].name];

                // Assign correct order, coordinates, and name to each atom added to temp_residue
                temp_residue.atoms[order].xyz[0] = pep.chains[c].residues[r].atoms[a].xyz[0];
                temp_residue.atoms[order].xyz[1] = pep.chains[c].residues[r].atoms[a].xyz[1];
                temp_residue.atoms[order].xyz[2] = pep.chains[c].residues[r].atoms[a].xyz[2];
                temp_residue.atoms[order].name = pep.chains[c].residues[r].atoms[a].name;
                temp_residue.atoms[order].movable = pep.chains[c].residues[r].atoms[a].movable;
                //cout<<"pep=>temp r="<<r<<" a="<<a<<" movable="<<pep.chains[c].residues[r].atoms[a].movable<<endl;
            }
            
            // Use temp_residue to update the actual peptide with the correct order
            for (int a=0;a<pep.chains[c].residues[r].atoms.size();a++)
            {
                pep.chains[c].residues[r].atoms[a].xyz[0] = temp_residue.atoms[a].xyz[0];
                pep.chains[c].residues[r].atoms[a].xyz[1] = temp_residue.atoms[a].xyz[1];
                pep.chains[c].residues[r].atoms[a].xyz[2] = temp_residue.atoms[a].xyz[2];
                pep.chains[c].residues[r].atoms[a].name = temp_residue.atoms[a].name; 
                pep.chains[c].residues[r].atoms[a].movable = temp_residue.atoms[a].movable; 
            
                //cout<<"temp=>pep r="<<r<<" a="<<a<<" movable="<<pep.chains[c].residues[r].atoms[a].movable<<endl;
            }
    
        }
    }  
} 

vector<string> getfilelist(string path)
{
    vector<string> strs;

    strs.clear();
    for (const auto& entry : fs::directory_iterator(path)) {
        if (fs::is_regular_file(entry) && entry.path().extension() == ".pdb") {
            strs.push_back(entry.path());
        }
    }

    return strs;
}

map<string, vector<ModelUnit>> loadtemplatestrs(string basepath, vector<string> subdirs)
{
    //string basepath="../param/3nthelix_template"
    map<string, vector<ModelUnit>> template_rna;

    int i, j;
    //vector<string> subdirs = {{"A-U"}, {"G-C"}, {"G-U"}};
    vector<string> strs;

    for(i=0;i<subdirs.size();i++){
        strs=getfilelist(basepath + "/" + subdirs[i]);

        template_rna[subdirs[i]].resize(strs.size());
        for(int j=0; j<strs.size();j++){
            map<int, map<int, int>> pdbidx2seqidx;
            map<string, int> chainstr2chainidx;
            template_rna[subdirs[i]][j] = read_pdb_structure(strs[j].c_str(), chainstr2chainidx, pdbidx2seqidx, 2, 0);

            standardize_atom_order(template_rna[subdirs[i]][j]);
        }
    }
    
    vector<string>().swap(strs);

    return template_rna;
}

void calc_interaction_distmap(const ModelUnit &rna_str, vector<vector<double>> &pmat){

    double dx, dy, dz, length;
    pmat.resize(rna_str.chains[0].residues[0].atoms.size());
    for(int i=0;i<rna_str.chains[0].residues[0].atoms.size();i++){
        pmat[i].resize(rna_str.chains[0].residues[1].atoms.size());
        for(int j=0;j<rna_str.chains[0].residues[1].atoms.size();j++){
            dx = rna_str.chains[0].residues[0].atoms[i].xyz[0] - rna_str.chains[0].residues[1].atoms[j].xyz[0];
            dy = rna_str.chains[0].residues[0].atoms[i].xyz[1] - rna_str.chains[0].residues[1].atoms[j].xyz[1];
            dz = rna_str.chains[0].residues[0].atoms[i].xyz[2] - rna_str.chains[0].residues[1].atoms[j].xyz[2];
            pmat[i][j] = sqrt(dx*dx + dy*dy + dz*dz+1e-7);
        }
    }
}

void get_distmap(const map<string, vector<ModelUnit>> &rna_str, map<string, vector<vector<vector<double>>>> &stackdistmap){
    for (const auto& pair1 : rna_str) {
        stackdistmap[pair1.first].resize(pair1.second.size());
        for(int i=0;i<pair1.second.size();i++){
            calc_interaction_distmap(pair1.second[i], stackdistmap[pair1.first][i]);
        }
    }
}

map<string, map<string, vector<ModelUnit>>>parse_ideal_rna(string paramdir)
{
    map<string, map<string, vector<ModelUnit>>> ideal_rna;

    ideal_rna["bp"] = loadtemplatestrs(paramdir+"/bp_template", {"A-U", "G-C", "G-U"});
    ideal_rna["loop"] = loadtemplatestrs(paramdir+"/loop_template", {"A", "C", "G", "U"});
    ideal_rna["stacking_adj"] = loadtemplatestrs(paramdir+"/stacking_template/1", {"A_A", "A_C", "A_G", "A_U", "C_A", "C_C", "C_G", "C_U", "G_A", "G_C", "G_G", "G_U", "U_A", "U_C", "U_G", "U_U"});
    ideal_rna["stacking_nonadj"] = loadtemplatestrs(paramdir+"/stacking_template/0", {"A_A", "A_C", "A_G", "A_U", "C_C", "C_G", "C_U", "G_G", "G_U", "U_U"});
    ideal_rna["nt"] = loadtemplatestrs(paramdir+"/single_nt_template", {"A", "C", "G", "U"});

    return ideal_rna;
}

bool load_pairwiseinfo(string filepath, string res1, string res2, vector<vector<double>> &pmat){

    map<string, map<int, int>> atomidx={
        {"A", {{0, 0}, {1, 3}, {2, 4}, {3, 5}, {4, 7}, {5, 8}, {6, 9}, {7, 6}, {8, 11}, {9, 10}, {10, 1}, {11, 2}, {12, 12}, {15, 18}, {13, 21}, {16, 19}, {14, 13}, {17, 16}, {18, 14}, {19, 15}, {20, 17}, {21, 20}}},
        {"C", {{0, 0}, {1, 3}, {2, 4}, {3, 5}, {4, 7}, {5, 8}, {6, 9}, {7, 6}, {8, 11}, {9, 10}, {10, 1}, {11, 2}, {12, 12}, {15, 16}, {13, 13}, {14, 19}, {16, 18}, {17, 15}, {18, 14}, {19, 17}}},
        {"G", {{0, 0}, {1, 3}, {2, 4}, {3, 5}, {4, 7}, {5, 8}, {6, 9}, {7, 6}, {8, 11}, {9, 10}, {10, 1}, {11, 2}, {12, 12}, {15, 18}, {13, 22}, {16, 19}, {14, 13}, {17, 16}, {18, 14}, {19, 15}, {20, 21}, {21, 17}, {22, 20}}},
        {"U", {{0, 0}, {1, 3}, {2, 4}, {3, 5}, {4, 7}, {5, 8}, {6, 9}, {7, 6}, {8, 11}, {9, 10}, {10, 1}, {11, 2}, {12, 12}, {15, 16}, {13, 13}, {14, 19}, {16, 18}, {17, 15}, {18, 14}, {19, 17}}},
    };

    // CSVファイルを開く
    std::ifstream file(filepath);  // ファイル名を適切に変更してください

    // ファイルが開けなかった場合のエラーチェック
    if (!file.is_open()) {
        std::cerr << "File cannnot be opened.\n";
        return 1;
    }

    std::string line;
    int idx1=0;
    while (std::getline(file, line)) {
        std::stringstream ss(line); 
        std::string val;
        std::vector<std::string> row;
        int idx2=0;
        while (std::getline(ss, val, ',')) {
            pmat[atomidx[res1][idx1]][atomidx[res2][idx2]]=stod(val);  // セルをベクターに追加
            //cout<<idx1 << "-" << idx2 << "-" << pmat[idx1][idx2]<<endl;
            idx2++;
        }
        idx1++;
    }

    // ファイルを閉じる
    file.close();

    return true;
}

map<string, vector<DNTTEMPLATEUnit>> loadpwstrs(string basepath, vector<string> subdirs)
{
    map<string, vector<DNTTEMPLATEUnit>> template_rna;

    int i, j;
    vector<string> strs;
    string res1, res2;
    map<string, pair<int, int>> NtAtomNum={{"A_A", {22,22}}, {"A_C", {22,20}}, {"A_G", {22,23}}, {"A_U", {22,20}}, {"C_A", {20,22}}, {"C_C", {20,20}}, {"C_G", {20,23}}, {"C_U", {20,20}}, {"G_A", {23,22}}, {"G_C", {23,20}}, {"G_G", {23,23}}, {"G_U", {23,20}}, {"U_A", {20,22}}, {"U_C", {20,20}}, {"U_G", {20,23}}, {"U_U",{20,20}}};

    for(i=0;i<subdirs.size();i++){
        strs=getfilelist(basepath + "/" + subdirs[i] + "/PDB");

        res1=subdirs[i][0];
        res2=subdirs[i][2];

        template_rna[subdirs[i]].resize(strs.size());
        for(int j=0; j<strs.size();j++){
            /*
            if(subdirs[i] == "A_G" && j == 42){
                cout<<strs[j].c_str()<<endl;
                exit(1);
            }
            */
            map<int, map<int, int>> pdbidx2seqidx;
            map<string, int> chainstr2chainidx;
            template_rna[subdirs[i]][j].model = read_pdb_structure(strs[j].c_str(), chainstr2chainidx, pdbidx2seqidx, 2, 0);
            standardize_atom_order(template_rna[subdirs[i]][j].model);

            fs::path filePath = strs[j].c_str();
            string filedir=filePath.parent_path();
            string filename=filePath.stem();
            string fileext=filePath.extension();

            template_rna[subdirs[i]][j].pwdist_mean.resize(NtAtomNum[subdirs[i]].first, vector<double>(NtAtomNum[subdirs[i]].second));
            template_rna[subdirs[i]][j].pwdist_std.resize(NtAtomNum[subdirs[i]].first, vector<double>(NtAtomNum[subdirs[i]].second));
            load_pairwiseinfo(basepath + "/" + subdirs[i] + "/DISTMAP_mean/" + filename + ".txt", res1, res2, template_rna[subdirs[i]][j].pwdist_mean);
            load_pairwiseinfo(basepath + "/" + subdirs[i] + "/DISTMAP_std/" + filename + ".txt", res1, res2, template_rna[subdirs[i]][j].pwdist_std);
        }
    }
    
    vector<string>().swap(strs);

    return template_rna;
}


map<string, map<string, vector<DNTTEMPLATEUnit>>>parse_ideal_rna_distmap(string paramdir)
{
    map<string, map<string, vector<DNTTEMPLATEUnit>>> ideal_rna;

    ideal_rna["loop"] = loadpwstrs(paramdir+"/loop_template/2ntloops", {"A_A", "A_C", "A_G", "A_U", "C_A", "C_C", "C_G", "C_U", "G_A", "G_C", "G_G", "G_U", "U_A", "U_C", "U_G", "U_U"});
    
    return ideal_rna;
}

void get_score_map(string seq, const vector<pair<double,vector<size_t>>> &bp_vec, vector<double> &score_vec){
    int idx1, idx2;
    int minbpntidx = -1000000;
    int maxbpntidx = 1000000;
    
    for(int i=0;i<bp_vec.size();i++){
        idx1 = bp_vec[i].second[1];
        idx2 = bp_vec[i].second[3];

        if(idx1 < minbpntidx) minbpntidx = idx1;
        if(idx2 < minbpntidx) minbpntidx = idx2;
        if(idx1 > maxbpntidx) maxbpntidx = idx1;
        if(idx2 > maxbpntidx) maxbpntidx = idx2;
    }

    minbpntidx-=1;
    maxbpntidx+=1;

    int seq_len = seq.length();
    score_vec.resize(seq_len, 1.0);

    if(minbpntidx >= 0){
        for(int i=0;i<minbpntidx+1;i++){
            score_vec[i] = 0.1;
        }
    }

    if(maxbpntidx < seq_len){
        for(int i=maxbpntidx;i<seq_len;i++){
            score_vec[i] = 0.1;
        }
    }
}




