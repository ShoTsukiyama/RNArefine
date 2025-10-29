
#include "MissingRNAatom.hpp"

/* return true if residue is standard nucleotide
 * return false otherwise */
bool checkMissingRNAatom(ResidueUnit &residue, 
    vector<string> &missing_atom_vec,
    vector<size_t> &nonstd_atom_vec, Logger &logger, const int option)
{
    missing_atom_vec.clear();
    nonstd_atom_vec.clear();
    if (coarse_grained_atoms_pdb.count(residue.resn)==0)
    {
        //cout<<residue.resn<<residue.resi<<residue.icode<<"\tnon-standard residue"<<endl;
        return false;
    }

    size_t a,i,j;
    for (a=0;a<residue.atoms.size();a++){
        if(find(coarse_grained_atoms_pdb[residue.resn].begin(), coarse_grained_atoms_pdb[residue.resn].end(), residue.atoms[a].name) == coarse_grained_atoms_pdb[residue.resn].end())
        {
            nonstd_atom_vec.push_back(a);
        }
    }

    if (nonstd_atom_vec.size())
    {
        if(!option)
        {
            logger.log(residue.resn+to_string(residue.resi)+string(1, residue.icode)+"\t"+to_string(nonstd_atom_vec.size())+" non-standard atoms\t");
            for (i=0;i<nonstd_atom_vec.size();i++) logger.log(residue.atoms[nonstd_atom_vec[i]].name+' ');
            logger.log("\n");
        }
        if (option)
        {
            for (i=0;i<nonstd_atom_vec.size();i++)
                residue.atoms[nonstd_atom_vec[i]].name.clear();
            for (i=0;i<residue.atoms.size();i++)
            {
                if (residue.atoms[i].name.size()) continue;
                for (j=i+1;j<residue.atoms.size();j++)
                {
                    if (residue.atoms[j].name.size()==0) continue;
                    residue.atoms[i].name  =residue.atoms[j].name;
                    residue.atoms[i].xyz[0]=residue.atoms[j].xyz[0];
                    residue.atoms[i].xyz[1]=residue.atoms[j].xyz[1];
                    residue.atoms[i].xyz[2]=residue.atoms[j].xyz[2];
                    residue.atoms[j].name.clear();
                    break;
                }
            }
            for (i=0;i<nonstd_atom_vec.size();i++) residue.atoms.pop_back();
        }
    }
    
    string atom;
    bool found;
    for (i=0; i!=coarse_grained_atoms_pdb[residue.resn].size(); i++)
    {
        atom=coarse_grained_atoms_pdb[residue.resn][i];
        found=false;
        for (a=0;a<residue.atoms.size();a++)
        {
            if (residue.atoms[a].name==atom)
            {
                found=true;
                break;
            }
        }
        if (found==false) missing_atom_vec.push_back(atom);
    }
    // terminal output that lists missing atoms
    if (missing_atom_vec.size())
    {
        logger.log(residue.resn+to_string(residue.resi)+string(1, residue.icode)+"\t"+to_string(missing_atom_vec.size())+" missing atoms\t");
        for (a=0;a<missing_atom_vec.size();a++)
        logger.log(missing_atom_vec[a]+' ');
        logger.log("\n");
    }

    return true;
}

size_t countUniqAtom(vector<vector<double> >&xyz_list2)
{
    size_t atomNum=xyz_list2.size();
    double Extra=1.0e-6;
    size_t dup;
    size_t a,b;
    for (a=1;a<xyz_list2.size();a++)
    {
        dup=0;
        for (b=0;b<a;b++){
            if (fabs(xyz_list2[a][0]-xyz_list2[b][0])+
                fabs(xyz_list2[a][1]-xyz_list2[b][1])+
                fabs(xyz_list2[a][2]-xyz_list2[b][2])<Extra) dup++;
        }
        if (dup) atomNum--;
    }
    return atomNum;
}


//bool fillMissingRNAatomInBasePair(const int c1, const int r1, const int c2, const int r2, 
    //const ModelUnit &pdb_entry, ChainUnit &fill_chain,
    //map<string, vector<ModelUnit>> &ideal_rna,
    //vector<string>missing_atom_vec_r1, vector<string>missing_atom_vec_r2, int &templateidx, bool fullatom,
    //map<string, vector<vector<double>>> &blengthparam, map<string, vector<vector<double>>> &adjblengthparam,
    //map<string, vector<vector<double>>> &ntangleparam, map<string, vector<vector<int>>> &ntangleparamidx, map<string, vector<vector<double>>> &adjangleparam, map<string, vector<vector<int>>> &adjangleparamidx)
bool AssembleBasePairFrag(const int c1, const int r1, const int c2, const int r2, 
    const ModelUnit &pdb_entry, map<string, vector<ModelUnit>> &ideal_rna,
    FragmentResidue &spresidue1, FragmentResidue &spresidue2, 
    int templateidx)
{
    int i;
    vector<double> tmp(3,0);
    vector<vector<double> > xyz_list1;
    vector<vector<double> > xyz_list2;
    vector<vector<double> > RotMatix, minRotMatix;  // U
    vector<double> TranVect, minTranVect;  // t
    string key, key1, key2;
    string base1, base2, base1_pdb, base2_pdb, bpkey;
    vector<vector<bool>> movableflag;
    vector<vector<int>> movableatomid;
    double rmsd, minrmsd;
    int mintempidx;
    int atom;
    int bp1, bp2;
    double startene, endene;

    base1=pdb_entry.chains[c1].residues[r1].resn[2];
    base1_pdb=pdb_entry.chains[c1].residues[r1].resn;
    base2=pdb_entry.chains[c2].residues[r2].resn[2];
    base2_pdb=pdb_entry.chains[c2].residues[r2].resn;
    bpkey=base1+"-"+base2;

    int idx1, idx2;
    if(bpkey == "A-U" || bpkey == "G-C" || bpkey == "G-U"){
        idx1=0;
        idx2=1;
    }else if(bpkey == "U-A" || bpkey == "C-G" || bpkey == "U-G"){
        idx1=1;
        idx2=0;
        bpkey=base2+"-"+base1;
    }else{
        cout<<bpkey<<endl;
        cout<<"Error"<<endl;
        exit(1);
    }

    size_t a,a1,a2;
    size_t atomNum=0;

    movableflag.resize(2, vector<bool>(23,false));
    movableatomid.resize(2, vector<int>(23,-1));
    for (a=0;a<pdb_entry.chains[c1].residues[r1].atoms.size();a++){
        auto it = find(coarse_grained_atoms_pdb[base1_pdb].begin(), coarse_grained_atoms_pdb[base1_pdb].end(), pdb_entry.chains[c1].residues[r1].atoms[a].name);
        
        if(pdb_entry.chains[c1].residues[r1].atoms[a].movable==0 && (it != coarse_grained_atoms_pdb[base1_pdb].end())){
            atomNum+=1;
            movableflag[0][a]=true;
        }
        movableatomid[0][a]=distance(coarse_grained_atoms_pdb[base1_pdb].begin(), it);
    }

    for (a=0;a<pdb_entry.chains[c2].residues[r2].atoms.size();a++){
        auto it = find(coarse_grained_atoms_pdb[base2_pdb].begin(), coarse_grained_atoms_pdb[base2_pdb].end(), pdb_entry.chains[c2].residues[r2].atoms[a].name);

        if(pdb_entry.chains[c2].residues[r2].atoms[a].movable==0 && (it != coarse_grained_atoms_pdb[base2_pdb].end())){
            atomNum+=1;
            movableflag[1][a]=true;
        }
        movableatomid[1][a]=distance(coarse_grained_atoms_pdb[base2_pdb].begin(), it);
    }

    if (ideal_rna.count(bpkey) == 0){
        for (const auto& pair : ideal_rna) cout<<pair.first << "-" <<endl;
        exit(1);
        return false;
    }

    minrmsd=10000000000; mintempidx=-1;
    for(i=0;i<ideal_rna[bpkey].size();i++)
    {
        xyz_list1.clear();
        xyz_list2.clear();

        xyz_list1.assign(atomNum,tmp);
        xyz_list2.assign(atomNum,tmp);
        a=0;

        for (a1=0;a1<pdb_entry.chains[c1].residues[r1].atoms.size();a1++)
        {
            //if (pdb_entry.chains[c1].residues[r1].atoms[a1].movable) continue;
            if (!movableflag[0][a1]) continue;

            xyz_list1[a][0]=ideal_rna[bpkey][i].chains[0].residues[idx1].atoms[movableatomid[0][a1]].xyz[0];
            xyz_list1[a][1]=ideal_rna[bpkey][i].chains[0].residues[idx1].atoms[movableatomid[0][a1]].xyz[1];
            xyz_list1[a][2]=ideal_rna[bpkey][i].chains[0].residues[idx1].atoms[movableatomid[0][a1]].xyz[2];

            xyz_list2[a][0]=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz[0];
            xyz_list2[a][1]=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz[1];
            xyz_list2[a][2]=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz[2];

            a++;
        }

        for (a2=0;a2<pdb_entry.chains[c2].residues[r2].atoms.size();a2++)
        {
            if (!movableflag[1][a2]) continue;

            xyz_list1[a][0]=ideal_rna[bpkey][i].chains[0].residues[idx2].atoms[movableatomid[1][a2]].xyz[0];
            xyz_list1[a][1]=ideal_rna[bpkey][i].chains[0].residues[idx2].atoms[movableatomid[1][a2]].xyz[1];
            xyz_list1[a][2]=ideal_rna[bpkey][i].chains[0].residues[idx2].atoms[movableatomid[1][a2]].xyz[2];

            xyz_list2[a][0]=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz[0];
            xyz_list2[a][1]=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz[1];
            xyz_list2[a][2]=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz[2];

            a++;
        }

        if(i == 0){
            atomNum=countUniqAtom(xyz_list2);

            if (atomNum<xyz_list2.size())
            cerr<<"duplicated coordinate when reconstructing "
                <<pdb_entry.chains[c1].residues[r1].resn
                <<pdb_entry.chains[c1].residues[r1].resi
                <<pdb_entry.chains[c1].residues[r1].icode<<" and "
                <<pdb_entry.chains[c2].residues[r2].resn
                <<pdb_entry.chains[c2].residues[r2].resi
                <<pdb_entry.chains[c2].residues[r2].icode<<endl;

            if (atomNum<3)
            {
                for (a=0;a<xyz_list1.size();a++)
                {
                    xyz_list1[a].clear();
                    xyz_list2[a].clear();
                }
                xyz_list1.clear();
                xyz_list2.clear();
                //key.clear();
                //key1.clear();
                //key2.clear();
                return false;
            }
        }
        
        RotMatix.clear();
        TranVect.clear();
        RotateCoor(xyz_list1,xyz_list2, RotMatix, TranVect);

        atomNum=0;
        rmsd=0;
        for (a=0;a<pdb_entry.chains[c1].residues[r1].atoms.size();a++)
        {
            if (!movableflag[0][a]) continue;
            ChangeCoor(ideal_rna[bpkey][i].chains[0].residues[idx1].atoms[movableatomid[0][a]].xyz,RotMatix, TranVect, xyz_list1[atomNum]);
            tmp[0]=(pdb_entry.chains[c1].residues[r1].atoms[a].xyz[0]-xyz_list1[atomNum][0]);
            tmp[1]=(pdb_entry.chains[c1].residues[r1].atoms[a].xyz[1]-xyz_list1[atomNum][1]);
            tmp[2]=(pdb_entry.chains[c1].residues[r1].atoms[a].xyz[2]-xyz_list1[atomNum][2]);

            atomNum++;
            rmsd+=sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
        }
        for (a=0;a<pdb_entry.chains[c2].residues[r2].atoms.size();a++)
        {
            if (!movableflag[1][a]) continue;
            ChangeCoor(ideal_rna[bpkey][i].chains[0].residues[idx2].atoms[movableatomid[1][a]].xyz,RotMatix,TranVect,xyz_list1[atomNum]);
            tmp[0]=(pdb_entry.chains[c2].residues[r2].atoms[a].xyz[0]-xyz_list1[atomNum][0]);
            tmp[1]=(pdb_entry.chains[c2].residues[r2].atoms[a].xyz[1]-xyz_list1[atomNum][1]);
            tmp[2]=(pdb_entry.chains[c2].residues[r2].atoms[a].xyz[2]-xyz_list1[atomNum][2]);

            atomNum++;
            rmsd+=sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
        }

        if(rmsd < minrmsd){
            minrmsd=rmsd;
            mintempidx=i;
            minRotMatix=RotMatix;
            minTranVect=TranVect;
        }
    }

    templateidx=mintempidx;

    spresidue1.restype = base1;
    spresidue2.restype = base2;

    FragmentAtoms coordtmp;
    coordtmp.coord.x = 0; coordtmp.coord.y = 0; coordtmp.coord.z = 0;
    spresidue1.atoms.resize(coarse_grained_atoms[base1].size(), coordtmp);
    spresidue2.atoms.resize(coarse_grained_atoms[base2].size(), coordtmp);

    for(i=0;i<coarse_grained_atoms[base1].size();i++){
        ChangeCoor(ideal_rna[bpkey][mintempidx].chains[0].residues[idx1].atoms[i].xyz,minRotMatix,minTranVect,spresidue1.atoms[i].coord);
    }
    for(i=0;i<coarse_grained_atoms[base2].size();i++){
        ChangeCoor(ideal_rna[bpkey][mintempidx].chains[0].residues[idx2].atoms[i].xyz,minRotMatix,minTranVect,spresidue2.atoms[i].coord);
    }
        
    /* clean up */
    for (a=0;a<xyz_list1.size();a++)
    {
        xyz_list1[a].clear();
        xyz_list2[a].clear();
    }
    xyz_list1.clear();
    xyz_list2.clear();
    TranVect.clear();
    for (a=0;a<3;a++){
        RotMatix[a].clear();
    }
    RotMatix.clear();
    movableflag.clear();
    movableatomid.clear();
    return true;
}

void UpdateStemStructure(const int c1, const int r1, const int c2, const int r2, 
    const ModelUnit &pdb_entry, ChainUnit &fill_chain,
    vector<string>missing_atom_vec_r1, vector<string>missing_atom_vec_r2,
    FragmentResidue &spresidue1, FragmentResidue &spresidue2, 
    bool fullatom)
{
    int a, a1, a2;
    string base1, base2, base1_pdb, base2_pdb;
    int atom;

    base1=pdb_entry.chains[c1].residues[r1].resn[2];
    base2=pdb_entry.chains[c2].residues[r2].resn[2];
    base1_pdb=pdb_entry.chains[c1].residues[r1].resn;
    base2_pdb=pdb_entry.chains[c2].residues[r2].resn;

    if(fullatom){
        for (a=0;a<pdb_entry.chains[c1].residues[r1].atoms.size();a++)
        {
            auto it = find(full_atoms_pdb[base1_pdb].begin(), full_atoms_pdb[base1_pdb].end(), pdb_entry.chains[c1].residues[r1].atoms[a].name);

            if(it != full_atoms_pdb[base1_pdb].end()){
                atom = distance(full_atoms_pdb[base1_pdb].begin(), it);
                fill_chain.residues[r1].atoms[a].xyz[0]=spresidue1.atoms[atom].coord.x;
                fill_chain.residues[r1].atoms[a].xyz[1]=spresidue1.atoms[atom].coord.y;
                fill_chain.residues[r1].atoms[a].xyz[2]=spresidue1.atoms[atom].coord.z;
            }
        }

        a=0;
        for (a1=pdb_entry.chains[c1].residues[r1].atoms.size();a1<pdb_entry.chains[c1].residues[r1].atoms.size()+missing_atom_vec_r1.size();a1++){
            auto it = find(full_atoms_pdb[base1_pdb].begin(), full_atoms_pdb[base1_pdb].end(), missing_atom_vec_r1[a]);

            if(it != full_atoms_pdb[base1_pdb].end()){
                atom = distance(full_atoms_pdb[base1_pdb].begin(), it);
                fill_chain.residues[r1].atoms[a1].xyz[0]=spresidue1.atoms[atom].coord.x;
                fill_chain.residues[r1].atoms[a1].xyz[1]=spresidue1.atoms[atom].coord.y;
                fill_chain.residues[r1].atoms[a1].xyz[2]=spresidue1.atoms[atom].coord.z;
            }
            a+=1;
        }


        for (a=0;a<pdb_entry.chains[c2].residues[r2].atoms.size();a++)
        {
            auto it = find(full_atoms_pdb[base2_pdb].begin(), full_atoms_pdb[base2_pdb].end(), pdb_entry.chains[c2].residues[r2].atoms[a].name);

            if(it != full_atoms_pdb[base2_pdb].end()){
                atom = distance(full_atoms_pdb[base2_pdb].begin(), it);
                fill_chain.residues[r2].atoms[a].xyz[0]=spresidue2.atoms[atom].coord.x;
                fill_chain.residues[r2].atoms[a].xyz[1]=spresidue2.atoms[atom].coord.y;
                fill_chain.residues[r2].atoms[a].xyz[2]=spresidue2.atoms[atom].coord.z;
            }
        }
        a=0;
        for (a2=pdb_entry.chains[c2].residues[r2].atoms.size();a2<pdb_entry.chains[c2].residues[r2].atoms.size()+missing_atom_vec_r2.size();a2++){
            auto it = find(full_atoms_pdb[base2_pdb].begin(), full_atoms_pdb[base2_pdb].end(), missing_atom_vec_r2[a]);

            if(it != full_atoms_pdb[base2_pdb].end()){
                atom = distance(full_atoms_pdb[base2_pdb].begin(), it);
                fill_chain.residues[r2].atoms[a2].xyz[0]=spresidue2.atoms[atom].coord.x;
                fill_chain.residues[r2].atoms[a2].xyz[1]=spresidue2.atoms[atom].coord.y;
                fill_chain.residues[r2].atoms[a2].xyz[2]=spresidue2.atoms[atom].coord.z;
            }
            a+=1;
        }
    }else{
        for (a=0;a<missing_atom_vec_r1.size();a++){
            auto it = find(full_atoms_pdb[base1_pdb].begin(), full_atoms_pdb[base1_pdb].end(), missing_atom_vec_r1[a]);

            if(it != full_atoms_pdb[base1_pdb].end()){
                atom = distance(full_atoms_pdb[base1_pdb].begin(), it);
                fill_chain.residues[r1].atoms[a].xyz[0]=spresidue1.atoms[atom].coord.x;
                fill_chain.residues[r1].atoms[a].xyz[1]=spresidue1.atoms[atom].coord.y;
                fill_chain.residues[r1].atoms[a].xyz[2]=spresidue1.atoms[atom].coord.z;
            }
        }

        for (a=0;a<missing_atom_vec_r2.size();a++){
            auto it = find(full_atoms_pdb[base2_pdb].begin(), full_atoms_pdb[base2_pdb].end(), missing_atom_vec_r2[a]);
            
            if(it != full_atoms_pdb[base2_pdb].end()){
                atom = distance(full_atoms_pdb[base2_pdb].begin(), it);
                fill_chain.residues[r2].atoms[a].xyz[0]=spresidue2.atoms[atom].coord.x;
                fill_chain.residues[r2].atoms[a].xyz[1]=spresidue2.atoms[atom].coord.y;
                fill_chain.residues[r2].atoms[a].xyz[2]=spresidue2.atoms[atom].coord.z;
            }
        }
    }
}

void ExtractNt(const int c, const int r, const ModelUnit &pdb_entry, FragmentResidue &spresidue)
{
    int a;
    int atom;
    string base, base_pdb;

    base=pdb_entry.chains[c].residues[r].resn[2];
    base_pdb=pdb_entry.chains[c].residues[r].resn;

    FragmentAtoms coordtmp;
    coordtmp.coord.x = 0; coordtmp.coord.y = 0; coordtmp.coord.z = 0;
    spresidue.restype = base;
    //initresidue.restype = base;
    if(full_atoms_pdb[base_pdb].size() == pdb_entry.chains[c].residues[r].atoms.size()){
        spresidue.atoms.resize(full_atoms_pdb[base_pdb].size(), coordtmp);
    }else{
        spresidue.atoms.resize(coarse_grained_atoms_pdb[base_pdb].size(), coordtmp);
    }
    //initresidue.atoms.resize(full_atoms[base].size(), coordtmp);
    //inputtedstrflag.resize(full_atoms[base].size(), false);
    for (a=0;a<pdb_entry.chains[c].residues[r].atoms.size();a++){
        auto it = find(full_atoms_pdb[base_pdb].begin(), full_atoms_pdb[base_pdb].end(), pdb_entry.chains[c].residues[r].atoms[a].name);
        if(it != full_atoms_pdb[base_pdb].end()){
            atom = distance(full_atoms_pdb[base_pdb].begin(), it);

            //cout<<atom << "-" << a  << "-" << full_atoms_pdb[base_pdb][atom] << "-" << pdb_entry.chains[c].residues[r].atoms[a].name << endl;

            spresidue.atoms[atom].coord.x = pdb_entry.chains[c].residues[r].atoms[a].xyz[0];
            spresidue.atoms[atom].coord.y = pdb_entry.chains[c].residues[r].atoms[a].xyz[1];
            spresidue.atoms[atom].coord.z = pdb_entry.chains[c].residues[r].atoms[a].xyz[2];

            //initresidue.atoms[atom].coord.x = pdb_entry.chains[c].residues[r].atoms[a].xyz[0];
            //initresidue.atoms[atom].coord.y = pdb_entry.chains[c].residues[r].atoms[a].xyz[1];
            //initresidue.atoms[atom].coord.z = pdb_entry.chains[c].residues[r].atoms[a].xyz[2];
            //inputtedstrflag[atom] = true;
        }
    }
}

void GetInitStr(int c, const ModelUnit &pdb_entry, FragmentStructure &initresidues, vector<vector<bool>> &inputtedstrflags)
{
    int r, a;
    int atom;
    string base, base_pdb;
    int reslen = pdb_entry.chains[c].residues.size();

    for(r = 0; r < reslen;r++){
        base_pdb=pdb_entry.chains[c].residues[r].resn;
        base=pdb_entry.chains[c].residues[r].resn[2];

        FragmentResidue initresidue;
        vector<bool> inputtedstrflag;
        FragmentAtoms coordtmp;
        coordtmp.coord.x = 0; coordtmp.coord.y = 0; coordtmp.coord.z = 0;

        initresidue.restype = base;
        initresidue.atoms.resize(full_atoms[base].size(), coordtmp);
        inputtedstrflag.resize(full_atoms[base].size(), false);

        for (a=0;a<pdb_entry.chains[c].residues[r].atoms.size();a++){
            auto it = find(full_atoms_pdb[base_pdb].begin(), full_atoms_pdb[base_pdb].end(), pdb_entry.chains[c].residues[r].atoms[a].name);
            atom = distance(full_atoms_pdb[base_pdb].begin(), it);
            initresidue.atoms[atom].coord.x = pdb_entry.chains[c].residues[r].atoms[a].xyz[0];
            initresidue.atoms[atom].coord.y = pdb_entry.chains[c].residues[r].atoms[a].xyz[1];
            initresidue.atoms[atom].coord.z = pdb_entry.chains[c].residues[r].atoms[a].xyz[2];
            inputtedstrflag[atom] = true;
        }
        initresidues.residues.push_back(initresidue);
        inputtedstrflags.push_back(inputtedstrflag);
    }  
}


//bool fillMissingRNAatomLocal(const int c, const int r,
    //const ModelUnit &pdb_entry, ChainUnit &fill_chain, map<string, vector<ModelUnit>> &ideal_rna, vector<string>missing_atom_vec, int &templateidx, bool fullatom,
    //map<string, vector<vector<double>>> &blengthparam, map<string, vector<vector<double>>> &adjblengthparam,
    //map<string, vector<vector<double>>> &ntangleparam, map<string, vector<vector<int>>> &ntangleparamidx, map<string, vector<vector<double>>> &adjangleparam, map<string, vector<vector<int>>> &adjangleparamidx)
bool AssembleNonBpNtFrag(const int c, const int r, const ModelUnit &pdb_entry, map<string, vector<ModelUnit>> &ideal_rna, FragmentResidue &spresidue, int templateidx)
{
    int i;
    vector<double> tmp(3,0);
    vector<vector<double> > xyz_list1;
    vector<vector<double> > xyz_list2;
    vector<vector<double> > RotMatix, minRotMatix;  // U
    vector<double> TranVect, minTranVect;  // t
    string key, key1, key2;
    bool stack_prev, stack_next;
    string base;
    vector<vector<bool>> movableflag;
    vector<vector<int>> movableatomid;
    double rmsd, minrmsd;
    int mintempidx;
    string base2, base2_pdb;
    int atom;
    int nt;
    double startene, endene;

    base2=pdb_entry.chains[c].residues[r].resn[2];
    base2_pdb=pdb_entry.chains[c].residues[r].resn;

    size_t a,a1,a2;
    size_t atomNum=0;

    movableflag.resize(1, vector<bool>(23,false));
    movableatomid.resize(1, vector<int>(23,-1));
    for (a=0;a<pdb_entry.chains[c].residues[r].atoms.size();a++){
        auto it = find(coarse_grained_atoms_pdb[base2_pdb].begin(), coarse_grained_atoms_pdb[base2_pdb].end(), pdb_entry.chains[c].residues[r].atoms[a].name);
        if(pdb_entry.chains[c].residues[r].atoms[a].movable==0 && (it != coarse_grained_atoms_pdb[base2_pdb].end())){
            atomNum+=1;
            movableflag[0][a]=true;
        }
        movableatomid[0][a]=distance(coarse_grained_atoms_pdb[base2_pdb].begin(), it);
    }

    if (ideal_rna.count(base2) == 0) return false;

    minrmsd=10000000000; mintempidx=-1;
    for(i=0;i<ideal_rna[base2].size();i++)
    {
        xyz_list1.clear();
        xyz_list2.clear();

        xyz_list1.assign(atomNum,tmp);
        xyz_list2.assign(atomNum,tmp);
        a=0;

        for (a1=0;a1<pdb_entry.chains[c].residues[r].atoms.size();a1++)
        {
            //if (pdb_entry.chains[c].residues[r].atoms[a1].movable) continue;
            if (!movableflag[0][a1]) continue;
            xyz_list1[a][0]=ideal_rna[base2][i].chains[0].residues[0].atoms[movableatomid[0][a1]].xyz[0];
            xyz_list1[a][1]=ideal_rna[base2][i].chains[0].residues[0].atoms[movableatomid[0][a1]].xyz[1];
            xyz_list1[a][2]=ideal_rna[base2][i].chains[0].residues[0].atoms[movableatomid[0][a1]].xyz[2];

            xyz_list2[a][0]=pdb_entry.chains[c].residues[r].atoms[a1].xyz[0];
            xyz_list2[a][1]=pdb_entry.chains[c].residues[r].atoms[a1].xyz[1];
            xyz_list2[a][2]=pdb_entry.chains[c].residues[r].atoms[a1].xyz[2];
            a++;
        }

        if(i == 0){
            atomNum=countUniqAtom(xyz_list2);
            if (atomNum<xyz_list2.size())
            cerr<<"duplicated coordinate when reconstructing "
                <<pdb_entry.chains[c].residues[r].resn
                <<pdb_entry.chains[c].residues[r].resi
                <<pdb_entry.chains[c].residues[r].icode<<endl;

            if (atomNum<3)
            {
                for (a=0;a<xyz_list1.size();a++)
                {
                    xyz_list1[a].clear();
                    xyz_list2[a].clear();
                }
                xyz_list1.clear();
                xyz_list2.clear();
                key.clear();
                key1.clear();
                key2.clear();
                return false;
            }
        }

        RotMatix.clear();
        TranVect.clear();
        RotateCoor(xyz_list1,xyz_list2, RotMatix, TranVect);

        atomNum=0;
        rmsd=0;
        for (a=0;a<pdb_entry.chains[c].residues[r].atoms.size();a++)
        {
            if (!movableflag[0][a]) continue;

            ChangeCoor(ideal_rna[base2][i].chains[0].residues[0].atoms[movableatomid[0][a]].xyz,RotMatix, TranVect, xyz_list1[atomNum]);

            tmp[0]=(pdb_entry.chains[c].residues[r].atoms[a].xyz[0]-xyz_list1[atomNum][0]);
            tmp[1]=(pdb_entry.chains[c].residues[r].atoms[a].xyz[1]-xyz_list1[atomNum][1]);
            tmp[2]=(pdb_entry.chains[c].residues[r].atoms[a].xyz[2]-xyz_list1[atomNum][2]);
            atomNum++;
            rmsd+=sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
        }

        if(rmsd < minrmsd){
            minrmsd=rmsd;
            mintempidx=i;
            minRotMatix=RotMatix;
            minTranVect=TranVect;
        }
    }

    templateidx=mintempidx;

    spresidue.restype = base2;

    FragmentAtoms coordtmp;
    coordtmp.coord.x = 0; coordtmp.coord.y = 0; coordtmp.coord.z = 0;
    spresidue.atoms.resize(coarse_grained_atoms_pdb[base2_pdb].size(), coordtmp);
    
    for(i=0;i<coarse_grained_atoms_pdb[base2_pdb].size();i++){
        ChangeCoor(ideal_rna[base2][mintempidx].chains[0].residues[0].atoms[i].xyz,minRotMatix,minTranVect,spresidue.atoms[i].coord);
    }
        
    for (a=0;a<xyz_list1.size();a++)
    {
        xyz_list1[a].clear();
        xyz_list2[a].clear();
    }
    xyz_list1.clear();
    xyz_list2.clear();
    TranVect.clear();
    for (a=0;a<3;a++){
        RotMatix[a].clear();
        minRotMatix[a].clear();
    }
    RotMatix.clear();
    minRotMatix.clear();
    key.clear();
    key1.clear();
    key2.clear();
    movableflag.clear();
    movableatomid.clear();
    return true;
}

void UpdateLoopStructure(const int c, const int r,
    const ModelUnit &pdb_entry, ChainUnit &fill_chain,
    vector<string>missing_atom_vec, FragmentResidue &spresidue,
    bool fullatom)
{
    int a, a1, a2;
    string base, base_pdb;
    int atom;

    base=pdb_entry.chains[c].residues[r].resn[2];
    base_pdb=pdb_entry.chains[c].residues[r].resn;

    if(fullatom){
        for (a=0;a<pdb_entry.chains[c].residues[r].atoms.size();a++)
        {
            auto it = find(full_atoms_pdb[base_pdb].begin(), full_atoms_pdb[base_pdb].end(), pdb_entry.chains[c].residues[r].atoms[a].name);

            if(it != full_atoms_pdb[base_pdb].end()){
                atom = distance(full_atoms_pdb[base_pdb].begin(), it);

                fill_chain.residues[r].atoms[a].xyz[0]=spresidue.atoms[atom].coord.x;
                fill_chain.residues[r].atoms[a].xyz[1]=spresidue.atoms[atom].coord.y;
                fill_chain.residues[r].atoms[a].xyz[2]=spresidue.atoms[atom].coord.z;
            }
        }
        a=0;
        for (a1=pdb_entry.chains[c].residues[r].atoms.size();a1<pdb_entry.chains[c].residues[r].atoms.size()+missing_atom_vec.size();a1++){
            auto it = find(full_atoms_pdb[base_pdb].begin(), full_atoms_pdb[base_pdb].end(), missing_atom_vec[a]);
            
            if(it != full_atoms_pdb[base_pdb].end()){
                atom = distance(full_atoms_pdb[base_pdb].begin(), it);

                fill_chain.residues[r].atoms[a1].xyz[0]=spresidue.atoms[atom].coord.x;
                fill_chain.residues[r].atoms[a1].xyz[1]=spresidue.atoms[atom].coord.y;
                fill_chain.residues[r].atoms[a1].xyz[2]=spresidue.atoms[atom].coord.z;
            }

            a+=1;
        }
    }else{
        for (a=0;a<missing_atom_vec.size();a++){
            auto it = find(full_atoms_pdb[base_pdb].begin(), full_atoms_pdb[base_pdb].end(), missing_atom_vec[a]);
            if(it != full_atoms_pdb[base_pdb].end()){
                atom = distance(full_atoms_pdb[base_pdb].begin(), it);

                fill_chain.residues[r].atoms[a].xyz[0]=spresidue.atoms[atom].coord.x;
                fill_chain.residues[r].atoms[a].xyz[1]=spresidue.atoms[atom].coord.y;
                fill_chain.residues[r].atoms[a].xyz[2]=spresidue.atoms[atom].coord.z;
            }
        }
    }
}

void AddMissingAtoms(ModelUnit &pdb_entry, int c, const ChainUnit &fill_chain, const map<size_t,vector<string> > &fill_residue_map, const vector<bool> &updateflag, int option)
{
    int r, a;
    if(option == 1){
        for (r=0;r<pdb_entry.chains[c].residues.size();r++)
        {
            if(!updateflag[r]) continue;
            for (a=0;a<fill_chain.residues[r].atoms.size();a++)
            {
                pdb_entry.chains[c].residues[r].atoms.push_back(fill_chain.residues[r].atoms[a]);
            }
        }
    }else if(option == 2){
        for (r=0;r<pdb_entry.chains[c].residues.size();r++)
        {
            if(!updateflag[r]) continue;
            for (a=0;a<pdb_entry.chains[c].residues[r].atoms.size();a++)
            {
                pdb_entry.chains[c].residues[r].atoms[a].xyz = fill_chain.residues[r].atoms[a].xyz;
            }
    
            int inputatomnum = pdb_entry.chains[c].residues[r].atoms.size();
            //for (a=pdb_entry.chains[c].residues[r].atoms.size();a<pdb_entry.chains[c].residues[r].atoms.size()+fill_residue_map[r].size();a++)
            for(a=inputatomnum;a<inputatomnum+fill_residue_map.at(r).size();a++)
            {
                pdb_entry.chains[c].residues[r].atoms.push_back(fill_chain.residues[r].atoms[a]);
            }
        }
    }
}

bool MissingRNAatom(ModelUnit &pdb_entry, const size_t c,
    map<string, map<string, vector<ModelUnit>>> &ideal_rna,
    const vector<vector<vector<int>>> &stembps, const map<int, vector<pair<int, int>>> &loopnts, const map<int, vector<vector<vector<int>>>> &connectnts,
    map<int, int> &bp_template_idx, map<int, int> &nonbp_template_idx,
    FragmentStructure &inputtedstrall, vector<vector<bool>> &inputtedstrallflag,
    ExtParams &extparameters, Logger &logger, const int option)
    //map<string, vector<vector<double>>> &blengthparam, map<string, vector<vector<double>>> &adjblengthparam, 
    //map<string, vector<vector<double>>> &ntangleparam, map<string, vector<vector<int>>> &ntangleparamidx, map<string, vector<vector<double>>> &adjangleparam, map<string, vector<vector<int>>> &adjangleparamidx, 
    //map<string, vector<vector<double>>> &hydparam, map<string, vector<vector<int>>> &hydparamidx, const int option=0)
{
    vector<string> missing_atom_vec;
    vector<size_t> nonstd_atom_vec;
    size_t r,i,j;
    vector<size_t> remove_residue_vec;
    map<size_t,vector<string> >fill_residue_map;
    vector<int> nonbpresidues;
    bool fullatom;
    int templateidx;
    bool bpprev, bpnext;
    string looptype;

    int missing_atom_num=0;
    for (r=0;r<pdb_entry.chains[c].residues.size();r++)
    {
        nonbpresidues.push_back(r);
        checkMissingRNAatom(pdb_entry.chains[c].residues[r], missing_atom_vec, nonstd_atom_vec, logger, option);

        if (pdb_entry.chains[c].residues[r].atoms.size()==0)
        {
            pdb_entry.chains[c].residues[r].resn.clear();
            remove_residue_vec.push_back(r);
        }
        else if (missing_atom_vec.size()) fill_residue_map[r]=missing_atom_vec;
        missing_atom_num+=missing_atom_vec.size();
    }

    if (remove_residue_vec.size())
    {
        for (i=0;i<pdb_entry.chains[c].residues.size();i++)
        {
            if (pdb_entry.chains[c].residues[i].resn.size()) continue;
            for (j=i+1;j<pdb_entry.chains[c].residues.size();j++)
            {
                if (pdb_entry.chains[c].residues[j].resn.size()==0) continue;
                pdb_entry.chains[c].residues[i].het  =pdb_entry.chains[c].residues[j].het;
                pdb_entry.chains[c].residues[i].resi =pdb_entry.chains[c].residues[j].resi;
                pdb_entry.chains[c].residues[i].icode=pdb_entry.chains[c].residues[j].icode;
                pdb_entry.chains[c].residues[i].resn =pdb_entry.chains[c].residues[j].resn;
                pdb_entry.chains[c].residues[i].atoms=pdb_entry.chains[c].residues[j].atoms;
                pdb_entry.chains[c].residues[j].resn.clear();
                break;
            }
        }
        for (i=0;i<remove_residue_vec.size();i++) pdb_entry.chains[c].residues.pop_back();
    }

    size_t bp;
    int r1,r2,c1,c2;

    ChainUnit fill_chain;
    ResidueUnit fill_residue;
    AtomUnit fill_atom;
    size_t a;
    vector<bool> updateflag(pdb_entry.chains[c].residues.size(), false);

    if(option > 0){
        if(missing_atom_num){
            cout<<"Number of missing atoms: "<< missing_atom_num<<endl;
            /*
            if(option == 1)
            {
                for (r=0;r<pdb_entry.chains[c].residues.size();r++)
                {
                    fill_chain.residues.push_back(fill_residue);
                    if (fill_residue_map.count(r)==0) continue;
                    for (a=0;a<fill_residue_map[r].size();a++)
                    {
                        fill_chain.residues[r].atoms.push_back(fill_atom);
                        fill_chain.residues[r].atoms[a].name=fill_residue_map[r][a];
                        fill_chain.residues[r].atoms[a].xyz.assign(3,0);
                        fill_chain.residues[r].atoms[a].movable=1;
                        fill_chain.residues[r].atoms[a].input=1;
                    }
                }
                fullatom=false;
            }else if(option == 2)
            {   
            */ 
            for (r=0;r<pdb_entry.chains[c].residues.size();r++)
            {
                fill_chain.residues.push_back(fill_residue);

                for (a=0;a<pdb_entry.chains[c].residues[r].atoms.size();a++)
                {
                    fill_chain.residues[r].atoms.push_back(fill_atom);
                    fill_chain.residues[r].atoms[a].name=pdb_entry.chains[c].residues[r].atoms[a].name;
                    fill_chain.residues[r].atoms[a].xyz.assign(3,0);
                    fill_chain.residues[r].atoms[a].movable=0;
                    fill_chain.residues[r].atoms[a].input=1;
                }
                i=0;
                for (a=pdb_entry.chains[c].residues[r].atoms.size();a<pdb_entry.chains[c].residues[r].atoms.size()+fill_residue_map[r].size();a++)
                {
                    fill_chain.residues[r].atoms.push_back(fill_atom);
                    fill_chain.residues[r].atoms[a].name=fill_residue_map[r][i];
                    fill_chain.residues[r].atoms[a].xyz.assign(3,0);
                    fill_chain.residues[r].atoms[a].movable=1;
                    fill_chain.residues[r].atoms[a].input=0;
                    i+=1;
                }
            }
            fullatom=true;
            //}

            //map<string, int> resnorder={{"  A  U", 0}, {"  U  A", 1}, {"  G  C", 0}, {"  C  G", 1}, {"  G  U", 0}, {"  U  G", 1}};
            double startene, endene;

            GetInitStr(c, pdb_entry, inputtedstrall, inputtedstrallflag);

            for (i=0;i<stembps.size();i++)
            {
                FragmentStructure superimposedstr;
                FragmentStructure inputtedstr;
                vector<vector<bool>> inputtedstrflag;
                int edgeleft_r = stembps[i][0][1];
                int bp_len = stembps[i].size();
                int bpnt_len = bp_len*2;
                vector<vector<int>> bindmap(bpnt_len, vector<int>(bpnt_len, 0));
                vector<bool> moveresidues(bpnt_len, false);

                int edgeright_r = stembps[i][bp_len-1][3];

                superimposedstr.residues.resize(bpnt_len);
                inputtedstr.residues.resize(bpnt_len);
                inputtedstrflag.resize(bpnt_len);

                for (j=0;j<stembps[i].size();j++)
                {
                    c1= stembps[i][j][0];
                    r1= stembps[i][j][1];
                    c2= stembps[i][j][2];
                    r2= stembps[i][j][3];
                    AssembleBasePairFrag(c1,r1,c2,r2, pdb_entry, ideal_rna["bp"], superimposedstr.residues[r1-edgeleft_r], superimposedstr.residues[r2-edgeright_r+bp_len], templateidx);
                    inputtedstr.residues[r1-edgeleft_r] = inputtedstrall.residues[r1];
                    inputtedstrflag[r1-edgeleft_r] = inputtedstrallflag[r1];
                    inputtedstr.residues[r2-edgeright_r+bp_len] = inputtedstrall.residues[r2];
                    inputtedstrflag[r2-edgeright_r+bp_len] = inputtedstrallflag[r2];

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
                }
                
                //if(option > 1){
                    //AttachHydrogen(superimposedstr, extparameters.hydparam, extparameters.hydparamidx);
                    
                    //fragmentoptimizationbp(superimposedstr, inputtedstr, inputtedstrflag, bp_len, startene, endene, extparameters, bindmap, stackmap, moveresidues, 5000);
                    //cout<<"Fragment optimization is finished (stem " << i+1 << "), " << "starting energy: " << startene << ", final energy: " << endene<<endl;
                    
                    //RemoveHydrogens(superimposedstr);
                //}

                for (j=0;j<stembps[i].size();j++)
                {
                    c1= stembps[i][j][0];
                    r1= stembps[i][j][1];
                    c2= stembps[i][j][2];
                    r2= stembps[i][j][3];

                    if(c == c1) updateflag[r1]=true;
                    if(c == c2) updateflag[r2]=true;

                    UpdateStemStructure(c1, r1, c2, r2, pdb_entry, fill_chain, fill_residue_map[r1], fill_residue_map[r2], superimposedstr.residues[r1-edgeleft_r], superimposedstr.residues[r2-edgeright_r+bp_len], fullatom);
                }
            }

            AddMissingAtoms(pdb_entry, c, fill_chain, fill_residue_map, updateflag, 2);

            updateflag.clear();
            updateflag.resize(pdb_entry.chains[c].residues.size(), false);

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
                    }

                    for(j=edgeleft_r;j<edgeright_r+1;j++)
                    {
                        AssembleNonBpNtFrag(c1, j, pdb_entry, ideal_rna["loop"], superimposedstr.residues[addres_f+j-edgeleft_r], templateidx);
                        inputtedstr.residues[addres_f+j-edgeleft_r] = inputtedstrall.residues[j];
                        inputtedstrflag[addres_f+j-edgeleft_r] = inputtedstrallflag[j];
                        
                        moveresidues[addres_f+j-edgeleft_r] = true;
                        if(j != edgeright_r){
                            bindmap[addres_f+j-edgeleft_r][addres_f+j-edgeleft_r+1] = 1;
                        }
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
                    }

                    //if(option > 1){
                        //AttachHydrogen(superimposedstr, extparameters.hydparam, extparameters.hydparamidx);

                        //fragmentoptimizationloop(superimposedstr, inputtedstr, inputtedstrflag, loop_len, startene, endene, extparameters, bindmap, stackmap, moveresidues, 5000);
                        //cout<<"Fragment optimization is finished (loop " << i+1 << "), " << "starting energy: " << startene << ", final energy: " << endene<<endl;

                        //RemoveHydrogens(superimposedstr);
                    //}

                    for(j=edgeleft_r;j<edgeright_r+1;j++){
                        if(c == c1) updateflag[j]=true;
                        UpdateLoopStructure(c1, j, pdb_entry, fill_chain, fill_residue_map[j], superimposedstr.residues[addres_f+j-edgeleft_r], fullatom);
                    }
                }
            }
            AddMissingAtoms(pdb_entry, c, fill_chain, fill_residue_map, updateflag, 2);
        }else{
            logger.log("Skip full-atomic structure genratation\n", true);
            GetInitStr(c, pdb_entry, inputtedstrall, inputtedstrallflag);

            double startene, endene;
        }
    }

    fill_atom.name.clear();
    fill_atom.xyz.clear();
    fill_residue.resn.clear();
    vector<AtomUnit>   ().swap(fill_residue.atoms);
    vector<ResidueUnit>().swap(fill_chain.residues);

    
    //clean up
    missing_atom_vec.clear();
    nonstd_atom_vec.clear();
    remove_residue_vec.clear();
    map<size_t,vector<string> >().swap(fill_residue_map);

    return true;
}

bool MissingRNAatom(ModelUnit &pdb_entry, 
    map<string, map<string, vector<ModelUnit>>> &ideal_rna,
    const vector<vector<vector<int>>> &stembps, const map<int, vector<pair<int, int>>> &loopnts, const map<int, vector<vector<vector<int>>>> &connectnts,
    map<int, int> &bp_template_idx, map<int, int> &nonbp_template_idx,
    FragmentStructure &inputtedstrall, vector<vector<bool>> &inputtedstrallflag,
    ExtParams &extparameters, const string &outdir, const bool &logflag, const int option)
    //map<string, vector<vector<double>>> &blengthparam, map<string, vector<vector<double>>> &adjblengthparam, 
    //map<string, vector<vector<double>>> &ntangleparam, map<string, vector<vector<int>>> &ntangleparamidx, map<string, vector<vector<double>>> &adjangleparam, map<string, vector<vector<int>>> &adjangleparamidx, 
    //map<string, vector<vector<double>>> &hydparam, map<string, vector<vector<int>>> &hydparamidx, const int option=0)
{

    Logger logger = Logger();
    logger.setPrint(false);
    if(logflag){
        logger.setVerbose(true);
        logger.setLogFile(outdir+"/logs/full_atom_generation.log");
    }

    size_t c;
    for (c=0;c<pdb_entry.chains.size();c++)
        MissingRNAatom(pdb_entry, c, ideal_rna, stembps, loopnts, connectnts, bp_template_idx, nonbp_template_idx, inputtedstrall, inputtedstrallflag, extparameters, logger, option);

    return true;
}


