
#include "BaseConformation.hpp"

bool SelectBPTemplate(ModelUnit &pdb_entry, map<string, vector<ModelUnit>> &ideal_rna, const vector<pair<double,vector<size_t>>> bp_vec, map<int, int> &bp_template_idx)
{
    bp_template_idx.clear();
    size_t bp,c1,c2,r1,r2;
    string resn_pair;
    size_t c,r;
    string base1, base2, bpkey;
    size_t a,a1,a2;
    double rmsd, minrmsd;
    int mintempidx;
    vector<double> tmp(3,0);
    vector<vector<double>> xyz_list1;
    vector<vector<double>> xyz_list2;
    vector<vector<double>> RotMatix, minRotMatix;
    vector<double> TranVect, minTranVect;
    
    for (bp=0;bp<bp_vec.size();bp++)
    {
        c1= bp_vec[bp].second[0];
        r1= bp_vec[bp].second[1];
        c2= bp_vec[bp].second[2];
        r2= bp_vec[bp].second[3];
        base1=pdb_entry.chains[c1].residues[r1].resn[2];
        base2=pdb_entry.chains[c2].residues[r2].resn[2];
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
            cout<<"Non-canonical bp!!"<<endl;
            return false;
        }

        size_t atomNum=pdb_entry.chains[c1].residues[r1].atoms.size()+pdb_entry.chains[c2].residues[r2].atoms.size()-22;
        minrmsd=10000000000; mintempidx=-1;
        for(int i=0;i<ideal_rna.at(bpkey).size();i++)
        {
            xyz_list1.clear();
            xyz_list2.clear();

            xyz_list1.assign(atomNum,tmp);
            xyz_list2.assign(atomNum,tmp);
            a=0;

            for (a1=11;a1<pdb_entry.chains[c1].residues[r1].atoms.size();a1++)
            {
                xyz_list1[a][0]=ideal_rna.at(bpkey)[i].chains[0].residues[idx1].atoms[a1].xyz[0];
                xyz_list1[a][1]=ideal_rna.at(bpkey)[i].chains[0].residues[idx1].atoms[a1].xyz[1];
                xyz_list1[a][2]=ideal_rna.at(bpkey)[i].chains[0].residues[idx1].atoms[a1].xyz[2];

                xyz_list2[a][0]=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz[0];
                xyz_list2[a][1]=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz[1];
                xyz_list2[a][2]=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz[2];

                a++;
            }

            for (a2=11;a2<pdb_entry.chains[c2].residues[r2].atoms.size();a2++)
            {
                xyz_list1[a][0]=ideal_rna.at(bpkey)[i].chains[0].residues[idx2].atoms[a2].xyz[0];
                xyz_list1[a][1]=ideal_rna.at(bpkey)[i].chains[0].residues[idx2].atoms[a2].xyz[1];
                xyz_list1[a][2]=ideal_rna.at(bpkey)[i].chains[0].residues[idx2].atoms[a2].xyz[2];

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
                    return false;
                }
            }

            RotMatix.clear();
            TranVect.clear();
            RotateCoor(xyz_list1,xyz_list2, RotMatix, TranVect);

            a=0;
            rmsd=0;
            for (a1=11;a1<pdb_entry.chains[c1].residues[r1].atoms.size();a1++)
            {
                ChangeCoor(ideal_rna.at(bpkey)[i].chains[c1].residues[idx1].atoms[a1].xyz,RotMatix, TranVect, xyz_list1[a]);
                tmp[0]=(pdb_entry.chains[c1].residues[r1].atoms[a1].xyz[0]-xyz_list1[a][0]);
                tmp[1]=(pdb_entry.chains[c1].residues[r1].atoms[a1].xyz[1]-xyz_list1[a][1]);
                tmp[2]=(pdb_entry.chains[c1].residues[r1].atoms[a1].xyz[2]-xyz_list1[a][2]);

                a++;
                rmsd+=tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2];
            }

            for (a2=11;a2<pdb_entry.chains[c2].residues[r2].atoms.size();a2++)
            {
                ChangeCoor(ideal_rna.at(bpkey)[i].chains[c2].residues[idx2].atoms[a2].xyz,RotMatix,TranVect,xyz_list1[a]);
                tmp[0]=(pdb_entry.chains[c2].residues[r2].atoms[a2].xyz[0]-xyz_list1[a][0]);
                tmp[1]=(pdb_entry.chains[c2].residues[r2].atoms[a2].xyz[1]-xyz_list1[a][1]);
                tmp[2]=(pdb_entry.chains[c2].residues[r2].atoms[a2].xyz[2]-xyz_list1[a][2]);

                a++;
                rmsd+=tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2];
            }
            //rmsd = sqrt(rmsd);

            if(rmsd < minrmsd){
                minrmsd=rmsd;
                mintempidx=i;
                minRotMatix=RotMatix;
                minTranVect=TranVect;
            }
        }

        bp_template_idx[bp]=mintempidx;
    }
    
    if(bp_vec.size() != 0){
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
        resn_pair.clear();
        return true;   
    }

    return true;
}

int fixBaseConformation(ResidueUnit &residue,
    ModelUnit&ideal_rna)
{
    vector<double> tmp(3,0);
    size_t atomNum=residue.atoms.size();
    vector<vector<double> > xyz_list1;
    vector<vector<double> > xyz_list2;
    vector<vector<double> > RotMatix;  // U
    vector<double> TranVect;  // t
    size_t a,b;

    int immovableCount=0;
    for (a=0;a<residue.atoms.size();a++)
        immovableCount+=(residue.atoms[a].movable==0);

    if (immovableCount==atomNum)
    {
        tmp.clear();
        return 0;
    }

    if (immovableCount>=3) atomNum=immovableCount;
    xyz_list1.assign(atomNum,tmp);
    xyz_list2.assign(atomNum,tmp);

    b=0;
    for (a=0;a<residue.atoms.size();a++)
    {
        if (immovableCount>=3 && residue.atoms[a].movable) continue;
        xyz_list1[b][0]=ideal_rna.chains[0].residues[0].atoms[a].xyz[0];
        xyz_list1[b][1]=ideal_rna.chains[0].residues[0].atoms[a].xyz[1];
        xyz_list1[b][2]=ideal_rna.chains[0].residues[0].atoms[a].xyz[2];

        xyz_list2[b][0]=residue.atoms[a].xyz[0];
        xyz_list2[b][1]=residue.atoms[a].xyz[1];
        xyz_list2[b][2]=residue.atoms[a].xyz[2];
        b++;
    }

    atomNum=countUniqAtom(xyz_list2);
    if (atomNum<xyz_list2.size())
        cerr<<"duplicated coordinate in "
            <<residue.resn<<residue.resi<<residue.icode<<endl;

    //if (residue.atoms.size()>=3)
    if (atomNum>=3)
        RotateCoor(xyz_list1,xyz_list2, RotMatix, TranVect);
    else
    {
        cerr<<"too few ("<<atomNum<<") unique atoms to fill "
            <<residue.resn<<residue.resi<<residue.icode<<endl;
        RotMatix.push_back(tmp);
        RotMatix.push_back(tmp);
        RotMatix.push_back(tmp);
        RotMatix[0][0]=RotMatix[1][1]=RotMatix[2][2]=1;
        
        TranVect=tmp;
        for (a=0;a<xyz_list1.size();a++)
        {
            TranVect[0]+=xyz_list2[a][0]-xyz_list1[a][0];
            TranVect[1]+=xyz_list2[a][1]-xyz_list1[a][1];
            TranVect[2]+=xyz_list2[a][2]-xyz_list1[a][2];
        }
        TranVect[0]/=xyz_list1.size();
        TranVect[1]/=xyz_list1.size();
        TranVect[2]/=xyz_list1.size();
    }

    double rmsd=0;
    for (a=0;a<xyz_list1.size();a++)
    {
        ChangeCoor(xyz_list1[a], RotMatix, TranVect, tmp);
        rmsd+=Points2Distance(tmp,xyz_list2[a]);
    }
    rmsd=sqrt(rmsd);
    
    for (a=0;a<xyz_list1.size();a++)
    {
        xyz_list1[a].clear();
        xyz_list2[a].clear();
    }
    xyz_list1.clear();
    xyz_list2.clear();
    
    int moved=0;
    if (rmsd>=0.1)
    {
        //cout<<residue.resn<<residue.resi<<residue.icode<<" rmsd="<<rmsd<<endl;
        for (a=0;a<residue.atoms.size();a++)
        {
            if (!residue.atoms[a].movable) continue;
            moved++;
            ChangeCoor(
                ideal_rna.chains[0].residues[0].atoms[a].xyz,
                RotMatix, TranVect, residue.atoms[a].xyz);
            if (immovableCount>=3) residue.atoms[a].movable=false;
        }
    }

    TranVect.clear();
    for (a=0;a<3;a++) RotMatix[a].clear();
    RotMatix.clear();
    return moved;
}

int fixBasePairConformation(ResidueUnit &residue1, ResidueUnit &residue2,
    ModelUnit&ideal_rna)
{
    vector<double> tmp(3,0);
    size_t atomNum=residue1.atoms.size()-11;
    atomNum      +=residue2.atoms.size()-11;
    vector<vector<double> > xyz_list1;
    vector<vector<double> > xyz_list2;
    string base, base1, base2, bpkey;
    vector<vector<double> > RotMatix;  // U
    vector<double> TranVect;  // t
    size_t a,b;

    //int immovableCount=0;
    //for (a=11;a<residue1.atoms.size();a++) immovableCount+=1;
        //immovableCount+=(residue1.atoms[a].movable==0);
    //for (a=11;a<residue2.atoms.size();a++) immovableCount+=1;
        //immovableCount+=(residue2.atoms[a].movable==0);

    //if (immovableCount==atomNum)
    //{
        //tmp.clear();
        //return 0;
    //}

    xyz_list1.assign(atomNum,vector<double>(3, 0));
    xyz_list2.assign(atomNum,vector<double>(3, 0));

    b=0;
    for (a=11;a<residue1.atoms.size();a++)
    {
        xyz_list1[b][0]=ideal_rna.chains[0].residues[0].atoms[a].xyz[0];
        xyz_list1[b][1]=ideal_rna.chains[0].residues[0].atoms[a].xyz[1];
        xyz_list1[b][2]=ideal_rna.chains[0].residues[0].atoms[a].xyz[2];

        xyz_list2[b][0]=residue1.atoms[a].xyz[0];
        xyz_list2[b][1]=residue1.atoms[a].xyz[1];
        xyz_list2[b][2]=residue1.atoms[a].xyz[2];
        b++;
    }
    for (a=11;a<residue2.atoms.size();a++)
    {
        xyz_list1[b][0]=ideal_rna.chains[0].residues[1].atoms[a].xyz[0];
        xyz_list1[b][1]=ideal_rna.chains[0].residues[1].atoms[a].xyz[1];
        xyz_list1[b][2]=ideal_rna.chains[0].residues[1].atoms[a].xyz[2];

        xyz_list2[b][0]=residue2.atoms[a].xyz[0];
        xyz_list2[b][1]=residue2.atoms[a].xyz[1];
        xyz_list2[b][2]=residue2.atoms[a].xyz[2];
        b++;
    }

    atomNum=countUniqAtom(xyz_list2);
    if (atomNum<3)
    {
        for (a=0;a<xyz_list1.size();a++)
        {
            xyz_list1[a].clear();
            xyz_list2[a].clear();
        }
        xyz_list1.clear();
        xyz_list2.clear();
        tmp.clear();
        return 0;
    }

    RotateCoor(xyz_list1,xyz_list2, RotMatix, TranVect);
    
    double N=0;
    double rmsd=0;
    double diff_x, diff_y, diff_z;
    for (a=11;a<xyz_list1.size();a++)
    {
        ChangeCoor(xyz_list1[a], RotMatix, TranVect, tmp);

        diff_x = tmp[0] - xyz_list2[a][0];
        diff_y = tmp[1] - xyz_list2[a][1];
        diff_z = tmp[2] - xyz_list2[a][2];

        rmsd+=diff_x*diff_x+diff_y*diff_y+diff_z*diff_z;
        N += 1;
    }
    rmsd=sqrt(rmsd/N);
    
    for (a=0;a<xyz_list1.size();a++)
    {
        xyz_list1[a].clear();
        xyz_list2[a].clear();
    }
    xyz_list1.clear();
    xyz_list2.clear();
    
    int moved=0;
    if (rmsd>=0.025)
    {
        //cout<<residue1.resn<<residue1.resi<<residue1.icode<<" :: "
        //    <<residue2.resn<<residue2.resi<<residue2.icode<<" rmsd="<<rmsd<<endl;
        b=0;
        for (a=11;a<residue1.atoms.size();a++)
        {
            //if (!residue1.atoms[a].movable) continue;
            moved++;
            ChangeCoor(
                ideal_rna.chains[0].residues[0].atoms[a].xyz,
                RotMatix, TranVect, residue1.atoms[a].xyz);
        }
        for (a=11;a<residue2.atoms.size();a++)
        {
            //if (!residue2.atoms[a].movable) continue;
            moved++;
            ChangeCoor(
                ideal_rna.chains[0].residues[1].atoms[a].xyz,
                RotMatix, TranVect, residue2.atoms[a].xyz);
        }
    }

    TranVect.clear();
    for (a=0;a<3;a++) RotMatix[a].clear();
    RotMatix.clear();
    return moved;
}

int fixBaseConformation(ModelUnit &pdb_entry, 
    map<string, vector<ModelUnit>> &ideal_rna,
    const vector<pair<double,vector<size_t> > >&bp_vec, map<int, int> &bp_template_idx)
{
    map<string, int> resnorder={{"  A  U", 0}, {"  U  A", 1}, {"  G  C", 0}, {"  C  G", 1}, {"  G  U", 0}, {"  U  G", 1}};
    size_t bp,c1,c2,r1,r2;
    string resn_pair;
    size_t c,r;
    string resn1, resn2, resn;
    int moved=0;
    
    for (bp=0;bp<bp_vec.size();bp++)
    {
        c1= bp_vec[bp].second[0];
        r1= bp_vec[bp].second[1];
        c2= bp_vec[bp].second[2];
        r2= bp_vec[bp].second[3];
        resn1=pdb_entry.chains[c1].residues[r1].resn[2];
        resn2=pdb_entry.chains[c2].residues[r2].resn[2];
        resn_pair=pdb_entry.chains[c1].residues[r1].resn+
                  pdb_entry.chains[c2].residues[r2].resn;
        if (resn_pair!="  A  U" && resn_pair!="  U  A" && 
            resn_pair!="  C  G" && resn_pair!="  G  C" &&
            resn_pair!="  G  U" && resn_pair!="  U  G" ) continue;

        if(bp_template_idx.count(bp)){
            if(resnorder[resn_pair] == 0){
                moved+=fixBasePairConformation(pdb_entry.chains[c1].residues[r1], pdb_entry.chains[c2].residues[r2], ideal_rna[resn1+"-"+resn2][bp_template_idx[bp]]);
            }else if(resnorder[resn_pair] == 1){
                moved+=fixBasePairConformation(pdb_entry.chains[c2].residues[r2], pdb_entry.chains[c1].residues[r1], ideal_rna[resn2+"-"+resn1][bp_template_idx[bp]]);
            }
        }
    }
    
    /*
    for (c=0;c<pdb_entry.chains.size();c++){
        for (r=0;r<pdb_entry.chains[c].residues.size();r++){
            if(nonbp_template_idx.count(r)){
                resn=pdb_entry.chains[c].residues[r].resn[2];
                moved+=fixBaseConformation(pdb_entry.chains[c].residues[r], ideal_rna["nt"][resn][nonbp_template_idx[r]]);
            }
        }
    }
    */
    
    resnorder.clear();
    resn_pair.clear();
    return moved;
}

