
#include "FragmentRefine.hpp"

bool fillFragmentInBasePair(const int c1, const int r1,
    const int c2, const int r2, const ModelUnit &pdb_entry, ChainUnit &fill_chain,
    map<string, vector<ModelUnit>> &ideal_rna, int &templateidx, bool fullatom)
{
    int i;
    vector<double> tmp(3,0);
    vector<vector<double> > xyz_list1;
    vector<vector<double> > xyz_list2;
    vector<vector<double> > RotMatix, minRotMatix;  // U
    vector<double> TranVect, minTranVect;  // t
    string key, key1, key2;
    bool stack_prev, stack_next;
    string base1, base2, bpkey;
    double rmsd, minrmsd;
    int mintempidx;

    map<string, vector<string>> full_atoms={
        {"A", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'", " N9 ", " C8 ", " N7 ", " C5 ", " C6 ", " N6 ", " N1 ", " C2 ", " N3 ", " C4 "}},
        {"C", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'", " N1 ", " C2 ", " O2 ", " N3 ", " C4 ", " N4 ", " C5 ", " C6 "}},
        {"G", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'", " N9 ", " C8 ", " N7 ", " C5 ", " C6 ", " O6 ", " N1 ", " C2 ", " N2 ", " N3 ", " C4 "}},
        {"U", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'", " N1 ", " C2 ", " O2 ", " N3 ", " C4 ", " O4 ", " C5 ", " C6 "}},
    };
    
    stack_prev=false;
    if(r1-1 >= 0 && r2+1 < pdb_entry.chains[c2].residues.size()) stack_prev=true;
    stack_next=false;
    if(r1+1 < pdb_entry.chains[c1].residues.size() && r2-1 >= 0) stack_next=true;

    base1=pdb_entry.chains[c1].residues[r1].resn[2];
    base2=pdb_entry.chains[c2].residues[r2].resn[2];
    bpkey=base1+"-"+base2;

    size_t a,a1,a2;
    size_t atomNum=0;

    //pdb_entry.chains[c1].residues[r1].atoms[a].movable==0

    atomNum+=pdb_entry.chains[c1].residues[r1].atoms.size();
    atomNum+=pdb_entry.chains[c2].residues[r2].atoms.size();
        
    if (stack_prev)
    {
        atomNum+=pdb_entry.chains[c1].residues[r1-1].atoms.size();
        atomNum+=pdb_entry.chains[c2].residues[r2+1].atoms.size();
    }
    if (stack_next)
    {
        atomNum+=pdb_entry.chains[c1].residues[r1+1].atoms.size();
        atomNum+=pdb_entry.chains[c2].residues[r2-1].atoms.size();
    }

    if (ideal_rna.count(bpkey) == 0){
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
            xyz_list1[a][0]=ideal_rna[bpkey][i].chains[0].residues[1].atoms[a1].xyz[0];
            xyz_list1[a][1]=ideal_rna[bpkey][i].chains[0].residues[1].atoms[a1].xyz[1];
            xyz_list1[a][2]=ideal_rna[bpkey][i].chains[0].residues[1].atoms[a1].xyz[2];
            xyz_list2[a][0]=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz[0];
            xyz_list2[a][1]=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz[1];
            xyz_list2[a][2]=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz[2];
            a++;
        }
    
        for (a2=0;a2<pdb_entry.chains[c2].residues[r2].atoms.size();a2++)
        {
            xyz_list1[a][0]=ideal_rna[bpkey][i].chains[0].residues[4].atoms[a2].xyz[0];
            xyz_list1[a][1]=ideal_rna[bpkey][i].chains[0].residues[4].atoms[a2].xyz[1];
            xyz_list1[a][2]=ideal_rna[bpkey][i].chains[0].residues[4].atoms[a2].xyz[2];
            xyz_list2[a][0]=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz[0];
            xyz_list2[a][1]=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz[1];
            xyz_list2[a][2]=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz[2];
            a++;
        }

        if (stack_prev)
        {
            for (a1=0;a1<ideal_rna[bpkey][i].chains[0].residues[0].atoms.size();a1++)
            {
                xyz_list1[a][0]=ideal_rna[bpkey][i].chains[0].residues[0].atoms[a1].xyz[0];
                xyz_list1[a][1]=ideal_rna[bpkey][i].chains[0].residues[0].atoms[a1].xyz[1];
                xyz_list1[a][2]=ideal_rna[bpkey][i].chains[0].residues[0].atoms[a1].xyz[2];
                xyz_list2[a][0]=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz[0];
                xyz_list2[a][1]=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz[1];
                xyz_list2[a][2]=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz[2];
                a++;
            }

            for (a2=0;a2<ideal_rna[bpkey][i].chains[0].residues[5].atoms.size();a2++)
            {
                xyz_list1[a][0]=ideal_rna[bpkey][i].chains[0].residues[5].atoms[a2].xyz[0];
                xyz_list1[a][1]=ideal_rna[bpkey][i].chains[0].residues[5].atoms[a2].xyz[1];
                xyz_list1[a][2]=ideal_rna[bpkey][i].chains[0].residues[5].atoms[a2].xyz[2];
                xyz_list2[a][0]=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz[0];
                xyz_list2[a][1]=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz[1];
                xyz_list2[a][2]=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz[2];
                a++;
            }
        }

        if (stack_next)
        {
            for (a1=0;a1<ideal_rna[bpkey][i].chains[0].residues[2].atoms.size();a1++)
            {
                xyz_list1[a][0]=ideal_rna[bpkey][i].chains[0].residues[2].atoms[a1].xyz[0];
                xyz_list1[a][1]=ideal_rna[bpkey][i].chains[0].residues[2].atoms[a1].xyz[1];
                xyz_list1[a][2]=ideal_rna[bpkey][i].chains[0].residues[2].atoms[a1].xyz[2];
                xyz_list2[a][0]=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz[0];
                xyz_list2[a][1]=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz[1];
                xyz_list2[a][2]=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz[2];
                a++;
            }

            for (a2=0;a2<ideal_rna[bpkey][i].chains[0].residues[3].atoms.size();a2++)
            {
                xyz_list1[a][0]=ideal_rna[bpkey][i].chains[0].residues[3].atoms[a2].xyz[0];
                xyz_list1[a][1]=ideal_rna[bpkey][i].chains[0].residues[3].atoms[a2].xyz[1];
                xyz_list1[a][2]=ideal_rna[bpkey][i].chains[0].residues[3].atoms[a2].xyz[2];
                xyz_list2[a][0]=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz[0];
                xyz_list2[a][1]=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz[1];
                xyz_list2[a][2]=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz[2];
                a++;
            }
        }

        if(i == 0){
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
                //key.clear();
                //key1.clear();
                //key2.clear();
                return false;
            }
        }
        
        RotMatix.clear();
        TranVect.clear();
        RotateCoor(xyz_list1,xyz_list2, RotMatix, TranVect);

        tmp[0]=0;tmp[1]=0;tmp[2]=0;
        atomNum=0;
        for (a=0;a<pdb_entry.chains[c1].residues[r1].atoms.size();a++)
        {   
            ChangeCoor(ideal_rna[bpkey][i].chains[0].residues[1].atoms[a].xyz,RotMatix,TranVect,xyz_list1[atomNum]);
            tmp[0]+=(pdb_entry.chains[c1].residues[r1].atoms[a].xyz[0]-xyz_list1[atomNum][0]);
            tmp[1]+=(pdb_entry.chains[c1].residues[r1].atoms[a].xyz[1]-xyz_list1[atomNum][1]);
            tmp[2]+=(pdb_entry.chains[c1].residues[r1].atoms[a].xyz[2]-xyz_list1[atomNum][2]);
            atomNum++;
        }
        for (a=0;a<pdb_entry.chains[c2].residues[r2].atoms.size();a++)
        {
            ChangeCoor(ideal_rna[bpkey][i].chains[0].residues[4].atoms[a].xyz,RotMatix,TranVect,xyz_list1[atomNum]);
            tmp[0]+=(pdb_entry.chains[c2].residues[r2].atoms[a].xyz[0]-xyz_list1[atomNum][0]);
            tmp[1]+=(pdb_entry.chains[c2].residues[r2].atoms[a].xyz[1]-xyz_list1[atomNum][1]);
            tmp[2]+=(pdb_entry.chains[c2].residues[r2].atoms[a].xyz[2]-xyz_list1[atomNum][2]);
            atomNum++;
        }

        TranVect[0]+=(tmp[0]/atomNum);
        TranVect[1]+=(tmp[1]/atomNum);
        TranVect[2]+=(tmp[2]/atomNum);

        atomNum=0;
        rmsd=0;
        for (a=0;a<pdb_entry.chains[c1].residues[r1].atoms.size();a++)
        {
            ChangeCoor(ideal_rna[bpkey][i].chains[0].residues[1].atoms[a].xyz,RotMatix, TranVect, xyz_list1[atomNum]);
            tmp[0]=(pdb_entry.chains[c1].residues[r1].atoms[a].xyz[0]-xyz_list1[atomNum][0]);
            tmp[1]=(pdb_entry.chains[c1].residues[r1].atoms[a].xyz[1]-xyz_list1[atomNum][1]);
            tmp[2]=(pdb_entry.chains[c1].residues[r1].atoms[a].xyz[2]-xyz_list1[atomNum][2]);
            atomNum++;
            rmsd+=sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
        }
        for (a=0;a<pdb_entry.chains[c2].residues[r2].atoms.size();a++)
        {
            ChangeCoor(ideal_rna[bpkey][i].chains[0].residues[4].atoms[a].xyz,RotMatix,TranVect,xyz_list1[atomNum]);
            tmp[0]=(pdb_entry.chains[c2].residues[r2].atoms[a].xyz[0]-xyz_list1[atomNum][0]);
            tmp[1]=(pdb_entry.chains[c2].residues[r2].atoms[a].xyz[1]-xyz_list1[atomNum][1]);
            tmp[2]=(pdb_entry.chains[c2].residues[r2].atoms[a].xyz[2]-xyz_list1[atomNum][2]);
            atomNum++;
            rmsd+=sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
        }

        if (stack_prev)
        {
            for (a=0;a<ideal_rna[bpkey][i].chains[0].residues[0].atoms.size();a++)
            {
                ChangeCoor(ideal_rna[bpkey][i].chains[0].residues[0].atoms[a].xyz,RotMatix, TranVect, xyz_list1[atomNum]);
                tmp[0]=(pdb_entry.chains[c1].residues[r1-1].atoms[a].xyz[0]-xyz_list1[atomNum][0]);
                tmp[1]=(pdb_entry.chains[c1].residues[r1-1].atoms[a].xyz[1]-xyz_list1[atomNum][1]);
                tmp[2]=(pdb_entry.chains[c1].residues[r1-1].atoms[a].xyz[2]-xyz_list1[atomNum][2]);
                atomNum++;
                rmsd+=sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
            }
            for (a=0;a<ideal_rna[bpkey][i].chains[0].residues[5].atoms.size();a++)
            {
                ChangeCoor(ideal_rna[bpkey][i].chains[0].residues[5].atoms[a].xyz,RotMatix,TranVect,xyz_list1[atomNum]);
                tmp[0]=(pdb_entry.chains[c2].residues[r2+1].atoms[a].xyz[0]-xyz_list1[atomNum][0]);
                tmp[1]=(pdb_entry.chains[c2].residues[r2+1].atoms[a].xyz[1]-xyz_list1[atomNum][1]);
                tmp[2]=(pdb_entry.chains[c2].residues[r2+1].atoms[a].xyz[2]-xyz_list1[atomNum][2]);
                atomNum++;
                rmsd+=sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
            }
        }

        if (stack_next)
        {
            for (a=0;a<ideal_rna[bpkey][i].chains[0].residues[2].atoms.size();a++)
            {
                ChangeCoor(ideal_rna[bpkey][i].chains[0].residues[2].atoms[a].xyz,RotMatix, TranVect, xyz_list1[atomNum]);
                tmp[0]=(pdb_entry.chains[c1].residues[r1+1].atoms[a].xyz[0]-xyz_list1[atomNum][0]);
                tmp[1]=(pdb_entry.chains[c1].residues[r1+1].atoms[a].xyz[1]-xyz_list1[atomNum][1]);
                tmp[2]=(pdb_entry.chains[c1].residues[r1+1].atoms[a].xyz[2]-xyz_list1[atomNum][2]);
                atomNum++;
                rmsd+=sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
            }
            for (a=0;a<ideal_rna[bpkey][i].chains[0].residues[3].atoms.size();a++)
            {
                ChangeCoor(ideal_rna[bpkey][i].chains[0].residues[3].atoms[a].xyz,RotMatix,TranVect,xyz_list1[atomNum]);
                tmp[0]=(pdb_entry.chains[c2].residues[r2-1].atoms[a].xyz[0]-xyz_list1[atomNum][0]);
                tmp[1]=(pdb_entry.chains[c2].residues[r2-1].atoms[a].xyz[1]-xyz_list1[atomNum][1]);
                tmp[2]=(pdb_entry.chains[c2].residues[r2-1].atoms[a].xyz[2]-xyz_list1[atomNum][2]);
                atomNum++;
                rmsd+=sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
            }
        }

        if(rmsd < minrmsd){
            minrmsd=rmsd;
            mintempidx=i;
            minRotMatix=RotMatix;
            minTranVect=TranVect;
        }
    }

    templateidx=mintempidx;

    if(fullatom){
        for (a=0;a<pdb_entry.chains[c1].residues[r1].atoms.size();a++)
        {
            ChangeCoor(ideal_rna[bpkey][mintempidx].chains[0].residues[1].atoms[a].xyz,
                minRotMatix, minTranVect, fill_chain.residues[r1].atoms[a].xyz);
        }
        for (a=0;a<pdb_entry.chains[c2].residues[r2].atoms.size();a++)
        {
            ChangeCoor(ideal_rna[bpkey][mintempidx].chains[0].residues[4].atoms[a].xyz,
                minRotMatix, minTranVect, fill_chain.residues[r2].atoms[a].xyz);
        }
    }else{
        for (a=0;a<pdb_entry.chains[c1].residues[r1].atoms.size();a++)
        {
            if(pdb_entry.chains[c1].residues[r1].atoms[a].movable==1){
                ChangeCoor(ideal_rna[bpkey][mintempidx].chains[0].residues[1].atoms[a].xyz,
                    minRotMatix, minTranVect, fill_chain.residues[r1].atoms[a].xyz);
            }else{
                fill_chain.residues[r1].atoms[a].xyz=pdb_entry.chains[c1].residues[r1].atoms[a].xyz;
            }
        }

        for (a=0;a<pdb_entry.chains[c2].residues[r2].atoms.size();a++)
        {
            if(pdb_entry.chains[c2].residues[r2].atoms[a].movable==1){
                ChangeCoor(ideal_rna[bpkey][mintempidx].chains[0].residues[4].atoms[a].xyz,
                    minRotMatix, minTranVect, fill_chain.residues[r2].atoms[a].xyz);
            }else{
                fill_chain.residues[r2].atoms[a].xyz=pdb_entry.chains[c2].residues[r2].atoms[a].xyz;
            }
        }
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
    minTranVect.clear();
    for (a=0;a<3;a++){
        RotMatix[a].clear();
        minRotMatix[a].clear();
    }
    RotMatix.clear();
    minRotMatix.clear();
    return true;
}

bool fillFragmentLocal(const int c, const int r,
    const ModelUnit &pdb_entry, ChainUnit &fill_chain, map<string, vector<ModelUnit>> &ideal_rna, int &templateidx, bool fullatom)
{
    int i;
    vector<double> tmp(3,0);
    vector<vector<double> > xyz_list1;
    vector<vector<double> > xyz_list2;
    vector<vector<double> > RotMatix, minRotMatix;  // U
    vector<double> TranVect, minTranVect;  // t
    string key, key1, key2;
    bool stack_prev, stack_next;
    double rmsd, minrmsd;
    int mintempidx;
    string base2;

    map<string, vector<string>> full_atoms={
        {"A", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'", " N9 ", " C8 ", " N7 ", " C5 ", " C6 ", " N6 ", " N1 ", " C2 ", " N3 ", " C4 "}},
        {"C", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'", " N1 ", " C2 ", " O2 ", " N3 ", " C4 ", " N4 ", " C5 ", " C6 "}},
        {"G", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'", " N9 ", " C8 ", " N7 ", " C5 ", " C6 ", " O6 ", " N1 ", " C2 ", " N2 ", " N3 ", " C4 "}},
        {"U", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'", " N1 ", " C2 ", " O2 ", " N3 ", " C4 ", " O4 ", " C5 ", " C6 "}},
    };
    
    stack_prev=false;
    stack_next=false;
    if(r-1 >= 0) stack_prev=true;
    if(r+1 < pdb_entry.chains[c].residues.size()) stack_next=true;

    base2=pdb_entry.chains[c].residues[r].resn[2];

    size_t a,a1,a2;
    size_t atomNum=0;

    atomNum+=pdb_entry.chains[c].residues[r].atoms.size();
    if (stack_prev) atomNum+=pdb_entry.chains[c].residues[r-1].atoms.size(); 
    if (stack_next) atomNum+=pdb_entry.chains[c].residues[r+1].atoms.size();

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
            xyz_list1[a][0]=ideal_rna[base2][i].chains[0].residues[1].atoms[a1].xyz[0];
            xyz_list1[a][1]=ideal_rna[base2][i].chains[0].residues[1].atoms[a1].xyz[1];
            xyz_list1[a][2]=ideal_rna[base2][i].chains[0].residues[1].atoms[a1].xyz[2];

            xyz_list2[a][0]=pdb_entry.chains[c].residues[r].atoms[a1].xyz[0];
            xyz_list2[a][1]=pdb_entry.chains[c].residues[r].atoms[a1].xyz[1];
            xyz_list2[a][2]=pdb_entry.chains[c].residues[r].atoms[a1].xyz[2];
            a++;
        }

        if (stack_prev)
        {
            for (a1=0;a1<ideal_rna[base2][i].chains[0].residues[0].atoms.size();a1++)
            {
                xyz_list1[a][0]=ideal_rna[base2][i].chains[0].residues[0].atoms[a1].xyz[0];
                xyz_list1[a][1]=ideal_rna[base2][i].chains[0].residues[0].atoms[a1].xyz[1];
                xyz_list1[a][2]=ideal_rna[base2][i].chains[0].residues[0].atoms[a1].xyz[2];

                xyz_list2[a][0]=pdb_entry.chains[c].residues[r-1].atoms[a1].xyz[0];
                xyz_list2[a][1]=pdb_entry.chains[c].residues[r-1].atoms[a1].xyz[1];
                xyz_list2[a][2]=pdb_entry.chains[c].residues[r-1].atoms[a1].xyz[2];
                a++;
            }
        }

        if (stack_next)
        {
            for (a1=0;a1<ideal_rna[base2][i].chains[0].residues[2].atoms.size();a1++)
            {
                xyz_list1[a][0]=ideal_rna[base2][i].chains[0].residues[2].atoms[a1].xyz[0];
                xyz_list1[a][1]=ideal_rna[base2][i].chains[0].residues[2].atoms[a1].xyz[1];
                xyz_list1[a][2]=ideal_rna[base2][i].chains[0].residues[2].atoms[a1].xyz[2];
                xyz_list2[a][0]=pdb_entry.chains[c].residues[r+1].atoms[a1].xyz[0];
                xyz_list2[a][1]=pdb_entry.chains[c].residues[r+1].atoms[a1].xyz[1];
                xyz_list2[a][2]=pdb_entry.chains[c].residues[r+1].atoms[a1].xyz[2];
                a++;
            }
        }

        if(i == 0){
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
                key.clear();
                key1.clear();
                key2.clear();
                return false;
            }
        }

        RotMatix.clear();
        TranVect.clear();
        RotateCoor(xyz_list1,xyz_list2, RotMatix, TranVect);

        tmp[0]=0; tmp[1]=0; tmp[2]=0;
        atomNum=0;
        for (a=0;a<pdb_entry.chains[c].residues[r].atoms.size();a++)
        {
            ChangeCoor(ideal_rna[base2][i].chains[0].residues[1].atoms[a].xyz,RotMatix,TranVect,xyz_list1[atomNum]);
            tmp[0]+=(pdb_entry.chains[c].residues[r].atoms[a].xyz[0]-xyz_list1[atomNum][0]);
            tmp[1]+=(pdb_entry.chains[c].residues[r].atoms[a].xyz[1]-xyz_list1[atomNum][1]);
            tmp[2]+=(pdb_entry.chains[c].residues[r].atoms[a].xyz[2]-xyz_list1[atomNum][2]);
            atomNum++;
        }
        TranVect[0]+=(tmp[0]/atomNum);
        TranVect[1]+=(tmp[1]/atomNum);
        TranVect[2]+=(tmp[2]/atomNum);

        atomNum=0;
        rmsd=0;
        for (a=0;a<pdb_entry.chains[c].residues[r].atoms.size();a++)
        {
            ChangeCoor(ideal_rna[base2][i].chains[0].residues[1].atoms[a].xyz,RotMatix, TranVect, xyz_list1[atomNum]);
            tmp[0]=(pdb_entry.chains[c].residues[r].atoms[a].xyz[0]-xyz_list1[atomNum][0]);
            tmp[1]=(pdb_entry.chains[c].residues[r].atoms[a].xyz[1]-xyz_list1[atomNum][1]);
            tmp[2]=(pdb_entry.chains[c].residues[r].atoms[a].xyz[2]-xyz_list1[atomNum][2]);
            atomNum++;
            rmsd+=sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
        }

        if (stack_prev)
        {
            for (a=0;a<ideal_rna[base2][i].chains[0].residues[0].atoms.size();a++)
            {
                ChangeCoor(ideal_rna[base2][i].chains[0].residues[0].atoms[a].xyz,RotMatix, TranVect, xyz_list1[atomNum]);
                tmp[0]=(pdb_entry.chains[c].residues[r-1].atoms[a].xyz[0]-xyz_list1[atomNum][0]);
                tmp[1]=(pdb_entry.chains[c].residues[r-1].atoms[a].xyz[1]-xyz_list1[atomNum][1]);
                tmp[2]=(pdb_entry.chains[c].residues[r-1].atoms[a].xyz[2]-xyz_list1[atomNum][2]);
                atomNum++;
                rmsd+=sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
            }
        }

        if (stack_next)
        {
            for (a=0;a<ideal_rna[base2][i].chains[0].residues[2].atoms.size();a++)
            {
                ChangeCoor(ideal_rna[base2][i].chains[0].residues[2].atoms[a].xyz,RotMatix, TranVect, xyz_list1[atomNum]);
                tmp[0]=(pdb_entry.chains[c].residues[r+1].atoms[a].xyz[0]-xyz_list1[atomNum][0]);
                tmp[1]=(pdb_entry.chains[c].residues[r+1].atoms[a].xyz[1]-xyz_list1[atomNum][1]);
                tmp[2]=(pdb_entry.chains[c].residues[r+1].atoms[a].xyz[2]-xyz_list1[atomNum][2]);
                atomNum++;
                rmsd+=sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
            }
        }

        if(rmsd < minrmsd){
            minrmsd=rmsd;
            mintempidx=i;
            minRotMatix=RotMatix;
            minTranVect=TranVect;
        }
    }

    templateidx=mintempidx;

    if(fullatom){
        for (a=0;a<pdb_entry.chains[c].residues[r].atoms.size();a++)
        {   
            ChangeCoor(ideal_rna[base2][mintempidx].chains[0].residues[1].atoms[a].xyz,
                    minRotMatix, minTranVect, fill_chain.residues[r].atoms[a].xyz);
        }
    }else{
        for (a=0;a<pdb_entry.chains[c].residues[r].atoms.size();a++)
        { 
            if(pdb_entry.chains[c].residues[r].atoms[a].movable==1){
                ChangeCoor(ideal_rna[base2][mintempidx].chains[0].residues[1].atoms[a].xyz,
                    minRotMatix, minTranVect, fill_chain.residues[r].atoms[a].xyz);
            }else{
                fill_chain.residues[r].atoms[a].xyz=pdb_entry.chains[c].residues[r].atoms[a].xyz;
            }
        }
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
        minRotMatix[a].clear();
    }
    RotMatix.clear();
    minRotMatix.clear();
    key.clear();
    key1.clear();
    key2.clear();
    return true;
}

bool FragRef(ModelUnit &pdb_entry, const size_t c,
    map<string, map<string, vector<ModelUnit>>> &ideal_rna,
    const vector<pair<double,vector<size_t>>>&bp_vec, map<int, int> &bp_template_idx, map<int, int> &nonbp_template_idx,
    const int option)
{
    size_t r,i,j;
    vector<int> nonbpresidues;
    bool fullatom;
    int templateidx;

    size_t bp;
    int r1,r2,c1,c2;

    ChainUnit fill_chain;
    ResidueUnit fill_residue;
    AtomUnit fill_atom;
    size_t a;

    for (r=0;r<pdb_entry.chains[c].residues.size();r++)
    {
        nonbpresidues.push_back(r);
        fill_chain.residues.push_back(fill_residue);

        for (a=0;a<pdb_entry.chains[c].residues[r].atoms.size();a++)
        {
            fill_chain.residues[r].atoms.push_back(fill_atom);
            fill_chain.residues[r].atoms[a].name=pdb_entry.chains[c].residues[r].atoms[a].name;
            fill_chain.residues[r].atoms[a].xyz.assign(3,0);
            fill_chain.residues[r].atoms[a].movable=0;
        }
    }

    if(option == 1) fullatom=false;
    else if(option == 2) fullatom=true;

    map<string, int> resnorder={{"  A  U", 0}, {"  U  A", 1}, {"  G  C", 0}, {"  C  G", 1}, {"  G  U", 0}, {"  U  G", 1}};
    for (bp=0;bp<bp_vec.size();bp++)
    {
        //if (bp_vec[bp].second[4]) continue; // Watson-Crick pair only
        c1= bp_vec[bp].second[0];
        r1= bp_vec[bp].second[1];
        c2= bp_vec[bp].second[2];
        r2= bp_vec[bp].second[3];

        if(c1==c || c2==c){
            if(resnorder[pdb_entry.chains[c1].residues[r1].resn+pdb_entry.chains[c2].residues[r2].resn] == 0)
                fillFragmentInBasePair(c1,r1,c2,r2, pdb_entry, fill_chain, ideal_rna["bp"], templateidx, fullatom);
                
            if(resnorder[pdb_entry.chains[c1].residues[r1].resn+pdb_entry.chains[c2].residues[r2].resn] == 1)
                fillFragmentInBasePair(c2,r2,c1,r1, pdb_entry, fill_chain, ideal_rna["bp"], templateidx, fullatom);
            
            bp_template_idx[bp]=templateidx;
        }

        if (find(nonbpresidues.begin(), nonbpresidues.end(), r1) != nonbpresidues.end()) {
            nonbpresidues.erase(remove(nonbpresidues.begin(), nonbpresidues.end(), r1), nonbpresidues.end());
        }
        if (find(nonbpresidues.begin(), nonbpresidues.end(), r2) != nonbpresidues.end()) {
            nonbpresidues.erase(remove(nonbpresidues.begin(), nonbpresidues.end(), r2), nonbpresidues.end());
        }
    }

    for (i=0;i<nonbpresidues.size();i++){
        r=nonbpresidues[i];
        fillFragmentLocal(c, r, pdb_entry, fill_chain, ideal_rna["nt"], templateidx, fullatom);
        nonbp_template_idx[r]=templateidx;
    }

    for (r=0;r<pdb_entry.chains[c].residues.size();r++)
    {
        for (a=0;a<pdb_entry.chains[c].residues[r].atoms.size();a++)
        {
            //cout<<pdb_entry.chains[c].residues[r].atoms[a].xyz[0] << "-" << pdb_entry.chains[c].residues[r].atoms[a].xyz[1] << "-" << pdb_entry.chains[c].residues[r].atoms[a].xyz[2]<<endl;
            pdb_entry.chains[c].residues[r].atoms[a].xyz = fill_chain.residues[r].atoms[a].xyz;
            //cout<<pdb_entry.chains[c].residues[r].atoms[a].xyz[0] << "-" << pdb_entry.chains[c].residues[r].atoms[a].xyz[1] << "-" << pdb_entry.chains[c].residues[r].atoms[a].xyz[2]<<endl;
        }
    }

    fill_atom.name.clear();
    fill_atom.xyz.clear();
    fill_residue.resn.clear();
    vector<AtomUnit>   ().swap(fill_residue.atoms);
    vector<ResidueUnit>().swap(fill_chain.residues);

    return true;
}

bool FragRef(ModelUnit &pdb_entry, 
    map<string, map<string, vector<ModelUnit>>> &ideal_rna,
    const vector<pair<double,vector<size_t> > >&bp_vec, map<int, int> &bp_template_idx, map<int, int> &nonbp_template_idx,
    const int option)
{
    size_t c;
    for (c=0;c<pdb_entry.chains.size();c++)
        FragRef(pdb_entry, c, ideal_rna, bp_vec, bp_template_idx, nonbp_template_idx, option);

    return true;
}


