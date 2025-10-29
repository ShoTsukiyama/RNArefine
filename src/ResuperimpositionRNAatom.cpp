
#include "ResuperimpositionRNAatom.hpp"

bool ResuperimpositionRNAatomInBasePair(const int c1, const int r1,
    const int c2, const int r2, const ModelUnit &pdb_entry, ChainUnit &fill_chain,
    map<string, vector<ModelUnit>> &ideal_rna, int &templateidx, double &energy, bool fullatom)
{
    int i, j;
    vector<double> tmp(3,0);
    vector<vector<double> > xyz_list1;
    vector<vector<double> > xyz_list2;
    vector<vector<double> > xyz_list3;
    vector<vector<double> > RotMatix, minRotMatix;  // U
    vector<double> TranVect, minTranVect;  // t
    string key, key1, key2;
    bool stack_prev, stack_next;
    string base, base1, base2, base3, base4, base5, base6, bpkey;
    vector<vector<int>> atomid;
    double rmsd, minrmsd;
    int mintempidx;
    int atom;

    map<string, vector<string>> full_atoms={
        {"A", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'", " N9 ", " C8 ", " N7 ", " C5 ", " C6 ", " N6 ", " N1 ", " C2 ", " N3 ", " C4 "}},
        {"C", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'", " N1 ", " C2 ", " O2 ", " N3 ", " C4 ", " N4 ", " C5 ", " C6 "}},
        {"G", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'", " N9 ", " C8 ", " N7 ", " C5 ", " C6 ", " O6 ", " N1 ", " C2 ", " N2 ", " N3 ", " C4 "}},
        {"U", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'", " N1 ", " C2 ", " O2 ", " N3 ", " C4 ", " O4 ", " C5 ", " C6 "}},
    };

    map<string, vector<string>> bb_atoms={
        {"A", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'"}},
        {"C", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'"}},
        {"G", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'"}},
        {"U", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'"}},
    };
    
    stack_prev=false;
    if(r1-1 >= 0 && r2+1 < pdb_entry.chains[c2].residues.size()) stack_prev=true;
    stack_next=false;
    if(r1+1 < pdb_entry.chains[c1].residues.size() && r2-1 >= 0) stack_next=true;
    
    /*
    cout<<"resuperimposition residue pair "
        <<pdb_entry.chains[c1].residues[r1].resn
        <<pdb_entry.chains[c1].residues[r1].resi
        <<pdb_entry.chains[c1].residues[r1].icode<<" and "
        <<pdb_entry.chains[c2].residues[r2].resn
        <<pdb_entry.chains[c2].residues[r2].resi
        <<pdb_entry.chains[c2].residues[r2].icode<<endl;
    */

    base1=pdb_entry.chains[c1].residues[r1].resn[2];
    base2=pdb_entry.chains[c2].residues[r2].resn[2];
    bpkey=base1+"-"+base2;

    size_t a,a1,a2;
    size_t atomNum=0;
    size_t atomNumall=0;


    //movableflag.resize(6, vector<bool>(23,false));
    for (a=0;a<full_atoms[base1].size();a++){
        if(pdb_entry.chains[c1].residues[r1].atoms[a].input == 1) atomNum+=1;
        atomNumall+=1;
    }
    for (a=0;a<full_atoms[base2].size();a++){
        if(pdb_entry.chains[c2].residues[r2].atoms[a].input == 1) atomNum+=1;
        atomNumall+=1;
    }
        
    if (stack_prev)
    {
        base3=pdb_entry.chains[c1].residues[r1-1].resn[2];
        for (a=0;a<bb_atoms[base3].size();a++){
            if(pdb_entry.chains[c1].residues[r1-1].atoms[a].input == 1) atomNum+=1;
            atomNumall+=1;
        }

        base4=pdb_entry.chains[c2].residues[r2+1].resn[2];
        for (a=0;a<bb_atoms[base4].size();a++){
            if(pdb_entry.chains[c2].residues[r2+1].atoms[a].input == 1) atomNum+=1;
            atomNumall+=1;
        }
    }
    if (stack_next)
    {
        base5=pdb_entry.chains[c1].residues[r1+1].resn[2];
        for (a=0;a<bb_atoms[base5].size();a++){
            if(pdb_entry.chains[c1].residues[r1+1].atoms[a].input == 1) atomNum+=1;
            atomNumall+=1;
        }

        base6=pdb_entry.chains[c2].residues[r2-1].resn[2];
        for (a=0;a<bb_atoms[base6].size();a++){
            if(pdb_entry.chains[c2].residues[r2-1].atoms[a].input == 1) atomNum+=1;
            atomNumall+=1;
        }
    }

    if (ideal_rna.count(bpkey) == 0){
        return false;
    }

    minrmsd=10000000000; mintempidx=-1;
    for(i=0;i<ideal_rna[bpkey].size();i++)
    {
        xyz_list1.clear();
        xyz_list2.clear();
        xyz_list3.clear();

        xyz_list1.assign(atomNum,tmp);
        xyz_list2.assign(atomNum,tmp);
        xyz_list3.assign(atomNumall,tmp);
        a=0;

        for (a1=0;a1<full_atoms[base1].size();a1++)
        {
            if(pdb_entry.chains[c1].residues[r1].atoms[a1].input != 1) continue; 

            xyz_list1[a][0]=ideal_rna[bpkey][i].chains[0].residues[1].atoms[a1].xyz[0];
            xyz_list1[a][1]=ideal_rna[bpkey][i].chains[0].residues[1].atoms[a1].xyz[1];
            xyz_list1[a][2]=ideal_rna[bpkey][i].chains[0].residues[1].atoms[a1].xyz[2];

            xyz_list2[a][0]=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz[0];
            xyz_list2[a][1]=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz[1];
            xyz_list2[a][2]=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz[2];
            a++;
        }

        for (a2=0;a2<full_atoms[base2].size();a2++)
        {
            if(pdb_entry.chains[c2].residues[r2].atoms[a2].input != 1) continue; 

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
            for (a1=0;a1<bb_atoms[base3].size();a1++)
            {
                if(pdb_entry.chains[c1].residues[r1-1].atoms[a1].input != 1) continue; 

                xyz_list1[a][0]=ideal_rna[bpkey][i].chains[0].residues[0].atoms[a1].xyz[0];
                xyz_list1[a][1]=ideal_rna[bpkey][i].chains[0].residues[0].atoms[a1].xyz[1];
                xyz_list1[a][2]=ideal_rna[bpkey][i].chains[0].residues[0].atoms[a1].xyz[2];

                xyz_list2[a][0]=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz[0];
                xyz_list2[a][1]=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz[1];
                xyz_list2[a][2]=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz[2];
                a++;
            }

            for (a2=0;a2<bb_atoms[base4].size();a2++)
            {
                if(pdb_entry.chains[c2].residues[r2+1].atoms[a2].input != 1) continue; 

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
            for (a1=0;a1<bb_atoms[base5].size();a1++)
            {
                if(pdb_entry.chains[c1].residues[r1+1].atoms[a1].input != 1) continue; 

                xyz_list1[a][0]=ideal_rna[bpkey][i].chains[0].residues[2].atoms[a1].xyz[0];
                xyz_list1[a][1]=ideal_rna[bpkey][i].chains[0].residues[2].atoms[a1].xyz[1];
                xyz_list1[a][2]=ideal_rna[bpkey][i].chains[0].residues[2].atoms[a1].xyz[2];

                xyz_list2[a][0]=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz[0];
                xyz_list2[a][1]=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz[1];
                xyz_list2[a][2]=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz[2];
                a++;
            }

            for (a2=0;a2<bb_atoms[base6].size();a2++)
            {
                if(pdb_entry.chains[c2].residues[r2-1].atoms[a2].input != 1) continue;

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
                for (a=0;a<xyz_list3.size();a++)
                {
                    xyz_list3[a].clear();
                }
                xyz_list1.clear();
                xyz_list2.clear();
                xyz_list3.clear();
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

        /*
        atomNum=0;
        for (a=0;a<full_atoms[base1].size();a++)
        {   
            ChangeCoor(ideal_rna[bpkey][i].chains[0].residues[1].atoms[a].xyz,RotMatix,TranVect,xyz_list1[atomNum]);
            tmp[0]+=(pdb_entry.chains[c1].residues[r1].atoms[a].xyz[0]-xyz_list1[atomNum][0]);
            tmp[1]+=(pdb_entry.chains[c1].residues[r1].atoms[a].xyz[1]-xyz_list1[atomNum][1]);
            tmp[2]+=(pdb_entry.chains[c1].residues[r1].atoms[a].xyz[2]-xyz_list1[atomNum][2]);
            atomNum++;
        }
        for (a=0;a<full_atoms[base2].size();a++)
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
        */

        j=0;
        rmsd=0;
        for (a=0;a<full_atoms[base1].size();a++)
        {
            //if(pdb_entry.chains[c1].residues[r1].atoms[a].input != 1) continue; 
            ChangeCoor(ideal_rna[bpkey][i].chains[0].residues[1].atoms[a].xyz,RotMatix, TranVect, xyz_list3[j]);
            tmp[0]=(pdb_entry.chains[c1].residues[r1].atoms[a].xyz[0]-xyz_list3[j][0]);
            tmp[1]=(pdb_entry.chains[c1].residues[r1].atoms[a].xyz[1]-xyz_list3[j][1]);
            tmp[2]=(pdb_entry.chains[c1].residues[r1].atoms[a].xyz[2]-xyz_list3[j][2]);

            j++;
            rmsd+=sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
        }
        for (a=0;a<full_atoms[base2].size();a++)
        {
            //if(pdb_entry.chains[c2].residues[r2].atoms[a].input != 1) continue; 
            ChangeCoor(ideal_rna[bpkey][i].chains[0].residues[4].atoms[a].xyz,RotMatix,TranVect,xyz_list3[j]);
            tmp[0]=(pdb_entry.chains[c2].residues[r2].atoms[a].xyz[0]-xyz_list3[j][0]);
            tmp[1]=(pdb_entry.chains[c2].residues[r2].atoms[a].xyz[1]-xyz_list3[j][1]);
            tmp[2]=(pdb_entry.chains[c2].residues[r2].atoms[a].xyz[2]-xyz_list3[j][2]);

            j++;
            rmsd+=sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
        }

        if (stack_prev)
        {
            for (a=0;a<bb_atoms[base3].size();a++)
            {
                //if(pdb_entry.chains[c1].residues[r1-1].atoms[a].input != 1) continue; 
                ChangeCoor(ideal_rna[bpkey][i].chains[0].residues[0].atoms[a].xyz,RotMatix, TranVect, xyz_list3[j]);
                tmp[0]=(pdb_entry.chains[c1].residues[r1-1].atoms[a].xyz[0]-xyz_list3[j][0]);
                tmp[1]=(pdb_entry.chains[c1].residues[r1-1].atoms[a].xyz[1]-xyz_list3[j][1]);
                tmp[2]=(pdb_entry.chains[c1].residues[r1-1].atoms[a].xyz[2]-xyz_list3[j][2]);

                j++;
                rmsd+=sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
            }
            for (a=0;a<bb_atoms[base4].size();a++)
            {
                //if(pdb_entry.chains[c2].residues[r2+1].atoms[a].input != 1) continue; 
                ChangeCoor(ideal_rna[bpkey][i].chains[0].residues[5].atoms[a].xyz,RotMatix,TranVect,xyz_list3[j]);
                tmp[0]=(pdb_entry.chains[c2].residues[r2+1].atoms[a].xyz[0]-xyz_list3[j][0]);
                tmp[1]=(pdb_entry.chains[c2].residues[r2+1].atoms[a].xyz[1]-xyz_list3[j][1]);
                tmp[2]=(pdb_entry.chains[c2].residues[r2+1].atoms[a].xyz[2]-xyz_list3[j][2]);

                j++;
                rmsd+=sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
            }
        }

        if (stack_next)
        {
            for (a=0;a<bb_atoms[base5].size();a++)
            {
                //if(pdb_entry.chains[c1].residues[r1+1].atoms[a].input != 1) continue; 
                ChangeCoor(ideal_rna[bpkey][i].chains[0].residues[2].atoms[a].xyz,RotMatix, TranVect, xyz_list3[j]);
                tmp[0]=(pdb_entry.chains[c1].residues[r1+1].atoms[a].xyz[0]-xyz_list3[j][0]);
                tmp[1]=(pdb_entry.chains[c1].residues[r1+1].atoms[a].xyz[1]-xyz_list3[j][1]);
                tmp[2]=(pdb_entry.chains[c1].residues[r1+1].atoms[a].xyz[2]-xyz_list3[j][2]);

                j++;
                rmsd+=sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
            }
            for (a=0;a<bb_atoms[base6].size();a++)
            {
                //if(pdb_entry.chains[c2].residues[r2-1].atoms[a].input != 1) continue; 
                ChangeCoor(ideal_rna[bpkey][i].chains[0].residues[3].atoms[a].xyz,RotMatix,TranVect,xyz_list3[j]);
                tmp[0]=(pdb_entry.chains[c2].residues[r2-1].atoms[a].xyz[0]-xyz_list3[j][0]);
                tmp[1]=(pdb_entry.chains[c2].residues[r2-1].atoms[a].xyz[1]-xyz_list3[j][1]);
                tmp[2]=(pdb_entry.chains[c2].residues[r2-1].atoms[a].xyz[2]-xyz_list3[j][2]);
                j++;
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
    energy+=minrmsd;

    if(fullatom){
        for (a=0;a<full_atoms[base1].size();a++)
        {
            ChangeCoor(ideal_rna[bpkey][mintempidx].chains[0].residues[1].atoms[a].xyz,
                minRotMatix, minTranVect, fill_chain.residues[r1].atoms[a].xyz);
        }

        for (a=0;a<full_atoms[base2].size();a++)
        {
            ChangeCoor(ideal_rna[bpkey][mintempidx].chains[0].residues[4].atoms[a].xyz,
                minRotMatix, minTranVect, fill_chain.residues[r2].atoms[a].xyz);
        }
    }else{
        /* eventually assign the coordinate */
        for (a=0;a<full_atoms[base1].size();a++)
        {
            if(pdb_entry.chains[c1].residues[r1].atoms[a].movable == 0) continue;

            ChangeCoor(ideal_rna[bpkey][mintempidx].chains[0].residues[1].atoms[a].xyz,
                minRotMatix, minTranVect, fill_chain.residues[r1].atoms[a].xyz);
        }

        for (a=0;a<full_atoms[base2].size();a++)
        {
            if(pdb_entry.chains[c2].residues[r2].atoms[a].movable == 0) continue;

            ChangeCoor(ideal_rna[bpkey][mintempidx].chains[0].residues[4].atoms[a].xyz,
                minRotMatix, minTranVect, fill_chain.residues[r2].atoms[a].xyz);
        }
    }
        
    /* clean up */
    for (a=0;a<xyz_list1.size();a++)
    {
        xyz_list1[a].clear();
        xyz_list2[a].clear();
    }
    for (a=0;a<xyz_list3.size();a++)
    {
        xyz_list3[a].clear();
    }
    xyz_list1.clear();
    xyz_list2.clear();
    xyz_list3.clear();
    TranVect.clear();
    for (a=0;a<3;a++){
        RotMatix[a].clear();
    }
    RotMatix.clear();
    return true;
}

bool ResuperimpositionRNAatomLocal(const int c, const int r,
    const ModelUnit &pdb_entry, ChainUnit &fill_chain, map<string, vector<ModelUnit>> &ideal_rna, int &templateidx, double &energy, bool fullatom)
{
    int i,j;
    vector<double> tmp(3,0);
    vector<vector<double> > xyz_list1;
    vector<vector<double> > xyz_list2;
    vector<vector<double> > xyz_list3;
    vector<vector<double> > RotMatix, minRotMatix;  // U
    vector<double> TranVect, minTranVect;  // t
    string key, key1, key2;
    bool stack_prev, stack_next;
    string base;
    double rmsd, minrmsd;
    int mintempidx;
    string base1, base2, base3;
    int atom;

    map<string, vector<string>> full_atoms={
        {"A", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'", " N9 ", " C8 ", " N7 ", " C5 ", " C6 ", " N6 ", " N1 ", " C2 ", " N3 ", " C4 "}},
        {"C", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'", " N1 ", " C2 ", " O2 ", " N3 ", " C4 ", " N4 ", " C5 ", " C6 "}},
        {"G", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'", " N9 ", " C8 ", " N7 ", " C5 ", " C6 ", " O6 ", " N1 ", " C2 ", " N2 ", " N3 ", " C4 "}},
        {"U", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'", " N1 ", " C2 ", " O2 ", " N3 ", " C4 ", " O4 ", " C5 ", " C6 "}},
    };

    map<string, vector<string>> bb_atoms={
        {"A", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'"}},
        {"C", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'"}},
        {"G", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'"}},
        {"U", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'"}},
    };
    
    stack_prev=false;
    stack_next=false;
    if(r-1 >= 0) stack_prev=true;
    if(r+1 < pdb_entry.chains[c].residues.size()) stack_next=true;
    
    /*
    cout<<"resuperimposing residue pair "
        <<pdb_entry.chains[c].residues[r].resn
        <<pdb_entry.chains[c].residues[r].resi
        <<pdb_entry.chains[c].residues[r].icode<<endl;
    */

    base2=pdb_entry.chains[c].residues[r].resn[2];

    size_t a,a1,a2;
    size_t atomNum=0;
    size_t atomNumall=0;

    for (a=0;a<full_atoms[base2].size();a++){
        if(pdb_entry.chains[c].residues[r].atoms[a].input == 1) atomNum+=1;
        atomNumall += 1;
    }
    if (stack_prev)
    {
        base1=pdb_entry.chains[c].residues[r-1].resn[2];
        for (a=0;a<bb_atoms[base1].size();a++){
            if(pdb_entry.chains[c].residues[r-1].atoms[a].input == 1) atomNum+=1;
            atomNumall += 1;
        }
    }

    if (stack_next)
    {
        base3=pdb_entry.chains[c].residues[r+1].resn[2];
        for (a=0;a<bb_atoms[base3].size();a++){
            if(pdb_entry.chains[c].residues[r+1].atoms[a].input == 1) atomNum+=1;
            atomNumall += 1;
        }
    }

    if (ideal_rna.count(base2) == 0) return false;

    minrmsd=10000000000; mintempidx=-1;
    for(i=0;i<ideal_rna[base2].size();i++)
    {
        xyz_list1.clear();
        xyz_list2.clear();
        xyz_list3.clear();

        xyz_list1.assign(atomNum,tmp);
        xyz_list2.assign(atomNum,tmp);
        xyz_list3.assign(atomNumall,tmp);
        a=0;
        for (a1=0;a1<full_atoms[base2].size();a1++)
        {
            if(pdb_entry.chains[c].residues[r].atoms[a1].input != 1)  continue;
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
            for (a1=0;a1<bb_atoms[base1].size();a1++)
            {
                if(pdb_entry.chains[c].residues[r-1].atoms[a1].input != 1)  continue;
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
            for (a1=0;a1<bb_atoms[base3].size();a1++)
            {
                if(pdb_entry.chains[c].residues[r+1].atoms[a1].input != 1)  continue;
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
                for (a=0;a<xyz_list3.size();a++)
                {
                    xyz_list3[a].clear();
                }
                xyz_list1.clear();
                xyz_list2.clear();
                xyz_list3.clear();
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
        /*
        atomNum=0;
        for (a=0;a<full_atoms[base2].size();a++)
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
        */

        j=0;
        rmsd=0;
        for (a=0;a<full_atoms[base2].size();a++)
        {
            //if(pdb_entry.chains[c].residues[r].atoms[a].input != 1)  continue;
            ChangeCoor(ideal_rna[base2][i].chains[0].residues[1].atoms[a].xyz,RotMatix, TranVect, xyz_list3[j]);
            tmp[0]=(pdb_entry.chains[c].residues[r].atoms[a].xyz[0]-xyz_list3[j][0]);
            tmp[1]=(pdb_entry.chains[c].residues[r].atoms[a].xyz[1]-xyz_list3[j][1]);
            tmp[2]=(pdb_entry.chains[c].residues[r].atoms[a].xyz[2]-xyz_list3[j][2]);
            j++;
            rmsd+=sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
        }

        if (stack_prev)
        {
            for (a=0;a<bb_atoms[base1].size();a++)
            {
                //if(pdb_entry.chains[c].residues[r-1].atoms[a].input != 1)  continue;
                ChangeCoor(ideal_rna[base2][i].chains[0].residues[0].atoms[a].xyz,RotMatix, TranVect, xyz_list3[j]);
                tmp[0]=(pdb_entry.chains[c].residues[r-1].atoms[a].xyz[0]-xyz_list3[j][0]);
                tmp[1]=(pdb_entry.chains[c].residues[r-1].atoms[a].xyz[1]-xyz_list3[j][1]);
                tmp[2]=(pdb_entry.chains[c].residues[r-1].atoms[a].xyz[2]-xyz_list3[j][2]);
                j++;
                rmsd+=sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
            }
        }

        if (stack_next)
        {
            for (a=0;a<bb_atoms[base3].size();a++)
            {
                //if(pdb_entry.chains[c].residues[r+1].atoms[a].input != 1)  continue;
                ChangeCoor(ideal_rna[base2][i].chains[0].residues[2].atoms[a].xyz,RotMatix, TranVect, xyz_list3[j]);
                tmp[0]=(pdb_entry.chains[c].residues[r+1].atoms[a].xyz[0]-xyz_list3[j][0]);
                tmp[1]=(pdb_entry.chains[c].residues[r+1].atoms[a].xyz[1]-xyz_list3[j][1]);
                tmp[2]=(pdb_entry.chains[c].residues[r+1].atoms[a].xyz[2]-xyz_list3[j][2]);
                j++;
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
    energy+=minrmsd;

    if(fullatom){
        for (a=0;a<full_atoms[base2].size();a++)
        {
            ChangeCoor(ideal_rna[base2][mintempidx].chains[0].residues[1].atoms[a].xyz,
                minRotMatix, minTranVect, fill_chain.residues[r].atoms[a].xyz);
        }
    }else{
        for (a=0;a<full_atoms[base2].size();a++)
        {
            if(pdb_entry.chains[c].residues[r].atoms[a].movable == 0) continue;

            ChangeCoor(ideal_rna[base2][mintempidx].chains[0].residues[1].atoms[a].xyz,
                minRotMatix, minTranVect, fill_chain.residues[r].atoms[a].xyz);
        }
    }
        
    /* clean up */
    for (a=0;a<xyz_list1.size();a++)
    {
        xyz_list1[a].clear();
        xyz_list2[a].clear();
    }
    for (a=0;a<xyz_list3.size();a++)
    {
        xyz_list3[a].clear();
    }
    xyz_list1.clear();
    xyz_list2.clear();
    xyz_list3.clear();
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

bool ResuperimpositionRNAatom(ModelUnit &pdb_entry, const size_t c,
    map<string, map<string, vector<ModelUnit>>> &ideal_rna,
    const vector<pair<double,vector<size_t>>>&bp_vec, map<int, int> &bp_template_idx, map<int, int> &nonbp_template_idx,
    double &energy, int &spnum,
    const int option)
{
    vector<string> missing_atom_vec;
    vector<size_t> nonstd_atom_vec;
    size_t r,i,j;
    vector<size_t> remove_residue_vec;
    vector<bool> target_residues(pdb_entry.chains[c].residues.size(),false);
    vector<int> nonbpresidues;
    bool fullatom;
    int templateidx;
    bool bpprev, bpnext;
    string looptype;

    size_t bp;
    int r1,r2,c1,c2;

    ChainUnit fill_chain;
    ResidueUnit fill_residue;
    AtomUnit fill_atom;
    size_t a;

    if(option >= 1){
        if(option == 1) fullatom=false;
        else if(option == 2) fullatom=true;

        for (r=0;r<pdb_entry.chains[c].residues.size();r++)
        {
            nonbpresidues.push_back(r);
            fill_chain.residues.push_back(fill_residue);
            for (a=0;a<pdb_entry.chains[c].residues[r].atoms.size();a++)
            {
                fill_chain.residues[r].atoms.push_back(fill_atom);
                fill_chain.residues[r].atoms[a].name=pdb_entry.chains[c].residues[r].atoms[a].name;

                fill_chain.residues[r].atoms[a].xyz.assign(3,0);
                fill_chain.residues[r].atoms[a].xyz[0]=pdb_entry.chains[c].residues[r].atoms[a].xyz[0];
                fill_chain.residues[r].atoms[a].xyz[1]=pdb_entry.chains[c].residues[r].atoms[a].xyz[1];
                fill_chain.residues[r].atoms[a].xyz[2]=pdb_entry.chains[c].residues[r].atoms[a].xyz[2];
                fill_chain.residues[r].atoms[a].movable=pdb_entry.chains[c].residues[r].atoms[a].movable;
                fill_chain.residues[r].atoms[a].input=pdb_entry.chains[c].residues[r].atoms[a].input;

                if(fill_chain.residues[r].atoms[a].movable == 1) target_residues[r]=true;
            }
        }
  
        map<string, int> resnorder={{"  A  U", 0}, {"  U  A", 1}, {"  G  C", 0}, {"  C  G", 1}, {"  G  U", 0}, {"  U  G", 1}};
        for (bp=0;bp<bp_vec.size();bp++)
        {
            //if (bp_vec[bp].second[4]) continue; // Watson-Crick pair only
            c1= bp_vec[bp].second[0];
            r1= bp_vec[bp].second[1];
            c2= bp_vec[bp].second[2];
            r2= bp_vec[bp].second[3];

            auto it1 = std::find(target_residues.begin(), target_residues.end(), r1);
            auto it2 = std::find(target_residues.begin(), target_residues.end(), r2);

            if((option == 2) || ((option == 1) && ((c1==c && it1 != target_residues.end()) || (c2==c && it2 != target_residues.end())))){
                if(resnorder[pdb_entry.chains[c1].residues[r1].resn+pdb_entry.chains[c2].residues[r2].resn] == 0){
                    ResuperimpositionRNAatomInBasePair(c1,r1,c2,r2, pdb_entry, fill_chain, ideal_rna["bp"], templateidx, energy, fullatom);
                    spnum+=1;
                    bp_template_idx[bp]=templateidx;
                }
                
                if(resnorder[pdb_entry.chains[c1].residues[r1].resn+pdb_entry.chains[c2].residues[r2].resn] == 1){
                    ResuperimpositionRNAatomInBasePair(c2,r2,c1,r1, pdb_entry, fill_chain, ideal_rna["bp"], templateidx, energy, fullatom);
                    spnum+=1;
                    bp_template_idx[bp]=templateidx;
                }
            }

            if (find(nonbpresidues.begin(), nonbpresidues.end(), r1) != nonbpresidues.end()) {
                nonbpresidues.erase(remove(nonbpresidues.begin(), nonbpresidues.end(), r1), nonbpresidues.end());
            }
            if (find(nonbpresidues.begin(), nonbpresidues.end(), r2) != nonbpresidues.end()) {
                nonbpresidues.erase(remove(nonbpresidues.begin(), nonbpresidues.end(), r2), nonbpresidues.end());
            }
        }

        for (i=0;i<nonbpresidues.size();i++)
        {
            r=nonbpresidues[i];

            //bpprev=true;
            //bpnext=true;
            //if (find(nonbpresidues.begin(), nonbpresidues.end(), r-1) != nonbpresidues.end()) bpprev=false;
            //if (find(nonbpresidues.begin(), nonbpresidues.end(), r+1) != nonbpresidues.end()) bpnext=false;

            auto it = std::find(target_residues.begin(), target_residues.end(), r);

            if((option == 2) || ((option == 1) && it != target_residues.end())){

                //if(bpprev && bpnext) looptype="oneloop";
                //else if((!bpprev) && bpnext) looptype="closing";
                //else if(bpprev && (!bpnext)) looptype="opening";
                //else looptype="loop";

                looptype="loop";
                
                ResuperimpositionRNAatomLocal(c, r, pdb_entry, fill_chain, ideal_rna[looptype], templateidx, energy, fullatom);
                spnum+=1;
                nonbp_template_idx[r]=templateidx;
            }
        }

        if(option == 1){
            for (r=0;r<pdb_entry.chains[c].residues.size();r++)
            {
                for (a=0;a<fill_chain.residues[r].atoms.size();a++)
                {
                    if(pdb_entry.chains[c].residues[r].atoms[a].movable) pdb_entry.chains[c].residues[r].atoms[a]=fill_chain.residues[r].atoms[a];
                }
            }
        }else if(option == 2){
            for (r=0;r<pdb_entry.chains[c].residues.size();r++)
            {
                for (a=0;a<pdb_entry.chains[c].residues[r].atoms.size();a++)
                {
                    pdb_entry.chains[c].residues[r].atoms[a]=fill_chain.residues[r].atoms[a];
                }
            }
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

    return true;
}

bool ResuperimpositionRNAatom(ModelUnit &pdb_entry, 
    map<string, map<string, vector<ModelUnit>>> &ideal_rna,
    const vector<pair<double,vector<size_t> > >&bp_vec, map<int, int> &bp_template_idx, map<int, int> &nonbp_template_idx,
    double &energy, int &spnum,
    const int option)
{
    energy=0; spnum=0;
    size_t c;
    for (c=0;c<pdb_entry.chains.size();c++)
        ResuperimpositionRNAatom(pdb_entry, c, ideal_rna, bp_vec, bp_template_idx, nonbp_template_idx, energy, spnum, option);

    return true;
}

