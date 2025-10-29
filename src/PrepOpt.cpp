
#include "PrepOpt.hpp"

void ChangeCoor(const vector<double> &v1, const vector<vector<double> > &rotation_matrix, const vector<double> &translate_vec,Point3s &v2)
{
	v2.x=v1[0]*rotation_matrix[0][0]+v1[1]*rotation_matrix[1][0]+v1[2]*rotation_matrix[2][0]+translate_vec[0];
	v2.y=v1[0]*rotation_matrix[0][1]+v1[1]*rotation_matrix[1][1]+v1[2]*rotation_matrix[2][1]+translate_vec[1];
	v2.z=v1[0]*rotation_matrix[0][2]+v1[1]*rotation_matrix[1][2]+v1[2]*rotation_matrix[2][2]+translate_vec[2];
}

void prepopt(const int c1, const int r1,
    const int c2, const int r2, const ModelUnit &pdb_entry,
    map<string, vector<ModelUnit>> &ideal_rna, bool stack_prev, bool stack_next, 
    vector<vector<double>> minRotMatix, vector<double> minTranVect, int mintempidx,
    FragmentStructure &superimposedstr, FragmentStructure &inputtedstr, vector<vector<bool>> &inputtedstrflag)
{
    int i;
    int atom;
    int idx;
    string base0, base1, base2, base3, base4, base5, bpkey;
    FragmentAtoms tmp;
    tmp.coord.x = 0; tmp.coord.y = 0; tmp.coord.z = 0;

    base1=pdb_entry.chains[c1].residues[r1].resn[2];
    base4=pdb_entry.chains[c2].residues[r2].resn[2];

    vector<int> residueidx;
    if (stack_prev && stack_next) residueidx={1,4,0,5,2,3};
    else if(!stack_prev && stack_next) residueidx={0,3,1,2};
    else if(stack_prev && !stack_next) residueidx={1,2,0,3};
    else residueidx={0,1};

    bpkey=base1+"-"+base4;

    if (stack_prev && stack_next){
        superimposedstr.residues.resize(6);
        inputtedstr.residues.resize(6);
        inputtedstrflag.resize(6);
    }else if((!stack_prev && stack_next) || (stack_prev && !stack_next)){
        superimposedstr.residues.resize(4);
        inputtedstr.residues.resize(4);
        inputtedstrflag.resize(4);
    }else{
        superimposedstr.residues.resize(2);
        inputtedstr.residues.resize(2);
        inputtedstrflag.resize(2);
    }

    int j=0;
    base1=pdb_entry.chains[c1].residues[r1].resn[2];
    idx=residueidx[j];
    superimposedstr.residues[idx].restype=base1;
    superimposedstr.residues[idx].atoms.resize(full_atoms[base1].size(), tmp);

    inputtedstr.residues[idx].restype=base1;
    inputtedstr.residues[idx].atoms.resize(full_atoms[base1].size(), tmp);

    inputtedstrflag[idx].resize(full_atoms[base1].size(), false);
    for(i=0;i<full_atoms[base1].size();i++){
        ChangeCoor(ideal_rna[bpkey][mintempidx].chains[0].residues[1].atoms[i].xyz,minRotMatix,minTranVect,superimposedstr.residues[idx].atoms[i].coord);
    }
    j += 1;

    base4=pdb_entry.chains[c2].residues[r2].resn[2];
    idx=residueidx[j];
    
    superimposedstr.residues[idx].restype=base4;
    superimposedstr.residues[idx].atoms.resize(full_atoms[base4].size(), tmp);

    inputtedstr.residues[idx].restype=base4;
    inputtedstr.residues[idx].atoms.resize(full_atoms[base4].size(), tmp);

    inputtedstrflag[idx].resize(full_atoms[base4].size(), false);
    for(i=0;i<full_atoms[base4].size();i++){
        ChangeCoor(ideal_rna[bpkey][mintempidx].chains[0].residues[4].atoms[i].xyz,minRotMatix,minTranVect,superimposedstr.residues[idx].atoms[i].coord);
    }
    j += 1;

    if(stack_prev){
        base0=pdb_entry.chains[c1].residues[r1-1].resn[2];
        idx=residueidx[j];
        
        superimposedstr.residues[idx].restype=base0;
        superimposedstr.residues[idx].atoms.resize(bb_atoms[base0].size(), tmp);

        inputtedstr.residues[idx].restype=base0;
        inputtedstr.residues[idx].atoms.resize(bb_atoms[base0].size(), tmp);

        inputtedstrflag[idx].resize(bb_atoms[base0].size(), false);
        for(i=0;i<bb_atoms[base0].size();i++){
            ChangeCoor(ideal_rna[bpkey][mintempidx].chains[0].residues[0].atoms[i].xyz,minRotMatix,minTranVect,superimposedstr.residues[idx].atoms[i].coord);    
        }
        j += 1;

        base5=pdb_entry.chains[c2].residues[r2+1].resn[2];
        idx=residueidx[j];
        superimposedstr.residues[idx].restype=base5;
        superimposedstr.residues[idx].atoms.resize(bb_atoms[base5].size(), tmp);

        inputtedstr.residues[idx].restype=base5;
        inputtedstr.residues[idx].atoms.resize(bb_atoms[base5].size(), tmp);

        inputtedstrflag[idx].resize(bb_atoms[base5].size(), false);
        for(i=0;i<bb_atoms[base5].size();i++){
            ChangeCoor(ideal_rna[bpkey][mintempidx].chains[0].residues[5].atoms[i].xyz,minRotMatix,minTranVect,superimposedstr.residues[idx].atoms[i].coord);
        }
        j += 1;
    }

    if(stack_next){
        base2=pdb_entry.chains[c1].residues[r1+1].resn[2];
        idx=residueidx[j];
        superimposedstr.residues[idx].restype=base2;
        superimposedstr.residues[idx].atoms.resize(bb_atoms[base2].size(), tmp);

        inputtedstr.residues[idx].restype=base2;
        inputtedstr.residues[idx].atoms.resize(bb_atoms[base2].size(), tmp);

        inputtedstrflag[idx].resize(bb_atoms[base2].size(), false);
        for(i=0;i<bb_atoms[base2].size();i++){
            ChangeCoor(ideal_rna[bpkey][mintempidx].chains[0].residues[2].atoms[i].xyz,minRotMatix,minTranVect, superimposedstr.residues[idx].atoms[i].coord);
        }
        j += 1;

        base3=pdb_entry.chains[c2].residues[r2-1].resn[2];
        idx=residueidx[j];
        superimposedstr.residues[idx].restype=base3;
        superimposedstr.residues[idx].atoms.resize(bb_atoms[base3].size(), tmp);

        inputtedstr.residues[idx].restype=base3;
        inputtedstr.residues[idx].atoms.resize(bb_atoms[base3].size(), tmp);

        inputtedstrflag[idx].resize(bb_atoms[base3].size(), false);
        for(i=0;i<bb_atoms[base3].size();i++){
            ChangeCoor(ideal_rna[bpkey][mintempidx].chains[0].residues[3].atoms[i].xyz,minRotMatix,minTranVect, superimposedstr.residues[idx].atoms[i].coord);
        }
        j += 1;
    }

    j = 0;
    for (int a=0;a<pdb_entry.chains[c1].residues[r1].atoms.size();a++)
    {
        auto it = find(full_atoms[base1].begin(), full_atoms[base1].end(), pdb_entry.chains[c1].residues[r1].atoms[a].name);
        
        if(it != full_atoms[base1].end()){
            atom = distance(full_atoms[base1].begin(), it);
            idx=residueidx[j];
            inputtedstr.residues[idx].atoms[atom].coord.x = pdb_entry.chains[c1].residues[r1].atoms[a].xyz[0];
            inputtedstr.residues[idx].atoms[atom].coord.y = pdb_entry.chains[c1].residues[r1].atoms[a].xyz[1];
            inputtedstr.residues[idx].atoms[atom].coord.z = pdb_entry.chains[c1].residues[r1].atoms[a].xyz[2];
            inputtedstrflag[idx][atom] = true;
        }
    }
    j += 1;

    for (int a=0;a<pdb_entry.chains[c2].residues[r2].atoms.size();a++)
    {
        auto it = find(full_atoms[base4].begin(), full_atoms[base4].end(), pdb_entry.chains[c2].residues[r2].atoms[a].name);

        if(it != full_atoms[base4].end()){
            atom = distance(full_atoms[base4].begin(), it);
            idx=residueidx[j];
            inputtedstr.residues[idx].atoms[atom].coord.x = pdb_entry.chains[c2].residues[r2].atoms[a].xyz[0];
            inputtedstr.residues[idx].atoms[atom].coord.y = pdb_entry.chains[c2].residues[r2].atoms[a].xyz[1];
            inputtedstr.residues[idx].atoms[atom].coord.z = pdb_entry.chains[c2].residues[r2].atoms[a].xyz[2];
            inputtedstrflag[idx][atom] = true;
            
        }
    }
    j += 1;

    if(stack_prev){
        for (int a=0;a<pdb_entry.chains[c1].residues[r1-1].atoms.size();a++)
        {
            auto it = find(bb_atoms[base0].begin(), bb_atoms[base0].end(), pdb_entry.chains[c1].residues[r1-1].atoms[a].name);
            
            if(it != bb_atoms[base0].end()){
                atom = distance(bb_atoms[base0].begin(), it);
                idx=residueidx[j];
                inputtedstr.residues[idx].atoms[atom].coord.x = pdb_entry.chains[c1].residues[r1-1].atoms[a].xyz[0];
                inputtedstr.residues[idx].atoms[atom].coord.y = pdb_entry.chains[c1].residues[r1-1].atoms[a].xyz[1];
                inputtedstr.residues[idx].atoms[atom].coord.z = pdb_entry.chains[c1].residues[r1-1].atoms[a].xyz[2];
                inputtedstrflag[idx][atom] = true;
            }
        }
        j += 1;

        for (int a=0;a<pdb_entry.chains[c2].residues[r2+1].atoms.size();a++)
        {
            auto it = find(bb_atoms[base5].begin(), bb_atoms[base5].end(), pdb_entry.chains[c2].residues[r2+1].atoms[a].name);
            
            if(it != bb_atoms[base5].end()){
                atom = distance(bb_atoms[base5].begin(), it);
                idx=residueidx[j];
                inputtedstr.residues[idx].atoms[atom].coord.x = pdb_entry.chains[c2].residues[r2+1].atoms[a].xyz[0];
                inputtedstr.residues[idx].atoms[atom].coord.y = pdb_entry.chains[c2].residues[r2+1].atoms[a].xyz[1];
                inputtedstr.residues[idx].atoms[atom].coord.z = pdb_entry.chains[c2].residues[r2+1].atoms[a].xyz[2];
                inputtedstrflag[idx][atom] = true;
            }
        }
        j += 1;
    }

    if(stack_next){
        base2=pdb_entry.chains[c1].residues[r1+1].resn[2];
        for (int a=0;a<pdb_entry.chains[c1].residues[r1+1].atoms.size();a++)
        {
            auto it = find(bb_atoms[base2].begin(), bb_atoms[base2].end(), pdb_entry.chains[c1].residues[r1+1].atoms[a].name);
            
            if(it != bb_atoms[base2].end()){
                atom = distance(bb_atoms[base2].begin(), it);
                idx=residueidx[j];
                inputtedstr.residues[idx].atoms[atom].coord.x = pdb_entry.chains[c1].residues[r1+1].atoms[a].xyz[0];
                inputtedstr.residues[idx].atoms[atom].coord.y = pdb_entry.chains[c1].residues[r1+1].atoms[a].xyz[1];
                inputtedstr.residues[idx].atoms[atom].coord.z = pdb_entry.chains[c1].residues[r1+1].atoms[a].xyz[2];
                inputtedstrflag[idx][atom] = true;
            }
        }
        j += 1;

        base3=pdb_entry.chains[c2].residues[r2-1].resn[2];
        for (int a=0;a<pdb_entry.chains[c2].residues[r2-1].atoms.size();a++)
        {
            auto it = find(bb_atoms[base3].begin(), bb_atoms[base3].end(), pdb_entry.chains[c2].residues[r2-1].atoms[a].name);
            
            if(it != bb_atoms[base3].end()){
                atom = distance(bb_atoms[base3].begin(), it);
                idx=residueidx[j];
                inputtedstr.residues[idx].atoms[atom].coord.x = pdb_entry.chains[c2].residues[r2-1].atoms[a].xyz[0];
                inputtedstr.residues[idx].atoms[atom].coord.y = pdb_entry.chains[c2].residues[r2-1].atoms[a].xyz[1];
                inputtedstr.residues[idx].atoms[atom].coord.z = pdb_entry.chains[c2].residues[r2-1].atoms[a].xyz[2];
                inputtedstrflag[idx][atom] = true;
            }
        }
        j += 1;
    }
}

void prepoptloop(const int c, const int r, const ModelUnit &pdb_entry,
    map<string, vector<ModelUnit>> &ideal_rna, bool stack_prev, bool stack_next, 
    vector<vector<double>> minRotMatix, vector<double> minTranVect, int mintempidx,
    FragmentStructure &superimposedstr, FragmentStructure &inputtedstr, vector<vector<bool>> &inputtedstrflag)
{
    int i;
    int atom;
    int idx;
    string base0, base1, base2;
    FragmentAtoms tmp;
    tmp.coord.x = 0; tmp.coord.y = 0; tmp.coord.z = 0;

    base1=pdb_entry.chains[c].residues[r].resn[2];

    vector<int> residueidx;
    if (stack_prev && stack_next) residueidx={1,0,2};
    else if(!stack_prev && stack_next) residueidx={0,1};
    else if(stack_prev && !stack_next) residueidx={1,0};
    else residueidx={0};

    if (stack_prev && stack_next){
        superimposedstr.residues.resize(3);
        inputtedstr.residues.resize(3);
        inputtedstrflag.resize(3);
    }else if((!stack_prev && stack_next) || (stack_prev && !stack_next)){
        superimposedstr.residues.resize(2);
        inputtedstr.residues.resize(2);
        inputtedstrflag.resize(2);
    }else{
        superimposedstr.residues.resize(1);
        inputtedstr.residues.resize(1);
        inputtedstrflag.resize(1);
    }

    int j=0;
    base1=pdb_entry.chains[c].residues[r].resn[2];
    idx=residueidx[j];
    superimposedstr.residues[idx].restype=base1;
    superimposedstr.residues[idx].atoms.resize(full_atoms[base1].size(), tmp);

    inputtedstr.residues[idx].restype=base1;
    inputtedstr.residues[idx].atoms.resize(full_atoms[base1].size(), tmp);

    inputtedstrflag[idx].resize(full_atoms[base1].size(), false);
    for(i=0;i<full_atoms[base1].size();i++){
        ChangeCoor(ideal_rna[base1][mintempidx].chains[0].residues[1].atoms[i].xyz,minRotMatix,minTranVect, superimposedstr.residues[idx].atoms[i].coord);
    }
    j += 1;

    if(stack_prev)
    {
        base0=pdb_entry.chains[c].residues[r-1].resn[2];
        idx=residueidx[j];
        
        superimposedstr.residues[idx].restype=base0;
        superimposedstr.residues[idx].atoms.resize(bb_atoms[base0].size(), tmp);

        inputtedstr.residues[idx].restype=base0;
        inputtedstr.residues[idx].atoms.resize(bb_atoms[base0].size(), tmp);

        inputtedstrflag[idx].resize(bb_atoms[base0].size(), false);
        for(i=0;i<bb_atoms[base0].size();i++){
            ChangeCoor(ideal_rna[base1][mintempidx].chains[0].residues[0].atoms[i].xyz,minRotMatix,minTranVect,superimposedstr.residues[idx].atoms[i].coord);    
        }
        j += 1;
    }

    if(stack_next)
    {
        base2=pdb_entry.chains[c].residues[r+1].resn[2];
        idx=residueidx[j];
        superimposedstr.residues[idx].restype=base2;
        superimposedstr.residues[idx].atoms.resize(bb_atoms[base2].size(), tmp);

        inputtedstr.residues[idx].restype=base2;
        inputtedstr.residues[idx].atoms.resize(bb_atoms[base2].size(), tmp);

        inputtedstrflag[idx].resize(bb_atoms[base2].size(), false);
        for(i=0;i<bb_atoms[base2].size();i++){
            ChangeCoor(ideal_rna[base1][mintempidx].chains[0].residues[2].atoms[i].xyz,minRotMatix,minTranVect, superimposedstr.residues[idx].atoms[i].coord);
        }
        j += 1;
    }

    j = 0;
    for (int a=0;a<pdb_entry.chains[c].residues[r].atoms.size();a++)
    {
        auto it = find(full_atoms[base1].begin(), full_atoms[base1].end(), pdb_entry.chains[c].residues[r].atoms[a].name);
        
        if(it != full_atoms[base1].end())
        {
            atom = distance(full_atoms[base1].begin(), it);
            idx=residueidx[j];
            inputtedstr.residues[idx].atoms[atom].coord.x = pdb_entry.chains[c].residues[r].atoms[a].xyz[0];
            inputtedstr.residues[idx].atoms[atom].coord.y = pdb_entry.chains[c].residues[r].atoms[a].xyz[1];
            inputtedstr.residues[idx].atoms[atom].coord.z = pdb_entry.chains[c].residues[r].atoms[a].xyz[2];
            inputtedstrflag[idx][atom] = true;
        }
    }
    j += 1;

    if(stack_prev)
    {
        for (int a=0;a<pdb_entry.chains[c].residues[r-1].atoms.size();a++)
        {
            auto it = find(bb_atoms[base0].begin(), bb_atoms[base0].end(), pdb_entry.chains[c].residues[r-1].atoms[a].name);
            
            if(it != bb_atoms[base0].end()){
                atom = distance(bb_atoms[base0].begin(), it);
                idx=residueidx[j];
                inputtedstr.residues[idx].atoms[atom].coord.x = pdb_entry.chains[c].residues[r-1].atoms[a].xyz[0];
                inputtedstr.residues[idx].atoms[atom].coord.y = pdb_entry.chains[c].residues[r-1].atoms[a].xyz[1];
                inputtedstr.residues[idx].atoms[atom].coord.z = pdb_entry.chains[c].residues[r-1].atoms[a].xyz[2];
                inputtedstrflag[idx][atom] = true;
            }
        }
        j += 1;
    }

    if(stack_next)
    {
        base2=pdb_entry.chains[c].residues[r+1].resn[2];
        for (int a=0;a<pdb_entry.chains[c].residues[r+1].atoms.size();a++)
        {
            auto it = find(bb_atoms[base2].begin(), bb_atoms[base2].end(), pdb_entry.chains[c].residues[r+1].atoms[a].name);
            
            if(it != bb_atoms[base2].end())
            {
                atom = distance(bb_atoms[base2].begin(), it);
                idx=residueidx[j];
                inputtedstr.residues[idx].atoms[atom].coord.x = pdb_entry.chains[c].residues[r+1].atoms[a].xyz[0];
                inputtedstr.residues[idx].atoms[atom].coord.y = pdb_entry.chains[c].residues[r+1].atoms[a].xyz[1];
                inputtedstr.residues[idx].atoms[atom].coord.z = pdb_entry.chains[c].residues[r+1].atoms[a].xyz[2];
                inputtedstrflag[idx][atom] = true;
            }
        }
        j += 1;
    }
}

void prepgrad(const FragmentStructure &stru, vector<vector<Point3s>> &gradptv){

    int i, mtype;
    gradptv.clear();
    gradptv.resize(stru.residues.size());
    for(i=0;i<stru.residues.size();i++){
        gradptv[i].resize(stru.residues[i].atoms.size());
        for(mtype=0;mtype<stru.residues[i].atoms.size();mtype++)
        {
            gradptv[i][mtype].x = 0;
            gradptv[i][mtype].y = 0;
            gradptv[i][mtype].z = 0;
        }
    }
}

void initgrad(const FragmentStructure &stru, vector<vector<Point3s>> &gradptv)
{
	int i, mtype;
	for(i=0;i<stru.residues.size();i++){
		for(mtype=0;mtype<stru.residues[i].atoms.size();mtype++)
		{
			gradptv[i][mtype].x = 0;
			gradptv[i][mtype].y = 0;
			gradptv[i][mtype].z = 0;
		}
	}
}

void setinitobj(FragmentStructure &stru, const FragmentStructure &targetstr, const vector<vector<bool>> &targetflag, vector<vector<Point3s>> &gradient, list<InitPosobj> &binitobj, double K, vector<bool> &moveresidues){
    string base;
    int atom;
    static map<string, vector<int>> targetatoms={
        {"A",{0}},
        {"C",{0}},
        {"G",{0}},
        {"U",{0}}
    };
    for(int i=0;i<stru.residues.size();i++){
        if(!moveresidues[i]) continue;
	base = stru.residues[i].restype;
        for(int j=0;j<targetatoms[base].size();j++){
		atom=targetatoms[base][j];
            if(targetflag[i][atom]){
                //cout<<i << "-" <<j << "-" << targetstr.residues[i].atoms[j].coord.x << "-" << targetstr.residues[i].atoms[j].coord.y << "-" << targetstr.residues[i].atoms[j].coord.z<<endl;
		InitPosobj binitobjone=InitPosobj::create(&stru.residues[i].atoms[atom].coord, &gradient[i][atom], targetstr.residues[i].atoms[atom].coord.x, targetstr.residues[i].atoms[atom].coord.y, targetstr.residues[i].atoms[atom].coord.z, K);
                binitobj.push_back(binitobjone);
            }
        }
    }
}

void setinitpairobj(FragmentStructure &stru, const FragmentStructure &targetstr, const vector<vector<bool>> &targetflag, vector<vector<Point3s>> &gradient, list<InitPairobj> &initpairobj, double K, vector<bool> &moveresidues){
    string base;
    int atom;
    vector<int> targetatoms= {0, 5};
    double x,y,z;
    double initlength;

    for(int i=0;i<stru.residues.size();i++){
        if(!moveresidues[i]) continue;

        for(int j=0;j<stru.residues.size();j++){
            if(!moveresidues[j]) continue;

            for(int k=0;k<targetatoms.size();k++){
                atom=targetatoms[k];
                if(targetflag[i][atom] && targetflag[j][atom]){
                    x = targetstr.residues[i].atoms[atom].coord.x - targetstr.residues[j].atoms[atom].coord.x;
                    y = targetstr.residues[i].atoms[atom].coord.y - targetstr.residues[j].atoms[atom].coord.y;
                    z = targetstr.residues[i].atoms[atom].coord.z - targetstr.residues[j].atoms[atom].coord.z;

                    initlength = sqrt(x*x+y*y+z*z);

                    InitPairobj initpairobjone= InitPairobj::create(&stru.residues[i].atoms[atom].coord, &stru.residues[j].atoms[atom].coord, &gradient[i][atom], &gradient[j][atom], initlength, K);
                    initpairobj.push_back(initpairobjone);
                }
            }
        }
    }
}

void settemplatepairobj(FragmentStructure &stru, const vector<FragmentStructure> &templatestrall, const vector<vector<vector<bool>>> &templatestrallflag, vector<vector<Point3s>> &gradient, list<InitPairobj> &templatepairobj, double K, vector<bool> &moveresidues){
    string base;
    int atom;
    double x,y,z;
    double templatelength;

    for(int n=0;n<templatestrall.size();n++){
        for(int i=0;i<stru.residues.size();i++){
            if(!moveresidues[i]) continue;

            for(int j=0;j<stru.residues.size();j++){
                if(!moveresidues[j]) continue;

                for(int m=0;m<template_target_atoms.size();m++){
                    atom=template_target_atoms_idx[m];
                    
                    if(atom < templatestrallflag[n][i].size() && atom < templatestrallflag[n][j].size()){
                        if(templatestrallflag[n][i][atom] && templatestrallflag[n][j][atom]){
                            x = templatestrall[n].residues[i].atoms[atom].coord.x - templatestrall[n].residues[j].atoms[atom].coord.x;
                            y = templatestrall[n].residues[i].atoms[atom].coord.y - templatestrall[n].residues[j].atoms[atom].coord.y;
                            z = templatestrall[n].residues[i].atoms[atom].coord.z - templatestrall[n].residues[j].atoms[atom].coord.z;

                            templatelength = sqrt(x*x+y*y+z*z);
                            InitPairobj templatepairobjone= InitPairobj::create(&stru.residues[i].atoms[atom].coord, &stru.residues[j].atoms[atom].coord, &gradient[i][atom], &gradient[j][atom], templatelength, K);
                            templatepairobj.push_back(templatepairobjone);
                        }
                    }
                }
            }
        }
    }
}

void setbasepairobj(FragmentStructure &stru, vector<vector<Point3s>> &gradient, vector<vector<int>> &bindmap, const vector<vector<double>> &cbpmap, const vector<vector<double>> &noncbpmap, const double &K, list<BasePairobj> &bpobj)
{
	int i,j,k;
	int res1, res2;
	string aat1, aat2;
    vector<Point3s*> AtomPtrs;
    vector<Point3s*> Atomgrad;

    //cout<<bindmap.size()<<endl;
    for(i=0;i<stru.residues.size();i++){
        //cout<<bindmap[i].size()<<endl;
        res1=i;
        aat1 =  stru.residues[res1].restype;
        for(j=i+1;j<stru.residues.size();j++){
            //cout<<i << "-" << j<<endl;
            res2=j;
            aat2 =  stru.residues[res2].restype;

            if(bindmap[i][j]==2){
                AtomPtrs.clear();
                Atomgrad.clear();
                for(k=12;k<heavily_atom_num[aat1];k++){
                    Point3s* respt = &stru.residues[res1].atoms[k].coord;
                    AtomPtrs.push_back(respt);

                    Point3s* resgradpt = &gradient[res1][k];
                    Atomgrad.push_back(resgradpt);
                }

                for(k=12;k<heavily_atom_num[aat2];k++){
                    Point3s* respt = &stru.residues[res2].atoms[k].coord;
                    AtomPtrs.push_back(respt);

                    Point3s* resgradpt = &gradient[res2][k];
                    Atomgrad.push_back(resgradpt);
                }

                BasePairobj bpobjone(AtomPtrs, Atomgrad, bindmap[i][j]);
                bpobj.push_back(bpobjone);
            }

            if(noncbpmap[i][j] > 0){
                AtomPtrs.clear();
                Atomgrad.clear();

                for(k=12;k<heavily_atom_num[aat1];k++){
                    Point3s* respt = &stru.residues[res1].atoms[k].coord;
                    AtomPtrs.push_back(respt);

                    Point3s* resgradpt = &gradient[res1][k];
                    Atomgrad.push_back(resgradpt);
                }

                for(k=12;k<heavily_atom_num[aat2];k++){
                    Point3s* respt = &stru.residues[res2].atoms[k].coord;
                    AtomPtrs.push_back(respt);

                    Point3s* resgradpt = &gradient[res2][k];
                    Atomgrad.push_back(resgradpt);
                }

                BasePairobj bpobjone(AtomPtrs, Atomgrad, noncbpmap[i][j]);
                bpobj.push_back(bpobjone);
            }
        }

        /*
        AtomPtrs.clear();
        Atomgrad.clear();
        for(k=12;k<heavily_atom_num[aat1];k++){
            Point3s* respt = &stru.residues[res1].atoms[k].coord;
            AtomPtrs.push_back(respt);

            Point3s* resgradpt = &gradient[res1][k];
            Atomgrad.push_back(resgradpt);
        }

        BasePairobj bpobjone1(AtomPtrs, Atomgrad, K);
        bpobj.push_back(bpobjone1);

        AtomPtrs.clear();
        Atomgrad.clear();
        for(k=12;k<heavily_atom_num[aat2];k++){
            Point3s* respt = &stru.residues[res2].atoms[k].coord;
            AtomPtrs.push_back(respt);

            Point3s* resgradpt = &gradient[res2][k];
            Atomgrad.push_back(resgradpt);
        }

        BasePairobj bpobjonent2(AtomPtrs, Atomgrad, K);
        bpobj.push_back(bpobjonent2);
        */
    }

    /*
    for(i=0;i<baselen;i++){
        res1=i;
        res2=baselen*2-i-1;

        vector<Point3s*> AtomPtrs;
        vector<Point3s*> Atomgrad;
        for(j=12;j<stru.residues[res1].atoms.size();j++){
            Point3s* respt = &stru.residues[res1].atoms[j].coord;
            AtomPtrs.push_back(respt);

            Point3s* resgradpt = &gradient[res1][j];
            Atomgrad.push_back(resgradpt);
        }

        for(j=12;j<stru.residues[res2].atoms.size();j++){
            Point3s* respt = &stru.residues[res2].atoms[j].coord;
            AtomPtrs.push_back(respt);

            Point3s* resgradpt = &gradient[res2][j];
            Atomgrad.push_back(resgradpt);
        }

        BasePairobj bpobjone(AtomPtrs, Atomgrad);
        bpobj.push_back(bpobjone);
    }
    */
}


void setbaseplaneobj(FragmentStructure &stru, vector<vector<Point3s>> &gradient, const double &K, list<BasePairobj> &baseplaneobj)
{
	int i,j,k;
	int res;
	string aat;
    vector<Point3s*> AtomPtrs;
    vector<Point3s*> Atomgrad;

    for(i=0;i<stru.residues.size();i++){
        aat =  stru.residues[i].restype;

        AtomPtrs.clear();
        Atomgrad.clear();

        for(k=12;k<heavily_atom_num[aat];k++){
            Point3s* respt = &stru.residues[res].atoms[k].coord;
            AtomPtrs.push_back(respt);

            Point3s* resgradpt = &gradient[res][k];
            Atomgrad.push_back(resgradpt);
        }

        BasePairobj bpobjone(AtomPtrs, Atomgrad, K);
        baseplaneobj.push_back(bpobjone);
    }
}


void stack_superimpose(FragmentStructure &stru, const ModelUnit &ideal_rna_one, int res1, int res2, int order, int ord, double &rmsd)
{
    int i,j;
    int count;
    string base1, base2;
    vector<vector<double>> RotMatix;
    vector<double> TranVect;
    vector<double> xyztmp(3,0);
    int atom;
    BasicFunc bf;
    Point3s vec;
    int base1_temp, base2_temp;

    base1 = stru.residues[res1].restype;
    base2 = stru.residues[res2].restype;

    map<string, vector<int>> targetatoms={
        {"A",{12, 13, 14, 15, 16, 17, 18, 19, 20, 21}},
        {"C",{12, 13, 14, 15, 16, 17, 18, 19}},
        {"G",{12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22}},
        {"U",{12, 13, 14, 15, 16, 17, 18, 19}}
    };

    vector<vector<double> > xyz_list1(targetatoms[base1].size(), vector<double>(3,0));
    vector<vector<double> > xyz_list2(targetatoms[base1].size(), vector<double>(3,0));
    vector<vector<double> > xyz_list3(targetatoms[base2].size(), vector<double>(3,0));

    if(ord == 0){
        if(order == 0){
            base1_temp = 0;
            base2_temp = 1;
        }else if(order == 1){
            base1_temp = 1;
            base2_temp = 0;
        }
    }else if(ord == 1){
        if(order == 0){
            base1_temp = 1;
            base2_temp = 0;
        }else if(order == 1){
            base1_temp = 0;
            base2_temp = 1;
        }
    }

    count = 0;
    for(j=0; j<targetatoms[base1].size();j++){
        atom=targetatoms[base1][j];

        xyz_list1[count][0]=stru.residues[res1].atoms[atom].coord.x;
        xyz_list1[count][1]=stru.residues[res1].atoms[atom].coord.y;
        xyz_list1[count][2]=stru.residues[res1].atoms[atom].coord.z;

        xyz_list2[count][0]=ideal_rna_one.chains[0].residues[base1_temp].atoms[atom].xyz[0];
        xyz_list2[count][1]=ideal_rna_one.chains[0].residues[base1_temp].atoms[atom].xyz[1];
        xyz_list2[count][2]=ideal_rna_one.chains[0].residues[base1_temp].atoms[atom].xyz[2];
        
        count+=1;
    }

    bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

    rmsd=0;
    count=0;
    for(j=0; j<targetatoms[base2].size();j++){
        atom=targetatoms[base2][j];
        xyz_list3[count][0]=ideal_rna_one.chains[0].residues[base2_temp].atoms[atom].xyz[0];
        xyz_list3[count][1]=ideal_rna_one.chains[0].residues[base2_temp].atoms[atom].xyz[1];
        xyz_list3[count][2]=ideal_rna_one.chains[0].residues[base2_temp].atoms[atom].xyz[2];

        bf.Transform(xyz_list3[count],RotMatix,TranVect,xyztmp);

        vec.x=stru.residues[res2].atoms[atom].coord.x-xyztmp[0];
        vec.y=stru.residues[res2].atoms[atom].coord.y-xyztmp[1];
        vec.z=stru.residues[res2].atoms[atom].coord.z-xyztmp[2];

        rmsd+=bf.norm(vec);
        count+=1;
    }
}

void selectstackfrag(FragmentStructure &stru, vector<vector<Point3s>> &gradient, int res1, int res2, const vector<ModelUnit> &ideal_rna, int stidx, vector<vector<int>> &indices, vector<vector<int>> &orders, int ord, const int &sttype, double K, list<StackTemplateobj2> &stackobj)
{
	int i,j;
	double rmsd;
    string base1, base2;

    base1 = stru.residues[res1].restype;
    base2 = stru.residues[res2].restype;

	//static map<string, vector<int>> targetatoms={
		//{"A",{6, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21}},
                //{"C",{6, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}},
                //{"G",{6, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22}},
                //{"U",{6, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}}
	//};

    map<string, vector<int>> targetatoms={
        {"A",{12, 13, 14, 15, 16, 17, 18, 19, 20, 21}},
        {"C",{12, 13, 14, 15, 16, 17, 18, 19}},
        {"G",{12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22}},
        {"U",{12, 13, 14, 15, 16, 17, 18, 19}}
    };
    
    vector<double> rmsds;
    if(sttype == 0 && base1 == base2) rmsds.resize(ideal_rna.size()*2, 0.0);
    else rmsds.resize(ideal_rna.size(), 0.0);

	for(i=0;i<ideal_rna.size();i++)
	{
        rmsd=0;
        stack_superimpose(stru, ideal_rna[i], res1, res2, 0, ord, rmsd);
        rmsds[i]=rmsd;
        stack_superimpose(stru, ideal_rna[i], res2, res1, 1, ord, rmsd);
        rmsds[i]+=rmsd;
	}

	if(sttype == 0 && base1 == base2)
	{
		for(i=0;i<ideal_rna.size();i++)
		{
            rmsd=0;
            stack_superimpose(stru, ideal_rna[i], res1, res2, 1, ord, rmsd);
            rmsds[i]=rmsd;
            stack_superimpose(stru, ideal_rna[i], res2, res1, 0, ord, rmsd);
            rmsds[i]+=rmsd;
		}
	}

    //for(i=0;i<rmsds.size();i++){
        //cout<<i << "-" << rmsds[i]<<endl;
    //}

    indices[stidx].resize(rmsds.size());
    orders[stidx].resize(rmsds.size());
    for (int i = 0; i < rmsds.size(); ++i) indices[stidx][i] = i;

    sort(indices[stidx].begin(), indices[stidx].end(), [&rmsds](int i1, int i2) {
        return rmsds[i1] < rmsds[i2];
    });

    //cout<<"============="<<endl;
    for(i=0;i<indices[stidx].size();i++){
        //cout<<i << "-" << indices[i] << "-" << rmsds[indices[i]]<<endl;

        if(indices[stidx][i] >= ideal_rna.size()){
            indices[stidx][i] = indices[stidx][i] - ideal_rna.size();
            orders[stidx][i] = 1;
        }else{
            orders[stidx][i] = ord;
        }
    }

    int template_idx = indices[stidx][0];
    int template_order = orders[stidx][0];

    StackTemplateobj2 stackobjone(&stru.residues[res1], &stru.residues[res2], &gradient[res1], &gradient[res2], &ideal_rna, &indices[stidx], &orders[stidx], 10, K);
    stackobj.push_back(stackobjone);
}


void setstackingobj(FragmentStructure &stru, vector<vector<Point3s>> &gradient, const map<string, map<string, vector<ModelUnit>>> &ideal_rna, const vector<vector<pair<int, double>>> &stackmap, const map<string, map<string, vector<vector<vector<double>>>>> &stackdistmap, vector<vector<int>> &indices, vector<vector<int>> &orders, double K, list<StackTemplateobj2> &stackobj)
{
	int i,j,k;
	int res1, res2;
	string base1, base2;
    int minfragidx; 
    double minrmsd;
    int baseorder;

    //cout<<bindmap.size()<<endl;

    int stnum=0;
    for(i=0;i<stru.residues.size();i++){
        for(j=i+1;j<stru.residues.size();j++){
            if(stackmap[i][j].first > 0) stnum+=1;
        }
    }

    indices.resize(stnum);
    orders.resize(stnum);

    int count=0;

    for(i=0;i<stru.residues.size();i++){
        //cout<<bindmap[i].size()<<endl;
        res1=i;
        base1 =  stru.residues[res1].restype;
        for(j=i+1;j<stru.residues.size();j++){
            //cout<<i << "-" << j<<endl;
            if(stackmap[i][j].first > 0){
                res2=j;
                base2 = stru.residues[res2].restype;

                if(stackmap[i][j].first == 2){
                    //selectstackfrag(stru, gradient, res1, res2, ideal_rna.at("stacking_adj").at(base1+"_"+base2), stackdistmap.at("stacking_adj").at(base1+"_"+base2), 1, minfragidx, minrmsd, baseorder, K, stackobj);
                    selectstackfrag(stru, gradient, res1, res2, ideal_rna.at("stacking_adj").at(base1+"_"+base2), count, indices, orders, 0, 1, stackmap[i][j].second*K, stackobj);
                
                }else{
                    map<string, string> basestackname = {{"A_A", "A_A"}, {"A_C", "A_C"}, {"A_G", "A_G"}, {"A_U", "A_U"}, {"C_A", "A_C"}, {"C_C", "C_C"}, {"C_G", "C_G"}, {"C_U", "C_U"}, {"G_A", "A_G"}, {"G_C", "C_G"}, {"G_G", "G_G"}, {"G_U", "G_U"}, {"U_A", "A_U"}, {"U_C", "C_U"}, {"U_G", "G_U"}, {"U_U", "U_U"}};
                    
                    if(base1+"_"+base2 == basestackname[base1+"_"+base2]){
                        //selectstackfrag(stru, gradient, res1, res2, ideal_rna.at("stacking_nonadj").at(basestackname[base1+"_"+base2]), stackdistmap.at("stacking_nonadj").at(basestackname[base1+"_"+base2]), 0, minfragidx, minrmsd, baseorder, K, stackobj);
                        selectstackfrag(stru, gradient, res1, res2, ideal_rna.at("stacking_nonadj").at(basestackname[base1+"_"+base2]), count, indices, orders, 0, 0, stackmap[i][j].second*K, stackobj);
                    
                    }else{
                        //selectstackfrag(stru, gradient, res2, res1, ideal_rna.at("stacking_nonadj").at(basestackname[base1+"_"+base2]), stackdistmap.at("stacking_nonadj").at(basestackname[base1+"_"+base2]), 0, minfragidx, minrmsd, baseorder, K, stackobj);
                        selectstackfrag(stru, gradient, res1, res2, ideal_rna.at("stacking_nonadj").at(basestackname[base1+"_"+base2]), count, indices, orders, 1, 0, stackmap[i][j].second*K, stackobj);
                    }
                    
                }

                count+=1;
            }
        }
    }
}


void selectstackfrag(FragmentStructure &stru, vector<vector<Point3s>> &gradient, int res1, int res2, const vector<ModelUnit> &ideal_rna, const vector<vector<vector<double>>> &distmap, const int &sttype, int &minfragidx, double &minrmsd, int &baseorder, double K, list<Stackingobj> &stackobj){
	int i,j;
	int atom, atom1, atom2;
	int count;
	vector<vector<double>> RotMatix;
    vector<double> TranVect;
	minrmsd=100000000000;
	double rmsd;
	vector<double> xyztmp(3,0);
	Point3s vec;
	BasicFunc bf;
    string base1, base2;

    base1 = stru.residues[res1].restype;
    base2 = stru.residues[res2].restype;
    baseorder=0;

	//static map<string, vector<int>> targetatoms={
		//{"A",{6, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21}},
                //{"C",{6, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}},
                //{"G",{6, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22}},
                //{"U",{6, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}}
	//};

        //static map<string, vector<int>> targetatoms={
	        //{"A",{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21}},
                //{"C",{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}},
                //{"G",{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22}},
                //{"U",{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}}
	//};

    //static map<string, vector<int>> targetatoms={
        //{"A",{12, 13, 14, 15, 16, 17, 18, 19, 20, 21}},
        //{"C",{12, 13, 14, 15, 16, 17, 18, 19}},
        //{"G",{12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22}},
        //{"U",{12, 13, 14, 15, 16, 17, 18, 19}}
    //};

    /*
    static map<string, vector<int>> targetatoms={
	{"A",{0,1,2,3,4,5,6,7,8,12,13,14,15,16,17,18,19,20,21}},
        {"C",{0,1,2,3,4,5,6,7,8,12,13,14,15,16,17,18,19}},
        {"G",{0,1,2,3,4,5,6,7,8,12,13,14,15,16,17,18,19,20,21,22}},
        {"U",{0,1,2,3,4,5,6,7,8,12,13,14,15,16,17,18,19}}
	};
    */

    
    static map<string, vector<int>> targetatoms={
	    {"A",{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21}},
        {"C",{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19}},
        {"G",{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}},
        {"U",{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19}}
	};
    
    static map<string, vector<int>> targetatoms_dist={
        {"A",{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21}},
        {"C",{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19}},
        {"G",{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}},
        {"U",{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19}}
    };

    /*
    static map<string, vector<int>> targetatoms={
	    {"A",{12,13,14,15,16,17,18,19,20,21}},
        {"C",{12,13,14,15,16,17,18,19}},
        {"G",{12,13,14,15,16,17,18,19,20,21,22}},
        {"U",{12,13,14,15,16,17,18,19}}
	};

    static map<string, vector<int>> targetatoms_dist={
        {"A",{12, 13, 14, 15, 16, 18, 19, 20, 21}},
        {"C",{12, 13, 15, 16, 18, 19}},
        {"G",{12, 13, 14, 15, 16, 18, 19, 21, 22}},
        {"U",{12, 13, 15, 16, 18, 19}}
    };
    */

	vector<vector<double> > xyz_list1(targetatoms[base1].size()+targetatoms[base2].size(), vector<double>(3,0));
    vector<vector<double> > xyz_list2(targetatoms[base1].size()+targetatoms[base2].size(), vector<double>(3,0));

	for(i=0;i<ideal_rna.size();i++)
	{
		count = 0;
		for(j=0; j<targetatoms[base1].size();j++){
			atom=targetatoms[base1][j];
			xyz_list1[count][0]=stru.residues[res1].atoms[atom].coord.x;
			xyz_list1[count][1]=stru.residues[res1].atoms[atom].coord.y;
			xyz_list1[count][2]=stru.residues[res1].atoms[atom].coord.z;

			xyz_list2[count][0]=ideal_rna[i].chains[0].residues[0].atoms[atom].xyz[0];
			xyz_list2[count][1]=ideal_rna[i].chains[0].residues[0].atoms[atom].xyz[1];
			xyz_list2[count][2]=ideal_rna[i].chains[0].residues[0].atoms[atom].xyz[2];

            count+=1;
		}

		for(j=0; j<targetatoms[base2].size();j++){
			atom=targetatoms[base2][j];
			xyz_list1[count][0]=stru.residues[res2].atoms[atom].coord.x;
			xyz_list1[count][1]=stru.residues[res2].atoms[atom].coord.y;
			xyz_list1[count][2]=stru.residues[res2].atoms[atom].coord.z;

			xyz_list2[count][0]=ideal_rna[i].chains[0].residues[1].atoms[atom].xyz[0];
			xyz_list2[count][1]=ideal_rna[i].chains[0].residues[1].atoms[atom].xyz[1];
			xyz_list2[count][2]=ideal_rna[i].chains[0].residues[1].atoms[atom].xyz[2];

			count+=1;
		}

		bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

		rmsd=0;
		count=0;
		for(j=0; j<targetatoms[base1].size();j++){
			atom=targetatoms[base1][j];
			bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);
			vec.x=stru.residues[res1].atoms[atom].coord.x-xyztmp[0];
			vec.y=stru.residues[res1].atoms[atom].coord.y-xyztmp[1];
			vec.z=stru.residues[res1].atoms[atom].coord.z-xyztmp[2];

			rmsd+=bf.norm(vec);
			count+=1;
		}

		for(j=0; j<targetatoms[base2].size();j++){
			atom=targetatoms[base2][j];
			bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);

			vec.x=stru.residues[res2].atoms[atom].coord.x-xyztmp[0];
			vec.y=stru.residues[res2].atoms[atom].coord.y-xyztmp[1];
			vec.z=stru.residues[res2].atoms[atom].coord.z-xyztmp[2];

			rmsd+=bf.norm(vec);
			count+=1;
		}

		if(minrmsd > rmsd){
            minrmsd=rmsd;
            minfragidx=i;
            
        }
	}
        

	if(sttype == 0 && base1 == base2)
	{
		for(i=0;i<ideal_rna.size();i++)
		{
			int count = 0;
			for(j=0; j<targetatoms[base2].size();j++){
				atom=targetatoms[base2][j];
				xyz_list1[count][0]=stru.residues[res2].atoms[atom].coord.x;
				xyz_list1[count][1]=stru.residues[res2].atoms[atom].coord.y;
				xyz_list1[count][2]=stru.residues[res2].atoms[atom].coord.z;

				xyz_list2[count][0]=ideal_rna[i].chains[0].residues[0].atoms[atom].xyz[0];
				xyz_list2[count][1]=ideal_rna[i].chains[0].residues[0].atoms[atom].xyz[1];
				xyz_list2[count][2]=ideal_rna[i].chains[0].residues[0].atoms[atom].xyz[2];
				count+=1;
			}

			for(j=0; j<targetatoms[base1].size();j++){
				atom=targetatoms[base1][j];
				xyz_list1[count][0]=stru.residues[res1].atoms[atom].coord.x;
				xyz_list1[count][1]=stru.residues[res1].atoms[atom].coord.y;
				xyz_list1[count][2]=stru.residues[res1].atoms[atom].coord.z;

				xyz_list2[count][0]=ideal_rna[i].chains[0].residues[1].atoms[atom].xyz[0];
				xyz_list2[count][1]=ideal_rna[i].chains[0].residues[1].atoms[atom].xyz[1];
				xyz_list2[count][2]=ideal_rna[i].chains[0].residues[1].atoms[atom].xyz[2];
				count+=1;
			}

			bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

			rmsd=0;
			count=0;
			//for(int j=0; j<SNTTEMPLATEATOMS[aat].size();j++){
			for(j=0; j<targetatoms[base2].size();j++){
				atom=targetatoms[base2][j];
				bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);
				vec.x=stru.residues[res2].atoms[atom].coord.x-xyztmp[0];
				vec.y=stru.residues[res2].atoms[atom].coord.y-xyztmp[1];
				vec.z=stru.residues[res2].atoms[atom].coord.z-xyztmp[2];
				rmsd+=bf.norm(vec);
				count+=1;
			}

			for(j=0; j<targetatoms[base1].size();j++){
				atom=targetatoms[base1][j];
				bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);
				vec.x=stru.residues[res1].atoms[atom].coord.x-xyztmp[0];
				vec.y=stru.residues[res1].atoms[atom].coord.y-xyztmp[1];
				vec.z=stru.residues[res1].atoms[atom].coord.z-xyztmp[2];
				rmsd+=bf.norm(vec);
				count+=1;
			}

			if(minrmsd > rmsd){
				minrmsd=rmsd;
                minfragidx=i;
                baseorder=1;
			}
		}
	}

    //cout<<res1 << "-" << res2 << "-" << minfragidx << "-" << minrmsd << "-" << sttype<<endl;
    //cout<<minrmsd << "-" << minfragidx << "-" << baseorder<<endl;

    double length, dx, dy, dz;
    if(baseorder == 0){
        for(int i=0;i<targetatoms_dist[base1].size();i++){
            atom1=targetatoms_dist[base1][i];
            for(int j=0;j<targetatoms_dist[base2].size();j++){
                atom2=targetatoms_dist[base2][j];
                
                /*
                dx = stru.residues[res1].atoms[atom1].coord.x - stru.residues[res2].atoms[atom2].coord.x;
                dy = stru.residues[res1].atoms[atom1].coord.y - stru.residues[res2].atoms[atom2].coord.y;
                dz = stru.residues[res1].atoms[atom1].coord.z - stru.residues[res2].atoms[atom2].coord.z;
                length = sqrt(dx*dx + dy*dy+ dz*dz);
                cout<< atom1 << "-" << atom2 << "-" << length << "-" << distmap[minfragidx][atom1][atom2]<<endl;
                */

                Stackingobj stackobjone(&stru.residues[res1].atoms[atom1].coord, &stru.residues[res2].atoms[atom2].coord, &gradient[res1][atom1], &gradient[res2][atom2], distmap[minfragidx][atom1][atom2], K);
                stackobj.push_back(stackobjone);
            }
        }
    }else{
        for(int i=0;i<targetatoms_dist[base1].size();i++){
            atom1=targetatoms_dist[base1][i];
            for(int j=0;j<targetatoms_dist[base2].size();j++){
                atom2=targetatoms_dist[base2][j];

                /*
                dx = stru.residues[res1].atoms[atom1].coord.x - stru.residues[res2].atoms[atom2].coord.x;
                dy = stru.residues[res1].atoms[atom1].coord.y - stru.residues[res2].atoms[atom2].coord.y;
                dz = stru.residues[res1].atoms[atom1].coord.z - stru.residues[res2].atoms[atom2].coord.z;
                length = sqrt(dx*dx + dy*dy+ dz*dz);
                cout<< atom1 << "-" << atom2 << "-" << length << "-" << distmap[minfragidx][atom2][atom1]<<endl;
                */

                Stackingobj stackobjone(&stru.residues[res1].atoms[atom1].coord, &stru.residues[res2].atoms[atom2].coord, &gradient[res1][atom1], &gradient[res2][atom2], distmap[minfragidx][atom2][atom1], K);
                stackobj.push_back(stackobjone);
            }
        }
    }
}

void setstackingobj(FragmentStructure &stru, vector<vector<Point3s>> &gradient, const map<string, map<string, vector<ModelUnit>>> &ideal_rna, const vector<vector<pair<int, double>>> &stackmap, const map<string, map<string, vector<vector<vector<double>>>>> &stackdistmap, double K, list<Stackingobj> &stackobj)
{
	int i,j,k;
	int res1, res2;
	string base1, base2;
    int minfragidx; 
    double minrmsd;
    int baseorder;

    //cout<<bindmap.size()<<endl;

    for(i=0;i<stru.residues.size();i++){
        //cout<<bindmap[i].size()<<endl;
        res1=i;
        base1 =  stru.residues[res1].restype;
        for(j=i+1;j<stru.residues.size();j++){
            //cout<<i << "-" << j<<endl;
            if(stackmap[i][j].first > 0){
                res2=j;
                base2 = stru.residues[res2].restype;
                //cout<<res1 << "-" << res2 << "-" << stackmap[i][j].first<<endl;

                if(stackmap[i][j].first == 2){
                    selectstackfrag(stru, gradient, res1, res2, ideal_rna.at("stacking_adj").at(base1+"_"+base2), stackdistmap.at("stacking_adj").at(base1+"_"+base2), 1, minfragidx, minrmsd, baseorder, K, stackobj);
                }else{
                    map<string, string> basestackname = {{"A_A", "A_A"}, {"A_C", "A_C"}, {"A_G", "A_G"}, {"A_U", "A_U"}, {"C_A", "A_C"}, {"C_C", "C_C"}, {"C_G", "C_G"}, {"C_U", "C_U"}, {"G_A", "A_G"}, {"G_C", "C_G"}, {"G_G", "G_G"}, {"G_U", "G_U"}, {"U_A", "A_U"}, {"U_C", "C_U"}, {"U_G", "G_U"}, {"U_U", "U_U"}};
                    
                    if(base1+"_"+base2 == basestackname[base1+"_"+base2]){
                        selectstackfrag(stru, gradient, res1, res2, ideal_rna.at("stacking_nonadj").at(basestackname[base1+"_"+base2]), stackdistmap.at("stacking_nonadj").at(basestackname[base1+"_"+base2]), 0, minfragidx, minrmsd, baseorder, stackmap[i][j].second*K, stackobj);
                    }else{
                        selectstackfrag(stru, gradient, res2, res1, ideal_rna.at("stacking_nonadj").at(basestackname[base1+"_"+base2]), stackdistmap.at("stacking_nonadj").at(basestackname[base1+"_"+base2]), 0, minfragidx, minrmsd, baseorder, stackmap[i][j].second*K, stackobj);
                    }
                    
                }
            }
        }
    }
}

void setstackinggeoobj(FragmentStructure &stru, vector<vector<Point3s>> &gradient, const map<int, map<string, map<string, vector<vector<double>>>>> &stackparam, const vector<vector<pair<int, double>>> &stackmap, vector<vector<Point3s*>> &AtomPtrslistA, vector<vector<Point3s*>> &AtomgradlistA, vector<vector<Point3s*>> &AtomPtrslistB, vector<vector<Point3s*>> &AtomgradlistB, double K, list<StackGeoobjMC> &stackgeoobj)
{
	int i,j,k;
	int res1, res2;
	string base1, base2;
    int minfragidx; 
    double minrmsd;
    int baseorder;

    //cout<<bindmap.size()<<endl;
    int stacknum=0;
    for(i=0;i<stru.residues.size();i++){
        for(j=i+1;j<stru.residues.size();j++){
            if(stackmap[i][j].first > 0) stacknum+=1;
        }
    }
    AtomPtrslistA.resize(stacknum);
    AtomgradlistA.resize(stacknum);
    AtomPtrslistB.resize(stacknum);
    AtomgradlistB.resize(stacknum);

    int count=0;
    for(i=0;i<stru.residues.size();i++){
        //cout<<bindmap[i].size()<<endl;
        res1=i;
        base1 =  stru.residues[res1].restype;
        for(j=i+1;j<stru.residues.size();j++){
            //cout<<i << "-" << j<<endl;
            if(stackmap[i][j].first > 0){
                res2=j;
                base2 = stru.residues[res2].restype;

                for(k=12;k<heavily_atom_num[base1];k++){
                    Point3s* respt = &stru.residues[res1].atoms[k].coord;
                    AtomPtrslistA[count].push_back(respt);

                    Point3s* resgradpt = &gradient[res1][k];
                    AtomgradlistA[count].push_back(resgradpt);
                }

                for(k=12;k<heavily_atom_num[base2];k++){
                    Point3s* respt = &stru.residues[res2].atoms[k].coord;
                    AtomPtrslistB[count].push_back(respt);

                    Point3s* resgradpt = &gradient[res2][k];
                    AtomgradlistB[count].push_back(resgradpt);
                }

                if(stackmap[i][j].first == 2){
                    StackGeoobjMC stackgeoobjone(AtomPtrslistA[count], AtomPtrslistB[count], AtomgradlistA[count], AtomgradlistB[count], stru.residues[res1].restype[0], stru.residues[res2].restype[0], stackparam.at(1).at(base1).at(base2)[0][0]+3*stackparam.at(1).at(base1).at(base2)[0][1], stackparam.at(1).at(base1).at(base2)[1][0]+3*stackparam.at(1).at(base1).at(base2)[1][1], stackmap[i][j].second*K, stackmap[i][j].second*K);
                    stackgeoobj.push_back(stackgeoobjone);
                }else{
                    StackGeoobjMC stackgeoobjone(AtomPtrslistA[count], AtomPtrslistB[count], AtomgradlistA[count], AtomgradlistB[count], stru.residues[res1].restype[0], stru.residues[res2].restype[0], stackparam.at(0).at(base1).at(base2)[0][0]+3*stackparam.at(0).at(base1).at(base2)[0][1], stackparam.at(0).at(base1).at(base2)[1][0]+3*stackparam.at(0).at(base1).at(base2)[1][1], stackmap[i][j].second*K, stackmap[i][j].second*K);
                    stackgeoobj.push_back(stackgeoobjone);
                }
                count+=1;
            }
        }
    }
}

void selectbpfrag(FragmentStructure &stru, vector<vector<Point3s>> &gradient, int res1, int res2, const vector<ModelUnit> &ideal_rna, const vector<vector<vector<double>>> &distmap, int &minfragidx, double &minrmsd, const double &K, list<BPdistobj> &bpdistobj){
	int i,j;
	int atom, atom1, atom2;
	int count;
	vector<vector<double>> RotMatix;
    vector<double> TranVect;
	minrmsd=100000000000;
	double rmsd;
	vector<double> xyztmp(3,0);
	Point3s vec;
	BasicFunc bf;
    string base1, base2;
    int base_idx_1, base_idx_2;

    base1 = stru.residues[res1].restype;
    base2 = stru.residues[res2].restype;

    static map<string, vector<int>> targetatoms={
	    {"A",{0,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19,20,21}},
        {"C",{0,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19}},
        {"G",{0,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19,20,21,22}},
        {"U",{0,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19}}
	};

    static map<string, vector<int>> targetatoms_dist={
        {"A",{12, 13, 14, 15, 16, 18, 19, 20, 21}},
        {"C",{12, 13, 15, 16, 18, 19}},
        {"G",{12, 13, 14, 15, 16, 18, 19, 21, 22}},
        {"U",{12, 13, 15, 16, 18, 19}}
    };

	vector<vector<double> > xyz_list1(targetatoms[base1].size()+targetatoms[base2].size(), vector<double>(3,0));
    vector<vector<double> > xyz_list2(targetatoms[base1].size()+targetatoms[base2].size(), vector<double>(3,0));

    /*
    if(order == 0){
        base_idx_1=0;
        base_idx_2=1;
    }else{
        base_idx_1=1;
        base_idx_2=0;
    }
    */
    base_idx_1=0;
    base_idx_2=1;
    
	for(i=0;i<ideal_rna.size();i++)
	{
		count = 0;
		for(j=0; j<targetatoms[base1].size();j++){
			atom=targetatoms[base1][j];
			xyz_list1[count][0]=stru.residues[res1].atoms[atom].coord.x;
			xyz_list1[count][1]=stru.residues[res1].atoms[atom].coord.y;
			xyz_list1[count][2]=stru.residues[res1].atoms[atom].coord.z;

			xyz_list2[count][0]=ideal_rna[i].chains[0].residues[base_idx_1].atoms[atom].xyz[0];
			xyz_list2[count][1]=ideal_rna[i].chains[0].residues[base_idx_1].atoms[atom].xyz[1];
			xyz_list2[count][2]=ideal_rna[i].chains[0].residues[base_idx_1].atoms[atom].xyz[2];

            count+=1;
		}

		for(j=0; j<targetatoms[base2].size();j++){
			atom=targetatoms[base2][j];
			xyz_list1[count][0]=stru.residues[res2].atoms[atom].coord.x;
			xyz_list1[count][1]=stru.residues[res2].atoms[atom].coord.y;
			xyz_list1[count][2]=stru.residues[res2].atoms[atom].coord.z;

			xyz_list2[count][0]=ideal_rna[i].chains[0].residues[base_idx_2].atoms[atom].xyz[0];
			xyz_list2[count][1]=ideal_rna[i].chains[0].residues[base_idx_2].atoms[atom].xyz[1];
			xyz_list2[count][2]=ideal_rna[i].chains[0].residues[base_idx_2].atoms[atom].xyz[2];

			count+=1;
		}

		bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

		rmsd=0;
		count=0;
		for(j=0; j<targetatoms[base1].size();j++){
			atom=targetatoms[base1][j];
			bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);
			vec.x=stru.residues[res1].atoms[atom].coord.x-xyztmp[0];
			vec.y=stru.residues[res1].atoms[atom].coord.y-xyztmp[1];
			vec.z=stru.residues[res1].atoms[atom].coord.z-xyztmp[2];

			rmsd+=bf.norm(vec);
			count+=1;
		}

		for(j=0; j<targetatoms[base2].size();j++){
			atom=targetatoms[base2][j];
			bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);

			vec.x=stru.residues[res2].atoms[atom].coord.x-xyztmp[0];
			vec.y=stru.residues[res2].atoms[atom].coord.y-xyztmp[1];
			vec.z=stru.residues[res2].atoms[atom].coord.z-xyztmp[2];

			rmsd+=bf.norm(vec);
			count+=1;
		}

		if(minrmsd > rmsd){
            minrmsd=rmsd;
            minfragidx=i;
            
        }
	}

    //cout<<res1 << "-" << res2 << "-" << minfragidx << "-" << minrmsd << "-" << sttype<<endl;
    //cout<<minrmsd << "-" << minfragidx << "-" << baseorder<<endl;

    double length, dx, dy, dz;
    for(int i=0;i<targetatoms_dist[base1].size();i++){
        atom1=targetatoms_dist[base1][i];
        for(int j=0;j<targetatoms_dist[base2].size();j++){
            atom2=targetatoms_dist[base2][j];
            
            dx = stru.residues[res1].atoms[atom1].coord.x - stru.residues[res2].atoms[atom2].coord.x;
            dy = stru.residues[res1].atoms[atom1].coord.y - stru.residues[res2].atoms[atom2].coord.y;
            dz = stru.residues[res1].atoms[atom1].coord.z - stru.residues[res2].atoms[atom2].coord.z;
            length = sqrt(dx*dx + dy*dy+ dz*dz);
            //cout<< atom1 << "-" << atom2 << "-" << length << "-" << distmap[minfragidx][atom1][atom2]<<endl;

            BPdistobj bpdistobjone(&stru.residues[res1].atoms[atom1].coord, &stru.residues[res2].atoms[atom2].coord, &gradient[res1][atom1], &gradient[res2][atom2], distmap[minfragidx][atom1][atom2], K);
            bpdistobj.push_back(bpdistobjone);
        }
    }
}


void setbpobj(FragmentStructure &stru, vector<vector<Point3s>> &gradient, const map<string, map<string, vector<ModelUnit>>> &ideal_rna, vector<vector<int>> &bpmap, const map<string, vector<vector<vector<double>>>> &bpdistmap, double K, list<BPdistobj> &bpdistobj)
{
	int i,j,k;
	int res1, res2;
	string base1, base2;
    int minfragidx; 
    double minrmsd;
    int baseorder;

    map<string, string> basebpname = {{"A-U", "A-U"}, {"U-A", "A-U"}, {"G-C", "G-C"}, {"C-G", "G-C"}, {"G-U", "G-U"}, {"U-G", "G-U"}};
                    
    //cout<<bindmap.size()<<endl;

    for(i=0;i<stru.residues.size();i++){
        res1=i;
        base1 =  stru.residues[res1].restype;
        for(j=i+1;j<stru.residues.size();j++){
            if(bpmap[i][j]==2){
                res2=j;
                base2 = stru.residues[res2].restype;
                if(base1+"-"+base2 == basebpname[base1+"-"+base2]){
                    selectbpfrag(stru, gradient, res1, res2, ideal_rna.at("bp").at(basebpname[base1+"-"+base2]), bpdistmap.at(basebpname[base1+"-"+base2]), minfragidx, minrmsd, K, bpdistobj);
                }else{
                    selectbpfrag(stru, gradient, res2, res1, ideal_rna.at("bp").at(basebpname[base1+"-"+base2]), bpdistmap.at(basebpname[base1+"-"+base2]), minfragidx, minrmsd, K, bpdistobj);
                }
            }
        }
    }
}

void nt_superimpose(FragmentStructure &stru, const ModelUnit &ideal_rna_one, int res, double &rmsd)
{
    int i,j;
    int count;
    string base;
    vector<vector<double>> RotMatix;
    vector<double> TranVect;
    vector<double> xyztmp(3,0);
    int atom;
    BasicFunc bf;
    Point3s vec;

    base = stru.residues[res].restype;

    map<string, vector<int>> targetatoms={
        {"A",{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21}},
        {"C",{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}},
        {"G",{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22}},
        {"U",{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}}
    };

    vector<vector<double> > xyz_list1(targetatoms[base].size(), vector<double>(3,0));
    vector<vector<double> > xyz_list2(targetatoms[base].size(), vector<double>(3,0));

    count = 0;
    for(j=0; j<targetatoms[base].size();j++){
        atom=targetatoms[base][j];

        xyz_list1[count][0]=stru.residues[res].atoms[atom].coord.x;
        xyz_list1[count][1]=stru.residues[res].atoms[atom].coord.y;
        xyz_list1[count][2]=stru.residues[res].atoms[atom].coord.z;

        xyz_list2[count][0]=ideal_rna_one.chains[0].residues[0].atoms[atom].xyz[0];
        xyz_list2[count][1]=ideal_rna_one.chains[0].residues[0].atoms[atom].xyz[1];
        xyz_list2[count][2]=ideal_rna_one.chains[0].residues[0].atoms[atom].xyz[2];
        
        count+=1;
    }

    bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

    rmsd=0;
    count=0;
    for(j=0; j<targetatoms[base].size();j++){
        atom=targetatoms[base][j];

        bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);

        vec.x=xyz_list1[count][0]-xyztmp[0];
        vec.y=xyz_list1[count][1]-xyztmp[1];
        vec.z=xyz_list1[count][2]-xyztmp[2];

        rmsd+=bf.norm(vec);
        count+=1;
    }
}


void selectntfrag(FragmentStructure &stru, vector<vector<Point3s>> &gradient, int res, const vector<ModelUnit> &ideal_rna, int stidx, vector<vector<int>> &indices, double K, list<SingleNtTemplateobj> &ntobj)
{
	int i,j;
	double rmsd;

    string base = stru.residues[res].restype;

    map<string, vector<int>> targetatoms={
        {"A",{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21}},
        {"C",{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}},
        {"G",{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22}},
        {"U",{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}}
    };
    
    vector<double> rmsds(ideal_rna.size(), 0.0);

	for(i=0;i<ideal_rna.size();i++)
	{
        rmsd=0;
        nt_superimpose(stru, ideal_rna[i], res, rmsd);
        rmsds[i]=rmsd;
	}

    indices[stidx].resize(rmsds.size());
    for (int i = 0; i < rmsds.size(); ++i) indices[stidx][i] = i;

    sort(indices[stidx].begin(), indices[stidx].end(), [&rmsds](int i1, int i2) {
        return rmsds[i1] < rmsds[i2];
    });

    int template_idx = indices[stidx][0];

    SingleNtTemplateobj nttempobjone(&stru.residues[res], &gradient[res], &ideal_rna, &indices[stidx], 10, K);
    ntobj.push_back(nttempobjone);
}


void setnttemplateobj(FragmentStructure &stru, vector<vector<Point3s>> &gradient, const map<string, map<string, vector<ModelUnit>>> &ideal_rna, vector<vector<int>> &indices, double K, list<SingleNtTemplateobj> &nttempobj)
{
	int i,j,k;
	int res;
	string base;
    int minfragidx; 
    double minrmsd;
    int baseorder;

    indices.resize(stru.residues.size());
    for(i=0;i<stru.residues.size();i++){
        res=i;
        base =  stru.residues[res].restype;
        selectntfrag(stru, gradient, res, ideal_rna.at("nt").at(base), i, indices, K, nttempobj);
    }
}


void setlengthobj(FragmentStructure &stru, vector<vector<Point3s>> &gradient, int trandposst, int trandposed, vector<vector<int>> &bindmap, const map<string, vector<vector<double>>> &blengthparam, const map<string, vector<vector<double>>> &blengthstdval, const map<string, vector<vector<double>>> &adjblengthparam, const map<string, vector<vector<double>>> &adjblengthstdval, const map<string, vector<vector<double>>> &hydparam, const map<string, vector<vector<int>>> &hydparamidx, const map<string, vector<vector<double>>> &addhydparam, const map<string, vector<vector<int>>> &addhydparamidx, double K, list<BondLengthobj> &blengthobj, vector<bool> &moveresidues)
{
    int i,j,k;
    int res1, res2;
    int atom1, atom2;
    double ideal,stdval;
    string base;

    for(i=trandposst;i<trandposed;i++){
        base = stru.residues[i].restype;
        if(i > 0 && bindmap[i-1][i] == 1){
            for(j=0;j<adjblengthparam.at(base).size();j++){
                for(k=0;k<adjblengthparam.at(base)[j].size();k++){
                    if(adjblengthparam.at(base)[j][k] < 0) continue;

                    res1=i-1; res2=i;
                    atom1=j; atom2=k;
                    ideal=adjblengthparam.at(base)[j][k];
                    stdval=adjblengthstdval.at(base)[j][k];

                    //cout<<ideal << "-" << stdval<<endl;

                    if(atom1 < stru.residues[res1].atoms.size() && atom2 < stru.residues[res2].atoms.size()){
                        BondLengthobj blengtheobjone = BondLengthobj::create(&stru.residues[res1].atoms[atom1].coord, &stru.residues[res2].atoms[atom2].coord, &gradient[res1][atom1], &gradient[res2][atom2], ideal, stdval, K);
                        blengthobj.push_back(blengtheobjone);
                    }
                }
            }
        }

        if(!moveresidues[i]) continue;
        for(j=0;j<blengthparam.at(base).size();j++){
            for(k=j+1;k<blengthparam.at(base)[j].size();k++){
                if(blengthparam.at(base)[j][k] < 0) continue;

                res1=i; res2=i;
                atom1=j; atom2=k;
                ideal=blengthparam.at(base)[j][k];
                stdval=blengthstdval.at(base)[j][k];

                //cout<<ideal << "-" << stdval<<endl;

                if(atom1 < stru.residues[res1].atoms.size() && atom2 < stru.residues[res2].atoms.size())
                {    
                    BondLengthobj blengtheobjone = BondLengthobj::create(&stru.residues[res1].atoms[atom1].coord, &stru.residues[res2].atoms[atom2].coord, &gradient[res1][atom1], &gradient[res2][atom2], ideal, stdval, K);
                    blengthobj.push_back(blengtheobjone);
                }
            }
        }

        for(j=0;j<hydparamidx.at(base).size();j++){
            res1=i; res2=i;
            atom1=hydparamidx.at(base)[j][2]; atom2=hydparamidx.at(base)[j][3];
            ideal=hydparam.at(base)[j][4];
            stdval=hydparam.at(base)[j][5];

            //cout<<ideal << "-" << stdval<<endl;

            if(atom1 < stru.residues[res1].atoms.size() && atom2 < stru.residues[res2].atoms.size())
            {
                BondLengthobj blengtheobjone = BondLengthobj::create(&stru.residues[res1].atoms[atom1].coord, &stru.residues[res2].atoms[atom2].coord, &gradient[res1][atom1], &gradient[res2][atom2], ideal, stdval, K);
                blengthobj.push_back(blengtheobjone);
            }
        }

        for(j=0;j<addhydparamidx.at(base).size();j++){
            res1=i; res2=i;
            atom1=addhydparamidx.at(base)[j][0]; atom2=addhydparamidx.at(base)[j][1];
            ideal=addhydparam.at(base)[j][0];
            stdval=addhydparam.at(base)[j][1];

            //cout<<ideal << "-" << stdval<<endl;

            if(atom1 < stru.residues[res1].atoms.size() && atom2 < stru.residues[res2].atoms.size())
            {
                BondLengthobj blengtheobjone = BondLengthobj::create(&stru.residues[res1].atoms[atom1].coord, &stru.residues[res2].atoms[atom2].coord, &gradient[res1][atom1], &gradient[res2][atom2], ideal, stdval, K);
                blengthobj.push_back(blengtheobjone);
            }
        }
    }
}

void setangleobj(FragmentStructure &stru, vector<vector<Point3s>> &gradient, int trandposst, int trandposed, vector<vector<int>> &bindmap, const map<string, vector<vector<double>>> &ntangleparam, const map<string, vector<vector<int>>> &ntangleparamidx, const map<string, vector<vector<double>>> &adjangleparam, const map<string, vector<vector<int>>> &adjangleparamidx, const map<string, vector<vector<double>>> &hydangleparam, const map<string, vector<vector<int>>> &hydangleparamidx, double K, list<BondAngleobj> &bangleobj, vector<bool> &moveresidues)
{
    int i,j,k;
    int res1, res2, res3;
    int atom1, atom2, atom3;
    double ideal, stdval;
    string base;
    for(i=trandposst;i<trandposed;i++){
        base = stru.residues[i].restype;

        if(i > 0 && bindmap[i-1][i] == 1){
            for(j=0;j<adjangleparamidx.at(base).size();j++){
                res1=i+adjangleparamidx.at(base)[j][0]; 
                res2=i+adjangleparamidx.at(base)[j][1]; 
                res3=i+adjangleparamidx.at(base)[j][2];

                atom1=adjangleparamidx.at(base)[j][3]; 
                atom2=adjangleparamidx.at(base)[j][4]; 
                atom3=adjangleparamidx.at(base)[j][5];

                ideal=adjangleparam.at(base)[j][0];
                stdval=adjangleparam.at(base)[j][1];

                //cout<<res1 << "-" << res2 << "-" << res3<<endl;
                //cout<<atom1 << "-" << atom2 << "-" << atom3<<endl;
                //cout<<ideal << "-" << stdval<<endl;

                if(atom1 < stru.residues[res1].atoms.size() && atom2 < stru.residues[res2].atoms.size() && atom3 < stru.residues[res3].atoms.size()){
                    BondAngleobj bangleobjone = BondAngleobj::create(&stru.residues[res1].atoms[atom1].coord, &stru.residues[res2].atoms[atom2].coord, &stru.residues[res3].atoms[atom3].coord, &gradient[res1][atom1], &gradient[res2][atom2], &gradient[res3][atom3], ideal, stdval, K);
                    bangleobj.push_back(bangleobjone);
                }
            }
        }

        if(!moveresidues[i]) continue;
        for(j=0;j<ntangleparamidx.at(base).size();j++){
            res1=i; res2=i; res3=i;
            atom1=ntangleparamidx.at(base)[j][0];
            atom2=ntangleparamidx.at(base)[j][1];
            atom3=ntangleparamidx.at(base)[j][2];

            ideal=ntangleparam.at(base)[j][0];
            stdval=ntangleparam.at(base)[j][1];

            //cout<<res1 << "-" << res2 << "-" << res3<<endl;
            //cout<<atom1 << "-" << atom2 << "-" << atom3<<endl;
            //cout<<ideal << "-" << stdval<<endl;

            if(atom1 < stru.residues[res1].atoms.size() && atom2 < stru.residues[res2].atoms.size() && atom3 < stru.residues[res3].atoms.size()){
                BondAngleobj bangleobjone = BondAngleobj::create(&stru.residues[res1].atoms[atom1].coord, &stru.residues[res2].atoms[atom2].coord, &stru.residues[res3].atoms[atom3].coord, &gradient[res1][atom1], &gradient[res2][atom2], &gradient[res3][atom3], ideal, stdval, K);
                bangleobj.push_back(bangleobjone);
            }
        }
        

        for(j=0;j<hydangleparamidx.at(base).size();j++){
            res1=i; res2=i; res3=i;
            atom1=hydangleparamidx.at(base)[j][0];
            atom2=hydangleparamidx.at(base)[j][1];
            atom3=hydangleparamidx.at(base)[j][2];

            ideal=hydangleparam.at(base)[j][0];
            stdval=hydangleparam.at(base)[j][1];

            //cout<<res1 << "-" << res2 << "-" << res3<<endl;
            //cout<<atom1 << "-" << atom2 << "-" << atom3<<endl;
            //cout<<ideal << "-" << stdval<<endl;

            if(atom1 < stru.residues[res1].atoms.size() && atom2 < stru.residues[res2].atoms.size() && atom3 < stru.residues[res3].atoms.size()){
                BondAngleobj bangleobjone = BondAngleobj::create(&stru.residues[res1].atoms[atom1].coord, &stru.residues[res2].atoms[atom2].coord, &stru.residues[res3].atoms[atom3].coord, &gradient[res1][atom1], &gradient[res2][atom2], &gradient[res3][atom3], ideal, stdval, K);
                bangleobj.push_back(bangleobjone);
            }
        }
    }
    //exit(1);
}

void setdihedralobj(FragmentStructure &stru, vector<vector<Point3s>> &gradient, int trandposst, int trandposed, const map<string, vector<pair<vector<int>, vector<double>>>> &hydtorsionparam, map<string, vector<pair<vector<int>, vector<double>>>> torsionparam, double K, list<BondDihedralobj> &dihedralobj, vector<bool> &moveresidues)
{
    int i,j,k;
    int res1, res2, res3, res4;
    int atom1, atom2, atom3, atom4;
    double idv, dki, dphase, dpn;
    string base;

    for(i=trandposst;i<trandposed;i++){
        base = stru.residues[i].restype;

        const vector<pair<vector<int>, vector<double>>> &torsionparam_temp = torsionparam.at(base);
        for(j=0;j<torsionparam_temp.size();j++){
            atom1 = torsionparam_temp[j].first[0];
            atom2 = torsionparam_temp[j].first[1];
            atom3 = torsionparam_temp[j].first[2];
            atom4 = torsionparam_temp[j].first[3];

            res1=i; res2=i; res3=i; res4=i;

            idv = torsionparam_temp[j].second[0];
            dki = torsionparam_temp[j].second[1];
            dphase = torsionparam_temp[j].second[2];
            dpn = torsionparam_temp[j].second[3];

            //cout<<atom1 << "-" << atom2 << "-" << atom3 << "-" << atom4<<endl;
            //cout<<idv << "-" << dki << "-" << dphase << "-" << dpn<<endl;
            
            BondDihedralobj dihedralobjone = BondDihedralobj::create(&stru.residues[res1].atoms[atom1].coord, &stru.residues[res2].atoms[atom2].coord, &stru.residues[res3].atoms[atom3].coord, &stru.residues[res4].atoms[atom4].coord, &gradient[res1][atom1], &gradient[res2][atom2], &gradient[res3][atom3], &gradient[res4][atom4], K, dpn, dphase, idv);
            dihedralobj.push_back(dihedralobjone);
        }


        const vector<pair<vector<int>, vector<double>>> &hydtorsionparam_temp = hydtorsionparam.at(base);

        for(j=0;j<hydtorsionparam_temp.size();j++){
            atom1 = hydtorsionparam_temp[j].first[0];
            atom2 = hydtorsionparam_temp[j].first[1];
            atom3 = hydtorsionparam_temp[j].first[2];
            atom4 = hydtorsionparam_temp[j].first[3];

            res1=i; res2=i; res3=i; res4=i;

            idv = hydtorsionparam_temp[j].second[0];
            dki = hydtorsionparam_temp[j].second[1];
            dphase = hydtorsionparam_temp[j].second[2];
            dpn = hydtorsionparam_temp[j].second[3];

            //cout<<atom1 << "-" << atom2 << "-" << atom3 << "-" << atom4<<endl;
            //cout<<idv << "-" << dki << "-" << dphase << "-" << dpn<<endl;
            
            BondDihedralobj dihedralobjone = BondDihedralobj::create(&stru.residues[res1].atoms[atom1].coord, &stru.residues[res2].atoms[atom2].coord, &stru.residues[res3].atoms[atom3].coord, &stru.residues[res4].atoms[atom4].coord, &gradient[res1][atom1], &gradient[res2][atom2], &gradient[res3][atom3], &gradient[res4][atom4], K, dpn, dphase, idv);
            dihedralobj.push_back(dihedralobjone);
        }

    }
}

void setimpropersobj(FragmentStructure &stru, vector<vector<Point3s>> &gradient, int trandposst, int trandposed, list<BondAngleobj> &bangleobj, list<BondImpropersobj> &improperobj)
{
}

void getgradientinit(list<InitPosobj> binitobj){
	for (std::list<InitPosobj>::iterator it = binitobj.begin(); it != binitobj.end(); ++it) {
		it->calcgradientinit();
	}
}

double getenergyinit(list<InitPosobj> binitobj){
	double energy=0.0;
	for (std::list<InitPosobj>::const_iterator it = binitobj.begin(); it != binitobj.end(); ++it) {
		energy+=it->calcenergyinit();
	}
	
	return energy;
}

void getgradientinitpair(list<InitPairobj> binitpairobj){
	for (std::list<InitPairobj>::iterator it = binitpairobj.begin(); it != binitpairobj.end(); ++it) {
		it->calcgradientinitpair();
	}
}

double getenergyinitpair(list<InitPairobj> binitpairobj){
	double energy=0.0;
	for (std::list<InitPairobj>::const_iterator it = binitpairobj.begin(); it != binitpairobj.end(); ++it) {
		energy+=it->calcenergyinitpair();
	}
	
	return energy;
}

void getgradienttemplatepair(list<InitPairobj> btemplatepairobj){
	for (std::list<InitPairobj>::iterator it = btemplatepairobj.begin(); it != btemplatepairobj.end(); ++it) {
		it->calcgradientinitpair();
	}
}

double getenergytemplatepair(list<InitPairobj> btemplatepairobj){
	double energy=0.0;
	for (std::list<InitPairobj>::const_iterator it = btemplatepairobj.begin(); it != btemplatepairobj.end(); ++it) {
		energy+=it->calcenergyinitpair();
	}
	
	return energy;
}

void getgradientbasepair(list<BasePairobj> &baseplaneobj){
	for (std::list<BasePairobj>::iterator it = baseplaneobj.begin(); it != baseplaneobj.end(); ++it) {
		it->calcgradientbasepair();
	}
}

double getenergybasepair(list<BasePairobj> &basepairobj){
	double energy=0.0;
	for (std::list<BasePairobj>::const_iterator it = basepairobj.begin(); it != basepairobj.end(); ++it) {
		energy+=it->calcenergybasepair();
	}
	
	return energy;
}

void getgradientbaseplane(list<BasePairobj> &basepairobj){
	for (std::list<BasePairobj>::iterator it = basepairobj.begin(); it != basepairobj.end(); ++it) {
		it->calcgradientbasepair();
	}
}

double getenergybaseplane(list<BasePairobj> &baseplaneobj){
	double energy=0.0;
	for (std::list<BasePairobj>::const_iterator it = baseplaneobj.begin(); it != baseplaneobj.end(); ++it) {
		energy+=it->calcenergybasepair();
	}
	
	return energy;
}

void getgradientbasepairdist(list<BPdistobj> &bpdistobj){
	for (std::list<BPdistobj>::iterator it = bpdistobj.begin(); it != bpdistobj.end(); ++it) {
		it->calcgradientbp();
	}
}

double getenergybasepairdist(list<BPdistobj> &bpdistobj){
	double energy=0.0;
	for (std::list<BPdistobj>::const_iterator it = bpdistobj.begin(); it != bpdistobj.end(); ++it) {
		energy+=it->calcenergybp();
	}
	
	return energy;
}

void getgradienthbond(list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj){
	for (std::list<Hbondobj>::const_iterator it = hbondobj.begin(); it != hbondobj.end(); ++it) {
		it->calcgradienthbond();
	}
    for (std::list<Hbond2hydobj>::const_iterator it = hbondbpobj.begin(); it != hbondbpobj.end(); ++it) {
		it->calcgradienthbond();
	}
}

double getenergyhbond(list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj){
	double energy=0.0;
	for (std::list<Hbondobj>::const_iterator it = hbondobj.begin(); it != hbondobj.end(); ++it) {
		energy+=it->calcenergyhbond();
	}
    for (std::list<Hbond2hydobj>::const_iterator it = hbondbpobj.begin(); it != hbondbpobj.end(); ++it) {
		energy+=it->calcenergyhbond();
	}
	
	return energy;
}

void getgradientblength(list<BondLengthobj> &blengthobj){
	for (std::list<BondLengthobj>::const_iterator it = blengthobj.begin(); it != blengthobj.end(); ++it) {
		it->calcgradientlength();
	}
}

double getenergyblength(const list<BondLengthobj> &blengthobj){
	double energy=0.0;
	for (std::list<BondLengthobj>::const_iterator it = blengthobj.begin(); it != blengthobj.end(); ++it) {
		energy+=it->calcenergylength();
	}
	
	return energy;
}

void getgradientstacking(list<Stackingobj> &stackobj){
	for (std::list<Stackingobj>::const_iterator it = stackobj.begin(); it != stackobj.end(); ++it) {
		it->calcgradientstack();
	}
}

double getenergystacking(list<Stackingobj> &stackobj){
	double energy=0.0;
	for (std::list<Stackingobj>::const_iterator it = stackobj.begin(); it != stackobj.end(); ++it) {
		energy+=it->calcenergystack();
        //cout<<it->calcenergystack()<<endl;
	}
	
	return energy;
}

void getgradientstackinggeo(list<StackGeoobjMC> &stackgeoobj){
	for (std::list<StackGeoobjMC>::iterator it = stackgeoobj.begin(); it != stackgeoobj.end(); ++it) {
		it->calcgradient();
	}
}

double getenergystackinggeo(list<StackGeoobjMC> &stackgeoobj){
	double energy=0.0;
	for (std::list<StackGeoobjMC>::const_iterator it = stackgeoobj.begin(); it != stackgeoobj.end(); ++it) {
		energy+=it->calcenergy();
        //cout<<it->calcenergystack()<<endl;
	}
	
	return energy;
}

void getgradientbangle(list<BondAngleobj> &bangleobj){
	for (std::list<BondAngleobj>::const_iterator it = bangleobj.begin(); it != bangleobj.end(); ++it) {
		it->calcgradientangle();
	}
}

void getgradienthydtorsion(list<BondDihedralobj> &dihedralobj){
	for (std::list<BondDihedralobj>::const_iterator it = dihedralobj.begin(); it != dihedralobj.end(); ++it) {
		it->calcgradientdihedral();
	}
}

double getenergybangle(list<BondAngleobj> &bangleobj){
	double energy=0.0;
	for (std::list<BondAngleobj>::const_iterator it = bangleobj.begin(); it != bangleobj.end(); ++it) {
		energy+=it->calcenergyangle();
	}
	
	return energy;
}

double getenergyhydtor(list<BondDihedralobj> &hydtorobj){
    double energy=0.0;
	for (std::list<BondDihedralobj>::const_iterator it = hydtorobj.begin(); it != hydtorobj.end(); ++it) {
		energy+=it->calcenergydihedral();
	}
	
	return energy;
}

void getgradientbimpropers(list<BondImpropersobj> &improperobj){
	for (std::list<BondImpropersobj>::const_iterator it = improperobj.begin(); it != improperobj.end(); ++it) {
		it->calcgradientimpropers();
	}
}

double getenergybimpropers(list<BondImpropersobj> &improperobj){
	double energy=0.0;
	for (std::list<BondImpropersobj>::const_iterator it = improperobj.begin(); it != improperobj.end(); ++it) {
		energy+=it->calcenergyimpropers();
	}
	
	return energy;
}

void getgradientbdihedral(list<BondDihedralobj> &dihedralobj){
	for (std::list<BondDihedralobj>::const_iterator it = dihedralobj.begin(); it != dihedralobj.end(); ++it) {
		it->calcgradientdihedral();
	}
}

double getenergybdihedral(list<BondDihedralobj> &dihedralobj){
	double energy=0.0;
	for (std::list<BondDihedralobj>::const_iterator it = dihedralobj.begin(); it != dihedralobj.end(); ++it) {
		energy+=it->calcenergydihedral();
	}
	
	return energy;
}

void getgradientelectros(list<TElectrobj> &telecteobj){
	for (std::list<TElectrobj>::const_iterator it = telecteobj.begin(); it != telecteobj.end(); ++it) {
		it->calcgradientelectros();
	}
}

double getenergyelectros(list<TElectrobj> &telecteobj){
	double energy=0.0;
	for (std::list<TElectrobj>::const_iterator it = telecteobj.begin(); it != telecteobj.end(); ++it) {
		energy+=it->calcenergyelectros();
	}
	
	return energy;
}

void getgradientlj(list<LennardsJobj> &ljobj){
	for (std::list<LennardsJobj>::const_iterator it = ljobj.begin(); it != ljobj.end(); ++it) {
		it->calcgradientlj();
	}
}

double getenergylj(list<LennardsJobj> &ljobj){
	double energy=0.0;
	for (std::list<LennardsJobj>::const_iterator it = ljobj.begin(); it != ljobj.end(); ++it) {
		energy+=it->calcenergylj();
	}
	
	return energy;
}

void getgradientclash(list<Clashobj> &clashobj){
	for (std::list<Clashobj>::const_iterator it = clashobj.begin(); it != clashobj.end(); ++it) {
		it->calcgradientclash();
	}
}

double getenergyclash(list<Clashobj> &clashobj){
	double energy=0.0;
	for (std::list<Clashobj>::const_iterator it = clashobj.begin(); it != clashobj.end(); ++it) {
		energy+=it->calcenergyclash();
	}
	
	return energy;
}

void getgradientnttemp(list<SingleNtTemplateobj> &nttempobj){
	for (std::list<SingleNtTemplateobj>::const_iterator it = nttempobj.begin(); it != nttempobj.end(); ++it) {
		it->calcgradientstacktemp();
	}
}

double getenergynttemp(list<SingleNtTemplateobj> &nttempobj){
	double energy=0.0;
	for (std::list<SingleNtTemplateobj>::iterator it = nttempobj.begin(); it != nttempobj.end(); ++it) {
		energy+=it->calcenergystacktemp();
	}
	
	return energy;
}


void getgradient(list<BondLengthobj> &blengthobj, list<BondAngleobj> &bangleobj, list<BondDihedralobj> &dihedralobj, list<BasePairobj> &baseplaneobj, list<BasePairobj> &basepairobj, list<BPdistobj> &bpdistobj, list<Stackingobj> &stackobj, list<StackGeoobjMC> &stackgeoobj, list<InitPosobj> &binitobj, list<InitPairobj> &binitpairobj, list<InitPairobj> &btemplatepairobj, list<TElectrobj> &telecteobj, list<LennardsJobj> &ljobj, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, list<Clashobj> &clashobj, list<SingleNtTemplateobj> &nttempobj)
{
	if(blengthobj.size() > 0) getgradientblength(blengthobj);
	if(bangleobj.size() > 0) getgradientbangle(bangleobj);
    if(dihedralobj.size() > 0) getgradienthydtorsion(dihedralobj);
	if(basepairobj.size() > 0) getgradientbasepair(basepairobj);
    if(bpdistobj.size() > 0) getgradientbasepairdist(bpdistobj);
    if(binitobj.size() > 0) getgradientinit(binitobj);
    if(binitpairobj.size() > 0) getgradientinitpair(binitpairobj);
    if(btemplatepairobj.size() > 0) getgradienttemplatepair(btemplatepairobj);
    if(telecteobj.size() > 0) getgradientelectros(telecteobj);
    if(ljobj.size() > 0) getgradientlj(ljobj);
    if(hbondobj.size() > 0) getgradienthbond(hbondobj, hbondbpobj);
    if(stackobj.size() > 0) getgradientstacking(stackobj);
    if(stackgeoobj.size() > 0) getgradientstackinggeo(stackgeoobj);
    if(clashobj.size() > 0) getgradientclash(clashobj);
    if(nttempobj.size() > 0) getgradientnttemp(nttempobj);
}

void getgradientloop(list<BondLengthobj> &blengthobj, list<BondAngleobj> &bangleobj, list<BondDihedralobj> &dihedralobj, list<InitPosobj> &binitobj, list<InitPairobj> &binitpairobj, list<InitPairobj> &btemplatepairobj, list<TElectrobj> &telecteobj, list<LennardsJobj> &ljobj, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, list<Clashobj> &clashobj, list<SingleNtTemplateobj> &nttempobj)
{
	if(blengthobj.size() > 0) getgradientblength(blengthobj);
	if(bangleobj.size() > 0) getgradientbangle(bangleobj);
    if(dihedralobj.size() > 0) getgradienthydtorsion(dihedralobj);
    if(binitobj.size() > 0) getgradientinit(binitobj);
    if(binitpairobj.size() > 0) getgradientinitpair(binitpairobj);
    if(btemplatepairobj.size() > 0) getgradienttemplatepair(btemplatepairobj);
    if(telecteobj.size() > 0) getgradientelectros(telecteobj);
    if(ljobj.size() > 0) getgradientlj(ljobj);
    if(hbondobj.size() > 0) getgradienthbond(hbondobj, hbondbpobj);
    if(clashobj.size() > 0) getgradientclash(clashobj);
    if(nttempobj.size() > 0) getgradientnttemp(nttempobj);
}

double getcovalentenergystem(list<BondLengthobj> &blengthobj, list<BondAngleobj> &bangleobj, list<BondDihedralobj> &dihedralobj, list<BasePairobj> &baseplaneobj, list<BasePairobj> &basepairobj, list<BPdistobj> &bpdistobj, list<Stackingobj> &stackobj, list<StackGeoobjMC> &stackgeoobj, list<InitPosobj> &binitobj, list<InitPairobj> &binitpairobj, list<InitPairobj> &btemplatepairobj, list<TElectrobj> &telecteobj, list<LennardsJobj> &ljobj, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, list<Clashobj> &clashobj, list<SingleNtTemplateobj> &nttempobj)
{
    double lengthene, angleene,hydtorene, bpene, bpdistene, initene, initpairene, templatepairene, telecteene, ljene, hbondene, stackene, clashene, nttempene, stackgeoene;
    lengthene=0; angleene=0; hydtorene=0; bpene=0; bpdistene=0; initene=0; initpairene=0; templatepairene=0; telecteene=0; ljene=0; hbondene=0; stackene=0; stackgeoene=0; clashene=0; nttempene=0;
    if(blengthobj.size() > 0) lengthene=getenergyblength(blengthobj);
    if(bangleobj.size() > 0) angleene=getenergybangle(bangleobj);
    if(dihedralobj.size() > 0) hydtorene=getenergyhydtor(dihedralobj);
    if(basepairobj.size() > 0) bpene=getenergybasepair(basepairobj);
    if(bpdistobj.size() > 0) bpdistene=getenergybasepairdist(bpdistobj);
    if(binitobj.size() > 0) initene=getenergyinit(binitobj);
    if(binitpairobj.size() > 0) initpairene=getenergyinitpair(binitpairobj);
    if(btemplatepairobj.size() > 0) templatepairene=getenergytemplatepair(btemplatepairobj);
    if(telecteobj.size() > 0) telecteene=getenergyelectros(telecteobj);
    //cout<<"telecteene"<<endl;
    //cout<<telecteene<<endl;
    if(ljobj.size() > 0) ljene=getenergylj(ljobj);
    //cout<<"ljene"<<endl;
    //cout<<ljene<<endl;
    if(hbondobj.size() > 0) hbondene=getenergyhbond(hbondobj, hbondbpobj);
    if(stackobj.size() > 0) stackene=getenergystacking(stackobj);
    if(stackgeoobj.size() > 0) stackgeoene=getenergystackinggeo(stackgeoobj);
    if(clashobj.size() > 0) clashene=getenergyclash(clashobj);
    if(nttempobj.size() > 0) nttempene=getenergynttemp(nttempobj);
    //cout<<"aaa-" << lengthene << "-" << angleene << "-" << initene << "-" << bpene << "-" << telecteene<< "-" << ljene << "-" << hbondene << "-" << stackene << "-" << clashene<<endl;
    //exit(1);

    //cout<<"init dist:" <<  initpairene<<endl;
    //cout<<"Stack:" <<  stackene<<endl;
    //cout<<lengthene+angleene+hydtorene+bpene+initene+initpairene+templatepairene+telecteene+ljene+hbondene+stackene+clashene<<endl;

    return lengthene+angleene+hydtorene+bpene+bpdistene+initene+initpairene+templatepairene+telecteene+ljene+hbondene+stackene+stackgeoene+clashene+nttempene;
}

double getcovalentenergystemloop(list<BondLengthobj> &blengthobj, list<BondAngleobj> &bangleobj, list<BondDihedralobj> &dihedralobj, list<BasePairobj> &baseplaneobj, list<BasePairobj> &basepairobj, list<Stackingobj> &stackobj, list<StackGeoobjMC> &stackgeoobj, list<InitPosobj> &binitobj, list<InitPairobj> &binitpairobj, list<InitPairobj> &btemplatepairobj, list<TElectrobj> &telecteobj, list<LennardsJobj> &ljobj, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, list<Clashobj> &clashobj, list<SingleNtTemplateobj> &nttempobj)
{
    double lengthene, angleene,hydtorene, bpene, initene, initpairene, templatepairene, telecteene, ljene, hbondene, stackene, clashene, nttempene, stackgeoene;
    lengthene=0; angleene=0; hydtorene=0; bpene=0; initene=0; initpairene=0; templatepairene=0; telecteene=0; ljene=0; hbondene=0; stackene=0; clashene=0; nttempene=0; stackgeoene=0;
    if(blengthobj.size() > 0) lengthene=getenergyblength(blengthobj);
    if(bangleobj.size() > 0) angleene=getenergybangle(bangleobj);
    if(dihedralobj.size() > 0) hydtorene=getenergyhydtor(dihedralobj);
    if(basepairobj.size() > 0) bpene=getenergybasepair(basepairobj);
    if(binitobj.size() > 0) initene=getenergyinit(binitobj);
    if(binitpairobj.size() > 0) initpairene=getenergyinitpair(binitpairobj);
    if(btemplatepairobj.size() > 0) templatepairene=getenergytemplatepair(btemplatepairobj);
    if(telecteobj.size() > 0) telecteene=getenergyelectros(telecteobj);
    if(ljobj.size() > 0) ljene=getenergylj(ljobj);
    if(hbondobj.size() > 0) hbondene=getenergyhbond(hbondobj, hbondbpobj);
    if(stackobj.size() > 0) stackene=getenergystacking(stackobj);
    if(stackgeoobj.size() > 0) stackgeoene=getenergystackinggeo(stackgeoobj);
    if(clashobj.size() > 0) clashene=getenergyclash(clashobj);
    if(nttempobj.size() > 0) nttempene=getenergynttemp(nttempobj);
    //cout<<"bbb-" << lengthene << "-" << angleene << "-" << bpene << "-" << initene << "-" << telecteene<< "-" << ljene << "-" << hbondene << "-" << stackene << "-" << clashene <<endl;

    return lengthene+angleene+hydtorene+bpene+initene+initpairene+templatepairene+telecteene+ljene+hbondene+stackene+stackgeoene+clashene+nttempene;
}

double calcaveragenorm(vector<vector<Point3s>> &gradient, double clip_val)
{
	int i, mtype;
	double dx, dy, dz;
	double totalnorm = 0;
  	unsigned int num = 0;
    double norm, scale;



	for(i=0;i<gradient.size();i++){
		for(mtype=0;mtype<gradient[i].size();mtype++)
		{
			dx=gradient[i][mtype].x;
			dy=gradient[i][mtype].y;
			dz=gradient[i][mtype].z;

			totalnorm+=sqrt(dx*dx+dy*dy+dz*dz);
			num+=1;
		}
	}

	scale = sqrt(totalnorm/double(num)) + 1e-8;

    for(i=0;i<gradient.size();i++){
		for(mtype=0;mtype<gradient[i].size();mtype++)
		{
            gradient[i][mtype].x/=scale;
			gradient[i][mtype].y/=scale;
			gradient[i][mtype].z/=scale;
		}
	}

    //double maxgrad = -1.0;

    if(clip_val > 0){
        for(i=0;i<gradient.size();i++){
		    for(mtype=0;mtype<gradient[i].size();mtype++)
            {
                dx=gradient[i][mtype].x;
                dy=gradient[i][mtype].y;
                dz=gradient[i][mtype].z;

                norm=sqrt(dx*dx+dy*dy+dz*dz);

                //if(norm > 10000) cout<<norm<<endl;

                //if(maxgrad < norm) maxgrad = norm;

                if(norm > clip_val){
                    scale = clip_val/norm;
                    gradient[i][mtype].x*=scale;
                    gradient[i][mtype].y*=scale;
                    gradient[i][mtype].z*=scale;
                }
            }
        }
    }

    return scale;
}

void updatestructure(FragmentStructure &stru, const vector<vector<Point3s>> &gradient, double step_size, vector<bool> &moveresidues){
	int i, j;
	char aat;
	for(i=0;i<stru.residues.size();i++){
        if(!moveresidues[i]) continue;
		for(j=0;j<stru.residues[i].atoms.size();j++)
		{
			//stru.residues[i].atoms[j].coord.x += step_size*gradient[i][j].x/averagenorm;
			//stru.residues[i].atoms[j].coord.y += step_size*gradient[i][j].y/averagenorm;
			//stru.residues[i].atoms[j].coord.z += step_size*gradient[i][j].z/averagenorm;
            //cout<<"-------------------"<<endl;
            //cout<<step_size*gradient[i][j].x << "-" << step_size*gradient[i][j].y << "-" << step_size*gradient[i][j].z<<endl;
            //cout<<gradient[i][j].x << "-" << gradient[i][j].y << "-" << gradient[i][j].z<<endl;
            stru.residues[i].atoms[j].coord.x += step_size*gradient[i][j].x;
			stru.residues[i].atoms[j].coord.y += step_size*gradient[i][j].y;
			stru.residues[i].atoms[j].coord.z += step_size*gradient[i][j].z;
		}
	}
}

double goldenelectionsearch(FragmentStructure &stru, const vector<vector<Point3s>> &gradient, list<BondLengthobj> &blengthobj, list<BondAngleobj> &bangleobj, list<BondDihedralobj> &dihedralobj, list<BasePairobj> &baseplaneobj, list<BasePairobj> &basepairobj, list<BPdistobj> &bpdistobj, list<Stackingobj> &stackobj, list<StackGeoobjMC> &stackgeoobj, list<InitPosobj> &binitobj, list<InitPairobj> &binitpairobj, list<InitPairobj> &btemplatepairobj, list<TElectrobj> &telecteobj, list<LennardsJobj> &ljobj, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, list<Clashobj> &clashobj, list<SingleNtTemplateobj> &nttempobj, vector<bool> &moveresidues)
{
	static const double goldenratio = 0.5*(sqrt(5.0)-1.);
	double step = 1.0e-4;
	double left = step*(1.0-goldenratio);
	double right = step*goldenratio;

	double energyst = getcovalentenergystem(blengthobj, bangleobj, dihedralobj, baseplaneobj, basepairobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);
	
	updatestructure(stru, gradient, -right, moveresidues);
	double energyR = getcovalentenergystem(blengthobj, bangleobj, dihedralobj, baseplaneobj, basepairobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);

	updatestructure(stru, gradient, right-left, moveresidues);
	double energyL = getcovalentenergystem(blengthobj, bangleobj, dihedralobj, baseplaneobj, basepairobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);

    //cout<<energyst<< "-" << energyR << "-" << energyL << "-" << step<<endl;

	int stepnum = 0;
	while(fabs(step) > 1e-8)
  	{
        if(energyL < energyR || (energyL > energyst && energyR > energyst))
    	{
			step = step - step*(1-goldenratio);
			right = left;
			energyR = energyL;
			left = step*(1.-goldenratio);
			updatestructure(stru, gradient, right-left, moveresidues);
			energyL = getcovalentenergystem(blengthobj, bangleobj, dihedralobj, baseplaneobj, basepairobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);
		
            //cout<<"L-" << energyst<< "-" << energyR << "-" << energyL << "-" << step<<endl;
        
        }else
        {
            //cout<<step << "-" << left<<endl;
            step = step - step*(1-goldenratio);
            left = right;
			energyL = energyR;
            right = step*goldenratio;
            updatestructure(stru, gradient, -right, moveresidues);
			energyR = getcovalentenergystem(blengthobj, bangleobj, dihedralobj, baseplaneobj, basepairobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);
			updatestructure(stru, gradient, right-left, moveresidues);

            //cout<<"R-" << energyst<< "-" << energyR << "-" << energyL << "-" << step<<endl;
		}
	}

    //exit(1);
	updatestructure(stru, gradient, left-step/2, moveresidues);
    return (energyL+energyR)/2.;
}

double linesearch(FragmentStructure &stru, const vector<vector<Point3s>> &gradient, list<BondLengthobj> &blengthobj, list<BondAngleobj> &bangleobj, list<BondDihedralobj> &dihedralobj, list<BasePairobj> &baseplaneobj, list<BasePairobj> &basepairobj, list<BPdistobj> &bpdistobj, list<Stackingobj> &stackobj, list<StackGeoobjMC> &stackgeoobj, list<InitPosobj> &binitobj, list<InitPairobj> &binitpairobj, list<InitPairobj> &btemplatepairobj, list<TElectrobj> &telecteobj, list<LennardsJobj> &ljobj, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, list<Clashobj> &clashobj, list<SingleNtTemplateobj> &nttempobj, vector<bool> &moveresidues, double alpha0, double rho, double beta)
{
    double alpha = alpha0;

    double new_energy;
    double init_energy =  getcovalentenergystem(blengthobj, bangleobj, dihedralobj, baseplaneobj, basepairobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);

    double grad_norm_sq = 0.0;
    for (int i=0;i<gradient.size();i++) {
        for (int j=0;j<gradient[i].size();j++) {
            grad_norm_sq += gradient[i][j].x * gradient[i][j].x + gradient[i][j].y * gradient[i][j].y + gradient[i][j].z * gradient[i][j].z;
        }
    }

    //cout<<"============"<<endl;
    while (true) {
        //cout<<"linesearch-" << alpha<<endl;
        updatestructure(stru, gradient, -alpha, moveresidues);
        new_energy =  getcovalentenergystem(blengthobj, bangleobj, dihedralobj, baseplaneobj, basepairobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);

        if (new_energy <= init_energy) {
            break;
        }else{
            updatestructure(stru, gradient, alpha, moveresidues);
        }

        alpha *= rho;
        if(alpha < 1e-8) break;
    }

    return new_energy;
}

double goldenelectionsearchloop(FragmentStructure &stru, const vector<vector<Point3s>> &gradient, list<BondLengthobj> &blengthobj, list<BondAngleobj> &bangleobj, list<BondDihedralobj> &dihedralobj, list<BasePairobj> &baseplaneobj, list<BasePairobj> &bpobj, list<Stackingobj> &stackobj, list<StackGeoobjMC> &stackgeoobj, list<InitPosobj> &binitobj, list<InitPairobj> &binitpairobj, list<InitPairobj> &btemplatepairobj, list<TElectrobj> &telecteobj, list<LennardsJobj> &ljobj, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, list<Clashobj> &clashobj, list<SingleNtTemplateobj> &nttempobj, vector<bool> &moveresidues)
{
	static const double goldenratio = .5*(sqrt(5.)-1.);
	double step = 1.e-4;
	double left = step*(1.-goldenratio);
	double right = step*goldenratio;

	double energyst = getcovalentenergystemloop(blengthobj, bangleobj, dihedralobj, baseplaneobj, bpobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);
	
	updatestructure(stru, gradient, -right, moveresidues);
	double energyR = getcovalentenergystemloop(blengthobj, bangleobj, dihedralobj, baseplaneobj, bpobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);

	updatestructure(stru, gradient, right-left, moveresidues);
	double energyL = getcovalentenergystemloop(blengthobj, bangleobj, dihedralobj, baseplaneobj, bpobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);

	int stepnum = 0;
	while(fabs(step) > 1.e-8)
  	{
		if(energyL < energyR || (energyL > energyst && energyR > energyst))
    	{
			step = right;
			right = left;
			energyR = energyL;
			left = step*(1.-goldenratio);
			updatestructure(stru, gradient, right-left, moveresidues);
			energyL = getcovalentenergystemloop(blengthobj, bangleobj, dihedralobj, baseplaneobj, bpobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);
		}else{
			step = step - left;
			left = right;
			energyL = energyR;
			right = step*goldenratio;
			updatestructure(stru, gradient, -right, moveresidues);
			energyR = getcovalentenergystemloop(blengthobj, bangleobj, dihedralobj, baseplaneobj, bpobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);
			updatestructure(stru, gradient, right-left, moveresidues);
		}
	}

	updatestructure(stru, gradient, left-step/2, moveresidues);
	return (energyL+energyR)/2.;
}

void getljpairs(FragmentStructure &superimposedstr, vector<vector<Point3s>> &gradient, const map<string, map<int, vector<double>>> &atomfeatures, int res1, int res2, vector<vector<int>> molConnect, list<LennardsJobj> &ljobj)
{
    int k, l, m;
	int atomnum1, atomnum2;
	string base1, base2;
    int atom1, atom2;
    Point3s vec12;
    double vdWInfodR01, vdWInfodR02, vdWInfodEpsilon1, vdWInfodEpsilon2, star, epsilon, LJ6, LJ12,Kval;
    
    base1 = superimposedstr.residues[res1].restype;
	base2 = superimposedstr.residues[res2].restype;
	atomnum1 = superimposedstr.residues[res1].atoms.size();
	atomnum2 = superimposedstr.residues[res2].atoms.size();

    if(res1 != res2){
        for(k=0;k<atomnum1;k++){
			for(l=0;l<atomnum2;l++){
                if(molConnect[k][l] >= 2){
                    const Point3s &pt1 = superimposedstr.residues[res1].atoms[k].coord;
					const Point3s &pt2 = superimposedstr.residues[res2].atoms[l].coord;

                    vec12.x = pt2.x - pt1.x;
					vec12.y = pt2.y - pt1.y;
					vec12.z = pt2.z - pt1.z;

                    double length12=vec12.x*vec12.x+vec12.y*vec12.y+vec12.z*vec12.z;

                    if(length12 < CUTOFF_SQR+10.0){
                        vdWInfodR01 = atomfeatures.at(base1).at(k)[1];
                        vdWInfodR02 = atomfeatures.at(base2).at(l)[1];
                        vdWInfodEpsilon1 = atomfeatures.at(base1).at(k)[2];
                        vdWInfodEpsilon2 = atomfeatures.at(base2).at(l)[2];
                        star = (vdWInfodR01 + vdWInfodR02);
                        epsilon =  sqrt(vdWInfodEpsilon1*vdWInfodEpsilon2);
                        LJ6 = -2*epsilon*pow(star,6);
                        LJ12 = epsilon*pow(star,12);

                        //cout<<full_atoms[base1][k] << "-" << full_atoms[base2][l] << "-" << res1 << "-" << res2 << "-" << LJ6 << "-" << LJ12<<endl;
                        
                        //if(molConnect[k][l] == 3){
							//LJ6 /=2.0;
							//LJ12 /=2.0;
						//}
                        
                        //if(stackmap[res1][res2] > 0){
                            //if(k > 12 && k < full_atoms_name[base1].size() && l > 12 && l < full_atoms_name[base2].size()) Kval = 1.2;
                            //else Kval = 1.0;
                        //}else{
                            //Kval = 1.0;
                        //}
                        Kval = 1.0;
                        LennardsJobj ljobjone = LennardsJobj::create(&superimposedstr.residues[res1].atoms[k].coord, &superimposedstr.residues[res2].atoms[l].coord, &gradient[res1][k], &gradient[res2][l], LJ6, LJ12, Kval);
			            ljobj.push_back(ljobjone);
                    }
                }
            }
        }
    }else{
        for(k=0;k<atomnum1;k++){
			for(l=k+1;l<atomnum2;l++){
                if(molConnect[k][l] >= 3){
                    const Point3s &pt1 = superimposedstr.residues[res1].atoms[k].coord;
					const Point3s &pt2 = superimposedstr.residues[res2].atoms[l].coord;

                    vec12.x = pt2.x - pt1.x;
					vec12.y = pt2.y - pt1.y;
					vec12.z = pt2.z - pt1.z;

                    double length12=vec12.x*vec12.x+vec12.y*vec12.y+vec12.z*vec12.z;

                    if(length12 < CUTOFF_SQR+5.0){
                        vdWInfodR01 = atomfeatures.at(base1).at(k)[1];
                        vdWInfodR02 = atomfeatures.at(base2).at(l)[1];
                        vdWInfodEpsilon1 = atomfeatures.at(base1).at(k)[2];
                        vdWInfodEpsilon2 = atomfeatures.at(base2).at(l)[2];
                        star = (vdWInfodR01 + vdWInfodR02);
                        epsilon =  sqrt(vdWInfodEpsilon1*vdWInfodEpsilon2);
                        LJ6 = -2.0*epsilon*pow(star,6);
                        LJ12 = epsilon*pow(star,12);

                        if(molConnect[k][l] == 3){
							LJ6 /=2.0;
							LJ12 /=2.0;
						}

                        //cout<<full_atoms[base1][k] << "-" << full_atoms[base2][l] << "-" << res1 << "-" << res2 << "-" << LJ6 << "-" << LJ12<<endl;

                        //if(stackmap[res1][res2] > 0){
                            //if(k > 12 && k < full_atoms_name[base1].size() && l > 12 && l < full_atoms_name[base2].size()) Kval = 1.2;
                            //else Kval = 1.0;
                        //}else{
                            //Kval = 1.0;
                        //}
                        Kval = 1.0;
                        LennardsJobj ljobjone = LennardsJobj::create(&superimposedstr.residues[res1].atoms[k].coord, &superimposedstr.residues[res2].atoms[l].coord, &gradient[res1][k], &gradient[res2][l], LJ6, LJ12, Kval);
			            ljobj.push_back(ljobjone);
                    }
                }
            }
        }
    }

}

void getclashpairs(FragmentStructure &superimposedstr, vector<vector<Point3s>> &gradient, const map<string, map<int, vector<double>>> &atomfeatures, int res1, int res2, vector<vector<int>> molConnect, double K, list<Clashobj> &clashobj)
//void getclashpairs(FragmentStructure &superimposedstr, vector<vector<Point3s>> &gradient, const map<string, vector<double>> &clashparam, int res1, int res2, vector<vector<int>> molConnect, double K, list<Clashobj> &clashobj)
{
    int k, l, m;
	int atomnum1, atomnum2;
	string base1, base2;
    int atom1, atom2;
    Point3s vec12;
    double vdWInfodR01, vdWInfodR02, vdW;
    
    base1 = superimposedstr.residues[res1].restype;
	base2 = superimposedstr.residues[res2].restype;
	atomnum1 = superimposedstr.residues[res1].atoms.size();
	atomnum2 = superimposedstr.residues[res2].atoms.size();

    if(res1 != res2){
        for(k=0;k<atomnum1;k++){
			for(l=0;l<atomnum2;l++){
                if(molConnect[k][l] >= 2){
                    const Point3s &pt1 = superimposedstr.residues[res1].atoms[k].coord;
					const Point3s &pt2 = superimposedstr.residues[res2].atoms[l].coord;

                    vec12.x = pt2.x - pt1.x;
					vec12.y = pt2.y - pt1.y;
					vec12.z = pt2.z - pt1.z;

                    double length12=vec12.x*vec12.x+vec12.y*vec12.y+vec12.z*vec12.z;

                    if(length12 < CUTOFF_SQR+10){
                        vdWInfodR01 = atomfeatures.at(base1).at(k)[1];
                        vdWInfodR02 = atomfeatures.at(base2).at(l)[1];
                        vdW = (vdWInfodR01 + vdWInfodR02);

                        //vdW = clashparam.at(base1).at(k) + clashparam.at(base2).at(l) - 0.4;
                        
                        Clashobj clashobjone = Clashobj::create(&superimposedstr.residues[res1].atoms[k].coord, &superimposedstr.residues[res2].atoms[l].coord, &gradient[res1][k], &gradient[res2][l], vdW, K);
			            clashobj.push_back(clashobjone);
                    }
                }
            }
        }
    }else{
        for(k=0;k<atomnum1;k++){
			for(l=k+1;l<atomnum2;l++){
                if(molConnect[k][l] >= 3){
                    const Point3s &pt1 = superimposedstr.residues[res1].atoms[k].coord;
					const Point3s &pt2 = superimposedstr.residues[res2].atoms[l].coord;

                    vec12.x = pt2.x - pt1.x;
					vec12.y = pt2.y - pt1.y;
					vec12.z = pt2.z - pt1.z;

                    double length12=vec12.x*vec12.x+vec12.y*vec12.y+vec12.z*vec12.z;

                    if(length12 < CUTOFF_SQR+5){
                        vdWInfodR01 = atomfeatures.at(base1).at(k)[1];
                        vdWInfodR02 = atomfeatures.at(base2).at(l)[1];
                        vdW = (vdWInfodR01 + vdWInfodR02);
                        
                        //vdW = clashparam.at(base1).at(k) + clashparam.at(base2).at(l) - 0.4;
                        
                        Clashobj clashobjone = Clashobj::create(&superimposedstr.residues[res1].atoms[k].coord, &superimposedstr.residues[res2].atoms[l].coord, &gradient[res1][k], &gradient[res2][l], vdW, K);
			            clashobj.push_back(clashobjone);
                    }
                }
            }
        }
    }

}


void getelectrospairs(FragmentStructure &superimposedstr, vector<vector<Point3s>> &gradient, const map<string, map<int, vector<double>>> &atomfeatures, int res1, int res2, vector<vector<int>> molConnect, list<TElectrobj> &telecteobj)
{
	int k, l;
	int atomnum1, atomnum2;
	string base1, base2;
    int atom1, atom2;
    double dQQ, dQQep;
	
	base1 = superimposedstr.residues[res1].restype;
	base2 = superimposedstr.residues[res2].restype;
	atomnum1 = superimposedstr.residues[res1].atoms.size();
	atomnum2 = superimposedstr.residues[res2].atoms.size();

	if(res1 != res2){
		for(k=0;k<atomnum1;k++){
			for(l=0;l<atomnum2;l++){
				if(molConnect[k][l] >= 2){
                    atom1 = k; atom2 = l;

                    //if(molConnect[k][l] == 3){
                        //dQQ = (COULOMB*atomfeatures.at(base1).at(k)[0]*atomfeatures.at(base2).at(l)[0])/1.2;
					//}else{
					dQQ = COULOMB*atomfeatures.at(base1).at(k)[0]*atomfeatures.at(base2).at(l)[0];
					//}
                    dQQep = dQQ/EpsilonRNA;

                    TElectrobj telecteobjone = TElectrobj::create(&superimposedstr.residues[res1].atoms[k].coord, &superimposedstr.residues[res2].atoms[l].coord, &gradient[res1][k], &gradient[res2][l], dQQep);
			        telecteobj.push_back(telecteobjone);
				}
			}
		}
	}else{
		for(k=0;k<atomnum1;k++){
			for(l=k+1;l<atomnum2;l++){
				if(molConnect[k][l] >= 3){
					atom1 = k; atom2 = l;

                    if(molConnect[k][l] == 3){
						dQQ = (COULOMB*atomfeatures.at(base1).at(k)[0]*atomfeatures.at(base2).at(l)[0])/1.2;
					}else{
						dQQ = COULOMB*atomfeatures.at(base1).at(k)[0]*atomfeatures.at(base2).at(l)[0];
					}
                    dQQep = dQQ/EpsilonRNA;

                    TElectrobj telecteobjone = TElectrobj::create(&superimposedstr.residues[res1].atoms[k].coord, &superimposedstr.residues[res2].atoms[l].coord, &gradient[res1][k], &gradient[res2][l], dQQep);
			        telecteobj.push_back(telecteobjone);
				}
			}
		}
	}
}

double calcenergyhbond(const Point3s &teAtomHyd, const Point3s &teAtomDon, const Point3s &teAtomAcc)
{
	double dREq=1.55;
	double dK=1.0;
	double dXiR=2.0;
	double dXiTheta=2.0;
	double dThreshold=-1.0;
	Point3s vec21, vec23, vec13, grad1, grad2, grad3, gradsum;

	vec21.x = teAtomDon.x - teAtomHyd.x;
	vec21.y = teAtomDon.y - teAtomHyd.y;
	vec21.z = teAtomDon.z - teAtomHyd.z;

	vec23.x = teAtomAcc.x - teAtomHyd.x;
	vec23.y = teAtomAcc.y - teAtomHyd.y;
	vec23.z = teAtomAcc.z - teAtomHyd.z;

	double length21=sqrt(vec21.x*vec21.x+vec21.y*vec21.y+vec21.z*vec21.z);
	double length23=sqrt(vec23.x*vec23.x+vec23.y*vec23.y+vec23.z*vec23.z);

	vec21.x /= length21; vec21.y /= length21; vec21.z /= length21;
	vec23.x /= length23; vec23.y /= length23; vec23.z /= length23;
	
	double innerprot = vec21.x*vec23.x+vec21.y*vec23.y+vec21.z*vec23.z;

	double dcosangle = innerprot;
	if     (dcosangle >  1.) dcosangle =  1.;
  	else if(dcosangle < -1.) dcosangle = -1.;

	return -dK*exp(-dXiR*(length23-dREq)*(length23-dREq) - dXiTheta*dcosangle);
}

/*
void gethbondpairs(FragmentStructure &superimposedstr, vector<vector<Point3s>> &gradient, int res1, int res2, vector<vector<int>> molConnect, vector<vector<bool>> intraatomConnect1, vector<vector<bool>> intraatomConnect2, list<Hbondobj> &hbondobj, const double &K)
{
	int k, l, m;
	int atomnum1, atomnum2;
	string base1, base2;
	Point3s vec12;
	double energy;

    //map<string, vector<int>> HBondHydrogenAtoms = {
        //{"A",{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,0}},
        //{"C",{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1}},
        //{"G",{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,1}},
        //{"U",{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1}}
    //};

    map<string, vector<int>> HBondHydrogenAtoms = {
        {"A",{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0}},
        {"C",{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0}},
        {"G",{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0}},
        {"U",{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0}}
    };
    
    map<string, vector<int>> HBondAcceptorAtoms = {
        {"A",{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0}},
        {"C",{0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}},
        {"G",{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0}},
        {"U",{0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0}}
    };
    

    base1 = superimposedstr.residues[res1].restype;
	base2 = superimposedstr.residues[res2].restype;
	atomnum1 = superimposedstr.residues[res1].atoms.size();
	atomnum2 = superimposedstr.residues[res2].atoms.size();

    for(k=0;k<atomnum1;k++){
        for(l=0;l<atomnum2;l++){
            //cout<<"ccc-" << k << "-" << l << "-" << HBondHydrogenAtoms[base1].size() << "-" << HBondAcceptorAtoms[base2].size()<<endl;
            //cout<<molConnect.size()<<endl;
            //cout<<molConnect[k].size()<<endl;
            if((HBondHydrogenAtoms[base1][k] == 1) && (HBondAcceptorAtoms[base2][l] == 1) && (molConnect[k][l] == 4)){
                for(m=0;m<atomnum1;m++){
                    //cout<<k << "-" << m << "-" << intraatomConnect.size()<<endl;
                    if(intraatomConnect1[k][m]){
                        //cout<<"bbb"<<endl;
                        const Point3s &ptH = superimposedstr.residues[res1].atoms[k].coord;
                        const Point3s &ptAcc = superimposedstr.residues[res2].atoms[l].coord;
                        const Point3s &ptDon = superimposedstr.residues[res1].atoms[m].coord;

                        //cout<<"aaa"<<endl;

                        energy=calcenergyhbond(ptH, ptDon, ptAcc);

                        //cout<<energy<<endl;
                        //cout<<base1 << "-"<< base2 << "-"<< res1 << "-"<< res2 << "-" << ptH.x << "-" << ptH.y << "-" << ptH.z << "-" << ptAcc.x << "-" << ptAcc.y << "-" << ptAcc.z << "-" << ptDon.x << "-" << ptDon.y << "-" << ptDon.z<<endl;
                        if(energy < HBONDThreshold){
                            Hbondobj hbondobjone = Hbondobj::create(&superimposedstr.residues[res1].atoms[k].coord, &superimposedstr.residues[res1].atoms[m].coord, &superimposedstr.residues[res2].atoms[l].coord, &gradient[res1][k], &gradient[res1][m], &gradient[res2][l], K);
			                hbondobj.push_back(hbondobjone);
                        }
                    }
                }
            }
        }
    }

    for(k=0;k<atomnum2;k++){
        for(l=0;l<atomnum1;l++){
            if((HBondHydrogenAtoms[base2][k] == 1) && (HBondAcceptorAtoms[base1][l] == 1) && (molConnect[l][k] == 4)){
                for(m=0;m<atomnum2;m++){
                    if(intraatomConnect2[k][m]){
                        const Point3s &ptH = superimposedstr.residues[res2].atoms[k].coord;
                        const Point3s &ptAcc = superimposedstr.residues[res1].atoms[l].coord;
                        const Point3s &ptDon = superimposedstr.residues[res2].atoms[m].coord;

                        energy=calcenergyhbond(ptH, ptDon, ptAcc);

                        if(energy < HBONDThreshold){
                            Hbondobj hbondobjone = Hbondobj::create(&superimposedstr.residues[res2].atoms[k].coord, &superimposedstr.residues[res2].atoms[m].coord, &superimposedstr.residues[res1].atoms[l].coord, &gradient[res2][k], &gradient[res2][m], &gradient[res1][l], K);
			                hbondobj.push_back(hbondobjone);
                        }
                    }
                }
            }
        }
    }
}
*/

bool checkhbond2(Point3s* teAtomHyd1, Point3s* teAtomHyd2, Point3s* teAtomDon, Point3s* teAtomAcc)
{
	Point3s vec23_1, vec23_2, vec23, vec13, vec21;

	vec23_1.x = teAtomAcc->x - teAtomHyd1->x;
	vec23_1.y = teAtomAcc->y - teAtomHyd1->y;
	vec23_1.z = teAtomAcc->z - teAtomHyd1->z;

	vec23_2.x = teAtomAcc->x - teAtomHyd2->x;
	vec23_2.y = teAtomAcc->y - teAtomHyd2->y;
	vec23_2.z = teAtomAcc->z - teAtomHyd2->z;

	double length23_1=sqrt(vec23_1.x*vec23_1.x+vec23_1.y*vec23_1.y+vec23_1.z*vec23_1.z+1e-8);
	double length23_2=sqrt(vec23_2.x*vec23_2.x+vec23_2.y*vec23_2.y+vec23_2.z*vec23_2.z+1e-8);

	int hydatom;
	if(length23_1 < length23_2){
		vec21.x = teAtomDon->x - teAtomHyd1->x;
		vec21.y = teAtomDon->y - teAtomHyd1->y;
		vec21.z = teAtomDon->z - teAtomHyd1->z;

		vec23.x = vec23_1.x; vec23.y = vec23_1.y; vec23.z = vec23_1.z;

		hydatom = -1;
	}else{
		vec21.x = teAtomDon->x - teAtomHyd2->x;
		vec21.y = teAtomDon->y - teAtomHyd2->y;
		vec21.z = teAtomDon->z - teAtomHyd2->z;

		vec23.x = vec23_2.x; vec23.y = vec23_2.y; vec23.z = vec23_2.z;

		hydatom = 1;
	}

	vec13.x = teAtomDon->x - teAtomAcc->x;
	vec13.y = teAtomDon->y - teAtomAcc->y;
	vec13.z = teAtomDon->z - teAtomAcc->z;

	double length21=sqrt(vec21.x*vec21.x+vec21.y*vec21.y+vec21.z*vec21.z+1e-8);
	double length23=sqrt(vec23.x*vec23.x+vec23.y*vec23.y+vec23.z*vec23.z+1e-8);
	double length13=sqrt(vec13.x*vec13.x+vec13.y*vec13.y+vec13.z*vec13.z+1e-8);

	vec21.x /= length21; vec21.y /= length21; vec21.z /= length21;
	vec23.x /= length23; vec23.y /= length23; vec23.z /= length23;
	vec13.x /= length13; vec13.y /= length13; vec13.z /= length13;

	double innerprot = vec21.x*vec23.x+vec21.y*vec23.y+vec21.z*vec23.z;
	double dcosangle = innerprot;
	if     (dcosangle >  1.) dcosangle =  1.;
	else if(dcosangle < -1.) dcosangle = -1.;

	double angle = acos(dcosangle)*DEGRAD;

	//cout<<angle << "-" << length13 << "-" << (length13 <= 4.5)<<endl;

	//return ((angle >= 95) && (length13 <= 4.5));
	return (length13 <= 4.0);
}

bool checkhbond(Point3s* teAtomHyd, Point3s* teAtomDon, Point3s* teAtomAcc)
{
	Point3s vec21, vec23, vec13;

	vec21.x = teAtomDon->x - teAtomHyd->x;
	vec21.y = teAtomDon->y - teAtomHyd->y;
	vec21.z = teAtomDon->z - teAtomHyd->z;

	vec23.x = teAtomAcc->x - teAtomHyd->x;
	vec23.y = teAtomAcc->y - teAtomHyd->y;
	vec23.z = teAtomAcc->z - teAtomHyd->z;

	vec13.x = teAtomDon->x - teAtomAcc->x;
	vec13.y = teAtomDon->y - teAtomAcc->y;
	vec13.z = teAtomDon->z - teAtomAcc->z;

	double length21=sqrt(vec21.x*vec21.x+vec21.y*vec21.y+vec21.z*vec21.z+1e-8);
	double length23=sqrt(vec23.x*vec23.x+vec23.y*vec23.y+vec23.z*vec23.z+1e-8);
	double length13=sqrt(vec13.x*vec13.x+vec13.y*vec13.y+vec13.z*vec13.z+1e-8);

	vec21.x /= length21; vec21.y /= length21; vec21.z /= length21;
	vec23.x /= length23; vec23.y /= length23; vec23.z /= length23;
	
	double innerprot = vec21.x*vec23.x+vec21.y*vec23.y+vec21.z*vec23.z;
	double dcosangle = innerprot;
	if     (dcosangle >  1.) dcosangle =  1.;
	else if(dcosangle < -1.) dcosangle = -1.;

	double angle = acos(dcosangle)*DEGRAD;

	//cout<<angle << "-" << length13 << "-" << (length13 <= 4.5)<<endl;

	//return ((angle >= 95) && (length13 <= 4.5));
	return (length13 <= 4.5);
}

void gethbondpairsbp(FragmentStructure &superimposedstr, vector<vector<Point3s>> &gradient, int res1, int res2, vector<vector<int>> molConnect, vector<vector<bool>> intraatomConnect1, vector<vector<bool>> intraatomConnect2, const map<char, map<char, map<int, map<int, HydroBondParameters>>>> &hbondparam, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, int bptype, const double &K)
{
	int k, l, m;
	int atomnum1, atomnum2;
	string base1, base2;
    char AccBase, DonBase;
	Point3s vec12;
	double energy;
    bool flag;

    //map<string, vector<int>> HBondHydrogenAtoms = {
        //{"A",{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,0}},
        //{"C",{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1}},
        //{"G",{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,1}},
        //{"U",{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1}}
    //};
    
    map<string, vector<int>> HBondHydrogenAtoms = {
        {"A",{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0}},
        {"C",{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0}},
        {"G",{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0}},
        {"U",{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0}}
    };
    
    map<string, vector<int>> HBondAcceptorAtoms = {
        {"A",{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0}},
        {"C",{0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}},
        {"G",{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0}},
        {"U",{0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0}}
    };

    //Acc idx, Hyd idx, Don idx, order
    static map<string, vector<vector<vector<int>>>> BpHbondInfo = {
        {"AA",{{{20,30,17,0}, {20,31,17,0}}, {{18,30,17,0}, {18,31,17,0}}, {{14,30,17,0}, {14,31,17,0}}}},
        {"AC",{{{15,30,17,1}, {15,31,17,1}}, {{14,29,17,0}, {14,30,17,0}}, {{14,30,17,1}, {14,31,17,1}}, {{20,29,17,0}, {20,30,17,0}}, {{18,29,17,0}, {18,30,17,0}}}},
        {"CA",{{{15,30,17,0}, {15,31,17,0}}, {{14,29,17,1}, {14,30,17,1}}, {{14,30,17,0}, {14,31,17,0}}, {{20,29,17,1}, {20,30,17,1}}, {{18,29,17,1}, {18,30,17,1}}}},
        {"AG",{{{14,30,17,1}, {14,31,17,1}}, {{18,31,20,0}, {18,32,20,0}}, {{14,31,20,0}, {14,32,20,0}}, {{20,31,20,0}, {20,32,20,0}}, {{18,33,18,0}}, {{17,30,17,1}, {17,31,17,1}}, {{14,33,18,0}}, {{20,33,18,0}}}},
        {"GA",{{{14,30,17,0}, {14,31,17,0}}, {{18,31,20,1}, {18,32,20,1}}, {{14,31,20,1}, {14,32,20,1}}, {{20,31,20,1}, {20,32,20,1}}, {{18,33,18,1}}, {{17,30,17,0}, {17,31,17,0}}, {{14,33,18,1}}, {{20,33,18,1}}}},
        {"AU",{{{17,30,17,1}, {17,31,17,1}}, {{18,29,15,0}}}},
        {"UA",{{{17,30,17,0}, {17,31,17,0}}, {{18,29,15,1}}}},        
        {"GG",{{{14,31,20,0}, {14,32,20,0}}, {{17,33,18,0}}, {{17,31,20,0},{17,32,20,0}}, {{14,33,18,0}}}},
        {"GC",{{{17,29,17,0}, {17,30,17,0}}, {{15,33,18,1}}, {{14,31,20,1},{14,32,20,1}}}},
        {"CG",{{{17,29,17,1}, {17,30,17,1}}, {{15,33,18,0}}, {{14,31,20,0},{14,32,20,0}}}},
        {"GU",{{{17,31,20,1}, {17,32,20,1}}, {{17,29,15,0}}}},
        {"UG",{{{17,31,20,0}, {17,32,20,0}}, {{17,29,15,1}}}},
        {"CC",{{{14,29,17,0}, {14,30,17,0}}, {{15,29,17,0},{15,30,17,0}}}},
        {"UC",{{{14,29,17,0}, {14,30,17,0}}, {{17,29,17,0},{17,30,17,0}}, {{15,29,15,1}}, {{14,29,15,1}}}},
        {"UC",{{{14,29,17,1}, {14,30,17,1}}, {{17,29,17,1},{17,30,17,1}}, {{15,29,15,0}}, {{14,29,15,0}}}},
        {"UU",{{{14,29,15,0}}, {{17,29,15,0}}}}
    };
    

    base1 = superimposedstr.residues[res1].restype;
	base2 = superimposedstr.residues[res2].restype;
	atomnum1 = superimposedstr.residues[res1].atoms.size();
	atomnum2 = superimposedstr.residues[res2].atoms.size();

    vector<vector<vector<int>>> hbondinfo = BpHbondInfo[base1+base2];

    for(k=0;k<hbondinfo.size();k++)
    {   
        if(hbondinfo[k].size() > 1)
        {
            if(hbondinfo[k][0][3] == 0){
                AccBase = base1[0]; DonBase = base2[0];

                flag = checkhbond2(&superimposedstr.residues[res2].atoms[hbondinfo[k][0][1]].coord, &superimposedstr.residues[res2].atoms[hbondinfo[k][1][1]].coord, &superimposedstr.residues[res2].atoms[hbondinfo[k][0][2]].coord, &superimposedstr.residues[res1].atoms[hbondinfo[k][0][0]].coord);
								
			    if(!flag && bptype==1) continue;

                Hbond2hydobj hbondobjone = Hbond2hydobj::create(&superimposedstr.residues[res2].atoms[hbondinfo[k][0][1]].coord, &superimposedstr.residues[res2].atoms[hbondinfo[k][1][1]].coord, &superimposedstr.residues[res2].atoms[hbondinfo[k][0][2]].coord, &superimposedstr.residues[res1].atoms[hbondinfo[k][0][0]].coord, &gradient[res2][hbondinfo[k][0][1]], &gradient[res2][hbondinfo[k][1][1]], &gradient[res2][hbondinfo[k][0][2]], &gradient[res1][hbondinfo[k][0][0]], hbondparam.at(AccBase).at(DonBase).at(hbondinfo[k][0][0]).at(hbondinfo[k][0][2]).adlength_mean, hbondparam.at(AccBase).at(DonBase).at(hbondinfo[k][0][0]).at(hbondinfo[k][0][2]).ahlength_mean, hbondparam.at(AccBase).at(DonBase).at(hbondinfo[k][0][0]).at(hbondinfo[k][0][2]).dhaangle_mean , K);
                hbondbpobj.push_back(hbondobjone);
            }else{
                AccBase = base2[0]; DonBase = base1[0];

                flag = checkhbond2(&superimposedstr.residues[res1].atoms[hbondinfo[k][0][1]].coord, &superimposedstr.residues[res1].atoms[hbondinfo[k][1][1]].coord, &superimposedstr.residues[res1].atoms[hbondinfo[k][0][2]].coord, &superimposedstr.residues[res2].atoms[hbondinfo[k][0][0]].coord);
								
			    if(!flag && bptype==1) continue;

                Hbond2hydobj hbondobjone = Hbond2hydobj::create(&superimposedstr.residues[res1].atoms[hbondinfo[k][0][1]].coord, &superimposedstr.residues[res1].atoms[hbondinfo[k][1][1]].coord, &superimposedstr.residues[res1].atoms[hbondinfo[k][0][2]].coord, &superimposedstr.residues[res2].atoms[hbondinfo[k][0][0]].coord, &gradient[res1][hbondinfo[k][0][1]], &gradient[res1][hbondinfo[k][1][1]], &gradient[res1][hbondinfo[k][0][2]], &gradient[res2][hbondinfo[k][0][0]], hbondparam.at(AccBase).at(DonBase).at(hbondinfo[k][0][0]).at(hbondinfo[k][0][2]).adlength_mean, hbondparam.at(AccBase).at(DonBase).at(hbondinfo[k][0][0]).at(hbondinfo[k][0][2]).ahlength_mean, hbondparam.at(AccBase).at(DonBase).at(hbondinfo[k][0][0]).at(hbondinfo[k][0][2]).dhaangle_mean, K);
                hbondbpobj.push_back(hbondobjone);
            }
        }else{
            if(hbondinfo[k][0][3] == 0){
                AccBase = base1[0]; DonBase = base2[0];

                flag = checkhbond(&superimposedstr.residues[res2].atoms[hbondinfo[k][0][1]].coord, &superimposedstr.residues[res2].atoms[hbondinfo[k][0][2]].coord, &superimposedstr.residues[res1].atoms[hbondinfo[k][0][0]].coord);
				
			    if(!flag && bptype==1) continue;

                Hbondobj hbondobjone = Hbondobj::create(&superimposedstr.residues[res2].atoms[hbondinfo[k][0][1]].coord, &superimposedstr.residues[res2].atoms[hbondinfo[k][0][2]].coord, &superimposedstr.residues[res1].atoms[hbondinfo[k][0][0]].coord, &gradient[res2][hbondinfo[k][0][1]], &gradient[res2][hbondinfo[k][0][2]], &gradient[res1][hbondinfo[k][0][0]], hbondparam.at(AccBase).at(DonBase).at(hbondinfo[k][0][0]).at(hbondinfo[k][0][2]).adlength_mean, hbondparam.at(AccBase).at(DonBase).at(hbondinfo[k][0][0]).at(hbondinfo[k][0][2]).ahlength_mean, hbondparam.at(AccBase).at(DonBase).at(hbondinfo[k][0][0]).at(hbondinfo[k][0][2]).dhaangle_mean, K);
                hbondobj.push_back(hbondobjone);
            }else{
                AccBase = base2[0]; DonBase = base1[0];

                flag = checkhbond(&superimposedstr.residues[res1].atoms[hbondinfo[k][0][1]].coord, &superimposedstr.residues[res1].atoms[hbondinfo[k][0][2]].coord, &superimposedstr.residues[res2].atoms[hbondinfo[k][0][0]].coord);
				
			    if(!flag && bptype==1) continue;

                Hbondobj hbondobjone = Hbondobj::create(&superimposedstr.residues[res1].atoms[hbondinfo[k][0][1]].coord, &superimposedstr.residues[res1].atoms[hbondinfo[k][0][2]].coord, &superimposedstr.residues[res2].atoms[hbondinfo[k][0][0]].coord, &gradient[res1][hbondinfo[k][0][1]], &gradient[res1][hbondinfo[k][0][2]], &gradient[res2][hbondinfo[k][0][0]], hbondparam.at(AccBase).at(DonBase).at(hbondinfo[k][0][0]).at(hbondinfo[k][0][2]).adlength_mean, hbondparam.at(AccBase).at(DonBase).at(hbondinfo[k][0][0]).at(hbondinfo[k][0][2]).ahlength_mean, hbondparam.at(AccBase).at(DonBase).at(hbondinfo[k][0][0]).at(hbondinfo[k][0][2]).dhaangle_mean, K);
                hbondobj.push_back(hbondobjone);
            }
        }
    }

    /*
    for(k=0;k<atomnum1;k++){
        for(l=0;l<atomnum2;l++){
            if((HBondHydrogenAtoms[base1][k] == 1) && (HBondAcceptorAtoms[base2][l] == 1) && (molConnect[k][l] == 4)){
                for(m=0;m<atomnum1;m++){
                    if(intraatomConnect1[k][m] && (k < 12 || m < 12)){
                        const Point3s &ptH = superimposedstr.residues[res1].atoms[k].coord;
                        const Point3s &ptAcc = superimposedstr.residues[res2].atoms[l].coord;
                        const Point3s &ptDon = superimposedstr.residues[res1].atoms[m].coord;

                        energy=calcenergyhbond(ptH, ptDon, ptAcc);

                        if(energy < HBONDThreshold){
                            Hbondobj hbondobjone = Hbondobj::create(&superimposedstr.residues[res1].atoms[k].coord, &superimposedstr.residues[res1].atoms[m].coord, &superimposedstr.residues[res2].atoms[l].coord, &gradient[res1][k], &gradient[res1][m], &gradient[res2][l], K);
			                hbondobj.push_back(hbondobjone);
                        }
                    }
                }
            }
        }
    }

    for(k=0;k<atomnum2;k++){
        for(l=0;l<atomnum1;l++){
            if((HBondHydrogenAtoms[base2][k] == 1) && (HBondAcceptorAtoms[base1][l] == 1) && (molConnect[l][k] == 4)){
                for(m=0;m<atomnum2;m++){
                    if(intraatomConnect2[k][m] && (k < 12 || m < 12)){
                        const Point3s &ptH = superimposedstr.residues[res2].atoms[k].coord;
                        const Point3s &ptAcc = superimposedstr.residues[res1].atoms[l].coord;
                        const Point3s &ptDon = superimposedstr.residues[res2].atoms[m].coord;

                        energy=calcenergyhbond(ptH, ptDon, ptAcc);

                        if(energy < HBONDThreshold){
                            Hbondobj hbondobjone = Hbondobj::create(&superimposedstr.residues[res2].atoms[k].coord, &superimposedstr.residues[res2].atoms[m].coord, &superimposedstr.residues[res1].atoms[l].coord, &gradient[res2][k], &gradient[res2][m], &gradient[res1][l], K);
			                hbondobj.push_back(hbondobjone);
                        }
                    }
                }
            }
        }
    }
    */
}

void getmolConnect(FragmentStructure &superimposedstr, vector<vector<Point3s>> &gradient, const vector<vector<pair<int, double>>> &stackmap, const map<string, map<int, vector<double>>> &atomfeatures, const map<string, vector<double>> &clashparam, const map<char, map<char, map<int, map<int, HydroBondParameters>>>> &hbondparam,
    const map<string, vector<vector<int>>> &molConnectone, const map<string, vector<vector<int>>> &molConnecttwo, const map<string, vector<vector<int>>> &molConnectthree, 
    const map<string, vector<vector<int>>> &molAdjConnectone, const map<string, vector<vector<int>>> &molAdjConnecttwo, const map<string, vector<vector<int>>> &molAdjConnectthree, 
    const map<string, vector<vector<int>>> &molAdjConnectoneidx, const map<string, vector<vector<int>>> &molAdjConnecttwoidx, const map<string, vector<vector<int>>> &molAdjConnectthreeidx, const vector<vector<int>> &bindmap, const vector<vector<double>> &cbpmap, const vector<vector<double>> &noncbpmap, vector<bool> &moveresidues,
    list<TElectrobj> &telecteobj, list<LennardsJobj> &ljobj, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, list<Clashobj> &clashobj, map<string, double> enecoeff)
{
    int i, j, k;
    string base1, base2, base;
    int atomnum1, atomnum2, atomnum;
    vector<string> nts={"A","C","G","U"};
    map<string, int> atomnumnts={{"A", 33}, {"C", 31}, {"G", 34}, {"U", 30}};
    vector<vector<int>> molConnect;
    map<string, vector<vector<bool>>> intraatomconnect;

    for(i=0;i<nts.size();i++){
        base=nts[i];
        atomnum = atomnumnts[base];
        intraatomconnect[base].resize(atomnum, vector<bool>(atomnum, false));
        for(k=0;k<molConnectone.at(base).size();k++){
            intraatomconnect[base][molConnectone.at(base)[k][0]][molConnectone.at(base)[k][1]] = true;
            intraatomconnect[base][molConnectone.at(base)[k][1]][molConnectone.at(base)[k][0]] = true;
        }
    }

    for(i=0;i<superimposedstr.residues.size();i++){
        base1 = superimposedstr.residues[i].restype;
        atomnum1 = superimposedstr.residues[i].atoms.size();
        for(j=i;j<superimposedstr.residues.size();j++){
            base2 = superimposedstr.residues[j].restype;
            atomnum2 = superimposedstr.residues[j].atoms.size();

            if(!moveresidues[i] && !moveresidues[j]) continue;

            molConnect.clear();
            molConnect.resize(atomnum1, vector<int>(atomnum2, 4));

            if(i == j){
                for(k=0;k<molConnectthree.at(base1).size();k++){
                    molConnect[molConnectthree.at(base1)[k][0]][molConnectthree.at(base1)[k][3]] = 3;
                    molConnect[molConnectthree.at(base1)[k][3]][molConnectthree.at(base1)[k][0]] = 3;
                }
                for(k=0;k<molConnecttwo.at(base1).size();k++){
                    molConnect[molConnecttwo.at(base1)[k][0]][molConnecttwo.at(base1)[k][2]] = 2;
                    molConnect[molConnecttwo.at(base1)[k][2]][molConnecttwo.at(base1)[k][0]] = 2;
                }
                for(k=0;k<molConnectone.at(base1).size();k++){
                    molConnect[molConnectone.at(base1)[k][0]][molConnectone.at(base1)[k][1]] = 1;
                    molConnect[molConnectone.at(base1)[k][1]][molConnectone.at(base1)[k][0]] = 1;
                }
                for(k=0;k<atomnum1;k++){
					molConnect[k][k] = 0;
				}
            }else if(bindmap[i][j] == 1){
                //cout<<i << "-" << j<<endl;
                for(k=0;k<molAdjConnectthree.at(base1).size();k++){
                    //cout<<full_atoms[base1][molAdjConnectthree.at(base1)[k][0]] << "-" << full_atoms[base2][molAdjConnectthree.at(base1)[k][3]] << "-" << 3<<endl;
					molConnect[molAdjConnectthree.at(base1)[k][0]][molAdjConnectthree.at(base1)[k][3]] = 3;
				}
                for(k=0;k<molAdjConnecttwo.at(base1).size();k++){
					molConnect[molAdjConnecttwo.at(base1)[k][0]][molAdjConnecttwo.at(base1)[k][2]] = 2;
				}
                for(k=0;k<molAdjConnectone.at(base1).size();k++){
					molConnect[molAdjConnectone.at(base1)[k][0]][molAdjConnectone.at(base1)[k][1]] = 1;
				}
            }

            if(enecoeff["electros"] > 0) getelectrospairs(superimposedstr, gradient, atomfeatures, i, j, molConnect, telecteobj);
            if(enecoeff["LJ"] > 0) getljpairs(superimposedstr, gradient, atomfeatures, i, j, molConnect, ljobj);
            if(enecoeff["HBond"] > 0){
                if(i != j){
                    //if(bindmap[i][j] == 2){
                    if(cbpmap[i][j] > 0){
                        gethbondpairsbp(superimposedstr, gradient, i, j, molConnect, intraatomconnect[base1], intraatomconnect[base2], hbondparam, hbondobj, hbondbpobj, 0, cbpmap[i][j]*enecoeff["HBond"]);
                    }

                    if(noncbpmap[i][j] > 0){
                        gethbondpairsbp(superimposedstr, gradient, i, j, molConnect, intraatomconnect[base1], intraatomconnect[base2], hbondparam, hbondobj, hbondbpobj, 1, noncbpmap[i][j]*enecoeff["HBond"]);
                    }
                    
                    //else{
                        //gethbondpairs(superimposedstr, gradient, i, j, molConnect, intraatomconnect[base1], intraatomconnect[base2], hbondobj, enecoeff["HBond"]);
                    //}
                }
            }
            if(enecoeff["Clash"] > 0){
                getclashpairs(superimposedstr, gradient, atomfeatures, i, j, molConnect, enecoeff["Clash"], clashobj);
                //getclashpairs(superimposedstr, gradient, clashparam, i, j, molConnect, enecoeff["Clash"], clashobj);
            }
        }
    }
}

void getmolConnectbp(FragmentStructure &superimposedstr, vector<vector<Point3s>> &gradient, int bp_len, const vector<vector<int>> &stackmap, const map<string, map<int, vector<double>>> &atomfeatures, const map<string, vector<double>> &clashparam, const map<char, map<char, map<int, map<int, HydroBondParameters>>>> &hbondparam,
    const map<string, vector<vector<int>>> &molConnectone, const map<string, vector<vector<int>>> &molConnecttwo, const map<string, vector<vector<int>>> &molConnectthree, 
    const map<string, vector<vector<int>>> &molAdjConnectone, const map<string, vector<vector<int>>> &molAdjConnecttwo, const map<string, vector<vector<int>>> &molAdjConnectthree, 
    const map<string, vector<vector<int>>> &molAdjConnectoneidx, const map<string, vector<vector<int>>> &molAdjConnecttwoidx, const map<string, vector<vector<int>>> &molAdjConnectthreeidx, vector<vector<int>> &bindmap, const vector<vector<double>> &cbpmap, const vector<vector<double>> &noncbpmap, vector<bool> &moveresidues,
    list<TElectrobj> &telecteobj, list<LennardsJobj> &ljobj, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, list<Clashobj> &clashobj, map<string, double> enecoeff)
{
    int i, j, k;
    string base1, base2, base;
    int atomnum1, atomnum2, atomnum;
    vector<string> nts={"A","C","G","U"};
    map<string, int> atomnumnts={{"A", 33}, {"C", 31}, {"G", 34}, {"U", 30}};
    vector<vector<int>> molConnect;
    map<string, vector<vector<bool>>> intraatomconnect;

    for(i=0;i<nts.size();i++){
        base=nts[i];
        atomnum = atomnumnts[base];
        intraatomconnect[base].resize(atomnum, vector<bool>(atomnum, false));
        for(k=0;k<molConnectone.at(base).size();k++){
            intraatomconnect[base][molConnectone.at(base)[k][0]][molConnectone.at(base)[k][1]] = true;
            intraatomconnect[base][molConnectone.at(base)[k][1]][molConnectone.at(base)[k][0]] = true;
        }
    }

    for(i=0;i<superimposedstr.residues.size();i++){
        base1 = superimposedstr.residues[i].restype;
        atomnum1 = superimposedstr.residues[i].atoms.size();
        for(j=i;j<superimposedstr.residues.size();j++){
            base2 = superimposedstr.residues[j].restype;
            atomnum2 = superimposedstr.residues[j].atoms.size();

            if(!moveresidues[i] && !moveresidues[j]) continue;

            molConnect.clear();
            molConnect.resize(atomnum1, vector<int>(atomnum2, 4));

            if(i == j){
                for(k=0;k<molConnectthree.at(base1).size();k++){
                    molConnect[molConnectthree.at(base1)[k][0]][molConnectthree.at(base1)[k][3]] = 3;
                    molConnect[molConnectthree.at(base1)[k][3]][molConnectthree.at(base1)[k][0]] = 3;
                }
                for(k=0;k<molConnecttwo.at(base1).size();k++){
                    molConnect[molConnecttwo.at(base1)[k][0]][molConnecttwo.at(base1)[k][2]] = 2;
                    molConnect[molConnecttwo.at(base1)[k][2]][molConnecttwo.at(base1)[k][0]] = 2;
                }
                for(k=0;k<molConnectone.at(base1).size();k++){
                    molConnect[molConnectone.at(base1)[k][0]][molConnectone.at(base1)[k][1]] = 1;
                    molConnect[molConnectone.at(base1)[k][1]][molConnectone.at(base1)[k][0]] = 1;
                }
                for(k=0;k<atomnum1;k++){
					molConnect[k][k] = 0;
				}
            }else if(j == i+1 && j!= bp_len){
                for(k=0;k<molAdjConnectthree.at(base1).size();k++){
                    molConnect[molAdjConnectthree.at(base1)[k][0]][molAdjConnectthree.at(base1)[k][3]] = 3;
				}
                for(k=0;k<molAdjConnecttwo.at(base1).size();k++){
                    molConnect[molAdjConnecttwo.at(base1)[k][0]][molAdjConnecttwo.at(base1)[k][2]] = 2;
				}
                for(k=0;k<molAdjConnectone.at(base1).size();k++){
                    molConnect[molAdjConnectone.at(base1)[k][0]][molAdjConnectone.at(base1)[k][1]] = 1;
				}
            }

            if(enecoeff["electros"] > 0) getelectrospairs(superimposedstr, gradient, atomfeatures, i, j, molConnect, telecteobj);
            if(enecoeff["LJ"] > 0) getljpairs(superimposedstr, gradient, atomfeatures, i, j, molConnect, ljobj);
            if(enecoeff["HBond"] > 0){
                //if(i != j){
                    //gethbondpairs(superimposedstr, gradient, i, j, molConnect, intraatomconnect[base1], intraatomconnect[base2], hbondobj);
                //}

                if(cbpmap[i][j] > 0){
                    gethbondpairsbp(superimposedstr, gradient, i, j, molConnect, intraatomconnect[base1], intraatomconnect[base2], hbondparam, hbondobj, hbondbpobj, 0, cbpmap[i][j]*enecoeff["HBond"]);
                }

                if(noncbpmap[i][j] > 0){
                    gethbondpairsbp(superimposedstr, gradient, i, j, molConnect, intraatomconnect[base1], intraatomconnect[base2], hbondparam, hbondobj, hbondbpobj, 1, noncbpmap[i][j]*enecoeff["HBond"]);
                }
                
                //else{
                    //gethbondpairs(superimposedstr, gradient, i, j, molConnect, intraatomconnect[base1], intraatomconnect[base2], hbondobj, enecoeff["HBond"]);
                //}
            }
            if(enecoeff["Clash"] > 0){
                getclashpairs(superimposedstr, gradient, atomfeatures, i, j, molConnect, enecoeff["Clash"], clashobj);
                //getclashpairs(superimposedstr, gradient, clashparam, i, j, molConnect, enecoeff["Clash"], clashobj);
            }
        }
    }
}

void lbfgsupdate(const vector<vector<Point3s>> &grad, const deque<vector<vector<Point3s>>> &s_list, const deque<vector<vector<Point3s>>> &y_list, const deque<double> &rho_list, const vector<bool> &movenodeinfo, vector<vector<Point3s>> &z){
	vector<vector<Point3s>> q=grad;

	size_t m = s_list.size();

	vector<double> alpha;
	alpha.resize(s_list.size());

	for (int i = m - 1; i >= 0; --i) {
		const vector<vector<Point3s>> &s = s_list[i];
        const vector<vector<Point3s>> &y = y_list[i];
        double rho = rho_list[i];

		double dot_s_q = 0.0;
		//for (auto& [key, s_vec] : s) {
        for(int k=0;k<s.size();k++){
			if(movenodeinfo[k]){
				for (size_t j = 0; j < s[k].size(); ++j) {
					dot_s_q += (s[k][j].x*q[k][j].x);
					dot_s_q += (s[k][j].y*q[k][j].y);
					dot_s_q += (s[k][j].z*q[k][j].z);

                    //cout<<s[k][j].x<<"-" <<s[k][j].y<<"-" <<s[k][j].z<<"-" <<q[k][j].x<<"-" <<q[k][j].y<<"-" <<q[k][j].z<<endl;
				}
			}
		}
		double alpha_i = rho * dot_s_q;

		alpha[i] = alpha_i;

		//for (auto& [key, y_vec] : y) {
        for(int k=0;k<y.size();k++){
			if(movenodeinfo[k]){
				for (size_t j = 0; j < y[k].size(); ++j) {
					q[k][j].x = q[k][j].x - alpha_i * y[k][j].x;
					q[k][j].y = q[k][j].y - alpha_i * y[k][j].y;
					q[k][j].z = q[k][j].z - alpha_i * y[k][j].z;
				}
			}
		}
	}

	double gamma_k = 1.0;
	if(m > 0){
		double gamma_num = 0.0;
		double gamma_den = 0.0;
		//for (auto& [key, s_vec] : s_list[m-1]) {
        for(int k=0;k<s_list[m-1].size();k++){
			if(movenodeinfo[k]){
				for (size_t j = 0; j < s_list[m-1][k].size(); ++j) {
					gamma_num += s_list[m-1][k][j].x * y_list[m-1][k][j].x;
					gamma_num += s_list[m-1][k][j].y * y_list[m-1][k][j].y;
					gamma_num += s_list[m-1][k][j].z * y_list[m-1][k][j].z;
					gamma_den += y_list[m-1][k][j].x * y_list[m-1][k][j].x;
					gamma_den += y_list[m-1][k][j].y * y_list[m-1][k][j].y;
					gamma_den += y_list[m-1][k][j].z * y_list[m-1][k][j].z;
				}
			}
		}
		gamma_k = gamma_num / (gamma_den + 1e-7);
	}

	//for (auto& [key, q_vec] : q) 
    for(int i=0;i<q.size();i++)
	{
		if(movenodeinfo[i]){
			for (size_t j = 0; j < q[i].size(); ++j)
			{
				q[i][j].x *= gamma_k;
				q[i][j].y *= gamma_k;
				q[i][j].z *= gamma_k;
			}
		}
	}

	z = q;

	for (size_t i = 0; i < m; ++i) {
		const vector<vector<Point3s>> &s = s_list[i];
    	const vector<vector<Point3s>> &y = y_list[i];
    	double rho = rho_list[i];

		double dot_y_z = 0.0;
    	//for (auto& [key, y_vec] : y) {
        for(int k=0;k<y.size();k++){
			if(movenodeinfo[k]){
				for (size_t j = 0; j < y[k].size(); ++j) {
					dot_y_z += (y[k][j].x * z[k][j].x);
					dot_y_z += (y[k][j].y * z[k][j].y);
					dot_y_z += (y[k][j].z * z[k][j].z);
				}
			}
    	}
        
		double beta_i = rho * dot_y_z;

		//for (auto& [key, s_vec] : s) {
        for(int k=0;k<s.size();k++){
			if(movenodeinfo[k]){
				for (size_t j = 0; j < s[k].size(); ++j) {
                    z[k][j].x = z[k][j].x + s[k][j].x * (alpha[i] - beta_i);
					z[k][j].y = z[k][j].y + s[k][j].y * (alpha[i] - beta_i);
					z[k][j].z = z[k][j].z + s[k][j].z * (alpha[i] - beta_i);
				}
			}
    	}
	}
}

void updatelbfgsfeatures(const FragmentStructure &superimposedstr, const FragmentStructure &superimposedstr_prev, const vector<vector<Point3s>> &grad, const vector<vector<Point3s>> &grad_prev, deque<vector<vector<Point3s>>> &s_list, deque<vector<vector<Point3s>>> &y_list, deque<double> &rho_list, const vector<bool> &movenodeinfo, int m){
    // 新しい s, y, ρ を保存
    double rho_new =0.0;
    vector<vector<Point3s>> s_new, y_new;
    s_new.resize(grad.size());
    y_new.resize(grad.size());
    //for (auto& [key, grad_vec] : gradptv) {
    for(int i=0;i<grad.size();i++){
        if(movenodeinfo[i]){
            s_new[i].resize(grad[i].size());
            y_new[i].resize(grad[i].size());
            for (size_t j = 0; j < grad[i].size(); ++j) {
                //cout<<key << "-" << j << "-" << gradptv[key][j].x << "-" << gradptvprev[key][j].x<<endl;

                s_new[i][j].x = superimposedstr.residues[i].atoms[j].coord.x - superimposedstr_prev.residues[i].atoms[j].coord.x;
                s_new[i][j].y = superimposedstr.residues[i].atoms[j].coord.y - superimposedstr_prev.residues[i].atoms[j].coord.y;
                s_new[i][j].z = superimposedstr.residues[i].atoms[j].coord.z - superimposedstr_prev.residues[i].atoms[j].coord.z;
                y_new[i][j].x = grad[i][j].x - grad_prev[i][j].x;
                y_new[i][j].y = grad[i][j].y - grad_prev[i][j].y;
                y_new[i][j].z = grad[i][j].z - grad_prev[i][j].z;

                //rho_new+=s_new[i][j].x*y_new[i][j].x;
                //rho_new+=s_new[i][j].y*y_new[i][j].y;
                //rho_new+=s_new[i][j].z*y_new[i][j].z;
                rho_new+=(s_new[i][j].x*y_new[i][j].x+s_new[i][j].y*y_new[i][j].y+s_new[i][j].z*y_new[i][j].z);
                //cout<<i << "-" << j << "-" << s_new[i][j].x << "-" << y_new[i][j].x << "-" << s_new[i][j].y << "-" << y_new[i][j].y <<"-"<< s_new[i][j].z << "-" << y_new[i][j].z<<endl;
            }
        }
    }
   
    rho_new = 1.0/(rho_new+1e-7);
    
    if(m > 0 ){
        if (s_list.size() == m) {
            s_list.pop_front();
            y_list.pop_front();
            rho_list.pop_front();
        }
        s_list.push_back(s_new);
        y_list.push_back(y_new);
        rho_list.push_back(rho_new);
    }
}

/*
void fragmentoptimizationbp(
    FragmentStructure &superimposedstr, const FragmentStructure &inputtedstr, const vector<vector<bool>> &inputtedstrflag, const vector<FragmentStructure> &templatestrall, const vector<vector<vector<bool>>> &templatestrallflag, const map<string, map<string, vector<ModelUnit>>> &ideal_rna, int bplen, 
    double &startene, double &endene, const ExtParams &extparameters, vector<vector<int>> &bindmap, vector<vector<int>> &stackmap, const map<string, map<string, vector<vector<vector<double>>>>> &stackdistmap, vector<bool> &moveresidues, const int &epochs)
    //map<string, vector<vector<double>>> &blengthparam, map<string, vector<vector<double>>> &adjblengthparam,
    //map<string, vector<vector<double>>> &ntangleparam, map<string, vector<vector<int>>> &ntangleparamidx, map<string, vector<vector<double>>> &adjangleparam, map<string, vector<vector<int>>> &adjangleparamidx)
{
    //map<string, double> enecoeff={{"length", 400}, {"angle", 4}, {"dihedral", 0}, {"init", 0.5}, {"improper", 0}, {"electros", 0}, {"LJ", 0}, {"HBond", 0}, {"BP", 1}};
    //map<string, double> enecoeff={{"length", 5000}, {"angle", 8}, {"hydtor", 0}, {"dihedral", 0}, {"init", 10}, {"improper", 0}, {"electros", 1}, {"LJ", 1}, {"HBond", 1}, {"BP", 1}};
    //map<string, double> enecoeff={{"length", 3000}, {"angle", 8}, {"hydtor", 0}, {"dihedral", 0}, {"init", 0}, {"initdist", 1}, {"templatedist", 1}, {"improper", 0}, {"electros", 0}, {"LJ", 0}, {"HBond", 1}, {"BP", 1}, {"Stack", 5}, {"Clash", 20}};
    //map<string, double> enecoeff={{"length", 500}, {"angle",160}, {"hydtor", 0}, {"dihedral", 0}, {"init", 0}, {"initdist", 0.1}, {"templatedist", 0.1}, {"improper", 0}, {"electros", 1}, {"LJ", 1}, {"HBond", 1}, {"BP", 1}, {"Stack", 2}, {"Clash", 4}};
    map<string, double> enecoeff={{"length", 400}, {"angle",140}, {"hydtor", 0}, {"dihedral", 0}, {"init", 0.05}, {"initdist", 0.0}, {"templatedist", 0.0}, {"improper", 0}, {"electros", 1}, {"LJ", 1}, {"HBond", 1}, {"BP", 3.0}, {"Stack", 3}, {"Clash", 5}};

    //FragmentStructure superimposedstr;
    int istart1, iend1, istart2, iend2;

    istart1 = 0;
    iend1 = bplen;
    istart2 = bplen;
    iend2 = bplen*2;

    //optimization
    vector<vector<Point3s>> gradptv;
    prepgrad(superimposedstr, gradptv);

    FragmentStructure superimposedstr_prev;
    vector<vector<Point3s>> gradptv_prev;
    deque<vector<vector<Point3s>>> s_list, y_list;
    deque<double> rho_list;
    vector<vector<Point3s>> z;
    prepgrad(superimposedstr, z);
    int m=15;

    list<InitPosobj> binitobj;
    list<InitPairobj> binitpairobj;
    list<InitPairobj> btemplatepairobj;
    list<BondLengthobj> blengthobj;
	list<BondAngleobj> bangleobj;
    list<BondDihedralobj> dihedralobj;
    list<BasePairobj> bpobj;
    list<BPdistobj> bpdistobj;
    list<TElectrobj> telecteobj;
	list<LennardsJobj> ljobj;
	list<Hbondobj> hbondobj;
    list<Hbond2hydobj> hbondbpobj;
    list<Stackingobj> stackobj;
    list<StackGeoobjMC> stackgeoobj;
    list<Clashobj> clashobj;
    list<SingleNtTemplateobj> nttempobj;

    vector<vector<int>> indices;
    vector<vector<int>> orders;

    if(enecoeff["init"] > 0) setinitobj(superimposedstr, inputtedstr, inputtedstrflag, gradptv, binitobj, enecoeff["init"], moveresidues);
    if(enecoeff["initdist"] > 0) setinitpairobj(superimposedstr, inputtedstr, inputtedstrflag, gradptv, binitpairobj, enecoeff["initdist"], moveresidues);
    if(enecoeff["templatedist"] > 0) settemplatepairobj(superimposedstr, templatestrall, templatestrallflag, gradptv, btemplatepairobj, enecoeff["templatedist"]/templatestrall.size(), moveresidues);

    if(enecoeff["length"] > 0){
        setlengthobj(superimposedstr, gradptv, istart1, iend1, bindmap, extparameters.blengthparam, extparameters.blengthstdval, extparameters.adjblengthparam, extparameters.adjblengthstdval, extparameters.hydparam, extparameters.hydparamidx, extparameters.hydaddlength, extparameters.hydaddlengthidx, enecoeff["length"], blengthobj, moveresidues);
        setlengthobj(superimposedstr, gradptv, istart2, iend2, bindmap, extparameters.blengthparam, extparameters.blengthstdval, extparameters.adjblengthparam, extparameters.adjblengthstdval, extparameters.hydparam, extparameters.hydparamidx, extparameters.hydaddlength, extparameters.hydaddlengthidx, enecoeff["length"], blengthobj, moveresidues);
    }
    if(enecoeff["angle"] > 0){
        setangleobj(superimposedstr, gradptv, istart1, iend1, bindmap, extparameters.ntangleparam, extparameters.ntangleparamidx, extparameters.adjangleparam, extparameters.adjangleparamidx, extparameters.hydaddangle, extparameters.hydaddangleidx, enecoeff["angle"], bangleobj, moveresidues);
        setangleobj(superimposedstr, gradptv, istart2, iend2, bindmap, extparameters.ntangleparam, extparameters.ntangleparamidx, extparameters.adjangleparam, extparameters.adjangleparamidx, extparameters.hydaddangle, extparameters.hydaddangleidx, enecoeff["angle"], bangleobj, moveresidues);
    }

    if(enecoeff["hydtor"] > 0){
        setdihedralobj(superimposedstr, gradptv, istart1, iend1, extparameters.hydtorsionparam, extparameters.btorsionparam, enecoeff["hydtor"], dihedralobj, moveresidues);
        setdihedralobj(superimposedstr, gradptv, istart2, iend2, extparameters.hydtorsionparam, extparameters.btorsionparam, enecoeff["hydtor"], dihedralobj, moveresidues);
    }

    if(enecoeff["BP"] > 0) setbasepairobj(superimposedstr, gradptv, bindmap, enecoeff["BP"], bpobj);
    if(enecoeff["Stack"] > 0){
        setstackingobj(superimposedstr, gradptv, ideal_rna, stackmap, stackdistmap, enecoeff["Stack"], stackobj);
        //setstackingobj(superimposedstr, gradptv, ideal_rna, stackmap, stackdistmap, indices, orders, enecoeff["Stack"], stackobj);
    }

    //if(enecoeff["Stack"] > 0) setstackingobj(superimposedstr, gradptv, ideal_rna, stackmap, stackdistmap, enecoeff["Stack"], stackobj);

    if(enecoeff["electros"] > 0 || enecoeff["LJ"] > 0 || enecoeff["HBond"] > 0 || enecoeff["Clash"] > 0) getmolConnectbp(superimposedstr, gradptv, bplen, stackmap, extparameters.atomfeatures, extparameters.clashparam, extparameters.bphbond, 
                                                                                              extparameters.molConnectone, extparameters.molConnecttwo, extparameters.molConnectthree, 
                                                                                              extparameters.molAdjConnectone, extparameters.molAdjConnecttwo, extparameters.molAdjConnectthree, 
                                                                                              extparameters.molAdjConnectoneidx, extparameters.molAdjConnecttwoidx, extparameters.molAdjConnectthreeidx, bindmap, moveresidues,
                                                                                              telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, enecoeff);

    int estop_counter, stop_count;
    double averagenorm, energy, minenergy;
    //int epochs=10000;
    stop_count = 50;
    minenergy = 100000000;
    estop_counter=0;  

    startene = getcovalentenergystem(blengthobj, bangleobj, dihedralobj, bpobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);

    for (int epoch = 0; epoch < epochs; ++epoch){
        initgrad(superimposedstr, gradptv);
        getgradient(blengthobj, bangleobj, dihedralobj, bpobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);
        //averagenorm = calcaveragenorm(gradptv);

        if(epoch > 1){
            updatelbfgsfeatures(superimposedstr, superimposedstr_prev, gradptv, gradptv_prev, s_list, y_list, rho_list, moveresidues, m); //lbfgs
            lbfgsupdate(gradptv, s_list, y_list, rho_list, moveresidues, z); //lbfgs
        }

        superimposedstr_prev=superimposedstr; //lbfgs
        gradptv_prev=gradptv; //lbfgs

        //energy = goldenelectionsearch(superimposedstr, gradptv, blengthobj, bangleobj, bpobj, binitobj, telecteobj, ljobj, hbondobj, moveresidues);
        if(epoch > 1){
            averagenorm = calcaveragenorm(z);
            //energy = goldenelectionsearch(superimposedstr, z, blengthobj, bangleobj, bpobj, binitobj, telecteobj, ljobj, hbondobj, moveresidues); //lbfgs
            
            double alpha0 = 1.0;
            double rho = 0.5;
            double beta = 0.1;
            energy = linesearch(superimposedstr, z, blengthobj, bangleobj, dihedralobj, bpobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj, moveresidues, alpha0, rho, beta);
        }
        else{
            averagenorm = calcaveragenorm(gradptv);
            energy = goldenelectionsearch(superimposedstr, gradptv, blengthobj, bangleobj, dihedralobj, bpobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj, moveresidues);
        }

        //cout<< "epoch "<<epoch+1 << ": "<< energy<<endl;

        if(energy < minenergy){
            minenergy = energy;
            estop_counter = 0;
        }else{
            estop_counter += 1;
        }

        if(stop_count == estop_counter) break;
    }

    endene = getcovalentenergystem(blengthobj, bangleobj, dihedralobj, bpobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);
}

void fragmentoptimizationloop(FragmentStructure &superimposedstr, const FragmentStructure &inputtedstr, const vector<vector<bool>> &inputtedstrflag, const vector<FragmentStructure> &templatestrall, const vector<vector<vector<bool>>> &templatestrallflag, const map<string, map<string, vector<ModelUnit>>> &ideal_rna, int looplen, 
    double &startene, double &endene, const ExtParams &extparameters, vector<vector<int>> &bindmap, vector<vector<int>> &stackmap, const map<string, map<string, vector<vector<vector<double>>>>> &stackdistmap, vector<bool> &moveresidues, const int &epochs)
    //, map<string, vector<vector<double>>> &blengthparam, map<string, vector<vector<double>>> &adjblengthparam,
    //map<string, vector<vector<double>>> &ntangleparam, map<string, vector<vector<int>>> &ntangleparamidx, map<string, vector<vector<double>>> &adjangleparam, map<string, vector<vector<int>>> &adjangleparamidx)
{
    //map<string, double> enecoeff={{"length", 400}, {"angle", 4}, {"dihedral", 0}, {"init", 0.5}, {"improper", 0}, {"electros", 0}, {"LJ", 0}, {"HBond", 0}, {"BP", 0}};
    //map<string, double> enecoeff={{"length", 500}, {"angle", 4}, {"dihedral", 0}, {"init", 100}, {"improper", 0}, {"electros", 1}, {"LJ", 1}, {"HBond", 1}, {"BP", 1}};
    //map<string, double> enecoeff={{"length", 5000}, {"angle", 8}, {"hydtor", 0}, {"dihedral", 0}, {"init", 10}, {"improper", 0}, {"electros", 1}, {"LJ", 1}, {"HBond", 1}, {"BP", 0}};
    map<string, double> enecoeff={{"length", 3000}, {"angle", 8}, {"hydtor", 0}, {"dihedral", 0}, {"init", 0}, {"initdist", 1}, {"templatedist", 1}, {"improper", 0}, {"electros", 0}, {"LJ", 0}, {"HBond", 1}, {"BP", 3.0}, {"Stack", 5}, {"Clash", 20}};

    int istart, iend;

    istart=0;
    iend=superimposedstr.residues.size();
    
    //optimization
    vector<vector<Point3s>> gradptv;
    prepgrad(superimposedstr, gradptv);

    FragmentStructure superimposedstr_prev;
    vector<vector<Point3s>> gradptv_prev;
    deque<vector<vector<Point3s>>> s_list, y_list;
    deque<double> rho_list;
    vector<vector<Point3s>> z;
    prepgrad(superimposedstr, z);
    int m=15;

    list<InitPosobj> binitobj;
    list<InitPairobj> binitpairobj;
    list<InitPairobj> btemplatepairobj;
    list<BondLengthobj> blengthobj;
	list<BondAngleobj> bangleobj;
    list<BondDihedralobj> dihedralobj;
    list<TElectrobj> telecteobj;
	list<LennardsJobj> ljobj;
	list<Hbondobj> hbondobj;
    list<BasePairobj> bpobj;
    list<BPdistobj> bpdistobj;
    list<Hbond2hydobj> hbondbpobj;
    list<Stackingobj> stackobj;
    list<StackGeoobjMC> stackgeoobj;
    list<Clashobj> clashobj;
    list<SingleNtTemplateobj> nttempobj;

    vector<vector<int>> indices;
    vector<vector<int>> orders;

    if(enecoeff["init"] > 0) setinitobj(superimposedstr, inputtedstr, inputtedstrflag, gradptv, binitobj, enecoeff["init"], moveresidues);
    if(enecoeff["initdist"] > 0) setinitpairobj(superimposedstr, inputtedstr, inputtedstrflag, gradptv, binitpairobj, enecoeff["initdist"], moveresidues);
    if(enecoeff["templatedist"] > 0) settemplatepairobj(superimposedstr, templatestrall, templatestrallflag, gradptv, btemplatepairobj, enecoeff["templatedist"]/templatestrall.size(), moveresidues);

    if(enecoeff["length"] > 0)
    {
        setlengthobj(superimposedstr, gradptv, istart, iend, bindmap, extparameters.blengthparam, extparameters.blengthstdval, extparameters.adjblengthparam, extparameters.adjblengthstdval, extparameters.hydparam, extparameters.hydparamidx, extparameters.hydaddlength, extparameters.hydaddlengthidx, enecoeff["length"], blengthobj, moveresidues);
    }
    if(enecoeff["angle"] > 0)
    {
        setangleobj(superimposedstr, gradptv, istart, iend, bindmap, extparameters.ntangleparam, extparameters.ntangleparamidx, extparameters.adjangleparam, extparameters.adjangleparamidx, extparameters.hydaddangle, extparameters.hydaddangleidx, enecoeff["angle"], bangleobj, moveresidues);
    }

    if(enecoeff["hydtor"] > 0){
        setdihedralobj(superimposedstr, gradptv, istart, iend, extparameters.hydtorsionparam, extparameters.btorsionparam, enecoeff["hydtor"], dihedralobj, moveresidues);
    }
    
    //if(enecoeff["BP"] > 0) setbasepairobj(superimposedstr, gradptv, bindmap, bpobj);
    if(enecoeff["Stack"] > 0){
        setstackingobj(superimposedstr, gradptv, ideal_rna, stackmap, stackdistmap, enecoeff["Stack"], stackobj);
        //setstackingobj(superimposedstr, gradptv, ideal_rna, stackmap, stackdistmap, indices, orders, enecoeff["Stack"], stackobj);
    }

    //if(enecoeff["Stack"] > 0) setstackingobj(superimposedstr, gradptv, ideal_rna, stackmap, stackdistmap, enecoeff["Stack"], stackobj);

    if(enecoeff["electros"] > 0 || enecoeff["LJ"] > 0 || enecoeff["HBond"] > 0 || enecoeff["Clash"] > 0) getmolConnect(superimposedstr, gradptv, stackmap, extparameters.atomfeatures, extparameters.clashparam, extparameters.bphbond,
                                                                                              extparameters.molConnectone, extparameters.molConnecttwo, extparameters.molConnectthree, 
                                                                                              extparameters.molAdjConnectone, extparameters.molAdjConnecttwo, extparameters.molAdjConnectthree, 
                                                                                              extparameters.molAdjConnectoneidx, extparameters.molAdjConnecttwoidx, extparameters.molAdjConnectthreeidx, bindmap, moveresidues,
                                                                                              telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, enecoeff);
                                                                                        
    startene = getcovalentenergystemloop(blengthobj, bangleobj, dihedralobj, bpobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);
    int estop_counter, stop_count;
    double averagenorm, energy, minenergy;
    //int epochs=50000;
    stop_count = 50;
    minenergy = 100000000;
    estop_counter=0;
    for (int epoch = 0; epoch < epochs; ++epoch){
        initgrad(superimposedstr, gradptv);
        getgradientloop(blengthobj, bangleobj, dihedralobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);

        if(epoch > 1){
            updatelbfgsfeatures(superimposedstr, superimposedstr_prev, gradptv, gradptv_prev, s_list, y_list, rho_list, moveresidues, m); //lbfgs
            lbfgsupdate(gradptv, s_list, y_list, rho_list, moveresidues, z); //lbfgs
        }

        superimposedstr_prev=superimposedstr; //lbfgs
        gradptv_prev=gradptv; //lbfgs
        
        if(epoch > 1){
            averagenorm = calcaveragenorm(z);
            //energy = goldenelectionsearch(superimposedstr, z, blengthobj, bangleobj, bpobj, binitobj, telecteobj, ljobj, hbondobj, moveresidues); //lbfgs
            double alpha0 = 3.0;
            double rho = 0.5;
            double beta = 0.0000001;
            energy = linesearch(superimposedstr, z, blengthobj, bangleobj, dihedralobj, bpobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj, moveresidues, alpha0, rho, beta);
        }else{
            averagenorm = calcaveragenorm(gradptv);
            energy = goldenelectionsearchloop(superimposedstr, gradptv, blengthobj, bangleobj, dihedralobj, bpobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj, moveresidues);
        }
        //cout<< "epoch "<<epoch+1 << ": "<< energy<<endl;

        if(energy < minenergy){
            minenergy = energy;
            estop_counter = 0;
        }else{
            estop_counter += 1;
        }

        if(stop_count == estop_counter) break;
    }

    endene = getcovalentenergystemloop(blengthobj, bangleobj, dihedralobj, bpobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);
}
*/

void fragmentoptimization(
    FragmentStructure &superimposedstr, const FragmentStructure &inputtedstr, const vector<vector<bool>> &inputtedstrflag, const vector<FragmentStructure> &templatestrall, const vector<vector<vector<bool>>> &templatestrallflag, const map<string, map<string, vector<ModelUnit>>> &ideal_rna, 
    double &startene, double &endene, const ExtParams &extparameters, vector<vector<int>> &bindmap, const vector<vector<pair<int, double>>> &stackmap, const vector<vector<double>> &cbpmap, const vector<vector<double>> &noncbpmap, const map<string, vector<vector<vector<double>>>> &bpdistmap, const map<string, map<string, vector<vector<vector<double>>>>> &stackdistmap, vector<bool> &moveresidues, const double maxtime, const int &epochs, Logger &logger)
    //map<string, vector<vector<double>>> &blengthparam, map<string, vector<vector<double>>> &adjblengthparam,
    //map<string, vector<vector<double>>> &ntangleparam, map<string, vector<vector<int>>> &ntangleparamidx, map<string, vector<vector<double>>> &adjangleparam, map<string, vector<vector<int>>> &adjangleparamidx)
{
    
    //map<string, double> enecoeff={{"length", 500}, {"angle",160}, {"hydtor", 0}, {"dihedral", 0}, {"init", 0}, {"initdist", 0.5}, {"templatedist", 0.5}, {"improper", 0}, {"electros", 1.0}, {"LJ", 1.0}, {"HBond", 20.0}, {"BP", 1.0}, {"Stack", 2.0}, {"StackGeo", 5.0}, {"Clash", 5.0}, {"Nt", 5.0}, {"BPdist", 0.0}, {"BasePlane", 5.0}}; 
    //map<string, double> enecoeff={{"length", 400}, {"angle",120}, {"hydtor", 0}, {"dihedral", 0}, {"init", 0}, {"initdist", 0.5}, {"templatedist", 0.5}, {"improper", 0}, {"electros", 1.0}, {"LJ", 1.0}, {"HBond", 20.0}, {"BP", 0.0}, {"Stack", 0.5}, {"StackGeo", 0.0}, {"Clash", 5.0}, {"Nt", 0.0}, {"BPdist", 0.0}, {"BasePlane", 5.0}}; 
    map<string, double> enecoeff={{"length", 400}, {"angle",120}, {"hydtor", 0}, {"dihedral", 0}, {"init", 0}, {"initdist", 0.5}, {"templatedist", 0.5}, {"improper", 0}, {"electros", 1.0}, {"LJ", 1.0}, {"HBond", 10.0}, {"BP", 0.0}, {"Stack", 0.5}, {"StackGeo", 0.0}, {"Clash", 5.0}, {"Nt", 0.5}, {"BPdist", 0.0}, {"BasePlane", 1.0}}; 

    //FragmentStructure superimposedstr;
    int istart, iend;

    istart = 0;
    iend=superimposedstr.residues.size();

    //optimization
    vector<vector<Point3s>> gradptv;
    prepgrad(superimposedstr, gradptv);
    std::ostringstream log_buffer;

    FragmentStructure superimposedstr_prev;
    vector<vector<Point3s>> gradptv_prev;
    deque<vector<vector<Point3s>>> s_list, y_list;
    deque<double> rho_list;
    vector<vector<Point3s>> z;
    prepgrad(superimposedstr, z);
    int m=15;

    vector<vector<int>> indices;
    vector<vector<int>> orders;
    vector<vector<Point3s*>> AtomPtrslistA;
    vector<vector<Point3s*>> AtomgradlistA;
    vector<vector<Point3s*>> AtomPtrslistB;
    vector<vector<Point3s*>> AtomgradlistB;

    list<InitPosobj> binitobj;
    list<InitPairobj> binitpairobj;
    list<InitPairobj> btemplatepairobj;
    list<BondLengthobj> blengthobj;
	list<BondAngleobj> bangleobj;
    list<BondDihedralobj> dihedralobj;
    list<BasePairobj> bpobj;
    list<BasePairobj> baseplaneobj;
    list<BPdistobj> bpdistobj;
    list<TElectrobj> telecteobj;
	list<LennardsJobj> ljobj;
	list<Hbondobj> hbondobj;
    list<Hbond2hydobj> hbondbpobj;
    list<Stackingobj> stackobj;
    list<StackGeoobjMC> stackgeoobj;
    list<Clashobj> clashobj;
    list<SingleNtTemplateobj> nttempobj;

    if(enecoeff["init"] > 0) setinitobj(superimposedstr, inputtedstr, inputtedstrflag, gradptv, binitobj, enecoeff["init"], moveresidues);
    if(enecoeff["initdist"] > 0) setinitpairobj(superimposedstr, inputtedstr, inputtedstrflag, gradptv, binitpairobj, enecoeff["initdist"]/(1+templatestrall.size()), moveresidues);
    if(enecoeff["templatedist"] > 0) settemplatepairobj(superimposedstr, templatestrall, templatestrallflag, gradptv, btemplatepairobj, enecoeff["initdist"]/(1+templatestrall.size()), moveresidues);

    if(enecoeff["length"] > 0){
        setlengthobj(superimposedstr, gradptv, istart, iend, bindmap, extparameters.blengthparam, extparameters.blengthstdval, extparameters.adjblengthparam, extparameters.adjblengthstdval, extparameters.hydparam, extparameters.hydparamidx, extparameters.hydaddlength, extparameters.hydaddlengthidx, enecoeff["length"], blengthobj, moveresidues);
    }
    if(enecoeff["angle"] > 0){
        setangleobj(superimposedstr, gradptv, istart, iend, bindmap, extparameters.ntangleparam, extparameters.ntangleparamidx, extparameters.adjangleparam, extparameters.adjangleparamidx, extparameters.hydaddangle, extparameters.hydaddangleidx, enecoeff["angle"], bangleobj, moveresidues);
    }

    if(enecoeff["hydtor"] > 0){
        setdihedralobj(superimposedstr, gradptv, istart, iend, extparameters.hydtorsionparam, extparameters.btorsionparam, enecoeff["hydtor"], dihedralobj, moveresidues);
    }

    if(enecoeff["BP"] > 0) setbasepairobj(superimposedstr, gradptv, bindmap, cbpmap, noncbpmap, enecoeff["BP"], bpobj);
    
    if(enecoeff["BasePlane"] > 0) setbaseplaneobj(superimposedstr, gradptv, enecoeff["BasePlane"], baseplaneobj);
    
    if(enecoeff["Stack"] > 0) setstackingobj(superimposedstr, gradptv, ideal_rna, stackmap, stackdistmap, enecoeff["Stack"], stackobj);

    if(enecoeff["StackGeo"] > 0) setstackinggeoobj(superimposedstr, gradptv, extparameters.stackparam, stackmap, AtomPtrslistA, AtomgradlistA, AtomPtrslistB, AtomgradlistB, enecoeff["StackGeo"], stackgeoobj);

    if(enecoeff["BPdist"] > 0) setbpobj(superimposedstr, gradptv, ideal_rna, bindmap, bpdistmap, enecoeff["BPdist"], bpdistobj);
    
        //setstackingobj(superimposedstr, gradptv, ideal_rna, stackmap, stackdistmap, indices, orders, enecoeff["Stack"], stackobj);
    if(enecoeff["Nt"] > 0) setnttemplateobj(superimposedstr, gradptv, ideal_rna, indices, enecoeff["Nt"], nttempobj);
    //if(enecoeff["Stack"] > 0) setstackingobj(superimposedstr, gradptv, ideal_rna, stackmap, stackdistmap, enecoeff["Stack"], stackobj);
    
    if(enecoeff["electros"] > 0 || enecoeff["LJ"] > 0 || enecoeff["HBond"] > 0 || enecoeff["Clash"] > 0) getmolConnect(superimposedstr, gradptv, stackmap, extparameters.atomfeatures, extparameters.clashparam, extparameters.bphbond,
                                                                                                extparameters.molConnectone, extparameters.molConnecttwo, extparameters.molConnectthree, 
                                                                                                extparameters.molAdjConnectone, extparameters.molAdjConnecttwo, extparameters.molAdjConnectthree, 
                                                                                                extparameters.molAdjConnectoneidx, extparameters.molAdjConnecttwoidx, extparameters.molAdjConnectthreeidx, bindmap, cbpmap, noncbpmap, moveresidues,
                                                                                                telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, enecoeff);
    
    if(enecoeff["electros"] <= 0) telecteobj.clear();
    if(enecoeff["LJ"] <= 0) ljobj.clear();
    
    if(enecoeff["HBond"] <= 0){
        hbondobj.clear();
        hbondbpobj.clear();
    }

    if(enecoeff["Clash"] <= 0){
        clashobj.clear();
    }

    int estop_counter, stop_count;
    double averagenorm, energy, minenergy;
    //int epochs=10000;
    stop_count = 30;
    minenergy = 100000000;
    estop_counter=0;  
    auto start_time = std::chrono::high_resolution_clock::now();

    startene = getcovalentenergystem(blengthobj, bangleobj, dihedralobj, baseplaneobj, bpobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);

    for (int epoch = 0; epoch < epochs; ++epoch){
        initgrad(superimposedstr, gradptv);
        getgradient(blengthobj, bangleobj, dihedralobj, baseplaneobj, bpobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);
        //averagenorm = calcaveragenorm(gradptv);

        if(epoch > 1){
            updatelbfgsfeatures(superimposedstr, superimposedstr_prev, gradptv, gradptv_prev, s_list, y_list, rho_list, moveresidues, m); //lbfgs
            lbfgsupdate(gradptv, s_list, y_list, rho_list, moveresidues, z); //lbfgs
        }

        superimposedstr_prev=superimposedstr; //lbfgs
        gradptv_prev=gradptv; //lbfgs

        //energy = goldenelectionsearch(superimposedstr, gradptv, blengthobj, bangleobj, bpobj, binitobj, telecteobj, ljobj, hbondobj, moveresidues);
        if(epoch > 1){
            averagenorm = calcaveragenorm(z, 10000);
            //energy = goldenelectionsearch(superimposedstr, z, blengthobj, bangleobj, bpobj, binitobj, telecteobj, ljobj, hbondobj, moveresidues); //lbfgs
            
            double alpha0 = 1.0;
            double rho = 0.5;
            double beta = 0.1;
            energy = linesearch(superimposedstr, z, blengthobj, bangleobj, dihedralobj, baseplaneobj, bpobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj, moveresidues, alpha0, rho, beta);
            //energy = goldenelectionsearch(superimposedstr, z, blengthobj, bangleobj, dihedralobj, bpobj, stackobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, moveresidues);
            
            //averagenorm = calcaveragenorm(gradptv);
            //energy = goldenelectionsearch(superimposedstr, gradptv, blengthobj, bangleobj, dihedralobj, bpobj, stackobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, moveresidues);
        
        }
        else{
            averagenorm = calcaveragenorm(gradptv, 10000);
            energy = goldenelectionsearch(superimposedstr, gradptv, blengthobj, bangleobj, dihedralobj, baseplaneobj, bpobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj, moveresidues);
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_time = end_time - start_time;
        //logger.log("epoch "+to_string(epoch+1)+", energy: "+to_string(energy)+" Time: "+to_string(elapsed_time.count()));

        log_buffer << "epoch (L-BFGS): " << epoch+1 << ", "
                   << "energy: " << energy << ", "
                   << "Time: " << elapsed_time.count() << "\n";

        if (epoch % 100 == 0 && !log_buffer.str().empty()) {
            logger.log(log_buffer.str());
            log_buffer.str("");
            log_buffer.clear();
        }

        if(energy < minenergy){
            minenergy = energy;
            estop_counter = 0;
        }else{
            estop_counter += 1;
        }

        if(stop_count == estop_counter) break;
        if(elapsed_time.count() > maxtime) break;
    }

    /*
    estop_counter = 0;
    for (int epoch = 0; epoch < epochs; ++epoch){
        initgrad(superimposedstr, gradptv);
        getgradient(blengthobj, bangleobj, dihedralobj, bpobj, stackobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj);
        
        averagenorm = calcaveragenorm(gradptv);
        energy = goldenelectionsearch(superimposedstr, gradptv, blengthobj, bangleobj, dihedralobj, bpobj, stackobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, moveresidues);
        
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_time = end_time - start_time;

        log_buffer << "epoch (Sweepy descent): " << epoch+1 << ", "
                   << "energy: " << energy << ", "
                   << "Time: " << elapsed_time.count() << "\n";

        if (epoch % 100 == 0 && !log_buffer.str().empty()) {
            logger.log(log_buffer.str());
            log_buffer.str("");
            log_buffer.clear();
        }

        if(energy < minenergy){
            minenergy = energy;
            estop_counter = 0;
        }else{
            estop_counter += 1;
        }

        if(stop_count == estop_counter) break;
        if(elapsed_time.count() > maxtime) break;
    }
    */

    if (log_buffer.tellp() > 0) { 
        logger.log(log_buffer.str());
    }

    endene = getcovalentenergystem(blengthobj, bangleobj, dihedralobj, baseplaneobj, bpobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);
}


/*
void fragmentoptimizationbpopt(
    FragmentStructure &superimposedstr, const FragmentStructure &inputtedstr, const vector<vector<bool>> &inputtedstrflag, const vector<FragmentStructure> &templatestrall, const vector<vector<vector<bool>>> &templatestrallflag, 
    double &startene, double &endene, const ExtParams &extparameters, vector<vector<int>> &bindmap, vector<vector<int>> &stackmap, vector<bool> &moveresidues, const int &epochs)
    //map<string, vector<vector<double>>> &blengthparam, map<string, vector<vector<double>>> &adjblengthparam,
    //map<string, vector<vector<double>>> &ntangleparam, map<string, vector<vector<int>>> &ntangleparamidx, map<string, vector<vector<double>>> &adjangleparam, map<string, vector<vector<int>>> &adjangleparamidx)
{
    //map<string, double> enecoeff={{"length", 400}, {"angle", 4}, {"dihedral", 0}, {"init", 0.5}, {"improper", 0}, {"electros", 0}, {"LJ", 0}, {"HBond", 0}, {"BP", 1}};
    map<string, double> enecoeff={{"length", 3000}, {"angle", 8}, {"dihedral", 0}, {"init", 0}, {"initdist", 1}, {"templatedist", 1}, {"improper", 0}, {"electros", 1}, {"LJ", 1}, {"HBond", 1}, {"BP", 3.0}, {"Clash", 0}};

    //FragmentStructure superimposedstr;

    int ntnum = superimposedstr.residues.size();
    
    //optimization
    vector<vector<Point3s>> gradptv;
    prepgrad(superimposedstr, gradptv);

    list<InitPosobj> binitobj;
    list<InitPairobj> binitpairobj;
    list<InitPairobj> btemplatepairobj;
    list<BondLengthobj> blengthobj;
	list<BondAngleobj> bangleobj;
    list<BondDihedralobj> dihedralobj;
    list<BasePairobj> bpobj;
    list<BPdistobj> bpdistobj;
    list<TElectrobj> telecteobj;
	list<LennardsJobj> ljobj;
	list<Hbondobj> hbondobj;
    list<Hbond2hydobj> hbondbpobj;
    list<Stackingobj> stackobj;
    list<StackGeoobjMC> stackgeoobj;
    list<Clashobj> clashobj;
    list<SingleNtTemplateobj> nttempobj;

    if(enecoeff["init"] > 0) setinitobj(superimposedstr, inputtedstr, inputtedstrflag, gradptv, binitobj, enecoeff["init"], moveresidues);
    if(enecoeff["initdist"] > 0) setinitpairobj(superimposedstr, inputtedstr, inputtedstrflag, gradptv, binitpairobj, enecoeff["initdist"], moveresidues);
    if(enecoeff["templatedist"] > 0) settemplatepairobj(superimposedstr, templatestrall, templatestrallflag, gradptv, btemplatepairobj, enecoeff["templatedist"]/templatestrall.size(), moveresidues);

    if(enecoeff["length"] > 0){
        setlengthobj(superimposedstr, gradptv, 0, ntnum, bindmap, extparameters.blengthparam, extparameters.blengthstdval, extparameters.adjblengthparam, extparameters.adjblengthstdval, extparameters.hydparam, extparameters.hydparamidx, extparameters.hydaddlength, extparameters.hydaddlengthidx, enecoeff["length"], blengthobj, moveresidues);
        //setlengthobj(superimposedstr, gradptv, 0, ntnum, bindmap, extparameters.blengthparam, extparameters.adjblengthparam, enecoeff["length"], blengthobj, moveresidues);
    }
    if(enecoeff["angle"] > 0){
        setangleobj(superimposedstr, gradptv, 0, ntnum, bindmap, extparameters.ntangleparam, extparameters.ntangleparamidx, extparameters.adjangleparam, extparameters.adjangleparamidx, extparameters.hydaddangle, extparameters.hydaddangleidx, enecoeff["angle"], bangleobj, moveresidues);
        //setangleobj(superimposedstr, gradptv, 0, ntnum, bindmap, extparameters.ntangleparam, extparameters.ntangleparamidx, extparameters.adjangleparam, extparameters.adjangleparamidx, enecoeff["angle"], bangleobj, moveresidues);
    }

    if(enecoeff["hydtor"] > 0){
        setdihedralobj(superimposedstr, gradptv, 0, ntnum, extparameters.hydtorsionparam, extparameters.btorsionparam, enecoeff["hydtor"], dihedralobj, moveresidues);
    }

    if(enecoeff["BP"] > 0) setbasepairobj(superimposedstr, gradptv, bindmap, enecoeff["BP"], bpobj);

    if(enecoeff["electros"] > 0 || enecoeff["LJ"] > 0 || enecoeff["HBond"] > 0 || enecoeff["Clash"] > 0) getmolConnect(superimposedstr, gradptv, stackmap, extparameters.atomfeatures, extparameters.clashparam, extparameters.bphbond,
                                                                                              extparameters.molConnectone, extparameters.molConnecttwo, extparameters.molConnectthree, 
                                                                                              extparameters.molAdjConnectone, extparameters.molAdjConnecttwo, extparameters.molAdjConnectthree, 
                                                                                              extparameters.molAdjConnectoneidx, extparameters.molAdjConnecttwoidx, extparameters.molAdjConnectthreeidx, bindmap, moveresidues,
                                                                                              telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, enecoeff);

    startene = getcovalentenergystem(blengthobj, bangleobj, dihedralobj, bpobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);
    
    int estop_counter, stop_count;
    double averagenorm, energy, minenergy;
    //int epochs=10000;
    stop_count = 50;
    minenergy = 100000000;
    estop_counter=0;
    for (int epoch = 0; epoch < epochs; ++epoch){
        initgrad(superimposedstr, gradptv);
        getgradient(blengthobj, bangleobj, dihedralobj, bpobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);

        averagenorm = calcaveragenorm(gradptv);
        energy = goldenelectionsearch(superimposedstr, gradptv, blengthobj, bangleobj, dihedralobj, bpobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj, moveresidues);
        //cout<< "epoch "<<epoch+1 << ": "<< energy<<endl;

        if(energy < minenergy){
            minenergy = energy;
            estop_counter = 0;
        }else{
            estop_counter += 1;
        }

        if(stop_count == estop_counter) break;
    }

    endene = getcovalentenergystem(blengthobj, bangleobj, dihedralobj, bpobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);
}

void fragmentoptimizationallopt(
    FragmentStructure &superimposedstr, const FragmentStructure &inputtedstr, const vector<vector<bool>> &inputtedstrflag, const vector<FragmentStructure> &templatestrall, const vector<vector<vector<bool>>> &templatestrallflag, 
    double &startene, double &endene, const ExtParams &extparameters, vector<vector<int>> &bindmap, vector<vector<int>> &stackmap, vector<bool> &moveresidues, const int &epochs)
    //map<string, vector<vector<double>>> &blengthparam, map<string, vector<vector<double>>> &adjblengthparam,
    //map<string, vector<vector<double>>> &ntangleparam, map<string, vector<vector<int>>> &ntangleparamidx, map<string, vector<vector<double>>> &adjangleparam, map<string, vector<vector<int>>> &adjangleparamidx)
{
    //map<string, double> enecoeff={{"length", 400}, {"angle", 4}, {"dihedral", 0}, {"init", 0.5}, {"improper", 0}, {"electros", 0}, {"LJ", 0}, {"HBond", 0}, {"BP", 1}};
    map<string, double> enecoeff={{"length", 3000}, {"angle", 8}, {"dihedral", 0}, {"init", 0}, {"initdist", 1}, {"templatedist", 1}, {"improper", 0}, {"electros", 1}, {"LJ", 1}, {"HBond", 1}, {"BP", 3.0}, {"Clash", 0}};

    //FragmentStructure superimposedstr;

    int ntnum = superimposedstr.residues.size();
    
    //optimization
    vector<vector<Point3s>> gradptv;
    prepgrad(superimposedstr, gradptv);

    list<InitPosobj> binitobj;
    list<InitPairobj> binitpairobj;
    list<InitPairobj> btemplatepairobj;
    list<BondLengthobj> blengthobj;
	list<BondAngleobj> bangleobj;
    list<BondDihedralobj> dihedralobj;
    list<BasePairobj> bpobj;
    list<BPdistobj> bpdistobj;
    list<TElectrobj> telecteobj;
	list<LennardsJobj> ljobj;
	list<Hbondobj> hbondobj;
    list<Hbond2hydobj> hbondbpobj;
    list<Stackingobj> stackobj;
    list<StackGeoobjMC> stackgeoobj;
    list<Clashobj> clashobj;
    list<SingleNtTemplateobj> nttempobj;

    if(enecoeff["init"] > 0) setinitobj(superimposedstr, inputtedstr, inputtedstrflag, gradptv, binitobj, enecoeff["init"], moveresidues);
    if(enecoeff["initdist"] > 0) setinitpairobj(superimposedstr, inputtedstr, inputtedstrflag, gradptv, binitpairobj, enecoeff["initdist"], moveresidues);
    if(enecoeff["templatedist"] > 0) settemplatepairobj(superimposedstr, templatestrall, templatestrallflag, gradptv, btemplatepairobj, enecoeff["templatedist"]/templatestrall.size(), moveresidues);

    if(enecoeff["length"] > 0){
        setlengthobj(superimposedstr, gradptv, 0, ntnum, bindmap, extparameters.blengthparam, extparameters.blengthstdval, extparameters.adjblengthparam, extparameters.adjblengthstdval, extparameters.hydparam, extparameters.hydparamidx, extparameters.hydaddlength, extparameters.hydaddlengthidx, enecoeff["length"], blengthobj, moveresidues);
        //setlengthobj(superimposedstr, gradptv, 0, ntnum, bindmap, extparameters.blengthparam, extparameters.adjblengthparam, enecoeff["length"], blengthobj, moveresidues);
    }
    if(enecoeff["angle"] > 0){
        setangleobj(superimposedstr, gradptv, 0, ntnum, bindmap, extparameters.ntangleparam, extparameters.ntangleparamidx, extparameters.adjangleparam, extparameters.adjangleparamidx, extparameters.hydaddangle, extparameters.hydaddangleidx, enecoeff["angle"], bangleobj, moveresidues);
        //setangleobj(superimposedstr, gradptv, 0, ntnum, bindmap, extparameters.ntangleparam, extparameters.ntangleparamidx, extparameters.adjangleparam, extparameters.adjangleparamidx, enecoeff["angle"], bangleobj, moveresidues);
    }

    if(enecoeff["hydtor"] > 0){
        setdihedralobj(superimposedstr, gradptv, 0, ntnum, extparameters.hydtorsionparam, extparameters.btorsionparam, enecoeff["hydtor"], dihedralobj, moveresidues);
    }

    if(enecoeff["BP"] > 0) setbasepairobj(superimposedstr, gradptv, bindmap, enecoeff["BP"], bpobj);

    if(enecoeff["electros"] > 0 || enecoeff["LJ"] > 0 || enecoeff["HBond"] > 0 || enecoeff["Clash"] > 0) getmolConnect(superimposedstr, gradptv, stackmap, extparameters.atomfeatures, extparameters.clashparam, extparameters.bphbond,
                                                                                              extparameters.molConnectone, extparameters.molConnecttwo, extparameters.molConnectthree, 
                                                                                              extparameters.molAdjConnectone, extparameters.molAdjConnecttwo, extparameters.molAdjConnectthree, 
                                                                                              extparameters.molAdjConnectoneidx, extparameters.molAdjConnecttwoidx, extparameters.molAdjConnectthreeidx, bindmap, moveresidues,
                                                                                              telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, enecoeff);
                                                                                        
    startene = getcovalentenergystem(blengthobj, bangleobj, dihedralobj, bpobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);

    int estop_counter, stop_count;
    double averagenorm, energy, minenergy;
    //int epochs=10000;
    stop_count = 50;
    minenergy = 100000000;
    estop_counter=0;
    for (int epoch = 0; epoch < epochs; ++epoch){
        initgrad(superimposedstr, gradptv);
        getgradient(blengthobj, bangleobj, dihedralobj, bpobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);
        averagenorm = calcaveragenorm(gradptv);

        energy = goldenelectionsearch(superimposedstr, gradptv, blengthobj, bangleobj, dihedralobj, bpobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj, moveresidues);
        //cout<< "epoch "<<epoch+1 << ": "<< energy<<endl;

        if(energy < minenergy){
            minenergy = energy;
            estop_counter = 0;
        }else{
            estop_counter += 1;
        }

        if(stop_count == estop_counter) break;
    }

    endene = getcovalentenergystem(blengthobj, bangleobj, dihedralobj, bpobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);
}
*/
