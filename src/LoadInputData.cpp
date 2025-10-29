
#include "LoadInputData.hpp"

bool check_targetpdb(ModelUnit pdb_entry, vector<string> neccatoms){
    if(pdb_entry.chains.size() != 1){
        cerr<<"Error: Input structure must be a single chain structure.";
        exit(1);
    }

    //vector<string> neccatoms = {"  P ", " C4'", " C3'"};
    int neccatoms_num = neccatoms.size();

    std::ostringstream neccatoms_oss;
    for (size_t i = 0; i < neccatoms.size(); ++i) {
        if (i != 0) neccatoms_oss << ", ";
        neccatoms_oss << neccatoms[i];
    }
    string neccatoms_str = neccatoms_oss.str();

    for(int i=0;i<pdb_entry.chains[0].residues.size();i++){
        int count=0;
        for(int j=0;j<pdb_entry.chains[0].residues[i].atoms.size();j++){
            string atomname=pdb_entry.chains[0].residues[i].atoms[j].name;
            auto it = find(neccatoms.begin(), neccatoms.end(), atomname);
            if(it != neccatoms.end()){
                count+=1;
            }
        }

        if(count >= neccatoms_num){
            cerr<<"Error: " << neccatoms_str << " atoms must be included in each residue (hint: residue "<< pdb_entry.chains[0].residues[i].resi<<").";
            exit(1);
        }
    }
    return true;
}

void gettargetfilesfiles(string tdir, vector<string> &targetfiles, string ext)
{
	//vector<string>refname;
	try {
        for (const auto& entry : fs::directory_iterator(tdir)) {
            if (fs::is_regular_file(entry) && entry.path().extension() == ext) {
				targetfiles.push_back(entry.path());
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Error:: template file cannot be loaded." << e.what() << std::endl;
    }
}

void get_multitemplate_info(const ModelUnit &pdb_entry, const string &tempstrdir, vector<FragmentStructure> &templatestrall, vector<vector<vector<bool>>> &templatestrallflag, int &atomic_detail, int &allowX, vector<string> neccatoms)
{
    vector<ModelUnit> template_pdb_entry;
    vector<string> templatestrs;

    if (!tempstrdir.empty()) {
        gettargetfilesfiles(tempstrdir, templatestrs,".pdb");
    }else{
        return;
    }

    //templatestrall.resize(templatestrs.size());
    //templatestrallflag.resize(templatestrs.size());

    int numbb=pdb_entry.chains[0].residues.size();

    for(int i=0;i<templatestrs.size();i++){
        cout<<"Loading " << templatestrs[i] << "..."<<endl;
        map<int, map<int, int>> pdbidx2seqidx;
        map<string, int> chainstr2chainidx;
        ModelUnit pdb_tmp=read_pdb_structure(templatestrs[i].c_str(), chainstr2chainidx, pdbidx2seqidx, atomic_detail,allowX);
        bool flag = check_targetpdb(pdb_tmp, neccatoms);

        if(!flag){
            cerr << "Warning: The strcuture does not include backbone atoms. [" << templatestrs[i] << " is ignored]" << endl;
            continue;
        }

        template_pdb_entry.push_back(pdb_tmp);

        FragmentStructure templatestr;
        vector<vector<bool>> templatestrflag;
        
        templatestr.residues.resize(numbb);
        templatestrflag.resize(numbb);

        if(template_pdb_entry[i].chains.size() != 1){
            cerr<<"Error: template structure must be single chain structures. (" << templatestrs[i] << ")";
            exit(1);
        }

        for(int k=0;k<template_pdb_entry[i].chains[0].residues.size();k++){
            if(template_pdb_entry[i].chains[0].residues[k].resi > pdb_entry.chains[0].residues[numbb-1].resi) continue;
            if(template_pdb_entry[i].chains[0].residues[k].resi < pdb_entry.chains[0].residues[0].resi) continue;

            int residx=template_pdb_entry[i].chains[0].residues[k].resi - pdb_entry.chains[0].residues[0].resi;
            templatestr.residues[residx].atoms.resize(template_target_atoms.size());
            templatestr.residues[residx].restype = template_pdb_entry[i].chains[0].residues[k].resn;
            templatestrflag[residx].resize(template_target_atoms.size(), false);

            for(int m=0;m<template_target_atoms.size();m++){
                templatestr.residues[residx].atoms[m].coord.x = 0.0;
                templatestr.residues[residx].atoms[m].coord.y = 0.0;
                templatestr.residues[residx].atoms[m].coord.z = 0.0;
            }

            for(int l=0;l<template_pdb_entry[i].chains[0].residues[k].atoms.size();l++){
                for(int m=0;m<template_target_atoms.size();m++){
                    if(template_target_atoms[m] == template_pdb_entry[i].chains[0].residues[k].atoms[l].name){
                        templatestr.residues[residx].atoms[m].coord.x=template_pdb_entry[i].chains[0].residues[k].atoms[l].xyz[0];
                        templatestr.residues[residx].atoms[m].coord.y=template_pdb_entry[i].chains[0].residues[k].atoms[l].xyz[1];
                        templatestr.residues[residx].atoms[m].coord.z=template_pdb_entry[i].chains[0].residues[k].atoms[l].xyz[2];
                        templatestrflag[residx][m]=true;
                    }
                }
            }
        }

        templatestrall.push_back(templatestr);
        templatestrallflag.push_back(templatestrflag);
    }

    if(templatestrall.size() > 0) cout<<"Template structures (" << templatestrall.size() << " numbers) were loaded."<<endl;
}

void get_multibp_info(const ModelUnit &pdb_entry, const string &bpdir, const string &seq, map<string, int> &chainstr2chainidx, map<int, map<int, int>> &pdbidx2seqidx, vector<vector<pair<double,vector<size_t>>>> &bps_list, vector<vector<pair<double,vector<size_t>>>> &noncbps_list, vector<vector<vector<double>>> &bps_scoremat, vector<vector<vector<bool>>> &bps_flagmat, vector<vector<int>> &bps_ranges)
{
     int i,j;
    vector<string> bpsfiles;
    double sc;
    int chainA, chainB, baseA, baseB;
    char basetype1, basetype2;
    int seq_len;
    double minscore, maxscore;

    if (!bpdir.empty()) {
        gettargetfilesfiles(bpdir, bpsfiles,".bpseq");
    }else{
        return;
    }

    seq_len = seq.length();

    for(i=0;i<bpsfiles.size();i++){
        //cout<<"[[" << i <<"]]"<<endl;
        //cout<<bpsfiles[i]<<endl;
        int bpseqst, bpseqend;
        vector<pair<double,vector<size_t>>> bp_vec;
        vector<pair<double,vector<size_t>>> noncbp_vec;

        vector<vector<double>> bps_scoremat_tmp(seq_len, vector<double>(seq_len, 0.0));
        vector<vector<bool>> bps_flagmat_tmp(seq_len, vector<bool>(seq_len, 0.0));
        
        loadBP(bpsfiles[i], seq, chainstr2chainidx, pdbidx2seqidx, bp_vec, noncbp_vec, bpseqst, bpseqend);
        bps_list.push_back(bp_vec);
        noncbps_list.push_back(noncbp_vec);

        for(int k=bpseqst;k<bpseqend+1;k++){
            for(int l=bpseqst;l<bpseqend+1;l++){
                bps_flagmat_tmp[k][l]=true;
                bps_flagmat_tmp[l][k]=true;
            }
        }

        for(j=0;j<bp_vec.size();j++){
            sc = bp_vec[j].first;
            chainA = bp_vec[j].second[0];
            baseA = bp_vec[j].second[1];
            chainA = bp_vec[j].second[2];
            baseB = bp_vec[j].second[3];

            if(baseA < bpseqst || baseA > bpseqend || baseB < bpseqst || baseB > bpseqend){
                cout<<"Warning: Nucleotide index in base-pair is out of range. The pair is ignored. (hint: " << bpsfiles[i]<<")"<<endl;
                continue;
            }

            bps_scoremat_tmp[baseA][baseB]+=sc;
            bps_scoremat_tmp[baseB][baseA]+=sc;
        }

        for(j=0;j<noncbp_vec.size();j++){
            sc = noncbp_vec[j].first;
            chainA = noncbp_vec[j].second[0];
            baseA = noncbp_vec[j].second[1];
            chainA = noncbp_vec[j].second[2];
            baseB = noncbp_vec[j].second[3];

            if(baseA < bpseqst || baseA > bpseqend || baseB < bpseqst || baseB > bpseqend){
                cout<<"Warning: Nucleotide index in base-pair is out of range. The pair is ignored. (hint: " << bpsfiles[i]<<")"<<endl;
                continue;
            }

            bps_scoremat_tmp[baseA][baseB]+=sc;
            bps_scoremat_tmp[baseB][baseA]+=sc;
        }

        bps_scoremat.push_back(bps_scoremat_tmp);
        bps_flagmat.push_back(bps_flagmat_tmp);
        bps_ranges.push_back({bpseqst, bpseqend});

        //for(j=0;j<seq_len;j++){
            //for(int k=0;k<seq_len;k++){
                //if(bps_flagmat_tmp[j][k]) cout<<j << "-" << k << "-" << bps_scoremat_tmp[j][k]<<endl;
            //}
        //}
    }
    
}

void get_multistack_info(const ModelUnit &pdb_entry, const string &stackdir, const string &seq, map<string, int> &chainstr2chainidx, map<int, map<int, int>> &pdbidx2seqidx, vector<vector<pair<double,vector<size_t>>>> &stack_list, vector<vector<vector<double>>> &stack_scoremat, vector<vector<vector<bool>>> &stack_flagmat, vector<vector<int>> &stack_ranges)
{
    int i,j;
    vector<string> stackfiles;
    double sc;
    int chainA, chainB, baseA, baseB;
    char basetype1, basetype2;
    int seq_len;
    double minscore, maxscore;
    bool flag;

    if (!stackdir.empty()) {
        gettargetfilesfiles(stackdir, stackfiles,".stack");
    }else{
        return;
    }

    seq_len = seq.length();
    for(i=0;i<stackfiles.size();i++){
        //cout<<"[[" << i <<"]]"<<endl;
        //cout<<stackfiles[i]<<endl;
        int stackseqst, stackseqend;
        vector<pair<double,vector<size_t>>> stack_vec;
        
        flag=loadST(stackfiles[i], seq, chainstr2chainidx, pdbidx2seqidx, stack_vec, stackseqst, stackseqend);
        if(!flag) continue;
        stack_list.push_back(stack_vec);

        vector<vector<double>> stack_scoremat_tmp(seq_len, vector<double>(seq_len, 0.0));
        vector<vector<bool>> stack_flagmat_tmp(seq_len, vector<bool>(seq_len, 0.0));

        for(int k=stackseqst;k<stackseqend+1;k++){
            for(int l=stackseqst;l<stackseqend+1;l++){
                stack_flagmat_tmp[k][l]=true;
                stack_flagmat_tmp[l][k]=true;
            }
        }

        for(j=0;j<stack_vec.size();j++){
            sc = stack_vec[j].first;
            chainA = stack_vec[j].second[0];
            baseA = stack_vec[j].second[1];
            chainA = stack_vec[j].second[2];
            baseB = stack_vec[j].second[3];

            if(baseA < stackseqst || baseA > stackseqend || baseB < stackseqst || baseB > stackseqend){
                cout<<"Warning: Nucleotide index in base-pair is out of range. The pair is ignored. (hint: " << stackfiles[i]<<")"<<endl;
                continue;
            }

            stack_scoremat_tmp[baseA][baseB]+=sc;
            stack_scoremat_tmp[baseB][baseA]+=sc;
        }

        stack_scoremat.push_back(stack_scoremat_tmp);
        stack_flagmat.push_back(stack_flagmat_tmp);
        stack_ranges.push_back({stackseqst, stackseqend});

        //for(j=0;j<seq_len;j++){
            //for(int k=0;k<seq_len;k++){
                //if(stack_flagmat_tmp[j][k]) cout<<j << "-" << k << "-" << stack_scoremat_tmp[j][k]<<endl;
            //}
        //}
    }
}

void getmultidistmap(const ModelUnit &pdb_entry, const string &distdir, const string &seq, vector<vector<vector<double>>> &distmaps, vector<vector<vector<bool>>> &distmapflags, vector<vector<int>> &distmapranges)
{
    int i,j;
    vector<string> distfiles;
    double sc;
    int chainA, chainB, baseA, baseB;
    char basetype1, basetype2;
    int seq_len;
    double minscore, maxscore;

    if (!distdir.empty()) {
        gettargetfilesfiles(distdir, distfiles,".dist");
    }else{
        return;
    }

    seq_len = seq.length();

    for(int i=0;i<distfiles.size();i++){
        int distseqst, distseqend;
        vector<int> distmaprange(2);
        vector<vector<double>> distmap;
        vector<vector<bool>> distmapflag;

        bool flag = loadDistmap(distfiles[i], seq, distmap, distmapflag, distmaprange[0], distmaprange[1]);
        if(flag){
            distmaps.push_back(distmap);
            distmapflags.push_back(distmapflag);
            distmapranges.push_back(distmaprange);
        }
    }
}


void getmulticontactmap(const ModelUnit &pdb_entry, const string &contactdir, const string &seq, vector<vector<vector<double>>> &contactmaps, vector<vector<vector<bool>>> &contactmapflags, vector<vector<int>> &contactmapranges)
{
    int i,j;
    vector<string> contactfiles;
    double sc;
    int chainA, chainB, baseA, baseB;
    char basetype1, basetype2;
    int seq_len;
    double minscore, maxscore;

    if (!contactdir.empty()) {
        gettargetfilesfiles(contactdir, contactfiles,".cont");
    }else{
        return;
    }

    seq_len = seq.length();

    for(int i=0;i<contactfiles.size();i++){
        int contactseqst, contactseqend;
        vector<int> contactmaprange(2);
        vector<vector<double>> contactmap;
        vector<vector<bool>> contactmapflag;

        bool flag = loadContactmap(contactfiles[i], seq, contactmap, contactmapflag, contactmaprange[0], contactmaprange[1]);
        if(flag){
            contactmaps.push_back(contactmap);
            contactmapflags.push_back(contactmapflag);
            contactmapranges.push_back(contactmaprange);
        }
    }
}

/*
void get_multibp_info(const ModelUnit &pdb_entry, const string &bpdir, const string &seq, vector<pair<double,vector<size_t>>> &bps_list, vector<vector<double>> &bps_scoremat, double threshold)
{
    int i,j;
    vector<string> bpsfiles;
    double sc;
    int chainA, chainB, baseA, baseB;
    char basetype1, basetype2;
    int seq_len;
    double minscore, maxscore;

    if (!bpdir.empty()) {
        gettargetfilesfiles(bpdir, bpsfiles,".bpseq");
    }else{
        return;
    }

    seq_len = seq.length();
    bps_scoremat.resize(seq_len, vector<double>(seq_len, 0.0));  
    vector<vector<double>> count(seq_len, vector<double>(seq_len, 0.0));

    for(i=0;i<bpsfiles.size();i++){
        int bpseqst, bpseqend;
        vector<pair<double,vector<size_t>>> bp_vec;
        
        loadBP(bpsfiles[i], seq, bp_vec, bpseqst, bpseqend);

        for(int k=bpseqst;k<bpseqend+1;k++){
            for(int l=bpseqst;l<bpseqend+1;l++){
                count[k][l]+=1;
                count[l][k]+=1;
            }
        }

        for(j=0;j<bp_vec.size();j++){
            sc = bp_vec[j].first;
            chainA = bp_vec[j].second[0];
            baseA = bp_vec[j].second[1];
            chainA = bp_vec[j].second[2];
            baseB = bp_vec[j].second[3];

            if(baseA < bpseqst || baseA > bpseqend || baseB < bpseqst || baseB > bpseqend){
                cout<<"Warning: Nucleotide index in base-pair is out of range. The pair is ignored. (hint: " << bpsfiles[i]<<")"<<endl;
                continue;
            }

            bps_scoremat[baseA][baseB]+=sc;
            bps_scoremat[baseB][baseA]+=sc;
        }
    }

    minscore = 100000000;
    maxscore = -1.0;
    for(int k=0;k<seq_len;k++){
        for(int l=k;l<seq_len;l++){
            bps_scoremat[k][l]/=count[k][l];
            bps_scoremat[l][k]=bps_scoremat[k][l];

            if(minscore > bps_scoremat[k][l]) minscore=bps_scoremat[k][l];
            if(maxscore < bps_scoremat[k][l]) maxscore=bps_scoremat[k][l];
        }
    }

    cout<<minscore << "-" << maxscore<<endl;

    for(int k=0;k<seq_len;k++){
        for(int l=k;l<seq_len;l++){
            bps_scoremat[k][l] = (bps_scoremat[k][l] - minscore)/(maxscore-minscore);
            bps_scoremat[l][k] = bps_scoremat[k][l];

            if(bps_scoremat[k][l] > threshold){
                pair<double,vector<size_t>> bp_vec_temp;
                bp_vec_temp.first = bps_scoremat[k][l];
                bp_vec_temp.second.resize(5,0);

                bp_vec_temp.second[0] = 0;
                bp_vec_temp.second[1] = k;
                bp_vec_temp.second[2] = 0;
                bp_vec_temp.second[3] = l;

                basetype1 = seq[bp_vec_temp.second[1]];
                basetype2 = seq[bp_vec_temp.second[3]];

                bool canonicaflag=false;
                if ((basetype1=='A' &&(basetype2=='U' || basetype2=='T')) ||
                    (basetype1=='C' && basetype2=='G')||
                    ((basetype1=='U' || basetype1=='T')&& basetype2=='A')  ||
                    (basetype1=='G' && basetype2=='C'))
                {
                        bp_vec_temp.second[4]=0;
                        canonicaflag=true;
                }
                else if ((basetype1=='G' && basetype2=='U') ||
                         (basetype1=='U' && basetype2=='G'))
                {
                            bp_vec_temp.second[4]=1;
                            canonicaflag=true;
                }

                if(canonicaflag){
                    bps_list.push_back(bp_vec_temp);
                }
            }
        }
    }

    for(int i=0;i<bps_list.size();i++){
        cout<<bps_list[i].first << "-" << bps_list[i].second[0] << "-" << bps_list[i].second[1] << "-" << bps_list[i].second[2] << "-" << bps_list[i].second[3] << "-" << bps_list[i].second[4]<<endl;
    }
}
*/

/*
void get_multistack_info(const ModelUnit &pdb_entry, const string &stackdir, const string &seq, vector<pair<double,vector<size_t>>> &stack_list, vector<vector<double>> &stack_scoremat, double threshold)
{
    int i,j;
    vector<string> stackfiles;
    double sc;
    int chainA, chainB, baseA, baseB;
    char basetype1, basetype2;
    int seq_len;
    double minscore, maxscore;

    if (!stackdir.empty()) {
        gettargetfilesfiles(stackdir, stackfiles,".stack");
    }else{
        return;
    }

    seq_len = seq.length();
    stack_scoremat.resize(seq_len, vector<double>(seq_len, 0.0));  
    vector<vector<double>> count(seq_len, vector<double>(seq_len, 0.0));

    for(i=0;i<stackfiles.size();i++){
        int stackseqst, stackseqend;
        vector<pair<double,vector<size_t>>> stack_vec;
        
        loadST(stackfiles[i], seq, stack_vec, stackseqst, stackseqend);

        for(int k=stackseqst;k<stackseqend+1;k++){
            for(int l=stackseqst;l<stackseqend+1;l++){
                count[k][l]+=1;
                count[l][k]+=1;
            }
        }

        for(j=0;j<stack_vec.size();j++){
            sc = stack_vec[j].first;
            chainA = stack_vec[j].second[0];
            baseA = stack_vec[j].second[1];
            chainA = stack_vec[j].second[2];
            baseB = stack_vec[j].second[3];

            if(baseA < stackseqst || baseA > stackseqend || baseB < stackseqst || baseB > stackseqend){
                cout<<"Warning: Nucleotide index in base-pair is out of range. The pair is ignored. (hint: " << stackfiles[i]<<")"<<endl;
                continue;
            }

            stack_scoremat[baseA][baseB]+=sc;
            stack_scoremat[baseB][baseA]+=sc;
        }
    }

    minscore = 100000000;
    maxscore = -1.0;
    for(int k=0;k<seq_len;k++){
        for(int l=k;l<seq_len;l++){
            stack_scoremat[k][l]/=count[k][l];
            stack_scoremat[l][k]=stack_scoremat[k][l];

            if(minscore > stack_scoremat[k][l]) minscore=stack_scoremat[k][l];
            if(maxscore < stack_scoremat[k][l]) maxscore=stack_scoremat[k][l];
        }
    }

    cout<<minscore << "-" << maxscore<<endl;

    for(int k=0;k<seq_len;k++){
        for(int l=k;l<seq_len;l++){
            stack_scoremat[k][l] = (stack_scoremat[k][l] - minscore)/(maxscore-minscore);
            stack_scoremat[l][k] = stack_scoremat[k][l];

            if(stack_scoremat[k][l] > threshold){
                pair<double,vector<size_t>> stack_vec_temp;
                stack_vec_temp.first = stack_scoremat[k][l];
                stack_vec_temp.second.resize(5,0);

                stack_vec_temp.second[0] = 0;
                stack_vec_temp.second[1] = k;
                stack_vec_temp.second[2] = 0;
                stack_vec_temp.second[3] = l;

                
                if(l == k + 1) stack_vec_temp.second[4]=1;

                stack_list.push_back(stack_vec_temp);
            }
        }
    }

    for(int i=0;i<stack_list.size();i++){
        cout<<stack_list[i].first << "-" << stack_list[i].second[0] << "-" << stack_list[i].second[1] << "-" << stack_list[i].second[2] << "-" << stack_list[i].second[3] << "-" << stack_list[i].second[4]<<endl;
    }
}
*/

























