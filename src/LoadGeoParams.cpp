
#include "LoadGeoParams.hpp"

bool loadlength(string basepath, map<string, vector<vector<double>>> &ntlengthparam, map<string, vector<vector<double>>> &ntlengthstdval, map<string, vector<vector<double>>> &adjlengthparam, map<string, vector<vector<double>>> &adjlengthstdval)
{
    int count;
    string	oneline, field;
    string atom1, atom2;
    int residx1, residx2;
    string base;
    int idx1, idx2;
    double mean, stdval;
    vector<string> nts = {"A","C","G","U"};
    for(int i=0;i<nts.size();i++)
	{
        base = nts[i];
        ifstream inputFile(basepath+"/length_"+base+".txt");
        if(!inputFile.is_open())
		{
			cout<<"Parameter file can not be opened (length-" << base << ")\n"<<endl;
			return false;
		}

        ntlengthparam[base].resize(full_atoms_name[base].size(), vector<double>(full_atoms_name[base].size(), -1.0));
        adjlengthparam[base].resize(full_atoms_name[base].size(), vector<double>(full_atoms_name[base].size(), -1.0));

        ntlengthstdval[base].resize(full_atoms_name[base].size(), vector<double>(full_atoms_name[base].size(), -1.0));
        adjlengthstdval[base].resize(full_atoms_name[base].size(), vector<double>(full_atoms_name[base].size(), -1.0));

        while(getline(inputFile, oneline))
		{
            stringstream ssbase(oneline);
			count=0;
			while (getline(ssbase, field, ','))
			{
                //if(count==0) ssidx=stoi(field); //ssidx
				if(count==0) atom1=field; 
				if(count==1) atom2=field; 
				if(count==2) residx1=stoi(field); 
				if(count==3) residx2=stoi(field); 
				if(count==4) mean=stof(field); //mean
				if(count==5) stdval=stof(field); //stdval
				count+=1;
            }

            if(count == 6){
                if(residx1 == 0 && residx2 == 0){
                    auto it1 = find(full_atoms_name[base].begin(), full_atoms_name[base].end(), atom1);
                    idx1 = -1;
                    if(it1 != full_atoms_name[base].end()){
                        idx1=distance(full_atoms_name[base].begin(), it1);
                    }else{
                        cout<<atom1<<endl;
                        cout<<"Error1"<<endl;
                        exit(1);
                    }

                    auto it2 = find(full_atoms_name[base].begin(), full_atoms_name[base].end(), atom2);
                    idx2 = -1;
                    if(it2 != full_atoms_name[base].end()){
                        idx2=distance(full_atoms_name[base].begin(), it2);
                    }else{
                        cout<<atom2<<endl;
                        cout<<"Error2"<<endl;
                        exit(1);
                    }

                    ntlengthparam[base][idx1][idx2] = mean;
                    ntlengthparam[base][idx2][idx1] = mean;

                    ntlengthstdval[base][idx1][idx2] = stdval;
                    ntlengthstdval[base][idx2][idx1] = stdval;

                }else if(residx1 == -1 && residx2 == 0){
                    auto it1 = find(full_atoms_name[base].begin(), full_atoms_name[base].end(), atom1);
                    idx1 = -1;
                    if(it1 != full_atoms_name[base].end()){
                        idx1=distance(full_atoms_name[base].begin(), it1);
                    }else{
                        cout<<atom1<<endl;
                        cout<<"Error1"<<endl;
                        exit(1);
                    }

                    auto it2 = find(full_atoms_name[base].begin(), full_atoms_name[base].end(), atom2);
                    idx2 = -1;
                    if(it2 != full_atoms_name[base].end()){
                        idx2=distance(full_atoms_name[base].begin(), it2);
                    }else{
                        cout<<atom2<<endl;
                        cout<<"Error2"<<endl;
                        exit(1);
                    }

                    adjlengthparam[base][idx1][idx2] = mean;
                    adjlengthstdval[base][idx1][idx2] = stdval;

                }
            }else{
                cout<<"Error3"<<endl;
                exit(1);
            }
        }
    }

    return true;
}

bool loadangle(string basepath, map<string, vector<vector<double>>> &ntangleparam, map<string, vector<vector<int>>> &ntangleparamidx, map<string, vector<vector<double>>> &adjangleparam, map<string, vector<vector<int>>> &adjangleparamidx)
{
    int count;
    string	oneline, field;
    string atom1, atom2, atom3;
    int residx1, residx2, residx3;
    string base;
    int idx1, idx2, idx3;
    double mean, stdval;
    vector<string> nts = {"A","C","G","U"};
    for(int i=0;i<nts.size();i++)
	{
        base = nts[i];
        ifstream inputFile(basepath+"/angle_"+base+".txt");
        if(!inputFile.is_open())
		{
			cout<<"Parameter file can not be opened (angle-" << base << ")\n"<<endl;
			return false;
		}

        vector<double> ntangleparamone(2);
        vector<int> ntangleparamidxone(3);
        vector<double> adjangleparamone(2);
        vector<int> adjangleparamidxone(6);

        while(getline(inputFile, oneline))
		{
            stringstream ssbase(oneline);
			count=0;
			while (getline(ssbase, field, ','))
			{
                //if(count==0) ssidx=stoi(field); //ssidx
				if(count==0) atom1=field; 
				if(count==1) atom2=field; 
                if(count==2) atom3=field; 
				if(count==3) residx1=stoi(field); 
				if(count==4) residx2=stoi(field); 
                if(count==5) residx3=stoi(field); 
				if(count==6) mean=stof(field); //mean
				if(count==7) stdval=stof(field); //stdval
				count+=1;
            }

            if(count == 8){
                if(residx1 != -1 && residx2 != -1 && residx3 != -1){
                    auto it1 = find(full_atoms_name[base].begin(), full_atoms_name[base].end(), atom1);
                    idx1 = -1;
                    if(it1 != full_atoms_name[base].end()){
                        idx1=distance(full_atoms_name[base].begin(), it1);
                    }else{
                        cout<<atom1<<endl;
                        cout<<"Error1"<<endl;
                        exit(1);
                    }

                    auto it2 = find(full_atoms_name[base].begin(), full_atoms_name[base].end(), atom2);
                    idx2 = -1;
                    if(it2 != full_atoms_name[base].end()){
                        idx2=distance(full_atoms_name[base].begin(), it2);
                    }else{
                        cout<<atom2<<endl;
                        cout<<"Error2"<<endl;
                        exit(1);
                    }

                    auto it3 = find(full_atoms_name[base].begin(), full_atoms_name[base].end(), atom3);
                    idx3 = -1;
                    if(it3 != full_atoms_name[base].end()){
                        idx3=distance(full_atoms_name[base].begin(), it3);
                    }else{
                        cout<<atom3<<endl;
                        cout<<"Error2"<<endl;
                        exit(1);
                    }

                    ntangleparamidxone = {idx1, idx2, idx3};
                    ntangleparamone = {mean, stdval};

                    ntangleparamidx[nts[i]].push_back(ntangleparamidxone);
                    ntangleparam[nts[i]].push_back(ntangleparamone);
                }else if(residx1 == -1 || residx2 == -1 || residx3 == -1){
                    auto it1 = find(full_atoms_name[base].begin(), full_atoms_name[base].end(), atom1);
                    idx1 = -1;
                    if(it1 != full_atoms_name[base].end()){
                        idx1=distance(full_atoms_name[base].begin(), it1);
                    }else{
                        cout<<atom1<<endl;
                        cout<<"Error1"<<endl;
                        exit(1);
                    }

                    auto it2 = find(full_atoms_name[base].begin(), full_atoms_name[base].end(), atom2);
                    idx2 = -1;
                    if(it2 != full_atoms_name[base].end()){
                        idx2=distance(full_atoms_name[base].begin(), it2);
                    }else{
                        cout<<atom2<<endl;
                        cout<<"Error2"<<endl;
                        exit(1);
                    }

                    auto it3 = find(full_atoms_name[base].begin(), full_atoms_name[base].end(), atom3);
                    idx3 = -1;
                    if(it3 != full_atoms_name[base].end()){
                        idx3=distance(full_atoms_name[base].begin(), it3);
                    }else{
                        cout<<atom3<<endl;
                        cout<<"Error2"<<endl;
                        exit(1);
                    }

                    adjangleparamidxone = {residx1, residx2, residx3, idx1, idx2, idx3};
                    adjangleparamone = {mean, stdval};

                    adjangleparamidx[nts[i]].push_back(adjangleparamidxone);
                    adjangleparam[nts[i]].push_back(adjangleparamone);
                }
            }else{
                cout<<"Error3"<<endl;
                exit(1);
            }
        }
    }

    return true;
}

bool loadhydtorsion(string basepath, map<string, vector<pair<vector<int>, vector<double>>>> &hydtorsionparam)
{
    string	oneline, field;
    string atom1, atom2, atom3, atom4;
    string base;
    double dki, idv, dpn, dphase;
    int idx1, idx2, idx3, idx4;
    double mean, stdval;
    vector<string> nts = {"A","C","G","U"};
    for(int i=0;i<nts.size();i++)
	{
        base = nts[i];
        ifstream inputFile(basepath+"/hyd_torsion_"+base+".txt");
        if(!inputFile.is_open())
		{
			cout<<"Parameter file can not be opened (hyd_torsion-" << base << ")\n"<<endl;
			return false;
		}

        bool loadflag=false;
        while(getline(inputFile, oneline))
		{
            stringstream ssbase(oneline);
			while (getline(ssbase, field))
			{
                istringstream iss(field);

                if (iss >> atom1 >> atom2 >> atom3 >> atom4 >> idv >> dki >> dphase >> dpn){
                    loadflag=true;
                };
            }

            pair<vector<int>, vector<double>> hydtorsionparam_tmp;

            if(loadflag){
                hydtorsionparam_tmp.first.resize(4, 0);
                hydtorsionparam_tmp.second.resize(4, 0.0);

                //cout<<atom1 << "-" << atom2 << "-" << atom3 << "-" << atom4<<endl;
                //cout<<idv << "-" << dki << "-" << dphase << "-" << dpn<<endl;

                auto it1 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom1);
                idx1 = -1;
                if(it1 != hyd_full_atoms_name[base].end()){
                    idx1=distance(hyd_full_atoms_name[base].begin(), it1);
                }else{
                    cout<<atom1<<endl;
                    cout<<"Error1"<<endl;
                    exit(1);
                }

                auto it2 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom2);
                idx2 = -1;
                if(it2 != hyd_full_atoms_name[base].end()){
                    idx2=distance(hyd_full_atoms_name[base].begin(), it2);
                }else{
                    cout<<atom2<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                auto it3 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom3);
                idx3 = -1;
                if(it3 != hyd_full_atoms_name[base].end()){
                    idx3=distance(hyd_full_atoms_name[base].begin(), it3);
                }else{
                    cout<<atom3<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                auto it4 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom4);
                idx4 = -1;
                if(it4 != hyd_full_atoms_name[base].end()){
                    idx4=distance(hyd_full_atoms_name[base].begin(), it4);
                }else{
                    cout<<atom4<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                hydtorsionparam_tmp.first[0] = idx1;
                hydtorsionparam_tmp.first[1] = idx2;
                hydtorsionparam_tmp.first[2] = idx3;
                hydtorsionparam_tmp.first[3] = idx4;

                //cout<<hydtorsionparam_tmp.first[0] << "-" << hydtorsionparam_tmp.first[1] << "-" << hydtorsionparam_tmp.first[2] << "-" << hydtorsionparam_tmp.first[3]<<endl;

                hydtorsionparam_tmp.second[0] = idv;
                hydtorsionparam_tmp.second[1] = dki;
                hydtorsionparam_tmp.second[2] = dphase;
                hydtorsionparam_tmp.second[3] = dpn;

                hydtorsionparam[nts[i]].push_back(hydtorsionparam_tmp);
            }else{
                cout<<"Error3"<<endl;
                exit(1);
            }
        }
    }

    return true;
}

bool loadtorsion(string basepath, map<string, vector<pair<vector<int>, vector<double>>>> &torsionparam)
{
    string	oneline, field;
    string atom1, atom2, atom3, atom4;
    string base;
    double dki, idv, dpn, dphase;
    int idx1, idx2, idx3, idx4;
    double mean, stdval;
    vector<string> nts = {"A","C","G","U"};
    for(int i=0;i<nts.size();i++)
	{
        base = nts[i];
        ifstream inputFile(basepath+"/torsion_"+base+".txt");
        if(!inputFile.is_open())
		{
			cout<<"Parameter file can not be opened (torsion-" << base << ")\n"<<endl;
			return false;
		}

        bool loadflag=false;
        while(getline(inputFile, oneline))
		{
            stringstream ssbase(oneline);
			while (getline(ssbase, field))
			{
                istringstream iss(field);

                if (iss >> atom1 >> atom2 >> atom3 >> atom4 >> idv >> dki >> dphase >> dpn){
                    loadflag=true;
                };
            }

            pair<vector<int>, vector<double>> torsionparam_tmp;

            if(loadflag)
            {
                torsionparam_tmp.first.resize(4, 0);
                torsionparam_tmp.second.resize(4, 0.0);

                //cout<<atom1 << "-" << atom2 << "-" << atom3 << "-" << atom4<<endl;
                //cout<<idv << "-" << dki << "-" << dphase << "-" << dpn<<endl;

                auto it1 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom1);
                idx1 = -1;
                if(it1 != hyd_full_atoms_name[base].end()){
                    idx1=distance(hyd_full_atoms_name[base].begin(), it1);
                }else{
                    cout<<atom1<<endl;
                    cout<<"Error1"<<endl;
                    exit(1);
                }

                auto it2 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom2);
                idx2 = -1;
                if(it2 != hyd_full_atoms_name[base].end()){
                    idx2=distance(hyd_full_atoms_name[base].begin(), it2);
                }else{
                    cout<<atom2<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                auto it3 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom3);
                idx3 = -1;
                if(it3 != hyd_full_atoms_name[base].end()){
                    idx3=distance(hyd_full_atoms_name[base].begin(), it3);
                }else{
                    cout<<atom3<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                auto it4 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom4);
                idx4 = -1;
                if(it4 != hyd_full_atoms_name[base].end()){
                    idx4=distance(hyd_full_atoms_name[base].begin(), it4);
                }else{
                    cout<<atom4<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                torsionparam_tmp.first[0] = idx1;
                torsionparam_tmp.first[1] = idx2;
                torsionparam_tmp.first[2] = idx3;
                torsionparam_tmp.first[3] = idx4;

                torsionparam_tmp.second[0] = idv;
                torsionparam_tmp.second[1] = dki;
                torsionparam_tmp.second[2] = dphase;
                torsionparam_tmp.second[3] = dpn;

                torsionparam[nts[i]].push_back(torsionparam_tmp);
            }else{
                cout<<"Error3"<<endl;
                exit(1);
            }
        }
    }

    return true;
}


bool loadhydrogenparam(string basepath, map<string, vector<vector<double>>> &hydparam, map<string, vector<vector<int>>> &hydparamidx)
{
    int count;
    string	oneline, field;
    string atom1, atom2, atom3, atom4;
    int residx1, residx2, residx3, residx4;
    string base;
    int idx1, idx2, idx3, idx4;
    double torsion_mean, torsion_stdval, angle_mean, angle_stdval, length_mean, length_stdval;
    vector<string> nts = {"A","C","G","U"};
    for(int i=0;i<nts.size();i++)
	{
        base = nts[i];
        ifstream inputFile(basepath+"/hyd_"+base+".txt");
        if(!inputFile.is_open())
		{
			cout<<"Parameter file can not be opened (hydrogen-" << base << ")\n"<<endl;
			return false;
		}

        vector<double> hydparamone(6);
        vector<int> hydparamidxone(8);

        while(getline(inputFile, oneline))
		{
            stringstream ssbase(oneline);
			count=0;
			while (getline(ssbase, field, ','))
			{
                //if(count==0) ssidx=stoi(field); //ssidx
				if(count==0) atom1=field; 
				if(count==1) atom2=field; 
                if(count==2) atom3=field; 
                if(count==3) atom4=field;
				if(count==4) residx1=stoi(field); 
				if(count==5) residx2=stoi(field); 
                if(count==6) residx3=stoi(field); 
                if(count==7) residx4=stoi(field); 
				if(count==8) torsion_mean=stof(field); //mean
				if(count==9) torsion_stdval=stof(field); //stdval
                if(count==10) angle_mean=stof(field); //mean
				if(count==11) angle_stdval=stof(field); //stdval
                if(count==12) length_mean=stof(field); //mean
				if(count==13) length_stdval=stof(field); //stdval
				count+=1;
            }

            if(count == 14){
                auto it1 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom1);
                idx1 = -1;
                if(it1 != hyd_full_atoms_name[base].end()){
                    idx1=distance(hyd_full_atoms_name[base].begin(), it1);
                }else{
                    cout<<atom1<<endl;
                    cout<<"Error1"<<endl;
                    exit(1);
                }

                auto it2 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom2);
                idx2 = -1;
                if(it2 != hyd_full_atoms_name[base].end()){
                    idx2=distance(hyd_full_atoms_name[base].begin(), it2);
                }else{
                    cout<<atom2<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                auto it3 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom3);
                idx3 = -1;
                if(it3 != hyd_full_atoms_name[base].end()){
                    idx3=distance(hyd_full_atoms_name[base].begin(), it3);
                }else{
                    cout<<atom3<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                auto it4 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom4);
                idx4 = -1;
                if(it4 != hyd_full_atoms_name[base].end()){
                    idx4=distance(hyd_full_atoms_name[base].begin(), it4);
                }else{
                    cout<<atom4<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                hydparamidxone = {idx1, idx2, idx3, idx4};
                hydparamone = {torsion_mean, torsion_stdval, angle_mean, angle_stdval, length_mean, length_stdval};

                hydparamidx[base].push_back(hydparamidxone);
                hydparam[base].push_back(hydparamone);
                
            }else{
                cout<<"Error3"<<endl;
                exit(1);
            }
        }
    }

    return true;
}

bool loadhydrogenaddangle(string basepath, map<string, vector<vector<double>>> &hydparam, map<string, vector<vector<int>>> &hydparamidx)
{
    int count;
    string	oneline, field;
    string atom1, atom2, atom3;
    string base;
    int idx1, idx2, idx3;
    double angle_mean, angle_stdval;
    vector<string> nts = {"A","C","G","U"};
    for(int i=0;i<nts.size();i++)
	{
        base = nts[i];
        ifstream inputFile(basepath+"/bin_info_inner_"+base+".txt");
        if(!inputFile.is_open())
		{
			cout<<"Parameter file can not be opened (hydrogen (add)-" << base << ")\n"<<endl;
			return false;
		}

        vector<double> hydparamone(2);
        vector<int> hydparamidxone(3);

        while(getline(inputFile, oneline))
		{
            stringstream ssbase(oneline);
			count=0;
			while (getline(ssbase, field, ','))
			{
                //if(count==0) ssidx=stoi(field); //ssidx
				if(count==0) atom1=field; 
				if(count==1) atom2=field; 
                if(count==2) atom3=field; 
                if(count==3) angle_mean=stof(field); //mean
				if(count==4) angle_stdval=stof(field); //stdval
				count+=1;
            }

            if(count == 5){
                auto it1 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom1);
                idx1 = -1;
                if(it1 != hyd_full_atoms_name[base].end()){
                    idx1=distance(hyd_full_atoms_name[base].begin(), it1);
                }else{
                    cout<<atom1<<endl;
                    cout<<"Error1"<<endl;
                    exit(1);
                }

                auto it2 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom2);
                idx2 = -1;
                if(it2 != hyd_full_atoms_name[base].end()){
                    idx2=distance(hyd_full_atoms_name[base].begin(), it2);
                }else{
                    cout<<atom2<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                auto it3 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom3);
                idx3 = -1;
                if(it3 != hyd_full_atoms_name[base].end()){
                    idx3=distance(hyd_full_atoms_name[base].begin(), it3);
                }else{
                    cout<<atom3<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                hydparamidxone = {idx1, idx2, idx3};
                hydparamone = {angle_mean, angle_stdval};

                hydparamidx[base].push_back(hydparamidxone);
                hydparam[base].push_back(hydparamone);
                
            }else{
                cout<<"Error3"<<endl;
                exit(1);
            }
        }
    }

    return true;
}

bool loadhydrogenaddlength(string basepath, map<string, vector<vector<double>>> &hydparam, map<string, vector<vector<int>>> &hydparamidx)
{
    int count;
    string	oneline, field;
    string atom1, atom2;
    string base;
    int idx1, idx2;
    double length_mean, length_stdval;
    vector<string> nts = {"A","C","G","U"};
    for(int i=0;i<nts.size();i++)
	{
        base = nts[i];
        ifstream inputFile(basepath+"/bin_info_length_"+base+".txt");
        if(!inputFile.is_open())
		{
			cout<<"Parameter file can not be opened (hydrogen (add)-" << base << ")\n"<<endl;
			return false;
		}

        vector<double> hydparamone(2);
        vector<int> hydparamidxone(2);

        while(getline(inputFile, oneline))
		{
            stringstream ssbase(oneline);
			count=0;
			while (getline(ssbase, field, ','))
			{
                //if(count==0) ssidx=stoi(field); //ssidx
				if(count==0) atom1=field; 
				if(count==1) atom2=field; 
                if(count==2) length_mean=stof(field); //mean
				if(count==3) length_stdval=stof(field); //stdval
				count+=1;
            }

            if(count == 4){
                auto it1 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom1);
                idx1 = -1;
                if(it1 != hyd_full_atoms_name[base].end()){
                    idx1=distance(hyd_full_atoms_name[base].begin(), it1);
                }else{
                    cout<<atom1<<endl;
                    cout<<"Error1"<<endl;
                    exit(1);
                }

                auto it2 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom2);
                idx2 = -1;
                if(it2 != hyd_full_atoms_name[base].end()){
                    idx2=distance(hyd_full_atoms_name[base].begin(), it2);
                }else{
                    cout<<atom2<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                hydparamidxone = {idx1, idx2};
                hydparamone = {length_mean, length_stdval};

                hydparamidx[base].push_back(hydparamidxone);
                hydparam[base].push_back(hydparamone);
                
            }else{
                cout<<"Error3"<<endl;
                exit(1);
            }
        }
    }

    return true;
}


bool getConnect(string basepath, map<string, vector<vector<int>>> &molConnectone, map<string, vector<vector<int>>> &molConnecttwo, map<string, vector<vector<int>>> &molConnectthree)
{
    vector<string> nts = {"A","C","G","U"};
    int count;
    string	oneline, field;
    string atom1, atom2, atom3, atom4;
    string base;
    int idx1, idx2, idx3, idx4;
    ifstream inputFile;

    for(int i=0;i<nts.size();i++)
	{
        base = nts[i];

        inputFile.open(basepath+"/"+base+"_one.connect");
        if(!inputFile)
		{
			cout<<"Parameter file can not be opened (Connect one-" << base << ")\n"<<endl;
			return false;
		}

        while(getline(inputFile, oneline))
		{
            stringstream ssbase(oneline);
			count=0;
			while (getline(ssbase, field, ','))
			{
                //if(count==0) ssidx=stoi(field); //ssidx
				if(count==0) atom1=field; 
				if(count==1) atom2=field;
                count+=1;
            }
            if(count == 2){
                auto it1 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom1);
                idx1 = -1;
                if(it1 != hyd_full_atoms_name[base].end()){
                    idx1=distance(hyd_full_atoms_name[base].begin(), it1);
                }else{
                    cout<<atom1<<endl;
                    cout<<"Error1"<<endl;
                    exit(1);
                }

                auto it2 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom2);
                idx2 = -1;
                if(it2 != hyd_full_atoms_name[base].end()){
                    idx2=distance(hyd_full_atoms_name[base].begin(), it2);
                }else{
                    cout<<atom2<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }
                if(idx1 >= 0 && idx2 >= 0) molConnectone[base].push_back({idx1, idx2});
            }
        }
        inputFile.close();
        inputFile.clear();

        inputFile.open(basepath+"/"+base+"_two.connect");
        if(!inputFile)
		{
			cout<<"Parameter file can not be opened (Connect two-" << base << ")\n"<<endl;
			return false;
		}

        while(getline(inputFile, oneline))
		{
            stringstream ssbase(oneline);
			count=0;
			while (getline(ssbase, field, ','))
			{
                //if(count==0) ssidx=stoi(field); //ssidx
				if(count==0) atom1=field; 
				if(count==1) atom2=field;
                if(count==2) atom3=field;
                count+=1;
            }
            if(count == 3){
                auto it1 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom1);
                idx1 = -1;
                if(it1 != hyd_full_atoms_name[base].end()){
                    idx1=distance(hyd_full_atoms_name[base].begin(), it1);
                }else{
                    cout<<atom1<<endl;
                    cout<<"Error1"<<endl;
                    exit(1);
                }

                auto it2 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom2);
                idx2 = -1;
                if(it2 != hyd_full_atoms_name[base].end()){
                    idx2=distance(hyd_full_atoms_name[base].begin(), it2);
                }else{
                    cout<<atom2<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                auto it3 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom3);
                idx3 = -1;
                if(it3 != hyd_full_atoms_name[base].end()){
                    idx3=distance(hyd_full_atoms_name[base].begin(), it3);
                }else{
                    cout<<atom3<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                if(idx1 >= 0 && idx2 >= 0 && idx3 >= 0) molConnecttwo[base].push_back({idx1, idx2, idx3});
            }
        }
        inputFile.close();
        inputFile.clear();

        inputFile.open(basepath+"/"+base+"_three.connect");
        if(!inputFile)
		{
			cout<<"Parameter file can not be opened (Connect three-" << base << ")\n"<<endl;
			return false;
		}

        while(getline(inputFile, oneline))
		{
            stringstream ssbase(oneline);
			count=0;
			while (getline(ssbase, field, ','))
			{
                //if(count==0) ssidx=stoi(field); //ssidx
				if(count==0) atom1=field; 
				if(count==1) atom2=field;
                if(count==2) atom3=field;
                if(count==3) atom4=field;
                count+=1;
            }
            if(count == 4){
                auto it1 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom1);
                idx1 = -1;
                if(it1 != hyd_full_atoms_name[base].end()){
                    idx1=distance(hyd_full_atoms_name[base].begin(), it1);
                }else{
                    cout<<atom1<<endl;
                    cout<<"Error1"<<endl;
                    exit(1);
                }

                auto it2 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom2);
                idx2 = -1;
                if(it2 != hyd_full_atoms_name[base].end()){
                    idx2=distance(hyd_full_atoms_name[base].begin(), it2);
                }else{
                    cout<<atom2<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                auto it3 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom3);
                idx3 = -1;
                if(it3 != hyd_full_atoms_name[base].end()){
                    idx3=distance(hyd_full_atoms_name[base].begin(), it3);
                }else{
                    cout<<atom3<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                auto it4 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom4);
                idx4 = -1;
                if(it4 != hyd_full_atoms_name[base].end()){
                    idx4=distance(hyd_full_atoms_name[base].begin(), it4);
                }else{
                    cout<<atom4<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                if(idx1 >= 0 && idx2 >= 0 && idx3 >= 0 && idx4 >=0) molConnectthree[base].push_back({idx1, idx2, idx3, idx4});
            }
        }
        inputFile.close();
        inputFile.clear();
    }
    return true;
}

bool getConnectAdj(string basepath, map<string, vector<vector<int>>> &molAdjConnectone, map<string, vector<vector<int>>> &molAdjConnecttwo, map<string, vector<vector<int>>> &molAdjConnectthree, map<string, vector<vector<int>>> &molAdjConnectoneidx, map<string, vector<vector<int>>> &molAdjConnecttwoidx, map<string, vector<vector<int>>> &molAdjConnectthreeidx)
{
    vector<string> nts = {"A","C","G","U"};
    int count;
    string	oneline, field;
    string atom1, atom2, atom3, atom4;
    string base;
    int idx1, idx2, idx3, idx4;
    int residx1, residx2, residx3, residx4;
    ifstream inputFile;

    for(int i=0;i<nts.size();i++)
	{
        base = nts[i];

        inputFile.open(basepath+"/"+base+"_one_adj.connect");
        if(!inputFile)
		{
			cout<<"Parameter file can not be opened (Adj connect one-" << base << ")\n"<<endl;
			return false;
		}

        while(getline(inputFile, oneline))
		{
            stringstream ssbase(oneline);
			count=0;
			while (getline(ssbase, field, ','))
			{
                //if(count==0) ssidx=stoi(field); //ssidx
				if(count==0) atom1=field; 
				if(count==1) atom2=field;
                if(count==2) residx1=stoi(field); 
				if(count==3) residx2=stoi(field);
                count+=1;
            }
            if(count == 4){
                auto it1 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom1);
                idx1 = -1;
                if(it1 != hyd_full_atoms_name[base].end()){
                    idx1=distance(hyd_full_atoms_name[base].begin(), it1);
                }else{
                    cout<<atom1<<endl;
                    cout<<"Error1"<<endl;
                    exit(1);
                }

                auto it2 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom2);
                idx2 = -1;
                if(it2 != hyd_full_atoms_name[base].end()){
                    idx2=distance(hyd_full_atoms_name[base].begin(), it2);
                }else{
                    cout<<atom2<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                molAdjConnectone[base].push_back({idx1, idx2});
                molAdjConnectoneidx[base].push_back({residx1, residx2});
            }
        }
        inputFile.close();
        inputFile.clear();

        inputFile.open(basepath+"/"+base+"_two_adj.connect");
        if(!inputFile)
		{
			cout<<"Parameter file can not be opened (Adj connect two-" << base << ")\n"<<endl;
			return false;
		}

        while(getline(inputFile, oneline))
		{
            stringstream ssbase(oneline);
			count=0;
			while (getline(ssbase, field, ','))
			{
                //if(count==0) ssidx=stoi(field); //ssidx
				if(count==0) atom1=field; 
				if(count==1) atom2=field;
                if(count==2) atom3=field;
                if(count==3) residx1=stoi(field); 
				if(count==4) residx2=stoi(field);
                if(count==5) residx3=stoi(field);
                count+=1;
            }
            if(count == 6){
                auto it1 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom1);
                idx1 = -1;
                if(it1 != hyd_full_atoms_name[base].end()){
                    idx1=distance(hyd_full_atoms_name[base].begin(), it1);
                }else{
                    cout<<atom1<<endl;
                    cout<<"Error1"<<endl;
                    exit(1);
                }

                auto it2 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom2);
                idx2 = -1;
                if(it2 != hyd_full_atoms_name[base].end()){
                    idx2=distance(hyd_full_atoms_name[base].begin(), it2);
                }else{
                    cout<<atom2<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                auto it3 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom3);
                idx3 = -1;
                if(it3 != hyd_full_atoms_name[base].end()){
                    idx3=distance(hyd_full_atoms_name[base].begin(), it3);
                }else{
                    cout<<atom3<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                molAdjConnecttwo[base].push_back({idx1, idx2, idx3});
                molAdjConnecttwoidx[base].push_back({residx1, residx2, residx3});
            }
        }
        inputFile.close();
        inputFile.clear();

        inputFile.open(basepath+"/"+base+"_three_adj.connect");
        if(!inputFile)
		{
			cout<<"Parameter file can not be opened (Adj connect three-" << base << ")\n"<<endl;
			return false;
		}

        while(getline(inputFile, oneline))
		{
            stringstream ssbase(oneline);
			count=0;
			while (getline(ssbase, field, ','))
			{
                //if(count==0) ssidx=stoi(field); //ssidx
				if(count==0) atom1=field; 
				if(count==1) atom2=field;
                if(count==2) atom3=field;
                if(count==3) atom4=field;
                if(count==4) residx1=stoi(field); 
				if(count==5) residx2=stoi(field);
                if(count==6) residx3=stoi(field);
                if(count==7) residx4=stoi(field);
                count+=1;
            }
            if(count == 8){
                auto it1 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom1);
                idx1 = -1;
                if(it1 != hyd_full_atoms_name[base].end()){
                    idx1=distance(hyd_full_atoms_name[base].begin(), it1);
                }else{
                    cout<<atom1<<endl;
                    cout<<"Error1"<<endl;
                    exit(1);
                }

                auto it2 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom2);
                idx2 = -1;
                if(it2 != hyd_full_atoms_name[base].end()){
                    idx2=distance(hyd_full_atoms_name[base].begin(), it2);
                }else{
                    cout<<atom2<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                auto it3 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom3);
                idx3 = -1;
                if(it3 != hyd_full_atoms_name[base].end()){
                    idx3=distance(hyd_full_atoms_name[base].begin(), it3);
                }else{
                    cout<<atom3<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                auto it4 = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom4);
                idx4 = -1;
                if(it4 != hyd_full_atoms_name[base].end()){
                    idx4=distance(hyd_full_atoms_name[base].begin(), it4);
                }else{
                    cout<<atom4<<endl;
                    cout<<"Error2"<<endl;
                    exit(1);
                }

                molAdjConnectthree[base].push_back({idx1, idx2, idx3, idx4});
                molAdjConnectthreeidx[base].push_back({residx1, residx2, residx3, residx4});
            }
        }
        inputFile.close();
        inputFile.clear();
    }
    return true;
}

bool loadfeatures(string basepath, map<string, map<int, vector<double>>> &atomfeatures)
{
    int count;
    string	oneline, field;
    string atom;
    double Charge, vdWR0, vdWEps;
    string base;
    int idx;
    double mean, stdval;
    vector<string> nts = {"A","C","G","U"};
    for(int i=0;i<nts.size();i++)
	{
        base = nts[i];
        ifstream inputFile(basepath+"/"+base+".fea");
        if(!inputFile.is_open())
		{
			cout<<"Parameter file can not be opened (feature-" << base << ")\n"<<endl;
			return false;
		}

        while(getline(inputFile, oneline))
		{
            stringstream ssbase(oneline);
			count=0;
			while (getline(ssbase, field, ','))
			{
                //if(count==0) ssidx=stoi(field); //ssidx
				if(count==0) atom=field; 
				if(count==1) Charge=stof(field); 
				if(count==2) vdWR0=stof(field); 
                if(count==3) vdWEps=stof(field); 
				count+=1;
            }

            if(count == 4){
                auto it = find(hyd_full_atoms_name[base].begin(), hyd_full_atoms_name[base].end(), atom);
                idx = -1;
                if(it != hyd_full_atoms_name[base].end()){
                    idx=distance(hyd_full_atoms_name[base].begin(), it);
                }else{
                    cout<<atom<<endl;
                    cout<<"Error1"<<endl;
                    exit(1);
                }
                atomfeatures[base][idx]={Charge, vdWR0, vdWEps};
            }
        }
    }

    return true;
}

bool loadvwd(string basepath, map<string, vector<double>> &parameter)
{
	int i, j;
	FILE *file;
	vector<string> nts = {"A","C","G","U"};
	char	oneline[600];
	char atomname[5];
	atomname[0]=' ';atomname[1]=' ';atomname[2]=' ';atomname[3]=' ';atomname[4]='\0';
	string atomnamest;
	double val;
	int atomid;
	bool flagps;

    map<char, map<int, int>> mcidx2atomid = {
        {'A', {{0,0}, {10,1}, {11,2}, {1,3}, {2,4}, {3,5}, {7,6}, {4,7}, {5,8}, {6,9}, {9,10}, {8,11}, {12,12}, {14,13}, {18,14}, {19,15}, {17,16}, {20,17}, {15,18}, {16,19}, {21,20}, {13,21}, {22,22}, {23,23}, {24,24}, {25,25}, {26,26}, {27,27}, {28,28}, {29,29}, {30,30}, {31,31}, {32,32}}},
        {'C', {{0,0}, {10,1}, {11,2}, {1,3}, {2,4}, {3,5}, {7,6}, {4,7}, {5,8}, {6,9}, {9,10}, {8,11}, {12,12}, {13,13}, {18,14}, {17,15}, {15,16}, {19,17}, {16,18}, {14,19}, {20,20}, {21,21}, {22,22}, {23,23}, {24,24}, {25,25}, {26,26}, {27,27}, {28,28}, {29,29}, {30,30}}},
        {'G', {{0,0}, {10,1}, {11,2}, {1,3}, {2,4}, {3,5}, {7,6}, {4,7}, {5,8}, {6,9}, {9,10}, {8,11}, {12,12}, {14,13}, {18,14}, {19,15}, {17,16}, {21,17}, {15,18}, {16,19}, {22,20}, {20,21}, {13,22}, {23,23}, {24,24}, {25,25}, {26,26}, {27,27}, {28,28}, {29,29}, {30,30}, {31,31}, {32,32}, {33,33}}},
        {'U', {{0,0}, {10,1}, {11,2}, {1,3}, {2,4}, {3,5}, {7,6}, {4,7}, {5,8}, {6,9}, {9,10}, {8,11}, {12,12}, {13,13}, {18,14}, {17,15}, {15,16}, {19,17}, {16,18}, {14,19}, {20,20}, {21,21}, {22,22}, {23,23}, {24,24}, {25,25}, {26,26}, {27,27}, {28,28}, {29,29}}}
    };

	for(i=0;i<nts.size();i++){
		string nt = nts[i];

		file = fopen((basepath+"/vdw_"+nt+".txt").c_str(), "rt");
		if(file==NULL)
		{
			printf( "Parameter file can not be opened %s\n", basepath.c_str());
			return false;
		}	
		parameter[nt].resize(full_atoms[nt].size());
		for(j=0; j<full_atoms[nt].size();j++) parameter[nt][j]=-1;	

		int idx=0;
		while(!feof(file))
		{
			flagps=fgets(oneline,600,file);
			if(!flagps) break;
			sscanf(oneline,"%c%c%c%c%lf", &atomname[0], &atomname[1], &atomname[2], &atomname[3], &val);
			parameter[nt][mcidx2atomid[nts[i][0]][idx]]=val;
			idx+=1;
		}
		fclose(file);
	}

	return true;
}

bool loadstackparameters(string basepath, map<int, map<string, map<string, vector<vector<double>>>>> &stparameter){
	vector<string> nts = {"A","U","G","C"};
	string tdir,pdbpath;
	vector<string>refname;
	FILE *filein;
	string oneline_str;
	string	oneline;	
	string st;
	int count;
	int paramtype, binnum;
	double minval, maxval, span, mean, stdval, mode;
	string field;
	bool flagf;

	map< char, int> Reschar2Idint={{'A',0},{'C',1},{'G',2},{'U',3}};

	for(int ss=0; ss<2; ss++){
		for(int nt1=0;nt1<4;nt1++)
		{
			for(int nt2=nt1;nt2<4;nt2++)
			{
				stparameter[ss][nts[nt1]][nts[nt2]].resize(2);
				stparameter[ss][nts[nt2]][nts[nt1]].resize(2);
				for(int i=0;i<2;i++){
					stparameter[ss][nts[nt1]][nts[nt2]][i].resize(2);
					stparameter[ss][nts[nt2]][nts[nt1]][i].resize(2);
				}
			}
		}
	}

	for(int ss=0; ss<2; ss++){
		for(int nt1=0;nt1<4;nt1++)
		{
			for(int nt2=0;nt2<4;nt2++)
			{
				if(ss == 0 && nt1 > nt2) continue;

				tdir=basepath+"/" + to_string(ss) + "/"+nts[nt1]+"_"+nts[nt2];
				
				ifstream inputFile_leng(tdir+"/length/bin_info.txt");
				if(!inputFile_leng.is_open())
				{
					cout<<"Bin info file can not be opened (Stacking parameters [length])\n"<<endl;
					return false;
				}

				while(getline(inputFile_leng, oneline))
				{
					stringstream ssbase(oneline);
					count=0;
					while (getline(ssbase, field, ','))
					{
						if(count==0) paramtype=stoi(field);
						if(count==1) minval=stof(field); //minval
						if(count==2) maxval=stof(field); //maxval
						if(count==3) span=stof(field); //span
						if(count==4) binnum=stoi(field); //binnum
						if(count==5) mean=stof(field); //mean
						if(count==6) stdval=stof(field); //stdval
						count+=1;
					}

					if(ss == 0){
						stparameter[ss][nts[nt1]][nts[nt2]][0][0] = mean;
						stparameter[ss][nts[nt1]][nts[nt2]][0][1] = stdval;
						stparameter[ss][nts[nt2]][nts[nt1]][0][0] = mean;
                        stparameter[ss][nts[nt2]][nts[nt1]][0][1] = stdval;

					}else{
						stparameter[ss][nts[nt1]][nts[nt2]][0][0] = mean;
						stparameter[ss][nts[nt1]][nts[nt2]][0][1] = stdval;
					}
				}

				ifstream inputFile_ang(tdir+"/angle/bin_info.txt");
				if(!inputFile_ang.is_open())
				{
					cout<<"Bin info file can not be opened (Stacking parameters [angle])\n"<<endl;
					return false;
				}

				while(getline(inputFile_ang, oneline))
				{
					stringstream ssbase(oneline);
					count=0;
					while (getline(ssbase, field, ','))
					{
						if(count==0) paramtype=stoi(field);
						if(count==1) minval=stof(field); //minval
						if(count==2) maxval=stof(field); //maxval
						if(count==3) span=stof(field); //span
						if(count==4) binnum=stoi(field); //binnum
						if(count==5) mean=stof(field); //mean
						if(count==6) stdval=stof(field); //stdval
						count+=1;
					}

					if(ss == 0){
						stparameter[ss][nts[nt1]][nts[nt2]][1][0] = mean;
						stparameter[ss][nts[nt1]][nts[nt2]][1][1] = stdval;
						stparameter[ss][nts[nt2]][nts[nt1]][1][0] = mean;
                        stparameter[ss][nts[nt2]][nts[nt1]][1][1] = stdval;

					}else{
						stparameter[ss][nts[nt1]][nts[nt2]][1][0] = mean;
						stparameter[ss][nts[nt1]][nts[nt2]][1][1] = stdval;
					}
				}
			}
		}
	}
	return true;
}

bool loadbpparameters(string basepath, map<string, map<string, vector<vector<double>>>> &stparameter){
	vector<string> nts = {"A","U","G","C"};
	string tdir,pdbpath;
	vector<string>refname;
	FILE *filein;
	string oneline_str;
	string	oneline;	
	string st;
	int count;
	double mean, stdval;
	string field;
	bool flagf;

	for(int nt1=0;nt1<4;nt1++)
	{
		for(int nt2=nt1;nt2<4;nt2++)
		{
			stparameter[nts[nt1]][nts[nt2]].resize(1);
			stparameter[nts[nt2]][nts[nt1]].resize(1);
			for(int i=0;i<1;i++){
				stparameter[nts[nt1]][nts[nt2]][i].resize(2);
				stparameter[nts[nt2]][nts[nt1]][i].resize(2);
			}
		}
	}

	for(int nt1=0;nt1<4;nt1++)
	{
		for(int nt2=nt1;nt2<4;nt2++)
		{

			tdir=basepath+"/" + "/"+nts[nt1]+"_"+nts[nt2];
			
			ifstream inputFile(tdir+"/bp_length.txt");
			if(!inputFile.is_open())
			{
				cout<<"Bin info file can not be opened (Basepair parameters [length])\n"<<endl;
				return false;
			}

			while(getline(inputFile, oneline))
			{
				stringstream ssbase(oneline);
				count=0;
				while (getline(ssbase, field, ','))
				{
					if(count==0) mean=stof(field); //mean
					if(count==1) stdval=stof(field); //stdval
					count+=1;
				}

				stparameter[nts[nt1]][nts[nt2]][0][0] = mean;
				stparameter[nts[nt1]][nts[nt2]][0][1] = stdval;
				stparameter[nts[nt2]][nts[nt1]][0][0] = mean;
				stparameter[nts[nt2]][nts[nt1]][0][1] = stdval;
			}
		}
	}

    /*
	for(int nt1=0;nt1<4;nt1++)
	{
		for(int nt2=nt1;nt2<4;nt2++)
		{
			cout<<nts[nt1] << "-" << nts[nt2] << "-" << stparameter[nts[nt1]][nts[nt2]][0][0] << "-" << stparameter[nts[nt1]][nts[nt2]][0][1]<<endl;
		}
	}
    */

	return true;
}

bool loadbphbondparameters(string basepath, map<char, map<char, map<int, map<int, HydroBondParameters>>>> &stparameter)
{
	vector<string> nts = {"A","U","G","C"};
	string tdir,pdbpath;
	vector<string>refname;
	FILE *filein;
	string oneline_str;
	string	oneline;	
	string st;
	int count;
	double mean, stdval;
	string field;
	string acceptor_base, donor_base;
	string acceptor_atom, donor_atom;
	double adlength_mean, adlength_std;
	double dhlength_mean, dhlength_std;
	double ahlength_mean, ahlength_std;
	double adhangle_mean, adhangle_std;
	double dhaangle_mean, dhaangle_std;
	double hadangle_mean, hadangle_std;
	bool flagf;

	for(int nt1=0;nt1<4;nt1++)
	{
		for(int nt2=nt1;nt2<4;nt2++)
		{
			tdir=basepath+"/"+nts[nt1]+"_"+nts[nt2];
			
			ifstream inputFile(tdir+"/bp_hbond.txt");
			if(!inputFile.is_open())
			{
				cout<<"Bin info file can not be opened (Basepair hbond parameters)\n"<<endl;
				return false;
			}

			while(getline(inputFile, oneline))
			{
				stringstream ssbase(oneline);
				count=0;
				while (getline(ssbase, field, ','))
				{
					if(count==0) acceptor_base=field; 
					if(count==1) donor_base=field; 
					if(count==2) acceptor_atom=field; 
					if(count==3) donor_atom=field; 
					if(count==4) adlength_mean=stof(field); 
					if(count==5) adlength_std=stof(field); 
					if(count==6) dhlength_mean=stof(field);
					if(count==7) dhlength_std=stof(field);
					if(count==8) ahlength_mean=stof(field);
					if(count==9) ahlength_std=stof(field);
					if(count==10) adhangle_mean=stof(field);
					if(count==11) adhangle_std=stof(field);
					if(count==12) dhaangle_mean=stof(field);
					if(count==13) dhaangle_std=stof(field);
					if(count==14) hadangle_mean=stof(field);
					if(count==15) hadangle_std=stof(field);
					count+=1;
				}

				stparameter[acceptor_base[0]][donor_base[0]][hydfullatomname2idx[acceptor_base][acceptor_atom]][hydfullatomname2idx[donor_base][donor_atom]].adlength_mean = adlength_mean;
				stparameter[acceptor_base[0]][donor_base[0]][hydfullatomname2idx[acceptor_base][acceptor_atom]][hydfullatomname2idx[donor_base][donor_atom]].adlength_std = adlength_std;
				stparameter[acceptor_base[0]][donor_base[0]][hydfullatomname2idx[acceptor_base][acceptor_atom]][hydfullatomname2idx[donor_base][donor_atom]].dhlength_mean = dhlength_mean;
				stparameter[acceptor_base[0]][donor_base[0]][hydfullatomname2idx[acceptor_base][acceptor_atom]][hydfullatomname2idx[donor_base][donor_atom]].dhlength_std = dhlength_std;
				stparameter[acceptor_base[0]][donor_base[0]][hydfullatomname2idx[acceptor_base][acceptor_atom]][hydfullatomname2idx[donor_base][donor_atom]].ahlength_mean = ahlength_mean;
				stparameter[acceptor_base[0]][donor_base[0]][hydfullatomname2idx[acceptor_base][acceptor_atom]][hydfullatomname2idx[donor_base][donor_atom]].ahlength_std = ahlength_std;
				stparameter[acceptor_base[0]][donor_base[0]][hydfullatomname2idx[acceptor_base][acceptor_atom]][hydfullatomname2idx[donor_base][donor_atom]].adhangle_mean = adhangle_mean;
				stparameter[acceptor_base[0]][donor_base[0]][hydfullatomname2idx[acceptor_base][acceptor_atom]][hydfullatomname2idx[donor_base][donor_atom]].adhangle_std = adhangle_std;
				stparameter[acceptor_base[0]][donor_base[0]][hydfullatomname2idx[acceptor_base][acceptor_atom]][hydfullatomname2idx[donor_base][donor_atom]].dhaangle_mean = dhaangle_mean;
				stparameter[acceptor_base[0]][donor_base[0]][hydfullatomname2idx[acceptor_base][acceptor_atom]][hydfullatomname2idx[donor_base][donor_atom]].dhaangle_std = dhaangle_std;
				stparameter[acceptor_base[0]][donor_base[0]][hydfullatomname2idx[acceptor_base][acceptor_atom]][hydfullatomname2idx[donor_base][donor_atom]].hadangle_mean = hadangle_mean;
				stparameter[acceptor_base[0]][donor_base[0]][hydfullatomname2idx[acceptor_base][acceptor_atom]][hydfullatomname2idx[donor_base][donor_atom]].hadangle_std = hadangle_std;
			}
		}
	}

	/*
	for(int nt1=0;nt1<4;nt1++)
	{
		for(int nt2=0;nt2<4;nt2++)
		{
			cout<<nts[nt1][0] << "-" << nts[nt2][0] << "-" << stparameter[nts[nt1][0]][nts[nt2][0]][0][0] << "-" << stparameter[nts[nt1][0]][nts[nt2][0]][0][1]<<endl;
		}
	}
	*/

	return true;
}

bool loadsugartorsion(string basepath, map<string, vector<vector<vector<double>>>> &torsionparam)
{
    string	oneline, field;
    double mean, stdval;
    int paramtype, count;
    vector<string> nts = {"A","C","G","U"};

    for(int nt=0;nt<4;nt++)
    {
        torsionparam[nts[nt]].resize(2);
        for(int j=0;j<2;j++){
            torsionparam[nts[nt]][j].resize(9);
            for(int k=0;k<9;k++){
                torsionparam[nts[nt]][j][k].resize(2);
            }
        }
    }

    for(int i=1;i<3;i++)
    {
        ifstream inputFile(basepath+"/class_"+ to_string(i) +".txt");
        if(!inputFile.is_open())
		{
			cout<<"Parameter file can not be opened (sugar torsion)\n"<<endl;
			return false;
		}

        paramtype = 0;
        while(getline(inputFile, oneline))
        {
            stringstream ssbase(oneline);
            count=0;
            while (getline(ssbase, field, ','))
            {
                if(count==0) mean=stof(field); //mean
                if(count==1) stdval=stof(field); //stdval
                count+=1;
            }
            
            for(int nt=0;nt<4;nt++){
                torsionparam[nts[nt]][i-1][paramtype][0] = mean;
                torsionparam[nts[nt]][i-1][paramtype][1] = stdval;
            }
            paramtype += 1;
        }
    }

    /*
    for(int nt=0;nt<4;nt++)
    {
        for(int j=0;j<2;j++){
            for(int k=0;k<6;k++){
                cout<<torsionparam[nts[nt]][j][k][0] << "-" << torsionparam[nts[nt]][j][k][1] << endl;
            }
        }
    }
    */

    return true;
}

bool loadbbtorsion(string basepath, map<string, vector<vector<double>>> &torsionparam)
{
    string	oneline, field;
    double mean, stdval;
    int paramtype, count;
    vector<string> nts = {"A","C","G","U"};

    for(int nt=0;nt<4;nt++)
    {
        torsionparam[nts[nt]].resize(5);
        for(int j=0;j<5;j++){
            torsionparam[nts[nt]][j].resize(2);
        }
    }

    ifstream inputFile(basepath+"/tor_bb.txt");
    if(!inputFile.is_open())
    {
        cout<<"Parameter file can not be opened (bb torsion)\n"<<endl;
        return false;
    }

    paramtype = 0;
    while(getline(inputFile, oneline))
    {
        stringstream ssbase(oneline);
        count=0;
        while (getline(ssbase, field, ','))
        {
            if(count==0) mean=stof(field); //mean
            if(count==1) stdval=stof(field); //stdval
            count+=1;
        }
        
        for(int nt=0;nt<4;nt++){
            torsionparam[nts[nt]][paramtype][0] = mean;
            torsionparam[nts[nt]][paramtype][1] = stdval;
        }
        paramtype += 1;
    }

    /*
    for(int nt=0;nt<4;nt++)
    {
        for(int j=0;j<5;j++){
            cout<<torsionparam[nts[nt]][j][0] << "-" << torsionparam[nts[nt]][j][1] << endl;
        }
    }
    */

    return true;
}

bool loadbbhbondparameters(string basepath, map<char, map<char, map<int, map<int, BBHydroBondParameters>>>> &stparameter)
{
	vector<string> nts = {"A","U","G","C"};
	string tdir,pdbpath;
	vector<string>refname;
	FILE *filein;
	string oneline_str;
	string	oneline;	
	string st;
	int count;
	double mean, stdval;
	string field;
	string atom1, atom2, atom3;
	double dhaangle_mean, dhaangle_std;
	double haaaangle_mean, haaaangle_std;
    double daaaangle_mean, daaaangle_std;
    double adlength_mean, adlength_std;
    double ahlength_mean, ahlength_std;
	bool flagf;

    for(int nt2=0;nt2<4;nt2++)
    {
        ifstream inputFile(basepath+"/hbond_O2'_base_"+nts[nt2]+".txt");
        if(!inputFile.is_open())
        {
            cout<<"Bin info file can not be opened (Backbone hbond parameters)\n"<<endl;
            return false;
        }

        while(getline(inputFile, oneline))
        {
            stringstream ssbase(oneline);
            count=0;
            while (getline(ssbase, field, ','))
            { 
                if(count==0) atom1=field; 
                if(count==1) atom2=field; 
                if(count==2) atom3=field; 
                if(count==3) dhaangle_mean=stof(field); 
                if(count==4) dhaangle_std=stof(field); 
                if(count==5) adlength_mean=stof(field); 
                if(count==6) adlength_std=stof(field); 
                if(count==7) ahlength_mean=stof(field); 
                if(count==8) ahlength_std=stof(field); 
                if(count==9) haaaangle_mean=stof(field); 
                if(count==10) haaaangle_std=stof(field); 
                if(count==11) daaaangle_mean=stof(field); 
                if(count==12) daaaangle_std=stof(field); 
                count+=1;
            }

            //cout<<nts[nt1][0] << "-" << hydfullatomname2idx[nts[nt1]][atom1] << "-" << hydfullatomname2idx[nts[nt1]][atom2] << "-" << adlength_mean << "-" << adlength_std << "-" << cadangle_mean << "-" << cadangle_std<<endl;
            for(int nt1=0;nt1<4;nt1++)
            {
                stparameter[nts[nt2][0]][nts[nt1][0]][hydfullatomname2idx[nts[nt2]][atom3]][hydfullatomname2idx[nts[nt1]][atom1]].dhaangle_mean = dhaangle_mean;
                stparameter[nts[nt2][0]][nts[nt1][0]][hydfullatomname2idx[nts[nt2]][atom3]][hydfullatomname2idx[nts[nt1]][atom1]].dhaangle_std = dhaangle_std;
                stparameter[nts[nt2][0]][nts[nt1][0]][hydfullatomname2idx[nts[nt2]][atom3]][hydfullatomname2idx[nts[nt1]][atom1]].adlength_mean = adlength_mean;
                stparameter[nts[nt2][0]][nts[nt1][0]][hydfullatomname2idx[nts[nt2]][atom3]][hydfullatomname2idx[nts[nt1]][atom1]].adlength_std = adlength_std;
                stparameter[nts[nt2][0]][nts[nt1][0]][hydfullatomname2idx[nts[nt2]][atom3]][hydfullatomname2idx[nts[nt1]][atom1]].ahlength_mean = ahlength_mean;
                stparameter[nts[nt2][0]][nts[nt1][0]][hydfullatomname2idx[nts[nt2]][atom3]][hydfullatomname2idx[nts[nt1]][atom1]].ahlength_std = ahlength_std;
                stparameter[nts[nt2][0]][nts[nt1][0]][hydfullatomname2idx[nts[nt2]][atom3]][hydfullatomname2idx[nts[nt1]][atom1]].haaaangle_mean = haaaangle_mean;
                stparameter[nts[nt2][0]][nts[nt1][0]][hydfullatomname2idx[nts[nt2]][atom3]][hydfullatomname2idx[nts[nt1]][atom1]].haaaangle_std = haaaangle_std;
                stparameter[nts[nt2][0]][nts[nt1][0]][hydfullatomname2idx[nts[nt2]][atom3]][hydfullatomname2idx[nts[nt1]][atom1]].daaaangle_mean = daaaangle_mean;
                stparameter[nts[nt2][0]][nts[nt1][0]][hydfullatomname2idx[nts[nt2]][atom3]][hydfullatomname2idx[nts[nt1]][atom1]].daaaangle_std = daaaangle_std;
            
            /*
            cout<<
                stparameter[nts[nt2][0]][nts[nt1][0]][hydfullatomname2idx[nts[nt2]][atom3]][hydfullatomname2idx[nts[nt1]][atom1]].dhaangle_mean << "-"<<
                stparameter[nts[nt2][0]][nts[nt1][0]][hydfullatomname2idx[nts[nt2]][atom3]][hydfullatomname2idx[nts[nt1]][atom1]].dhaangle_std << "-"<<
                stparameter[nts[nt2][0]][nts[nt1][0]][hydfullatomname2idx[nts[nt2]][atom3]][hydfullatomname2idx[nts[nt1]][atom1]].adlength_mean << "-"<<
                stparameter[nts[nt2][0]][nts[nt1][0]][hydfullatomname2idx[nts[nt2]][atom3]][hydfullatomname2idx[nts[nt1]][atom1]].adlength_std << "-"<<
                stparameter[nts[nt2][0]][nts[nt1][0]][hydfullatomname2idx[nts[nt2]][atom3]][hydfullatomname2idx[nts[nt1]][atom1]].ahlength_mean << "-"<<
                stparameter[nts[nt2][0]][nts[nt1][0]][hydfullatomname2idx[nts[nt2]][atom3]][hydfullatomname2idx[nts[nt1]][atom1]].ahlength_std << "-"<<
                stparameter[nts[nt2][0]][nts[nt1][0]][hydfullatomname2idx[nts[nt2]][atom3]][hydfullatomname2idx[nts[nt1]][atom1]].haaaangle_mean << "-"<<
                stparameter[nts[nt2][0]][nts[nt1][0]][hydfullatomname2idx[nts[nt2]][atom3]][hydfullatomname2idx[nts[nt1]][atom1]].haaaangle_std << "-"<<
                stparameter[nts[nt2][0]][nts[nt1][0]][hydfullatomname2idx[nts[nt2]][atom3]][hydfullatomname2idx[nts[nt1]][atom1]].daaaangle_mean << "-"<<
                stparameter[nts[nt2][0]][nts[nt1][0]][hydfullatomname2idx[nts[nt2]][atom3]][hydfullatomname2idx[nts[nt1]][atom1]].daaaangle_std << "-"<<endl;
            */
            }
        }
        //cout<<"================="<<endl;
    }
    //exit(1);

	return true;
}


void loadparameters(string paramdir, ExtParams &extparameters){
    loadlength(paramdir+"/length", extparameters.blengthparam, extparameters.blengthstdval, extparameters.adjblengthparam, extparameters.adjblengthstdval);
    loadangle(paramdir+"/angle", extparameters.ntangleparam, extparameters.ntangleparamidx, extparameters.adjangleparam, extparameters.adjangleparamidx);
    loadhydtorsion(paramdir+"/hyd_torsion", extparameters.hydtorsionparam);
    loadtorsion(paramdir+"/torsion", extparameters.btorsionparam);
    loadhydrogenparam(paramdir+"/hyd", extparameters.hydparam, extparameters.hydparamidx);
    loadhydrogenaddangle(paramdir+"/hyd_add", extparameters.hydaddangle, extparameters.hydaddangleidx);
    loadhydrogenaddlength(paramdir+"/hyd_add", extparameters.hydaddlength, extparameters.hydaddlengthidx);
    getConnect(paramdir+"/connect", extparameters.molConnectone, extparameters.molConnecttwo, extparameters.molConnectthree);
    getConnectAdj(paramdir+"/connect", extparameters.molAdjConnectone, extparameters.molAdjConnecttwo, extparameters.molAdjConnectthree, extparameters.molAdjConnectoneidx, extparameters.molAdjConnecttwoidx, extparameters.molAdjConnectthreeidx);
    loadfeatures(paramdir+"/features", extparameters.atomfeatures);
    loadvwd(paramdir+"/clash", extparameters.clashparam);
    loadstackparameters(paramdir+"/stacking", extparameters.stackparam);
    loadbpparameters(paramdir+"/basepair", extparameters.bpparam);
    loadbphbondparameters(paramdir+"/basepair_hbond", extparameters.bphbond);
    loadsugartorsion(paramdir+"/sugar_torsion", extparameters.sugartorsion);
    loadbbtorsion(paramdir+"/bb_torsion", extparameters.bbtorsion);
    loadbbhbondparameters(paramdir+"/hbond_o2o4", extparameters.bbhbond);
}

