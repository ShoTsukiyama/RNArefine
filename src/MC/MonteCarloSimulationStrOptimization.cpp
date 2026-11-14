
#include "MonteCarloSimulation.hpp"
#include "Gradfuncobj.hpp"


void MonteCarloSimulation::MCprepgrad(vector<vector<Point3s>> &gradptv){

    int i, mtype;
    gradptv.clear();
    gradptv.resize(numbb);
    for(i=0;i<numbb;i++){
        char base = refinestr_temp[i].aaa;
        gradptv[i].resize(AccAtomParamNum[base][4]);
        for(mtype=0;mtype<AccAtomParamNum[base][4];mtype++)
        {
            gradptv[i][mtype].x = 0;
            gradptv[i][mtype].y = 0;
            gradptv[i][mtype].z = 0;
        }
    }
}

void MonteCarloSimulation::MCsetinitobj(const map<int, vector<Point3s>> &initstructure, vector<vector<Point3s>> &gradient, list<InitPosobj> &binitobj, double K){
    
    for (auto it = initstructure.begin(); it != initstructure.end(); ++it) {
        const pair<const int, vector<Point3s>>* pair = &(*it);

        const int key_ptr = pair->first;
        const std::vector<Point3s>* value_ptr = &(pair->second);

        for (size_t i = 0; i < value_ptr->size(); ++i) {
            const Point3s& pt = (*value_ptr)[i];
	    //cout<<refinestr_temp[key_ptr].ptv[i].x << "-" << refinestr_temp[key_ptr].ptv[i].y << "-" << refinestr_temp[key_ptr].ptv[i].z<<endl;
	    //cout<<pt.x << "-" << pt.y << "-" << pt.z<<endl;
            InitPosobj binitobjone=InitPosobj::create(&refinestr_temp[key_ptr].ptv[i], &gradient[key_ptr][i], pt.x, pt.y, pt.z, K);
            binitobj.push_back(binitobjone);
        }       
    }
}

void MonteCarloSimulation::MCsettemplatepairobj(const vector<bool> &moveresidues, vector<vector<Point3s>> &gradient, list<InitPairobj> &btemplatedistobj, double K){
    int atom1, atom2;

    for(int n=0;n<(*MCtemplatedistmap).size();n++){
        for(int i=0;i<numbb;i++){
            for(int j=i+1;j<numbb;j++){
                if(moveresidues[i] || moveresidues[j]){
                    for(int m=0;m<template_target_atoms.size();m++){
                        int k = template_target_atoms_idx[m];
                        atom1=atomid2mcidx[refinestr_temp[i].aaa][k];
                        atom2=atomid2mcidx[refinestr_temp[j].aaa][k];

                        if(m < (*MCtemplatedistflag)[n][i][j].size()){
                            if((*MCtemplatedistflag)[n][i][j][k]){
                                InitPairobj templatedistobjone= InitPairobj::create(&refinestr_temp[i].ptv[atom1], &refinestr_temp[j].ptv[atom2], &gradient[i][atom1], &gradient[j][atom2], (*MCtemplatedistmap)[n][i][j][k], K);
                                btemplatedistobj.push_back(templatedistobjone);
                            }
                        }
                    }
                }
            }
        }
    }
}

void MonteCarloSimulation::MCsetinitpairobj(const map<int, vector<Point3s>> &initstructure, const vector<bool> &moveresidues, vector<vector<Point3s>> &gradient, list<InitPairobj> &binitdistobj, double K){
    int atom1, atom2;

    for(int i=0;i<numbb;i++){
        for(int j=i+1;j<numbb;j++){
            if(moveresidues[i] || moveresidues[j]){
                for(int m=0;m<template_target_atoms.size();m++){
                    int k = template_target_atoms_idx[m];
                    atom1=atomid2mcidx[refinestr_temp[i].aaa][k];
                    atom2=atomid2mcidx[refinestr_temp[j].aaa][k];

                    if(m < (*MCinitdistflag)[i][j].size()){
                        if((*MCinitdistflag)[i][j][k]){
                            //InitPairobj initdistobjone= InitPairobj::create(&refinestr_temp[i].ptv[atom1], &refinestr_temp[j].ptv[atom2], &gradient[i][atom1], &gradient[j][atom2], (*MCinitdistmap)[i][j][m], (*MCinitdistscore)[i][j]*K);
                            InitPairobj initdistobjone= InitPairobj::create(&refinestr_temp[i].ptv[atom1], &refinestr_temp[j].ptv[atom2], &gradient[i][atom1], &gradient[j][atom2], (*MCinitdistmap)[i][j][k],  (*MCinitdistscore)[i][j]*K);
                            binitdistobj.push_back(initdistobjone);
                        }
                    }
                }
            }
        }
    }

    /*
    for (const auto& pair1 : initstructure) {
        int key_ptr1 = pair1.first;
        const vector<Point3s> value1 = pair1.second;

        for (const auto& pair2 : initstructure) {
            int key_ptr2 = pair1.first;
            const vector<Point3s> value2 = pair1.second;

            if(key_ptr1 >= key_ptr2) continue;

            for(int k=0;k<targetatoms.size();k++){
                atom=targetatoms[k];
                x = value1[atom].x - value2[atom].x;
                y = value1[atom].y - value2[atom].y;
                z = value1[atom].z - value2[atom].z;

                initlength = sqrt(x*x+y*y+z*z);

                InitPairobj initdistobjone= InitPairobj::create(&refinestr_temp[key_ptr1].ptv[atom], &refinestr_temp[key_ptr2].ptv[atom], &gradient[key_ptr1][atom], &gradient[key_ptr2][atom], initlength, K);
                binitdistobj.push_back(initdistobjone);
            }
        }
    }
    */
}

void MonteCarloSimulation::MCsetlengthobj(const vector<vector<int>> &bindingmap, vector<vector<Point3s>> &gradient, double K, list<BondLengthobj> &blengthobj)
{
    int i, j, idx;
    int res3, res4;
	int atom3, atom4;
    char base;

    for(i=0;i<numbb;i++)
	{
		base = refinestr_temp[i].aaa;
		
		for(j=0;j<15;j++)
		{
			
			if(i == 0 && j == 10){
				res3 = i; res4 = i;
				atom3=FIRSTOPParamsList[base][j][2]; atom4=FIRSTOPParamsList[base][j][3];
                
                if(bindingmap[res3][res4] == 2){
                    //cout<<"aa-" << refinestr_temp[res3].ptv[atom3].x<<endl;
			        BondLengthobj blengtheobjone = BondLengthobj::create(&refinestr_temp[res3].ptv[atom3], &refinestr_temp[res4].ptv[atom4], &gradient[res3][atom3], &gradient[res4][atom4], MCOPLENGTHPARAM[base][0][0], MCOPLENGTHPARAM[base][0][1], K);
                    blengthobj.push_back(blengtheobjone);
                }
			}else if(i == 0 && j == 11){
				res3 = i; res4 = i;
				atom3=FIRSTOPParamsList[base][j][2]; atom4=FIRSTOPParamsList[base][j][3];
				
                if(bindingmap[res3][res4] == 2){
			        //cout<<"bb-" << refinestr_temp[res3].ptv[atom3].x<<endl;
                    BondLengthobj blengtheobjone = BondLengthobj::create(&refinestr_temp[res3].ptv[atom3], &refinestr_temp[res4].ptv[atom4], &gradient[res3][atom3], &gradient[res4][atom4], MCOPLENGTHPARAM[base][1][0], MCOPLENGTHPARAM[base][1][1], K);
                    blengthobj.push_back(blengtheobjone);
                }
			}else{
				atom3 = ParamsList[base][j][2]; atom4 = ParamsList[base][j][3];
                if((i != 0) || (i == 0 && j != 0)){
                    if(j == 0){
                        res3 = i-1; res4 = i;
                        if(bindingmap[res3][res4] == 1){
				            //cout<<"vv-" << refinestr_temp[res3].ptv[atom3].x<<endl;
                            BondLengthobj blengtheobjone = BondLengthobj::create(&refinestr_temp[res3].ptv[atom3], &refinestr_temp[res4].ptv[atom4], &gradient[res3][atom3], &gradient[res4][atom4], MCLENGTHPARAM[base][j][0], MCLENGTHPARAM[base][j][1], K);
                            blengthobj.push_back(blengtheobjone);
                        }
                    }else{
                        res3 = i; res4 = i;
                        if(bindingmap[res3][res4] == 2){
				            //cout<<"bb-" << refinestr_temp[res3].ptv[atom3].x<<endl;
                            BondLengthobj blengtheobjone = BondLengthobj::create(&refinestr_temp[res3].ptv[atom3], &refinestr_temp[res4].ptv[atom4], &gradient[res3][atom3], &gradient[res4][atom4], MCLENGTHPARAM[base][j][0], MCLENGTHPARAM[base][j][1], K);
                            blengthobj.push_back(blengtheobjone);
                        }
                    }
                }
			//}
		    }
        }
	        
	        
		res3 = i; res4=i;
        if(bindingmap[res3][res4] == 2){
            for(j=15;j<15+BaseParamsList[base].size();j++)
            {
                idx=j-15;
                atom3=BaseParamsList[base][idx][2]; atom4=BaseParamsList[base][idx][3];
                
		        //cout<<"ff-" << refinestr_temp[res3].ptv[atom3].x<<endl;
                BondLengthobj blengtheobjone = BondLengthobj::create(&refinestr_temp[res3].ptv[atom3], &refinestr_temp[res4].ptv[atom4], &gradient[res3][atom3], &gradient[res4][atom4], MCBASELENGTHPARAM[base][idx][0], MCBASELENGTHPARAM[base][idx][1], K);
                blengthobj.push_back(blengtheobjone);
            }
		}
		
        res3 = i; res4=i;
        if(bindingmap[res3][res4] == 2){
            for(j=0;j<HydParamsList[base].size();j++)
            {
                atom3=HydParamsList[base][j][2]; atom4=HydParamsList[base][j][3];
                
				//cout<<"ww-" << refinestr_temp[res3].ptv[atom3].x<<endl;
                BondLengthobj blengtheobjone = BondLengthobj::create(&refinestr_temp[res3].ptv[atom3], &refinestr_temp[res4].ptv[atom4], &gradient[res3][atom3], &gradient[res4][atom4], MCHYDLENGTHPARAM[base][j][4], MCHYDLENGTHPARAM[base][j][1], K);
                blengthobj.push_back(blengtheobjone);
            }
		}
		
        
        for(j=0;j<AddLengthParamsList[base].size();j++)
        {
            if(i == 0 && j < 4) continue;
            if(j < 4){
                res3 = i-1; res4 = i;
                if(bindingmap[res3][res4] == 1){
			//cout<<"qq-"<<refinestr_temp[res3].ptv[atom3].x<<endl;
                    atom3=AddLengthParamsList[base][j][0]; atom4=AddLengthParamsList[base][j][1];
                    BondLengthobj blengtheobjone = BondLengthobj::create(&refinestr_temp[res3].ptv[atom3], &refinestr_temp[res4].ptv[atom4], &gradient[res3][atom3], &gradient[res4][atom4], MCADDLENGTHPARAM[base][j][0], MCADDLENGTHPARAM[base][j][1], K);
                    blengthobj.push_back(blengtheobjone);
                }
            }else{
                res3 = i; res4 = i;
                if(bindingmap[res3][res4] == 2){
			//cout<<"pp-"<<refinestr_temp[res3].ptv[atom3].x<<endl;
                    atom3=AddLengthParamsList[base][j][0]; atom4=AddLengthParamsList[base][j][1];
                    BondLengthobj blengtheobjone = BondLengthobj::create(&refinestr_temp[res3].ptv[atom3], &refinestr_temp[res4].ptv[atom4], &gradient[res3][atom3], &gradient[res4][atom4], MCADDLENGTHPARAM[base][j][0], MCADDLENGTHPARAM[base][j][1], K);
                    blengthobj.push_back(blengtheobjone);
                }
            }
        }
        
        res3 = i; res4 = i;
        if(bindingmap[res3][res4] == 2){
            for(j=0;j<AddHydLengthParamsList[base].size();j++)
            {
                atom3=AddHydLengthParamsList[base][j][0]; atom4=AddHydLengthParamsList[base][j][1];

                //cout<<res3 << "-" << res4 << "-" << atom3 << "-" << atom4<<endl;

                BondLengthobj blengtheobjone = BondLengthobj::create(&refinestr_temp[res3].ptv[atom3], &refinestr_temp[res4].ptv[atom4], &gradient[res3][atom3], &gradient[res4][atom4], MCHYDADDLENGTHPARAM[base][j][0], MCHYDADDLENGTHPARAM[base][j][1], K);
                blengthobj.push_back(blengtheobjone);

		//double dx = refinestr_temp[res3].ptv[atom3].x - refinestr_temp[res4].ptv[atom4].x;
		//double dy = refinestr_temp[res3].ptv[atom3].y - refinestr_temp[res4].ptv[atom4].y;
		//double dz = refinestr_temp[res3].ptv[atom3].z - refinestr_temp[res4].ptv[atom4].z;

		//double leng_tt = sqrt(dx*dx+dy*dy+dz*dz);
	        //cout<<res3 << "-" << res4 << "-" << atom3 << "-" << atom4 << "-" << base << "-" << leng_tt << "-" << MCHYDADDLENGTHPARAM[base][j][0]<<endl;
                //cout<<res3 << "-" << atom3 << "-" << refinestr_temp[res3].ptv[atom3].x << "-" <<  refinestr_temp[res3].ptv[atom3].y << "-" << refinestr_temp[res3].ptv[atom3].z<<endl;
		//cout<<res4 << "-" << atom4 << "-" << refinestr_temp[res4].ptv[atom4].x << "-" << refinestr_temp[res4].ptv[atom4].y << "-" << refinestr_temp[res4].ptv[atom4].z<<endl;
	    }
        }
        
	}
}

void MonteCarloSimulation::MCsetangleobj(const vector<vector<int>> &bindingmap, vector<vector<Point3s>> &gradient, double K, list<BondAngleobj> &bangleobj)
{
    int i, j, idx;
    int res2, res3, res4;
	int atom2, atom3, atom4;
    char base;

    for(i=0;i<numbb;i++)
	{
		base = refinestr_temp[i].aaa;
		
		for(j=0;j<15;j++)
		{
			if(i == 0 && j == 10){
				res2 = i; res3 = i; res4 = i;
                if(bindingmap[res3][res4] == 2){
                    atom2=FIRSTOPParamsList[base][j][1]; atom3=FIRSTOPParamsList[base][j][2]; atom4=FIRSTOPParamsList[base][j][3];

                    BondAngleobj bangleobjone = BondAngleobj::create(&refinestr_temp[res2].ptv[atom2], &refinestr_temp[res3].ptv[atom3], &refinestr_temp[res4].ptv[atom4], &gradient[res2][atom2], &gradient[res3][atom3], &gradient[res4][atom4], MCOPINNERANGPARAM[base][0][0], MCOPINNERANGPARAM[base][0][1], K);
                    bangleobj.push_back(bangleobjone);
                }
			}else if(i == 0 && j == 11){
				res2 = i; res3 = i; res4 = i;
                if(bindingmap[res3][res4] == 2){
                    atom2=FIRSTOPParamsList[base][j][1]; atom3=FIRSTOPParamsList[base][j][2]; atom4=FIRSTOPParamsList[base][j][3];

                    BondAngleobj bangleobjone = BondAngleobj::create(&refinestr_temp[res2].ptv[atom2], &refinestr_temp[res3].ptv[atom3], &refinestr_temp[res4].ptv[atom4], &gradient[res2][atom2], &gradient[res3][atom3], &gradient[res4][atom4], MCOPINNERANGPARAM[base][1][0], MCOPINNERANGPARAM[base][1][1], K);
                    bangleobj.push_back(bangleobjone);
                }
			}else{
                if((i != 0) || (i == 0 && j != 0 && j != 1)){
                    atom2 = ParamsList[base][j][1]; atom3 = ParamsList[base][j][2]; atom4 = ParamsList[base][j][3];
                    if(j == 0){
                        res2 = i-1; res3 = i-1; res4 = i;
                        if(bindingmap[res3][res4] == 1){
                            BondAngleobj bangleobjone = BondAngleobj::create(&refinestr_temp[res2].ptv[atom2], &refinestr_temp[res3].ptv[atom3], &refinestr_temp[res4].ptv[atom4], &gradient[res2][atom2], &gradient[res3][atom3], &gradient[res4][atom4], MCINNERANGPARAM[base][j][0], MCINNERANGPARAM[base][j][1], K);
                            bangleobj.push_back(bangleobjone);
                        }
                    }else if(j == 1){
                        res2 = i-1; res3 = i; res4 = i;
                        if(bindingmap[res2][res3] == 1){
                            BondAngleobj bangleobjone = BondAngleobj::create(&refinestr_temp[res2].ptv[atom2], &refinestr_temp[res3].ptv[atom3], &refinestr_temp[res4].ptv[atom4], &gradient[res2][atom2], &gradient[res3][atom3], &gradient[res4][atom4], MCINNERANGPARAM[base][j][0], MCINNERANGPARAM[base][j][1], K);
                            bangleobj.push_back(bangleobjone);
                        }
                    }else if(j == 11){
                        res2 = i-1; res3 = i; res4 = i;
                        if(bindingmap[res2][res3] == 1){
                            BondAngleobj bangleobjone = BondAngleobj::create(&refinestr_temp[res2].ptv[atom2], &refinestr_temp[res3].ptv[atom3], &refinestr_temp[res4].ptv[atom4], &gradient[res2][atom2], &gradient[res3][atom3], &gradient[res4][atom4], MCINNERANGPARAM[base][j][0], MCINNERANGPARAM[base][j][1], K);
                            bangleobj.push_back(bangleobjone);
                        }
                    }else{
                        res2 = i; res3 = i; res4 = i;
                        if(bindingmap[res2][res3] == 2){
                            BondAngleobj bangleobjone = BondAngleobj::create(&refinestr_temp[res2].ptv[atom2], &refinestr_temp[res3].ptv[atom3], &refinestr_temp[res4].ptv[atom4], &gradient[res2][atom2], &gradient[res3][atom3], &gradient[res4][atom4], MCINNERANGPARAM[base][j][0], MCINNERANGPARAM[base][j][1], K);
                            bangleobj.push_back(bangleobjone);
                        }
                    }
                }
			}
		}
		
		res2 = i; res3 = i; res4=i;
		for(j=15;j<15+BaseParamsList[base].size();j++)
		{
            if(bindingmap[res2][res3] == 2){
                idx=j-15;
			    atom2=BaseParamsList[base][idx][1]; atom3=BaseParamsList[base][idx][2]; atom4=BaseParamsList[base][idx][3];
		
                BondAngleobj bangleobjone = BondAngleobj::create(&refinestr_temp[res2].ptv[atom2], &refinestr_temp[res3].ptv[atom3], &refinestr_temp[res4].ptv[atom4], &gradient[res2][atom2], &gradient[res3][atom3], &gradient[res4][atom4], MCBASEINNERANGPARAM[base][idx][0], MCBASEINNERANGPARAM[base][idx][1], K);
                bangleobj.push_back(bangleobjone);
            }
		}
		
        res2 = i; res3 = i; res4=i;
		for(j=0;j<MCHYDLENGTHPARAM[base].size();j++)
		{
            if(bindingmap[res2][res3] == 2){
                atom2=HydParamsList[base][j][1]; atom3=HydParamsList[base][j][2]; atom4=HydParamsList[base][j][3];
		
                BondAngleobj bangleobjone = BondAngleobj::create(&refinestr_temp[res2].ptv[atom2], &refinestr_temp[res3].ptv[atom3], &refinestr_temp[res4].ptv[atom4], &gradient[res2][atom2], &gradient[res3][atom3], &gradient[res4][atom4], MCHYDINNERANGPARAM[base][j][4], MCHYDINNERANGPARAM[base][j][1], K);
                bangleobj.push_back(bangleobjone);
            }
		}

        for(j=0;j<MCADDINNERANGPARAM[base].size();j++)
        {
            if(i == 0 && j < 1) continue;
            atom2=AddAngleParamsList[base][j][0]; atom3=AddAngleParamsList[base][j][1]; atom4=AddAngleParamsList[base][j][2];
            
            if(j < 1){
                res2 = i-1; res3 = i; res4 = i;
                if(bindingmap[res2][res3] == 1){
                    BondAngleobj bangleobjone = BondAngleobj::create(&refinestr_temp[res2].ptv[atom2], &refinestr_temp[res3].ptv[atom3], &refinestr_temp[res4].ptv[atom4], &gradient[res2][atom2], &gradient[res3][atom3], &gradient[res4][atom4], MCADDINNERANGPARAM[base][j][0], MCADDINNERANGPARAM[base][j][1], K);
                    bangleobj.push_back(bangleobjone);
                }
            }else{
                res2 = i; res3 = i; res4 = i;
                if(bindingmap[res2][res3] == 2){
                    BondAngleobj bangleobjone = BondAngleobj::create(&refinestr_temp[res2].ptv[atom2], &refinestr_temp[res3].ptv[atom3], &refinestr_temp[res4].ptv[atom4], &gradient[res2][atom2], &gradient[res3][atom3], &gradient[res4][atom4], MCADDINNERANGPARAM[base][j][0], MCADDINNERANGPARAM[base][j][1], K);
                    bangleobj.push_back(bangleobjone);
                }
            }
        }

        
        res2 = i; res3 = i; res4 = i;
        if(bindingmap[res2][res3] == 2){
            for(j=0;j<AddHydAngleParamsList[base].size();j++)
            {
                atom2=AddHydAngleParamsList[base][j][0]; atom3=AddHydAngleParamsList[base][j][1]; atom4=AddHydAngleParamsList[base][j][2];

                BondAngleobj bangleobjone = BondAngleobj::create(&refinestr_temp[res2].ptv[atom2], &refinestr_temp[res3].ptv[atom3], &refinestr_temp[res4].ptv[atom4], &gradient[res2][atom2], &gradient[res3][atom3], &gradient[res4][atom4], MCHYDADDANGLEPARAM[base][j][0], MCHYDADDANGLEPARAM[base][j][1], K);
                bangleobj.push_back(bangleobjone);

                //cout<<res3 << "-" << res4 << "-" << atom3 << "-" << atom4 << "-" << MCHYDADDANGLEPARAM[base][j][0]<<endl;
            }
        }
	}
    //exit(1);
}

void MonteCarloSimulation::MCsetdangleobj(const vector<vector<int>> &bindingmap, vector<vector<Point3s>> &gradient, double K, list<BondDihedralobj> &bangleobj)
{
    int i, j, idx;
    int res1, res2, res3, res4;
	int atom1, atom2, atom3, atom4;
    int pos1, pos2, pos3, pos4;
    double ideal, rotn;
    char base;
    char base1, base2;
    int minfragidx;
    double minrmsd;

    vector<vector<int>> torsion_idx_list = {{3,4,5,0}, {4,5,0,1}, {5,0,1,2}, {0,1,2,3}, {1,2,3,4}};
    vector<vector<int>> torsion_rel_pos = {{0,0,0,1}, {0,0,1,1}, {0,1,1,1}, {1,1,1,1}, {1,1,1,1}};

    for(i=0;i<numbb-1;i++)
	{
        for(j=0;j<torsion_idx_list.size();j++){
            pos1 = i + torsion_rel_pos[j][0]; pos2 = i + torsion_rel_pos[j][1]; pos3 = i + torsion_rel_pos[j][2]; pos4 = i + torsion_rel_pos[j][3];
            base =  refinestr_temp[i].aaa;

            if((bindingmap[pos1][pos2] == 1) || (bindingmap[pos2][pos3] == 1) || (bindingmap[pos3][pos4] == 1) || (bindingmap[pos1][pos2] == 2) || (bindingmap[pos2][pos3] == 2) || (bindingmap[pos3][pos4] == 2)){
                atom1 = torsion_idx_list[j][0]; atom2 = torsion_idx_list[j][1]; atom3 = torsion_idx_list[j][2]; atom4 = torsion_idx_list[j][3];

                //cout<<pos1 << "-" << pos2 << "-" << pos3 << "-" << pos4<<endl;
                //cout<<atom1 << "-" << atom2 << "-" << atom3 << "-" << atom4<<endl;

                BondDihedralobj dihedralobjone1 = BondDihedralobj::create(&refinestr_temp[pos1].ptv[atom1], &refinestr_temp[pos2].ptv[atom2], &refinestr_temp[pos3].ptv[atom3], &refinestr_temp[pos4].ptv[atom4], &gradient[pos1][atom1], &gradient[pos2][atom2], &gradient[pos3][atom3], &gradient[pos4][atom4], K, MCBBTORSION[base][j][0], MCBBTORSION[base][j][1]);
                bangleobj.push_back(dihedralobjone1);

            }
        }
    }
}

void MonteCarloSimulation::MCsetdsugarangleobj(const vector<vector<int>> &bindingmap, vector<vector<Point3s>> &gradient, double K, list<MultiBondDihedralobj> &multibangleobj)
{
    vector<int> sugaratomidx = {2,3,4,5,7,8,6,1,12,13,0};
    vector<int> resrelpos = {0,0,0,0,0,0,0,0,0,0,1};

    for(int i=0;i<numbb-1;i++)
	{
		char base = refinestr_temp[i].aaa;

        vector<Point3s*> AtomPtrs;
        vector<Point3s*> Atomgrad;
        vector<vector<vector<double*>>> ideals;

        for(int j=0;j<sugaratomidx.size();j++){
            int k = sugaratomidx[j];
            int pos = i + resrelpos[j];
            Point3s* respt = &refinestr_temp[pos].ptv[k];
            AtomPtrs.push_back(respt);

            Point3s* resgradpt = &gradient[pos][k];
            Atomgrad.push_back(resgradpt);
        }

        ideals.resize(2);
        for(int j=0;j<2;j++){
            for(int k=0;k<9;k++){
                double* idealpt1 = &MCSUGARTORSION[base][j][k][0];
                double* idealpt2 = &MCSUGARTORSION[base][j][k][1];
                ideals[j].push_back({idealpt1, idealpt2});
            }
        }

        MultiBondDihedralobj multibangleobjone = MultiBondDihedralobj(AtomPtrs, Atomgrad, ideals, K);
        multibangleobj.push_back(multibangleobjone);

    }
}

void MonteCarloSimulation::MCstack_superimpose(int res1, int res2, char base1, char base2, const int &sttype, double &rmsd, const int &order, const int &tmpidx)
{
    int i,j;
    int count;
    vector<vector<double>> RotMatix;
    vector<double> TranVect;
    vector<double> xyztmp(3,0);
    int atom;
    BasicFunc bf;
    Point3s vec;
    int base1_temp, base2_temp;

    /*
    map<char, vector<int>> targetatoms={
        {'A',{12, 13, 14, 15, 16, 17, 18, 19, 20, 21}},
        {'C',{12, 13, 14, 15, 16, 17, 18, 19}},
        {'G',{12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22}},
        {'U',{12, 13, 14, 15, 16, 17, 18, 19}}
    };
    */
    
    static map<char, vector<int>> targetatoms={
        {'A',{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21}},
        {'C',{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}},
        {'G',{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22}},
        {'U',{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}}
    };

    if(order == 0){
        base1_temp = 0;
        base2_temp = 1;
    }else if(order == 1){
        base1_temp = 1;
        base2_temp = 0;
    }

    vector<vector<double> > xyz_list1(targetatoms[base1].size()+targetatoms[base2].size(), vector<double>(3,0));
    vector<vector<double> > xyz_list2(targetatoms[base1].size()+targetatoms[base2].size(), vector<double>(3,0));

    count = 0;
    for(j=0; j<targetatoms[base1].size();j++){
        atom=targetatoms[base1][j];
        xyz_list1[count][0]=refinestr_temp[res1].ptv[atom].x;
        xyz_list1[count][1]=refinestr_temp[res1].ptv[atom].y;
        xyz_list1[count][2]=refinestr_temp[res1].ptv[atom].z;

        xyz_list2[count][0]=STACKTEMPLATE[sttype][base1][base2][tmpidx][base1_temp].ptv[atom].x;
        xyz_list2[count][1]=STACKTEMPLATE[sttype][base1][base2][tmpidx][base1_temp].ptv[atom].y;
        xyz_list2[count][2]=STACKTEMPLATE[sttype][base1][base2][tmpidx][base1_temp].ptv[atom].z;

        count+=1;
    }

    for(j=0; j<targetatoms[base2].size();j++){
        atom=targetatoms[base2][j];
        xyz_list1[count][0]=refinestr_temp[res2].ptv[atom].x;
        xyz_list1[count][1]=refinestr_temp[res2].ptv[atom].y;
        xyz_list1[count][2]=refinestr_temp[res2].ptv[atom].z;

        xyz_list2[count][0]=STACKTEMPLATE[sttype][base1][base2][tmpidx][base2_temp].ptv[atom].x;
        xyz_list2[count][1]=STACKTEMPLATE[sttype][base1][base2][tmpidx][base2_temp].ptv[atom].y;
        xyz_list2[count][2]=STACKTEMPLATE[sttype][base1][base2][tmpidx][base2_temp].ptv[atom].z;

        count+=1;
    }

    bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

    rmsd=0;
    count=0;
    for(j=0; j<targetatoms[base1].size();j++){
        atom=targetatoms[base1][j];
        bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);
        vec.x=refinestr_temp[res1].ptv[atom].x-xyztmp[0];
        vec.y=refinestr_temp[res1].ptv[atom].y-xyztmp[1];
        vec.z=refinestr_temp[res1].ptv[atom].z-xyztmp[2];

        rmsd+=(vec.x*vec.x+vec.y*vec.y+vec.z*vec.z);
        count+=1;
    }

    for(j=0; j<targetatoms[base2].size();j++){
        atom=targetatoms[base2][j];
        bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);

        vec.x=refinestr_temp[res2].ptv[atom].x-xyztmp[0];
        vec.y=refinestr_temp[res2].ptv[atom].y-xyztmp[1];
        vec.z=refinestr_temp[res2].ptv[atom].z-xyztmp[2];

        rmsd+=(vec.x*vec.x+vec.y*vec.y+vec.z*vec.z);
        count+=1;
    }
    rmsd = sqrt(rmsd/count);
}

void MonteCarloSimulation::MCselectstackfrag(vector<vector<Point3s>> &gradient, int res1, int res2, const int &sttype, int stidx, double K, list<Stackingobj> &stackobj)
{
	int i,j;
	int atom, atom1, atom2;
	int count;
	vector<vector<double>> RotMatix;
    vector<double> TranVect;
	double rmsd;
	vector<double> xyztmp(3,0);
	Point3s vec;
	BasicFunc bf;
    char base1, base2;

    base1 = refinestr_temp[res1].aaa;
    base2 = refinestr_temp[res2].aaa;
    vector<int> indices;
    vector<int> orders;
    
    /*
	static map<char, vector<int>> targetatoms={
		{'A',{12, 13, 14, 15, 16, 17, 18, 19, 20, 21}},
        {'C',{12, 13, 14, 15, 16, 17, 18, 19}},
        {'G',{12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22}},
        {'U',{12, 13, 14, 15, 16, 17, 18, 19}}
	};
    
    
    static map<char, vector<int>> targetatoms_dist={
		{'A',{12, 16, 19}},
        {'C',{12, 16, 17}},
        {'G',{12, 16, 19}},
        {'U',{12, 16, 17}}
	};
    */

    /*
    static map<char, vector<int>> targetatoms={
        {'A',{0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21}},
        {'C',{0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 15, 16, 17, 18, 19}},
        {'G',{0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22}},
        {'U',{0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 15, 16, 17, 18, 19}}
    };
    */

    //static map<char, vector<int>> targetatoms_dist={
        //{'A',{12, 13, 14, 15, 16, 17, 18, 19, 21}},
        //{'C',{12, 13, 14, 15, 16, 17}},
        //{'G',{12, 13, 14, 15, 16, 17, 18, 19, 20}},
        //{'U',{12, 13, 14, 15, 16, 17}}
    //};

    /*
    static map<char, vector<int>> targetatoms={
		{'A',{0,1,2,3,4,5,6,7,8,12,13,14,15,16,17,18,19,20,21}},
        {'C',{0,1,2,3,4,5,6,7,8,12,13,14,15,16,17,18,19}},
        {'G',{0,1,2,3,4,5,6,7,8,12,13,14,15,16,17,18,19,20,21,22}},
        {'U',{0,1,2,3,4,5,6,7,8,12,13,14,15,16,17,18,19}}
	};
    */

    /*
    static map<char, vector<int>> targetatoms_dist={
	    {'A',{12, 13, 14, 15, 16, 17, 18, 19, 21}},
        {'C',{12, 13, 14, 15, 16, 17}},
        {'G',{12, 13, 14, 15, 16, 17, 18, 19, 20}},
        {'U',{12, 13, 14, 15, 16, 17}}
    };
    */
    
    static map<char, vector<int>> targetatoms_dist={
	    {'A',{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21}},
        {'C',{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19}},
        {'G',{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}},
        {'U',{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19}}
    };
    

    /*
    static map<char, vector<int>> targetatoms={
		{'A',{0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21}},
        {'C',{0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 15, 16, 17, 18, 19}},
        {'G',{0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22}},
        {'U',{0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 15, 16, 17, 18, 19}}
	};
    */

    vector<vector<Residueatoms>> &stacktemplate=STACKTEMPLATE[sttype][base1][base2];
    vector<vector<vector<double>>> &stackdistancemap=STACKDISTMAP[sttype][base1][base2];

    vector<double> rmsds;
    if(sttype == 0 && base1 == base2) rmsds.resize(stacktemplate.size()*2, 0.0);
    else rmsds.resize(stacktemplate.size(), 0.0);

	for(i=0;i<stacktemplate.size();i++)
	{
        rmsd=0;
		count=0;
        MCstack_superimpose(res1, res2, base1, base2, sttype, rmsd, 0, i);
        rmsds[i]=rmsd;
        //MCstack_superimpose(res2, res1, base1, base2, sttype, rmsd, 1, i);
        //rmsds[i]+=rmsd;
	}

	if(sttype == 0 && base1 == base2)
	{
		for(i=0;i<stacktemplate.size();i++)
		{
			rmsd=0;
		    count=0;
            MCstack_superimpose(res1, res2, base1, base2, sttype, rmsd, 1, i);
            rmsds[i+stacktemplate.size()]=rmsd;
            //MCstack_superimpose(res2, res1, base1, base2, sttype, rmsd, 0, i);
            //rmsds[i]+=rmsd;
			
		}
	}

    indices.resize(rmsds.size());
    orders.resize(rmsds.size());
    for (int i = 0; i < rmsds.size(); ++i) indices[i] = i;

    sort(indices.begin(), indices.end(), [&rmsds](int i1, int i2) {
        return rmsds[i1] < rmsds[i2];
    });

    for(i=0;i<indices.size();i++){
        //cout<<i << "-" << indices[i] << "-" << rmsds[indices[i]]<<endl;

        if(indices[i] >= STACKTEMPLATE[sttype][base1][base2].size()){
            indices[i] = indices[i] - STACKTEMPLATE[sttype][base1][base2].size();
            orders[i] = 1;
        }else{
            orders[i] = 0;
        }
    }

    //StackTemplateobj2MC stackobjone(refinestr_temp[res1].ptv, refinestr_temp[res2].ptv, &gradient[res1], &gradient[res2], refinestr_temp[res1].aaa, refinestr_temp[res2].aaa, &STACKTEMPLATE[sttype][base1][base2], &indices[stidx], &orders[stidx], 5, K);
    //stackobj.push_back(stackobjone);

    int template_idx = indices[0];
    int template_order = orders[0];

    //cout<<res1 << "-" << res2 << "-" << template_idx << "-" << rmsds[template_idx] << "-" << sttype << "-" << base1 << "-" << base2<<endl;

    if(template_order == 0){
        for(int i=0;i<targetatoms_dist[base1].size();i++){
            atom1=targetatoms_dist[base1][i];
            for(int j=0;j<targetatoms_dist[base2].size();j++){
                atom2=targetatoms_dist[base2][j];
                
                /*
                double dx = refinestr_temp[res1].ptv[atom1].x - refinestr_temp[res2].ptv[atom2].x;
                double dy = refinestr_temp[res1].ptv[atom1].y - refinestr_temp[res2].ptv[atom2].y;
                double dz = refinestr_temp[res1].ptv[atom1].z - refinestr_temp[res2].ptv[atom2].z;
                double length = sqrt(dx*dx + dy*dy+ dz*dz);
                */
                
                Stackingobj stackobjone(&refinestr_temp[res1].ptv[atom1], &refinestr_temp[res2].ptv[atom2], &gradient[res1][atom1], &gradient[res2][atom2], stackdistancemap[template_idx][atom1][atom2], K);
                stackobj.push_back(stackobjone);
            }
        }
    }else{
        for(int i=0;i<targetatoms_dist[base1].size();i++){
            atom1=targetatoms_dist[base1][i];
            for(int j=0;j<targetatoms_dist[base2].size();j++){
                atom2=targetatoms_dist[base2][j];

                /*
                double dx = refinestr_temp[res1].ptv[atom1].x - refinestr_temp[res2].ptv[atom2].x;
                double dy = refinestr_temp[res1].ptv[atom1].y - refinestr_temp[res2].ptv[atom2].y;
                double dz = refinestr_temp[res1].ptv[atom1].z - refinestr_temp[res2].ptv[atom2].z;
                double length = sqrt(dx*dx + dy*dy+ dz*dz);
                */

                Stackingobj stackobjone(&refinestr_temp[res1].ptv[atom1], &refinestr_temp[res2].ptv[atom2], &gradient[res1][atom1], &gradient[res2][atom2], stackdistancemap[template_idx][atom2][atom1], K);
                stackobj.push_back(stackobjone);
            }
        }
    }
}

void MonteCarloSimulation::MCsetstackingobj(vector<vector<Point3s>> &gradient, const vector<vector<double>> &stackmap, double K, list<Stackingobj> &stackobj)
{
	int i,j,k;
	int res1, res2;
	char base1, base2;
    int minfragidx; 
    double minrmsd;
    int baseorder;

    //cout<<bindmap.size()<<endl;

    int stnum=0;
    for(i=0;i<refinestr_temp.size();i++){
        for(j=i+1;j<refinestr_temp.size();j++){
            if(stackmap[i][j] > 0) stnum+=1;
        }
    }

    int count=0;

    for(i=0;i<refinestr_temp.size();i++){
        //cout<<bindmap[i].size()<<endl;
        res1=i;
        base1 =  refinestr_temp[res1].aaa;
        for(j=i+1;j<refinestr_temp.size();j++){
            //cout<<i << "-" << j<<endl;
            res2=j;
            base2 = refinestr_temp[res2].aaa;
            if(stackmap[res1][res2] > 0){
                if(i+1 == j){
                    MCselectstackfrag(gradient, res1, res2, 1, count, stackmap[res1][res2]*K, stackobj);
                }else{
                    MCselectstackfrag(gradient, res1, res2, 0, count, stackmap[res1][res2]*K, stackobj);
                    
                }

                count+=1;
            }
        }
    }
}

void MonteCarloSimulation::MCsetstackinggeoobj(vector<vector<Point3s>> &gradient, const vector<vector<double>> &stackmap, vector<vector<Point3s*>> &basestackptrsA, vector<vector<Point3s*>> &basestackptrsB, vector<vector<Point3s*>> basestackgradA, vector<vector<Point3s*>> basestackgradB, double K, list<StackGeoobjMC> &stackgeoobj)
{
	int i,j,k,l;
	int res1, res2;
	char base1, base2;
    int minfragidx; 
    double minrmsd;
    int baseorder;

    map<char, vector<int>> targetatoms={
        {'A',{12, 13, 14, 15, 16, 17, 18, 19, 20, 21}},
        {'C',{12, 13, 14, 15, 16, 17, 18, 19}},
        {'G',{12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22}},
        {'U',{12, 13, 14, 15, 16, 17, 18, 19}}
    };

    int stnum=0;
    for(i=0;i<refinestr_temp.size();i++){
        for(j=i+1;j<refinestr_temp.size();j++){
            if(stackmap[i][j] > 0) stnum+=1;
        }
    }

    basestackptrsA.resize(stnum);
    basestackptrsB.resize(stnum);
    basestackgradA.resize(stnum);
    basestackgradB.resize(stnum);

    int count=0;

    for(i=0;i<refinestr_temp.size();i++){
        res1=i;
        base1 =  refinestr_temp[res1].aaa;
        for(j=i+1;j<refinestr_temp.size();j++){
            res2=j;
            base2 = refinestr_temp[res2].aaa;
            if(stackmap[res1][res2] > 0){
                basestackptrsA[count].clear();
                basestackptrsB[count].clear();
                basestackgradA[count].clear();
                basestackgradB[count].clear();

                for(l=0;l<targetatoms[base1].size();l++){
                    k = targetatoms[base1][l];
                    Point3s* respt = &refinestr_temp[res1].ptv[k];
                    basestackptrsA[count].push_back(respt);

                    Point3s* resgradpt = &gradient[res1][k];
                    basestackgradA[count].push_back(resgradpt);
                }

                //for(k=12;k<AccAtomParamNum[base2][3];k++){
                for(l=0;l<targetatoms[base2].size();l++){
                    k = targetatoms[base2][l];
                    Point3s* respt = &refinestr_temp[res2].ptv[k];
                    basestackptrsB[count].push_back(respt);

                    Point3s* resgradpt = &gradient[res2][k];
                    basestackgradB[count].push_back(resgradpt);
                }

                if(i+1 == j){
                    StackGeoobjMC stackgeoobjone(basestackptrsA[count], basestackptrsB[count], basestackgradA[count], basestackgradB[count], refinestr_temp[res1].aaa, refinestr_temp[res2].aaa, STACKPARAM[1][base1][base2][0][0]+2*STACKPARAM[1][base1][base2][0][1], STACKPARAM[1][base1][base2][1][0]+2*STACKPARAM[1][base1][base2][1][0], stackmap[res1][res2]*K, stackmap[res1][res2]*K);
                    stackgeoobj.push_back(stackgeoobjone);
                }else{
                    StackGeoobjMC stackgeoobjone(basestackptrsA[count], basestackptrsB[count], basestackgradA[count], basestackgradB[count], refinestr_temp[res1].aaa, refinestr_temp[res2].aaa, STACKPARAM[0][base1][base2][0][0]+2*STACKPARAM[0][base1][base2][0][0], STACKPARAM[0][base1][base2][1][0]+2*STACKPARAM[0][base1][base2][1][0], stackmap[res1][res2]*K, stackmap[res1][res2]*K);
                    stackgeoobj.push_back(stackgeoobjone);
                }

                count+=1;
            }
        }
    }
}

void MonteCarloSimulation::MCbp_superimpose(int res1, int res2, char base1, char base2, double &rmsd, const int &tmpidx)
{
    int j;
    int count;
    vector<vector<double>> RotMatix;
    vector<double> TranVect;
    vector<double> xyztmp(3,0);
    int atom;
    BasicFunc bf;
    Point3s vec;

    static map<char, vector<int>> targetatoms={
	    {'A',{0,1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19,20,21}},
        {'C',{0,1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19}},
        {'G',{0,1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19,20,21,22}},
        {'U',{0,1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19}}
	};

    vector<vector<double> > xyz_list1(targetatoms[base1].size()+targetatoms[base2].size(), vector<double>(3,0));
    vector<vector<double> > xyz_list2(targetatoms[base1].size()+targetatoms[base2].size(), vector<double>(3,0));

    count = 0;
    for(j=0; j<targetatoms[base1].size();j++){
        atom=targetatoms[base1][j];
        xyz_list1[count][0]=refinestr_temp[res1].ptv[atom].x;
        xyz_list1[count][1]=refinestr_temp[res1].ptv[atom].y;
        xyz_list1[count][2]=refinestr_temp[res1].ptv[atom].z;
        count+=1;
    }

    for(j=0; j<targetatoms[base2].size();j++){
        atom=targetatoms[base2][j];
        xyz_list1[count][0]=refinestr_temp[res2].ptv[atom].x;
        xyz_list1[count][1]=refinestr_temp[res2].ptv[atom].y;
        xyz_list1[count][2]=refinestr_temp[res2].ptv[atom].z;
        count+=1;
    }

    count = 0;
    for(j=0; j<targetatoms[base1].size();j++){
        atom=targetatoms[base1][j];
        //xyz_list1[count][0]=refinestr_temp[res1].ptv[atom].x;
        //xyz_list1[count][1]=refinestr_temp[res1].ptv[atom].y;
        //xyz_list1[count][2]=refinestr_temp[res1].ptv[atom].z;

        xyz_list2[count][0]=BPTEMPLATE[base1][base2][tmpidx][0].ptv[atom].x;
        xyz_list2[count][1]=BPTEMPLATE[base1][base2][tmpidx][0].ptv[atom].y;
        xyz_list2[count][2]=BPTEMPLATE[base1][base2][tmpidx][0].ptv[atom].z;

        count+=1;
    }

    for(j=0; j<targetatoms[base2].size();j++){
        atom=targetatoms[base2][j];
        //xyz_list1[count][0]=refinestr_temp[res2].ptv[atom].x;
        //xyz_list1[count][1]=refinestr_temp[res2].ptv[atom].y;
        //xyz_list1[count][2]=refinestr_temp[res2].ptv[atom].z;

        xyz_list2[count][0]=BPTEMPLATE[base1][base2][tmpidx][1].ptv[atom].x;
        xyz_list2[count][1]=BPTEMPLATE[base1][base2][tmpidx][1].ptv[atom].y;
        xyz_list2[count][2]=BPTEMPLATE[base1][base2][tmpidx][1].ptv[atom].z;

        count+=1;
    }

    bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

    rmsd=0;
    count=0;
    for(j=0; j<targetatoms[base1].size();j++){
        //atom=targetatoms[base1][j];
        bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);
        vec.x=xyz_list1[count][0]-xyztmp[0];
        vec.y=xyz_list1[count][1]-xyztmp[1];
        vec.z=xyz_list1[count][2]-xyztmp[2];

        rmsd+=(vec.x*vec.x+vec.y*vec.y+vec.z*vec.z);
        count+=1;
    }

    for(j=0; j<targetatoms[base2].size();j++){
        //atom=targetatoms[base2][j];
        bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);

        vec.x=xyz_list1[count][0]-xyztmp[0];
        vec.y=xyz_list1[count][1]-xyztmp[1];
        vec.z=xyz_list1[count][2]-xyztmp[2];

        rmsd+=(vec.x*vec.x+vec.y*vec.y+vec.z*vec.z);
        count+=1;
    }
    rmsd = sqrt(rmsd/count);
}

void MonteCarloSimulation::MCselectbpfrag(vector<vector<Point3s>> &gradient, int res1, int res2, double K, list<BPdistobj> &bpdistobj)
{
	int i,j;
	int atom, atom1, atom2;
	int count;
	vector<vector<double>> RotMatix;
    vector<double> TranVect;
	double rmsd;
	vector<double> xyztmp(3,0);
	Point3s vec;
	BasicFunc bf;
    char base1, base2;
    int template_idx=-1;
    double minrmsd=1e8;

    base1 = refinestr_temp[res1].aaa;
    base2 = refinestr_temp[res2].aaa;

    /*
    static map<char, vector<int>> targetatoms_dist={
	    {'A',{12, 13, 14, 15, 16, 17, 18, 19, 21}},
        {'C',{12, 13, 14, 15, 16, 17}},
        {'G',{12, 13, 14, 15, 16, 17, 18, 19, 20}},
        {'U',{12, 13, 14, 15, 16, 17}}
    };
    */

    static map<char, vector<int>> targetatoms_dist={
	    {'A',{12, 13, 14}},
        {'C',{12, 13, 14}},
        {'G',{12, 13, 14}},
        {'U',{12, 13, 14}}
    };

    vector<vector<Residueatoms>> &bptemplate=BPTEMPLATE[base1][base2];
    vector<vector<vector<double>>> &bpdistancemap=BPDISTMAP[base1][base2];

    vector<double> rmsds;
    rmsds.resize(bptemplate.size(), 0.0);

	for(i=0;i<bptemplate.size();i++)
	{
        rmsd=0;
		count=0;
        MCbp_superimpose(res1, res2, base1, base2, rmsd, i);
        //rmsds[i]=rmsd;

        if(minrmsd > rmsd){
            minrmsd = rmsd;
            template_idx = i;
        }
	}

    if(template_idx == -1) {
        return;
    }

    //indices[bpidx].resize(rmsds.size());
    //for (int i = 0; i < rmsds.size(); ++i) indices[bpidx][i] = i;

    //sort(indices[bpidx].begin(), indices[bpidx].end(), [&rmsds](int i1, int i2) {
        //return rmsds[i1] < rmsds[i2];
    //});

    //int template_idx = indices[bpidx][0];

    for(int i=0;i<targetatoms_dist[base1].size();i++){
        atom1=targetatoms_dist[base1][i];
        for(int j=0;j<targetatoms_dist[base2].size();j++){
            atom2=targetatoms_dist[base2][j];
            
            //double dx = refinestr_temp[res1].ptv[atom1].x - refinestr_temp[res2].ptv[atom2].x;
            //double dy = refinestr_temp[res1].ptv[atom1].y - refinestr_temp[res2].ptv[atom2].y;
            //double dz = refinestr_temp[res1].ptv[atom1].z - refinestr_temp[res2].ptv[atom2].z;
            //double length = sqrt(dx*dx + dy*dy+ dz*dz);
            BPdistobj bpdistobjone(&refinestr_temp[res1].ptv[atom1], &refinestr_temp[res2].ptv[atom2], &gradient[res1][atom1], &gradient[res2][atom2], bpdistancemap[template_idx][atom1][atom2], K);
            bpdistobj.push_back(bpdistobjone);
        }
    }
}

void MonteCarloSimulation::MCsetbpdistobj(vector<vector<Point3s>> &gradient, const vector<vector<int>> &basepairlist, const vector<double> &basepairscorelist, double K, list<BPdistobj> &bpdistobj)
{
	int i,j,k;
	int res1, res2;
	char base1, base2;
    int minfragidx; 

    //int bpnum=0;
    //for(int i=0;i<basepairlist.size();i++) bpnum+=1;
    //indices.resize(bpnum);

    for(int i=0;i<basepairlist.size();i++){
        res1=basepairlist[i][0];
        res2=basepairlist[i][1];

        base1 = refinestr_temp[res1].aaa;
        base2 = refinestr_temp[res2].aaa;

        MCselectbpfrag(gradient, res1, res2, K*basepairscorelist[i], bpdistobj);
    }
}


void MonteCarloSimulation::MCnt_superimpose(int res, char base, double &rmsd, const int &tmpidx)
{
    int i,j;
    int count;
    vector<vector<double>> RotMatix;
    vector<double> TranVect;
    vector<double> xyztmp(3,0);
    int atom;
    BasicFunc bf;
    Point3s vec;

    map<char, vector<int>> targetatoms={
        {'A',{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21}},
        {'C',{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}},
        {'G',{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22}},
        {'U',{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}}
    };

    vector<vector<double> > xyz_list1(targetatoms[base].size(), vector<double>(3,0));
    vector<vector<double> > xyz_list2(targetatoms[base].size(), vector<double>(3,0));

    count = 0;
    for(j=0; j<targetatoms[base].size();j++){
        atom=targetatoms[base][j];
        xyz_list1[count][0]=refinestr_temp[res].ptv[atom].x;
        xyz_list1[count][1]=refinestr_temp[res].ptv[atom].y;
        xyz_list1[count][2]=refinestr_temp[res].ptv[atom].z;

        xyz_list2[count][0]=NTTEMPLATE[base][tmpidx][0].ptv[atom].x;
        xyz_list2[count][1]=NTTEMPLATE[base][tmpidx][0].ptv[atom].y;
        xyz_list2[count][2]=NTTEMPLATE[base][tmpidx][0].ptv[atom].z;

        count+=1;
    }

    bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

    rmsd=0;
    count=0;
    for(j=0; j<targetatoms[base].size();j++){
        atom=targetatoms[base][j];
        bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);
        vec.x=refinestr_temp[res].ptv[atom].x-xyztmp[0];
        vec.y=refinestr_temp[res].ptv[atom].y-xyztmp[1];
        vec.z=refinestr_temp[res].ptv[atom].z-xyztmp[2];

        rmsd+=(vec.x*vec.x+vec.y*vec.y+vec.z*vec.z);
        count+=1;
    }
    rmsd = sqrt(rmsd/count);
}


void MonteCarloSimulation::MCselectntfrag(vector<vector<Point3s>> &gradient, int res, int ntidx, vector<vector<int>> &indices, double K, list<SingleNtTemplateobjMC> &nttempobj)
{
	int i,j;
	int atom;
	int count;
	vector<vector<double>> RotMatix;
    vector<double> TranVect;
	double rmsd;
	vector<double> xyztmp(3,0);
	Point3s vec;
	BasicFunc bf;
    char base;

    base = refinestr_temp[res].aaa;

    map<string, vector<int>> targetatoms={
        {"A",{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21}},
        {"C",{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}},
        {"G",{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22}},
        {"U",{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}}
    };

    vector<vector<Residueatoms>> nttemplate=NTTEMPLATE[base];

    vector<double> rmsds;
    rmsds.resize(nttemplate.size(), 0.0);

	for(i=0;i<nttemplate.size();i++)
	{
        rmsd=0;
		count=0;
        MCnt_superimpose(res, base, rmsd, i);
        rmsds[i]=rmsd;
	}

    indices[ntidx].resize(rmsds.size());
    for (int i = 0; i < rmsds.size(); ++i) indices[ntidx][i] = i;

    sort(indices[ntidx].begin(), indices[ntidx].end(), [&rmsds](int i1, int i2) {
        return rmsds[i1] < rmsds[i2];
    });

    int template_idx = indices[ntidx][0];
    SingleNtTemplateobjMC nttempobjone(refinestr_temp[res].ptv, &gradient[res], refinestr_temp[res].aaa, &NTTEMPLATE[base], &indices[ntidx], 1, K);
    nttempobj.push_back(nttempobjone);
}


void MonteCarloSimulation::MCsetnttemplateobj(vector<vector<Point3s>> &gradient, vector<vector<int>> &indices, double K, list<SingleNtTemplateobjMC> &nttempobj)
{
	int i,j,k;
	char base;
    int minfragidx; 

    indices.resize(refinestr_temp.size());
    for(int i=0;i<refinestr_temp.size();i++){
        MCselectntfrag(gradient, i, i, indices, K, nttempobj);
    }
}

void MonteCarloSimulation::MCsetbaseplaneobj(vector<vector<Point3s>> &gradient, vector<vector<Point3s*>> &BaseAtomPtrslist, vector<vector<Point3s*>> &BaseAtomgradlist, const vector<bool> &moveresidues, const double &K, list<BasePairobj> &bplaneobj)
{
	int i,j,k;
	int res1, res2;
	int atom1, atom2;
	string aat1, aat2;
    char base;

    static map<char, vector<int>> targetatoms={
		{'A',{12, 13, 14, 15, 16, 17, 18, 19, 20, 21}},
        {'C',{12, 13, 14, 15, 16, 17, 18, 19}},
        {'G',{12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22}},
        {'U',{12, 13, 14, 15, 16, 17, 18, 19}}
	};

    int ntnum = refinestr_temp.size();

    BaseAtomPtrslist.resize(ntnum);
    BaseAtomgradlist.resize(ntnum);

    for(int i=0;i<refinestr_temp.size();i++){

        if(moveresidues[i]){

            base = refinestr_temp[i].aaa;

            for(j=0;j<targetatoms[base].size();j++){
                k = targetatoms[base][j];
                Point3s* respt = &refinestr_temp[i].ptv[k];
                BaseAtomPtrslist[i].push_back(respt);

                Point3s* resgradpt = &gradient[i][k];
                BaseAtomgradlist[i].push_back(resgradpt);
            }

            BasePairobj baseobjone(BaseAtomPtrslist[i], BaseAtomgradlist[i], K);
            bplaneobj.push_back(baseobjone);
        }
    }
}


void MonteCarloSimulation::MCsetbasepairobj(vector<vector<Point3s>> &gradient, const vector<vector<int>> &basepairlist, const vector<vector<int>> &noncbasepairlist, const vector<double> &basepairscorelist, const vector<double> &noncbasepairscorelist, vector<vector<Point3s*>> &AtomPtrslist, vector<vector<Point3s*>> &Atomgradlist, const double &K, list<BasePairobj> &bpobj)
{
	int i,j,k;
	int res1, res2;
	int atom1, atom2;
	string aat1, aat2;
    char base1, base2;

    static map<char, vector<int>> targetatoms={
		{'A',{12, 13, 14, 15, 16, 17, 18, 19, 20, 21}},
        {'C',{12, 13, 14, 15, 16, 17, 18, 19}},
        {'G',{12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22}},
        {'U',{12, 13, 14, 15, 16, 17, 18, 19}}
	};

    /*
    static map<char, vector<int>> targetatoms={
		{'A',{0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21}},
        {'C',{0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 15, 16, 17, 18, 19}},
        {'G',{0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22}},
        {'U',{0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 15, 16, 17, 18, 19}}
	};
    */

    int bpnum=0;
    //for(int i=0;i<basepairlist.size();i++) bpnum+=1;
    bpnum+=basepairlist.size();
    bpnum+=noncbasepairlist.size();

    AtomPtrslist.resize(bpnum);
    Atomgradlist.resize(bpnum);

    int count=0;
    for(int i=0;i<basepairlist.size();i++)
    {
        res1=basepairlist[i][0];
        res2=basepairlist[i][1];

        base1 = refinestr_temp[res1].aaa;
        base2 = refinestr_temp[res2].aaa;

        //if((base1 == 'A' && base2 == 'U') || (base1 == 'U' && base2 == 'A') || (base1 == 'C' && base2 == 'G') || (base1 == 'G' && base2 == 'C'))
        //{
            //vector<Point3s*> AtomPtrs;
            //vector<Point3s*> Atomgrad;

            //for(k=12;k<AccAtomParamNum[base1][3];k++){
        for(j=0;j<targetatoms[base1].size();j++){
            k = targetatoms[base1][j];
            Point3s* respt = &refinestr_temp[res1].ptv[k];
            AtomPtrslist[count].push_back(respt);

            Point3s* resgradpt = &gradient[res1][k];
            Atomgradlist[count].push_back(resgradpt);
        }

            //for(k=12;k<AccAtomParamNum[base2][3];k++){
        for(j=0;j<targetatoms[base2].size();j++){
            k = targetatoms[base2][j];
            Point3s* respt = &refinestr_temp[res2].ptv[k];
            AtomPtrslist[count].push_back(respt);

            Point3s* resgradpt = &gradient[res2][k];
            Atomgradlist[count].push_back(resgradpt);
        }

        //cout<<res1 << "-" << res2 << "-" <<  base1 << "-" << base2 << "-" << basepairscorelist[i] * K<<endl;

        BasePairobj bpobjone(AtomPtrslist[count], Atomgradlist[count], basepairscorelist[i] * K);
        bpobj.push_back(bpobjone);

        count+=1;
        //}
    }

    /*
    for(int i=0;i<noncbasepairlist.size();i++){
        res1=noncbasepairlist[i][0];
        res2=noncbasepairlist[i][1];

        base1 = refinestr_temp[res1].aaa;
        base2 = refinestr_temp[res2].aaa;

        for(j=0;j<targetatoms[base1].size();j++){
            k = targetatoms[base1][j];
            Point3s* respt = &refinestr_temp[res1].ptv[k];
            AtomPtrslist[count].push_back(respt);

            Point3s* resgradpt = &gradient[res1][k];
            Atomgradlist[count].push_back(resgradpt);
        }

        //for(k=12;k<AccAtomParamNum[base2][3];k++){
        for(j=0;j<targetatoms[base2].size();j++){
            k = targetatoms[base2][j];
            Point3s* respt = &refinestr_temp[res2].ptv[k];
            AtomPtrslist[count].push_back(respt);

            Point3s* resgradpt = &gradient[res2][k];
            Atomgradlist[count].push_back(resgradpt);
        }

        //cout<<base1 << "-" << base2 << "-" << noncbasepairscorelist[i]<<endl;

        BasePairobj bpobjone(AtomPtrslist[count], Atomgradlist[count], noncbasepairscorelist[i] * K);
        bpobj.push_back(bpobjone);

        count+=1;
    }
    */

}

void MonteCarloSimulation::MCsetNonbondedobj(vector<vector<Point3s>> &gradient, const vector<vector<int>> &contactmap, const vector<vector<int>> &bindingmap, const vector<vector<int>> &cbasepairlist, const vector<vector<int>> &noncbasepairlist, const vector<double> &cbasepairscorelist, const vector<double> &noncbasepairscorelist,list<TElectrobj> &telecteobj, list<LennardsJobj> &ljobj, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, list<Clashobj> &clashobj, map<string, double> &enecoeff)
{
    int i, j, k;
    char base1, base2;
    string base1_str, base2_str, base;
    double dQQ, dQQep;
    Point3d p12;
    char AccBase, DonBase;
    int atomnum1, atomnum2, atomnum;
    int res1, res2;
    double length_thresh;
    bool flag;
    vector<vector<vector<int>>> hbondinfo;
    vector<string> nts={"A","C","G","U"};
    map<string, int> atomnumnts={{"A", 33}, {"C", 31}, {"G", 34}, {"U", 30}};

    for(i=0;i<numbb;i++)
	{
		base1=refinestr_temp[i].aaa;
		base1_str=string(1, base1);
        atomnum1=AccAtomParamNum[base1][4];
		for(j=i;j<numbb;j++)
		{
			base2=refinestr_temp[j].aaa;
			base2_str=string(1, base2);
            atomnum2=AccAtomParamNum[base2][4];
            
            if(i != j)
            {
                /*
                if(bpmap[i][j] == 2)
                {
                    if((base1_str == "G" && base2_str == "C") || (base1_str == "C" && base2_str == "G") || (base1_str == "A" && base2_str == "U") || (base1_str == "U" && base2_str == "A"))
                    {
                        vector<vector<vector<int>>> hbondinfo = BpHbondInfo[base1_str+base2_str];
                        for(k=0;k<hbondinfo.size();k++)
                        {
                            if(hbondinfo[k].size() > 1)
                            {
                                if(hbondinfo[k][0][3] == 0){
                                    const Point3s &ptH1 = refinestr_temp[j].ptv[hbondinfo[k][0][1]];
                                    const Point3s &ptH2 = refinestr_temp[j].ptv[hbondinfo[k][1][1]];
                                    const Point3s &ptAcc = refinestr_temp[i].ptv[hbondinfo[k][0][0]];
                                    const Point3s &ptDon = refinestr_temp[j].ptv[hbondinfo[k][0][2]];

                                    AccBase = base1; DonBase = base2;
                                    Hbond2hydobj hbondobjone = Hbond2hydobj::create(&refinestr_temp[j].ptv[hbondinfo[k][0][1]], &refinestr_temp[j].ptv[hbondinfo[k][1][1]], &refinestr_temp[j].ptv[hbondinfo[k][0][2]], &refinestr_temp[i].ptv[hbondinfo[k][0][0]], &gradient[j][hbondinfo[k][0][1]], &gradient[j][hbondinfo[k][1][1]], &gradient[j][hbondinfo[k][0][2]], &gradient[i][hbondinfo[k][0][0]], BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean, enecoeff["HBond"]);
                                    hbondbpobj.push_back(hbondobjone);
                                }else{
                                    const Point3s &ptH1 = refinestr_temp[i].ptv[hbondinfo[k][0][1]];
                                    const Point3s &ptH2 = refinestr_temp[i].ptv[hbondinfo[k][1][1]];
                                    const Point3s &ptAcc = refinestr_temp[j].ptv[hbondinfo[k][0][0]];
                                    const Point3s &ptDon = refinestr_temp[i].ptv[hbondinfo[k][0][2]];

                                    AccBase = base2; DonBase = base1;
                                    Hbond2hydobj hbondobjone = Hbond2hydobj::create(&refinestr_temp[i].ptv[hbondinfo[k][0][1]], &refinestr_temp[i].ptv[hbondinfo[k][1][1]], &refinestr_temp[i].ptv[hbondinfo[k][0][2]], &refinestr_temp[j].ptv[hbondinfo[k][0][0]], &gradient[i][hbondinfo[k][0][1]], &gradient[i][hbondinfo[k][1][1]], &gradient[i][hbondinfo[k][0][2]], &gradient[j][hbondinfo[k][0][0]], BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean, enecoeff["HBond"]);
                                    hbondbpobj.push_back(hbondobjone);
                                }
                            }else{
                                if(hbondinfo[k][0][3] == 0){
                                    const Point3s &ptH = refinestr_temp[j].ptv[hbondinfo[k][0][1]];
                                    const Point3s &ptAcc = refinestr_temp[i].ptv[hbondinfo[k][0][0]];
                                    const Point3s &ptDon = refinestr_temp[j].ptv[hbondinfo[k][0][2]];

                                    AccBase = base1; DonBase = base2;
                                    Hbondobj hbondobjone = Hbondobj::create(&refinestr_temp[j].ptv[hbondinfo[k][0][1]], &refinestr_temp[j].ptv[hbondinfo[k][0][2]], &refinestr_temp[i].ptv[hbondinfo[k][0][0]], &gradient[j][hbondinfo[k][0][1]], &gradient[j][hbondinfo[k][0][2]], &gradient[i][hbondinfo[k][0][0]], BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean, enecoeff["HBond"]);
                                    hbondobj.push_back(hbondobjone);
                                }else{
                                    const Point3s &ptH = refinestr_temp[i].ptv[hbondinfo[k][0][1]];
                                    const Point3s &ptAcc = refinestr_temp[j].ptv[hbondinfo[k][0][0]];
                                    const Point3s &ptDon = refinestr_temp[i].ptv[hbondinfo[k][0][2]];

                                    AccBase = base2; DonBase = base1;
                                    Hbondobj hbondobjone = Hbondobj::create(&refinestr_temp[i].ptv[hbondinfo[k][0][1]], &refinestr_temp[i].ptv[hbondinfo[k][0][2]], &refinestr_temp[j].ptv[hbondinfo[k][0][0]], &gradient[i][hbondinfo[k][0][1]], &gradient[i][hbondinfo[k][0][2]], &gradient[j][hbondinfo[k][0][0]], BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean, enecoeff["HBond"]);
                                    hbondobj.push_back(hbondobjone);
                                }
                            }
                        }
                    }
                }
                */

                hbondinfo = BBHbondInfo[base1_str+base2_str];

                k = 0;
                AccBase = base1; DonBase = base2;
                length_thresh = BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean + HBOND_CONST*BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_std;

                //cout<<hbondinfo[k][0][0] << "-" << hbondinfo[k][0][1] << "-" << hbondinfo[k][0][2]<<endl;
                //cout << BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean << "-" << length_thresh<<endl;

                flag = checkhbond(&refinestr_temp[j].ptv[hbondinfo[k][0][1]], &refinestr_temp[j].ptv[hbondinfo[k][0][2]], &refinestr_temp[i].ptv[hbondinfo[k][0][0]], length_thresh);

                if(flag){
                    Hbondobj hbondobjone = Hbondobj::create(&refinestr_temp[j].ptv[hbondinfo[k][0][1]], &refinestr_temp[j].ptv[hbondinfo[k][0][2]], &refinestr_temp[i].ptv[hbondinfo[k][0][0]], &gradient[j][hbondinfo[k][0][1]], &gradient[j][hbondinfo[k][0][2]], &gradient[i][hbondinfo[k][0][0]], BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean - BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_std, enecoeff["BBHBond"]);
                    hbondobj.push_back(hbondobjone);
                }

                for(k=1;k<hbondinfo.size();k++)
                {
                    AccBase = base1; DonBase = base2;

                    length_thresh = BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean + HBOND_CONST*BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_std;

                    //cout<<hbondinfo[k][0][0] << "-" << hbondinfo[k][0][1] << "-" << hbondinfo[k][0][2]<<endl;
                    //cout << BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean << "-" << length_thresh<<endl;

                    flag = checkhbond(&refinestr_temp[j].ptv[hbondinfo[k][0][1]], &refinestr_temp[j].ptv[hbondinfo[k][0][2]], &refinestr_temp[i].ptv[hbondinfo[k][0][0]], length_thresh);
					if(!flag) continue;

                    Hbondobj hbondobjone = Hbondobj::create(&refinestr_temp[j].ptv[hbondinfo[k][0][1]], &refinestr_temp[j].ptv[hbondinfo[k][0][2]], &refinestr_temp[i].ptv[hbondinfo[k][0][0]], &gradient[j][hbondinfo[k][0][1]], &gradient[j][hbondinfo[k][0][2]], &gradient[i][hbondinfo[k][0][0]], BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean - BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_std, enecoeff["BBBaseHBond"]);
                    hbondobj.push_back(hbondobjone);
                }

                hbondinfo = BBHbondInfo[base2_str+base1_str];

                k = 0;
                AccBase = base2; DonBase = base1;
                length_thresh = BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean + HBOND_CONST*BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_std;

                //cout<<hbondinfo[k][0][0] << "-" << hbondinfo[k][0][1] << "-" << hbondinfo[k][0][2]<<endl;
                //cout<< BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean << "-" << length_thresh<<endl;

                flag = checkhbond(&refinestr_temp[i].ptv[hbondinfo[k][0][1]], &refinestr_temp[i].ptv[hbondinfo[k][0][2]], &refinestr_temp[j].ptv[hbondinfo[k][0][0]], length_thresh);

                if(flag){
                    Hbondobj hbondobjone = Hbondobj::create(&refinestr_temp[i].ptv[hbondinfo[k][0][1]], &refinestr_temp[i].ptv[hbondinfo[k][0][2]], &refinestr_temp[j].ptv[hbondinfo[k][0][0]], &gradient[i][hbondinfo[k][0][1]], &gradient[i][hbondinfo[k][0][2]], &gradient[j][hbondinfo[k][0][0]], BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean - BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_std, enecoeff["BBHBond"]);
                    hbondobj.push_back(hbondobjone);
                }

                for(k=1;k<hbondinfo.size();k++)
                {
                    AccBase = base2; DonBase = base1;

                    length_thresh = BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean + HBOND_CONST*BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_std;

                    //cout<<hbondinfo[k][0][0] << "-" << hbondinfo[k][0][1] << "-" << hbondinfo[k][0][2]<<endl;
                    //cout<< BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean << "-" << length_thresh<<endl;

                    flag = checkhbond(&refinestr_temp[i].ptv[hbondinfo[k][0][1]], &refinestr_temp[i].ptv[hbondinfo[k][0][2]], &refinestr_temp[j].ptv[hbondinfo[k][0][0]], length_thresh);
					if(!flag) continue;

                    Hbondobj hbondobjone = Hbondobj::create(&refinestr_temp[i].ptv[hbondinfo[k][0][1]], &refinestr_temp[i].ptv[hbondinfo[k][0][2]], &refinestr_temp[j].ptv[hbondinfo[k][0][0]], &gradient[i][hbondinfo[k][0][1]], &gradient[i][hbondinfo[k][0][2]], &gradient[j][hbondinfo[k][0][0]], BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean - BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_std, enecoeff["BBBaseHBond"]);
                    hbondobj.push_back(hbondobjone);
                }

                if(contactmap[i][j] >= 1)
                {
                    for(int ii=0;ii<atomid2mcidx[base1].size();ii++)
                    {
                        for(int jj=0;jj<atomid2mcidx[base2].size();jj++)
                        {  
                            if((bindingmap[i][j] != 1) || (bindingmap[i][j] == 1 && intramolConnect[base1_str][base2_str][ii][jj] >= 2)){
                                if(bindingmap[i][j] == 1 && intramolConnect[base1_str][base2_str][ii][jj] == 3){
                                    dQQ = (COULOMB*AtomFeatures[base1_str][ii][0]*AtomFeatures[base2_str][jj][0])/1.2;
                                }else{
                                    dQQ = COULOMB*AtomFeatures[base1_str][ii][0]*AtomFeatures[base2_str][jj][0];
                                }
                                dQQep = dQQ/EpsilonRNA;

                                TElectrobj telecteobjone = TElectrobj::create(&refinestr_temp[i].ptv[ii], &refinestr_temp[j].ptv[jj], &gradient[i][ii], &gradient[j][jj], dQQep);
                                telecteobj.push_back(telecteobjone);
                            
                                if(contactmap[i][j] >= 2)
                                {
                                    double vdWInfodR01 = AtomFeatures.at(base1_str)[ii][1];
                                    double vdWInfodR02 = AtomFeatures.at(base2_str)[jj][1];
                                    double vdWInfodEpsilon1 = AtomFeatures.at(base1_str)[ii][2];
                                    double vdWInfodEpsilon2 = AtomFeatures.at(base2_str)[jj][2];
                                    double star = (vdWInfodR01 + vdWInfodR02);
                                    double epsilon =  sqrt(vdWInfodEpsilon1*vdWInfodEpsilon2);
                                    double LJ6 = -2*epsilon*pow(star,6);
                                    double LJ12 = epsilon*pow(star,12);

                                    if(bindingmap[i][j] == 1 && intramolConnect[base1_str][base2_str][ii][jj] == 3){
                                        LJ6 /=2.0;
                                        LJ12 /=2.0;
                                    }

                                    double Kval = 1.0;
                                    LennardsJobj ljobjone = LennardsJobj::create(&refinestr_temp[i].ptv[ii], &refinestr_temp[j].ptv[jj], &gradient[i][ii], &gradient[j][jj], LJ6, LJ12, Kval);
                                    ljobj.push_back(ljobjone);

                                    double vdW = (vdWInfodR01 + vdWInfodR02);
                                    //double vdW = MCVWPARAM[base1][ii]+MCVWPARAM[base2][jj]-0.4;

                                    //cout<<ii << "-" << jj << "-" << vdW<<endl;
                                    Clashobj clashobjone = Clashobj::create(&refinestr_temp[i].ptv[ii], &refinestr_temp[j].ptv[jj], &gradient[i][ii], &gradient[j][jj], vdW, enecoeff["Clash"]);
                                    clashobj.push_back(clashobjone);
                                }
                            }

                            /*
                            if((bindingmap[i][j] != 1) || (bindingmap[i][j] == 1 && intramolConnect[base1_str][base2_str][ii][jj] == 4)){
                                if((HBondHydrogenAtoms[base1_str][ii] == 1) && (HBondAcceptorAtoms[base2_str][jj] == 1))
                                {
                                    for(int m=0;m<atomid2mcidx[base1].size();m++){
                                        if(intraAtomConnect[base1_str][ii][m]){
                                            const Point3s &ptH = refinestr_temp[i].ptv[ii];
                                            const Point3s &ptAcc = refinestr_temp[j].ptv[jj];
                                            const Point3s &ptDon = refinestr_temp[i].ptv[m];
        
                                            p12=bf.setv(ptAcc.x-ptDon.x,ptAcc.y-ptDon.y,ptAcc.z-ptDon.z);
                                            double hleng=bf.norm(p12);
        
                                            //double energy=calcenergyhydrogenbond(ptH, ptDon, ptAcc);
                                            //if(energy < HBONDThreshold){
                                            if(hleng < HBONDADLengThreshold){
                                                Hbondobj hbondobjone = Hbondobj::create(&refinestr_temp[i].ptv[ii], &refinestr_temp[i].ptv[m], &refinestr_temp[j].ptv[jj], &gradient[i][ii], &gradient[i][m], &gradient[j][jj], enecoeff["HBond"]);
                                                hbondobj.push_back(hbondobjone);
                                            }
                                        }
                                    }
                                }
        
                                
                                if((HBondHydrogenAtoms[base2_str][jj] == 1) && (HBondAcceptorAtoms[base1_str][ii] == 1))
                                {
                                    for(int m=0;m<atomid2mcidx[base2].size();m++){
                                        if(intraAtomConnect[base2_str][jj][m]){
                                            const Point3s &ptH = refinestr_temp[j].ptv[jj];
                                            const Point3s &ptAcc = refinestr_temp[i].ptv[ii];
                                            const Point3s &ptDon = refinestr_temp[j].ptv[m];
        
                                            //double energy=calcenergyhydrogenbond(ptH, ptDon, ptAcc);
                                            p12=bf.setv(ptAcc.x-ptDon.x,ptAcc.y-ptDon.y,ptAcc.z-ptDon.z);
                                            double hleng=bf.norm(p12);
        
                                            //if(energy < HBONDThreshold){
                                            if(hleng < HBONDADLengThreshold){
                                                Hbondobj hbondobjone = Hbondobj::create(&refinestr_temp[j].ptv[jj], &refinestr_temp[j].ptv[m], &refinestr_temp[i].ptv[ii], &gradient[j][jj], &gradient[j][m], &gradient[i][ii], enecoeff["HBond"]);
                                                hbondobj.push_back(hbondobjone);
                                            }
                                        }
                                    }
                                }
                            }
                            */
                        }
                    }
                }
            }else{
                hbondinfo = BBHbondInfo[base1_str+base2_str];
                for(k=1;k<hbondinfo.size();k++)
                {
                    AccBase = base1; DonBase = base2;

                    length_thresh = BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean + HBOND_CONST*BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_std;

                    //cout<<hbondinfo[k][0][0] << "-" << hbondinfo[k][0][1] << "-" << hbondinfo[k][0][2]<<endl;
                    //cout << BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean << "-" << length_thresh<<endl;

                    flag = checkhbond(&refinestr_temp[j].ptv[hbondinfo[k][0][1]], &refinestr_temp[j].ptv[hbondinfo[k][0][2]], &refinestr_temp[i].ptv[hbondinfo[k][0][0]], length_thresh);
					if(!flag) continue;

                    Hbondobj hbondobjone = Hbondobj::create(&refinestr_temp[j].ptv[hbondinfo[k][0][1]], &refinestr_temp[j].ptv[hbondinfo[k][0][2]], &refinestr_temp[i].ptv[hbondinfo[k][0][0]], &gradient[j][hbondinfo[k][0][1]], &gradient[j][hbondinfo[k][0][2]], &gradient[i][hbondinfo[k][0][0]], BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean - BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_std, enecoeff["BBBaseHBond"]);
                    hbondobj.push_back(hbondobjone);
                }

                hbondinfo = BBHbondInfo[base2_str+base1_str];
                for(k=1;k<hbondinfo.size();k++)
                {
                    AccBase = base2; DonBase = base1;

                    length_thresh = BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean + HBOND_CONST*BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_std;

                    //cout<<hbondinfo[k][0][0] << "-" << hbondinfo[k][0][1] << "-" << hbondinfo[k][0][2]<<endl;
                    //cout<< BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean << "-" << length_thresh<<endl;

                    flag = checkhbond(&refinestr_temp[i].ptv[hbondinfo[k][0][1]], &refinestr_temp[i].ptv[hbondinfo[k][0][2]], &refinestr_temp[j].ptv[hbondinfo[k][0][0]], length_thresh);
					if(!flag) continue;

                    Hbondobj hbondobjone = Hbondobj::create(&refinestr_temp[i].ptv[hbondinfo[k][0][1]], &refinestr_temp[i].ptv[hbondinfo[k][0][2]], &refinestr_temp[j].ptv[hbondinfo[k][0][0]], &gradient[i][hbondinfo[k][0][1]], &gradient[i][hbondinfo[k][0][2]], &gradient[j][hbondinfo[k][0][0]], BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean - BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_std, enecoeff["BBBaseHBond"]);
                    hbondobj.push_back(hbondobjone);
                }

                if(contactmap[i][j] >= 1)
                {
                    for(int ii=0;ii<atomid2mcidx[base1].size();ii++)
                    {
                        for(int jj=ii+1;jj<atomid2mcidx[base2].size();jj++)
                        {
                            if(ntmolConnect[base1_str][ii][jj] >= 3){
                                if(ntmolConnect[base1_str][ii][jj] == 3){
                                    dQQ = (COULOMB*AtomFeatures[base1_str][ii][0]*AtomFeatures[base2_str][jj][0])/1.2;
                                }else{
                                    dQQ = COULOMB*AtomFeatures[base1_str][ii][0]*AtomFeatures[base2_str][jj][0];
                                }
                                dQQep = dQQ/EpsilonRNA;

                                TElectrobj telecteobjone = TElectrobj::create(&refinestr_temp[i].ptv[ii], &refinestr_temp[j].ptv[jj], &gradient[i][ii], &gradient[j][jj], dQQep);
                                telecteobj.push_back(telecteobjone);

                                if(contactmap[i][j] >= 2)
                                {
                                    double vdWInfodR01 = AtomFeatures.at(base1_str)[ii][1];
                                    double vdWInfodR02 = AtomFeatures.at(base2_str)[jj][1];
                                    double vdWInfodEpsilon1 = AtomFeatures.at(base1_str)[ii][2];
                                    double vdWInfodEpsilon2 = AtomFeatures.at(base2_str)[jj][2];
                                    double star = (vdWInfodR01 + vdWInfodR02);
                                    double epsilon =  sqrt(vdWInfodEpsilon1*vdWInfodEpsilon2);
                                    double LJ6 = -2*epsilon*pow(star,6);
                                    double LJ12 = epsilon*pow(star,12);

                                    if(ntmolConnect[base1_str][ii][jj] == 3){
                                        LJ6 /=2.0;
                                        LJ12 /=2.0;
                                    }

                                    double Kval = 1.0;
                                    LennardsJobj ljobjone = LennardsJobj::create(&refinestr_temp[i].ptv[ii], &refinestr_temp[j].ptv[jj], &gradient[i][ii], &gradient[j][jj], LJ6, LJ12, Kval);
                                    ljobj.push_back(ljobjone);

                                    double vdW = (vdWInfodR01 + vdWInfodR02);
                                    //double vdW = MCVWPARAM[base1][ii]+MCVWPARAM[base2][jj]-0.4;
                                    Clashobj clashobjone = Clashobj::create(&refinestr_temp[i].ptv[ii], &refinestr_temp[j].ptv[jj], &gradient[i][ii], &gradient[j][jj], vdW, enecoeff["Clash"]);
                                    clashobj.push_back(clashobjone);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    for(int i=0;i<cbasepairlist.size();i++)
    {
        res1=cbasepairlist[i][0];
        res2=cbasepairlist[i][1];

        base1 = refinestr_temp[res1].aaa;
        base2 = refinestr_temp[res2].aaa;

        base1_str=string(1, base1);
        base2_str=string(1, base2);

        //cout<<"aaa"<<endl;
        //cout<<base1 << "-" << base2 <<"-" << cbasepairscorelist[i]<< endl;

        hbondinfo = BpHbondInfo[base1_str+base2_str];
        for(k=0;k<hbondinfo.size();k++)
        {
            if(hbondinfo[k].size() > 1)
            {
                if(hbondinfo[k][0][3] == 0){
                    AccBase = base1; DonBase = base2;

                    //cout<<hbondinfo[k][0][0] << "-" << hbondinfo[k][0][1] << "-" << hbondinfo[k][0][2]<<endl;

                    //cout<<BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean << "-" << BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean << "-" << BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean<<endl;
                    Hbond2hydobj hbondobjone = Hbond2hydobj::create(&refinestr_temp[res2].ptv[hbondinfo[k][0][1]], &refinestr_temp[res2].ptv[hbondinfo[k][1][1]], &refinestr_temp[res2].ptv[hbondinfo[k][0][2]], &refinestr_temp[res1].ptv[hbondinfo[k][0][0]], &gradient[res2][hbondinfo[k][0][1]], &gradient[res2][hbondinfo[k][1][1]], &gradient[res2][hbondinfo[k][0][2]], &gradient[res1][hbondinfo[k][0][0]], BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean - BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_std, cbasepairscorelist[i]*enecoeff["HBond"]);
                    hbondbpobj.push_back(hbondobjone);
                }else{
                    AccBase = base2; DonBase = base1;

                    //cout<<hbondinfo[k][0][0] << "-" << hbondinfo[k][0][1] << "-" << hbondinfo[k][0][2]<<endl;

                    //cout<<BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean << "-" << BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean << "-" << BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean<<endl;
                    Hbond2hydobj hbondobjone = Hbond2hydobj::create(&refinestr_temp[res1].ptv[hbondinfo[k][0][1]], &refinestr_temp[res1].ptv[hbondinfo[k][1][1]], &refinestr_temp[res1].ptv[hbondinfo[k][0][2]], &refinestr_temp[res2].ptv[hbondinfo[k][0][0]], &gradient[res1][hbondinfo[k][0][1]], &gradient[res1][hbondinfo[k][1][1]], &gradient[res1][hbondinfo[k][0][2]], &gradient[res2][hbondinfo[k][0][0]], BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean - BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_std, cbasepairscorelist[i]*enecoeff["HBond"]);
                    hbondbpobj.push_back(hbondobjone);
                }
            }else{
                if(hbondinfo[k][0][3] == 0){
                    AccBase = base1; DonBase = base2;

                    //cout<<hbondinfo[k][0][0] << "-" << hbondinfo[k][0][1] << "-" << hbondinfo[k][0][2]<<endl;

                    //cout<<BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean << "-" << BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean << "-" << BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean<<endl;
                    Hbondobj hbondobjone = Hbondobj::create(&refinestr_temp[res2].ptv[hbondinfo[k][0][1]], &refinestr_temp[res2].ptv[hbondinfo[k][0][2]], &refinestr_temp[res1].ptv[hbondinfo[k][0][0]], &gradient[res2][hbondinfo[k][0][1]], &gradient[res2][hbondinfo[k][0][2]], &gradient[res1][hbondinfo[k][0][0]], BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean - BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_std, cbasepairscorelist[i]*enecoeff["HBond"]);
                    hbondobj.push_back(hbondobjone);
                }else{
                    AccBase = base2; DonBase = base1;

                    //cout<<hbondinfo[k][0][0] << "-" << hbondinfo[k][0][1] << "-" << hbondinfo[k][0][2]<<endl;

                    //cout<<BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean << "-" << BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean << "-" << BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean<<endl;
                    Hbondobj hbondobjone = Hbondobj::create(&refinestr_temp[res1].ptv[hbondinfo[k][0][1]], &refinestr_temp[res1].ptv[hbondinfo[k][0][2]], &refinestr_temp[res2].ptv[hbondinfo[k][0][0]], &gradient[res1][hbondinfo[k][0][1]], &gradient[res1][hbondinfo[k][0][2]], &gradient[res2][hbondinfo[k][0][0]], BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean - BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_std, cbasepairscorelist[i]*enecoeff["HBond"]);
                    hbondobj.push_back(hbondobjone);
                }
            }
        }

    }
    

    for(int i=0;i<noncbasepairlist.size();i++){
        res1=noncbasepairlist[i][0];
        res2=noncbasepairlist[i][1];

        base1 = refinestr_temp[res1].aaa;
        base2 = refinestr_temp[res2].aaa;

        base1_str=string(1, base1);
        base2_str=string(1, base2);

        //cout<<"bbb"<<endl;
        //cout<<base1 << "-" << base2 <<"-" << noncbasepairscorelist[i] << endl;

        hbondinfo = BpHbondInfo[base1_str+base2_str];
        for(k=0;k<hbondinfo.size();k++)
        {
            if(hbondinfo[k].size() > 1)
            {
                if(hbondinfo[k][0][3] == 0){
                    AccBase = base1; DonBase = base2;

                    //cout<<BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean << "-" << BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean << "-" << BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean<<endl;
                    
                    length_thresh = BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean + HBOND_CONST*BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_std;

                    //cout<<hbondinfo[k][0][0] << "-" << hbondinfo[k][0][1] << "-" << hbondinfo[k][0][2]<<endl;
                    //cout<<BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean << "-" << length_thresh<<endl;

                    flag = checkhbond2(&refinestr_temp[res2].ptv[hbondinfo[k][0][1]], &refinestr_temp[res2].ptv[hbondinfo[k][1][1]], &refinestr_temp[res2].ptv[hbondinfo[k][0][2]], &refinestr_temp[res1].ptv[hbondinfo[k][0][0]], length_thresh);
                    
                    if(!flag) continue;

                    Hbond2hydobj hbondobjone = Hbond2hydobj::create(&refinestr_temp[res2].ptv[hbondinfo[k][0][1]], &refinestr_temp[res2].ptv[hbondinfo[k][1][1]], &refinestr_temp[res2].ptv[hbondinfo[k][0][2]], &refinestr_temp[res1].ptv[hbondinfo[k][0][0]], &gradient[res2][hbondinfo[k][0][1]], &gradient[res2][hbondinfo[k][1][1]], &gradient[res2][hbondinfo[k][0][2]], &gradient[res1][hbondinfo[k][0][0]], BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean - BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_std, noncbasepairscorelist[i]*enecoeff["HBond"]);
                    hbondbpobj.push_back(hbondobjone);
                }else{
                    AccBase = base2; DonBase = base1;

                    length_thresh = BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean + HBOND_CONST*BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_std;

                    //cout<<hbondinfo[k][0][0] << "-" << hbondinfo[k][0][1] << "-" << hbondinfo[k][0][2]<<endl;
                    //cout<<BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean << "-" << length_thresh<<endl;

                    flag = checkhbond2(&refinestr_temp[res1].ptv[hbondinfo[k][0][1]], &refinestr_temp[res1].ptv[hbondinfo[k][1][1]], &refinestr_temp[res1].ptv[hbondinfo[k][0][2]], &refinestr_temp[res2].ptv[hbondinfo[k][0][0]], length_thresh);
                    
                    if(!flag) continue;

                    //cout<<BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean << "-" << BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean << "-" << BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean<<endl;
                    Hbond2hydobj hbondobjone = Hbond2hydobj::create(&refinestr_temp[res1].ptv[hbondinfo[k][0][1]], &refinestr_temp[res1].ptv[hbondinfo[k][1][1]], &refinestr_temp[res1].ptv[hbondinfo[k][0][2]], &refinestr_temp[res2].ptv[hbondinfo[k][0][0]], &gradient[res1][hbondinfo[k][0][1]], &gradient[res1][hbondinfo[k][1][1]], &gradient[res1][hbondinfo[k][0][2]], &gradient[res2][hbondinfo[k][0][0]], BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean - BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_std, noncbasepairscorelist[i]*enecoeff["HBond"]);
                    hbondbpobj.push_back(hbondobjone);
                }
            }else{
                if(hbondinfo[k][0][3] == 0){
                    AccBase = base1; DonBase = base2;

                    length_thresh = BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean + HBOND_CONST*BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_std;

                    //cout<<hbondinfo[k][0][0] << "-" << hbondinfo[k][0][1] << "-" << hbondinfo[k][0][2]<<endl;
                    //cout<<BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean << "-" << length_thresh<<endl;

                    flag = checkhbond(&refinestr_temp[res2].ptv[hbondinfo[k][0][1]], &refinestr_temp[res2].ptv[hbondinfo[k][0][2]], &refinestr_temp[res1].ptv[hbondinfo[k][0][0]], length_thresh);

					if(!flag) continue;

                    //cout<<BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean << "-" << BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean << "-" << BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean<<endl;
                    Hbondobj hbondobjone = Hbondobj::create(&refinestr_temp[res2].ptv[hbondinfo[k][0][1]], &refinestr_temp[res2].ptv[hbondinfo[k][0][2]], &refinestr_temp[res1].ptv[hbondinfo[k][0][0]], &gradient[res2][hbondinfo[k][0][1]], &gradient[res2][hbondinfo[k][0][2]], &gradient[res1][hbondinfo[k][0][0]], BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean - BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_std, noncbasepairscorelist[i]*enecoeff["HBond"]);
                    hbondobj.push_back(hbondobjone);
                }else{
                    AccBase = base2; DonBase = base1;

                    length_thresh = BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean + HBOND_CONST*BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_std;

                    //cout<<hbondinfo[k][0][0] << "-" << hbondinfo[k][0][1] << "-" << hbondinfo[k][0][2]<<endl;
                    //cout<<BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean << "-" << length_thresh<<endl;

                    flag = checkhbond(&refinestr_temp[res1].ptv[hbondinfo[k][0][1]], &refinestr_temp[res1].ptv[hbondinfo[k][0][2]], &refinestr_temp[res2].ptv[hbondinfo[k][0][0]], length_thresh);

					if(!flag) continue;

                    //cout<<BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean << "-" << BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean << "-" << BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean<<endl;
                    Hbondobj hbondobjone = Hbondobj::create(&refinestr_temp[res1].ptv[hbondinfo[k][0][1]], &refinestr_temp[res1].ptv[hbondinfo[k][0][2]], &refinestr_temp[res2].ptv[hbondinfo[k][0][0]], &gradient[res1][hbondinfo[k][0][1]], &gradient[res1][hbondinfo[k][0][2]], &gradient[res2][hbondinfo[k][0][0]], BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean - BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_std, noncbasepairscorelist[i]*enecoeff["HBond"]);
                    hbondobj.push_back(hbondobjone);
                }
            }
        }
    }

}


void MonteCarloSimulation::MCresetNonbondedobj(vector<vector<Point3s>> &gradient, const vector<vector<int>> &contactmap, 
    const vector<vector<int>> &bindingmap, const vector<vector<int>> &cbasepairlist, const vector<vector<int>> &noncbasepairlist, 
    const vector<double> &cbasepairscorelist, const vector<double> &noncbasepairscorelist,  
    list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, map<string, double> &enecoeff)
{
    hbondobj.clear();
    hbondbpobj.clear();

    int i, j, k;
    char base1, base2;
    string base1_str, base2_str, base;
    double dQQ, dQQep;
    Point3d p12;
    char AccBase, DonBase;
    int atomnum1, atomnum2, atomnum;
    int res1, res2;
    double length_thresh;
    bool flag;
    vector<vector<vector<int>>> hbondinfo;
    vector<string> nts={"A","C","G","U"};
    map<string, int> atomnumnts={{"A", 33}, {"C", 31}, {"G", 34}, {"U", 30}};

    for(i=0;i<numbb;i++)
	{
		base1=refinestr_temp[i].aaa;
		base1_str=string(1, base1);
        atomnum1=AccAtomParamNum[base1][4];
		for(j=i;j<numbb;j++)
		{
			base2=refinestr_temp[j].aaa;
			base2_str=string(1, base2);
            atomnum2=AccAtomParamNum[base2][4];
            
            if(i != j)
            {
                hbondinfo = BBHbondInfo[base1_str+base2_str];

                k = 0;
                AccBase = base1; DonBase = base2;
                length_thresh = BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean + HBOND_CONST*BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_std;

                //cout<<hbondinfo[k][0][0] << "-" << hbondinfo[k][0][1] << "-" << hbondinfo[k][0][2]<<endl;
                //cout << BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean << "-" << length_thresh<<endl;

                flag = checkhbond(&refinestr_temp[j].ptv[hbondinfo[k][0][1]], &refinestr_temp[j].ptv[hbondinfo[k][0][2]], &refinestr_temp[i].ptv[hbondinfo[k][0][0]], length_thresh);

                if(enecoeff["BBHBond"] > 0){
                    if(flag){
                        Hbondobj hbondobjone = Hbondobj::create(&refinestr_temp[j].ptv[hbondinfo[k][0][1]], &refinestr_temp[j].ptv[hbondinfo[k][0][2]], &refinestr_temp[i].ptv[hbondinfo[k][0][0]], &gradient[j][hbondinfo[k][0][1]], &gradient[j][hbondinfo[k][0][2]], &gradient[i][hbondinfo[k][0][0]], BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean - BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_std, enecoeff["BBHBond"]);
                        hbondobj.push_back(hbondobjone);
                    }
                }

                if(enecoeff["BBBaseHBond"] > 0){
                    for(k=1;k<hbondinfo.size();k++)
                    {
                        AccBase = base1; DonBase = base2;

                        length_thresh = BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean + HBOND_CONST*BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_std;

                        //cout<<hbondinfo[k][0][0] << "-" << hbondinfo[k][0][1] << "-" << hbondinfo[k][0][2]<<endl;
                        //cout << BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean << "-" << length_thresh<<endl;

                        flag = checkhbond(&refinestr_temp[j].ptv[hbondinfo[k][0][1]], &refinestr_temp[j].ptv[hbondinfo[k][0][2]], &refinestr_temp[i].ptv[hbondinfo[k][0][0]], length_thresh);
                        if(!flag) continue;

                        Hbondobj hbondobjone = Hbondobj::create(&refinestr_temp[j].ptv[hbondinfo[k][0][1]], &refinestr_temp[j].ptv[hbondinfo[k][0][2]], &refinestr_temp[i].ptv[hbondinfo[k][0][0]], &gradient[j][hbondinfo[k][0][1]], &gradient[j][hbondinfo[k][0][2]], &gradient[i][hbondinfo[k][0][0]], BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean - BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_std, enecoeff["BBBaseHBond"]);
                        hbondobj.push_back(hbondobjone);
                    }
                }

                hbondinfo = BBHbondInfo[base2_str+base1_str];

                k = 0;
                AccBase = base2; DonBase = base1;
                length_thresh = BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean + HBOND_CONST*BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_std;

                //cout<<hbondinfo[k][0][0] << "-" << hbondinfo[k][0][1] << "-" << hbondinfo[k][0][2]<<endl;
                //cout<< BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean << "-" << length_thresh<<endl;

                flag = checkhbond(&refinestr_temp[i].ptv[hbondinfo[k][0][1]], &refinestr_temp[i].ptv[hbondinfo[k][0][2]], &refinestr_temp[j].ptv[hbondinfo[k][0][0]], length_thresh);

                if(enecoeff["BBHBond"] > 0){
                    if(flag){
                        Hbondobj hbondobjone = Hbondobj::create(&refinestr_temp[i].ptv[hbondinfo[k][0][1]], &refinestr_temp[i].ptv[hbondinfo[k][0][2]], &refinestr_temp[j].ptv[hbondinfo[k][0][0]], &gradient[i][hbondinfo[k][0][1]], &gradient[i][hbondinfo[k][0][2]], &gradient[j][hbondinfo[k][0][0]], BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean - BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_std, enecoeff["BBHBond"]);
                        hbondobj.push_back(hbondobjone);
                    }
                }

                if(enecoeff["BBBaseHBond"] > 0){
                    for(k=1;k<hbondinfo.size();k++)
                    {
                        AccBase = base2; DonBase = base1;

                        length_thresh = BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean + HBOND_CONST*BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_std;

                        //cout<<hbondinfo[k][0][0] << "-" << hbondinfo[k][0][1] << "-" << hbondinfo[k][0][2]<<endl;
                        //cout<< BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean << "-" << length_thresh<<endl;

                        flag = checkhbond(&refinestr_temp[i].ptv[hbondinfo[k][0][1]], &refinestr_temp[i].ptv[hbondinfo[k][0][2]], &refinestr_temp[j].ptv[hbondinfo[k][0][0]], length_thresh);
                        if(!flag) continue;

                        Hbondobj hbondobjone = Hbondobj::create(&refinestr_temp[i].ptv[hbondinfo[k][0][1]], &refinestr_temp[i].ptv[hbondinfo[k][0][2]], &refinestr_temp[j].ptv[hbondinfo[k][0][0]], &gradient[i][hbondinfo[k][0][1]], &gradient[i][hbondinfo[k][0][2]], &gradient[j][hbondinfo[k][0][0]], BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean - BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_std, enecoeff["BBBaseHBond"]);
                        hbondobj.push_back(hbondobjone);
                    }
                }
            }else{
                if(enecoeff["BBBaseHBond"] > 0){
                    hbondinfo = BBHbondInfo[base1_str+base2_str];

                    for(k=1;k<hbondinfo.size();k++)
                    {
                        AccBase = base1; DonBase = base2;

                        length_thresh = BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean + HBOND_CONST*BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_std;

                        //cout<<hbondinfo[k][0][0] << "-" << hbondinfo[k][0][1] << "-" << hbondinfo[k][0][2]<<endl;
                        //cout << BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean << "-" << length_thresh<<endl;

                        flag = checkhbond(&refinestr_temp[j].ptv[hbondinfo[k][0][1]], &refinestr_temp[j].ptv[hbondinfo[k][0][2]], &refinestr_temp[i].ptv[hbondinfo[k][0][0]], length_thresh);
                        if(!flag) continue;

                        Hbondobj hbondobjone = Hbondobj::create(&refinestr_temp[j].ptv[hbondinfo[k][0][1]], &refinestr_temp[j].ptv[hbondinfo[k][0][2]], &refinestr_temp[i].ptv[hbondinfo[k][0][0]], &gradient[j][hbondinfo[k][0][1]], &gradient[j][hbondinfo[k][0][2]], &gradient[i][hbondinfo[k][0][0]], BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean - BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_std, enecoeff["BBBaseHBond"]);
                        hbondobj.push_back(hbondobjone);
                    }
                }

                if(enecoeff["BBBaseHBond"] > 0){
                    hbondinfo = BBHbondInfo[base2_str+base1_str];

                    for(k=1;k<hbondinfo.size();k++)
                    {
                        AccBase = base2; DonBase = base1;

                        length_thresh = BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean + HBOND_CONST*BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_std;

                        //cout<<hbondinfo[k][0][0] << "-" << hbondinfo[k][0][1] << "-" << hbondinfo[k][0][2]<<endl;
                        //cout<< BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean << "-" << length_thresh<<endl;

                        flag = checkhbond(&refinestr_temp[i].ptv[hbondinfo[k][0][1]], &refinestr_temp[i].ptv[hbondinfo[k][0][2]], &refinestr_temp[j].ptv[hbondinfo[k][0][0]], length_thresh);
                        if(!flag) continue;

                        Hbondobj hbondobjone = Hbondobj::create(&refinestr_temp[i].ptv[hbondinfo[k][0][1]], &refinestr_temp[i].ptv[hbondinfo[k][0][2]], &refinestr_temp[j].ptv[hbondinfo[k][0][0]], &gradient[i][hbondinfo[k][0][1]], &gradient[i][hbondinfo[k][0][2]], &gradient[j][hbondinfo[k][0][0]], BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean - BBHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_std, enecoeff["BBBaseHBond"]);
                        hbondobj.push_back(hbondobjone);
                    }
                }
            }
        }
    }

    if(enecoeff["HBond"] > 0){
        for(int i=0;i<cbasepairlist.size();i++)
        {
            res1=cbasepairlist[i][0];
            res2=cbasepairlist[i][1];

            base1 = refinestr_temp[res1].aaa;
            base2 = refinestr_temp[res2].aaa;

            base1_str=string(1, base1);
            base2_str=string(1, base2);
            
            hbondinfo = BpHbondInfo[base1_str+base2_str];
            for(k=0;k<hbondinfo.size();k++)
            {
                if(hbondinfo[k].size() > 1)
                {
                    if(hbondinfo[k][0][3] == 0){
                        AccBase = base1; DonBase = base2;

                        Hbond2hydobj hbondobjone = Hbond2hydobj::create(&refinestr_temp[res2].ptv[hbondinfo[k][0][1]], &refinestr_temp[res2].ptv[hbondinfo[k][1][1]], &refinestr_temp[res2].ptv[hbondinfo[k][0][2]], &refinestr_temp[res1].ptv[hbondinfo[k][0][0]], &gradient[res2][hbondinfo[k][0][1]], &gradient[res2][hbondinfo[k][1][1]], &gradient[res2][hbondinfo[k][0][2]], &gradient[res1][hbondinfo[k][0][0]], BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean - BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_std, cbasepairscorelist[i]*enecoeff["HBond"]);
                        hbondbpobj.push_back(hbondobjone);
                    }else{
                        AccBase = base2; DonBase = base1;
                        
                        Hbond2hydobj hbondobjone = Hbond2hydobj::create(&refinestr_temp[res1].ptv[hbondinfo[k][0][1]], &refinestr_temp[res1].ptv[hbondinfo[k][1][1]], &refinestr_temp[res1].ptv[hbondinfo[k][0][2]], &refinestr_temp[res2].ptv[hbondinfo[k][0][0]], &gradient[res1][hbondinfo[k][0][1]], &gradient[res1][hbondinfo[k][1][1]], &gradient[res1][hbondinfo[k][0][2]], &gradient[res2][hbondinfo[k][0][0]], BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean - BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_std, cbasepairscorelist[i]*enecoeff["HBond"]);
                        hbondbpobj.push_back(hbondobjone);
                    }
                }else{
                    if(hbondinfo[k][0][3] == 0){
                        AccBase = base1; DonBase = base2;
                        
                        Hbondobj hbondobjone = Hbondobj::create(&refinestr_temp[res2].ptv[hbondinfo[k][0][1]], &refinestr_temp[res2].ptv[hbondinfo[k][0][2]], &refinestr_temp[res1].ptv[hbondinfo[k][0][0]], &gradient[res2][hbondinfo[k][0][1]], &gradient[res2][hbondinfo[k][0][2]], &gradient[res1][hbondinfo[k][0][0]], BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean - BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_std, cbasepairscorelist[i]*enecoeff["HBond"]);
                        hbondobj.push_back(hbondobjone);
                    }else{
                        AccBase = base2; DonBase = base1;
                        
                        Hbondobj hbondobjone = Hbondobj::create(&refinestr_temp[res1].ptv[hbondinfo[k][0][1]], &refinestr_temp[res1].ptv[hbondinfo[k][0][2]], &refinestr_temp[res2].ptv[hbondinfo[k][0][0]], &gradient[res1][hbondinfo[k][0][1]], &gradient[res1][hbondinfo[k][0][2]], &gradient[res2][hbondinfo[k][0][0]], BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean - BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_std, cbasepairscorelist[i]*enecoeff["HBond"]);
                        hbondobj.push_back(hbondobjone);
                    }
                }
            }
        }
        

        for(int i=0;i<noncbasepairlist.size();i++){
            res1=noncbasepairlist[i][0];
            res2=noncbasepairlist[i][1];

            base1 = refinestr_temp[res1].aaa;
            base2 = refinestr_temp[res2].aaa;

            base1_str=string(1, base1);
            base2_str=string(1, base2);

            //cout<<"bbb"<<endl;
            //cout<<base1 << "-" << base2 <<"-" << noncbasepairscorelist[i] << endl;

            hbondinfo = BpHbondInfo[base1_str+base2_str];
            for(k=0;k<hbondinfo.size();k++)
            {
                if(hbondinfo[k].size() > 1)
                {
                    if(hbondinfo[k][0][3] == 0){
                        AccBase = base1; DonBase = base2;

                        length_thresh = BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean + HBOND_CONST*BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_std;

                        flag = checkhbond2(&refinestr_temp[res2].ptv[hbondinfo[k][0][1]], &refinestr_temp[res2].ptv[hbondinfo[k][1][1]], &refinestr_temp[res2].ptv[hbondinfo[k][0][2]], &refinestr_temp[res1].ptv[hbondinfo[k][0][0]], length_thresh);
                        
                        if(!flag) continue;

                        Hbond2hydobj hbondobjone = Hbond2hydobj::create(&refinestr_temp[res2].ptv[hbondinfo[k][0][1]], &refinestr_temp[res2].ptv[hbondinfo[k][1][1]], &refinestr_temp[res2].ptv[hbondinfo[k][0][2]], &refinestr_temp[res1].ptv[hbondinfo[k][0][0]], &gradient[res2][hbondinfo[k][0][1]], &gradient[res2][hbondinfo[k][1][1]], &gradient[res2][hbondinfo[k][0][2]], &gradient[res1][hbondinfo[k][0][0]], BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean - BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_std, noncbasepairscorelist[i]*enecoeff["HBond"]);
                        hbondbpobj.push_back(hbondobjone);
                    }else{
                        AccBase = base2; DonBase = base1;

                        length_thresh = BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean + HBOND_CONST*BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_std;

                        flag = checkhbond2(&refinestr_temp[res1].ptv[hbondinfo[k][0][1]], &refinestr_temp[res1].ptv[hbondinfo[k][1][1]], &refinestr_temp[res1].ptv[hbondinfo[k][0][2]], &refinestr_temp[res2].ptv[hbondinfo[k][0][0]], length_thresh);
                        
                        if(!flag) continue;
                        
                        Hbond2hydobj hbondobjone = Hbond2hydobj::create(&refinestr_temp[res1].ptv[hbondinfo[k][0][1]], &refinestr_temp[res1].ptv[hbondinfo[k][1][1]], &refinestr_temp[res1].ptv[hbondinfo[k][0][2]], &refinestr_temp[res2].ptv[hbondinfo[k][0][0]], &gradient[res1][hbondinfo[k][0][1]], &gradient[res1][hbondinfo[k][1][1]], &gradient[res1][hbondinfo[k][0][2]], &gradient[res2][hbondinfo[k][0][0]], BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean - BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_std, noncbasepairscorelist[i]*enecoeff["HBond"]);
                        hbondbpobj.push_back(hbondobjone);
                    }
                }else{
                    if(hbondinfo[k][0][3] == 0){
                        AccBase = base1; DonBase = base2;

                        length_thresh = BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean + HBOND_CONST*BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_std;

                        flag = checkhbond(&refinestr_temp[res2].ptv[hbondinfo[k][0][1]], &refinestr_temp[res2].ptv[hbondinfo[k][0][2]], &refinestr_temp[res1].ptv[hbondinfo[k][0][0]], length_thresh);

                        if(!flag) continue;

                        Hbondobj hbondobjone = Hbondobj::create(&refinestr_temp[res2].ptv[hbondinfo[k][0][1]], &refinestr_temp[res2].ptv[hbondinfo[k][0][2]], &refinestr_temp[res1].ptv[hbondinfo[k][0][0]], &gradient[res2][hbondinfo[k][0][1]], &gradient[res2][hbondinfo[k][0][2]], &gradient[res1][hbondinfo[k][0][0]], BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean - BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_std, noncbasepairscorelist[i]*enecoeff["HBond"]);
                        hbondobj.push_back(hbondobjone);
                    }else{
                        AccBase = base2; DonBase = base1;

                        length_thresh = BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean + HBOND_CONST*BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_std;

                        flag = checkhbond(&refinestr_temp[res1].ptv[hbondinfo[k][0][1]], &refinestr_temp[res1].ptv[hbondinfo[k][0][2]], &refinestr_temp[res2].ptv[hbondinfo[k][0][0]], length_thresh);

                        if(!flag) continue;
                        
                        Hbondobj hbondobjone = Hbondobj::create(&refinestr_temp[res1].ptv[hbondinfo[k][0][1]], &refinestr_temp[res1].ptv[hbondinfo[k][0][2]], &refinestr_temp[res2].ptv[hbondinfo[k][0][0]], &gradient[res1][hbondinfo[k][0][1]], &gradient[res1][hbondinfo[k][0][2]], &gradient[res2][hbondinfo[k][0][0]], BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].adlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].ahlength_mean, BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_mean - BPHBOND[AccBase][DonBase][hbondinfo[k][0][0]][hbondinfo[k][0][2]].dhaangle_std, noncbasepairscorelist[i]*enecoeff["HBond"]);
                        hbondobj.push_back(hbondobjone);
                    }
                }
            }
        }
    }
}


void MCgetgradientinit(list<InitPosobj> binitobj){
	for (std::list<InitPosobj>::iterator it = binitobj.begin(); it != binitobj.end(); ++it) {
		it->calcgradientinit();
	}
}

double MCgetenergyinit(list<InitPosobj> binitobj){
	double energy=0.0;
	for (std::list<InitPosobj>::const_iterator it = binitobj.begin(); it != binitobj.end(); ++it) {
		energy+=it->calcenergyinit();
	}
	
	return energy;
}

void MCgetgradientinitdist(list<InitPairobj> binitdistobj){
	for (std::list<InitPairobj>::iterator it = binitdistobj.begin(); it != binitdistobj.end(); ++it) {
        it->calcgradientinitpair();
	}
}

void MCgetgradienttemplatedist(list<InitPairobj> btemplatedistobj){
	for (std::list<InitPairobj>::iterator it = btemplatedistobj.begin(); it != btemplatedistobj.end(); ++it) {
        it->calcgradientinitpair();
	}
}

double MCgetenergyinitdist(list<InitPairobj> binitdistobj){
	double energy=0.0;
	for (std::list<InitPairobj>::const_iterator it = binitdistobj.begin(); it != binitdistobj.end(); ++it) {
		energy+=it->calcenergyinitpair();
	}
	
	return energy;
}

double MCgetenergytemplatedist(list<InitPairobj> btemplatedistobj){
	double energy=0.0;
	for (std::list<InitPairobj>::const_iterator it = btemplatedistobj.begin(); it != btemplatedistobj.end(); ++it) {
		energy+=it->calcenergyinitpair();
	}
	
	return energy;
}

void MCgetgradientbasepair(list<BasePairobj> &basepairobj){
	for (std::list<BasePairobj>::iterator it = basepairobj.begin(); it != basepairobj.end(); ++it) {
		it->calcgradientbasepair();
	}
}

double MCgetenergybasepair(list<BasePairobj> &basepairobj){
	double energy=0.0;
	for (std::list<BasePairobj>::iterator it = basepairobj.begin(); it != basepairobj.end(); ++it) {
		energy+=it->calcenergybasepair();
	}
	
	return energy;
}

void MCgetgradientbaseplane(list<BasePairobj> &baseplaneobj){
	for (std::list<BasePairobj>::iterator it = baseplaneobj.begin(); it != baseplaneobj.end(); ++it) {
		it->calcgradientbasepair();
	}
}

double MCgetenergybaseplane(list<BasePairobj> &baseplaneobj){
	double energy=0.0;
	for (std::list<BasePairobj>::iterator it = baseplaneobj.begin(); it != baseplaneobj.end(); ++it) {
		energy+=it->calcenergybasepair();
	}
	
	return energy;
}

void MCgetgradientbasepairdist(list<BPdistobj> &bpdistobj){
	for (std::list<BPdistobj>::const_iterator it = bpdistobj.begin(); it != bpdistobj.end(); ++it) {
		it->calcgradientbp();
	}
}

double MCgetenergybasepairdist(list<BPdistobj> &bpdistobj){
	double energy=0.0;
	for (std::list<BPdistobj>::const_iterator it = bpdistobj.begin(); it != bpdistobj.end(); ++it) {
		energy+=it->calcenergybp();
	}
	
	return energy;
}


void MCgetgradienthbond(list<Hbondobj> &hbondobj){
	for (std::list<Hbondobj>::const_iterator it = hbondobj.begin(); it != hbondobj.end(); ++it) {
		it->calcgradienthbond();
	}
}

double MCgetenergyhbond(list<Hbondobj> &hbondobj){
	double energy=0.0;
	for (std::list<Hbondobj>::const_iterator it = hbondobj.begin(); it != hbondobj.end(); ++it) {
		energy+=it->calcenergyhbond();
	}
	
	return energy;
}

void MCgetgradienthbond(list<Hbond2hydobj> &hbondobj){
	for (std::list<Hbond2hydobj>::const_iterator it = hbondobj.begin(); it != hbondobj.end(); ++it) {
		it->calcgradienthbond();
	}
}

void MCgetgradientstacking(list<Stackingobj> &stackobj){
	for (std::list<Stackingobj>::const_iterator it = stackobj.begin(); it != stackobj.end(); ++it) {
		it->calcgradientstack();
	}
}

double MCgetenergystacking(list<Stackingobj> &stackobj){
	double energy=0.0;
	for (std::list<Stackingobj>::const_iterator it = stackobj.begin(); it != stackobj.end(); ++it) {
		energy+=it->calcenergystack();
	}
	
	return energy;
}

void MCgetgradientstackinggeo(list<StackGeoobjMC> &stackgeoobj){
	for (std::list<StackGeoobjMC>::iterator it = stackgeoobj.begin(); it != stackgeoobj.end(); ++it) {
		it->calcgradient();
	}
}

double MCgetenergystackinggeo(list<StackGeoobjMC> &stackgeoobj){
	double energy=0.0;
	for (std::list<StackGeoobjMC>::const_iterator it = stackgeoobj.begin(); it != stackgeoobj.end(); ++it) {
		energy+=it->calcenergy();
	}
	
	return energy;
}



void MCgetgradientclash(list<Clashobj> &clashobj){
	for (std::list<Clashobj>::const_iterator it = clashobj.begin(); it != clashobj.end(); ++it) {
		it->calcgradientclash();
	}
}

double MCgetenergyclash(list<Clashobj> &clashobj){
	double energy=0.0;
	for (std::list<Clashobj>::const_iterator it = clashobj.begin(); it != clashobj.end(); ++it) {
		energy+=it->calcenergyclash();
	}
	
	return energy;
}

void MCgradientnttemp(list<SingleNtTemplateobjMC> &nttempobj){
	for (std::list<SingleNtTemplateobjMC>::const_iterator it = nttempobj.begin(); it != nttempobj.end(); ++it) {
		it->calcgradientstacktemp();
	}
}

double MCenergynttemp(list<SingleNtTemplateobjMC> &nttempobj){
	double energy=0.0;
	for (std::list<SingleNtTemplateobjMC>::iterator it = nttempobj.begin(); it != nttempobj.end(); ++it) {
		energy+=it->calcenergystacktemp();
	}
	
	return energy;
}

double MCgetenergyhbond(list<Hbond2hydobj> &hbondobj){
	double energy=0.0;
	for (std::list<Hbond2hydobj>::const_iterator it = hbondobj.begin(); it != hbondobj.end(); ++it) {
		energy+=it->calcenergyhbond();
	}
	
	return energy;
}

void MCgetgradientblength(list<BondLengthobj> &blengthobj){
	for (std::list<BondLengthobj>::const_iterator it = blengthobj.begin(); it != blengthobj.end(); ++it) {
		it->calcgradientlength();
	}
}

double MCgetenergyblength(const list<BondLengthobj> &blengthobj){
	double energy=0.0;
	for (std::list<BondLengthobj>::const_iterator it = blengthobj.begin(); it != blengthobj.end(); ++it) {
		energy+=it->calcenergylength();
	}
	
	return energy;
}

void MCgetgradientbangle(list<BondAngleobj> &bangleobj){
	for (std::list<BondAngleobj>::const_iterator it = bangleobj.begin(); it != bangleobj.end(); ++it) {
		it->calcgradientangle();
	}
}

void MCgetgradienthydtorsion(list<BondDihedralobj> &dihedralobj){
	for (std::list<BondDihedralobj>::const_iterator it = dihedralobj.begin(); it != dihedralobj.end(); ++it) {
		it->calcgradientdihedral();
	}
}

void MCgetgradientsugarringtorsion(list<MultiBondDihedralobj> &dihedralobj){
	for (std::list<MultiBondDihedralobj>::iterator it = dihedralobj.begin(); it != dihedralobj.end(); ++it) {
		it->calcgradientmultidihedral();
	}
}

double MCgetenergybangle(list<BondAngleobj> &bangleobj){
	double energy=0.0;
	for (std::list<BondAngleobj>::const_iterator it = bangleobj.begin(); it != bangleobj.end(); ++it) {
		energy+=it->calcenergyangle();
	}
	
	return energy;
}

double MCgetenergyhydtor(list<BondDihedralobj> &hydtorobj){
    double energy=0.0;
	for (std::list<BondDihedralobj>::const_iterator it = hydtorobj.begin(); it != hydtorobj.end(); ++it) {
		energy+=it->calcenergydihedral();
	}
	
	return energy;
}

double MCgetenergysugarringtor(list<MultiBondDihedralobj> &torobj){
    double energy=0.0;
	for (std::list<MultiBondDihedralobj>::const_iterator it = torobj.begin(); it != torobj.end(); ++it) {
		energy+=it->calcenergymultidihedral();
	}
	
	return energy;
}

void MCgetgradientbimpropers(list<BondImpropersobj> &improperobj){
	for (std::list<BondImpropersobj>::const_iterator it = improperobj.begin(); it != improperobj.end(); ++it) {
		it->calcgradientimpropers();
	}
}

double MCgetenergybimpropers(list<BondImpropersobj> &improperobj){
	double energy=0.0;
	for (std::list<BondImpropersobj>::const_iterator it = improperobj.begin(); it != improperobj.end(); ++it) {
		energy+=it->calcenergyimpropers();
	}
	
	return energy;
}

void MCgetgradientbdihedral(list<BondDihedralobj> &dihedralobj){
	for (std::list<BondDihedralobj>::const_iterator it = dihedralobj.begin(); it != dihedralobj.end(); ++it) {
		it->calcgradientdihedral();
	}
}

double MCgetenergybdihedral(list<BondDihedralobj> &dihedralobj){
	double energy=0.0;
	for (std::list<BondDihedralobj>::const_iterator it = dihedralobj.begin(); it != dihedralobj.end(); ++it) {
		energy+=it->calcenergydihedral();
	}
	
	return energy;
}

void MCgetgradientelectros(list<TElectrobj> &telecteobj){
	for (std::list<TElectrobj>::const_iterator it = telecteobj.begin(); it != telecteobj.end(); ++it) {
		it->calcgradientelectros();
	}
}

double MCgetenergyelectros(list<TElectrobj> &telecteobj){
	double energy=0.0;
	for (std::list<TElectrobj>::const_iterator it = telecteobj.begin(); it != telecteobj.end(); ++it) {
		energy+=it->calcenergyelectros();
	}
	
	return energy;
}

void MCgetgradientlj(list<LennardsJobj> &ljobj){
	for (std::list<LennardsJobj>::const_iterator it = ljobj.begin(); it != ljobj.end(); ++it) {
		it->calcgradientlj();
	}
}

double MCgetenergylj(list<LennardsJobj> &ljobj){
	double energy=0.0;
	for (std::list<LennardsJobj>::const_iterator it = ljobj.begin(); it != ljobj.end(); ++it) {
		energy+=it->calcenergylj();
	}
	
	return energy;
}


void MCgetgradient(list<BondLengthobj> &blengthobj, list<BondAngleobj> &bangleobj, list<BondDihedralobj> &dihedralobj, list<MultiBondDihedralobj> &sugardihedralobj, list<BasePairobj> &bplaneobj, list<BasePairobj> &basepairobj, list<BPdistobj> &bpdistobj, list<Stackingobj> &stackobj, list<StackGeoobjMC> &stackgeoobj, list<InitPosobj> &binitobj, list<InitPairobj> &binitpairobj, list<InitPairobj> &btemplatepairobj, list<TElectrobj> &telecteobj, list<LennardsJobj> &ljobj, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, list<Clashobj> &clashobj, list<SingleNtTemplateobjMC> &nttempobj)
{
	if(blengthobj.size() > 0)MCgetgradientblength(blengthobj);
	if(bangleobj.size() > 0) MCgetgradientbangle(bangleobj);
    if(dihedralobj.size() > 0) MCgetgradienthydtorsion(dihedralobj);
    if(sugardihedralobj.size() > 0) MCgetgradientsugarringtorsion(sugardihedralobj);
	if(basepairobj.size() > 0) MCgetgradientbasepair(basepairobj);
    if(bplaneobj.size() > 0) MCgetgradientbaseplane(bplaneobj);
    if(bpdistobj.size() > 0) MCgetgradientbasepairdist(bpdistobj);
    if(binitobj.size() > 0) MCgetgradientinit(binitobj);
    if(binitpairobj.size() > 0) MCgetgradientinitdist(binitpairobj);
    if(btemplatepairobj.size() > 0) MCgetgradienttemplatedist(btemplatepairobj);
    if(telecteobj.size() > 0) MCgetgradientelectros(telecteobj);
    if(ljobj.size() > 0) MCgetgradientlj(ljobj);
    if(hbondobj.size() > 0) MCgetgradienthbond(hbondobj);
    if(hbondbpobj.size() > 0) MCgetgradienthbond(hbondbpobj);
    if(stackobj.size() > 0) MCgetgradientstacking(stackobj);
    if(stackgeoobj.size() > 0) MCgetgradientstackinggeo(stackgeoobj);
    if(clashobj.size() > 0) MCgetgradientclash(clashobj);
    if(nttempobj.size() > 0) MCgradientnttemp(nttempobj);
}

double MCgetcovalentenergystem(list<BondLengthobj> &blengthobj, list<BondAngleobj> &bangleobj, list<BondDihedralobj> &dihedralobj, list<MultiBondDihedralobj> &sugardihedralobj, list<BasePairobj> &bplaneobj, list<BasePairobj> &basepairobj, list<BPdistobj> &bpdistobj, list<Stackingobj> &stackobj, list<StackGeoobjMC> &stackgeoobj, list<InitPosobj> &binitobj, list<InitPairobj> &binitpairobj, list<InitPairobj> &btemplatepairobj, list<TElectrobj> &telecteobj, list<LennardsJobj> &ljobj, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, list<Clashobj> &clashobj, list<SingleNtTemplateobjMC> &nttempobj)
{
    bool disp = false;
    double lengthene, angleene, torene, sugartorene, bplaneene, bpene, bpdistene, initene, initdistene, templatedistene, telecteene, ljene, hbondene, stackene, stackgeoene, clashene, nttempene;
    lengthene=0; angleene=0; torene=0; sugartorene=0; bplaneene=0; bpene=0; bpdistene=0; initene=0;  initdistene=0; templatedistene=0; telecteene=0; ljene=0; hbondene=0; stackene=0; stackgeoene=0; clashene=0; nttempene=0;
    if(blengthobj.size() > 0){
            //cout<<"a1"<<endl;
	    lengthene=MCgetenergyblength(blengthobj);
        if(disp) cout<<"Length: "<<lengthene<<endl;
    }
    if(bangleobj.size() > 0){
            //cout<<"a2"<<endl;
	    angleene=MCgetenergybangle(bangleobj);
        if(disp) cout<<"Angle: "<<angleene<<endl;
    }
    if(dihedralobj.size() > 0){
            //cout<<"a3"<<endl;
	    torene=MCgetenergyhydtor(dihedralobj);
        if(disp) cout<<"Hyd: "<<torene<<endl;
	   }
    
    if(sugardihedralobj.size() > 0){
            //cout<<"a3"<<endl;
	    sugartorene=MCgetenergysugarringtor(sugardihedralobj);
        if(disp) cout<<"Sugar tor: "<<sugartorene<<endl;
	   }

    if(basepairobj.size() > 0){
            //cout<<"a4"<<endl;
	    bpene=MCgetenergybasepair(basepairobj);
        if(disp) cout<<"Base pair:" <<  bpene<<endl;
    }

    if(bplaneobj.size() > 0){
            //cout<<"a4"<<endl;
	    bplaneene=MCgetenergybaseplane(bplaneobj);
        if(disp) cout<<"Base plane:" <<  bplaneene<<endl;
    }
    if(bpdistobj.size() > 0){
            //cout<<"a4"<<endl;
	    bpdistene=MCgetenergybasepairdist(bpdistobj);
        if(disp) cout<<"Base pair dist:" <<  bpdistene<<endl;
    }
    if(binitobj.size() > 0){
            //cout<<"a5"<<endl;
	    initene=MCgetenergyinit(binitobj);
        if(disp) cout<<"Init:" <<  initene<<endl;
    }
    if(binitpairobj.size() > 0){
            //cout<<"a5"<<endl;
	    initdistene=MCgetenergyinitdist(binitpairobj);
        if(disp) cout<<"Init Dist:" <<  initdistene<<endl;
    }
    if(btemplatepairobj.size() > 0){
            //cout<<"a5"<<endl;
	    templatedistene=MCgetenergytemplatedist(btemplatepairobj);
        if(disp) cout<<"Template Dist:" <<  templatedistene<<endl;
    }
    if(telecteobj.size() > 0){
	    //cout<<"a6"<<endl;
	   telecteene=MCgetenergyelectros(telecteobj);
       if(disp) cout<<"Electe:" <<  telecteene<<endl;
    }
    if(ljobj.size() > 0){
	   //cout<<"a7"<<endl;
	   ljene=MCgetenergylj(ljobj);
       if(disp) cout<<"LJ:" <<  ljene<<endl;
    }
    if(hbondobj.size() > 0){
            //cout<<"a8"<<endl;
	    hbondene+=MCgetenergyhbond(hbondobj);
        if(disp) cout<<"Hbond:" <<  hbondene<<endl;
    }
    if(hbondbpobj.size() > 0){
	    //cout<<"a9"<<endl;
	    hbondene+=MCgetenergyhbond(hbondbpobj);
        if(disp) cout<<"Hbond Bp:" <<  hbondene<<endl;
    }
    if(stackobj.size() > 0){
        stackene=MCgetenergystacking(stackobj);
        if(disp) cout<<"Stack:" <<  stackene<<endl;
    }
    if(stackgeoobj.size() > 0){
        stackgeoene=MCgetenergystackinggeo(stackgeoobj);
        if(disp) cout<<"Stack Geo:" <<  stackgeoene<<endl;
    }
    if(clashobj.size() > 0){
        clashene=MCgetenergyclash(clashobj);
        if(disp) cout<<"Clash:" <<  clashene<<endl;
    }

    if(nttempobj.size() > 0){
        nttempene=MCenergynttemp(nttempobj);
        if(disp) cout<<"Nt:" <<  nttempene<<endl;
    }

    return lengthene+angleene+torene+sugartorene+bplaneene+bpene+bpdistene+initene+initdistene+templatedistene+telecteene+ljene+hbondene+stackene+stackgeoene+clashene+nttempene;
}

void MCinitgrad(vector<vector<Point3s>> &gradptv)
{
	int i, mtype;
	for(i=0;i<gradptv.size();i++){
		for(mtype=0;mtype<gradptv[i].size();mtype++)
		{
			gradptv[i][mtype].x = 0;
			gradptv[i][mtype].y = 0;
			gradptv[i][mtype].z = 0;
		}
	}
}

void MonteCarloSimulation::MCupdatelbfgsfeatures(const vector<vector<Point3s>> &refinestr_temp_prev, const vector<vector<Point3s>> &grad, const vector<vector<Point3s>> &grad_prev, deque<vector<vector<Point3s>>> &s_list, deque<vector<vector<Point3s>>> &y_list, deque<double> &rho_list, const vector<bool> &movenodeinfo, int m){
    double rho_new =0.0;
    vector<vector<Point3s>> s_new, y_new;
    s_new.resize(grad.size());
    y_new.resize(grad.size());
    for(int i=0;i<grad.size();i++){
        if(movenodeinfo[i]){
            s_new[i].resize(grad[i].size());
            y_new[i].resize(grad[i].size());
            for (size_t j = 0; j < grad[i].size(); ++j) {
                s_new[i][j].x = refinestr_temp[i].ptv[j].x - refinestr_temp_prev[i][j].x;
                s_new[i][j].y = refinestr_temp[i].ptv[j].y - refinestr_temp_prev[i][j].y;
                s_new[i][j].z = refinestr_temp[i].ptv[j].z - refinestr_temp_prev[i][j].z;
                y_new[i][j].x = grad[i][j].x - grad_prev[i][j].x;
                y_new[i][j].y = grad[i][j].y - grad_prev[i][j].y;
                y_new[i][j].z = grad[i][j].z - grad_prev[i][j].z;

                rho_new+=(s_new[i][j].x*y_new[i][j].x+s_new[i][j].y*y_new[i][j].y+s_new[i][j].z*y_new[i][j].z);
            }
        }
    }
   
    if (rho_new <= 0.0) return;
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

void MClbfgsupdate(const vector<vector<Point3s>> &grad, const deque<vector<vector<Point3s>>> &s_list, const deque<vector<vector<Point3s>>> &y_list, const deque<double> &rho_list, const vector<bool> &movenodeinfo, vector<vector<Point3s>> &z){
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

double MCcalcaveragenorm(vector<vector<Point3s>> &gradient, double clip_val=-1.0)
{
	int i, mtype;
	double dx, dy, dz;
	double totalnorm = 0;
  	unsigned int num = 0;
    double norm, scale;

    /*
    for(i=0;i<gradient.size();i++){
		for(mtype=0;mtype<gradient[i].size();mtype++)
		{
			dx=gradient[i][mtype].x;
			dy=gradient[i][mtype].y;
			dz=gradient[i][mtype].z;

            norm=sqrt(dx*dx+dy*dy+dz*dz);

            if(norm > 1000){
                scale = 1000/norm;
                gradient[i][mtype].x*=scale;
                gradient[i][mtype].y*=scale;
                gradient[i][mtype].z*=scale;
            }
		}
	}
    */

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

void MonteCarloSimulation::MCupdatestructure(const vector<vector<Point3s>> &gradient, double step_size, const vector<bool> &moveresidues){
	int i, j;
	char base;
	for(i=0;i<refinestr_temp.size();i++){
        if(!moveresidues[i]) continue;
        base = refinestr_temp[i].aaa;
		for(j=0;j<AccAtomParamNum[base][4];j++)
		{
            refinestr_temp[i].ptv[j].x += step_size*gradient[i][j].x;
			refinestr_temp[i].ptv[j].y += step_size*gradient[i][j].y;
			refinestr_temp[i].ptv[j].z += step_size*gradient[i][j].z;
		}
	}
}

double MonteCarloSimulation::MClinesearch(const vector<vector<Point3s>> &gradient, list<BondLengthobj> &blengthobj, list<BondAngleobj> &bangleobj, list<BondDihedralobj> &dihedralobj, list<MultiBondDihedralobj> &sugardihedralobj, list<BasePairobj> &bplaneobj, list<BasePairobj> &basepairobj, list<BPdistobj> &bpdistobj, list<Stackingobj> &stackobj, list<StackGeoobjMC> &stackgeoobj, list<InitPosobj> &binitobj, list<InitPairobj> &binitdistobj, list<InitPairobj> &btemplatedistobj, list<TElectrobj> &telecteobj, list<LennardsJobj> &ljobj, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, list<Clashobj> &clashobj, list<SingleNtTemplateobjMC> &nttempobj, const vector<bool> &moveresidues, double alpha0, double rho, double beta)
{
    double alpha = alpha0;

    double new_energy;
    double init_energy =  MCgetcovalentenergystem(blengthobj, bangleobj, dihedralobj, sugardihedralobj, bplaneobj, basepairobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitdistobj, btemplatedistobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);

    double grad_norm_sq = 0.0;
    for (int i=0;i<gradient.size();i++) {
        for (int j=0;j<gradient[i].size();j++) {
            grad_norm_sq += gradient[i][j].x * gradient[i][j].x + gradient[i][j].y * gradient[i][j].y + gradient[i][j].z * gradient[i][j].z;
        }
    }

    //cout<<"============"<<endl;
    while (true) {
        //cout<<"linesearch-" << alpha<<endl;
        MCupdatestructure(gradient, -alpha, moveresidues);
        new_energy =  MCgetcovalentenergystem(blengthobj, bangleobj, dihedralobj, sugardihedralobj, bplaneobj, basepairobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitdistobj, btemplatedistobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);

        if (new_energy <= init_energy) {
            break;
        }else{
            MCupdatestructure(gradient, alpha, moveresidues);
        }

        alpha *= rho;
        if(alpha < 1e-8) break;
    }

    return new_energy;
}

double MonteCarloSimulation::MCgoldenelectionsearch(const vector<vector<Point3s>> &gradient, list<BondLengthobj> &blengthobj, list<BondAngleobj> &bangleobj, list<BondDihedralobj> &dihedralobj, list<MultiBondDihedralobj> &sugardihedralobj, list<BasePairobj> &bplaneobj, list<BasePairobj> &basepairobj, list<BPdistobj> &bpdistobj, list<Stackingobj> &stackobj, list<StackGeoobjMC> &stackgeoobj, list<InitPosobj> &binitobj, list<InitPairobj> &binitdistobj, list<InitPairobj> &btemplatedistobj, list<TElectrobj> &telecteobj, list<LennardsJobj> &ljobj, list<Hbondobj> &hbondobj, list<Hbond2hydobj> &hbondbpobj, list<Clashobj> &clashobj, list<SingleNtTemplateobjMC> &nttempobj, const vector<bool> &moveresidues)
{
	static const double goldenratio = .5*(sqrt(5.)-1.);
	double step = 1.e-4;
	double left = step*(1.-goldenratio);
	double right = step*goldenratio;
        
	double energyst = MCgetcovalentenergystem(blengthobj, bangleobj, dihedralobj, sugardihedralobj, bplaneobj, basepairobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitdistobj, btemplatedistobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);
	
	MCupdatestructure(gradient, -right, moveresidues);
	double energyR = MCgetcovalentenergystem(blengthobj, bangleobj, dihedralobj, sugardihedralobj, bplaneobj, basepairobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitdistobj, btemplatedistobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);

	MCupdatestructure(gradient, right-left, moveresidues);
	double energyL = MCgetcovalentenergystem(blengthobj, bangleobj, dihedralobj, sugardihedralobj, bplaneobj, basepairobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitdistobj, btemplatedistobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);

	int stepnum = 0;
	while(fabs(step) > 1e-10)
  	{
        if(energyL < energyR || (energyL > energyst && energyR > energyst))
    	{
			step = step - step*(1-goldenratio);
			right = left;
			energyR = energyL;
			left = step*(1.-goldenratio);
			MCupdatestructure(gradient, right-left, moveresidues);
			energyL = MCgetcovalentenergystem(blengthobj, bangleobj, dihedralobj, sugardihedralobj, bplaneobj, basepairobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitdistobj, btemplatedistobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);
		}else
        {
            step = step - step*(1-goldenratio);
            //left = right;
            left = right - left;
			energyL = energyR;
            right = step*goldenratio;
	        MCupdatestructure(gradient, -right, moveresidues);
	        energyR = MCgetcovalentenergystem(blengthobj, bangleobj, dihedralobj, sugardihedralobj, bplaneobj, basepairobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitdistobj, btemplatedistobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);
	        MCupdatestructure(gradient, right-left, moveresidues);
	    }
	}

	MCupdatestructure(gradient, left-step/2, moveresidues);
    //return (energyL+energyR)/2.;
    return MCgetcovalentenergystem(
        blengthobj, bangleobj, dihedralobj, sugardihedralobj, 
        bplaneobj, basepairobj, bpdistobj, stackobj, stackgeoobj, 
        binitobj, binitdistobj, btemplatedistobj, telecteobj, 
        ljobj, hbondobj, hbondbpobj, clashobj, nttempobj
    );
}

void MonteCarloSimulation::supbpfrag(int res1, int res2, int &minfragidx, double &minrmsd){
	int i,j;
	int atom;
	int count;
	vector<vector<double>> RotMatix;
    vector<double> TranVect;
    //vector<vector<double>> minRotMatix;
    //vector<double> minTranVect;
	minrmsd=100000000000;
	double rmsd;
    vector<double> xyz(3,0);
	vector<double> xyztmp(3,0);
	Point3s vec;
	BasicFunc bf;
    char base1, base2;
    int base_idx_1, base_idx_2;

    base1 = refinestr_temp[res1].aaa;
    base2 = refinestr_temp[res2].aaa;

    static map<char, vector<int>> targetatoms={
	    {'A',{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21}},
        {'C',{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19}},
        {'G',{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}},
        {'U',{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19}}
	};

    vector<vector<Residueatoms>> bptemplate=BPTEMPLATE[base1][base2];

	vector<vector<double> > xyz_list1(targetatoms[base1].size()+targetatoms[base2].size(), vector<double>(3,0));
    vector<vector<double> > xyz_list2(targetatoms[base1].size()+targetatoms[base2].size(), vector<double>(3,0));

    base_idx_1=0;
    base_idx_2=1;
    
	for(i=0;i<bptemplate.size();i++)
	{
		count = 0;
		for(j=0; j<targetatoms[base1].size();j++){
            atom=targetatoms[base1][j];
            xyz_list1[count][0]=refinestr_temp[res1].ptv[atom].x;
            xyz_list1[count][1]=refinestr_temp[res1].ptv[atom].y;
            xyz_list1[count][2]=refinestr_temp[res1].ptv[atom].z;

            xyz_list2[count][0]=BPTEMPLATE[base1][base2][i][0].ptv[atom].x;
            xyz_list2[count][1]=BPTEMPLATE[base1][base2][i][0].ptv[atom].y;
            xyz_list2[count][2]=BPTEMPLATE[base1][base2][i][0].ptv[atom].z;

            count+=1;
		}

		for(j=0; j<targetatoms[base2].size();j++){
            atom=targetatoms[base2][j];
            xyz_list1[count][0]=refinestr_temp[res2].ptv[atom].x;
            xyz_list1[count][1]=refinestr_temp[res2].ptv[atom].y;
            xyz_list1[count][2]=refinestr_temp[res2].ptv[atom].z;

            xyz_list2[count][0]=BPTEMPLATE[base1][base2][i][1].ptv[atom].x;
            xyz_list2[count][1]=BPTEMPLATE[base1][base2][i][1].ptv[atom].y;
            xyz_list2[count][2]=BPTEMPLATE[base1][base2][i][1].ptv[atom].z;

            count+=1;
		}

		bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

		rmsd=0;
		count=0;
		for(j=0; j<targetatoms[base1].size();j++){
            atom=targetatoms[base1][j];
            bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);
            vec.x=refinestr_temp[res1].ptv[atom].x-xyztmp[0];
            vec.y=refinestr_temp[res1].ptv[atom].y-xyztmp[1];
            vec.z=refinestr_temp[res1].ptv[atom].z-xyztmp[2];

            rmsd+=(vec.x*vec.x+vec.y*vec.y+vec.z*vec.z);
            count+=1;
		}

		for(j=0; j<targetatoms[base2].size();j++){
            atom=targetatoms[base2][j];
            bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);

            vec.x=refinestr_temp[res2].ptv[atom].x-xyztmp[0];
            vec.y=refinestr_temp[res2].ptv[atom].y-xyztmp[1];
            vec.z=refinestr_temp[res2].ptv[atom].z-xyztmp[2];

            rmsd+=(vec.x*vec.x+vec.y*vec.y+vec.z*vec.z);
            count+=1;
		}
        rmsd = sqrt(rmsd/count);

		if(minrmsd > rmsd){
            minrmsd=rmsd;
            minfragidx=i;
            //minRotMatix = RotMatix;
            //minTranVect = TranVect;
        }
	}

    static map<char, vector<int>> targetatoms_rep={
	    {'A',{12,13,14,15,16,17,18,19,20,21}},
        {'C',{12,13,14,15,16,17,18,19}},
        {'G',{12,13,14,15,16,17,18,19,20,21,22}},
        {'U',{12,13,14,15,16,17,18,19}}
	};

    vector<vector<double> > xyz_list3(targetatoms_rep[base1].size()+targetatoms_rep[base2].size(), vector<double>(3,0));
    vector<vector<double> > xyz_list4(targetatoms_rep[base1].size()+targetatoms_rep[base2].size(), vector<double>(3,0));

    if(minrmsd < 2.0){
        count = 0;
		for(j=0; j<targetatoms_rep[base1].size();j++){
            atom=targetatoms[base1][j];
            xyz_list3[count][0]=refinestr_temp[res1].ptv[atom].x;
            xyz_list3[count][1]=refinestr_temp[res1].ptv[atom].y;
            xyz_list3[count][2]=refinestr_temp[res1].ptv[atom].z;

            xyz_list4[count][0]=BPTEMPLATE[base1][base2][minfragidx][base_idx_1].ptv[atom].x;
            xyz_list4[count][1]=BPTEMPLATE[base1][base2][minfragidx][base_idx_1].ptv[atom].y;
            xyz_list4[count][2]=BPTEMPLATE[base1][base2][minfragidx][base_idx_1].ptv[atom].z;

            count+=1;
		}

		for(j=0; j<targetatoms_rep[base2].size();j++){
            atom=targetatoms[base2][j];
            xyz_list3[count][0]=refinestr_temp[res2].ptv[atom].x;
            xyz_list3[count][1]=refinestr_temp[res2].ptv[atom].y;
            xyz_list3[count][2]=refinestr_temp[res2].ptv[atom].z;

            xyz_list4[count][0]=BPTEMPLATE[base1][base2][minfragidx][base_idx_2].ptv[atom].x;
            xyz_list4[count][1]=BPTEMPLATE[base1][base2][minfragidx][base_idx_2].ptv[atom].y;
            xyz_list4[count][2]=BPTEMPLATE[base1][base2][minfragidx][base_idx_2].ptv[atom].z;

            count+=1;
		}

        bf.RotateCoor(xyz_list3,xyz_list4,RotMatix,TranVect);

        for(j=0; j<targetatoms_rep[base1].size();j++)
        {
            atom=targetatoms_rep[base1][j];

            xyz[0]=BPTEMPLATE[base1][base2][minfragidx][base_idx_1].ptv[atom].x;
            xyz[1]=BPTEMPLATE[base1][base2][minfragidx][base_idx_1].ptv[atom].y;
            xyz[2]=BPTEMPLATE[base1][base2][minfragidx][base_idx_1].ptv[atom].z;

            bf.Transform(xyz,RotMatix,TranVect,xyztmp);
            refinestr_temp[res1].ptv[atom].x = xyztmp[0];
            refinestr_temp[res1].ptv[atom].y = xyztmp[1];
            refinestr_temp[res1].ptv[atom].z = xyztmp[2];
        }

        for(j=0; j<targetatoms_rep[base2].size();j++)
        {
            atom=targetatoms_rep[base2][j];

            xyz[0]=BPTEMPLATE[base1][base2][minfragidx][base_idx_2].ptv[atom].x;
            xyz[1]=BPTEMPLATE[base1][base2][minfragidx][base_idx_2].ptv[atom].y;
            xyz[2]=BPTEMPLATE[base1][base2][minfragidx][base_idx_2].ptv[atom].z;

            bf.Transform(xyz,RotMatix,TranVect,xyztmp);
            refinestr_temp[res2].ptv[atom].x = xyztmp[0];
            refinestr_temp[res2].ptv[atom].y = xyztmp[1];
            refinestr_temp[res2].ptv[atom].z = xyztmp[2];
        }
    }
}

void MonteCarloSimulation::supbpfragmain(const vector<vector<int>> &basepairlist)
{
	int i,j,k;
	int res1, res2;
	string base1, base2;
    int minfragidx; 
    double minrmsd;
    int baseorder;

    map<string, string> basebpname = {{"A-U", "A-U"}, {"U-A", "A-U"}, {"G-C", "G-C"}, {"C-G", "G-C"}, {"G-U", "G-U"}, {"U-G", "G-U"}};
    
    for(int i=0;i<basepairlist.size();i++)
    {
        res1=basepairlist[i][0];
        res2=basepairlist[i][1];
        supbpfrag(res1, res2, minfragidx, minrmsd);
    }
}

void MonteCarloSimulation::MCfragmentoptimization(
    const map<int, vector<Point3s>> &initstructure, const vector<vector<int>> &bindingmap,
    const vector<bool> &moveresidues, const vector<vector<int>> &contactmap, const vector<vector<double>> &stackingmap,
    const vector<vector<int>> &basepairlist, const vector<vector<int>> &noncbasepairlist,
    const vector<double> &basepairscorelist, const vector<double> &noncbasepairscorelist, 
    const bool bpoption){
    double weight_rand=Random();

    //map<string, double> enecoeff={{"length", 400}, {"angle",120}, {"hydtor", 0}, {"dihedral", 0}, {"init", 0}, {"init_dist", 1.0}, {"template_dist", 1.0}, {"improper", 0}, {"electros", 1.0}, {"LJ", 1.0}, {"HBond", 20.0}, {"BP", 1.0}, {"BPdist", 0.0}, {"Stack", 2.0}, {"StackGeo", 5.0}, {"Clash", 5.0}, {"Nt", 5.0}, {"BasePlane", 5.0}};
    //map<string, double> enecoeff={{"length", 400}, {"angle",120}, {"hydtor", 0}, {"dihedral", 0}, {"init", 0}, {"init_dist", 0.5}, {"template_dist", 0.5}, {"improper", 0}, {"electros", 1.0}, {"LJ", 1.0}, {"HBond", 20.0}, {"BP", 0.0}, {"BPdist", 0.0}, {"Stack", 0.5}, {"StackGeo", 0.0}, {"Clash", 5.0}, {"Nt", 0.0}, {"BasePlane", 5.0}};
    map<string, double> enecoeff={{"length", 400}, {"angle",120}, {"tor", 0.0}, {"sugar", 0.0}, {"dihedral", 0}, {"init", 0}, {"init_dist", 3.0}, {"template_dist", 3.0}, {"improper", 0}, {"electros", 1.0}, {"LJ", 1.0}, {"HBond", 25.0}, {"BBHBond", 7.5}, {"BBBaseHBond", 0.0}, {"BP", 7.5}, {"BPdist", 3.0}, {"Stack", 0.5}, {"StackGeo", 0.0}, {"Clash", 5.0}, {"Nt", 0.0}, {"BasePlane", 17.5}};
    //map<string, double> enecoeff={{"length", 0}, {"angle",0}, {"tor", 0.0}, {"sugar", 0.0}, {"init", 0}, {"init_dist", 0}, {"template_dist", 0}, {"improper", 0}, {"electros", 0}, {"LJ", 0}, {"HBond", 0.0}, {"BBHBond", 5.0}, {"BP", 0.0}, {"BPdist", 0.0}, {"Stack", 0.0}, {"StackGeo", 0.0}, {"Clash", 0}, {"Nt", 0}, {"BasePlane", 0.0}};

    //optimization
    vector<vector<Point3s>> gradptv;
    MCprepgrad(gradptv);

    vector<vector<Point3s>> refinestr_temp_prev(numbb);
    for(int i=0;i<numbb;i++){
        if(moveresidues[i]){
            char base=refinestr_temp[i].aaa;
            refinestr_temp_prev[i].resize(AccAtomParamNum[base][4]);
        }
    }

    vector<vector<Point3s>> gradptv_prev;
    deque<vector<vector<Point3s>>> s_list, y_list;
    deque<double> rho_list;
    vector<vector<Point3s>> z;
    MCprepgrad(z);
    int m=15;

    list<InitPosobj> binitobj;
    list<InitPairobj> binitpairobj;
    list<InitPairobj> btemplatepairobj;
    list<BondLengthobj> blengthobj;
	list<BondAngleobj> bangleobj;
    list<BondDihedralobj> dihedralobj;
    list<MultiBondDihedralobj> sugardihedralobj;
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
    list<SingleNtTemplateobjMC> nttempobj;

    vector<vector<int>> nt_indices;
    vector<vector<int>> bp_indices;
    vector<vector<Point3s*>> basestackptrsA;
    vector<vector<Point3s*>> basestackptrsB;
    vector<vector<Point3s*>> basestackgradA;
    vector<vector<Point3s*>> basestackgradB;
    vector<vector<Point3s*>> AtomPtrslist;
    vector<vector<Point3s*>> Atomgradlist;

    vector<vector<Point3s*>> BaseAtomPtrslist;
    vector<vector<Point3s*>> BaseAtomgradlist;

    if(enecoeff["init"] > 0) MCsetinitobj(initstructure, gradptv, binitobj, enecoeff["init"]);
    if(enecoeff["init_dist"] > 0) MCsetinitpairobj(initstructure, moveresidues, gradptv, binitpairobj, enecoeff["init_dist"]/(1+(*MCtemplatedistmap).size()));
    if(enecoeff["template_dist"] > 0) MCsettemplatepairobj(moveresidues, gradptv, btemplatepairobj, enecoeff["init_dist"]/(1+(*MCtemplatedistmap).size()));

    if(enecoeff["length"] > 0) MCsetlengthobj(bindingmap, gradptv, enecoeff["length"], blengthobj);

    if(enecoeff["angle"] > 0) MCsetangleobj(bindingmap, gradptv, enecoeff["angle"], bangleobj);

    if(enecoeff["tor"] > 0) MCsetdangleobj(bindingmap, gradptv, enecoeff["tor"], dihedralobj);

    if(enecoeff["sugar"] > 0) MCsetdsugarangleobj(bindingmap, gradptv, enecoeff["sugar"], sugardihedralobj);
    
    if(enecoeff["BP"] > 0) MCsetbasepairobj(gradptv, basepairlist, noncbasepairlist, basepairscorelist, noncbasepairscorelist, AtomPtrslist, Atomgradlist, enecoeff["BP"], bpobj);

    if(enecoeff["BasePlane"] > 0) MCsetbaseplaneobj(gradptv, BaseAtomPtrslist, BaseAtomgradlist, moveresidues, enecoeff["BasePlane"], baseplaneobj);

    if(enecoeff["BPdist"] > 0) MCsetbpdistobj(gradptv, basepairlist, basepairscorelist, enecoeff["BPdist"], bpdistobj);

    if(enecoeff["Stack"] > 0) MCsetstackingobj(gradptv, stackingmap, enecoeff["Stack"], stackobj);

    if(enecoeff["StackGeo"] > 0) MCsetstackinggeoobj(gradptv, stackingmap, basestackptrsA, basestackptrsB, basestackgradA, basestackgradB, enecoeff["StackGeo"], stackgeoobj);

    if(enecoeff["Nt"] > 0) MCsetnttemplateobj(gradptv, nt_indices, enecoeff["Nt"], nttempobj);

    if(enecoeff["electros"] > 0 || enecoeff["LJ"] > 0 || enecoeff["HBond"] > 0 || enecoeff["BBHBond"] > 0 || enecoeff["BBBaseHBond"] > 0 || enecoeff["Clash"] > 0) MCsetNonbondedobj(gradptv, contactmap, bindingmap, basepairlist, noncbasepairlist, basepairscorelist, noncbasepairscorelist, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, enecoeff);

    if(enecoeff["electros"] <= 0) telecteobj.clear();
    if(enecoeff["LJ"] <= 0) ljobj.clear();
    
    if(enecoeff["HBond"] <= 0 && enecoeff["BBHBond"] <= 0 && enecoeff["BBBaseHBond"] <= 0){
        hbondobj.clear();
        hbondbpobj.clear();
    }

    if(enecoeff["Clash"] <= 0){
        clashobj.clear();
    }

    int estop_counter, stop_count;
    double averagenorm, energy, minenergy;
    int steps=200;
    stop_count = 3;
    minenergy = 1.0 / 0.0;
    estop_counter=0;  

    bool reset_frag = true;
    int proc_step = 4;
    int frag_rest_iter = 15;

    double startene = MCgetcovalentenergystem(blengthobj, bangleobj, dihedralobj, sugardihedralobj, baseplaneobj, bpobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);

    for (int step = 0; step < steps; ++step){

        if(step > 0 && (step % frag_rest_iter) == 0)
        {
            if(bpoption) supbpfragmain(basepairlist);

            //cout<<"===================="<<endl;
            //cout<<hbondobj.size() << "-" << hbondbpobj.size()<<endl;
            MCresetNonbondedobj(gradptv, contactmap, bindingmap, basepairlist, noncbasepairlist, basepairscorelist, noncbasepairscorelist, hbondobj, hbondbpobj, enecoeff);
            //cout<<hbondobj.size() << "-" << hbondbpobj.size()<<endl;

            //cout<<"===================="<<endl;
            //cout<<stackobj.size() << "-" << bpdistobj.size()<<endl;
            if(enecoeff["Stack"] > 0){
                stackobj.clear();
                MCsetstackingobj(gradptv, stackingmap, enecoeff["Stack"], stackobj);
            }

            if(enecoeff["BPdist"] > 0){
                bpdistobj.clear();
                MCsetbpdistobj(gradptv, basepairlist, basepairscorelist, enecoeff["BPdist"], bpdistobj);
            }
            //cout<<stackobj.size() << "-" << bpdistobj.size()<<endl;

            reset_frag = true;
            s_list.clear();
            y_list.clear();
            rho_list.clear();
        }

        MCinitgrad(gradptv);
        MCgetgradient(blengthobj, bangleobj, dihedralobj, sugardihedralobj, baseplaneobj, bpobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);
        //averagenorm = calcaveragenorm(gradptv);
        
        if(!reset_frag){
            MCupdatelbfgsfeatures(refinestr_temp_prev, gradptv, gradptv_prev, s_list, y_list, rho_list, moveresidues, m);
            MClbfgsupdate(gradptv, s_list, y_list, rho_list, moveresidues, z); //lbfgs
        }


        for(int i=0;i<numbb;i++){
            if(moveresidues[i]){
                char base=refinestr_temp[i].aaa;
                for(int j=0;j<AccAtomParamNum[base][4];j++){
                    refinestr_temp_prev[i][j].x = refinestr_temp[i].ptv[j].x;
                    refinestr_temp_prev[i][j].y = refinestr_temp[i].ptv[j].y;
                    refinestr_temp_prev[i][j].z = refinestr_temp[i].ptv[j].z;
                }
            }
        }
        gradptv_prev=gradptv; //lbfgs

        //energy = goldenelectionsearch(superimposedstr, gradptv, blengthobj, bangleobj, bpobj, binitobj, telecteobj, ljobj, hbondobj, moveresidues);
        if(!reset_frag){
            averagenorm = MCcalcaveragenorm(z, 10000);
            //energy = goldenelectionsearch(superimposedstr, z, blengthobj, bangleobj, bpobj, binitobj, telecteobj, ljobj, hbondobj, moveresidues); //lbfgs
            double alpha0 = 1.0;
            double rho = 0.5;
            double beta = 0.1;
            energy = MClinesearch(z, blengthobj, bangleobj, dihedralobj, sugardihedralobj, baseplaneobj, bpobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj, moveresidues, alpha0, rho, beta);
        }
        else{
            averagenorm = MCcalcaveragenorm(gradptv, 10000);
            //cout<<"fff"<<endl;
            energy = MCgoldenelectionsearch(gradptv, blengthobj, bangleobj, dihedralobj, sugardihedralobj, baseplaneobj, bpobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj, moveresidues);
            //cout<<"ggg"<<endl;
            reset_frag = false;
        }

        //cout<< "epoch "<<epoch+1 << ": "<< energy << "=================="<<endl;

        if(((step+1) % frag_rest_iter) == 0){
            if(energy < minenergy){
                minenergy = energy;
                estop_counter = 0;
            }else{
                estop_counter += 1;
            }

            if(stop_count == estop_counter) break;
        }

        /*
        if(energy < minenergy){
            minenergy = energy;
            estop_counter = 0;
        }else{
            estop_counter += 1;
        }
        if(stop_count == estop_counter) break;
        */
    }

    double endene = MCgetcovalentenergystem(blengthobj, bangleobj, dihedralobj, sugardihedralobj, bpobj, baseplaneobj, bpdistobj, stackobj, stackgeoobj, binitobj, binitpairobj, btemplatepairobj, telecteobj, ljobj, hbondobj, hbondbpobj, clashobj, nttempobj);
    //exit(1);

}























