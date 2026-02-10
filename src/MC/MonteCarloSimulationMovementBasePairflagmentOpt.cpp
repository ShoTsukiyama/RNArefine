
#include "MonteCarloSimulation.hpp"

void MonteCarloSimulation::setbpfragflap(int istart1, int iend1, int istart2, int iend2)
{
	int i,j;
	int atom1, atom2, atom3;

	for(i=istart1;i<iend1;i++){
		char base=refinestr_temp[i].aaa;

        for(j=0;j<MCLENGTHPARAM[base].size();j++){
            lengthflag[i][j] = 1;
        }

        for(j=0;j<MCINNERANGPARAM[base].size();j++){
            angleflag[i][j] = 1;
        }

        for(j=0;j<MCTORSIONPARAM[base].size();j++){
            torsionflag[i][j] = 1;
        }

        for(j=0;j<MCHYDLENGTHPARAM[base].size();j++){
            hydlengthflag[i][j] = 1;
        }

        for(j=0;j<MCHYDINNERANGPARAM[base].size();j++){
            hydangleflag[i][j] = 1;
        }

        for(j=0;j<MCHYDTORSIONPARAM[base].size();j++){
            hydtorsionflag[i][j] = 1;
        }

        for(j=0;j<AddLengthParamsList[base].size();j++){
            addlengthflag[i][j] = 1;
        }

        for(j=0;j<AddAngleParamsList[base].size();j++){
            addangleflag[i][j] = 1;
        }

		for(j=0;j<AddHydLengthParamsList[base].size();j++)
		{
			hydaddlengthflag[i][j]=1;
		}

		for(j=0;j<AddHydAngleParamsList[base].size();j++)
		{
			hydaddangleflag[i][j]=1;
		}
	}

	for(i=istart2;i<iend2;i++){
		char base=refinestr_temp[i].aaa;

        for(j=0;j<MCLENGTHPARAM[base].size();j++){
            lengthflag[i][j] = 1;
        }

        for(j=0;j<MCINNERANGPARAM[base].size();j++){
            angleflag[i][j] = 1;
        }

        for(j=0;j<MCTORSIONPARAM[base].size();j++){
            torsionflag[i][j] = 1;
        }

        for(j=0;j<MCHYDLENGTHPARAM[base].size();j++){
            hydlengthflag[i][j] = 1;
        }

        for(j=0;j<MCHYDINNERANGPARAM[base].size();j++){
            hydangleflag[i][j] = 1;
        }

        for(j=0;j<MCHYDTORSIONPARAM[base].size();j++){
            hydtorsionflag[i][j] = 1;
        }

        for(j=0;j<AddLengthParamsList[base].size();j++){
            addlengthflag[i][j] = 1;
        }

        for(j=0;j<AddAngleParamsList[base].size();j++){
            addangleflag[i][j] = 1;
        }

		for(j=0;j<AddHydLengthParamsList[base].size();j++)
		{
			hydaddlengthflag[i][j]=1;
		}

		for(j=0;j<AddHydAngleParamsList[base].size();j++)
		{
			hydaddangleflag[i][j]=1;
		}
	}

	int istart = min(istart1, istart2);
	int iend = max(iend1, iend2);

	for(i=0;i<numbb;i++)
	{
		ntflag[i] = 1;
		for(j=i;j<numbb;j++)
		{
			if(j<istart || i>iend){}
			else
			{
				mapbbflag[i][j]=1;
				mapbbflag[j][i]=1;
				mapflag[i][j]=1;
				mapflag[j][i]=1;
			}
		}
	}
}

//double MonteCarloSimulation::getbasepairenergy(const int &bp_base1, const int &bp_base2, vector<double> &templateenergies, int &templateidx){
double MonteCarloSimulation::getbasepairenergy(const int &bp_base1, const int &bp_base2, int &templateidx){
	int i,j;
	int atom;
	int count;
	vector<vector<double>> RotMatix;
    vector<double> TranVect;
	double minrmsd=100000000000;
	double rmsd;
	vector<double> xyztmp(3,0);
	vector<double> atom_coord;
	Point3s vec;
	BasicFunc bf;

	char base1, base2;
	base1 = refinestr_temp[bp_base1].aaa;
	base2 = refinestr_temp[bp_base2].aaa;

	/*
	static map<char, vector<int>> targetatoms_sup={
		{'A',{0,1,2,3,4,5}},
        {'C',{0,1,2,3,4,5}},
        {'G',{0,1,2,3,4,5}},
        {'U',{0,1,2,3,4,5}}
	};
	*/
	
	static map<char, vector<int>> targetatoms_sup={
		{'A',{12,13,14,15,16,17,18,19,20,21}},
		{'C',{12,13,14,15,16,17,18,19}},
		{'G',{12,13,14,15,16,17,18,19,20,21,22}},
		{'U',{12,13,14,15,16,17,18,19}}
	};
	

	vector<vector<double> > xyz_list1(targetatoms_sup[base1].size()+targetatoms_sup[base2].size(), vector<double>(3,0));
    vector<vector<double> > xyz_list2(targetatoms_sup[base1].size()+targetatoms_sup[base2].size(), vector<double>(3,0));

	for(i=0;i<BPTEMPLATE[base1][base2].size();i++)
	{
		count = 0;
		for(j=0; j<targetatoms_sup[base1].size();j++){
			atom=targetatoms_sup[base1][j];

			xyz_list1[count][0]=refinestr_temp[bp_base1].ptv[atom].x;
			xyz_list1[count][1]=refinestr_temp[bp_base1].ptv[atom].y;
			xyz_list1[count][2]=refinestr_temp[bp_base1].ptv[atom].z;

			xyz_list2[count][0]=BPTEMPLATE[base1][base2][i][0].ptv[atom].x;
			xyz_list2[count][1]=BPTEMPLATE[base1][base2][i][0].ptv[atom].y;
			xyz_list2[count][2]=BPTEMPLATE[base1][base2][i][0].ptv[atom].z;

			count+=1;
		}

		for(j=0; j<targetatoms_sup[base2].size();j++){
			atom=targetatoms_sup[base2][j];

			xyz_list1[count][0]=refinestr_temp[bp_base2].ptv[atom].x;
			xyz_list1[count][1]=refinestr_temp[bp_base2].ptv[atom].y;
			xyz_list1[count][2]=refinestr_temp[bp_base2].ptv[atom].z;

			xyz_list2[count][0]=BPTEMPLATE[base1][base2][i][1].ptv[atom].x;
			xyz_list2[count][1]=BPTEMPLATE[base1][base2][i][1].ptv[atom].y;
			xyz_list2[count][2]=BPTEMPLATE[base1][base2][i][1].ptv[atom].z;

			count+=1;
		}

		bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

		rmsd=0;
		count=0;
		for(j=0; j<targetatoms_sup[base1].size();j++){
			atom=targetatoms_sup[base1][j];

			atom_coord.clear();
			atom_coord.resize(3, 0.0);
			atom_coord[0]=BPTEMPLATE[base1][base2][i][0].ptv[atom].x;
			atom_coord[1]=BPTEMPLATE[base1][base2][i][0].ptv[atom].y;
			atom_coord[2]=BPTEMPLATE[base1][base2][i][0].ptv[atom].z;

			bf.Transform(atom_coord,RotMatix,TranVect,xyztmp);
			vec.x=refinestr_temp[bp_base1].ptv[atom].x-xyztmp[0];
			vec.y=refinestr_temp[bp_base1].ptv[atom].y-xyztmp[1];
			vec.z=refinestr_temp[bp_base1].ptv[atom].z-xyztmp[2];

			rmsd+=(vec.x*vec.x+vec.y*vec.y+vec.z*vec.z);
			count+=1;
		}

		for(j=0; j<targetatoms_sup[base2].size();j++){
			atom=targetatoms_sup[base2][j];

			atom_coord.clear();
			atom_coord.resize(3, 0.0);
			atom_coord[0]=BPTEMPLATE[base1][base2][i][1].ptv[atom].x;
			atom_coord[1]=BPTEMPLATE[base1][base2][i][1].ptv[atom].y;
			atom_coord[2]=BPTEMPLATE[base1][base2][i][1].ptv[atom].z;

			bf.Transform(atom_coord,RotMatix,TranVect,xyztmp);
			vec.x=refinestr_temp[bp_base2].ptv[atom].x-xyztmp[0];
			vec.y=refinestr_temp[bp_base2].ptv[atom].y-xyztmp[1];
			vec.z=refinestr_temp[bp_base2].ptv[atom].z-xyztmp[2];

			rmsd+=(vec.x*vec.x+vec.y*vec.y+vec.z*vec.z);
			count+=1;
		}
		rmsd = sqrt(rmsd/count);

		if(rmsd < minrmsd){
			templateidx=i;
			minrmsd=rmsd;
		}
	}

	return minrmsd;
}


bool MonteCarloSimulation::replacebpfrag(int &bp_base1, int &bp_base2, int &selected_fragidx)
{
	int i,j;
	int atom;
	int count;
	vector<vector<double>> RotMatix;
    vector<double> TranVect;
	double minrmsd=100000000000;
	double rmsd;
	vector<double> xyztmp(3,0);
	vector<double> atom_coord;
	BasicFunc bf;

	char base1, base2;
	base1 = refinestr_temp[bp_base1].aaa;
	base2 = refinestr_temp[bp_base2].aaa;

	/*	
	static map<char, vector<int>> targetatoms_sup={
		{'A',{0,1,2,3,4,5}},
		{'C',{0,1,2,3,4,5}},
		{'G',{0,1,2,3,4,5}},
		{'U',{0,1,2,3,4,5}}
	};
	
	static map<char, vector<int>> targetatoms={
		{'A',{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21}},
		{'C',{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19}},
		{'G',{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}},
		{'U',{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19}}
	};
	*/
	
	static map<char, vector<int>> targetatoms_sup={
		{'A',{12,13,14,15,16,17,18,19,20,21}},
		{'C',{12,13,14,15,16,17,18,19}},
		{'G',{12,13,14,15,16,17,18,19,20,21,22}},
		{'U',{12,13,14,15,16,17,18,19}}
	};
	
	static map<char, vector<int>> targetatoms={
		{'A',{12,13,14,15,16,17,18,19,20,21}},
		{'C',{12,13,14,15,16,17,18,19}},
		{'G',{12,13,14,15,16,17,18,19,20,21,22}},
		{'U',{12,13,14,15,16,17,18,19}}
	};
	
	
	vector<vector<double>> xyz_list1(targetatoms_sup[base1].size()+targetatoms_sup[base2].size(), vector<double>(3,0));
    vector<vector<double>> xyz_list2(targetatoms_sup[base1].size()+targetatoms_sup[base2].size(), vector<double>(3,0));

	count = 0;
	for(j=0; j<targetatoms_sup[base1].size();j++){
		atom = targetatoms_sup[base1][j];
		xyz_list1[count][0]=refinestr_temp[bp_base1].ptv[atom].x;
		xyz_list1[count][1]=refinestr_temp[bp_base1].ptv[atom].y;
		xyz_list1[count][2]=refinestr_temp[bp_base1].ptv[atom].z;

		xyz_list2[count][0]=BPTEMPLATE[base1][base2][selected_fragidx][0].ptv[atom].x;
		xyz_list2[count][1]=BPTEMPLATE[base1][base2][selected_fragidx][0].ptv[atom].y;
		xyz_list2[count][2]=BPTEMPLATE[base1][base2][selected_fragidx][0].ptv[atom].z;

		count+=1;
	}

	for(j=0; j<targetatoms_sup[base2].size();j++){
		atom = targetatoms_sup[base2][j];
		xyz_list1[count][0]=refinestr_temp[bp_base2].ptv[atom].x;
		xyz_list1[count][1]=refinestr_temp[bp_base2].ptv[atom].y;
		xyz_list1[count][2]=refinestr_temp[bp_base2].ptv[atom].z;

		xyz_list2[count][0]=BPTEMPLATE[base1][base2][selected_fragidx][1].ptv[atom].x;
		xyz_list2[count][1]=BPTEMPLATE[base1][base2][selected_fragidx][1].ptv[atom].y;
		xyz_list2[count][2]=BPTEMPLATE[base1][base2][selected_fragidx][1].ptv[atom].z;
		count+=1;
	}

	bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

	for(j=0; j<targetatoms[base1].size();j++){
		atom = targetatoms[base1][j];
		atom_coord.clear();
		atom_coord.resize(3, 0.0);
		atom_coord[0]=BPTEMPLATE[base1][base2][selected_fragidx][0].ptv[atom].x;
		atom_coord[1]=BPTEMPLATE[base1][base2][selected_fragidx][0].ptv[atom].y;
		atom_coord[2]=BPTEMPLATE[base1][base2][selected_fragidx][0].ptv[atom].z;

		bf.Transform(atom_coord,RotMatix,TranVect,xyztmp);

		refinestr_temp[bp_base1].ptv[atom].x = xyztmp[0];
		refinestr_temp[bp_base1].ptv[atom].y = xyztmp[1];
		refinestr_temp[bp_base1].ptv[atom].z = xyztmp[2];
	}

	for(j=0; j<targetatoms[base2].size();j++){
		atom = targetatoms[base2][j];
		atom_coord.clear();
		atom_coord.resize(3, 0.0);
		atom_coord[0]=BPTEMPLATE[base1][base2][selected_fragidx][1].ptv[atom].x;
		atom_coord[1]=BPTEMPLATE[base1][base2][selected_fragidx][1].ptv[atom].y;
		atom_coord[2]=BPTEMPLATE[base1][base2][selected_fragidx][1].ptv[atom].z;

		bf.Transform(atom_coord,RotMatix,TranVect,xyztmp);

		refinestr_temp[bp_base2].ptv[atom].x = xyztmp[0];
		refinestr_temp[bp_base2].ptv[atom].y = xyztmp[1];
		refinestr_temp[bp_base2].ptv[atom].z = xyztmp[2];
	}

	return true;
}


int MonteCarloSimulation::selecttargetbp(int &basei, int &basej, double &bpscore)
{
	if (bpscoresum <= 0.0) return 0;
	int i,j;

	double totalval = 0.0;
	double randval = Random()*bpscoresum;

	for(int i=0;i<numbb;i++)
	{
		for(int j=i+1;j<numbb;j++)
		{
			if(cbpmap[i][j] > 0)
			{
				if(totalval <= randval && randval < totalval+cbpmap[i][j]){
					basei = i;
					basej = j;
					bpscore = cbpmap[i][j];
					return 1;
				}
				totalval += cbpmap[i][j];
			}

			/*
			if(noncbpmap[i][j] > 0)
			{
				if(totalval <= randval && randval < totalval+noncbpmap[i][j]){
					basei = i;
					basej = j;
					bpscore = noncbpmap[i][j];
					return 2;
				}
				totalval += noncbpmap[i][j];
			}
			*/
		}
	}
	return 0;
}

void MonteCarloSimulation::selectbptemplateidx(const int &bp_base1, const int &bp_base2, const vector<double> &templateenergies, int &templateidx)
{
	char base1, base2;
	base1 = refinestr_temp[bp_base1].aaa;
	base2 = refinestr_temp[bp_base2].aaa;

	vector<double> templateenergies_norm(templateenergies.size());

	int maxval=20;
	int minval=0;

	double totalval = 0;
	double tempenergy;
	for(int i=0;i<templateenergies.size();i++){
		tempenergy = min(20.0, templateenergies[i]);
		templateenergies_norm[i] = 1 - ((tempenergy - minval)/(maxval-minval));
		totalval += templateenergies_norm[i];
	}

	if(totalval <= 0){
		templateidx=int(Random()*BPTEMPLATE[base1][base2].size());
		return;
	}

	double randval = Random()*totalval;
	double val=0;
	
	templateidx = -1;
	for(int i=0;i<templateenergies_norm.size();i++)
	{
		if(val <= randval && randval < val+templateenergies_norm[i]){
			templateidx = i;
			break;
		}
		val += templateenergies_norm[i];
	}

	if(templateidx < 0) templateidx=templateenergies_norm.size()-1;

	if(templateidx > BPTEMPLATE[base1][base2].size()){
		templateidx=templateidx-BPTEMPLATE[base1][base2].size();
	}
}

bool MonteCarloSimulation::fmovementbpfrag(double oldenergy, double &newenergy)
{
	int i, j;
	char base1, base2;

	if(bplist.size() == 0){
		newenergy=oldenergy+1000000000;
		return false;
	}

	//int target_bp_idx = selecttargetbp();
	//int bp_base1 = bplist[target_bp_idx][0];
	//int bp_base2 = bplist[target_bp_idx][1];
	vector<vector<int>> basepairlist;
	vector<vector<int>> noncbasepairlist;
	vector<double> basepairscorelist;
	vector<double> noncbasepairscorelist;

	int bp_base1, bp_base2;
	double bpscore=-1;
	int bpselflag = selecttargetbp(bp_base1, bp_base2, bpscore);

	if(bpselflag == 0){
		newenergy=oldenergy+1000000000;
		return false;
	}else if(bpselflag == 1){
		basepairlist.push_back({bp_base1, bp_base2});
		basepairscorelist.push_back(bpscore);
	}else if(bpselflag == 2){
		noncbasepairlist.push_back({bp_base1, bp_base2});
		noncbasepairscorelist.push_back(bpscore);
	}

	StrCopy(refinestr, refinestr_temp);
	resetfrag(0);

	base1 = refinestr_temp[bp_base1].aaa;
	base2 = refinestr_temp[bp_base2].aaa;

	//cout<<"~~~~~~~~~~~~~~"<<endl;
	//cout<<bp_base1 << "-" << bp_base2 << "-" << base1 << "-" << base2<<endl;
	//cout<<"~~~~~~~~~~~~~~"<<endl;

	//vector<double> templateenergies;
	int templateidx;
	//double minrmsd = getbasepairenergy(bp_base1, bp_base2, templateenergies, templateidx);
	double minrmsd = getbasepairenergy(bp_base1, bp_base2, templateidx);

	cout<<"Base pair-" << minrmsd << "-" << bpscore <<endl;

	if(minrmsd > 20.0){
		newenergy=oldenergy+1000000000;
		return false;
	}

	//cout<<bpscore << "-" << minrmsd<<endl;

	//if(minrmsd < 0.1){
	//selectbptemplateidx(bp_base1, bp_base2, templateenergies, templateidx);
	replacebpfrag(bp_base1, bp_base2, templateidx);
	param2coordhydrogen(refinestr_temp,bp_base1,bp_base1+1);
	param2coordhydrogen(refinestr_temp,bp_base2,bp_base2+1);
	//}

	map<int, vector<Point3s>> initstructure;
	initstructure[bp_base1].resize(AccAtomParamNum[base1][3]);
	initstructure[bp_base2].resize(AccAtomParamNum[base2][3]);

	for(j=0; j<AccAtomParamNum[base1][3];j++){
		initstructure[bp_base1][j].x=refinestr_temp[bp_base1].ptv[j].x;
		initstructure[bp_base1][j].y=refinestr_temp[bp_base1].ptv[j].y;
		initstructure[bp_base1][j].z=refinestr_temp[bp_base1].ptv[j].z;
	}
	for(j=0; j<AccAtomParamNum[base2][3];j++){
		initstructure[bp_base2][j].x=refinestr_temp[bp_base2].ptv[j].x;
		initstructure[bp_base2][j].y=refinestr_temp[bp_base2].ptv[j].y;
		initstructure[bp_base2][j].z=refinestr_temp[bp_base2].ptv[j].z;
	}

	Point3s coordA, coordB;
	coordA.x = refinestr_temp[bp_base1].ptv[8].x;
	coordA.y = refinestr_temp[bp_base1].ptv[8].y;
	coordA.z = refinestr_temp[bp_base1].ptv[8].z;
	coordB.x = refinestr_temp[bp_base2].ptv[8].x;
	coordB.y = refinestr_temp[bp_base2].ptv[8].y;
	coordB.z = refinestr_temp[bp_base2].ptv[8].z;

	vector<vector<int>> bindingmap(numbb, vector<int>(numbb, 0));
	if(bp_base1 - 1 >= 0) bindingmap[bp_base1 - 1][bp_base1] = 1;
	if(bp_base2 - 1 >= 0) bindingmap[bp_base2 - 1][bp_base2] = 1;
	if(bp_base1 + 1 < numbb) bindingmap[bp_base1][bp_base1 + 1] = 1;
	if(bp_base2 + 1 < numbb) bindingmap[bp_base2][bp_base2 + 1] = 1;
	bindingmap[bp_base1][bp_base1] = 2;
	bindingmap[bp_base2][bp_base2] = 2;

	vector<bool> moveresidues(numbb,false);
	moveresidues[bp_base1]=true;
	moveresidues[bp_base2]=true;

	int leng;
	vector<vector<int>> contactmap(numbb, vector<int>(numbb, 0));
	vector<vector<double>> stackingmap(numbb, vector<double>(numbb, 0));
	for(int i=0; i<numbb;i++){
		Point3s rp;
		if(i == bp_base1 || i == bp_base2) continue;
		rp.x = refinestr_temp[i].ptv[8].x - coordA.x;
		rp.y = refinestr_temp[i].ptv[8].y - coordA.y;
		rp.z = refinestr_temp[i].ptv[8].z - coordA.z;

		leng = rp.x*rp.x+rp.y*rp.y+rp.z*rp.z;

		//if(leng < GLOBALCUTOFF_SQR){
		contactmap[bp_base1][i] = 1;
		contactmap[i][bp_base1] = 1;
		//}

		if(bp_base1 == i){
			if(leng < CUTOFF_SQR+5.0){
				contactmap[bp_base1][i] = 2;
				contactmap[i][bp_base1] = 2;
			}
		}else{
			if(leng < CUTOFF_SQR+10.0){
				contactmap[bp_base1][i] = 2;
				contactmap[i][bp_base1] = 2;
			}
		}
		

		rp.x = refinestr_temp[i].ptv[8].x - coordB.x;
		rp.y = refinestr_temp[i].ptv[8].y - coordB.y;
		rp.z = refinestr_temp[i].ptv[8].z - coordB.z;

		leng = rp.x*rp.x+rp.y*rp.y+rp.z*rp.z;

		//if(leng < GLOBALCUTOFF_SQR){
		contactmap[bp_base2][i] = 1;
		contactmap[i][bp_base2] = 1;
		//}

		if(bp_base2 == i){
			if(leng < CUTOFF_SQR+5.0){
				contactmap[bp_base2][i] = 2;
				contactmap[i][bp_base2] = 2;
			}
		}else{
			if(leng < CUTOFF_SQR+10.0){
				contactmap[bp_base2][i] = 2;
				contactmap[i][bp_base2] = 2;
			}
		}

		stackingmap[bp_base1][i] = stmap[bp_base1][i];
		stackingmap[i][bp_base1] = stmap[i][bp_base1];
		stackingmap[bp_base2][i] = stmap[bp_base2][i];
		stackingmap[i][bp_base2] = stmap[i][bp_base2];
	}

	//basepairlist.push_back({bp_base1, bp_base2});

	MonteCarloSimulation::MCfragmentoptimization(initstructure, bindingmap, moveresidues, contactmap, stackingmap, basepairlist, noncbasepairlist, basepairscorelist, noncbasepairscorelist, true);

	
	int istart1=bp_base1-1;
	int iend1=bp_base1+2;
	if(istart1<0) istart1=0;
	if(iend1>=numbb) iend1=numbb;

	int istart2=bp_base2-1;
	int iend2=bp_base2+2;
	if(istart2<0) istart2=0;
	if(iend2>=numbb) iend2=numbb;

	coord2parambb(refinestr_temp,istart1,iend1);
	coord2paramsugar(refinestr_temp,istart1,iend1);
	coord2parambase(refinestr_temp,istart1,iend1);
	coord2paramhydrogen(refinestr_temp,istart1,iend1);
	coord2addparamhydrogen(refinestr_temp,istart1,iend1);
	coord2paramadd(refinestr_temp,istart1,iend1);

	coord2parambb(refinestr_temp,istart2,iend2);
	coord2paramsugar(refinestr_temp,istart2,iend2);
	coord2parambase(refinestr_temp,istart2,iend2);
	coord2paramhydrogen(refinestr_temp,istart2,iend2);
	coord2addparamhydrogen(refinestr_temp,istart2,iend2);
	coord2paramadd(refinestr_temp,istart2,iend2);

	setbpfragflap(istart1, iend1, istart2, iend2);


	newenergy = calcenergymain();
	bool upflag=metropolisUpdate(oldenergy, newenergy);
	if(upflag){
		StrCopy(refinestr_temp, refinestr);
		initializeatomenergy();
		updateenergies();
	}

	/*
	cout<<"===================="<<endl;
	cout<<"QqqqqqqqqqqqqqqqqQ"<<endl;
	checkupdate();
	cout<<"PppppppppppppppppP"<<endl;
	resetfrag(1);
	double newenergy_tmp = calcenergymain();

	if(abs(newenergy - newenergy_tmp)/max(abs(newenergy), abs(newenergy_tmp)) > 0.002){
		cout<<"Energy Error"<<endl;
		cout<<newenergy << "-" << newenergy_tmp<<endl;
		exit(1);
	}else{
		cout<<"OK"<<endl;
	}
	*/

	return upflag;
}

