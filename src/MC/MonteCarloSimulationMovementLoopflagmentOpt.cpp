
#include "MonteCarloSimulation.hpp"

void MonteCarloSimulation::setloopntfragflap(int istart, int iend)
{
	int i,j;
	int atom1, atom2, atom3;

	for(i=istart;i<iend;i++){
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

		for(j=0;j<AddHydLengthParamsList[base].size();j++){
			hydaddlengthflag[i][j]=1;
		}

		for(j=0;j<AddHydAngleParamsList[base].size();j++)
		{
			hydaddangleflag[i][j]=1;
		}
	}

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

double MonteCarloSimulation::getbaseenergy(const int &nt_base, vector<double> &templateenergies, int &templateidx){
	int i,j;
	int atom;
	int count;
	vector<vector<double>> RotMatix;
    vector<double> TranVect;
	double minrmsd=100000000000;
	double rmsd;
	vector<double> xyztmp(3,0);
	Point3s vec;
	BasicFunc bf;

	char base;
	base = refinestr_temp[nt_base].aaa;

    static map<char, vector<int>> targetatoms={
		{'A',{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21}},
        {'C',{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19}},
        {'G',{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}},
        {'U',{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19}}
	};

	vector<vector<double> > xyz_list1(targetatoms[base].size(), vector<double>(3,0));
    vector<vector<double> > xyz_list2(targetatoms[base].size(), vector<double>(3,0));

	for(i=0;i<LOOPNTTEMPLATE[base].size();i++)
	{
		count = 0;
		for(j=0; j<targetatoms[base].size();j++){
			atom=targetatoms[base][j];

			xyz_list1[count][0]=refinestr_temp[nt_base].ptv[atom].x;
			xyz_list1[count][1]=refinestr_temp[nt_base].ptv[atom].y;
			xyz_list1[count][2]=refinestr_temp[nt_base].ptv[atom].z;
			xyz_list2[count][0]=LOOPNTTEMPLATE[base][i][0].ptv[atom].x;
			xyz_list2[count][1]=LOOPNTTEMPLATE[base][i][0].ptv[atom].y;
			xyz_list2[count][2]=LOOPNTTEMPLATE[base][i][0].ptv[atom].z;
			count+=1;
		}

		bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

		rmsd=0;
		count=0;
		for(j=0; j<targetatoms[base].size();j++){
			atom=targetatoms[base][j];
			bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);
			vec.x=refinestr_temp[nt_base].ptv[atom].x-xyztmp[0];
			vec.y=refinestr_temp[nt_base].ptv[atom].y-xyztmp[1];
			vec.z=refinestr_temp[nt_base].ptv[atom].z-xyztmp[2];
			rmsd+=bf.norm(vec);
			count+=1;
		}

		templateenergies.push_back(rmsd);

		if(rmsd < minrmsd){
			minrmsd=rmsd;
			templateidx=i;
		}
	}

	//cout<<minrmsd<<endl;

	return minrmsd/count;
}


bool MonteCarloSimulation::replacentfrag(int &base, int &selected_fragidx)
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

	char base1;
	base1 = refinestr_temp[base].aaa;

	//int selected_fragidx = int(Random()*LOOPNTTEMPLATE[base1].size());

	vector<vector<double>> xyz_list1(AccAtomParamNum[base1][3], vector<double>(3,0));
    vector<vector<double>> xyz_list2(AccAtomParamNum[base1][3], vector<double>(3,0));

	count = 0;
	for(j=0; j<AccAtomParamNum[base1][3];j++){
		xyz_list1[count][0]=refinestr_temp[base].ptv[j].x;
		xyz_list1[count][1]=refinestr_temp[base].ptv[j].y;
		xyz_list1[count][2]=refinestr_temp[base].ptv[j].z;

		xyz_list2[count][0]=LOOPNTTEMPLATE[base1][selected_fragidx][0].ptv[j].x;
		xyz_list2[count][1]=LOOPNTTEMPLATE[base1][selected_fragidx][0].ptv[j].y;
		xyz_list2[count][2]=LOOPNTTEMPLATE[base1][selected_fragidx][0].ptv[j].z;

		count+=1;
	}

	bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

	for(j=0; j<AccAtomParamNum[base1][3];j++){
		atom_coord.clear();
		atom_coord.resize(3, 0.0);
		atom_coord[0]=LOOPNTTEMPLATE[base1][selected_fragidx][0].ptv[j].x;
		atom_coord[1]=LOOPNTTEMPLATE[base1][selected_fragidx][0].ptv[j].y;
		atom_coord[2]=LOOPNTTEMPLATE[base1][selected_fragidx][0].ptv[j].z;

		bf.Transform(atom_coord,RotMatix,TranVect,xyztmp);

		refinestr_temp[base].ptv[j].x = xyztmp[0];
		refinestr_temp[base].ptv[j].y = xyztmp[1];
		refinestr_temp[base].ptv[j].z = xyztmp[2];

	}

	return true;
}

void MonteCarloSimulation::selectnttemplateidx(const int &nt_base, const vector<double> &templateenergies, int &templateidx)
{
	char base;
	base = refinestr_temp[nt_base].aaa;

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
		templateidx=int(Random()*LOOPNTTEMPLATE[base].size());
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

	if(templateidx > LOOPNTTEMPLATE[base].size()){
		templateidx=templateidx-LOOPNTTEMPLATE[base].size();
	}
}

int MonteCarloSimulation::selecttargetloopnt()
{
	int i,j;

	return int(Random()*looplist.size());
}

bool MonteCarloSimulation::fmovementloopfrag(double oldenergy, double &newenergy)
{
	int i, j;
	char base;

	if(looplist.size() == 0){
		newenergy=oldenergy+1000000000;
		return false;
	}

	int target_bp_idx = selecttargetloopnt();
	int nt_base = looplist[target_bp_idx];

	//cout<<"~~~~~~~~~~~~~~"<<endl;
	//cout<<nt_base<<endl;
	//cout<<"~~~~~~~~~~~~~~"<<endl;

	base = refinestr_temp[nt_base].aaa;

	vector<double> templateenergies;
	int templateidx;
	double randval=Random();

	StrCopy(refinestr, refinestr_temp);
	resetfrag(0);
	
	//if(randval < 0.5){	
	double minrmsd = getbaseenergy(nt_base, templateenergies, templateidx);
	
	//if(minrmsd < 0.1){
	//selectnttemplateidx(nt_base, templateenergies, templateidx);
		//}else{
			//templateidx = int(Random()*LOOPNTTEMPLATE[base].size());
		//}
	replacentfrag(nt_base, templateidx);
	param2coordhydrogen(refinestr_temp,nt_base,nt_base+1);
	//}

	map<int, vector<Point3s>> initstructure;
	initstructure[nt_base].resize(AccAtomParamNum[base][3]);

	for(j=0; j<AccAtomParamNum[base][3];j++){
		initstructure[nt_base][j].x=refinestr_temp[nt_base].ptv[j].x;
		initstructure[nt_base][j].y=refinestr_temp[nt_base].ptv[j].y;
		initstructure[nt_base][j].z=refinestr_temp[nt_base].ptv[j].z;
	}

	Point3s coordA;
	coordA.x = refinestr_temp[nt_base].ptv[8].x;
	coordA.y = refinestr_temp[nt_base].ptv[8].y;
	coordA.z = refinestr_temp[nt_base].ptv[8].z;

	vector<vector<int>> bindingmap(numbb, vector<int>(numbb, 0));
	if(nt_base - 1 >= 0) bindingmap[nt_base - 1][nt_base] = 1;
	if(nt_base + 1 < numbb) bindingmap[nt_base][nt_base + 1] = 1;
	bindingmap[nt_base][nt_base] = 2;

	vector<bool> moveresidues(numbb,false);
	moveresidues[nt_base]=true;

	int leng;
	vector<vector<int>> basepairlist;
	vector<vector<int>> noncbasepairlist;
	vector<double> basepairscorelist;
	vector<double> noncbasepairscorelist;

	vector<vector<int>> contactmap(numbb, vector<int>(numbb, 0));
	vector<vector<double>> stackingmap(numbb, vector<double>(numbb, 0));
	for(int i=0; i<numbb;i++){
		Point3s rp;
		if(i == nt_base) continue;
		rp.x = refinestr_temp[i].ptv[8].x - coordA.x;
		rp.y = refinestr_temp[i].ptv[8].y - coordA.y;
		rp.z = refinestr_temp[i].ptv[8].z - coordA.z;

		leng = rp.x*rp.x+rp.y*rp.y+rp.z*rp.z;

		//if(leng < GLOBALCUTOFF_SQR){
		contactmap[nt_base][i] = 1;
		contactmap[i][nt_base] = 1;
		//}

		if(nt_base == i){
			if(leng < CUTOFF_SQR+5.0){
				contactmap[nt_base][i] = 2;
				contactmap[i][nt_base] = 2;
			}
		}else{
			if(leng < CUTOFF_SQR+10.0){
				contactmap[nt_base][i] = 2;
				contactmap[i][nt_base] = 2;
			}
		}

		if(cbpmap[nt_base][i] > 0){ basepairlist.push_back({nt_base, i}); basepairscorelist.push_back(cbpmap[nt_base][i]);}
		if(noncbpmap[nt_base][i] > 0){ noncbasepairlist.push_back({nt_base, i}); noncbasepairscorelist.push_back(noncbpmap[nt_base][i]);}

		stackingmap[target_bp_idx][i] = stmap[target_bp_idx][i];
		stackingmap[i][target_bp_idx] = stmap[i][target_bp_idx];
	}

	MonteCarloSimulation::MCfragmentoptimization(initstructure, bindingmap, moveresidues, contactmap, stackingmap, basepairlist, noncbasepairlist, basepairscorelist, noncbasepairscorelist);

	
	int istart=nt_base-1;
	int iend=nt_base+2;
	if(istart<0) istart=0;
	if(iend>=numbb) iend=numbb;

	coord2parambb(refinestr_temp,istart,iend);
	coord2paramsugar(refinestr_temp,istart,iend);
	coord2parambase(refinestr_temp,istart,iend);
	coord2paramhydrogen(refinestr_temp,istart,iend);
	coord2addparamhydrogen(refinestr_temp,istart,iend);
	coord2paramadd(refinestr_temp,istart,iend);

	setloopntfragflap(istart, iend);

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

