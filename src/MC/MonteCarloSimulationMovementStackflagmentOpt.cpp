
#include "MonteCarloSimulation.hpp"

void MonteCarloSimulation::setadjstackfragflap(int istart, int iend)
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

void MonteCarloSimulation::setnonadjstackfragflap(int istart1, int iend1, int istart2, int iend2)
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

		for(j=0;j<AddHydLengthParamsList[base].size();j++){
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

		for(j=0;j<AddHydLengthParamsList[base].size();j++){
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

double MonteCarloSimulation::getstackenergy(const int &stack_base1, const int &stack_base2, const int &sttype, vector<double> &templateenergies, int &mintempidx, int &minorder){
	int i,j;
	int atom;
	int count;
	vector<vector<double>> RotMatix;
    vector<double> TranVect;
	double minrmsd=1.0/0.0;
	double rmsd;
	vector<double> xyztmp(3,0);
	vector<double> atom_coord;
	Point3s vec;
	BasicFunc bf;

	char base1, base2;
	base1 = refinestr_temp[stack_base1].aaa;
	base2 = refinestr_temp[stack_base2].aaa;

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
	
	mintempidx = -1;

	/*
	static map<char, vector<int>> targetatoms={
		{'A',{12,13,14,15,16,17,18,19,20,21}},
        {'C',{12,13,14,15,16,17,18,19}},
        {'G',{12,13,14,15,16,17,18,19,20,21,22}},
        {'U',{12,13,14,15,16,17,18,19}}
	};
	*/
	vector<vector<double> > xyz_list1(targetatoms_sup[base1].size()+targetatoms_sup[base2].size(), vector<double>(3,0));
    vector<vector<double> > xyz_list2(targetatoms_sup[base1].size()+targetatoms_sup[base2].size(), vector<double>(3,0));

	for(i=0;i<STACKTEMPLATE[sttype][base1][base2].size();i++)
	{
		count = 0;
		for(j=0; j<targetatoms_sup[base1].size();j++){
			atom=targetatoms_sup[base1][j];

			xyz_list1[count][0]=refinestr_temp[stack_base1].ptv[atom].x;
			xyz_list1[count][1]=refinestr_temp[stack_base1].ptv[atom].y;
			xyz_list1[count][2]=refinestr_temp[stack_base1].ptv[atom].z;

			xyz_list2[count][0]=STACKTEMPLATE[sttype][base1][base2][i][0].ptv[atom].x;
			xyz_list2[count][1]=STACKTEMPLATE[sttype][base1][base2][i][0].ptv[atom].y;
			xyz_list2[count][2]=STACKTEMPLATE[sttype][base1][base2][i][0].ptv[atom].z;

			count+=1;
		}

		for(j=0; j<targetatoms_sup[base2].size();j++){
			atom=targetatoms_sup[base2][j];

			xyz_list1[count][0]=refinestr_temp[stack_base2].ptv[atom].x;
			xyz_list1[count][1]=refinestr_temp[stack_base2].ptv[atom].y;
			xyz_list1[count][2]=refinestr_temp[stack_base2].ptv[atom].z;

			xyz_list2[count][0]=STACKTEMPLATE[sttype][base1][base2][i][1].ptv[atom].x;
			xyz_list2[count][1]=STACKTEMPLATE[sttype][base1][base2][i][1].ptv[atom].y;
			xyz_list2[count][2]=STACKTEMPLATE[sttype][base1][base2][i][1].ptv[atom].z;

			count+=1;
		}

		bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

		rmsd=0;
		count=0;
		for(j=0; j<targetatoms[base1].size();j++){
			atom=targetatoms[base1][j];

			atom_coord.clear();
			atom_coord.resize(3, 0.0);
			atom_coord[0]=STACKTEMPLATE[sttype][base1][base2][i][0].ptv[atom].x;
			atom_coord[1]=STACKTEMPLATE[sttype][base1][base2][i][0].ptv[atom].y;
			atom_coord[2]=STACKTEMPLATE[sttype][base1][base2][i][0].ptv[atom].z;

			bf.Transform(atom_coord,RotMatix,TranVect,xyztmp);

			vec.x=refinestr_temp[stack_base1].ptv[atom].x-xyztmp[0];
			vec.y=refinestr_temp[stack_base1].ptv[atom].y-xyztmp[1];
			vec.z=refinestr_temp[stack_base1].ptv[atom].z-xyztmp[2];

			rmsd+=bf.norm(vec);
			count+=1;
		}

		for(j=0; j<targetatoms[base2].size();j++){
			atom=targetatoms[base2][j];

			atom_coord.clear();
			atom_coord.resize(3, 0.0);
			atom_coord[0]=STACKTEMPLATE[sttype][base1][base2][i][1].ptv[atom].x;
			atom_coord[1]=STACKTEMPLATE[sttype][base1][base2][i][1].ptv[atom].y;
			atom_coord[2]=STACKTEMPLATE[sttype][base1][base2][i][1].ptv[atom].z;

			bf.Transform(atom_coord,RotMatix,TranVect,xyztmp);
			vec.x=refinestr_temp[stack_base2].ptv[atom].x-xyztmp[0];
			vec.y=refinestr_temp[stack_base2].ptv[atom].y-xyztmp[1];
			vec.z=refinestr_temp[stack_base2].ptv[atom].z-xyztmp[2];

			rmsd+=bf.norm(vec);
			count+=1;
		}
		templateenergies.push_back(rmsd);

		if(rmsd < minrmsd){
			mintempidx=i;
			minorder=0;
			minrmsd=rmsd;
		}
	}

	if(sttype == 0 && base1 == base2)
	{
		for(i=0;i<STACKTEMPLATE[sttype][base1][base2].size();i++)
		{
			int count = 0;
			for(j=0; j<targetatoms_sup[base2].size();j++){
				atom=targetatoms_sup[base2][j];

				xyz_list1[count][0]=refinestr_temp[stack_base2].ptv[atom].x;
				xyz_list1[count][1]=refinestr_temp[stack_base2].ptv[atom].y;
				xyz_list1[count][2]=refinestr_temp[stack_base2].ptv[atom].z;

				xyz_list2[count][0]=STACKTEMPLATE[sttype][base1][base2][i][0].ptv[atom].x;
				xyz_list2[count][1]=STACKTEMPLATE[sttype][base1][base2][i][0].ptv[atom].y;
				xyz_list2[count][2]=STACKTEMPLATE[sttype][base1][base2][i][0].ptv[atom].z;

				count+=1;
			}

			for(j=0; j<targetatoms_sup[base1].size();j++){
				atom=targetatoms_sup[base1][j];

				xyz_list1[count][0]=refinestr_temp[stack_base1].ptv[atom].x;
				xyz_list1[count][1]=refinestr_temp[stack_base1].ptv[atom].y;
				xyz_list1[count][2]=refinestr_temp[stack_base1].ptv[atom].z;

				xyz_list2[count][0]=STACKTEMPLATE[sttype][base1][base2][i][1].ptv[atom].x;
				xyz_list2[count][1]=STACKTEMPLATE[sttype][base1][base2][i][1].ptv[atom].y;
				xyz_list2[count][2]=STACKTEMPLATE[sttype][base1][base2][i][1].ptv[atom].z;

				count+=1;
			}

			bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

			rmsd=0;
			count=0;
			//for(int j=0; j<SNTTEMPLATEATOMS[aat].size();j++){
			for(j=0; j<targetatoms[base2].size();j++){
				atom=targetatoms[base2][j];

				atom_coord.clear();
				atom_coord.resize(3, 0.0);
				atom_coord[0]=STACKTEMPLATE[sttype][base1][base2][i][0].ptv[atom].x;
				atom_coord[1]=STACKTEMPLATE[sttype][base1][base2][i][0].ptv[atom].y;
				atom_coord[2]=STACKTEMPLATE[sttype][base1][base2][i][0].ptv[atom].z;

				bf.Transform(atom_coord,RotMatix,TranVect,xyztmp);
				vec.x=refinestr_temp[stack_base2].ptv[atom].x-xyztmp[0];
				vec.y=refinestr_temp[stack_base2].ptv[atom].y-xyztmp[1];
				vec.z=refinestr_temp[stack_base2].ptv[atom].z-xyztmp[2];

				rmsd+=bf.norm(vec);
				count+=1;
			}

			for(j=0; j<targetatoms[base1].size();j++){
				atom=targetatoms[base1][j];

				atom_coord.clear();
				atom_coord.resize(3, 0.0);
				atom_coord[0]=STACKTEMPLATE[sttype][base1][base2][i][1].ptv[atom].x;
				atom_coord[1]=STACKTEMPLATE[sttype][base1][base2][i][1].ptv[atom].y;
				atom_coord[2]=STACKTEMPLATE[sttype][base1][base2][i][1].ptv[atom].z;

				bf.Transform(atom_coord,RotMatix,TranVect,xyztmp);
				vec.x=refinestr_temp[stack_base1].ptv[atom].x-xyztmp[0];
				vec.y=refinestr_temp[stack_base1].ptv[atom].y-xyztmp[1];
				vec.z=refinestr_temp[stack_base1].ptv[atom].z-xyztmp[2];

				rmsd+=bf.norm(vec);
				count+=1;
			}

			templateenergies.push_back(rmsd);

			if(rmsd < minrmsd){
				mintempidx=i;
				minorder=1;
				minrmsd=rmsd;
			}
		}
	}

	return minrmsd;
}

void MonteCarloSimulation::selectstacktemplateidx(const int &stack_base1, const int &stack_base2, const int &sttype, const vector<double> &templateenergies, int &templateidx, int &temporder)
{
	char base1, base2;
	base1 = refinestr_temp[stack_base1].aaa;
	base2 = refinestr_temp[stack_base2].aaa;

	vector<double> templateenergies_norm(templateenergies.size());

	/*
	int maxval=-10000;
	int minval=100000;
	for(int i=0;i<templateenergies.size();i++){
		if(maxval < templateenergies[i]) maxval = templateenergies[i];
		if(minval > templateenergies[i]) minval = templateenergies[i];
	}
	*/
	int maxval=20;
	int minval=0;
	temporder=0;

	double totalval = 0;
	double tempenergy;
	for(int i=0;i<templateenergies.size();i++){
		tempenergy = min(20.0, templateenergies[i]);
		templateenergies_norm[i] = 1 - ((tempenergy - minval)/(maxval-minval));
		totalval += templateenergies_norm[i];
	}

	if(totalval <= 0){
		temporder=0;
		templateidx=int(Random()*STACKTEMPLATE[sttype][base1][base2].size());
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

	if(templateidx >= STACKTEMPLATE[sttype][base1][base2].size()){
		temporder=1;
		templateidx=templateidx-STACKTEMPLATE[sttype][base1][base2].size();
	}

}

bool MonteCarloSimulation::replacestackfrag(int &stack_base1, int &stack_base2, int &sttype, int &templateidx, int &temporder)
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
	
	/*
	static map<char, vector<int>> targetatoms={
		{'A',{12,13,14,15,16,17,18,19,20,21}},
		{'C',{12,13,14,15,16,17,18,19}},
		{'G',{12,13,14,15,16,17,18,19,20,21,22}},
		{'U',{12,13,14,15,16,17,18,19}}
	};
	*/

	char base1, base2;
	base1 = refinestr_temp[stack_base1].aaa;
	base2 = refinestr_temp[stack_base2].aaa;

	//int selected_fragidx = int(Random()*STACKTEMPLATE[sttype][base1][base2].size());
	int selected_fragidx = templateidx;
	int tempresidx1, tempresidx2;
	if(temporder == 0){
		tempresidx1 = 0;
		tempresidx2 = 1;
	}else if(temporder == 1){
		tempresidx1 = 1;
		tempresidx2 = 0;
	}

	vector<vector<double>> xyz_list1(targetatoms_sup[base1].size()+targetatoms_sup[base2].size(), vector<double>(3,0));
    vector<vector<double>> xyz_list2(targetatoms_sup[base1].size()+targetatoms_sup[base2].size(), vector<double>(3,0));

	if(base1 == base2 && sttype == 0){
		if(Random() > 0.5){
			int temp = stack_base1;
			stack_base1 = stack_base2;
			stack_base2 = temp;
		}
	}

	count = 0;
	for(j=0; j<targetatoms_sup[base1].size();j++)
	{
		atom=targetatoms_sup[base1][j];
		xyz_list1[count][0]=refinestr_temp[stack_base1].ptv[atom].x;
		xyz_list1[count][1]=refinestr_temp[stack_base1].ptv[atom].y;
		xyz_list1[count][2]=refinestr_temp[stack_base1].ptv[atom].z;

		xyz_list2[count][0]=STACKTEMPLATE[sttype][base1][base2][selected_fragidx][tempresidx1].ptv[atom].x;
		xyz_list2[count][1]=STACKTEMPLATE[sttype][base1][base2][selected_fragidx][tempresidx1].ptv[atom].y;
		xyz_list2[count][2]=STACKTEMPLATE[sttype][base1][base2][selected_fragidx][tempresidx1].ptv[atom].z;
		count+=1;
	}

	for(j=0; j<targetatoms_sup[base2].size();j++)
	{
		atom=targetatoms_sup[base2][j];
		xyz_list1[count][0]=refinestr_temp[stack_base2].ptv[atom].x;
		xyz_list1[count][1]=refinestr_temp[stack_base2].ptv[atom].y;
		xyz_list1[count][2]=refinestr_temp[stack_base2].ptv[atom].z;

		xyz_list2[count][0]=STACKTEMPLATE[sttype][base1][base2][selected_fragidx][tempresidx2].ptv[atom].x;
		xyz_list2[count][1]=STACKTEMPLATE[sttype][base1][base2][selected_fragidx][tempresidx2].ptv[atom].y;
		xyz_list2[count][2]=STACKTEMPLATE[sttype][base1][base2][selected_fragidx][tempresidx2].ptv[atom].z;
		count+=1;
	}


	bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

	for(j=0; j<targetatoms[base1].size();j++){
		atom=targetatoms[base1][j];
		atom_coord.clear();
		atom_coord.resize(3, 0.0);
		atom_coord[0]=STACKTEMPLATE[sttype][base1][base2][selected_fragidx][tempresidx1].ptv[atom].x;
		atom_coord[1]=STACKTEMPLATE[sttype][base1][base2][selected_fragidx][tempresidx1].ptv[atom].y;
		atom_coord[2]=STACKTEMPLATE[sttype][base1][base2][selected_fragidx][tempresidx1].ptv[atom].z;

		bf.Transform(atom_coord,RotMatix,TranVect,xyztmp);

		refinestr_temp[stack_base1].ptv[atom].x = xyztmp[0];
		refinestr_temp[stack_base1].ptv[atom].y = xyztmp[1];
		refinestr_temp[stack_base1].ptv[atom].z = xyztmp[2];
	}

	for(j=0; j<targetatoms[base2].size();j++){
        atom=targetatoms[base2][j];
		atom_coord.clear();
		atom_coord.resize(3, 0.0);
		atom_coord[0]=STACKTEMPLATE[sttype][base1][base2][selected_fragidx][tempresidx2].ptv[atom].x;
		atom_coord[1]=STACKTEMPLATE[sttype][base1][base2][selected_fragidx][tempresidx2].ptv[atom].y;
		atom_coord[2]=STACKTEMPLATE[sttype][base1][base2][selected_fragidx][tempresidx2].ptv[atom].z;

		bf.Transform(atom_coord,RotMatix,TranVect,xyztmp);

		refinestr_temp[stack_base2].ptv[atom].x = xyztmp[0];
		refinestr_temp[stack_base2].ptv[atom].y = xyztmp[1];
		refinestr_temp[stack_base2].ptv[atom].z = xyztmp[2];
	}

	return true;
}


int MonteCarloSimulation::selecttargetstack(int &basei, int &basej, double &bpscore)
{
	if (stackscoresum <= 0.0) return 0;
	int i,j;

	double totalval = 0.0;
	double randval = Random()*stackscoresum;

	for(int i=0;i<numbb;i++)
	{
		for(int j=i+1;j<numbb;j++)
		{
			if(stmap[i][j] > 0)
			{
				if(totalval <= randval && randval < totalval+stmap[i][j]){
					basei = i;
					basej = j;
					bpscore = stmap[i][j];
					return 1;
				}
				totalval += stmap[i][j];
			}
		}
	}
	return 0;
}

bool MonteCarloSimulation::fmovementstackfrag(double oldenergy, double &newenergy)
{
	int i, j, k;
	bool flagtor;
	int trandpos, trandposst, trandposed,span;
	double mean, std, upper, lower;
	int idx;
	double rmax, rnorm;
	int mtype1, mtype2;
	char base, base1, base2;
	bool flagdone;
	bool exitOuterLoop = false;

	if(stacklist.size() == 0){
		newenergy=oldenergy+1000000000;
		return false;
	}

	int stack_base1, stack_base2;
	double stackscore=-1;
	int stacktype;
	int stackselflag = selecttargetstack(stack_base1, stack_base2, stackscore);

	if(stackselflag == 0){
		newenergy=oldenergy+1000000000;
		return false;
	}
	
	if((stack_base1+1 == stack_base2) || (stack_base1 == stack_base2+1)){
		stacktype = 1;
	}else{
		stacktype = 0;
	}

	//cout<<"=============="<<endl;
	//cout<<target_stack_idx << "-" << stack_base1 << "-" << stack_base2 << "-" << stacktype<<endl;
	//cout<<stack_base1 << "-" << stack_base2<<endl;

	StrCopy(refinestr, refinestr_temp);
	resetfrag(0);

	base1 = refinestr_temp[stack_base1].aaa;
	base2 = refinestr_temp[stack_base2].aaa;

	vector<double> templateenergies;
	int templateidx;
	int temporder;
	double minrmsd = getstackenergy(stack_base1, stack_base2, stacktype, templateenergies, templateidx, temporder);

	cout<<"Base stacking-" << minrmsd << "-" << stackscore <<endl;
	
	if(minrmsd > 20.0){
		newenergy=oldenergy+1000000000;
		return false;
	}
	
	//if(minrmsd < 0.1){
	//selectstacktemplateidx(stack_base1, stack_base2, stacktype, templateenergies, templateidx, temporder);
	replacestackfrag(stack_base1, stack_base2, stacktype, templateidx, temporder);
	param2coordhydrogen(refinestr_temp,stack_base1,stack_base1+1);    
	param2coordhydrogen(refinestr_temp,stack_base2,stack_base2+1);
	//}

	map<int, vector<Point3s>> initstructure;
	initstructure[stack_base1].resize(AccAtomParamNum[base1][3]);
	initstructure[stack_base2].resize(AccAtomParamNum[base2][3]);

	for(j=0; j<AccAtomParamNum[base1][3];j++){
		initstructure[stack_base1][j].x=refinestr_temp[stack_base1].ptv[j].x;
		initstructure[stack_base1][j].y=refinestr_temp[stack_base1].ptv[j].y;
		initstructure[stack_base1][j].z=refinestr_temp[stack_base1].ptv[j].z;
	}

	for(j=0; j<AccAtomParamNum[base2][3];j++){
		initstructure[stack_base2][j].x=refinestr_temp[stack_base2].ptv[j].x;
		initstructure[stack_base2][j].y=refinestr_temp[stack_base2].ptv[j].y;
		initstructure[stack_base2][j].z=refinestr_temp[stack_base2].ptv[j].z;
	}

	Point3s coordA, coordB;
	coordA.x = refinestr_temp[stack_base1].ptv[8].x;
	coordA.y = refinestr_temp[stack_base1].ptv[8].y;
	coordA.z = refinestr_temp[stack_base1].ptv[8].z;
	coordB.x = refinestr_temp[stack_base2].ptv[8].x;
	coordB.y = refinestr_temp[stack_base2].ptv[8].y;
	coordB.z = refinestr_temp[stack_base2].ptv[8].z;

	vector<vector<int>> bindingmap(numbb, vector<int>(numbb, 0));
	if(stack_base1 - 1 >= 0) bindingmap[stack_base1 - 1][stack_base1] = 1;
	if(stack_base2 - 1 >= 0) bindingmap[stack_base2 - 1][stack_base2] = 1;
	if(stack_base1 + 1 < numbb) bindingmap[stack_base1][stack_base1 + 1] = 1;
	if(stack_base2 + 1 < numbb) bindingmap[stack_base2][stack_base2 + 1] = 1;
	bindingmap[stack_base1][stack_base1] = 2;
	bindingmap[stack_base2][stack_base2] = 2;

	vector<bool> moveresidues(numbb,false);
	moveresidues[stack_base1]=true;
	moveresidues[stack_base2]=true;

	int leng;
	vector<vector<int>> basepairlist;
	vector<vector<int>> noncbasepairlist;
	vector<double> basepairscorelist;
	vector<double> noncbasepairscorelist;

	vector<vector<int>> contactmap(numbb, vector<int>(numbb, 0));
	vector<vector<double>> stackingmap(numbb, vector<double>(numbb, 0));
	for(int i=0; i<numbb;i++){
		Point3s rp;
		if(i == stack_base1 || i == stack_base2) continue;
		rp.x = refinestr_temp[i].ptv[8].x - coordA.x;
		rp.y = refinestr_temp[i].ptv[8].y - coordA.y;
		rp.z = refinestr_temp[i].ptv[8].z - coordA.z;

		leng = rp.x*rp.x+rp.y*rp.y+rp.z*rp.z;

		//if(leng < GLOBALCUTOFF_SQR){
		contactmap[stack_base1][i] = 1;
		contactmap[i][stack_base1] = 1;
		//}

		if(stack_base1 == i){
			if(leng < CUTOFF_SQR+5.0){
				contactmap[stack_base1][i] = 2;
				contactmap[i][stack_base1] = 2;
			}
		}else{
			if(leng < CUTOFF_SQR+10.0){
				contactmap[stack_base1][i] = 2;
				contactmap[i][stack_base1] = 2;
			}
		}

		rp.x = refinestr_temp[i].ptv[8].x - coordB.x;
		rp.y = refinestr_temp[i].ptv[8].y - coordB.y;
		rp.z = refinestr_temp[i].ptv[8].z - coordB.z;

		leng = rp.x*rp.x+rp.y*rp.y+rp.z*rp.z;

		//if(leng < GLOBALCUTOFF_SQR){
		contactmap[stack_base2][i] = 1;
		contactmap[i][stack_base2] = 1;
		//}

		if(stack_base2 == i){
			if(leng < CUTOFF_SQR+5.0){
				contactmap[stack_base2][i] = 2;
				contactmap[i][stack_base2] = 2;
			}
		}else{
			if(leng < CUTOFF_SQR+10.0){
				contactmap[stack_base2][i] = 2;
				contactmap[i][stack_base2] = 2;
			}
		}

		if(cbpmap[stack_base1][i] > 0){ basepairlist.push_back({stack_base1, i}); basepairscorelist.push_back(cbpmap[stack_base1][i]);}
		if(cbpmap[stack_base2][i] > 0){ basepairlist.push_back({stack_base2, i}); basepairscorelist.push_back(cbpmap[stack_base2][i]);}
		if(noncbpmap[stack_base1][i] > 0){ noncbasepairlist.push_back({stack_base1, i}); noncbasepairscorelist.push_back(noncbpmap[stack_base1][i]);}
		if(noncbpmap[stack_base2][i] > 0){ noncbasepairlist.push_back({stack_base2, i}); noncbasepairscorelist.push_back(noncbpmap[stack_base2][i]);}

		stackingmap[stack_base1][i] = stmap[stack_base1][i];
		stackingmap[i][stack_base1] = stmap[i][stack_base1];
		stackingmap[stack_base2][i] = stmap[stack_base2][i];
		stackingmap[i][stack_base2] = stmap[i][stack_base2];
	}

	MonteCarloSimulation::MCfragmentoptimization(initstructure, bindingmap, moveresidues, contactmap, stackingmap, basepairlist, noncbasepairlist, basepairscorelist, noncbasepairscorelist, true);

	if(stacktype == 1){
		int istart=min(stack_base1-1, stack_base2-1);
		int iend=max(stack_base1+2, stack_base2+2);
		if(istart<0) istart=0;
		if(iend>=numbb) iend=numbb;

		coord2parambb(refinestr_temp,istart,iend);
		coord2paramsugar(refinestr_temp,istart,iend);
		coord2parambase(refinestr_temp,istart,iend);
		coord2paramhydrogen(refinestr_temp,istart,iend);
		coord2addparamhydrogen(refinestr_temp,istart,iend);
		coord2paramadd(refinestr_temp,istart,iend);

		setadjstackfragflap(istart,iend);
	}else{
		int istart1=stack_base1-1;
		int iend1=stack_base1+2;
		if(istart1<0) istart1=0;
		if(iend1>=numbb) iend1=numbb;

		int istart2=stack_base2-1;
		int iend2=stack_base2+2;
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

		setnonadjstackfragflap(istart1, iend1, istart2, iend2);
	}

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

