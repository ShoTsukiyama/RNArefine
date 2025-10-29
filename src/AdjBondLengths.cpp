


#include "AdjBondLengths.hpp"

using namespace std;

void get_template_idx(const ModelUnit &pdb_entry, const int c1, const int c2, const int r1, const int r2, map<string, map<string, vector<DNTTEMPLATEUnit>>> ideal_rna, int &templateidx){
	int i, j;
	double rmsd, minrmsd;
	int mintempidx;
	vector<double> tmp(3,0);
	vector<vector<double> > xyz_list1;
    vector<vector<double> > xyz_list2;
	vector<vector<double> > xyz_list3;
	vector<vector<double> > RotMatix, minRotMatix;  // U
	vector<double> TranVect, minTranVect;  // t
	int a,a1,a2;

	map<string, vector<string>> full_atoms={
        {"A", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'", " N9 ", " C8 ", " N7 ", " C5 ", " C6 ", " N6 ", " N1 ", " C2 ", " N3 ", " C4 "}},
        {"C", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'", " N1 ", " C2 ", " O2 ", " N3 ", " C4 ", " N4 ", " C5 ", " C6 "}},
        {"G", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'", " N9 ", " C8 ", " N7 ", " C5 ", " C6 ", " O6 ", " N1 ", " C2 ", " N2 ", " N3 ", " C4 "}},
        {"U", {" P  ", " OP1", " OP2", " O5'", " C5'", " C4'", " O4'", " C3'", " O3'", " C2'", " O2'", " C1'", " N1 ", " C2 ", " O2 ", " N3 ", " C4 ", " O4 ", " C5 ", " C6 "}},
    };

	string resname1(1, pdb_entry.chains[c1].residues[r1].resn[2]);
	string resname2(1, pdb_entry.chains[c2].residues[r2].resn[2]);
	string adj=resname1+"_"+resname2;

	int atomnum1 = full_atoms[resname1].size();
	int atomnum2 = full_atoms[resname2].size();

	//ideal_rna["loop"][adj].model

	minrmsd=10000000000; mintempidx=-1;
	for(i=0;i<ideal_rna["loop"][adj].size();i++){
		xyz_list1.clear();
        xyz_list2.clear();
		xyz_list3.clear();

		xyz_list1.assign(atomnum1+atomnum2,tmp);
        xyz_list2.assign(atomnum1+atomnum2,tmp);
		xyz_list3.assign(atomnum1+atomnum2,tmp);

		a=0;
        for (a1=0;a1<full_atoms[resname1].size();a1++)
        {
            xyz_list1[a][0]=ideal_rna["loop"][adj][i].model.chains[0].residues[0].atoms[a1].xyz[0];
            xyz_list1[a][1]=ideal_rna["loop"][adj][i].model.chains[0].residues[0].atoms[a1].xyz[1];
            xyz_list1[a][2]=ideal_rna["loop"][adj][i].model.chains[0].residues[0].atoms[a1].xyz[2];

            xyz_list2[a][0]=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz[0];
            xyz_list2[a][1]=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz[1];
            xyz_list2[a][2]=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz[2];
            a++;
        }

		for (a2=0;a2<full_atoms[resname2].size();a2++)
        {
			xyz_list1[a][0]=ideal_rna["loop"][adj][i].model.chains[0].residues[1].atoms[a2].xyz[0];
            xyz_list1[a][1]=ideal_rna["loop"][adj][i].model.chains[0].residues[1].atoms[a2].xyz[1];
            xyz_list1[a][2]=ideal_rna["loop"][adj][i].model.chains[0].residues[1].atoms[a2].xyz[2];
            xyz_list2[a][0]=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz[0];
            xyz_list2[a][1]=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz[1];
            xyz_list2[a][2]=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz[2];
            a++;
        }

		RotMatix.clear();
        TranVect.clear();
        RotateCoor(xyz_list1,xyz_list2, RotMatix, TranVect);

		j=0;
        rmsd=0;
        for (a=0;a<full_atoms[resname1].size();a++)
        {
            //if(pdb_entry.chains[c].residues[r].atoms[a].input != 1)  continue;
            ChangeCoor(ideal_rna["loop"][adj][i].model.chains[0].residues[0].atoms[a].xyz,RotMatix, TranVect, xyz_list3[j]);
            tmp[0]=(pdb_entry.chains[c1].residues[r1].atoms[a].xyz[0]-xyz_list3[j][0]);
            tmp[1]=(pdb_entry.chains[c1].residues[r1].atoms[a].xyz[1]-xyz_list3[j][1]);
            tmp[2]=(pdb_entry.chains[c1].residues[r1].atoms[a].xyz[2]-xyz_list3[j][2]);
            j++;
            rmsd+=sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
        }

		for (a=0;a<full_atoms[resname2].size();a++)
        {
            //if(pdb_entry.chains[c].residues[r].atoms[a].input != 1)  continue;
            ChangeCoor(ideal_rna["loop"][adj][i].model.chains[0].residues[1].atoms[a].xyz,RotMatix, TranVect, xyz_list3[j]);
            tmp[0]=(pdb_entry.chains[c2].residues[r2].atoms[a].xyz[0]-xyz_list3[j][0]);
            tmp[1]=(pdb_entry.chains[c2].residues[r2].atoms[a].xyz[1]-xyz_list3[j][1]);
            tmp[2]=(pdb_entry.chains[c2].residues[r2].atoms[a].xyz[2]-xyz_list3[j][2]);
            j++;
            rmsd+=sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
        }

		if(rmsd < minrmsd){
            minrmsd=rmsd;
            mintempidx=i;
            minRotMatix=RotMatix;
            minTranVect=TranVect;
        }
	}

	templateidx=mintempidx;
}

int fix_adj_bond_lengths(ModelUnit &pep, map<size_t, vector<int>> loop_vec, map<string, map<string, vector<DNTTEMPLATEUnit>>> ideal_rna, double tolerance){
	
	// Initialize variables
	double b_ideal = 0;
	double b_std = 0;
	double b = 0;
	double deltaB = 0;
	double diff = 0;
	int a2 = 0;

	double a1x = 0;
	double a1y = 0;
	double a1z = 0;
	double a2x = 0;
	double a2y = 0;
	double a2z = 0;

	int moved = 0;
	int templateidx;

	// Iterate through chains
	for (int c=0; c<pep.chains.size(); c++){
		// Iterate through residues

		for(int r_idx=0; r_idx<loop_vec[c].size();r_idx++){
			if(loop_vec[c][r_idx] == 0) continue;

			int r2=loop_vec[c][r_idx];
			int r1=r2-1;

			if (r1 >= 0 && isConnected(pep,c,r1)) {
				get_template_idx(pep, c, c, r1, r2, ideal_rna, templateidx);

				for (int a1=0; a1<pep.chains[c].residues[r1].atoms.size(); a1++){
					string atom1 = pep.chains[c].residues[r1].atoms[a1].name;
					for (a2=0; a2<pep.chains[c].residues[r2].atoms.size(); a2++){
						string atom2 = pep.chains[c].residues[r2].atoms[a2].name;
						bool pass = false;
						// if neither atom is movable, immediately continue
						if (pep.chains[c].residues[r1].atoms[a1].movable==0 && pep.chains[c].residues[r2].atoms[a2].movable==0){
							//cout << "none movable " << resname << resnumber << atom1 << atom2 << endl;
							continue;
						}

						string resname1(1, pep.chains[c].residues[r1].resn[2]);
						string resname2(1, pep.chains[c].residues[r2].resn[2]);
						string adj=resname1+"_"+resname2;

						b_ideal = ideal_rna["loop"][adj][templateidx].pwdist_mean[a1][a2];
						b_std = ideal_rna["loop"][adj][templateidx].pwdist_std[a1][a2];

						a1x = pep.chains[c].residues[r1].atoms[a1].xyz[0];
						a1y = pep.chains[c].residues[r1].atoms[a1].xyz[1];
						a1z = pep.chains[c].residues[r1].atoms[a1].xyz[2];
						a2x = pep.chains[c].residues[r2].atoms[a2].xyz[0];
						a2y = pep.chains[c].residues[r2].atoms[a2].xyz[1];
						a2z = pep.chains[c].residues[r2].atoms[a2].xyz[2];
						// Calculate the current bond length
						b = sqrt((a2x-a1x)*(a2x-a1x) + (a2y-a1y)*(a2y-a1y) + (a2z-a1z)*(a2z-a1z));
						// Calculate fractional difference between the current and ideal bond lengths
						diff = fabs(b - b_ideal) / b_ideal;

						// Check if difference is greater than or equal to the allowed tolerance
						if (diff >= tolerance && b>0.0001){
							//cout<<b << "-" << b_ideal << "-" << diff << "-" << tolerance<<endl;
							deltaB = b - b_ideal;

							if (pep.chains[c].residues[r1].atoms[a1].movable==1 && pep.chains[c].residues[r2].atoms[a2].movable==0){
								//cout << "a1 movable " << resname << resnumber << atom1 << atom2 << "b: " << b << "b_ideal: " << b_ideal << "diff: " << diff << endl;
								moved++;
								pep.chains[c].residues[r1].atoms[a1].xyz[0] = 0.01*(deltaB/b)*(a2x-a1x) + a1x;
								pep.chains[c].residues[r1].atoms[a1].xyz[1] = 0.01*(deltaB/b)*(a2y-a1y) + a1y;
								pep.chains[c].residues[r1].atoms[a1].xyz[2] = 0.01*(deltaB/b)*(a2z-a1z) + a1z;
							} else if (pep.chains[c].residues[r1].atoms[a1].movable==0 && pep.chains[c].residues[r2].atoms[a2].movable==1){
								//cout << "a2 movable " << resname << resnumber << atom1 << atom2 << "b: " << b << "b_ideal: " << b_ideal << "diff: " << diff << endl;
								moved++;
								pep.chains[c].residues[r2].atoms[a2].xyz[0] = 0.01*(deltaB/b)*(a1x-a2x) + a2x;
								pep.chains[c].residues[r2].atoms[a2].xyz[1] = 0.01*(deltaB/b)*(a1y-a2y) + a2y;
								pep.chains[c].residues[r2].atoms[a2].xyz[2] = 0.01*(deltaB/b)*(a1z-a2z) + a2z;
							} else if (pep.chains[c].residues[r1].atoms[a1].movable==1 && pep.chains[c].residues[r2].atoms[a2].movable==1){
								//cout << "both movable " << resname << resnumber << atom1 << atom2 << "b: " << b << "b_ideal: " << b_ideal << "diff: " << diff << endl;
								moved++;
								// move each atom half the required distance
								pep.chains[c].residues[r2].atoms[a2].xyz[0] = 0.01*0.5*((deltaB/b)*(a1x-a2x)) + a2x;
								pep.chains[c].residues[r2].atoms[a2].xyz[1] = 0.01*0.5*((deltaB/b)*(a1y-a2y)) + a2y;
								pep.chains[c].residues[r2].atoms[a2].xyz[2] = 0.01*0.5*((deltaB/b)*(a1z-a2z)) + a2z;
								pep.chains[c].residues[r1].atoms[a1].xyz[0] = 0.01*0.5*((deltaB/b)*(a2x-a1x)) + a1x;
								pep.chains[c].residues[r1].atoms[a1].xyz[1] = 0.01*0.5*((deltaB/b)*(a2y-a1y)) + a1y;
								pep.chains[c].residues[r1].atoms[a1].xyz[2] = 0.01*0.5*((deltaB/b)*(a2z-a1z)) + a1z;
							}
						}	
					}
				}
			}
		}
	}

	return moved;
}



