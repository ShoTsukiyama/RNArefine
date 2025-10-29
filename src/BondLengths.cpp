
#include "BondLengths.hpp"

/* Check if residue r and residue r+1 are connected */
bool isConnected(const ModelUnit &pep, const int c, const int r)
{
    for (int a=0;a<=12;a++) if (Points2Distance(
        pep.chains[c].residues[r].atoms[a].xyz,
        pep.chains[c].residues[r+1].atoms[a].xyz)<=8) return true;
    return false;
}

int fix_bond_lengths (ModelUnit &pep, double tolerance){
	
	// Initialize variables
	double b_ideal = 0;
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

	// Iterate through chains
	for (int c=0; c<pep.chains.size(); c++){
		// Iterate through residues
		for (int r=0; r<pep.chains[c].residues.size(); r++){
			string resname = pep.chains[c].residues[r].resn;
			int resnumber = pep.chains[c].residues[r].resi;
			
			/* Special case to check the O3'-P bonds */
			if (r != (pep.chains[c].residues.size()-1) && isConnected(pep,c,r)) {
                for (int a1=7;a1<=8;a1++) {
			        // Get (x,y,z) coordinates
			        a1x = pep.chains[c].residues[r].atoms[a1].xyz[0];
			        a1y = pep.chains[c].residues[r].atoms[a1].xyz[1];
			        a1z = pep.chains[c].residues[r].atoms[a1].xyz[2];
                    for (a2=0;a2<=3;a2++) {
                        b_ideal = OSP_bonds[a1][a2];
                        if (b_ideal<=0) continue;
					    a2x = pep.chains[c].residues[r+1].atoms[a2].xyz[0];
					    a2y = pep.chains[c].residues[r+1].atoms[a2].xyz[1];
					    a2z = pep.chains[c].residues[r+1].atoms[a2].xyz[2];
					    // Calculate the current bond length
					    b = sqrt((a2x-a1x)*(a2x-a1x) + (a2y-a1y)*(a2y-a1y) + (a2z-a1z)*(a2z-a1z));
					    // Calculate fractional difference between the current and ideal bond lengths
					    diff = fabs(b - b_ideal) / b_ideal;
					    // Check if difference is greater than or equal to the allowed tolerance
					    if (diff >= tolerance && b>0.0001){
						    //cout << resname << resnumber << atom1 << atom2 << b << ' ' << b_ideal << ' ' << diff <<endl;
						    deltaB = b - b_ideal; // difference in bond length from ideal
						    // Check which atoms are movable and assign new coordinates to the movable atom
						    if (pep.chains[c].residues[r].atoms[a1].movable==1 && pep.chains[c].residues[r].atoms[a2].movable==0){
							    //cout << "a1 movable " << resname << resnumber << atom1 << atom2 << "b: " << b << "b_ideal: " << b_ideal << "diff: " << diff << endl;
							    moved++;
							    pep.chains[c].residues[r].atoms[a1].xyz[0] = (deltaB/b)*(a2x-a1x) + a1x;
							    pep.chains[c].residues[r].atoms[a1].xyz[1] = (deltaB/b)*(a2y-a1y) + a1y;
							    pep.chains[c].residues[r].atoms[a1].xyz[2] = (deltaB/b)*(a2z-a1z) + a1z;
						    } else if (pep.chains[c].residues[r].atoms[a1].movable==0 && pep.chains[c].residues[r].atoms[a2].movable==1){
							    //cout << "a2 movable " << resname << resnumber << atom1 << atom2 << "b: " << b << "b_ideal: " << b_ideal << "diff: " << diff << endl;
							    moved++;
							    pep.chains[c].residues[r+1].atoms[a2].xyz[0] = (deltaB/b)*(a1x-a2x) + a2x;
							    pep.chains[c].residues[r+1].atoms[a2].xyz[1] = (deltaB/b)*(a1y-a2y) + a2y;
							    pep.chains[c].residues[r+1].atoms[a2].xyz[2] = (deltaB/b)*(a1z-a2z) + a2z;
						    } else if (pep.chains[c].residues[r].atoms[a1].movable==1 && pep.chains[c].residues[r].atoms[a2].movable==1){
							    //cout << "both movable " << resname << resnumber << atom1 << atom2 << "b: " << b << "b_ideal: " << b_ideal << "diff: " << diff << endl;
							    moved++;
							    // move each atom half the required distance
							    pep.chains[c].residues[r+1].atoms[a2].xyz[0] = 0.5*((deltaB/b)*(a1x-a2x)) + a2x;
							    pep.chains[c].residues[r+1].atoms[a2].xyz[1] = 0.5*((deltaB/b)*(a1y-a2y)) + a2y;
							    pep.chains[c].residues[r+1].atoms[a2].xyz[2] = 0.5*((deltaB/b)*(a1z-a2z)) + a2z;
							    pep.chains[c].residues[r].atoms[a1].xyz[0] = 0.5*((deltaB/b)*(a2x-a1x)) + a1x;
							    pep.chains[c].residues[r].atoms[a1].xyz[1] = 0.5*((deltaB/b)*(a2y-a1y)) + a1y;
							    pep.chains[c].residues[r].atoms[a1].xyz[2] = 0.5*((deltaB/b)*(a2z-a1z)) + a1z;
						    }
					    }	
					}
				}
			}

			// Iterate through all pairs of backbone atoms
			for (int a1=0; a1<pep.chains[c].residues[r].atoms.size(); a1++){
				string atom1 = pep.chains[c].residues[r].atoms[a1].name;
				for (a2=a1+1; a2<pep.chains[c].residues[r].atoms.size(); a2++){
					string atom2 = pep.chains[c].residues[r].atoms[a2].name;
					bool pass = false;
					// if neither atom is movable, immediately continue
					if (pep.chains[c].residues[r].atoms[a1].movable==0 && pep.chains[c].residues[r].atoms[a2].movable==0){
						//cout << "none movable " << resname << resnumber << atom1 << atom2 << endl;
						continue;
					}
					//cout << atom1 << pep.chains[c].residues[r].atoms[a1].movable << endl;
					//cout << "find atom in array" << endl;
					
					// Get ideal bond lengths
					if (resname == "  A"){
						if (A_bonds[a1][a2]>0){
							pass = true;
							b_ideal = A_bonds[a1][a2];
						}
					} else if (resname == "  C"){
						if (C_bonds[a1][a2]>0){
							pass = true;
							b_ideal = C_bonds[a1][a2];
						}
					} else if (resname == "  G"){
						if (G_bonds[a1][a2]>0){
							pass = true;
							b_ideal = G_bonds[a1][a2];
						}
					} else if (resname == "  U"){
						if (U_bonds[a1][a2]>0){
							pass = true;
							b_ideal = U_bonds[a1][a2];
						}
					}

					if (! pass) continue;
					// Get (x,y,z) coordinates
					a1x = pep.chains[c].residues[r].atoms[a1].xyz[0];
					a1y = pep.chains[c].residues[r].atoms[a1].xyz[1];
					a1z = pep.chains[c].residues[r].atoms[a1].xyz[2];
					a2x = pep.chains[c].residues[r].atoms[a2].xyz[0];
					a2y = pep.chains[c].residues[r].atoms[a2].xyz[1];
					a2z = pep.chains[c].residues[r].atoms[a2].xyz[2];
					// Calculate the current bond length
					b = sqrt((a2x-a1x)*(a2x-a1x) + (a2y-a1y)*(a2y-a1y) + (a2z-a1z)*(a2z-a1z));
					// Calculate fractional difference between the current and ideal bond lengths
					diff = fabs(b - b_ideal) / b_ideal;
					// Check if difference is greater than or equal to the allowed tolerance
					if (diff >= tolerance && b>0.0001){
						//cout << resname << resnumber << atom1 << atom2 << b << ' ' << b_ideal << ' ' << diff <<endl;
						deltaB = b - b_ideal; // difference in bond length from ideal
						// Check which atoms are movable and assign new coordinates to the movable atom
						if (pep.chains[c].residues[r].atoms[a1].movable==1 && pep.chains[c].residues[r].atoms[a2].movable==0){
							//cout << "a1 movable " << resname << resnumber << atom1 << atom2 << "b: " << b << "b_ideal: " << b_ideal << "diff: " << diff << endl;
							moved++;
							pep.chains[c].residues[r].atoms[a1].xyz[0] = (deltaB/b)*(a2x-a1x) + a1x;
							pep.chains[c].residues[r].atoms[a1].xyz[1] = (deltaB/b)*(a2y-a1y) + a1y;
							pep.chains[c].residues[r].atoms[a1].xyz[2] = (deltaB/b)*(a2z-a1z) + a1z;
						} else if (pep.chains[c].residues[r].atoms[a1].movable==0 && pep.chains[c].residues[r].atoms[a2].movable==1){
							//cout << "a2 movable " << resname << resnumber << atom1 << atom2 << "b: " << b << "b_ideal: " << b_ideal << "diff: " << diff << endl;
							moved++;
							pep.chains[c].residues[r].atoms[a2].xyz[0] = (deltaB/b)*(a1x-a2x) + a2x;
							pep.chains[c].residues[r].atoms[a2].xyz[1] = (deltaB/b)*(a1y-a2y) + a2y;
							pep.chains[c].residues[r].atoms[a2].xyz[2] = (deltaB/b)*(a1z-a2z) + a2z;
						} else if (pep.chains[c].residues[r].atoms[a1].movable==1 && pep.chains[c].residues[r].atoms[a2].movable==1){
							//cout << "both movable " << resname << resnumber << atom1 << atom2 << "b: " << b << "b_ideal: " << b_ideal << "diff: " << diff << endl;
							moved++;
							// move each atom half the required distance
							pep.chains[c].residues[r].atoms[a2].xyz[0] = 0.5*((deltaB/b)*(a1x-a2x)) + a2x;
							pep.chains[c].residues[r].atoms[a2].xyz[1] = 0.5*((deltaB/b)*(a1y-a2y)) + a2y;
							pep.chains[c].residues[r].atoms[a2].xyz[2] = 0.5*((deltaB/b)*(a1z-a2z)) + a2z;
							pep.chains[c].residues[r].atoms[a1].xyz[0] = 0.5*((deltaB/b)*(a2x-a1x)) + a1x;
							pep.chains[c].residues[r].atoms[a1].xyz[1] = 0.5*((deltaB/b)*(a2y-a1y)) + a1y;
							pep.chains[c].residues[r].atoms[a1].xyz[2] = 0.5*((deltaB/b)*(a2z-a1z)) + a1z;
						}
					}	
				}
			}
		}
	}

	return moved;
}




