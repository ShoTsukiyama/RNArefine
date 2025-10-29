
#include "RefineRNA.hpp"

void getContactResidueMap(ModelUnit &pdb_entry, int c, vector<vector<bool>> &rescontacts){
    double a1x, a1y, a1z, a2x, a2y, a2z, b;
    double contactthresh=18.161;
    int L=pdb_entry.chains[c].residues.size();
    rescontacts.resize(L, vector<bool>(L, false));

    for(int i=0;i<L;i++){
        a1x = pdb_entry.chains[c].residues[i].atoms[11].xyz[0];
		a1y = pdb_entry.chains[c].residues[i].atoms[11].xyz[1];
		a1z = pdb_entry.chains[c].residues[i].atoms[11].xyz[2];
        for(int j=0;j<L;j++){
            a2x = pdb_entry.chains[c].residues[j].atoms[11].xyz[0];
			a2y = pdb_entry.chains[c].residues[j].atoms[11].xyz[1];
			a2z = pdb_entry.chains[c].residues[j].atoms[11].xyz[2];

            b = sqrt((a2x-a1x)*(a2x-a1x) + (a2y-a1y)*(a2y-a1y) + (a2z-a1z)*(a2z-a1z));

            if(b < contactthresh){
                rescontacts[i][j] = true;
            }
        }
    }
}

int correct_geometry(ModelUnit &pdb_entry, const vector<vector<bool>> &rescontacts, double alpha, double tolerance){
    double a1x, a1y, a1z, a2x, a2y, a2z, b;
    int c, r, a, a2;
    double b_ideal, diff, deltaB;
    double dx, dy, dz;
    int moved = 0;
    //double alpha=0.2;

    for (c=0; c<pdb_entry.chains.size(); c++){
        for (r=0; r<pdb_entry.chains[c].residues.size(); r++){
            for(a=0;a<pdb_entry.chains[c].residues[r].atoms.size();a++){
                dx=0; dy=0;dz=0;

                string atom = pdb_entry.chains[c].residues[r].atoms[a].name;
                string resname = pdb_entry.chains[c].residues[r].resn;
			    int resnumber = pdb_entry.chains[c].residues[r].resi;

                a1x = pdb_entry.chains[c].residues[r].atoms[a].xyz[0];
			    a1y = pdb_entry.chains[c].residues[r].atoms[a].xyz[1];
			    a1z = pdb_entry.chains[c].residues[r].atoms[a].xyz[2];

                if (r != (pdb_entry.chains[c].residues.size()-1) && isConnected(pdb_entry,c,r))
                {
                    if(a == 7 || a == 8){
                        for (a2=0;a2<=3;a2++) {
                            b_ideal = OSP_bonds[a][a2];
                            if (b_ideal<=0) continue;
                            a2x = pdb_entry.chains[c].residues[r+1].atoms[a2].xyz[0];
					        a2y = pdb_entry.chains[c].residues[r+1].atoms[a2].xyz[1];
					        a2z = pdb_entry.chains[c].residues[r+1].atoms[a2].xyz[2];

                            b = sqrt((a2x-a1x)*(a2x-a1x) + (a2y-a1y)*(a2y-a1y) + (a2z-a1z)*(a2z-a1z));
                            deltaB = b - b_ideal;
                            diff = fabs(deltaB) / b_ideal;

                            if (diff >= tolerance && b>0.0001){
                                moved++;
                                dx += alpha*((deltaB/b)*(a2x-a1x));
                                dy += alpha*((deltaB/b)*(a2y-a1y));
                                dz += alpha*((deltaB/b)*(a2z-a1z));
                            }
                        }
                    }
                }

                if (r-1 > 0 && isConnected(pdb_entry,c,r-1))
                {
                    if(a == 0 || a == 1 || a == 2 || a == 3)
                    {
                        for (int a2=7;a2<=8;a2++) {
                            b_ideal = OSP_bonds[a2][a];
                            if (b_ideal<=0) continue;
                            a2x = pdb_entry.chains[c].residues[r-1].atoms[a2].xyz[0];
                            a2y = pdb_entry.chains[c].residues[r-1].atoms[a2].xyz[1];
                            a2z = pdb_entry.chains[c].residues[r-1].atoms[a2].xyz[2];

                            b = sqrt((a2x-a1x)*(a2x-a1x) + (a2y-a1y)*(a2y-a1y) + (a2z-a1z)*(a2z-a1z));
                            deltaB = b - b_ideal;
                            diff = fabs(deltaB) / b_ideal;

                            if (diff >= tolerance && b>0.0001){
                                moved++;
                                dx += alpha*((deltaB/b)*(a2x-a1x));
                                dy += alpha*((deltaB/b)*(a2y-a1y));
                                dz += alpha*((deltaB/b)*(a2z-a1z));
                            }
                        }
                    }
                }

                for (int a2=0; a2<pdb_entry.chains[c].residues[r].atoms.size(); a2++){
                    string atom2 = pdb_entry.chains[c].residues[r].atoms[a2].name;

                    bool pass = false;
                    if (resname == "  A"){
						if (A_bonds[a][a2]>0){
							pass = true;
							b_ideal = A_bonds[a][a2];
						}
					} else if (resname == "  C"){
						if (C_bonds[a][a2]>0){
							pass = true;
							b_ideal = C_bonds[a][a2];
						}
					} else if (resname == "  G"){
						if (G_bonds[a][a2]>0){
							pass = true;
							b_ideal = G_bonds[a][a2];
						}
					} else if (resname == "  U"){
						if (U_bonds[a][a2]>0){
							pass = true;
							b_ideal = U_bonds[a][a2];
						}
					}

                   
                    if (! pass) continue;
                    a2x = pdb_entry.chains[c].residues[r].atoms[a2].xyz[0];
					a2y = pdb_entry.chains[c].residues[r].atoms[a2].xyz[1];
					a2z = pdb_entry.chains[c].residues[r].atoms[a2].xyz[2];

                    b = sqrt((a2x-a1x)*(a2x-a1x) + (a2y-a1y)*(a2y-a1y) + (a2z-a1z)*(a2z-a1z));
                    deltaB = b - b_ideal;
                    diff = fabs(deltaB) / b_ideal;

                    if (diff >= tolerance && b>0.0001){
                        moved++;
                        
                        dx += alpha*((deltaB/b)*(a2x-a1x));
						dy += alpha*((deltaB/b)*(a2y-a1y));
						dz += alpha*((deltaB/b)*(a2z-a1z));
                    }
                }

                pdb_entry.chains[c].residues[r].atoms[a].xyz[0] +=dx;
				pdb_entry.chains[c].residues[r].atoms[a].xyz[1] +=dy;
				pdb_entry.chains[c].residues[r].atoms[a].xyz[2] +=dz;
            }
        }
    }
    return moved;
}


int correct_clash(ModelUnit &pdb_entry, const vector<vector<bool>> &rescontacts, double alpha){
    double a1x, a1y, a1z, a2x, a2y, a2z, b;
    int c, r, a, r2, a2;
    double b_ideal, diff, deltaB;
    double dx, dy, dz;
    double vdw1 = 0;
	double vdw2 = 0;
    int moved = 0;
    //double alpha=0.25;

    for (c=0; c<pdb_entry.chains.size(); c++){
        for (r=0; r<pdb_entry.chains[c].residues.size(); r++){
            for(a=0;a<pdb_entry.chains[c].residues[r].atoms.size();a++){
                dx=0; dy=0;dz=0;

                string atom = pdb_entry.chains[c].residues[r].atoms[a].name;
                string resname = pdb_entry.chains[c].residues[r].resn;
			    int resnumber = pdb_entry.chains[c].residues[r].resi;
                //string atom1 = pdb_entry.chains[c].residues[r].atoms[a].name;

                if (resname == "  A") vdw1 = A_vdw[a];
                else if (resname == "  C") vdw1 = C_vdw[a];
                else if (resname == "  G") vdw1 = G_vdw[a];
                else if (resname == "  U") vdw1 = U_vdw[a];

                a1x = pdb_entry.chains[c].residues[r].atoms[a].xyz[0];
			    a1y = pdb_entry.chains[c].residues[r].atoms[a].xyz[1];
			    a1z = pdb_entry.chains[c].residues[r].atoms[a].xyz[2];

                //for (c2=0; c2<pdb_entry.chains.size(); c2++){
                for (int r2=0; r2<pdb_entry.chains[c].residues.size(); r2++){
                    if(rescontacts[r][r2]){
                        string resname2 = pdb_entry.chains[c].residues[r2].resn;
					    int resnumber2 = pdb_entry.chains[c].residues[r2].resi;

                        for (int a2=0; a2<pdb_entry.chains[c].residues[r2].atoms.size(); a2++){
                            //string atom2 = pdb_entry.chains[c].residues[r2].atoms[a2].name;

                            if(r+1==r2 && a < 9 && a2 < 4 && OSP_bonds[a][a2] > 0) continue;
                            if(r2+1==r && a2 < 9 && a < 4 && OSP_bonds[a2][a] > 0) continue;
                            if(r == r2) continue;
                            if(r == r2 && a == a2) continue;
                            
                            if (resname2 == "  A") vdw2 = A_vdw[a2];
                            else if (resname2 == "  C") vdw2 = C_vdw[a2];
                            else if (resname2 == "  G") vdw2 = G_vdw[a2];
                            else if (resname2 == "  U") vdw2 = U_vdw[a2];

                            double vdw = vdw1 + vdw2;
                            double vdw_cutoff =  vdw - 0.4;

                            double a2x = pdb_entry.chains[c].residues[r2].atoms[a2].xyz[0];
							double a2y = pdb_entry.chains[c].residues[r2].atoms[a2].xyz[1];
							double a2z = pdb_entry.chains[c].residues[r2].atoms[a2].xyz[2];

                            double d = sqrt((a2x-a1x)*(a2x-a1x) + (a2y-a1y)*(a2y-a1y) + (a2z-a1z)*(a2z-a1z));
					
                            if (d < vdw_cutoff && d > 0.0001){
                                double delta = d - vdw_cutoff; 

                                dx += alpha*((delta/d)*(a2x-a1x));
								dy += alpha*((delta/d)*(a2y-a1y));
							    dz += alpha*((delta/d)*(a2z-a1z));

                                moved++;
                            }
                        }
                    }
                }
                //}
                pdb_entry.chains[c].residues[r].atoms[a].xyz[0] +=dx;
				pdb_entry.chains[c].residues[r].atoms[a].xyz[1] +=dy;
				pdb_entry.chains[c].residues[r].atoms[a].xyz[2] +=dz;
            }
        }
    }
    return moved;
}

int correct_combine(ModelUnit &pdb_entry, const vector<vector<bool>> &rescontacts, double alpha, double beta, double tolerance){
    double a1x, a1y, a1z, a2x, a2y, a2z, b;
    int c, r, a, a2;
    double b_ideal, diff, deltaB;
    double dx, dy, dz;
    double vdw1 = 0;
	double vdw2 = 0;
    int moved = 0;
    //double alpha=0.1;
    //double beta=0.1;

    //cout<<alpha << "-" << beta<<endl;

    for (c=0; c<pdb_entry.chains.size(); c++){
        for (r=0; r<pdb_entry.chains[c].residues.size(); r++){
            for(a=0;a<pdb_entry.chains[c].residues[r].atoms.size();a++){
                dx=0; dy=0;dz=0;

                string atom = pdb_entry.chains[c].residues[r].atoms[a].name;
                string resname = pdb_entry.chains[c].residues[r].resn;
			    int resnumber = pdb_entry.chains[c].residues[r].resi;

                a1x = pdb_entry.chains[c].residues[r].atoms[a].xyz[0];
			    a1y = pdb_entry.chains[c].residues[r].atoms[a].xyz[1];
			    a1z = pdb_entry.chains[c].residues[r].atoms[a].xyz[2];

                if (r != (pdb_entry.chains[c].residues.size()-1) && isConnected(pdb_entry,c,r))
                {
                    if(a == 7 || a == 8){
                        for (a2=0;a2<=3;a2++) {
                            b_ideal = OSP_bonds[a][a2];
                            if (b_ideal<=0) continue;
                            a2x = pdb_entry.chains[c].residues[r+1].atoms[a2].xyz[0];
					        a2y = pdb_entry.chains[c].residues[r+1].atoms[a2].xyz[1];
					        a2z = pdb_entry.chains[c].residues[r+1].atoms[a2].xyz[2];

                            b = sqrt((a2x-a1x)*(a2x-a1x) + (a2y-a1y)*(a2y-a1y) + (a2z-a1z)*(a2z-a1z));
                            deltaB = b - b_ideal;
                            diff = fabs(deltaB) / b_ideal;

                            if (diff >= tolerance && b>0.0001){
                                moved++;
                                dx += alpha*((deltaB/b)*(a2x-a1x));
                                dy += alpha*((deltaB/b)*(a2y-a1y));
                                dz += alpha*((deltaB/b)*(a2z-a1z));
                            }
                        }
                    }
                }

                if (r-1 > 0 && isConnected(pdb_entry,c,r-1))
                {
                    if(a == 0 || a == 1 || a == 2 || a == 3)
                    {
                        for (int a2=7;a2<=8;a2++) {
                            b_ideal = OSP_bonds[a2][a];
                            if (b_ideal<=0) continue;
                            a2x = pdb_entry.chains[c].residues[r-1].atoms[a2].xyz[0];
                            a2y = pdb_entry.chains[c].residues[r-1].atoms[a2].xyz[1];
                            a2z = pdb_entry.chains[c].residues[r-1].atoms[a2].xyz[2];

                            b = sqrt((a2x-a1x)*(a2x-a1x) + (a2y-a1y)*(a2y-a1y) + (a2z-a1z)*(a2z-a1z));
                            deltaB = b - b_ideal;
                            diff = fabs(deltaB) / b_ideal;

                            if (diff >= tolerance && b>0.0001){
                                moved++;
                                dx += alpha*((deltaB/b)*(a2x-a1x));
                                dy += alpha*((deltaB/b)*(a2y-a1y));
                                dz += alpha*((deltaB/b)*(a2z-a1z));
                            }
                        }
                    }
                }

                for (int a2=0; a2<pdb_entry.chains[c].residues[r].atoms.size(); a2++){
                    string atom2 = pdb_entry.chains[c].residues[r].atoms[a2].name;

                    bool pass = false;
                    if (resname == "  A"){
						if (A_bonds[a][a2]>0){
							pass = true;
							b_ideal = A_bonds[a][a2];
						}
					} else if (resname == "  C"){
						if (C_bonds[a][a2]>0){
							pass = true;
							b_ideal = C_bonds[a][a2];
						}
					} else if (resname == "  G"){
						if (G_bonds[a][a2]>0){
							pass = true;
							b_ideal = G_bonds[a][a2];
						}
					} else if (resname == "  U"){
						if (U_bonds[a][a2]>0){
							pass = true;
							b_ideal = U_bonds[a][a2];
						}
					}

                   
                    if (! pass) continue;
                    a2x = pdb_entry.chains[c].residues[r].atoms[a2].xyz[0];
					a2y = pdb_entry.chains[c].residues[r].atoms[a2].xyz[1];
					a2z = pdb_entry.chains[c].residues[r].atoms[a2].xyz[2];

                    b = sqrt((a2x-a1x)*(a2x-a1x) + (a2y-a1y)*(a2y-a1y) + (a2z-a1z)*(a2z-a1z));
                    deltaB = b - b_ideal;
                    diff = fabs(deltaB) / b_ideal;

                    if (diff >= tolerance && b>0.0001){
                        moved++;
                        
                        dx += alpha*((deltaB/b)*(a2x-a1x));
						dy += alpha*((deltaB/b)*(a2y-a1y));
						dz += alpha*((deltaB/b)*(a2z-a1z));
                    }
                }

                if (resname == "  A") vdw1 = A_vdw[a];
                else if (resname == "  C") vdw1 = C_vdw[a];
                else if (resname == "  G") vdw1 = G_vdw[a];
                else if (resname == "  U") vdw1 = U_vdw[a];

                a1x = pdb_entry.chains[c].residues[r].atoms[a].xyz[0];
			    a1y = pdb_entry.chains[c].residues[r].atoms[a].xyz[1];
			    a1z = pdb_entry.chains[c].residues[r].atoms[a].xyz[2];

                //for (c2=0; c2<pdb_entry.chains.size(); c2++){
                for (int r2=0; r2<pdb_entry.chains[c].residues.size(); r2++){
                    if(rescontacts[r][r2]){
                        string resname2 = pdb_entry.chains[c].residues[r2].resn;
					    int resnumber2 = pdb_entry.chains[c].residues[r2].resi;

                        for (int a2=0; a2<pdb_entry.chains[c].residues[r2].atoms.size(); a2++){
                            //string atom2 = pdb_entry.chains[c].residues[r2].atoms[a2].name;

                            if(r+1==r2 && a < 9 && a2 < 4 && OSP_bonds[a][a2] > 0) continue;
                            if(r2+1==r && a2 < 9 && a < 4 && OSP_bonds[a2][a] > 0) continue;
                            if(r == r2) continue;
                            if(r == r2 && a == a2) continue;
                            
                            if (resname2 == "  A") vdw2 = A_vdw[a2];
                            else if (resname2 == "  C") vdw2 = C_vdw[a2];
                            else if (resname2 == "  G") vdw2 = G_vdw[a2];
                            else if (resname2 == "  U") vdw2 = U_vdw[a2];

                            double vdw = vdw1 + vdw2;
                            double vdw_cutoff =  vdw - 0.4;

                            double a2x = pdb_entry.chains[c].residues[r2].atoms[a2].xyz[0];
							double a2y = pdb_entry.chains[c].residues[r2].atoms[a2].xyz[1];
							double a2z = pdb_entry.chains[c].residues[r2].atoms[a2].xyz[2];

                            double d = sqrt((a2x-a1x)*(a2x-a1x) + (a2y-a1y)*(a2y-a1y) + (a2z-a1z)*(a2z-a1z));
					
                            if (d < vdw_cutoff && d > 0.0001){
                                double delta = d - vdw_cutoff; 

                                dx += beta*((delta/d)*(a2x-a1x));
								dy += beta*((delta/d)*(a2y-a1y));
							    dz += beta*((delta/d)*(a2z-a1z));

                                moved++;
                            }
                        }
                    }
                }

                pdb_entry.chains[c].residues[r].atoms[a].xyz[0] +=dx;
				pdb_entry.chains[c].residues[r].atoms[a].xyz[1] +=dy;
				pdb_entry.chains[c].residues[r].atoms[a].xyz[2] +=dz;
            }
        }
    }
    return moved;
}


void RefineRNA(ModelUnit &pdb_entry, map<string, map<string, vector<ModelUnit>>> &ideal_rna,
    const vector<pair<double,vector<size_t>>>&bp_vec, map<int, int> &bp_template_idx, map<int, int> &nonbp_template_idx, double alpha, double beta, double tolerance, const string outdir, const bool &logflag)
{
    int moved;

    Logger logger = Logger();
    logger.setPrint(false);
    if(logflag){
        logger.setVerbose(true);
        logger.setLogFile(outdir+"/logs/non_energy_refinement.log");
    }

    standardize_pdb_order(pdb_entry);

    SelectBPTemplate(pdb_entry, ideal_rna["bp"], bp_vec, bp_template_idx);

    for (int c=0; c<pdb_entry.chains.size(); c++){
        vector<vector<bool>> rescontacts;
        getContactResidueMap(pdb_entry, c, rescontacts);

        int length = pdb_entry.chains[c].residues.size();
        int iterations = round(length+3000);

        for (int t=0; t<iterations; t++){
            //moved = fix_bond_lengths(pdb_entry, tolerance);
            //moved += fix_adj_bond_lengths(pdb_entry, loop_vec, ideal_adj_rna, tolerance);
            //moved+=fixBaseConformation(pdb_entry, ideal_rna, bp_vec, bp_template_idx, nonbp_template_idx);
            //moved+=fix_clashes(pdb_entry);

            //moved=correct_geometry(pdb_entry, rescontacts, alpha, tolerance);
            //moved+=fixBaseConformation(pdb_entry, ideal_rna["bp"], bp_vec, bp_template_idx);
            //moved+=correct_clash(pdb_entry, rescontacts, beta);
            //moved+=fixBaseConformation(pdb_entry, ideal_rna["bp"], bp_vec, bp_template_idx);
            moved=correct_combine(pdb_entry, rescontacts, alpha, beta, tolerance);
            moved+=fixBaseConformation(pdb_entry, ideal_rna["bp"], bp_vec, bp_template_idx);

            logger.log("t="+to_string(t)+" moved="+to_string(moved));

            // if no atoms are moved, immediately break out of this for loop
            if (moved==0) break;

        }
    }

    logger.log("Final moved time="+to_string(moved)+"\n", true);
}


