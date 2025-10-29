#include "MonteCarloSimulation.hpp"

double BondLengthMCeneobj::calcenergyblength() const
{
	if(*flag){
		double newene;
		double deltaene;

		if(*bleng < *lower || *upper < *bleng){
			double tmp = *bleng - *ideal;
			
			//if(sqrt(tmp*tmp) > 1.0) cout<< remark << "-" << sqrt(tmp*tmp)<<endl;

			return tmp*tmp;
		}else{
			return 0;
		}
	}else{
		return prevene;
	}
}

void BondLengthMCeneobj::setprevene(const double &_prevene) const
{
	prevene=_prevene;
	*atomenergy1+=_prevene;
	*atomenergy2+=_prevene;
}

double BondAngleMCeneobj::calcenergybangle() const
{
	if(*flag){
		double newene;
		double deltaene;

		if(*bangle < *lower || *upper < *bangle){
			double tmp = *bangle - *ideal;

			//if(sqrt(tmp*tmp) > 3.0) cout<< remark << "-" << *bangle << "-" << *ideal << "-" << sqrt(tmp*tmp) << "-" << *lower << "-" << *upper << endl;

			return tmp*tmp;
		}else{
			return 0;
		}
	}else{
		return prevene;
	}
}

void BondAngleMCeneobj::setprevene(const double &_prevene) const
{
	prevene=_prevene;
	*atomenergy1+=_prevene;
	*atomenergy2+=_prevene;
	*atomenergy3+=_prevene;
}

int getbin(double val, const vector<double> &boundary){
	for(int i=0;i<boundary.size();i++){
		if(val < boundary[i]) return i;
	}
	return boundary.size()-1;
}

double BondTorsioneneobj::calcenergybtorsion() const
{
	if(*flag){
		int binidx=getbin(*btorsion, *boundary);
		double diff=-(*energy)[binidx];
		return diff;
	}else{
		return prevene;
	}
}

void BondTorsioneneobj::setprevene(const double &_prevene) const
{
	prevene=_prevene;
	*atomenergy1+=_prevene;
	*atomenergy2+=_prevene;
	*atomenergy3+=_prevene;
	*atomenergy4+=_prevene;
}

double CombTorsionseneobj::calcenergycombtorsion() const
{
	if(*flag1 || *flag2){
		int binidx1=getbin(*btorsion1, *boundary1);
		int binidx2=getbin(*btorsion2, *boundary2);

		double diff=-(*energy)[binidx1*(*binnum2)+binidx2];
		return diff;
	}else{
		return prevene;
	}
}

void CombTorsionseneobj::setprevene(const double &_prevene) const
{
	prevene=_prevene;
	*atomenergy11+=_prevene;
	*atomenergy12+=_prevene;
	*atomenergy13+=_prevene;
	*atomenergy14+=_prevene;
	*atomenergy21+=_prevene;
	*atomenergy22+=_prevene;
	*atomenergy23+=_prevene;
	*atomenergy24+=_prevene;
}

double BondDihedralMCeneobj::calcenergybdihedral() const
{
	if(*flag){
		double newene;
		double deltaene;

		if(*bdangle < *lower || *upper < *bdangle){
			double tmp = *bdangle - *ideal;
			return sqrt(tmp*tmp);
		}else{
			return 0;
		}
	}else{
		return prevene;
	}
}

void BondDihedralMCeneobj::setprevene(const double &_prevene) const
{
	prevene=_prevene;
	*atomenergy1+=_prevene;
	*atomenergy2+=_prevene;
	*atomenergy3+=_prevene;
	*atomenergy4+=_prevene;
}

void MonteCarloSimulation::setgeometricaleneobj(){
	int i,j, idx;
	char base;
	int res1, res2, res3, res4;
	int atom1, atom2, atom3, atom4;

	BasePtrsA.resize(numbb);
	BaseAtomenergiesA.resize(numbb);

	for(i=0;i<numbb;i++)
	{
		base = refinestr_temp[i].aaa;

		if(weightenergy["Baseplane energy"])
		{
			for(int ii=AccAtomParamNum[base][1];ii<AccAtomParamNum[base][3];ii++){
				BasePtrsA[i].push_back(&refinestr_temp[i].ptv[ii]);
			}

			for(int ii=AccAtomParamNum[base][1];ii<AccAtomParamNum[base][3];ii++){
				BaseAtomenergiesA[i].push_back(&pairwiseenergyatoms[i][ii]);
			}

			BasePairMCobj2 basep = BasePairMCobj2(BasePtrsA[i], BaseAtomenergiesA[i], 1.0, &ntflag[i]);
			baseplane.push_back(basep);
		}
		
		for(j=0;j<15;j++)
		{
			if(i == 0 && j == 10){
				res1 = i; res2 = i; res3 = i; res4 = i;
				atom1=FIRSTOPParamsList[base][j][0]; atom2=FIRSTOPParamsList[base][j][1]; atom3=FIRSTOPParamsList[base][j][2]; atom4=FIRSTOPParamsList[base][j][3];
				
				if(weightenergy["Length energy"] > 0){
					BondLengthMCeneobj blmc = BondLengthMCeneobj::create(&refinestr_temp[i].len[j], &lengthenergyatoms[res3][atom3], &lengthenergyatoms[res4][atom4], &lengthflag[i][j], &MCOPLENGTHPARAM[base][0][0], &MCOPLENGTHPARAM[base][0][3], &MCOPLENGTHPARAM[base][0][2], to_string(i)+"_"+to_string(j)+"_"+to_string(res3)+"_"+to_string(res4)+"_"+to_string(atom3)+"_"+to_string(atom4));
					bondlengthmcene.push_back(blmc);
				}

				if(weightenergy["Angle energy"] > 0){
					BondAngleMCeneobj bamc = BondAngleMCeneobj::create(&refinestr_temp[i].ang[j], &angleenergyatoms[res2][atom2], &angleenergyatoms[res3][atom3], &angleenergyatoms[res4][atom4], &angleflag[i][j], &MCOPINNERANGPARAM[base][0][0], &MCOPINNERANGPARAM[base][0][3], &MCOPINNERANGPARAM[base][0][2], to_string(i)+"_"+to_string(j)+"_"+to_string(res2)+"_"+to_string(res3)+"_"+to_string(res4)+"_"+to_string(atom2)+"_"+to_string(atom3)+"_"+to_string(atom4));
					bondanglemcene.push_back(bamc);
				}
				
				if(weightenergy["Torsion energy"] > 0){
					BondTorsioneneobj btmc = BondTorsioneneobj::create(&refinestr_temp[i].dih[j], &torsionenergyatoms[res1][atom1], &torsionenergyatoms[res2][atom2], &torsionenergyatoms[res3][atom3], &torsionenergyatoms[res4][atom4], &torsionflag[i][j], &MCOPTORSIONDISTRIBUTION[base][0].boundary, &MCOPTORSIONDISTRIBUTION[base][0].energy);
					bondtorsionmcene.push_back(btmc);
				}
			}else if(i == 0 && j == 11){
				res1 = i; res2 = i; res3 = i; res4 = i;
				atom1=FIRSTOPParamsList[base][j][0]; atom2=FIRSTOPParamsList[base][j][1]; atom3=FIRSTOPParamsList[base][j][2]; atom4=FIRSTOPParamsList[base][j][3];
				
				if(weightenergy["Length energy"] > 0){
					BondLengthMCeneobj blmc = BondLengthMCeneobj::create(&refinestr_temp[i].len[j], &lengthenergyatoms[res3][atom3], &lengthenergyatoms[res4][atom4], &lengthflag[i][j], &MCOPLENGTHPARAM[base][1][0], &MCOPLENGTHPARAM[base][1][3], &MCOPLENGTHPARAM[base][1][2], to_string(i)+"_"+to_string(j)+"_"+to_string(res3)+"_"+to_string(res4)+"_"+to_string(atom3)+"_"+to_string(atom4));
					bondlengthmcene.push_back(blmc);
				}

				if(weightenergy["Angle energy"] > 0){
					BondAngleMCeneobj bamc = BondAngleMCeneobj::create(&refinestr_temp[i].ang[j], &angleenergyatoms[res2][atom2], &angleenergyatoms[res3][atom3], &angleenergyatoms[res4][atom4], &angleflag[i][j], &MCOPINNERANGPARAM[base][1][0], &MCOPINNERANGPARAM[base][1][3], &MCOPINNERANGPARAM[base][1][2], to_string(i)+"_"+to_string(j)+"_"+to_string(res2)+"_"+to_string(res3)+"_"+to_string(res4)+"_"+to_string(atom2)+"_"+to_string(atom3)+"_"+to_string(atom4));
					bondanglemcene.push_back(bamc);
				}

				if(weightenergy["Torsion energy"] > 0){
					BondTorsioneneobj btmc = BondTorsioneneobj::create(&refinestr_temp[i].dih[j], &torsionenergyatoms[res1][atom1], &torsionenergyatoms[res2][atom2], &torsionenergyatoms[res3][atom3], &torsionenergyatoms[res4][atom4], &torsionflag[i][j], &MCOPTORSIONDISTRIBUTION[base][1].boundary, &MCOPTORSIONDISTRIBUTION[base][1].energy);
					bondtorsionmcene.push_back(btmc);
				}
			}else{
				atom1 = ParamsList[base][j][0]; atom2 = ParamsList[base][j][1]; atom3 = ParamsList[base][j][2]; atom4 = ParamsList[base][j][3];
				if(weightenergy["Length energy"] > 0){
					if((i != 0) || (i == 0 && j != 0)){
						if(j == 0){
							res3 = i-1; res4 = i;
						}else{
							res3 = i; res4 = i;
						}
						BondLengthMCeneobj blmc = BondLengthMCeneobj::create(&refinestr_temp[i].len[j], &lengthenergyatoms[res3][atom3], &lengthenergyatoms[res4][atom4], &lengthflag[i][j], &MCLENGTHPARAM[base][j][0], &MCLENGTHPARAM[base][j][3], &MCLENGTHPARAM[base][j][2], to_string(i)+"_"+to_string(j)+"_"+to_string(res3)+"_"+to_string(res4)+"_"+to_string(atom3)+"_"+to_string(atom4));
						bondlengthmcene.push_back(blmc);
					}
				}
				
				if(weightenergy["Angle energy"] > 0){
					if((i != 0) || (i == 0 && j != 0 && j != 1)){
						if(j == 0){
							res2 = i-1; res3 = i-1; res4 = i;
						}else if(j == 1){
							res2 = i-1; res3 = i; res4 = i;
						}else if(j == 11){
							res2 = i-1; res3 = i; res4 = i;
						}else{
							res2 = i; res3 = i; res4 = i;
						}
						BondAngleMCeneobj bamc = BondAngleMCeneobj::create(&refinestr_temp[i].ang[j], &angleenergyatoms[res2][atom2], &angleenergyatoms[res3][atom3], &angleenergyatoms[res4][atom4], &angleflag[i][j], &MCINNERANGPARAM[base][j][0], &MCINNERANGPARAM[base][j][3], &MCINNERANGPARAM[base][j][2], to_string(i)+"_"+to_string(j)+"_"+to_string(res2)+"_"+to_string(res3)+"_"+to_string(res4)+"_"+to_string(atom2)+"_"+to_string(atom3)+"_"+to_string(atom4));
						bondanglemcene.push_back(bamc);
					}
				}

				if(weightenergy["Torsion energy"] > 0){
					if((i != 0) || (i == 0 && j != 0 && j != 1 && j != 2)){
						if(j == 0){
							res1 = i-1; res2 = i-1; res3 = i-1; res4 = i;
						}else if(j == 1){
							res1 = i-1; res2 = i-1; res3 = i; res4 = i;
						}else if(j == 2){
							res1 = i-1; res2 = i; res3 = i; res4 = i;
						}else if(j == 10){
							res1 = i-1; res2 = i; res3 = i; res4 = i;
						}else if(j == 11){
							res1 = i; res2 = i-1; res3 = i; res4 = i;
						}else{
							res1 = i; res2 = i; res3 = i; res4 = i;
						}
						BondTorsioneneobj btmc = BondTorsioneneobj::create(&refinestr_temp[i].dih[j], &torsionenergyatoms[res1][atom1], &torsionenergyatoms[res2][atom2], &torsionenergyatoms[res3][atom3], &torsionenergyatoms[res4][atom4], &torsionflag[i][j], &MCTORSIONDISTRIBUTION[base][j].boundary, &MCTORSIONDISTRIBUTION[base][j].energy);
						bondtorsionmcene.push_back(btmc);
					}
				}
			}
		}
		
		res1 = i; res2 = i; res3 = i; res4=i;
		for(j=15;j<15+BaseParamsList[base].size();j++)
		{
			idx=j-15;
			atom1=BaseParamsList[base][idx][0]; atom2=BaseParamsList[base][idx][1]; atom3=BaseParamsList[base][idx][2]; atom4=BaseParamsList[base][idx][3];
			
			if(weightenergy["Length energy"] > 0){
				BondLengthMCeneobj blmc = BondLengthMCeneobj::create(&refinestr_temp[i].len[j], &lengthenergyatoms[res3][atom3], &lengthenergyatoms[res4][atom4], &lengthflag[i][j], &MCBASELENGTHPARAM[base][idx][0], &MCBASELENGTHPARAM[base][idx][3], &MCBASELENGTHPARAM[base][idx][2], to_string(i)+"_"+to_string(j)+"_"+to_string(res3)+"_"+to_string(res4)+"_"+to_string(atom3)+"_"+to_string(atom4));
				bondlengthmcene.push_back(blmc);
			}
		
			if(weightenergy["Angle energy"] > 0){
				BondAngleMCeneobj bamc = BondAngleMCeneobj::create(&refinestr_temp[i].ang[j], &angleenergyatoms[res2][atom2], &angleenergyatoms[res3][atom3], &angleenergyatoms[res4][atom4], &angleflag[i][j], &MCBASEINNERANGPARAM[base][idx][0], &MCBASEINNERANGPARAM[base][idx][3], &MCBASEINNERANGPARAM[base][idx][2], "base_"+to_string(i)+"_"+to_string(j)+"_"+to_string(res2)+"_"+to_string(res3)+"_"+to_string(res4)+"_"+to_string(atom2)+"_"+to_string(atom3)+"_"+to_string(atom4));
				bondanglemcene.push_back(bamc);
			}
		
			if(weightenergy["Torsion energy"] > 0){
				BondTorsioneneobj btmc = BondTorsioneneobj::create(&refinestr_temp[i].dih[j], &torsionenergyatoms[res1][atom1], &torsionenergyatoms[res2][atom2], &torsionenergyatoms[res3][atom3], &torsionenergyatoms[res4][atom4], &torsionflag[i][j], &MCBASETORSIONDISTRIBUTION[base][idx].boundary, &MCBASETORSIONDISTRIBUTION[base][idx].energy);
				bondtorsionmcene.push_back(btmc);
			}
		}
		
		for(j=0;j<HydParamsList[base].size();j++)
		{
			atom1=HydParamsList[base][j][0]; atom2=HydParamsList[base][j][1]; atom3=HydParamsList[base][j][2]; atom4=HydParamsList[base][j][3];

			if(weightenergy["Length energy"] > 0){
				BondLengthMCeneobj blmc = BondLengthMCeneobj::create(&refinestr_temp[i].hydlen[j], &lengthenergyatoms[res3][atom3], &lengthenergyatoms[res4][atom4], &hydlengthflag[i][j], &MCHYDLENGTHPARAM[base][j][4], &MCHYDLENGTHPARAM[base][j][3], &MCHYDLENGTHPARAM[base][j][2], to_string(i)+"_"+to_string(j)+"_"+to_string(res3)+"_"+to_string(res4)+"_"+to_string(atom3)+"_"+to_string(atom4));
				bondlengthmcene.push_back(blmc);
			}
		
			if(weightenergy["Angle energy"] > 0){
				BondAngleMCeneobj bamc = BondAngleMCeneobj::create(&refinestr_temp[i].hydang[j], &angleenergyatoms[res2][atom2], &angleenergyatoms[res3][atom3], &angleenergyatoms[res4][atom4], &hydangleflag[i][j], &MCHYDINNERANGPARAM[base][j][4], &MCHYDINNERANGPARAM[base][j][3], &MCHYDINNERANGPARAM[base][j][2], "hyd_"+to_string(i)+"_"+to_string(j)+"_"+to_string(res2)+"_"+to_string(res3)+"_"+to_string(res4)+"_"+to_string(atom2)+"_"+to_string(atom3)+"_"+to_string(atom4));
				bondanglemcene.push_back(bamc);
			}
			
			if(weightenergy["Torsion energy"] > 0){
				BondTorsioneneobj btmc = BondTorsioneneobj::create(&refinestr_temp[i].hyddih[j], &torsionenergyatoms[res1][atom1], &torsionenergyatoms[res2][atom2], &torsionenergyatoms[res3][atom3], &torsionenergyatoms[res4][atom4], &hydtorsionflag[i][j], &MCHYDTORSIONDISTRIBUTION[base][j].boundary, &MCHYDTORSIONDISTRIBUTION[base][j].energy);
				bondtorsionmcene.push_back(btmc);
			}
		}
		
		if(weightenergy["Length energy"] > 0){
			for(j=0;j<AddLengthParamsList[base].size();j++)
			{
				if(i == 0 && j < 4) continue;
				if(j < 4){
					res2 = i-1; res3 = i;
				}else{
					res2 = i; res3 = i;
				}

				atom2=AddLengthParamsList[base][j][0]; atom3=AddLengthParamsList[base][j][1];

				BondLengthMCeneobj blmc = BondLengthMCeneobj::create(&refinestr_temp[i].addlen[j], &lengthenergyatoms[res2][atom2], &lengthenergyatoms[res3][atom3], &addlengthflag[i][j], &MCADDLENGTHPARAM[base][j][0], &MCADDLENGTHPARAM[base][j][3], &MCADDLENGTHPARAM[base][j][2], to_string(i)+"_"+to_string(j)+"_"+to_string(res2)+"_"+to_string(res3)+"_"+to_string(atom2)+"_"+to_string(atom3));
				addbondlengthmcene.push_back(blmc);
			}

			
			for(j=0;j<AddHydLengthParamsList[base].size();j++)
			{
				res2 = i; res3 = i;
				atom2=AddHydLengthParamsList[base][j][0]; atom3=AddHydLengthParamsList[base][j][1];

				BondLengthMCeneobj blmc = BondLengthMCeneobj::create(&refinestr_temp[i].addhydlen[j], &lengthenergyatoms[res2][atom2], &lengthenergyatoms[res3][atom3], &hydaddlengthflag[i][j], &MCHYDADDLENGTHPARAM[base][j][0], &MCHYDADDLENGTHPARAM[base][j][3], &MCHYDADDLENGTHPARAM[base][j][2], to_string(i)+"_"+to_string(j)+"_"+to_string(res2)+"_"+to_string(res3)+"_"+to_string(atom2)+"_"+to_string(atom3));
				addhydbondlengthmcene.push_back(blmc);
			}
			

			for(j=0;j<AddHydAngleParamsList[base].size();j++)
			{
				res2 = i; res3 = i; res4 = i;
				atom2=AddHydAngleParamsList[base][j][0]; atom3=AddHydAngleParamsList[base][j][1]; atom4=AddHydAngleParamsList[base][j][2];

				BondAngleMCeneobj bamc = BondAngleMCeneobj::create(&refinestr_temp[i].addhydang[j], &angleenergyatoms[res2][atom2], &angleenergyatoms[res3][atom3], &angleenergyatoms[res4][atom4], &hydaddangleflag[i][j], &MCHYDADDANGLEPARAM[base][j][0], &MCHYDADDANGLEPARAM[base][j][3], &MCHYDADDANGLEPARAM[base][j][2], "addhyd_"+to_string(i)+"_"+to_string(j)+"_"+to_string(res2)+"_"+to_string(res3)+"_"+to_string(res4)+"_"+to_string(atom2)+"_"+to_string(atom3)+"_"+to_string(atom4));
				addhydbondanglemcene.push_back(bamc);
			}
		}

		if(weightenergy["Angle energy"] > 0){
			for(j=0;j<AddAngleParamsList[base].size();j++)
			{
				if(i == 0 && j < 1) continue;
				if(j < 1){
					res1 = i-1; res2 = i; res3 = i;
				}else{
					res1 = i; res2 = i; res3 = i;
				}
				atom1=AddAngleParamsList[base][j][0]; atom2=AddAngleParamsList[base][j][1]; atom3=AddAngleParamsList[base][j][2];
				BondAngleMCeneobj bamc = BondAngleMCeneobj::create(&refinestr_temp[i].addang[j], &angleenergyatoms[res1][atom1], &angleenergyatoms[res2][atom2], &angleenergyatoms[res3][atom3], &addangleflag[i][j], &MCADDINNERANGPARAM[base][j][0], &MCADDINNERANGPARAM[base][j][3], &MCADDINNERANGPARAM[base][j][2], "add_"+to_string(i)+"_"+to_string(j)+"_"+to_string(res1)+"_"+to_string(res2)+"_"+to_string(res3)+"_"+to_string(atom1)+"_"+to_string(atom2)+"_"+to_string(atom3));
				addbondanglemcene.push_back(bamc);
			}
		}

		if(weightenergy["Comb torsion energy"] > 0){
			for(j=0;j<combtorsion[base].size();j++){
				int tor1 = combtorsion[base][j][0];
				int tor2 = combtorsion[base][j][1];

				int res11=i; int res12=i; int res13=i; int res14=i;
				int res21=i; int res22=i; int res23=i; int res24=i;

				int atom11=ParamsList[base][tor1][0];
				int atom12=ParamsList[base][tor1][1];
				int atom13=ParamsList[base][tor1][2];
				int atom14=ParamsList[base][tor1][3];
				int atom21=ParamsList[base][tor2][0];
				int atom22=ParamsList[base][tor2][1];
				int atom23=ParamsList[base][tor2][2];
				int atom24=ParamsList[base][tor2][3];

				CombTorsionseneobj cbtmc = CombTorsionseneobj::create(&refinestr_temp[i].dih[tor1], &refinestr_temp[i].dih[tor2], 
					&torsionenergyatoms[res11][atom11], &torsionenergyatoms[res12][atom12], &torsionenergyatoms[res13][atom13], &torsionenergyatoms[res14][atom14], 
					&torsionenergyatoms[res21][atom21], &torsionenergyatoms[res22][atom22], &torsionenergyatoms[res23][atom23], &torsionenergyatoms[res24][atom24], 
					&torsionflag[i][tor1], &torsionflag[i][tor2], &TORSIONPAIRSPARAM[base][j].boundary1, &TORSIONPAIRSPARAM[base][j].boundary2, &TORSIONPAIRSPARAM[base][j].energy, &TORSIONPAIRSPARAM[base][j].binnum1, &TORSIONPAIRSPARAM[base][j].binnum2);
				combtorsionmcene.push_back(cbtmc);
			}
		}	

		if(weightenergy["InitPos"] > 0){
			for(int m=0;m<template_target_atoms.size();m++){
				int k = template_target_atoms_idx[m];
				int ii = atomid2mcidx[base][k];
				if(m < (*MCinputtedstrallflagPtr)[i].size()){
					if((*MCinputtedstrallflagPtr)[i][m]){
						InitMCPosobj initposmc = InitMCPosobj::create(&refinestr_temp[i].ptv[ii], 20, &pairwiseenergyatoms[i][ii], (*MCinputtedstrallPtr).residues[i].atoms[m].coord.x, (*MCinputtedstrallPtr).residues[i].atoms[m].coord.y, (*MCinputtedstrallPtr).residues[i].atoms[m].coord.z, 1.0, &ntflag[i]);
						initposene.push_back(initposmc);
					}
				}
			}
		}
	}
}

double MonteCarloSimulation::calcenelength(){
	double totalene=0;
	double ene;
	lengthenergies.clear();

	for (std::list<BondLengthMCeneobj>::const_iterator it = bondlengthmcene.begin(); it != bondlengthmcene.end(); ++it) {
		ene=it->calcenergyblength();
		lengthenergies.push_back(ene);
		totalene+=ene;
	}

	for (std::list<BondLengthMCeneobj>::const_iterator it = addbondlengthmcene.begin(); it != addbondlengthmcene.end(); ++it) {
		ene=it->calcenergyblength();
		lengthenergies.push_back(ene);
		totalene+=ene;
	}

	for (std::list<BondLengthMCeneobj>::const_iterator it = addhydbondlengthmcene.begin(); it != addhydbondlengthmcene.end(); ++it) {
		ene=it->calcenergyblength();
		lengthenergies.push_back(ene);
		totalene+=ene;
	}

	return totalene;
}

void MonteCarloSimulation::updateenelength(){
	int count = 0;
	for (std::list<BondLengthMCeneobj>::const_iterator it = bondlengthmcene.begin(); it != bondlengthmcene.end(); ++it){
		it->setprevene(lengthenergies[count]);
		count += 1;
	}

	for (std::list<BondLengthMCeneobj>::const_iterator it = addbondlengthmcene.begin(); it != addbondlengthmcene.end(); ++it) {
		it->setprevene(lengthenergies[count]);
		count += 1;
	}

	for (std::list<BondLengthMCeneobj>::const_iterator it = addhydbondlengthmcene.begin(); it != addhydbondlengthmcene.end(); ++it) {
		it->setprevene(lengthenergies[count]);
		count += 1;
	}
}

double MonteCarloSimulation::calceneangle(){
	double totalene = 0;
	double ene;
	angleenergies.clear();

	for (std::list<BondAngleMCeneobj>::const_iterator it = bondanglemcene.begin(); it != bondanglemcene.end(); ++it) {
		ene=it->calcenergybangle();
		angleenergies.push_back(ene);
		totalene+=ene;
	}

	for (std::list<BondAngleMCeneobj>::const_iterator it = addbondanglemcene.begin(); it != addbondanglemcene.end(); ++it) {
		ene=it->calcenergybangle();
		angleenergies.push_back(ene);
		totalene+=ene;
	}

	for (std::list<BondAngleMCeneobj>::const_iterator it = addhydbondanglemcene.begin(); it != addhydbondanglemcene.end(); ++it) {
		ene=it->calcenergybangle();
		angleenergies.push_back(ene);
		totalene+=ene;
	}

	return totalene;
}

void MonteCarloSimulation::updateeneangle(){
	int count = 0;
	for (std::list<BondAngleMCeneobj>::const_iterator it = bondanglemcene.begin(); it != bondanglemcene.end(); ++it) {
		it->setprevene(angleenergies[count]);
		count += 1;
	}

	for (std::list<BondAngleMCeneobj>::const_iterator it = addbondanglemcene.begin(); it != addbondanglemcene.end(); ++it) {
		it->setprevene(angleenergies[count]);
		count += 1;
	}

	for (std::list<BondAngleMCeneobj>::const_iterator it = addhydbondanglemcene.begin(); it != addhydbondanglemcene.end(); ++it) {
		it->setprevene(angleenergies[count]);
		count += 1;
	}
}

double MonteCarloSimulation::calcenetorsion(){
	double totalene = 0;
	double ene;
	torsionenergies.clear();

	for (std::list<BondTorsioneneobj>::const_iterator it = bondtorsionmcene.begin(); it != bondtorsionmcene.end(); ++it) {
		ene=it->calcenergybtorsion();
		torsionenergies.push_back(ene);
		totalene+=ene;
	}

	return totalene;
}

void MonteCarloSimulation::updateenetorsion(){
	int count = 0;
	for (std::list<BondTorsioneneobj>::const_iterator it = bondtorsionmcene.begin(); it != bondtorsionmcene.end(); ++it) {
		it->setprevene(torsionenergies[count]);
		count += 1;
	}
}

double MonteCarloSimulation::calcenecombtorsion(){
	double totalene = 0;
	double ene;
	combtorsionenergies.clear();

	for (std::list<CombTorsionseneobj>::const_iterator it = combtorsionmcene.begin(); it != combtorsionmcene.end(); ++it) {
		ene=it->calcenergycombtorsion();
		combtorsionenergies.push_back(ene);
		totalene+=ene;
	}

	return totalene;
}

void MonteCarloSimulation::updateenecombtorsion(){
	int count = 0;
	for (std::list<CombTorsionseneobj>::const_iterator it = combtorsionmcene.begin(); it != combtorsionmcene.end(); ++it) {
		it->setprevene(combtorsionenergies[count]);
		count += 1;
	}
}




















