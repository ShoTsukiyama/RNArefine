#include "MonteCarloSimulation.hpp"

void MonteCarloSimulation::updateenergies()
{
	if(weightenergy["InitPos"] > 0) updateeneinitpos();
	if(weightenergy["Init"] > 0) updateeneinit();
	if(weightenergy["Template"] > 0) updateenetemplate();
	if(weightenergy["Length energy"] > 0) updateenelength();
	if(weightenergy["Angle energy"] > 0) updateeneangle();
	if(weightenergy["Torsion energy"] > 0) updateenetorsion();
	if(weightenergy["Comb torsion energy"] > 0) updateenecombtorsion();
	if(weightenergy["Clash energy"] > 0) updateeneclash();
	if(weightenergy["Electros energy"] > 0) updateeneelectros();
	if(weightenergy["LJ energy"] > 0) updateenelj();
	if(weightenergy["Hbond energy"] > 0) updateenehbond();
	if(weightenergy["BP energy"] > 0) updateenebp();
	if(weightenergy["Stack energy"] > 0) updateenestack();
	if(weightenergy["Stack geometrical energy"] > 0) updateenestackgeo();
	if(weightenergy["Baseplane energy"] > 0) updateenebaseplane();
}

void MonteCarloSimulation::initializeatomenergy()
{
	int i,j;
	for(i=0;i<numbb;i++){
		char base=refinestr_temp[i].aaa;
		for(j=0;j<AccAtomParamNum[base][4];j++){
			lengthenergyatoms[i][j]=0;
			torsionenergyatoms[i][j]=0;
			angleenergyatoms[i][j]=0;
			pairwiseenergyatoms[i][j]=0;
		}
	}
}

double MonteCarloSimulation::calcenergymain()
{
	double totalene=0;
	bool dispflag=false;
	
	if(dispflag) cout<<"Calculate energy======================="<<endl;

	if(weightenergy["InitPos"] > 0){
		double initpos_ene = weightenergy["InitPos"] * calceneinitpos();
		totalene+=initpos_ene;
		if(dispflag) cout<<"InitPos: "<< initpos_ene <<endl;
	}

	if(weightenergy["Init"] > 0){
		double init_ene = weightenergy["Init"] * calceneinit();
		totalene+=init_ene;
		if(dispflag) cout<<"Init: "<< init_ene <<endl;
	}

	if(weightenergy["Template"] > 0){
		double template_ene = weightenergy["Template"] * calcenetemplate();
		totalene+=template_ene;
		if(dispflag) cout<<"Template: "<< template_ene <<endl;
	}
	
	if(weightenergy["Length energy"] > 0){
		double length_ene = weightenergy["Length energy"] * calcenelength();
		totalene+=length_ene;
		if(dispflag) cout<<"Length energy: "<< length_ene <<endl;
	}
	
	if(weightenergy["Angle energy"] > 0){
		double angle_ene = weightenergy["Angle energy"] * calceneangle();
		totalene+=angle_ene;
		if(dispflag) cout<<"Angle energy: "<< angle_ene <<endl;
	}
	
	if(weightenergy["Torsion energy"] > 0){
		double torsion_ene = weightenergy["Torsion energy"] * calcenetorsion();
		totalene+=torsion_ene;
		if(dispflag) cout<<"Torsion energy: "<< torsion_ene <<endl;
	}
	
	if(weightenergy["Comb torsion energy"] > 0){
		double combtorsion_ene = weightenergy["Comb torsion energy"] * calcenecombtorsion();
		totalene+=combtorsion_ene;
		if(dispflag) cout<<"Combination torsion energy: "<< combtorsion_ene <<endl;
	}
	
	if(weightenergy["Clash energy"] > 0){
		double clash_ene = weightenergy["Clash energy"] * calceneclash();
		totalene+=clash_ene;
		if(dispflag) cout<<"Clash energy: "<< clash_ene <<endl;
	}
	
	if(weightenergy["Electros energy"] > 0){
		double electros_ene = weightenergy["Electros energy"] * calceneelectros();
		totalene+=electros_ene;
		if(dispflag) cout<<"Electros energy: "<< electros_ene <<endl;
	}
	
	if(weightenergy["LJ energy"] > 0){
		double lj_ene = weightenergy["LJ energy"] * calcenelj();
		totalene+=lj_ene;
		if(dispflag) cout<<"LJ energy: "<< lj_ene <<endl;
	}
	
	if(weightenergy["Hbond energy"] > 0){
		double hbond_ene = weightenergy["Hbond energy"] * calcenehbond();
		totalene+=hbond_ene;
		if(dispflag) cout<<"Hbond energy: "<< hbond_ene <<endl;
	}
	
	if(weightenergy["BP energy"] > 0){
		double bp_ene = weightenergy["BP energy"] * calcenebp();
		totalene+=bp_ene;
		if(dispflag) cout<<"BP energy: "<< bp_ene <<endl;
	}

	if(weightenergy["Stack energy"] > 0){
		double stack_ene = weightenergy["Stack energy"] * calcenestack();
		totalene+=stack_ene;
		if(dispflag) cout<<"Stack energy: "<< stack_ene <<endl;
	}

	if(weightenergy["Stack geometrical energy"] > 0){
		double stackgeo_ene = weightenergy["Stack geometrical energy"] * calcenestackgeo();
		totalene+=stackgeo_ene;
		if(dispflag) cout<<"Stack geometrical energy: "<< stackgeo_ene <<endl;
	}

	if(weightenergy["Baseplane energy"] > 0){
		double baseplane_ene = weightenergy["Baseplane energy"] * calcenebaseplane();
		totalene+=baseplane_ene;
		if(dispflag) cout<<"Baseplane energy: "<< baseplane_ene <<endl;
	}
	
	return totalene;
}



























































