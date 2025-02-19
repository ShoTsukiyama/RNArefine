// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

void RNASimulation::selectbptemplate2(Residueatoms restr1, Residueatoms restr2, vector<Residueatoms> &dnttemplate, TemplateDistMap &bptemplatelength){
	char aat1=restr1.aaa;
	char aat2=restr2.aaa;

	dnttemplate.resize(2);

	int i, j;
	Residueatoms templatetemp1, templatetemp2;
	int atom;
	vector<vector<double>> RotMatix;
    vector<double> TranVect;
	vector<vector<double>> RotMatixtmp;
    vector<double> TranVecttmp;
	vector<double> xyztmp(3,0);
	Point3s vec;
	double minrmsd=100000000000;
	int minidx=-1;

	vector<double> tmp(3,0);
	vector<vector<double> > xyz_list1(AccAtomParamNum[restr1.aaa][2]+AccAtomParamNum[restr2.aaa][2],tmp);
    vector<vector<double> > xyz_list2(AccAtomParamNum[restr1.aaa][2]+AccAtomParamNum[restr2.aaa][2],tmp);
	vector<vector<double> > xyz_listtmp(15,tmp);

	//cout<<"bp template"<<endl;
	//for(int j=0; j<SNTTEMPLATEATOMS[aat1].size();j++){
		//atom=SNTTEMPLATEATOMS[aat1][j];
		//cout<< "0" << "-" << j << "-" << restr1.ptv[atom].x << "-" << restr1.ptv[atom].y << "-" << restr1.ptv[atom].z<<endl;
	//}
	//for(int j=0; j<SNTTEMPLATEATOMS[aat2].size();j++){
		//atom=SNTTEMPLATEATOMS[aat2][j];
		//cout<< "1" << "-" << j << "-" << restr2.ptv[atom].x << "-" << restr2.ptv[atom].y << "-" << restr2.ptv[atom].z<<endl;
	//}
	
	for(i=0;i<BPLENTTEMPLATE[aat1][aat2].size();i++){
		//for(int j=0; j<SNTTEMPLATEATOMS[aat].size();j++){
		int count = 0;
		for(j=0; j<AccAtomParamNum[restr1.aaa][2];j++){
			atom=SNTTEMPLATEATOMS[aat1][j];
			xyz_list1[count][0]=restr1.ptv[atom].x;
			xyz_list1[count][1]=restr1.ptv[atom].y;
			xyz_list1[count][2]=restr1.ptv[atom].z;

			xyz_list2[count][0]=BPLENTTEMPLATE[aat1][aat2][i][0].ptv[atom].x;
			xyz_list2[count][1]=BPLENTTEMPLATE[aat1][aat2][i][0].ptv[atom].y;
			xyz_list2[count][2]=BPLENTTEMPLATE[aat1][aat2][i][0].ptv[atom].z;
			count+=1;
		}

		for(j=0; j<AccAtomParamNum[restr2.aaa][2];j++){
			atom=SNTTEMPLATEATOMS[aat2][j];
			xyz_list1[count][0]=restr2.ptv[atom].x;
			xyz_list1[count][1]=restr2.ptv[atom].y;
			xyz_list1[count][2]=restr2.ptv[atom].z;

			xyz_list2[count][0]=BPLENTTEMPLATE[aat1][aat2][i][1].ptv[atom].x;
			xyz_list2[count][1]=BPLENTTEMPLATE[aat1][aat2][i][1].ptv[atom].y;
			xyz_list2[count][2]=BPLENTTEMPLATE[aat1][aat2][i][1].ptv[atom].z;
			count+=1;
		}

		bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

		count=0;
		//for(int j=0; j<SNTTEMPLATEATOMS[aat].size();j++){
		for(j=0; j<AccAtomParamNum[restr1.aaa][2];j++){
			atom=SNTTEMPLATEATOMS[aat1][j];
			bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);
			templatetemp1.ptv[atom].x=xyztmp[0];
			templatetemp1.ptv[atom].y=xyztmp[1];
			templatetemp1.ptv[atom].z=xyztmp[2];
			count+=1;
		}

		for(j=0; j<AccAtomParamNum[restr2.aaa][2];j++){
			atom=SNTTEMPLATEATOMS[aat2][j];
			bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);
			templatetemp2.ptv[atom].x=xyztmp[0];
			templatetemp2.ptv[atom].y=xyztmp[1];
			templatetemp2.ptv[atom].z=xyztmp[2];
			count+=1;
		}

		double rmsd=0;
		//for(int j=0; j<SNTTEMPLATEATOMS[aat].size();j++){
		for(j=0; j<AccAtomParamNum[restr1.aaa][2];j++){
			atom=SNTTEMPLATEATOMS[aat1][j];
			vec.x=restr1.ptv[atom].x-templatetemp1.ptv[atom].x;
			vec.y=restr1.ptv[atom].y-templatetemp1.ptv[atom].y;
			vec.z=restr1.ptv[atom].z-templatetemp1.ptv[atom].z;
			rmsd+=bf.norm(vec);
		}

		for(j=0; j<AccAtomParamNum[restr2.aaa][2];j++){
			atom=SNTTEMPLATEATOMS[aat2][j];
			vec.x=restr2.ptv[atom].x-templatetemp2.ptv[atom].x;
			vec.y=restr2.ptv[atom].y-templatetemp2.ptv[atom].y;
			vec.z=restr2.ptv[atom].z-templatetemp2.ptv[atom].z;
			rmsd+=bf.norm(vec);
		}

		if(minrmsd > rmsd){
			RotMatixtmp=RotMatix;
			TranVecttmp=TranVect;
			minidx=i;
			minrmsd=rmsd;
		}
	}

	//cout<<"============"<<endl;

	//for(int j=0; j<SNTTEMPLATEATOMS[aat1].size();j++){
	for(j=0; j<AccAtomParamNum[restr1.aaa][2];j++){
		atom=SNTTEMPLATEATOMS[aat1][j];
		xyz_listtmp[j][0]=BPLENTTEMPLATE[aat1][aat2][minidx][0].ptv[atom].x;
		xyz_listtmp[j][1]=BPLENTTEMPLATE[aat1][aat2][minidx][0].ptv[atom].y;
		xyz_listtmp[j][2]=BPLENTTEMPLATE[aat1][aat2][minidx][0].ptv[atom].z;

		bf.Transform(xyz_listtmp[j],RotMatixtmp,TranVecttmp,xyztmp);
		dnttemplate[0].ptv[atom].x=xyztmp[0];
		dnttemplate[0].ptv[atom].y=xyztmp[1];
		dnttemplate[0].ptv[atom].z=xyztmp[2];
		//cout<< "0" << "-" << j << "-" << dnttemplate[0].ptv[atom].x << "-" << dnttemplate[0].ptv[atom].y << "-" << dnttemplate[0].ptv[atom].z<<endl;
	}

	//for(int j=0; j<SNTTEMPLATEATOMS[aat2].size();j++){
	for(j=0; j<AccAtomParamNum[restr2.aaa][2];j++){
		atom=SNTTEMPLATEATOMS[aat2][j];
		xyz_listtmp[j][0]=BPLENTTEMPLATE[aat1][aat2][minidx][1].ptv[atom].x;
		xyz_listtmp[j][1]=BPLENTTEMPLATE[aat1][aat2][minidx][1].ptv[atom].y;
		xyz_listtmp[j][2]=BPLENTTEMPLATE[aat1][aat2][minidx][1].ptv[atom].z;

		bf.Transform(xyz_listtmp[j],RotMatixtmp,TranVecttmp,xyztmp);
		dnttemplate[1].ptv[atom].x=xyztmp[0];
		dnttemplate[1].ptv[atom].y=xyztmp[1];
		dnttemplate[1].ptv[atom].z=xyztmp[2];
		//cout<< "1" << "-" << j << "-" << dnttemplate[1].ptv[atom].x << "-" << dnttemplate[1].ptv[atom].y << "-" << dnttemplate[1].ptv[atom].z<<endl;
	}

	bptemplatelength=BPLENTTEMPLATELENGTH[aat1][aat2][minidx];

	//exit(1);
}

void RNASimulation::selectdnttemplate2(Residueatoms restr1, Residueatoms restr2, int basei, int basej, vector<Residueatoms> &dnttemplate, TemplateDistMap &dnttemplatelength){
	char aat1=restr1.aaa;
	char aat2=restr2.aaa;

	dnttemplate.resize(2);

	int i, j;
	Residueatoms templatetemp1, templatetemp2;
	int atom;
	vector<vector<double>> RotMatix;
    vector<double> TranVect;
	vector<vector<double>> RotMatixtmp;
    vector<double> TranVecttmp;
	vector<double> xyztmp(3,0);
	Point3s vec;
	double minrmsd=100000000000;
	int minidx=-1;

	vector<double> tmp(3,0);
	//vector<vector<double> > xyz_list1(SNTTEMPLATEATOMS[aat].size(),tmp);
    //vector<vector<double> > xyz_list2(SNTTEMPLATEATOMS[aat].size(),tmp);
	//vector<vector<double> > xyz_listtmp(SNTTEMPLATEATOMS[aat].size(),tmp);
	vector<vector<double> > xyz_list1(AccAtomParamNum[restr1.aaa][2]+AccAtomParamNum[restr2.aaa][2],tmp);
    vector<vector<double> > xyz_list2(AccAtomParamNum[restr1.aaa][2]+AccAtomParamNum[restr2.aaa][2],tmp);
	vector<vector<double> > xyz_listtmp(15,tmp);
	vector<vector<Residueatoms>> templatestructures;

	if(stpairednt[basei] != basej){
		templatestructures = DILENTTEMPLATE[aat1][aat2];
		//templatestructures = STACKTEMPLATE[2][aat1][aat2];
	}else{
		if(basei+1 == basej){
			templatestructures = STACKTEMPLATE[1][aat1][aat2];
		}else if(basej+1 == basei){
			templatestructures.resize(2);
			templatestructures[0] = STACKTEMPLATE[1][aat2][aat1][1];
			templatestructures[1] = STACKTEMPLATE[1][aat2][aat1][0];
		}else{
			templatestructures = STACKTEMPLATE[0][aat1][aat2];
		}
	}

	//cout<<"dnt template"<<endl;
	//for(int j=0; j<SNTTEMPLATEATOMS[aat1].size();j++){
		//atom=SNTTEMPLATEATOMS[aat1][j];
		//cout<< "0" << "-" << j << "-"  << restr1.ptv[atom].x << "-" << restr1.ptv[atom].y << "-" << restr1.ptv[atom].z<<endl;
	//}
	//for(int j=0; j<SNTTEMPLATEATOMS[aat2].size();j++){
		//atom=SNTTEMPLATEATOMS[aat1][j];
		//cout<< "1" << "-" << j << "-"  << restr2.ptv[atom].x << "-" << restr2.ptv[atom].y << "-" << restr2.ptv[atom].z<<endl;
	//}
	for(i=0;i<templatestructures.size();i++){
		int count = 0;
		//for(int j=0; j<SNTTEMPLATEATOMS[aat].size();j++){
		for(j=0; j<AccAtomParamNum[restr1.aaa][2];j++){
			atom=SNTTEMPLATEATOMS[aat1][j];
			xyz_list1[count][0]=restr1.ptv[atom].x;
			xyz_list1[count][1]=restr1.ptv[atom].y;
			xyz_list1[count][2]=restr1.ptv[atom].z;

			xyz_list2[count][0]=templatestructures[i][0].ptv[atom].x;
			xyz_list2[count][1]=templatestructures[i][0].ptv[atom].y;
			xyz_list2[count][2]=templatestructures[i][0].ptv[atom].z;
			count+=1;
		}

		for(j=0; j<AccAtomParamNum[restr2.aaa][2];j++){
			atom=SNTTEMPLATEATOMS[aat2][j];
			xyz_list1[count][0]=restr2.ptv[atom].x;
			xyz_list1[count][1]=restr2.ptv[atom].y;
			xyz_list1[count][2]=restr2.ptv[atom].z;

			xyz_list2[count][0]=templatestructures[i][1].ptv[atom].x;
			xyz_list2[count][1]=templatestructures[i][1].ptv[atom].y;
			xyz_list2[count][2]=templatestructures[i][1].ptv[atom].z;
			count+=1;
		}

		bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

		count = 0;
		//for(int j=0; j<SNTTEMPLATEATOMS[aat].size();j++){
		for(j=0; j<AccAtomParamNum[restr1.aaa][2];j++){
			atom=SNTTEMPLATEATOMS[aat1][j];
			bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);
			templatetemp1.ptv[atom].x=xyztmp[0];
			templatetemp1.ptv[atom].y=xyztmp[1];
			templatetemp1.ptv[atom].z=xyztmp[2];
			count+=1;
		}

		for(j=0; j<AccAtomParamNum[restr2.aaa][2];j++){
			atom=SNTTEMPLATEATOMS[aat2][j];
			bf.Transform(xyz_list2[count],RotMatix,TranVect,xyztmp);
			templatetemp2.ptv[atom].x=xyztmp[0];
			templatetemp2.ptv[atom].y=xyztmp[1];
			templatetemp2.ptv[atom].z=xyztmp[2];
			count+=1;
		}

		double rmsd=0;
		//for(int j=0; j<SNTTEMPLATEATOMS[aat].size();j++){
		for(j=0; j<AccAtomParamNum[restr1.aaa][2];j++){
			atom=SNTTEMPLATEATOMS[aat1][j];
			vec.x=restr1.ptv[atom].x-templatetemp1.ptv[atom].x;
			vec.y=restr1.ptv[atom].y-templatetemp1.ptv[atom].y;
			vec.z=restr1.ptv[atom].z-templatetemp1.ptv[atom].z;
			rmsd+=bf.norm(vec);
		}

		for(j=0; j<AccAtomParamNum[restr2.aaa][2];j++){
			atom=SNTTEMPLATEATOMS[aat2][j];
			vec.x=restr2.ptv[atom].x-templatetemp2.ptv[atom].x;
			vec.y=restr2.ptv[atom].y-templatetemp2.ptv[atom].y;
			vec.z=restr2.ptv[atom].z-templatetemp2.ptv[atom].z;
			rmsd+=bf.norm(vec);
		}

		if(minrmsd > rmsd){
			RotMatixtmp=RotMatix;
			TranVecttmp=TranVect;
			minidx=i;
			minrmsd=rmsd;
		}
	}

	//cout<<"==========="<<endl;

	for(int j=0; j<AccAtomParamNum[restr1.aaa][2];j++){
	//for(j=0; j<6;j++){
		atom=SNTTEMPLATEATOMS[aat1][j];
		xyz_listtmp[j][0]=templatestructures[minidx][0].ptv[atom].x;
		xyz_listtmp[j][1]=templatestructures[minidx][0].ptv[atom].y;
		xyz_listtmp[j][2]=templatestructures[minidx][0].ptv[atom].z;

		bf.Transform(xyz_listtmp[j],RotMatixtmp,TranVecttmp,xyztmp);
		dnttemplate[0].ptv[atom].x=xyztmp[0];
		dnttemplate[0].ptv[atom].y=xyztmp[1];
		dnttemplate[0].ptv[atom].z=xyztmp[2];
		//cout<< "0" << "-" << j << "-"  << dnttemplate[0].ptv[atom].x << "-" << dnttemplate[0].ptv[atom].y << "-" << dnttemplate[0].ptv[atom].z<<endl;
	}

	for(int j=0; j<AccAtomParamNum[restr1.aaa][2];j++){
	//for(j=0; j<6;j++){
		atom=SNTTEMPLATEATOMS[aat2][j];
		xyz_listtmp[j][0]=templatestructures[minidx][1].ptv[atom].x;
		xyz_listtmp[j][1]=templatestructures[minidx][1].ptv[atom].y;
		xyz_listtmp[j][2]=templatestructures[minidx][1].ptv[atom].z;

		bf.Transform(xyz_listtmp[j],RotMatixtmp,TranVecttmp,xyztmp);
		dnttemplate[1].ptv[atom].x=xyztmp[0];
		dnttemplate[1].ptv[atom].y=xyztmp[1];
		dnttemplate[1].ptv[atom].z=xyztmp[2];
		//cout<< "1" << "-" << j << "-" << dnttemplate[1].ptv[atom].x << "-" << dnttemplate[1].ptv[atom].y << "-" << dnttemplate[1].ptv[atom].z<<endl;
	}

	dnttemplatelength=DILENTTEMPLATELENGTH[aat1][aat2][minidx];
}


void RNASimulation::getcenterpoints(vector<Residueatoms> decstr, map<int, Point3s> &centerpoints, int st, int ed, vector<int> atomlist){
	int i,j;
	for(i=st;i<ed;i++){
		centerpoints[i].x=0; centerpoints[i].y=0; centerpoints[i].z=0;
		for(j=0;j<atomlist.size();j++){
			centerpoints[i].x+=decstr[i].ptv[atomlist[j]].x;
			centerpoints[i].y+=decstr[i].ptv[atomlist[j]].y;
			centerpoints[i].z+=decstr[i].ptv[atomlist[j]].z;
		}
		centerpoints[i].x/=6;
		centerpoints[i].x/=6;
		centerpoints[i].x/=6;
	}
}

double RNASimulation::eneinitpoints(vector<Residueatoms> recrrstr, map<int, Point3s> initpoints1, map<int, Point3s> initpoints2, int st, int ed){
	Point3s vec;
	int i,j;
	double rnorm;
	int mtype1, mtype2;
	double rmsd=0;
	int resnumst, resnum;
	int basei, basej;
	double diff;
	map<int, Point3s> centerpoints1, centerpoints2;
	vector<int> atomlist1 = {0,1,2,3,4,5};
	vector<int> atomlist2 = {3,4,7,8,9};

	getcenterpoints(recrrstr, centerpoints1, st, ed, atomlist1);
	getcenterpoints(recrrstr, centerpoints2, st, ed, atomlist2);
	
	for(i=st;i<ed;i++){
		vec.x = initpoints1[i].x-centerpoints1[i].x; 
		vec.y = initpoints1[i].y-centerpoints1[i].y; 
		vec.z = initpoints1[i].z-centerpoints1[i].z;
		rmsd+=bf.norm(vec);

		vec.x = initpoints2[i].x-centerpoints2[i].x; 
		vec.y = initpoints2[i].y-centerpoints2[i].y; 
		vec.z = initpoints2[i].z-centerpoints2[i].z;
		rmsd+=bf.norm(vec);
	}

	return rmsd;
}

void RNASimulation::eneedgegeoparams(vector<Residueatoms> recrrstr, int ed, int numseq, double &feneangle, double &fenelength){
	bool flagtor;
	char aat;
	double feneangle1, feneangle2;
	int ieneangle1, ieneangle2, ienelength;

	if(ed < numseq){
		coord2parambb(recrrstr,ed,ed+1,1);
		aat=recrrstr[ed].aaa;
		
		feneangle1=0; ieneangle1=0; feneangle2=0; ieneangle2=0; fenelength=0; ienelength=0;
		getouteneangle(recrrstr[ed].ang[INNERANGLETARGET[0]], INNERANGPARAM[aat][INNERANGLETARGET[0]].mean, INNERANGPARAM[aat][INNERANGLETARGET[0]].lower, INNERANGPARAM[aat][INNERANGLETARGET[0]].upper, INNERANGPARAM[aat][INNERANGLETARGET[0]].slower, INNERANGPARAM[aat][INNERANGLETARGET[0]].supper, 1, feneangle1, ieneangle1);
		getouteneangle(recrrstr[ed].ang[INNERANGLETARGET[1]], INNERANGPARAM[aat][INNERANGLETARGET[1]].mean, INNERANGPARAM[aat][INNERANGLETARGET[1]].lower, INNERANGPARAM[aat][INNERANGLETARGET[1]].upper, INNERANGPARAM[aat][INNERANGLETARGET[1]].slower, INNERANGPARAM[aat][INNERANGLETARGET[1]].supper, 1, feneangle2, ieneangle2);

		getoutenelength(recrrstr[ed].len[BLENGTHTARGET[0]], LENGTHPARAM[aat][BLENGTHTARGET[0]].mean, LENGTHPARAM[aat][BLENGTHTARGET[0]].lower, LENGTHPARAM[aat][BLENGTHTARGET[0]].upper, LENGTHPARAM[aat][BLENGTHTARGET[0]].slower, LENGTHPARAM[aat][BLENGTHTARGET[0]].supper, 1, fenelength, ienelength);
	
		feneangle=feneangle1+feneangle2;
		fenelength=fenelength;
	}else{
		feneangle=0;
		fenelength=0;
	}
}

void RNASimulation::eneaddgeoparams(vector<Residueatoms> recrrstr, int st, int ed, double &feneangle, double &fenelength)
{
	int i, j;
	char aat;
	double feneangletmp, fenelengthtmp;
	int ieneangletmp, ienelengthtmp;

	feneangle=0.0; fenelength=0.0;
	coord2paramadd(recrrstr, st, ed);

	for(i=st;i<ed;i++){
		aat=recrrstr[i].aaa;
		for(j=0;j<ADDBLENGTHTARGET.size();j++)
		{
			getoutenelength(recrrstr[i].addlen[ADDBLENGTHTARGET[j]], ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[j]].mean, ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[j]].lower, ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[j]].upper, ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[j]].slower, ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[j]].supper, 1, fenelengthtmp, ienelengthtmp);
			fenelength+=fenelengthtmp;
		}

		for(j=0;j<ADDINNERANGLETARGET.size();j++)
		{
			getouteneangle(recrrstr[i].addang[ADDINNERANGLETARGET[j]], ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].mean, ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].lower, ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].upper, ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].slower, ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].supper, 1, feneangletmp, ieneangletmp);
			feneangle+=feneangletmp;
		}
	}
}

void RNASimulation::calbplength(Residueatoms res1, Residueatoms res2, vector<vector<double>> &interdistmap){
	Point3s vec;
	for(int i=0;i<PAIRLENGTHTARGETMOVNUM;i++){
		for(int j=0;j<PAIRLENGTHTARGETMOVNUM;j++){
			vec.x=res1.ptv[i].x-res2.ptv[j].x;
			vec.y=res1.ptv[i].y-res2.ptv[j].y;
			vec.z=res1.ptv[i].z-res2.ptv[j].z;

			interdistmap[i][j]=bf.norm(vec);
		}
	}
}

double RNASimulation::compbplengthmap(vector<vector<double>> map1, vector<vector<double>> map2){
	double tmp;
	double rmsd = 0.0;
	for(int i=0;i<map1.size();i++){
		for(int j=0;j<map1.size();j++){
			tmp=map1[i][j]-map2[i][j];
			rmsd+=sqrt(tmp*tmp);
		}
	}
	return rmsd;
}

double RNASimulation::enebpplane(vector<Residueatoms> decstr, int basei, int basej){
	double a, b, c, d;
	Point3s pointA, pointB, pointC;
	char aat1, aat2;
	double ene=0.0;
	aat1 = decstr[basei].aaa;
	aat2 = decstr[basej].aaa;

	vector<int> baseatoms1=BaseHeavyAtom[aat1];
	vector<int> baseatoms2=BaseHeavyAtom[aat2];

	pointA=decstr[basej].ptv[baseatoms2[0]]; pointB=decstr[basej].ptv[baseatoms2[1]]; pointC=decstr[basej].ptv[baseatoms2[2]];
	getplanebp(pointA, pointB, pointC, a, b, c, d);
	
	for(int i=0; i<3;i++){
		ene+=getdistfromplane(a, b, c, d, decstr[basei].ptv[baseatoms1[i]]);
	}

	return ene;
}

double RNASimulation::enebptor(vector<Residueatoms> decstr, int basei, int basej){
	int atom1, atom2, atom3, atom4;
	char aat1, aat2;
	double dih;
	Point3s tmp1, tmp2, tmp3, tmp4;
	double indftmp=0.0; 
	int indntmp=0.0;
	aat1=decstr[basei].aaa;
	aat2=decstr[basej].aaa;

	atom1=BasePairTorList[aat1][aat2][0][0];
	atom2=BasePairTorList[aat1][aat2][0][1];
	atom3=BasePairTorList[aat1][aat2][0][2];
	atom4=BasePairTorList[aat1][aat2][0][3];
	tmp1=decstr[basei].ptv[atom1];tmp2=decstr[basei].ptv[atom2];tmp3=decstr[basej].ptv[atom3];tmp4=decstr[basej].ptv[atom4];
	dih=bf.phi(tmp1.x,tmp1.y,tmp1.z,tmp2.x,tmp2.y,tmp2.z,tmp3.x,tmp3.y,tmp3.z,tmp4.x,tmp4.y,tmp4.z);
	if(dih>360) dih=dih-360;
	else if(dih<0) dih=dih+360;
	getoutenetorangle(dih, BASEPAIRPARAM[1][0].mean, BASEPAIRPARAM[1][0].lower, BASEPAIRPARAM[1][0].upper, 1, indftmp, indntmp);

	return indftmp;
}

double RNASimulation::enebptemplatelength(vector<Residueatoms> recrrstr, vector<pair<int,int>> bpairs, map<int, TemplateDistMap> bptemplateslength){

	int i, basei, basej;
	double rmsd=0.0;
	vector<vector<double>> interdistmap;
	interdistmap.clear();
	interdistmap.resize(PAIRLENGTHTARGETMOVNUM, vector<double>(PAIRLENGTHTARGETMOVNUM, 0.0));
	for(i=0; i<bpairs.size();i++){
		basei=bpairs[i].first;
		basej=bpairs[i].second;
		if(basej >= 0){
			calbplength(recrrstr[basei], recrrstr[basej], interdistmap);
			rmsd+=compbplengthmap(interdistmap, bptemplateslength[basei].interdistmap);
			//rmsd+=enebpplane(recrrstr, basei, basej);
			//rmsd+=enebptor(recrrstr, basei, basej);
		}
	}
	return rmsd;
}

double RNASimulation::eneadjtemplatelength(vector<Residueatoms> recrrstr, int st, int ed, int numseq, map<int, TemplateDistMap> dnttemplateslength){

	int i, basei, basej;
	char aat1, aat2;
	double rmsd=0.0;
	vector<vector<double>> interdistmap;
	interdistmap.clear();
	interdistmap.resize(PAIRLENGTHTARGETMOVNUM, vector<double>(PAIRLENGTHTARGETMOVNUM, 0.0));

	if(st > 0){
		basei=st-1; basej=st;
		aat1 = recrrstr[basei].aaa; aat2 = recrrstr[basej].aaa;
		calbplength(recrrstr[basei], recrrstr[basej], interdistmap);
		rmsd+=compbplengthmap(interdistmap, dnttemplateslength[basei].interdistmap);
	}
	for(i=st;i<ed-1;i++){
		basei=i; basej=i+1;
		aat1 = restr[basei].aaa; aat2 = restr[basej].aaa;
		calbplength(recrrstr[basei], recrrstr[basej], interdistmap);
		rmsd+=compbplengthmap(interdistmap, dnttemplateslength[basei].interdistmap);
	}
	if(ed < numseq){
		basei=ed-1; basej=ed;
		aat1 = restr[basei].aaa; aat2 = restr[basej].aaa;
		calbplength(recrrstr[basei], recrrstr[basej], interdistmap);
		rmsd+=compbplengthmap(interdistmap, dnttemplateslength[basei].interdistmap);
	}

	return rmsd;
}

double RNASimulation::distfromene2(vector<Residueatoms> recrrstr, vector<Residueatoms> decstr, map<int, Point3s> initpoints1, map<int, Point3s> initpoints2, vector<pair<int,int>> bpairs, map<int, TemplateDistMap> bptemplateslength, map<int, TemplateDistMap> dnttemplateslength, int st, int ed, int numseq){
	double ene1, ene2, ene3, ene4, ene5, ene6, ene7;
	ene1=eneinitpoints(recrrstr, initpoints1, initpoints2, st, ed);

	//ene1=distfromene(recrrstr, decstr, st, ed, 1);
	eneedgegeoparams(recrrstr, ed, numseq, ene2, ene3);
	eneaddgeoparams(recrrstr, st, ed, ene4, ene5);
	ene6=enebptemplatelength(recrrstr, bpairs, bptemplateslength);
	ene7=eneadjtemplatelength(recrrstr, st, ed, numseq, dnttemplateslength);

	//cout<<"============"<<endl;
	//cout<<ene1<<endl;
	//cout<<10*ene2<<endl;
	//cout<<100*ene3<<endl;
	//cout<<5*ene4<<endl;
	//cout<<100*ene5<<endl;
	//cout<<ene6<<endl;
	//cout<<ene7<<endl;

	//return ene1+10*ene2+100*ene3+5*ene4+100*ene5+100*ene6;

	return ene1+10*ene2+100*ene3+5*ene4+100*ene5+ene6+ene7;
}


double RNASimulation::goldenSearchcalcene2(vector<Residueatoms> recrrstr, vector<Residueatoms> decstr, map<int, Point3s> initpoints1, map<int, Point3s> initpoints2, vector<pair<int,int>> bpairs, map<int, TemplateDistMap> bptemplateslength, map<int, TemplateDistMap> dnttemplateslength, double torval, int dihtype, int trandpos, int trandposst, int trandposed, int numseq){
	vector<Residueatoms> restrtmp;
	bool flagtor;
	restrtmp.resize(recrrstr.size());
	restrtmp=recrrstr;

	restrtmp[trandpos].dih[dihtype]=torval;
	flagtor=param2coordbb(restrtmp,trandpos,trandposed+1);
	param2coordsugar(restrtmp,trandpos,trandposed+1);
	param2coordbase(restrtmp,trandpos,trandposed+1);

	double ene=distfromene2(restrtmp, decstr, initpoints1, initpoints2, bpairs, bptemplateslength, dnttemplateslength, trandposst, trandposed+1, numseq);

	return ene;
}

void RNASimulation::initgeo2(vector<Residueatoms> &recrrstr, int st, int ed, int type){
	int i,j;
	char aat;
	bool flagtor;
	double randoval;
	int ist;

	if(st == 0){
		aat=recrrstr[0].aaa;
		for(j=0;j<INNERANGLETARGET.size();j++){
			if(lp1.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][0] > 0){
				//randoval=Random();
				if(j == 10){
					if(checkoutlinerlgm(recrrstr[0].ang[j], OPINNERANGPARAM[aat][0].lower, OPINNERANGPARAM[aat][0].upper)){
					//if(randoval < 0.1){
						//cout<<"aaa"<<endl;
						//recrrstr[0].ang[j]=getrandomnewparam(OPINNERANGPARAM[aat][0].lower, OPINNERANGPARAM[aat][0].upper);
						recrrstr[0].ang[j]=getmovangval(recrrstr[0].ang[j], OPINNERANGPARAM[aat][0].mean, OPINNERANGPARAM[aat][0].std, 0.5);
					//}
					}
				}else if(j == 11){
					if(checkoutlinerlgm(recrrstr[0].ang[j], OPINNERANGPARAM[aat][1].lower, OPINNERANGPARAM[aat][1].upper)){
					//if(randoval < 0.1){
						//cout<<"aaa"<<endl;
						//recrrstr[0].ang[j]=getrandomnewparam(OPINNERANGPARAM[aat][1].lower, OPINNERANGPARAM[aat][1].upper);
						recrrstr[0].ang[j]=getmovangval(recrrstr[0].ang[j], OPINNERANGPARAM[aat][1].mean, OPINNERANGPARAM[aat][1].std, 0.5);
					//}
					}
				}else{
					if(checkoutlinerlgm(recrrstr[0].ang[j], INNERANGPARAM[aat][j].lower, INNERANGPARAM[aat][j].upper)){
					//if(randoval < 0.1){
						//cout<<"aaa"<<endl;
						//recrrstr[0].ang[j]=getrandomnewparam(INNERANGPARAM[aat][j].lower, INNERANGPARAM[aat][j].upper);
						recrrstr[0].ang[j]=getmovangval(recrrstr[0].ang[j], INNERANGPARAM[aat][j].mean, INNERANGPARAM[aat][j].std, 0.5);
					//}
					}
				}
			}
		}

		for(j=0;j<BLENGTHTARGET.size();j++){
			//if(lp1.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][0] > 0){
				//randoval=Random();
				if(j == 10){
					if(checkoutlinerlgm(recrrstr[0].len[j], OPLENGTHPARAM[aat][0].lower, OPLENGTHPARAM[aat][0].upper)){
					//if(randoval < 0.1){
						//cout<<"aaa"<<endl;
						//recrrstr[0].len[j]=getrandomnewparam(OPLENGTHPARAM[aat][0].lower, OPLENGTHPARAM[aat][0].upper);
						recrrstr[0].len[j]=getmovlenval(recrrstr[0].len[j], OPLENGTHPARAM[aat][0].mean, OPLENGTHPARAM[aat][0].std, 1);
					//}
					}
				}else if(j == 11){
					if(checkoutlinerlgm(recrrstr[0].len[j], OPLENGTHPARAM[aat][1].lower, OPLENGTHPARAM[aat][1].upper)){
					//if(randoval < 0.1){
						//cout<<"aaa"<<endl;
						//recrrstr[0].len[j]=getrandomnewparam(OPLENGTHPARAM[aat][1].lower, OPLENGTHPARAM[aat][1].upper);
						recrrstr[0].len[j]=getmovlenval(recrrstr[0].len[j], OPLENGTHPARAM[aat][1].mean, OPLENGTHPARAM[aat][1].std, 1);
					//}
					}
				}else{
					if(checkoutlinerlgm(recrrstr[0].len[j], LENGTHPARAM[aat][j].lower, LENGTHPARAM[aat][j].upper)){
					//if(randoval < 0.1){
						//cout<<"aaa"<<endl;
						//recrrstr[0].len[j]=getrandomnewparam(LENGTHPARAM[aat][j].lower, LENGTHPARAM[aat][j].upper);
						recrrstr[0].len[j]=getmovlenval(recrrstr[0].len[j], LENGTHPARAM[aat][j].mean, LENGTHPARAM[aat][j].std, 1);
					//}
					}
				}
			//}
		}



		//flagtor=param2coordbb(recrrstr,0,1);
		//if(type > 0){
			//param2coordsugar(recrrstr,0,1);
			//param2coordbase(recrrstr,0,1);
			//param2coordhydrogen(recrrstr,0,1);
		//}
		ist=1;
	}else{
		ist=st;
	}

	for(int i=ist;i<ed;i++){
		aat=recrrstr[i].aaa;
		for(j=0;j<INNERANGLETARGET.size();j++){
			randoval=Random();
			//if(lp1.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][i] > 0){
				if(checkoutlinerlgm(recrrstr[i].ang[j], INNERANGPARAM[aat][j].lower, INNERANGPARAM[aat][j].upper)){
				//if(randoval < 0.1){
					//cout<<"aaa"<<endl;
					//recrrstr[i].ang[j]=getrandomnewparam(INNERANGPARAM[aat][j].lower, INNERANGPARAM[aat][j].upper);
					recrrstr[i].ang[j]=getmovangval(recrrstr[i].ang[j], INNERANGPARAM[aat][j].mean, INNERANGPARAM[aat][j].std, 0.5);
				//}
				}
			//}
		}

		for(j=0;j<BLENGTHTARGET.size();j++){
			randoval=Random();
			//if(lp1.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][i] > 0){
				if(checkoutlinerlgm(recrrstr[i].len[j], LENGTHPARAM[aat][j].lower, LENGTHPARAM[aat][j].upper)){
				//if(randoval < 0.1){
					//cout<<"aaa"<<endl;
					//recrrstr[i].len[j]=getrandomnewparam(LENGTHPARAM[aat][j].lower, LENGTHPARAM[aat][j].upper);
					recrrstr[i].len[j]=getmovlenval(recrrstr[i].len[j], LENGTHPARAM[aat][j].mean, LENGTHPARAM[aat][j].std, 1);
				//}
				}
			//}
		}
	}

	flagtor=param2coordbb(recrrstr,st,ed);
	if(type > 0){
		param2coordsugar(recrrstr,st,ed);
		param2coordbase(recrrstr,st,ed);
		param2coordhydrogen(recrrstr,st,ed);
	}
}

bool RNASimulation::fmovementlgm2(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc)
{
	int i,j,k,inext;
	int istart,iend;
	int trandpos, trandposed, trandposst;
	int wpos;
	char aat;
	double randval=Random();
	bool exitOuterLoop = false;
	map<int, vector<Point3s>> snttemplatestr;
	int lentype, angtype, dihtype;
	int basei, basej, tmp;
	bool flagtor;
	int dihcls;
	//double tolerance=1;
	int pairwiseidx;
	int span;

	map<string, double> enecoeff = {
    	//{"inittemp", 10},
		//{"adjtemplength", 1},
		//{"bptemplength", 1},
		//{"angle", 2},
		//{"length", 20}

		{"inittemp", 0.1},
		{"adjtemplength", 0.1},
		{"bptemplength", 0},
		{"angle", 1},
		{"length", 100}
	};

	randval=Random();

	//cout<<"aaa"<<endl;

	int totwrong1=lp1.nn[CLASH_ENERGY_ST]*2;
	int totwrong2=lp1.nn[BASE_ENERGY_ST+0] + lp1.nn[BASE_ENERGY_ST+1];
	int totwrong3=0;
	for(i=0;i<BLENGTHTARGET.size();i++){
		totwrong3+=lp1.nn[BLENGTH_ENERGY_ST+BLENGTHTARGET[i]];
		//totwrong3-=lp1.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[i]][0];
	}
	for(i=0;i<INNERANGLETARGET.size();i++){
		totwrong3+=lp1.nn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[i]];
		//totwrong3-=lp1.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[i]][0];
	}

	//cout<<"bbb"<<endl;

	if(totwrong3>0 && randval<0.5){
		//cout<<"a"<<endl;
		wpos=int(Random()*totwrong3-EPS);
		i=0;
		for(trandpos=0;trandpos<numseq;trandpos++)
		{
			for(j=0;j<BLENGTHTARGET.size();j++){
				inext=i+lp1.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][trandpos];
				if(i<=wpos && wpos<inext){
					exitOuterLoop = true;
					break;
				}
				i=inext;
			}

			for(j=0;j<INNERANGLETARGET.size();j++){
				inext=i+lp1.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][trandpos];
				if(i<=wpos && wpos<inext){
					exitOuterLoop = true;
					break;
				}
				i=inext;
			}
			if(exitOuterLoop) break;
		}
	}else if(totwrong1>0 && randval>=0.5 && randval<0.7){
		//cout<<"b"<<endl;
		wpos=int(Random()*totwrong1-EPS);

		i=0;
		for(trandpos=0;trandpos<numseq;trandpos++)
		{
			inext=i+lp1.indn[CLASH_ENERGY_ST][trandpos];
			if(i<=wpos && wpos<inext){
				break;
			}
			i=inext;
		}
	}else if(totwrong2>0 && randval>=0.7 && randval<0.9){
		//cout<<"c"<<endl;
		wpos=int(Random()*totwrong2-EPS);

		i=0;
		for(trandpos=0;trandpos<BASEPAIRS.size();trandpos++)
		{
			basei=BASEPAIRS[trandpos].basepidx.first;
			basej=BASEPAIRS[trandpos].basepidx.second;

			if(basei > basej){
				tmp=basei;
				basei=basej;
				basej=tmp;
			}
			pairwiseidx=bpmapid2idx[basei*numseq+basej];

			inext=i+bpenergy[pairwiseidx].pcw[0]+bpenergy[pairwiseidx].pcw[1];
			if(i<=wpos && wpos<inext){
				break;
			}
			i=inext;
		}

		randval=Random();
		if(randval > 0.5) trandpos=basei;
		else trandpos=basej;
	}else{
		//cout<<"d"<<endl;
		trandpos=int((numseq-2)*Random())+1;
	}

	//cout<<"eee"<<endl;

	span=2+int(2*Random());
	trandposst=trandpos-span/2;
	trandposed=trandposst+span;
	if(trandposst<1) trandposst=1;
	if(trandposed>numseq-1) trandposed=numseq-1;
	span=trandposed-trandposst;
	if(span < 2){
		newenergy=oldenergy;
		return false;
	}

	//cout<<"fff"<<endl;
	//cout<<trandposst << "-" << trandposed<<endl;

	//始端または終端ノードがbase pairに含まれないようにする
	vector<pair<int,int>> pairsf; 
	for(i=trandposst;i<trandposed;i++){
		if(bppairednt[i] == -1){
			pairsf.push_back({i,-1});
		}else if(bppairednt[i] < 1){
			break;
		}else if(bppairednt[i] >= numseq-1){
			break;
		}else{
			pairsf.push_back({i,bppairednt[i]});
		}
	}

	//cout<<"ggg"<<endl;

	vector<pair<int,int>> pairsb;
	for(i=trandposed-1;i>=trandposst;i--){
		if(bppairednt[i] == -1){
			pairsb.push_back({i,-1});
		}else if(bppairednt[i] < 1){
			break;
		}else if(bppairednt[i] >= numseq-1){
			break;
		}else{
			pairsb.push_back({i,bppairednt[i]});
		}
	}

	//cout<<"hhh"<<endl;

	vector<pair<int,int>> pairs;
	if(pairsf.size() >= pairsb.size()){
		pairs=pairsf;
	}else{
		for (i=pairsb.size()-1; i>=0; i--) {
			pairs.push_back(pairsb[i]);
		}
	}

	if(pairs.size() == 0){
            newenergy=oldenergy;
            return false;
	}

	//cout<<"iii"<<endl;

	copystr(decstr, numseq);

        //cout<<"jjj"<<endl;

	char aat1, aat2;
	vector<vector<int>> lmpparamidx;
	vector<bool> nodemflag;
	vector<pair<int, int>> edges;
	vector<int> edgetype;
	vector<int> paraminfo;
	vector<int> mnodes, nodes;
	map<int, vector<Residueatoms>> dnttemplates;
	map<int, TemplateDistMap> dnttemplateslength;
	map<int, vector<Residueatoms>> bptemplates;
	map<int, TemplateDistMap> bptemplateslength;

	//cout<<"kkk"<<endl;

	setlmpbpparam(restr, pairs, lmpparamidx);
	getgraphinfoLMP(pairs, numseq, lmpparamidx, nodemflag, edges, edgetype, paraminfo);

	//cout<<"lll"<<endl;
	//cout<<pairs.size()<<endl;

	for(i=0;i<pairs.size();i++)
	{
		basei = pairs[i].first; basej = pairs[i].second;
		if(basej == -1){
			nodes.push_back(basei);
		}else{
			nodes.push_back(basei);
			nodes.push_back(basej);
		}

		if(nodemflag[basei] & find(mnodes.begin(), mnodes.end(), basei) == mnodes.end()){
			mnodes.push_back(basei);
		}

		if((basej != -1) & nodemflag[basej] & find(mnodes.begin(), mnodes.end(), basej) == mnodes.end()){
			mnodes.push_back(basej);
		}
	}

	//cout<<"mmm"<<endl;


	for(i=0;i<edges.size();i++)
	{
		basei=edges[i].first; 
		basej=edges[i].second;

		if(edgetype[i] == 1){
			//if(basei < basej){
			selectdnttemplate(restr[basei], restr[basej], basei, basej, dnttemplates[i], dnttemplateslength[i]);
			//}else{
				//selectdnttemplate(restr[basej], restr[basei], basei, basej, dnttemplates[i], dnttemplateslength[i]);
			//}
		}
	}

	//cout<<"nnn"<<endl;


	for(i=0;i<edges.size();i++){
		basei=edges[i].first; basej=edges[i].second;
		aat1 = restr[basei].aaa; aat2 = restr[basej].aaa;

		if(edgetype[i] == 2){
			selectbptemplate(restr[basei], restr[basej], bptemplates[i], bptemplateslength[i]);
		}
	}

	//cout<<"ooo"<<endl;

	/*
	for(i=0;i<mnodes.size();i++){
		initrandommovegss(restr[mnodes[i]]);
		if(mnodes[i]+2 > numseq) iend=numseq;
		else iend=mnodes[i]+2;
		coord2parambb(restr,mnodes[i],iend,1);
		coord2paramsugar(restr,mnodes[i],iend,1);
	}
	*/

	//cout<<"ppp"<<endl;

	int maxiter=1000;
	int iter;
	int estopiter=20;
	int nonupiter=0;
	double energy, enetmp;

	vector<vector<bool>> atomflag;
	vector<bool> residueflag;
	map<string, vector<double>> edgeenergy, edgeenergytmp;
	map<string, vector<double>> nodeenergy, nodeenergytmp;
	bool initflag=true;
	vector<Residueatoms> restrtmp;
	restrtmp.resize(restr.size());

	//cout<<"qqq"<<endl;

	resetflaggss(residueflag, atomflag, numseq);
	//cout<<"rrr"<<endl;
	initenergygss(edgeenergy, nodeenergy, edges, mnodes, numseq);
	//cout<<"sss"<<endl;
	energy=calcenegss(restr, decstr, numseq, edges, edgetype, paraminfo, nodes, dnttemplates, bptemplates, dnttemplateslength, bptemplateslength, enecoeff, residueflag, atomflag, edgeenergy, nodeenergy, 0);

	//cout<<"ttt"<<endl;

	for(iter=0;iter<maxiter;iter++){
                //cout<<"uuu"<<endl;
		wpos=mnodes[int(Random()*mnodes.size()-EPS)];
		//cout<<"vvv"<<endl;
		restrtmp=restr;
		edgeenergytmp=edgeenergy;
		nodeenergytmp=nodeenergy;

		//cout<<"www"<<endl;
		aat=restrtmp[wpos].aaa;
                
		//cout<<"xxx"<<endl;

		randval=Random();
		if(randval < 0.34){
			//cout<<"y1"<<endl;
			dihtype=TORSIONANGLETARGET[int(Random()*TORSIONANGLETARGET.size())];
			enetmp=goldenSearchdih(restrtmp, decstr, wpos, dihtype, 1, numseq, edges, edgetype, paraminfo, nodes, dnttemplates, bptemplates, dnttemplateslength, bptemplateslength, enecoeff, edgeenergytmp, nodeenergytmp);
		}else if((0.34 <= randval) & (randval < 0.67)){
			//cout<<"y2"<<endl;
			angtype=INNERANGLETARGET[int(Random()*INNERANGLETARGET.size())];
			enetmp=goldenSearchang(restrtmp, decstr, wpos, angtype, 1, numseq, edges, edgetype, paraminfo, nodes, dnttemplates, bptemplates, dnttemplateslength, bptemplateslength, enecoeff, edgeenergytmp, nodeenergytmp);
		}else{
			//cout<<"y3"<<endl;
			lentype=BLENGTHTARGET[int(Random()*BLENGTHTARGET.size())];
			enetmp=goldenSearchlen(restrtmp, decstr, wpos, lentype, 0.1, numseq, edges, edgetype, paraminfo, nodes, dnttemplates, bptemplates, dnttemplateslength, bptemplateslength, enecoeff, edgeenergytmp, nodeenergytmp);
		}

		//cout<<"zzz"<<endl;

		
		if(energy>enetmp){
			//cout<<"a0"<<endl;
			energy=enetmp;
			restr=restrtmp;
			edgeenergy=edgeenergytmp;
			nodeenergy=nodeenergytmp;
			nonupiter=0;
			//cout<<ene<<endl;
		}else{
			//cout<<"a1"<<endl;
			nonupiter+=1;
		}

		//cout<<"a2"<<endl;

		if(nonupiter == estopiter) break;
	}

	//cout<<"a3"<<endl;

	vector<bool> targetidxtmp(numseq, false);
	for(i=0;i<pairs.size();i++)
	{
		basei = pairs[i].first;
		basej = pairs[i].second;

		targetidxtmp[basei] = 1;
		if(basej != -1){
			targetidxtmp[basej] = 1;
		}
	}

	//cout<<"a4"<<endl;
	vector<pair<int,int>> continuousidx;
	int start = -1;

	for(i=0;i<numseq;i++){
		if (targetidxtmp[i] == 1) {
			if (start == -1) {
				start = i;
			}
		}else {
            if (start != -1) {
                continuousidx.push_back({start, i - 1});
                start = -1;
            }
        }
	}
	//cout<<"a5"<<endl;
	if (start != -1) {
        continuousidx.push_back({start, numseq - 1});
    }
	for(i = 0; i<continuousidx.size();i++){
		istart=continuousidx[i].first-1;
		iend=continuousidx[i].second+2;
		if(iend>=numseq) iend=numseq;
		//cout<<to_string(istart) + "-" + to_string(iend)<<endl;
		coord2parambb(restr,continuousidx[i].first,iend,1);
		coord2paramsugar(restr,continuousidx[i].first,iend,1);
		param2coordbase(restr,continuousidx[i].first,iend);
		param2coordhydrogen(restr,continuousidx[i].first,iend);

		coord2paramadd(restr,istart,iend);
		settmplmpflap(numseq, istart, iend);
	}

	//cout<<"a6"<<endl;

	newenergy=fcalcrmsdenergy(restr, numseq, lp2, rppmatbk, bpenergybk, bsenergybk, nmapbk, distmapbk, ntpairdistbk, adjpairdistbk, distenebk, multidistenebk, contactenebk, mapflagbk, mapbbflagbk, flagacc);

	//cout<<"a7"<<endl;
	//exit(1);
	
	/*
	cout<<"QqqqqqqqqqqqqqqqqQ"<<endl;
	checkupdate(restr, numseq);
	cout<<"PppppppppppppppppP"<<endl;
	fsetinitial(lp2,rppmatbk, bpenergybk, bsenergybk, nmapbk, distmapbk, ntpairdistbk, adjpairdistbk, distenebk, multidistenebk, contactenebk, mapflagbk, mapbbflagbk, numseq);
	double newenergy1=fcalcrmsdenergy(restr, numseq, lp2, rppmatbk, bpenergybk, bsenergybk, nmapbk, distmapbk, ntpairdistbk, adjpairdistbk, distenebk, multidistenebk, contactenebk, mapflagbk, mapbbflagbk, false);
	if(fabs(newenergy1-newenergy) > 10){
		cout<<"Error"<<endl;
		cout<<"2:: " + to_string(newenergy1)<<endl;
		cout<<"3:: " + to_string(newenergy)<<endl;
		exit(1);
	}
	cout<<"OK"<<endl;
	*/

	//exit(1);
	
	//cout<<"a8"<<endl;
	
	bool upidx=metropolisUpdate(oldenergy, newenergy, Bt);
	//bool upidx=metropolisUpdate(oldenergy, newenergy, 30*Random());

        //cout<<"a9"<<endl;

	if(upidx)
	{
		rcopystr(decstr, restr, numseq);
		return true;
	}
	else
	{
		return false;
	}
}












