
#include "MonteCarloSimulation.hpp"

void MonteCarloSimulation::setrotateflap(double trandpos1, double trandpos2)
{
	int i,j;
	int atom1, atom2, atom3;
	char base10,base11,base20,base21;

	base10=refinestr_temp[trandpos1].aaa;
	base11=refinestr_temp[trandpos1+1].aaa;
	base20=refinestr_temp[trandpos2].aaa;

	lengthflag[trandpos1][4] = 1;
	angleflag[trandpos1][4] = 1;
	torsionflag[trandpos1][4] = 1;
	lengthflag[trandpos1][5] = 1;
	angleflag[trandpos1][5] = 1;
	torsionflag[trandpos1][5] = 1;
	lengthflag[trandpos1+1][0] = 1;
	angleflag[trandpos1+1][0] = 1;
	torsionflag[trandpos1+1][0] = 1;
	angleflag[trandpos1+1][1] = 1;
	torsionflag[trandpos1+1][1] = 1;
	torsionflag[trandpos1+1][2] = 1;

	for(j=0;j<AddLengthParamsList[base10].size();j++){
		addlengthflag[trandpos1][j] = 1;
	}

	for(j=0;j<AddLengthParamsList[base11].size();j++){
		addlengthflag[trandpos1+1][j] = 1;
	}

	for(j=0;j<AddAngleParamsList[base10].size();j++){
		addangleflag[trandpos1][j] = 1;
	}

	for(j=0;j<AddAngleParamsList[base11].size();j++){
		addangleflag[trandpos1+1][j] = 1;
	}

	for(j=0;j<AddHydLengthParamsList[base10].size();j++)
	{
		hydaddlengthflag[trandpos1][j]=1;
	}

	for(j=0;j<AddHydLengthParamsList[base11].size();j++)
	{
		hydaddlengthflag[trandpos1+1][j]=1;
	}

	for(j=0;j<AddHydAngleParamsList[base10].size();j++)
	{
		hydaddangleflag[trandpos1][j]=1;
	}

	for(j=0;j<AddHydAngleParamsList[base11].size();j++)
	{
		hydaddangleflag[trandpos1+1][j]=1;
	}

	lengthflag[trandpos2][0] = 1;
	angleflag[trandpos2][0] = 1;
	torsionflag[trandpos2][0] = 1;
	lengthflag[trandpos2][1] = 1;
	angleflag[trandpos2][1] = 1;
	torsionflag[trandpos2][1] = 1;
	lengthflag[trandpos2][2] = 1;
	angleflag[trandpos2][2] = 1;
	torsionflag[trandpos2][2] = 1;
	lengthflag[trandpos2][3] = 1;
	angleflag[trandpos2][3] = 1;
	torsionflag[trandpos2][3] = 1;
	angleflag[trandpos2][4] = 1;
	torsionflag[trandpos2][4] = 1;
	torsionflag[trandpos2][5] = 1;

	for(j=0;j<AddLengthParamsList[base20].size();j++){
		addlengthflag[trandpos2][j] = 1;
	}

	if(trandpos2 < numbb-1){
		addlengthflag[trandpos2+1][0] = 1;
	}

	for(j=0;j<AddAngleParamsList[base20].size();j++){
		addangleflag[trandpos2][j] = 1;
	}

	for(j=0;j<AddHydLengthParamsList[base20].size();j++){
		hydaddlengthflag[trandpos2][j]=1;
	}

	for(j=0;j<AddHydAngleParamsList[base20].size();j++)
	{
		hydaddangleflag[trandpos2][j]=1;
	}

	for(i=0;i<numbb;i++)
	{
		ntflag[i] = 1;
		for(j=i;j<numbb;j++)
		{
			if(j<trandpos1 || i>trandpos2){}
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



bool MonteCarloSimulation::fmovementrot(double oldenergy, double &newenergy)
{
	int i, j, k;
	bool flagtor;
	int trandpos, trandposst, trandposed,span;
	double mean, std, upper, lower;
	int idx;
	double rmax, rnorm;
	Point3s v1;
	double mati[9];
	int mtype1, mtype2;
	char base, base1, base2;
	double trandang;
	bool flagdone;
	bool exitOuterLoop = false;

	vector<int> BBATOMIDX={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};

	selecttargetposLMP(trandpos, 0.8);
	base=refinestr[trandpos].aaa;

	if(Random()<0.5)//leftend
	{
		trandposed=int(5*Random())+1;
		trandposed+=trandpos;
	}else//rightend
	{
		trandposed=trandpos;
		trandpos=int(5*Random())+1;
		trandpos=trandposed-trandpos;
	}

	if(trandpos<2){
		trandpos=2;
		if(trandposed<=2){
			trandposed=trandpos+int(5*Random())+1;
		}
	}

	if(trandposed>numbb-2){
		trandposed=numbb-2;
		if(trandpos>=numbb-2){
			trandpos=trandposed-int(5*Random())+1;
			if(trandpos<2) trandpos=2;
		}
	}

	//cout<<"=============="<<endl;
	//cout<<trandpos << "-" << trandposed<<endl;

	StrCopy(refinestr, refinestr_temp);
	resetfrag(0);

	base1=refinestr_temp[trandpos-1].aaa;
	base2=refinestr_temp[trandposed+1].aaa;

	trandang=180.0*Random()-90.0;
	trandang=trandang*RADDEG;//rotangle

	v1.x=refinestr_temp[trandposed+1].ptv[5].x-refinestr_temp[trandpos-1].ptv[5].x;
	v1.y=refinestr_temp[trandposed+1].ptv[5].y-refinestr_temp[trandpos-1].ptv[5].y;
	v1.z=refinestr_temp[trandposed+1].ptv[5].z-refinestr_temp[trandpos-1].ptv[5].z;
	v1=bf.units(v1);//direction

	double vcos=cos(trandang);
	double vsin=sin(trandang);
	double vxx=v1.x*v1.x;
	double vxy=v1.x*v1.y;
	double vxz=v1.x*v1.z;
	double vyy=v1.y*v1.y;
	double vyz=v1.y*v1.z;
	double vzz=v1.z*v1.z;
	mati[0]=vxx+(1-vxx)*vcos;
	mati[1]=vxy*(1-vcos)-v1.z*vsin;
	mati[2]=vxz*(1-vcos)+v1.y*vsin;
	mati[3]=vxy*(1-vcos)+v1.z*vsin;
	mati[4]=vyy+(1-vyy)*vcos;
	mati[5]=vyz*(1-vcos)-v1.x*vsin;
	mati[6]=vxz*(1-vcos)-v1.y*vsin;
	mati[7]=vyz*(1-vcos)+v1.x*vsin;
	mati[8]=vzz+(1-vzz)*vcos;

	for(i=trandpos;i<=trandposed;i++)
	{
		for(j=0;j<BBATOMIDX.size();j++)
		{
			v1.x=(refinestr_temp[i].ptv[BBATOMIDX[j]].x-refinestr_temp[trandpos-1].ptv[5].x)*mati[0]
				+(refinestr_temp[i].ptv[BBATOMIDX[j]].y-refinestr_temp[trandpos-1].ptv[5].y)*mati[1]
				+(refinestr_temp[i].ptv[BBATOMIDX[j]].z-refinestr_temp[trandpos-1].ptv[5].z)*mati[2]
				+refinestr_temp[trandpos-1].ptv[5].x;
			v1.y=(refinestr_temp[i].ptv[BBATOMIDX[j]].x-refinestr_temp[trandpos-1].ptv[5].x)*mati[3]
				+(refinestr_temp[i].ptv[BBATOMIDX[j]].y-refinestr_temp[trandpos-1].ptv[5].y)*mati[4]
				+(refinestr_temp[i].ptv[BBATOMIDX[j]].z-refinestr_temp[trandpos-1].ptv[5].z)*mati[5]
				+refinestr_temp[trandpos-1].ptv[5].y;
			v1.z=(refinestr_temp[i].ptv[BBATOMIDX[j]].x-refinestr_temp[trandpos-1].ptv[5].x)*mati[6]
				+(refinestr_temp[i].ptv[BBATOMIDX[j]].y-refinestr_temp[trandpos-1].ptv[5].y)*mati[7]
				+(refinestr_temp[i].ptv[BBATOMIDX[j]].z-refinestr_temp[trandpos-1].ptv[5].z)*mati[8]
				+refinestr_temp[trandpos-1].ptv[5].z;
			refinestr_temp[i].ptv[BBATOMIDX[j]].x=v1.x;
			refinestr_temp[i].ptv[BBATOMIDX[j]].y=v1.y;
			refinestr_temp[i].ptv[BBATOMIDX[j]].z=v1.z;
		}
	}
	i=trandpos-1;
	for(j=4;j<BBATOMIDX.size();j++)
	{
		v1.x=(refinestr_temp[i].ptv[BBATOMIDX[j]].x-refinestr_temp[trandpos-1].ptv[5].x)*mati[0]
			+(refinestr_temp[i].ptv[BBATOMIDX[j]].y-refinestr_temp[trandpos-1].ptv[5].y)*mati[1]
			+(refinestr_temp[i].ptv[BBATOMIDX[j]].z-refinestr_temp[trandpos-1].ptv[5].z)*mati[2]
			+refinestr_temp[trandpos-1].ptv[5].x;
		v1.y=(refinestr_temp[i].ptv[BBATOMIDX[j]].x-refinestr_temp[trandpos-1].ptv[5].x)*mati[3]
			+(refinestr_temp[i].ptv[BBATOMIDX[j]].y-refinestr_temp[trandpos-1].ptv[5].y)*mati[4]
			+(refinestr_temp[i].ptv[BBATOMIDX[j]].z-refinestr_temp[trandpos-1].ptv[5].z)*mati[5]
			+refinestr_temp[trandpos-1].ptv[5].y;
		v1.z=(refinestr_temp[i].ptv[BBATOMIDX[j]].x-refinestr_temp[trandpos-1].ptv[5].x)*mati[6]
			+(refinestr_temp[i].ptv[BBATOMIDX[j]].y-refinestr_temp[trandpos-1].ptv[5].y)*mati[7]
			+(refinestr_temp[i].ptv[BBATOMIDX[j]].z-refinestr_temp[trandpos-1].ptv[5].z)*mati[8]
			+refinestr_temp[trandpos-1].ptv[5].z;
		refinestr_temp[i].ptv[BBATOMIDX[j]].x=v1.x;
		refinestr_temp[i].ptv[BBATOMIDX[j]].y=v1.y;
		refinestr_temp[i].ptv[BBATOMIDX[j]].z=v1.z;
	}
	i=trandposed+1;
	for(j=0;j<3;j++)
	{
		v1.x=(refinestr_temp[i].ptv[BBATOMIDX[j]].x-refinestr_temp[trandpos-1].ptv[5].x)*mati[0]
			+(refinestr_temp[i].ptv[BBATOMIDX[j]].y-refinestr_temp[trandpos-1].ptv[5].y)*mati[1]
			+(refinestr_temp[i].ptv[BBATOMIDX[j]].z-refinestr_temp[trandpos-1].ptv[5].z)*mati[2]
			+refinestr_temp[trandpos-1].ptv[5].x;
		v1.y=(refinestr_temp[i].ptv[BBATOMIDX[j]].x-refinestr_temp[trandpos-1].ptv[5].x)*mati[3]
			+(refinestr_temp[i].ptv[BBATOMIDX[j]].y-refinestr_temp[trandpos-1].ptv[5].y)*mati[4]
			+(refinestr_temp[i].ptv[BBATOMIDX[j]].z-refinestr_temp[trandpos-1].ptv[5].z)*mati[5]
			+refinestr_temp[trandpos-1].ptv[5].y;
		v1.z=(refinestr_temp[i].ptv[BBATOMIDX[j]].x-refinestr_temp[trandpos-1].ptv[5].x)*mati[6]
			+(refinestr_temp[i].ptv[BBATOMIDX[j]].y-refinestr_temp[trandpos-1].ptv[5].y)*mati[7]
			+(refinestr_temp[i].ptv[BBATOMIDX[j]].z-refinestr_temp[trandpos-1].ptv[5].z)*mati[8]
			+refinestr_temp[trandpos-1].ptv[5].z;
		refinestr_temp[i].ptv[BBATOMIDX[j]].x=v1.x;refinestr_temp[i].ptv[BBATOMIDX[j]].y=v1.y;refinestr_temp[i].ptv[BBATOMIDX[j]].z=v1.z;
	}
	int istart=trandpos-1;
	int iend=trandposed+2;
	if(iend>=numbb) iend=numbb;

	param2coordsugar(refinestr_temp,istart,iend);
	param2coordbase(refinestr_temp,istart,iend);
	param2coordhydrogen(refinestr_temp,istart,iend);

	coord2parambb(refinestr_temp,istart,iend);
	//coord2paramsugar(refinestr_temp,istart,iend);

	int end = min(iend+2, numbb);
	coord2paramadd(refinestr_temp,istart,end);
	
	coord2addparamhydrogen(refinestr_temp,iend-1,end);
	coord2addparamhydrogen(refinestr_temp,istart,trandpos+2);

	setrotateflap(trandpos-1, trandposed+1);
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

