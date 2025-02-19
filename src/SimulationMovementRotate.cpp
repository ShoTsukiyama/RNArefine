// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

void RNASimulation::settmprotateflap(double numseq, double trandpos1, double trandpos2, char aat1, char aat2)
{
	int i,j;
	int st, ed;

	lp2.flap[BLENGTH_ENERGY_ST+4][trandpos1]=true;
	lp2.flap[INNERANGLE_ENERGY_ST+4][trandpos1]=true;
	lp2.flap[TORSIONANGLE_ENERGY_ST+4][trandpos1]=true;
	lp2.flap[BLENGTH_ENERGY_ST+5][trandpos1]=true;
	lp2.flap[INNERANGLE_ENERGY_ST+5][trandpos1]=true;
	lp2.flap[TORSIONANGLE_ENERGY_ST+5][trandpos1]=true;
	lp2.flap[BLENGTH_ENERGY_ST+0][trandpos1+1]=true;
	lp2.flap[INNERANGLE_ENERGY_ST+0][trandpos1+1]=true;
	lp2.flap[TORSIONANGLE_ENERGY_ST+0][trandpos1+1]=true;
	lp2.flap[INNERANGLE_ENERGY_ST+1][trandpos1+1]=true;
	lp2.flap[TORSIONANGLE_ENERGY_ST+1][trandpos1+1]=true;
	lp2.flap[TORSIONANGLE_ENERGY_ST+2][trandpos1+1]=true;

	for(j=0;j<ADDBLENGTHTARGET.size();j++)
	{
		lp2.flap[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]][trandpos1]=true;
		lp2.flap[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]][trandpos1+1]=true;
	}
	if(trandpos1+1 < numseq-1){
		lp2.flap[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[0]][trandpos1+2]=true;
	}
	for(j=0;j<ADDINNERANGLETARGET.size();j++)
	{
		lp2.flap[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]][trandpos1]=true;
		lp2.flap[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]][trandpos1+1]=true;
	}

	lp2.flap[BLENGTH_ENERGY_ST+0][trandpos2]=true;
	lp2.flap[INNERANGLE_ENERGY_ST+0][trandpos2]=true;
	lp2.flap[TORSIONANGLE_ENERGY_ST+0][trandpos2]=true;
	lp2.flap[BLENGTH_ENERGY_ST+1][trandpos2]=true;
	lp2.flap[INNERANGLE_ENERGY_ST+1][trandpos2]=true;
	lp2.flap[TORSIONANGLE_ENERGY_ST+1][trandpos2]=true;
	lp2.flap[BLENGTH_ENERGY_ST+2][trandpos2]=true;
	lp2.flap[INNERANGLE_ENERGY_ST+2][trandpos2]=true;
	lp2.flap[TORSIONANGLE_ENERGY_ST+2][trandpos2]=true;
	lp2.flap[BLENGTH_ENERGY_ST+3][trandpos2]=true;
	lp2.flap[INNERANGLE_ENERGY_ST+3][trandpos2]=true;
	lp2.flap[TORSIONANGLE_ENERGY_ST+3][trandpos2]=true;
	lp2.flap[INNERANGLE_ENERGY_ST+4][trandpos2]=true;
	lp2.flap[TORSIONANGLE_ENERGY_ST+4][trandpos2]=true;
	lp2.flap[TORSIONANGLE_ENERGY_ST+5][trandpos2]=true;

	for(j=0;j<ADDBLENGTHTARGET.size();j++)
	{
		lp2.flap[ADDBLENGTH_ENERGY_ST+ADDBLENGTHTARGET[j]][trandpos2]=true;
	}
	if(trandpos2 < numseq-1){
		lp2.flap[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[0]][trandpos2+1]=true;
	}
	for(j=0;j<ADDINNERANGLETARGET.size();j++)
	{
		lp2.flap[ADDINNERANGLE_ENERGY_ST+ADDINNERANGLETARGET[j]][trandpos2]=true;
	}
	
	for(i=0;i<numseq;i++)
	{
		for(j=i;j<numseq;j++)
		//for(j=i+1;j<numseq;j++)
		{
			if(j<trandpos1 || i>trandpos2){}
			else
			{
				mapbbflagbk[i][j]=true;
				mapbbflagbk[j][i]=true;
				mapflagbk[i][j]=true;
				mapflagbk[j][i]=true;
			}
		}
	}
}


bool RNASimulation::fmovementrot(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc)
{
	int i,j;
	int istart,iend;
	int trandpos,trandpos2,tmp;
	double trandang;
	double mati[9];
	Point3d v1;
	int wpos;
	char aat1, aat2;
	bool idxtmp1, idxtmp2;
	vector<int> scpos(numseq,0);
	vector<int> accsc(numseq,0);
	vector<int> accscr(numseq,0);
	int transbase;
	int total;
	int totwrong=0;
	for(i=0;i<BASE_ENERGY_NUM;i++) totwrong+=lp1.nn[BASE_ENERGY_ST+i];

	if(numseq < 5){
	    newenergy=oldenergy;
		return false;
	}

	double randval=Random();
	bool breakfrag=false;
	if(totwrong>0 && randval<0.1+totwrong/BASEPAIRS.size()){
		wpos=int(Random()*totwrong);

		i=-1;
		for(trandpos=0;trandpos<BASEPAIRS.size();trandpos++)
		{
			for(j=0;j<BASE_ENERGY_NUM;j++){
				if(lp1.indn[BASE_ENERGY_ST+j][BASEPAIRS[trandpos].basepidx.first]==1) i++;
				if(i==wpos){
					breakfrag=true;
					break;
				}
			}
			if(breakfrag) break;
		}
		if(Random()<0.5)//leftend
		{
			trandpos2=int(19*Random())+1;
			trandpos2+=trandpos;
		}
		else//rightend
		{
			trandpos2=trandpos;
			trandpos=int(19*Random())+1;
			trandpos=trandpos2-trandpos;
		}
		if(trandpos<2){
			trandpos=2;
			if(trandpos2<=2){
				trandpos2=trandpos+int(19*Random())+1;
				if(trandpos2>numseq-2) trandpos2=numseq-2;
			}
		}
		if(trandpos2>numseq-2){
			trandpos2=numseq-2;
			if(trandpos>=numseq-2){
				trandpos=trandpos2-int(19*Random())+1;
				if(trandpos<2) trandpos=2;
			}
		}
	}else{
		trandpos=int(numseq*Random());
		if(Random()<0.5)//leftend
		{
			trandpos2=int(19*Random())+1;
			trandpos2+=trandpos;
		}
		else//rightend
		{
			trandpos2=trandpos;
			trandpos=int(19*Random())+1;
			trandpos=trandpos2-trandpos;
		}
		if(trandpos<2){
			trandpos=2;
			if(trandpos2<=2){
				trandpos2=trandpos+int(19*Random())+1;
				if(trandpos2>numseq-2) trandpos2=numseq-2;
			}
		}
		if(trandpos2>numseq-2){
			trandpos2=numseq-2;
			if(trandpos>=numseq-2){
				trandpos=trandpos2-int(19*Random())+1;
				if(trandpos<2) trandpos=2;
			}
		}
		if(trandpos < 2) trandpos=2;
		if(trandpos2>numseq-2) trandpos2=numseq-2;
	}

	copystr(decstr, numseq);

	aat1=restr[trandpos-1].aaa;
	aat2=restr[trandpos2+1].aaa;

	trandang=180.0*Random()-90.0;
	trandang=trandang*RADDEG;//rotangle
	v1.x=restr[trandpos2+1].ptv[5].x-restr[trandpos-1].ptv[5].x;
	v1.y=restr[trandpos2+1].ptv[5].y-restr[trandpos-1].ptv[5].y;
	v1.z=restr[trandpos2+1].ptv[5].z-restr[trandpos-1].ptv[5].z;
	v1=bf.unit(v1);//direction

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

	for(i=trandpos;i<=trandpos2;i++)
	{
		for(j=0;j<BBATOMIDX.size();j++)
		{
			v1.x=(restr[i].ptv[BBATOMIDX[j]].x-restr[trandpos-1].ptv[5].x)*mati[0]
				+(restr[i].ptv[BBATOMIDX[j]].y-restr[trandpos-1].ptv[5].y)*mati[1]
				+(restr[i].ptv[BBATOMIDX[j]].z-restr[trandpos-1].ptv[5].z)*mati[2]
				+restr[trandpos-1].ptv[5].x;
			v1.y=(restr[i].ptv[BBATOMIDX[j]].x-restr[trandpos-1].ptv[5].x)*mati[3]
				+(restr[i].ptv[BBATOMIDX[j]].y-restr[trandpos-1].ptv[5].y)*mati[4]
				+(restr[i].ptv[BBATOMIDX[j]].z-restr[trandpos-1].ptv[5].z)*mati[5]
				+restr[trandpos-1].ptv[5].y;
			v1.z=(restr[i].ptv[BBATOMIDX[j]].x-restr[trandpos-1].ptv[5].x)*mati[6]
				+(restr[i].ptv[BBATOMIDX[j]].y-restr[trandpos-1].ptv[5].y)*mati[7]
				+(restr[i].ptv[BBATOMIDX[j]].z-restr[trandpos-1].ptv[5].z)*mati[8]
				+restr[trandpos-1].ptv[5].z;
			restr[i].ptv[BBATOMIDX[j]].x=v1.x;
			restr[i].ptv[BBATOMIDX[j]].y=v1.y;
			restr[i].ptv[BBATOMIDX[j]].z=v1.z;
		}
	}
	i=trandpos-1;
	for(j=4;j<BBATOMIDX.size();j++)
	{
		v1.x=(restr[i].ptv[BBATOMIDX[j]].x-restr[trandpos-1].ptv[5].x)*mati[0]
			+(restr[i].ptv[BBATOMIDX[j]].y-restr[trandpos-1].ptv[5].y)*mati[1]
			+(restr[i].ptv[BBATOMIDX[j]].z-restr[trandpos-1].ptv[5].z)*mati[2]
			+restr[trandpos-1].ptv[5].x;
		v1.y=(restr[i].ptv[BBATOMIDX[j]].x-restr[trandpos-1].ptv[5].x)*mati[3]
			+(restr[i].ptv[BBATOMIDX[j]].y-restr[trandpos-1].ptv[5].y)*mati[4]
			+(restr[i].ptv[BBATOMIDX[j]].z-restr[trandpos-1].ptv[5].z)*mati[5]
			+restr[trandpos-1].ptv[5].y;
		v1.z=(restr[i].ptv[BBATOMIDX[j]].x-restr[trandpos-1].ptv[5].x)*mati[6]
			+(restr[i].ptv[BBATOMIDX[j]].y-restr[trandpos-1].ptv[5].y)*mati[7]
			+(restr[i].ptv[BBATOMIDX[j]].z-restr[trandpos-1].ptv[5].z)*mati[8]
			+restr[trandpos-1].ptv[5].z;
		restr[i].ptv[BBATOMIDX[j]].x=v1.x;
		restr[i].ptv[BBATOMIDX[j]].y=v1.y;
		restr[i].ptv[BBATOMIDX[j]].z=v1.z;
	}
	i=trandpos2+1;
	for(j=0;j<3;j++)
	{
		v1.x=(restr[i].ptv[BBATOMIDX[j]].x-restr[trandpos-1].ptv[5].x)*mati[0]
			+(restr[i].ptv[BBATOMIDX[j]].y-restr[trandpos-1].ptv[5].y)*mati[1]
			+(restr[i].ptv[BBATOMIDX[j]].z-restr[trandpos-1].ptv[5].z)*mati[2]
			+restr[trandpos-1].ptv[5].x;
		v1.y=(restr[i].ptv[BBATOMIDX[j]].x-restr[trandpos-1].ptv[5].x)*mati[3]
			+(restr[i].ptv[BBATOMIDX[j]].y-restr[trandpos-1].ptv[5].y)*mati[4]
			+(restr[i].ptv[BBATOMIDX[j]].z-restr[trandpos-1].ptv[5].z)*mati[5]
			+restr[trandpos-1].ptv[5].y;
		v1.z=(restr[i].ptv[BBATOMIDX[j]].x-restr[trandpos-1].ptv[5].x)*mati[6]
			+(restr[i].ptv[BBATOMIDX[j]].y-restr[trandpos-1].ptv[5].y)*mati[7]
			+(restr[i].ptv[BBATOMIDX[j]].z-restr[trandpos-1].ptv[5].z)*mati[8]
			+restr[trandpos-1].ptv[5].z;
		restr[i].ptv[BBATOMIDX[j]].x=v1.x;restr[i].ptv[BBATOMIDX[j]].y=v1.y;restr[i].ptv[BBATOMIDX[j]].z=v1.z;
	}
	istart=trandpos-1;
	iend=trandpos2+2;
	if(iend>=numseq) iend=numseq;

	param2coordsugar(restr,istart,iend);
	param2coordbase(restr,istart,iend);
	param2coordhydrogen(restr,istart,numseq);

	coord2parambb(restr,istart,iend,1);

	int end = min(iend+2, numseq);
	coord2paramadd(restr,istart, end);

	//coord2paramsugar(restr,istart,iend,1);

	if(fabs(restr[trandpos-1].ang[4]-INNERANGPARAM[aat1][4].mean)<=10.0) idxtmp1=true;
	else idxtmp1=false;
	if(fabs(restr[trandpos2+1].ang[4]-INNERANGPARAM[aat2][4].mean)<=10.0) idxtmp2=true;
	else idxtmp2=false;

	if(idxtmp1 && idxtmp2){

		settmprotateflap(numseq, trandpos-1, trandpos2+1, aat1, aat2);
		newenergy=fcalcrmsdenergy(restr, numseq, lp2, rppmatbk, bpenergybk, bsenergybk, nmapbk, distmapbk, ntpairdistbk, adjpairdistbk, distenebk, multidistenebk, contactenebk, mapflagbk, mapbbflagbk, flagacc);

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
		
		bool upidx=metropolisUpdate(oldenergy, newenergy, Bt);

		if(upidx)
		{
			rcopystr(decstr, restr, numseq);
			return true;
		}
		else
		{
			return false;
		}
	}else{
		newenergy=oldenergy;
		return false;
	}
}












