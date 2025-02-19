// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

void RNASimulation::setlenflaplocal(double numseq, int trandpos, double aat, int lentype)
{
	int i,j;

	lp2.flap[BLENGTH_ENERGY_ST+lentype][trandpos]=true;

	for(i=0;i<numseq;i++)
	{
		for(j=i;j<numseq;j++)
		//for(j=i+1;j<numseq;j++)
		{
			if(j<trandpos || i>trandpos+1){}
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


void RNASimulation::coord2paramlocal(vector<Residueatoms> &decstr, int trandpos, int moveatom, int start, int end)
{
	int i, j;
	char aat;
	char indnt;
	Point3d p12,p23;
	int atom1, atom2, atom3, atom4;
	Point3s tmp1,tmp2,tmp3,tmp4;
	BasicFunc bf;
	int realstart;

	if((start == 0) & (trandpos == 0))
	{
		indnt=decstr[0].aaa;
		aat=decstr[0].aaa;
		decstr[0].dih[0]=TORSIONPARAM[aat][0].mean;
		decstr[0].dih[1]=TORSIONPARAM[aat][1].mean;
		decstr[0].dih[2]=TORSIONPARAM[aat][2].mean;
		decstr[0].ang[0]=INNERANGPARAM[aat][0].mean;
		decstr[0].ang[1]=INNERANGPARAM[aat][1].mean;
		decstr[0].len[0]=LENGTHPARAM[aat][0].mean;
	
		atom3=ParamsList[aat][1][2];atom4=ParamsList[aat][1][3];
		tmp3=decstr[0].ptv[atom3];tmp4=decstr[0].ptv[atom4];
		if((trandpos == 0) && ((moveatom == atom3) || (moveatom == atom4))){
			coord2len(tmp3, tmp4, decstr[0].len[1]);

			lp2.flap[BLENGTH_ENERGY_ST+1][0]=true;
		}
		
		atom2=ParamsList[aat][2][1];atom3=ParamsList[aat][2][2];atom4=ParamsList[aat][2][3];
		tmp2=decstr[0].ptv[atom2];tmp3=decstr[0].ptv[atom3];tmp4=decstr[0].ptv[atom4];

		if((trandpos == 0) && ((moveatom == atom2) || (moveatom == atom3) || (moveatom == atom4))){
			coord2anglen(tmp2, tmp3, tmp4, decstr[0].len[2], decstr[0].ang[2]);

			lp2.flap[BLENGTH_ENERGY_ST+2][0]=true;
			lp2.flap[INNERANGLE_ENERGY_ST+2][0]=true;
		}
		
		for(j=3;j<AccAtomParamNum[aat][0];j++){
			atom1=ParamsList[aat][j][0];atom2=ParamsList[aat][j][1];atom3=ParamsList[aat][j][2];atom4=ParamsList[aat][j][3];
			tmp1=decstr[0].ptv[atom1];tmp2=decstr[0].ptv[atom2];tmp3=decstr[0].ptv[atom3];tmp4=decstr[0].ptv[atom4];

			if((trandpos == 0) && ((moveatom == atom1) || (moveatom == atom2) || (moveatom == atom3) || (moveatom == atom4))){
				coord2params(tmp1, tmp2, tmp3, tmp4, decstr[0].len[j], decstr[0].ang[j], decstr[0].dih[j]);

				lp2.flap[BLENGTH_ENERGY_ST+j][0]=true;
				lp2.flap[INNERANGLE_ENERGY_ST+j][0]=true;
				lp2.flap[TORSIONANGLE_ENERGY_ST+j][0]=true;
			}
		}
		realstart=1;
	}else{
		realstart=start;
	}

	for(i=realstart;i<end;i++)
	{
		aat=decstr[i].aaa;
		atom1=ParamsList[aat][0][0];atom2=ParamsList[aat][0][1];atom3=ParamsList[aat][0][2];atom4=ParamsList[aat][0][3];
		tmp1=decstr[i-1].ptv[atom1];tmp2=decstr[i-1].ptv[atom2];tmp3=decstr[i-1].ptv[atom3];tmp4=decstr[i].ptv[atom4];

		if(((trandpos == i-1) && ((moveatom == atom1) || (moveatom == atom2) || (moveatom == atom3))) || ((trandpos == i) && (moveatom == atom4)))
		{
			coord2params(tmp1, tmp2, tmp3, tmp4, decstr[i].len[0], decstr[i].ang[0], decstr[i].dih[0]);

			lp2.flap[BLENGTH_ENERGY_ST][i]=true;
			lp2.flap[INNERANGLE_ENERGY_ST][i]=true;
			lp2.flap[TORSIONANGLE_ENERGY_ST][i]=true;
		}
		
		atom1=ParamsList[aat][1][0];atom2=ParamsList[aat][1][1];atom3=ParamsList[aat][1][2];atom4=ParamsList[aat][1][3];
		tmp1=decstr[i-1].ptv[atom1];tmp2=decstr[i-1].ptv[atom2];tmp3=decstr[i].ptv[atom3];tmp4=decstr[i].ptv[atom4];

		if(((trandpos == i-1) && ((moveatom == atom1) || (moveatom == atom2))) || ((trandpos == i) && ((moveatom == atom3) || (moveatom == atom4))))
		{
			coord2params(tmp1, tmp2, tmp3, tmp4, decstr[i].len[1], decstr[i].ang[1], decstr[i].dih[1]);

			lp2.flap[BLENGTH_ENERGY_ST+1][i]=true;
			lp2.flap[INNERANGLE_ENERGY_ST+1][i]=true;
			lp2.flap[TORSIONANGLE_ENERGY_ST+1][i]=true;
		}

		atom1=ParamsList[aat][2][0];atom2=ParamsList[aat][2][1];atom3=ParamsList[aat][2][2];atom4=ParamsList[aat][2][3];
		tmp1=decstr[i-1].ptv[atom1];tmp2=decstr[i].ptv[atom2];tmp3=decstr[i].ptv[atom3];tmp4=decstr[i].ptv[atom4];
		

		if(((trandpos == i-1) && (moveatom == atom1)) || ((trandpos == i) && ((moveatom == atom2) || (moveatom == atom3) || (moveatom == atom4))))
		{
			coord2params(tmp1, tmp2, tmp3, tmp4, decstr[i].len[2], decstr[i].ang[2], decstr[i].dih[2]);

			lp2.flap[BLENGTH_ENERGY_ST+2][i]=true;
			lp2.flap[INNERANGLE_ENERGY_ST+2][i]=true;
			lp2.flap[TORSIONANGLE_ENERGY_ST+2][i]=true;
		}
		
		for(j=3;j<AccAtomParamNum[aat][0];j++){

			atom1=ParamsList[aat][j][0];atom2=ParamsList[aat][j][1];atom3=ParamsList[aat][j][2];atom4=ParamsList[aat][j][3];
			tmp1=decstr[i].ptv[atom1];tmp2=decstr[i].ptv[atom2];tmp3=decstr[i].ptv[atom3];tmp4=decstr[i].ptv[atom4];

			if(((trandpos == i) && ((moveatom == atom1) || (moveatom == atom2) || (moveatom == atom3) || (moveatom == atom4))))
			{
				coord2params(tmp1, tmp2, tmp3, tmp4, decstr[i].len[j], decstr[i].ang[j], decstr[i].dih[j]);

				lp2.flap[BLENGTH_ENERGY_ST+j][i]=true;
				lp2.flap[INNERANGLE_ENERGY_ST+j][i]=true;
				lp2.flap[TORSIONANGLE_ENERGY_ST+j][i]=true;
			}
		}
	}
}

void RNASimulation::coord2paramsugarlocal(vector<Residueatoms> &decstr, int trandpos, int moveatom, int start, int end)
{
	int i,j;
	char aat;
	int atom1, atom2, atom3, atom4;
	int realstart;
	BasicFunc bf;
	Point3s tmp1,tmp2,tmp3,tmp4;
	float len, ang, dih;
	
	if(start == 0 && trandpos == 0){
		aat=decstr[0].aaa;
		for(j=AccAtomParamNum[aat][0];j<AccAtomParamNum[aat][2];j++)
		{
			if(j == 10 || j == 11){
				atom1=FIRSTOPParamsList[aat][j][0]; atom2=FIRSTOPParamsList[aat][j][1]; atom3=FIRSTOPParamsList[aat][j][2]; atom4=FIRSTOPParamsList[aat][j][3];
			}else{
				atom1=ParamsList[aat][j][0]; atom2=ParamsList[aat][j][1]; atom3=ParamsList[aat][j][2]; atom4=ParamsList[aat][j][3];
			}

			if(atom1 == moveatom | atom2 == moveatom | atom3 == moveatom | atom3 == moveatom){
				tmp1=decstr[0].ptv[atom1]; tmp2=decstr[0].ptv[atom2]; tmp3=decstr[0].ptv[atom3]; tmp4=decstr[0].ptv[atom4];
				coord2params(tmp1, tmp2, tmp3, tmp4, len, ang, dih);
				decstr[0].dih[j]=dih; decstr[0].len[j]=len; decstr[0].ang[j]=ang;

				lp2.flap[BLENGTH_ENERGY_ST+j][0]=true;
				lp2.flap[INNERANGLE_ENERGY_ST+j][0]=true;
				lp2.flap[TORSIONANGLE_ENERGY_ST+j][0]=true;
			}
		}
		realstart=1;
	}else{
		realstart=start;
	}
	for(i=realstart;i<end;i++)
	{
		aat=decstr[i].aaa;
		for(j=AccAtomParamNum[aat][0];j<AccAtomParamNum[aat][2];j++){
			atom1=ParamsList[aat][j][0]; atom2=ParamsList[aat][j][1]; atom3=ParamsList[aat][j][2]; atom4=ParamsList[aat][j][3];
			if(j == 10){
				tmp1=decstr[i-1].ptv[atom1]; tmp2=decstr[i].ptv[atom2]; tmp3=decstr[i].ptv[atom3]; tmp4=decstr[i].ptv[atom4];

				if(((trandpos == i-1) && (moveatom == atom1)) || ((trandpos == i) && ((moveatom == atom2) || (moveatom == atom3) || (moveatom == atom4))))
				{
					coord2params(tmp1, tmp2, tmp3, tmp4, decstr[i].len[j], decstr[i].ang[j], decstr[i].dih[j]);

					lp2.flap[BLENGTH_ENERGY_ST+j][i]=true;
					lp2.flap[INNERANGLE_ENERGY_ST+j][i]=true;
					lp2.flap[TORSIONANGLE_ENERGY_ST+j][i]=true;
				}

			}else if(j == 11){
				tmp1=decstr[i].ptv[atom1]; tmp2=decstr[i-1].ptv[atom2]; tmp3=decstr[i].ptv[atom3]; tmp4=decstr[i].ptv[atom4];

				if(((trandpos == i-1) && (moveatom == atom2)) || ((trandpos == i) && ((moveatom == atom1) || (moveatom == atom3) || (moveatom == atom4))))
				{
					coord2params(tmp1, tmp2, tmp3, tmp4, decstr[i].len[j], decstr[i].ang[j], decstr[i].dih[j]);

					lp2.flap[BLENGTH_ENERGY_ST+j][i]=true;
					lp2.flap[INNERANGLE_ENERGY_ST+j][i]=true;
					lp2.flap[TORSIONANGLE_ENERGY_ST+j][i]=true;
				}

			}else{
				tmp1=decstr[i].ptv[atom1]; tmp2=decstr[i].ptv[atom2]; tmp3=decstr[i].ptv[atom3]; tmp4=decstr[i].ptv[atom4];

				if(((trandpos == i) && ((moveatom == atom1) || (moveatom == atom2) || (moveatom == atom3) || (moveatom == atom4))))
				{
					coord2params(tmp1, tmp2, tmp3, tmp4, decstr[i].len[j], decstr[i].ang[j], decstr[i].dih[j]);

					lp2.flap[BLENGTH_ENERGY_ST+j][i]=true;
					lp2.flap[INNERANGLE_ENERGY_ST+j][i]=true;
					lp2.flap[TORSIONANGLE_ENERGY_ST+j][i]=true;
				}
			}
		}
	}
}

void RNASimulation::coord2paramaddlocal(vector<Residueatoms> &decstr, int trandpos, int moveatom,int start,int end)
{
	int i,j;
	char aat;
	int atom1, atom2, atom3, atom4;
	int realstart;
	int batomnum;
	int aaind[4];
	BasicFunc bf;
	Point3s tmp1,tmp2,tmp3,tmp4;
	float len, ang, dih;
	
	if(start == 0 && trandpos == 0){
		aat=decstr[0].aaa;
		for(j=1;j<AddAngleParamsList[aat].size();j++)
		{
			atom1=AddAngleParamsList[aat][j][0]; atom2=AddAngleParamsList[aat][j][1]; atom3=AddAngleParamsList[aat][j][2];

			if((atom1 == moveatom) | (atom2 == moveatom) | (atom3 == moveatom))
			{
				tmp1=decstr[0].ptv[atom1]; tmp2=decstr[0].ptv[atom2]; tmp3=decstr[0].ptv[atom3];
				coord2ang(tmp1, tmp2, tmp3, ang);
				decstr[0].addang[j]=ang;

				lp2.flap[ADDINNERANGLE_ENERGY_ST+j][0]=true;
			}
		}

		for(j=0;j<AddLengthParamsList[aat].size();j++)
		{
			atom1=AddLengthParamsList[aat][j][0]; atom2=AddLengthParamsList[aat][j][1];

			if((atom1 == moveatom) | (atom2 == moveatom))
			{
				tmp1=decstr[0].ptv[atom1]; tmp2=decstr[0].ptv[atom2];
				coord2len(tmp1, tmp2, len);
				decstr[0].addlen[j]=len;

				lp2.flap[ADDBLENGTH_ENERGY_ST+j][0]=true;
			}
		}
		realstart=1;
	}else{
		realstart=start;
	}
	for(i=realstart;i<end;i++)
	{
		aat=decstr[i].aaa;

		//cout<<aat<<endl;
		//cout<<to_string(i) + "-" + to_string(0) + "-" + to_string(decstr[i].addang[0])<<endl;

		atom1=AddAngleParamsList[aat][0][0]; atom2=AddAngleParamsList[aat][0][1]; atom3=AddAngleParamsList[aat][0][2];

		if(((trandpos == i-1) && (moveatom == atom1)) || ((trandpos == i) && ((moveatom == atom2) || (moveatom == atom3))))
		{
			tmp1=decstr[i-1].ptv[atom1]; tmp2=decstr[i].ptv[atom2]; tmp3=decstr[i].ptv[atom3];
			coord2ang(tmp1, tmp2, tmp3, ang);
			decstr[i].addang[0]=ang;

			lp2.flap[ADDINNERANGLE_ENERGY_ST][i]=true;
		}

		//cout<<to_string(i) + "-" + to_string(0) + "-" + to_string(ang)<<endl;

		for(j=1;j<AddAngleParamsList[aat].size();j++)
		{
			//cout<<aat<<endl;

			atom1=AddAngleParamsList[aat][j][0]; atom2=AddAngleParamsList[aat][j][1]; atom3=AddAngleParamsList[aat][j][2];
			if(((trandpos == i) && ((moveatom == atom1) || (moveatom == atom2) || (moveatom == atom3))))
			{
				tmp1=decstr[i].ptv[atom1]; tmp2=decstr[i].ptv[atom2]; tmp3=decstr[i].ptv[atom3];
				coord2ang(tmp1, tmp2, tmp3, ang);
				decstr[i].addang[j]=ang;

				lp2.flap[ADDINNERANGLE_ENERGY_ST+j][i]=true;
			}
			
			//cout<<to_string(i) + "-" + to_string(j) + "-" + to_string(ang)<<endl;
		}

		for(j=0;j<AddLengthParamsList[aat].size();j++)
		{
			atom1=AddLengthParamsList[aat][j][0]; atom2=AddLengthParamsList[aat][j][1];

			if(((trandpos == i) && ((moveatom == atom1) || (moveatom == atom2))))
			{
				tmp1=decstr[i].ptv[atom1]; tmp2=decstr[i].ptv[atom2];
				coord2len(tmp1, tmp2, len);
				decstr[i].addlen[j]=len;

				lp2.flap[ADDBLENGTH_ENERGY_ST+j][i]=true;
			}
		}
	}
}


void RNASimulation::changelocallen(vector<Residueatoms> &decstr, int trandpos, int lentype, int numseq)
{
	char aat=decstr[trandpos].aaa;
	Point3s pn;
	bool flagts;
	int atom1, atom2, atom3, atom4;
	int moveatom;

	if(trandpos == 0){
		if(lentype == 10 || lentype == 11){
			atom1=FIRSTOPParamsList[aat][lentype][0];atom2=FIRSTOPParamsList[aat][lentype][1];atom3=FIRSTOPParamsList[aat][lentype][2];moveatom=FIRSTOPParamsList[aat][lentype][3];
		}else{
			atom1=ParamsList[aat][lentype][0];atom2=ParamsList[aat][lentype][1];atom3=ParamsList[aat][lentype][2];moveatom=ParamsList[aat][lentype][3];
		}

		flagts=bf.tor2pos22(decstr[trandpos].ptv[atom1].x,decstr[trandpos].ptv[atom1].y,decstr[trandpos].ptv[atom1].z,
							decstr[trandpos].ptv[atom2].x,decstr[trandpos].ptv[atom2].y,decstr[trandpos].ptv[atom2].z,
							decstr[trandpos].ptv[atom3].x,decstr[trandpos].ptv[atom3].y,decstr[trandpos].ptv[atom3].z,
							decstr[trandpos].dih[lentype]*RADDEG,decstr[trandpos].len[lentype],decstr[trandpos].ang[lentype]*RADDEG,
							&pn.x,&pn.y,&pn.z);
	}else{
		atom1=ParamsList[aat][lentype][0];atom2=ParamsList[aat][lentype][1];atom3=ParamsList[aat][lentype][2];moveatom=ParamsList[aat][lentype][3];
		
		flagts=bf.tor2pos22(decstr[trandpos].ptv[atom1].x,decstr[trandpos].ptv[atom1].y,decstr[trandpos].ptv[atom1].z,
							decstr[trandpos].ptv[atom2].x,decstr[trandpos].ptv[atom2].y,decstr[trandpos].ptv[atom2].z,
							decstr[trandpos].ptv[atom3].x,decstr[trandpos].ptv[atom3].y,decstr[trandpos].ptv[atom3].z,
							decstr[trandpos].dih[lentype]*RADDEG,decstr[trandpos].len[lentype],decstr[trandpos].ang[lentype]*RADDEG,
							&pn.x,&pn.y,&pn.z);
	}

	decstr[trandpos].ptv[moveatom].x=pn.x;
	decstr[trandpos].ptv[moveatom].y=pn.y;
	decstr[trandpos].ptv[moveatom].z=pn.z;

	int iend=trandpos+2;
	if(iend > numseq) iend=numseq;
	coord2paramlocal(decstr, trandpos, moveatom, trandpos, iend);
	coord2paramsugarlocal(decstr, trandpos, moveatom, trandpos, iend);
	coord2paramaddlocal(decstr, trandpos, moveatom, trandpos, iend);
}

bool RNASimulation::fmovementlenlocal(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc)//length 
{

	int i,j;
	double r,l;
	double nr,nl;
	int trandpos;
	bool flagtor;
	int wpos;
	char aat;
	double pval;
	bool upidx;

	int lentype=BLENGTHTARGET[int(Random()*BLENGTHTARGET.size())];
	int totwrong=lp1.nn[BLENGTH_ENERGY_ST+lentype];
	int wrong=0;
	double randval=Random();
	int uptype;
	//double Btang;

	//cout<<to_string(0) + "-" + to_string(basepairs[0].basepidx.first) + "-" + to_string(basepairs[0].basepidx.second)<<endl;


	if(totwrong>0){
		uptype=0;
		wpos=int(Random()*totwrong);

		i=-1;
		for(trandpos=0;trandpos<numseq;trandpos++)
		{
			if(lp1.indn[BLENGTH_ENERGY_ST+lentype][trandpos]==1) i++;
			if(i==wpos) break;
		}
	}
	else{
		newenergy=oldenergy;
        return false;
	}

	copystr(decstr, numseq);
	aat=restr[trandpos].aaa;

	if(trandpos==0 && lentype == 10){
		if(uptype == 0){
			restr[trandpos].len[lentype]=getmovlenval(restr[trandpos].len[lentype], OPLENGTHPARAM[aat][0].mean, OPLENGTHPARAM[aat][0].std, 0.1);
			//Btang=(1.4+0.02*numseq);
		}else{
			restr[trandpos].len[lentype]=getrandomnewparam(OPLENGTHPARAM[aat][0].lower, OPLENGTHPARAM[aat][0].upper);
			//Btang=-1;
		}
	}else if(trandpos==0 && lentype == 11){
		if(uptype == 0){
			restr[trandpos].len[lentype]=getmovlenval(restr[trandpos].len[lentype], OPLENGTHPARAM[aat][1].mean, OPLENGTHPARAM[aat][1].std, 0.1);
			//Btang=(1.4+0.02*numseq);
		}else{
			restr[trandpos].len[lentype]=getrandomnewparam(OPLENGTHPARAM[aat][1].lower, OPLENGTHPARAM[aat][1].upper);
			//Btang=-1;
		}
	}else{
		if(uptype == 0){
			restr[trandpos].len[lentype]=getmovlenval(restr[trandpos].len[lentype], LENGTHPARAM[aat][lentype].mean, LENGTHPARAM[aat][lentype].std, 0.1);
			//Btang=(1.4+0.02*numseq);
		}else{
			restr[trandpos].len[lentype]=getrandomnewparam(LENGTHPARAM[aat][lentype].lower, LENGTHPARAM[aat][lentype].upper);
			//Btang=-1;
		}
	}

	changelocallen(restr, trandpos, lentype, numseq);

	param2coordbase(restr,trandpos,numseq);
	param2coordhydrogen(restr,trandpos,numseq);

	setlenflaplocal(numseq, trandpos, aat, lentype);

	newenergy=fcalcrmsdenergy(restr, numseq, lp2, rppmatbk, bpenergybk, bsenergybk, nmapbk, distmapbk, ntpairdistbk, adjpairdistbk, distenebk, multidistenebk, contactenebk, mapflagbk, mapbbflagbk, flagacc);
	upidx=metropolisUpdate(oldenergy, newenergy, 20);
	//upidx=metropolisUpdate(oldenergy, newenergy, Btang);

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

	if(upidx){
		rcopystr(decstr, restr, numseq);
		return true;
	}
	else
	{
		return false;
	}
}












