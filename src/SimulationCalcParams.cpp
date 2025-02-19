// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


void RNASimulation::coord2params(Point3s tmp1, Point3s tmp2, Point3s tmp3, Point3s tmp4, float &len, float &ang, float &dih)
{
	BasicFunc bf;
	Point3d p12,p23;
	p12=bf.setv(tmp2.x-tmp3.x,tmp2.y-tmp3.y,tmp2.z-tmp3.z);
	p23=bf.setv(tmp4.x-tmp3.x,tmp4.y-tmp3.y,tmp4.z-tmp3.z);
	len=bf.norm(p23);
	ang=bf.angv(p12,p23);
	dih=bf.phi(tmp1.x,tmp1.y,tmp1.z,
			   tmp2.x,tmp2.y,tmp2.z,
			   tmp3.x,tmp3.y,tmp3.z,
			   tmp4.x,tmp4.y,tmp4.z);
}

void RNASimulation::coord2dangle(Point3s tmp1, Point3s tmp2, Point3s tmp3, Point3s tmp4, float &dih)
{
	BasicFunc bf;
	Point3d p12,p23;
	dih=bf.phi(tmp1.x,tmp1.y,tmp1.z,
			   tmp2.x,tmp2.y,tmp2.z,
			   tmp3.x,tmp3.y,tmp3.z,
			   tmp4.x,tmp4.y,tmp4.z);
}

void RNASimulation::coord2anglen(Point3s tmp2, Point3s tmp3, Point3s tmp4, float &len, float &ang)
{
	BasicFunc bf;
	Point3d p12,p23;
	p12=bf.setv(tmp2.x-tmp3.x,tmp2.y-tmp3.y,tmp2.z-tmp3.z);
	p23=bf.setv(tmp4.x-tmp3.x,tmp4.y-tmp3.y,tmp4.z-tmp3.z);
	len=bf.norm(p23);
	ang=bf.angv(p12,p23);
}

void RNASimulation::coord2ang(Point3s tmp2, Point3s tmp3, Point3s tmp4, float &ang)
{
	BasicFunc bf;
	Point3d p12,p23;
	p12=bf.setv(tmp2.x-tmp3.x,tmp2.y-tmp3.y,tmp2.z-tmp3.z);
	p23=bf.setv(tmp4.x-tmp3.x,tmp4.y-tmp3.y,tmp4.z-tmp3.z);
	ang=bf.angv(p12,p23);
}

void RNASimulation::coord2len(Point3s tmp3, Point3s tmp4, float &len)
{
	BasicFunc bf;
	Point3d p23;
	p23=bf.setv(tmp4.x-tmp3.x,tmp4.y-tmp3.y,tmp4.z-tmp3.z);
	len=bf.norm(p23);
}


void RNASimulation::coord2parambb(vector<Residueatoms> &decstr,int start, int end,int simmode)
{
	int i, j;
	char aat;
	int aaind[4];
	char indnt;
	Point3d p12,p23;
	int atom1, atom2, atom3, atom4;
	Point3s tmp1,tmp2,tmp3,tmp4;
	BasicFunc bf;
	int realstart;

	if(start == 0)
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
		coord2len(tmp3, tmp4, decstr[0].len[1]);
		
		atom2=ParamsList[aat][2][1];atom3=ParamsList[aat][2][2];atom4=ParamsList[aat][2][3];
		tmp2=decstr[0].ptv[atom2];tmp3=decstr[0].ptv[atom3];tmp4=decstr[0].ptv[atom4];
		coord2anglen(tmp2, tmp3, tmp4, decstr[0].len[2], decstr[0].ang[2]);
		
		for(j=3;j<AccAtomParamNum[aat][0];j++){
			atom1=ParamsList[aat][j][0];atom2=ParamsList[aat][j][1];atom3=ParamsList[aat][j][2];atom4=ParamsList[aat][j][3];
			tmp1=decstr[0].ptv[atom1];tmp2=decstr[0].ptv[atom2];tmp3=decstr[0].ptv[atom3];tmp4=decstr[0].ptv[atom4];
			coord2params(tmp1, tmp2, tmp3, tmp4, decstr[0].len[j], decstr[0].ang[j], decstr[0].dih[j]);
		}
		realstart=1;
	}else{
		realstart=start;
	}

	if(simmode == 0 || simmode == 1)
	{
		//}
		
		for(i=realstart;i<end;i++)
		{
			aat=decstr[i].aaa;
			atom1=ParamsList[aat][0][0];atom2=ParamsList[aat][0][1];atom3=ParamsList[aat][0][2];atom4=ParamsList[aat][0][3];
			tmp1=decstr[i-1].ptv[atom1];tmp2=decstr[i-1].ptv[atom2];tmp3=decstr[i-1].ptv[atom3];tmp4=decstr[i].ptv[atom4];
			coord2params(tmp1, tmp2, tmp3, tmp4, decstr[i].len[0], decstr[i].ang[0], decstr[i].dih[0]);
			
			atom1=ParamsList[aat][1][0];atom2=ParamsList[aat][1][1];atom3=ParamsList[aat][1][2];atom4=ParamsList[aat][1][3];
			tmp1=decstr[i-1].ptv[atom1];tmp2=decstr[i-1].ptv[atom2];tmp3=decstr[i].ptv[atom3];tmp4=decstr[i].ptv[atom4];
			coord2params(tmp1, tmp2, tmp3, tmp4, decstr[i].len[1], decstr[i].ang[1], decstr[i].dih[1]);
			
			atom1=ParamsList[aat][2][0];atom2=ParamsList[aat][2][1];atom3=ParamsList[aat][2][2];atom4=ParamsList[aat][2][3];
			tmp1=decstr[i-1].ptv[atom1];tmp2=decstr[i].ptv[atom2];tmp3=decstr[i].ptv[atom3];tmp4=decstr[i].ptv[atom4];
			coord2params(tmp1, tmp2, tmp3, tmp4, decstr[i].len[2], decstr[i].ang[2], decstr[i].dih[2]);
			
			for(j=3;j<AccAtomParamNum[aat][0];j++){
				atom1=ParamsList[aat][j][0];atom2=ParamsList[aat][j][1];atom3=ParamsList[aat][j][2];atom4=ParamsList[aat][j][3];
				tmp1=decstr[i].ptv[atom1];tmp2=decstr[i].ptv[atom2];tmp3=decstr[i].ptv[atom3];tmp4=decstr[i].ptv[atom4];
				
				coord2params(tmp1, tmp2, tmp3, tmp4, decstr[i].len[j], decstr[i].ang[j], decstr[i].dih[j]);
			}
		}
	}else if(simmode == 2)
	{
		for(i=realstart;i<end;i++)
		{
			aat=decstr[i].aaa;

			decstr[i].len[0]=LENGTHPARAM[aat][0].mean;
			decstr[i].ang[0]=INNERANGPARAM[aat][0].mean;
			atom1=ParamsList[aat][0][0];atom2=ParamsList[aat][0][1];atom3=ParamsList[aat][0][2];atom4=ParamsList[aat][0][3];
			tmp1=decstr[i-1].ptv[atom1];tmp2=decstr[i-1].ptv[atom2];tmp3=decstr[i-1].ptv[atom3];tmp4=decstr[i].ptv[atom4];
			coord2dangle(tmp1, tmp2, tmp3, tmp4, decstr[i].dih[0]);

			decstr[i].len[1]=LENGTHPARAM[aat][1].mean;
			decstr[i].ang[1]=INNERANGPARAM[aat][1].mean;
			atom1=ParamsList[aat][1][0];atom2=ParamsList[aat][1][1];atom3=ParamsList[aat][1][2];atom4=ParamsList[aat][1][3];
			tmp1=decstr[i-1].ptv[atom1];tmp2=decstr[i-1].ptv[atom2];tmp3=decstr[i].ptv[atom3];tmp4=decstr[i].ptv[atom4];
			coord2dangle(tmp1, tmp2, tmp3, tmp4, decstr[i].dih[1]);
			
			decstr[i].len[2]=LENGTHPARAM[aat][2].mean;
			decstr[i].ang[2]=INNERANGPARAM[aat][2].mean;
			atom1=ParamsList[aat][2][0];atom2=ParamsList[aat][2][1];atom3=ParamsList[aat][2][2];atom4=ParamsList[aat][2][3];
			tmp1=decstr[i-1].ptv[atom1];tmp2=decstr[i].ptv[atom2];tmp3=decstr[i].ptv[atom3];tmp4=decstr[i].ptv[atom4];
			coord2dangle(tmp1, tmp2, tmp3, tmp4, decstr[i].dih[2]);
			
			for(j=3;j<AccAtomParamNum[aat][0];j++)
			{
				decstr[i].len[j]=LENGTHPARAM[aat][j].mean;
				decstr[i].ang[j]=INNERANGPARAM[aat][j].mean;
				atom1=ParamsList[aat][j][0];atom2=ParamsList[aat][j][1];atom3=ParamsList[aat][j][2];atom4=ParamsList[aat][j][3];
				tmp1=decstr[i].ptv[atom1];tmp2=decstr[i].ptv[atom2];tmp3=decstr[i].ptv[atom3];tmp4=decstr[i].ptv[atom4];
				coord2dangle(tmp1, tmp2, tmp3, tmp4, decstr[i].dih[j]);
			}
		}
	}
}


void RNASimulation::coord2paramsugar(vector<Residueatoms> &decstr,int start,int end,int simmode)
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
	
	if(simmode == 0)
	{
		if(start == 0){
			aat=decstr[0].aaa;
			for(j=AccAtomParamNum[aat][0];j<AccAtomParamNum[aat][2];j++)
			{
				if(j == 10 || j == 11){
					atom1=FIRSTOPParamsList[aat][j][0]; atom2=FIRSTOPParamsList[aat][j][1]; atom3=FIRSTOPParamsList[aat][j][2]; atom4=FIRSTOPParamsList[aat][j][3];
				}else{
					atom1=ParamsList[aat][j][0]; atom2=ParamsList[aat][j][1]; atom3=ParamsList[aat][j][2]; atom4=ParamsList[aat][j][3];
				}
				tmp1=decstr[0].ptv[atom1]; tmp2=decstr[0].ptv[atom2]; tmp3=decstr[0].ptv[atom3]; tmp4=decstr[0].ptv[atom4];
				coord2params(tmp1, tmp2, tmp3, tmp4, len, ang, dih);
				decstr[0].dih[j]=dih; decstr[0].len[j]=len; decstr[0].ang[j]=ang;
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
				}else if(j == 11){
					tmp1=decstr[i].ptv[atom1]; tmp2=decstr[i-1].ptv[atom2]; tmp3=decstr[i].ptv[atom3]; tmp4=decstr[i].ptv[atom4];
				}else{
					tmp1=decstr[i].ptv[atom1]; tmp2=decstr[i].ptv[atom2]; tmp3=decstr[i].ptv[atom3]; tmp4=decstr[i].ptv[atom4];
				}
				coord2params(tmp1, tmp2, tmp3, tmp4, decstr[i].len[j], decstr[i].ang[j], decstr[i].dih[j]);
			}
		}
	}
	else if(simmode == 1)
	{
		if(start == 0){
			aat=decstr[0].aaa;
			for(j=AccAtomParamNum[aat][0];j<AccAtomParamNum[aat][2];j++)
			{		
				auto it = std::find(MODTARGET.begin(), MODTARGET.end(), j);
				if(it != MODTARGET.end())
				{
					if(j == 10 || j == 11){
						atom1=FIRSTOPParamsList[aat][j][0]; atom2=FIRSTOPParamsList[aat][j][1]; atom3=FIRSTOPParamsList[aat][j][2]; atom4=FIRSTOPParamsList[aat][j][3];
					}else{
						atom1=ParamsList[aat][j][0]; atom2=ParamsList[aat][j][1]; atom3=ParamsList[aat][j][2]; atom4=ParamsList[aat][j][3];
					}
					tmp1=decstr[0].ptv[atom1]; tmp2=decstr[0].ptv[atom2]; tmp3=decstr[0].ptv[atom3]; tmp4=decstr[0].ptv[atom4];
					coord2params(tmp1, tmp2, tmp3, tmp4, len, ang, dih);
					decstr[0].dih[j]=dih; decstr[0].len[j]=len; decstr[0].ang[j]=ang;
				}else{
					if(j == 10){
						decstr[0].len[j] = OPLENGTHPARAM[aat][0].mean;
						decstr[0].ang[j] = OPINNERANGPARAM[aat][0].mean;
						decstr[0].dih[j] = OPTORSIONPARAM[aat][0].mean;
					}else if(j == 11){
						decstr[0].len[j] = OPLENGTHPARAM[aat][1].mean;
						decstr[0].ang[j] = OPINNERANGPARAM[aat][1].mean;
						decstr[0].dih[j] = OPTORSIONPARAM[aat][1].mean;
					}else{
						decstr[0].len[j]=LENGTHPARAM[aat][j].mean;
						decstr[0].ang[j]=INNERANGPARAM[aat][j].mean;
						decstr[0].dih[j]=TORSIONPARAM[aat][j].mean;
					}
				}
			}
			realstart=1;
		}else{
			realstart=start;
		}
		for(i=realstart;i<end;i++)
		{
			aat=decstr[i].aaa;
			for(j=AccAtomParamNum[aat][0];j<AccAtomParamNum[aat][2];j++)
			{
				auto it = std::find(MODTARGET.begin(), MODTARGET.end(), j);
				if(it != MODTARGET.end())
				{
					atom1=ParamsList[aat][j][0]; atom2=ParamsList[aat][j][1]; atom3=ParamsList[aat][j][2]; atom4=ParamsList[aat][j][3];
					if(j == 10){	
						tmp1=decstr[i-1].ptv[atom1]; tmp2=decstr[i].ptv[atom2]; tmp3=decstr[i].ptv[atom3]; tmp4=decstr[i].ptv[atom4];
					}else if(j == 11){
						tmp1=decstr[i].ptv[atom1]; tmp2=decstr[i-1].ptv[atom2]; tmp3=decstr[i].ptv[atom3]; tmp4=decstr[i].ptv[atom4];
					}else{
						tmp1=decstr[i].ptv[atom1]; tmp2=decstr[i].ptv[atom2]; tmp3=decstr[i].ptv[atom3]; tmp4=decstr[i].ptv[atom4];
					}
					coord2params(tmp1, tmp2, tmp3, tmp4, len, ang, dih);
					decstr[i].dih[j]=dih; decstr[i].len[j]=len; decstr[i].ang[j]=ang;
				}else{
					decstr[i].len[j]=LENGTHPARAM[aat][j].mean;
					decstr[i].ang[j]=INNERANGPARAM[aat][j].mean;
					decstr[i].dih[j]=TORSIONPARAM[aat][j].mean;
				}
			}
		}
	}else if(simmode == 2)
	{
		if(start == 0){
			aat=decstr[0].aaa;
			for(j=AccAtomParamNum[aat][0];j<AccAtomParamNum[aat][2];j++)
			{
				if(j == 10){
					decstr[0].len[j] = OPLENGTHPARAM[aat][0].mean;
					decstr[0].ang[j] = OPINNERANGPARAM[aat][0].mean;
					decstr[0].dih[j] = OPTORSIONPARAM[aat][0].mean;
				}else if(j == 11){
					decstr[0].len[j] = OPLENGTHPARAM[aat][1].mean;
					decstr[0].ang[j] = OPINNERANGPARAM[aat][1].mean;
					decstr[0].dih[j] = OPTORSIONPARAM[aat][1].mean;
				}else{
					decstr[0].len[j]=LENGTHPARAM[aat][j].mean;
					decstr[0].ang[j]=INNERANGPARAM[aat][j].mean;
					decstr[0].dih[j]=TORSIONPARAM[aat][j].mean;
				}
			}
			realstart=1;
		}else{
			realstart=start;
		}

		for(i=realstart;i<end;i++)
		{
			aat=decstr[i].aaa;
			for(j=AccAtomParamNum[aat][0];j<AccAtomParamNum[aat][2];j++)
			{
				decstr[i].len[j]=LENGTHPARAM[aat][j].mean;
				decstr[i].ang[j]=INNERANGPARAM[aat][j].mean;
				decstr[i].dih[j]=TORSIONPARAM[aat][j].mean;
			}
		}
	}
}


void RNASimulation::coord2paramadd(vector<Residueatoms> &decstr,int start,int end)
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
	
	if(start == 0){
		aat=decstr[0].aaa;
		for(j=1;j<AddAngleParamsList[aat].size();j++)
		{
			atom1=AddAngleParamsList[aat][j][0]; atom2=AddAngleParamsList[aat][j][1]; atom3=AddAngleParamsList[aat][j][2];
			tmp1=decstr[0].ptv[atom1]; tmp2=decstr[0].ptv[atom2]; tmp3=decstr[0].ptv[atom3];
			coord2ang(tmp1, tmp2, tmp3, ang);
			decstr[0].addang[j]=ang;
		}

		for(j=0;j<AddLengthParamsList[aat].size();j++)
		{
			atom1=AddLengthParamsList[aat][j][0]; atom2=AddLengthParamsList[aat][j][1];
			tmp1=decstr[0].ptv[atom1]; tmp2=decstr[0].ptv[atom2];
			coord2len(tmp1, tmp2, len);
			decstr[0].addlen[j]=len;
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
		tmp1=decstr[i-1].ptv[atom1]; tmp2=decstr[i].ptv[atom2]; tmp3=decstr[i].ptv[atom3];
		coord2ang(tmp1, tmp2, tmp3, ang);
		decstr[i].addang[0]=ang;

		//cout<<to_string(i) + "-" + to_string(0) + "-" + to_string(ang)<<endl;

		for(j=1;j<AddAngleParamsList[aat].size();j++)
		{
			//cout<<aat<<endl;

			atom1=AddAngleParamsList[aat][j][0]; atom2=AddAngleParamsList[aat][j][1]; atom3=AddAngleParamsList[aat][j][2];
			tmp1=decstr[i].ptv[atom1]; tmp2=decstr[i].ptv[atom2]; tmp3=decstr[i].ptv[atom3];
			coord2ang(tmp1, tmp2, tmp3, ang);
			decstr[i].addang[j]=ang;
			
			//cout<<to_string(i) + "-" + to_string(j) + "-" + to_string(ang)<<endl;
		}

		for(j=0;j<AddLengthParamsList[aat].size();j++)
		{
			atom1=AddLengthParamsList[aat][j][0]; atom2=AddLengthParamsList[aat][j][1];
			tmp1=decstr[i].ptv[atom1]; tmp2=decstr[i].ptv[atom2];
			coord2len(tmp1, tmp2, len);
			decstr[i].addlen[j]=len;
		}
	}
}

void RNASimulation::coord2paramhydrogen(vector<Residueatoms> &decstr,int start,int end,int simmode)
{
	int i,j;
	char aat;
	int atom1, atom2, atom3, atom4;
	int realstart;
	int batomnum;
	int aaind[4];
	int idx;
	Point3s tmp1,tmp2,tmp3,tmp4;
	float len, ang, dih;
	
	if(simmode == 0 | simmode == 1)
	{
		for(i=start;i<end;i++)
		{
			aat=decstr[i].aaa;
			idx=0;
			for(j=AccAtomParamNum[aat][3];j<AccAtomParamNum[aat][4];j++)
			{
				atom1=HydParamsList[aat][idx][0]; atom2=HydParamsList[aat][idx][1]; atom3=HydParamsList[aat][idx][2]; atom4=HydParamsList[aat][idx][3];
				tmp1=decstr[i].ptv[atom1]; tmp2=decstr[i].ptv[atom2]; tmp3=decstr[i].ptv[atom3]; tmp4=decstr[i].ptv[atom4];
				coord2params(tmp1, tmp2, tmp3, tmp4, len, ang, dih);
				decstr[i].hyddih[idx]=dih; decstr[i].hydlen[idx]=len; decstr[i].hydang[idx]=ang;
				idx+=1;
			}
		}
	}
	if(simmode == 2)
	{
		for(i=start;i<end;i++)
		{
			aat=decstr[i].aaa;
			idx=0;
			for(j=AccAtomParamNum[aat][3];j<AccAtomParamNum[aat][4];j++)
			{
				decstr[i].hyddih[idx]=HYDTORSIONPARAM[aat][idx].mean; 
				decstr[i].hydlen[idx]=HYDLENGTHPARAM[aat][idx].mean;
				decstr[i].hydang[idx]=HYDINNERANGPARAM[aat][idx].mean; 
				//cout<<to_string(decstr[i].hyddih[idx]) + "-" + to_string(decstr[i].hydlen[idx]) + "-" + to_string(decstr[i].hydang[idx])<<endl;
				idx+=1;
			}
		}
	}
}

void RNASimulation::coord2paramhydrogenBase(vector<Residueatoms> &decstr,int start,int end,int simmode)
{
	int i,j;
	char aat;
	int atom1, atom2, atom3, atom4;
	int realstart;
	int batomnum;
	int aaind[4];
	Point3s tmp1,tmp2,tmp3,tmp4;
	float len, ang, dih;
	
	if(simmode == 0 | simmode == 1)
	{
		for(i=start;i<end;i++)
		{
			aat=decstr[i].aaa;
			for(auto idx: BaseHydParamIdx[aat])
			{
				atom1=HydParamsList[aat][idx][0]; atom2=HydParamsList[aat][idx][1]; atom3=HydParamsList[aat][idx][2]; atom4=HydParamsList[aat][idx][3];
				tmp1=decstr[i].ptv[atom1]; tmp2=decstr[i].ptv[atom2]; tmp3=decstr[i].ptv[atom3]; tmp4=decstr[i].ptv[atom4];
				coord2params(tmp1, tmp2, tmp3, tmp4, len, ang, dih);
				decstr[i].hyddih[idx]=dih; decstr[i].hydlen[idx]=len; decstr[i].hydang[idx]=ang;
			}
		}
	}
	if(simmode == 2)
	{
		for(i=start;i<end;i++)
		{
			aat=decstr[i].aaa;
			for(auto idx: BaseHydParamIdx[aat])
			{
				decstr[i].hyddih[idx]=HYDTORSIONPARAM[aat][idx].mean; 
				decstr[i].hydlen[idx]=HYDLENGTHPARAM[aat][idx].mean;
				decstr[i].hydang[idx]=HYDINNERANGPARAM[aat][idx].mean; 
			}
		}
	}
}












