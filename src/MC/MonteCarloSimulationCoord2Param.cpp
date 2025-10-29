#include "MonteCarloSimulation.hpp"

void MonteCarloSimulation::coord2params(const Point3s &tmp1, const Point3s &tmp2, const Point3s &tmp3, const Point3s &tmp4, float &len, float &ang, float &dih)
{
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

void MonteCarloSimulation::coord2dangle(const Point3s &tmp1, const Point3s &tmp2, const Point3s &tmp3, const Point3s &tmp4, float &dih)
{
	BasicFunc bf;
	Point3d p12,p23;
	dih=bf.phi(tmp1.x,tmp1.y,tmp1.z,
			   tmp2.x,tmp2.y,tmp2.z,
			   tmp3.x,tmp3.y,tmp3.z,
			   tmp4.x,tmp4.y,tmp4.z);
}

void MonteCarloSimulation::coord2anglen(const Point3s &tmp2, const Point3s &tmp3, const Point3s &tmp4, float &len, float &ang)
{
	BasicFunc bf;
	Point3d p12,p23;
	p12=bf.setv(tmp2.x-tmp3.x,tmp2.y-tmp3.y,tmp2.z-tmp3.z);
	p23=bf.setv(tmp4.x-tmp3.x,tmp4.y-tmp3.y,tmp4.z-tmp3.z);
	len=bf.norm(p23);
	ang=bf.angv(p12,p23);
}

void MonteCarloSimulation::coord2ang(const Point3s &tmp2, const Point3s &tmp3, const Point3s &tmp4, float &ang)
{
	BasicFunc bf;
	Point3d p12,p23;
	p12=bf.setv(tmp2.x-tmp3.x,tmp2.y-tmp3.y,tmp2.z-tmp3.z);
	p23=bf.setv(tmp4.x-tmp3.x,tmp4.y-tmp3.y,tmp4.z-tmp3.z);
	ang=bf.angv(p12,p23);
}

void MonteCarloSimulation::coord2len(const Point3s &tmp3, const Point3s &tmp4, float &len)
{
	BasicFunc bf;
	Point3d p23;
	p23=bf.setv(tmp4.x-tmp3.x,tmp4.y-tmp3.y,tmp4.z-tmp3.z);
	len=bf.norm(p23);
}


void MonteCarloSimulation::coord2parambb(vector<Residueatoms> &decstr,int start, int end)
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
		decstr[0].dih[0]=MCTORSIONPARAM[aat][0][0];
		decstr[0].dih[1]=MCTORSIONPARAM[aat][1][0];
		decstr[0].dih[2]=MCTORSIONPARAM[aat][2][0];
		decstr[0].ang[0]=MCINNERANGPARAM[aat][0][0];
		decstr[0].ang[1]=MCINNERANGPARAM[aat][1][0];
		decstr[0].len[0]=MCLENGTHPARAM[aat][0][0];
	
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
}


void MonteCarloSimulation::coord2paramsugar(vector<Residueatoms> &decstr,int start,int end)
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


void MonteCarloSimulation::coord2paramadd(vector<Residueatoms> &decstr,int start,int end)
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

		for(j=4;j<AddLengthParamsList[aat].size();j++)
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

		atom1=AddAngleParamsList[aat][0][0]; atom2=AddAngleParamsList[aat][0][1]; atom3=AddAngleParamsList[aat][0][2];
		tmp1=decstr[i-1].ptv[atom1]; tmp2=decstr[i].ptv[atom2]; tmp3=decstr[i].ptv[atom3];
		coord2ang(tmp1, tmp2, tmp3, ang);
		decstr[i].addang[0]=ang;

		for(j=1;j<AddAngleParamsList[aat].size();j++)
		{
			atom1=AddAngleParamsList[aat][j][0]; atom2=AddAngleParamsList[aat][j][1]; atom3=AddAngleParamsList[aat][j][2];
			tmp1=decstr[i].ptv[atom1]; tmp2=decstr[i].ptv[atom2]; tmp3=decstr[i].ptv[atom3];
			coord2ang(tmp1, tmp2, tmp3, ang);
			decstr[i].addang[j]=ang;
		}

		for(j=0;j<4;j++){
			atom1=AddLengthParamsList[aat][j][0]; atom2=AddLengthParamsList[aat][j][1];
			tmp1=decstr[i-1].ptv[atom1]; tmp2=decstr[i].ptv[atom2];
			coord2len(tmp1, tmp2, len);
			decstr[i].addlen[j]=len;
		}

		for(j=4;j<AddLengthParamsList[aat].size();j++)
		{
			atom1=AddLengthParamsList[aat][j][0]; atom2=AddLengthParamsList[aat][j][1];
			tmp1=decstr[i].ptv[atom1]; tmp2=decstr[i].ptv[atom2];
			coord2len(tmp1, tmp2, len);
			decstr[i].addlen[j]=len;
		}
	}
}

void MonteCarloSimulation::coord2paramhydrogen(vector<Residueatoms> &decstr,int start,int end)
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
	
    for(i=start;i<end;i++)
    {
		//cout<<"[[" << i << "]]"<<endl;
        aat=decstr[i].aaa;
        idx=0;
        for(j=AccAtomParamNum[aat][3];j<AccAtomParamNum[aat][4];j++)
        {
			//cout<<j<<endl;
            atom1=HydParamsList[aat][idx][0]; atom2=HydParamsList[aat][idx][1]; atom3=HydParamsList[aat][idx][2]; atom4=HydParamsList[aat][idx][3];
            tmp1=decstr[i].ptv[atom1]; tmp2=decstr[i].ptv[atom2]; tmp3=decstr[i].ptv[atom3]; tmp4=decstr[i].ptv[atom4];
            coord2params(tmp1, tmp2, tmp3, tmp4, len, ang, dih);
            decstr[i].hyddih[idx]=dih; decstr[i].hydlen[idx]=len; decstr[i].hydang[idx]=ang;
            idx+=1;
        }
	}
}

void MonteCarloSimulation::coord2addparamhydrogen(vector<Residueatoms> &decstr,int start,int end)
{
	int i,j;
	char aat;
	int atom2, atom3, atom4;
	Point3s tmp2, tmp3,tmp4;
	float len, ang;
	for(i=start;i<end;i++)
    {
		//cout<<"[[" << i << "]]"<<endl;
        aat=decstr[i].aaa;
		for(j=0;j<AddHydLengthParamsList[aat].size();j++)
		{
			atom3=AddHydLengthParamsList[aat][j][0]; atom4=AddHydLengthParamsList[aat][j][1];
			tmp3=decstr[i].ptv[atom3]; tmp4=decstr[i].ptv[atom4];
			coord2len(tmp3, tmp4, len);
			decstr[i].addhydlen[j]=len;
			//cout<< i << "-" << atom3 << "-" << atom4 << "-" << decstr[i].addhydlen[j]<<endl;
			//cout<< i << "-" << atom3 << "-" << decstr[i].ptv[atom3].x << "-" << decstr[i].ptv[atom3].y << "-" << decstr[i].ptv[atom3].z<<endl;
			//cout<< i << "-" << atom4 << "-" << decstr[i].ptv[atom4].x << "-" << decstr[i].ptv[atom4].y << "-" << decstr[i].ptv[atom4].z<<endl;
			//cout<<decstr[i].addhydlen[j]<<endl;
		}

		for(j=0;j<AddHydAngleParamsList[aat].size();j++)
		{
			atom2=AddHydAngleParamsList[aat][j][0]; atom3=AddHydAngleParamsList[aat][j][1]; atom4=AddHydAngleParamsList[aat][j][2];
			tmp2=decstr[i].ptv[atom2]; tmp3=decstr[i].ptv[atom3]; tmp4=decstr[i].ptv[atom4];
			coord2ang(tmp2, tmp3, tmp4, ang);
			decstr[i].addhydang[j]=ang;
			//cout<<i << "-" << atom2 << "-" << atom3 << "-" << atom4 << "-" << ang<<endl;
			//cout<<decstr[i].ptv[atom2].x << "-" << decstr[i].ptv[atom2].y << "-" << decstr[i].ptv[atom2].z<<endl;
			//cout<<decstr[i].ptv[atom3].x << "-" << decstr[i].ptv[atom3].y << "-" << decstr[i].ptv[atom3].z<<endl;
			//cout<<decstr[i].ptv[atom4].x << "-" << decstr[i].ptv[atom4].y << "-" << decstr[i].ptv[atom4].z<<endl;
		}
	}
	//exit(1);
}

void MonteCarloSimulation::coord2paramhydrogenBase(vector<Residueatoms> &decstr,int start,int end)
{
	int i,j;
	char aat;
	int atom1, atom2, atom3, atom4;
	int realstart;
	int batomnum;
	int aaind[4];
	Point3s tmp1,tmp2,tmp3,tmp4;
	float len, ang, dih;
	
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

void MonteCarloSimulation::coord2parambase(vector<Residueatoms> &decstr,int start,int end)
{
	int i,j;
	char aat;
	int atom1, atom2, atom3, atom4;
	Point3s tmp1,tmp2,tmp3,tmp4;
	float len, ang, dih;
	int idx;

	for(i=start;i<end;i++)
	{
		aat=decstr[i].aaa;
		idx = 0;
		for(j=AccAtomParamNum[aat][2];j<AccAtomParamNum[aat][3];j++)
		{
			atom1=BaseParamsList[aat][idx][0]; atom2=BaseParamsList[aat][idx][1]; atom3=BaseParamsList[aat][idx][2]; atom4=BaseParamsList[aat][idx][3];
			
			tmp1=decstr[i].ptv[atom1]; tmp2=decstr[i].ptv[atom2]; tmp3=decstr[i].ptv[atom3]; tmp4=decstr[i].ptv[atom4];
			coord2params(tmp1, tmp2, tmp3, tmp4, len, ang, dih);
			decstr[i].dih[j]=dih; 
			decstr[i].len[j]=len; 
			decstr[i].ang[j]=ang;
			idx+=1;
		}
	}
}























