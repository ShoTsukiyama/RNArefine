
#include "AttachHydrogens.hpp"

double norm(Point3s p1)
{
	double t=sqrt(p1.x*p1.x+p1.y*p1.y+p1.z*p1.z+EPS);
	return t;
}

int maxnormal(Point3s p1)
{
	int i;
	double tem,temp[3];	
	temp[0]=fabs(p1.x);
	temp[1]=fabs(p1.y);
	temp[2]=fabs(p1.z);
	i=0;
	tem=temp[0];
	if(temp[1]>tem)
	{
		i=1;
		tem=temp[1];
	}
	if(temp[2]>tem)
	{
		i=2;
		tem=temp[2];
	}
	return i;
}

Point3s unit(Point3s p1)
{
	Point3s temp;
	double t=sqrt(p1.x*p1.x+p1.y*p1.y+p1.z*p1.z);
	if(t<1e-40) 
	{
		temp.x=0;
		temp.y=0;
		temp.z=0;
		return temp;
	}
	temp.x=p1.x/t;
	temp.y=p1.y/t;
	temp.z=p1.z/t;
	return temp;
}

Point3s prod(Point3s p1,Point3s p2)//from p1 to p2
{
	Point3s temp;
	temp.x=p1.y*p2.z-p1.z*p2.y;
	temp.y=p1.z*p2.x-p1.x*p2.z;
	temp.z=p1.x*p2.y-p1.y*p2.x;
	return temp;

}

Point3s scal(Point3s p1,double f)
{
	Point3s temp;
	temp.x=f*p1.x;
	temp.y=f*p1.y;
	temp.z=f*p1.z;
	return temp;
}

void q2rot(double *q,double rmax[])
{
	int i,j;
	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			rmax[i*3+j]=0;
		}
	}
	rmax[0]=q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3];
	rmax[1]=2*(q[1]*q[2]-q[0]*q[3]);
	rmax[2]=2*(q[1]*q[3]+q[0]*q[2]);
	rmax[3]=2*(q[1]*q[2]+q[0]*q[3]);
	rmax[4]=q[0]*q[0]-q[1]*q[1]+q[2]*q[2]-q[3]*q[3];
	rmax[5]=2*(q[2]*q[3]-q[0]*q[1]);
	rmax[6]=2*(q[1]*q[3]-q[0]*q[2]);
	rmax[7]=2*(q[2]*q[3]+q[0]*q[1]);
	rmax[8]=q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3];
}

Point3s mmat(double mat[9],Point3s p1)
{
	Point3s temp;
	temp.x=p1.x*mat[0]+p1.y*mat[1]+p1.z*mat[2];
	temp.y=p1.x*mat[3]+p1.y*mat[4]+p1.z*mat[5];
	temp.z=p1.x*mat[6]+p1.y*mat[7]+p1.z*mat[8];
	return temp;
}

bool tor2pos22(float xi,float yi,float zi,float xj,float yj,float zj,float xk,
    float yk,float zk,float tang,float tleng,float tinner,float *xl,float *yl,float *zl)
{
    Point3s p12,p23,e1,e2,e3;
    Point3s p1,p2,p3;
    double tpangle,tcos,tsin,rmat[9];
    double q[4];
    p12.x=xi-xj;p12.y=yi-yj;p12.z=zi-zj;
    p23.x=xk-xj;p23.y=yk-yj;p23.z=zk-zj;
    ///*
    //if(norm(p12)<epsilon*0.00001 || norm(p23)<epsilon*0.00001 || angv(p12,p23)/(180.0 / M_PI)<epsilon*0.001 || (PI-angv(p12,p23))/(180.0 / M_PI)<epsilon*0.001)
    if(norm(p12)<EPS || norm(p23)<EPS)
    {
    int imax=maxnormal(p23);
    if(imax==0)
    {
        yj+=0.01f;
    }
    else if(imax==1)
    {
        zj+=0.01f;
    }
    else
    {
        xj+=0.01f;
    }
    p12.x=xi-xj;p12.y=yi-yj;p12.z=zi-zj;
    p23.x=xk-xj;p23.y=yk-yj;p23.z=zk-zj;
    printf("make adjustment tor2pos22\n");
    }
    //*/
    e2=unit(p23);
    e3=unit(p12);
    e1=prod(e2,e3);
    e1=unit(e1);
    if(norm(e1)<EPS || norm(e2)<EPS)
    {
    printf("wrong in tor2pos22 [%f %f %f] [%f %f %f] [%f %f %f]\n",xi,yi,zi,xj,yj,zj,xk,yk,zk);
    *xl=xk;*yl=yk;*zl=zk;
    return false;
    }
    p1=scal(e2,tleng);
    tpangle=(PI_MACRO-tinner)/2.0;
    tcos=cos(tpangle);
    tsin=sin(tpangle);
    q[0]=tcos;q[1]=tsin*e1.x;q[2]=tsin*e1.y;q[3]=tsin*e1.z;
    q2rot(q,rmat);
    p2=mmat(rmat,p1);

    tpangle=tang/2.0;
    tcos=cos(tpangle);
    tsin=sin(tpangle);
    q[0]=tcos;q[1]=tsin*e2.x;q[2]=tsin*e2.y;q[3]=tsin*e2.z;
    q2rot(q,rmat);
    p3=mmat(rmat,p2);

    *xl=p3.x+xk;*yl=p3.y+yk;*zl=p3.z+zk;
    return true;
}


void AttachHydrogen(FragmentStructure &stru, const map<string, vector<vector<double>>> &hydparam, const map<string, vector<vector<int>>> &hydparamidx)
{
	Point3s pn;
	bool flagts;
	int idx;
	int atom1, atom2, atom3, atom4;
    string base;
    double tor, ang, leng;

    map<string, int> all_atom_num ={{"A", 33}, {"C", 31}, {"G", 34}, {"U", 30}};

    for(int i=0;i<stru.residues.size();i++){
        base = stru.residues[i].restype;
        if(stru.residues[i].atoms.size() == all_atom_num[base]) continue;
        for(int j=0;j<hydparamidx.at(base).size();j++){
            tor = hydparam.at(base).at(j)[0];
            ang = hydparam.at(base).at(j)[2];
            leng = hydparam.at(base).at(j)[4];

            atom1 = hydparamidx.at(base).at(j)[0];
            atom2 = hydparamidx.at(base).at(j)[1];
            atom3 = hydparamidx.at(base).at(j)[2];
            atom4 = hydparamidx.at(base).at(j)[3];

            FragmentAtoms hydatom;
            flagts=tor2pos22(stru.residues[i].atoms[atom1].coord.x,stru.residues[i].atoms[atom1].coord.y,stru.residues[i].atoms[atom1].coord.z,
                             stru.residues[i].atoms[atom2].coord.x,stru.residues[i].atoms[atom2].coord.y,stru.residues[i].atoms[atom2].coord.z,
                             stru.residues[i].atoms[atom3].coord.x,stru.residues[i].atoms[atom3].coord.y,stru.residues[i].atoms[atom3].coord.z,
                             tor*RADDEG,leng,ang*RADDEG,&hydatom.coord.x,&hydatom.coord.y,&hydatom.coord.z);
            stru.residues[i].atoms.push_back(hydatom);
        }
        //cout<<base << "-" << stru.residues[i].atoms.size()<<endl;
        //for(int j=0;j<stru.residues[i].atoms.size();j++){
            //cout<<j << "-" << stru.residues[i].atoms[j].coord.x << "-" << stru.residues[i].atoms[j].coord.y << "-" << stru.residues[i].atoms[j].coord.z<<endl;
        //}
    }
}

void AttachHydrogen(FragmentStructure &stru, const vector<vector<vector<double>>> &hydrogencoords, const map<string, vector<vector<double>>> &hydparam, const map<string, vector<vector<int>>> &hydparamidx)
{
	Point3s pn;
	bool flagts;
	int idx;
	int atom1, atom2, atom3, atom4;
    string base;
    double tor, ang, leng;

    map<string, int> all_atom_num ={{"A", 33}, {"C", 31}, {"G", 34}, {"U", 30}};

    for(int i=0;i<stru.residues.size();i++){
        base = stru.residues[i].restype;
        if(stru.residues[i].atoms.size() == all_atom_num[base]) continue;
        for(int j=0;j<hydparamidx.at(base).size();j++){
            tor = hydparam.at(base).at(j)[0];
            ang = hydparam.at(base).at(j)[2];
            leng = hydparam.at(base).at(j)[4];

            atom1 = hydparamidx.at(base).at(j)[0];
            atom2 = hydparamidx.at(base).at(j)[1];
            atom3 = hydparamidx.at(base).at(j)[2];
            atom4 = hydparamidx.at(base).at(j)[3];

            FragmentAtoms hydatom;
            flagts=tor2pos22(stru.residues[i].atoms[atom1].coord.x,stru.residues[i].atoms[atom1].coord.y,stru.residues[i].atoms[atom1].coord.z,
                             stru.residues[i].atoms[atom2].coord.x,stru.residues[i].atoms[atom2].coord.y,stru.residues[i].atoms[atom2].coord.z,
                             stru.residues[i].atoms[atom3].coord.x,stru.residues[i].atoms[atom3].coord.y,stru.residues[i].atoms[atom3].coord.z,
                             tor*RADDEG,leng,ang*RADDEG,&hydatom.coord.x,&hydatom.coord.y,&hydatom.coord.z);

            //cout<<hydatom.coord.x << "-" << hydatom.coord.y << "-" << hydatom.coord.z<<endl;

            hydatom.coord.x = hydrogencoords[i][j][0];
            hydatom.coord.y = hydrogencoords[i][j][1];
            hydatom.coord.z = hydrogencoords[i][j][2];

            //cout<<hydatom.coord.x << "-" << hydatom.coord.y << "-" << hydatom.coord.z<<endl;

            stru.residues[i].atoms.push_back(hydatom);
        }
        //cout<<base << "-" << stru.residues[i].atoms.size()<<endl;
        //for(int j=0;j<stru.residues[i].atoms.size();j++){
            //cout<<j << "-" << stru.residues[i].atoms[j].coord.x << "-" << stru.residues[i].atoms[j].coord.y << "-" << stru.residues[i].atoms[j].coord.z<<endl;
        //}
    }
    //exit(1);
}

void UpdateHydrogen(FragmentStructure &stru, const map<string, vector<vector<double>>> &hydparam, const map<string, vector<vector<int>>> &hydparamidx)
{
	Point3s pn;
	bool flagts;
	int idx;
	int atom1, atom2, atom3, atom4;
    string base;
    double tor, ang, leng;

    map<string, int> all_atom_num ={{"A", 33}, {"C", 31}, {"G", 34}, {"U", 30}};
    map<string, int> coarse_grained_atom_num ={{"A", 22}, {"C", 20}, {"G", 23}, {"U", 20}};

    for(int i=0;i<stru.residues.size();i++){
        base = stru.residues[i].restype;
        for(int j=0;j<hydparamidx.at(base).size();j++){
            tor = hydparam.at(base).at(j)[0];
            ang = hydparam.at(base).at(j)[2];
            leng = hydparam.at(base).at(j)[4];

            atom1 = hydparamidx.at(base).at(j)[0];
            atom2 = hydparamidx.at(base).at(j)[1];
            atom3 = hydparamidx.at(base).at(j)[2];
            atom4 = hydparamidx.at(base).at(j)[3];

            FragmentAtoms hydatom;
            flagts=tor2pos22(stru.residues[i].atoms[atom1].coord.x,stru.residues[i].atoms[atom1].coord.y,stru.residues[i].atoms[atom1].coord.z,
                             stru.residues[i].atoms[atom2].coord.x,stru.residues[i].atoms[atom2].coord.y,stru.residues[i].atoms[atom2].coord.z,
                             stru.residues[i].atoms[atom3].coord.x,stru.residues[i].atoms[atom3].coord.y,stru.residues[i].atoms[atom3].coord.z,
                             tor*RADDEG,leng,ang*RADDEG,&hydatom.coord.x,&hydatom.coord.y,&hydatom.coord.z);
            
            stru.residues[i].atoms[coarse_grained_atom_num[base] + j].coord.x = hydatom.coord.x;
            stru.residues[i].atoms[coarse_grained_atom_num[base] + j].coord.y = hydatom.coord.y;
            stru.residues[i].atoms[coarse_grained_atom_num[base] + j].coord.z = hydatom.coord.z;
        }
        //cout<<base << "-" << stru.residues[i].atoms.size()<<endl;
        //for(int j=0;j<stru.residues[i].atoms.size();j++){
            //cout<<j << "-" << stru.residues[i].atoms[j].coord.x << "-" << stru.residues[i].atoms[j].coord.y << "-" << stru.residues[i].atoms[j].coord.z<<endl;
        //}
    }
}

void RemoveHydrogens(FragmentStructure &stru)
{
    string base;
    int atom_num;

    for(int i=0;i<stru.residues.size();i++){
        base = stru.residues[i].restype;
        atom_num = heavily_atom_num[base];
        stru.residues[i].atoms.erase(stru.residues[i].atoms.begin() + atom_num, stru.residues[i].atoms.end());
    }
}
