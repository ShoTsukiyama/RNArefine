// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

bool RNASimulation::param2coordbb(vector<Residueatoms> &decstr,int istart,int iend){
	int i, j;
	Point3s pt,pn,pc;
	BasicFunc bf;
	int atom1, atom2, atom3, atom4;
	bool flagts;
	bool flagwhole=true;
	int realstart;
	int ss;
	char aat;

	if(istart==0)
	{
		aat=decstr[0].aaa;
		
		decstr[0].ptv[0].x=0;
		decstr[0].ptv[0].y=0;
		decstr[0].ptv[0].z=0;
		decstr[0].ptv[1].x=LENGTHPARAM[aat][1].mean;
		decstr[0].ptv[1].y=0;
		decstr[0].ptv[1].z=0;
		decstr[0].ptv[2].x=LENGTHPARAM[aat][1].mean-LENGTHPARAM[aat][2].mean*cos(INNERANGPARAM[aat][2].mean*RADDEG);
		decstr[0].ptv[2].y=LENGTHPARAM[aat][2].mean*sin(INNERANGPARAM[aat][2].mean*RADDEG);
		decstr[0].ptv[2].z=0;

		for(i=3;i<AccAtomParamNum[aat][0];i++){
			atom1=ParamsList[aat][i][0];
			atom2=ParamsList[aat][i][1];
			atom3=ParamsList[aat][i][2];
			atom4=ParamsList[aat][i][3];
			
			flagts=bf.tor2pos22(decstr[0].ptv[atom1].x,decstr[0].ptv[atom1].y,decstr[0].ptv[atom1].z,
							    decstr[0].ptv[atom2].x,decstr[0].ptv[atom2].y,decstr[0].ptv[atom2].z,
							    decstr[0].ptv[atom3].x,decstr[0].ptv[atom3].y,decstr[0].ptv[atom3].z,
							    decstr[0].dih[i]*RADDEG,decstr[0].len[i],decstr[0].ang[i]*RADDEG,
								&pn.x,&pn.y,&pn.z);
			decstr[0].ptv[atom4].x=pn.x;
			decstr[0].ptv[atom4].y=pn.y;
			decstr[0].ptv[atom4].z=pn.z;
		}
		realstart=1;

	}
	else realstart=istart;

	for(i=realstart;i<iend;i++)
	{
		aat=decstr[i].aaa;
		ss=decstr[i].ss;

		atom1=ParamsList[aat][0][0];
		atom2=ParamsList[aat][0][1];
		atom3=ParamsList[aat][0][2];
		atom4=ParamsList[aat][0][3];
		flagts=bf.tor2pos22(decstr[i-1].ptv[atom1].x,decstr[i-1].ptv[atom1].y,decstr[i-1].ptv[atom1].z,
							decstr[i-1].ptv[atom2].x,decstr[i-1].ptv[atom2].y,decstr[i-1].ptv[atom2].z,
							decstr[i-1].ptv[atom3].x,decstr[i-1].ptv[atom3].y,decstr[i-1].ptv[atom3].z,
							decstr[i].dih[0]*RADDEG,decstr[i].len[0],decstr[i].ang[0]*RADDEG,
							&pn.x,&pn.y,&pn.z);
		decstr[i].ptv[atom4].x=pn.x;
		decstr[i].ptv[atom4].y=pn.y;
		decstr[i].ptv[atom4].z=pn.z;

		atom1=ParamsList[aat][1][0];
		atom2=ParamsList[aat][1][1];
		atom3=ParamsList[aat][1][2];
		atom4=ParamsList[aat][1][3];
		flagts=bf.tor2pos22(decstr[i-1].ptv[atom1].x,decstr[i-1].ptv[atom1].y,decstr[i-1].ptv[atom1].z,
							decstr[i-1].ptv[atom2].x,decstr[i-1].ptv[atom2].y,decstr[i-1].ptv[atom2].z,
							decstr[i].ptv[atom3].x,decstr[i].ptv[atom3].y,decstr[i].ptv[atom3].z,
							decstr[i].dih[1]*RADDEG,decstr[i].len[1],decstr[i].ang[1]*RADDEG,
							&pn.x,&pn.y,&pn.z);
		decstr[i].ptv[atom4].x=pn.x;
		decstr[i].ptv[atom4].y=pn.y;
		decstr[i].ptv[atom4].z=pn.z;

		atom1=ParamsList[aat][2][0];
		atom2=ParamsList[aat][2][1];
		atom3=ParamsList[aat][2][2];
		atom4=ParamsList[aat][2][3];
		flagts=bf.tor2pos22(decstr[i-1].ptv[atom1].x,decstr[i-1].ptv[atom1].y,decstr[i-1].ptv[atom1].z,
							decstr[i].ptv[atom2].x,decstr[i].ptv[atom2].y,decstr[i].ptv[atom2].z,
							decstr[i].ptv[atom3].x,decstr[i].ptv[atom3].y,decstr[i].ptv[atom3].z,
							decstr[i].dih[2]*RADDEG,decstr[i].len[2],decstr[i].ang[2]*RADDEG,
							&pn.x,&pn.y,&pn.z);
		decstr[i].ptv[atom4].x=pn.x;
		decstr[i].ptv[atom4].y=pn.y;
		decstr[i].ptv[atom4].z=pn.z;

		for(j=3;j<AccAtomParamNum[aat][0];j++){
			atom1=ParamsList[aat][j][0];
			atom2=ParamsList[aat][j][1];
			atom3=ParamsList[aat][j][2];
			atom4=ParamsList[aat][j][3];
			flagts=bf.tor2pos22(decstr[i].ptv[atom1].x,decstr[i].ptv[atom1].y,decstr[i].ptv[atom1].z,
								decstr[i].ptv[atom2].x,decstr[i].ptv[atom2].y,decstr[i].ptv[atom2].z,
								decstr[i].ptv[atom3].x,decstr[i].ptv[atom3].y,decstr[i].ptv[atom3].z,
								decstr[i].dih[j]*RADDEG,decstr[i].len[j],decstr[i].ang[j]*RADDEG,
								&pn.x,&pn.y,&pn.z);
			decstr[i].ptv[atom4].x=pn.x;
			decstr[i].ptv[atom4].y=pn.y;
			decstr[i].ptv[atom4].z=pn.z;
		}
	}
	return true;
}

void RNASimulation::param2coordsugar(vector<Residueatoms> &decstr,int lstart,int lend)
{
	Point3s pn;
	int realstart;
	bool flagts;
	int i,j;
	char aat;
	int atom1, atom2, atom3, atom4;

	if(lstart == 0){
		aat=decstr[0].aaa;
		for(j=AccAtomParamNum[aat][0];j<AccAtomParamNum[aat][2];j++)
		{
			if(j == 10 || j == 11){
				atom1=FIRSTOPParamsList[aat][j][0];atom2=FIRSTOPParamsList[aat][j][1];atom3=FIRSTOPParamsList[aat][j][2];atom4=FIRSTOPParamsList[aat][j][3];
			}else{
				atom1=ParamsList[aat][j][0];atom2=ParamsList[aat][j][1];atom3=ParamsList[aat][j][2];atom4=ParamsList[aat][j][3];
			}
			flagts=bf.tor2pos22(decstr[0].ptv[atom1].x,decstr[0].ptv[atom1].y,decstr[0].ptv[atom1].z,
								decstr[0].ptv[atom2].x,decstr[0].ptv[atom2].y,decstr[0].ptv[atom2].z,
								decstr[0].ptv[atom3].x,decstr[0].ptv[atom3].y,decstr[0].ptv[atom3].z,
								decstr[0].dih[j]*RADDEG,decstr[0].len[j],decstr[0].ang[j]*RADDEG,
								&pn.x,&pn.y,&pn.z);
			decstr[0].ptv[atom4].x=pn.x;decstr[0].ptv[atom4].y=pn.y;decstr[0].ptv[atom4].z=pn.z;
		}
		realstart = 1;
	}else{
		realstart=lstart;
	}
	
	for(i=realstart;i<lend;i++){
		aat=decstr[i].aaa;
		for(j=AccAtomParamNum[aat][0];j<AccAtomParamNum[aat][2];j++)
		{
			atom1=ParamsList[aat][j][0];atom2=ParamsList[aat][j][1];atom3=ParamsList[aat][j][2];atom4=ParamsList[aat][j][3];
			
			if(j == 10){
				flagts=bf.tor2pos22(decstr[i-1].ptv[atom1].x,decstr[i-1].ptv[atom1].y,decstr[i-1].ptv[atom1].z,
									decstr[i].ptv[atom2].x,decstr[i].ptv[atom2].y,decstr[i].ptv[atom2].z,
									decstr[i].ptv[atom3].x,decstr[i].ptv[atom3].y,decstr[i].ptv[atom3].z,
									decstr[i].dih[j]*RADDEG,decstr[i].len[j],decstr[i].ang[j]*RADDEG,
									&pn.x,&pn.y,&pn.z);
			}else if(j == 11){
				flagts=bf.tor2pos22(decstr[i].ptv[atom1].x,decstr[i].ptv[atom1].y,decstr[i].ptv[atom1].z,
									decstr[i-1].ptv[atom2].x,decstr[i-1].ptv[atom2].y,decstr[i-1].ptv[atom2].z,
									decstr[i].ptv[atom3].x,decstr[i].ptv[atom3].y,decstr[i].ptv[atom3].z,
									decstr[i].dih[j]*RADDEG,decstr[i].len[j],decstr[i].ang[j]*RADDEG,
									&pn.x,&pn.y,&pn.z);
			}else{
				flagts=bf.tor2pos22(decstr[i].ptv[atom1].x,decstr[i].ptv[atom1].y,decstr[i].ptv[atom1].z,
									decstr[i].ptv[atom2].x,decstr[i].ptv[atom2].y,decstr[i].ptv[atom2].z,
									decstr[i].ptv[atom3].x,decstr[i].ptv[atom3].y,decstr[i].ptv[atom3].z,
									decstr[i].dih[j]*RADDEG,decstr[i].len[j],decstr[i].ang[j]*RADDEG,
									&pn.x,&pn.y,&pn.z);
			}
			decstr[i].ptv[atom4].x=pn.x;decstr[i].ptv[atom4].y=pn.y;decstr[i].ptv[atom4].z=pn.z;
		}
	}
}

/*
Point3d basesetv(double tx,double ty,double tz)
{
	Point3d temp;
	temp.x=tx;
	temp.y=ty;
	temp.z=tz;
	return temp;
}

Point3d calcsubtract(Point3d p1, Point3d p2)
{
  Point3d p3;
  p3=basesetv(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
  return p3;
}

void basenormalize(Point3d &p1)
{
  double  b = 1/sqrt((double)(p1.x*p1.x + p1.y*p1.y + p1.z*p1.z));
  p1.x *= b;
  p1.y *= b;
  p1.z *= b;
}

Point3d calcadd(Point3d p1, Point3d p2)
{
  Point3d p3;
  p3=basesetv(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z);
  return p3;
}

Point3d basecross(Point3d p1, Point3d p2)
{
  return basesetv(
            p1.y * p2.z - p1.z * p2.y,
            p1.z * p2.x - p1.x * p2.z,
            p1.x * p2.y - p1.y * p2.x
            );
}

void getcoordsystem(Point3d res1, Point3d res2, Point3d res3, Point3d &x, Point3d &y, Point3d &z)
{
  Point3d u,v;

  u=calcsubtract(res2, res1);
  v=calcsubtract(res3, res1);

  basenormalize(u);
  basenormalize(v);  

  y=calcadd(u,v);
  basenormalize(y);

  z=basecross(u,v);
  basenormalize(z);

  x=basecross(y,z);
  basenormalize(x);
}
*/

void getrotl2g(Point3s ptx, Point3s pty, Point3s ptz, vector<vector<double>> &rot){
	rot[0][0] = ptx.x; rot[0][1] = ptx.y; rot[0][2] = ptx.z;
	rot[1][0] = pty.x; rot[1][1] = pty.y; rot[1][2] = pty.z;
	rot[2][0] = ptz.x; rot[2][1] = ptz.y; rot[2][2] = ptz.z;
}

void local2global(vector<vector<double>> rot, vector<double> tl, vector<vector<double>> localcoords, vector<vector<double>> &globalcoords){
	vector<vector<double>> globalcoordstmp(globalcoords.size(), vector<double>(globalcoords[1].size()));

	for (int i = 0; i < localcoords.size(); ++i) {
		for (int j = 0; j < localcoords[i].size(); ++j) {
			globalcoords[i][j]=0;
			for (int k = 0; k < 3; ++k) {
				globalcoords[i][j]+=localcoords[i][k]*rot[k][j];
			}
		}

		for (int j = 0; j < localcoords[i].size(); ++j) {
			globalcoords[i][j] = globalcoords[i][j] + tl[j];
		}
	}
}

bool RNASimulation::param2coordbase(vector<Residueatoms> &decstr, int istart, int iend){
	int i,j;
	char aat;
	Point3s bcoord1, bcoord2, bcoord3;
	Point3s ptx, pty, ptz;
	int atom1, atom2, atom3;
	vector<vector<double>> rot(3, vector<double>(3));
	vector<double> tl(3);
	vector<vector<double>> baseatomcoordstmp;

	for(i=istart;i<iend;i++){
		aat=decstr[i].aaa;

		atom1=BASESATOMS[aat][0];atom2=BASESATOMS[aat][1];atom3=BASESATOMS[aat][2];
			
		bcoord1.x=decstr[i].ptv[atom1].x; bcoord1.y=decstr[i].ptv[atom1].y; bcoord1.z=decstr[i].ptv[atom1].z;
		bcoord2.x=decstr[i].ptv[atom2].x; bcoord2.y=decstr[i].ptv[atom2].y; bcoord2.z=decstr[i].ptv[atom2].z;
		bcoord3.x=decstr[i].ptv[atom3].x; bcoord3.y=decstr[i].ptv[atom3].y; bcoord3.z=decstr[i].ptv[atom3].z;

		bf.getcoordsystem(bcoord2, bcoord1, bcoord3, ptx, pty, ptz);
		getrotl2g(ptx, pty, ptz, rot);
		tl[0] = bcoord2.x; tl[1] = bcoord2.y; tl[2] = bcoord2.z;

		baseatomcoordstmp.resize(BASEATOM[aat].size(), vector<double>(3));
		local2global(rot, tl, BASETEMPLATE[aat], baseatomcoordstmp);

		for(j=0;j<BASEATOM[aat].size();j++){
			decstr[i].ptv[BASEATOM[aat][j]].x=baseatomcoordstmp[j][0];
			decstr[i].ptv[BASEATOM[aat][j]].y=baseatomcoordstmp[j][1];
			decstr[i].ptv[BASEATOM[aat][j]].z=baseatomcoordstmp[j][2];
		}
	}
	return true;
}

void RNASimulation::param2coordhydrogen(vector<Residueatoms> &decstr,int lstart,int lend)
{
	Point3s pn;
	bool flagts;
	int i,j;
	char aat;
	int idx;
	int atom1, atom2, atom3, atom4;
	
	for(i=lstart;i<lend;i++){
		aat=decstr[i].aaa;
		idx=0;
		for(j=AccAtomParamNum[aat][3];j<AccAtomParamNum[aat][4];j++)
		{
			atom1=HydParamsList[aat][idx][0]; atom2=HydParamsList[aat][idx][1]; atom3=HydParamsList[aat][idx][2]; atom4=HydParamsList[aat][idx][3];
			
			flagts=bf.tor2pos22(decstr[i].ptv[atom1].x,decstr[i].ptv[atom1].y,decstr[i].ptv[atom1].z,
								decstr[i].ptv[atom2].x,decstr[i].ptv[atom2].y,decstr[i].ptv[atom2].z,
								decstr[i].ptv[atom3].x,decstr[i].ptv[atom3].y,decstr[i].ptv[atom3].z,
								decstr[i].hyddih[idx]*RADDEG,decstr[i].hydlen[idx],decstr[i].hydang[idx]*RADDEG,
								&pn.x,&pn.y,&pn.z);
			
			//cout<<to_string(i) + "-" +  to_string(j) + "-" +  to_string(atom4) + "-" +  to_string(decstr[i].ptv[atom4].x) + "-" + to_string(decstr[i].ptv[atom4].y) + "-" + to_string(decstr[i].ptv[atom4].z)<<endl;

			decstr[i].ptv[atom4].x=pn.x;decstr[i].ptv[atom4].y=pn.y;decstr[i].ptv[atom4].z=pn.z;

			//cout<<to_string(i) + "-" +  to_string(j) + "-" +  to_string(atom4) + "-" +  to_string(decstr[i].ptv[atom4].x) + "-" + to_string(decstr[i].ptv[atom4].y) + "-" + to_string(decstr[i].ptv[atom4].z)<<endl;

			//cout<<to_string(i) + "-" + to_string(j) + "-" + to_string(decstr[i].ptv[atom4].x) + "-" + to_string(decstr[i].ptv[atom4].y) + "-" + to_string(decstr[i].ptv[atom4].z)<<endl;
			idx+=1;
		}
	}
}



























