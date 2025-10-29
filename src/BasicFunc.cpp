// BasicFunc.cpp: implementation of the BasicFunc class.
//
//////////////////////////////////////////////////////////////////////
#include "BasicFunc.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

#define ROTATE(a,i,j,k,l) { g = a[i][j]; \
                            h = a[k][l]; \
                            a[i][j] = g-s*(h+g*tau); \
                            a[k][l] = h+s*(g-h*tau); }

BasicFunc::BasicFunc()
{

}

BasicFunc::~BasicFunc()
{

}

//ノルムを算出する関数
double BasicFunc::norm(Point3d p1)
{
	double t=sqrt(p1.x*p1.x+p1.y*p1.y+p1.z*p1.z+EPS);
	return t;
}

double BasicFunc::norm(Point3s p1)
{
	double t=sqrt(p1.x*p1.x+p1.y*p1.y+p1.z*p1.z+EPS);
	return t;
}

//内積を算出するための関数
double BasicFunc::dotv(Point3s p1,Point3s p2)
{
	double temp;
	temp=p1.x*p2.x+p1.y*p2.y+p1.z*p2.z;
	return temp;
}

double BasicFunc::dotv(Point3d p1,Point3d p2)
{
	double temp;
	temp=p1.x*p2.x+p1.y*p2.y+p1.z*p2.z;
	return temp;
}

Point3d BasicFunc::setv(double tx,double ty,double tz)
{
	Point3d temp;
	temp.x=tx;
	temp.y=ty;
	temp.z=tz;
	return temp;
}

Point3s BasicFunc::setvs(double tx,double ty,double tz)
{
	Point3s temp;
	temp.x=tx;
	temp.y=ty;
	temp.z=tz;
	return temp;
}

void BasicFunc::basenormalize(Point3s &p1)
{
  double  b = 1/sqrt((double)(p1.x*p1.x + p1.y*p1.y + p1.z*p1.z));
  p1.x *= b;
  p1.y *= b;
  p1.z *= b;
}

Point3s BasicFunc::basesetv(double tx,double ty,double tz)
{
	Point3s temp;
	temp.x=tx;
	temp.y=ty;
	temp.z=tz;
	return temp;
}

Point3s BasicFunc::calcsubtract(Point3s p1, Point3s p2)
{
  Point3s p3;
  p3=setvs(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
  return p3;
}

Point3s BasicFunc::calcadd(Point3s p1, Point3s p2)
{
  Point3s p3;
  p3=setvs(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z);
  return p3;
}

Point3s BasicFunc::basecross(Point3s p1, Point3s p2)
{
  return basesetv(
            p1.y * p2.z - p1.z * p2.y,
            p1.z * p2.x - p1.x * p2.z,
            p1.x * p2.y - p1.y * p2.x
            );
}

void BasicFunc::getcoordsystem(Point3s res1N, Point3s res1C1, Point3s res1C2, Point3s &x, Point3s &y, Point3s &z)
{
  Point3s u,v;

  u=calcsubtract(res1C1, res1N);
  v=calcsubtract(res1C2, res1N);

  basenormalize(u);
  basenormalize(v);  

  y=calcadd(u,v);
  basenormalize(y);

  z=basecross(u,v);
  basenormalize(z);

  x=basecross(y,z);
  basenormalize(x);
}

//dihedral angleを算出する関数
double BasicFunc::phi(double xi,double yi,double zi,double xj,double yj,double zj,double xk,
		   double yk,double zk,double xl,double yl,double zl)
{
	double xij,yij,zij,
           xkj,ykj,zkj,
           xkl,ykl,zkl,
           dxi,dyi,dzi,
           gxi,gyi,gzi,
           bi,bk,ct,
           boi2,boj2,
           z1,z2,ap,s,
           bioj,bjoi;
	
	
	/* Calculate the vectors C,B,C                                       */
	xij = xi - xj;
	yij = yi - yj;
	zij = zi - zj;
	xkj = xk - xj;
	ykj = yk - yj;
	zkj = zk - zj;
	xkl = xk - xl;
	ykl = yk - yl;
	zkl = zk - zl;
	
	/* Calculate the normals to the two planes n1 and n2
	this is given as the cross products:
	AB x BC
	--------- = n1
	|AB x BC|
	
	  BC x CD
      --------- = n2
      |BC x CD|
	*/
	dxi = yij * zkj - zij * ykj;     /* Normal to plane 1  面1の法線              */
	dyi = zij * xkj - xij * zkj;    
	dzi = xij * ykj - yij * xkj;
	gxi = zkj * ykl - ykj * zkl;     /* Mormal to plane 2  面2の法線              */
	gyi = xkj * zkl - zkj * xkl;
	gzi = ykj * xkl - xkj * ykl;
	
	/* Calculate the length of the two normals                           */
	bi = dxi * dxi + dyi * dyi + dzi * dzi;
	bk = gxi * gxi + gyi * gyi + gzi * gzi;
	ct = dxi * gxi + dyi * gyi + dzi * gzi;
	
	//boi2 = 1./bi;
	//boj2 = 1./bk;
	bi   = (double)sqrt((double)bi);
	bk   = (double)sqrt((double)bk);
	if(bi<EPS*0.01 || bk<EPS*0.01) return 180;
	z1   = 1./bi;
	z2   = 1./bk;
	//bioj = bi * z2;
	//bjoi = bk * z1;
	ct   = ct * z1 * z2;
	if (ct >  1.0)   ct = 1.0;
	if (ct < (-1.0)) ct = -1.0;
	ap   = acos(ct);
	
	s = xkj * (dzi * gyi - dyi * gzi)
	  + ykj * (dxi * gzi - dzi * gxi)
	  + zkj * (dyi * gxi - dxi * gyi);
	
	if (s < 0.0) ap = -ap;
	
	ap = (ap > 0.0) ? PI-ap : -(PI+ap);

	// angle
	ap *= (double)180.0/PI;
	if(ap<0) ap+=360;
	return(ap);
}


void BasicFunc::transposeRotationMatrix(vector<vector<double> > &RotMatix){
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            // 元の行列の (i, j) と (j, i) の要素を交換する
            std::swap(RotMatix[i][j], RotMatix[j][i]);
        }
    }
}

void BasicFunc::Transform(vector<double> xyz, vector<vector<double> > rotationMatrix, vector<double> TranVect, vector<double> &xyzres){
    // 逆回転を行う
	xyzres[0]=rotationMatrix[0][0] * (xyz[0] - TranVect[0]) + rotationMatrix[0][1] * (xyz[1] - TranVect[1]) + rotationMatrix[0][2] * (xyz[2] - TranVect[2]);
	xyzres[1]=rotationMatrix[1][0] * (xyz[0] - TranVect[0]) + rotationMatrix[1][1] * (xyz[1] - TranVect[1]) + rotationMatrix[1][2] * (xyz[2] - TranVect[2]);
	xyzres[2]=rotationMatrix[2][0] * (xyz[0] - TranVect[0]) + rotationMatrix[2][1] * (xyz[1] - TranVect[1]) + rotationMatrix[2][2] * (xyz[2] - TranVect[2]);
}

void BasicFunc::Transform(vector<double> xyz_tmp, vector<vector<double> > rotationMatrix, vector<double> TranVect, Point3d &tp){
    // 逆回転を行う
    tp.x=rotationMatrix[0][0] * (xyz_tmp[0] - TranVect[0]) + rotationMatrix[0][1] * (xyz_tmp[1] - TranVect[1]) + rotationMatrix[0][2] * (xyz_tmp[2] - TranVect[2]);
	tp.y=rotationMatrix[1][0] * (xyz_tmp[0] - TranVect[0]) + rotationMatrix[1][1] * (xyz_tmp[1] - TranVect[1]) + rotationMatrix[1][2] * (xyz_tmp[2] - TranVect[2]);
    tp.z=rotationMatrix[2][0] * (xyz_tmp[0] - TranVect[0]) + rotationMatrix[2][1] * (xyz_tmp[1] - TranVect[1]) + rotationMatrix[2][2] * (xyz_tmp[2] - TranVect[2]);
}

void BasicFunc::cross(double a[3], double b[3], double c[3])
{
  a[0] = b[1]*c[2] - b[2]*c[1];
  a[1] = b[2]*c[0] - b[0]*c[2];
  a[2] = b[0]*c[1] - b[1]*c[0];
}

Point3s BasicFunc::cross(Point3s p1, Point3s p2)
{
  return setvs(
            p1.y * p2.z - p1.z * p2.y,
            p1.z * p2.x - p1.x * p2.z,
            p1.x * p2.y - p1.y * p2.x
            );
}

Point3d BasicFunc::crossang(Point3d p1,Point3d p2)
{
  	Point3d out;
	out.x = p1.y * p2.z - p1.z * p2.y;
	out.y = p1.z * p2.x - p1.x * p2.z;
	out.z = p1.x * p2.y - p1.y * p2.x;
	return out;
}

void BasicFunc::normalize(double a[3])
{
  double  b = 1/sqrt((double)(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]));
  a[0] *= b;
  a[1] *= b;
  a[2] *= b;
}

void BasicFunc::normalize(Point3s &p1)
{
  double  b = 1/sqrt((double)(p1.x*p1.x + p1.y*p1.y + p1.z*p1.z));
  p1.x *= b;
  p1.y *= b;
  p1.z *= b;
}

double dot_prod(double a[3], double b[3])
{
  return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

int jacobi3(double a[3][3], double d[3], double v[3][3], int* n_rot)
{
  int count, k, i, j;
  double tresh, theta, tau, t, sum, s, h, g, c, b[3], z[3];

  /*Initialize v to the identity matrix.*/
  for (i=0; i<3; i++) 
  { 
    for (j=0; j<3; j++) 
      v[i][j] = 0.0;
    v[i][i] = 1.0;
  }

  /* Initialize b and d to the diagonal of a */
  for (i=0; i<3; i++) 
    b[i] = d[i] = a[i][i];

  /* z will accumulate terms */
  for (i=0; i<3; i++) 
    z[i] = 0.0; 
  
  *n_rot = 0;

  /* 50 tries */
  for (count=0; count<50; count++)     
  {

    /* sum off-diagonal elements */
    sum = 0.0;
    for (i=0; i<2; i++) 
    {
      for (j=i+1; j<3; j++)
         sum += fabs(a[i][j]);
    }

    /* if converged to machine underflow */
    if (sum == 0.0) 
      return(1);

    /* on 1st three sweeps... */
    if (count < 3) 
      tresh = sum * 0.2 / 9.0;    
    else       
      tresh = 0.0;      

    for (i=0; i<2; i++) 
    {
      for (j=i+1; j<3; j++) 
      {
        g = 100.0 * fabs(a[i][j]);

        /*  after four sweeps, skip the rotation if
         *   the off-diagonal element is small 
         */
        if ( count > 3  &&  fabs(d[i])+g == fabs(d[i])
              &&  fabs(d[j])+g == fabs(d[j]) ) 
        {
          a[i][j] = 0.0;
        } 
        else if (fabs(a[i][j]) > tresh) 
        {
          h = d[j] - d[i];
          
          if (fabs(h)+g == fabs(h))
          {
            t = a[i][j] / h;
          }
          else 
          {
            theta = 0.5 * h / (a[i][j]);
            t = 1.0 / ( fabs(theta) +
                        (double)sqrt(1.0 + theta*theta) );
            if (theta < 0.0) 
              t = -t;
          }
          
          c = 1.0 / (double) sqrt(1 + t*t);
          s = t * c;
          tau = s / (1.0 + c);
          h = t * a[i][j];

          z[i] -= h;
          z[j] += h;
          d[i] -= h;
          d[j] += h;

          a[i][j] = 0.0;

          for (k=0; k<=i-1; k++) 
            ROTATE(a, k, i, k, j)

          for (k=i+1; k<=j-1; k++) 
            ROTATE(a, i, k, k, j)

          for (k=j+1; k<3; k++) 
            ROTATE(a, i, k, j, k)

          for (k=0; k<3; k++) 
            ROTATE(v, k, i, k, j)

          ++(*n_rot);
        }
      }
    }

    for (i=0; i<3; i++) 
    {
      b[i] += z[i];
      d[i] = b[i];
      z[i] = 0.0;
    }
  }

  printf("Too many iterations in jacobi3\n");
  return (0);
}  

int diagonalize_symmetric(double matrix[3][3], 
                          double eigen_vec[3][3], 
                          double eigenval[3])
{
  int n_rot, i, j, k;
  double vec[3][3];
  double val; 
  
  if (!jacobi3(matrix, eigenval, vec, &n_rot)) 
  {
    printf("convergence failed\n");
    return (0);
  }

  /* sort solutions by eigenval */
  for (i=0; i<3; i++) 
  {
    k = i;
    val = eigenval[i];
    
    for (j=i+1; j<3; j++)
      if (eigenval[j] >= val)
      { 
        k = j;
        val = eigenval[k];
      }
       
    if (k != i) 
    {
      eigenval[k] = eigenval[i];
      eigenval[i] = val;
      for (j=0; j<3; j++) 
      {
        val = vec[j][i];
        vec[j][i] = vec[j][k];
        vec[j][k] = val;
      }
    }
  }

  /* transpose such that first index refers to solution index */
  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      eigen_vec[i][j] = vec[j][i];

  return (1);
}

void BasicFunc::VectorTimesMatrix(double *v1, const vector<vector<double> > &rotation_matrix,  double *v2)
{
	v2[0]=v1[0]*rotation_matrix[0][0]+v1[1]*rotation_matrix[1][0]+v1[2]*rotation_matrix[2][0];
	v2[1]=v1[0]*rotation_matrix[0][1]+v1[1]*rotation_matrix[1][1]+v1[2]*rotation_matrix[2][1];
	v2[2]=v1[0]*rotation_matrix[0][2]+v1[1]*rotation_matrix[1][2]+v1[2]*rotation_matrix[2][2];
}

void BasicFunc::ChangeCoor(vector<double> &v1, const vector<vector<double> > &rotation_matrix, const vector<double> &translate_vec,Point3d &v2)
{
	v2.x=v1[0]*rotation_matrix[0][0]+v1[1]*rotation_matrix[1][0]+v1[2]*rotation_matrix[2][0]+translate_vec[0];
	v2.y=v1[0]*rotation_matrix[0][1]+v1[1]*rotation_matrix[1][1]+v1[2]*rotation_matrix[2][1]+translate_vec[1];
	v2.z=v1[0]*rotation_matrix[0][2]+v1[1]*rotation_matrix[1][2]+v1[2]*rotation_matrix[2][2]+translate_vec[2];
}

void BasicFunc::RotateCoor(const vector<vector<double> > target_stru,const vector<vector<double> > refrence_stru, vector<vector<double> > &RotMatix, vector<double> &TranVect)
{
	int i, j, k, n;
	
	const int n_list=target_stru.size();
	
	//save two structure
	double ref_xlist[n_list][3],mov_xlist[n_list][3];
	for(i=0;i<n_list;i++) 
	{
		for(j=0;j<3;j++)
		{
			ref_xlist[i][j]=target_stru[i][j];
			mov_xlist[i][j]=refrence_stru[i][j];
		}
	}
	
	//translate two structure to the centre of mass, calculate correlation matrix
	double ref_com[3];
	double mov_com[3];
	double mov_to_ref[3];
      double R[3][3];
	
	/* calculate the centre of mass */
	for (i=0; i<3; i++)
	{ 
		mov_com[i] = 0.0;
		ref_com[i] = 0.0;
	}
	
	for (n=0; n<n_list; n++) 
	{
		for (i=0; i<3; i++)
		{ 
			mov_com[i] += mov_xlist[n][i];
			ref_com[i] += ref_xlist[n][i];
		}
	}
	
	for (i=0; i<3; i++)
	{
		mov_com[i] /= n_list;
		ref_com[i] /= n_list;
		mov_to_ref[i] = ref_com[i] - mov_com[i];
	}
	
	/* shift mov_xlist and ref_xlist to centre of mass */
	for (n=0; n<n_list; n++) 
	{
		for (i=0; i<3; i++)
		{ 
		  mov_xlist[n][i] -= mov_com[i];
		  ref_xlist[n][i] -= ref_com[i];
		}
	}
	
	/* initialize */
	for (i=0; i<3; i++)
	{
		for (j=0; j<3; j++) R[i][j] = 0.0;
	}
	
	for (n=0; n<n_list; n++) 
	{
	/*
	 * correlation matrix R:   
	 *   R[i,j) = sum(over n): y(n,i) * x(n,j)  
	 *   where x(n) and y(n) are two vector sets   
	 */
		for (i=0; i<3; i++)
		{
		  	for (j=0; j<3; j++) R[i][j] += mov_xlist[n][i] * ref_xlist[n][j];
		}
	}
	
	//calculate rotation matrix and translate vector
	double Rt[3][3], RtR[3][3];
	double left_eigenvec[3][3], right_eigenvec[3][3], eigenval[3];
	double v[3];
	
	/* build Rt, transpose of R  */
	for (i=0; i<3; i++)
	{
		for (j=0; j<3; j++) Rt[i][j] = R[j][i];
	}
	
	/* make symmetric RtR = Rt X R */
	for (i=0; i<3; i++) 
	{
		for (j=0; j<3; j++)
		{
		  	RtR[i][j] = 0.0;
		  	for (k = 0; k<3; k++) RtR[i][j] += Rt[k][i] * R[j][k];
		}
	}
	
	diagonalize_symmetric(RtR, right_eigenvec, eigenval);
	
	
	/* right_eigenvec's should be an orthogonal system but could be left
	* or right-handed. Let's force into right-handed system.
	*/
	cross(&right_eigenvec[2][0], &right_eigenvec[0][0], &right_eigenvec[1][0]);
	
	/* From the Kabsch algorithm, the eigenvec's of RtR
	* are identical to the right_eigenvec's of R.
	* This means that left_eigenvec = R x right_eigenvec 
	*/
	for (i=0; i<3; i++) 
	{
		for (j=0; j<3; j++) left_eigenvec[i][j] = dot_prod(&right_eigenvec[i][0], &Rt[j][0]);
	}
	
	for (i=0; i<3; i++) normalize(&left_eigenvec[i][0]);
	
	/* 
	* Force left_eigenvec[2] to be orthogonal to the other vectors.
	* First check if the rotational matrices generated from the 
	* orthogonal eigenvectors are in a right-handed or left-handed
	* co-ordinate system - given by sigma. Sigma is needed to
	* resolve this ambiguity in calculating the RMSD.
	*/
	cross(v, &left_eigenvec[0][0], &left_eigenvec[1][0]);
	
	for (i=0; i<3; i++) left_eigenvec[2][i] = v[i]; 
	
	TranVect.clear();
	TranVect.resize(3, 0.0);
	RotMatix.clear();
	RotMatix.resize(3, vector<double>(3, 0.0));

	/* calc optimal rotation matrix that minimises residual */
	for (i=0;i<3; i++)
	{
		for (j=0; j<3; j++) 
		{
			for (k=0; k<3; k++) 
			      RotMatix[i][j] += left_eigenvec[k][i] * right_eigenvec[k][j];
		}
	}
	double temp[3];
	VectorTimesMatrix(ref_com, RotMatix, temp);
	TranVect[0]=refrence_stru[0][0]-mov_xlist[0][0]-temp[0];
	TranVect[1]=refrence_stru[0][1]-mov_xlist[0][1]-temp[1];
	TranVect[2]=refrence_stru[0][2]-mov_xlist[0][2]-temp[2];
	
}

//inner angleを算出するための関数
double BasicFunc::angv(Point3s p1,Point3s p2)
{
	double t1=norm(p1);
	double t2=norm(p2);
	if(t1<EPS || t2<EPS)
	{
		return 0;
	}
	double t3=dotv(p1,p2)/t1/t2;
	if(t3<-1.0)
	{
		t3=-1.0;
	}
	else if(t3>1.0)
	{
		t3=1.0;
	}
	return(acos(t3)*180.0 / M_PI);
	
	/*
	double dotval = dotv(p1,p2);
	point3d crossval = crossang(p1,p2);

	double angle = atan2(norm(crossval), dotval);
	//double angle = atan(crossval.z/dotval);
	angle = angle * 180.0 / M_PI;

	if(angle>180){
		angle=angle-180;
	}
	else if(angle<0) {
		angle=angle+180;
	}
	return angle;
	*/
}

//inner angleを算出するための関数
double BasicFunc::angv(Point3d p1,Point3d p2)
{
	double t1=norm(p1);
	double t2=norm(p2);
	if(t1<EPS || t2<EPS)
	{
		return 0;
	}
	double t3=dotv(p1,p2)/t1/t2;
	if(t3<-1.0)
	{
		t3=-1.0;
	}
	else if(t3>1.0)
	{
		t3=1.0;
	}
	return(acos(t3)*180.0 / M_PI);
	
	/*
	double dotval = dotv(p1,p2);
	point3d crossval = crossang(p1,p2);

	double angle = atan2(norm(crossval), dotval);
	//double angle = atan(crossval.z/dotval);
	angle = angle * 180.0 / M_PI;

	if(angle>180){
		angle=angle-180;
	}
	else if(angle<0) {
		angle=angle+180;
	}
	return angle;
	*/
}

void  BasicFunc::q2rot(double *q,double rmax[])
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

Point3d BasicFunc::scal(Point3d p1,double f)
{
	Point3d temp;
	temp.x=f*p1.x;
	temp.y=f*p1.y;
	temp.z=f*p1.z;
	return temp;
}

Point3d BasicFunc::prod(Point3d p1,Point3d p2)//from p1 to p2
{
	Point3d temp;
	temp.x=p1.y*p2.z-p1.z*p2.y;
	temp.y=p1.z*p2.x-p1.x*p2.z;
	temp.z=p1.x*p2.y-p1.y*p2.x;
	return temp;

}

Point3d BasicFunc::mmat(double mat[9],Point3d p1)
{
	Point3d temp;
	temp.x=p1.x*mat[0]+p1.y*mat[1]+p1.z*mat[2];
	temp.y=p1.x*mat[3]+p1.y*mat[4]+p1.z*mat[5];
	temp.z=p1.x*mat[6]+p1.y*mat[7]+p1.z*mat[8];
	return temp;
}

Point3d BasicFunc::unit(Point3d p1)
{
	Point3d temp;
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

Point3s BasicFunc::units(Point3s p1)
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

int BasicFunc::maxnormal(Point3d p1)
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
int BasicFunc::minnormal(Point3d p1)
{
	int i;
	double tem,temp[3];	
	temp[0]=fabs(p1.x);
	temp[1]=fabs(p1.y);
	temp[2]=fabs(p1.z);
	i=0;
	tem=temp[0];
	if(temp[1]<tem)
	{
		i=1;
		tem=temp[1];
	}
	if(temp[2]<tem)
	{
		i=2;
		tem=temp[2];
	}
	return i;
}

Point3d BasicFunc::ranv(double fac)
{
	Point3d tp;
	double dtheta,dphi,dr;
//	dr=0.00001+fac*rand()/double(RAND_MAX+1.0);
//	dtheta=(1.0*rand()/double(RAND_MAX+1.0))*PI;
//	dphi=(2.0*rand()/double(RAND_MAX+1.0)-1.0)*PI;
	dr=0.00001+fac*Random();
	dtheta=(1.0*Random())*PI;
	dphi=(2.0*Random()-1.0)*PI;
	tp.x=dr*sin(dtheta)*cos(dphi);
	tp.y=dr*sin(dtheta)*sin(dphi);
	tp.z=dr*cos(dtheta);
	return(tp);
}

bool BasicFunc::tor2pos22(float xi,float yi,float zi,float xj,float yj,float zj,float xk,
		   float yk,float zk,float tang,float tleng,float tinner,float *xl,float *yl,float *zl)
{
	Point3d p12,p23,e1,e2,e3;
	Point3d p1,p2,p3;
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
	tpangle=(PI-tinner)/2.0;
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

void BasicFunc::getnormalvec(Point3s center1, Point3s res11, Point3s res12, Point3s &z){
  Point3s u,v;

  u=calcsubtract(res11, center1);
  v=calcsubtract(res12, center1);

  z=basecross(u,v);
  basenormalize(z);

  //baseangle = angleBetweenVectors(z1, z2);
}

double BasicFunc::angleBetweenVectors(Point3s a, Point3s b){
    BasicFunc bf;
    double dot = bf.dotv(a, b);
    double normA = bf.norm(a);
    double normB = bf.norm(b);
    double cosTheta = dot / (normA * normB);

	return cosTheta;
    //double angle = acos(cosTheta);
    //double angleInDegrees = angle * (180.0 / M_PI);
    //if(angleInDegrees > 180) angleInDegrees=angleInDegrees-180;
    //if(angleInDegrees < 0) angleInDegrees=angleInDegrees+180;
    //return angleInDegrees;
}



