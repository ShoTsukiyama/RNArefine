
#include "Superpose.hpp"

void normalize(double a[3])
{
  double  b = 1/sqrt((double)(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]));
  a[0] *= b;
  a[1] *= b;
  a[2] *= b;
}

void cross(double a[3], double b[3], double c[3])
{
  a[0] = b[1]*c[2] - b[2]*c[1];
  a[1] = b[2]*c[0] - b[0]*c[2];
  a[2] = b[0]*c[1] - b[1]*c[0];
}

void VectorTimesMatrix(double *v1, const vector<vector<double> > &rotation_matrix,  double *v2)
{
	v2[0]=v1[0]*rotation_matrix[0][0]+v1[1]*rotation_matrix[1][0]+v1[2]*rotation_matrix[2][0];
	v2[1]=v1[0]*rotation_matrix[0][1]+v1[1]*rotation_matrix[1][1]+v1[2]*rotation_matrix[2][1];
	v2[2]=v1[0]*rotation_matrix[0][2]+v1[1]*rotation_matrix[1][2]+v1[2]*rotation_matrix[2][2];
}

void ChangeCoor(vector<double> &v1, const vector<vector<double> > &rotation_matrix, const vector<double> &translate_vec,vector<double> &v2)
{
	v2[0]=v1[0]*rotation_matrix[0][0]+v1[1]*rotation_matrix[1][0]+v1[2]*rotation_matrix[2][0]+translate_vec[0];
	v2[1]=v1[0]*rotation_matrix[0][1]+v1[1]*rotation_matrix[1][1]+v1[2]*rotation_matrix[2][1]+translate_vec[1];
	v2[2]=v1[0]*rotation_matrix[0][2]+v1[1]*rotation_matrix[1][2]+v1[2]*rotation_matrix[2][2]+translate_vec[2];
}

void RotateCoor(const vector<vector<double> > &target_stru,const vector<vector<double> > &refrence_stru, vector<vector<double> > &RotMatix, vector<double> &TranVect)
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
	
	
	TranVect.assign(3, 0.);
	RotMatix.assign(3, TranVect);
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

