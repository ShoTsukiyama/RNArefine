// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

void RNASimulation::initrandommovegss(Residueatoms &recrrstr){
	double rmax, rnorm;
	Point3d rp;
	double randval=Random();
	
	for(int j=0;j<PAIRLENGTHTARGETMOVNUM;j++){
		rmax=0.5+0.5*Random();
		rp=bf.ranv(rmax);
		recrrstr.ptv[j].x+=rp.x;
		recrrstr.ptv[j].y+=rp.y;
		recrrstr.ptv[j].z+=rp.z;
	}
}

void RNASimulation::resetflaggss(vector<bool> &residueflag, vector<vector<bool>> &atomflag, int numseq)
{
	atomflag.clear();
	atomflag.resize(numseq, vector<bool>(PAIRLENGTHTARGETMOVNUM, false));

	residueflag.clear();
	residueflag.resize(numseq, false);
}

double RNASimulation::calcinitstrgss(vector<Residueatoms> restr, vector<Residueatoms> decstr, vector<int> nodes, vector<double> &energy, vector<bool> residueflag, int uptype){
	int i, j;
	Point3s vec;
	double rnorm;
	double totalenergy=0.0;

	if(uptype == 0){
		for(i=0;i<nodes.size();i++)
		{
			energy[i]=0.0;
			for(j=0;j<PAIRLENGTHTARGETMOVNUM;j++){
				vec.x = restr[nodes[i]].ptv[j].x - decstr[nodes[i]].ptv[j].x; 
				vec.y = restr[nodes[i]].ptv[j].x - decstr[nodes[i]].ptv[j].x; 
				vec.z = restr[nodes[i]].ptv[j].x - decstr[nodes[i]].ptv[j].x; 
				energy[i]+=bf.norm(vec);
			}
			totalenergy+=energy[i];
		}
	}else if(uptype == 1){
		double iold, inew;
		for(i=0;i<nodes.size();i++)
		{
			if(residueflag[nodes[i]]){
				inew=0.0;
				iold=energy[i];
				for(j=0;j<PAIRLENGTHTARGETMOVNUM;j++){
					vec.x = restr[nodes[i]].ptv[j].x - decstr[nodes[i]].ptv[j].x; 
					vec.y = restr[nodes[i]].ptv[j].x - decstr[nodes[i]].ptv[j].x; 
					vec.z = restr[nodes[i]].ptv[j].x - decstr[nodes[i]].ptv[j].x; 
					inew+=bf.norm(vec);
				}
				totalenergy+=(inew-iold);
			}else{
				totalenergy+=energy[i];
			}
		}
	}else if(uptype == 2){
		double iold, inew;
		for(i=0;i<nodes.size();i++)
		{
			if(residueflag[nodes[i]]){
				inew=0.0;
				iold=energy[i];
				for(j=0;j<PAIRLENGTHTARGETMOVNUM;j++){
					vec.x = restr[nodes[i]].ptv[j].x - decstr[nodes[i]].ptv[j].x; 
					vec.y = restr[nodes[i]].ptv[j].x - decstr[nodes[i]].ptv[j].x; 
					vec.z = restr[nodes[i]].ptv[j].x - decstr[nodes[i]].ptv[j].x; 
					inew+=bf.norm(vec);
				}
				energy[i]+=(inew-iold);
				totalenergy+=energy[i];
			}else{
				totalenergy+=energy[i];
			}
		}
	}

	//cout<<"init-" << totalenergy<<endl;

	return totalenergy;
}


double RNASimulation::calcmsesuptempgss(Point3s *restr1, Point3s *restr2, char aat1, char aat2, vector<Residueatoms> strtemplate){
	int mtype1, mtype2;
	double rnorm;
	Point3s vec;

	int i, j, k;
	vector<vector<double>> RotMatix;
    vector<double> TranVect;
	vector<double> xyztmp(3,0);
	double rmsd, length;

	int parammun = AccAtomParamNum[aat1][2] + AccAtomParamNum[aat2][2];
	vector<double> tmp(3,0);
	vector<vector<double> > xyz_list1(parammun,tmp);
    vector<vector<double> > xyz_list2(parammun,tmp);

	k = 0;
	for(j=0; j<AccAtomParamNum[aat1][2];j++){
		xyz_list1[k][0]=restr1[j].x;
		xyz_list1[k][1]=restr1[j].y;
		xyz_list1[k][2]=restr1[j].z;

		xyz_list2[k][0]=strtemplate[0].ptv[j].x;
		xyz_list2[k][1]=strtemplate[0].ptv[j].y;
		xyz_list2[k][2]=strtemplate[0].ptv[j].z;

		k+=1;
	}

	for(j=0; j<AccAtomParamNum[aat2][2];j++){
		xyz_list1[k][0]=restr2[j].x;
		xyz_list1[k][1]=restr2[j].y;
		xyz_list1[k][2]=restr2[j].z;

		xyz_list2[k][0]=strtemplate[1].ptv[j].x;
		xyz_list2[k][1]=strtemplate[1].ptv[j].y;
		xyz_list2[k][2]=strtemplate[1].ptv[j].z;
		k+=1;
	}

	bf.RotateCoor(xyz_list1,xyz_list2,RotMatix,TranVect);

	rmsd=0;
	k=0;
	for(j=0; j<AccAtomParamNum[aat1][2];j++){
		bf.Transform(xyz_list2[k],RotMatix,TranVect,xyztmp);
		vec.x=restr1[j].x-xyztmp[0];
		vec.y=restr1[j].y-xyztmp[1];
		vec.z=restr1[j].z-xyztmp[2];
		length=bf.norm(vec);
		rmsd+=length;
		k+=1;
	}
	for(j=0; j<AccAtomParamNum[aat2][2];j++){
		bf.Transform(xyz_list2[k],RotMatix,TranVect,xyztmp);
		vec.x=restr2[j].x-xyztmp[0];
		vec.y=restr2[j].y-xyztmp[1];
		vec.z=restr2[j].z-xyztmp[2];
		length=bf.norm(vec);
		rmsd+=length;
		k+=1;
	}
	
	return rmsd;
}

double RNASimulation::calcmseditmpgss(vector<Residueatoms> restr, vector<pair<int, int>> edgeinfo, vector<int> edgetype, map<int, vector<Residueatoms>> dnttemplates, vector<double> &energy, vector<bool> residueflag, int uptype)
{
	int i;
	int basei, basej;
	int mtype1, mtype2;
	char aat1, aat2;
	double rnorm;
	Point3s vec;
	double tmp;
	double iold, inew;
	double totalenergy=0.0;

	if(uptype == 0){
		for(i=0;i<edgeinfo.size();i++){
			if(edgetype[i] == 1){
				basei=edgeinfo[i].first; basej=edgeinfo[i].second;
				aat1=restr[basei].aaa; aat2=restr[basej].aaa;
				energy[i] = 0.0;

				energy[i] = calcmsesuptempgss(restr[basei].ptv, restr[basej].ptv, aat1, aat2, dnttemplates[i]);

				totalenergy+=energy[i];
			}
		}
	}else if(uptype == 1){
		for(i=0;i<edgeinfo.size();i++){
			if(edgetype[i] == 1){
				basei=edgeinfo[i].first; basej=edgeinfo[i].second;
				aat1=restr[basei].aaa; aat2=restr[basej].aaa;
				if(residueflag[basei] | residueflag[basej]){
					inew=0.0;
					iold=energy[i];

					inew = calcmsesuptempgss(restr[basei].ptv, restr[basej].ptv, aat1, aat2, dnttemplates[i]);

					totalenergy+=(inew-iold);
				}
			}
		}
	}else if(uptype == 2){
		for(i=0;i<edgeinfo.size();i++){
			if(edgetype[i] == 1){
				basei=edgeinfo[i].first; basej=edgeinfo[i].second;
				aat1=restr[basei].aaa; aat2=restr[basej].aaa;
				if(residueflag[basei] | residueflag[basej]){
					inew=0.0;
					iold=energy[i];

					inew = calcmsesuptempgss(restr[basei].ptv, restr[basej].ptv, aat1, aat2, dnttemplates[i]);

					energy[i]+=(inew-iold);
					totalenergy+=energy[i];
				}
			}
		}
	}
	return totalenergy;
}


double RNASimulation::calcmseditmplengthgss(vector<Residueatoms> restr, vector<pair<int, int>> edgeinfo, vector<int> edgetype, map<int, TemplateDistMap> templateslength, vector<double> &energy, vector<bool> residueflag, int uptype)
{
	int i;
	int basei, basej;
	int mtype1, mtype2;
	double rnorm;
	Point3s vec;
	double tmp;
	double iold, inew;
	double totalenergy=0.0;

	if(uptype == 0){
		for(i=0;i<edgeinfo.size();i++){
			if(edgetype[i] == 1){
				basei=edgeinfo[i].first; basej=edgeinfo[i].second;

				energy[i] = 0.0;
				//if(basei < basej){
					for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
					{
						for(mtype2=0;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
						{
							vec.x = restr[basei].ptv[mtype1].x - restr[basej].ptv[mtype2].x; 
							vec.y = restr[basei].ptv[mtype1].y - restr[basej].ptv[mtype2].y; 
							vec.z = restr[basei].ptv[mtype1].z - restr[basej].ptv[mtype2].z;

							rnorm=bf.norm(vec);
							
							tmp = rnorm-templateslength[i].interdistmap[mtype1][mtype2];
							energy[i] += (tmp*tmp);
						}
					}
					for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
					{
						for(mtype2=mtype1+1;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
						{
							vec.x = restr[basei].ptv[mtype1].x - restr[basei].ptv[mtype2].x; 
							vec.y = restr[basei].ptv[mtype1].y - restr[basei].ptv[mtype2].y; 
							vec.z = restr[basei].ptv[mtype1].z - restr[basei].ptv[mtype2].z;

							rnorm=bf.norm(vec);

							tmp = rnorm-templateslength[i].intradistmap[0][mtype1][mtype2];
							energy[i] += (tmp*tmp);
						}
					}
					for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
					{
						for(mtype2=mtype1+1;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
						{
							vec.x = restr[basej].ptv[mtype1].x - restr[basej].ptv[mtype2].x; 
							vec.y = restr[basej].ptv[mtype1].y - restr[basej].ptv[mtype2].y; 
							vec.z = restr[basej].ptv[mtype1].z - restr[basej].ptv[mtype2].z;

							rnorm=bf.norm(vec);

							tmp = rnorm-templateslength[i].intradistmap[1][mtype1][mtype2];
							energy[i] += (tmp*tmp);
						}
					}
				/*
				}else{
					for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
					{
						for(mtype2=0;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
						{
							vec.x = restr[basej].ptv[mtype2].x - restr[basei].ptv[mtype1].x; 
							vec.y = restr[basej].ptv[mtype2].y - restr[basei].ptv[mtype1].y; 
							vec.z = restr[basej].ptv[mtype2].z - restr[basei].ptv[mtype1].z;

							rnorm=bf.norm(vec);
							
							tmp = rnorm-templateslength[i].interdistmap[mtype2][mtype1];
							energy[i] += (tmp*tmp);
						}
					}
					for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
					{
						for(mtype2=mtype1+1;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
						{
							vec.x = restr[basej].ptv[mtype1].x - restr[basej].ptv[mtype2].x; 
							vec.y = restr[basej].ptv[mtype1].y - restr[basej].ptv[mtype2].y; 
							vec.z = restr[basej].ptv[mtype1].z - restr[basej].ptv[mtype2].z;

							rnorm=bf.norm(vec);

							tmp = rnorm-templateslength[i].intradistmap[0][mtype1][mtype2];
							energy[i] += (tmp*tmp);
						}
					}
					for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
					{
						for(mtype2=mtype1+1;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
						{
							vec.x = restr[basei].ptv[mtype1].x - restr[basei].ptv[mtype2].x; 
							vec.y = restr[basei].ptv[mtype1].y - restr[basei].ptv[mtype2].y; 
							vec.z = restr[basei].ptv[mtype1].z - restr[basei].ptv[mtype2].z;

							rnorm=bf.norm(vec);

							tmp = rnorm-templateslength[i].intradistmap[1][mtype1][mtype2];
							energy[i] += (tmp*tmp);
						}
					}
				}
				*/
				totalenergy+=energy[i];
			}
		}
	}else if(uptype == 1){
		for(i=0;i<edgeinfo.size();i++){
			if(edgetype[i] == 1){
				basei=edgeinfo[i].first; basej=edgeinfo[i].second;
				if(residueflag[basei] | residueflag[basej]){
					inew=0.0;
					iold=energy[i];

					//if(basei < basej){
						for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
						{
							for(mtype2=0;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
							{
								vec.x = restr[basei].ptv[mtype1].x - restr[basej].ptv[mtype2].x; 
								vec.y = restr[basei].ptv[mtype1].y - restr[basej].ptv[mtype2].y; 
								vec.z = restr[basei].ptv[mtype1].z - restr[basej].ptv[mtype2].z;

								rnorm=bf.norm(vec);

								tmp = rnorm-templateslength[i].interdistmap[mtype1][mtype2];
								inew += (tmp*tmp);
							}
						}
						for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
						{
							for(mtype2=mtype1+1;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
							{
								vec.x = restr[basei].ptv[mtype1].x - restr[basei].ptv[mtype2].x; 
								vec.y = restr[basei].ptv[mtype1].y - restr[basei].ptv[mtype2].y; 
								vec.z = restr[basei].ptv[mtype1].z - restr[basei].ptv[mtype2].z;

								rnorm=bf.norm(vec);

								tmp = rnorm-templateslength[i].intradistmap[0][mtype1][mtype2];
								inew += (tmp*tmp);
							}
						}
						for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
						{
							for(mtype2=mtype1+1;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
							{
								vec.x = restr[basej].ptv[mtype1].x - restr[basej].ptv[mtype2].x; 
								vec.y = restr[basej].ptv[mtype1].y - restr[basej].ptv[mtype2].y; 
								vec.z = restr[basej].ptv[mtype1].z - restr[basej].ptv[mtype2].z;

								rnorm=bf.norm(vec);

								tmp = rnorm-templateslength[i].intradistmap[1][mtype1][mtype2];
								inew += (tmp*tmp);
							}
						}
					/*
					}else{
						for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
						{
							for(mtype2=0;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
							{
								vec.x = restr[basej].ptv[mtype2].x - restr[basei].ptv[mtype1].x; 
								vec.y = restr[basej].ptv[mtype2].y - restr[basei].ptv[mtype1].y; 
								vec.z = restr[basej].ptv[mtype2].z - restr[basei].ptv[mtype1].z;

								rnorm=bf.norm(vec);
								
								tmp = rnorm-templateslength[i].interdistmap[mtype2][mtype1];
								inew += (tmp*tmp);

								//cout<<"========"<<endl;
								//cout<<rnorm<<endl;
								//cout<<templateslength[i].interdistmap[mtype2][mtype1]<<endl;
							}
						}
						for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
						{
							for(mtype2=mtype1+1;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
							{
								vec.x = restr[basej].ptv[mtype1].x - restr[basej].ptv[mtype2].x; 
								vec.y = restr[basej].ptv[mtype1].y - restr[basej].ptv[mtype2].y; 
								vec.z = restr[basej].ptv[mtype1].z - restr[basej].ptv[mtype2].z;

								rnorm=bf.norm(vec);

								tmp = rnorm-templateslength[i].intradistmap[0][mtype1][mtype2];
								inew += (tmp*tmp);
							}
						}
						for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
						{
							for(mtype2=mtype1+1;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
							{
								vec.x = restr[basei].ptv[mtype1].x - restr[basei].ptv[mtype2].x; 
								vec.y = restr[basei].ptv[mtype1].y - restr[basei].ptv[mtype2].y; 
								vec.z = restr[basei].ptv[mtype1].z - restr[basei].ptv[mtype2].z;

								rnorm=bf.norm(vec);

								tmp = rnorm-templateslength[i].intradistmap[1][mtype1][mtype2];
								inew += (tmp*tmp);
							}
						}
					}
					*/
					totalenergy+=(inew-iold);
				}else{
					totalenergy+=energy[i];
				}
			}
		}
	}else if(uptype == 2){
		for(i=0;i<edgeinfo.size();i++){
			if(edgetype[i] == 1){
				basei=edgeinfo[i].first; basej=edgeinfo[i].second;

				if(residueflag[basei] | residueflag[basej]){
					inew=0.0;
					iold=energy[i];
					//if(basei < basej){
						for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
						{
							for(mtype2=0;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
							{
								vec.x = restr[basei].ptv[mtype1].x - restr[basej].ptv[mtype2].x; 
								vec.y = restr[basei].ptv[mtype1].y - restr[basej].ptv[mtype2].y; 
								vec.z = restr[basei].ptv[mtype1].z - restr[basej].ptv[mtype2].z;

								rnorm=bf.norm(vec);
								
								tmp = rnorm-templateslength[i].interdistmap[mtype1][mtype2];
								inew += (tmp*tmp);
							}
						}
						for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
						{
							for(mtype2=mtype1+1;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
							{
								vec.x = restr[basei].ptv[mtype1].x - restr[basei].ptv[mtype2].x; 
								vec.y = restr[basei].ptv[mtype1].y - restr[basei].ptv[mtype2].y; 
								vec.z = restr[basei].ptv[mtype1].z - restr[basei].ptv[mtype2].z;

								rnorm=bf.norm(vec);

								tmp = rnorm-templateslength[i].intradistmap[0][mtype1][mtype2];
								inew += (tmp*tmp);
							}
						}
						for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
						{
							for(mtype2=mtype1+1;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
							{
								vec.x = restr[basej].ptv[mtype1].x - restr[basej].ptv[mtype2].x; 
								vec.y = restr[basej].ptv[mtype1].y - restr[basej].ptv[mtype2].y; 
								vec.z = restr[basej].ptv[mtype1].z - restr[basej].ptv[mtype2].z;

								rnorm=bf.norm(vec);

								tmp = rnorm-templateslength[i].intradistmap[1][mtype1][mtype2];
								inew += (tmp*tmp);
							}
						}
					/*
					}else{
						for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
						{
							for(mtype2=0;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
							{
								vec.x = restr[basej].ptv[mtype2].x - restr[basei].ptv[mtype1].x; 
								vec.y = restr[basej].ptv[mtype2].y - restr[basei].ptv[mtype1].y; 
								vec.z = restr[basej].ptv[mtype2].z - restr[basei].ptv[mtype1].z;

								rnorm=bf.norm(vec);
								
								tmp = rnorm-templateslength[i].interdistmap[mtype2][mtype1];
								inew += (tmp*tmp);
							}
						}
						for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
						{
							for(mtype2=mtype1+1;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
							{
								vec.x = restr[basej].ptv[mtype1].x - restr[basej].ptv[mtype2].x; 
								vec.y = restr[basej].ptv[mtype1].y - restr[basej].ptv[mtype2].y; 
								vec.z = restr[basej].ptv[mtype1].z - restr[basej].ptv[mtype2].z;

								rnorm=bf.norm(vec);

								tmp = rnorm-templateslength[i].intradistmap[0][mtype1][mtype2];
								inew += (tmp*tmp);
							}
						}
						for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
						{
							for(mtype2=mtype1+1;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
							{
								vec.x = restr[basei].ptv[mtype1].x - restr[basei].ptv[mtype2].x; 
								vec.y = restr[basei].ptv[mtype1].y - restr[basei].ptv[mtype2].y; 
								vec.z = restr[basei].ptv[mtype1].z - restr[basei].ptv[mtype2].z;

								rnorm=bf.norm(vec);

								tmp = rnorm-templateslength[i].intradistmap[1][mtype1][mtype2];
								inew += (tmp*tmp);
							}
						}
					}
					*/
					energy[i]+=(inew-iold);
					totalenergy+=energy[i];
				}else{
					totalenergy+=energy[i];
				}
			}
		}
	}

	//cout<<"adjlength-" << totalenergy<<endl;
	
	return totalenergy;
}

double RNASimulation::calcmsebptmplengthgss(vector<Residueatoms> restr, vector<pair<int, int>> edgeinfo, vector<int> edgetype, map<int, TemplateDistMap> templateslength, vector<double> &energy, vector<bool> residueflag, int uptype)
{
	int i;
	int basei, basej;
	int mtype1, mtype2;
	double rnorm;
	Point3s vec;
	double tmp;
	double iold, inew;
	double totalenergy=0.0;

	if(uptype == 0){
		for(i=0;i<edgeinfo.size();i++){
			if(edgetype[i] == 2){
				basei=edgeinfo[i].first; basej=edgeinfo[i].second;

				energy[i] = 0.0;
				for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
				{
					for(mtype2=0;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
					{
						vec.x = restr[basei].ptv[mtype1].x - restr[basej].ptv[mtype2].x; 
						vec.y = restr[basei].ptv[mtype1].y - restr[basej].ptv[mtype2].y; 
						vec.z = restr[basei].ptv[mtype1].z - restr[basej].ptv[mtype2].z;

						rnorm=bf.norm(vec);
						
						tmp = rnorm-templateslength[i].interdistmap[mtype1][mtype2];
						energy[i] += (tmp*tmp);
					}
				}
				for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
				{
					for(mtype2=mtype1+1;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
					{
						vec.x = restr[basei].ptv[mtype1].x - restr[basei].ptv[mtype2].x; 
						vec.y = restr[basei].ptv[mtype1].y - restr[basei].ptv[mtype2].y; 
						vec.z = restr[basei].ptv[mtype1].z - restr[basei].ptv[mtype2].z;

						rnorm=bf.norm(vec);

						tmp = rnorm-templateslength[i].intradistmap[0][mtype1][mtype2];
						energy[i] += (tmp*tmp);
					}
				}
				for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
				{
					for(mtype2=mtype1+1;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
					{
						vec.x = restr[basej].ptv[mtype1].x - restr[basej].ptv[mtype2].x; 
						vec.y = restr[basej].ptv[mtype1].y - restr[basej].ptv[mtype2].y; 
						vec.z = restr[basej].ptv[mtype1].z - restr[basej].ptv[mtype2].z;

						rnorm=bf.norm(vec);

						tmp = rnorm-templateslength[i].intradistmap[1][mtype1][mtype2];
						energy[i] += (tmp*tmp);
					}
				}
				totalenergy+=energy[i];
			}
		}
	}else if(uptype == 1){
		for(i=0;i<edgeinfo.size();i++){
			if(edgetype[i] == 2){
				basei=edgeinfo[i].first; basej=edgeinfo[i].second;
				if(residueflag[basei] | residueflag[basej]){
					inew=0.0;
					iold=energy[i];

					for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
					{
						for(mtype2=0;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
						{
							vec.x = restr[basei].ptv[mtype1].x - restr[basej].ptv[mtype2].x; 
							vec.y = restr[basei].ptv[mtype1].y - restr[basej].ptv[mtype2].y; 
							vec.z = restr[basei].ptv[mtype1].z - restr[basej].ptv[mtype2].z;

							rnorm=bf.norm(vec);
							
							tmp = rnorm-templateslength[i].interdistmap[mtype1][mtype2];
							inew += (tmp*tmp);
						}
					}
					for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
					{
						for(mtype2=mtype1+1;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
						{
							vec.x = restr[basei].ptv[mtype1].x - restr[basei].ptv[mtype2].x; 
							vec.y = restr[basei].ptv[mtype1].y - restr[basei].ptv[mtype2].y; 
							vec.z = restr[basei].ptv[mtype1].z - restr[basei].ptv[mtype2].z;

							rnorm=bf.norm(vec);

							tmp = rnorm-templateslength[i].intradistmap[0][mtype1][mtype2];
							inew += (tmp*tmp);
						}
					}
					for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
					{
						for(mtype2=mtype1+1;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
						{
							vec.x = restr[basej].ptv[mtype1].x - restr[basej].ptv[mtype2].x; 
							vec.y = restr[basej].ptv[mtype1].y - restr[basej].ptv[mtype2].y; 
							vec.z = restr[basej].ptv[mtype1].z - restr[basej].ptv[mtype2].z;

							rnorm=bf.norm(vec);

							tmp = rnorm-templateslength[i].intradistmap[1][mtype1][mtype2];
							inew += (tmp*tmp);
						}
					}
					
					totalenergy+=(inew-iold);
				}else{
					totalenergy+=energy[i];
				}
			}
		}
	}else if(uptype == 2){
		for(i=0;i<edgeinfo.size();i++){
			if(edgetype[i] == 2){
				basei=edgeinfo[i].first; basej=edgeinfo[i].second;

				if(residueflag[basei] | residueflag[basej]){
					inew=0.0;
					iold=energy[i];
					
					for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
					{
						for(mtype2=0;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
						{
							vec.x = restr[basei].ptv[mtype1].x - restr[basej].ptv[mtype2].x; 
							vec.y = restr[basei].ptv[mtype1].y - restr[basej].ptv[mtype2].y; 
							vec.z = restr[basei].ptv[mtype1].z - restr[basej].ptv[mtype2].z;

							rnorm=bf.norm(vec);
							
							tmp = rnorm-templateslength[i].interdistmap[mtype1][mtype2];
							inew += (tmp*tmp);
						}
					}
					for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
					{
						for(mtype2=mtype1+1;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
						{
							vec.x = restr[basei].ptv[mtype1].x - restr[basei].ptv[mtype2].x; 
							vec.y = restr[basei].ptv[mtype1].y - restr[basei].ptv[mtype2].y; 
							vec.z = restr[basei].ptv[mtype1].z - restr[basei].ptv[mtype2].z;

							rnorm=bf.norm(vec);

							tmp = rnorm-templateslength[i].intradistmap[0][mtype1][mtype2];
							inew += (tmp*tmp);
						}
					}
					for(mtype1=0;mtype1<PAIRLENGTHTARGETMOVNUM;mtype1++)
					{
						for(mtype2=mtype1+1;mtype2<PAIRLENGTHTARGETMOVNUM;mtype2++)
						{
							vec.x = restr[basej].ptv[mtype1].x - restr[basej].ptv[mtype2].x; 
							vec.y = restr[basej].ptv[mtype1].y - restr[basej].ptv[mtype2].y; 
							vec.z = restr[basej].ptv[mtype1].z - restr[basej].ptv[mtype2].z;

							rnorm=bf.norm(vec);

							tmp = rnorm-templateslength[i].intradistmap[1][mtype1][mtype2];
							inew += (tmp*tmp);
						}
					}
					
					energy[i]+=(inew-iold);
					totalenergy+=energy[i];
				}else{
					totalenergy+=energy[i];
				}
			}
		}
	}

	//cout<<"adjlength-" << totalenergy<<endl;
	
	return totalenergy;
}



double RNASimulation::getouteneanglegss(double val, double mean, double lower, double upper, double slower, double supper, bool &outflag){
	double tmp;

	if(lower <= val && val <= upper){
		outflag=false;
		return 0;
	}else if(lower > val){
		tmp=abs(mean-val);
		outflag=true;
		return tmp;
	}else if(val > upper){
		tmp=abs(val-mean);
		outflag=true;
		return tmp;
	}
	return 0;
}

double RNASimulation::calcenebondanglegss(vector<Residueatoms> restr, vector<int> nodes, vector<double> &energy, vector<vector<bool>> &angleoutflag, vector<bool> residueflag, int numseq, int uptype){
	int i, j;
	Point3s vec;
	double rnorm;
	char aat;
	int atom1, atom2, atom3, atom4;
	int residx1, residx2, residx3, residx4;
	double totalenergy=0.0;
	bool outflag;
	if(uptype == 0){
		for(i=0;i<nodes.size();i++)
		{
			energy[i]=0.0;
			aat=restr[nodes[i]].aaa;

			for(j=0;j<INNERANGLETARGET.size();j++){
				if(nodes[i] == 0){
					if(INNERANGLETARGET[j] == 10){
						energy[i]+=getouteneanglegss(restr[0].ang[INNERANGLETARGET[j]], OPINNERANGPARAM[aat][0].mean, OPINNERANGPARAM[aat][0].lower, OPINNERANGPARAM[aat][0].upper, OPINNERANGPARAM[aat][0].slower, OPINNERANGPARAM[aat][0].supper, outflag);
						angleoutflag[i][j]=outflag;
					}else if(INNERANGLETARGET[j] == 11){
						energy[i]+=getouteneanglegss(restr[0].ang[INNERANGLETARGET[j]], OPINNERANGPARAM[aat][1].mean, OPINNERANGPARAM[aat][1].lower, OPINNERANGPARAM[aat][1].upper, OPINNERANGPARAM[aat][1].slower, OPINNERANGPARAM[aat][1].supper, outflag);
						angleoutflag[i][j]=outflag;
					}else{
						energy[i]+=getouteneanglegss(restr[0].ang[INNERANGLETARGET[j]], INNERANGPARAM[aat][INNERANGLETARGET[j]].mean, INNERANGPARAM[aat][INNERANGLETARGET[j]].lower, INNERANGPARAM[aat][INNERANGLETARGET[j]].upper, INNERANGPARAM[aat][INNERANGLETARGET[j]].slower, INNERANGPARAM[aat][INNERANGLETARGET[j]].supper, outflag);
						angleoutflag[i][j]=outflag;
					}
				}else{
					energy[i]+=getouteneanglegss(restr[nodes[i]].ang[INNERANGLETARGET[j]], INNERANGPARAM[aat][INNERANGLETARGET[j]].mean, INNERANGPARAM[aat][INNERANGLETARGET[j]].lower, INNERANGPARAM[aat][INNERANGLETARGET[j]].upper, INNERANGPARAM[aat][INNERANGLETARGET[j]].slower, INNERANGPARAM[aat][INNERANGLETARGET[j]].supper, outflag);
					angleoutflag[i][j]=outflag;
				}
			}

			for(j=0;j<ADDINNERANGLETARGET.size();j++)
			{
				energy[i]+=getouteneanglegss(restr[nodes[i]].addang[ADDINNERANGLETARGET[j]], ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].mean, ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].lower, ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].upper, ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].slower, ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].supper, outflag);
				angleoutflag[i][j]=outflag;
			}

			totalenergy+=energy[i];
		}
	}else if(uptype == 1){
		double iold, inew;		

		for(i=0;i<nodes.size();i++)
		{
			if(residueflag[nodes[i]] || ((nodes[i]-1 > 0) && (residueflag[nodes[i]-1]))){
				inew=0.0;
				iold=energy[i];
				aat=restr[nodes[i]].aaa;

				for(j=0;j<INNERANGLETARGET.size();j++){
					if(nodes[i] == 0){
						if(INNERANGLETARGET[j] == 10){
							inew+=getouteneanglegss(restr[0].ang[INNERANGLETARGET[j]], OPINNERANGPARAM[aat][0].mean, OPINNERANGPARAM[aat][0].lower, OPINNERANGPARAM[aat][0].upper, OPINNERANGPARAM[aat][0].slower, OPINNERANGPARAM[aat][0].supper, outflag);
						}else if(INNERANGLETARGET[j] == 11){
							inew+=getouteneanglegss(restr[0].ang[INNERANGLETARGET[j]], OPINNERANGPARAM[aat][1].mean, OPINNERANGPARAM[aat][1].lower, OPINNERANGPARAM[aat][1].upper, OPINNERANGPARAM[aat][1].slower, OPINNERANGPARAM[aat][1].supper, outflag);
						}else{
							inew+=getouteneanglegss(restr[0].ang[INNERANGLETARGET[j]], INNERANGPARAM[aat][INNERANGLETARGET[j]].mean, INNERANGPARAM[aat][INNERANGLETARGET[j]].lower, INNERANGPARAM[aat][INNERANGLETARGET[j]].upper, INNERANGPARAM[aat][INNERANGLETARGET[j]].slower, INNERANGPARAM[aat][INNERANGLETARGET[j]].supper, outflag);
						}
					}else{
						inew+=getouteneanglegss(restr[nodes[i]].ang[INNERANGLETARGET[j]], INNERANGPARAM[aat][INNERANGLETARGET[j]].mean, INNERANGPARAM[aat][INNERANGLETARGET[j]].lower, INNERANGPARAM[aat][INNERANGLETARGET[j]].upper, INNERANGPARAM[aat][INNERANGLETARGET[j]].slower, INNERANGPARAM[aat][INNERANGLETARGET[j]].supper, outflag);
					}
				}

				for(j=0;j<ADDINNERANGLETARGET.size();j++)
				{
					inew+=getouteneanglegss(restr[nodes[i]].addang[ADDINNERANGLETARGET[j]], ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].mean, ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].lower, ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].upper, ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].slower, ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].supper, outflag);
				}

				totalenergy+=(inew-iold);
			}else{
				totalenergy+=energy[i];
			}
		}
	}else if(uptype == 2){
		double iold, inew;
		for(i=0;i<nodes.size();i++)
		{
			if(residueflag[nodes[i]] || ((nodes[i]-1 > 0) && (residueflag[nodes[i]-1]))){
				inew=0.0;
				iold=energy[i];
				aat=restr[nodes[i]].aaa;

				for(j=0;j<INNERANGLETARGET.size();j++){
					if(nodes[i] == 0){
						if(INNERANGLETARGET[j] == 10){
							inew+=getouteneanglegss(restr[0].ang[INNERANGLETARGET[j]], OPINNERANGPARAM[aat][0].mean, OPINNERANGPARAM[aat][0].lower, OPINNERANGPARAM[aat][0].upper, OPINNERANGPARAM[aat][0].slower, OPINNERANGPARAM[aat][0].supper, outflag);
							angleoutflag[i][j]=outflag;
						}else if(INNERANGLETARGET[j] == 11){
							inew+=getouteneanglegss(restr[0].ang[INNERANGLETARGET[j]], OPINNERANGPARAM[aat][1].mean, OPINNERANGPARAM[aat][1].lower, OPINNERANGPARAM[aat][1].upper, OPINNERANGPARAM[aat][1].slower, OPINNERANGPARAM[aat][1].supper, outflag);
							angleoutflag[i][j]=outflag;
						}else{
							inew+=getouteneanglegss(restr[0].ang[INNERANGLETARGET[j]], INNERANGPARAM[aat][INNERANGLETARGET[j]].mean, INNERANGPARAM[aat][INNERANGLETARGET[j]].lower, INNERANGPARAM[aat][INNERANGLETARGET[j]].upper, INNERANGPARAM[aat][INNERANGLETARGET[j]].slower, INNERANGPARAM[aat][INNERANGLETARGET[j]].supper, outflag);
							angleoutflag[i][j]=outflag;
						}
					}else{
						inew+=getouteneanglegss(restr[nodes[i]].ang[INNERANGLETARGET[j]], INNERANGPARAM[aat][INNERANGLETARGET[j]].mean, INNERANGPARAM[aat][INNERANGLETARGET[j]].lower, INNERANGPARAM[aat][INNERANGLETARGET[j]].upper, INNERANGPARAM[aat][INNERANGLETARGET[j]].slower, INNERANGPARAM[aat][INNERANGLETARGET[j]].supper, outflag);
						angleoutflag[i][j]=outflag;
					}
				}

				for(j=0;j<ADDINNERANGLETARGET.size();j++)
				{
					inew+=getouteneanglegss(restr[nodes[i]].addang[ADDINNERANGLETARGET[j]], ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].mean, ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].lower, ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].upper, ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].slower, ADDINNERANGPARAM[aat][ADDINNERANGLETARGET[j]].supper, outflag);
					angleoutflag[i][j]=outflag;
				}

				energy[i]+=(inew-iold);
				totalenergy+=energy[i];
			}else{
				totalenergy+=energy[i];
			}
		}
	}

	return totalenergy;
}

double RNASimulation::getoutenelength(double val, double mean, double lower, double upper, double slower, double supper, bool &outflag){
	double tmp;

	if(lower <= val && val <= upper){
		outflag=false;
		return 0;
	}else if(slower > val){
		tmp=abs(mean-val)+1;
		outflag=true;
		return tmp*tmp;
	}else if(val > supper){
		tmp=abs(mean-val)+1;
		outflag=true;
		return tmp*tmp;
	}else if(lower > val){
		tmp=abs(mean-val);
		outflag=true;
		return tmp;
	}else if(val > upper){
		tmp=abs(val-mean);
		outflag=true;
		return tmp;
	}
	return 0;
}

double RNASimulation::calcenebondlengthgss(vector<Residueatoms> restr, vector<int> nodes, vector<double> &energy, vector<vector<bool>> &lengthoutflag, vector<bool> residueflag, int numseq, int uptype){
	int i, j;
	Point3s vec;
	double rnorm;
	char aat;
	int atom1, atom2, atom3, atom4;
	int residx1, residx2, residx3, residx4;
	double totalenergy=0.0;
	bool outflag;
	if(uptype == 0){
		for(i=0;i<nodes.size();i++)
		{
			energy[i]=0.0;
			aat=restr[nodes[i]].aaa;

			for(j=0;j<BLENGTHTARGET.size();j++){
				if(nodes[i] == 0){
					if(BLENGTHTARGET[j] == 10){
						energy[i]+=getoutenelength(restr[0].len[BLENGTHTARGET[j]], OPLENGTHPARAM[aat][0].mean, OPLENGTHPARAM[aat][0].lower, OPLENGTHPARAM[aat][0].upper, OPLENGTHPARAM[aat][0].slower, OPLENGTHPARAM[aat][0].supper, outflag);
						lengthoutflag[i][j]=outflag;
					}else if(BLENGTHTARGET[j] == 11){
						energy[i]+=getoutenelength(restr[0].len[BLENGTHTARGET[j]], OPLENGTHPARAM[aat][1].mean, OPLENGTHPARAM[aat][1].lower, OPLENGTHPARAM[aat][1].upper, OPLENGTHPARAM[aat][1].slower, OPLENGTHPARAM[aat][1].supper, outflag);
						lengthoutflag[i][j]=outflag;
					}else{
						energy[i]+=getoutenelength(restr[0].len[BLENGTHTARGET[j]], LENGTHPARAM[aat][BLENGTHTARGET[j]].mean, LENGTHPARAM[aat][BLENGTHTARGET[j]].lower, LENGTHPARAM[aat][BLENGTHTARGET[j]].upper, LENGTHPARAM[aat][BLENGTHTARGET[j]].slower, LENGTHPARAM[aat][BLENGTHTARGET[j]].supper, outflag);
						lengthoutflag[i][j]=outflag;
					}
				}else{
					energy[i]+=getoutenelength(restr[nodes[i]].len[BLENGTHTARGET[j]], LENGTHPARAM[aat][BLENGTHTARGET[j]].mean, LENGTHPARAM[aat][BLENGTHTARGET[j]].lower, LENGTHPARAM[aat][BLENGTHTARGET[j]].upper, LENGTHPARAM[aat][BLENGTHTARGET[j]].slower, LENGTHPARAM[aat][BLENGTHTARGET[j]].supper, outflag);
					lengthoutflag[i][j]=outflag;
				}
			}

			for(j=0;j<ADDBLENGTHTARGET.size();j++)
			{
				energy[i]+=getoutenelength(restr[nodes[i]].addlen[ADDBLENGTHTARGET[j]], ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[j]].mean, ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[j]].lower, ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[j]].upper, ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[j]].slower, ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[j]].supper, outflag);
				lengthoutflag[i][j]=outflag;
			}

			totalenergy+=energy[i];
		}
	}else if(uptype == 1){
		double iold, inew;	
		for(i=0;i<nodes.size();i++)
		{
			if(residueflag[nodes[i]] || ((nodes[i]-1 > 0) && (residueflag[nodes[i]-1]))){
				inew=0.0;
				iold=energy[i];
				aat=restr[nodes[i]].aaa;

				for(j=0;j<BLENGTHTARGET.size();j++){
					if(nodes[i] == 0){
						if(BLENGTHTARGET[j] == 10){
							inew+=getoutenelength(restr[0].len[BLENGTHTARGET[j]], OPLENGTHPARAM[aat][0].mean, OPLENGTHPARAM[aat][0].lower, OPLENGTHPARAM[aat][0].upper, OPLENGTHPARAM[aat][0].slower, OPLENGTHPARAM[aat][0].supper, outflag);
						}else if(BLENGTHTARGET[j] == 11){
							inew+=getoutenelength(restr[0].len[BLENGTHTARGET[j]], OPLENGTHPARAM[aat][1].mean, OPLENGTHPARAM[aat][1].lower, OPLENGTHPARAM[aat][1].upper, OPLENGTHPARAM[aat][1].slower, OPLENGTHPARAM[aat][1].supper, outflag);
						}else{
							inew+=getoutenelength(restr[0].len[BLENGTHTARGET[j]], LENGTHPARAM[aat][BLENGTHTARGET[j]].mean, LENGTHPARAM[aat][BLENGTHTARGET[j]].lower, LENGTHPARAM[aat][BLENGTHTARGET[j]].upper, LENGTHPARAM[aat][BLENGTHTARGET[j]].slower, LENGTHPARAM[aat][BLENGTHTARGET[j]].supper, outflag);
						}
					}else{
						inew+=getoutenelength(restr[nodes[i]].len[BLENGTHTARGET[j]], LENGTHPARAM[aat][BLENGTHTARGET[j]].mean, LENGTHPARAM[aat][BLENGTHTARGET[j]].lower, LENGTHPARAM[aat][BLENGTHTARGET[j]].upper, LENGTHPARAM[aat][BLENGTHTARGET[j]].slower, LENGTHPARAM[aat][BLENGTHTARGET[j]].supper, outflag);
					}
				}

				for(j=0;j<ADDBLENGTHTARGET.size();j++)
				{
					inew+=getoutenelength(restr[nodes[i]].addlen[ADDBLENGTHTARGET[j]], ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[j]].mean, ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[j]].lower, ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[j]].upper, ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[j]].slower, ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[j]].supper, outflag);
				}
				totalenergy+=(inew-iold);
			}else{
				totalenergy+=energy[i];
			}
		}	
	}else if(uptype == 2){
		double iold, inew;
		for(i=0;i<nodes.size();i++)
		{
			if(residueflag[nodes[i]] || ((nodes[i]-1 > 0) && (residueflag[nodes[i]-1]))){
				inew=0.0;
				iold=energy[i];
				aat=restr[nodes[i]].aaa;

				for(j=0;j<BLENGTHTARGET.size();j++){
					if(nodes[i] == 0){
						if(BLENGTHTARGET[j] == 10){
							inew+=getoutenelength(restr[0].len[BLENGTHTARGET[j]], OPLENGTHPARAM[aat][0].mean, OPLENGTHPARAM[aat][0].lower, OPLENGTHPARAM[aat][0].upper, OPLENGTHPARAM[aat][0].slower, OPLENGTHPARAM[aat][0].supper, outflag);
							lengthoutflag[i][j]=outflag;
						}else if(BLENGTHTARGET[j] == 11){
							inew+=getoutenelength(restr[0].len[BLENGTHTARGET[j]], OPLENGTHPARAM[aat][1].mean, OPLENGTHPARAM[aat][1].lower, OPLENGTHPARAM[aat][1].upper, OPLENGTHPARAM[aat][1].slower, OPLENGTHPARAM[aat][1].supper, outflag);
							lengthoutflag[i][j]=outflag;
						}else{
							inew+=getoutenelength(restr[0].len[BLENGTHTARGET[j]], LENGTHPARAM[aat][BLENGTHTARGET[j]].mean, LENGTHPARAM[aat][BLENGTHTARGET[j]].lower, LENGTHPARAM[aat][BLENGTHTARGET[j]].upper, LENGTHPARAM[aat][BLENGTHTARGET[j]].slower, LENGTHPARAM[aat][BLENGTHTARGET[j]].supper, outflag);
							lengthoutflag[i][j]=outflag;
						}
					}else{
						inew+=getoutenelength(restr[nodes[i]].len[BLENGTHTARGET[j]], LENGTHPARAM[aat][BLENGTHTARGET[j]].mean, LENGTHPARAM[aat][BLENGTHTARGET[j]].lower, LENGTHPARAM[aat][BLENGTHTARGET[j]].upper, LENGTHPARAM[aat][BLENGTHTARGET[j]].slower, LENGTHPARAM[aat][BLENGTHTARGET[j]].supper, outflag);
						lengthoutflag[i][j]=outflag;
					}
				}

				for(j=0;j<ADDBLENGTHTARGET.size();j++)
				{
					inew+=getoutenelength(restr[nodes[i]].addlen[ADDBLENGTHTARGET[j]], ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[j]].mean, ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[j]].lower, ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[j]].upper, ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[j]].slower, ADDLENGTHPARAM[aat][ADDBLENGTHTARGET[j]].supper, outflag);
					lengthoutflag[i][j]=outflag;
				}
				energy[i]+=(inew-iold);
				totalenergy+=energy[i];
			}else{
				totalenergy+=energy[i];
			}
		}	
	}

	return totalenergy;
}

void RNASimulation::initopleng(vector<Residueatoms> &recrrstr, vector<int> nodes){
	int i,j;
	char aat;
	int atom2, atom4;
	int residx2, residx4;
	Point3s tmp1,tmp2,tmp3,tmp4;
	int iend;
	int base;

	for(i=0;i<nodes.size();i++)
	{
		base=nodes[i];
		aat=recrrstr[base].aaa;

		for(j=0;j<AccAtomParamNum[aat][2];j++)
		{
			if(SUBBLENGTHTARGET[j] == 0){
				if(base == 0) continue;
				residx2=base-1; residx4=base;
			}else if(SUBBLENGTHTARGET[j] == 1){
				if(base == 0) continue;
				residx2=base-1; residx4=base;
			}else if(SUBBLENGTHTARGET[j] == 11){
				if(base == 0) continue;
				residx2=base-1; residx4=base;
			}else{
				residx2=base; residx4=base;
			}

			atom2=ParamsList[aat][j][1]; atom4=ParamsList[aat][j][3];
			tmp2=recrrstr[residx2].ptv[atom2];tmp4=recrrstr[residx4].ptv[atom4];

			coord2len(tmp2, tmp4, recrrstr[base].oplen[SUBBLENGTHTARGET[j]]);
		}
	}
}


double RNASimulation::calcenenbondlengthgss(vector<Residueatoms> restr, vector<int> nodes, vector<double> &energy, vector<bool> residueflag, int uptype)
{
	int i, j;
	int basei, basej, base;
	char aat1, aat2, aat;
	int mtype1, mtype2;
	double rnorm;
	Point3s vec;
	double tmp;
	double iold, inew;
	double totalenergy=0.0;
	bool outflag;

	if(uptype == 0){
		for(i=0;i<nodes.size();i++)
		{
			base=nodes[i];
			energy[i]=0.0;
			aat=restr[base].aaa;

			if(base == 0){
				for(j=0;j<AccAtomParamNum[aat][2];j++){
					if(SUBBLENGTHTARGET[j]==0 | SUBBLENGTHTARGET[j]==1 | SUBBLENGTHTARGET[j]==11) continue;
					else{
						//cout<<restr[base].oplen[SUBBLENGTHTARGET[j]] << "-" << SUBLENGTHPARAM[aat][BLENGTHTARGET[j]].mean<<endl;
						energy[i]+=getoutenelength(restr[base].oplen[SUBBLENGTHTARGET[j]], SUBLENGTHPARAM[aat][SUBBLENGTHTARGET[j]].mean, SUBLENGTHPARAM[aat][SUBBLENGTHTARGET[j]].lower, SUBLENGTHPARAM[aat][SUBBLENGTHTARGET[j]].upper, SUBLENGTHPARAM[aat][SUBBLENGTHTARGET[j]].slower, SUBLENGTHPARAM[aat][SUBBLENGTHTARGET[j]].supper, outflag);
					}
				}
			}else{
				for(j=0;j<AccAtomParamNum[aat][2];j++){
					//cout<<restr[base].oplen[SUBBLENGTHTARGET[j]] << "-" << SUBLENGTHPARAM[aat][BLENGTHTARGET[j]].mean<<endl;
					energy[i]+=getoutenelength(restr[base].oplen[SUBBLENGTHTARGET[j]], SUBLENGTHPARAM[aat][SUBBLENGTHTARGET[j]].mean, SUBLENGTHPARAM[aat][SUBBLENGTHTARGET[j]].lower, SUBLENGTHPARAM[aat][SUBBLENGTHTARGET[j]].upper, SUBLENGTHPARAM[aat][SUBBLENGTHTARGET[j]].slower, SUBLENGTHPARAM[aat][SUBBLENGTHTARGET[j]].supper, outflag);
				}
			}
			totalenergy+=energy[i];
		}
	}else if(uptype == 1){
		double iold, inew;	
		for(i=0;i<nodes.size();i++)
		{
			base=nodes[i];

			if(residueflag[base] || ((base-1 > 0) && (residueflag[base-1]))){
				inew=0.0;
				iold=energy[i];
				aat=restr[base].aaa;

				if(base == 0){
					for(j=0;j<AccAtomParamNum[aat][2];j++){
						if(SUBBLENGTHTARGET[j]==0 | SUBBLENGTHTARGET[j]==1 | SUBBLENGTHTARGET[j]==11) continue;
						else{
							if(SUBBLENGTHTARGET[j] == 8) continue;
							//cout<<restr[base].oplen[SUBBLENGTHTARGET[j]] << "-" << SUBLENGTHPARAM[aat][BLENGTHTARGET[j]].mean<<endl;
							inew+=getoutenelength(restr[base].oplen[SUBBLENGTHTARGET[j]], SUBLENGTHPARAM[aat][SUBBLENGTHTARGET[j]].mean, SUBLENGTHPARAM[aat][SUBBLENGTHTARGET[j]].lower, SUBLENGTHPARAM[aat][SUBBLENGTHTARGET[j]].upper, SUBLENGTHPARAM[aat][SUBBLENGTHTARGET[j]].slower, SUBLENGTHPARAM[aat][SUBBLENGTHTARGET[j]].supper, outflag);
						}
					}
				}else{
					for(j=0;j<AccAtomParamNum[aat][2];j++){
						if(SUBBLENGTHTARGET[j] == 8) continue;
						//cout<<restr[base].oplen[SUBBLENGTHTARGET[j]] << "-" << SUBLENGTHPARAM[aat][BLENGTHTARGET[j]].mean<<endl;
						inew+=getoutenelength(restr[base].oplen[SUBBLENGTHTARGET[j]], SUBLENGTHPARAM[aat][SUBBLENGTHTARGET[j]].mean, SUBLENGTHPARAM[aat][SUBBLENGTHTARGET[j]].lower, SUBLENGTHPARAM[aat][SUBBLENGTHTARGET[j]].upper, SUBLENGTHPARAM[aat][SUBBLENGTHTARGET[j]].slower, SUBLENGTHPARAM[aat][SUBBLENGTHTARGET[j]].supper, outflag);
					}
				}
				totalenergy+=(inew-iold);
			}else{
				totalenergy+=energy[i];
			}
		}
	}else if(uptype == 2){
		double iold, inew;	
		for(i=0;i<nodes.size();i++)
		{
			base=nodes[i];

			if(residueflag[base] || ((base-1 > 0) && (residueflag[base-1]))){
				inew=0.0;
				iold=energy[i];
				aat=restr[base].aaa;

				if(base == 0){
					for(j=0;j<AccAtomParamNum[aat][2];j++){
						if(SUBBLENGTHTARGET[j]==0 | SUBBLENGTHTARGET[j]==1 | SUBBLENGTHTARGET[j]==11) continue;
						else{
							if(SUBBLENGTHTARGET[j] == 8) continue;
							//cout<<restr[base].oplen[SUBBLENGTHTARGET[j]] << "-" << SUBLENGTHPARAM[aat][BLENGTHTARGET[j]].mean<<endl;
							inew+=getoutenelength(restr[base].oplen[SUBBLENGTHTARGET[j]], SUBLENGTHPARAM[aat][SUBBLENGTHTARGET[j]].mean, SUBLENGTHPARAM[aat][SUBBLENGTHTARGET[j]].lower, SUBLENGTHPARAM[aat][SUBBLENGTHTARGET[j]].upper, SUBLENGTHPARAM[aat][SUBBLENGTHTARGET[j]].slower, SUBLENGTHPARAM[aat][SUBBLENGTHTARGET[j]].supper, outflag);
						}
					}
				}else{
					for(j=0;j<AccAtomParamNum[aat][2];j++){
						if(SUBBLENGTHTARGET[j] == 8) continue;
						//cout<<restr[base].oplen[SUBBLENGTHTARGET[j]] << "-" << SUBLENGTHPARAM[aat][BLENGTHTARGET[j]].mean<<endl;
						inew+=getoutenelength(restr[base].oplen[SUBBLENGTHTARGET[j]], SUBLENGTHPARAM[aat][SUBBLENGTHTARGET[j]].mean, SUBLENGTHPARAM[aat][SUBBLENGTHTARGET[j]].lower, SUBLENGTHPARAM[aat][SUBBLENGTHTARGET[j]].upper, SUBLENGTHPARAM[aat][SUBBLENGTHTARGET[j]].slower, SUBLENGTHPARAM[aat][SUBBLENGTHTARGET[j]].supper, outflag);
					}
				}
				energy[i]+=(inew-iold);
				totalenergy+=energy[i];
			}else{
				totalenergy+=energy[i];
			}
		}
	}
	
	return totalenergy;
}


void RNASimulation::initenergygss(map<string, vector<double>> &edgeenergy, map<string, vector<double>> &nodeenergy, vector<vector<bool>> &lengthoutflag, vector<vector<bool>> &angleoutflag, vector<pair<int, int>> edges, vector<int> nodes, int numseq)
{
	nodeenergy["inittemp"].clear();
	nodeenergy["inittemp"].resize(nodes.size(), 0.0);
	nodeenergy["angle"].clear();
	nodeenergy["angle"].resize(nodes.size(), 0.0);
	nodeenergy["length"].clear();
	nodeenergy["length"].resize(nodes.size(), 0.0);
	nodeenergy["cnnlength"].clear();
	nodeenergy["cnnlength"].resize(nodes.size(), 0.0);
	lengthoutflag.resize(nodes.size(), vector<bool>(BLENGTHTARGET.size()));
	angleoutflag.resize(nodes.size(), vector<bool>(INNERANGLETARGET.size()));

	edgeenergy["bptemplength"].clear();
	edgeenergy["bptemplength"].resize(edges.size(), 0.0);
	edgeenergy["adjtemplength"].clear();
	edgeenergy["adjtemplength"].resize(edges.size(), 0.0);
}

double RNASimulation::calcenegss(vector<Residueatoms> restr, vector<Residueatoms> decstr, int numseq, vector<pair<int, int>> edgeinfo, vector<int> edgetype, vector<int> paraminfo, vector<int> nodes, map<int, vector<Residueatoms>> dnttemplates, map<int, vector<Residueatoms>> bptemplates, map<int, TemplateDistMap> dnttemplateslength, map<int, TemplateDistMap> bptemplateslength, map<string, double> enecoeff, vector<bool> residueflag, vector<vector<bool>> atomflag, map<string, vector<double>> &edgeenergy, map<string, vector<double>> &nodeenergy, vector<vector<bool>> &lengthoutflag, vector<vector<bool>> &angleoutflag, int uptype)
{	
	double totalene;
	totalene=0.0;
	double ene1, ene2, ene3, ene4, ene5, ene6;
	ene1=enecoeff["inittemp"]*calcinitstrgss(restr, decstr, nodes, nodeenergy["inittemp"], residueflag, uptype);
	totalene+=ene1;
	//cout<<"inittemp-" <<tempene<<endl;
	//ene2=enecoeff["adjtemplength"]*calcmseditmplengthgss(restr, edgeinfo, edgetype, dnttemplateslength, edgeenergy["adjtemplength"], residueflag, uptype);
	ene2=enecoeff["adjtemplength"]*calcmseditmpgss(restr, edgeinfo, edgetype, dnttemplates, edgeenergy["adjtemplength"], residueflag, uptype);
	//calcmseditmplengthgss(restr, edgeinfo, edgetype, dnttemplateslength, edgeenergy["adjtemplength"], residueflag, uptype);
	totalene+=ene2;
	//cout<<"adjtemplength-" <<tempene<<endl;

	ene3=enecoeff["bptemplength"]*calcmsebptmplengthgss(restr, edgeinfo, edgetype, bptemplateslength, edgeenergy["bptemplength"], residueflag, uptype);
	totalene+=ene3;

	ene4=enecoeff["angle"]*calcenebondanglegss(restr, nodes, nodeenergy["angle"], angleoutflag, residueflag, numseq, uptype);
	//cout<<"angle-" <<tempene<<endl;
	totalene+=ene4;

	ene5=enecoeff["length"]*calcenebondlengthgss(restr, nodes, nodeenergy["length"], lengthoutflag, residueflag, numseq, uptype);
	totalene+=ene5;

	ene6=enecoeff["cnnlength"]*calcenenbondlengthgss(restr, nodes, nodeenergy["cnnlength"], residueflag, uptype);
	totalene+=ene6;

	//cout<<" inittemp-" <<ene1 << " adjtemplength-" <<ene2 << " bptemplength-" << ene3 << " angle-" <<ene4 << " length-" << ene5 << " cnnlength-" << ene6 <<endl;

	return totalene;
}

void RNASimulation::getflaggss(vector<bool> &residueflag, vector<vector<bool>> &atomflag, int trandpos, char aat, int paramtype){
	int atom=ParamsList[aat][paramtype][3];
	atomflag[trandpos][atom]=true;
	residueflag[trandpos]=true;
	
}

void RNASimulation::changeparamgss(vector<Residueatoms> &recrrstr, int paramtype, int trandpos)
{
	char aat=recrrstr[trandpos].aaa;
	int atom1, atom2, atom3, atom4;
	Point3s pn;
	bool flagts;

	atom1=ParamsList[aat][paramtype][0];
	atom2=ParamsList[aat][paramtype][1];
	atom3=ParamsList[aat][paramtype][2];
	atom4=ParamsList[aat][paramtype][3];

	if(paramtype == 0){
		flagts=bf.tor2pos22(recrrstr[trandpos-1].ptv[atom1].x,recrrstr[trandpos-1].ptv[atom1].y,recrrstr[trandpos-1].ptv[atom1].z,
							recrrstr[trandpos-1].ptv[atom2].x,recrrstr[trandpos-1].ptv[atom2].y,recrrstr[trandpos-1].ptv[atom2].z,
							recrrstr[trandpos-1].ptv[atom3].x,recrrstr[trandpos-1].ptv[atom3].y,recrrstr[trandpos-1].ptv[atom3].z,
							recrrstr[trandpos].dih[paramtype]*RADDEG,recrrstr[trandpos].len[paramtype],recrrstr[trandpos].ang[paramtype]*RADDEG,
							&pn.x,&pn.y,&pn.z);
	}else if(paramtype == 1){
		flagts=bf.tor2pos22(recrrstr[trandpos-1].ptv[atom1].x,recrrstr[trandpos-1].ptv[atom1].y,recrrstr[trandpos-1].ptv[atom1].z,
							recrrstr[trandpos-1].ptv[atom2].x,recrrstr[trandpos-1].ptv[atom2].y,recrrstr[trandpos-1].ptv[atom2].z,
							recrrstr[trandpos].ptv[atom3].x,recrrstr[trandpos].ptv[atom3].y,recrrstr[trandpos].ptv[atom3].z,
							recrrstr[trandpos].dih[paramtype]*RADDEG,recrrstr[trandpos].len[paramtype],recrrstr[trandpos].ang[paramtype]*RADDEG,
							&pn.x,&pn.y,&pn.z);
	}else if(paramtype == 2){
		flagts=bf.tor2pos22(recrrstr[trandpos-1].ptv[atom1].x,recrrstr[trandpos-1].ptv[atom1].y,recrrstr[trandpos-1].ptv[atom1].z,
							recrrstr[trandpos].ptv[atom2].x,recrrstr[trandpos].ptv[atom2].y,recrrstr[trandpos].ptv[atom2].z,
							recrrstr[trandpos].ptv[atom3].x,recrrstr[trandpos].ptv[atom3].y,recrrstr[trandpos].ptv[atom3].z,
							recrrstr[trandpos].dih[paramtype]*RADDEG,recrrstr[trandpos].len[paramtype],recrrstr[trandpos].ang[paramtype]*RADDEG,
							&pn.x,&pn.y,&pn.z);
	}else if(paramtype == 10){
		flagts=bf.tor2pos22(recrrstr[trandpos-1].ptv[atom1].x,recrrstr[trandpos-1].ptv[atom1].y,recrrstr[trandpos-1].ptv[atom1].z,
							recrrstr[trandpos].ptv[atom2].x,recrrstr[trandpos].ptv[atom2].y,recrrstr[trandpos].ptv[atom2].z,
							recrrstr[trandpos].ptv[atom3].x,recrrstr[trandpos].ptv[atom3].y,recrrstr[trandpos].ptv[atom3].z,
							recrrstr[trandpos].dih[paramtype]*RADDEG,recrrstr[trandpos].len[paramtype],recrrstr[trandpos].ang[paramtype]*RADDEG,
							&pn.x,&pn.y,&pn.z);
	}else if(paramtype == 11){
		flagts=bf.tor2pos22(recrrstr[trandpos].ptv[atom1].x,recrrstr[trandpos].ptv[atom1].y,recrrstr[trandpos].ptv[atom1].z,
							recrrstr[trandpos-1].ptv[atom2].x,recrrstr[trandpos-1].ptv[atom2].y,recrrstr[trandpos-1].ptv[atom2].z,
							recrrstr[trandpos].ptv[atom3].x,recrrstr[trandpos].ptv[atom3].y,recrrstr[trandpos].ptv[atom3].z,
							recrrstr[trandpos].dih[paramtype]*RADDEG,recrrstr[trandpos].len[paramtype],recrrstr[trandpos].ang[paramtype]*RADDEG,
							&pn.x,&pn.y,&pn.z);
	}else{
		flagts=bf.tor2pos22(recrrstr[trandpos].ptv[atom1].x,recrrstr[trandpos].ptv[atom1].y,recrrstr[trandpos].ptv[atom1].z,
							recrrstr[trandpos].ptv[atom2].x,recrrstr[trandpos].ptv[atom2].y,recrrstr[trandpos].ptv[atom2].z,
							recrrstr[trandpos].ptv[atom3].x,recrrstr[trandpos].ptv[atom3].y,recrrstr[trandpos].ptv[atom3].z,
							recrrstr[trandpos].dih[paramtype]*RADDEG,recrrstr[trandpos].len[paramtype],recrrstr[trandpos].ang[paramtype]*RADDEG,
							&pn.x,&pn.y,&pn.z);
	}

	recrrstr[trandpos].ptv[atom4].x=pn.x;
	recrrstr[trandpos].ptv[atom4].y=pn.y;
	recrrstr[trandpos].ptv[atom4].z=pn.z;
	
}

void RNASimulation::updateparamgss(vector<Residueatoms> &recrrstr, vector<vector<bool>> atomflag, int trandpos, int numseq){
	int i,j;
	char aat;
	int atom1, atom2, atom3, atom4;
	int residx1, residx2, residx3, residx4;
	Point3s tmp1,tmp2,tmp3,tmp4;
	int iend;

	if(trandpos+2 > numseq) iend=numseq;
	else iend=trandpos+2;

	for(i=trandpos;i<iend;i++){
		aat=recrrstr[i].aaa;
		for(j=0;j<AccAtomParamNum[aat][2];j++)
		{
			if(j == 0){
				residx1=i-1; residx2=i-1; residx3=i-1; residx4=i;
			}else if(j == 1){
				residx1=i-1; residx2=i-1; residx3=i; residx4=i;
			}else if(j == 2){
				residx1=i-1; residx2=i; residx3=i; residx4=i;
			}else if(j == 10){
				residx1=i-1; residx2=i; residx3=i; residx4=i;
			}else if(j == 11){
				residx1=i; residx2=i-1; residx3=i; residx4=i;
			}else{
				residx1=i; residx2=i; residx3=i; residx4=i;
			}

			atom1=ParamsList[aat][j][0]; atom2=ParamsList[aat][j][1]; atom3=ParamsList[aat][j][2]; atom4=ParamsList[aat][j][3];
			tmp1=recrrstr[residx1].ptv[atom1];tmp2=recrrstr[residx2].ptv[atom2];tmp3=recrrstr[residx3].ptv[atom3];tmp4=recrrstr[residx4].ptv[atom4];
			
			if(atomflag[residx4][atom4] | atomflag[residx3][atom3]){
				coord2params(tmp1, tmp2, tmp3, tmp4, recrrstr[i].len[j], recrrstr[i].ang[j], recrrstr[i].dih[j]);
				coord2len(tmp2, tmp4, recrrstr[i].oplen[j]);
			}else if(atomflag[residx2][atom2]){
				coord2ang(tmp2, tmp3, tmp4, recrrstr[i].ang[j]);
				coord2dangle(tmp1, tmp2, tmp3, tmp4, recrrstr[i].dih[j]);
				coord2len(tmp2, tmp4, recrrstr[i].oplen[j]);
			}else if(atomflag[residx1][atom1]){
				coord2dangle(tmp1, tmp2, tmp3, tmp4, recrrstr[i].dih[j]);
			}
		}

		for(j=0;j<AddAngleParamsList[aat].size();j++)
		{
			if(j == 0){
				residx1=i-1; residx2=i; residx3=i;
			}else{
				residx1=i; residx2=i; residx3=i;
			}
			atom1=AddAngleParamsList[aat][j][0]; atom2=AddAngleParamsList[aat][j][1]; atom3=AddAngleParamsList[aat][j][2];
			tmp1=recrrstr[residx1].ptv[atom1]; tmp2=recrrstr[residx2].ptv[atom2]; tmp3=recrrstr[residx3].ptv[atom3];
			if(atomflag[residx1][atom1] | atomflag[residx2][atom2] | atomflag[residx3][atom3]){
				coord2ang(tmp1, tmp2, tmp3, recrrstr[i].addang[j]);
			}
		}

		for(j=0;j<AddLengthParamsList[aat].size();j++)
		{
			residx1=i; residx2=i; residx3=i;
			atom1=AddLengthParamsList[aat][j][0]; atom2=AddLengthParamsList[aat][j][1];
			tmp1=recrrstr[i].ptv[atom1]; tmp2=recrrstr[i].ptv[atom2];
			if(atomflag[residx1][atom1] | atomflag[residx2][atom2]){
				coord2len(tmp1, tmp2, recrrstr[i].addlen[j]);
			}
		}
	}
	
}

double RNASimulation::goldenSearchcalcene(vector<Residueatoms> recrrstr, vector<Residueatoms> decstr, double torval, int paramtype, int paramidx, int trandpos, int numseq, vector<pair<int, int>> edgeinfo, vector<int> edgetype, vector<int> paraminfo, vector<int> nodes, map<int, vector<Residueatoms>> dnttemplates, map<int, vector<Residueatoms>> bptemplates, map<int, TemplateDistMap> dnttemplateslength, map<int, TemplateDistMap> bptemplateslength, map<string, double> enecoeff, vector<bool> residueflag, vector<vector<bool>> atomflag, map<string, vector<double>> edgeenergy, map<string, vector<double>> nodeenergy, vector<vector<bool>> lengthoutflag, vector<vector<bool>> angleoutflag){
	bool flagtor;
	double energy;

	if(paramtype == 0){
		recrrstr[trandpos].dih[paramidx]=torval;
	}else if(paramtype == 1){
		recrrstr[trandpos].ang[paramidx]=torval;
	}else if(paramtype == 2){
		recrrstr[trandpos].len[paramidx]=torval;
	}

	changeparamgss(recrrstr, paramidx, trandpos);

	updateparamgss(recrrstr, atomflag, trandpos, numseq);

	energy=calcenegss(recrrstr, decstr, numseq, edgeinfo, edgetype, paraminfo, nodes, dnttemplates, bptemplates, dnttemplateslength, bptemplateslength, enecoeff, residueflag, atomflag, edgeenergy, nodeenergy, lengthoutflag, angleoutflag, 1);

	return energy;
}

const double PHI = (1 + std::sqrt(5)) / 2;
double RNASimulation::goldenSearchdih(vector<Residueatoms> &recrrstr, vector<Residueatoms> decstr, int trandpos, int dihtype, double tolerance, int numseq, vector<pair<int, int>> edgeinfo, vector<int> edgetype, vector<int> paraminfo, vector<int> nodes, map<int, vector<Residueatoms>> dnttemplates, map<int, vector<Residueatoms>> bptemplates, map<int, TemplateDistMap> dnttemplateslength, map<int, TemplateDistMap> bptemplateslength, map<string, double> enecoeff, map<string, vector<double>> &edgeenergy, map<string, vector<double>> &nodeenergy, vector<vector<bool>> &lengthoutflag, vector<vector<bool>> &angleoutflag){
	char aat=recrrstr[trandpos].aaa;
	double a=TORSIONPARAM[aat][dihtype].lower; 
	double b=TORSIONPARAM[aat][dihtype].upper;
	double energy;
	
	double c = b - (b - a) / PHI;
    double d = a + (b - a) / PHI;
	double cene, dene;
	int count=0;

	vector<bool> residueflag;
	vector<vector<bool>> atomflag;
	resetflaggss(residueflag, atomflag, numseq);
	getflaggss(residueflag, atomflag, trandpos, aat, dihtype);

	while (std::abs(b - a) > tolerance) {
		cene=goldenSearchcalcene(recrrstr, decstr, c, 0, dihtype, trandpos, numseq, edgeinfo, edgetype, paraminfo, nodes, dnttemplates, bptemplates, dnttemplateslength, bptemplateslength, enecoeff,residueflag, atomflag, edgeenergy, nodeenergy, lengthoutflag, angleoutflag);
		dene=goldenSearchcalcene(recrrstr, decstr, d, 0, dihtype, trandpos, numseq, edgeinfo, edgetype, paraminfo, nodes, dnttemplates, bptemplates, dnttemplateslength, bptemplateslength, enecoeff,residueflag, atomflag, edgeenergy, nodeenergy, lengthoutflag, angleoutflag);

		if (cene < dene) {
            b = d;
        } else {
            a = c;
        }

		c = b - (b - a) / PHI;
        d = a + (b - a) / PHI;
		count+=1;
	}

	//cout<<dihtype << "-" << a << "-" << b << "-" << count<<endl;

	recrrstr[trandpos].dih[dihtype]=(a + b)/2;

	changeparamgss(recrrstr, dihtype, trandpos);
	updateparamgss(recrrstr, atomflag, trandpos, numseq);
	energy=calcenegss(recrrstr, decstr, numseq, edgeinfo, edgetype, paraminfo, nodes, dnttemplates, bptemplates, dnttemplateslength, bptemplateslength, enecoeff, residueflag, atomflag, edgeenergy, nodeenergy, lengthoutflag, angleoutflag, 2);

	double energy2=calcenegss(recrrstr, decstr, numseq, edgeinfo, edgetype, paraminfo, nodes, dnttemplates, bptemplates, dnttemplateslength, bptemplateslength, enecoeff, residueflag, atomflag, edgeenergy, nodeenergy, lengthoutflag, angleoutflag, 0);

	//cout<<energy << "-" << energy2<<endl;

	//if(abs(energy-energy2) > 0.0001){
		//cout<<"error"<<endl;
		//exit(1);
	//}

	return energy;
}

double RNASimulation::goldenSearchang(vector<Residueatoms> &recrrstr, vector<Residueatoms> decstr, int trandpos, int angtype, double tolerance, int numseq, vector<pair<int, int>> edgeinfo, vector<int> edgetype, vector<int> paraminfo, vector<int> nodes, map<int, vector<Residueatoms>> dnttemplates, map<int, vector<Residueatoms>> bptemplates, map<int, TemplateDistMap> dnttemplateslength, map<int, TemplateDistMap> bptemplateslength, map<string, double> enecoeff, map<string, vector<double>> &edgeenergy, map<string, vector<double>> &nodeenergy, vector<vector<bool>> &lengthoutflag, vector<vector<bool>> &angleoutflag){
	char aat=recrrstr[trandpos].aaa;
	double a=INNERANGPARAM[aat][angtype].lower; 
	double b=INNERANGPARAM[aat][angtype].upper;
	double energy;
	
	double c = b - (b - a) / PHI;
    double d = a + (b - a) / PHI;
	double cene, dene;
	int count=0;

	vector<bool> residueflag;
	vector<vector<bool>> atomflag;
	resetflaggss(residueflag, atomflag, numseq);
	getflaggss(residueflag, atomflag, trandpos, aat, angtype);

	while (std::abs(b - a) > tolerance) {
		cene=goldenSearchcalcene(recrrstr, decstr, c, 1, angtype, trandpos, numseq, edgeinfo, edgetype, paraminfo, nodes, dnttemplates, bptemplates, dnttemplateslength, bptemplateslength, enecoeff,residueflag, atomflag, edgeenergy, nodeenergy, lengthoutflag, angleoutflag);
		dene=goldenSearchcalcene(recrrstr, decstr, d, 1, angtype, trandpos, numseq, edgeinfo, edgetype, paraminfo, nodes, dnttemplates, bptemplates, dnttemplateslength, bptemplateslength, enecoeff,residueflag, atomflag, edgeenergy, nodeenergy, lengthoutflag, angleoutflag);

		if (cene < dene) {
            b = d;
        } else {
            a = c;
        }

		c = b - (b - a) / PHI;
        d = a + (b - a) / PHI;
		count+=1;
	}

	recrrstr[trandpos].ang[angtype]=(a + b)/2;
	changeparamgss(recrrstr, angtype, trandpos);
	updateparamgss(recrrstr, atomflag, trandpos, numseq);
	energy=calcenegss(recrrstr, decstr, numseq, edgeinfo, edgetype, paraminfo, nodes, dnttemplates, bptemplates, dnttemplateslength, bptemplateslength, enecoeff, residueflag, atomflag, edgeenergy, nodeenergy, lengthoutflag, angleoutflag, 2);

	double energy2=calcenegss(recrrstr, decstr, numseq, edgeinfo, edgetype, paraminfo, nodes, dnttemplates, bptemplates, dnttemplateslength, bptemplateslength, enecoeff, residueflag, atomflag, edgeenergy, nodeenergy, lengthoutflag, angleoutflag, 0);

	//cout<<energy << "-" << energy2<<endl;

	//if(abs(energy-energy2) > 0.0001){
		//cout<<"error"<<endl;
		//exit(1);
	//}

	return energy;
}

double RNASimulation::goldenSearchlen(vector<Residueatoms> &recrrstr, vector<Residueatoms> decstr, int trandpos, int lentype, double tolerance, int numseq, vector<pair<int, int>> edgeinfo, vector<int> edgetype, vector<int> paraminfo, vector<int> nodes, map<int, vector<Residueatoms>> dnttemplates, map<int, vector<Residueatoms>> bptemplates, map<int, TemplateDistMap> dnttemplateslength, map<int, TemplateDistMap> bptemplateslength, map<string, double> enecoeff, map<string, vector<double>> &edgeenergy, map<string, vector<double>> &nodeenergy, vector<vector<bool>> &lengthoutflag, vector<vector<bool>> &angleoutflag){
	char aat=recrrstr[trandpos].aaa;
	double a=LENGTHPARAM[aat][lentype].lower; 
	double b=LENGTHPARAM[aat][lentype].upper;
	double energy;
	
	double c = b - (b - a) / PHI;
    double d = a + (b - a) / PHI;
	double cene, dene;
	int count=0;

	vector<bool> residueflag;
	vector<vector<bool>> atomflag;
	resetflaggss(residueflag, atomflag, numseq);
	getflaggss(residueflag, atomflag, trandpos, aat, lentype);

	while (std::abs(b - a) > tolerance) {
		cene=goldenSearchcalcene(recrrstr, decstr, c, 2, lentype, trandpos, numseq, edgeinfo, edgetype, paraminfo, nodes, dnttemplates, bptemplates, dnttemplateslength, bptemplateslength, enecoeff,residueflag, atomflag, edgeenergy, nodeenergy, lengthoutflag, angleoutflag);
		dene=goldenSearchcalcene(recrrstr, decstr, d, 2, lentype, trandpos, numseq, edgeinfo, edgetype, paraminfo, nodes, dnttemplates, bptemplates, dnttemplateslength, bptemplateslength, enecoeff,residueflag, atomflag, edgeenergy, nodeenergy, lengthoutflag, angleoutflag);

		if (cene < dene) {
            b = d;
        } else {
            a = c;
        }

		c = b - (b - a) / PHI;
        d = a + (b - a) / PHI;
		count+=1;
	}

	recrrstr[trandpos].len[lentype]=(a + b)/2;
	changeparamgss(recrrstr, lentype, trandpos);
	updateparamgss(recrrstr, atomflag, trandpos, numseq);
	energy=calcenegss(recrrstr, decstr, numseq, edgeinfo, edgetype, paraminfo, nodes, dnttemplates, bptemplates, dnttemplateslength, bptemplateslength, enecoeff, residueflag, atomflag, edgeenergy, nodeenergy, lengthoutflag, angleoutflag, 2);

	double energy2=calcenegss(recrrstr, decstr, numseq, edgeinfo, edgetype, paraminfo, nodes, dnttemplates, bptemplates, dnttemplateslength, bptemplateslength, enecoeff, residueflag, atomflag, edgeenergy, nodeenergy, lengthoutflag, angleoutflag, 0);

	//cout<<energy << "-" << energy2<<endl;

	//if(abs(energy-energy2) > 0.0001){
		//cout<<"error"<<endl;
		//exit(1);
	//}

	return energy;
}





