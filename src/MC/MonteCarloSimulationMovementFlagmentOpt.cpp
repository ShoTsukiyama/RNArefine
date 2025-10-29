
#include "MonteCarloSimulation.hpp"

void MonteCarloSimulation::setadjfragflap(int istart, int iend)
{
	int i,j;
	int atom1, atom2, atom3;

	for(i=istart;i<iend;i++){
		char base=refinestr_temp[i].aaa;

        for(j=0;j<MCLENGTHPARAM[base].size();j++){
            lengthflag[i][j] = 1;
        }

        for(j=0;j<MCINNERANGPARAM[base].size();j++){
            angleflag[i][j] = 1;
        }

        for(j=0;j<MCTORSIONPARAM[base].size();j++){
            torsionflag[i][j] = 1;
        }

        for(j=0;j<MCHYDLENGTHPARAM[base].size();j++){
            hydlengthflag[i][j] = 1;
        }

        for(j=0;j<MCHYDINNERANGPARAM[base].size();j++){
            hydangleflag[i][j] = 1;
        }

        for(j=0;j<MCHYDTORSIONPARAM[base].size();j++){
            hydtorsionflag[i][j] = 1;
        }

        for(j=0;j<AddLengthParamsList[base].size();j++){
            addlengthflag[i][j] = 1;
        }

        for(j=0;j<AddAngleParamsList[base].size();j++){
            addangleflag[i][j] = 1;
        }

		for(j=0;j<AddHydLengthParamsList[base].size();j++){
			hydaddlengthflag[i][j]=1;
		}

		for(j=0;j<AddHydAngleParamsList[base].size();j++)
		{
			hydaddangleflag[i][j]=1;
		}
	}

	for(i=0;i<numbb;i++)
	{
		ntflag[i] = 1;
		for(j=i;j<numbb;j++)
		{
			if(j<istart || i>iend){}
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

void MonteCarloSimulation::setnonadjfragflap(int istart1, int iend1, int istart2, int iend2)
{
	int i,j;
	int atom1, atom2, atom3;

	for(i=istart1;i<iend1;i++){
		char base=refinestr_temp[i].aaa;

        for(j=0;j<MCLENGTHPARAM[base].size();j++){
            lengthflag[i][j] = 1;
        }

        for(j=0;j<MCINNERANGPARAM[base].size();j++){
            angleflag[i][j] = 1;
        }

        for(j=0;j<MCTORSIONPARAM[base].size();j++){
            torsionflag[i][j] = 1;
        }

        for(j=0;j<MCHYDLENGTHPARAM[base].size();j++){
            hydlengthflag[i][j] = 1;
        }

        for(j=0;j<MCHYDINNERANGPARAM[base].size();j++){
            hydangleflag[i][j] = 1;
        }

        for(j=0;j<MCHYDTORSIONPARAM[base].size();j++){
            hydtorsionflag[i][j] = 1;
        }

        for(j=0;j<AddLengthParamsList[base].size();j++){
            addlengthflag[i][j] = 1;
        }

        for(j=0;j<AddAngleParamsList[base].size();j++){
            addangleflag[i][j] = 1;
        }

		for(j=0;j<AddHydLengthParamsList[base].size();j++){
			hydaddlengthflag[i][j]=1;
		}

		for(j=0;j<AddHydAngleParamsList[base].size();j++)
		{
			hydaddangleflag[i][j]=1;
		}
	}

	for(i=istart2;i<iend2;i++){
		char base=refinestr_temp[i].aaa;

        for(j=0;j<MCLENGTHPARAM[base].size();j++){
            lengthflag[i][j] = 1;
        }

        for(j=0;j<MCINNERANGPARAM[base].size();j++){
            angleflag[i][j] = 1;
        }

        for(j=0;j<MCTORSIONPARAM[base].size();j++){
            torsionflag[i][j] = 1;
        }

        for(j=0;j<MCHYDLENGTHPARAM[base].size();j++){
            hydlengthflag[i][j] = 1;
        }

        for(j=0;j<MCHYDINNERANGPARAM[base].size();j++){
            hydangleflag[i][j] = 1;
        }

        for(j=0;j<MCHYDTORSIONPARAM[base].size();j++){
            hydtorsionflag[i][j] = 1;
        }

        for(j=0;j<AddLengthParamsList[base].size();j++){
            addlengthflag[i][j] = 1;
        }

        for(j=0;j<AddAngleParamsList[base].size();j++){
            addangleflag[i][j] = 1;
        }

		for(j=0;j<AddHydLengthParamsList[base].size();j++){
			hydaddlengthflag[i][j]=1;
		}

		for(j=0;j<AddHydAngleParamsList[base].size();j++)
		{
			hydaddangleflag[i][j]=1;
		}
	}

	int istart = min(istart1, istart2);
	int iend = max(iend1, iend2);

	for(i=0;i<numbb;i++)
	{
		ntflag[i] = 1;
		for(j=i;j<numbb;j++)
		{
			if(j<istart || i>iend){}
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

void MonteCarloSimulation::setfragflap(int istart, int iend)
{
	int i,j;
	int atom1, atom2, atom3;

	for(i=istart;i<iend;i++){
		char base=refinestr_temp[i].aaa;

        for(j=0;j<MCLENGTHPARAM[base].size();j++){
            lengthflag[i][j] = 1;
        }

        for(j=0;j<MCINNERANGPARAM[base].size();j++){
            angleflag[i][j] = 1;
        }

        for(j=0;j<MCTORSIONPARAM[base].size();j++){
            torsionflag[i][j] = 1;
        }

        for(j=0;j<MCHYDLENGTHPARAM[base].size();j++){
            hydlengthflag[i][j] = 1;
        }

        for(j=0;j<MCHYDINNERANGPARAM[base].size();j++){
            hydangleflag[i][j] = 1;
        }

        for(j=0;j<MCHYDTORSIONPARAM[base].size();j++){
            hydtorsionflag[i][j] = 1;
        }

        for(j=0;j<AddLengthParamsList[base].size();j++){
            addlengthflag[i][j] = 1;
        }

        for(j=0;j<AddAngleParamsList[base].size();j++){
            addangleflag[i][j] = 1;
        }

		for(j=0;j<AddHydLengthParamsList[base].size();j++){
			hydaddlengthflag[i][j]=1;
		}

		for(j=0;j<AddHydAngleParamsList[base].size();j++)
		{
			hydaddangleflag[i][j]=1;
		}
	}

	ntflag[istart] = 1;
	for(j=0;j<numbb;j++)
	{
		mapbbflag[istart][j]=1;
		mapbbflag[j][istart]=1;
		mapflag[istart][j]=1;
		mapflag[j][istart]=1;
	}
}



int MonteCarloSimulation::selecttargetfrag(const int &fraglen)
{
	int i,j;
	return int(Random()*(numbb - fraglen+1));
}

bool MonteCarloSimulation::fmovementfrag(double oldenergy, double &newenergy)
{
	int i, j, l;
	bool flagtor;
	int trandpos, trandposst, trandposed,span;
	double mean, std, upper, lower;
	int idx;
	double rmax, rnorm;
	int mtype1, mtype2;
	char base, pbase;
	bool flagdone;
	bool exitOuterLoop = false;
	int leng;
	Point3s coordA, coordB;

	int minflaglen=2;
	int maxflaglen=min(5, numbb);

	int fraglen=Random()*(maxflaglen-minflaglen+1)+minflaglen;

	if(fraglen > numbb) fraglen=numbb;

	StrCopy(refinestr, refinestr_temp);
	resetfrag(0);

	int stpos = selecttargetfrag(fraglen);
	int edpos = stpos+fraglen;

	map<int, vector<Point3s>> initstructure;
	vector<vector<int>> bindingmap(numbb, vector<int>(numbb, 0));
	vector<bool> moveresidues(numbb,false);
	vector<bool> updateresidues(numbb,false);
	vector<vector<int>> basepairlist;
	vector<vector<int>> noncbasepairlist;
	vector<double> basepairscorelist;
	vector<double> noncbasepairscorelist;

	vector<vector<int>> contactmap(numbb, vector<int>(numbb, 0));
	vector<vector<double>> stackingmap(numbb, vector<double>(numbb, 0));
	vector<vector<int>> basepairmap(numbb, vector<int>(numbb, 0));

	for(i=stpos;i<edpos;i++){
		base = refinestr_temp[i].aaa;

		initstructure[i].resize(AccAtomParamNum[base][3]);
		for(j=0; j<AccAtomParamNum[base][3];j++){
			initstructure[i][j].x=refinestr_temp[i].ptv[j].x;
			initstructure[i][j].y=refinestr_temp[i].ptv[j].y;
			initstructure[i][j].z=refinestr_temp[i].ptv[j].z;
		}

		if(i - 1 >= 0) bindingmap[i - 1][i] = 1;
		if(i + 1 < numbb) bindingmap[i][i + 1] = 1;
		bindingmap[i][i] = 2;
		moveresidues[i]=true;
		updateresidues[i]=true;
		if(i+1 < numbb) updateresidues[i+1]=true;

		coordA.x = refinestr_temp[i].ptv[8].x;
		coordA.y = refinestr_temp[i].ptv[8].y;
		coordA.z = refinestr_temp[i].ptv[8].z;
		
		vector<int> noncbp_ks;
		vector<int> bp_ks;
		int noncbp_k = -1;
		int bp_k = -1;
		for(int j=0; j<numbb;j++)
		{
			Point3s rp;
			//if(i == j) continue;
			rp.x = refinestr_temp[j].ptv[8].x - coordA.x;
			rp.y = refinestr_temp[j].ptv[8].y - coordA.y;
			rp.z = refinestr_temp[j].ptv[8].z - coordA.z;

			leng = rp.x*rp.x+rp.y*rp.y+rp.z*rp.z;

			//if(leng < GLOBALCUTOFF_SQR){
			contactmap[i][j] = 1;
			contactmap[j][i] = 1;
			//}

			if(i == j){
				if(leng < CUTOFF_SQR+5.0){
					contactmap[i][j] = 2;
					contactmap[j][i] = 2;
				}
			}else{
				if(leng < CUTOFF_SQR+10.0){
					contactmap[i][j] = 2;
					contactmap[j][i] = 2;
				}
			}

			stackingmap[i][j] = stmap[i][j];
			stackingmap[j][i] = stmap[j][i];

			if(noncbpmap[i][j] > 0){
				vector<int> temp{min(i, j), max(i, j)};

				if (find(noncbasepairlist.begin(), noncbasepairlist.end(), temp) == noncbasepairlist.end()) {
					noncbasepairlist.push_back({i, j});
					noncbasepairscorelist.push_back(noncbpmap[i][j]);
					//noncbp_k = j;
					noncbp_ks.push_back(j);
				}
			}

			if(cbpmap[i][j] > 0){
				vector<int> temp{min(i, j), max(i, j) };

				if (find(basepairlist.begin(), basepairlist.end(), temp) == basepairlist.end()) {
					basepairlist.push_back({i, j});
					basepairscorelist.push_back(cbpmap[i][j]);
					//bp_k = j;
					bp_ks.push_back(j);
				}

			}
		}

		if(bp_ks.size() > 0)
		{
			for(j=0;j<bp_ks.size();j++)
			{
				bp_k=bp_ks[j];
				pbase=refinestr_temp[bp_k].aaa;
				initstructure[bp_k].resize(AccAtomParamNum[pbase][3]);
				for(l=0; l<AccAtomParamNum[pbase][3];l++){
					initstructure[bp_k][l].x=refinestr_temp[bp_k].ptv[l].x;
					initstructure[bp_k][l].y=refinestr_temp[bp_k].ptv[l].y;
					initstructure[bp_k][l].z=refinestr_temp[bp_k].ptv[l].z;
				}

				if(bp_k - 1 >= 0) bindingmap[bp_k - 1][bp_k] = 1;
				if(bp_k + 1 < numbb) bindingmap[bp_k][bp_k + 1] = 1;
				bindingmap[bp_k][bp_k] = 2;
				moveresidues[bp_k]=true;
				updateresidues[bp_k]=true;
				if(bp_k+1 < numbb) updateresidues[bp_k+1]=true;

				coordA.x = refinestr_temp[bp_k].ptv[8].x;
				coordA.y = refinestr_temp[bp_k].ptv[8].y;
				coordA.z = refinestr_temp[bp_k].ptv[8].z;
				for(int j=0; j<numbb;j++){
					//if(k == j) continue;
					Point3s rp;
					rp.x = refinestr_temp[j].ptv[8].x - coordA.x;
					rp.y = refinestr_temp[j].ptv[8].y - coordA.y;
					rp.z = refinestr_temp[j].ptv[8].z - coordA.z;

					leng = rp.x*rp.x+rp.y*rp.y+rp.z*rp.z;

					//if(leng < GLOBALCUTOFF_SQR){
					contactmap[bp_k][j] = 1;
					contactmap[j][bp_k] = 1;
					//}
					
					if(bp_k == j){
						if(leng < CUTOFF_SQR+5.0){
							contactmap[bp_k][j] = 2;
							contactmap[j][bp_k] = 2;
						}
					}else{
						if(leng < CUTOFF_SQR+10.0){
							contactmap[bp_k][j] = 2;
							contactmap[j][bp_k] = 2;
						}
					}

					stackingmap[bp_k][j] = stmap[bp_k][j];
					stackingmap[j][bp_k] = stmap[j][bp_k];
				}
			}
		}

		if(noncbp_ks.size() > 0)
		{
			for(j=0;j<noncbp_ks.size();j++)
			{
				noncbp_k=noncbp_ks[j];
				pbase=refinestr_temp[noncbp_k].aaa;
				initstructure[noncbp_k].resize(AccAtomParamNum[pbase][3]);
				for(l=0; l<AccAtomParamNum[pbase][3];l++){
					initstructure[noncbp_k][l].x=refinestr_temp[noncbp_k].ptv[l].x;
					initstructure[noncbp_k][l].y=refinestr_temp[noncbp_k].ptv[l].y;
					initstructure[noncbp_k][l].z=refinestr_temp[noncbp_k].ptv[l].z;
				}

				if(noncbp_k - 1 >= 0) bindingmap[noncbp_k - 1][noncbp_k] = 1;
				if(noncbp_k + 1 < numbb) bindingmap[noncbp_k][noncbp_k + 1] = 1;
				bindingmap[noncbp_k][noncbp_k] = 2;
				moveresidues[noncbp_k]=true;
				updateresidues[noncbp_k]=true;
				if(noncbp_k+1 < numbb) updateresidues[noncbp_k+1]=true;

				coordA.x = refinestr_temp[noncbp_k].ptv[8].x;
				coordA.y = refinestr_temp[noncbp_k].ptv[8].y;
				coordA.z = refinestr_temp[noncbp_k].ptv[8].z;
				for(int j=0; j<numbb;j++){
					//if(k == j) continue;
					Point3s rp;
					rp.x = refinestr_temp[j].ptv[8].x - coordA.x;
					rp.y = refinestr_temp[j].ptv[8].y - coordA.y;
					rp.z = refinestr_temp[j].ptv[8].z - coordA.z;

					leng = rp.x*rp.x+rp.y*rp.y+rp.z*rp.z;

					//if(leng < GLOBALCUTOFF_SQR){
					contactmap[noncbp_k][j] = 1;
					contactmap[j][noncbp_k] = 1;
					//}
					
					if(noncbp_k == j){
						if(leng < CUTOFF_SQR+5.0){
							contactmap[noncbp_k][j] = 2;
							contactmap[j][noncbp_k] = 2;
						}
					}else{
						if(leng < CUTOFF_SQR+10.0){
							contactmap[noncbp_k][j] = 2;
							contactmap[j][noncbp_k] = 2;
						}
					}

					stackingmap[noncbp_k][j] = stmap[noncbp_k][j];
					stackingmap[j][noncbp_k] = stmap[j][noncbp_k];
				}
			}
		}

		/*
		vector<vector<int>> ks;
		for(j=0;j<stacklist.size();j++){
			if(i == stacklist[j][0]){
				k=stacklist[j][1];
				ks.push_back({i,k});
			}
			else if(i == stacklist[j][1]){
				k=stacklist[j][0];
				ks.push_back({i,k});
			}
		}

		for(int l=0;l<ks.size();l++){
			k=ks[l][1];

			pbase=refinestr_temp[k].aaa;
			initstructure[k].resize(AccAtomParamNum[pbase][3]);
			for(int m=0; m<AccAtomParamNum[pbase][3];m++){
				initstructure[k][m].x=refinestr_temp[k].ptv[m].x;
				initstructure[k][m].y=refinestr_temp[k].ptv[m].y;
				initstructure[k][m].z=refinestr_temp[k].ptv[m].z;
			}

			if(k - 1 >= 0) bindingmap[k - 1][k] = 1;
			if(k + 1 < numbb) bindingmap[k][k + 1] = 1;
			bindingmap[k][k] = 2;
			moveresidues[k]=true;
			updateresidues[k]=true;
			if(k+1 < numbb) updateresidues[k+1]=true;

			coordA.x = refinestr_temp[k].ptv[8].x;
			coordA.y = refinestr_temp[k].ptv[8].y;
			coordA.z = refinestr_temp[k].ptv[8].z;
			for(int j=0; j<numbb;j++){
				//if(k == j) continue;
				Point3s rp;
				rp.x = refinestr_temp[j].ptv[8].x - coordA.x;
				rp.y = refinestr_temp[j].ptv[8].y - coordA.y;
				rp.z = refinestr_temp[j].ptv[8].z - coordA.z;

				leng = rp.x*rp.x+rp.y*rp.y+rp.z*rp.z;

				if(leng < GLOBALCUTOFF_SQR){
					contactmap[k][j] = 1;
					contactmap[j][k] = 1;
				}
				
				if(k == j){
					if(leng < CUTOFF_SQR+5.0){
						contactmap[k][j] = 2;
						contactmap[j][k] = 2;
					}
				}else{
					if(leng < CUTOFF_SQR+10.0){
						contactmap[k][j] = 2;
						contactmap[j][k] = 2;
					}
				}

				stackingmap[k][j] = stmap[k][j];
				stackingmap[j][k] = stmap[j][k];
			}
		}
		*/
	}

	//cout<<"=============="<<endl;
	//cout<<target_stack_idx << "-" << stack_base1 << "-" << stack_base2 << "-" << stacktype<<endl;

	MonteCarloSimulation::MCfragmentoptimization(initstructure, bindingmap, moveresidues, contactmap, stackingmap, basepairlist, noncbasepairlist, basepairscorelist, noncbasepairscorelist);


	int istart=stpos-1;
	int iend=edpos+1;
	if(istart<0) istart=0;
	if(iend>=numbb) iend=numbb;

	for(int k=0;k<numbb;k++){
		if(updateresidues[k]){
			coord2parambb(refinestr_temp,k,k+1);
			coord2paramsugar(refinestr_temp,k,k+1);
			coord2parambase(refinestr_temp,k,k+1);
			coord2paramhydrogen(refinestr_temp,k,k+1);
			coord2addparamhydrogen(refinestr_temp,k,k+1);
			coord2paramadd(refinestr_temp,k,k+1);
			//setadjstackfragflap(k,k+1);
			setfragflap(k,k+1);
		}
	}

	/*
	coord2parambb(refinestr_temp,istart,iend);
	coord2paramsugar(refinestr_temp,istart,iend);
	coord2parambase(refinestr_temp,istart,iend);
	coord2paramhydrogen(refinestr_temp,istart,iend);
	coord2addparamhydrogen(refinestr_temp,istart,iend);
	coord2paramadd(refinestr_temp,istart,iend);
	setadjstackfragflap(istart,iend);

	vector<int> processnt; 
	for(i=0;i<basepairlist.size();i++){
		k=basepairlist[i][1];
		coord2parambb(refinestr_temp,k,k+1);
		coord2paramsugar(refinestr_temp,k,k+1);
		coord2parambase(refinestr_temp,k,k+1);
		coord2paramhydrogen(refinestr_temp,k,k+1);
		coord2addparamhydrogen(refinestr_temp,k,k+1);
		coord2paramadd(refinestr_temp,k,k+1);
		setadjstackfragflap(k,k+1);
		processnt.push_back(k);


		if(k-1 >= 0){
			auto it = std::find(processnt.begin(), processnt.end(), k-1);
			if (it == processnt.end()) {
				j=k-1;
				coord2parambb(refinestr_temp,j,j+1);
				coord2paramsugar(refinestr_temp,j,j+1);
				coord2parambase(refinestr_temp,j,j+1);
				coord2paramhydrogen(refinestr_temp,j,j+1);
				coord2addparamhydrogen(refinestr_temp,j,j+1);
				coord2paramadd(refinestr_temp,j,j+1);
				setadjstackfragflap(j,j+1);
				processnt.push_back(j);
			}
		}

		if(k+1 < numbb){
			auto it = std::find(processnt.begin(), processnt.end(), k+1);
			if (it == processnt.end()) {
				j=k+1;
				coord2parambb(refinestr_temp,j,j+1);
				coord2paramsugar(refinestr_temp,j,j+1);
				coord2parambase(refinestr_temp,j,j+1);
				coord2paramhydrogen(refinestr_temp,j,j+1);
				coord2addparamhydrogen(refinestr_temp,j,j+1);
				coord2paramadd(refinestr_temp,j,j+1);
				setadjstackfragflap(j,j+1);
				processnt.push_back(j);
			}
		}
	}
	*/
	

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

