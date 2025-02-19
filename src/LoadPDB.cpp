// ParsePDB.cpp: implementation of the ParsePDB class.
//
//////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>
#include <stdlib.h>
#include <sstream>
#include "LoadPDB.hpp"

LoadPDB::LoadPDB()
{
	allocpromod=1;
	allocproseq=800;
	flagter=false;
}


int LoadPDB::dealline(string st,vector<Atom> &proseq, vector<Promodel> &promod)
{
	//cout<<st<<endl;
	
	int i,j;

	if(st.substr(0,5) == "MODEL"){

		int tempnumi;
		istringstream iss(st);
		vector<std::string> tokens;
		std::string token;
		while (iss >> token) {
			tokens.push_back(token);
		}
		//cout<<tokens[1]<<endl;
		//exit(1);
		//cout<<st.substr(10+i,1).c_str()<<endl;
		//tempnumi=atoi(st.substr(10+i,1).c_str());

		//tempnumi=atoi(tokens[1].c_str()[0]);
		if(numpromod==0 && flagfirstmodel)
		{
			flagfirstmodel=false;
			return 42;
		}
		
		if(numpromod==allocpromod-1)
		{
			allocpromod*=2;
		    promod.resize(allocpromod);
			for(i=numpromod+1;i<allocpromod;i++)
			{
				promod[i].nchain=0;
				for(int j=0;j<MAX_CHAIN_NUM;j++)
				{
					promod[i].procha[j].chainseg.init=-1;
					promod[i].procha[j].chainseg.term=-1;
				}
			}
		}
		// end of model
		if(promod[numpromod].procha[promod[numpromod].nchain].chainseg.term==-1)
		{
			promod[numpromod].procha[promod[numpromod].nchain].chainseg.term=numproseq-1;
		}
		flagbeginmodel=true;
		// only here change numpromod
		numpromod++;
		return 14;

	}
	
	else if(st.substr(0,6) == "ENDMDL")
	{
		// end of model
		if(promod[numpromod].procha[promod[numpromod].nchain].chainseg.term==-1)
		{
			promod[numpromod].procha[promod[numpromod].nchain].chainseg.term=numproseq-1;
		}
		return 15;
	}
	else if (st.substr(0,3) == "TER")
	{
		//end a chain
		flagter=true;
		return 16;
	}
	//else if(st.substr(0,5)=="ATOM " || st.substr(0,6)=="HETATM")
	else if(st.substr(0,4)=="ATOM" || st.substr(0,6)=="HETATM")
	{
		if(st.substr(16,1)[0]!=' ' && st.substr(16,1)[0]!='A' && numproseq>0)
		{
			if(proseq[numproseq-1].detailtype==st.substr(12,4))
			{
				return 41;		
			}	
		}
		
		if(numproseq==0)
		{
			proseq.resize(allocproseq);
		}
		
		else if(numproseq==allocproseq-1)
		{
			allocproseq*=2;
			proseq.resize(allocproseq);	
		}
		if(st.substr(0,4)=="ATOM")
		{
			proseq[numproseq].simpletype=1;
		}
		else if(st.substr(0,6)=="HETATM")
		{
			proseq[numproseq].simpletype=2;
		}

		proseq[numproseq].seqno=st.substr(4,7).c_str();
		proseq[numproseq].detailtype=st.substr(12,4);

		auto it = Atom2id.find(st.substr(12,4));
		if(it != Atom2id.end()){
			proseq[numproseq].detail=Atom2id[st.substr(12,4)];
		}else{
			proseq[numproseq].detail=ATOM_TYPE_NUM;
		}

		//proseq_tmp.alt=st.substr(16,1)[0];
		proseq[numproseq].alt=st.substr(16,1)[0];
		
		//amino

		//proseq_tmp.residue=st.substr(17,3);
		proseq[numproseq].residue=st.substr(17,3);

		if(proseq[numproseq].residue!="  A" && proseq[numproseq].residue!="  C" && proseq[numproseq].residue!="  G" && proseq[numproseq].residue!="  U"){
			if(proseq[numproseq].residue=="RA3" || proseq[numproseq].residue=="RC3" || proseq[numproseq].residue=="RG3" || proseq[numproseq].residue=="RU3" || proseq[numproseq].residue=="DA" || proseq[numproseq].residue=="DC" || proseq[numproseq].residue=="DG" || proseq[numproseq].residue=="DT"){
				proseq[numproseq].residueid=PdbRes2ResIdchar[proseq[numproseq].residue];
				proseq[numproseq].residue=nonstPdbRes2stPdbRes[proseq[numproseq].residue];
			}
			else{
				proseq[numproseq].residueid='4';
				proseq[numproseq].residue="  X";
			}
		}
		//if(proseq_tmp.residue!="  A" && proseq_tmp.residue!="  C" && proseq_tmp.residue!="  G" && proseq_tmp.residue!="  U"){
			//proseq_tmp.residueid=4;
		//}

		else{
			proseq[numproseq].residueid=PdbRes2ResIdchar[proseq[numproseq].residue];
			//proseq_tmp.residueid=pres2id[proseq_tmp.residue];
		}

		//modelid
		proseq[numproseq].modelid=numpromod;

		//chainid
		proseq[numproseq].chainid=st.substr(21,1)[0];
		//proseq_tmp.chainid=st.substr(21,1)[0];

		proseq[numproseq].resno=atoi(st.substr(22,4).c_str());
		//proseq_tmp.resno=atoi(st.substr(22,4).c_str());

		proseq[numproseq].ins=st.substr(26,1)[0];
		//proseq_tmp.ins=st.substr(26,1)[0];

		proseq[numproseq].x=float(atof(st.substr(30,8).c_str())); 
		proseq[numproseq].y=float(atof(st.substr(38,8).c_str())); 
		proseq[numproseq].z=float(atof(st.substr(46,8).c_str())); 
		//proseq_tmp.x=float(atof(st.substr(30,8).c_str())); 
		//proseq_tmp.y=float(atof(st.substr(38,8).c_str())); 
		//proseq_tmp.z=float(atof(st.substr(46,8).c_str()));   

		//一番最初のMODEL、chainの場合
		if(promod[numpromod].nchain==0 && flagbeginmodel)
		{
			promod[numpromod].procha[promod[numpromod].nchain].chainid=st.substr(21,1)[0];
			promod[numpromod].procha[promod[numpromod].nchain].chainseg.init=numproseq;	
			flagbeginmodel=false;
		}
		
		//TERの次のレコードである||現在と1つ前の残基がどちらもHETATMで、1つ前の残基とchainが異なる or 現在のATOMと1つ前の残基とが異なるchainである場合
		else if((flagter) 
		|| (st.substr(21,1)[0] != proseq[numproseq-1].chainid && !(st.substr(0,6)!="HETATM" && proseq[numproseq-1].simpletype==2))
		|| (st.substr(21,1)[0] != proseq[numproseq-1].chainid && st.substr(0,4)=="ATOM"))
		{
			//numproseq-1:終了位置のatomのid term:: terminate
			promod[numpromod].procha[promod[numpromod].nchain].chainseg.term=numproseq-1;
			//新たなchainに変更 nchain=nchain+1
			promod[numpromod].nchain+=1;
			//chainのidを格納
			promod[numpromod].procha[promod[numpromod].nchain].chainid=st.substr(21,1)[0];
			//numproseq: 開始位置のatomのid init:: initial
			promod[numpromod].procha[promod[numpromod].nchain].chainseg.init=numproseq;
		}

		if(flagter)
		{
			flagter=false;
		}
		numproseq++;
		
		
		return 17;
	}
	return 0;
	
}

//42:最初のモデルのMODELレコード
//15:ENDMDLレコード
//16:TERレコード
//17:ATOMのレコード
//41:Alternate location識別子のレコード


bool LoadPDB::loadpdb(string pdbfile,vector<Atom> &proseq, vector<Promodel> &promod)
{
	FILE	*filein;
	string oneline_str;
	char	oneline[600];	
	bool flagf;
	bool modelend=false;
	int headnumchain;//number of chains in header file
	int indexchain;// num of aminos in each chain

	numpromod=0;allocpromod=1;
	numproseq=0;allocproseq=800;

    filein = fopen(pdbfile.c_str(), "r");
	if(filein==NULL)
	{
		printf( "PDB file can not be opened %s\n", pdbfile.c_str());
		return false;
	}

	int i;
	pdbname=pdbfile;
	
	//init
	headnumchain=-1;
	indexchain=0;
	
	promod.resize(allocpromod);

	flagfirstmodel=true;
	flagbeginmodel=true;
	for(i=0;i<allocpromod;i++)
	{
		promod[i].nchain=0;
		for(int j=0;j<MAX_CHAIN_NUM;j++)
		{
			promod[i].procha[j].chainseg.init=-1;
			promod[i].procha[j].chainseg.term=-1;
		}
	}
	
	//read
	char oldline[600];oldline[0]='\0';
	while(!feof(filein) && !modelend)
	{
		flagf=fgets(oneline,600,filein);
		if(strcmp(oneline,oldline)!=0 || oldline[0]=='\0')
		{
			oneline_str=oneline;
			dealline(oneline_str,proseq, promod);
		}
		strcpy(oldline,oneline);
    }

	fclose(filein);
	promod[numpromod].nchain+=1;
	//realloc
	if(promod[numpromod].procha[promod[numpromod].nchain-1].chainseg.term==-1)
	{
		promod[numpromod].procha[promod[numpromod].nchain-1].chainseg.term=numproseq-1;	
	}

	if(numproseq!=0){
		proseq.resize(numproseq);
		proseq.shrink_to_fit();
	}
	if(numpromod!=-1){
		promod.resize(numpromod+1);
		promod.shrink_to_fit();
	}

	return true;
}


void LoadPDB::extractbb(vector<Atom> &proseq, vector<Promodel> &promod, vector<Boneinfo> &bb, int model, int chain)
{
	if(numproseq==0 || promod.size()==0)
	{
		return;
	}
	if(model<-1 || model>numpromod)
	{
		printf("no this model\n");
		return;
	}
	
	if( (model>=0 && chain> promod[model].nchain)  || chain<-1)
	{
		printf("no this chain\n");
		return;
	}

	numbb=0;
	if(numproseq==0) return;
	int seqinit,seqterm;
	if(model==-1)
	{
		seqinit=promod[0].procha[0].chainseg.init;
		seqterm=promod[numpromod].procha[promod[numpromod].nchain-1].chainseg.term;
	}
	else if(chain==-1)
	{
		seqinit=promod[model].procha[0].chainseg.init;
		seqterm=promod[model].procha[promod[model].nchain-1].chainseg.term;
	}	
	else
	{
		seqinit=promod[model].procha[chain].chainseg.init;
		seqterm=promod[model].procha[chain].chainseg.term;
	}

	if(seqinit==-1 || seqterm==-1)
	{
		printf("no atoms in the model\n");
		return;
	}

	int allocbb=1+numproseq/4;
	bb.resize(allocbb);
	int i,j,k,l;
	string tpname;
	int indnt;
	
	for(i=seqinit;i<=seqterm;i++)
	{
		j=i;
		while(j<=seqterm && proseq[i].resno==proseq[j].resno && proseq[i].chainid==proseq[j].chainid && proseq[i].ins==proseq[j].ins && proseq[i].chainid==proseq[j].chainid)
		{
			j++;
		}
		//for(k=0;k<13;k++)
		//{
			//bb[numbb].indall[k]=-1;
		//}
		
		for(l=0;l<Atom2id.size();l++) bb[numbb].indatom[l]=-1;

		for(k=i;k<j;k++)
		{
			bb[numbb].indatom[proseq[k].detail] = k;
		}
		bb[numbb].istart=i;
		bb[numbb].iend=j-1;
		bb[numbb].resid=ResIdchar2Reschar[proseq[i].residueid];

		bb[numbb].resind=proseq[i].resno;

		numbb++;

		if(numbb==allocbb)
		{
			allocbb*=2;
			bb.resize(allocbb);
		}
		i=j-1;
	}
	if(numbb!=0){
		bb.resize(numbb);
		bb.shrink_to_fit();
	}
}

LoadPDB::~LoadPDB()
{	
}













