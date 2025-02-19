// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

RNASimulation::RNASimulation()
{
}

RNASimulation::~RNASimulation()
{
}

void RNASimulation::checkupdate(vector<Residueatoms> decstr, int numseq){
	int i,j;
	char aat;
	int atom1, atom2, atom3, atom4;
	BasicFunc bf;
	Point3s tmp1,tmp2,tmp3,tmp4;
	float len, ang, dih;
	int idx;

	for(i=1;i<numseq;i++)
	{
		aat=decstr[i].aaa;
		atom1=ParamsList[aat][0][0];atom2=ParamsList[aat][0][1];atom3=ParamsList[aat][0][2];atom4=ParamsList[aat][0][3];
		tmp1=decstr[i-1].ptv[atom1];tmp2=decstr[i-1].ptv[atom2];tmp3=decstr[i-1].ptv[atom3];tmp4=decstr[i].ptv[atom4];
		coord2params(tmp1, tmp2, tmp3, tmp4, len, ang, dih);
		if(abs(decstr[i].len[0]-len) > 0.1 | abs(decstr[i].ang[0]-ang) > 1 | abs(decstr[i].dih[0]-dih) > 1){
			cout<<i<<endl;
			cout<<to_string(decstr[i].len[0]) + "-" + to_string(decstr[i].ang[0]) + "-" + to_string(decstr[i].dih[0])<<endl;
			cout<<to_string(len) + "-" + to_string(ang) + "-" + to_string(dih)<<endl;
			cout<<"ErrorA1"<<endl;
			exit(1);
		}
		
		atom1=ParamsList[aat][1][0];atom2=ParamsList[aat][1][1];atom3=ParamsList[aat][1][2];atom4=ParamsList[aat][1][3];
		tmp1=decstr[i-1].ptv[atom1];tmp2=decstr[i-1].ptv[atom2];tmp3=decstr[i].ptv[atom3];tmp4=decstr[i].ptv[atom4];
		coord2params(tmp1, tmp2, tmp3, tmp4, len, ang, dih);
		if(abs(decstr[i].len[1]-len) > 0.1 | abs(decstr[i].ang[1]-ang) > 1 | abs(decstr[i].dih[1]-dih) > 1){
			cout<<i<<endl;
			cout<<to_string(decstr[i].len[1]) + "-" + to_string(decstr[i].ang[1]) + "-" + to_string(decstr[i].dih[1])<<endl;
			cout<<to_string(len) + "-" + to_string(ang) + "-" + to_string(dih)<<endl;
			cout<<"ErrorA2"<<endl;
			exit(1);
		}
		
		atom1=ParamsList[aat][2][0];atom2=ParamsList[aat][2][1];atom3=ParamsList[aat][2][2];atom4=ParamsList[aat][2][3];
		tmp1=decstr[i-1].ptv[atom1];tmp2=decstr[i].ptv[atom2];tmp3=decstr[i].ptv[atom3];tmp4=decstr[i].ptv[atom4];
		coord2params(tmp1, tmp2, tmp3, tmp4, len, ang, dih);
		if(abs(decstr[i].len[2]-len) > 0.1 | abs(decstr[i].ang[2]-ang) > 1 | abs(decstr[i].dih[2]-dih) > 1){
			cout<<i<<endl;
			cout<<to_string(decstr[i].len[2]) + "-" + to_string(decstr[i].ang[2]) + "-" + to_string(decstr[i].dih[2])<<endl;
			cout<<to_string(len) + "-" + to_string(ang) + "-" + to_string(dih)<<endl;
			cout<<"ErrorA3"<<endl;
			exit(1);
		}
		
		for(j=3;j<AccAtomParamNum[aat][0];j++){
			atom1=ParamsList[aat][j][0];atom2=ParamsList[aat][j][1];atom3=ParamsList[aat][j][2];atom4=ParamsList[aat][j][3];
			tmp1=decstr[i].ptv[atom1];tmp2=decstr[i].ptv[atom2];tmp3=decstr[i].ptv[atom3];tmp4=decstr[i].ptv[atom4];
			
			coord2params(tmp1, tmp2, tmp3, tmp4, len, ang, dih);
			if(abs(decstr[i].len[j]-len) > 0.1 | abs(decstr[i].ang[j]-ang) > 1 | abs(decstr[i].dih[j]-dih) > 1){
				cout<<i<<endl;
				cout<<to_string(decstr[i].len[j]) + "-" + to_string(decstr[i].ang[j]) + "-" + to_string(decstr[i].dih[j])<<endl;
				cout<<to_string(len) + "-" + to_string(ang) + "-" + to_string(dih)<<endl;
				cout<<"ErrorA"+to_string(j)<<endl;
				exit(1);
			}
		}
	}


	aat=decstr[0].aaa;
	for(j=AccAtomParamNum[aat][0];j<AccAtomParamNum[aat][2];j++)
	{
		auto it = std::find(MODTARGET.begin(), MODTARGET.end(), j);
		if(it != MODTARGET.end())
		{
			if(j == 10 || j == 11){
				atom1=FIRSTOPParamsList[aat][j][0];atom2=FIRSTOPParamsList[aat][j][1];atom3=FIRSTOPParamsList[aat][j][2];atom4=FIRSTOPParamsList[aat][j][3];
			}
			//else if(j == 11){
				//atom1=EDGEOPParamsList[aat][j][0];atom2=EDGEOPParamsList[aat][j][1];atom3=EDGEOPParamsList[aat][j][2];atom4=EDGEOPParamsList[aat][j][3];
			//}
			else{
				atom1=ParamsList[aat][j][0];atom2=ParamsList[aat][j][1];atom3=ParamsList[aat][j][2];atom4=ParamsList[aat][j][3];
			}
			tmp1=decstr[0].ptv[atom1];tmp2=decstr[0].ptv[atom2];tmp3=decstr[0].ptv[atom3];tmp4=decstr[0].ptv[atom4];
			coord2params(tmp1, tmp2, tmp3, tmp4, len, ang, dih);
			if(abs(decstr[0].len[j]-len) > 0.1 | abs(decstr[0].ang[j]-ang) > 1 | abs(decstr[0].dih[j]-dih) > 1){
				cout<<"0"<<endl;
				cout<<to_string(decstr[0].len[j]) + "-" + to_string(decstr[0].ang[j]) + "-" + to_string(decstr[0].dih[j])<<endl;
				cout<<to_string(len) + "-" + to_string(ang) + "-" + to_string(dih)<<endl;
				cout<<"ErrorB"+to_string(j)<<endl;
				exit(1);
			}
		}
	}

	
	for(i=1;i<numseq;i++)
	{
		aat=decstr[i].aaa;
		for(j=AccAtomParamNum[aat][0];j<AccAtomParamNum[aat][2];j++)
		{
			auto it = std::find(MODTARGET.begin(), MODTARGET.end(), j);
			if(it != MODTARGET.end())
			{
				atom1=ParamsList[aat][j][0];atom2=ParamsList[aat][j][1];atom3=ParamsList[aat][j][2];atom4=ParamsList[aat][j][3];
				if(j == 10){
					tmp1=decstr[i-1].ptv[atom1];tmp2=decstr[i].ptv[atom2];tmp3=decstr[i].ptv[atom3];tmp4=decstr[i].ptv[atom4];
				}else if(j == 11){
					tmp1=decstr[i].ptv[atom1];tmp2=decstr[i-1].ptv[atom2];tmp3=decstr[i].ptv[atom3];tmp4=decstr[i].ptv[atom4];
				}else{
					tmp1=decstr[i].ptv[atom1];tmp2=decstr[i].ptv[atom2];tmp3=decstr[i].ptv[atom3];tmp4=decstr[i].ptv[atom4];
				}
				coord2params(tmp1, tmp2, tmp3, tmp4, len, ang, dih);

				if(abs(decstr[i].len[j]-len) > 0.1 | abs(decstr[i].ang[j]-ang) > 1 | abs(decstr[i].dih[j]-dih) > 1){
					cout<<i<<endl;
					cout<<to_string(decstr[i].len[j]) + "-" + to_string(decstr[i].ang[j]) + "-" + to_string(decstr[i].dih[j])<<endl;
					cout<<to_string(len) + "-" + to_string(ang) + "-" + to_string(dih)<<endl;
					cout<<"ErrorB"+to_string(j)<<endl;
					exit(1);
				}
			}
		}


		/*
		for(j=AccAtomParamNum[aat][2];j<AccAtomParamNum[aat][3];j++)
		{
			cout<<j<<endl;
			atom1=ParamsList[aat][j][0];
			atom2=ParamsList[aat][j][1];
			atom3=ParamsList[aat][j][2];
			atom4=ParamsList[aat][j][3];
			
			tmp1=decstr[i].ptv[atom1];
			tmp2=decstr[i].ptv[atom2];
			tmp3=decstr[i].ptv[atom3];
			tmp4=decstr[i].ptv[atom4];
			
			coord2params(tmp1, tmp2, tmp3, tmp4, len, ang, dih);

			if(dih >=350 & decstr[i].dih[j] < 350) dih=dih-360;

			if(abs(decstr[i].len[j]-len) > 0.1 | abs(decstr[i].ang[j]-ang) > 1 | abs(decstr[i].dih[j]-dih) > 1){
				cout<<i<<endl;
				cout<<to_string(decstr[i].len[j]) + "-" + to_string(decstr[i].ang[j]) + "-" + to_string(decstr[i].dih[j])<<endl;
				cout<<to_string(len) + "-" + to_string(ang) + "-" + to_string(dih)<<endl;
				cout<<"ErrorC"+to_string(j)<<endl;
				exit(1);
			}
		}
		*/
	}

	for(i=0;i<numseq;i++)
	{
		aat=decstr[i].aaa;
		idx=0;
		for(j=AccAtomParamNum[aat][3];j<AccAtomParamNum[aat][4];j++)
		{
			
			atom1=HydParamsList[aat][idx][0]; atom2=HydParamsList[aat][idx][1]; atom3=HydParamsList[aat][idx][2]; atom4=HydParamsList[aat][idx][3];
			tmp1=decstr[i].ptv[atom1];tmp2=decstr[i].ptv[atom2];tmp3=decstr[i].ptv[atom3];tmp4=decstr[i].ptv[atom4];
			
			coord2params(tmp1, tmp2, tmp3, tmp4, len, ang, dih);

			if(abs(decstr[i].hydlen[idx]-len) > 0.1 | abs(decstr[i].hydang[idx]-ang) > 1 | abs(decstr[i].hyddih[idx]-dih) > 1){
				cout<<i<<endl;
				cout<<to_string(decstr[i].hydlen[idx]) + "-" + to_string(decstr[i].hydang[idx]) + "-" + to_string(decstr[i].hyddih[idx])<<endl;
				cout<<to_string(len) + "-" + to_string(ang) + "-" + to_string(dih)<<endl;
				cout<<"ErrorH"+to_string(j)<<endl;
				exit(1);
			}

			idx+=1;
		}
	}

	
	for(i=1;i<numseq;i++)
	{
		aat=decstr[i].aaa;
		atom1=AddAngleParamsList[aat][0][0]; atom2=AddAngleParamsList[aat][0][1]; atom3=AddAngleParamsList[aat][0][2];
		tmp1=decstr[i-1].ptv[atom1]; tmp2=decstr[i].ptv[atom2]; tmp3=decstr[i].ptv[atom3];
		coord2ang(tmp1, tmp2, tmp3, ang);
		if(abs(decstr[i].addang[0]-ang) > 1){
			cout<<i<<endl;
			cout<<to_string(decstr[i].addang[0]) + "-" + to_string(ang)<<endl;
			cout<<"ErrorD"+to_string(j)<<endl;
			exit(1);
		}

		for(j=1;j<AddAngleParamsList[aat].size();j++)
		{
			atom1=AddAngleParamsList[aat][j][0]; atom2=AddAngleParamsList[aat][j][1]; atom3=AddAngleParamsList[aat][j][2];
			
			tmp1=decstr[i].ptv[atom1]; tmp2=decstr[i].ptv[atom2]; tmp3=decstr[i].ptv[atom3];
			coord2ang(tmp1, tmp2, tmp3, ang);
			if(abs(decstr[i].addang[j]-ang) > 1){
				cout<<to_string(i) + "-" + to_string(j) + "-" + to_string(decstr[i].addang[j]) + "-" + to_string(ang)<<endl;
				cout<<"ErrorD"+to_string(j)<<endl;
				exit(1);
			}

		}

		for(j=0;j<AddLengthParamsList[aat].size();j++)
		{
			atom1=AddLengthParamsList[aat][j][0]; atom2=AddLengthParamsList[aat][j][1];
			tmp1=decstr[i].ptv[atom1]; tmp2=decstr[i].ptv[atom2];
			coord2len(tmp1, tmp2, len);
			if(abs(decstr[i].addlen[j]-len) > 0.1){
				cout<<i<<endl;
				cout<<to_string(decstr[i].addlen[j]) + "-" + to_string(len)<<endl;
				cout<<"ErrorD"+to_string(j)<<endl;
				exit(1);
			}
		}
	}
}

void RNASimulation::setStr(vector<Residueatoms> &decstr, vector<Atom> proseq, vector<Boneinfo> bb, string seq)
{
	int i,j,k;
	BPBasicInfo basepair;
	int nti,ntj;
	int stbaseidx,edbaseidx;
	int pairidx;

	//toroutflag.resize(bb.size(), vector<bool>(TorsionClsList[0].size(), false));
	//torclsidx.resize(bb.size());
	
	for(i=0;i<bb.size();i++)
	{
		decstr[i].aaa=bb[i].resid;
		decstr[i].resind=Reschar2Idint[decstr[i].aaa];
		decstr[i].residx=i;
		
		for(k=0;k<NNheavy[decstr[i].resind][0];k++)
		{
			decstr[i].ptv[k].x=proseq[bb[i].indatom[NNheavy[decstr[i].resind][1+k]]].x;
			decstr[i].ptv[k].y=proseq[bb[i].indatom[NNheavy[decstr[i].resind][1+k]]].y;
			decstr[i].ptv[k].z=proseq[bb[i].indatom[NNheavy[decstr[i].resind][1+k]]].z;
		}
	}

	for(i=0;i<bb.size();i++)
	{
		decstr[i].ss=0;
		decstr[i].base.baseidx=-1;
		decstr[i].base.basetype=-1;
	}

	for(i=0; i<bb.size(); i++){
		bppairednt[i]=-1;
		stpairednt[i]=-1;
	}

	for(i=0;i<BASEPAIRS.size();i++){
		nti=BASEPAIRS[i].basepidx.first;
		ntj=BASEPAIRS[i].basepidx.second;
		bpmapid2idx[nti*bb.size()+ntj]=i;
		bpmapid2idx[ntj*bb.size()+nti]=i;

		bppairednt[nti]=ntj;
		bppairednt[ntj]=nti;
		
		decstr[nti].ss=1;
		decstr[ntj].ss=1;
		decstr[nti].base.baseidx=ntj;
		decstr[ntj].base.baseidx=nti;
		decstr[nti].base.basetype=BaseTypeList[seq[nti]][seq[ntj]];
		decstr[ntj].base.basetype=BaseTypeList[seq[ntj]][seq[nti]];
	}

	for(i=0;i<STACKINGPAIRS.size();i++){
		nti=STACKINGPAIRS[i].basepidx.first;
                ntj=STACKINGPAIRS[i].basepidx.second;
		stpairednt[nti]=ntj;
		stpairednt[ntj]=nti;

		bsmapid2idx[nti*bb.size()+ntj]=i;
		bsmapid2idx[ntj*bb.size()+nti]=i;
	}
}

void RNASimulation::rcopystr(vector<Residueatoms> &decstr, vector<Residueatoms> restr, int numseq)
{
	int i,j,k;
	//decstr.clear();
	//decstr.shrink_to_fit();
	decstr.resize(restr.size());
	std::copy(restr.begin(), restr.end(), decstr.begin());
	memcpy(lp1.nn,lp2.nn,TOTAL_ENERGY_NUM*sizeof(int));
	memcpy(lp1.fn,lp2.fn,TOTAL_ENERGY_NUM*sizeof(double));

	rppmat.resize(rppmatbk.size());
	bpenergy.resize(bpenergybk.size());
	bsenergy.resize(bsenergybk.size());
	nmap.resize(nmapbk.size());
	bpenergy=bpenergybk;
	bsenergy=bsenergybk;
	for(i=0;i<numseq;i++){
		rppmat[i].resize(rppmatbk[i].size());
		nmap[i].resize(nmapbk[i].size());
		rppmat[i]=rppmatbk[i];
		nmap[i]=nmapbk[i];
	}

	for(i=0;i<numseq*numseq;i++){
		contactene[i]=contactenebk[i];
	}

	mapflag=mapflagbk;
	mapbbflag=mapbbflagbk;
	distmap=distmapbk;
	//contactene=contactenebk;
	ntpairdist=ntpairdistbk;
	adjpairdist=adjpairdistbk;
	distene=distenebk;
	multidistene=multidistenebk;

	lp1.indn.resize(lp2.indn.size());
	lp1.indf.resize(lp2.indf.size());
	lp1.flap.resize(lp2.flap.size());
	for(i=0;i<TOTAL_ENERGY_NUM;i++)
	{
		copy(lp2.indn[i].begin(), lp2.indn[i].end(), lp1.indn[i].begin());
		copy(lp2.indf[i].begin(), lp2.indf[i].end(), lp1.indf[i].begin());
		copy(lp2.flap[i].begin(), lp2.flap[i].end(), lp1.flap[i].begin());
	}
}

void RNASimulation::copystr(vector<Residueatoms> decstr, int numseq)
{
	int i, j, k;
	//restr.clear();
	//restr.shrink_to_fit();
	//restr.resize(decstr.size());
	std::copy(decstr.begin(), decstr.end(), restr.begin());

	memcpy(lp2.nn,lp1.nn,TOTAL_ENERGY_NUM*sizeof(int));
	memcpy(lp2.fn,lp1.fn,TOTAL_ENERGY_NUM*sizeof(double));

	rppmatbk.resize(rppmat.size());
	bpenergybk.resize(bpenergy.size());
	bsenergybk.resize(bsenergy.size());
	nmapbk.resize(nmap.size());
	bpenergybk=bpenergy;
	bsenergybk=bsenergy;
	for(i=0;i<numseq;i++){
		rppmatbk[i].resize(rppmat[i].size());
		nmapbk[i].resize(nmap[i].size());
		rppmatbk[i]=rppmat[i];
		nmapbk[i]=nmap[i];
	}
	for(i=0;i<numseq*numseq;i++){
		contactenebk[i]=contactene[i];
	}

	mapflagbk=mapflag;
	mapbbflagbk=mapbbflag;
	distmapbk=distmap;
	ntpairdistbk=ntpairdist;
	adjpairdistbk=adjpairdist;
	//contactenebk=contactene;
	distenebk=distene;
	multidistenebk=multidistene;

	lp2.indn.resize(TOTAL_ENERGY_NUM);
	lp2.indf.resize(TOTAL_ENERGY_NUM);
	lp2.flap.resize(TOTAL_ENERGY_NUM);
	for(i=0;i<TOTAL_ENERGY_NUM;i++)
	{
		lp2.indn[i].resize(lp1.indn[i].size());
		lp2.indf[i].resize(lp1.indf[i].size());
		lp2.flap[i].resize(lp1.flap[i].size());

		std::copy(lp1.indn[i].begin(), lp1.indn[i].end(), lp2.indn[i].begin());
		std::copy(lp1.indf[i].begin(), lp1.indf[i].end(), lp2.indf[i].begin());
		std::copy(lp1.flap[i].begin(), lp1.flap[i].end(), lp2.flap[i].begin());
	}

}

void RNASimulation::fsetinitial(Energyterms &lp, vector<vector<Rescontact>> &rpp, vector<Basecontact> &bprpp, vector<Stackcontact> &bsrpp, vector<vector<int>> &nmaprpp, vector<vector<vector<double>>> &distmaprpp, vector<vector<vector<double>>> &ntpairdistrpp, vector<vector<vector<double>>> &adjpairdistrpp, vector<double> &distenerpp, vector<vector<double>> &multidistenerpp, vector<Phycontact> &contactenerpp, vector<vector<bool>> &mapflagrpp, vector<vector<bool>> &mapbbflagrpp, int numseq){
	int i,j,k,l;
	//lp.indn.clear(); lp.indn.shrink_to_fit();
	//lp.indf.clear(); lp.indf.shrink_to_fit();
	//lp.flap.clear(); lp.flap.shrink_to_fit();
	//bprpp.clear(); bprpp.shrink_to_fit();
	//bsrpp.clear(); bsrpp.shrink_to_fit();
	//rpp.clear(); rpp.shrink_to_fit();
	//mapflagrpp.clear(); mapflagrpp.shrink_to_fit();
	//nmaprpp.clear(); nmaprpp.shrink_to_fit();
	//distmaprpp.clear(); distmaprpp.shrink_to_fit();
	//distenerpp.clear(); distenerpp.shrink_to_fit();
	//contactenerpp.clear(); contactenerpp.shrink_to_fit();
	//multidistenerpp.clear(); multidistenerpp.shrink_to_fit();

	//cout<<"bbb"<<endl;

	//lp.indn.resize(TOTAL_ENERGY_NUM);
	//lp.indf.resize(TOTAL_ENERGY_NUM);
	//lp.flap.resize(TOTAL_ENERGY_NUM);
	
	//for(i=0;i<TOTAL_ENERGY_NUM;i++)
	//{
		//lp.indn[i].resize(numseq);
		//lp.indf[i].resize(numseq);
		//lp.flap[i].resize(numseq);
	//}

	//cout<<"ccc"<<endl;

	//rpp.resize(numseq, vector<Rescontact>(numseq));
	//nmaprpp.resize(numseq, vector<int>(numseq, 0));
	//mapflagrpp.resize(numseq, vector<bool>(numseq, false));
	
	//distmaprpp.resize(numseq, vector<vector<double>>(numseq, vector<double>(DISTTARGETATOMS.size(), 0.0)));
	//distenerpp.resize(numseq*numseq,0.0);
	//contactenerpp.resize(numseq*numseq,0.0);
	//multidistenerpp.resize(npaa3.size(), vector<double>(numseq*numseq, 0.0));

	//bprpp.resize(BASEPAIRS.size());
	//bsrpp.resize(STACKINGPAIRS.size());

	for (auto& row : mapflagrpp) {
    	fill(row.begin(), row.end(), false);
	}

	for (auto& row : mapbbflagrpp) {
    	fill(row.begin(), row.end(), false);
	}

	for (auto& row : nmaprpp) {
    	fill(row.begin(), row.end(), 0);
	}

	for (auto& row : distmaprpp) {
    	for (auto& subrow : row) {
        	fill(subrow.begin(), subrow.end(), 0.0);
    	}
	}

	for (auto& row : ntpairdistrpp) {
    	for (auto& subrow : row) {
        	fill(subrow.begin(), subrow.end(), 0.0);
    	}
	}

	for (auto& row : adjpairdistrpp) {
    	for (auto& subrow : row) {
        	fill(subrow.begin(), subrow.end(), 0.0);
    	}
	}

	fill(distenerpp.begin(), distenerpp.end(), 0.0);
	//fill(contactenerpp.begin(), contactenerpp.end(), 0.0);

	for (auto& row : multidistenerpp) {
    	fill(row.begin(), row.end(), 0.0);
	}

	for(i=0;i<BASEPAIRS.size();i++){
		for(k=0;k<BASE_ENERGY_NUM+BASEHYD_ENERGY_NUM;k++)
		{
			bprpp[i].pval[k]=0;
			bprpp[i].pcw[k]=0;
		}
		for(k=0;k<BASEPAIRTEMP_ENERGY_NUM;k++)
		{
			bprpp[i].tval[k]=0;
		}
	}

	for(i=0;i<STACKINGPAIRS.size();i++){
		for(k=0;k<STACKING_ENERGY_NUM;k++)
		{
			bsrpp[i].pval[k]=0;
			bsrpp[i].pcw[k]=0;
		}
		for(k=0;k<STACKINGTEMP_ENERGY_NUM;k++)
		{
			bsrpp[i].tval[k]=0;
		}
	}

	for(i=0;i<TOTAL_ENERGY_NUM;i++)
	{
		lp.nn[i]=0;lp.fn[i]=0;
		fill(lp.indn[i].begin(), lp.indn[i].end(), 0);
		fill(lp.indf[i].begin(), lp.indf[i].end(), 0);
		fill(lp.flap[i].begin(), lp.flap[i].end(), 0);
	}

	for(i=0;i<numseq;i++) 
	{
		for(j=0;j<numseq;j++)
		{

			for(k=0;k<CLASH_ENERGY_NUM;k++)
			{
				//rpp[i][j].flag[k]=false;
				rpp[i][j].pval[k]=0;
				rpp[i][j].pcw[k]=0;
			}

			l=j*numseq+i;
			for(k=0;k<CONTACT_ENERGY_NUM;k++)
			{
				contactenerpp[l].pval[k]=0;
				contactenerpp[l].pcw[k]=0;
			}
		}
	}
}

void RNASimulation::fsetinitialarray(int numseq)
{
	int i,j,k,l;
	if(restr.size()!=0){
		restr.clear(); restr.shrink_to_fit();

		mapflag.clear(); mapflag.shrink_to_fit();
		mapflagbk.clear(); mapflagbk.shrink_to_fit();

		mapbbflag.clear(); mapbbflag.shrink_to_fit();
		mapbbflagbk.clear(); mapbbflagbk.shrink_to_fit();

		distmap.clear(); distmap.shrink_to_fit();
		distmapbk.clear(); distmapbk.shrink_to_fit();

		ntpairdist.clear(); ntpairdist.shrink_to_fit();
		ntpairdistbk.clear(); ntpairdistbk.shrink_to_fit();

		adjpairdist.clear(); adjpairdist.shrink_to_fit();
		adjpairdistbk.clear(); adjpairdistbk.shrink_to_fit();

		distene.clear(); distene.shrink_to_fit();
		distenebk.clear(); distenebk.shrink_to_fit();

		contactene.clear(); contactene.shrink_to_fit();
		contactenebk.clear(); contactenebk.shrink_to_fit();

		multidistene.clear(); multidistene.shrink_to_fit();
		multidistenebk.clear(); multidistenebk.shrink_to_fit();

		nmap.clear(); nmap.shrink_to_fit();
		nmapbk.clear(); nmapbk.shrink_to_fit();

		rppmat.clear(); rppmat.shrink_to_fit();
		rppmatbk.clear(); rppmatbk.shrink_to_fit();

		bpenergy.clear(); bpenergy.shrink_to_fit();
		bpenergybk.clear(); bpenergybk.shrink_to_fit();

		bsenergy.clear(); bsenergy.shrink_to_fit();
		bsenergybk.clear(); bsenergybk.shrink_to_fit();
	}

	restr.resize(numseq);
	lp1.indn.resize(TOTAL_ENERGY_NUM);lp1.indf.resize(TOTAL_ENERGY_NUM);lp1.flap.resize(TOTAL_ENERGY_NUM);
	lp2.indn.resize(TOTAL_ENERGY_NUM);lp2.indf.resize(TOTAL_ENERGY_NUM);lp2.flap.resize(TOTAL_ENERGY_NUM);
	for(i=0;i<TOTAL_ENERGY_NUM;i++)
	{
		lp1.nn[i]=0;
		lp1.indn[i].resize(numseq);lp1.indf[i].resize(numseq);lp1.flap[i].resize(numseq);

		lp2.nn[i]=0;
		lp2.indn[i].resize(numseq);lp2.indf[i].resize(numseq);lp2.flap[i].resize(numseq);
	}

	nmap.resize(numseq, vector<int>(numseq, 0));
	nmapbk.resize(numseq, vector<int>(numseq, 0));

	mapflag.resize(numseq, vector<bool>(numseq, false));
	mapflagbk.resize(numseq, vector<bool>(numseq, false));

	mapbbflag.resize(numseq, vector<bool>(numseq, false));
	mapbbflagbk.resize(numseq, vector<bool>(numseq, false));

	rppmat.resize(numseq, vector<Rescontact>(numseq));
	rppmatbk.resize(numseq, vector<Rescontact>(numseq));

	distmap.resize(numseq, vector<vector<double>>(numseq, vector<double>(DISTTARGETATOMS.size(), 0.0)));
	distmapbk.resize(numseq, vector<vector<double>>(numseq, vector<double>(DISTTARGETATOMS.size(), 0.0)));

	ntpairdist.resize(numseq, vector<vector<double>>(PAIRLENGTHTARGETNUM, vector<double>(PAIRLENGTHTARGETNUM, 0.0)));
	ntpairdistbk.resize(numseq, vector<vector<double>>(PAIRLENGTHTARGETNUM, vector<double>(PAIRLENGTHTARGETNUM, 0.0)));

	adjpairdist.resize(numseq, vector<vector<double>>(PAIRLENGTHTARGETNUM, vector<double>(PAIRLENGTHTARGETNUM, 0.0)));
	adjpairdistbk.resize(numseq, vector<vector<double>>(PAIRLENGTHTARGETNUM, vector<double>(PAIRLENGTHTARGETNUM, 0.0)));

	distene.resize(numseq*numseq, 0.0);
	distenebk.resize(numseq*numseq, 0.0);

	//contactene.resize(numseq*numseq, 0.0);
	//contactenebk.resize(numseq*numseq, 0.0);
	contactene.resize(numseq*numseq);
	contactenebk.resize(numseq*numseq);

	multidistene.resize(npaa3.size(), vector<double>(numseq*numseq, 0.0));
	multidistenebk.resize(npaa3.size(), vector<double>(numseq*numseq, 0.0));

	bpenergy.resize(BASEPAIRS.size());
	bpenergybk.resize(BASEPAIRS.size());

	bsenergy.resize(STACKINGPAIRS.size());
	bsenergybk.resize(STACKINGPAIRS.size());

	for(i=0;i<TOTAL_ENERGY_NUM;i++)
	{
		lp1.nn[i]=0;lp1.fn[i]=0;
		lp2.nn[i]=0;lp2.fn[i]=0;

		fill(lp1.indn[i].begin(), lp1.indn[i].end(), 0);
		fill(lp1.indf[i].begin(), lp1.indf[i].end(), 0);
		fill(lp1.flap[i].begin(), lp1.flap[i].end(), 0);
		fill(lp2.indn[i].begin(), lp2.indn[i].end(), 0);
		fill(lp2.indf[i].begin(), lp2.indf[i].end(), 0);
		fill(lp2.flap[i].begin(), lp2.flap[i].end(), 0);
	}

	for(i=0;i<numseq;i++) 
	{
		for(j=0;j<numseq;j++)
		{
			for(k=0;k<CLASH_ENERGY_NUM;k++)
			{
				rppmat[i][j].pval[k]=0;rppmat[i][j].pcw[k]=0;
				rppmatbk[i][j].pval[k]=0;rppmatbk[i][j].pcw[k]=0;
			}

			l=j*numseq+i;
			for(k=0;k<CONTACT_ENERGY_NUM;k++)
			{
				contactene[l].pval[k]=0;
				contactene[l].pcw[k]=0;
				contactenebk[l].pval[k]=0;
				contactenebk[l].pcw[k]=0;
			}
		}
	}
}


























