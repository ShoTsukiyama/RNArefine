
#include "MonteCarloSimulation.hpp"

void MonteCarloSimulation::checkupdate(){
	int i,j;
	char aat;
	int atom1, atom2, atom3, atom4;
	BasicFunc bf;
	Point3s tmp1,tmp2,tmp3,tmp4;
	float len, ang, dih;
	int idx;

	for(i=1;i<numbb;i++)
	{
		aat=refinestr_temp[i].aaa;
		atom1=ParamsList[aat][0][0];atom2=ParamsList[aat][0][1];atom3=ParamsList[aat][0][2];atom4=ParamsList[aat][0][3];
		tmp1=refinestr_temp[i-1].ptv[atom1];tmp2=refinestr_temp[i-1].ptv[atom2];tmp3=refinestr_temp[i-1].ptv[atom3];tmp4=refinestr_temp[i].ptv[atom4];
		coord2params(tmp1, tmp2, tmp3, tmp4, len, ang, dih);
		if(abs(refinestr_temp[i].len[0]-len) > 0.1 | abs(refinestr_temp[i].ang[0]-ang) > 1 | abs(refinestr_temp[i].dih[0]-dih) > 1){
			cout<<i<<endl;
			cout<<to_string(refinestr_temp[i].len[0]) + "-" + to_string(refinestr_temp[i].ang[0]) + "-" + to_string(refinestr_temp[i].dih[0])<<endl;
			cout<<to_string(len) + "-" + to_string(ang) + "-" + to_string(dih)<<endl;
			cout<<"ErrorA1"<<endl;
			exit(1);
		}
		
		atom1=ParamsList[aat][1][0];atom2=ParamsList[aat][1][1];atom3=ParamsList[aat][1][2];atom4=ParamsList[aat][1][3];
		tmp1=refinestr_temp[i-1].ptv[atom1];tmp2=refinestr_temp[i-1].ptv[atom2];tmp3=refinestr_temp[i].ptv[atom3];tmp4=refinestr_temp[i].ptv[atom4];
		coord2params(tmp1, tmp2, tmp3, tmp4, len, ang, dih);
		if(abs(refinestr_temp[i].len[1]-len) > 0.1 | abs(refinestr_temp[i].ang[1]-ang) > 1 | abs(refinestr_temp[i].dih[1]-dih) > 1){
			cout<<i<<endl;
			cout<<to_string(refinestr_temp[i].len[1]) + "-" + to_string(refinestr_temp[i].ang[1]) + "-" + to_string(refinestr_temp[i].dih[1])<<endl;
			cout<<to_string(len) + "-" + to_string(ang) + "-" + to_string(dih)<<endl;
			cout<<"ErrorA2"<<endl;
			exit(1);
		}
		
		atom1=ParamsList[aat][2][0];atom2=ParamsList[aat][2][1];atom3=ParamsList[aat][2][2];atom4=ParamsList[aat][2][3];
		tmp1=refinestr_temp[i-1].ptv[atom1];tmp2=refinestr_temp[i].ptv[atom2];tmp3=refinestr_temp[i].ptv[atom3];tmp4=refinestr_temp[i].ptv[atom4];
		coord2params(tmp1, tmp2, tmp3, tmp4, len, ang, dih);
		if(abs(refinestr_temp[i].len[2]-len) > 0.1 | abs(refinestr_temp[i].ang[2]-ang) > 1 | abs(refinestr_temp[i].dih[2]-dih) > 1){
			cout<<i<<endl;
			cout<<to_string(refinestr_temp[i].len[2]) + "-" + to_string(refinestr_temp[i].ang[2]) + "-" + to_string(refinestr_temp[i].dih[2])<<endl;
			cout<<to_string(len) + "-" + to_string(ang) + "-" + to_string(dih)<<endl;
			cout<<"ErrorA3"<<endl;
			exit(1);
		}
		
		for(j=3;j<AccAtomParamNum[aat][0];j++){
			atom1=ParamsList[aat][j][0];atom2=ParamsList[aat][j][1];atom3=ParamsList[aat][j][2];atom4=ParamsList[aat][j][3];
			tmp1=refinestr_temp[i].ptv[atom1];tmp2=refinestr_temp[i].ptv[atom2];tmp3=refinestr_temp[i].ptv[atom3];tmp4=refinestr_temp[i].ptv[atom4];
			
			coord2params(tmp1, tmp2, tmp3, tmp4, len, ang, dih);
			if(abs(refinestr_temp[i].len[j]-len) > 0.1 | abs(refinestr_temp[i].ang[j]-ang) > 1 | abs(refinestr_temp[i].dih[j]-dih) > 1){
				cout<<i<<endl;
				cout<<to_string(refinestr_temp[i].len[j]) + "-" + to_string(refinestr_temp[i].ang[j]) + "-" + to_string(refinestr_temp[i].dih[j])<<endl;
				cout<<to_string(len) + "-" + to_string(ang) + "-" + to_string(dih)<<endl;
				cout<<"ErrorA"+to_string(j)<<endl;
				exit(1);
			}
		}
	}


	aat=refinestr_temp[0].aaa;
	for(j=AccAtomParamNum[aat][0];j<AccAtomParamNum[aat][2];j++)
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
		tmp1=refinestr_temp[0].ptv[atom1];tmp2=refinestr_temp[0].ptv[atom2];tmp3=refinestr_temp[0].ptv[atom3];tmp4=refinestr_temp[0].ptv[atom4];
		coord2params(tmp1, tmp2, tmp3, tmp4, len, ang, dih);
		if(abs(refinestr_temp[0].len[j]-len) > 0.1 | abs(refinestr_temp[0].ang[j]-ang) > 1 | abs(refinestr_temp[0].dih[j]-dih) > 1){
			cout<<"0"<<endl;
			cout<<to_string(refinestr_temp[0].len[j]) + "-" + to_string(refinestr_temp[0].ang[j]) + "-" + to_string(refinestr_temp[0].dih[j])<<endl;
			cout<<to_string(len) + "-" + to_string(ang) + "-" + to_string(dih)<<endl;
			cout<<"ErrorB"+to_string(j)<<endl;
			exit(1);
		}
	}

	
	for(i=1;i<numbb;i++)
	{
		aat=refinestr_temp[i].aaa;
		for(j=AccAtomParamNum[aat][0];j<AccAtomParamNum[aat][2];j++)
		{
			atom1=ParamsList[aat][j][0];atom2=ParamsList[aat][j][1];atom3=ParamsList[aat][j][2];atom4=ParamsList[aat][j][3];
			if(j == 10){
				tmp1=refinestr_temp[i-1].ptv[atom1];tmp2=refinestr_temp[i].ptv[atom2];tmp3=refinestr_temp[i].ptv[atom3];tmp4=refinestr_temp[i].ptv[atom4];
			}else if(j == 11){
				tmp1=refinestr_temp[i].ptv[atom1];tmp2=refinestr_temp[i-1].ptv[atom2];tmp3=refinestr_temp[i].ptv[atom3];tmp4=refinestr_temp[i].ptv[atom4];
			}else{
				tmp1=refinestr_temp[i].ptv[atom1];tmp2=refinestr_temp[i].ptv[atom2];tmp3=refinestr_temp[i].ptv[atom3];tmp4=refinestr_temp[i].ptv[atom4];
			}
			coord2params(tmp1, tmp2, tmp3, tmp4, len, ang, dih);

			if(abs(refinestr_temp[i].len[j]-len) > 0.1 | abs(refinestr_temp[i].ang[j]-ang) > 1 | abs(refinestr_temp[i].dih[j]-dih) > 1){
				cout<<i<<endl;
				cout<<to_string(refinestr_temp[i].len[j]) + "-" + to_string(refinestr_temp[i].ang[j]) + "-" + to_string(refinestr_temp[i].dih[j])<<endl;
				cout<<to_string(len) + "-" + to_string(ang) + "-" + to_string(dih)<<endl;
				cout<<"ErrorB"+to_string(j)<<endl;
				exit(1);
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
			
			tmp1=refinestr_temp[i].ptv[atom1];
			tmp2=refinestr_temp[i].ptv[atom2];
			tmp3=refinestr_temp[i].ptv[atom3];
			tmp4=refinestr_temp[i].ptv[atom4];
			
			coord2params(tmp1, tmp2, tmp3, tmp4, len, ang, dih);

			if(dih >=350 & refinestr_temp[i].dih[j] < 350) dih=dih-360;

			if(abs(refinestr_temp[i].len[j]-len) > 0.1 | abs(refinestr_temp[i].ang[j]-ang) > 1 | abs(refinestr_temp[i].dih[j]-dih) > 1){
				cout<<i<<endl;
				cout<<to_string(refinestr_temp[i].len[j]) + "-" + to_string(refinestr_temp[i].ang[j]) + "-" + to_string(refinestr_temp[i].dih[j])<<endl;
				cout<<to_string(len) + "-" + to_string(ang) + "-" + to_string(dih)<<endl;
				cout<<"ErrorC"+to_string(j)<<endl;
				exit(1);
			}
		}
		*/
	}

	for(i=0;i<numbb;i++)
	{
		aat=refinestr_temp[i].aaa;
		idx=0;
		for(j=AccAtomParamNum[aat][3];j<AccAtomParamNum[aat][4];j++)
		{
			
			atom1=HydParamsList[aat][idx][0]; atom2=HydParamsList[aat][idx][1]; atom3=HydParamsList[aat][idx][2]; atom4=HydParamsList[aat][idx][3];
			tmp1=refinestr_temp[i].ptv[atom1];tmp2=refinestr_temp[i].ptv[atom2];tmp3=refinestr_temp[i].ptv[atom3];tmp4=refinestr_temp[i].ptv[atom4];
			
			coord2params(tmp1, tmp2, tmp3, tmp4, len, ang, dih);

			if(abs(refinestr_temp[i].hydlen[idx]-len) > 0.1 | abs(refinestr_temp[i].hydang[idx]-ang) > 1 | abs(refinestr_temp[i].hyddih[idx]-dih) > 1){
				cout<<i << "-" << j << "-" << idx << "-" << aat<<endl;
				cout<<to_string(refinestr_temp[i].hydlen[idx]) + "-" + to_string(refinestr_temp[i].hydang[idx]) + "-" + to_string(refinestr_temp[i].hyddih[idx])<<endl;
				cout<<to_string(len) + "-" + to_string(ang) + "-" + to_string(dih)<<endl;

				cout<<to_string(tmp1.x) + "-" + to_string(tmp1.y) + "-" + to_string(tmp1.z)<<endl;
				cout<<to_string(tmp2.x) + "-" + to_string(tmp2.y) + "-" + to_string(tmp2.z)<<endl;
				cout<<to_string(tmp3.x) + "-" + to_string(tmp3.y) + "-" + to_string(tmp3.z)<<endl;
				cout<<to_string(tmp4.x) + "-" + to_string(tmp4.y) + "-" + to_string(tmp4.z)<<endl;
				cout<<"ErrorH"+to_string(j)<<endl;
				exit(1);
			}

			idx+=1;
		}
	}

	
	for(i=1;i<numbb;i++)
	{
		aat=refinestr_temp[i].aaa;
		atom1=AddAngleParamsList[aat][0][0]; atom2=AddAngleParamsList[aat][0][1]; atom3=AddAngleParamsList[aat][0][2];
		tmp1=refinestr_temp[i-1].ptv[atom1]; tmp2=refinestr_temp[i].ptv[atom2]; tmp3=refinestr_temp[i].ptv[atom3];
		coord2ang(tmp1, tmp2, tmp3, ang);
		if(abs(refinestr_temp[i].addang[0]-ang) > 1){
			cout<<i<<endl;
			cout<<to_string(refinestr_temp[i].addang[0]) + "-" + to_string(ang)<<endl;
			cout<<"ErrorD"+to_string(j)<<endl;
			exit(1);
		}

		for(j=1;j<AddAngleParamsList[aat].size();j++)
		{
			atom1=AddAngleParamsList[aat][j][0]; atom2=AddAngleParamsList[aat][j][1]; atom3=AddAngleParamsList[aat][j][2];
			
			tmp1=refinestr_temp[i].ptv[atom1]; tmp2=refinestr_temp[i].ptv[atom2]; tmp3=refinestr_temp[i].ptv[atom3];
			coord2ang(tmp1, tmp2, tmp3, ang);
			if(abs(refinestr_temp[i].addang[j]-ang) > 1){
				cout<<to_string(i) + "-" + to_string(j) + "-" + to_string(refinestr_temp[i].addang[j]) + "-" + to_string(ang)<<endl;
				cout<<"ErrorD"+to_string(j)<<endl;
				exit(1);
			}

		}

		for(j=0;j<4;j++)
		{
			atom1=AddLengthParamsList[aat][j][0]; atom2=AddLengthParamsList[aat][j][1];
			tmp1=refinestr_temp[i-1].ptv[atom1]; tmp2=refinestr_temp[i].ptv[atom2];
			coord2len(tmp1, tmp2, len);
			if(abs(refinestr_temp[i].addlen[j]-len) > 0.1){
				cout<<i<<endl;
				cout<<to_string(refinestr_temp[i].addlen[j]) + "-" + to_string(len)<<endl;
				cout<<"ErrorD"+to_string(j)<<endl;
				exit(1);
			}
		}

		for(j=4;j<AddLengthParamsList[aat].size();j++)
		{
			atom1=AddLengthParamsList[aat][j][0]; atom2=AddLengthParamsList[aat][j][1];
			tmp1=refinestr_temp[i].ptv[atom1]; tmp2=refinestr_temp[i].ptv[atom2];
			coord2len(tmp1, tmp2, len);
			if(abs(refinestr_temp[i].addlen[j]-len) > 0.1){
				cout<<i<<endl;
				cout<<to_string(refinestr_temp[i].addlen[j]) + "-" + to_string(len)<<endl;
				cout<<"ErrorD"+to_string(j)<<endl;
				exit(1);
			}
		}

		for(j=0;j<AddHydLengthParamsList[aat].size();j++)
		{
			atom1=AddHydLengthParamsList[aat][j][0]; atom2=AddHydLengthParamsList[aat][j][1];
			tmp1=refinestr_temp[i].ptv[atom1]; tmp2=refinestr_temp[i].ptv[atom2];
			coord2len(tmp1, tmp2, len);

			if(abs(refinestr_temp[i].addhydlen[j]-len) > 0.1){
				cout<<i<<endl;
				cout<<to_string(refinestr_temp[i].addhydlen[j]) + "-" + to_string(len)<<endl;
				cout<<"ErrorHADD_length"+to_string(j)<<endl;
				exit(1);
			}
		}

		for(j=0;j<AddHydAngleParamsList[aat].size();j++)
		{
			atom1=AddHydAngleParamsList[aat][j][0]; atom2=AddHydAngleParamsList[aat][j][1]; atom3=AddHydAngleParamsList[aat][j][2];
			tmp1=refinestr_temp[i].ptv[atom1]; tmp2=refinestr_temp[i].ptv[atom2]; tmp3=refinestr_temp[i].ptv[atom3];
			coord2ang(tmp1, tmp2, tmp3, ang);

			if(abs(refinestr_temp[i].addhydang[j]-ang) > 1){
				cout<<to_string(i) + "-" + to_string(j) + "-" + to_string(refinestr_temp[i].addhydang[j]) + "-" + to_string(ang)<<endl;
				cout<<"ErrorHADD_angle"+to_string(j)<<endl;
				exit(1);
			}
		}
	}
}

























































