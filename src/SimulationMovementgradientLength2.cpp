// EMRefinement.cpp: implementation of the EMRefinement class.
//
//////////////////////////////////////////////////////////////////////

//#include "MCRefinement.hpp"
#include "SimulationBase.hpp"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

void RNASimulation::calcene2(vector<Residueatoms> restr, vector<Residueatoms> decstr, Gradmain2 &gradobj, vector<pair<int, int>> edgeinfo, vector<int> edgetype, vector<int> paraminfo, vector<int> movementresidues, map<int, vector<Residueatoms>> dnttemplates, map<int, vector<Residueatoms>> bptemplates, map<int, TemplateDistMap> dnttemplateslength, map<int, TemplateDistMap> bptemplateslength, map<string, double> enecoeff){
	int i;
	int basei, basej;
	char aat1, aat2;
	int paramnum=AccAtomParamNum['A'][2]+AddAngleParamsList['A'].size();

	if(enecoeff["adjlength"] > 0 || enecoeff["ntlength"] > 0 || enecoeff["bplength"] > 0){
		gradobj.lengthobj.lengths.clear();
		gradobj.lengthobj.diffs.clear();
		gradobj.lengthobj.dists.clear();
		gradobj.lengthobj.outlineridx.clear();
		gradobj.lengthobj.lengths.resize(edgeinfo.size(), vector<vector<double>>(PAIRLENGTHTARGETMOVNUM, vector<double>(PAIRLENGTHTARGETMOVNUM, 0.0)));
		gradobj.lengthobj.diffs.resize(edgeinfo.size(), vector<vector<double>>(PAIRLENGTHTARGETMOVNUM, vector<double>(PAIRLENGTHTARGETMOVNUM, 0.0)));
		gradobj.lengthobj.dists.resize(edgeinfo.size(), vector<vector<Point3s>>(PAIRLENGTHTARGETMOVNUM, vector<Point3s>(PAIRLENGTHTARGETMOVNUM)));
		gradobj.lengthobj.outlineridx.resize(edgeinfo.size(), vector<vector<bool>>(PAIRLENGTHTARGETMOVNUM, vector<bool>(PAIRLENGTHTARGETMOVNUM, false)));
	}

	gradobj.lengthobj.countadj=0; gradobj.lengthobj.countnt=0; gradobj.lengthobj.countbp=0;
	gradobj.adjtempobj.count=0; gradobj.bptempobj.count=0; gradobj.inittempobj.count=0;
	gradobj.angleobj.count=0; gradobj.adjtemplengthobj.count=0;

	for(i=0;i<edgeinfo.size();i++){
		basei=edgeinfo[i].first; basej=edgeinfo[i].second;
		aat1 = restr[basei].aaa; aat2 = restr[basej].aaa;

		if(edgetype[i] == 1){
			if(enecoeff["adjlength"] > 0){
				calcmseadjlength(restr[basei].ptv, restr[basej].ptv, basei, basej, aat1, aat2, restr[basei].resind, restr[basej].resind, paraminfo[i], gradobj.lengthobj.lengths[i], gradobj.lengthobj.diffs[i], gradobj.lengthobj.dists[i], gradobj.lengthobj.mseadj, gradobj.lengthobj.countadj, gradobj.lengthobj.outlineridx[i], enecoeff["adjlength"], false);
			}

			/*
			if(enecoeff["adjtemplength"] > 0){
				gradobj.adjtemplengthobj.intralengths1[i].clear();
				gradobj.adjtemplengthobj.intralengths1[i].resize(PAIRLENGTHTARGETMOVNUM, vector<double>(PAIRLENGTHTARGETMOVNUM, 0.0));
				gradobj.adjtemplengthobj.intradiffs1[i].clear();
				gradobj.adjtemplengthobj.intradiffs1[i].resize(PAIRLENGTHTARGETMOVNUM, vector<double>(PAIRLENGTHTARGETMOVNUM, 0.0));
				gradobj.adjtemplengthobj.intradists1[i].clear();
				gradobj.adjtemplengthobj.intradists1[i].resize(PAIRLENGTHTARGETMOVNUM, vector<Point3s>(PAIRLENGTHTARGETMOVNUM));

				gradobj.adjtemplengthobj.intralengths2[i].clear();
				gradobj.adjtemplengthobj.intralengths2[i].resize(PAIRLENGTHTARGETMOVNUM, vector<double>(PAIRLENGTHTARGETMOVNUM, 0.0));
				gradobj.adjtemplengthobj.intradiffs2[i].clear();
				gradobj.adjtemplengthobj.intradiffs2[i].resize(PAIRLENGTHTARGETMOVNUM, vector<double>(PAIRLENGTHTARGETMOVNUM, 0.0));
				gradobj.adjtemplengthobj.intradists2[i].clear();
				gradobj.adjtemplengthobj.intradists2[i].resize(PAIRLENGTHTARGETMOVNUM, vector<Point3s>(PAIRLENGTHTARGETMOVNUM));

				calcmsecnnditmpintralength(restr, basei, dnttemplateslength[i].intradistmap[0], dnttemplateslength[i].interdistmap, gradobj.adjtemplengthobj.intralengths1[i], gradobj.adjtemplengthobj.intradiffs1[i], gradobj.adjtemplengthobj.intradists1[i], gradobj.adjtemplengthobj.mse, gradobj.adjtemplengthobj.count, enecoeff["adjtemplength"]);
				calcmsecnnditmpintralength(restr, basej, dnttemplateslength[i].intradistmap[1], dnttemplateslength[i].interdistmap, gradobj.adjtemplengthobj.intralengths2[i], gradobj.adjtemplengthobj.intradiffs2[i], gradobj.adjtemplengthobj.intradists2[i], gradobj.adjtemplengthobj.mse, gradobj.adjtemplengthobj.count, enecoeff["adjtemplength"]);
			}
			*/

			if(enecoeff["adjtemp"] > 0){
				//gradobj.adjtempobj.tmpdiff[i].resize(PAIRLENGTHTARGETMOVNUM*2, 0.0);
				//gradobj.adjtempobj.tmpdists[i].resize(PAIRLENGTHTARGETMOVNUM*2);

				//cout<<"aaa0"<<endl;
				gradobj.adjtempobj.lengths[i].clear();
				gradobj.adjtempobj.lengths[i].resize(AccAtomParamNum[aat1][2] + AccAtomParamNum[aat2][2], 0.0);
				
				gradobj.adjtempobj.dists[i].clear();
				gradobj.adjtempobj.dists[i].resize(AccAtomParamNum[aat1][2] + AccAtomParamNum[aat2][2]);

				//gradobj.adjtempobj.lengths[i].clear();
				//gradobj.adjtempobj.lengths[i].resize(PAIRLENGTHTARGETMOVNUM, vector<double>(PAIRLENGTHTARGETMOVNUM, 0.0));
				//gradobj.adjtempobj.diffs[i].clear();
				//gradobj.adjtempobj.diffs[i].resize(PAIRLENGTHTARGETMOVNUM, vector<double>(PAIRLENGTHTARGETMOVNUM, 0.0));
				//gradobj.adjtempobj.dists[i].clear();
				//gradobj.adjtempobj.dists[i].resize(PAIRLENGTHTARGETMOVNUM, vector<Point3s>(PAIRLENGTHTARGETMOVNUM));

				//gradobj.adjtempobj.intralengths1[i].resize(PAIRLENGTHTARGETMOVNUM, vector<double>(PAIRLENGTHTARGETMOVNUM, 0.0));
				//gradobj.adjtempobj.intradiffs1[i].resize(PAIRLENGTHTARGETMOVNUM, vector<double>(PAIRLENGTHTARGETMOVNUM, 0.0));
				//gradobj.adjtempobj.intradists1[i].resize(PAIRLENGTHTARGETMOVNUM, vector<Point3s>(PAIRLENGTHTARGETMOVNUM));

				//gradobj.adjtempobj.intralengths2[i].resize(PAIRLENGTHTARGETMOVNUM, vector<double>(PAIRLENGTHTARGETMOVNUM, 0.0));
				//gradobj.adjtempobj.intradiffs2[i].resize(PAIRLENGTHTARGETMOVNUM, vector<double>(PAIRLENGTHTARGETMOVNUM, 0.0));
				//gradobj.adjtempobj.intradists2[i].resize(PAIRLENGTHTARGETMOVNUM, vector<Point3s>(PAIRLENGTHTARGETMOVNUM));

				//cout<<"aaa"<<endl;
				//if(basei < basej){
				calcmsetempadj(restr[basei].ptv, restr[basej].ptv, aat1, aat2, dnttemplates[i], gradobj.adjtempobj.lengths[i], gradobj.adjtempobj.dists[i], gradobj.adjtempobj.mse, gradobj.adjtempobj.count, enecoeff["adjtemp"]);
					//calcmseditmplength(restr[basei].ptv, dnttemplates[i][0].ptv, restr[basej].ptv, dnttemplates[i][1].ptv, dnttemplateslength[i], gradobj.adjtempobj.tmpdiff[i], gradobj.adjtempobj.tmpdists[i], gradobj.adjtempobj.mse, gradobj.adjtempobj.count, enecoeff["adjtemp"]);
					//calcmseditmplength(restr[basei].ptv, dnttemplates[i][0].ptv, restr[basej].ptv, dnttemplates[i][1].ptv, dnttemplateslength[i], gradobj.adjtempobj.lengths[i], gradobj.adjtempobj.diffs[i], gradobj.adjtempobj.dists[i], gradobj.adjtempobj.mse, gradobj.adjtempobj.count, enecoeff["adjtemp"]);

					//calcmseditmpintralength(restr[basei].ptv, dnttemplateslength[i].intradistmap[0], gradobj.adjtempobj.intralengths1[i], gradobj.adjtempobj.intradiffs1[i], gradobj.adjtempobj.intradists1[i], gradobj.adjtempobj.mse, gradobj.adjtempobj.count, enecoeff["adjtemp"]);
					//calcmseditmpintralength(restr[basej].ptv, dnttemplateslength[i].intradistmap[1], gradobj.adjtempobj.intralengths2[i], gradobj.adjtempobj.intradiffs2[i], gradobj.adjtempobj.intradists2[i], gradobj.adjtempobj.mse, gradobj.adjtempobj.count, enecoeff["adjtemp"]);
				//}else{
					//calcmsetempadj(restr[basej].ptv, restr[basei].ptv, aat2, aat1, dnttemplates[i], gradobj.adjtempobj.lengths[i], gradobj.adjtempobj.dists[i], gradobj.adjtempobj.mse, gradobj.adjtempobj.count, enecoeff["adjtemp"]);
					//calcmseditmplength(restr[basej].ptv, dnttemplates[i][0].ptv, restr[basei].ptv, dnttemplates[i][1].ptv, dnttemplateslength[i], gradobj.adjtempobj.tmpdiff[i], gradobj.adjtempobj.tmpdists[i], gradobj.adjtempobj.mse, gradobj.adjtempobj.count, enecoeff["adjtemp"]);
					//calcmseditmplength(restr[basej].ptv, dnttemplates[i][0].ptv, restr[basei].ptv, dnttemplates[i][1].ptv, dnttemplateslength[i], gradobj.adjtempobj.lengths[i], gradobj.adjtempobj.diffs[i], gradobj.adjtempobj.dists[i], gradobj.adjtempobj.mse, gradobj.adjtempobj.count, enecoeff["adjtemp"]);

					//calcmseditmpintralength(restr[basei].ptv, dnttemplateslength[i].intradistmap[1], gradobj.adjtempobj.intralengths2[i], gradobj.adjtempobj.intradiffs2[i], gradobj.adjtempobj.intradists2[i], gradobj.adjtempobj.mse, gradobj.adjtempobj.count, enecoeff["adjtemp"]);
					//calcmseditmpintralength(restr[basej].ptv, dnttemplateslength[i].intradistmap[0], gradobj.adjtempobj.intralengths1[i], gradobj.adjtempobj.intradiffs1[i], gradobj.adjtempobj.intradists1[i], gradobj.adjtempobj.mse, gradobj.adjtempobj.count, enecoeff["adjtemp"]);
				//}
			}
		}else if(edgetype[i] == 0){
			if(enecoeff["ntlength"] > 0){
				calcmsentlength(restr[basei].ptv, restr[basej].ptv, aat1, paraminfo[i], gradobj.lengthobj.lengths[i], gradobj.lengthobj.diffs[i], gradobj.lengthobj.dists[i], gradobj.lengthobj.msent, gradobj.lengthobj.countnt, gradobj.lengthobj.outlineridx[i], enecoeff["ntlength"], false);
			}

			if(enecoeff["adjtemplength"] > 0){
				gradobj.adjtemplengthobj.intralengths1[i].clear();
				gradobj.adjtemplengthobj.intralengths1[i].resize(PAIRLENGTHTARGETMOVNUM, vector<double>(PAIRLENGTHTARGETMOVNUM, 0.0));
				gradobj.adjtemplengthobj.intradiffs1[i].clear();
				gradobj.adjtemplengthobj.intradiffs1[i].resize(PAIRLENGTHTARGETMOVNUM, vector<double>(PAIRLENGTHTARGETMOVNUM, 0.0));
				gradobj.adjtemplengthobj.intradists1[i].clear();
				gradobj.adjtemplengthobj.intradists1[i].resize(PAIRLENGTHTARGETMOVNUM, vector<Point3s>(PAIRLENGTHTARGETMOVNUM));
				gradobj.adjtemplengthobj.outlineridx1[i].clear();
				gradobj.adjtemplengthobj.outlineridx1[i].resize(PAIRLENGTHTARGETMOVNUM, vector<bool>(PAIRLENGTHTARGETMOVNUM, false));

				//gradobj.adjtemplengthobj.intralengths2[i].clear();
				//gradobj.adjtemplengthobj.intralengths2[i].resize(PAIRLENGTHTARGETMOVNUM, vector<double>(PAIRLENGTHTARGETMOVNUM, 0.0));
				//gradobj.adjtemplengthobj.intradiffs2[i].clear();
				//gradobj.adjtemplengthobj.intradiffs2[i].resize(PAIRLENGTHTARGETMOVNUM, vector<double>(PAIRLENGTHTARGETMOVNUM, 0.0));
				//gradobj.adjtemplengthobj.intradists2[i].clear();
				//gradobj.adjtemplengthobj.intradists2[i].resize(PAIRLENGTHTARGETMOVNUM, vector<Point3s>(PAIRLENGTHTARGETMOVNUM));

				calcmsecnnditmpintralength(restr, basei, gradobj.adjtemplengthobj.intralengths1[i], gradobj.adjtemplengthobj.intradiffs1[i], gradobj.adjtemplengthobj.intradists1[i], gradobj.adjtemplengthobj.outlineridx1[i], gradobj.adjtemplengthobj.mse, gradobj.adjtemplengthobj.count, enecoeff["adjtemplength"]);
				//calcmsecnnditmpintralength(restr, basej, dnttemplateslength[i].intradistmap[1], dnttemplateslength[i].interdistmap, gradobj.adjtemplengthobj.intralengths2[i], gradobj.adjtemplengthobj.intradiffs2[i], gradobj.adjtemplengthobj.intradists2[i], gradobj.adjtemplengthobj.mse, gradobj.adjtemplengthobj.count, enecoeff["adjtemplength"]);
			}
		}else if(edgetype[i] == 2){
			if(enecoeff["bplength"] > 0){
				calcmsebplength(restr[basei].ptv, restr[basej].ptv, aat1, aat2, gradobj.lengthobj.lengths[i], gradobj.lengthobj.diffs[i], gradobj.lengthobj.dists[i], gradobj.lengthobj.msebp, gradobj.lengthobj.countbp, gradobj.lengthobj.outlineridx[i], enecoeff["bplength"], false);
			}

			if(enecoeff["bptemp"] > 0){
				//gradobj.bptempobj.tmpdiff[i].resize(PAIRLENGTHTARGETMOVNUM*2, 0.0);
				//gradobj.bptempobj.tmpdists[i].resize(PAIRLENGTHTARGETMOVNUM*2);
				//calcmseditmplength(restr[basei].ptv, bptemplates[i][0].ptv, restr[basej].ptv, bptemplates[i][1].ptv, bptemplateslength[i], gradobj.bptempobj.tmpdiff[i], gradobj.bptempobj.tmpdists[i], gradobj.bptempobj.mse, gradobj.bptempobj.count, enecoeff["bptemp"]);

				gradobj.bptempobj.lengths[i].clear();
				gradobj.bptempobj.lengths[i].resize(AccAtomParamNum[aat1][2] + AccAtomParamNum[aat2][2], 0.0);
				gradobj.bptempobj.dists[i].clear();
				gradobj.bptempobj.dists[i].resize(AccAtomParamNum[aat1][2] + AccAtomParamNum[aat2][2]);

				//gradobj.bptempobj.lengths[i].clear();
				//gradobj.bptempobj.lengths[i].resize(PAIRLENGTHTARGETMOVNUM, vector<double>(PAIRLENGTHTARGETMOVNUM, 0.0));
				//gradobj.bptempobj.diffs[i].clear();
				//gradobj.bptempobj.diffs[i].resize(PAIRLENGTHTARGETMOVNUM, vector<double>(PAIRLENGTHTARGETMOVNUM, 0.0));
				//gradobj.bptempobj.dists[i].clear();
				//gradobj.bptempobj.dists[i].resize(PAIRLENGTHTARGETMOVNUM, vector<Point3s>(PAIRLENGTHTARGETMOVNUM));

				//gradobj.bptempobj.intralengths1[i].resize(PAIRLENGTHTARGETMOVNUM, vector<double>(PAIRLENGTHTARGETMOVNUM, 0.0));
				//gradobj.bptempobj.intradiffs1[i].resize(PAIRLENGTHTARGETMOVNUM, vector<double>(PAIRLENGTHTARGETMOVNUM, 0.0));
				//gradobj.bptempobj.intradists1[i].resize(PAIRLENGTHTARGETMOVNUM, vector<Point3s>(PAIRLENGTHTARGETMOVNUM));

				//gradobj.bptempobj.intralengths2[i].resize(PAIRLENGTHTARGETMOVNUM, vector<double>(PAIRLENGTHTARGETMOVNUM, 0.0));
				//gradobj.bptempobj.intradiffs2[i].resize(PAIRLENGTHTARGETMOVNUM, vector<double>(PAIRLENGTHTARGETMOVNUM, 0.0));
				//gradobj.bptempobj.intradists2[i].resize(PAIRLENGTHTARGETMOVNUM, vector<Point3s>(PAIRLENGTHTARGETMOVNUM));

				calcmsetempbp(restr[basei].ptv, restr[basej].ptv, aat1, aat2, bptemplates[i], gradobj.bptempobj.lengths[i], gradobj.bptempobj.dists[i], gradobj.bptempobj.mse, gradobj.bptempobj.count, enecoeff["bptemp"]);
					
				//calcmseditmplength(restr[basei].ptv, bptemplates[i][0].ptv, restr[basej].ptv, bptemplates[i][1].ptv, bptemplateslength[i], gradobj.bptempobj.lengths[i], gradobj.bptempobj.diffs[i], gradobj.bptempobj.dists[i], gradobj.bptempobj.mse, gradobj.bptempobj.count, enecoeff["bptemp"]);
			
				//calcmseditmpintralength(restr[basei].ptv, bptemplateslength[i].intradistmap[0], gradobj.bptempobj.intralengths1[i], gradobj.bptempobj.intradiffs1[i], gradobj.bptempobj.intradists1[i], gradobj.bptempobj.mse, gradobj.bptempobj.count, enecoeff["bptemp"]);
				//calcmseditmpintralength(restr[basej].ptv, bptemplateslength[i].intradistmap[1], gradobj.bptempobj.intralengths2[i], gradobj.bptempobj.intradiffs2[i], gradobj.bptempobj.intradists2[i], gradobj.bptempobj.mse, gradobj.bptempobj.count, enecoeff["bptemp"]);
			}
		}
	}
	//exit(1);

	if(enecoeff["inittemp"] > 0) gradobj.inittempobj.count=0; 
	if(enecoeff["angle"] > 0) gradobj.angleobj.count=0;
	for(i=0;i<movementresidues.size();i++)
	{
		basei = movementresidues[i];

		if(enecoeff["inittemp"] > 0)
		{
			gradobj.inittempobj.tmpdiff[basei].clear();
			gradobj.inittempobj.tmpdiff[basei].resize(PAIRLENGTHTARGETMOVNUM, 0.0);
			gradobj.inittempobj.tmpdists[basei].clear();
			gradobj.inittempobj.tmpdists[basei].resize(PAIRLENGTHTARGETMOVNUM);
			calcmseinitlength(restr[basei].ptv, decstr[basei].ptv, gradobj.inittempobj.tmpdiff[basei], gradobj.inittempobj.tmpdists[basei], gradobj.inittempobj.mse, gradobj.inittempobj.count, enecoeff["inittemp"]);
		
		}

		if(enecoeff["angle"] > 0)
		{
			gradobj.angleobj.angles[basei].clear();
			gradobj.angleobj.angles[basei].resize(paramnum, 0.0);
			gradobj.angleobj.anglediffs[basei].clear();
			gradobj.angleobj.anglediffs[basei].resize(paramnum, 0.0);
			gradobj.angleobj.angleouterflags[basei].clear();
			gradobj.angleobj.angleouterflags[basei].resize(paramnum, false);
			gradobj.angleobj.theta[basei].clear();
			gradobj.angleobj.theta[basei].resize(paramnum, 0.0);
			gradobj.angleobj.normPQs[basei].clear();
			gradobj.angleobj.normPQs[basei].resize(paramnum, 0.0);
			gradobj.angleobj.normRQs[basei].clear();
			gradobj.angleobj.normRQs[basei].resize(paramnum, 0.0);
			gradobj.angleobj.dotps[basei].clear();
			gradobj.angleobj.dotps[basei].resize(paramnum, 0.0);
			gradobj.angleobj.vecPQs[basei].clear();
			gradobj.angleobj.vecPQs[basei].resize(paramnum);
			gradobj.angleobj.vecRQs[basei].clear();
			gradobj.angleobj.vecRQs[basei].resize(paramnum);
			calcmseangles(restr, basei, gradobj.angleobj.angles[basei], gradobj.angleobj.anglediffs[basei], gradobj.angleobj.angleouterflags[basei], gradobj.angleobj.theta[basei], gradobj.angleobj.normPQs[basei], gradobj.angleobj.normRQs[basei], gradobj.angleobj.dotps[basei], gradobj.angleobj.vecPQs[basei], gradobj.angleobj.vecRQs[basei], gradobj.angleobj.mse, gradobj.angleobj.count, enecoeff["angle"]);
		}
	}

	if(gradobj.lengthobj.countadj > 0 && enecoeff["adjlength"] > 0) gradobj.lengthobj.mseadj=gradobj.lengthobj.mseadj/gradobj.lengthobj.countadj;
	if(gradobj.lengthobj.countnt > 0 && enecoeff["ntlength"] > 0) gradobj.lengthobj.msent=gradobj.lengthobj.msent/gradobj.lengthobj.countnt;
	if(gradobj.lengthobj.countbp > 0 && enecoeff["bplength"] > 0) gradobj.lengthobj.msebp=gradobj.lengthobj.msebp/gradobj.lengthobj.countbp;
	if(gradobj.inittempobj.count > 0 && enecoeff["inittemp"] > 0) gradobj.inittempobj.mse=gradobj.inittempobj.mse/gradobj.inittempobj.count;
	if(gradobj.adjtempobj.count > 0 && enecoeff["adjtemp"] > 0) gradobj.adjtempobj.mse=gradobj.adjtempobj.mse/gradobj.adjtempobj.count;
	if(gradobj.adjtemplengthobj.count > 0 && enecoeff["adjtemplength"] > 0) gradobj.adjtemplengthobj.mse=gradobj.adjtemplengthobj.mse/gradobj.adjtemplengthobj.count;
	if(gradobj.bptempobj.count > 0 && enecoeff["bptemp"] > 0) gradobj.bptempobj.mse=gradobj.bptempobj.mse/gradobj.bptempobj.count;
	if(gradobj.angleobj.count > 0 && enecoeff["angle"] > 0) gradobj.angleobj.mse=gradobj.angleobj.mse/gradobj.angleobj.count;
	gradobj.mse=gradobj.lengthobj.mseadj+gradobj.lengthobj.msent+gradobj.lengthobj.msebp+gradobj.inittempobj.mse+gradobj.adjtempobj.mse+gradobj.bptempobj.mse+gradobj.angleobj.mse;
}

void RNASimulation::calcgrad2(vector<Residueatoms> restr, vector<Residueatoms> decstr, map<int, vector<Point3s>> &gradptv, Gradmain2 gradobj , vector<pair<int, int>> edgeinfo, vector<int> edgetype, vector<bool> movenodeinfo, vector<int> movementresidues, map<string, double> enecoeff)
{
	int i;
	int basei, basej;
	char aat1, aat2;

	for(i=0;i<edgeinfo.size();i++){
		basei=edgeinfo[i].first; basej=edgeinfo[i].second;
		aat1 = restr[basei].aaa; aat2 = restr[basej].aaa;
		if(edgetype[i] == 1){

			if(enecoeff["adjlength"] > 0 && gradobj.lengthobj.countadj > 0){
				calcgradientlength(restr[basei].ptv, restr[basej].ptv, basei, basej, gradptv, movenodeinfo, gradobj.lengthobj.lengths[i], gradobj.lengthobj.diffs[i], gradobj.lengthobj.dists[i], gradobj.lengthobj.countadj, gradobj.lengthobj.outlineridx[i], enecoeff["adjlength"]);
			}

			/*
			if(enecoeff["adjtemplength"] > 0 && gradobj.adjtemplengthobj.count > 0){
				//if(basei < basej){
					//calcgradientditmplength(basei, basej, gradptv, movenodeinfo, gradobj.adjtemplengthobj.lengths[i], gradobj.adjtemplengthobj.diffs[i], gradobj.adjtemplengthobj.dists[i], gradobj.adjtemplengthobj.count, enecoeff["adjtemplength"]);
					//calcgradientditmplength(basei, gradptv, movenodeinfo, gradobj.adjtemplengthobj.intralengths1[i], gradobj.adjtemplengthobj.intradiffs1[i], gradobj.adjtemplengthobj.intradists1[i], gradobj.adjtemplengthobj.count, enecoeff["adjtemplength"]);
					//calcgradientditmplength(basej, gradptv, movenodeinfo, gradobj.adjtemplengthobj.intralengths2[i], gradobj.adjtemplengthobj.intradiffs2[i], gradobj.adjtemplengthobj.intradists2[i], gradobj.adjtemplengthobj.count, enecoeff["adjtemplength"]);
				
					calcgracnnentditmpintra(basei, gradptv, movenodeinfo, gradobj.adjtemplengthobj.intralengths1[i], gradobj.adjtemplengthobj.intradiffs1[i], gradobj.adjtemplengthobj.intradists1[i], gradobj.adjtemplengthobj.count, enecoeff["adjtemplength"]);
					calcgracnnentditmpintra(basej, gradptv, movenodeinfo, gradobj.adjtemplengthobj.intralengths2[i], gradobj.adjtemplengthobj.intradiffs2[i], gradobj.adjtemplengthobj.intradists2[i], gradobj.adjtemplengthobj.count, enecoeff["adjtemplength"]);

				//}else{
					//calcgradientditmplength(basej, basei, gradptv, movenodeinfo, gradobj.adjtemplengthobj.lengths[i], gradobj.adjtemplengthobj.diffs[i], gradobj.adjtemplengthobj.dists[i], gradobj.adjtemplengthobj.count, enecoeff["adjtemplength"]);
					//calcgradientditmplength(basei, gradptv, movenodeinfo, gradobj.adjtemplengthobj.intralengths2[i], gradobj.adjtemplengthobj.intradiffs2[i], gradobj.adjtemplengthobj.intradists2[i], gradobj.adjtemplengthobj.count, enecoeff["adjtemplength"]);
					//calcgradientditmplength(basej, gradptv, movenodeinfo, gradobj.adjtemplengthobj.intralengths1[i], gradobj.adjtemplengthobj.intradiffs1[i], gradobj.adjtemplengthobj.intradists1[i], gradobj.adjtemplengthobj.count, enecoeff["adjtemplength"]);
				
				//}
			}
			*/
			
			if(enecoeff["adjtemp"] > 0 && gradobj.adjtempobj.count > 0){
				//if(basei < basej){
				calcgradienttempadj(basei, basej, aat1, aat2, gradptv, gradobj.adjtempobj.lengths[i], gradobj.adjtempobj.dists[i], gradobj.adjtempobj.count, enecoeff["adjtemp"]);
					//calcgradientditmplength(basei, basej, gradptv, gradobj.adjtempobj.tmpdiff[i], gradobj.adjtempobj.tmpdists[i], gradobj.adjtempobj.count, enecoeff["adjtemp"]);
					//calcgradientditmplength(basei, basej, gradptv, movenodeinfo, gradobj.adjtempobj.lengths[i], gradobj.adjtempobj.diffs[i], gradobj.adjtempobj.dists[i], gradobj.adjtempobj.count, enecoeff["adjtemp"]);
					//calcgradientditmplength(basei, gradptv, movenodeinfo, gradobj.adjtempobj.intralengths1[i], gradobj.adjtempobj.intradiffs1[i], gradobj.adjtempobj.intradists1[i], gradobj.adjtempobj.count, enecoeff["adjtemp"]);
					//calcgradientditmplength(basej, gradptv, movenodeinfo, gradobj.adjtempobj.intralengths2[i], gradobj.adjtempobj.intradiffs2[i], gradobj.adjtempobj.intradists2[i], gradobj.adjtempobj.count, enecoeff["adjtemp"]);
				//}else{
					//calcgradienttempadj(basej, basei, aat2, aat1, gradptv, gradobj.adjtempobj.lengths[i], gradobj.adjtempobj.dists[i], gradobj.adjtempobj.count, enecoeff["adjtemp"]);
					//calcgradientditmplength(basej, basei, gradptv, gradobj.adjtempobj.tmpdiff[i], gradobj.adjtempobj.tmpdists[i], gradobj.adjtempobj.count, enecoeff["adjtemp"]);
					//calcgradientditmplength(basej, basei, gradptv, movenodeinfo, gradobj.adjtempobj.lengths[i], gradobj.adjtempobj.diffs[i], gradobj.adjtempobj.dists[i], gradobj.adjtempobj.count, enecoeff["adjtemp"]);
					//calcgradientditmplength(basei, gradptv, movenodeinfo, gradobj.adjtempobj.intralengths2[i], gradobj.adjtempobj.intradiffs2[i], gradobj.adjtempobj.intradists2[i], gradobj.adjtempobj.count, enecoeff["adjtemp"]);
					//calcgradientditmplength(basej, gradptv, movenodeinfo, gradobj.adjtempobj.intralengths1[i], gradobj.adjtempobj.intradiffs1[i], gradobj.adjtempobj.intradists1[i], gradobj.adjtempobj.count, enecoeff["adjtemp"]);
				//}
			}
		}else if(edgetype[i] == 0){
			if(enecoeff["ntlength"] > 0 && gradobj.lengthobj.countnt > 0){
				calcgradientlength(restr[basei].ptv, restr[basej].ptv, basei, basej, gradptv, movenodeinfo, gradobj.lengthobj.lengths[i], gradobj.lengthobj.diffs[i], gradobj.lengthobj.dists[i], gradobj.lengthobj.countnt, gradobj.lengthobj.outlineridx[i], enecoeff["ntlength"]);
			}

			if(enecoeff["adjtemplength"] > 0 && gradobj.adjtemplengthobj.count > 0){
				calcgracnnentditmpintra(basei, gradptv, movenodeinfo, gradobj.adjtemplengthobj.intralengths1[i], gradobj.adjtemplengthobj.intradiffs1[i], gradobj.adjtemplengthobj.intradists1[i], gradobj.adjtemplengthobj.outlineridx1[i], gradobj.adjtemplengthobj.count, enecoeff["adjtemplength"]);
				//calcgracnnentditmpintra(basej, gradptv, movenodeinfo, gradobj.adjtemplengthobj.intralengths2[i], gradobj.adjtemplengthobj.intradiffs2[i], gradobj.adjtemplengthobj.intradists2[i], gradobj.adjtemplengthobj.count, enecoeff["adjtemplength"]);
			}

		}else if(edgetype[i] == 2){
			if(enecoeff["bplength"] > 0 && gradobj.lengthobj.countbp > 0){
				calcgradientlength(restr[basei].ptv, restr[basej].ptv, basei, basej, gradptv, movenodeinfo, gradobj.lengthobj.lengths[i], gradobj.lengthobj.diffs[i], gradobj.lengthobj.dists[i], gradobj.lengthobj.countbp, gradobj.lengthobj.outlineridx[i], enecoeff["bplength"]);
			}

			if(enecoeff["bptemp"] > 0 && gradobj.bptempobj.count > 0){
				calcgradienttempbp(basei, basej, aat1, aat2, gradptv, gradobj.bptempobj.lengths[i], gradobj.bptempobj.dists[i], gradobj.bptempobj.count, enecoeff["bptemp"]);
					
				//calcgradientditmplength(basei, basej, gradptv, gradobj.bptempobj.tmpdiff[i], gradobj.bptempobj.tmpdists[i], gradobj.bptempobj.count, enecoeff["bptemp"]);
				//calcgradientditmplength(basei, basej, gradptv, movenodeinfo, gradobj.bptempobj.lengths[i], gradobj.bptempobj.diffs[i], gradobj.bptempobj.dists[i], gradobj.bptempobj.count, enecoeff["bptemp"]);
			
				//calcgradientditmplength(basei, gradptv, movenodeinfo, gradobj.bptempobj.intralengths1[i], gradobj.bptempobj.intradiffs1[i], gradobj.bptempobj.intradists1[i], gradobj.bptempobj.count, enecoeff["bptemp"]);
				//calcgradientditmplength(basej, gradptv, movenodeinfo, gradobj.bptempobj.intralengths2[i], gradobj.bptempobj.intradiffs2[i], gradobj.bptempobj.intradists2[i], gradobj.bptempobj.count, enecoeff["bptemp"]);
			}
		}
	}
	for(i=0;i<movementresidues.size();i++)
	{
		basei = movementresidues[i];

		if(enecoeff["inittemp"] > 0 && gradobj.inittempobj.count > 0){ 
			calcgradientinitlength(restr[basei].ptv, decstr[basei].ptv, basei, gradptv, gradobj.inittempobj.tmpdiff[basei], gradobj.inittempobj.tmpdists[basei], gradobj.inittempobj.count, enecoeff["inittemp"]);
		}

		if(enecoeff["angle"] > 0 && gradobj.angleobj.count > 0){
			calcgradientangles(restr, basei, gradptv, movenodeinfo, gradobj.angleobj.anglediffs[basei], gradobj.angleobj.angleouterflags[basei], gradobj.angleobj.theta[basei], gradobj.angleobj.normPQs[basei], gradobj.angleobj.normRQs[basei], gradobj.angleobj.dotps[basei], gradobj.angleobj.vecPQs[basei], gradobj.angleobj.vecRQs[basei], gradobj.angleobj.count, enecoeff["angle"]);
		}
	}
}

void RNASimulation::initgrad2obj(Gradmain2 &gradobj, int index){
	if(index == 0){
		gradobj.mse=0;
		gradobj.lengthobj.mseadj=0;
		gradobj.lengthobj.msent=0;
		gradobj.lengthobj.msebp=0;
		gradobj.adjtempobj.mse=0;
		gradobj.adjtemplengthobj.mse=0;
		gradobj.inittempobj.mse=0;
		gradobj.bptempobj.mse=0;
		gradobj.angleobj.mse=0;

		gradobj.mseprev=0;
		gradobj.lengthobj.mseadjprev=0;
		gradobj.lengthobj.msentprev=0;
		gradobj.lengthobj.msebpprev=0;
		gradobj.adjtempobj.mseprev=0;
		gradobj.adjtemplengthobj.mseprev=0;
		gradobj.inittempobj.mseprev=0;
		gradobj.bptempobj.mseprev=0;
		gradobj.angleobj.mseprev=0;

		gradobj.lengthobj.countadj=0;
		gradobj.lengthobj.countnt=0;
		gradobj.lengthobj.countbp=0;
		gradobj.adjtempobj.count=0;
		gradobj.adjtemplengthobj.count=0;
		gradobj.inittempobj.count=0;
		gradobj.bptempobj.count=0;
		gradobj.angleobj.count=0;
	}else if(index == 1){
		gradobj.mse=0;
		gradobj.lengthobj.mseadj=0;
		gradobj.lengthobj.msent=0;
		gradobj.lengthobj.msebp=0;
		gradobj.adjtempobj.mse=0;
		gradobj.adjtemplengthobj.mse=0;
		gradobj.inittempobj.mse=0;
		gradobj.bptempobj.mse=0;
		gradobj.angleobj.mse=0;
	}else if(index == 2){
		gradobj.mseprev=0;
		gradobj.lengthobj.mseadjprev=0;
		gradobj.lengthobj.msentprev=0;
		gradobj.lengthobj.msebpprev=0;
		gradobj.adjtempobj.mseprev=0;
		gradobj.adjtemplengthobj.mseprev=0;
		gradobj.inittempobj.mseprev=0;
		gradobj.bptempobj.mseprev=0;
		gradobj.angleobj.mseprev=0;
	}else if(index == 3){
		gradobj.mseprev=gradobj.mse;
		gradobj.lengthobj.mseadjprev=gradobj.lengthobj.mseadj;
		gradobj.lengthobj.msentprev=gradobj.lengthobj.msent;
		gradobj.lengthobj.msebpprev=gradobj.lengthobj.msebp;
		gradobj.inittempobj.mseprev=gradobj.inittempobj.mse;
		gradobj.adjtempobj.mseprev=gradobj.adjtempobj.mse;
		gradobj.adjtemplengthobj.mseprev=gradobj.adjtemplengthobj.mse;
		gradobj.bptempobj.mseprev=gradobj.bptempobj.mse;
		gradobj.angleobj.mseprev=gradobj.angleobj.mse;
	}
}

double RNASimulation::lbfgsgradnorm(map<int, vector<Point3s>> grad, vector<bool> movenodeinfo){
	double normval=0;

	for(int i=0;i<movenodeinfo.size();i++){
		if(movenodeinfo[i]){
			for(int mtype=0;mtype<PAIRLENGTHTARGETMOVNUM;mtype++)
			{
				normval+=sqrt(grad[i][mtype].x*grad[i][mtype].x);
				normval+=sqrt(grad[i][mtype].y*grad[i][mtype].y);
				normval+=sqrt(grad[i][mtype].z*grad[i][mtype].z);
			}
		}
	}

	return normval;
}

void RNASimulation::lbfgsupdate(map<int, vector<Point3s>> grad, deque<map<int, vector<Point3s>>> s_list, deque<map<int, vector<Point3s>>> y_list, deque<double> rho_list, vector<bool> movenodeinfo, map<int, vector<Point3s>> &z){
	map<int, vector<Point3s>> q=grad;
	size_t m = s_list.size();

	vector<double> alpha;
	alpha.resize(s_list.size());

	for (int i = m - 1; i >= 0; --i) {
		map<int, vector<Point3s>> s = s_list[i];
        map<int, vector<Point3s>> y = y_list[i];
        double rho = rho_list[i];

		double alpha_i = 0.0;
		for (auto& [key, s_vec] : s) {
			if(movenodeinfo[key]){
				for (size_t j = 0; j < s_vec.size(); ++j) {
					alpha_i += rho * (s_vec[j].x*q[key][j].x);
					alpha_i += rho * (s_vec[j].y*q[key][j].y);
					alpha_i += rho * (s_vec[j].z*q[key][j].z);
				}
			}
		}

		alpha[i] = alpha_i;

		for (auto& [key, y_vec] : y) {
			if(movenodeinfo[key]){
				for (size_t j = 0; j < y_vec.size(); ++j) {
					q[key][j].x = q[key][j].x - alpha_i * y_vec[j].x;
					q[key][j].y = q[key][j].y - alpha_i * y_vec[j].y;
					q[key][j].z = q[key][j].z - alpha_i * y_vec[j].z;
				}
			}
		}
	}

	if(s_list.size() > 0){
		double gamma_num = 0.0;
		double gamma_den = 0.0;
		for (auto& [key, s_vec] : s_list[m-1]) {
			if(movenodeinfo[key]){
				for (size_t j = 0; j < s_vec.size(); ++j) {
					gamma_num += s_vec[j].x * y_list[m-1][key][j].x;
					gamma_num += s_vec[j].y * y_list[m-1][key][j].y;
					gamma_num += s_vec[j].z * y_list[m-1][key][j].z;
					gamma_den += y_list[m-1][key][j].x * y_list[m-1][key][j].x;
					gamma_den += y_list[m-1][key][j].y * y_list[m-1][key][j].y;
					gamma_den += y_list[m-1][key][j].z * y_list[m-1][key][j].z;
				}
			}
		}
		double gamma_k = gamma_num / (gamma_den+EPS);

		for (auto& [key, q_vec] : q) 
		{
			if(movenodeinfo[key]){
				for (size_t j = 0; j < q_vec.size(); ++j)
				{
					q[key][j].x = gamma_k * q[key][j].x;
					q[key][j].y = gamma_k * q[key][j].y;
					q[key][j].z = gamma_k * q[key][j].z;
				}
			}
		}
	}
	z = q;

	for (size_t i = 0; i < s_list.size(); ++i) {
		map<int, vector<Point3s>> s = s_list[i];
    	map<int, vector<Point3s>> y = y_list[i];
    	double rho = rho_list[i];

		double beta_i = 0.0;
    	for (auto& [key, y_vec] : y) {
			if(movenodeinfo[key]){
				for (size_t j = 0; j < y_vec.size(); ++j) {
					beta_i += rho * (y_vec[j].x * z[key][j].x);
					beta_i += rho * (y_vec[j].y * z[key][j].y);
					beta_i += rho * (y_vec[j].z * z[key][j].z);
				}
			}
    	}

		for (auto& [key, s_vec] : s) {
			if(movenodeinfo[key]){
				for (size_t j = 0; j < s_vec.size(); ++j) {
					z[key][j].x = z[key][j].x + s_vec[j].x * (alpha[i] - beta_i);
					z[key][j].y = z[key][j].y + s_vec[j].y * (alpha[i] - beta_i);
					z[key][j].z = z[key][j].z + s_vec[j].z * (alpha[i] - beta_i);
				}
			}
    	}
	}

	for (auto& [key, z_vec] : z) {
		if(movenodeinfo[key]){
			for (size_t j = 0; j < z_vec.size(); ++j) {
				z[key][j].x = -z[key][j].x;
				z[key][j].y = -z[key][j].y;
				z[key][j].z = -z[key][j].z;
			}
		}
	}
}

double RNASimulation::getgraddotproduct(map<int, vector<Point3s>> grad, map<int, vector<Point3s>> z, vector<bool> movenodeinfo)
{
	int i, j;
	double dot_product=0.0;
	for(i=0;i<movenodeinfo.size();i++){
		if(movenodeinfo[i]){
			for(j=0;j<PAIRLENGTHTARGETMOVNUM;j++)
			{
				dot_product+= grad[i][j].x*z[i][j].x;
				dot_product+= grad[i][j].y*z[i][j].y;
				dot_product+= grad[i][j].z*z[i][j].z;
			}
		}
	}

	return dot_product;
}

void RNASimulation::strupdate(vector<Residueatoms> &restr, double step_size, map<int, vector<Point3s>> z, vector<bool> movenodeinfo)
{
	int i, j;
	for(i=0;i<movenodeinfo.size();i++){
		if(movenodeinfo[i]){
			for(j=0;j<PAIRLENGTHTARGETMOVNUM;j++)
			{
				restr[i].ptv[j].x += step_size*z[i][j].x;
				restr[i].ptv[j].y += step_size*z[i][j].y;
				restr[i].ptv[j].z += step_size*z[i][j].z;
			}
		}
	}
}

double RNASimulation::linesearch(vector<Residueatoms>& restr, vector<Residueatoms> decstr, map<int, vector<Point3s>>& grad, map<int, vector<Point3s>> z, Gradmain2 &gradobj, vector<pair<int, int>> edgeinfo, vector<int> edgetype, vector<int> paraminfo, vector<bool> movenodeinfo, vector<int> movementresidues, map<int, vector<Residueatoms>> dnttemplates, map<int, vector<Residueatoms>> bptemplates, map<int, TemplateDistMap> dnttemplateslength, map<int, TemplateDistMap> bptemplateslength, map<string, double> enecoeff, double alpha, double beta, double gamma)
{
	int i, mtype;
	double step_size = alpha;
    double old_energy, new_energy, grad_norm;
	const double tolerance = 1e-7;
	vector<Residueatoms> restr_copy = restr;
	map<int, vector<Point3s>> grad_copy = grad;
	double armijo_condition;
	double dot_product, new_dot_product;

	initgrad2obj(gradobj, 1);
	calcene2(restr, decstr, gradobj, edgeinfo, edgetype, paraminfo, movementresidues, dnttemplates, bptemplates, dnttemplateslength, bptemplateslength, enecoeff);
	old_energy = gradobj.mse;
	initgrad(grad, movenodeinfo);
	calcgrad2(restr, decstr, grad, gradobj, edgeinfo, edgetype, movenodeinfo, movementresidues, enecoeff);

	dot_product=getgraddotproduct(grad, z, movenodeinfo);

	//const int max_iterations = 100; // 例: 最大50回
	//int iter_count = 0;

	double grad_norm_init = lbfgsgradnorm(grad, movenodeinfo);
	double curvature_condition = beta * grad_norm_init;

	while (true) {
		restr_copy = restr;
		grad_copy = grad;

		strupdate(restr_copy, step_size, z, movenodeinfo);

		initgrad2obj(gradobj, 1);
		calcene2(restr_copy, decstr, gradobj, edgeinfo, edgetype, paraminfo, movementresidues, dnttemplates, bptemplates, dnttemplateslength, bptemplateslength, enecoeff);
		new_energy = gradobj.mse;
		initgrad(grad_copy, movenodeinfo);
		calcgrad2(restr_copy, decstr, grad_copy, gradobj, edgeinfo, edgetype, movenodeinfo, movementresidues, enecoeff);
		grad_norm = lbfgsgradnorm(grad_copy, movenodeinfo);

		new_dot_product=getgraddotproduct(grad_copy, z, movenodeinfo);

		armijo_condition = old_energy + gamma * step_size * dot_product;

		if (new_energy < armijo_condition){
			break;
		}

		step_size *= beta;
		if (step_size < tolerance) {
            break;
        }

		//iter_count+=1;
		//if (iter_count >= max_iterations) break;
	}

	restr=restr_copy;
	grad=grad_copy;
	return step_size;
}

bool RNASimulation::fmovementgradientlength2(vector<Residueatoms> &decstr,int numseq,double oldenergy,double &newenergy,double Bt,bool flagacc)
{
	int i,j,k,inext;
	int istart,iend;
	int trandpos, trandposed, trandposst, span;
	int basei, basej;
	double trandang;
	Point3d rp;
	double mati[9];
	int wpos;
	char aat1, aat2, aat;
	bool idxtmp1, idxtmp2;
	int transbase;
	int total;
	int numiter;
	double rmax, rnorm;
	int mtype1, mtype2;
	vector<int> ntparamidx;
	vector<int> adjparamidx;
	bool flagdone;
	bool exitOuterLoop = false;
	int mtype;
	int rtype;
	int pairwiseidx;
	vector<bool> movenodeinfo;
	vector<pair<int, int>> edgeinfo;
	vector<int> edgetype;
	vector<int> paraminfo;
	map<int, vector<Point3s>> snttemplatestr;

	int epochs;
	double tolerance;
	int tmp;
	double m;

	double randval=Random();

	
	map<string, double> enecoeff;
	
	/*
	if(randval < 0.25){
		enecoeff={
			{"inittemp", 0},
			{"bptemp", 0.01},
			{"adjtemp", 0.5},
			{"adjlength", 100.0},
			{"bplength", 0},
			{"ntlength", 100.0},
			{"angle", 0},
			{"adjtemplength", 10}
		};
	}else if(0.25<=randval && 0.50<randval){
		enecoeff={
			{"inittemp", 0},
			{"bptemp", 0.1},
			{"adjtemp", 5.0},
			{"adjlength", 100.0},
			{"bplength", 0},
			{"ntlength", 100.0},
			{"angle", 0},
			{"adjtemplength", 100}
		};
	}else if(0.50<=randval && 0.75<randval){
		enecoeff={
			{"inittemp", 0},
			{"bptemp", 0.1},
			{"adjtemp", 10.0},
			{"adjlength", 100.0},
			{"bplength", 0},
			{"ntlength", 100.0},
			{"angle", 0},
			{"adjtemplength", 300}
		};
	}else{
		enecoeff={
			{"inittemp", 0},
			{"bptemp", 0},
			{"adjtemp", 4.0},
			{"adjlength", 10.0},
			{"bplength", 0},
			{"ntlength", 10.0},
			{"angle", 0},
			{"adjtemplength", 10.0}
		};
	}
	*/

	enecoeff={
			{"inittemp", 0},
			{"bptemp", 1.0},
			{"adjtemp", 4.0},
			{"adjlength", 100.0},
			{"bplength", 0},
			{"ntlength", 100.0},
			{"angle", 0},
			{"adjtemplength", 100.0}
		};

	//epochs=100+100*Random();
	epochs=50;
	tolerance=0.000001;
	m = 15;
	double alpha = 5.0;
	double beta = 0.8;
	double gamma = 0.0001;

	randval=Random();

	int totwrong1=lp1.nn[CLASH_ENERGY_ST]*2;
	int totwrong2=lp1.nn[BASE_ENERGY_ST+0] + lp1.nn[BASE_ENERGY_ST+1];
	int totwrong3=0;
	for(i=0;i<BLENGTHTARGET.size();i++){
		totwrong3+=lp1.nn[BLENGTH_ENERGY_ST+BLENGTHTARGET[i]];
		//totwrong3-=lp1.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[i]][0];
	}
	for(i=0;i<INNERANGLETARGET.size();i++){
		totwrong3+=lp1.nn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[i]];
		//totwrong3-=lp1.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[i]][0];
	}

	if(totwrong3>0 && randval<0.5){
		wpos=int(Random()*totwrong3-EPS);

		i=0;
		for(trandpos=0;trandpos<numseq;trandpos++)
		{
			for(j=0;j<BLENGTHTARGET.size();j++){
				inext=i+lp1.indn[BLENGTH_ENERGY_ST+BLENGTHTARGET[j]][trandpos];
				if(i<=wpos && wpos<inext){
					exitOuterLoop = true;
					break;
				}
				i=inext;
			}

			for(j=0;j<INNERANGLETARGET.size();j++){
				inext=i+lp1.indn[INNERANGLE_ENERGY_ST+INNERANGLETARGET[j]][trandpos];
				if(i<=wpos && wpos<inext){
					exitOuterLoop = true;
					break;
				}
				i=inext;
			}
			if(exitOuterLoop) break;
		}
	}else if(totwrong1>0 && randval>=0.5 && randval<0.7){
		wpos=int(Random()*totwrong1-EPS);

		i=0;
		for(trandpos=0;trandpos<numseq;trandpos++)
		{
			inext=i+lp1.indn[CLASH_ENERGY_ST][trandpos];
			if(i<=wpos && wpos<inext){
				break;
			}
			i=inext;
		}
	}else if(totwrong2>0 && randval>=0.7 && randval<0.9){
		wpos=int(Random()*totwrong2-EPS);

		i=0;
		for(trandpos=0;trandpos<BASEPAIRS.size();trandpos++)
		{
			basei=BASEPAIRS[trandpos].basepidx.first;
			basej=BASEPAIRS[trandpos].basepidx.second;

			if(basei > basej){
				tmp=basei;
				basei=basej;
				basej=tmp;
			}
			pairwiseidx=bpmapid2idx[basei*numseq+basej];

			inext=i+bpenergy[pairwiseidx].pcw[0]+bpenergy[pairwiseidx].pcw[1];
			if(i<=wpos && wpos<inext){
				break;
			}
			i=inext;
		}

		randval=Random();
		if(randval > 0.5) trandpos=basei;
		else trandpos=basej;
	}else{
		trandpos=int((numseq-2)*Random())+1;
	}

	span=2+int(4*Random());
	trandposst=trandpos-span/2;
	trandposed=trandposst+span;
	if(trandposst<1) trandposst=1;
	if(trandposed>numseq-1) trandposed=numseq-1;
	span=trandposed-trandposst;
	if(span < 2){
		newenergy=oldenergy;
		return false;
	}

	//始端または終端ノードがbase pairに含まれないようにする
	vector<pair<int,int>> pairsf; 
	for(i=trandposst;i<trandposed;i++){
		if(bppairednt[i] == -1){
			pairsf.push_back({i,-1});
		}else if(bppairednt[i] < 1){
			break;
		}else if(bppairednt[i] >= numseq-1){
			break;
		}else{
			pairsf.push_back({i,bppairednt[i]});
		}
	}

	vector<pair<int,int>> pairsb;
	for(i=trandposed-1;i>=trandposst;i--){
		if(bppairednt[i] == -1){
			pairsb.push_back({i,-1});
		}else if(bppairednt[i] < 1){
			break;
		}else if(bppairednt[i] >= numseq-1){
			break;
		}else{
			pairsb.push_back({i,bppairednt[i]});
		}
	}

	vector<pair<int,int>> pairs;
	if(pairsf.size() >= pairsb.size()){
		pairs=pairsf;
	}else{
		for (i=pairsb.size()-1; i>=0; i--) {
			pairs.push_back(pairsb[i]);
		}
	}

	copystr(decstr, numseq);

	vector<vector<int>> lmpparamidx;
	setlmpbpparam(restr, pairs, lmpparamidx);

	getgraphinfoLMP(pairs, numseq, lmpparamidx, movenodeinfo, edgeinfo, edgetype, paraminfo);

	map<int, vector<Residueatoms>> dnttemplates;
	map<int, TemplateDistMap> dnttemplateslength;
	for(i=0;i<edgeinfo.size();i++)
	{
		basei=edgeinfo[i].first; 
		basej=edgeinfo[i].second;
		
		if(edgetype[i] == 1){
			//if(basei < basej){
			selectdnttemplate(restr[basei], restr[basej], basei, basej, dnttemplates[i], dnttemplateslength[i]);
			//}else{
			//selectdnttemplate(restr[basej], restr[basei], basei, basej, dnttemplates[i], dnttemplateslength[i]);
			//}
		}
	}

	map<int, vector<Residueatoms>> bptemplates;
	map<int, TemplateDistMap> bptemplateslength;
	for(i=0;i<edgeinfo.size();i++){
		basei=edgeinfo[i].first; basej=edgeinfo[i].second;
		aat1 = restr[basei].aaa; aat2 = restr[basej].aaa;

		if(edgetype[i] == 2){
			selectbptemplate(restr[basei], restr[basej], bptemplates[i], bptemplateslength[i]);
		}
	}

	vector<int> movementresidues;
	for(i=0;i<pairs.size();i++)
	{
		basei = pairs[i].first; basej = pairs[i].second;
		if(basej == -1){
			movementresidues.push_back(basei);
		}else{
			movementresidues.push_back(basei);
			movementresidues.push_back(basej);
		}
	}

	Gradmain2 gradobj;
	vector<Residueatoms> restrprev, finalrestr;
	map<int, vector<Point3s>> gradptv, gradptvprev;
	deque<map<int, vector<Point3s>>> s_list, y_list;
	map<int, vector<Point3s>> z;
    deque<double> rho_list;
	double estop;
	int estopcount, estopnum;

	//initgradgeo(restr, movenodeinfo, numseq);
	randval = Random();
	if(randval < 0.5){
		initrandommove(restr, movenodeinfo, numseq);
	}

	initgrad2obj(gradobj, 0);
	calcene2(restr, decstr, gradobj, edgeinfo, edgetype, paraminfo, movementresidues, dnttemplates, bptemplates, dnttemplateslength, bptemplateslength, enecoeff);
	initgrad2obj(gradobj, 3);

	
	initgrad(gradptv, movenodeinfo);
	calcgrad2(restr, decstr, gradptv, gradobj, edgeinfo, edgetype, movenodeinfo, movementresidues, enecoeff);
	
	restrprev = restr;
	gradptvprev = gradptv;
	
	estopnum = 15;
	estopcount = 0;
	estop = 10000000000;
	finalrestr=restr;
	for (int epoch = 0; epoch < epochs; ++epoch){
		//cout<<"a-1"<<endl;
		if(epoch%100 == 0){
			for(i=0;i<edgeinfo.size();i++)
			{
				basei=edgeinfo[i].first; basej=edgeinfo[i].second;
				if(edgetype[i] == 1){
					//if(basei < basej){
					selectdnttemplate(restr[basei], restr[basej], basei, basej, dnttemplates[i], dnttemplateslength[i]);
					//}else{
						//selectdnttemplate(restr[basej], restr[basei], basei, basej, dnttemplates[i], dnttemplateslength[i]);
					//}
				}
			}

			for(i=0;i<edgeinfo.size();i++){
				basei=edgeinfo[i].first; basej=edgeinfo[i].second;
				aat1 = restr[basei].aaa; aat2 = restr[basej].aaa;

				if(edgetype[i] == 2){
					selectbptemplate(restr[basei], restr[basej], bptemplates[i], bptemplateslength[i]);
				}
			}
		}

		
		//cout << "epoch:" << epoch+1 << " mse: " << gradobj.mse << " mse(init): " << gradobj.inittempobj.mse << " mse(adjtmp): " << gradobj.adjtempobj.mse  << " mse(bptmp): " << gradobj.bptempobj.mse << " mse(length nt): " << gradobj.lengthobj.msent << " mse(adj temp length): " << gradobj.adjtemplengthobj.mse << " mse(adj length): " << gradobj.lengthobj.mseadj<<endl;
		
		//cout<<"a0"<<endl;
		lbfgsupdate(gradptv, s_list, y_list, rho_list, movenodeinfo, z);

		linesearch(restr, decstr, gradptv, z, gradobj, edgeinfo, edgetype, paraminfo, movenodeinfo, movementresidues, dnttemplates, bptemplates, dnttemplateslength, bptemplateslength, enecoeff, alpha, beta, gamma);
		// 新しい s, y, ρ を保存
		double rho_new =0.0;
		map<int, vector<Point3s>> s_new, y_new;
		for (auto& [key, grad_vec] : gradptv) {
			if(movenodeinfo[key]){
				s_new[key].resize(grad_vec.size());
				y_new[key].resize(grad_vec.size());
				for (size_t j = 0; j < grad_vec.size(); ++j) {
					//cout<<key << "-" << j << "-" << gradptv[key][j].x << "-" << gradptvprev[key][j].x<<endl;

					s_new[key][j].x = restr[key].ptv[j].x - restrprev[key].ptv[j].x;
					s_new[key][j].y = restr[key].ptv[j].y - restrprev[key].ptv[j].y;
					s_new[key][j].z = restr[key].ptv[j].z - restrprev[key].ptv[j].z;
					y_new[key][j].x = gradptv[key][j].x - gradptvprev[key][j].x;
					y_new[key][j].y = gradptv[key][j].y - gradptvprev[key][j].y;
					y_new[key][j].z = gradptv[key][j].z - gradptvprev[key][j].z;

					rho_new+=s_new[key][j].x*y_new[key][j].x;
					rho_new+=s_new[key][j].y*y_new[key][j].y;
					rho_new+=s_new[key][j].z*y_new[key][j].z;

					//cout<<key << "-" << j << "-" << s_new[key][j].x << "-" << y_new[key][j].x << "-" << s_new[key][j].y << "-" << y_new[key][j].y <<"-"<< s_new[key][j].z << "-" << y_new[key][j].z<<endl;
				}
			}
		}
		rho_new = 1.0/(rho_new+EPS);

		//cout<<rho_new<<endl;

		if(rho_new < 100000){
			if(m != 0 ){
				if (s_list.size() == m) {
					s_list.pop_front();
					y_list.pop_front();
					rho_list.pop_front();
				}
				s_list.push_back(s_new);
				y_list.push_back(y_new);
				rho_list.push_back(rho_new);
			}
		}

		if(abs(gradobj.mseprev - gradobj.mse) < tolerance){
			//cout << "Converged in " << epoch + 1 << " iterations." << endl;
            break;
		}

		if(estop > gradobj.mse){
			finalrestr=restr;
			estopcount=0;
			estop=gradobj.mse;
		}else{
			estopcount+=1;
		}

		if(estopnum == estopcount){
			break;
		}

		initgrad2obj(gradobj, 3); 
		restrprev = restr;
		gradptvprev = gradptv;

	}
	//cout << "Final:" << " mse: " << mse << " mse(length adj): " << mselengthadj << " mse(length nt): " << mselengthnt << " mse (length bp): " << mselengthbp << " mse(init): " << mseinit<<endl;
	restr=finalrestr;

	vector<bool> targetidxtmp(numseq, false);
	for(i=0;i<pairs.size();i++)
	{
		basei = pairs[i].first;
		basej = pairs[i].second;

		targetidxtmp[basei] = 1;
		if(basej != -1){
			targetidxtmp[basej] = 1;
		}
	}

	vector<pair<int,int>> continuousidx;
	int start = -1;

	for(i=0;i<numseq;i++){
		if (targetidxtmp[i] == 1) {
			if (start == -1) {
				start = i;
			}
		}else {
            if (start != -1) {
                continuousidx.push_back({start, i - 1});
                start = -1;
            }
        }
	}
	if (start != -1) {
        continuousidx.push_back({start, numseq - 1});
    }

	for(i = 0; i<continuousidx.size();i++){
		istart=continuousidx[i].first-1;
		iend=continuousidx[i].second+2;
		if(iend>=numseq) iend=numseq;
		//cout<<to_string(istart) + "-" + to_string(iend)<<endl;
		coord2parambb(restr,continuousidx[i].first,iend,1);
		coord2paramsugar(restr,continuousidx[i].first,iend,1);
		param2coordbase(restr,continuousidx[i].first,iend);
		param2coordhydrogen(restr,continuousidx[i].first,iend);

		coord2paramadd(restr,istart,iend);
		settmplmpflap(numseq, istart, iend);
	}

	newenergy=fcalcrmsdenergy(restr, numseq, lp2, rppmatbk, bpenergybk, bsenergybk, nmapbk, distmapbk, ntpairdistbk, adjpairdistbk, distenebk, multidistenebk, contactenebk, mapflagbk, mapbbflagbk, flagacc);

	/*¥
	cout<<"QqqqqqqqqqqqqqqqqQ"<<endl;
	checkupdate(restr, numseq);
	cout<<"PppppppppppppppppP"<<endl;
	fsetinitial(lp2,rppmatbk, bpenergybk, bsenergybk, nmapbk, distmapbk, ntpairdistbk, adjpairdistbk, distenebk, multidistenebk, contactenebk, mapflagbk, mapbbflagbk, numseq);
	double newenergy1=fcalcrmsdenergy(restr, numseq, lp2, rppmatbk, bpenergybk, bsenergybk, nmapbk, distmapbk, ntpairdistbk, adjpairdistbk, distenebk, multidistenebk, contactenebk, mapflagbk, mapbbflagbk, false);
	if(fabs(newenergy1-newenergy) > 10){
		cout<<"Error"<<endl;
		cout<<"2:: " + to_string(newenergy1)<<endl;
		cout<<"3:: " + to_string(newenergy)<<endl;
		exit(1);
	}
	cout<<"OK"<<endl;
	*/

	//exit(1);
	
	bool upidx=metropolisUpdate(oldenergy, newenergy, Bt);

	if(upidx)
	{
		rcopystr(decstr, restr, numseq);
		return true;
	}
	else
	{
		return false;
	}
}












