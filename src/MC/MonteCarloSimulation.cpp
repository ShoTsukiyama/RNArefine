
#include "MonteCarloSimulation.hpp"

MonteCarloSimulation::MonteCarloSimulation(const FragmentStructure &pdb_entry, const int &c, const ExtParams &extparameters, const FragmentStructure &inputtedstrall, const vector<vector<bool>> &inputtedstrallflag, const vector<FragmentStructure> &templatestrall, const vector<vector<vector<bool>>> &templatestrallflag, const vector<vector<vector<double>>> &initdistmap, const vector<vector<vector<bool>>> &initdistflag, const vector<vector<double>> &initdistscores, const vector<vector<vector<vector<double>>>> &templatedistmap, const vector<vector<vector<vector<bool>>>> &templatedistflag, const vector<vector<int>> &bindmap_input, const vector<vector<pair<int, double>>> &stmap_input, const vector<vector<double>> &cbpmap_input, const vector<vector<double>> &noncbpmap_input){
    
    setStr(pdb_entry, c);
    setOrigin();

    //LoadParams(pramdir);
    coord2parambase(refinestr,0,numbb);
    coord2paramsugar(refinestr,0,numbb);
    coord2parambb(refinestr,0,numbb);
    coord2paramadd(refinestr,0,numbb);
    coord2paramhydrogen(refinestr,0,numbb);
    
    param2coordbb(refinestr,0,numbb);
    param2coordsugar(refinestr,0,numbb);
    param2coordbase(refinestr,0,numbb);
    param2coordhydrogen(refinestr,0,numbb);

    coord2addparamhydrogen(refinestr,0,numbb);
   
    StrCopy(refinestr, refinestr_temp);
    Init(extparameters, inputtedstrall, inputtedstrallflag, templatestrall, templatestrallflag, initdistmap, initdistflag, initdistscores, templatedistmap, templatedistflag, bindmap_input, stmap_input, cbpmap_input, noncbpmap_input);
    //Simulation();
}

//MonteCarloSimulation::~MonteCarloSimulation()
//{
//}

void MonteCarloSimulation::StrCopy(vector<Residueatoms> strA, vector<Residueatoms> &strB)
{
    char base;
    int idx, j;
    map<char, int> atomnum = {{'A', 33}, {'C', 31}, {'G', 34}, {'U', 30}};
    for(int r=0;r<numbb;r++)
    {
        base = refinestr[r].aaa;
        for (int a=0;a<atomnum[base];a++)
        {
            strB[r].ptv[a].x=strA[r].ptv[a].x;
            strB[r].ptv[a].y=strA[r].ptv[a].y;
            strB[r].ptv[a].z=strA[r].ptv[a].z;
        }

        for(j=0;j<AccAtomParamNum[base][3];j++)
		{
            strB[r].dih[j]=strA[r].dih[j]; 
			strB[r].len[j]=strA[r].len[j]; 
			strB[r].ang[j]=strA[r].ang[j];
        }

        idx=0;
        for(j=AccAtomParamNum[base][3];j<AccAtomParamNum[base][4];j++)
        {
            strB[r].hyddih[idx]=strA[r].hyddih[idx]; 
			strB[r].hydlen[idx]=strA[r].hydlen[idx]; 
			strB[r].hydang[idx]=strA[r].hydang[idx];
            idx+=1;
        }

        for(j=0;j<AddAngleParamsList[base].size();j++) strB[r].addang[j]=strA[r].addang[j];
        for(j=0;j<AddLengthParamsList[base].size();j++) strB[r].addlen[j]=strA[r].addlen[j]; 
        for(j=0;j<AddHydLengthParamsList[base].size();j++){
            strB[r].addhydlen[j]=strA[r].addhydlen[j]; 
        }
        for(j=0;j<AddHydAngleParamsList[base].size();j++){
            strB[r].addhydang[j]=strA[r].addhydang[j]; 
        }
    }
}

void MonteCarloSimulation::setOrigin()
{
    originpoint.resize(3, vector<double>(3, 0.0));
    originpoint[0][0]=refinestr[0].ptv[0].x;
    originpoint[0][1]=refinestr[0].ptv[0].y;
    originpoint[0][2]=refinestr[0].ptv[0].z;
    originpoint[1][0]=refinestr[0].ptv[1].x;
    originpoint[1][1]=refinestr[0].ptv[1].y;
    originpoint[1][2]=refinestr[0].ptv[1].z;
    originpoint[2][0]=refinestr[0].ptv[2].x;
    originpoint[2][1]=refinestr[0].ptv[2].y;
    originpoint[2][2]=refinestr[0].ptv[2].z;
}

void MonteCarloSimulation::resetfrag(int index)
{
    int i,j;
    for(i=0;i<numbb;i++){
        ntflag[i]=index;

        for(j=0;j<numbb;j++){
            mapflag[i][j]=index;
            mapbbflag[i][j]=index;
        }

        char base=refinestr_temp[i].aaa;
        for(j=0;j<MCLENGTHPARAM[base].size();j++){
            lengthflag[i][j] = index;
        }

        for(j=0;j<MCINNERANGPARAM[base].size();j++){
            angleflag[i][j] = index;
        }

        for(j=0;j<MCTORSIONPARAM[base].size();j++){
            torsionflag[i][j] = index;
        }

        for(j=0;j<MCHYDLENGTHPARAM[base].size();j++){
            hydlengthflag[i][j] = index;
        }

        for(j=0;j<MCHYDINNERANGPARAM[base].size();j++){
            hydangleflag[i][j] = index;
        }

        for(j=0;j<MCHYDTORSIONPARAM[base].size();j++){
            hydtorsionflag[i][j] = index;
        }

        for(j=0;j<AddLengthParamsList[base].size();j++){
            addlengthflag[i][j] = index;
        }

        for(j=0;j<AddHydLengthParamsList[base].size();j++){
            hydaddlengthflag[i][j] = index;
        }

        for(j=0;j<AddHydAngleParamsList[base].size();j++){
            hydaddangleflag[i][j] = index;
        }

        for(j=0;j<AddAngleParamsList[base].size();j++){
            addangleflag[i][j] = index;
        }
    }
}

void MonteCarloSimulation::setStr(const FragmentStructure &pdb_entry, int c)
{
    int r, a, idx;
    map<char, int> nt2ntidx = {{'A', 0}, {'C', 1}, {'G', 2}, {'U', 3}};

    numbb = pdb_entry.residues.size();
    refinestr.resize(numbb);
    refinestr_temp.resize(numbb);
    for(r=0;r<numbb;r++)
    {
        refinestr[r].aaa = pdb_entry.residues[r].restype[0];
        refinestr[r].resind = nt2ntidx[refinestr[r].aaa];
        refinestr[r].residx = r;

        refinestr_temp[r].aaa = pdb_entry.residues[r].restype[0];
        refinestr_temp[r].resind = nt2ntidx[refinestr[r].aaa];
        refinestr_temp[r].residx = r;

        for (a=0;a<pdb_entry.residues[r].atoms.size();a++)
        {
            idx=atomid2mcidx[refinestr[r].aaa][a];
            refinestr[r].ptv[idx].x = pdb_entry.residues[r].atoms[a].coord.x;
            refinestr[r].ptv[idx].y = pdb_entry.residues[r].atoms[a].coord.y;
            refinestr[r].ptv[idx].z = pdb_entry.residues[r].atoms[a].coord.z;
        
            refinestr_temp[r].ptv[idx].x = pdb_entry.residues[r].atoms[a].coord.x;
            refinestr_temp[r].ptv[idx].y = pdb_entry.residues[r].atoms[a].coord.y;
            refinestr_temp[r].ptv[idx].z = pdb_entry.residues[r].atoms[a].coord.z;
        }
    }
}

void MonteCarloSimulation::Init(const ExtParams &extparameters, const FragmentStructure &inputtedstrall, const vector<vector<bool>> &inputtedstrallflag, const vector<FragmentStructure> &templatestrall, const vector<vector<vector<bool>>> &templatestrallflag, const vector<vector<vector<double>>> &initdistmap, const vector<vector<vector<bool>>> &initdistflag, const vector<vector<double>> &initdistscores, const vector<vector<vector<vector<double>>>> &templatedistmap, const vector<vector<vector<vector<bool>>>> &templatedistflag, const vector<vector<int>> &bindmap_input, const vector<vector<pair<int, double>>> &stmap_input, const vector<vector<double>> &cbpmap_input, const vector<vector<double>> &noncbpmap_input)
{   
    char base;
    temperature=-1;
    lsweepnum=5000;

    MCinputtedstrallPtr=&inputtedstrall;
    MCinputtedstrallflagPtr=&inputtedstrallflag;
    MCtemplatestrallPtr=&templatestrall;
    MCtemplatestrallflagPtr=&templatestrallflag;

    MCinitdistmap=&initdistmap;
    MCinitdistflag=&initdistflag;
    MCinitdistscore=&initdistscores;
    MCtemplatedistmap=&templatedistmap;
    MCtemplatedistflag=&templatedistflag;


    mapflag.resize(numbb, vector<int>(numbb, 1));
    mapbbflag.resize(numbb, vector<int>(numbb, 1));

    lengthenergyatoms.resize(numbb);
    angleenergyatoms.resize(numbb);
    torsionenergyatoms.resize(numbb);
    pairwiseenergyatoms.resize(numbb);

    for(int j=0;j<numbb;j++){
        base=refinestr_temp[j].aaa;
        pairwiseenergyatoms[j].resize(AccAtomParamNum[base][4],0);
    }

    stmap.resize(numbb, vector<double>(numbb));
    for(int i=0;i<numbb;i++){
        for(int j=0;j<numbb;j++){
            stmap[i][j] = stmap_input[i][j].second;
        }
    }

    bpmap.resize(numbb, vector<int>(numbb, 0));
    noncbpmap.resize(numbb, vector<double>(numbb, 0));
    cbpmap.resize(numbb, vector<double>(numbb, 0));
    for(int i=0;i<numbb;i++){
        for(int j=0;j<numbb;j++){
            if(bindmap_input[i][j] == 2) bpmap[i][j] = bindmap_input[i][j];
            if(cbpmap_input[i][j] > 0) cbpmap[i][j] = cbpmap_input[i][j];
            if(noncbpmap_input[i][j] > 0) noncbpmap[i][j] = noncbpmap_input[i][j];
        }
    }

    setpairwiseeneobj(extparameters, inputtedstrall, inputtedstrallflag, templatestrall, templatestrallflag);

    lengthflag.resize(numbb);
    angleflag.resize(numbb);
    torsionflag.resize(numbb);

    hydlengthflag.resize(numbb);
    hydangleflag.resize(numbb);
    hydtorsionflag.resize(numbb);

    addlengthflag.resize(numbb);
    addangleflag.resize(numbb);

    hydaddlengthflag.resize(numbb);
    hydaddangleflag.resize(numbb);

    ntflag.resize(numbb,1);

    for(int j=0;j<numbb;j++){
        base=refinestr_temp[j].aaa;
        lengthenergyatoms[j].resize(AccAtomParamNum[base][4],0);
        angleenergyatoms[j].resize(AccAtomParamNum[base][4],0);
        torsionenergyatoms[j].resize(AccAtomParamNum[base][4],0);
        
        lengthflag[j].resize(MCLENGTHPARAM[base].size()+MCBASELENGTHPARAM[base].size(),1);
        angleflag[j].resize(MCINNERANGPARAM[base].size()+MCBASEINNERANGPARAM[base].size(),1);
        torsionflag[j].resize(MCTORSIONPARAM[base].size()+MCBASETORSIONPARAM[base].size(),1);

        hydlengthflag[j].resize(MCHYDLENGTHPARAM[base].size(),1);
        hydangleflag[j].resize(MCHYDINNERANGPARAM[base].size(),1);
        hydtorsionflag[j].resize(MCHYDTORSIONPARAM[base].size(),1);

        hydaddlengthflag[j].resize(MCHYDADDLENGTHPARAM[base].size(),1);
        hydaddangleflag[j].resize(MCHYDADDANGLEPARAM[base].size(),1);

        addlengthflag[j].resize(AddLengthParamsList[base].size(),1);
        addangleflag[j].resize(AddAngleParamsList[base].size(),1);
    }

    setgeometricaleneobj();
}

void MonteCarloSimulation::setLogger(bool logflag, string logfile){
    logger.setPrint(false);
    if(logflag){
        logger.setVerbose(true);
        logger.setLogFile(logfile);
    }
}

double MonteCarloSimulation::Simulation(double maxtime, int startepochidx)
{
    double current_energy = calcenergymain();
    int moveidx;
    updateenergies();
    double newenergy;
    bool upflag;
    std::ostringstream log_buffer;
    //cout<<"Initial energy:: "<< current_energy<<endl;
    //logger.log("Initial energy: " + std::to_string(current_energy));

    resetfrag();
    
    double timelimit = (maxtime > 0) ? maxtime : std::numeric_limits<double>::infinity();

    //cout<<timelimit<<endl;
    //cout<<maxtime<<endl;

    //cout<<"Time: " << timelimit<<endl;
    //logger.log("Time limit: " + std::to_string(timelimit));
    auto start_time = std::chrono::high_resolution_clock::now();
    
    for(int k=0;k<lsweepnum;k++)
	{
        upflag=calcmovementmain(current_energy, newenergy, moveidx);

        log_buffer << startepochidx + k + 1 << ", "
                   << moveidx << ", "
                   << int(upflag) << ", "
                   << temperature << ", "
                   << current_energy << ", "
                   << newenergy << "\n";

        cout<<"Epoch "<<startepochidx+k+1<<", Current energy: "<<current_energy << ", New energy:" << newenergy<<endl;
        if (k % 100 == 0 && !log_buffer.str().empty()) {
            cout<<"Epoch "<<startepochidx+k+1<<", Current energy: "<<current_energy << ", New energy:" << newenergy<<endl;
            logger.log(log_buffer.str());
            log_buffer.str("");
            log_buffer.clear();
        }

        if(upflag){
            current_energy = newenergy;
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_time = end_time - start_time;
    
        if(elapsed_time.count() > timelimit) break;
    }

    if (log_buffer.tellp() > 0) { 
        logger.log(log_buffer.str());
    }

    return current_energy;
}

void MonteCarloSimulation::returnStr(FragmentStructure &pdb_entry, int c)
{
    int r, a;
    map<char, int> nt2ntidx = {{'A', 0}, {'C', 1}, {'G', 2}, {'U', 3}};

    for(r=0;r<pdb_entry.residues.size();r++)
    {
        for (a=0;a<pdb_entry.residues[r].atoms.size();a++)
        {
            pdb_entry.residues[r].atoms[a].coord.x = refinestr[r].ptv[atomid2mcidx[refinestr[r].aaa][a]].x;
            pdb_entry.residues[r].atoms[a].coord.y = refinestr[r].ptv[atomid2mcidx[refinestr[r].aaa][a]].y;
            pdb_entry.residues[r].atoms[a].coord.z = refinestr[r].ptv[atomid2mcidx[refinestr[r].aaa][a]].z;
        }
    }
}

void MonteCarloSimulation::printStrParams(){
    int r,a;

    for(r=0;r<refinestr.size();r++)
    {
        cout<<"[[" << refinestr[r].aaa << "]]"<<endl;
        for(a=0;a<ATOM_RES_NUM;a++)
        {
            cout<<r << "-" << a << "-" << refinestr[r].ptv[a].x << "-" << refinestr[r].ptv[a].y << "-" << refinestr[r].ptv[a].z<<endl;
        }
    }
}

void MonteCarloSimulation::setInitStr(const FragmentStructure &pdb_entry, int c){
    int r, a;
    initstr.resize(pdb_entry.residues.size());
    for(r=0;r<pdb_entry.residues.size();r++)
    {
        initstr[r].aaa = pdb_entry.residues[r].restype[0];
        for (a=0;a<pdb_entry.residues[r].atoms.size();a++)
        {
            initstr[r].ptv[a].x = pdb_entry.residues[r].atoms[a].coord.x;
            initstr[r].ptv[a].y = pdb_entry.residues[r].atoms[a].coord.y;
            initstr[r].ptv[a].z = pdb_entry.residues[r].atoms[a].coord.z;
        }
    }
}

void Swap(double &Ti, double &Tj){
    double temp;
    temp = Ti;
    Ti = Tj;
    Tj = temp;
}

bool SwapCheck(double Ei, double Ej, double Ti, double Tj)
{
    double v = exp((Ej-Ei)/((1/Tj)-(1/Ti)));
    return (v > 1.0) || (Random() < v);
}

void REMC_swapping(const vector<double> &energies, vector<double> &temperatures)
{
    int i;
    double Ei, Ej, Ti, Tj;
    double randval = Random();
    bool swapflag;

    //cout<<temperatures.size()<<endl;

    if(randval > 0.5){
        for(i=0;i<temperatures.size()-1;i+=2){
            //cout<<"========================="<<endl;
            //cout<<i << "-" << i+1 << "-" << temperatures[i] << "-" << temperatures[i+1]<<endl;
            Ei=energies[i];
            Ej=energies[i+1];
            Ti=temperatures[i];
            Tj=temperatures[i+1];
            swapflag = SwapCheck(Ei, Ej, Ti, Tj);
            if (swapflag) Swap(temperatures[i], temperatures[i+1]);
            //cout<<i << "-" << i+1 << "-" << temperatures[i] << "-" << temperatures[i+1]<<endl;
        }
    }else{
        for(i=1;i<temperatures.size()-1;i+=2){
            //cout<<"========================="<<endl;
            //cout<<i << "-" << i+1 << "-" << temperatures[i] << "-" << temperatures[i+1]<<endl;
            Ei=energies[i];
            Ej=energies[i+1];
            Ti=temperatures[i];
            Tj=temperatures[i+1];
            swapflag = SwapCheck(Ei, Ej, Ti, Tj);
            if (swapflag) Swap(temperatures[i], temperatures[i+1]);
            //cout<<i << "-" << i+1 << "-" << temperatures[i] << "-" << temperatures[i+1]<<endl;
        }
    }
}

void copy_pdb_entry(const ModelUnit &pdb_entry, ModelUnit &pdb_entry_temp)
{
    pdb_entry_temp.chains.resize(pdb_entry.chains.size());
    for(int i=0;i<pdb_entry.chains.size();i++){
        pdb_entry_temp.chains[i].chainID=pdb_entry.chains[i].chainID;
        pdb_entry_temp.chains[i].chainID_full=pdb_entry.chains[i].chainID_full;
        pdb_entry_temp.chains[i].sequence=pdb_entry.chains[i].sequence;
        pdb_entry_temp.chains[i].sarst=pdb_entry.chains[i].sarst;
        pdb_entry_temp.chains[i].residues.resize(pdb_entry.chains[i].residues.size());
        for(int j=0;j<pdb_entry.chains[i].residues.size();j++){
            pdb_entry_temp.chains[i].residues[j].het=pdb_entry.chains[i].residues[j].het;
            pdb_entry_temp.chains[i].residues[j].resi=pdb_entry.chains[i].residues[j].resi;
            pdb_entry_temp.chains[i].residues[j].icode=pdb_entry.chains[i].residues[j].icode;
            pdb_entry_temp.chains[i].residues[j].resn=pdb_entry.chains[i].residues[j].resn;
            pdb_entry_temp.chains[i].residues[j].atoms.resize(heavily_atom_num[string(1,pdb_entry_temp.chains[i].residues[j].resn[2])]);
            for(int k=0;k<heavily_atom_num[string(1,pdb_entry_temp.chains[i].residues[j].resn[2])];k++){
                pdb_entry_temp.chains[i].residues[j].atoms[k].name=pdb_entry.chains[i].residues[j].atoms[k].name;
                pdb_entry_temp.chains[i].residues[j].atoms[k].xyz.resize(3,0);
                pdb_entry_temp.chains[i].residues[j].atoms[k].xyz[0] = pdb_entry.chains[i].residues[j].atoms[k].xyz[0];
                pdb_entry_temp.chains[i].residues[j].atoms[k].xyz[1] = pdb_entry.chains[i].residues[j].atoms[k].xyz[1];
                pdb_entry_temp.chains[i].residues[j].atoms[k].xyz[2] = pdb_entry.chains[i].residues[j].atoms[k].xyz[2];
                pdb_entry_temp.chains[i].residues[j].atoms[k].movable=pdb_entry.chains[i].residues[j].atoms[k].movable;
                pdb_entry_temp.chains[i].residues[j].atoms[k].input=pdb_entry.chains[i].residues[j].atoms[k].input;
            }
        }
    }
}

void copy_fragment_str(const FragmentStructure &str, FragmentStructure &str_temp)
{
    str_temp.residues.resize(str.residues.size());
    for(int i=0;i<str.residues.size();i++){
        str_temp.residues[i].restype=str.residues[i].restype;
        str_temp.residues[i].atoms.resize(heavily_atom_num[str_temp.residues[i].restype]);
        for(int j=0;j<heavily_atom_num[str_temp.residues[i].restype];j++){
            str_temp.residues[i].atoms[j].coord.x=str.residues[i].atoms[j].coord.x;
            str_temp.residues[i].atoms[j].coord.y=str.residues[i].atoms[j].coord.y;
            str_temp.residues[i].atoms[j].coord.z=str.residues[i].atoms[j].coord.z;
        }
    }
}

void outReplica(const vector<std::unique_ptr<MonteCarloSimulation>> &mcsims, int c, FragmentStructure &fragstr, ModelUnit &pdb_str, const string &outdir){
    int i, j;

    for(int i=0;i<mcsims.size();i++){
        mcsims[i]->returnStr(fragstr, c);
        for(int j=0;j<pdb_str.chains[c].residues.size();j++){
            UpdateStructure(c, j, pdb_str, fragstr.residues[j]);
        }
        string remark="";
        write_pdb_structure((outdir+"/MC/replica_"+to_string(i+1)+".pdb").c_str(),pdb_str,remark);
    }
}

bool checkREMCtermination(const std::vector<double>& energies){
    
    const int num_low_temp = 10;
    if (energies.size() < num_low_temp) return false;
    std::vector<double> low_temp_energies(energies.begin(), energies.begin() + num_low_temp);

    double avg = std::accumulate(low_temp_energies.begin(), low_temp_energies.end(), 0.0) / num_low_temp;

    double variance = 0.0;
    for (double e : low_temp_energies) {
        variance += (e - avg) * (e - avg);
    }
    variance /= num_low_temp;

    return variance < (10e-4 * avg);
}

void REMC(const ModelUnit &pdb_entry, FragmentStructure &superimposedstr, int c, const string &pramdir, const ExtParams &extparameters, const FragmentStructure &inputtedstrall, const vector<vector<bool>> &inputtedstrallflag, const vector<FragmentStructure> &templatestrall, const vector<vector<vector<bool>>> &templatestrallflag, const vector<vector<vector<double>>> &initdistmap, const vector<vector<vector<bool>>> &initdistflag, const vector<vector<double>> &initdistscores, const vector<vector<vector<vector<double>>>> &templatedistmap, const vector<vector<vector<vector<bool>>>> &templatedistflag, const vector<vector<int>> &bindmap, const vector<vector<pair<int, double>>> &stmap, const vector<vector<double>> &cbpmap, const vector<vector<double>> &noncbpmap, double maxtime, const string &outdir, const bool &logflag){
    int globaliter=300;
    int swapspan=200;
    int repnum=40;
    double mintemperature=0;
    double maxtemperature=80;
    double tempspan = (maxtemperature - mintemperature) / (repnum - 1);
    bool term;

    fs::create_directories(outdir+"/MC");

    LoadParams(pramdir);

    FragmentStructure mediatestr;
    copy_fragment_str(superimposedstr, mediatestr);
    //RemoveHydrogens(superimposedstr);

    ModelUnit pdb_entry_temp;
    copy_pdb_entry(pdb_entry, pdb_entry_temp);

    //vector<FragmentStructure> superimposedstrs(repnum, superimposedstr);
    vector<std::unique_ptr<MonteCarloSimulation>> mcsims;

    vector<double> temperatures(repnum);
    for (int i = 0; i < repnum; ++i){
        temperatures[i] = mintemperature + i * tempspan;
        //cout<<superimposedstrs[i].residues.size()<<endl;
        //MonteCarloSimulation mcsim(superimposedstrs[i], c, pramdir, extparameters, inputtedstrall, inputtedstrallflag, bindmap, stmap);
        //cout<<"a2"<<endl;
        //mcsim.setsimulationparams(temperatures[i], swapspan);
        mcsims.push_back(std::make_unique<MonteCarloSimulation>(superimposedstr, c, extparameters, inputtedstrall, inputtedstrallflag, templatestrall, templatestrallflag, initdistmap, initdistflag, initdistscores, templatedistmap, templatedistflag, bindmap, stmap, cbpmap, noncbpmap));
        mcsims.back()->setsimulationparams(temperatures[i], swapspan);
    }

    if(logflag){
        fs::create_directories(outdir+"/logs/REMC");
        for (int i = 0; i < repnum; ++i){
            mcsims[i]->setLogger(logflag, outdir+"/logs/REMC/replica_"+to_string(i+1)+".log");
        }
    }

    auto start_time = std::chrono::high_resolution_clock::now();

    vector<double> energies(repnum);
    for(int i=0;i<globaliter;i++){
        
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(int j=0;j<mcsims.size();j++){
            //cout<<"[[" << j << "]]"<<endl;
            energies[j] = mcsims[j]->Simulation(-1, i*swapspan);
        }

        outReplica(mcsims, c, mediatestr, pdb_entry_temp, outdir);

        term=checkREMCtermination(energies);
        if(term) break;

        REMC_swapping(energies, temperatures);

        for(int j=0;j<mcsims.size();j++){
            mcsims[j]->setsimulationparams(temperatures[j], swapspan);
            //cout<<j << "-" << mcsims[j]->gettemperature()<<endl;
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_time = end_time - start_time;
        if(elapsed_time.count() > maxtime) break;
    }

    auto min_it = min_element(energies.begin(), energies.end());
    int min_str_index = std::distance(energies.begin(), min_it);

    mcsims[min_str_index]->returnStr(superimposedstr, c);
}

double exeEM(FragmentStructure &superimposedstr, int c, string pramdir, const ExtParams &extparameters, const FragmentStructure &inputtedstrall, const vector<vector<bool>> &inputtedstrallflag, const vector<FragmentStructure> &templatestrall, const vector<vector<vector<bool>>> &templatestrallflag, const vector<vector<vector<double>>> &initdistmap, const vector<vector<vector<bool>>> &initdistflag, const vector<vector<double>> &initdistscores, const vector<vector<vector<vector<double>>>> &templatedistmap, const vector<vector<vector<vector<bool>>>> &templatedistflag, const vector<vector<int>> &bindmap, const vector<vector<pair<int, double>>> &stmap, const vector<vector<double>> &cbpmap, const vector<vector<double>> &noncbpmap, double maxtime, const string &outdir, const bool &logflag)
{
    LoadParams(pramdir);
    MonteCarloSimulation mcsim(superimposedstr, c, extparameters, inputtedstrall, inputtedstrallflag, templatestrall, templatestrallflag, initdistmap, initdistflag, initdistscores, templatedistmap, templatedistflag, bindmap, stmap, cbpmap, noncbpmap);
    mcsim.setsimulationparams(3.0, 50000);
    if(logflag) mcsim.setLogger(logflag, outdir+"/logs/mc_simulation.log");
    double ene = mcsim.Simulation(maxtime, 0);
    mcsim.returnStr(superimposedstr, c);
    return ene;
}

double MonteCarloSimulationProcess(ModelUnit &pdb_entry, const size_t c, const string &pramdir,
    const vector<pair<double,vector<size_t>>> &bp_vec, const map<int, vector<pair<int, int>>> &loopnts, const vector<pair<double,vector<size_t>>> &noncbp_vec,
    const vector<vector<int>> &stemconnectnts, const map<int, vector<vector<vector<int>>>> &connectnts, const vector<vector<pair<int, double>>> &stmap,
    const FragmentStructure &inputtedstrall, const vector<vector<bool>> &inputtedstrallflag, const vector<double> &init_weight,
    const vector<FragmentStructure> &templatestrall, const vector<vector<vector<bool>>> &templatestrallflag,
    const ExtParams &extparameters, vector<vector<vector<double>>> &hydrogencoords, double maxtime, const string &outdir, const bool &logflag, const int option)
{
    size_t r,i,j;
    int r1,r2,c1,c2;
    size_t a;
    
    double startene, endene;
    int resnum=pdb_entry.chains[c].residues.size();
    FragmentStructure superimposedstr;
    FragmentStructure inputtedstr;
    vector<vector<bool>> inputtedstrflag;
    vector<vector<int>> bindmap(resnum, vector<int>(resnum, 0));
    vector<vector<double>> cbpmap(resnum, vector<double>(resnum, 0));
    vector<vector<double>> noncbpmap(resnum, vector<double>(resnum, 0));
    vector<bool> moveresidues(resnum, false);
    double score;
    
    superimposedstr.residues.resize(resnum);
    inputtedstr.residues.resize(resnum);
    inputtedstrflag.resize(resnum);
    for(int j=0;j<pdb_entry.chains[c].residues.size()-1;j++){
        ExtractNt(c, j, pdb_entry, superimposedstr.residues[j]);
        inputtedstr.residues[j] = inputtedstrall.residues[j];
        inputtedstrflag[j] = inputtedstrallflag[j];
        bindmap[j][j+1] = 1;
        moveresidues[j] = true;
    }
    j=pdb_entry.chains[c].residues.size()-1;
    ExtractNt(c, j, pdb_entry, superimposedstr.residues[j]);
    inputtedstr.residues[j] = inputtedstrall.residues[j];
    inputtedstrflag[j] = inputtedstrallflag[j];
    moveresidues[j] = true;

    for(i=0;i<bp_vec.size();i++){
        c1 = bp_vec[i].second[0];
        r1 = bp_vec[i].second[1];
        c2 = bp_vec[i].second[2];
        r2 = bp_vec[i].second[3];

        if(c1 == c && c2 == c){
            bindmap[r1][r2]=2;
            bindmap[r2][r1]=2;
        }
    }

    for(i=0;i<bp_vec.size();i++){
        score = bp_vec[i].first;
        c1 = bp_vec[i].second[0];
        r1 = bp_vec[i].second[1];
        c2 = bp_vec[i].second[2];
        r2 = bp_vec[i].second[3];

        if(c1 == c && c2 == c){
            cbpmap[r1][r2]=score;
            cbpmap[r2][r1]=score;
        }
    }

    for(i=0;i<noncbp_vec.size();i++){
        score = noncbp_vec[i].first;
        c1 = noncbp_vec[i].second[0];
        r1 = noncbp_vec[i].second[1];
        c2 = noncbp_vec[i].second[2];
        r2 = noncbp_vec[i].second[3];

        if(c1 == c && c2 == c){
            noncbpmap[r1][r2]=score;
            noncbpmap[r2][r1]=score;
        }
    }
    
    AttachHydrogen(superimposedstr, extparameters.hydparam, extparameters.hydparamidx);

    vector<vector<vector<double>>> initdistmap;
    vector<vector<vector<bool>>> initdistflag;
    vector<vector<double>> initdistscores;
    vector<vector<vector<vector<double>>>> templatedistmap;
    vector<vector<vector<vector<bool>>>> templatedistflag;
    string base1, base2;
    initdistmap.resize(resnum);
    initdistflag.resize(resnum);
    initdistscores.resize(resnum);
    for(int i=0;i<resnum;i++){
        base1=superimposedstr.residues[i].restype;
        initdistmap[i].resize(resnum);
        initdistflag[i].resize(resnum);
        initdistscores[i].resize(resnum);
        for(int j=i+1;j<resnum;j++){
            base2=superimposedstr.residues[j].restype;
            initdistmap[i][j].resize(AccAtomParamNum[base1[0]][2], 0);
            initdistflag[i][j].resize(AccAtomParamNum[base1[0]][2], false);
            initdistscores[i][j] = min(init_weight[i], init_weight[j]);
            for(int k=0;k<template_target_atoms_idx.size();k++){
                if(k < inputtedstrallflag[i].size() && k < inputtedstrallflag[j].size()){
                    int l =  template_target_atoms_idx[k];
					if(inputtedstrallflag[i][l] && inputtedstrallflag[j][l]){
                        Point3s initatomA = inputtedstrall.residues[i].atoms[l].coord;
                        Point3s initatomB = inputtedstrall.residues[j].atoms[l].coord;

                        Point3s initvec;
                        initvec.x = initatomA.x - initatomB.x;
                        initvec.y = initatomA.y - initatomB.y;
                        initvec.z = initatomA.z - initatomB.z;

                        double initdist = sqrt(initvec.x*initvec.x+initvec.y*initvec.y+initvec.z*initvec.z);
                        initdistmap[i][j][l] = initdist;
                        initdistflag[i][j][l] = true;
                    }
                }
            }
        }
        initdistmap[i][i].resize(template_target_atoms_idx.size(), 0);
    }

    templatedistmap.resize(templatestrall.size());
    templatedistflag.resize(templatestrall.size());
    for(int n=0;n<templatestrall.size();n++){
        templatedistmap[n].resize(resnum);
        templatedistflag[n].resize(resnum);
        for(int i=0;i<resnum;i++){
            base1=superimposedstr.residues[i].restype;
            templatedistmap[n][i].resize(resnum);
            templatedistflag[n][i].resize(resnum);
            for(int j=i+1;j<resnum;j++){
                base2=superimposedstr.residues[j].restype;
                templatedistmap[n][i][j].resize(AccAtomParamNum[base1[0]][2], 0);
                templatedistflag[n][i][j].resize(AccAtomParamNum[base1[0]][2], false);
                for(int k=0;k<template_target_atoms_idx.size();k++){
                    if(k < templatestrallflag[n][i].size() && k < templatestrallflag[n][j].size()){
                        int l =  template_target_atoms_idx[k];
                        if(templatestrallflag[n][i][l] && templatestrallflag[n][j][l]){
                            Point3s initatomA = templatestrall[n].residues[i].atoms[l].coord;
                            Point3s initatomB = templatestrall[n].residues[j].atoms[l].coord;

                            Point3s initvec;
                            initvec.x = initatomA.x - initatomB.x;
                            initvec.y = initatomA.y - initatomB.y;
                            initvec.z = initatomA.z - initatomB.z;

                            double initdist = sqrt(initvec.x*initvec.x+initvec.y*initvec.y+initvec.z*initvec.z);
                            templatedistmap[n][i][j][l] = initdist;
                            templatedistflag[n][i][j][l] = true;
                        }
                    }
                }
            }
        }
    }

    //MonteCarloSimulation mcsim(superimposedstr, c, pramdir, extparameters, inputtedstrall, inputtedstrallflag, bindmap, noncbpmap);
    //mcsim.returnStr(superimposedstr, c);
    //mcsim.printStrParams();
    double final_ene=exeEM(superimposedstr, c, pramdir, extparameters, inputtedstrall, inputtedstrallflag, templatestrall, templatestrallflag, initdistmap, initdistflag, initdistscores, templatedistmap, templatedistflag, bindmap, stmap, cbpmap, noncbpmap, maxtime, outdir, logflag);
    //REMC(pdb_entry, superimposedstr, c, pramdir, extparameters, inputtedstrall, inputtedstrallflag, bindmap, stmap, maxtime, outdir, logflag);

    hydrogencoords.resize(pdb_entry.chains[c].residues.size());
    for(int j=0;j<pdb_entry.chains[c].residues.size();j++)
    {
        GetHydrogenInfo(c, j, pdb_entry, superimposedstr.residues[j], hydrogencoords[j]);
    }
    
    RemoveHydrogens(superimposedstr);

    for(int j=0;j<pdb_entry.chains[c].residues.size();j++)
    {
        UpdateStructure(c, j, pdb_entry, superimposedstr.residues[j]);
        //UpdateStructure(c, j, pdb_entry, superimposedstr.residues[j]);
    }
    string remark="";
    write_pdb_structure((outdir+"/mc.pdb").c_str(),pdb_entry,remark);

    //return true;
    return final_ene;
}

double MonteCarloSimulationProcess(ModelUnit &pdb_entry, const string &pramdir, 
    const vector<pair<double,vector<size_t>>> &bp_vec, const map<int, vector<pair<int, int>>> &loopnts, const vector<pair<double,vector<size_t>>> &noncbp_vec,
    const vector<vector<int>> &stemconnectnts, const map<int, vector<vector<vector<int>>>> &connectnts, const vector<vector<pair<int, double>>> &stmap,
    const FragmentStructure &inputtedstrall, const vector<vector<bool>> &inputtedstrallflag, const vector<double> &init_weight,
    const vector<FragmentStructure> &templatestrall, const vector<vector<vector<bool>>> &templatestrallflag,
    const ExtParams &extparameters, vector<vector<vector<double>>> &hydrogencoords, double maxtime, const string &outdir, const bool &logflag, const int option)
{
    double final_ene;
    size_t c;
    for (c=0;c<pdb_entry.chains.size();c++){
        final_ene = MonteCarloSimulationProcess(pdb_entry, c, pramdir, bp_vec, loopnts, noncbp_vec, stemconnectnts, connectnts, stmap, inputtedstrall, inputtedstrallflag, init_weight, templatestrall, templatestrallflag, extparameters, hydrogencoords, maxtime, outdir, logflag, option);
    }
    
    return final_ene;
}




