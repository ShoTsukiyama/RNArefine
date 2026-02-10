
#include "GetInteraction.hpp"

const double Pm_mu    = 53.755; const double Pm_sd    =42.773; // P[i-1]-P[i]-P[j]-P[j-1]
const double Pp_mu    = 80.233; const double Pp_sd    =20.466; // P[i+1]-P[i]-P[j]-P[j+1]
const double O5m_mu   = 60.047; const double O5m_sd   =38.731; // O5'[i-1]-O5'[i]-O5'[j]-O5'[j-1]
const double O5p_mu   = 96.645; const double O5p_sd   =20.185; // O5'[i+1]-O5'[i]-O5'[j]-O5'[j+1]
const double C5m_mu   = 58.526; const double C5m_sd   =36.222; // C5'[i-1]-C5'[i]-C5'[j]-C5'[j-1]
const double C5p_mu   =100.129; const double C5p_sd   =23.505; // C5'[i+1]-C5'[i]-C5'[j]-C5'[j+1]
const double C4m_mu   = 64.830; const double C4m_sd   =31.640; // C4'[i-1]-C4'[i]-C4'[j]-C4'[j-1]
const double C4p_mu   =115.546; const double C4p_sd   =23.989; // C4'[i+1]-C4'[i]-C4'[j]-C4'[j+1]
const double C3m_mu   = 80.193; const double C3m_sd   =26.156; // C3'[i-1]-C3'[i]-C3'[j]-C3'[j-1]
const double C3p_mu   =133.291; const double C3p_sd   =26.435; // C3'[i+1]-C3'[i]-C3'[j]-C3'[j+1]
const double C2m_mu   = 87.048; const double C2m_sd   =31.557; // C2'[i-1]-C2'[i]-C2'[j]-C2'[j-1]
const double C2p_mu   =147.550; const double C2p_sd   =29.374; // C2'[i+1]-C2'[i]-C2'[j]-C2'[j+1]
const double C1m_mu   = 74.283; const double C1m_sd   =38.789; // C1'[i-1]-C1'[i]-C1'[j]-C1'[j-1]
const double C1p_mu   =134.728; const double C1p_sd   =27.557; // C1'[i+1]-C1'[i]-C1'[j]-C1'[j+1]
const double O4m_mu   = 62.955; const double O4m_sd   =38.254; // O4'[i-1]-O4'[i]-O4'[j]-O4'[j-1]
const double O4p_mu   =113.863; const double O4p_sd   =24.303; // O4'[i+1]-O4'[i]-O4'[j]-O4'[j+1]
const double O3m_mu   = 80.432; const double O3m_sd   =22.356; // O3'[i-1]-O3'[i]-O3'[j]-O3'[j-1]
const double O3p_mu   =137.697; const double O3p_sd   =28.655; // O3'[i+1]-O3'[i]-O3'[j]-O3'[j+1]
const double P44P_mu  = -0.166; const double P44P_sd  =28.901; // P[i]-C4'[i]-C4'[j]-P[j]
const double C4114C_mu=172.313; const double C4114C_sd=43.508; // C4'[i]-C1'[i]-C1'[j]-C4'[j]
const double PP_mu    = 18.564; const double PP_sd    = 0.896; // P[i]-P[j]
const double O5O5_mu  = 16.808; const double O5O5_sd  = 0.774; // O5'[i]-O5'[j]
const double C5C5_mu  = 17.160; const double C5C5_sd  = 0.513; // C5'[i]-C5'[j]
const double C4C4_mu  = 15.018; const double C4C4_sd  = 0.419; // C4'[i]-C4'[j]
const double C3C3_mu  = 13.676; const double C3C3_sd  = 0.498; // C3'[i]-C3'[j]
const double C2C2_mu  = 11.046; const double C2C2_sd  = 0.515; // C2'[i]-C2'[j]
const double C1C1_mu  = 10.659; const double C1C1_sd  = 0.336; // C1'[i]-C1'[j]
const double O4O4_mu  = 13.232; const double O4O4_sd  = 0.399; // O4'[i]-O4'[j]
const double O3O3_mu  = 15.266; const double O3O3_sd  = 0.641; // O3'[i]-O3'[j]
const double NN_mu    =  8.929; const double NN_sd    = 0.250; // N[i]-N[j]
const double aPm_mu   =109.192; const double aPm_sd   =21.448; // <P[i-1]P[i],P[j+1]P[j]>
const double aPp_mu   =109.326; const double aPp_sd   =20.496; // <P[i+1]P[i],P[j-1]P[j]>
const double aO5m_mu  = 97.545; const double aO5m_sd  =18.629; // <O5'[i-1]O5'[i],O5'[j+1]O5'[j]>
const double aO5p_mu  = 97.816; const double aO5p_sd  =18.227; // <O5'[i+1]O5'[i],O5'[j-1]O5'[j]>
const double aC5m_mu  = 94.309; const double aC5m_sd  =19.724; // <C5'[i-1]C5'[i],C5'[j+1]C5'[j]>
const double aC5p_mu  = 94.640; const double aC5p_sd  =19.843; // <C5'[i+1]C5'[i],C5'[j-1]C5'[j]>
const double aC4m_mu  = 81.590; const double aC4m_sd  =17.368; // <C4'[i-1]C4'[i],C4'[j+1]C4'[j]>
const double aC4p_mu  = 81.352; const double aC4p_sd  =18.082; // <C4'[i+1]C4'[i],C4'[j-1]C4'[j]>
const double aC3m_mu  = 71.338; const double aC3m_sd  =17.754; // <C3'[i-1]C3'[i],C3'[j+1]C3'[j]>
const double aC3p_mu  = 70.590; const double aC3p_sd  =18.100; // <C3'[i+1]C3'[i],C3'[j-1]C3'[j]>
const double aC2m_mu  = 59.411; const double aC2m_sd  =20.249; // <C2'[i-1]C2'[i],C2'[j+1]C2'[j]>
const double aC2p_mu  = 58.352; const double aC2p_sd  =20.141; // <C2'[i+1]C2'[i],C2'[j-1]C2'[j]>
const double aC1m_mu  = 64.542; const double aC1m_sd  =19.817; // <C1'[i-1]C1'[i],C1'[j+1]C1'[j]>
const double aC1p_mu  = 63.616; const double aC1p_sd  =20.319; // <C1'[i+1]C1'[i],C1'[j-1]C1'[j]>
const double aO4m_mu  = 78.198; const double aO4m_sd  =17.830; // <O4'[i-1]O4'[i],O4'[j+1]O4'[j]>
const double aO4p_mu  = 77.672; const double aO4p_sd  =18.859; // <O4'[i+1]O4'[i],O4'[j-1]O4'[j]>
const double aO3m_mu  = 70.536; const double aO3m_sd  =17.393; // <O3'[i-1]O3'[i],O3'[j+1]O3'[j]>
const double aO3p_mu  = 69.919; const double aO3p_sd  =17.786; // <O3'[i+1]O3'[i],O3'[j-1]O3'[j]>
const double aPC_mu   = 59.764; const double aPC_sd   =18.397; // <P[i]C4'[i],P[j]C4'[j]>
const double aCC_mu   =165.213; const double aCC_sd   =13.607; // <C4'[i]C1'[i],C4'[j]C1'[j]>


/* if a residue is paired with multiple other residues, just keep the
 * highest scoring residue pair */
void filter_bp(vector<pair<double,vector<size_t> > >&bp_vec)
{
    vector<pair<size_t,size_t> > paired_res_vec;
    vector<pair<double,size_t> > bp_idx_vec;
    vector<size_t> filtered_bp_vec;
    size_t bp;
    for (bp=0;bp<bp_vec.size();bp++)
        bp_idx_vec.push_back(pair<double,size_t>(-bp_vec[bp].first,bp));
    sort(bp_idx_vec.begin(), bp_idx_vec.end());
    
    for (bp=0;bp<bp_idx_vec.size();bp++)
    {
        //if (-bp_idx_vec[bp].first<0.5) break;
        if (find(paired_res_vec.begin(), paired_res_vec.end(), 
            pair<size_t,size_t>(bp_vec[bp_idx_vec[bp].second].second[0],
            bp_vec[bp_idx_vec[bp].second].second[1]))!=paired_res_vec.end())
            continue;
        if (find(paired_res_vec.begin(), paired_res_vec.end(),
            pair<size_t,size_t>(bp_vec[bp_idx_vec[bp].second].second[2],
            bp_vec[bp_idx_vec[bp].second].second[3]))!=paired_res_vec.end())
            continue;
        
        filtered_bp_vec.push_back(bp_idx_vec[bp].second);
        paired_res_vec.push_back(pair<size_t,size_t>(
                    bp_vec[bp_idx_vec[bp].second].second[0],
                    bp_vec[bp_idx_vec[bp].second].second[1]));
        paired_res_vec.push_back(pair<size_t,size_t>(
                    bp_vec[bp_idx_vec[bp].second].second[2],
                    bp_vec[bp_idx_vec[bp].second].second[3]));
    }
    vector<pair<size_t,size_t> >().swap(paired_res_vec);
    vector<pair<double,size_t> >().swap(bp_idx_vec);

    vector<pair<double,vector<size_t> > >bp_vec_tmp;
    for (bp=0;bp<bp_vec.size();bp++)
    {
        if (find(filtered_bp_vec.begin(), filtered_bp_vec.end(),bp)!=filtered_bp_vec.end()) bp_vec_tmp.push_back(bp_vec[bp]);
        bp_vec[bp].second.clear();
    }
    bp_vec.clear();

    for (bp=0;bp<bp_vec_tmp.size();bp++) bp_vec.push_back(bp_vec_tmp[bp]);

    vector<pair<double,vector<size_t> > >().swap(bp_vec_tmp);
    vector<size_t>().swap(filtered_bp_vec);
    return;
}

bool bp_tor_score(
    const vector<double> &c1, const vector<double> &c2,
    const vector<double> &c3, const vector<double> &c4,
    const double mu, const double sd, const double tol, const double weight,
    double &nominator, double &denominator)
{
    double tor=rad2deg(Points2Dihedral(c1,c2,c3,c4));
    if(tor<-180) return false;
    denominator+=weight;
    double diff=fabs(tor-mu); 
    if (diff>180) diff-=180;
    nominator+=weight*(1-diff/(tol* sd));
    return true;
}

bool bp_ang_score(
    const vector<double> &c1, const vector<double> &c2,
    const vector<double> &c3, const vector<double> &c4,
    const double mu, const double sd, const double tol, const double weight,
    double &nominator, double &denominator)
{
    double ang=rad2deg(Points4Angle(c1,c2,c3,c4));
    if(ang<-180) return false;
    denominator+=weight;
    nominator+=weight*(1-fabs(ang-mu)/(tol* sd));
    return true;
}

bool bp_len_score(const vector<double>&c1, const vector<double>&c2,
    const double mu, const double sd, const double tol, const double weight,
    double &nominator, double &denominator)
{
    denominator+=weight;
    nominator+=weight*(1-fabs(Points2Distance(c1,c2)-mu)/(tol*sd));
    return true;
}

bool bp_nn_score(const bool previnextj, const bool nextiprevj,
    const bool has_prev_ci,       const bool has_next_ci,
    const bool has_prev_cj,       const bool has_next_cj,
    const vector<double> &prev_ci, const vector<double> &next_ci,
    const vector<double> &prev_cj, const vector<double> &next_cj,
    const double mu, const double sd, const double tol, const double weight,
    double &nominator, double &denominator)
{
    if (!(has_prev_ci && has_next_cj) && !(has_next_ci && has_prev_cj))
        return false;
    denominator+=weight;
    double nominator_previnextj=0;
    double nominator_nextiprevj=0;
    if (has_prev_ci  && has_next_cj )
    {
        nominator_previnextj=1-fabs(Points2Distance(prev_ci,next_cj)-mu)/(tol*sd);
        if (!previnextj && nominator_previnextj>0) nominator_previnextj=0;
    }
    if (has_next_ci  && has_prev_cj )
    {
        nominator_nextiprevj=1-fabs(Points2Distance(next_ci,prev_cj)-mu)/(tol*sd);
        if (!nextiprevj && nominator_nextiprevj>0) nominator_nextiprevj=0;
    }

    if (nominator_previnextj>nominator_nextiprevj)
         nominator+=weight*nominator_previnextj;
    else nominator+=weight*nominator_nextiprevj;
    return true;
}

double mean_vec(const vector<double>&vec)
{
    double sum=0;
    int i;
    for (i=0;i<vec.size();i++) sum+=vec[i];
    return sum/vec.size();
}

void cssr(const ModelUnit &pdb_entry, vector<vector<size_t> >&res_str_vec,
    vector<pair<double,vector<size_t> > >&bp_vec, 
    const bool interchain)
{
    /* pre-trained parameters */
    double weight_len=1;
    double weight_nn =1;
    double weight_tor=1;
    double weight_ang=1;
    double tol=3; // tolerance, in the unit of standard deviation
    double adjust1=0.1;  // adjust for varying number of tests
    double adjust2=0;    // adjust for varying std
    double adjust3=-0.07;// baseline value
    double totaltest=(10-1)*(weight_len>0)+
                    (10-1)*(weight_nn >0)+
                    ( 9-1)*(weight_tor>0)+
                    (10-1)*(weight_ang>0);
    
    /* other variables */
    vector<size_t> tmp_bp(5,0);
    vector<size_t> tmp_nt(2,0);
    double ang;
    size_t c1, c2, r1, r2, a1, a2;
    char base1,base2;
    char base1prev,base1next;
    char base2prev,base2next;
    char icode,chainID;
    stringstream ss;
    double nominator, denominator;
    double nominator_previnextj;
    double nominator_nextiprevj;
    int   successtest;
    bool previnextj,nextiprevj;

    /* coordinates of previous residue */
    vector<double>prev_Pi(3,0);  bool has_prev_Pi;  vector<double>prev_Pj(3,0);  bool has_prev_Pj; 
    vector<double>prev_O5i(3,0); bool has_prev_O5i; vector<double>prev_O5j(3,0); bool has_prev_O5j;
    vector<double>prev_C5i(3,0); bool has_prev_C5i; vector<double>prev_C5j(3,0); bool has_prev_C5j;
    vector<double>prev_C4i(3,0); bool has_prev_C4i; vector<double>prev_C4j(3,0); bool has_prev_C4j;
    vector<double>prev_C3i(3,0); bool has_prev_C3i; vector<double>prev_C3j(3,0); bool has_prev_C3j;
    vector<double>prev_C2i(3,0); bool has_prev_C2i; vector<double>prev_C2j(3,0); bool has_prev_C2j;
    vector<double>prev_C1i(3,0); bool has_prev_C1i; vector<double>prev_C1j(3,0); bool has_prev_C1j;
    vector<double>prev_O4i(3,0); bool has_prev_O4i; vector<double>prev_O4j(3,0); bool has_prev_O4j;
    vector<double>prev_O3i(3,0); bool has_prev_O3i; vector<double>prev_O3j(3,0); bool has_prev_O3j;
    vector<double>prev_Nxi(3,0); bool has_prev_Nxi; vector<double>prev_Nxj(3,0); bool has_prev_Nxj;
    /* coordinates of current residue */
    vector<double>Pi(3,0);  bool has_Pi;  vector<double>Pj(3,0);  bool has_Pj; 
    vector<double>O5i(3,0); bool has_O5i; vector<double>O5j(3,0); bool has_O5j;
    vector<double>C5i(3,0); bool has_C5i; vector<double>C5j(3,0); bool has_C5j;
    vector<double>C4i(3,0); bool has_C4i; vector<double>C4j(3,0); bool has_C4j;
    vector<double>C3i(3,0); bool has_C3i; vector<double>C3j(3,0); bool has_C3j;
    vector<double>C2i(3,0); bool has_C2i; vector<double>C2j(3,0); bool has_C2j;
    vector<double>C1i(3,0); bool has_C1i; vector<double>C1j(3,0); bool has_C1j;
    vector<double>O4i(3,0); bool has_O4i; vector<double>O4j(3,0); bool has_O4j;
    vector<double>O3i(3,0); bool has_O3i; vector<double>O3j(3,0); bool has_O3j;
    vector<double>Nxi(3,0); bool has_Nxi; vector<double>Nxj(3,0); bool has_Nxj;
    /* coordinates of next residue */
    vector<double>next_Pi(3,0);  bool has_next_Pi;  vector<double>next_Pj(3,0);  bool has_next_Pj; 
    vector<double>next_O5i(3,0); bool has_next_O5i; vector<double>next_O5j(3,0); bool has_next_O5j;
    vector<double>next_C5i(3,0); bool has_next_C5i; vector<double>next_C5j(3,0); bool has_next_C5j;
    vector<double>next_C4i(3,0); bool has_next_C4i; vector<double>next_C4j(3,0); bool has_next_C4j;
    vector<double>next_C3i(3,0); bool has_next_C3i; vector<double>next_C3j(3,0); bool has_next_C3j;
    vector<double>next_C2i(3,0); bool has_next_C2i; vector<double>next_C2j(3,0); bool has_next_C2j;
    vector<double>next_C1i(3,0); bool has_next_C1i; vector<double>next_C1j(3,0); bool has_next_C1j;
    vector<double>next_O4i(3,0); bool has_next_O4i; vector<double>next_O4j(3,0); bool has_next_O4j;
    vector<double>next_O3i(3,0); bool has_next_O3i; vector<double>next_O3j(3,0); bool has_next_O3j;
    vector<double>next_Nxi(3,0); bool has_next_Nxi; vector<double>next_Nxj(3,0); bool has_next_Nxj;

    /* loop over base pairs */
    for (c1=0; c1<pdb_entry.chains.size(); c1++)
    {
        tmp_nt[0]=tmp_bp[0]=c1;
        for (r1=0; r1<pdb_entry.chains[c1].residues.size(); r1++)
        {
            base1=pdb_entry.chains[c1].residues[r1].resn[2];
            base1prev=base1next=0;
            if (!pdb_entry.chains[c1].residues[r1].atoms.size()) continue;
            //if ((chainID=pdb_entry.chains[c1].chainID)!=' ') 
                //ss<<chainID<<".";
            //ss<<base1<<pdb_entry.chains[c1].residues[r1].resi;
            //if ((icode=pdb_entry.chains[c1].residues[r1].icode)!=' ')
                //ss<<'^'<<icode;
            //tmp_bp[0]=ss.str();
            tmp_nt[1]=tmp_bp[1]=r1;
            res_str_vec.push_back(tmp_nt);
            ss.str(string());

            has_prev_Pi  = false;
            has_prev_O5i = false;
            has_prev_C5i = false;
            has_prev_C4i = false;
            has_prev_C3i = false;
            has_prev_C2i = false;
            has_prev_C1i = false;
            has_prev_O4i = false;
            has_prev_O3i = false;
            has_prev_Nxi = false;
            has_Pi       = false;
            has_O5i      = false;
            has_C5i      = false;
            has_C4i      = false;
            has_C3i      = false;
            has_C2i      = false;
            has_C1i      = false;
            has_O4i      = false;
            has_O3i      = false;
            has_Nxi      = false;
            has_next_Pi  = false;
            has_next_O5i = false;
            has_next_C5i = false;
            has_next_C4i = false;
            has_next_C3i = false;
            has_next_C2i = false;
            has_next_C1i = false;
            has_next_O4i = false;
            has_next_O3i = false;
            has_next_Nxi = false;

            for (a1=0;a1<pdb_entry.chains[c1].residues[r1].atoms.size();a1++)
            {
                if (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" P  ")
                {
                    has_Pi=true;
                    Pi=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
                }
                else if (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" O5'")
                {
                    has_O5i=true;
                    O5i=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
                }
                else if (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" C5'")
                {
                    has_C5i=true;
                    C5i=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
                }
                else if (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" C4'")
                {
                    has_C4i=true;
                    C4i=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
                }
                else if (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" C3'")
                {
                    has_C3i=true;
                    C3i=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
                }
                else if (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" C2'")
                {
                    has_C2i=true;
                    C2i=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
                }
                else if (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" C1'")
                {
                    has_C1i=true;
                    C1i=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
                }
                else if (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" O4'")
                {
                    has_O4i=true;
                    O4i=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
                }
                else if (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" O3'")
                {
                    has_O3i=true;
                    O3i=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
                }
                else if ((pdb_entry.chains[c1].residues[r1].atoms[a1].name==" N9 "
                        && (base1=='A' || base1=='G')) ||
                        (pdb_entry.chains[c1].residues[r1].atoms[a1].name==" N1 "
                        && (base1=='C' || base1=='T' || base1=='U')))
                {
                    has_Nxi=true;
                    Nxi=pdb_entry.chains[c1].residues[r1].atoms[a1].xyz;
                }
            }

            if (r1>0 && (weight_tor>0 || weight_ang>0 || weight_nn>0))
            {
                base1prev=pdb_entry.chains[c1].residues[r1-1].resn[2];
                for (a1=0;a1<pdb_entry.chains[c1].residues[r1-1].atoms.size();a1++)
                {
                    if (pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" P  ")
                    {
                        has_prev_Pi=true;
                        prev_Pi=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz;
                    }
                    else if (pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" O5'")
                    {
                        has_prev_O5i=true;
                        prev_O5i=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz;
                    }
                    else if (pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" C5'")
                    {
                        has_prev_C5i=true;
                        prev_C5i=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz;
                    }
                    else if (pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" C4'")
                    {
                        has_prev_C4i=true;
                        prev_C4i=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz;
                    }
                    else if (pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" C3'")
                    {
                        has_prev_C3i=true;
                        prev_C3i=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz;
                    }
                    else if (pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" C2'")
                    {
                        has_prev_C2i=true;
                        prev_C2i=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz;
                    }
                    else if (pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" C1'")
                    {
                        has_prev_C1i=true;
                        prev_C1i=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz;
                    }
                    else if (pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" O4'")
                    {
                        has_prev_O4i=true;
                        prev_O4i=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz;
                    }
                    else if (pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" O3'")
                    {
                        has_prev_O3i=true;
                        prev_O3i=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz;
                    }
                    else if ((pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" N9 "
                            && (base1prev=='A' || base1prev=='G')) ||
                             (pdb_entry.chains[c1].residues[r1-1].atoms[a1].name==" N1 "
                            && (base1prev=='C' || base1prev=='T' || base1prev=='U')))
                    {
                        has_prev_Nxi=true;
                        prev_Nxi=pdb_entry.chains[c1].residues[r1-1].atoms[a1].xyz;
                    }
                }
            }

            if (r1<pdb_entry.chains[c1].residues.size()-1 && (weight_tor>0 || weight_ang>0 || weight_nn>0))
            {
                base1next=pdb_entry.chains[c1].residues[r1+1].resn[2];
                for (a1=0;a1<pdb_entry.chains[c1].residues[r1+1].atoms.size();a1++)
                {
                    if (pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" P  ")
                    {
                        has_next_Pi=true;
                        next_Pi=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz;
                    }
                    else if (pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" O5'")
                    {
                        has_next_O5i=true;
                        next_O5i=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz;
                    }
                    else if (pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" C5'")
                    {
                        has_next_C5i=true;
                        next_C5i=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz;
                    }
                    else if (pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" C4'")
                    {
                        has_next_C4i=true;
                        next_C4i=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz;
                    }
                    else if (pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" C3'")
                    {
                        has_next_C3i=true;
                        next_C3i=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz;
                    }
                    else if (pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" C2'")
                    {
                        has_next_C2i=true;
                        next_C2i=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz;
                    }
                    else if (pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" C1'")
                    {
                        has_next_C1i=true;
                        next_C1i=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz;
                    }
                    else if (pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" O4'")
                    {
                        has_next_O4i=true;
                        next_O4i=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz;
                    }
                    else if (pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" O3'")
                    {
                        has_next_O3i=true;
                        next_O3i=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz;
                    }
                    else if ((pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" N9 "
                            && (base1next=='A' || base1next=='G')) ||
                             (pdb_entry.chains[c1].residues[r1+1].atoms[a1].name==" N1 "
                            && (base1next=='C' || base1next=='T' || base1next=='U')))
                    {
                        has_next_Nxi=true;
                        next_Nxi=pdb_entry.chains[c1].residues[r1+1].atoms[a1].xyz;
                    }
                }
            }


            for (c2=c1+interchain; c2<pdb_entry.chains.size(); c2++)
            {
                tmp_bp[2]=c2;
                for (r2=(c1==c2)*(r1+1); r2<pdb_entry.chains[c2].residues.size(); r2++)
                {
                    if (!pdb_entry.chains[c2].residues[r2].atoms.size() ||
                        Points2Distance2(
                        pdb_entry.chains[c1].residues[r1].atoms[0].xyz,
                        pdb_entry.chains[c2].residues[r2].atoms[0].xyz)>530)
                        continue;
                    base2=pdb_entry.chains[c2].residues[r2].resn[2];
                    base2prev=base2next=0;
                    if ((base1=='A' &&(base2=='U' || base2=='T')) ||
                        (base1=='C' && base2=='G')||
                       ((base1=='U' || base1=='T')&& base2=='A')  ||
                        (base1=='G' && base2=='C'))
                        tmp_bp[4]=0;
                        //tmp_bp[2]=base1+string("-")+base2+" WC";
                    else if ((base1=='G' && base2=='U') ||
                             (base1=='U' && base2=='G')) 
                        tmp_bp[4]=1;
                        //tmp_bp[2]=base1+string("-")+base2+" Wobble";
                    else continue;

                    //if ((chainID=pdb_entry.chains[c2].chainID)!=' ') 
                        //ss<<chainID<<".";
                    //ss<<base2<<pdb_entry.chains[c2].residues[r2].resi;
                    //if ((icode=pdb_entry.chains[c2].residues[r2].icode)!=' ')
                        //ss<<'^'<<icode;
                    //tmp_bp[1]=ss.str();
                    tmp_bp[3]=r2;
                    //ss.str(string());

                    has_prev_Pj  = false;
                    has_prev_O5j = false;
                    has_prev_C5j = false;
                    has_prev_C4j = false;
                    has_prev_C3j = false;
                    has_prev_C2j = false;
                    has_prev_C1j = false;
                    has_prev_O4j = false;
                    has_prev_O3j = false;
                    has_Pj       = false;
                    has_O5j      = false;
                    has_C5j      = false;
                    has_C4j      = false;
                    has_C3j      = false;
                    has_C2j      = false;
                    has_C1j      = false;
                    has_O4j      = false;
                    has_O3j      = false;
                    has_Nxj      = false;
                    has_next_Pj  = false;
                    has_next_O5j = false;
                    has_next_C5j = false;
                    has_next_C4j = false;
                    has_next_C3j = false;
                    has_next_C2j = false;
                    has_next_C1j = false;
                    has_next_O4j = false;
                    has_next_O3j = false;

                    for (a2=0;a2<pdb_entry.chains[c2].residues[r2].atoms.size();a2++)
                    {
                        if (pdb_entry.chains[c2].residues[r2].atoms[a2].name==" P  ")
                        {
                            has_Pj=true;
                            Pj=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz;
                        }
                        else if (pdb_entry.chains[c2].residues[r2].atoms[a2].name==" O5'")
                        {
                            has_O5j=true;
                            O5j=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz;
                        }
                        else if (pdb_entry.chains[c2].residues[r2].atoms[a2].name==" C5'")
                        {
                            has_C5j=true;
                            C5j=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz;
                        }
                        else if (pdb_entry.chains[c2].residues[r2].atoms[a2].name==" C4'")
                        {
                            has_C4j=true;
                            C4j=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz;
                        }
                        else if (pdb_entry.chains[c2].residues[r2].atoms[a2].name==" C3'")
                        {
                            has_C3j=true;
                            C3j=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz;
                        }
                        else if (pdb_entry.chains[c2].residues[r2].atoms[a2].name==" C2'")
                        {
                            has_C2j=true;
                            C2j=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz;
                        }
                        else if (pdb_entry.chains[c2].residues[r2].atoms[a2].name==" C1'")
                        {
                            has_C1j=true;
                            C1j=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz;
                        }
                        else if (pdb_entry.chains[c2].residues[r2].atoms[a2].name==" O4'")
                        {
                            has_O4j=true;
                            O4j=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz;
                        }
                        else if (pdb_entry.chains[c2].residues[r2].atoms[a2].name==" O3'")
                        {
                            has_O3j=true;
                            O3j=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz;
                        }
                        else if ((pdb_entry.chains[c2].residues[r2].atoms[a2].name==" N9 "
                                && (base2=='A' || base2=='G')) ||
                                (pdb_entry.chains[c2].residues[r2].atoms[a2].name==" N1 "
                                && (base2=='C' || base2=='T' || base2=='U')))
                        {
                            has_Nxj=true;
                            Nxj=pdb_entry.chains[c2].residues[r2].atoms[a2].xyz;
                        }
                    }

                    if (r2>0 && (weight_tor>0 || weight_ang>0 || weight_nn>0) && c2-c1+r2-r1>1)
                    {
                        base2prev=pdb_entry.chains[c2].residues[r2-1].resn[2];
                        for (a2=0;a2<pdb_entry.chains[c2].residues[r2-1].atoms.size();a2++)
                        {
                            if (pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" P  ")
                            {
                                has_prev_Pj=true;
                                prev_Pj=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz;
                            }
                            else if (pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" O5'")
                            {
                                has_prev_O5j=true;
                                prev_O5j=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz;
                            }
                            else if (pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" C5'")
                            {
                                has_prev_C5j=true;
                                prev_C5j=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz;
                            }
                            else if (pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" C4'")
                            {
                                has_prev_C4j=true;
                                prev_C4j=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz;
                            }
                            else if (pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" C3'")
                            {
                                has_prev_C3j=true;
                                prev_C3j=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz;
                            }
                            else if (pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" C2'")
                            {
                                has_prev_C2j=true;
                                prev_C2j=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz;
                            }
                            else if (pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" C1'")
                            {
                                has_prev_C1j=true;
                                prev_C1j=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz;
                            }
                            else if (pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" O4'")
                            {
                                has_prev_O4j=true;
                                prev_O4j=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz;
                            }
                            else if (pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" O3'")
                            {
                                has_prev_O3j=true;
                                prev_O3j=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz;
                            }
                            else if ((pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" N9 "
                                  && (base2prev=='A' || base2prev=='G')) ||
                                     (pdb_entry.chains[c2].residues[r2-1].atoms[a2].name==" N1 "
                                  && (base2prev=='C' || base2prev=='T' || base2prev=='U')))
                            {
                                has_prev_Nxj=true;
                                prev_Nxj=pdb_entry.chains[c2].residues[r2-1].atoms[a2].xyz;
                            }
                        }
                    }

                    if (r2<pdb_entry.chains[c2].residues.size()-1 && 
                       (weight_tor>0 || weight_ang>0 || weight_nn>0) && c2-c1+r2-r1>1)
                    {
                        base2next=pdb_entry.chains[c2].residues[r2+1].resn[2];
                        for (a2=0;a2<pdb_entry.chains[c2].residues[r2+1].atoms.size();a2++)
                        {
                            if (pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" P  ")
                            {
                                has_next_Pj=true;
                                next_Pj=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz;
                            }
                            else if (pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" O5'")
                            {
                                has_next_O5j=true;
                                next_O5j=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz;
                            }
                            else if (pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" C5'")
                            {
                                has_next_C5j=true;
                                next_C5j=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz;
                            }
                            else if (pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" C4'")
                            {
                                has_next_C4j=true;
                                next_C4j=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz;
                            }
                            else if (pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" C3'")
                            {
                                has_next_C3j=true;
                                next_C3j=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz;
                            }
                            else if (pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" C2'")
                            {
                                has_next_C2j=true;
                                next_C2j=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz;
                            }
                            else if (pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" C1'")
                            {
                                has_next_C1j=true;
                                next_C1j=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz;
                            }
                            else if (pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" O4'")
                            {
                                has_next_O4j=true;
                                next_O4j=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz;
                            }
                            else if (pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" O3'")
                            {
                                has_next_O3j=true;
                                next_O3j=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz;
                            }
                            else if ((pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" N9 "
                                  && (base2next=='A' || base2next=='G')) ||
                                     (pdb_entry.chains[c2].residues[r2+1].atoms[a2].name==" N1 "
                                  && (base2next=='C' || base2next=='T' || base2next=='U')))
                            {
                                has_next_Nxj=true;
                                next_Nxj=pdb_entry.chains[c2].residues[r2+1].atoms[a2].xyz;
                            }
                        }
                    }

                    nominator=denominator=successtest=0;
                    if (weight_len>0)
                    {
                        successtest--;
                        if (has_Pi  && has_Pj ) successtest+=bp_len_score( Pi, Pj,  PP_mu,  PP_sd,tol,weight_len,nominator,denominator); //   PP: P[i]-P[j]
                        if (has_O5i && has_O5j) successtest+=bp_len_score(O5i,O5j,O5O5_mu,O5O5_sd,tol,weight_len,nominator,denominator); // O5O5: O5'[i]-O5'[j]
                        if (has_C5i && has_C5j) successtest+=bp_len_score(C5i,C5j,C5C5_mu,C5C5_sd,tol,weight_len,nominator,denominator); // C5C5: C5'[i]-C5'[j]
                        if (has_C4i && has_C4j) successtest+=bp_len_score(C4i,C4j,C4C4_mu,C4C4_sd,tol,weight_len,nominator,denominator); // C4C4: C4'[i]-C4'[j]
                        if (has_C3i && has_C3j) successtest+=bp_len_score(C3i,C3j,C3C3_mu,C3C3_sd,tol,weight_len,nominator,denominator); // C3C3: C3'[i]-C3'[j]
                        if (has_C2i && has_C2j) successtest+=bp_len_score(C2i,C2j,C2C2_mu,C2C2_sd,tol,weight_len,nominator,denominator); // C2C2: C2'[i]-C2'[j]
                        if (has_C1i && has_C1j) successtest+=bp_len_score(C1i,C1j,C1C1_mu,C1C1_sd,tol,weight_len,nominator,denominator); // C1C1: C1'[i]-C1'[j]
                        if (has_O4i && has_O4j) successtest+=bp_len_score(O4i,O4j,O4O4_mu,O4O4_sd,tol,weight_len,nominator,denominator); // O4O4: O4'[i]-O4'[j]
                        if (has_O3i && has_O3j) successtest+=bp_len_score(O3i,O3j,O3O3_mu,O3O3_sd,tol,weight_len,nominator,denominator); // O3O3: O3'[i]-O3'[j]
                        if (has_Nxi && has_Nxj) successtest+=bp_len_score(Nxi,Nxj,  NN_mu,  NN_sd,tol,weight_len,nominator,denominator); //   NN: N[i]-N[j]
                    }
                    if (weight_nn>0)
                    {
                        successtest--;
                        previnextj=((base1prev=='A' &&(base2next=='U' || base2next=='T')) ||
                                    (base1prev=='C' && base2next=='G')||
                                   ((base1prev=='U' || base1prev=='T')&& base2next=='A')  ||
                                    (base1prev=='G' && base2next=='C'));
                        nextiprevj=((base1next=='A' &&(base2prev=='U' || base2prev=='T')) ||
                                    (base1next=='C' && base2prev=='G')||
                                   ((base1next=='U' || base1next=='T')&& base2prev=='A')  ||
                                    (base1next=='G' && base2prev=='C'));
                        if (has_Pi  && has_Pj) //   PP: P[i]-P[j]
                            successtest+=bp_nn_score(previnextj,nextiprevj,has_prev_Pi, has_next_Pi, has_prev_Pj, has_next_Pj,
                                prev_Pi, next_Pi, prev_Pj, next_Pj,   PP_mu,  PP_sd,tol,weight_nn,nominator,denominator);
                        if (has_O5i && has_O5j) // O5O5: O5'[i]-O5'[j]
                            successtest+=bp_nn_score(previnextj,nextiprevj,has_prev_O5i,has_next_O5i,has_prev_O5j,has_next_O5j,
                                prev_O5i,next_O5i,prev_O5j,next_O5j,O5O5_mu,O5O5_sd,tol,weight_nn,nominator,denominator);
                        if (has_C5i && has_C5j) // C5C5: C5'[i]-C5'[j]
                            successtest+=bp_nn_score(previnextj,nextiprevj,has_prev_C5i,has_next_C5i,has_prev_C5j,has_next_C5j,
                                prev_C5i,next_C5i,prev_C5j,next_C5j,C5C5_mu,C5C5_sd,tol,weight_nn,nominator,denominator);
                        if (has_C4i && has_C4j) // C4C4: C4'[i]-C4'[j]
                            successtest+=bp_nn_score(previnextj,nextiprevj,has_prev_C4i,has_next_C4i,has_prev_C4j,has_next_C4j,
                                prev_C4i,next_C4i,prev_C4j,next_C4j,C4C4_mu,C4C4_sd,tol,weight_nn,nominator,denominator);
                        if (has_C3i && has_C3j) // C3C3: C3'[i]-C3'[j]
                            successtest+=bp_nn_score(previnextj,nextiprevj,has_prev_C3i,has_next_C3i,has_prev_C3j,has_next_C3j,
                                prev_C3i,next_C3i,prev_C3j,next_C3j,C3C3_mu,C3C3_sd,tol,weight_nn,nominator,denominator);
                        if (has_C2i && has_C2j) // C2C2: C2'[i]-C2'[j]
                            successtest+=bp_nn_score(previnextj,nextiprevj,has_prev_C2i,has_next_C2i,has_prev_C2j,has_next_C2j,
                                prev_C2i,next_C2i,prev_C2j,next_C2j,C2C2_mu,C2C2_sd,tol,weight_nn,nominator,denominator);
                        if (has_C1i && has_C1j) // C1C1: C1'[i]-C1'[j]
                            successtest+=bp_nn_score(previnextj,nextiprevj,has_prev_C1i,has_next_C1i,has_prev_C1j,has_next_C1j,
                                prev_C1i,next_C1i,prev_C1j,next_C1j,C1C1_mu,C1C1_sd,tol,weight_nn,nominator,denominator);
                        if (has_O4i && has_O4j) // O4O4: O4'[i]-O4'[j]
                            successtest+=bp_nn_score(previnextj,nextiprevj,has_prev_O4i,has_next_O4i,has_prev_O4j,has_next_O4j,
                                prev_O4i,next_O4i,prev_O4j,next_O4j,O4O4_mu,O4O4_sd,tol,weight_nn,nominator,denominator);
                        if (has_O3i && has_O3j) // O3O3: O3'[i]-O3'[j]
                            successtest+=bp_nn_score(previnextj,nextiprevj,has_prev_O3i,has_next_O3i,has_prev_O3j,has_next_O3j,
                                prev_O3i,next_O3i,prev_O3j,next_O3j,O3O3_mu,O3O3_sd,tol,weight_nn,nominator,denominator);
                        if (has_Nxi && has_Nxj) //   NN: N[i]-N[j]
                            successtest+=bp_nn_score(previnextj,nextiprevj,has_prev_Nxi,has_next_Nxi,has_prev_Nxj,has_next_Nxj,
                                prev_Nxi,next_Nxi,prev_Nxj,next_Nxj,  NN_mu,  NN_sd,tol,weight_nn,nominator,denominator);
                    }
                    if (weight_tor>0)
                    {
                        successtest--;
                        if (has_next_Pi  && has_Pi  && has_Pj  && has_next_Pj ) successtest+=bp_tor_score( next_Pi, Pi, Pj, next_Pj, Pp_mu, Pp_sd, tol,weight_tor,nominator,denominator); //  Pp:  P[i+1]-P[i]-P[j]-P[j-1]
                        if (has_next_O5i && has_O5i && has_O5j && has_next_O5j) successtest+=bp_tor_score(next_O5i,O5i,O5j,next_O5j,O5p_mu,O5p_sd, tol,weight_tor,nominator,denominator); // O5p: O5'[i+1]-O5'[i]-O5'[j]-O5'[j-1]
                        if (has_next_C5i && has_C5i && has_C5j && has_next_C5j) successtest+=bp_tor_score(next_C5i,C5i,C5j,next_C5j,C5p_mu,C5p_sd, tol,weight_tor,nominator,denominator); // C5p: C5'[i+1]-C5'[i]-C5'[j]-C5'[j-1]
                        if (has_next_C4i && has_C4i && has_C4j && has_next_C4j) successtest+=bp_tor_score(next_C4i,C4i,C4j,next_C4j,C4p_mu,C4p_sd, tol,weight_tor,nominator,denominator); // C4p: C4'[i+1]-C4'[i]-C4'[j]-C4'[j-1]
                        if (has_next_C3i && has_C3i && has_C3j && has_next_C3j) successtest+=bp_tor_score(next_C3i,C3i,C3j,next_C3j,C3p_mu,C3p_sd, tol,weight_tor,nominator,denominator); // C3p: C3'[i+1]-C3'[i]-C3'[j]-C3'[j-1]
                        if (has_next_C2i && has_C2i && has_C2j && has_next_C2j) successtest+=bp_tor_score(next_C2i,C2i,C2j,next_C2j,C2p_mu,C2p_sd, tol,weight_tor,nominator,denominator); // C2p: C2'[i+1]-C2'[i]-C2'[j]-C2'[j-1]
                        if (has_next_C1i && has_C1i && has_C1j && has_next_C1j) successtest+=bp_tor_score(next_C1i,C1i,C1j,next_C1j,C1p_mu,C1p_sd, tol,weight_tor,nominator,denominator); // C1p: C1'[i+1]-C1'[i]-C1'[j]-C1'[j-1]
                        if (has_next_O4i && has_O4i && has_O4j && has_next_O4j) successtest+=bp_tor_score(next_O4i,O4i,O4j,next_O4j,O4p_mu,O4p_sd, tol,weight_tor,nominator,denominator); // O4p: O4'[i+1]-O4'[i]-O4'[j]-O4'[j-1]
                        if (has_prev_O3i && has_O3i && has_O3j && has_prev_O3j) successtest+=bp_tor_score(prev_O3i,O3i,O3j,prev_O3j,O3m_mu,O3m_sd, tol,weight_tor,nominator,denominator); // O3m: O3'[i-1]-O3'[i]-O3'[j]-O3'[j+1]
                    }
                    if (weight_ang>0)
                    {
                        successtest--;
                        if (has_next_Pi  && has_Pi  && has_prev_Pj  && has_Pj ) successtest+=bp_ang_score( next_Pi, Pi, prev_Pj, Pj, aPp_mu, aPp_sd, tol,weight_ang,nominator,denominator); //  aPp:      <P[i+1]P[i],P[j-1]P[j]>
                        if (has_next_O5i && has_O5i && has_prev_O5j && has_O5j) successtest+=bp_ang_score(next_O5i,O5i,prev_O5j,O5j,aO5p_mu,aO5p_sd, tol,weight_ang,nominator,denominator); // aO5p: <O5'[i+1]O5'[i],O5'[j-1]O5'[j]> 
                        if (has_next_C5i && has_C5i && has_prev_C5j && has_C5j) successtest+=bp_ang_score(next_C5i,C5i,prev_C5j,C5j,aC5p_mu,aC5p_sd, tol,weight_ang,nominator,denominator); // aC5p: <C5'[i+1]C5'[i],C5'[j-1]C5'[j]> 
                        if (has_next_C4i && has_C4i && has_prev_C4j && has_C4j) successtest+=bp_ang_score(next_C4i,C4i,prev_C4j,C4j,aC4p_mu,aC4p_sd, tol,weight_ang,nominator,denominator); // aC4p: <C4'[i+1]C4'[i],C4'[j-1]C4'[j]> 
                        if (has_next_C3i && has_C3i && has_prev_C3j && has_C3j) successtest+=bp_ang_score(next_C3i,C3i,prev_C3j,C3j,aC3p_mu,aC3p_sd, tol,weight_ang,nominator,denominator); // aC3p: <C3'[i+1]C3'[i],C3'[j-1]C3'[j]> 
                        if (has_next_C2i && has_C2i && has_prev_C2j && has_C2j) successtest+=bp_ang_score(next_C2i,C2i,prev_C2j,C2j,aC2p_mu,aC2p_sd, tol,weight_ang,nominator,denominator); // aC2p: <C2'[i+1]C2'[i],C2'[j-1]C2'[j]> 
                        if (has_next_C1i && has_C1i && has_prev_C1j && has_C1j) successtest+=bp_ang_score(next_C1i,C1i,prev_C1j,C1j,aC1p_mu,aC1p_sd, tol,weight_ang,nominator,denominator); // aC1p: <C1'[i+1]C1'[i],C1'[j-1]C1'[j]> 
                        if (has_next_O4i && has_O4i && has_prev_O4j && has_O4j) successtest+=bp_ang_score(next_O4i,O4i,prev_O4j,O4j,aO4p_mu,aO4p_sd, tol,weight_ang,nominator,denominator); // aO4p: <O4'[i+1]O4'[i],O4'[j-1]O4'[j]> 
                        if (has_next_O3i && has_O3i && has_prev_O3j && has_O3j) successtest+=bp_ang_score(next_O3i,O3i,prev_O3j,O3j,aO3p_mu,aO3p_sd, tol,weight_ang,nominator,denominator); // aO3p: <O3'[i+1]O3'[i],O3'[j-1]O3'[j]> 
                        if (has_C4i      && has_C1i && has_C4j      && has_C1j)//successtest+=bp_ang_score(    C4i,C1i,     C4j,C1j, aCC_mu, aCC_sd, tol,weight_ang,nominator,denominator); //  aCC:    <C4'[i]C1'[i],C4'[j]C1'[j]>
                        {
                            ang=rad2deg(Points4Angle(C4i,C1i,C4j,C1j));
                            if(ang>=-180)
                            {
                                successtest++;
                                denominator+=weight_ang; 
                                if (ang>aCC_mu) nominator+=weight_ang;
                                else nominator+=weight_ang*(1-fabs(ang-aCC_mu)/(tol*aCC_sd));
                            }
                        }
                    }
                    
                    if (denominator>0) bp_vec.push_back(
                        pair<double,vector<size_t> >(nominator/denominator+
                        adjust1*successtest/totaltest+adjust3, tmp_bp));
                }
            }
        }
        //cerr<<endl;
    }

    /* clean up */
    vector<size_t>().swap(tmp_bp);

    prev_Pi.clear();  prev_Pj.clear();
    prev_O5i.clear(); prev_O5j.clear();
    prev_C5i.clear(); prev_C5j.clear();
    prev_C4i.clear(); prev_C4j.clear();
    prev_C3i.clear(); prev_C3j.clear();
    prev_C2i.clear(); prev_C2j.clear();
    prev_C1i.clear(); prev_C1j.clear();
    prev_O4i.clear(); prev_O4j.clear();
    prev_O3i.clear(); prev_O3j.clear();
    prev_Nxi.clear(); prev_Nxj.clear();

    Pi.clear();  Pj.clear();
    O5i.clear(); O5j.clear();
    C5i.clear(); C5j.clear();
    C4i.clear(); C4j.clear();
    C3i.clear(); C3j.clear();
    C2i.clear(); C2j.clear();
    C1i.clear(); C1j.clear();
    O4i.clear(); O4j.clear();
    O3i.clear(); O3j.clear();
    Nxi.clear(); Nxj.clear();

    next_Pi.clear();  next_Pj.clear();  
    next_O5i.clear(); next_O5j.clear();
    next_C5i.clear(); next_C5j.clear();
    next_C4i.clear(); next_C4j.clear();
    next_C3i.clear(); next_C3j.clear();
    next_C2i.clear(); next_C2j.clear();
    next_C1i.clear(); next_C1j.clear();
    next_O4i.clear(); next_O4j.clear();
    next_O3i.clear(); next_O3j.clear();
    next_Nxi.clear(); next_Nxj.clear();
    return;
}

string getseq(ModelUnit pdb_entry){
    string seq="";
    for(int i=0;i<pdb_entry.chains[0].residues.size();i++){
        seq+=pdb_entry.chains[0].residues[i].resn[2];
    }
    return seq;
}

bool loadDistmap(string distpath, string seq, vector<vector<double>> &distmap, vector<vector<bool>> &distmapflag, int &distseqst, int &distseqend){
	int i,j;
    int count, count_row, count_col;
	string line, nt;
    bool flag1, flag2;
    int idx1, idx2;
	double val;
	double sc;
    bool flag;
    cout<<"Loading " << distpath << "..."<<endl;
	ifstream file(distpath);
	if (!file.is_open()) {
        cerr << "File of distance map cannot be opened [" << distpath << "]" << endl;
        return false;
    }

    int seq_len = seq.length();
    distmap.resize(seq_len, vector<double>(seq_len, 0));
    distmapflag.resize(seq_len, vector<bool>(seq_len, false));

    distseqst=0;
    distseqend=seq_len-1;
    int distrange=seq_len;

    int line_idx=0;
    count_row = 0;
	while (getline(file, line)) {
        istringstream stream(line);

        if (line_idx == 0 && line.find(",") != std::string::npos) {
            flag1=false; 
            flag2=false; 
            count = 0;
            while (getline(stream, nt, ',')) {
                if(count==0){
                    idx1=stoi(nt);
                    flag1=true;
                }
                if(count==1){
                    idx2=stoi(nt);
                    flag2=true;
                }
                count+=1;
            }   
            if(count == 2){
                distseqst=idx1-1;
                distseqend=idx2-1;

                if(distseqst < 0){
                    cerr << "Warning: The specified range is uncorrect. [" << distpath << " is ignored]" << endl;
                    distmap.clear();
                    distmapflag.clear();
                    return false;
                }
                if(distseqend+1 > seq_len){
                    cerr << "Warning: The specified range is uncorrect. [" << distpath << " is ignored]" << endl;
                    distmap.clear();
                    distmapflag.clear();
                    return false;
                }

                distrange=distseqend-distseqst+1;

                continue;
            }
        }else{
            if(distseqst+count_row > distseqend){
                distmap.clear();
                distmapflag.clear();
                cerr << "Warning: The range of the distance map is out of range [" << distpath << " is ignored]" << endl;
                return false;
            }

            count_col=0;
            while (stream >> nt) {
                val=stof(nt);

                if(distseqst+count_col > distseqend){
                    distmap.clear();
                    distmapflag.clear();
                    cerr << "Warning: The range of the distance map is out of range [" << distpath << " is ignored]" << endl;
                    return false;
                }

                distmap[distseqst+count_row][distseqst+count_col]=val;
                distmapflag[distseqst+count_row][distseqst+count_col]=true;

                //cout<<distseqst+count_row << "-" << distseqst+count_col << "-" << val<<endl;
                count_col+=1;
            }
        }
        count_row+=1;
        line_idx+=1;
    }

    if(count_col != distrange || count_row != distrange){
        cerr << "Warning: The range of the distance map does not match to the range [" << distpath << " is ignored]" << endl;
        distmap.clear();
        distmapflag.clear();
        return false;
    }

    for(i=0;i<distmap.size();i++){
        if(int(distmap[i][i]) != 0){
            cerr << "The diagonal elements of the distance map must be 0." << distpath << " is ignored]" << endl;
            distmap.clear();
            distmapflag.clear();
            return false;
        }

        for(j=i+1;j<distmap[i].size();j++){
            if(fabs(distmap[i][j] - distmap[j][i]) > 0.0001){
                cerr << "The distance map must be a symmetric matrix." << distpath << " is ignored]" << endl;
                distmap.clear();
                distmapflag.clear();
                return false;
            }
        }
    }

    file.close();

	return true;
}

bool loadContactmap(string contpath, string seq, vector<vector<double>> &contactmap, vector<vector<bool>> &contactmapflag, int &contactseqst, int &contactseqend){
	int i,j;
    int count, count_row, count_col;
	string line, nt;
    bool flag1, flag2;
    int idx1, idx2;
	double val;
	double sc;
    bool flag;
    cout<<"Loading " << contpath << "..."<<endl;
	ifstream file(contpath);
	if (!file.is_open()) {
        cerr << "File of contact map cannot be opened [" << contpath << "]" << endl;
        return false;
    }

    int seq_len = seq.length();
    contactmap.resize(seq_len, vector<double>(seq_len, 0));
    contactmapflag.resize(seq_len, vector<bool>(seq_len, false));

    contactseqst=0;
    contactseqend=seq_len-1;
    int contactrange=seq_len;

    int line_idx=0;
    count_row = 0;
	while (getline(file, line)) {
        istringstream stream(line);

        if (line_idx == 0 && line.find(",") != std::string::npos) {
            flag1=false; 
            flag2=false; 
            count = 0;
            while (getline(stream, nt, ',')) {
                if(count==0){
                    idx1=stoi(nt);
                    flag1=true;
                }
                if(count==1){
                    idx2=stoi(nt);
                    flag2=true;
                }
                count+=1;
            }   
            if(count == 2){
                contactseqst=idx1-1;
                contactseqend=idx2-1;

                if(contactseqst < 0){
                    cerr << "Warning: The specified range is uncorrect. [" << contpath << " is ignored]" << endl;
                    contactmap.clear();
                    contactmapflag.clear();
                    return false;
                }
                if(contactseqend+1 > seq_len){
                    cerr << "Warning: The specified range is uncorrect. [" << contpath << " is ignored]" << endl;
                    contactmap.clear();
                    contactmapflag.clear();
                    return false;
                }

                contactrange=contactseqend-contactseqst+1;

                continue;
            }
        }else{
            if(contactseqst+count_row > contactseqend){
                contactmap.clear();
                contactmapflag.clear();
                cerr << "Warning: The range of the contact map is out of range [" << contpath << " is ignored]" << endl;
                return false;
            }

            count_col=0;
            while (stream >> nt) {
                val=stof(nt);

                if(contactseqst+count_col > contactseqend){
                    contactmap.clear();
                    contactmapflag.clear();
                    cerr << "Warning: The range of the contact map is out of range [" << contpath << " is ignored]" << endl;
                    return false;
                }

                contactmap[contactseqst+count_row][contactseqst+count_col]=val;
                contactmapflag[contactseqst+count_row][contactseqst+count_col]=true;

                //cout<<contactseqst+count_row << "-" << contactseqst+count_col << "-" << val<<endl;
                count_col+=1;
            }
        }
        count_row+=1;
        line_idx+=1;
    }

    if(count_col != contactrange || count_row != contactrange){
        cerr << "Warning: The range of the contact map does not match to the range [" << contpath << " is ignored]" << endl;
        contactmap.clear();
        contactmapflag.clear();
        return false;
    }

    for(i=0;i<contactmap.size();i++){
        for(j=i+1;j<contactmap[i].size();j++){
            if(fabs(contactmap[i][j] - contactmap[j][i]) > 0.0001){
                cerr << "The contact map must be a symmetric matrix." << contpath << " is ignored]" << endl;
                contactmap.clear();
                contactmapflag.clear();
                return false;
            }

            //if(contactmapflag[i][j]) cout<< i << "-" << j << "-" << contactmap[i][j]<<endl;
        }
    }

    file.close();

	return true;
}

/*
bool loadST(string stpath, string seq, vector<pair<double,vector<size_t> > > &stack_vec, int &stackseqst, int &stackseqend){
	int count;
	string line, nt;
	bool flag1, flag2, flag3;
	int idx1, idx2;
	double sc;
    cout<<"Loading " << stpath << "..."<<endl;
	ifstream file(stpath);
	if (!file.is_open()) {
        cerr << "stacking information cannot be opened [" << stpath << "]" << endl;
        return false;
    }

    vector<vector<bool>> stackmat(seq.length(), vector<bool>(seq.length(), false));

    stackseqst=0;
    stackseqend=seq.length()-1;

    int line_idx = 0;
	while (getline(file, line)) {
        istringstream stream(line);

        if (line_idx == 0 && line.find(",") != std::string::npos) {
            flag1=false; 
            flag2=false; 
            count = 0;
            while (getline(stream, nt, ',')) {
                if(count==0){
                    idx1=stoi(nt);
                    flag1=true;
                }
                if(count==1){
                    idx2=stoi(nt);
                    flag2=true;
                }
                count+=1;
            }   

            if(count == 2){
                stackseqst=idx1-1;
                stackseqend=idx2-1;
		continue;
            }
        }else{

            flag1=false;
            flag2=false;
            flag3=false;        
            count = 0;
            sc=1;
            while (stream >> nt) {
                if(count==0){
                    idx1=stoi(nt);
                    flag1=true;
                }
                if(count==1){
                    idx2=stoi(nt);
                    flag2=true;
                }
                if(count==2){
                    sc=stof(nt);
                    flag3=true;
                }
                count+=1;
            }

            if(flag1 & flag2 & flag3){
                pair<double,vector<size_t> > stack_vec_temp;
                //if(idx1 < idx2){
                stack_vec_temp.first = sc;
                stack_vec_temp.second.resize(5,0);

                stack_vec_temp.second[0]=0;
                stack_vec_temp.second[1]=idx1-1;
                stack_vec_temp.second[2]=0;
                stack_vec_temp.second[3]=idx2-1;

                if((idx1 == idx2 + 1) || (idx2 == idx1 + 1)) stack_vec_temp.second[4]=1;

                if(!stackmat[stack_vec_temp.second[1]][stack_vec_temp.second[3]]){
                    stack_vec.push_back(stack_vec_temp);
                    stackmat[stack_vec_temp.second[1]][stack_vec_temp.second[3]]=true;
                    stackmat[stack_vec_temp.second[3]][stack_vec_temp.second[1]]=true;
                }
            }
        }
        line_idx+=1;
    }

    //for(int i=0;i<stack_vec.size();i++){
        //cout<<stack_vec[i].first << "-" << stack_vec[i].second[0] << "-" << stack_vec[i].second[1] << "-" << stack_vec[i].second[2] << "-" << stack_vec[i].second[3]<<endl;
    //}

    if(stackseqst < 0){
        cerr << "Warning: The specified range is uncorrect. [" << stpath << " is ignored]" << endl;
        stack_vec.clear();
        return false;
    }
    if(stackseqend > seq.length()-1){
        cerr << "Warning: The specified range is uncorrect. [" << stpath << " is ignored]" << endl;
        stack_vec.clear();
        return false;
    }
    
    // ファイルを閉じる
    file.close();

	return true;
}


bool loadBP(string bppath, string seq, vector<pair<double,vector<size_t>>> &bp_vec, vector<pair<double,vector<size_t>>> &noncbp_vec, int &bpseqst, int &bpseqend){
	int count;
	string line, nt;
	bool flag1, flag2, flag3, flag4;
	int idx1, idx2;
    string nt1;
    char base1, base2;
	double sc;
	//BPBasicInfo basepair;
    cout<<"Loading " << bppath << "..."<<endl;
	ifstream file(bppath);
	if (!file.is_open()) {
        cerr << "base paring information cannot be opened [" << bppath << "]" << endl;
        return false;
    }

    bpseqst=0;
    bpseqend=seq.length()-1;
    
    vector<vector<bool>> bpmat(seq.length(), vector<bool>(seq.length(), false));

    int line_idx = 0;
	while (getline(file, line)) {
        istringstream stream(line);

        if (line_idx == 0 && line.find(",") != std::string::npos) {
            flag1=false; 
            flag2=false; 
            count = 0;
            while (getline(stream, nt, ',')) {
                if(count==0){
                    idx1=stoi(nt);
                    flag1=true;
                }
                if(count==1){
                    idx2=stoi(nt);
                    flag2=true;
                }
                count+=1;
            }   

            if(count == 2){
                bpseqst=idx1-1;
                bpseqend=idx2-1;

                continue;
            }
        }else{
            flag1=false; 
            flag2=false; 
            flag3=false;
            flag4=false;
            count = 0;
            sc=1;
            while (stream >> nt) 
            {
                if(count==0){
                    idx1=stoi(nt);
                    flag1=true;
                }
                if(count==1){
                    nt1=nt;
                    flag2=true;
                }
                if(count==2){
                    idx2=stoi(nt);
                    flag3=true;
                }
                if(count==3){
                    sc=stof(nt);
                    flag4=true;
                }
                count+=1;
            }

            if(idx2 <= 0){
                line_idx+=1;
                continue;
            }

            if(flag1 & flag2 & flag3){
                pair<double,vector<size_t>> bp_vec_temp;
                bool canonicaflag=false;
                if(idx1 < idx2){
                    bp_vec_temp.first = sc;
                    bp_vec_temp.second.resize(5, 0.0);

                    bp_vec_temp.second[0] = 0;
                    bp_vec_temp.second[1] = idx1-1;
                    bp_vec_temp.second[2] = 0;
                    bp_vec_temp.second[3] = idx2-1;

                    base1 = seq[bp_vec_temp.second[1]];
                    base2 = seq[bp_vec_temp.second[3]];

                    if ((base1=='A' &&(base2=='U' || base2=='T')) || (base1=='C' && base2=='G')|| ((base1=='U' || base1=='T')&& base2=='A')  || (base1=='G' && base2=='C')){
                        bp_vec_temp.second[4]=0;
                        canonicaflag=true;
                    }
                    else if ((base1=='G' && base2=='U') || (base1=='U' && base2=='G')){
                        bp_vec_temp.second[4]=1;
                        canonicaflag=true;
                    }else{
                        bp_vec_temp.second[4]=2;
                        canonicaflag=false;
                    }
                    //else{
                        //cout<<bp_vec_temp.second[1] << "-" << bp_vec_temp.second[3]<<endl;
                    //}
                }else{
                    bp_vec_temp.first = sc;
                    bp_vec_temp.second.resize(5, 0.0);

                    bp_vec_temp.second[0] = 0;
                    bp_vec_temp.second[1] = idx2-1;
                    bp_vec_temp.second[2] = 0;
                    bp_vec_temp.second[3] = idx1-1;

                    base1 = seq[bp_vec_temp.second[1]];
                    base2 = seq[bp_vec_temp.second[3]];

                    if ((base1=='A' &&(base2=='U' || base2=='T')) || (base1=='C' && base2=='G')|| ((base1=='U' || base1=='T')&& base2=='A') || (base1=='G' && base2=='C')){
                        bp_vec_temp.second[4]=0;
                        canonicaflag=true;
                    }
                    else if ((base1=='G' && base2=='U') || (base1=='U' && base2=='G')){
                        bp_vec_temp.second[4]=1;
                        canonicaflag=true;
                    }else{
                        bp_vec_temp.second[4]=2;
                        canonicaflag=false;
                    }
                }
                if(canonicaflag){
                    if(!bpmat[bp_vec_temp.second[1]][bp_vec_temp.second[3]]){
                        bp_vec.push_back(bp_vec_temp);
                        bpmat[bp_vec_temp.second[1]][bp_vec_temp.second[3]]=true;
                        bpmat[bp_vec_temp.second[3]][bp_vec_temp.second[1]]=true;
                    }
                }else{
                    if(!bpmat[bp_vec_temp.second[1]][bp_vec_temp.second[3]]){
                        noncbp_vec.push_back(bp_vec_temp);
                        bpmat[bp_vec_temp.second[1]][bp_vec_temp.second[3]]=true;
                        bpmat[bp_vec_temp.second[3]][bp_vec_temp.second[1]]=true;
                    }
                }
            }
        }

        line_idx+=1;
    }

    //for(int i=0;i<bp_vec.size();i++){
        //cout<<bp_vec[i].first << "-" << bp_vec[i].second[0] << "-" << bp_vec[i].second[1] << "-" << bp_vec[i].second[2] << "-" << bp_vec[i].second[3]<<endl;
    //}

    if(bpseqst < 0){
        cerr << "Warning: The specified range is uncorrect. [" << bppath << " is ignored]" << endl;
        bp_vec.clear();
        noncbp_vec.clear();
        return false;
    }
    if(bpseqend > seq.length()-1){
        cerr << "Warning: The specified range is uncorrect. [" << bppath << " is ignored]" << endl;
        bp_vec.clear();
        noncbp_vec.clear();
        return false;
    }

    // ファイルを閉じる
    file.close();

	return true;
}
*/

double exp01(double x, double k, bool clamp_x = true) {
    if (clamp_x) {
        x = std::min(1.0, std::max(0.0, x));
    }
    // 極小 k の処理（分母が 0 に近くなるのを回避）
    const double eps = 1e-12;
    if (std::abs(k) < eps) {
        // lim_{k->0} (exp(k*x)-1)/(exp(k)-1) = x
        return x;
    }
    // expm1(z) = exp(z) - 1 : 小さい z で精度良く計算
    const double num = std::expm1(k * x);
    const double den = std::expm1(k);

    //cout<< x << "-" << num / den<<endl;

    return num / den;
}

bool loadST(string stpath, string seq, map<string, int> &chainstr2chainidx, map<int, map<int, int>> &pdbidx2seqidx, vector<pair<double,vector<size_t> > > &stack_vec, int &stackseqst, int &stackseqend){
	int count;
	string line, nt;
	bool flag1, flag2, flag3, flag4, flag5;
    string chain1, chain2;
    int chain1idx, chain2idx;
    int pdbidx1, pdbidx2;
	int idx1, idx2;
    char base1, base2;
	double sc;
    cout<<"Loading " << stpath << "..."<<endl;
	ifstream file(stpath);
	if (!file.is_open()) {
        cerr << "stacking information cannot be opened [" << stpath << "]" << endl;
        return false;
    }

    vector<vector<bool>> stackmat(seq.length(), vector<bool>(seq.length(), false));

    stackseqst=0;
    stackseqend=seq.length()-1;

	while (getline(file, line)) {
        istringstream stream(line);

        flag1=false; flag2=false; flag3=false; flag4=false; flag5=false;
        count = 0;
        sc=1.0;
        while (std::getline(stream, nt, ','))
        {
            if(count==0){
                chain1 = nt;
                if(chainstr2chainidx.find(chain1) != chainstr2chainidx.end()) {
                    chain1idx = chainstr2chainidx[nt];
                    flag1=true;
                }else{
                    cout << "Warning: chain " << chain1 << " is not included in the inputted pdb structure." << endl;
                    exit(1);
                }
            }
            if(count==1){
                chain2 = nt;
                if(chainstr2chainidx.find(chain2) != chainstr2chainidx.end()) {
                    chain2idx = chainstr2chainidx[nt];
                    flag2=true;
                }else{
                    cout << "Warning: chain " << chain2 << " is not included in the inputted pdb structure." << endl;
                    exit(1);
                }
            }
            if(count==2){
                pdbidx1 = stoi(nt);
                if(pdbidx2seqidx[chain1idx].find(pdbidx1) != pdbidx2seqidx[chain1idx].end()) {
                    idx1 = pdbidx2seqidx[chain1idx][pdbidx1];
                    flag3=true;
                }else{
                    cout << "Warning: Sequence idx " << pdbidx1 << " is not included in the inputted pdb structure." << endl;
                    exit(1);
                }
            }
            if(count==3){
                pdbidx2 = stoi(nt);
                if(pdbidx2seqidx[chain2idx].find(pdbidx2) != pdbidx2seqidx[chain2idx].end()) {
                    idx2 = pdbidx2seqidx[chain2idx][pdbidx2];
                    flag4=true;
                }else{
                    cout << "Warning: Sequence idx " << pdbidx2 << " is not included in the inputted pdb structure." << endl;
                    exit(1);
                }
            }
            if(count==4){
                sc=stof(nt);
                flag5=true;
            }
            count+=1;
        }

        if(flag1 & flag2 & flag3 & flag4){
            pair<double,vector<size_t> > stack_vec_temp;
            
            if(idx1 < idx2){
                stack_vec_temp.first = exp01(sc, 2.0);
                stack_vec_temp.second.resize(5,0);

                stack_vec_temp.second[0]=chain1idx;
                stack_vec_temp.second[1]=idx1;
                stack_vec_temp.second[2]=chain2idx;
                stack_vec_temp.second[3]=idx2;
            }else{
                stack_vec_temp.first = exp01(sc, 2.0);
                stack_vec_temp.second.resize(5,0);

                stack_vec_temp.second[0]=chain2idx;
                stack_vec_temp.second[1]=idx2;
                stack_vec_temp.second[2]=chain1idx;
                stack_vec_temp.second[3]=idx1;
            }

            //cout<<stack_vec_temp.second[0] << "-" << stack_vec_temp.second[1] << "-" << stack_vec_temp.second[2] << "-" << stack_vec_temp.second[3] << "-" << stack_vec_temp.first<<endl;

            if((idx1 == idx2 + 1) || (idx2 == idx1 + 1)) stack_vec_temp.second[4]=1;

            if(!stackmat[stack_vec_temp.second[1]][stack_vec_temp.second[3]]){
                stack_vec.push_back(stack_vec_temp);
                stackmat[stack_vec_temp.second[1]][stack_vec_temp.second[3]]=true;
                stackmat[stack_vec_temp.second[3]][stack_vec_temp.second[1]]=true;
            }
        }
    }

    //for(int i=0;i<stack_vec.size();i++){
        //cout<<stack_vec[i].first << "-" << stack_vec[i].second[0] << "-" << stack_vec[i].second[1] << "-" << stack_vec[i].second[2] << "-" << stack_vec[i].second[3]<<endl;
    //}

    if(stackseqst < 0){
        cerr << "Warning: The specified range is uncorrect. [" << stpath << " is ignored]" << endl;
        stack_vec.clear();
        return false;
    }
    if(stackseqend > seq.length()-1){
        cerr << "Warning: The specified range is uncorrect. [" << stpath << " is ignored]" << endl;
        stack_vec.clear();
        return false;
    }
    
    // ファイルを閉じる
    file.close();

	return true;
}

bool loadBP(string bppath, string seq, map<string, int> &chainstr2chainidx, map<int, map<int, int>> &pdbidx2seqidx, vector<pair<double,vector<size_t>>> &bp_vec, vector<pair<double,vector<size_t>>> &noncbp_vec, int &bpseqst, int &bpseqend){
	int count;
	string line, nt;
	bool flag1, flag2, flag3, flag4, flag5;
    string chain1, chain2;
    int chain1idx, chain2idx;
    int pdbidx1, pdbidx2;
	int idx1, idx2;
    char base1, base2;
	double sc;
	//BPBasicInfo basepair;
    cout<<"Loading " << bppath << "..."<<endl;
	ifstream file(bppath);
	if (!file.is_open()) {
        cerr << "base paring information cannot be opened [" << bppath << "]" << endl;
        return false;
    }

    bpseqst=0;
    bpseqend=seq.length()-1;
    
    vector<vector<bool>> bpmat(seq.length(), vector<bool>(seq.length(), false));

	while (getline(file, line)) {
        istringstream stream(line);

        flag1=false; flag2=false; flag3=false; flag4=false; flag5=false;
        count = 0;
        sc=1.0;
        while (std::getline(stream, nt, ','))
        {
            if(count==0){
                chain1 = nt;
                if(chainstr2chainidx.find(chain1) != chainstr2chainidx.end()) {
                    chain1idx = chainstr2chainidx[nt];
                    flag1=true;
                }else{
                    cout << "Warning: chain " << chain1 << " is not included in the inputted pdb structure." << endl;
                    exit(1);
                }
            }
            if(count==1){
                chain2 = nt;
                if(chainstr2chainidx.find(chain2) != chainstr2chainidx.end()) {
                    chain2idx = chainstr2chainidx[nt];
                    flag2=true;
                }else{
                    cout << "Warning: chain " << chain2 << " is not included in the inputted pdb structure." << endl;
                    exit(1);
                }
            }
            if(count==2){
                pdbidx1 = stoi(nt);
                if(pdbidx2seqidx[chain1idx].find(pdbidx1) != pdbidx2seqidx[chain1idx].end()) {
                    idx1 = pdbidx2seqidx[chain1idx][pdbidx1];
                    flag3=true;
                }else{
                    cout << "Warning: Sequence idx " << pdbidx1 << " is not included in the inputted pdb structure." << endl;
                    exit(1);
                }
            }
            if(count==3){
                pdbidx2 = stoi(nt);
                if(pdbidx2seqidx[chain2idx].find(pdbidx2) != pdbidx2seqidx[chain2idx].end()) {
                    idx2 = pdbidx2seqidx[chain2idx][pdbidx2];
                    flag4=true;
                }else{
                    cout << "Warning: Sequence idx " << pdbidx2 << " is not included in the inputted pdb structure." << endl;
                    exit(1);
                }
            }
            if(count==4){
                sc=stof(nt);
                flag5=true;
            }
            count+=1;
        }

        if(flag1 & flag2 & flag3 & flag4){
            pair<double,vector<size_t>> bp_vec_temp;
            bool canonicaflag=false;
            if(idx1 < idx2){
                bp_vec_temp.first = exp01(sc, 2.0);
                bp_vec_temp.second.resize(5, 0.0);

                bp_vec_temp.second[0] = chain1idx;
                bp_vec_temp.second[1] = idx1;
                bp_vec_temp.second[2] = chain2idx;
                bp_vec_temp.second[3] = idx2;

                base1 = seq[bp_vec_temp.second[1]];
                base2 = seq[bp_vec_temp.second[3]];

                //cout<<idx1 << "-" << idx2 << "-" << base1 << "-" << base2<<endl;

                //cout<<bp_vec_temp.second[0] << "-" << bp_vec_temp.second[1] << "-" << bp_vec_temp.second[2] << "-" << bp_vec_temp.second[3] << "-" << base1 << "-" << base2 << "-" << bp_vec_temp.first<<endl;

                if ((base1=='A' &&(base2=='U' || base2=='T')) || (base1=='C' && base2=='G')|| ((base1=='U' || base1=='T')&& base2=='A')  || (base1=='G' && base2=='C')){
                    bp_vec_temp.second[4]=0;
                    canonicaflag=true;
                }
                else if ((base1=='G' && base2=='U') || (base1=='U' && base2=='G')){
                    bp_vec_temp.second[4]=1;
                    canonicaflag=true;
                }else{
                    bp_vec_temp.second[4]=2;
                    canonicaflag=false;
                }
                //else{
                    //cout<<bp_vec_temp.second[1] << "-" << bp_vec_temp.second[3]<<endl;
                //}
            }else{
                bp_vec_temp.first = exp01(sc, 2.0);
                bp_vec_temp.second.resize(5, 0.0);

                bp_vec_temp.second[0] = chain2idx;
                bp_vec_temp.second[1] = idx2;
                bp_vec_temp.second[2] = chain1idx;
                bp_vec_temp.second[3] = idx1;

                base1 = seq[bp_vec_temp.second[1]];
                base2 = seq[bp_vec_temp.second[3]];

                //cout<<idx2 << "-" << idx1 << "-" << base1 << "-" << base2<<endl;

                //cout<<bp_vec_temp.second[0] << "-" << bp_vec_temp.second[1] << "-" << bp_vec_temp.second[2] << "-" << bp_vec_temp.second[3] << "-" << base1 << "-" << base2 << "-" << bp_vec_temp.first<<endl;

                if ((base1=='A' &&(base2=='U' || base2=='T')) || (base1=='C' && base2=='G')|| ((base1=='U' || base1=='T')&& base2=='A') || (base1=='G' && base2=='C')){
                    bp_vec_temp.second[4]=0;
                    canonicaflag=true;
                }
                else if ((base1=='G' && base2=='U') || (base1=='U' && base2=='G')){
                    bp_vec_temp.second[4]=1;
                    canonicaflag=true;
                }else{
                    bp_vec_temp.second[4]=2;
                    canonicaflag=false;
                }
            }
            if(canonicaflag){
                if(!bpmat[bp_vec_temp.second[1]][bp_vec_temp.second[3]]){
                    bp_vec.push_back(bp_vec_temp);
                    bpmat[bp_vec_temp.second[1]][bp_vec_temp.second[3]]=true;
                    bpmat[bp_vec_temp.second[3]][bp_vec_temp.second[1]]=true;
                }
            }else{
                if(!bpmat[bp_vec_temp.second[1]][bp_vec_temp.second[3]]){
                    noncbp_vec.push_back(bp_vec_temp);
                    bpmat[bp_vec_temp.second[1]][bp_vec_temp.second[3]]=true;
                    bpmat[bp_vec_temp.second[3]][bp_vec_temp.second[1]]=true;
                }
            }
        }
    }

    //for(int i=0;i<bp_vec.size();i++){
        //cout<<bp_vec[i].first << "-" << bp_vec[i].second[0] << "-" << bp_vec[i].second[1] << "-" << bp_vec[i].second[2] << "-" << bp_vec[i].second[3]<<endl;
    //}

    if(bpseqst < 0){
        cerr << "Warning: The specified range is uncorrect. [" << bppath << " is ignored]" << endl;
        bp_vec.clear();
        noncbp_vec.clear();
        return false;
    }
    if(bpseqend > seq.length()-1){
        cerr << "Warning: The specified range is uncorrect. [" << bppath << " is ignored]" << endl;
        bp_vec.clear();
        noncbp_vec.clear();
        return false;
    }

    // ファイルを閉じる
    file.close();

	return true;
}


void getsselement(const ModelUnit &pdb_entry, vector<pair<double,vector<size_t> > > bp_vec, vector<vector<vector<int>>> &stembps, map<int, vector<pair<int, int>>> &loopnts, vector<vector<int>> &stemconnectnts, map<int, vector<vector<vector<int>>>> &connectnts){
    vector<vector<int>> stembp;
    int prev_paired;
    int i;
    int ci, cj, c1, c2, r1, r2, bp;
    int istart, iend, istart_pair, iend_pair, istart_pair_c, iend_pair_c;
    string spt;

    for(ci=0;ci<pdb_entry.chains.size();ci++)
    {
        for(cj=ci;cj<pdb_entry.chains.size();cj++)
        {
            map<int, int> bppairednt;
            //vector<int> loopresidues, loopnt;
            if(ci == cj){
                for (r1=0; r1<pdb_entry.chains[ci].residues.size(); r1++) bppairednt[r1] = -1;
                for (bp=0;bp<bp_vec.size();bp++)
                {
                    c1= bp_vec[bp].second[0];
                    r1= bp_vec[bp].second[1];
                    c2= bp_vec[bp].second[2];
                    r2= bp_vec[bp].second[3];

                    if(ci == c1 && cj == c2){
                        bppairednt[r1]=r2;
                        bppairednt[r2]=r1;
                    }
                }

                prev_paired = -1;
                for(i=0;i<pdb_entry.chains[ci].residues.size();i++)
                {
                    if(bppairednt[i] >= 0 && i < bppairednt[i]){
                        if(prev_paired < 0){
                            stembp.clear();
                            stembp.push_back({ci, i, cj, bppairednt[i]});
                            prev_paired = bppairednt[i];
                        }else if(prev_paired >= 0 && prev_paired-1 != bppairednt[i]){
                            stembps.push_back(stembp);
                            stembp.clear();
                            stembp.push_back({ci, i, cj, bppairednt[i]});
                            prev_paired = bppairednt[i];
                        }else if(prev_paired >= 0 && prev_paired-1 == bppairednt[i]){
                            stembp.push_back({ci, i, cj, bppairednt[i]});
                            prev_paired = bppairednt[i];
                        }
                    }else{
                        if(prev_paired >= 0){
                            stembps.push_back(stembp);
                            stembp.clear();
                            prev_paired = -1;
                        }
                    }
                }
                if(prev_paired >= 0)
                {
                    stembps.push_back(stembp);
                    stembp.clear();
                }

            }else{
                for (r1=0; r1<pdb_entry.chains[ci].residues.size()+pdb_entry.chains[cj].residues.size(); r1++) bppairednt[r1] = -1;
                int seqlen_c1 = pdb_entry.chains[ci].residues.size();

                for (bp=0;bp<bp_vec.size();bp++)
                {
                    c1= bp_vec[bp].second[0];
                    r1= bp_vec[bp].second[1];
                    c2= bp_vec[bp].second[2];
                    r2= bp_vec[bp].second[3];

                    if(ci == c1 && cj == c2){
                        bppairednt[r1]=seqlen_c1+r2;
                        bppairednt[seqlen_c1+r2]=r1;
                    }

                    if(ci == c2 && cj == c1){
                        bppairednt[seqlen_c1+r1]=r2;
                        bppairednt[r2]=seqlen_c1+r1;
                    }
                }

                prev_paired = -1;
                for(i=0;i<pdb_entry.chains[ci].residues.size();i++)
                {
                    if(bppairednt[i] >= 0)
                    {
                        if(prev_paired < 0){
                            stembp.clear();
                            stembp.push_back({ci, i, cj, bppairednt[i]-seqlen_c1});
                            prev_paired = bppairednt[i] - seqlen_c1;
                        }else if(prev_paired >= 0 && prev_paired-1 != bppairednt[i]-seqlen_c1){
                            stembps.push_back(stembp);
                            stembp.clear();
                            stembp.push_back({ci, i, cj, bppairednt[i]-seqlen_c1});
                            prev_paired = bppairednt[i];
                        }else if(prev_paired >= 0 && prev_paired-1 == bppairednt[i]-seqlen_c1){
                            stembp.push_back({ci, i, cj, bppairednt[i]-seqlen_c1});
                            prev_paired = bppairednt[i]-seqlen_c1;
                        }
                    }else{
                        if(prev_paired >= 0){
                            stembps.push_back(stembp);
                            stembp.clear();
                            prev_paired = -1;
                        }
                    }
                }
                if(prev_paired >= 0)
                {
                    stembps.push_back(stembp);
                    stembp.clear();
                }
            }
        }
    }

    for(int i=0;i<stembps.size();i++)
    {
        int bp_len = stembps[i].size()-1;
        int baseleftprev_c = stembps[i][0][0];
        int baseleftprev_r = stembps[i][0][1]-1;
        int baseleftnext_c = stembps[i][bp_len][0];
        int baseleftnext_r = stembps[i][bp_len][1]+1;

        int baserightprev_c = stembps[i][bp_len][2];
        int baserightprev_r = stembps[i][bp_len][3]-1;
        int baserightnext_c = stembps[i][0][2];
        int baserightnext_r = stembps[i][0][3]+1;

        vector<int> stemconnectnt(8, -1);
        if(0 < baseleftprev_r && baseleftprev_r < pdb_entry.chains[baseleftprev_c].residues.size())
        {
            stemconnectnt[0] = baseleftprev_c;
            stemconnectnt[1] = baseleftprev_r;
        }
        if(0 < baseleftnext_r && baseleftnext_r < pdb_entry.chains[baseleftnext_c].residues.size())
        {
            stemconnectnt[2] = baseleftnext_c;
            stemconnectnt[3] = baseleftnext_r;
        }
        if(0 < baserightprev_r && baserightprev_r < pdb_entry.chains[baserightprev_c].residues.size())
        {
            stemconnectnt[4] = baserightprev_c;
            stemconnectnt[5] = baserightprev_r;
        }
        if(0 < baserightnext_r && baserightnext_r < pdb_entry.chains[baserightnext_c].residues.size())
        {
            stemconnectnt[6] = baserightnext_c;
            stemconnectnt[7] = baserightnext_r;
        }
        stemconnectnts.push_back(stemconnectnt);
    }

    for(ci=0;ci<pdb_entry.chains.size();ci++)
    {
        vector<int> bppairednt(pdb_entry.chains[ci].residues.size(), -100);
        vector<int> bppairedchain(pdb_entry.chains[ci].residues.size(), -100);
        vector<int> loopnt(pdb_entry.chains[ci].residues.size(),true);
        for (bp=0;bp<bp_vec.size();bp++)
        {
            c1= bp_vec[bp].second[0];
            r1= bp_vec[bp].second[1];
            c2= bp_vec[bp].second[2];
            r2= bp_vec[bp].second[3];

            //if(ci == c1) bppairednt[r1]=r2;
            //if(cj == c2) bppairednt[r2]=r1;
            if(ci == c1){
                bppairedchain[r1] = c2;
                bppairednt[r1] = r2;
                loopnt[r1] = false;
            }
            if(ci == c2){
                bppairedchain[r2] = c1;
                bppairednt[r2] = r1;
                loopnt[r2] = false;
            }
        }
    
        istart = -1; iend = -1;
        for(i=0; i<pdb_entry.chains[ci].residues.size(); i++){
            if(loopnt[i] && istart == -1){
                istart = i;
            }
            else if(!loopnt[i] && istart != -1){
                iend=i-1;
                loopnts[ci].push_back(make_pair(istart, iend));
                istart = -1;
            }
        }
        if(istart != -1){
            loopnts[ci].push_back(make_pair(istart, pdb_entry.chains[ci].residues.size()-1));
        }

        connectnts[ci].resize(loopnts[ci].size());
        for(i=0;i<loopnts[ci].size();i++){
            istart = loopnts[ci][i].first-1;
            iend = loopnts[ci][i].second+1;

            auto it1 = find(bppairednt.begin(), bppairednt.end(), istart);
            istart_pair = -1;
            istart_pair_c = -1;
            if(it1 != bppairednt.end()){
                istart_pair = bppairednt[istart];
                istart_pair_c = bppairedchain[istart];
            }

            auto it2 = find(bppairednt.begin(), bppairednt.end(), iend);
            iend_pair = -1;
            iend_pair_c = -1;
            if(it2 != bppairednt.end()){
                iend_pair = bppairednt[iend];
                iend_pair_c = bppairedchain[iend];
            }

            if(istart_pair != -1 && istart == iend_pair){
                connectnts[ci][i].push_back({istart, iend, istart_pair_c, iend_pair_c, 2});
                connectnts[ci][i].push_back({-1, -1, -1, -1, 0});
            }else if(istart_pair != -1 && iend_pair != -1){
                connectnts[ci][i].push_back({istart, istart_pair, ci, istart_pair_c, 1});
                connectnts[ci][i].push_back({iend, iend_pair, ci, iend_pair_c, 1});
            }else if(istart_pair != -1){
                connectnts[ci][i].push_back({istart, istart_pair, ci, istart_pair_c,1});
                connectnts[ci][i].push_back({-1, -1, -1, -1, 0});
            }else if(iend_pair != -1){
                connectnts[ci][i].push_back({-1, -1, -1, -1, 0});
                connectnts[ci][i].push_back({iend, iend_pair, ci, iend_pair_c, 1});
            }else{
                connectnts[ci][i].push_back({-1, -1, -1, -1, 0});
                connectnts[ci][i].push_back({-1, -1, -1, -1, 0});
            }
        }
    }
}   

void getstmap(const vector<pair<double,vector<size_t>>> &stack_vec, vector<vector<pair<int, double>>> &stackmap, const int &length){
    int i,j;
    size_t base1, base2, sttype;
    
    stackmap.resize(length, vector<pair<int, double>>(length));

    for(i=0;i<length;i++){
        for(j=0;j<length;j++){
            stackmap[i][j].first = 0;
            stackmap[i][j].second = 0.0;
        }
    }

    for(i=0;i<stack_vec.size();i++){
        base1=stack_vec[i].second[1];
        base2=stack_vec[i].second[3];
        sttype=stack_vec[i].second[4]+1;

        stackmap[base1][base2].first = sttype;
        stackmap[base2][base1].first = sttype;
        stackmap[base1][base2].second = stack_vec[i].first;
        stackmap[base2][base1].second = stack_vec[i].first;
    }
}

void get_edge_loop(vector<pair<double,vector<size_t>>> &bp_vec, vector<pair<double,vector<size_t>>> &noncbp_vec, int numbb, vector<double> &nt_weight)
{
    int st, ed;
    int idx1, idx2;
    st = numbb;
    ed = 0;

    nt_weight.resize(numbb, 1.0);
    if(bp_vec.size() > 0){
        for(int i=0;i<bp_vec.size();i++){
            idx1 = bp_vec[i].second[1];
            idx2 = bp_vec[i].second[3];

            //nt_weight[idx1] = 1;
            //nt_weight[idx2] = 1;

            //if(bp_vec[i].first > 0.5)
            //{
            if(st > idx1) st = idx1;
            if(st > idx2) st = idx2;
            if(ed < idx1) ed = idx1;
            if(ed < idx2) ed = idx2;
            //}else{
                //cout<<idx1 << "-" << idx2<<endl;
            //}
        }

        for(int i=0;i<st;i++){
            //cout<<i<<endl;
            nt_weight[i]=0.1;
        }

        for(int i=numbb-1;i>ed;i--){
            //cout<<i<<endl;
            nt_weight[i]=0.1;
        }

    }

    //for(int i=0;i<nt_weight.size();i++){
        //cout<<i << "-" << nt_weight[i]<<endl;
    //}

    //exit(1);
}






























