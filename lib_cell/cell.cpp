#include "cell.hpp"
#include "const.hpp"

const double CCell::vc           = -80;
const double CCell::stim         = 12.5;
const double CCell::stimduration = 5;

CCell::CCell()  {

    y = new double[101];
    ydot = new double[101];
    // initial conditions
    dt    = 0.1;
    
    double tmp[] = {0, 0.0202919327386704, 0.936490564532038, 0.936304480468966, 5.42243465425516e-06, 0.999659741420995, 0.0387439718487863, 0.0286342994782812,
                    0.00405157400000000, 0.994551100000000, 0.000648056012992726, 0.973061283727168, 0.00193864192924699, 0.00419293036123511, 0.809901538287505,
                    1.32068131238543e-06, 3.09984986234518e-07, 3.76351356810290, 0.821291689576353, 0.0163997256637698, 0.124384775412331, 0.00737380343482592,
                    0.000611244439499118, 0.00302488782583208, 0.136472047162854, 0.00399721388615559, 0.0117402146171874, 0.0188546757489768, 0.0962772087252904,
                    0.176838896529319, 1.05807454671198, 0.446032394313757, 9.90313914499844, 9.90459022726698, 9.90483249546865, 120, 0.000281071550870317,
                    0.000204821658962907, 0.000182893525702204, -81.7495598043969, 0.994600000000000, 1, 0.00150000000000000, 0.0244000000000000, 0.149400000000000,
                    0.407100000000000, 0.416100000000000, 0, 0.000100000000000000, 0.000600000000000000, 0.000800000000000000, 0, 0, 0, 0, 0, 0, 0, 0.000149496595416145,
                    0.988960482230953, 0.00264830219853136, 0.173614838476469, -2164.01487897279, 0.802586369629388, 0.0286375258151476, 0.000416030212867675,
                    1.16241361165146e-06, 8.11698968573295e-07, 2.89626479052790e-08, 4.20742890084225e-10, 1.17510183512479e-12, 0.105921364476028, 0.00378635125822485,
                    .03834795404533e-05, 0.0585360319738570, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
                   };

    for (int i = 0; i < 101; i++) {
        y[i] = tmp[i];
        ydot[i] = 0;
    }

    // All these flag parameters are not changed in all the simulation
    ISO         = 0;
    flagMina    = 1;
    AF          = 0;
    RA          = 0;
    epi         = 1;
    bGal        = 1;
    EAD         = 0;
    drug_index  = 0;
    drug        = 0;
    Ach         = 0.0;

    // channel scaling factors
    gnabar  = 1;    // Na Current conductance
    gnabbar = 1;    // Na Background Current conductance
    gnakbar = 1;    // Na/K Pump Current conductance
    gkrbar  = 1;    // Rapidly Activating K Current conductance
    gksbar  = 1;    // Slowly Activating K Current conductance
    gkpbar  = 1;    // Plateau K Current conductance
    gkachbar = 1;   // Muscarinic-Receptor-Activated K Current conductance
    gtobar  = 1;    // Transient Outward K Current conductance
    gkurbar = 1;    // Ultra rapid delayed rectifier Outward K Current conductance
    gkibar  = 1;    // Time-Independent K Current conductance
    vupbar  = 1;    // SERCA strength
    gcabbar = 1;    // Ca Background conductance
    gpcabar = 1;    // Sarcolemmal Ca Pump Current conductance
    gncxbar = 1;    // Na/Ca Exchanger flux conductance
    gcabar  = 1;    // L-type Calcium Current conductance

    v = y[39];
    ci = y[38];
    cnsr = y[31];
    vold  = v;
}

CCell::~CCell() {
    delete[] y;
    delete[] ydot;
}

CCell& CCell::operator=(const CCell& cell) {
    if (&cell != this) {
        for (int i = 0; i < N; i++) {
            y[i] = cell.y[i];
        }

        vold    = cell.vold;
        dt      = cell.dt;
        gnabar  = cell.gnabar;  // Na Current conductance
        gnabbar = cell.gnabbar; // Na Background Current conductance
        gnakbar = cell.gnakbar; // Na/K Pump Current conductance
        gkrbar  = cell.gkrbar;  // Rapidly Activating K Current conductance
        gksbar  = cell.gksbar;  // Slowly Activating K Current conductance
        gkpbar  = cell.gkpbar;  // Plateau K Current conductance
        gkachbar = cell.gkachbar; // Muscarinic-Receptor-Activated K Current conductance
        gtobar  = cell.gtobar;  // Transient Outward K Current conductance
        gkurbar = cell.gkurbar; // Ultra rapid delayed rectifier Outward K Current conductance
        gkibar  = cell.gkibar;  // Time-Independent K Current conductance
        vupbar  = cell.vupbar;  // SERCA strength
        gcabbar = cell.gcabbar; // Ca Background conductance
        gpcabar = cell.gpcabar; // Sarcolemmal Ca Pump Current conductance
        gncxbar = cell.gncxbar; // Na/Ca Exchanger flux conductance
        gcabar  = cell.gcabar;  // L-type Calcium Current conductance
    }
    return (*this);
}

void CCell::pace(double st) {
    y[39] = v;
    ci = y[38];
    cnsr = y[31];
    double dv = (vold - v) / dt;

    vold      = v;


    ddt = dt;
    pacex(st);

    v  = y[39];
    return;
}

void CCell::pacex(double st) {
    // Nernst Potentials
    ek              = (1 / FoRT) * log(Ko / y[35]); // [mV]
    ecl             = (1 / FoRT) * log(Cli / Clo);  // [mV]

    double I_Na_junc        = I_Na_junc_3d;
    double I_Na_sl          = I_Na_sl_3d;

    double I_nabk_junc      = I_Nabk_junc_3d;
    double I_nabk_sl        = I_Nabk_sl_3d;

    double I_ncx_junc       = incx_junc_ave;
    double I_ncx_sl         = incx_sl_ave;

    double I_ncx            = I_ncx_junc + I_ncx_sl;

    double I_nak_junc       = I_NaK_junc_3d;
    double I_nak_sl         = I_NaK_sl_3d;
    double I_nak            = I_nak_junc + I_nak_sl;

    double I_CaNa_junc      = I_CaNa_junc_3d;
    double I_CaNa_sl        = 0.0;

    double I_sk             = I_SK_junc_3d + I_SK_sl_3d;

    double I_Na_tot_junc    = I_Na_junc + I_nabk_junc + 3 * I_ncx_junc + 3 * I_nak_junc + I_CaNa_junc;  // [uA/uF]
    double I_Na_tot_sl      = I_Na_sl + I_nabk_sl + 3 * I_ncx_sl + 3 * I_nak_sl + I_CaNa_sl;            // [uA/uF]

    // K Concentration
    double I_to     = comp_ito();
    double I_kr     = comp_ikr();
    double I_ks     = comp_iks();
    double I_ki     = comp_ik1();
    double I_kp     = comp_ikp();
    double I_kur    = comp_ikur();
    double I_KAch   = comp_ikach();
    double I_CaK    = I_CaK_junc_3d;

    double I_K_tot  = I_to + I_kr + I_ks + I_ki - 2 * I_nak + I_CaK + I_kp + I_kur + I_KAch + I_sk;     // [uA/uF]
    ydot[35] = 0; // [mM/msec] Clamped [K+]

    double I_cabk_junc      = 0;
    double I_cabk_sl        = I_cabk_ave;

    double I_pca_junc       = 0;
    double I_pca_sl         = I_pca_ave;

    double I_Ca_junc        = ica_ave;
    double I_Ca_sl          = 0;
    // #endif


    double I_Ca_tot_junc    = I_Ca_junc + I_cabk_junc + I_pca_junc - 2 * I_ncx_junc; // [uA/uF]
    double I_Ca_tot_sl      = I_Ca_sl + I_cabk_sl + I_pca_sl - 2 * I_ncx_sl;         // [uA/uF]

    double I_ClCa_junc    = I_ClCa_junc_3d;
    double I_ClCa_sl      = I_ClCa_sl_3d;
    double I_ClCa         = I_ClCa_junc + I_ClCa_sl;
    double I_Clbk         = comp_I_Clbk();

    //// Membrane Potential
    double I_Na_tot       = I_Na_tot_junc + I_Na_tot_sl;
    double I_Cl_tot       = I_ClCa + I_Clbk;
    double I_Ca_tot       = I_Ca_tot_junc + I_Ca_tot_sl;
    double I_tot          = I_Na_tot + I_Cl_tot + I_Ca_tot + I_K_tot; // [pA/pF]
    ydot[39]              = -(I_tot - st);  // voltage

    //update variables
    for (int i = 0; i < N; i++) {
        y[i]    += ydot[i] * ddt;
    }

    _I_to   = I_to;
    _I_kr   = I_kr;
    _I_ks   = I_ks;
    _I_ki   = I_ki;
    _I_kp   = I_kp;
    _I_kur  = I_kur;
    _I_KAch = I_KAch;
    _I_CaK  = I_CaK;
    _I_ncx  = I_ncx;

    _I_Clbk = I_Clbk;

    double I_Ca     = I_Ca_junc + I_Ca_sl;
    double I_CaNa   = I_CaNa_junc + I_CaNa_sl;
    double I_Catot  = I_Ca + I_CaK + I_CaNa;
    _I_Catot        = I_Catot;
}

// I_kr: Rapidly Activating K Current
double CCell::comp_ikr(void) {
    double factor_rano_kr;
    double IC50_kr;

    if ( (flagMina == 1) && (drug_index == 1) ) {
        IC50_kr         = 35 * (1e-6);
        factor_rano_kr  = 1 / (1 + (drug / IC50_kr) );
    } else {
        factor_rano_kr  = 1;
    }

    double gkr      = 0.035 * sqrt(Ko / 5.4) * factor_rano_kr;
    // gkr = 0.035*sqrt(Ko/5.4);
    double xrss     = 1 / (1 + exp(-(y[39] + 10) / 5) );
    double tauxr    = 550 / (1 + exp((-22 - y[39]) / 9) ) * 6 / (1 + exp((y[39] - (-11)) / 9) ) + 230 / (1 + exp((y[39] - (-40)) / 20) );

    ydot[12]        = (xrss - y[12]) / tauxr;

    double rkr      = 1 / (1 + exp((y[39] + 74) / 24) );
    double I_kr     = gkrbar * (gkr * y[12] * rkr * (y[39] - ek) );

    return I_kr;
}

// I_ks: Slowly Activating K Current
double CCell::comp_iks(void) {
    double eks      = (1 / FoRT) * log( (Ko + pNaK * Nao) / (y[35] + pNaK * y[34]) );
    double gks_junc = (1 + 1 * AF + 2 * ISO) * 0.0035 * 1;
    double gks_sl   = (1 + 1 * AF + 2 * ISO) * 0.0035 * 1;
    double  xsss    = 1 / (1 + exp(-(y[39] + 40 * ISO + 3.8) / 14.25) );
    double tauxs    = 990.1 / (1 + exp(-(y[39] + 40 * ISO + 2.436) / 14.12) );

    ydot[13]        = (xsss - y[13]) / tauxs;

    double I_ks_junc = gksbar * (Fjunc * gks_junc * pow((y[13]), 2) * (y[39] - eks) );
    double I_ks_sl  = gksbar * (Fsl * gks_sl * pow((y[13]), 2) * (y[39] - eks) );
    double I_ks     = I_ks_junc + I_ks_sl;

    return I_ks;
}

// I_kp: Plateau K Current
double CCell::comp_ikp(void) {
    double kp_kp    = 1 / (1 + exp(7.488 - y[39] / 5.98) );
    double I_kp_junc = gkpbar * (Fjunc * gkp * kp_kp * (y[39] - ek) );
    double I_kp_sl  = gkpbar * (Fsl * gkp * kp_kp * (y[39] - ek) );
    double I_kp     = I_kp_junc + I_kp_sl;

    return I_kp;
}

// I_k,ach: Muscarinic-Receptor-Activated K Current
double CCell::comp_ikach(void) {
    double I_KAch   = gkachbar * (1 / (1 + pow((0.03 / Ach), 2.1) ) * (0.08 + 0.4 / (1 + exp((y[39] + 91) / 12) ) ) * (y[39] - ek) );

    return I_KAch;
}

// I_to: Transient Outward K Current (slow and fast components)
double CCell::comp_ito(void) {
    double GtoFast  = (1.0 - 0.7 * AF) * 0.165 * 1.0; //nS/pF

    // equations for activation;
    double xtoss    = ( (1) / ( 1 + exp( -(y[39] + 1.0) / 11.0 ) ) );
    double tauxtof  = 3.5 * exp( -(pow((y[39] / 30.0), 2.0)) ) + 1.5;

    // equations for inactivation;
    double ytoss    = ( (1.0) / ( 1 + exp( (y[39] + 40.5) / 11.5) ) ) ;
    double tauytof  = 25.635 * exp( -(pow(((y[39] + 52.45) / 15.8827), 2.0)) ) + 24.14;

    ydot[10]        = (xtoss - y[10]) / tauxtof;
    ydot[11]        = (ytoss - y[11]) / tauytof;
    double I_tof    = gtobar * (1.0 * GtoFast * y[10] * y[11] * (y[39] - ek) );

    double I_to     = I_tof;

    return I_to;
}

// I_kur: Ultra rapid delayed rectifier Outward K Current
double CCell::comp_ikur(void) {
    double Gkur     = 0.8 * (1.0 - 0.5 * AF) * (1 + 2 * ISO) * 0.045 * (1 + 0.2 * RA); //nS/pF
    double xkurss   = ( (1) / ( 1 + exp( (y[39] + 6) / -8.6 ) ) );
    double tauxkur  = 9 / (1 + exp((y[39] + 5) / 12.0)) + 0.5;

    double ykurss   = ( (1) / ( 1 + exp( (y[39] + 7.5) / 10 ) ) );
    double tauykur  = 590 / (1 + exp((y[39] + 60) / 10.0) ) + 3050;

    ydot[58]        = (xkurss - y[58]) / tauxkur;
    ydot[59]        = (ykurss - y[59]) / tauykur;
    double I_kur    = gkurbar * (1 * Gkur * y[58] * y[59] * (y[39] - ek));

    return I_kur;
}

// I_ki: Time-Independent K Current
double CCell::comp_ik1(void) {
    double aki      = 1.02 / (1 + exp(0.2385 * (y[39] - ek - 59.215)));
    double bki      = (0.49124 * exp(0.08032 * (y[39] + 5.476 - ek)) + exp(0.06175 * (y[39] - ek - 594.31)) ) / (1 + exp(-0.5143 * (y[39] - ek + 4.753)));
    double kiss     = aki / (aki + bki);

    double I_ki     = gkibar * ((1 + 1 * AF) * 0.0525 * sqrt(Ko / 5.4) * kiss * (y[39] - ek) );

    return I_ki;
}

double CCell::comp_I_Clbk(void) {
    double I_Clbk   = GClB * (y[39] - ecl);
    return I_Clbk;
}

void CCell::get_steady_state() 
{
        char            filename105 [1000];
        sprintf         (filename105, "./steady_state_output/steady_state_Ik.txt" );   

        output_double_array_txt(filename105, y, 101);
}

void CCell::set_steady_state()
{
        char            filename109 [1000];
        sprintf         (filename109, "./steady_state_init/steady_state_Ik.txt" );   

        read_double_from_txt( filename109, y, 101 );

        v = y[39];
        set_vold(v);
}