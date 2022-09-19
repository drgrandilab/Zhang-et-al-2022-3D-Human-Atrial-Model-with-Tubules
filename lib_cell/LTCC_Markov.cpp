#include <cmath>
#include "LTCC_Markov.hpp"

LTCC_Markov::LTCC_Markov(int NCaL_set, int ID)
    : random_num_LTCC(1, ID)
{
    open_NCaL_num   = 0;
    NCaL = NCaL_set;

    // Starting from C2 state
    for (int i = 0; i < NCaL; ++i) {
        LTCC_State.push_back(2);
    }
}

LTCC_Markov::~LTCC_Markov() {

}

double LTCC_Markov::update_states_v7(double dt, double Caj, double Vm) {

    double cpt     = 3.75e-3 * 0.55;    // [mM]
    double tmp     = cpt / Caj;
    double fcp     = 0.5 / (1 + (tmp / 3.0) * (tmp / 3.0))  + 0.50 / (1 + (tmp) * (tmp) * (tmp) * (tmp));

    double R2      = aR2 / (1 + exp(-(Vm - sR2) * pR2));

    double T2      = 0.6 * 1.5 * (0.59 + (0.8 * exp(0.052 * (Vm + 13.0))) / (1 + exp(0.132 * (Vm + 13.0))));

    double PT      = 1 - (1 / (1 + exp(-(Vm + sp2) / sp3)));
    double R1      = aR1 / (1 + exp(-(Vm - sR1) * pR1));
    double T1      = aT1 / (1 + exp(-(Vm - sT1) * pT1));
    double Rv1     = aRv1 / (1 + exp(-(Vm - sRv1) * pRv1));

    double Tv1_K5  = 20 + 30 * exp(-(Vm + 50) * (Vm + 50) / 150) +  30.0 / (1 + exp((Vm + 25) / 5));

    double alphaLCC     = R2 / T2;
    double betaLCC      = (1 - R2) / T2;

    double r1           = R1 / T1;
    double r2           = (1 - R1) / T1;

    double Vm_test      = 32.0;

    double SLOPE        = 8;
    double tau_k1k2     = 1.0 / (k1o * fcp);
    double Rv1_k1k2     = 0.85 / (1 + exp((Vm + Vm_test - 14 + 3 * fcp) / SLOPE)) + 0.15;

    double k1     = (1 - Rv1_k1k2) / tau_k1k2;
    double k2     = Rv1_k1k2  / tau_k1k2 ;

    double k3     = PT / sp1;

    double fcp_k5 =  fcp;
    double Is_Vss = (1.0 / (1 + exp( (Vm  + Vm_test + -14 + 3 * fcp) / SLOPE )));
    double k5     =  Is_Vss / Tv1_K5 / (fcp_k5 * 1.0 + 1.0);
    double k6     = (1 - Is_Vss) / Tv1_K5 * (fcp_k5 * 1.0 + 1.0) * 1;

    double k4     = k3 * (alphaLCC / betaLCC) * (k1 / k2) * (k5 / k6);

    double Tv1    =  45 + 3000 * exp(-(Vm + 50) * (Vm + 50) / 150) + 200 * 1.0 / (1 + exp((Vm + 60) / -3));

    double k1p    = Rv1 / Tv1;
    double k2p    = (1 - Rv1) / Tv1;
    double k3p    = k3;

    double Is_Vss_V        = (1 / (1 + exp( (Vm + Vm_test) / SLOPE )) + 0);

    double Tv1_K5p         =  45 + 2000 * exp(-(Vm + 50) * (Vm + 50) / 150) + 200 * 1.0 / (1 + exp((Vm + 60) / -3));
    double k5p    = Is_Vss_V / Tv1_K5p ;
    double k6p    =  (1 - Is_Vss_V) / Tv1_K5p;

    double k4p    = k3p * k5p * alphaLCC * k1p / (k2p * betaLCC * k6p);
    double s1p    = k1p;

    double s2p    = s1p * k2p * r1 / (r2 * k1p);

    double s1     = k1;
    double s2     = 0.2 * k2;

    double k13    = r2;
    double k14    =  k13 * k2 * r1 * s1 / ( s2 * r2 * k1);

    Rv1_k1k2      = 1 / (1 + exp((Vm + Vm_test - 14 + 8 * fcp) / SLOPE));
    double k7     = (1 - Rv1_k1k2) / 50;
    double k8     = Rv1_k1k2 / 50;

    double Tv1_K9 = 60 + 30 * exp(-(Vm + 55) * (Vm + 55) / 150);

    Is_Vss        = Rv1_k1k2;

    double k9     = 6 * Is_Vss / Tv1_K9;
    double k10    = 6 * (1 - Is_Vss) / Tv1_K9;

    double k12    = k3;
    double k11    = k12 * (k9 * k4 * k7) / (k8 * k3 * k10 );
    double fcp_k15         = 0.5 / (1 + (cpt / Caj / 2.5) * (cpt / Caj / 2.5) )  + 0.5 / (1 + (cpt / Caj) * (cpt / Caj) * (cpt / Caj) * (cpt / Caj));
    double k15_factor      = 1;
    double tau_k15 = 1 / (0.5 * fcp_k15 + 0.01);
    double k15     = k15_factor * (1 - fcp_k15) / tau_k15;
    double k15p    = k15_factor * fcp_k15 / tau_k15;
    double tau_k16 = 1 / (0.5 * fcp_k15 + 0.01);

    double k16     = k15_factor * (1 - fcp_k15) / tau_k16;
    double k16p    = k15_factor * fcp_k15 / tau_k16;


    //          ---- I4 -- I3 --------
    //          |    |     |         |
    //          |    I2 -- I1 -- I0  |
    //          |    |     |     |   |
    //          |    C2 -- C1 -- O   |
    //          |    |     |     |   |
    //          ---- I2'-- I1'--------

    //          ---- 7  -- 6  --------
    //          |    |     |         |
    //          |    3  -- 2  -- 8   |
    //          |    |     |     |   |
    //          |    0  -- 1  -- 9   |
    //          |    |     |     |   |
    //          ---- 5  -- 4  --------

    int open_NCaL_num = 0;

    for (int j = 0; j < NCaL; j++)
    {
        double rand_value = random_num_LTCC.gen_rand();
        double rr = rand_value / dt;

        switch (LTCC_State[j])
        {
        case 0:
            if (rr < alphaLCC)
                LTCC_State[j]   =   1;
            else if (rr < (alphaLCC + k6) )
                LTCC_State[j]   =   3;
            else if (rr < (alphaLCC + k6 + k6p) )
                LTCC_State[j]   =   5;
            break;

        case 1:
            if (rr < betaLCC)
                LTCC_State[j]   =   0;
            else if (rr < (betaLCC + k1) )
                LTCC_State[j]   =   2;
            else if (rr < (betaLCC + k1 + k1p) )
                LTCC_State[j]   =   4;
            else if (rr < (betaLCC + k1 + k1p + r1) )
                LTCC_State[j]   =   9;
            break;

        case 2:
            if (rr < k2)
                LTCC_State[j]   =   1;
            else if (rr < (k2 + k3) )
                LTCC_State[j]   =   3;
            else if (rr < (k2 + k3 + k7) )
                LTCC_State[j]   =   6;
            else if (rr < (k2 + k3 + k7 + k14) )
                LTCC_State[j]   =   8;
            break;

        case 3:
            if (rr < k5)
                LTCC_State[j]   =   0;
            else if (rr < (k5 + k4) )
                LTCC_State[j]   =   2;
            else if (rr < (k5 + k4 + k10) )
                LTCC_State[j]   =   7;
            break;

        case 4:
            if (rr < k2p)
                LTCC_State[j]   =   1;
            else if (rr < (k2p + k3p) )
                LTCC_State[j]   =   5;
            else if (rr < (k2p + k3p + k16p) )
                LTCC_State[j]   =   6;
            else if (rr < (k2p + k3p + k16p + s2p) )
                LTCC_State[j]   =   9;
            break;

        case 5:
            if (rr < k5p)
                LTCC_State[j]   =   0;
            else if (rr < (k5p + k4p) )
                LTCC_State[j]   =   4;
            else if (rr < (k5p + k4p + k15p) )
                LTCC_State[j]   =   7;
            break;

        case 6:
            if (rr < k8)
                LTCC_State[j]   =   2;
            else if (rr < (k8 + k16) )
                LTCC_State[j]   =   4;
            else if (rr < (k8 + k16 + k12) )
                LTCC_State[j]   =   7;
            break;

        case 7:
            if (rr < k9)
                LTCC_State[j]   =   3;
            else if (rr < (k9 + k15) )
                LTCC_State[j]   =   5;
            else if (rr < (k9 + k15 + k11) )
                LTCC_State[j]   =   6;
            break;

        case 8:
            if (rr < k13)
                LTCC_State[j]   =   2;
            else if (rr < (k13 + s2) )
                LTCC_State[j]   =   9;
            break;

        case 9:
            open_NCaL_num ++;
            if (rr < r2)
                LTCC_State[j]   =   1;
            else if (rr < (r2 + s1p) )
                LTCC_State[j]   =   4;
            else if (rr < (r2 + s1p + s1) )
                LTCC_State[j]   =   8;
            break;
        }
    }

    return open_NCaL_num;
}