#include "subcell.hpp"
#include <stdio.h>
#include <string.h>

#ifdef ___DEBUG
#include <iomanip>
#endif

inline unsigned int xorshift(unsigned int *xx, unsigned int *yy, unsigned int *zz, unsigned int *ww)
{
        unsigned int t = (*xx ^ (*xx << 11)); *xx = *yy; *yy = *zz; *zz = *ww;
        return ( *ww = (*ww ^ (*ww >> 19)) ^ (t ^ (t >> 8)) );
}

CSubcell::CSubcell(int sizex, int sizey, int sizez, int fmesh, double xii, int _seed)
{
        if (fmesh == 1)
                dt      = 0.1;
        else if (fmesh == 5)
                dt      = 0.005;
        else if (fmesh == 3)
                dt      = 0.01;
        else {
                cerr << "fine mesh incorrect !\n";
                exit(1);
        }

        nx      = sizex;
        ny      = sizey;
        nz      = sizez;

        num_tubule      = nx * ny * nz;
        finemesh        = fmesh;
        xi              = xii;
        
        n = nn = 0;
        
        layer           = 2;

        cao     = 1.8;                  // [mM]
        vup     = 1.2 * 0.3;            // [uM/ms]
        
        #ifdef ___USE_ORG_PARAM
                Jmax    = 0.0147;
                Ku      = 3.8 * 0.0001;
                Kb      = 5 * 0.00001;
                tauu    = 125.00;
                taub    = 5.0;
                tauc1   = 1.0;
                tauc2   = 1.0;
                hh      = 23;
                KK      = 850;
                gcabk   = 0;
                qslcap  = 0;
        #else
                Jmax    = 0.0147 * 18;
                Ku      = 5.0;
                Kb      = 0.005;
                tauu    = 1250.0;
                taub    = 0.5;
                tauc1   = 2.0;
                tauc2   = 0.3;
                gcabk   = 0.5 * 0.002513;
                qslcap  = 2.35 * 4;
        #endif


        kk_pow_10 = pow((KK / 1000.0), (hh));

        MaxSR   = 15;
        MinSR   = 1;
        ec50SR  = 450;
        hkosrca = 2.5;

        kup     = 2.5 * 0.123 * 2.0;
        KNSR    = 1700;

        NoCaL   = 4;
        gleak   = 1.035 * 0.00001;

        taup    = 0.022;
        tautr   = 5.0;
        BCSQN0  = 400;

        Kcp     = 100;
        pedk12  = 0.000001;
        pedk43  = 0.000001;

        seed    = _seed;
        bc      = 0;
        initialized     = false;

        #ifdef ___NCX
                NCXalpha        = 0.11;
                Fjunc           = 0.11;
        #endif

        #ifdef ___DEBUG
                bSTOP   = false;
        #endif

        #ifndef __UNIFORM_TUBULE
                hp_ryr  = 1.0;
        #endif

        num_open_ICaL   = 0;

        irave           = 0;    // sumir;
        iupave          = 0;    // sumjup;
        icaave          = 0;    // sumica;
        ileakave        = 0;    // sumjleak;
        icabkave        = 0;    // sumjcabk;
        islcapave       = 0;    // sumjslcap;

        incx_sl_ave     = 0;    // sum_jnaca_j_flux;
        incx_junc_ave   = 0;    // sum_jnaca_sl * vs;
        incxave         = 0;    // incx_junc_ave + incx_sl_ave;


        i_clca_junc_ave = 0;    // sum_iclca_j;
        i_clca_sl_ave   = 0;    // sum_iclca_sl;
        i_cana_ave      = 0;    // sumicana;
        i_cak_ave       = 0;    // sumicak;

        i_nak_junc_ave  = 0;    // sum_inak_junc;
        i_nak_sl_ave    = 0;    // sum_inak_sl;
        i_nabk_junc_ave = 0;    // sum_inabk_junc;
        i_nabk_sl_ave   = 0;    // sum_inabk_sl;
        i_na_junc_ave   = 0;    // sum_ina_junc;
        i_na_sl_ave     = 0;    // sum_ina_sl;

        i_sk_junc_ave   = 0;    // sum_isk_junc;
        i_sk_sl_ave     = 0;    // sum_isk_sl;
        
        double par_sk_tmp[] = {0, 
                0.0506381114404388,0.273335569451572,2.96381060498817,0.199981221802789,0.279328126521496,-86.9289059836381,0.00636311816933264,5.22915055145375};
        par_sk  = new double [9];
        int len_par_sk      = sizeof(par_sk_tmp)/sizeof(*par_sk_tmp);
        for(int i = 0; i < (len_par_sk); i++) 
        {
                par_sk[i]   = par_sk_tmp[i];
        }        
        
        gsk     = 0.8 * 1.0 / 8113.0 * 0.459661 * par_sk[1] * 0.96245 * (1e-10);        // [A]
        kdsk    = 0.794288 * ( pow(10, (-3.45)) );                                      // (mM)
 
        steady_state_CaNa       = nullptr;
        steady_state_currents   = nullptr;
        steady_state_buffers    = nullptr;
        
        gca     = 1.475;;
        pNa     = 1.0 / 8113.0 * 0.75e-8 * 0.96245 * (1e-10); // [cm/sec*F]
        pK      = 1.0 / 8113.0 * 1.35e-7 * 0.96245 * (1e-10); // [cm/sec*F]
}


void CSubcell::init(double initci, double initcj, double initnai)
{
        nxny            = nx * ny;
        n               = nx * ny * nz;

        finemesh3       = finemesh * finemesh * finemesh;
        nnx             = finemesh * nx;
        nny             = finemesh * ny;
        nnz             = finemesh * nz;
        nnxnny          = nnx * nny;
        nn              = nnx * nny * nnz;

        //cell parameters
        vi      = 0.969 / finemesh3;
        vs      = 0.025 / finemesh3;
        vnsr    = 0.025 / finemesh3;
     
        tausT   = 0.9 * xi * 1.42 / (finemesh * finemesh);
        tausL   = 0.9 * xi * 3.4 / (finemesh * finemesh);
        taumninv        = 2 / tausL + 4 / tausT;

        tauiT   = 0.9 * xi * 0.81 / 0.68 / (finemesh * finemesh);
        tauiL   = 0.9 * xi * 3.39 / 0.68 * 0.39 / (finemesh * finemesh);
        taunsrT         = xi * 7.2 / (finemesh * finemesh);
        taunsrL         = xi * 24.0 / (finemesh * finemesh);

        tausi = 0.17 / (finemesh * finemesh);      

        ci      = new double [nn];
        cs      = new double [nn];
        cp      = new double [n];
        cjsr    = new double [n];
        cnsr    = new double [nn];

        nai     = new double [nn];
        nas     = new double [nn];
        nap     = new double [n];

        // Buffers components
        cats    = new double [nn];

        Cyto_MMg        =   new double [nn];
        Cyto_MCa        =   new double [nn];
        Cyto_TnCHc      =   new double [nn];
        Cyto_TnCHm      =   new double [nn];
        Cyto_SR         =   new double [nn];
        Cyto_cam        =   new double [nn];
        cati    = new double [nn];
        
        cleft_cam       = new double [n];
        cleft_SR        = new double [n];
        cleft_SLL       = new double [n];
        cleft_SLH       = new double [n];
        
        SL_cam          = new double [nn];
        SL_SR           = new double [nn];
        SL_SLL          = new double [nn];
        SL_SLH          = new double [nn];
        
        JSR_CSQN        = new double [n];
        
        cleft_NaB       = new double [n];

        SL_NaB          = new double [nn];

        #ifdef ___DETERMINISTIC
                c1      = new double [n];
                c2      = new double [n];
                i1ca    = new double [n];
                i1ba    = new double [n];
                i2ca    = new double [n];
                i2ba    = new double [n];

                fryr1   = new double [n];
                fryr2   = new double [n];
                fryr3   = new double [n];
        #else
                y       = new int [n * NoCaL];
                ryr1    = new int [n];
                ryr2    = new int [n];
                ryr3    = new int [n];
                nryr    = new int [n];
        #endif

        // Ina.ina();
        ICaNa_array       = new double [n];
        ICaK_array        = new double [n];
        
        I_ClCa_junc_array = new double [n];
        I_nak_junc_array  = new double [n];
        I_nabk_junc_array = new double [n];
        I_Na_junc_array   = new double [n];


        // Xianwei: output currents
        JSR_buffer_factors = new double [n];
        j_serca         = new double [nn];
        j_ryr           = new double [n];
        ica_array       = new double [n];

        memset(ica_array, 0, n*sizeof(double));

        LTCC_Open_num = new int[n];
        ncx_array       = new double [nn];
        icabk_array     = new double [nn];
        jpca_array      = new double [nn];
        tubule_flag     = new int [n];
        I_ClCa_sl_array = new double[nn];
        I_nak_sl_array = new double[nn];
        I_nabk_sl_array = new double[nn];
        I_Na_sl_array = new double[nn];

        I_SK_junc_array         = new double[nn];
        I_SK_sl_array           = new double[nn];

        vp      = new double [n];

        PoSpark = new double [n];

        Jmaxx   = new double [n];
        cscp1   = new double [nn];
        cscp2   = new double [nn];

        Ici     = new double [nn];
        Icnsr   = new double [nn];

        Inai    = new double [nn];

        #ifndef __RYR_UNIFORM
                sigma   = new double [n];
        #endif

        #ifdef ___NO_CS_BUFFER
                csmn    = new double [nn];
        #else
                Ics     = new double [nn];
                Idps    = new double [nn];

                Inas    = new double [nn];
        #endif

        #ifdef ___EGTA
                caEGTAi = new double [nn];
                caEGTAs = new double [nn];
        #endif

        // Xianwei: cleft NCX flux
        #ifdef ___NCX
                ncx_array_p     = new double [n];
        #endif

        map_serca       = new double [nn];
        map_nryr        = new int [n];
        map_nlcc        = new int [n];
        map_ncx         = new double [nn]; 

        Itr     = new double [nn];

        crupos  = new int [n];

        //random number
        xsx     = new unsigned int [n];
        xsy     = new unsigned int [n];
        xsz     = new unsigned int [n];
        xsw     = new unsigned int [n];
        #pragma omp parallel for
        #pragma ivdep
        #pragma vector always
        for (int id = 0; id < n; id++)
        {
                xsx[id]         = 123456789 + id + seed;
                xsy[id]         = 362436069 + id * 100 + seed * 10;
                xsz[id]         = 521288629 + id * 1000 + seed * 100;
                xsw[id]         = 88675123 + id * 10000 + seed * 1000;
                
                for (int i = 0; i < 1000; i++)          xorshift(&xsx[id], &xsy[id], &xsz[id], &xsw[id]);
        }

        //initial conditions
        #pragma omp parallel for
        #pragma ivdep
        #pragma vector always
        for (int id = 0; id < nn; id++)
        {
                ci[id]          = initci;
                cs[id]          = initci;
                cnsr[id]        = initcj;
                Icnsr[id]       = 0;

                Inai[id]        = 0;

                #ifdef ___NO_CS_BUFFER
                        csmn[id]        = 0;
                #else
                        Ics[id]         = 0;
                        Idps[id]        = 0;

                        Inas[id]        = 0;
                #endif
                
                j_serca[id]     = 0;

                Ici[id]         = 0;

                cscp1[id]       = 0;
                cscp2[id]       = 0;
                Itr[id]         = 0;

                nai[id]         = initnai;
                nas[id]         = initnai;

                Inas[id]        = 0;

                ncx_array_p[id] = 0;

                ICaNa_array[id]       =0;
                ICaK_array[id]        =0;
                I_ClCa_junc_array[id] =0;
                I_nak_junc_array[id]  =0;
                I_nabk_junc_array[id] =0;
                I_Na_junc_array[id]   =0;
                ncx_array[id]         =0;
                icabk_array[id]       =0;
                jpca_array[id]        =0;
                I_ClCa_sl_array[id] = 0;
                I_nak_sl_array[id]=0;
                I_nabk_sl_array[id]=0;
                I_Na_sl_array[id]=0;

                I_SK_junc_array[id]     = 0;
                I_SK_sl_array[id]       = 0;
        }

        resetBuffer();
        
        #pragma omp parallel for
        #pragma ivdep
        #pragma vector always
        for (int id = 0; id < n; id++)
        {
                ICaNa_array [id] = 0;
                I_ClCa_junc_array[id]= 0;
                ICaK_array [id] = 0;
                I_nak_junc_array [id] = 0;
                I_nabk_junc_array [id] = 0;
                I_Na_junc_array [id] = 0;
                Jmaxx[id]       = Jmax;

                PoSpark[id]     = 0.0;

                cp[id]          = initci;
                cjsr[id]        = initcj;

                nap[id]         = initnai;
                
                #ifdef ___DETERMINISTIC
                        c1[id]          = 1;
                        c2[id]          = 0;
                        i1ca[id]        = 0;
                        i1ba[id]        = 0;
                        i2ca[id]        = 0;
                        i2ba[id]        = 0;

                        fryr1[id]       = 0.03;
                        fryr2[id]       = 0;
                        fryr3[id]       = 0;

                #else
                        ryr1[id]        = 0 + int(5.0 * xorshift(&xsx[id], &xsy[id], &xsz[id], &xsw[id]) / (double)(UINT_MAX));
                        ryr2[id]        = 0;
                        ryr3[id]        = 0;
                        nryr[id]        = 100;
                
                        for (int j = 0; j < NoCaL; j++)
                                y[id * NoCaL + j] = 2;
                #endif

                j_ryr[id] = 0;

                JSR_buffer_factors[id] = 0;

                #ifdef ___UNIFORM
                        double r = 0;
                #else
                        double r = calcvp(0, 0.3, -0.8, 0.8, id); //Gaussian distribution (0,0.3) range(-0.8~0.8)
                #endif

                vp[id] = vp_ave * (1 + r);

                #ifndef __RYR_UNIFORM
                        sigma[id] = 1 + calcvp(0, 0.2, -0.8, 0.8, id);
                #endif

                tubule_flag[id] = 1;

                int kk = id / nxny;
                kk = kk * finemesh + finemesh / 2;

                int modi = id % nxny;
                
                int jj = modi / nx;
                jj = jj * finemesh + finemesh / 2;
                
                int ii = modi % nx;
                ii = ii * finemesh + finemesh / 2;
                
                crupos[id] = ii + jj * nnx + kk * nnxnny;
                cscp2[crupos[id]] = 1 / taup * vp[id] / vs;
        }

        #pragma omp parallel for
        #pragma ivdep
        #pragma vector always
        for (int id = 0; id < n; id++)
        {
                map_nryr[id]     = 41;
                map_nlcc[id]     = 4;
        }

        #pragma omp parallel for
        #pragma ivdep
        #pragma vector always
        for (int id = 0; id < nn; id++)
        {
                map_serca[id]   = 1.0;
                map_ncx[id]     = 1.0;
        }

        for (int i = 0; i < n; ++i)
        {
                LTCC_Markov_vec.push_back( LTCC_Markov(map_nlcc[i], i+1000*seed) );
                random_num.push_back(xor_rand(seed, i+1000*seed));
                RyR_vec.push_back(RyR(map_nryr[i], i+1000*seed));
        }

        iupave = icaave = incxave = irave = ileakave = icabkave = islcapave = 0;
        i_clca_junc_ave = i_clca_sl_ave = i_cana_ave = i_cak_ave = 0;
        i_nak_junc_ave = i_nak_sl_ave = i_nabk_junc_ave = i_nabk_sl_ave = i_na_junc_ave = i_na_sl_ave = 0;
        i_sk_junc_ave = i_sk_sl_ave = 0;

        NCX_fit_scale   = 1.0;
        INafast_gating  = 0.0;

        initialized = true;
}

/* initialise random generator with a different sed */
void CSubcell::srand(int sed)
{
        seed = sed;

        #pragma ivdep
        #pragma vector always
        for (int id = 0; id < n; id++)
        {
                xsx[id] = 123456789 + id + seed;
                xsy[id] = 362436069 + id * 100 + seed * 10;
                xsz[id] = 521288629 + id * 1000 + seed * 100;
                xsw[id] = 88675123 + id * 10000 + seed * 1000;
        }
        
        for (int id = 0; id < n; id++)
                for (int i = 0; i < 1000; i++)
                        xorshift(&xsx[id], &xsy[id], &xsz[id], &xsw[id]);
}

CSubcell::~CSubcell()
{
        if (initialized)
        {
                delarray();
        }
}
void CSubcell::delarray(void)
{
        delete [] ci;
        delete [] cs;
        delete [] cp;
        delete [] cjsr;
        delete [] cnsr;

        delete [] nai;
        delete [] nas;
        delete [] nap;

        #ifdef ___DETERMINISTIC
                delete [] c1;
                delete [] c2;
                delete [] i1ca;
                delete [] i1ba;
                delete [] i2ca;
                delete [] i2ba;

                delete [] fryr1;
                delete [] fryr2;
                delete [] fryr3;
        #else
                delete [] ryr1;
                delete [] ryr2;
                delete [] ryr3;
                delete [] nryr;
                delete [] y;
        #endif

        delete [] j_serca;
        delete [] j_ryr;
        delete [] ica_array;
        delete [] LTCC_Open_num;
        delete [] ncx_array;
        delete [] icabk_array;
        delete [] jpca_array;
        delete [] tubule_flag;


        delete [] I_ClCa_junc_array;
        delete [] I_nak_junc_array ;
        delete [] I_nabk_junc_array;
        delete [] I_Na_junc_array  ;

        delete [] ICaNa_array;
        delete [] ICaK_array;
        #ifdef ___NCX
                delete [] ncx_array_p;
        #endif

        delete [] cati;

        delete [] cats;

        delete [] Cyto_MMg;
        delete [] Cyto_MCa;
        delete [] Cyto_TnCHc;
        delete [] Cyto_TnCHm;
        delete [] Cyto_SR;
        delete [] Cyto_cam;

        delete [] cleft_cam;
        delete [] cleft_SR;
        delete [] cleft_SLL;
        delete [] cleft_SLH;
        
        delete [] SL_cam;
        delete [] SL_SR;
        delete [] SL_SLL;
        delete [] SL_SLH;
        
        delete [] JSR_CSQN;

        delete [] SL_NaB;

        delete [] cleft_NaB;

        delete [] Jmaxx;

        delete [] PoSpark;

        #ifndef __RYR_UNIFORM
                delete[] sigma;
        #endif

        delete [] vp;
        delete [] cscp1;
        delete [] cscp2;
        delete [] Itr;

        delete [] Ici;
        delete [] Icnsr;

        delete [] Inai;

        #ifdef ___NO_CS_BUFFER
                delete [] csmn;
        #else

        delete [] Ics;
        delete [] Idps;

        delete [] Inas;
        #endif
        delete [] crupos;

        delete [] xsx;
        delete [] xsy;
        delete [] xsz;
        delete [] xsw;
        delete [] Tubular_map;

        delete [] par_sk;

        delete [] I_SK_junc_array;
        delete [] I_SK_sl_array;
}

CSubcell& CSubcell::operator=(const CSubcell& sc)
{
        if (&sc == this)        return (*this);
        
        if (initialized)
        {
                delarray();
        }

        //constructor
        dt = sc.dt;
        nx = sc.nx;
        ny = sc.ny;
        nz = sc.nz;
        finemesh = sc.finemesh;
        xi = sc.xi;

        layer = sc.layer;

        cao = sc.cao;
        vup = sc.vup;
        kup = sc.kup;
        KNSR = sc.KNSR;
        Jmax = sc.Jmax;
        gca = sc.gca;

        gcabk = sc.gcabk;
        qslcap = sc.qslcap;
        gleak = sc.gleak;

        Ku = sc.Ku;
        Kb = sc.Kb;
        tauu = sc.tauu;
        taub = sc.taub;
        tauc1 = sc.tauc1;
        tauc2 = sc.tauc2;
        BCSQN0 = sc.BCSQN0;

        //#ifdef ___SIGMOID
        Kcp = sc.Kcp;
        pedk12 = sc.pedk12;
        pedk43 = sc.pedk43;
        //#endif
        NoCaL = sc.NoCaL;
        #ifdef ___NCX
                NCXalpha        = sc.NCXalpha;
                Fjunc         = sc.Fjunc;
        #endif

        MaxSR = sc.MaxSR;
        MinSR = sc.MinSR;
        ec50SR = sc.ec50SR;
        hkosrca = sc.hkosrca;

        NCX_fit_scale = sc.NCX_fit_scale;

        init();

        //initial conditions
        #pragma ivdep
        #pragma vector always
        for (int id = 0; id < nn; id++)
        {
                ci[id] = sc.ci[id];
                cs[id] = sc.cs[id];
                cati[id] = sc.cati[id];

                cats[id]=sc.cats[id];

                Cyto_MMg[id]     =   sc.Cyto_MMg[id];
                Cyto_MCa[id]     =   sc.Cyto_MCa[id];
                Cyto_TnCHc[id]   =   sc.Cyto_TnCHc[id];
                Cyto_TnCHm[id]   =   sc.Cyto_TnCHm[id];
                Cyto_SR[id]      =   sc.Cyto_SR[id];
                Cyto_cam[id]     =   sc.Cyto_cam[id];

                SL_cam[id]       = sc.SL_cam[id];
                SL_SR[id]        = sc.SL_SR[id];
                SL_SLL[id]       = sc.SL_SLL[id];
                SL_SLH[id]       = sc.SL_SLH[id];

                cnsr[id] = sc.cnsr[id];
                Icnsr[id] = sc.Icnsr[id];

                #ifdef ___NO_CS_BUFFER
                        csmn[id] = sc.csmn[id];
                #else
                        Ics[id] = sc.Ics[id];
                        Idps[id] = sc.Idps[id];
                #endif

                j_serca[id] = sc.j_serca[id];
                ncx_array[id] = sc.ncx_array[id];
                icabk_array[id] = sc.icabk_array[id];
                jpca_array[id] = sc.jpca_array[id];

                cscp1[id] = sc.cscp1[id];
                cscp2[id] = sc.cscp2[id];
                Itr[id] = sc.Itr[id];
        }
        
        #pragma ivdep
        #pragma vector always
        for (int id = 0; id < n; id++)
        {
                Jmaxx[id] = sc.Jmaxx[id];

                PoSpark[id] = sc.PoSpark[id];

                #ifndef __RYR_UNIFORM
                        sigma[id] = sc.sigma[id];
                #endif

                cp[id] = sc.cp[id];
                cjsr[id] = sc.cjsr[id];


                cleft_cam[id]   = sc.cleft_cam[id];
                cleft_SR[id]    = sc.cleft_SR[id];
                cleft_SLL[id]   = sc.cleft_SLL[id];
                cleft_SLH[id]   = sc.cleft_SLH[id];
                
                JSR_CSQN[id]    = sc.JSR_CSQN[id];

                #ifdef ___DETERMINISTIC
                        c1[id] = sc.c1[id];
                        c2[id] = sc.c2[id];
                        i1ca[id] = sc.i1ca[id];
                        i1ba[id] = sc.i1ba[id];
                        i2ca[id] = sc.i2ca[id];
                        i2ba[id] = sc.i2ba[id];

                        fryr1[id] = sc.fryr1[id];
                        fryr2[id] = sc.fryr2[id];
                        fryr3[id] = sc.fryr3[id];
                #else
                        ryr1[id] = sc.ryr1[id];
                        ryr2[id] = sc.ryr2[id];
                        ryr3[id] = sc.ryr3[id];
                        nryr[id] = sc.nryr[id];
                        
                        for (int j = 0; j < NoCaL; j++)
                                y[id * NoCaL + j] = sc.y[id * NoCaL + j];
                #endif

                j_ryr[id] = sc.j_ryr[id];
                ica_array[id] = sc.ica_array[id];
                LTCC_Open_num[id] = sc.LTCC_Open_num[id];

                #ifndef __UNIFORM_TUBULE
                        tubule_flag[id] = sc.tubule_flag[id];
                        hp_ryr = sc.hp_ryr;
                #endif

                #ifdef ___NCX
                        ncx_array_p[id] = sc.ncx_array_p[n];
                #endif

                vp[id] = sc.vp[id];
                cscp2[crupos[id]] = sc.cscp2[crupos[id]];
                Ici[id] = sc.Ici[id];
        }

        #pragma ivdep
        #pragma vector always
        for (int id = 0; id < n; id++)
        {
                xsx[id] = sc.xsx[id];
                xsy[id] = sc.xsy[id];
                xsz[id] = sc.xsz[id];
                xsw[id] = sc.xsw[id];
        }

        return (*this);
}

// gaussian distribution generator???
double CSubcell::calcvp(double mean, double std, double lim1, double lim2, int ii)
{
        double res;

        do {
                //Gaussian
                double x1, x2, w;
                do
                {
                        x1 = 2.0 * xorshift(&xsx[ii], &xsy[ii], &xsz[ii], &xsw[ii]) / (double)(UINT_MAX) - 1.0;
                        x2 = 2.0 * xorshift(&xsx[ii], &xsy[ii], &xsz[ii], &xsw[ii]) / (double)(UINT_MAX) - 1.0;
                        w = x1 * x1 + x2 * x2;
                } while (w >= 1.0);

                w       = sqrt((-2.0 * log(w)) / w );
                double y1 = x1 * w;
                res     = y1 * std + mean;

        } while (res < lim1 || res > lim2);

        return res;
}


double CSubcell::compute_avg_ci(void)
{
        double sum = 0;

        #pragma omp parallel for schedule(static) reduction(+:sum)
        for (int id = 0; id < nn; id++)
        sum += ci[id];
        return (sum / nn);
}
double CSubcell::compute_avg_cs(void)
{
        double sum = 0;
        #pragma omp parallel for schedule(static) reduction(+:sum)
        for (int id = 0; id < nn; id++)
        sum += cs[id];
        return (sum / nn);
}
double CSubcell::compute_avg_cnsr(void)
{
        double sum = 0;
        #pragma omp parallel for schedule(static) reduction(+:sum)
        for (int id = 0; id < nn; id++)
        sum += cnsr[id];
        return (sum / nn);
}

double CSubcell::compute_avg_cp(void)
{
        double sum = 0;
        #pragma omp parallel for schedule(static) reduction(+:sum)
        for (int id = 0; id < nn; id++)
        sum += cp[id];
        return (sum / n);
}

double CSubcell::compute_avg_cjsr(void)
{
        double sum = 0;
        #pragma omp parallel for schedule(static) reduction(+:sum)
        for (int id = 0; id < nn; id++)
        sum += cjsr[id];
        return (sum / n);
}

double  CSubcell::compute_ica(void) 
{
        double ica_stan         = 0.001 * icaave * (1e-15) * 2 * F * 1000 / Cmem;               // [pA/pF]
        return ica_stan;
}

double  CSubcell::compute_incx(void) 
{
        
        double incx_stan        = incxave * 0.001 * (1e-15) * F * 1000 / Cmem;                  // [pA/pF]
        return incx_stan;
}
double CSubcell::compute_incx_junc(void)
{
        double incx_stan_junc   = incx_junc_ave * 0.001 * (1e-15) * F * 1000 / Cmem;            // [pA/pF]
        return incx_stan_junc;
}
double CSubcell::compute_incx_sl(void)
{
        double incx_stan_sl     = incx_sl_ave * 0.001 * (1e-15) * F * 1000 / Cmem;              // [pA/pF]
        return incx_stan_sl;
}
double  CSubcell::compute_icabk(void) 
{
        double icabk_stan       = icabkave * 0.001 * vs * (1e-15) * 2 * F * 1000 / Cmem;        // [pA/pF]
        return icabk_stan;
}
double  CSubcell::compute_ipca(void) 
{
        double ipca_stan        = islcapave * 0.001 * vs * (1e-15) * 2 * F * 1000 / Cmem;       // [pA/pF]
        return ipca_stan;
}

double  CSubcell::compute_serca(void)
{
        double iup_stan         = iupave * 0.001 * vi * (1e-15) * 2 * F * 1000 / Cmem;          // [pA/pF]
        return iup_stan;
}

double  CSubcell::compute_ir(void)
{
        double ir_stan          = irave * 0.001 * vp_ave * (1e-15) * 2 * F * 1000 / Cmem;       // [pA/pF]
        return ir_stan;
}

double  CSubcell::compute_ileak(void)
{
        double ileak_stan       = ileakave * 0.001 * vi * (1e-15) * 2 * F * 1000 / Cmem;        // [pA/pF]
        return ileak_stan;
}

double CSubcell::compute_I_ClCa_junc(void) {
        // return i_clca_junc_ave;
        return (i_clca_junc_ave / Cmem);  // [A/F]
}

double CSubcell::compute_I_ClCa_sl(void) {
        return (i_clca_sl_ave / Cmem);    // [A/F]
}

double CSubcell::update_ICaNa(double Nai, double NL) {
        // double pNa              = 0.75e-8 / 10285;       // [cm/sec]
        // double pNa              = 0.75e-8 / 8113.0 * 0.96245 * (1e-10);       // [cm/sec*F]
        double ibarna           = pNa * (v * F * 1e3 / rtf) *(0.75 * Nai * exp(v/rtf) - 0.75 * nao)  / (exp(v/rtf) - 1);
        
        double I_CaNa           = gca * (ibarna * NL * 0.45); // [A]
        
        return I_CaNa;
}

double CSubcell::update_ICaK(double NL) {
        double ibark    = pK * (v * F * 1e3 / rtf) * (0.75 * Ki * exp(v/rtf) - 0.75 * Ko) / (exp(v/rtf) - 1);
        
        double I_CaK    = gca * (ibark * NL * 0.45); // [A]
        
        return I_CaK;
}

double CSubcell::compute_ICaNa(void) {
        return (i_cana_ave / Cmem);
}
double CSubcell::compute_ICaK(void) {
        return (i_cak_ave / Cmem);
}

// output Na-related current (other than NCX)
double CSubcell::get_inak_junc(void) {
        return (i_nak_junc_ave / Cmem);
}

double CSubcell::get_inak_sl(void) {
        return (i_nak_sl_ave / Cmem);
}

double CSubcell::get_inabk_junc(void) {
        return (i_nabk_junc_ave / Cmem);
}

double CSubcell::get_inabk_sl(void) {
        return (i_nabk_sl_ave / Cmem);
}

double CSubcell::get_isk_junc(void) {
        return (i_sk_junc_ave / Cmem);
}
double CSubcell::get_isk_sl(void) {
        return (i_sk_sl_ave / Cmem);
}

double CSubcell::get_ina_junc(void) {
        return (i_na_junc_ave / Cmem);
}

double CSubcell::get_ina_sl(void) {
        return (i_na_sl_ave / Cmem);
}


double CSubcell::get_avg_nai(void)
{
        double sum = 0;
        for (int id = 0; id < nn; id++)
        sum += nai[id];
        return (sum / nn);
}

double CSubcell::get_avg_nas(void)
{
        double sum = 0;
        for (int id = 0; id < nn; id++)
        sum += nas[id];
        return (sum / nn);
}

double CSubcell::get_avg_nap(void)
{
        double sum = 0;
        for (int id = 0; id < nn; id++)
        sum += nap[id];
        return (sum / n);
}
void CSubcell::setboundary(int bcc)
{
        if (bc > 0) //corner
        {
                Jmaxx[0 + 0 * nx + 0 * nxny]            = 0;
                Jmaxx[0 + (ny - 1)*nx + (nz - 1)*nxny]  = 0;
                Jmaxx[0 + (ny - 1)*nx + 0 * nxny]       = 0;
                Jmaxx[0 + 0 * nx + (nz - 1)*nxny]       = 0;
                Jmaxx[(nx - 1) + (ny - 1)*nx + 0 * nxny]        = 0;
                Jmaxx[(nx - 1) + 0 * nx + (nz - 1)*nxny]        = 0;
                Jmaxx[(nx - 1) + 0 * nx + 0 * nxny]             = 0;
                Jmaxx[(nx - 1) + (ny - 1)*nx + (nz - 1)*nxny]   = 0;
        }

        if (bc > 1) //edge
        {
                //x fixed
                // #ifdef _OPENMP
                #pragma omp parallel for
                // #endif
                for (int j = 1; j < ny - 1; j++)
                {
                        Jmaxx[0 + j * nx + 0 * nxny]            = 0;
                        Jmaxx[(nx - 1) + j * nx + 0 * nxny]     = 0;
                        Jmaxx[0 + j * nx + (nz - 1)*nxny]       = 0;
                        Jmaxx[(nx - 1) + j * nx + (nz - 1)*nxny]= 0;
                }
                
                //y fixed
                // #ifdef _OPENMP
                #pragma omp parallel for
                // #endif
    
                for (int i = 1; i < (nx - 1); i++)
                {
                        Jmaxx[i + 0 * nx + 0 * nxny]            = 0;
                        Jmaxx[i + (ny - 1)*nx + 0 * nxny]       = 0;
                        Jmaxx[i + 0 * nx + (nz - 1)*nxny]       = 0;
                        Jmaxx[i + (ny - 1)*nx + (nz - 1)*nxny]  = 0;
                }
                
                #pragma ivdep
                #pragma vector always
                for (int k = 1; k < nz - 1; k++)
                {
                        Jmaxx[0 + 0 * nx + k * nxny]            = 0;
                        Jmaxx[0 + (ny - 1)*nx + k * nxny]       = 0;
                        Jmaxx[(nx - 1) + 0 * nx + k * nxny]     = 0;
                        Jmaxx[(nx - 1) + (ny - 1)*nx + k * nxny]= 0;
                }

        }
        
        if (bc > 2) //surface
        {
                //x fixed
                // #ifdef _OPENMP
                #pragma omp parallel for
                // #endif
                for (int j = 1; j < ny - 1; j++)
                {
                        #pragma ivdep
                        #pragma vector always
                        for (int k = 1; k < nz - 1; k++)
                        {
                                Jmaxx[0 + j * nx + k * nxny] = 0;
                                Jmaxx[(nx - 1) + j * nx + k * nxny] = 0;
                        }
                }
                
                //y fixed
                // #ifdef _OPENMP
                #pragma omp parallel for
                // #endif
                for (int i = 1; i < (nx - 1); i++)
                {
                        #pragma ivdep
                        #pragma vector always
                        for (int k = 1; k < nz - 1; k++)
                        {
                                Jmaxx[i + 0 * nx + k * nxny] = 0;
                                Jmaxx[i + (ny - 1)*nx + k * nxny] = 0;
                        }

                        //z fixed
                        #pragma ivdep
                        #pragma vector always
                        for (int j = 1; j < ny - 1; j++)
                        {
                                Jmaxx[i + j * nx + 0 * nxny] = 0;
                                Jmaxx[i + j * nx + (nz - 1)*nxny] = 0;
                        }
                }
        }

        if (bc > 3) //corner more
        {
                Jmaxx[1 + 1 * nx + 1 * nxny] = 0;
                Jmaxx[1 + (ny - 2)*nx + (nz - 2)*nxny] = 0;
                Jmaxx[1 + (ny - 2)*nx + 1 * nxny] = 0;
                Jmaxx[1 + 1 * nx + (nz - 2)*nxny] = 0;
                Jmaxx[(nx - 2) + (ny - 2)*nx + 1 * nxny] = 0;
                Jmaxx[(nx - 2) + 1 * nx + (nz - 2)*nxny] = 0;
                Jmaxx[(nx - 2) + 1 * nx + 1 * nxny] = 0;
                Jmaxx[(nx - 2) + (ny - 2)*nx + (nz - 2)*nxny] = 0;
        }
}

void CSubcell::setJmax(double newJmax) 
{
        Jmax = newJmax;
        for (int id = 0; id < n; id++)Jmaxx[id] = Jmax;
        setboundary(bc);
}
void CSubcell::resetBuffer(void) 
{
        #pragma omp parallel for
        #pragma ivdep
        #pragma vector always
        for (int id = 0; id < nn; id++)
        {
                cati[id] = kon * ci[id] * BT / (kon * ci[id] + koff);

                Cyto_MCa[id]     =   1.9;
                Cyto_MMg[id]     =   135;
                Cyto_TnCHc[id]   =   0.112 * 1e3;
                Cyto_TnCHm[id]   =   0.01 * 1e3;
                Cyto_SR[id]      =   3.3;
                Cyto_cam[id]     =   0.4;
                
                SL_cam[id]       = 0.1*0.4;
                SL_SR[id]        = 0.1*3.3;
                SL_SLL[id]       = 0.1*12;
                SL_SLH[id]       = 0.1*130;
                
                #ifdef ___EGTA
                        caEGTAi[id] = konEGTA * ci[id] * BEGTA / (konEGTA * ci[id] + koffEGTA);
                #endif

                SL_NaB[id]      = 0.8;
        }
        #pragma omp parallel for

        #pragma ivdep
        #pragma vector always
        for (int id = 0; id < n; id++)
        {
                cleft_cam[id]   = 0.1*0.4;          // [uM]
                cleft_SR[id]    = 0.1*3.3;          // [uM]
                cleft_SLL[id]   = 0.1*12;           // [uM]
                cleft_SLH[id]   = 0.1*130;          // [uM]     
                cleft_NaB[id]   = 1.8;              
        }

}

void CSubcell::set_lateral_Ttubule(void) {
        int nxny = nx * ny;
        int nxyz = nx * ny * nz;

        int counter = 0;

        for (int id = 0; id < (nxyz); id++) {
                if ((id % nx < layer) || (id % nx > (nx - 1 - layer))) {
                        set_Ttubule(id, 1.0); counter++;
                }
                else if ((((id % nxny - id % nx) / nx) < layer) || (((id % nxny - id % nx) / nx) > (ny - 1 - layer))) {
                        set_Ttubule(id, 1.0); counter++;
                }
                else if ((((id - id % nxny) / nxny) < layer) || (((id - id % nxny) / nxny) > (nz - 1 - layer))) {
                        set_Ttubule(id, 1.0); counter++;
                }
                else    set_Ttubule(id, 0.0);
        }

        num_tubule = counter;
}


void CSubcell::set_TubuleScale(void) {
        int nxyz = nx * ny * nz;

        int counter = 0;

        for (int id = 0; id < (nxyz); id++) {
                if (tubule_flag[id]==1)      counter++;
        }

        TubuleScale =  1.0;
}

void CSubcell::generate_tubular_map() {
        int nxyz = nx * ny * nz;

        int counter = 0;

        for (int id = 0; id < nxyz; id++) {
                if (tubule_flag[id]==1)      counter++;
        }
        num_tubule = counter;
        Tubular_map = new int [counter];
        counter = 0;
        for (int id = 0; id < nxyz; id++) {
                if (tubule_flag[id]==1)      {
                        Tubular_map[counter] = id;
                        counter++;
                }
        }

}

void CSubcell::set_new_Cmem(void) {
        double specific_cmem_surface    = 1e-2;         // [pF/um^2]
        double specific_cmem_tubule     = 0.56e-2;      // [pF/um^2]

        double cru_length       = 1.84;
        double cru_width        = 0.9;
        double pi               = 3.1415926;
        double diameter_tubule  = 0.3;
        double at_tt_ratio      = 6.0;
        double length   = cru_length * nx;
        double width    = cru_width * ny;
        double depth    = cru_width * nz;

        double area_surf        = 2.0 * (length * width + length * depth + width * depth);
        double area_single_tub  = pi * diameter_tubule * ( at_tt_ratio/(at_tt_ratio + 2.0) * cru_length + 2.0/(at_tt_ratio + 2.0) * cru_width );
        double num_surf_cru     = (nx * ny ) * 2 + ( (nz - 2) * nx ) * 2 + ( (ny - 2) * (nz - 2) )* 2;  // 3130 (one layer) rather than 5644 (two layers)

        // the most densely tubulated cell has 8998 coupled CRUs - Assumption Cmem = 110 [pF]
        // e.g., the de-tubulated cell has 3130 surface CRUs and (5644-3130) second-layer inner CRUs
        // double num_inner_cru_dense_tub  = num_dense_tubulated_CRU * 1.0 - num_surf_cru;
        // double area_sum_dense_tub       = num_inner_cru_dense_tub * area_single_tub + area_surf;
        double num_inner_cru    = num_tubule - num_surf_cru;
        Cmem    = (area_surf * specific_cmem_surface + num_inner_cru * area_single_tub * specific_cmem_tubule) * 1e-12; // [F]
}

void CSubcell::get_steady_state(){
        // ci, cs, cp, cnsr, cjsr; nai, nas, nap
        steady_state_id         = 8 * nn;
        steady_state_CaNa       = new double [steady_state_id];
        memset(steady_state_CaNa, 0, steady_state_id*sizeof(double));
        // RyR_1, RyR_2, RyR_3; state1_ina, state2_ina, state3_ina
        steady_state_id         = 7 * nn + 3;
        steady_state_currents   = new double [steady_state_id];
        memset(steady_state_currents, 0, steady_state_id*sizeof(double));
        // cleft_cam, cleft_SR, cleft_SLL, cleft_SLH
        // SL_cam, SL_SR, SL_SLL, SL_SLH
        // cati, Cyto_MCa, Cyto_MMg, Cyto_SR, Cyto_cam, Cyto_TnCHc, Cyto_TnCHm
        // caEGTAi, caEGTAs
        steady_state_id         = 17 * nn;
        steady_state_buffers    = new double [steady_state_id];
        memset(steady_state_buffers, 0, steady_state_id*sizeof(double));

        // Steady State files
        // file102 : ion concentration  // ci, cs, cp, cnsr, cjsr; nai, nas, nap
        char            filename102 [1000];
        sprintf         (filename102, "./steady_state_output/steady_state_Ca_Na.txt" );
        // file103 : currents  // RyR_1, RyR_2, RyR_3; state1_ina, state2_ina, state3_ina // LTCC 4 channels
        char            filename103 [1000];
        sprintf         (filename103, "./steady_state_output/steady_state_currents.txt" );
        // file104 : buffers
        // cleft_cam, cleft_SR, cleft_SLL, cleft_SLH
        // SL_cam, SL_SR, SL_SLL, SL_SLH
        // cati, Cyto_MCa, Cyto_MMg, Cyto_SR, Cyto_cam, Cyto_TnCHc, Cyto_TnCHm
        // caEGTAi, caEGTAs
        char            filename104 [1000];
        sprintf         (filename104, "./steady_state_output/steady_state_buffers.txt" );

        int single_state_CaNa_id        = 0;
        int single_state_currents_id    = 0;
        int single_state_buffers_id     = 0;
        
        for(int i = 0; i < nn; ++i){

                // ci, cs, cp, cnsr, cjsr; nai, nas, nap
                single_state_CaNa_id    = 0 + i * 8;
                steady_state_CaNa[single_state_CaNa_id] = ci[i];
                
                single_state_CaNa_id    = 1 + i * 8;
                steady_state_CaNa[single_state_CaNa_id] = cs[i];

                single_state_CaNa_id    = 2 + i * 8;
                steady_state_CaNa[single_state_CaNa_id] = cp[i];

                single_state_CaNa_id    = 3 + i * 8;
                steady_state_CaNa[single_state_CaNa_id] = cnsr[i];

                single_state_CaNa_id    = 4 + i * 8;
                steady_state_CaNa[single_state_CaNa_id] = cjsr[i];

                single_state_CaNa_id    = 5 + i * 8;
                steady_state_CaNa[single_state_CaNa_id] = nai[i];

                single_state_CaNa_id    = 6 + i * 8;
                steady_state_CaNa[single_state_CaNa_id] = nas[i];

                single_state_CaNa_id    = 7 + i * 8;
                steady_state_CaNa[single_state_CaNa_id] = nap[i];

                // RyR_1, RyR_2, RyR_3; state1_ina, state2_ina, state3_ina; LTCC_State[0] ... LTCC_State[4-1]
                single_state_currents_id = 0 + i * 3;
                steady_state_currents[single_state_currents_id] = RyR_vec[i].RyR_1;
                
                single_state_currents_id = 1 + i * 3;
                steady_state_currents[single_state_currents_id] = RyR_vec[i].RyR_2;

                single_state_currents_id = 2 + i * 3;
                steady_state_currents[single_state_currents_id] = RyR_vec[i].RyR_3;

                // cleft_cam, cleft_SR, cleft_SLL, cleft_SLH
                // SL_cam, SL_SR, SL_SLL, SL_SLH
                // cati, Cyto_MCa, Cyto_MMg, Cyto_SR, Cyto_cam, Cyto_TnCHc, Cyto_TnCHm
                // caEGTAi, caEGTAs
                single_state_buffers_id = 0 + i * 17;
                steady_state_buffers[single_state_buffers_id]   = cleft_cam[i]; 
                
                single_state_buffers_id = 1 + i * 17;
                steady_state_buffers[single_state_buffers_id]   = cleft_SR[i];  
                
                single_state_buffers_id = 2 + i * 17;
                steady_state_buffers[single_state_buffers_id]   = cleft_SLL[i]; 
                
                single_state_buffers_id = 3 + i * 17;
                steady_state_buffers[single_state_buffers_id]   = cleft_SLH[i]; 
                
                single_state_buffers_id = 4 + i * 17;
                steady_state_buffers[single_state_buffers_id]   = SL_cam[i];    
                
                single_state_buffers_id = 5 + i * 17;
                steady_state_buffers[single_state_buffers_id]   = SL_SR[i];     
                
                single_state_buffers_id = 6 + i * 17;
                steady_state_buffers[single_state_buffers_id]   = SL_SLL[i];    
                
                single_state_buffers_id = 7 + i * 17;
                steady_state_buffers[single_state_buffers_id]   = SL_SLH[i];    
                
                single_state_buffers_id = 8 + i * 17;
                steady_state_buffers[single_state_buffers_id]   = cati[i];      
                
                single_state_buffers_id = 9 + i * 17;
                steady_state_buffers[single_state_buffers_id]   = Cyto_MCa[i];  
                
                single_state_buffers_id = 10 + i * 17;
                steady_state_buffers[single_state_buffers_id]   = Cyto_MMg[i];  
                
                single_state_buffers_id = 11 + i * 17;
                steady_state_buffers[single_state_buffers_id]   = Cyto_SR[i];   
                
                single_state_buffers_id = 12 + i * 17;
                steady_state_buffers[single_state_buffers_id]   = Cyto_cam[i];  
                
                single_state_buffers_id = 13 + i * 17;
                steady_state_buffers[single_state_buffers_id]   = Cyto_TnCHc[i];
                
                single_state_buffers_id = 14 + i * 17;
                steady_state_buffers[single_state_buffers_id]   = Cyto_TnCHm[i];
                
                #ifdef ___EGTA                
                        single_state_buffers_id = 15 + i * 17;
                        steady_state_buffers[single_state_buffers_id]   = caEGTAi[i];
                        
                        single_state_buffers_id = 16 + i * 17;
                        steady_state_buffers[single_state_buffers_id]   = caEGTAs[i];

                #else
                        single_state_buffers_id = 15 + i * 17;
                        steady_state_buffers[single_state_buffers_id]   = 0.0;
                        
                        single_state_buffers_id = 16 + i * 17;
                        steady_state_buffers[single_state_buffers_id]   = 0.0;
                #endif
        }
        // RyR_1, RyR_2, RyR_3; state1_ina, state2_ina, state3_ina; LTCC_State[0] ... LTCC_State[4-1]
        single_state_currents_id        = 0 + nn * 3; 
        steady_state_currents[single_state_currents_id]         = Ina.state1_ina;
        single_state_currents_id        = 1 + nn * 3; 
        steady_state_currents[single_state_currents_id]         = Ina.state2_ina;
        single_state_currents_id        = 2 + nn * 3; 
        steady_state_currents[single_state_currents_id]         = Ina.state3_ina;

        for(int i = 0; i < nn; ++i){
                single_state_currents_id        = 0 + i * 4 + (nn+1) * 3; 
                steady_state_currents[single_state_currents_id]         = LTCC_Markov_vec[i].LTCC_State[0];
                single_state_currents_id        = 1 + i * 4 + (nn+1) * 3; 
                steady_state_currents[single_state_currents_id]         = LTCC_Markov_vec[i].LTCC_State[1];
                single_state_currents_id        = 2 + i * 4 + (nn+1) * 3; 
                steady_state_currents[single_state_currents_id]         = LTCC_Markov_vec[i].LTCC_State[2];
                single_state_currents_id        = 3 + i * 4 + (nn+1) * 3; 
                steady_state_currents[single_state_currents_id]         = LTCC_Markov_vec[i].LTCC_State[3];
        }

        output_double_array_txt(filename102, steady_state_CaNa, 8 * n);
        output_double_array_txt(filename103, steady_state_currents, 7 * n + 3);
        output_double_array_txt(filename104, steady_state_buffers, 17 * n);

        if(steady_state_CaNa) delete [] steady_state_CaNa;
        if(steady_state_currents) delete [] steady_state_currents;
        if(steady_state_buffers) delete [] steady_state_buffers;

        steady_state_CaNa = nullptr;
        steady_state_currents = nullptr;
        steady_state_buffers = nullptr;
}

void CSubcell::set_steady_state(){
        // Steady State files
        // file106 : ion concentration  // ci, cs, cp, cnsr, cjsr; nai, nas, nap
        char            filename106 [1000];
        sprintf         (filename106, "./steady_state_init/steady_state_Ca_Na.txt" );

        // file107 : currents  // RyR_1, RyR_2, RyR_3; state1_ina, state2_ina, state3_ina // LTCC 4 channels
        char            filename107 [1000];
        sprintf         (filename107, "./steady_state_init/steady_state_currents.txt" );

        // file108 : buffers
        // cleft_cam, cleft_SR, cleft_SLL, cleft_SLH
        // SL_cam, SL_SR, SL_SLL, SL_SLH
        // cati, Cyto_MCa, Cyto_MMg, Cyto_SR, Cyto_cam, Cyto_TnCHc, Cyto_TnCHm
        // caEGTAi, caEGTAs
        char            filename108 [1000];
        sprintf         (filename108, "./steady_state_init/steady_state_buffers.txt" );

        // ci, cs, cp, cnsr, cjsr; nai, nas, nap
        steady_state_id         = 8 * nn;
        steady_state_CaNa       = new double [steady_state_id];
        memset(steady_state_CaNa, 0, steady_state_id*sizeof(double));
        // RyR_1, RyR_2, RyR_3; state1_ina, state2_ina, state3_ina
        steady_state_id         = 7 * nn + 3;
        steady_state_currents   = new double [steady_state_id];
        memset(steady_state_currents, 0, steady_state_id*sizeof(double));
        // cleft_cam, cleft_SR, cleft_SLL, cleft_SLH
        // SL_cam, SL_SR, SL_SLL, SL_SLH
        // cati, Cyto_MCa, Cyto_MMg, Cyto_SR, Cyto_cam, Cyto_TnCHc, Cyto_TnCHm
        // caEGTAi, caEGTAs
        steady_state_id         = 17 * nn;
        steady_state_buffers    = new double [steady_state_id];
        memset(steady_state_buffers, 0, steady_state_id*sizeof(double));

       // // Set Steady State
        read_double_from_txt( filename106, steady_state_CaNa, 8 * n );
        read_double_from_txt( filename107, steady_state_currents, 7 * n + 3 );
        read_double_from_txt( filename108, steady_state_buffers, 17 * n );    


        int single_state_CaNa_id        = 0;
        int single_state_currents_id    = 0;
        int single_state_buffers_id     = 0;


        #pragma omp parallel for
        for(int i = 0; i < nn; ++i){

                // ci, cs, cp, cnsr, cjsr; nai, nas, nap
                single_state_CaNa_id    = 0 + i * 8;
                ci[i]   = steady_state_CaNa[single_state_CaNa_id];
                
                single_state_CaNa_id    = 1 + i * 8;
                cs[i]   = steady_state_CaNa[single_state_CaNa_id];

                single_state_CaNa_id    = 2 + i * 8;
                cp[i]   = steady_state_CaNa[single_state_CaNa_id];

                single_state_CaNa_id    = 3 + i * 8;
                cnsr[i] = steady_state_CaNa[single_state_CaNa_id];

                single_state_CaNa_id    = 4 + i * 8;
                cjsr[i] = steady_state_CaNa[single_state_CaNa_id];

                // RyR_1, RyR_2, RyR_3; state1_ina, state2_ina, state3_ina; LTCC_State[0] ... LTCC_State[4-1]
                single_state_currents_id = 0 + i * 3;
                RyR_vec[i].RyR_1        = steady_state_currents[single_state_currents_id];
                
                single_state_currents_id = 1 + i * 3;
                RyR_vec[i].RyR_2        = steady_state_currents[single_state_currents_id];

                single_state_currents_id = 2 + i * 3;
                RyR_vec[i].RyR_3        = steady_state_currents[single_state_currents_id];

                // cleft_cam, cleft_SR, cleft_SLL, cleft_SLH
                // SL_cam, SL_SR, SL_SLL, SL_SLH
                // cati, Cyto_MCa, Cyto_MMg, Cyto_SR, Cyto_cam, Cyto_TnCHc, Cyto_TnCHm
                // caEGTAi, caEGTAs
                single_state_buffers_id = 0 + i * 17;
                cleft_cam[i]    = steady_state_buffers[single_state_buffers_id];
                
                single_state_buffers_id = 1 + i * 17;
                cleft_SR[i]     = steady_state_buffers[single_state_buffers_id];
                
                single_state_buffers_id = 2 + i * 17;
                cleft_SLL[i]    = steady_state_buffers[single_state_buffers_id];
                
                single_state_buffers_id = 3 + i * 17;
                cleft_SLH[i]    = steady_state_buffers[single_state_buffers_id];
                
                single_state_buffers_id = 4 + i * 17;
                SL_cam[i]       = steady_state_buffers[single_state_buffers_id];
                
                single_state_buffers_id = 5 + i * 17;
                SL_SR[i]        = steady_state_buffers[single_state_buffers_id];
                
                single_state_buffers_id = 6 + i * 17;
                SL_SLL[i]       = steady_state_buffers[single_state_buffers_id];
                
                single_state_buffers_id = 7 + i * 17;
                SL_SLH[i]       = steady_state_buffers[single_state_buffers_id];
                
                single_state_buffers_id = 8 + i * 17;
                cati[i]         = steady_state_buffers[single_state_buffers_id];
                
                single_state_buffers_id = 9 + i * 17;
                Cyto_MCa[i]     = steady_state_buffers[single_state_buffers_id];
                
                single_state_buffers_id = 10 + i * 17;
                Cyto_MMg[i]     = steady_state_buffers[single_state_buffers_id];
                
                single_state_buffers_id = 11 + i * 17;
                Cyto_SR[i]      = steady_state_buffers[single_state_buffers_id];
                
                single_state_buffers_id = 12 + i * 17;
                Cyto_cam[i]     = steady_state_buffers[single_state_buffers_id];
                
                single_state_buffers_id = 13 + i * 17;
                Cyto_TnCHc[i]   = steady_state_buffers[single_state_buffers_id];
                
                single_state_buffers_id = 14 + i * 17;
                Cyto_TnCHm[i]   = steady_state_buffers[single_state_buffers_id];
                
                #ifdef ___EGTA                
                        single_state_buffers_id = 15 + i * 17;
                        caEGTAi[i]      = steady_state_buffers[single_state_buffers_id];                     
                        single_state_buffers_id = 16 + i * 17;
                        caEGTAs[i]      = steady_state_buffers[single_state_buffers_id];                   
                #endif
        }

        // RyR_1, RyR_2, RyR_3; state1_ina, state2_ina, state3_ina; LTCC_State[0] ... LTCC_State[4-1]
        single_state_currents_id        = 0 + nn * 3; 
        Ina.state1_ina  = steady_state_currents[single_state_currents_id];
        single_state_currents_id        = 1 + nn * 3; 
        Ina.state2_ina  = steady_state_currents[single_state_currents_id];
        single_state_currents_id        = 2 + nn * 3; 
        Ina.state3_ina  = steady_state_currents[single_state_currents_id];

        for(int i = 0; i < nn; ++i){
                single_state_currents_id        = 0 + i * 4 + (nn+1) * 3; 
                LTCC_Markov_vec[i].LTCC_State[0]        = steady_state_currents[single_state_currents_id];
                single_state_currents_id        = 1 + i * 4 + (nn+1) * 3; 
                LTCC_Markov_vec[i].LTCC_State[1]        = steady_state_currents[single_state_currents_id];
                single_state_currents_id        = 2 + i * 4 + (nn+1) * 3; 
                LTCC_Markov_vec[i].LTCC_State[2]        = steady_state_currents[single_state_currents_id];
                single_state_currents_id        = 3 + i * 4 + (nn+1) * 3; 
                LTCC_Markov_vec[i].LTCC_State[3]        = steady_state_currents[single_state_currents_id];
        }

        if(steady_state_CaNa) delete [] steady_state_CaNa;
        if(steady_state_currents) delete [] steady_state_currents;
        if(steady_state_buffers) delete [] steady_state_buffers;


        steady_state_CaNa = nullptr;
        steady_state_currents = nullptr;
        steady_state_buffers = nullptr;
}

void CSubcell::set_na_after_pausing(double time_after_pausing, double bcl)
{
        double a        = 15;   // [mM]
        double b        = 0.6;  // [Hz]
        // When bcl = 0.33s, [Na]i = 11.14 mM; When bcl = 1s, [Na]i = 9.375 mM
        double na_clamped  = a / (1 + b * sqrt(bcl * 1e-3) );

        double decay_percent_rate_dependance = -1.529e-9;       // 1/[ms]
        double fastest_decay_percent = 2.631e-6;
        double decay_percent =  decay_percent_rate_dependance * bcl + fastest_decay_percent; // [100%]/[ms]
        double na_dynamic = (1 - decay_percent * time_after_pausing) * na_clamped; // [mM]

        for (int id = 0; id < nn; id++)
        {
                nai[id]         = na_dynamic;
                nas[id]         = na_dynamic;
        }
        for (int id = 0; id < n; id++)
        {
                nap[id]         = na_dynamic;
        }
}

void CSubcell::set_LCC(int id, double nLCC_state) 
{
        map_nlcc[id]    = nLCC_state;
}

void CSubcell::set_RyR(int id, double nRyR_state) 
{
        map_nryr[id]    = nRyR_state;
}

void CSubcell::set_ncx(int id, double ncx_scale)
{
        map_ncx[id]     = ncx_scale;
}

void CSubcell::set_serca(int id, double serca_scale)
{
        map_serca[id]   = serca_scale;
}