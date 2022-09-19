#define ___NCX
#define __RYR_UNIFORM

#define ___CPDIFF

#define ___UNIFORM

#ifndef ___SUBCELL_H
#define ___SUBCELL_H
#define _USE_MATH_DEFINES

#include <iostream>
using namespace std;
#include <cmath>
#include <math.h>
#include <fstream>
#include <cfloat>
#include <climits>
#include <cstdlib>
#include "LTCC_Markov.hpp"
#include "RyR.hpp"
#include "xor_rand.hpp"
#include <vector>
#include "inaca.hpp"
#include "ina.hpp"
#include "LTCC_unitary.hpp"
#include "cell.hpp"
#include "input_output.h"

class CSubcell {
        private:
                int nxny;

                int finemesh3;
                int nnx;        //the number of finemesh (x)
                int nny;        //the number of finemesh (y)
                int nnz;        //the number of finemesh (z)
                int nnxnny;

                double vi;
                double vs;

                // static const double vp_ave;
                double tausT;
                double tausL;
                double taumninv;
                double tauiT;
                double tauiL;
                double taunsrT;
                double taunsrL;

                double tausT_na;
                double tausL_na;
                double taumninv_na;
                double tauiT_na;
                double tauiL_na;

                void computeIci(void);
                void computeIcnsr(void);

                double *Ici, *Icnsr;
                double *Inai, *Inas;

                #ifdef ___NO_CS_BUFFER
                        double *csmn;
                        void computecsmn(void);
                        double *nasmn;
                        void computenasmn(void);
                #else
                        void computeIcs(void);
                        double *Ics;
                        double *Idps;
                #endif

                void computeInas(void);

                int bino(double num, double p, int ii);
                double calcvp(double mean, double std, double lim1, double lim2, int ii);

                double dt;
                double vup, kup, KNSR;
                double vnaca;
                double Jmax;
                double taup;
                double tausi;
                double tautr;
                double gca;

                double taup_na;
                double tausi_na;
                double tautr_na;

                double cao;

                double gleak;
                double gcabk;
                double qslcap;

                double kk_pow_10;

                double Ku;      // CSQN-unbound opening rate
                double Kb;      // CSQN-bound opening rate
                double tauu;    // CSQN unbinding timescale
                double taub;    // CSQN binding timescale
                double tauc1;
                double tauc2;
                double BCSQN0;  // CSQN concentration

                int seed;
                bool initialized;
                int bc;

                double Kcp;
                double pedk12;
                double pedk43;
                int NoCaL;      //# of CaL channel per CaRU

                double MaxSR;
                double MinSR;
                double ec50SR;
                double hkosrca;

                void delarray(void);

                #ifdef ___EGTA
                        double BEGTA;
                #endif
        public:
                int layer;      // # layers with surface channels on each side

                #ifndef __RYR_UNIFORM
                        double *sigma;
                #endif

                unsigned int *xsx, *xsy, *xsz, *xsw;
                int nx;         // the number of CRU (x)
                int ny;         // the number of CRU (y)
                int nz;         // the number of CRU (z)
                int finemesh;
                double xi;

                int n;          // the number of CRU x * y * z
                int nn;         // the number of fine-mesh

                double vnsr;

                std::vector<LTCC_Markov> LTCC_Markov_vec;
                std::vector<xor_rand> random_num;

                std::vector<RyR> RyR_vec;

                int * LTCC_Open_num;

                int Total_tubule_num;
                int *Tubular_map;

                ina Ina;

                CSubcell(int sizex = 55, int sizey = 17, int sizez = 11, int fmesh = 1, double xii = 1, int seed = 0);

                virtual ~CSubcell();

                double *ci, *cs, *cp, *cjsr, *cnsr, *cati, *cats;
                double *cleft_cam, *cleft_SR, *cleft_SLL, *cleft_SLH;
                double *SL_cam, *SL_SR, *SL_SLL, *SL_SLH;
                double *Cyto_MCa, *Cyto_MMg, *Cyto_SR, *Cyto_cam, *Cyto_TnCHc, *Cyto_TnCHm;     // Cytosolic Tropnin C is "cati"
                double *JSR_CSQN;

                double *cscp1, *cscp2, *Itr;

                double *nai, *nas, *nap;
                double *SL_NaB, *cleft_NaB;
                double NCX_fit_scale;

                #ifdef ___DETERMINISTIC
                        double *c1, *c2, *i1ca, *i1ba, *i2ca, *i2ba, *fryr1, *fryr2, *fryr3;
                #else
                        int *y, *ryr1, *ryr2, *ryr3, *nryr;
                #endif

                double *j_serca, *j_ryr;
                double *ica_array, *ncx_array;

                #ifdef ___NCX
                        double *ncx_array_p;
                #endif

                double *icabk_array, *jpca_array, *I_ClCa_sl_array;
                double *I_nak_sl_array, *I_nabk_sl_array, *I_Na_sl_array;
                double *I_SK_sl_array, *I_SK_junc_array;

                double *ICaNa_array,  *ICaK_array;
                double *I_ClCa_junc_array;
                double *I_nak_junc_array, *I_nabk_junc_array, *I_Na_junc_array;
                double *JSR_buffer_factors;

                double *vp;
                double *Jmaxx;

                int *tubule_flag;
                double hp_ryr;
                double TubuleScale;

                int *crupos;
                double AF, ISO, gCa_iso_af_scale;
                double pNa, pK;
                double v;

                void pace(double voltage);
                void pace2(double voltage, double local_dt);
                void  update_LTCC_RyR(double voltage, double local_dt);
                
                void set_na_after_pausing(double time_after_pausing, double bcl);
                
                double get_jSR_inst_buffering_simp(double const & CaSR)
                {
                        double term1    = std::pow(double(CaSR / 1000.0), double(hh-1));
                        double term2    = term1 * (CaSR / 1000.0);
                        double term3    = std::pow(double(KK / 1000.0), double(hh));
                        
                        double rho              = rhoinf / (term3 / term2 + 1.0);
                        if (rho<0.0000001)      rho = 0.0000001;        // to avoid div small (for CaSR<200)
                        double rho_BCSQN        = rho * BCSQN;
                        
                        double MM               = (sqrt(1 + 8 * rho_BCSQN) - 1) / (4 * rho_BCSQN);
                        double ncjsr            = MM*nM + (1-MM)*nD;

                        double rhopri           = 0.001 * (hh * rhoinf * term1 * (term3 + term2) - rhoinf * term2 * hh * term1) / (term3 + term2) / (term3 + term2);
                        double rhopri_BCSQN     = rhopri * BCSQN;

                        double dMdc     = ((0.5 / sqrt(1 + 8 * rho_BCSQN) * 8 * rhopri_BCSQN * 4 * rho_BCSQN) 
                                        - (sqrt(1 + 8 * rho_BCSQN) - 1) * 4 * rhopri_BCSQN) / (4 * rho_BCSQN) / (4 * rho_BCSQN);
                        double dndc     = dMdc * (nM-nD);

                        return 1.0 / (1.0 + (kd_csqn * BCSQN*ncjsr + dndc*(CaSR*kd_csqn + CaSR*CaSR)) / ((kd_csqn + CaSR) * (kd_csqn + CaSR)));
                }  

                void init(double initci = 0.1, double initcj = 500, double initnai = 10);

                //set Nerst parameters
                double iupave, icaave, incxave, irave, ileakave, icabkave, islcapave;
                double incx_sl_ave, incx_junc_ave;
                double ir_ss, ir_ct;
                double i_clca_junc_ave, i_clca_sl_ave;
                double i_cana_ave, i_cak_ave;

                double i_nak_junc_ave, i_nak_sl_ave, i_nabk_junc_ave, i_nabk_sl_ave, i_na_junc_ave, i_na_sl_ave;
                double i_sk_junc_ave, i_sk_sl_ave;

                double count_tubule;

                double *PoSpark;
                double INafast_gating;

                double *map_serca, *map_ncx;
                int *map_nryr, *map_nlcc;

                double compute_avg_ci(void);
                double compute_avg_cs(void);
                double compute_avg_cnsr(void);
                double compute_avg_cp(void);
                double compute_avg_cjsr(void);

                double compute_ica(void);
                double compute_incx(void);
                double compute_icabk(void);
                double compute_ipca(void);
                double compute_incx_junc(void);
                double compute_incx_sl(void);
                double compute_ir(void);
                double compute_serca(void);
                double compute_ileak(void);
                double compute_I_ClCa_junc(void);
                double compute_I_ClCa_sl(void);
                double update_ICaNa(double Nai, double NL);
                double update_ICaK(double NL);
                // void update_ina_parameters(void);
                double compute_ICaNa(void);
                double compute_ICaK(void);

                double get_inak_junc(void);
                double get_inak_sl(void);
                double get_inabk_junc(void);
                double get_inabk_sl(void);
                double get_isk_junc(void);
                double get_isk_sl(void);
                double get_ina_junc(void);
                double get_ina_sl(void);
                double get_avg_nai(void);
                double get_avg_nas(void);
                double get_avg_nap(void);

                int num_open_ICaL;
                int steady_state_id;
                double *steady_state_CaNa, *steady_state_currents, *steady_state_buffers;
                void get_steady_state();
                void set_steady_state();
                CSubcell& operator=(const CSubcell& sc);//  the type of operator function and sc are both CSubcell

                #ifndef __UNIFORM_TUBULE
                        void sethyryr(int newhp_ryr) {hp_ryr = newhp_ryr;};
                        double gethyryr(void) {return hp_ryr;};
                #endif

                void setlayer(int newlayer) {layer = newlayer;};
                int getlayer(void) {return layer;}; // *) xw: add the layer setting functon
                void setdt(double newdt) {dt = newdt;}; //set the time step as 0.01, in pace.cc

                void set_Ttubule(int id, int tubule) {tubule_flag[id] = tubule;};

                void set_LCC(int id, double nLCC_state);
                void set_RyR(int id, double nRyR_state);
                void set_ncx(int id, double ncx_scale);
                void set_serca(int id, double serca_scale);

                double getdt(void) {return dt;};
                void setgca(double newgca) {gca = newgca;};
                double getgca(void) {return gca;};
                void setgncx(double newgncx) {vnaca = newgncx;};
                double getgncx(void) {return vnaca;};
                void setJmax(double newJmax);
                double getJmax(void) {return Jmax;};
                void setvup(double newvup) {vup = newvup;};
                double getvup(void) {return vup;};
                void setkup(double newkup) {kup = newkup;};
                double getkup(void) {return kup;};
                void setKNSR(double newKNSR) {KNSR = newKNSR;};
                double getKNSR(void) {return KNSR;};

                void setgleak(double newgleak) {gleak = newgleak;};
                double getgleak(void) {return gleak;};
                void setcao(double newcao) {cao = newcao;};
                double getcao(void) {return cao;};
                double getBCSQN(void) {return BCSQN;};

                void setgcabk(double newgcabk) {gcabk = newgcabk;};
                double getgcabk(void) {return gcabk;};
                void setqslcap(double newqslcap) {qslcap = newqslcap;};
                double getqslcap(void) {return qslcap;};

                double getKc(void) {return Kc;};
                double getnM(void) {return nM;};
                double getnD(void) {return nD;};
                double gethh(void) {return hh;};
                double getKK(void) {return KK;};
                double getrhoinf(void) {return rhoinf;};
                void setKu(double newKu) {Ku = newKu;};
                double getKu(void) {return Ku;};
                void setKb(double newKb) {Kb = newKb;};
                double getKb(void) {return Kb;};
                void settauu(double newtauu) {tauu = newtauu;};
                double gettauu(void) {return tauu;};
                void settaub(double newtaub) {taub = newtaub;};
                double gettaub(void) {return taub;};
                void settauc1(double newtauc1) {tauc1 = newtauc1;};
                double gettauc1(void) {return tauc1;};
                void settauc2(double newtauc2) {tauc2 = newtauc2;};
                double gettauc2(void) {return tauc2;};
                void setBCSQN0(double newBCSQN0) {BCSQN0 = newBCSQN0;};
                double getBCSQN0(void) {return BCSQN0;};

                void setvi(double newvi) {vi = newvi;}
                double getvi(void) {return vi;}

                void setvs(double newvs) 
                {
                        vs = newvs;
                        for (int id = 0; id < n; id++) 
                        {
                                cscp2[crupos[id]] = 1 / taup * vp[id] / vs;
                        }
                }

                double getvs(void) {return vs;}

                double gettausi(void) {return tausi;}
                double settausi(double  newval) {tausi = newval; return tausi;}

                double gettautr(void) {return tautr;}
                double settautr(double  newval) {tautr = newval; return tautr;}

                void setKcp(double newKcp) {Kcp = newKcp;};
                double getKcp(void) {return Kcp;};
                void setpedk12(double newpedk12) {pedk12 = newpedk12;};
                double getpedk12(void) {return pedk12;};
                void setpedk43(double newpedk43) {pedk43 = newpedk43;};
                double getpedk43(void) {return pedk43;};

                double getMaxSR(void) {return MaxSR;};
                void setMaxSR(double newMaxSR) {MaxSR = newMaxSR;};
                double getMinSR(void) {return MinSR;};
                void setMinSR(double newMinSR) {MinSR = newMinSR;};
                double getec50SR(void) {return ec50SR;};
                void setec50SR(double newec50SR) {ec50SR = newec50SR;};
                double gethkosrca(void) {return hkosrca;};
                void sethkosrca(double newhkosrca) {hkosrca = newhkosrca;};


                // Diffusion
                double gettausT(void) {return tausT;}
                double settausT(double  newval) {tausT = newval; taumninv = 2 / tausL + 4 / tausT; return tausT;}
                double gettausL(void) {return tausL;}
                double settausL(double  newval) {tausL = newval; taumninv = 2 / tausL + 4 / tausT; return tausL;}
                double gettauiT(void) {return tauiT;}
                double settauiT(double  newval) {tauiT = newval; return tauiT;}
                double gettauiL(void) {return tauiL;}
                double settauiL(double  newval) {tauiL = newval; return tauiL;}
                double gettaunsrT(void) {return taunsrT;}
                double settaunsrT(double  newval) {taunsrT = newval; return taunsrT;}
                double gettaunsrL(void) {return taunsrL;}
                double settaunsrL(double  newval) {taunsrL = newval; return taunsrL;}

                void srand(int sed = 0);
                void setboundary(int bcc = 3);
                void resetBuffer(void);

                void check_crupos() {
                        for (int i = 0; i < n; ++i)
                        {
                                std::cout << i << " " << crupos[i] << std::endl;
                        }
                }

                void set_lateral_Ttubule();
                void set_TubuleScale();
                void output_Ttubule_map();
                void output_nRyR_map();

                int num_tubule;

                void update_single_CRU(double V);

                void generate_tubular_map();
                void set_new_Cmem();

                const double F          = 96.5;         // [C/mmol]
                const double R          = 8.314;        // [J/mol*K]
                const double T          = 310.0;        // [K]
                const double rtf        = R * T / F;    // ~26.5
                const double rtf2       = R * T / (2.0 * F);

                double GClCa      = 0.0548 / 8113.0 * 0.96245 * (1e-10);        // [mS/uF*F]
                const double KdClCa     = 100;                                  // [mM] -> [uM]
                const double Cli        = 15.0;                                 // Intracellular Cl  [mM]
                const double Clo        = 150.0;                                // Extracellular Cl  [mM]
                const double ecl        = (1.0 * rtf) * log(Cli / Clo);         // [mV]
                const double Ki         = 120.0;                                // [mM]
                const double Ko         = 5.4;                                  // Extracellular K   [mM]
                const double nao        = 140.0;                                // 140 Extracellular Na  [mM]
                const double ek         = (1.0 * rtf) * log(Ko / Ki);           // [mV]
                
                double *par_sk;
                double gsk, kdsk;

                // Buffering parameters
                // Troponin C dynamic buffering current ITCi and ITCs
                const double BT         =   70.0;           // Tnc low affinity
                const double kon        =   0.0327;
                const double koff       =   0.0196;

                #ifdef ___EGTA_
                        const double BEGTA      = 350.0;
                        const double konEGTA    =   4E-3;
                        const double koffEGTA   =   2E-3;
                #endif

                // Buffer parameters from Grandi-Bers model
                const double Bmax_myosin        = 140e-3 * 1e3;         // [uM]

                const double koff_myoca         = 0.46e-3;              // [1/ms]
                const double kon_myoca          = 13.8 * 1e-3;          // [1/uM/ms]
                const double kd_myoca           = koff_myoca / kon_myoca;

                const double koff_myomg         = 0.057e-3;             // [1/ms]
                const double kon_myomg          = 0.0157 * 1e-3;        // [1/m [1/uM/ms]
                const double kd_myomg           = koff_myomg / kon_myomg;

                const double Mgi                = 1 * 1e3;              // [uM]

                const double Bmax_TnChigh       = 140e-3 * 1e3;         // [uM]

                const double koff_tnchca        = 0.032e-3;             // [1/ms]
                const double kon_tnchca         = 2.37 * 1e-3;          // [1/uM/ms]
                const double kd_tnchca          = koff_tnchca / kon_tnchca;

                const double koff_tnchmg        = 3.33e-3;              // [1/ms]
                const double kon_tnchmg         = 3e-3 * 1e-3;          // [1/uM/ms]
                const double kd_tnchmg          = koff_tnchmg / kon_tnchmg;

                const double Bmax_sr            = 47e-3 * 1e3;          // [mM]

                const double koff_sr            = 60e-3;                // [1/ms]
                const double kon_sr             = 100 * 1e-3;           // [1/uM/ms]
                const double kd_sr              = koff_sr / kon_sr;

                const double Bmax_cam           = 24e-3 * 1e3;          // [uM]

                const double koff_cam           = 238e-3;               // [1/ms]
                const double kon_cam            = 34 * 1e-3;            //[1/uM/ms]
                const double kd_cam             = koff_cam / kon_cam;

                const double Bmax_SLlowsl       = 37.4e-3 * 1e3;        // [uM]
                const double Bmax_SLlowj        = 4.6e-3 * 1e3;         // [uM]
                const double koff_sll           = 1300e-3;              // [1/ms]
                const double kon_sll            = 100 * 1e-3;           // [1/uM/ms]
                const double kd_sll             = koff_sll / kon_sll;

                const double Bmax_SLhighsl      = 13.4e-3 * 1e3;        // [uM]
                const double Bmax_SLhighj       = 1.65e-3 * 1e3;        // [uM]
                const double koff_slh           = 30e-3;                // [1/ms]
                const double kon_slh            = 100 * 1e-3;           // [1/m [1/uM/ms]
                const double kd_slh             = koff_slh / kon_slh;

                // RyR gating realted
                double static const  Bmax_Csqn  = 400;                  // [uM]
                double  BCSQN                   = Bmax_Csqn;
                double  koff_csqn               = 65;                   // [1/ms]
                double  kon_csqn                = 100 * 1e-3;           // [1/uM/ms]
                double  kd_csqn                 = koff_csqn / kon_csqn;
                double  Kc                      = kd_csqn;

                double  nM                      = 15;
                double  nD                      = 35;
                double  rhoinf                  = 5000;
                double  hh                      = 23.0;
                double  KK                      = 900;

                // Na buffering parameters (are not used in the code)
                const double Bmax_Naj           = 7.561;                // [mM]
                const double Bmax_Nasl          = 1.65;                 // [mM]
                const double koff_na            = 1e-3;                 // [1/ms]
                const double kon_na             = 0.1e-3;               // [1/mM/ms]
                const double kd_na              = koff_na / kon_na;

                // Na transport parameters (are not used in the code)
                double GNaB               = 1 * 0.597e-3 / 8113.0 * 0.96245 * (1e-10);  // [mS/uF*F]
                double IbarNaK            = 1 * 1.26 / 8113.0 * 0.96245 * (1e-10);      // [A]
                const double KmKo               = 1.5;                                  // [mM]
                const double Q10NaK             = 1.63;
                const double Q10KmNai           = 1.39;

                #ifdef ___EGTA
                        const double konEGTA        =   4E-3;
                        const double koffEGTA       =   2E-3;

                        double *caEGTAi, *caEGTAs;
                        double setBEGTA(double newBEGTA) {BEGTA = newBEGTA; resetBuffer(); return BEGTA;}
                        double getBEGTA(void) {return BEGTA;}
                #endif

                double NCXalpha;
                double Fjunc;
                double ena_junc, ena_sl;

                // // Ina parameters
                double GNa_hh           = 9.0 / 8113.0 * 0.96245 * (1e-10);     // [mS/uF*F]
                double GNaL_hh          = 0.0025 / 8113.0 * 0.96245 * (1e-10);  // [mS/uF*F]

                //  Na/K pump current
                double KmNaip           = 11.0;                                 // [mM]
                double sigma_nak        = (exp(nao / 67.3) - 1) / 7.0;
                double fnak;

                //cell parameters
                const double vjsr       = 0.02 * 0.75;
                const double vp_ave     = 0.00126;                              // [um^3]
                double Cmem             = 1.1 * (1e-10);                        // [F]

                int caffeine_app = 0;
};

#endif