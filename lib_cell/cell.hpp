#ifndef ___CELL_H
#define ___CELL_H

using namespace std;
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
#include <iostream>
#include "input_output.h"

class CCell {
        private:
                void pacex(double stim=0);
                static const int N = 101;
                static const double vc;
                static const double stim;
                static const double stimduration;
        
                double comp_ito(void);
                double comp_ikr(void);
                double comp_iks(void);
                double comp_ikp(void);
                double comp_ikach(void);
                double comp_ik1(void);
                double comp_ikur(void);
                double comp_I_Clbk(void);
                
                double ek;            // [mV]
                double ecl;           // [mV]
                double dt, ddt;
                double vold;

        public:
                void pace(double stim=0);
                double setdt(double dtt){dt=dtt; return dt;}
                double getdt(void){return dt;}
                CCell();
                virtual ~CCell();
                CCell& operator=(const CCell& cell);
                double *y;
                double *ydot;
                double v, ci, cnsr;

                double ica_ave, incx_ave, I_cabk_ave, I_pca_ave;
                void set_ica(double ica_diff) {ica_ave          = ica_diff;}
                void set_incx(double incx_diff) {incx_ave       = incx_diff;}
                void set_icabk(double I_cabk_diff) {I_cabk_ave  = I_cabk_diff;}
                void set_ipca(double I_pca_diff) {I_pca_ave     = I_pca_diff;}

                double incx_junc_diff, incx_junc_ave;
                double incx_sl_diff, incx_sl_ave;
                void set_incx_junc(double incx_junc_diff) {incx_junc_ave= incx_junc_diff;}
                void set_incx_sl(double incx_sl_diff) {incx_sl_ave      = incx_sl_diff;}

                double I_CaNa_junc_3d, I_CaK_junc_3d;
                void set_I_CaNa_junc(double new_I_CaNa_junc) {I_CaNa_junc_3d = new_I_CaNa_junc;}
                void set_I_CaK_junc(double new_I_CaK_junc) {I_CaK_junc_3d = new_I_CaK_junc;}

                double I_ClCa_junc_3d, I_ClCa_sl_3d;
                void set_I_ClCa_junc(double new_I_ClCa_junc) {I_ClCa_junc_3d = new_I_ClCa_junc;}
                void set_I_ClCa_sl(double new_I_ClCa_sl) {I_ClCa_sl_3d = new_I_ClCa_sl;}

                double I_NaK_junc_3d, I_NaK_sl_3d;
                void set_I_NaK_junc(double new_I_NaK_junc) {I_NaK_junc_3d = new_I_NaK_junc;}
                void set_I_NaK_sl(double new_I_NaK_sl) {I_NaK_sl_3d = new_I_NaK_sl;}

                double I_Nabk_junc_3d, I_Nabk_sl_3d;
                void set_I_Nabk_junc(double new_I_Nabk_junc) {I_Nabk_junc_3d = new_I_Nabk_junc;}
                void set_I_Nabk_sl(double new_I_Nabk_sl) {I_Nabk_sl_3d = new_I_Nabk_sl;}

                double I_Na_junc_3d, I_Na_sl_3d;
                void set_I_Na_junc(double new_I_Na_junc) {I_Na_junc_3d = new_I_Na_junc;}
                void set_I_Na_sl(double new_I_Na_sl) {I_Na_sl_3d = new_I_Na_sl;}
                void set_vold(double v_set) {vold = v_set;}

                double I_SK_junc_3d, I_SK_sl_3d;
                void set_I_SK_junc(double new_I_SK_junc ) {I_SK_junc_3d = new_I_SK_junc;}
                void set_I_SK_sl(double new_I_SK_sl ) {I_SK_sl_3d = new_I_SK_sl;}
                
                void get_steady_state();
                void set_steady_state();

                // conditions (just keep at the baseline for all the simulation)
                double ISO;
                double flagMina;
                double drug_index;
                double drug;
                double AF;
                double RA;
                double epi;
                double bGal;
                double EAD;
                double Ach;

                double _I_to;
                double _I_kr;
                double _I_ks;
                double _I_ki;
                double _I_kp;
                double _I_kur;
                double _I_KAch;
                double _I_CaK;

                double _I_Catot;
                double _I_ncx;

                double _I_Clbk;

                double gnabar;        // Na Current conductance
                double gnabbar;       // Na Background Current conductance
                double gnakbar;       // Na/K Pump Current conductance
                double gkrbar;        // Rapidly Activating K Current conductance
                double gksbar;        // Slowly Activating K Current conductance
                double gkpbar;        // Plateau K Current conductance
                double gkachbar;      // Muscarinic-Receptor-Activated K Current conductance
                double gtobar;        // Transient Outward K Current conductance
                double gkurbar;       // Ultra rapid delayed rectifier Outward K Current conductance
                double gkibar;        // Time-Independent K Current conductance
                double vupbar;        // SERCA strength
                double gcabbar;       // Ca Background conductance
                double gpcabar;       // Sarcolemmal Ca Pump Current conductance
                double gncxbar;       // Na/Ca Exchanger flux conductance
                double gcabar;        // L-type Calcium Current conductance
};

#endif
