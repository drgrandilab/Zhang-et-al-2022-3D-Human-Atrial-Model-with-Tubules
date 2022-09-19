#include "subcell.hpp"

void CSubcell::update_LTCC_RyR(double voltage, double local_dt) {
        v       = voltage;
        fnak    = 1 / (1 + 0.1245 * exp(-0.1 * v / rtf) + 0.0365 * sigma_nak * exp(-v / rtf) );
        
        double gsk_vm           = par_sk[2]/(1+exp((v-ek+par_sk[3])*par_sk[4])) + par_sk[5]/(1+exp((-(v-ek+par_sk[6])*par_sk[7])));
        gsk_vm                  = gsk_vm * gsk * (v - ek);

        LTCC_unitary  ICaL_unitary(voltage);
        inaca NCX(nap[0], nao, voltage, cao);

        double sum_inak_junc        = 0;
        double sum_inabk_junc       = 0;
        double sum_ina_junc         = 0;
        double sumica               = 0;
        double sumicana             = 0;
        double sumicak              = 0;
        double sum_jnaca_j_flux     = 0;
        double sum_iclca_j          = 0;
        double sum_jnaca_sl         = 0;
        double sumjcabk             = 0;
        double sumjslcap            = 0;
        double sum_iclca_sl         = 0;
        double sum_inak_sl          = 0;
        double sum_inabk_sl         = 0;
        double sum_ina_sl           = 0;
        double sumir                = 0;
        int total_ICaL_open         = 0;
        double sum_isk_junc         = 0;
        double sum_isk_sl           = 0;

        INafast_gating   = Ina.Update_INa(local_dt, v);

        #pragma omp parallel
        {
                #pragma omp for schedule(auto) reduction(+: sum_inak_junc, sum_inabk_junc, sum_ina_junc, sumica, sumicana, sumicak, sum_jnaca_j_flux, sum_iclca_j, \
                sum_jnaca_sl, sumjcabk, sumjslcap, sum_iclca_sl, sum_inak_sl, sum_inabk_sl, sum_ina_sl, total_ICaL_open, sumir, sum_isk_junc, sum_isk_sl)
                #pragma ivdep
                #pragma vector always
                for (int index = 0; index < num_tubule + n; index++) 
                {
                        if (index < num_tubule) 
                        {
                                int NL  = 0;
                                int id  = Tubular_map[index];

                                NL      = LTCC_Markov_vec[id].update_states_v7(local_dt, 0.001 * cp[id], voltage);      // [A/F]

                                ica_array[id]   = gca * ICaL_unitary.compute_LTCC_unitary(cp[id], cao) * NL;
                                ICaNa_array[id] = update_ICaNa(nap[id], NL);      // [A]
                                ICaK_array[id]  = update_ICaK(NL);                // [A]

                                LTCC_Open_num[id]       = NL;
                                total_ICaL_open         += tubule_flag[id] * LTCC_Open_num[id];

                                // ncx in cleft area
                                ncx_array_p[id] = map_ncx[id] * tubule_flag[id] * NCXalpha * NCX.compute_NCX_soltis(nap[id], nao, cp[id], cao)  * vs / vp[id] * NCX_fit_scale;
                                ncx_array[id]   = map_ncx[id] * tubule_flag[id] * (1 - NCXalpha) * NCX.compute_NCX_soltis(nas[id], nao, cs[id], cao) * NCX_fit_scale;
                                I_ClCa_junc_array[id]   = Fjunc * (GClCa) / (1 + KdClCa / cp[id]) * (v - ecl);   // [A]

                                double ena_junc         = (1 * rtf) * log(nao / nap[id]); // [mV]
                                I_nak_junc_array[id]    = (1 * Fjunc * (IbarNaK) * fnak * Ko / (1 + std::pow(KmNaip / (nap[id]), 4.0) ) / (Ko + KmKo) ); // [A]
                                I_nabk_junc_array[id]   = (Fjunc * (GNaB) * (v - ena_junc) ); // [A]
                                I_Na_junc_array[id]     = Fjunc * GNa_hh * INafast_gating * (v - ena_junc); // [A]
                                
                                I_ClCa_sl_array[id]     = tubule_flag[id] * (1 - Fjunc) * GClCa / (1 + KdClCa / cs[id]) * (v - ecl);      // [A]

                                double eca      = rtf2 * log(cao * 1000 / cs[id]); // [mV]
                                icabk_array[id] = tubule_flag[id] * gcabk * (v - eca) * 1.1; 
                                const double vmax       = 2.2 * 0.01;
                                const double km         = 0.5;
                                const double h          = 1.6;
                                
                                jpca_array[id]  = tubule_flag[id] * qslcap * vmax / (1 + std::pow(km / cs[id], h)) *  1.1;

                                I_nak_sl_array[id]         = tubule_flag[id] * (1 * (1 - Fjunc) * (IbarNaK) * fnak * Ko / (1 + std::pow(KmNaip / (nas[id]), 4.0) ) / (Ko + KmKo) ); // [A]
                                
                                double ena_sl           = (1 * rtf) * log(nao / nas[id]);         // [mV]
                                I_nabk_sl_array[id]        = tubule_flag[id] * (1 - Fjunc) * GNaB * (v - ena_sl); // [A]
                                I_Na_sl_array[id]          = tubule_flag[id] * (1 - Fjunc) * GNa_hh * INafast_gating * (v - ena_sl); // [A]
                                
                                double gsk_ca_junc      = 1.0 / (1.0 + exp((log10(kdsk)-log10(cp[id] * 1e-3))/0.3));
                                I_SK_junc_array[id]     = tubule_flag[id] * Fjunc * gsk_ca_junc * gsk_vm; // [A]

                                double gsk_ca_sl        = 1.0 / (1.0 + exp((log10(kdsk)-log10(cs[id] * 1e-3))/0.3));
                                I_SK_sl_array[id]       = tubule_flag[id] * (1.0 - Fjunc) * gsk_ca_sl * gsk_vm; // [A]
                                
                                sum_inak_junc    += I_nak_junc_array[id];       // [A*#]
                                sum_inabk_junc   += I_nabk_junc_array[id];      // [A*#]
                                sum_ina_junc     += I_Na_junc_array[id];        // [A*#]
                                sumica           += ica_array[id] * vp[id];
                                sumicana         += ICaNa_array[id];            // [A*#]
                                sumicak          += ICaK_array[id];             // [A*#]
                                sum_jnaca_j_flux += ncx_array_p[id] * vp[id];
                                sum_iclca_j      += I_ClCa_junc_array[id];      // [A*#]
                                sum_jnaca_sl     += ncx_array[id];
                                sumjcabk         += icabk_array[id];
                                sumjslcap        += jpca_array[id];
                                sum_iclca_sl     += I_ClCa_sl_array[id];        // [A*#]
                                sum_inak_sl      += I_nak_sl_array[id];         // [A*#]
                                sum_inabk_sl     += I_nabk_sl_array[id];        // [A*#]
                                sum_ina_sl       += I_Na_sl_array[id];          // [A*#]

                                sum_isk_junc     += I_SK_junc_array[id];        // [A*#]
                                sum_isk_sl       += I_SK_sl_array[id];          // [A*#]
                        } else 
                        {
                                int id = index - num_tubule;
                                if (!caffeine_app) 
                                {
                                        PoSpark[id]  = RyR_vec[id].Update_RyR_stochastic(local_dt, cp[id], cjsr[id]);
                                } else 
                                {
                                        PoSpark[id]     = 0.2;
                                }

                                JSR_buffer_factors[id] = get_jSR_inst_buffering_simp(cjsr[id]);

                                double Po       = PoSpark[id];
                                double Ir       = RyR_vec[id].Jmax * Po * (cjsr[id] - cp[id]) / vp[id];
                                sumir           += Ir;
                                j_ryr[id]       = Ir; // [uM/ms]
                        }
                }
        }
        icabkave        = sumjcabk;
        icaave          = sumica;
        islcapave       = sumjslcap;
        incx_sl_ave     = sum_jnaca_j_flux;
        incx_junc_ave   = sum_jnaca_sl * vs;
        incxave         = incx_junc_ave + incx_sl_ave;

        i_clca_junc_ave = sum_iclca_j;
        i_clca_sl_ave   = sum_iclca_sl;
        i_cana_ave      = sumicana;
        i_cak_ave       = sumicak;

        i_nak_junc_ave  = sum_inak_junc;
        i_nak_sl_ave    = sum_inak_sl;
        i_nabk_junc_ave = sum_inabk_junc;
        i_nabk_sl_ave   = sum_inabk_sl;
        i_na_junc_ave   = sum_ina_junc;
        i_na_sl_ave     = sum_ina_sl;
        
        i_sk_junc_ave   = sum_isk_junc;
        i_sk_sl_ave     = sum_isk_sl;   

        num_open_ICaL = total_ICaL_open;

        irave           = sumir;
}

void CSubcell::pace2(double voltage, double local_dt)
{
        v       = voltage;

        #ifndef ___NO_DIFFUSION
                //set diffusion terms
                computeIci();           // diffusion ci
                computeIcnsr();         // diffusion cnsr

                #ifdef ___NO_CS_BUFFER
                        computecsmn();  // diffusion cs
                #else
                        computeIcs();   // diffusion cs
                #endif
        #endif

        double sumjup           = 0;
        double sumjleak         = 0;

        #pragma omp parallel
        {
                // ------------Calculate cleft and jSR calcium concentration--------------------------------
                #pragma omp for schedule(auto)
                #pragma ivdep
                #pragma vector always
                for (int id = 0; id < n; id++)
                {
                        double Ir       = j_ryr[id];

                        Itr[crupos[id]] = (cnsr[crupos[id]] - cjsr[id]) / tautr;
                        double tmp      = JSR_buffer_factors[id] * (Itr[crupos[id]] - Ir * (vp[id] / vjsr));
                        cjsr[id]        =  cjsr[id] + local_dt *  tmp;

                        Idps[crupos[id]]        = (cp[id] - cs[crupos[id]]) / taup;

                        // Cleft buffers:  Cam, SR, SLL, SLH
                        double Icam_on          = kon_cam * ( Bmax_cam - cleft_cam[id] );
                        double Isr_on           = kon_sr * ( Bmax_sr - cleft_SR[id] );
                        double ISLL_cleft_on    = kon_sll * ( Bmax_SLlowj - cleft_SLL[id] );
                        double ISLH_cleft_on    = kon_slh * ( Bmax_SLhighj - cleft_SLH[id] );

                        double Icam_off         = koff_cam * cleft_cam[id];
                        double Isr_off          = koff_sr * cleft_SR[id];
                        double ISLL_cleft_off   = koff_sll * cleft_SLL[id];
                        double ISLH_cleft_off   = koff_slh * cleft_SLH[id];

                        // Update buffers
                        double Icam_cleft       = cp[id] * Icam_on - Icam_off;
                        double Isr_cleft        = cp[id] * Isr_on - Isr_off;
                        double ISLL_cleft       = cp[id] * ISLL_cleft_on - ISLL_cleft_off;
                        double ISLH_cleft       = cp[id] * ISLH_cleft_on - ISLH_cleft_off;

                        cleft_cam[id]   += Icam_cleft * local_dt;
                        cleft_SR[id]    += Isr_cleft * local_dt;
                        cleft_SLL[id]   += ISLL_cleft * local_dt;
                        cleft_SLH[id]   += ISLH_cleft * local_dt;

                        double Cleft_buffer     = Icam_cleft + Isr_cleft + ISLL_cleft + ISLH_cleft;

                        double Ica      = tubule_flag[id] * ica_array[id];;
                        double ICaNa    = tubule_flag[id] * ICaNa_array[id];
                        double ICaK     = tubule_flag[id] * ICaK_array[id];
                        double jnaca    = tubule_flag[id] * ncx_array_p[id];

                        double dcp      = Ir - Ica + jnaca - Idps[crupos[id]] - Cleft_buffer;

                        cp[id]          += dcp * local_dt;

                        if (cp[id] <= 0.00001) cp[id] = 0.00001;

                        double I_nak_junc       = tubule_flag[id] * I_nak_junc_array[id];
                        double I_nabk_junc      = tubule_flag[id] * I_nabk_junc_array[id];
                        double I_Na_junc        = tubule_flag[id] * I_Na_junc_array[id];
                }

                // -------------Calculate submembrane, cytosol and nSR calcium concentration----------------
                #pragma omp for reduction(+: sumjup, sumjleak) schedule(auto) 
                #pragma ivdep
                #pragma vector always
                for (int id = 0; id < nn; id++)
                {
                        double cs_i = cs[id];
                        double ci_i = ci[id];

                        double Iup = 0;
                        if (!caffeine_app) {
                                const double H  = 1.787;
                                double tmp1     = pow(ci_i / kup, H);
                                double tmp2     = pow(cnsr[id] / KNSR, H);
                                Iup     = map_serca[id] * vup * (tmp1 - tmp2) / (1 + tmp1 + tmp2);
                        }

                        double jnaca    = tubule_flag[id] * ncx_array[id]; // output funciton
                        double jcabk    = tubule_flag[id] * icabk_array[id];
                        double jslcap   = tubule_flag[id] * jpca_array[id];

                        // calculate EGTA buffer changes
                        #ifdef ___EGTA
                                double IEGTAs   = konEGTA * cs[id] * (BEGTA - caEGTAs[id]) - koffEGTA * caEGTAs[id];
                                double IEGTAi   = konEGTA * ci_i * (BEGTA - caEGTAi[id]) - koffEGTA * caEGTAi[id];
                        #endif

                        // Submembrane buffers:  EGTA (optional, dynamic), Cam, SR, SLL, SLH
                        double Icam_submem      = kon_cam * cs_i * ( Bmax_cam - SL_cam[id] ) - koff_cam * SL_cam[id];
                        double Isr_submem       = kon_sr * cs_i * ( Bmax_sr - SL_SR[id] ) - koff_sr * SL_SR[id];
                        double ISLL_submem      = kon_sll * cs_i * ( Bmax_SLlowsl - SL_SLL[id] ) - koff_sll * SL_SLL[id];
                        double ISLH_submem      = kon_slh * cs_i * ( Bmax_SLhighsl - SL_SLH[id] ) - koff_slh * SL_SLH[id];
                        SL_cam[id]      += Icam_submem * local_dt;
                        SL_SR[id]       += Isr_submem * local_dt;
                        SL_SLL[id]      += ISLL_submem * local_dt;
                        SL_SLH[id]      += ISLH_submem * local_dt;
                        double Submem_buffer   = Icam_submem + Isr_submem + ISLL_submem + ISLH_submem;

                        //Diffusion from submembrane to myoplasm Idsi
                        double Idsi     = (cs_i - ci_i) / tausi;

                        #ifdef ___DEBUG
                                if (isnan(Idsi)) //(Idsi != Idsi)
                                {
                                        cout << setprecision(10) << id << "\t" << Idsi << "\t cs=" << cs_i << "\t ci=" << ci_i << endl;
                                        bSTOP = true;
                                }
                        #endif

                        #ifdef ___EGTA
                                double dcs      = Idps[id] * vp[id] / vs + jnaca - Idsi - IEGTAs + Ics[id] - jcabk - jslcap - Submem_buffer;
                        #else
                                double dcs      = Idps[id] * vp[id] / vs + jnaca - Idsi + Ics[id] - jcabk - jslcap - Submem_buffer;
                        #endif

                        cs[id]  += dcs * local_dt;

                        // Cytosolic buffers:  EGTA (optional, dynamic), MCa, MMg, SR, Cam, Troponin C binding (dynamic), TnChc, Cyto_TnCHm
                        double ITCi     =   kon * ci_i * ( BT - cati[id] ) - koff * cati[id];
                        double IMCa_Cyto        =   kon_myoca * ci_i * ( Bmax_myosin - Cyto_MCa[id] - Cyto_MMg[id] ) - koff_myoca * Cyto_MCa[id];       // Myosin_ca [mM/ms] -> [uM/ms]
                        double IMMg_Cyto        =   kon_myomg * Mgi * ( Bmax_myosin - Cyto_MCa[id] - Cyto_MMg[id] ) - koff_myomg * Cyto_MMg[id];        // Myosin_mg [mM/ms] -> [uM/ms]
                        double Isr_Cyto         =   kon_sr * ci_i * ( Bmax_sr - Cyto_SR[id] ) - koff_sr * Cyto_SR[id];
                        double Icam_Cyto        =   kon_cam * ci_i * ( Bmax_cam - Cyto_cam[id] ) - koff_cam * Cyto_cam[id];
                        double ITnCHc_Cyto      =   kon_tnchca * ci_i * ( Bmax_TnChigh - Cyto_TnCHc[id] - Cyto_TnCHm[id] ) - koff_tnchca * Cyto_TnCHc[id];      // Cyto_TnCHc     [mM/ms] -> [uM/ms]
                        double ITnCHm_Cyto      =   kon_tnchmg * Mgi * ( Bmax_TnChigh - Cyto_TnCHc[id] - Cyto_TnCHm[id] ) - koff_tnchmg * Cyto_TnCHm[id];       // Cyto_TnCHm     [mM/ms] -> [uM/ms]
                        double Cyto_buffer      =   IMCa_Cyto + IMMg_Cyto + Isr_Cyto + Icam_Cyto + ITCi + ITnCHc_Cyto + ITnCHm_Cyto;
                        
                        cati[id]        +=  ITCi * local_dt;
                        Cyto_MCa[id]    +=  IMCa_Cyto * local_dt;
                        Cyto_MMg[id]    +=  IMMg_Cyto * local_dt;
                        Cyto_SR[id]     +=  Isr_Cyto * local_dt;
                        Cyto_cam[id]    +=  Icam_Cyto * local_dt;
                        Cyto_TnCHc[id]  +=  ITnCHc_Cyto * local_dt;
                        Cyto_TnCHm[id]  +=  ITnCHm_Cyto * local_dt;

                        const double KJSR       = 500;
                        double cjsr2    = cnsr[id] * cnsr[id];
                        double Ileak    = gleak * (cjsr2 * (cnsr[id] - ci_i)) / (cjsr2 + KJSR * KJSR);

                        #ifdef ___EGTA
                                double dci      = Idsi * (vs / vi) - Iup + Ileak + Ici[id] - IEGTAi - Cyto_buffer;
                        #else
                                double dci      = Idsi * (vs / vi) - Iup + Ileak + Ici[id] - Cyto_buffer;
                        #endif

                        // Update calcium concentration
                        ci[id]          +=  dci * local_dt;
                        double dcnsr    = ((Iup - Ileak) * (vi / vnsr) - Itr[id] * (vjsr / vnsr) + Icnsr[id]);
                        cnsr[id]        += dcnsr * local_dt;
                        sumjup          += Iup;
                        sumjleak        += Ileak;

                        // update EGTA buffer in submembrane and cytosol
                        #ifdef ___EGTA
                                caEGTAi[id]     += IEGTAi * local_dt;
                                caEGTAs[id]     += IEGTAs * local_dt;
                        #endif
                }
        }
        iupave          = sumjup;
        ileakave        = sumjleak;
}
