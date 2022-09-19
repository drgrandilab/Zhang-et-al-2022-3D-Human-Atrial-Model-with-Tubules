#ifndef LTCC_Markov_HPP
#define LTCC_Markov_HPP

#include <fstream>
#include <vector>
#include "xor_rand.hpp"

class LTCC_Markov
{
public:
	LTCC_Markov(int NCaL_set = 4, int ID = 0);
	~LTCC_Markov();

	std::vector<int> LTCC_State;
	xor_rand random_num_LTCC;

        // LTCC Current - Fixed Parameters
        double static const cpt     = 3.75e-3 * 0.55;     // [mM]
        double static const cat     = 7.617e-3;           // [mM]
        double static const s1o     = 4.2 * 0.0182688;    // [1/ms]
        double static const k1o     = 4.2 * 0.024168;     // [1/ms]
        double static const k2o     = 0.000103615;  // [1/ms]
        double static const sp0     = 1.5;
        double static const sp1     = 3;            // [ms]
        double static const sp2     = 27;           // [mV]
        double static const sp3     = 3;            // [mV]
        double static const sp4     = 4;            // [mV]
        double static const sp5     = 7.1 ;         // [mV]
        double static const sp6     = 15.6;         // [mV]
        double static const sp7     = 10;           // [ms]
        double static const sp8     = 4954;         // [ms]
        double static const sp9     = 78.0329;      // [ms]
        double static const sp10    = 0.1;          // [ms]
        double static const aR2     = 1;
        double static const sR2     = -9;           // [mV]
        double static const pR2     = 1.0 / 6.0;    //[1/mV]
        double static const aT2     = 1;            // [1/ms]
        double static const sT2     = -1000;        // [mV]
        double static const pT2     = 0.100;        // [1/mV]
        double static const aR1     = 0.09091;
        double static const sR1     = -1000;        // [mV]
        double static const pR1     = 0.100;        // [1/mV]
        double static const aT1     = 0.30303;      // [1/ms]
        double static const sT1     = -1000;        // [mV]
        double static const pT1     = 0.100;        // [1/mV]
        double static const aRv2    = 0.9;
        double static const sRv2    = -29;          // [mV]
        double static const pRv2    = 0.135;        // [1/mV]
        double static const aTv2    = 500;          // [1/ms]
        double static const sTv2    = -25;          // [mV]
        double static const pTv2    = 0.050;        // [1/mV]
        double static const aRv1    = 0.85;
        double static const sRv1    = -27;          // [mV]
        double static const pRv1    = 0.090;        // [1/mV]
        double static const aTv1    = 270;          // [1/ms]
        double static const sTv1    = -180;         // [mV]
        double static const pTv1    = 0.100;        // [1/mV]

	int NCaL;
	int open_NCaL_num;
	double Po_LCCj_m1;
	int flag_7_state, ICa_speed;

	double update_states_v7(double dt, double Caj, double Vm);
	void print_to_file(double t, std::ofstream & output_file);

        // Constants
	const static double R    = 8314;    // [J/kmol*K]
	const static double Frdy = 96485;   // [C/mol]
	const static double Temp = 310;     // [K]
	const static double FoRT = Frdy / R / Temp;
	const static double Qpow = (Temp - 310) / 10;

        // Fixed ion concentrations
	double static const Ko  = 5.4;         // Extracellular K   [mM]
	double static const Nao = 140;         // Extracellular Na  [mM]
	double static const Cao = 1.8;         // Extracellular Ca  [mM]
};
#endif
