#ifndef RyR_HPP
#define RyR_HPP
#include "xor_rand.hpp"

class RyR
{
public:
	RyR(int RYR_Num = 41, int ID = 0, double Jmax_set = 0.0147 * 18 * 0.5);
	~RyR();

	int RyR_1, RyR_2, RyR_3;
	xor_rand random_num_RyR;

	double Po;
        double Jmax;
	int N_RyR;

	double Update_RyR_stochastic(double dt, double Caj, double CaSRj);
	int Update_RyR_rates_stochastic(double num, double p);
	void set_Jmax(double Jmax_set) {
		Jmax = Jmax_set;
	}

	double static const MaxSR      = 15;
	double static const MinSR      = 1;
	double static const ec50SR     = 450;
	double static const hkosrca    = 2.5;

	double static const Ku         = 5.0;
	double static const Kb         = 0.005;
        double static const tauu       = 125.0;
	double static const taub       = 0.5;
	double static const tauc1      = 2.0;
	double static const tauc2      = 0.3;

	double static const taup       = 0.022;
	double static const tautr      = 5.0;
	double static const BCSQN0     = 400;
        double static const Kcp        = 10.0;

	double static const pedk12     = 0.000001;
	double static const pedk43     = 0.000001;
        
        double static const  Bmax_Csqn          = 400;          // [uM]
        double static const  BCSQN              = Bmax_Csqn;
        double static const  koff_csqn          = 65;           // [1/ms]
        double static const  kon_csqn           = 100 * 1e-3;   // [1/uM/ms]
        double static const  kd_csqn            = koff_csqn / kon_csqn;
        double static const  Kc                 = kd_csqn;

        double static const  nM                 = 15;
        double static const  nD                 = 35;
        double static const  rhoinf             = 5000;
        double static const  hh                 = 23.0;
        double static const  KK                 = 900;
};

#endif
