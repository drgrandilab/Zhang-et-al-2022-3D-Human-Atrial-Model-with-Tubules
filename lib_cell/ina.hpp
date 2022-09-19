#ifndef INA_HPP
#define INA_HPP

#include <math.h>

class ina
{
public:
	ina();

	double Update_INa(double dt, double v);
	~ina() {};

        double state1_ina, state2_ina, state3_ina;

        double am, bm;
        double ah, bh;
        double aj, bj;

        double gating;
};
#endif