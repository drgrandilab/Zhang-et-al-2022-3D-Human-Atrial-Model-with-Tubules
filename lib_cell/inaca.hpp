#ifndef INACA_HPP
#define INACA_HPP

#include <math.h>

class inaca
{
public:
	inaca(double nai, double nao, double V, double cao);

        double compute_NCX_soltis(double nai,  double nao, double Caj, double cao);

        ~inaca() {};

	//NCX (ca independent part)
	const double Kmnai = 12.3;             // [mM]
	const double Kmcai = 3.59;             // [uM]
	const double ksat = 0.27;
	const double eta = 0.35;

	const static double F = 96.5;          // [C/mmol]
	const static double R = 8.314;         // [J/mol*K]
	const static double T = 310;           // [K]
	const static double rtf = R * T / F;
	const static double rtf2 = R * T / (2 * F);
	const static  double vnaca = 7.452 * 3 * 0.8;

	double t1 ;
	double x1a;
	double x1b;
	double x2;
};


#endif