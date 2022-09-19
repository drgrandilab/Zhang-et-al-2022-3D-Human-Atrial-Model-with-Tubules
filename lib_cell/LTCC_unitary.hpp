#ifndef LTCC_UNITARY_HPP
#define LTCC_UNITARY_HPP

class LTCC_unitary
{
public:
	LTCC_unitary(double v) {
                const double pca = 1.4 * 6.1 * 0.65;
		z = v * F / (R * T);
		za = z * 2.0;
		//L- Ca current (ca independent part)
		factor1 = 4.0 * pca * F * F / (R * T);
		factor = v * factor1;
	};

	~LTCC_unitary() {};
        
	double compute_LTCC_unitary(double Caj, double cao) {
		double ica = (factor / (exp(za) - 1.0)) * (Caj * 0.001 * exp(za) - cao);
		if (fabs(za) < 0.1)
		{
                        ica = (factor1 / (2.0 * F / (R * T))) * (Caj * 0.001 * exp(za) - cao); // [uM/ms]
		}
		if (ica > 0) ica = 0;
		return ica;
	}

	double z, za, factor1, factor;

	const double F = 96.5;         // [C/mmol]
	const double R = 8.314;        // [J/mol*K]
	const double T = 310;          // [K]
};

#endif