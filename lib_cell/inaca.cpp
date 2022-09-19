#include "inaca.hpp"

inaca::inaca(double nai, double nao, double V, double cao) {
	double z = V * F / (R * T);
	t1 = (Kmcai * 0.001) * nao * nao * nao * (1 + (nai * nai * nai / (Kmnai * Kmnai * Kmnai)));
	x1a = exp(eta * z) * nai * nai * nai * cao;
	x1b = exp((eta - 1) * z) * nao * nao * nao;
	x2 = (1 + ksat * exp((eta - 1) * z));
}

double inaca::compute_NCX_soltis(double nai,  double nao, double Caj, double cao) {
        double Kmcao = 1.3;
        double Kmnao = 87.5;
        double Kda = 0.256;

        double x3 = vnaca;
        double Ka = 1 / (1 + (Kda * Kda * Kda / (Caj * Caj * Caj)));
        double t2 = Kmnao * Kmnao * Kmnao * (Caj * 0.001) +  Kmnai * Kmnai * Kmnai * cao * (1 + Caj / Kmcai);
        double t3 = Kmcao * nai * nai * nai + nai * nai * nai * cao + nao * nao * nao * (Caj * 0.001);
        return (Ka * x3 * (x1a - x1b * (Caj * 0.001)) / ((t1 + t2 + t3) * x2));
}