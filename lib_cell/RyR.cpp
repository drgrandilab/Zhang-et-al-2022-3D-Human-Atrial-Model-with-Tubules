#include "RyR.hpp"


RyR::RyR(int RYR_Num, int ID, double Jmax_set)
	: random_num_RyR(0, ID), Jmax(Jmax_set)
{

	N_RyR = RYR_Num;
	RyR_1 = 0 + int (5.0 * random_num_RyR.gen_rand_uint() / (double)(UINT_MAX));
	RyR_2 = 0;
	RyR_3 = 0;
}

RyR::~RyR() {}

double RyR::Update_RyR_stochastic(double dt, double Caj, double CaSRj) {
	double rho = rhoinf * 1.0 / (1.0 + pow(KK / CaSRj, hh));
	if (rho < 0.0000001)rho = 0.0000001;
	double MM = (sqrt(1 + 8 * rho * BCSQN) - 1) / (4 * rho * BCSQN);
	
        Po = (RyR_2 + RyR_3) / 100.0;

	double cp2 = Caj * Caj;
	double sgmd = cp2 * Caj / (Kcp * Kcp * Kcp + cp2 * Caj);
	double k12 = Ku * sgmd + pedk12;
	double k43 = Kb * sgmd + pedk43;

	double k14 = MM / taub * BCSQN / BCSQN0;
	double k21 = 1 / tauc1;
	double k23 = MM / taub * BCSQN / BCSQN0;
	double k41 = 1 / tauu;
	double k34 = 1 / tauc2;
	double k32 = k41 * k12 * k23 * k34 / k43 / k21 / k14;

	int RyR_4 = N_RyR - RyR_1 - RyR_2 - RyR_3;
	int ryr12 = Update_RyR_rates_stochastic(RyR_1, k12 * dt);
	int ryr14 = Update_RyR_rates_stochastic(RyR_1, k14 * dt);
	int ryr21 = Update_RyR_rates_stochastic(RyR_2, k21 * dt);
	int ryr23 = Update_RyR_rates_stochastic(RyR_2, k23 * dt);
	int ryr43 = Update_RyR_rates_stochastic(RyR_4, k43 * dt);
	int ryr41 = Update_RyR_rates_stochastic(RyR_4, k41 * dt);
	int ryr34 = Update_RyR_rates_stochastic(RyR_3, k34 * dt);
	int ryr32 = Update_RyR_rates_stochastic(RyR_3, k32 * dt);
	RyR_1 = RyR_1 - (ryr12 + ryr14) + (ryr21 + ryr41);
	RyR_2 = RyR_2 - (ryr21 + ryr23) + (ryr12 + ryr32);
	RyR_3 = RyR_3 - (ryr34 + ryr32) + (ryr43 + ryr23);

	if (RyR_1 < 0 || RyR_2 < 0 || RyR_3 < 0 || RyR_1 + RyR_2 + RyR_3 > N_RyR)
	{
		if (RyR_1 < 0)
		{
			if (random_num_RyR.gen_rand_uint() % 2)
				RyR_2 += RyR_1;
			RyR_1 = 0;
		}
		if (RyR_2 < 0)
		{
			if (random_num_RyR.gen_rand_uint() % 2)
			{
				RyR_1 += RyR_2;
				if (RyR_1 < 0)RyR_1 = 0;
			}
			else
				RyR_3 += RyR_2;
			RyR_2 = 0;
		}
		if (RyR_3 < 0)
		{
			if (random_num_RyR.gen_rand_uint() % 2)
			{
				RyR_2 += RyR_3;
				if (RyR_2 < 0)RyR_2 = 0;
			}
			RyR_3 = 0;
		}
		if (RyR_1 + RyR_2 + RyR_3 > N_RyR)
		{
			RyR_4 = N_RyR - (RyR_1 + RyR_2 + RyR_3);
			if (random_num_RyR.gen_rand_uint() % 2)
			{
				RyR_3 += RyR_4;
				if (RyR_3 < 0)
				{
					RyR_3 -= RyR_4;
					RyR_1 += RyR_4;
					if (RyR_1 < 0)
					{
						RyR_1 -= RyR_4;
						RyR_2 += RyR_4;
					}
				}
			}
			else
			{
				RyR_1 += RyR_4;
				if (RyR_1 < 0)
				{
					RyR_1 -= RyR_4;
					RyR_3 += RyR_4;
					if (RyR_3 < 0)
					{
						RyR_3 -= RyR_4;
						RyR_2 += RyR_4;
					}
				}
			}
		}
	}

	return Po;
}

// this function simulates the transition of RyR in the model/ stochastically .
int RyR::Update_RyR_rates_stochastic(double num, double p)
{
	int res;
	double lambda = num * p;

	if (lambda > 12)
	{
		//Gaussian
		double x1, x2, w;
		do
		{
			x1 = 2.0 * random_num_RyR.gen_rand() - 1.0;
			x2 = 2.0 * random_num_RyR.gen_rand() - 1.0;
			w = x1 * x1 + x2 * x2;
		} while (w >= 1.0);
		w = sqrt((-2.0 * log(w)) / w);
		double y1 = x1 * w;
		//double y2=x2*w;
		res = y1 * sqrt(num * p * (1 - p)) + num * p; // *** ave=num*p , rho^2=num*p*(1-p)
		res = int(res + 0.5); //round
	}
	else if (100 * p < 6.6 + 52 * pow(num, double(-0.5)))
	{
		//Poisson
		double L = exp(-lambda);
		double k = 0;
		double pp = 1;
		do
		{
			k++;
			double u = random_num_RyR.gen_rand();
			pp *= u;
		} while (pp >= L);
		res = k - 1;
	}
	else
	{
		//Gaussian
		double x1, x2, w;
		do
		{
			x1 = 2.0 * random_num_RyR.gen_rand() - 1.0;
			x2 = 2.0 * random_num_RyR.gen_rand() - 1.0;
			w = x1 * x1 + x2 * x2;
		} while (w >= 1.0);
		w = sqrt((-2.0 * log(w)) / w);
		double y1 = x1 * w;
		//double y2=x2*w;
		res = y1 * sqrt(num * p * (1 - p)) + num * p; // *** ave=num*p , rho^2=num*p*(1-p)
		res = int(res + 0.5); //round
	}

	if (res < 0)
		res = 0;

	return res;
}
