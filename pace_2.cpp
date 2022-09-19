#include "Spatial_Cell.h"
#include "stimulus.h"

int main(int argc, char const *argv[])
{
        double dt      = 0.01;
        int count      = 0;

        double stim_current     = 12.5;
        double stim_duration    = 5;
        double pre_rest         = 4000;
        double bcl              = atof(argv[1]);
        int iter                = atoi(argv[2]);
        double delay_rest       = 5000;
        double T_total          = pre_rest + bcl * iter + delay_rest;

        Spatial_Cell cell(dt, 0, 0, bcl, argv[3]);

	for (double t = 0; t < T_total; t += dt)
	{
		cell.dt = dt;
		double stim   = 0;
                if ( t >= pre_rest && t < T_total - delay_rest )
                {
                        stim    = S1(pre_rest, stim_current, bcl, t, stim_duration);
                }

		std::string nnn;
		if (count % 100 == 0)
                {
			cell.output_Cai(nnn);
		}

		if (count % 1000 == 0)
                {
			cell.output_binary(nnn);
		}

                cell.output_ina_dvdt(nnn, bcl);

                cell.pace(stim);

                if( (t > T_total - delay_rest - 0.5 * dt) && (t < T_total - delay_rest + 0.5 * dt) )
                {
                        cell.sc.get_steady_state();
                        cell.cc.get_steady_state();  
                }

		count ++;
	}
        
	return 0;
}