#include <omp.h>

#include <iostream>
#include <fstream>
/*using namespace std;*/
#include <stdio.h>
#include <string>
#include <sstream>
#include "subcell.hpp"
#include "input_output.h"
#include <cstdlib>
#include <string>

class Spatial_Cell
{
public:
	Spatial_Cell(double dt = 0.01, int _seed = 0, int _cell_ID = 0, double bcl = 1000, std::string tubular_map_file = "pool_tubule/tub_input_ver2_1.txt"); // seed for random number generators
	~Spatial_Cell();

	CSubcell sc;
	CCell cc;

	double t;
	double dt;
        double voltage_tmp;
        double dvdt;

	int m_cell_ID;
	int seed;
	float *CiArray   ;
	float *CpArray   ;
	float *CsArray   ;
	float *CnsrArray ;
	float *CjsrArray ;
	float *RyRArray  ;
	void pace(double stim);
	void pace_no_cc(double stim);

	FILE            *output8;
	FILE            *output9;
        FILE            *output100;
        FILE            *output101;

	void output_Cai(std::string folder);
	void output_binary(std::string folder);
        void output_ina_dvdt(std::string folder, double bcl);
};