#include "Spatial_Cell.h"

Spatial_Cell::Spatial_Cell (double _dt, int seed, int cell_ID, double bcl, std::string tubular_map_file) : sc(CSubcell(55, 17,  11,  1,  1, cell_ID))
	, cc (CCell()), m_cell_ID(cell_ID)
{
	dt              = _dt;
	sc.finemesh     = 1;

	sc.setdt(dt);
	cc.setdt(dt);

        // When bcl = 0.33s, [Na]i = 11.14 mM; When bcl = 1s, [Na]i = 9.375 mM
	double a        = 15;           // [mM]
	double b        = 0.6;          // [Hz]
        double nai_res  = 7.757289 + (15.13319 - 7.757289) / pow((1.0 + pow((bcl/129429500.0), 0.7283862)), 8754.445);

	sc.init(0.2, 500, nai_res);

	int n           = sc.nx * sc.ny * sc.nz;

	sc.set_lateral_Ttubule();       // detubulate the myocyte

	std::ifstream inTub(tubular_map_file.c_str());
	std::string Tubline;

	double TubPos[sc.nx * sc.ny * sc.nz];
	int tub_id = 0;
	int tub_count   = 0;

	while (std::getline(inTub, Tubline))
	{
		double Tubvalue;
		std::stringstream ssTub(Tubline);
		ssTub >> Tubvalue;
		sc.set_Ttubule( tub_id, Tubvalue );
		if (Tubvalue == 1)       tub_count++;
		++tub_id;
	}

	sc.generate_tubular_map();
        sc.set_new_Cmem();
        cout << sc.num_tubule << '\t' << sc.Cmem << endl;

	CiArray          = new float [n];
	CpArray          = new float [n];
	CsArray          = new float [n];
	CnsrArray        = new float [n];
	CjsrArray        = new float [n];
	RyRArray         = new float [n];

	double voltage_tmp      = -74.0;
	double dvdt_tmp         = 0.0;
	double dvdt             = 0.0;

	t = 0;

	output8        = fopen("./global_result/global_cai.txt", "w+");
	output9        = fopen("./global_result/global_current.txt", "w+");
        output100      = fopen("./global_result/ina.txt", "w+");
        output101      = fopen("./global_result/voltage.txt", "w+");

        cc.set_steady_state();
	sc.set_steady_state();
}

Spatial_Cell::~Spatial_Cell() {
	if (CiArray  ) delete [] CiArray  ;
	if (CpArray  ) delete [] CpArray  ;
	if (CsArray  ) delete [] CsArray  ;
	if (CnsrArray) delete [] CnsrArray;
	if (CjsrArray) delete [] CjsrArray;
	if (RyRArray ) delete [] RyRArray ;

	fclose(output8);
	fclose(output9);
        fclose(output100);
        fclose(output101);
}

void Spatial_Cell::pace(double stim) {
	// miss setting functions in cell.hpp and cell.cpp
	cc.set_ica( sc.compute_ica() );
	cc.set_icabk( sc.compute_icabk() );
	cc.set_ipca( sc.compute_ipca() );

	cc.set_incx_junc( sc.compute_incx_junc() );
	cc.set_incx_sl( sc.compute_incx_sl() );

	cc.set_I_ClCa_junc( sc.compute_I_ClCa_junc() );
	cc.set_I_ClCa_sl( sc.compute_I_ClCa_sl() );

	cc.set_I_CaNa_junc( sc.compute_ICaNa() );
	cc.set_I_CaK_junc( sc.compute_ICaK() );

	cc.set_I_NaK_junc( sc.get_inak_junc() );
	cc.set_I_NaK_sl( sc.get_inak_sl() );

	cc.set_I_Nabk_junc( sc.get_inabk_junc() );
	cc.set_I_Nabk_sl( sc.get_inabk_sl() );

	cc.set_I_Na_junc( sc.get_ina_junc() );
	cc.set_I_Na_sl( sc.get_ina_sl() );

        cc.set_I_SK_junc( sc.get_isk_junc() );
        cc.set_I_SK_sl( sc.get_isk_sl() );

	cc.pace( stim );

        sc.pace2( cc.v , dt);
	sc.update_LTCC_RyR(cc.v, dt);
        

	if (sc.ci[0] != sc.ci[0]) {
		std::cerr << sc.ci[0] <<  " NaNs at Cell ID " << m_cell_ID << " and t = " << t << std::endl << "Code existing ... !!!!! " << std::endl;
                std::exit(0);
	}

	t += dt;
}

void Spatial_Cell::pace_no_cc(double stim) {
	sc.pace2( cc.v , dt / 5.0 );
	sc.pace2( cc.v , dt / 5.0 );
	sc.update_LTCC_RyR( cc.v, dt);

	sc.pace2( cc.v , dt / 5.0 );
	sc.pace2( cc.v , dt / 5.0 );
	sc.pace2( cc.v , dt / 5.0 );

	if (sc.ci[0] != sc.ci[0]) 
        {
		std::cerr << sc.ci[0] <<  " NaNs at Cell ID " << m_cell_ID << " and t = " << t << std::endl << "Code existing ... !!!!! " << std::endl;;
		std::exit(0);
	}

	t += dt;
}

void Spatial_Cell::output_ina_dvdt(std::string folder, double bcl) {

        fprintf ( output100, "%f\t%f\n", t, sc.get_ina_junc()+sc.get_ina_sl() );

        // output dvdt
        double dvdt_tmp         = ( cc.v - voltage_tmp ) / dt;
        
        if ( dvdt_tmp > dvdt  ) 
        {
                dvdt    = dvdt_tmp;
        }
        
        voltage_tmp     = cc.v;
        
        if ( int(t/dt-50/dt)%int(bcl/dt) == 0.0 )
        {
                fprintf ( output101, "%f\t%f\n", t, dvdt );
                cout << (t-50) << "\t" << dvdt << endl;
                dvdt    = 0.0;
        }  
}

void Spatial_Cell::output_Cai(std::string folder) {
	fprintf ( output9, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", t, cc.v, sc.compute_ica(), sc.compute_ir(), sc.compute_serca(), sc.compute_incx(), sc.compute_icabk(), sc.compute_ipca(), sc.compute_ileak(), sc.get_ina_sl(), sc.get_ina_junc(), sc.get_inabk_sl(), sc.get_inabk_junc(), sc.get_inak_sl(), sc.get_inak_junc(), sc.compute_ICaK(), sc.compute_ICaNa(), cc._I_to, cc._I_kr, cc._I_ks, cc._I_ki, cc._I_kp, cc._I_kur, cc._I_KAch, cc._I_Clbk, sc.compute_I_ClCa_sl(), sc.compute_I_ClCa_junc(), sc.get_isk_sl() + sc.get_isk_junc() );

	fprintf ( output8, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", t, cc.v, sc.compute_avg_ci(), sc.compute_avg_cs(), sc.compute_avg_cp(), sc.compute_avg_cjsr(), sc.compute_avg_cnsr(), sc.get_avg_nai(), sc.get_avg_nas(), sc.get_avg_nap(), sc.compute_incx() );
}

void Spatial_Cell::output_binary(std::string folder) {

	// file3 : map of global calcium transient(ci)
	char            filename3 [1000];
	sprintf         (filename3, "./BinaryFiles/map_ci_%d.ID.%d.bin", int(nearbyint(t)), m_cell_ID);

	for ( int id = 0; id < sc.n; id ++ )       CiArray[id]     = float(sc.ci[id]);

	output_float_array_bin( filename3, CiArray, sc.n);

	// file5 : map of cleft area(cp)
	char            filename5 [1000];
	sprintf         (filename5, "./BinaryFiles/map_cp_%d.ID.%d.bin", int(nearbyint(t)), m_cell_ID);

	for ( int id = 0; id < sc.n; id ++ )       CpArray[id]     = float(sc.cp[id]);

	output_float_array_bin( filename5, CpArray, sc.n);

	// file7 : map of network sr(cnsr)
	char            filename7 [1000];
	sprintf         (filename7, "./BinaryFiles/map_cnsr_%d.ID.%d.bin", int(nearbyint(t)), m_cell_ID);

	for ( int id = 0; id < sc.nn; id ++ )       CnsrArray[id]     = float(sc.cnsr[id]);

	output_float_array_bin( filename7, CnsrArray, sc.nn);
}