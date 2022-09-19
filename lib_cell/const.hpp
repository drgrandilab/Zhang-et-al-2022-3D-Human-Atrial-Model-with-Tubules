// Constants for 'cell.cpp'
double Temp = 310;     // [K]
double R = 8314;       // [J/kmol*K]
double Frdy = 96485;   // [C/mol]
double FoRT = Frdy/R/Temp;

// Fractional currents in compartments
double Fjunc = 0.11;
double Fsl = 1-Fjunc;

// // Fixed ion concentrations
double Cli = 15;   // Intracellular Cl  [mM]
double Clo = 150;  // Extracellular Cl  [mM]
double Ko = 5.4;   // Extracellular K   [mM]
double Nao = 140;  // Extracellular Na  [mM]

// K current parameters
double pNaK = 0.01833;
double gkp = 0.002;

// Cl current parameters
double GClB = 9e-3;        // [mS/uF]