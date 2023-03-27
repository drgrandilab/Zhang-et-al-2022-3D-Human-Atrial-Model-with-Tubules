### What is this code?
This is a random-walk-based 3D tubular generator.

### Assumptions made in this tubular generator:
0. Experiment-based parameters shown in line 9-15 in "tubule_population_generator.m".
1. The smallest step is 1.84[um] in x direction and 0.9[um] in y / z direction
2. The surface 2 layers of CRUs are treated as surface CRUs
3. The density of y-direction tubule is as same as that in z-direction
4. The ratio of axial and transversal tubular density is 6
5. The default cell contains 55*17*11 CRUs

### Adjustable parameters:
1. size of the cell (CELL_LEN, CELL_WID, CELL_DEP) in [# CRU]
2. size of the population (TUBULE_POPULATION_NUM)

### How to use this code?
0. The default tubular length of 20 cells (19 cells with tubule + 1 detubulated cell) used in the publication are shown in line 62-65 in "tubule_population_generator.m". If you need to generate new tubular structures, please comment line 62-65 and uncomment line 49-59.
1. Use Matlab to run "tubule_population_generator.m" (This could run without step 0)
2. Validation figures will be generated in the same folder
3. Specific data and paraview files will be output in 'TRIAL' folder