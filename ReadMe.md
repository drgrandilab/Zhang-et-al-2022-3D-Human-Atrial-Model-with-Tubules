C++ code for Zhang et al. three-dimensional model of human atrial myocyte coupling electrophysiology and calcium signaling, with replicated transversal-axial tubular system.
_____________________________________________________________________________________________________
### Contents:

* ```ReadMe.md```:					this file  

* ```Makefile```:					the default make file  
* ```build```:						folder of object files  
* ```run.slurm```:				the command file in slurm-managed system   

* ```pace_2.cpp```:					the main function with inputs: basic-cycle-length in ms, number of pacing beats, input file of tubular structure  
* ```lib_cell```: 					folder with all the functions  
* ```pool_tubule```: 				folder of input fles of tubular structures  

* ```global_result```:				folder of output files (whole-cell averaged ion concentrations and ion currents)  
* ```BinaryFiles```: 				folder of output binary files (local CRU calcium concentration)  

* ```steady_state_init```: 			folder of files with initial steady state variables  
* ```steady_state_output```: 		folder of files with varible values at the end of pacing period   

_____________________________________________________________________________________________________
### Sample command in the terminal command line:
```make; ./pace_2 1000 28 "pool_tubule/tub_input_ver2_10.txt"```
_____________________________________________________________________________________________________

### References:

Zhang, Xianwei, Haibo Ni, Stefano Morotti, Charlotte E. R. Smith, Daisuke Sato, William E. Louch, Andrew G. Edwards, and Eleonora Grandi. Mechanisms of Spontaneous Ca2+ Release-Mediated Arrhythmia in a Novel 3D Human Atrial Myocyte Model: I. Transverse-Axial Tubule Variation. The Journal of Physiology. https://doi.org/10.1113/JP283363.

Zhang, Xianwei, Charlotte E. R. Smith, Stefano Morotti, Andrew G. Edwards, Daisuke Sato, William E. Louch, Haibo Ni, and Eleonora Grandi. Mechanisms of Spontaneous Ca2+ Release-Mediated Arrhythmia in a Novel 3D Human Atrial Myocyte Model: II. Ca2+-Handling Protein Variation. The Journal of Physiology.  https://doi.org/10.1113/JP283602.

Please cite the above papers when using this model.
