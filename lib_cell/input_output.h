#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <cstdlib>


#ifndef INPUT_OUTPUT_H
#define INPUT_OUTPUT_H

char read_char_from_gz(const char *filename, char *array, long int num);
char read_float_from_txt(const char *filename, float *array, long int num);
char read_double_from_txt(const char *filename, double *array, long int num);
char read_float_from_bin(const char *filename, float *array, long int num);
char read_double_from_bin(const char *filename, double *array, long int num);
void print_error_info_file_not_found(const char *filename);
void print_error_info_array_not_malloced(const char *filename);
void output_matrix(const char * filename, const double **matrix, int row, int col);
char output_float_array_bin(const char *filename, const float *voltages, long int num);
char output_double_array_bin(const char *filename, const double *voltages, long int num);
void print_error_info_filename_empty(const char *filename);
char output_double_array_txt(FILE *file, const double *array, long int num);
char output_double_array_txt(const char *filename, const double *array, long int num);
bool OutPutArrayToTxt(std::ofstream &output, const double *array, long int num );
template <typename T>
void print_error_info_file_open_failure(T filename) ;
#endif