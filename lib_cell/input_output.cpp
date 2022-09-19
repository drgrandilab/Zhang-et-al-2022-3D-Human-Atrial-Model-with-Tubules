#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <cstdlib>

#ifndef INPUT_OUTPUT_CPP
#define INPUT_OUTPUT_CPP

#include "input_output.h"


char read_char_from_gz(const char *filename, char *array, long int num) {
    long int i = 0;

    if (!array)
    {
        print_error_info_array_not_malloced(filename);
        exit(EXIT_FAILURE);
    }
    gzFile gz = gzopen(filename, "r");
    if (!gz) {
        print_error_info_file_not_found(filename);
        for (i = 0; i < num; ++i)
        {
            array[i] = 0;
        }
        return 0;
    }
    gzread(gz, array, num);  // because we are using the old version
    gzclose(gz);
    return 1;
}

char read_float_from_txt(const char *filename, float *array, long int num) {
    long int i = 0;

    if (!array)
    {
        print_error_info_array_not_malloced(filename);
        exit(EXIT_FAILURE);
    }

    FILE *in = fopen(filename, "r");
    if (!in) {
        print_error_info_file_not_found(filename);
        for (i = 0; i < num; ++i)
        {
            array[i] = 0;
        }
        return 0;
    }

    for (i = 0; i < num; ++i)
    {
        fscanf(in, "%f\n", & array[i]);
    }
    fclose(in);
    return 1;
}



char read_double_from_txt(const char *filename, double *array, long int num) {
    long int i = 0;

    if (!array)
    {
        print_error_info_array_not_malloced(filename);
        exit(EXIT_FAILURE);
    }

    FILE *in = fopen(filename, "r");
    if (!in) {
        print_error_info_file_not_found(filename);
        for (i = 0; i < num; ++i)
        {
            array[i] = 0;
        }
        return 0;
    }

    for (i = 0; i < num; ++i)
    {
        fscanf(in, "%lf\n", & array[i]);
    }
    fclose(in);
    return 1;
}



char read_float_from_bin(const char *filename, float *array, long int num) {
    long int i = 0;

    if (!array)
    {
        print_error_info_array_not_malloced(filename);
        exit(EXIT_FAILURE);
    }
    FILE *in = fopen(filename, "rb");
    if (!in) {
        print_error_info_file_not_found(filename);
        for (i = 0; i < num; ++i)
        {
            array[i] = 0;
        }
        return 0;
    }

    long int read_num = fread(array, sizeof(float), num, in);
    if (read_num != num) {
        fprintf(stderr, "%ld out of %ld read, does not match! in file %s\n", read_num, num, filename);
        exit(EXIT_FAILURE);
        return 0;
    }
    fclose(in);
    return 1;
}

char read_double_from_bin(const char *filename, double *array, long int num) {
    long int i = 0;

    if (!array)
    {
        print_error_info_array_not_malloced(filename);
        exit(EXIT_FAILURE);
    }
    FILE *in = fopen(filename, "rb");
    if (!in) {
        print_error_info_array_not_malloced(filename);
        for (i = 0; i < num; ++i)
        {
            array[i] = 0;
        }
        return 0;
    }

    long int read_num = fread(array, sizeof(double), num, in);
    if (read_num != num) {
        fprintf(stderr, "%ld out of %ld read, does not match! in file %s\n", read_num, num, filename);
        exit(EXIT_FAILURE);
    }
    fclose(in);
    return 1;
}


void print_error_info_file_not_found(const char *filename) {
    fprintf(stderr, "\n\n\n");
    fprintf(stderr, "**************************************\n");
    fprintf(stderr, "can not open file %s\n", filename);
    fprintf(stderr, "**************************************\n");
    fprintf(stderr, "assign array to be 0 \n");
    fprintf(stderr, "\n\n\n");
}

void print_error_info_array_not_malloced(const char *filename) {
    fprintf(stderr, "\n\n\n");
    fprintf(stderr, "**************************************\n");
    fprintf(stderr, "The array was not allocated when trying to read %s\n", filename);
    fprintf(stderr, "**************************************\n");
    fprintf(stderr, "\n\n\n");
}

void print_error_info_filename_empty(const char *filename) {
    fprintf(stderr, "\n\n\n");
    fprintf(stderr, "**************************************\n");
    fprintf(stderr, "filename was not assigned %s\n", filename);
    fprintf(stderr, "**************************************\n");
    fprintf(stderr, "\n\n\n");
}


void output_matrix(const char *filename, const double **matrix, int row, int col) {

    FILE *out;
    out = fopen( filename, "w+");
    int i, j;
    for (i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            fprintf(out, "%f ", matrix[i][j]);
        }
        fprintf(out, "\n");
    }
    fclose(out);
}


char output_float_array_bin(const char *filename, const float *array, long int num) {
    if (!filename) {
        print_error_info_filename_empty(filename);
        exit(EXIT_FAILURE);
    }

    FILE *file;
    // open the file
    file = fopen(filename, "wb");
    if (!file) {
        perror(filename);
        print_error_info_array_not_malloced(filename);
        exit(EXIT_FAILURE);
    }

    // output
    long int rw = fwrite(array, sizeof(float), num, file);
    if (rw != num) {
        fprintf(stderr, "%ld/%ld floats written to '%s', exiting!\n", rw, num, filename);
        exit(EXIT_FAILURE);
    }
    if (file)
        fclose(file);

    return 1;
}


char output_double_array_bin(const char *filename, const double *array, long int num) {
    if (!filename) {
        print_error_info_filename_empty(filename);
        exit(EXIT_FAILURE);
    }

    FILE *file;
    // open the file
    file = fopen(filename, "wb");
    if (!file) {
        perror(filename);
        print_error_info_file_open_failure(filename);
        exit(EXIT_FAILURE);
    }

    // output
    long int rw = fwrite(array, sizeof(double), num, file);
    if (rw != num) {
        fprintf(stderr, "%ld/%ld floats written to '%s', exiting!\n", rw, num, filename);
        exit(EXIT_FAILURE);
    }
    if (file)
        fclose(file);
    return 1;
}


char output_double_array_txt(FILE *file, const double *array, long int num) {

    if (!file) {
        perror("file is not open, please double check\n");
        // print_error_info_array_not_malloced(filename);
        exit(EXIT_FAILURE);
    }

    // output
    // long int rw = fwrite(array, sizeof(double), num, file);
    if (!array)
    {
        perror("The array to output is NULL!\n");
        exit(EXIT_FAILURE);
    }
    int i = 0;
    for (i = 0; i < num; i++) {
        fprintf(file, "%f ", array[i]);
    }
    fprintf(file, "\n");  // add line break at the end of line;

    return 1;
}


char output_double_array_txt(const char *filename, const double *array, long int num) {
    if (!filename) {
        print_error_info_filename_empty(filename);
        exit(EXIT_FAILURE);
    }

    FILE *file;
    // open the file
    file = fopen(filename, "w+");

    if (!file) {
        perror("file is not open, please double check\n");
        // print_error_info_array_not_malloced(filename);
        exit(EXIT_FAILURE);
    }

    // output
    // long int rw = fwrite(array, sizeof(double), num, file);
    if (!array)
    {
        perror("The array to output is NULL!\n");
        exit(EXIT_FAILURE);
    }
    int i = 0;
    for (i = 0; i < num; i++) {
        fprintf(file, "%.10f ", array[i]);
    }
    fprintf(file, "\n");  // add line break at the end of line;
    if (file)
        fclose(file);

    return 1;
}


bool OutPutArrayToTxt(std::ofstream &output, const double *array, long int num )
 {
    if (output.is_open()) {
        if (array) {
            for (int i = 0; i < num; ++i)
            {
                output << array[i] << " ";
            }
            output << std::endl;
            return true;
        } else {
            std::cerr << "The array to output is NULL!\n";
            std::exit(0);
        }


    } else {
        std::cerr << "file is not open, please double check\n";
        std::exit(0);
    }
}


template <typename T>
void print_error_info_file_open_failure(T filename) {

    fprintf(stderr, "\n\n\n");
    fprintf(stderr, "**************************************\n");
    fprintf(stderr, "can not open file %s\n", filename);
    fprintf(stderr, "**************************************\n");
    fprintf(stderr, "assign array to be 0 \n");
    fprintf(stderr, "\n\n\n");
}

#endif // end of INPUT_OUTPUT_C