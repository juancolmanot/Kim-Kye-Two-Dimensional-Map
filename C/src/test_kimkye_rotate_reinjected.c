// ===============================================================================
// Include libraries.
// ===============================================================================

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "../include/progress_handle.h"
#include "../include/file_handle.h"
#include "../include/kimkye.h"
#include "../include/dynamical_systems.h"
#include "../include/intermittency.h"
#include "../include/linear_algebra.h"
#include "../include/stats.h"
#include "../include/parameters.h"
#include "../include/ini.h"
#include "../include/random.h"


int main(int argc, char *argv[]) {

    // ===============================================================================
    // Check for right passed arguments.
    // ===============================================================================
    if (argc != 5){
        perror("Wrong amount of arguments.");
        printf("Usage: ./run-script.sh script_to_run datafiles/write_to_file.dat read_from_file.dat theta\n");
        printf("script_to_run must be a .c file in scr/ that is not a functions file.");
        printf("write_to_file.dat must be any .dat located int datafiles/.");
        printf("theta must be an angle from in grades from 0° - 360°");
        exit(EXIT_FAILURE);
    }

    // ===============================================================================
    // Load arguments into variables.
    // ===============================================================================
    const char *write_filename = argv[2];
    const char *read_file = argv[3];
    long double theta = (long double)atof(argv[4]);
    // ===============================================================================
    // File to write to
    // ===============================================================================
    FILE *f = open_file(write_filename);

    // ===============================================================================
    // Store data from read file.
    // ===============================================================================
    unsigned int rows = 0, cols = 0;
    long double **data;
    read_data_file_unsigned(read_file, &data, &rows, &cols);
    // ===============================================================================
    // Declare and initialize arrays to store states.
    // ===============================================================================
    long double *xr, *yr, *xr_prev, *yr_prev;
    xr = calloc(rows, sizeof(long double));
    yr = calloc(rows, sizeof(long double));
    xr_prev = calloc(rows, sizeof(long double));
    yr_prev = calloc(rows, sizeof(long double));

    // ===============================================================================
    // Load data into arrays.
    // ===============================================================================

    for (unsigned int i = 0; i < rows; i++){
    	xr[i] = data[i][1];
    	yr[i] = data[i][2];
    	xr_prev[i] = data[i][3];
    	yr_prev[i] = data[i][4];
    }
	
    // ===============================================================================
    // Rotate vectors.
    // ===============================================================================
    long double *xr_rot, *yr_rot, *xr_prev_rot, *yr_prev_rot;
    
    xr_rot = calloc(rows, sizeof(long double));
    yr_rot = calloc(rows, sizeof(long double));
    xr_prev_rot = calloc(rows, sizeof(long double));
    yr_prev_rot = calloc(rows, sizeof(long double));
    
    rotate_vectors(
	    xr,
	    yr,
	    xr_rot,
	    yr_rot,
	    theta,
	    rows
	);

	rotate_vectors(
	    xr_prev,
	    yr_prev,
	    xr_prev_rot,
	    yr_prev_rot,
	    theta,
	    rows
	);
	
	// ===============================================================================
    // Write data to file.
    // ===============================================================================
    for (unsigned int i = 0; i < rows; i++){
        fprintf(f, "%d %3.7Lf %3.7Lf %3.7Lf %3.7Lf\n", i, xr_rot[i], yr_rot[i], xr_prev_rot[i], yr_prev_rot[i]);
    }

    // ===============================================================================
    // Free allocated memory, close opened files.
    // ===============================================================================
	for (unsigned int i = 0; i < rows; i++){
        free(data[i]);
    }
	
    free(data);
    free(xr);
    free(yr);
    free(xr_prev);
    free(yr_prev);
    free(xr_rot);
    free(yr_rot);
    free(xr_prev_rot);
    free(yr_prev_rot);
    fclose(f);
	
    // ===============================================================================
    // Finish statement.
    // ===============================================================================
    //printf("Process finished. Results stored in %s\n", write_filename);
    return 0;
}