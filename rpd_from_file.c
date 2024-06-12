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
    if (argc != 4	){
        perror("Wrong amount of arguments.");
        printf("Usage: ./run-script.sh script_to_run write_to_file.dat read_from_file.dat\n");
        exit(EXIT_FAILURE);
    }

    // ===============================================================================
    // Load arguments into variables.
    // ===============================================================================
    const char *write_filename = argv[2];
    const char *read_file = argv[3];
    
    // ===============================================================================
    // File to write to
    // ===============================================================================
    FILE *f = open_file(write_filename);

    // ===============================================================================
    // Store data from read file.
    // ===============================================================================
    int rows = 0, cols = 0;
    long double **data;
    read_data_file(read_file, &data, &rows, &cols);

    // ===============================================================================
    // Arrays for histogram of RPD.
    // ===============================================================================
    unsigned int nbins = 100;
    long double *xi_r, *yi_r, **rpdi_r;
    long double *xi_rprev, *yi_rprev, **rpdi_rprev;
    xi_r = calloc(nbins, sizeof(long double));
    yi_r = calloc(nbins, sizeof(long double));
    xi_rprev = calloc(nbins, sizeof(long double));
    yi_rprev = calloc(nbins, sizeof(long double));
    rpdi_r = (long double **)calloc(nbins, sizeof(long double));
    rpdi_rprev = (long double **)calloc(nbins, sizeof(long double));
    for (unsigned int i = 0; i < nbins; i++){
        rpdi_r[i] = (long double *)calloc(nbins, sizeof(long double));
        rpdi_rprev[i] = (long double *)calloc(nbins, sizeof(long double));
    }

    // ===============================================================================
    // Fill arrays of data.
    // ===============================================================================
    unsigned int size_arr = (unsigned int)rows;
    long double *xr, *yr, *xr_prev, *yr_prev;
    xr = calloc(size_arr, sizeof(long double));
    yr = calloc(size_arr, sizeof(long double));
    xr_prev = calloc(size_arr, sizeof(long double));
    yr_prev = calloc(size_arr, sizeof(long double));
    for (int i = 0; i < rows; i++){
    	xr[i] = data[i][1];
    	yr[i] = data[i][2];
    	xr_prev[i] = data[i][3];
    	yr_prev[i] = data[i][4];
    }

    // ===============================================================================
    // Compute RPD.
    // ===============================================================================
    
    histogram_3d(
		rpdi_r,
		xi_r,
		yi_r,
		xr,
		yr,
		size_arr,
		nbins
	);

	histogram_3d(
		rpdi_rprev,
		xi_rprev,
		yi_rprev,
		xr_prev,
		yr_prev,
		size_arr,
		nbins
	);
    
    // ===============================================================================
    // Write to write file.
    // ===============================================================================
    
    for (unsigned int i = 0; i < nbins; i++){
        for (unsigned int j = 0; j < nbins; j++){
            fprintf(f, "%3.7Lf %3.7Lf %3.7Lf %3.7Lf %3.7Lf %3.7Lf\n",
            	xi_r[i],
            	yi_r[j],
            	rpdi_r[i][j],
            	xi_rprev[i],
            	yi_rprev[j],
            	rpdi_rprev[i][j]
            );
        }
    }

    // ===============================================================================
    // Free allocated memory, close opened files.
    // ===============================================================================
    
    for (unsigned int i = 0; i < nbins; i++){
        free(rpdi_r[i]);
        free(rpdi_rprev[i]);
    }

    for (int i = 0; i < rows; i++){
        free(data[i]);
    }

    free(data);
    free(rpdi_r);
    free(rpdi_rprev);
    free(xi_r);
    free(yi_r);
    free(xi_rprev);
    free(yi_rprev);
    fclose(f);

    // ===============================================================================
    // Finish statement.
    // ===============================================================================
    printf("Process finished. Results stored in %s\n", write_filename);
    return 0;
}

