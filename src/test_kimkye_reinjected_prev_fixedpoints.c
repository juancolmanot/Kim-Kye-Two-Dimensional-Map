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
    if (argc != 6){
        perror("Wrong amount of arguments.");
        printf("Usage: ./run-script.sh script_to_run write_to_file.dat read_from_file.dat params_file.ini\n");
        printf("write_to_file usually has names like 'test_kimkye_rpd_fixed_points_region(i).dat'\n");
        printf("read_from_file is the file where the list of fixed points is stored, normally in test_kimkye_evol_fixed.dat'\n");
        exit(EXIT_FAILURE);
    }

    // ===============================================================================
    // Load arguments into variables.
    // ===============================================================================
    const char *write_filename = argv[2];
    const char *read_file = argv[3];
    int line_read = atoi(argv[4]);
    const char *params_file = argv[5];

    // ===============================================================================
    // We load parameters in params1 variable
    // ===============================================================================
    Parameters1 params1;
    load_parameters_from_file(params_file, &params1, handler1);    
    
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
    // Create arrays of states and initialize them.
    // ===============================================================================
    long double *xn;
    long double *xn1 = calloc(2, sizeof(long double));
    random_uniform_n(&xn, 2, (long unsigned int)time(NULL), 0, 1);
    xn1[0] = xn1[1] = 0;

    // ===============================================================================
    // Instanciate parameters for system.
    // ===============================================================================
    Parameters_kimkye *params = malloc(sizeof(Parameters_kimkye));
    params->alpha = params1.p_alpha - expl(-18);
    params->beta = params1.p_beta;

    // ===============================================================================
    // Load fixed points from data array, using the inicated line.
    // ===============================================================================

    long double *fixedpoint = calloc(2, sizeof(long double));
    fixedpoint[0] = data[line_read][1];
    fixedpoint[1] = data[line_read][2];
    long double c = 5e-3; // Size of laminar radius.
    unsigned int finalsize = (unsigned int)params1.ntarget;

    // ===============================================================================
    // Declare and initialize arrays to store states.
    // ===============================================================================
    long double *xr, *yr, *xr_prev, *yr_prev;
    xr = calloc(finalsize, sizeof(long double));
    yr = calloc(finalsize, sizeof(long double));
    xr_prev = calloc(finalsize, sizeof(long double));
    yr_prev = calloc(finalsize, sizeof(long double));

    // ===============================================================================
    // Compute reinjected and previous to reinjected points.
    // ===============================================================================
    
    printf("%Lf %Lf\n", fixedpoint[0], fixedpoint[1]);

    reinject_states_prev_fixedpoints(
        map_kimkye,
        xn,
        xr,
        yr,
        xr_prev,
        yr_prev,
        params1.N,
        params,
        params1.ntarget,
        params1.transient,
        params1.n_map,
        &finalsize,
        fixedpoint,
        c
    );

    // ===============================================================================
    // Write data to file.
    // ===============================================================================
    for (unsigned int i = 0; i < finalsize; i++){
        fprintf(f, "%d %3.7Lf %3.7Lf %3.7Lf %3.7Lf\n", i, xr[i], yr[i], xr_prev[i], yr_prev[i]);
    }

    // ===============================================================================
    // Free allocated memory, close opened files.
    // ===============================================================================

    for (int i = 0; i < rows; i++){
        free(data[i]);
    }

    free(data);
    free(fixedpoint);
    free(xr);
    free(yr);
    free(xr_prev);
    free(yr_prev);
    free(xn);
    free(xn1);
    free(params);
    fclose(f);

    // ===============================================================================
    // Finish statement.
    // ===============================================================================
    printf("Process finished. Results stored in %s\n", write_filename);
    return 0;
}