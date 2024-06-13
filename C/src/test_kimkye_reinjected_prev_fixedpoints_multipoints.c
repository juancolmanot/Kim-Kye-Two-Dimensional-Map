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
    const char *params_file = argv[4];

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
    // Instanciate parameters for system.
    // ===============================================================================
    Parameters_kimkye *params = malloc(sizeof(Parameters_kimkye));
    params->alpha = params1.p_alpha - expl(-18);
    params->beta = params1.p_beta;

    // ===============================================================================
    // Create arrays of states and initialize them.
    // ===============================================================================
    long double *xn;
    long double *xn1 = calloc(2, sizeof(long double));
    random_uniform_n(&xn, 2, (long unsigned int)time(NULL), 0, 1);
    xn1[0] = xn1[1] = 0;
    unsigned int finalsize = (unsigned int)params1.ntarget;
    
    // ===============================================================================
    // Arrays for store values of reinjection
    // ===============================================================================
    long double *xr, *yr, *xr_prev, *yr_prev;

    xr = calloc((size_t)params1.ntarget, sizeof(long double));
    yr = calloc((size_t)params1.ntarget, sizeof(long double));
    xr_prev = calloc((size_t)params1.ntarget, sizeof(long double));
    yr_prev = calloc((size_t)params1.ntarget, sizeof(long double));

    // ===============================================================================
    // Load fixed points from data array
    // ===============================================================================
    unsigned int npoints = (unsigned int)rows;
    long double **fixed_points = (long double **)calloc(npoints, sizeof(long double));
    for (unsigned int i = 0; i < npoints; i++){
        fixed_points[i] = (long double *)calloc(2, sizeof(long double));
    }
    for (unsigned int i = 0; i < npoints; i ++){
        for (unsigned int j = 0; j < 2; j++){
            fixed_points[i][j] = data[i][j + 1];
        }
    }
    long double c = 1e-2;

    // ===============================================================================
    // Compute reinjected points.
    // ===============================================================================
    reinject_states_prev_fixedpoints_multipoints(
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
        fixed_points,
        npoints,
        c
    );

    // ===============================================================================
    // Write to file.
    // ===============================================================================
    printf("Writing to file: %s\n", write_filename);
    for (unsigned int i = 0; i < finalsize; i++){
        fprintf(f, "%3.7Lf %3.7Lf %3.7Lf %3.7Lf\n", xr[i], yr[i], xr_prev[i], yr_prev[i]);
    }

        for (int i = 0; i < rows; i++){
        free(data[i]);
    }

    free(data);
    fclose(f);
    free(fixed_points);
    free(xn);
    free(xn1);
    free(xr);
    free(yr);
    free(xr_prev);
    free(yr_prev);
    free(params);
    // ===============================================================================
    // Finish statement.
    // ===============================================================================
    printf("Process finished. Results stored in %s\n", write_filename);
    return 0;
}

