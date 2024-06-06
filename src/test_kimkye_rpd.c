#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "../include/progress_handle.h"
#include "../include/file_handle.h"
#include "../include/kimkye.h"
#include "../include/linear_algebra.h"
#include "../include/stats.h"
#include "../include/dynamical_systems.h"
#include "../include/intermittency.h"
#include "../include/parameters.h"
#include "../include/ini.h"

int main(int argc, char *argv[]) {

    if (argc != 6){
        perror("Wrong amount of arguments.");
        printf("Usage: ./run-script.sh filescript.c write_to_file read_from_file line_to_read config.ini");
        exit(EXIT_FAILURE);
    }

    // We load parameters in params1 variable
    Parameters1 params1;
    load_parameters_from_file(argv[5], &params1, handler1);

    // Read file and open file to write data
    const char *read_file = argv[3];
    const char *write_filename = argv[2];
    int readline = atoi(argv[4]);
    FILE *f = open_file(write_filename);
    int rows = 0, cols = 0;
    long double **data;
    read_data_file(read_file, &data, &rows, &cols);
    // Get data needed
    long double *fileline = calloc((size_t)cols, sizeof(long double));
    for (int i = 0; i < cols; i++){
        fileline[i] = data[readline][i];
    }
    for (int i = 0; i < rows; i++){
        free(data[i]);
    }
    free(data);

    // State variables for evolution
    long double *xn = calloc(2, sizeof(long double));
    long double *xn1 = calloc(2, sizeof(long double));
    const gsl_rng_type *Tr;
    gsl_rng *r;

    gsl_rng_env_setup();

    Tr = gsl_rng_default;
    r = gsl_rng_alloc(Tr);

    gsl_rng_set(r, (long unsigned int)time(NULL));

    xn[0] = (long double)gsl_ran_flat(r, 0.0, 1.0);
    xn[1] = (long double)gsl_ran_flat(r, 0.0, 1.0);
    xn1[0] = xn1[1] = 0;

    Parameters_kimkye *params = malloc(sizeof(Parameters_kimkye));
    params->alpha = params1.p_alpha - expl(-18);
    params->beta = params1.p_beta;

    long double *xi, *yi;

    unsigned int nbins = 100;

    xi = calloc((size_t)nbins, sizeof(long double));
    yi = calloc((size_t)nbins, sizeof(long double));

    long double *lbound, *ubound;
    lbound = calloc(2, sizeof(long double));
    ubound = calloc(2, sizeof(long double));

    lbound[0] = fileline[0];
    lbound[1] = fileline[1];
    ubound[0] = fileline[2];
    ubound[1] = fileline[3];

    long double **rpdi = (long double **)calloc(nbins, sizeof(long double));
    for (unsigned int i = 0; i < nbins; i++){
        rpdi[i] = (long double *)calloc(nbins, sizeof(long double));
    }
    
    unsigned int finalsize = (unsigned int)params1.N;

    rpd_funct_3d(
        map_kimkye,
        xn,
        xi,
        yi,
        rpdi,
        nbins,
        params1.N,
        params,
        params1.threshold,
        params1.ntarget,
        params1.transient,
        params1.n_map,
        &finalsize,
        lbound,
        ubound
    );

    for (unsigned int i = 0; i < nbins - 1; i++) {
        for (unsigned int j = 0; j < nbins - 1; j++){
            fprintf(f, "%5.10Lf %5.10Lf %5.10Lf\n", xi[i], yi[j], rpdi[i][j]);
        }
    }

    fclose(f);
    free(xn);
    free(xn1);
    free(xi);
    free(yi);
    free(rpdi);
    free(params);
    
    printf("Process finished. Results stored in %s\n", write_filename);
    return 0;
}

