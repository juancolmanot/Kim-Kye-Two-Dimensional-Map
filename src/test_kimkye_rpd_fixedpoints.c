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

int main(int argc, char *argv[]) {
    
    if (argc != 5){
        perror("Wrong amount of arguments.");
        printf("Usage: ./run-script.sh filescript.c write_to_file read_from_file config.ini");
        exit(EXIT_FAILURE);
    }

    // We load parameters in params1 variable
    Parameters1 params1;
    load_parameters_from_file(argv[4], &params1, handler1);

    // Read file and open file to write data
    const char *read_file = argv[3];
    const char *write_filename = argv[2];
    FILE *f = open_file(write_filename);
    int rows = 0, cols = 0;
    long double **data;
    read_data_file(read_file, &data, &rows, &cols);
    // Get data needed
    
    //const char *filename = "../datafiles/test_kimkye_evol_fixed.dat";
    //char[200] datafilename = "../datafiles/test_kimkye_rpd_fixed_points_region"
    //f = fopen("../datafiles/test_kimkye_rpd_fixed_points_region1.dat", "w");

    /*long unsigned int N = 100000000000000;
    unsigned int ntarget = 10000;
    unsigned int transient = 100000;
    unsigned int n_map = 10;
    long double p_alpha = 0.6890015659550103 - expl(-18);
    double p_beta = 0.5;
    //long double* xn = calloc(2, sizeof(long double));
    long double* xn1 = calloc(2, sizeof(long double));
    */

    //xn[0] = (long double)gsl_ran_flat(r, 0.0, 1.0);
    //xn[1] = (long double)gsl_ran_flat(r, 0.0, 1.0);
    
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

    unsigned int nbins = 100;

    long double *xi, *yi, **rpdi;
    xi = calloc(nbins, sizeof(long double));
    yi = calloc(nbins, sizeof(long double));
    rpdi = (long double **)calloc(nbins, sizeof(long double));
    for (unsigned int i = 0; i < nbins; i++){
        rpdi[i] = (long double *)calloc(nbins, sizeof(long double));
    }

    unsigned int region_line = 1;
    long double *fixedpoint = calloc(2, sizeof(long double));
    fixedpoint[0] = data[region_line][1];
    fixedpoint[1] = data[region_line][2];
    long double c = 5e-3;
    unsigned int finalsize = (unsigned int)params1.N;

    rpd_funct_3d_fixedpoints(
        map_kimkye,
        xn,
        xi,
        yi,
        rpdi,
        nbins,
        params1.N,
        params,
        params1.ntarget,
        params1.transient,
        params1.n_map,
        &finalsize,
        fixedpoint,
        c
    );

    for (unsigned int i = 0; i < nbins; i++){
        for (unsigned int j = 0; j < nbins; j++){
            fprintf(f, "%3.7Lf %3.7Lf %3.7Lf\n", xi[i], yi[j], rpdi[i][j]);
        }
    }

    for (unsigned int i = 0; i < nbins; i++){
        free(rpdi[i]);
    }

    for (int i = 0; i < rows; i++){
        free(data[i]);
    }

    free(data);
    fclose(f);
    free(fixedpoint);
    free(rpdi);
    free(xn);
    free(xn1);
    free(xi);
    free(yi);
    gsl_rng_free(r);
    free(params);
    return 0;
}

