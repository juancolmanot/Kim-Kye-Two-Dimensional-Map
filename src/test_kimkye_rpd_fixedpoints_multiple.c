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
    unsigned int finalsize = (unsigned int)params1.N;
    
    // ===============================================================================
    // Arrays for histogram of RPD.
    // ===============================================================================
    unsigned int nbins = 100;
    long double *xi, *yi, **rpdi;
    xi = calloc(nbins, sizeof(long double));
    yi = calloc(nbins, sizeof(long double));
    rpdi = (long double **)calloc(nbins, sizeof(long double));
    for (unsigned int i = 0; i < nbins; i++){
        rpdi[i] = (long double *)calloc(nbins, sizeof(long double));
    }

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
    long double c = 5e-3;

    // ===============================================================================
    // Compute RPD.
    // ===============================================================================
    rpd_funct_3d_fixedpoints_multiple(
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
        fixed_points,
        npoints,
        c
    );

    // ===============================================================================
    // Write to file.
    // ===============================================================================
    printf("Writing to file: %s\n", write_filename);
    for (unsigned int i = 0; i < nbins; i++){
        for (unsigned int j = 0; j < nbins; j++){
            fprintf(f, "%3.7Lf %3.7Lf %3.7Lf\n", xi[i], yi[j], rpdi[i][j]);
        }
    }

    for (unsigned int i = 0; i < nbins; i++){
        free(rpdi[i]);
    }

    for (int = 0; i < rows; i++){
        free(data[i]);
    }

    free(data);
    fclose(f);
    free(fixed_points);
    free(rpdi);
    free(xn);
    free(xn1);
    free(xi);
    free(yi);
    free(params);
    // ===============================================================================
    // Finish statement.
    // ===============================================================================
    printf("Process finished. Results stored in %s\n", write_filename);
    return 0;
}

