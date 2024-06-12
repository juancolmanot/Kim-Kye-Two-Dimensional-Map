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
    if (argc != 4){
        perror("Wrong amount of arguments.");
        printf("Usage: ./run-script.sh script_to_run write_to_file.dat params_file.ini\n");
        printf("write_to_file usually has names like 'test_kimkye_rpd_fixed_points_region(i).dat'\n");
        exit(EXIT_FAILURE);
    }

    // ===============================================================================
    // Load arguments into variables.
    // ===============================================================================
    const char *write_filename = argv[2];
    const char *params_file = argv[3];

    // ===============================================================================
    // We load parameters in params1 variable
    // ===============================================================================
    Parameters2 params2;
    load_parameters_from_file(params_file, &params2, handler2);

    // ===============================================================================
    // File to write to
    // ===============================================================================
    FILE *f = open_file(write_filename);

    // ===============================================================================
    // Instanciate parameters for system.
    // ===============================================================================
    Parameters_kimkye *params = malloc(sizeof(Parameters_kimkye));
    params->alpha = params2.p_alpha;
    params->beta = params2.p_beta;

    // ===============================================================================
    // Create arrays of states and initialize them.
    // ===============================================================================
    long double *xn;
    long double *xn1 = calloc(2, sizeof(long double));
    random_uniform_n(&xn, 2, (long unsigned int)time(NULL), 0, 1);
    xn1[0] = xn1[1] = 0;

    // ===============================================================================
    // Evolve transient period
    // ===============================================================================

    relax_map_n(
        map_kimkye,
        xn1,
        xn,
        (int)params2.n_map,
        params2.transient,
        params
    );
    
    // ===============================================================================
    // Declare rotated arrays
    // ===============================================================================
    long double *xnrot = calloc(2, sizeof(long double));
    long double *xn1rot = calloc(2, sizeof(long double));
    long double theta = 45;

    for (unsigned int i = 0; i < params2.N; i++){
        map_n(map_kimkye, xn1, xn, (int)params2.n_map, params);
        rotate_state(xn, xnrot, theta);
        rotate_state(xn1, xn1rot, theta);
        fprintf(f, "%d %Lf %Lf %Lf %Lf\n", i, xnrot[0], xnrot[1], xn1rot[0], xn1rot[1]);

        xn[0] = xn1[0];
        xn[1] = xn1[1];
    }

    free(params);
    free(xn);
    free(xn1);
    free(xnrot);
    free(xn1rot);
    fclose(f);
    // ===============================================================================
    // Finish statement.
    // ===============================================================================
    printf("Process finished. Results stored in %s\n", write_filename);
    return 0;
}

    