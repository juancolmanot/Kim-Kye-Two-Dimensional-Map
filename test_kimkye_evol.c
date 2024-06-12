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
    Parameters1 params1;
    load_parameters_from_file(params_file, &params1, handler2);

    // ===============================================================================
    // File to write to
    // ===============================================================================
    FILE *f = open_file(write_filename);

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

    // ===============================================================================
    // Evolve transient period
    // ===============================================================================

    relax_map_n(
        map_kimkye,
        xn1,
        xn,
        (int)params1.n_map,
        params1.transient,
        params
    );
    
    // ===============================================================================
    // Evolve stationary and write to file
    // ===============================================================================
    for (unsigned int i = 0; i < params1.N; i++){
        map_n(map_kimkye, xn1, xn, (int)params1.n_map, params);

        fprintf(f, "%d %Lf %Lf\n", i, xn[0], xn[1]);
3
        xn[0] = xn1[0];
        xn[1] = xn1[1];
    }

    free(params);
    free(xn);
    free(xn1);
    fclose(f);
    // ===============================================================================
    // Finish statement.
    // ===============================================================================
    printf("Process finished. Results stored in %s\n", write_filename);
    return 0;
}

