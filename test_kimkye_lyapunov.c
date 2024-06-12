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
    if (argc != 3){
        perror("Wrong amount of arguments.");
        printf("Usage: ./run-script.sh script_to_run params_file.ini\n");
        exit(EXIT_FAILURE);
    }

    // ===============================================================================
    // Load arguments into variables.
    // ===============================================================================
    const char *params_file = argv[2];

    // ===============================================================================
    // We load parameters in params1 variable
    // ===============================================================================
    Parameters3 params3;
    load_parameters_from_file(params_file, &params3, handler3);

	// ===============================================================================
    // Instanciate parameters for system.
    // ===============================================================================
    Parameters_kimkye *params = malloc(sizeof(Parameters_kimkye));
    params->alpha = params3.p_alpha;
    params->beta = params3.p_beta;
    int n_map = params3.n_map;
    int n = params3.n;
    // ===============================================================================
    // Create arrays of states and initialize them.
    // ===============================================================================
    long double *xn;
    long double *xn1 = calloc(2, sizeof(long double));
    random_uniform_n(&xn, 2, (long unsigned int)time(NULL), 0, 1);
    xn1[0] = xn1[1] = 0;

    // ===============================================================================
    // Evolve transient.
    // ===============================================================================
    relax_map_n(
	    map_kimkye,
	    xn1,
	    xn,
	    n_map,
	    transient,
	    params
	);

    // ===============================================================================
    // Compute trajectory.
    // ===============================================================================
    long double **x_trajectory;
    unsigned int n_steps = params3.n_steps;
    
    x_trajectory = trajectory(
    	map_kimkye,
    	xn,
    	n,
    	n_steps,
    	params,
    	n_map
    );

    // ===============================================================================
    // Compute lyapunov spectrum.
    // ===============================================================================
    long double *lyap;
    long double eps = params3.epsilon;
    lyap = lyapunov_spectrum(
    	map_kimkye,
    	x_trajectory,
    	n_steps,
    	params,
    	n,
    	n_map,
    	eps
    );

    for (int i = 0; i < n; i++) {
    	printf("L%d: %Lf\n", i + 1, lyap[i]);
    }
	
    // ===============================================================================
    // Print results.
    // ===============================================================================
    for (unsigned int i = 0; i < n_steps; i++) {
    	free(x_trajectory[i]);
    }

    free(x_trajectory);

	return 0;
}