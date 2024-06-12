#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "random.h"

// Generates random numbers array with uniform distribution from a to b, with size n.
void random_uniform_n(
    long double **xn,
    unsigned int n,
    long unsigned int seed,
    double a,
    double b
)
{
    *xn = calloc(2, sizeof(long double));

    const gsl_rng_type *Tr;
    gsl_rng *r;

    gsl_rng_env_setup();

    Tr = gsl_rng_default;
    r = gsl_rng_alloc(Tr);

    gsl_rng_set(r, seed);

    for (unsigned int i = 0; i < n; i++){
        (*xn)[i] = (long double)gsl_ran_flat(r, a, b);
    }
    gsl_rng_free(r);
}
