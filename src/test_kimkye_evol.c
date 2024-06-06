#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

typedef struct {
    long double alpha;
    long double beta;
} Parameters;

long double dif(
    long double *xn,
    long double *xn1
){
    long double difference;

    difference = sqrtl(powl(xn1[0] - xn[0], 2) + powl(xn1[1] - xn[1], 2));

    return difference;
}

long double distance(long double *xn, long double *xn1) {

    long double *xs, *x1s;
    xs = calloc(2, sizeof(long double));
    x1s = calloc(2, sizeof(long double));

    xs[0] = (xn[0] + xn1[0]) / 2;
    xs[1] = (xn[1] + xn1[1]) / 2;

    x1s[0] = xs[0];
    x1s[1] = xs[1];

    long double d;

    d = sqrtl(
        powl(xn[0] - xs[0], 2)
        + powl(xn[1] - xs[1], 2)
        + powl(xn1[0] - x1s[0], 2)
        + powl(xn1[1] - x1s[1], 2)
        );

    return d;
}

void map(long double *xn1, long double *xn, void *params) {
    Parameters *p = (Parameters*)params;
    long double alpha = p->alpha;
    long double beta = p->beta;

    xn1[0] = 4 * alpha * xn[0] * (1 - xn[0]) + beta * xn[1] * (1 - xn[0]);
    xn1[1] = 4 * alpha * xn[1] * (1 - xn[1]) + beta * xn[0] * (1 - xn[1]);

}

void map_n(long double *xn1, long double *xn, int n, void *params) {
    long double *x1_aux = calloc(2, sizeof(long double));
    long double *x_aux = calloc(2, sizeof(long double));

    x1_aux[0] = xn1[0];
    x1_aux[1] = xn1[1];
    x_aux[0] = xn[0];
    x_aux[1] = xn[1];

    for (int i = 0; i < n; i++){
        map(x1_aux, x_aux, params);
        x_aux[0] = x1_aux[0];
        x_aux[1] = x1_aux[1];
    }

    xn1[0] = x1_aux[0];
    xn1[1] = x1_aux[1];

    free(x_aux);
    free(x1_aux);
}


int main() {

    FILE *f = fopen("../datafiles/test_kimkye_evol.dat", "w");

    unsigned int N = 1000;
    unsigned int transient = 20000;
    int n_mapa = 1;
    Parameters *p = malloc(sizeof(Parameters));
    p->alpha = 0.6891;
    p->beta = 0.5;

    const gsl_rng_type *Tr;
    gsl_rng *r;
    gsl_rng_env_setup();

    Tr = gsl_rng_default;
    r = gsl_rng_alloc(Tr);

    gsl_rng_set(r, (unsigned int)time(NULL));

    long double *x0 = calloc(2, sizeof(long double));
    x0[0] = (long double)gsl_ran_flat(r, 0.0, 1.0);
    x0[1] = (long double)gsl_ran_flat(r, 0.0, 1.0);

    long double *xn = calloc(2, sizeof(long double));
    long double *xn1 = calloc(2, sizeof(long double));

    xn[0] = x0[0];
    xn[1] = x0[1];

    long double d;

    for (unsigned int i = 0; i < transient; i++){
        map_n(xn1, xn, n_mapa, p);
        xn[0] = xn1[0];
        xn[1] = xn1[1];
    }

    for (unsigned int i = 0; i < N; i++){
        map_n(xn1, xn, n_mapa, p);

        d = distance(xn, xn1);
        fprintf(f, "%d %Lf %Lf %Lf\n", i, xn[0], xn[1], d);

        xn[0] = xn1[0];
        xn[1] = xn1[1];
    }

    free(p);
    gsl_rng_free(r);
    free(x0);
    free(xn);
    free(xn1);
    fclose(f);
    return 0;
}

