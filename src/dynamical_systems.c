#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dynamical_systems.h"

void map_n(
    void (*func)(
        long double *,
        long double *,
        void *
    ),
    long double *xn1,
    long double *xn,
    int n,
    void *params
){
    long double *x1_aux = calloc(2, sizeof(long double));
    long double *x_aux = calloc(2, sizeof(long double));

    x1_aux[0] = xn1[0];
    x1_aux[1] = xn1[1];
    x_aux[0] = xn[0];
    x_aux[1] = xn[1];

    for (int i = 0; i < n; i++){
        func(x1_aux, x_aux, params);
        x_aux[0] = x1_aux[0];
        x_aux[1] = x1_aux[1];
    }

    xn1[0] = x1_aux[0];
    xn1[1] = x1_aux[1];

    free(x_aux);
    free(x1_aux);
}

void relax_map_n(
    void (*func)(
        long double *,
        long double *,
        void *
    ),
    long double *xn1,
    long double *x0,
    int n_return,
    unsigned int transient,
    void *params
){
    long double *x, *x1;
    x = calloc(2, sizeof(long double));
    x1 = calloc(2, sizeof(long double));
    x[0] = x0[0];
    x[1] = x0[1];

    for (unsigned int i = 0; i < transient; i++){
        map_n(func, x1, x, n_return, params);
        x[0] = x1[0];
        x[1] = x1[1];
    }

    xn1[0] = x1[0];
    xn1[1] = x1[1];

    free(x);
    free(x1);
}