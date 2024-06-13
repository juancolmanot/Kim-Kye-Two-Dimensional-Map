#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kimkye.h"
#include "linear_algebra.h"
#include "dynamical_systems.h"


void print_matrix(
    long double **mat,
    int rows,
    int cols
){
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%Lf ", mat[i][j]);
        }
        printf("\n");
    }
}

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

long double **trajectory(
    void (*func)(
        long double *,
        long double *,
        void *
    ),
    long double *x0,
    int n,
    unsigned int n_steps,
    void *params,
    int n_map
) {
    long double **x = (long double **)calloc(n_steps, sizeof(long double));
    for (unsigned int i = 0; i < n_steps; i++) {
        x[i] = (long double *)calloc((size_t)n, sizeof(long double));
    }

    long double *xn = calloc((size_t)n, sizeof(long double));
    long double *xn1 = calloc((size_t)n, sizeof(long double));

    xn[0] = x0[0];
    xn[1] = x0[1];

    for (unsigned int i = 0; i < n_steps; i++) {
        map_n(func, xn1, xn, n_map, params);
        x[i][0] = xn1[0];
        x[i][1] = xn1[1];
        xn[0] = xn1[0];
        xn[1] = xn1[1];
    }

    free(xn1);
    free(xn);
    return x;
}

long double *lyapunov_spectrum(
    void (*func)(
        long double *,
        long double *,
        void *
    ),
    long double **trajectory,
    unsigned int size_trajectory,
    void *params,
    int n,
    int n_map,
    long double epsilon
){
    // Lyapunov array
    long double *lyapunov = calloc((size_t)n, sizeof(long double));

    // Declare jacobian
    long double **jac = (long double **)malloc((size_t)n * sizeof(long double));

    // Declare matrices for Q and R.
    long double **Q, **R;
    Q = (long double **)calloc((size_t)n, sizeof(long double));
    R = (long double **)calloc((size_t)n, sizeof(long double));
    for (int i = 0; i < n; i++) {
        Q[i] = (long double *)calloc((size_t)n, sizeof(long double));
        R[i] = (long double *)calloc((size_t)n, sizeof(long double));
    }

    for (unsigned int i = 0; i < size_trajectory; i++) {
        
        jac = jacobian_numerical(
            map_n,
            func,
            trajectory[i],
            params,
            n,
            n_map,
            epsilon
        );

        qr_decomposition(jac, Q, R, n, n);

        for (int j = 0; j < n; j++){
            lyapunov[j] += logl(fabsl(R[j][j]));
        }
    }
    
    for (int i = 0; i < n; i++) {
        lyapunov[i] /= (long double)size_trajectory;
    }

    for (int i = 0; i < n; i++) {
        free(jac[i]);
        free(Q[i]);
        free(R[i]);
    }
    free(jac);
    free(Q);
    free(R);
    
    return lyapunov;
}