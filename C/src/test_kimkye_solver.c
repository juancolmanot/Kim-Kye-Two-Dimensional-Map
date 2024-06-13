#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define EPSILON 1e-6

struct kimkyeparams {
    long double alpha;
    long double beta;
    unsigned int n_map;
};

struct solveparams {
    long double threshold;
    long double atol;
    unsigned int n_max;
};

void matrix_to_array(
    long double **xpoint,
    long double xarr[],
    unsigned int size,
    unsigned int m,
    unsigned int n
){
    for (unsigned int i = 0; i < m; i++){
        for (unsigned int j = 0; j < n; j++){
            xarr[i * n + j] = xpoint[i][j];
        }
    }
}

void kimkye_n(
    long double *x,
    struct kimkyeparams *params,
    long double *fx
){
    long double alpha = params->alpha;
    long double beta = params->beta;
    unsigned int n_map = params->n_map;

    long double xn = x[0];
    long double yn = x[1];

    long double xn1 = 0, yn1 = 0;

    for (unsigned int i = 0; i < n_map; i++){
        xn1 = (4 * alpha * xn + beta * yn) * (1 - xn);
        yn1 = (4 * alpha * yn + beta * xn) * (1 - yn);
        xn = xn1;
        yn = yn1;
    }

    xn1 = xn1 - x[0];
    yn1 = yn1 - x[1];

    fx[0] = xn1;
    fx[1] = yn1;
}

void kimkye_n_map(
    long double *x,
    struct kimkyeparams *params,
    long double *fx
){
    long double alpha = params->alpha;
    long double beta = params->beta;
    unsigned int n_map = params->n_map;

    long double xn = x[0];
    long double yn = x[1];

    long double xn1 = 0, yn1 = 0;

    for (unsigned int i = 0; i < n_map; i++){
        xn1 = (4 * alpha * xn + beta * yn) * (1 - xn);
        yn1 = (4 * alpha * yn + beta * xn) * (1 - yn);
        xn = xn1;
        yn = yn1;
    }

    fx[0] = xn1;
    fx[1] = yn1;
}

void kimkye_jac(
    void (*func)(
        long double *,
        struct kimkyeparams *,
        long double *
    ),
    long double *x,
    struct kimkyeparams *params,
    long double **jac
){

    long double *fx1 = calloc(2, sizeof(long double));
    long double *fx2 = calloc(2, sizeof(long double));
    long double *fx0 = calloc(2, sizeof(long double));

    long double *xh1 = calloc(2, sizeof(long double));
    long double *xh2 = calloc(2, sizeof(long double));

    long double dfx, dfy, dgx, dgy;

    xh1[0] = x[0] + EPSILON;
    xh1[1] = x[1];
    xh2[0] = x[0];
    xh2[1] = x[1] + EPSILON;

    func(x, params, fx0);
    func(xh1, params, fx1);
    func(xh2, params, fx2);

    dfx = (fx1[0] - fx0[0]) / EPSILON;
    dfy = (fx1[1] - fx0[1]) / EPSILON;
    dgx = (fx2[0] - fx0[0]) / EPSILON;
    dgy = (fx2[1] - fx0[1]) / EPSILON;

    jac[0][0] = dfx;
    jac[0][1] = dfy;
    jac[1][0] = dgx;
    jac[1][1] = dgy;

    free(xh1);
    free(xh2);
    free(fx0);
    free(fx1);
    free(fx2);
}

void Newton_Rapshon_step(
    void (*func)(
        long double *,
        struct kimkyeparams *,
        long double *
    ),
    void (*jac_func)(
        void (*)(
            long double *,
            struct kimkyeparams *,
            long double *
        ),
        long double *,
        struct kimkyeparams *,
        long double **
    ),
    long double *x,
    struct kimkyeparams *params,
    long double *x1
){
    long double *f = calloc(2, sizeof(long double));
    long double det;
    long double **jac = (long double **)calloc(2, sizeof(long double));
    for (unsigned int i = 0; i < 2; i++){
        jac[i] = (long double *)calloc(2, sizeof(long double));
    }

    func(x, params, f);
    jac_func(func, x, params, jac);

    det = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];

    long double dx[2];

    dx[0] = -1 * (jac[1][1] * f[0] - jac[0][1] * f[1]) / det;
    dx[1] = -1 * (jac[0][0] * f[1] - jac[1][0] * f[0]) / det;

    x1[0] = x[0] + dx[0];
    x1[1] = x[1] + dx[1];

    free(f);
    free(jac);
}

void fsolve(
    void (*func)(
        long double *,
        struct kimkyeparams *,
        long double *
    ),
    void (*jac_func)(
        void (*)(
            long double *,
            struct kimkyeparams *,
            long double *
        ),
        long double *,
        struct kimkyeparams *,
        long double **
    ),
    long double *x,
    struct kimkyeparams *params,
    struct solveparams *solvep,
    long double *xroot,
    int *status
)
{
    long double threshold = solvep->threshold;
    long double atol = solvep->atol;
    unsigned int max_steps = solvep->n_max;

    const gsl_rng_type *Tr;
    gsl_rng *r;

    gsl_rng_env_setup();

    Tr = gsl_rng_default;
    r = gsl_rng_alloc(Tr);

    gsl_rng_set(r, (long unsigned int)time(NULL));

    long double *xprev = calloc(2, sizeof(long double));
    long double *xprev_2 = calloc(2, sizeof(long double));
    long double *x1 = calloc(2, sizeof(long double));
    long double *fx = calloc(2, sizeof(long double));
    xprev[0] = x[0];
    xprev[1] = x[1];

    *status = 0;

    for (unsigned int i = 0; i < max_steps; i++) {
        Newton_Rapshon_step(func, jac_func, xprev, params, x1);
        xprev[0] = x1[0];
        xprev[1] = x1[1];
        func(x1, params, fx);
        if (fabsl(fx[0]) > threshold || fabsl(fx[1]) > threshold) {
            xprev[0] = xprev_2[0] + (long double)gsl_ran_flat(r, -0.01, 0.01);
            xprev[1] = xprev_2[1] + (long double)gsl_ran_flat(r, -0.01, 0.01);
            xprev_2[0] = xprev[0];
            xprev_2[1] = xprev[1];
        }
        if (fabsl(fx[0]) < atol && fabsl(fx[1]) < atol) {
            *status = 1;
            break;
        }
    }

    xroot[0] = x1[0];
    xroot[1] = x1[1];

    free(xprev);
    free(xprev_2);
    free(x1);
    gsl_rng_free(r);
}


// void stability(
//     void (*func)(
//         long double *,
//         struct kimkyeparams *,
//         long double *
//     ),
//     void (*jac_func)(
//         void (*)(
//             long double *,
//             struct kimkyeparams *,
//             long double *
//         ),
//         long double *,
//         struct kimkyeparams *,
//         long double **
//     ),
//     long double *x,
//     unsigned int *condition
// ){
//     long double *f = calloc(2, sizeof(long double));
//     long double det;
//     long double **jac = (long double **)calloc(2, sizeof(long double));
//     for (unsigned int i = 0; i < 2; i++){
//         jac[i] = (long double *)calloc(2, sizeof(long double));
//     }
//
//     func(x, params, f);
//     jac_func(func, x, params, jac);
//
//     long double jacarr[2 * 2];
//
//     matrix_to_array(jac, jacarr, 2 * 2, 2, 2);
//
//     int n = 2;
//     int lda = n;
//     long double w[n];
//     long double vl[n * n];
//     long double vr[n * n];
//     int info;
//
//     info LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'V', 'V', n, jacarr, lda, w, v1, n, vr, n);
//
//
//
// }

int main() {

    FILE *f = fopen("../datafiles/test_solver.dat", "w");

    long double *x0 = calloc(2, sizeof(long double));
    struct kimkyeparams p = {0.6891, 0.5, 10};

    const gsl_rng_type *Tr;
    gsl_rng *r;

    gsl_rng_env_setup();

    Tr = gsl_rng_default;
    r = gsl_rng_alloc(Tr);

    gsl_rng_set(r, (long unsigned int)time(NULL));

    x0[0] = (long double)gsl_ran_flat(r, 0.0, 1.0);
    x0[1] = (long double)gsl_ran_flat(r, 0.0, 1.0);

    long double *fx = calloc(2, sizeof(long double));

    int status = 0;
    struct solveparams solvepar = {1, 1e-4, 5000};
    long double *xzero = calloc(2, sizeof(long double));

    long double a, b, n, h;
    a = 0;
    b = 1;
    n = 100;
    h = (b - a) / n;

    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = 0; j < n; j++) {
            x0[0] = h * (long double)i;
            x0[1] = h * (long double)j;
            fsolve(kimkye_n, kimkye_jac, x0, &p, &solvepar, xzero, &status);
            if (xzero[0] > 0 && xzero[0] < 1 && xzero[1] > 0 && xzero[1] < 1 && status == 1) {
                fprintf(f, "%5.10Lf %5.10Lf %5.10Lf %5.10Lf\n", xzero[0], xzero[1], fx[0], fx[1]);
            }
        }
    }

    free(x0);
    free(fx);
    free(xzero);
    gsl_rng_free(r);
    fclose(f);
    return 0;
}

