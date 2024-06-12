#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kimkye.h"

void map_kimkye(
    long double *xn1,
    long double *xn,
    void *params
){
    Parameters_kimkye *p = (Parameters_kimkye*)params;
    long double alpha = p->alpha;
    long double beta = p->beta;

    xn1[0] = 4 * alpha * xn[0] * (1 - xn[0]) + beta * xn[1] * (1 - xn[0]);
    xn1[1] = 4 * alpha * xn[1] * (1 - xn[1]) + beta * xn[0] * (1 - xn[1]);
}

void jac_kimkye(
    long double **J,
    long double *x,
    void *params
){
    Parameters_kimkye *p = (Parameters_kimkye*)params;
    long double alpha = p->alpha;
    long double beta = p->beta;

    J[0][0] = 4 * alpha * (1 - 2 * x[0]) - beta * x[1];
    J[0][1] = beta * (1 - x[0]);
    J[1][0] = beta * (1 - x[1]);
    J[1][1] = 4 * alpha * (1 - 2 * x[1]) - beta * x[0];
}