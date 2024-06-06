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
