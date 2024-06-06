#ifndef KIMKYE_H
#define KIMKYE_H

typedef struct {
    long double alpha;
    double beta;
} Parameters_kimkye;

void map_kimkye(
    long double *xn1,
    long double *xn,
    void *params
);

#endif