#ifndef DYNAMICAL_SYSTEMS_H
#define DYNAMICAL_SYSTEMS_H

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
);

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
);

#endif