#ifndef RANDOM_H
#define RANDOM_H

// Generates random numbers array with uniform distribution from a to b, with size n.
void random_uniform_n(
    long double **xn,
    unsigned int n,
    long unsigned int seed,
    double a,
    double b
);

#endif
