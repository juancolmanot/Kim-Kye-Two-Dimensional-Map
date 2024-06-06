#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "linear_algebra.h"
#include "stats.h"

// Get min value from array
long double la_min(
	long double* x,
	unsigned int length
){
    long double minimum = x[0];
    for (unsigned int i = 0; i < length; i++) {
        if (x[i] < minimum) {
            minimum = x[i];
        }
    }
    return minimum;
}

// Get max value from array
long double la_max(
	long double* x,
	unsigned int length
){
    long double maximum = x[0];
    for (unsigned int i = 0; i < length; i++) {
        if (x[i] > maximum) {
            maximum = x[i];
        }
    }
    return maximum;
}

// Rotate two dimensional state an angle theta.
void rotate_state(
    long double *x,
    long double *xrot,
    long double theta
){
    long double theta_rad = theta * M_PI / 180;
    
    xrot[0] = x[0] * cosl(theta_rad) - x[1] * sinl(theta_rad);
    xrot[1] = x[0] * sinl(theta_rad) + x[1] * cosl(theta_rad);
}

// Rotate two vectors of two dimensional states an angle theta.
void rotate_vectors(
    long double *x,
    long double *y,
    long double *xr,
    long double *yr,
    long double theta,
    unsigned int size
){
    long double *state, *state_rot;
    state = calloc(2, sizeof(long double));
    state_rot = calloc(2, sizeof(long double));
    
    for (unsigned int i = 0; i < size; i++){
        state[0] = x[i];
        state[1] = y[i];
        rotate_state(state, state_rot, theta);
        xr[i] = state_rot[0];
        yr[i] = state_rot[1];
        
    }
    free(state);
    free(state_rot);
}

// Find angle that aligns x data with xaxis.
void align_axes(
    long double *x,
    long double *y,
    long double *xal,
    long double *yal,
    unsigned int size,
    long double theta0,
    long double tol,
    unsigned int max_iter,
    long double *err,
    unsigned int *errsize
){
    long double theta = theta0;
    long double dtheta = 4;
    unsigned int count = 0;
    long double m, b;
    
    for (unsigned int i = 0; i < max_iter; i++){
        rotate_vectors(x, y, xal, yal, theta, size);
        
        linear_regression(xal, yal, size, &m, &b);
        
        if (fabsl(m) < tol) {
            printf("Axes aligned, theta: %Lf, m: %Lf\n", theta, m);
            break;
        }
        theta -= m * dtheta;
        err[i] = m;
        count++;
    }
    
    long double *resizeerr = realloc(err, count * sizeof(long double));
    err = resizeerr;
    *errsize = count;
}

// Compute distance from 4 dimensional state to hyper-surface.
long double distance(
	long double *xn,
	long double *xn1
){
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