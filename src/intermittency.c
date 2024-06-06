#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include "intermittency.h"
#include "stats.h"
#include "progress_handle.h"
#include "linear_algebra.h"
#include "dynamical_systems.h"

bool check_reinjection_circular(
    long double *fixed_point,
    long double *xn,
    long double *xn1,
    unsigned int size,
    long double c
){

    bool reinject = 0;
    long double drn = 0, drn1 = 0;

    for (unsigned int i = 0; i < size; i++){
        drn += (xn[i] - fixed_point[i]) * (xn[i] - fixed_point[i]);
        drn1 += (xn1[i] - fixed_point[i]) * (xn1[i] - fixed_point[i]);
    }

    drn = sqrtl(drn);
    drn1 = sqrtl(drn1);

    if (drn > c && drn1 < c){
        reinject = 1;
    }

    return reinject;
}

bool check_ejection_circular(
    long double *fixed_point,
    long double *xn,
    long double *xn1,
    unsigned int size,
    long double c
){

    bool eject = 0;
    long double drn = 0, drn1 = 0;

    for (unsigned int i = 0; i < size; i++){
        drn += (xn[i] - fixed_point[i]) * (xn[i] - fixed_point[i]);
        drn1 += (xn1[i] - fixed_point[i]) * (xn1[i] - fixed_point[i]);
    }

    drn = sqrtl(drn);
    drn1 = sqrtl(drn1);

    if (drn < c && drn1 > c){
        eject = 1;
    }

    return eject;
}

bool check_reinjection(
    long double *fixed_point,
    long double *xn,
    long double *xn1,
    unsigned int size,
    long double c
){

    bool reinject = 0;
    unsigned int count = 0;
    for (unsigned int i = 0; i < size; i++){
        if (xn[i] < fixed_point[i] - c || xn[i] > fixed_point[i] + c){
            count++;
        }
        if (xn1[i] > fixed_point[i] - c && xn1[i] < fixed_point[i] + c){
            count++;
        }
    }

    if (count == size * 2){
        reinject = 1;
    }

    return reinject;
}

bool check_ejection(
    long double *fixed_point,
    long double *xn,
    long double *xn1,
    unsigned int size,
    long double c
){

    bool eject = 0;
    unsigned int count = 0;
    for (unsigned int i = 0; i < size; i++){
        if (xn[i] > fixed_point[i] - c && xn[i] < fixed_point[i] + c){
            count++;
        }
        if (xn1[i] < fixed_point[i] - c || xn1[i] > fixed_point[i] + c){
            count++;
        }
    }

    if (count == size * 2){
        eject = 1;
    }

    return eject;
}


void reinject_states(
    void (*func)(
        long double *,
        long double *,
        void *
    ),
    long double *x0,
    long double *xr,
    long double *yr,
    long unsigned int nevol,
    void *params,
    double threshold,
    unsigned int ntarget,
    unsigned int ntransient,
    unsigned int n_return,
    unsigned int *finalsize
){
    long double *xn, *xn1;
    xn = calloc(2, sizeof(long double));
    xn1 = calloc(2, sizeof(long double));
    xn[0] = x0[0];
    xn[1] = x0[1];

    relax_map_n(func, xn1, xn, (int)n_return, ntransient, params);
    xn[0] = xn1[0];
    xn[1] = xn1[1];

    long double d;

    unsigned int rcount = 0, laminar = 0;

    time_t prev_time;
    
    double frec = 1;
    
    time(&prev_time);
    
    for (long unsigned int i = 0; i < nevol; i++){
        map_n(func, xn1, xn, (int)n_return, params);
        d = distance(xn, xn1);
        if (d < threshold && laminar == 0){
            xr[rcount] = xn1[0];
            yr[rcount] = xn1[1];
            laminar = 1;
            rcount++;
        }
        else if (d > threshold && laminar == 1)
        {
            laminar = 0;
        }
        
        print_progress(rcount, ntarget, &prev_time, frec);

        xn[0] = xn1[0];
        xn[1] = xn1[1];

        if (rcount >= ntarget) {
            break;
        }
    }

    long double *resizexr = realloc(xr, rcount * sizeof(long double));
    long double *resizeyr = realloc(yr, rcount * sizeof(long double));
    xr = resizexr;
    yr = resizeyr;
    *finalsize = rcount;
    
    free(xn);
    free(xn1);
}

void reinject_states_prev(
    void (*func)(
        long double *,
        long double *,
        void *
    ),
    long double *x0,
    long double *xr,
    long double *yr,
    long double *xr_prev,
    long double *yr_prev,
    long unsigned int nevol,
    void *params,
    double threshold,
    unsigned int ntarget,
    unsigned int ntransient,
    unsigned int n_return,
    unsigned int *finalsize
){
    long double *xn, *xn1;
    xn = calloc(2, sizeof(long double));
    xn1 = calloc(2, sizeof(long double));
    xn[0] = x0[0];
    xn[1] = x0[1];

    relax_map_n(func, xn1, xn, (int)n_return, ntransient, params);
    xn[0] = xn1[0];
    xn[1] = xn1[1];

    long double d;

    unsigned int rcount = 0, laminar = 0;

    time_t prev_time;
    
    double frec = 1;
    
    time(&prev_time);
    
    for (long unsigned int i = 0; i < nevol; i++){
        map_n(func, xn1, xn, (int)n_return, params);
        d = distance(xn, xn1);
        if (d < threshold && laminar == 0){
            xr[rcount] = xn1[0];
            yr[rcount] = xn1[1];
            xr_prev[rcount] = xn[0];
            yr_prev[rcount] = xn[1];
            laminar = 1;
            rcount++;
        }
        else if (d > threshold && laminar == 1){
            laminar = 0;
        }
        
        print_progress(rcount, ntarget, &prev_time, frec);

        xn[0] = xn1[0];
        xn[1] = xn1[1];

        if (rcount >= ntarget) {
            break;
        }
    }

    long double *resizexr = realloc(xr, rcount * sizeof(long double));
    long double *resizeyr = realloc(yr, rcount * sizeof(long double));
    long double *resizexr_prev = realloc(xr_prev, rcount * sizeof(long double));
    long double *resizeyr_prev = realloc(yr_prev, rcount * sizeof(long double));
    xr = resizexr;
    yr = resizeyr;
    xr_prev = resizexr_prev;
    yr_prev = resizeyr_prev;
    *finalsize = rcount;
    
    free(xn);
    free(xn1);

}

void reinject_states_fixedpoints(
    void (*func)(
        long double *,
        long double *,
        void *
    ),
    long double *x0,
    long double *xr,
    long double *yr,
    long unsigned int nevol,
    void *params,
    double threshold,
    unsigned int ntarget,
    unsigned int ntransient,
    unsigned int n_return,
    unsigned int *finalsize,
    long double *fixed_point,
    long double c
) {
    long double *xn, *xn1;
    xn = calloc(2, sizeof(long double));
    xn1 = calloc(2, sizeof(long double));
    xn[0] = x0[0];
    xn[1] = x0[1];

    relax_map_n(func, xn1, xn, (int)n_return, ntransient, params);
    xn[0] = xn1[0];
    xn[1] = xn1[1];

    unsigned int rcount = 0;

    time_t prev_time;
    
    double frec = 1;
    
    time(&prev_time);

    for (long unsigned int i = 0; i < nevol; i++){
        map_n(func, xn1, xn, (int)n_return, params);

        if (check_reinjection_circular(fixed_point, xn, xn1, 2, c) == 1){
            xr[rcount] = xn1[0];
            yr[rcount] = xn1[1];
            rcount++;
        }

        print_progress(rcount, ntarget, &prev_time, frec);

        xn[0] = xn1[0];
        xn[1] = xn1[1];

        if (rcount >= ntarget) {
            break;
        }
    }

    long double *resizexr = realloc(xr, rcount * sizeof(long double));
    long double *resizeyr = realloc(yr, rcount * sizeof(long double));
    xr = resizexr;
    yr = resizeyr;
    *finalsize = rcount;
    
    free(xn);
    free(xn1);
}

void reinject_states_prev_fixedpoints(
    void (*func)(
        long double *,
        long double *,
        void *
    ),
    long double *x0,
    long double *xr,
    long double *yr,
    long double *xr_prev,
    long double *yr_prev,
    long unsigned int nevol,
    void *params,
    double threshold,
    unsigned int ntarget,
    unsigned int ntransient,
    unsigned int n_return,
    unsigned int *finalsize,
    long double *fixed_point,
    long double c
) {
    long double *xn, *xn1;
    xn = calloc(2, sizeof(long double));
    xn1 = calloc(2, sizeof(long double));
    xn[0] = x0[0];
    xn[1] = x0[1];

    relax_map_n(func, xn1, xn, (int)n_return, ntransient, params);
    xn[0] = xn1[0];
    xn[1] = xn1[1];

    unsigned int rcount = 0;

    time_t prev_time;
    
    double frec = 1;
    
    time(&prev_time);

    for (long unsigned int i = 0; i < nevol; i++){
        map_n(func, xn1, xn, (int)n_return, params);

        if (check_reinjection_circular(fixed_point, xn, xn1, 2, c) == 1){
            xr[rcount] = xn1[0];
            yr[rcount] = xn1[1];
            xr_prev[rcount] = xn[0];
            yr_prev[rcount] = xn[1];
            rcount++;
        }

        print_progress(rcount, ntarget, &prev_time, frec);

        xn[0] = xn1[0];
        xn[1] = xn1[1];

        if (rcount >= ntarget) {
            break;
        }
    }

    long double *resizexr = realloc(xr, rcount * sizeof(long double));
    long double *resizeyr = realloc(yr, rcount * sizeof(long double));
    long double *resizexr_prev = realloc(xr_prev, rcount * sizeof(long double));
    long double *resizeyr_prev = realloc(yr_prev, rcount * sizeof(long double));
    xr = resizexr;
    yr = resizeyr;
    xr_prev = resizexr_prev;
    yr_prev = resizeyr_prev;
    *finalsize = rcount;
    
    free(xn);
    free(xn1);
}

void rpd_funct(
    void (*func)(
        long double *,
        long double *,
        void *
    ),
    long double *x0,
    long double *xi,
    long double *rpdxi,
    unsigned int nbins,
    long unsigned int nevol,
    void *params,
    double threshold,
    unsigned int ntarget,
    unsigned int ntransient,
    unsigned int n_return,
    unsigned int *finalsize,
    long double lbound,
    long double ubound,
    unsigned int variable
){
    long double *xn, *xn1;
    xn = calloc(2, sizeof(long double));
    xn1 = calloc(2, sizeof(long double));
    xn[0] = x0[0];
    xn[1] = x0[1];

    relax_map_n(func, xn1, xn, (int)n_return, ntransient, params);
    xn[0] = xn1[0];
    xn[1] = xn1[1];

    long double d;

    long double *xreinjected = calloc(ntarget, sizeof(long double));

    unsigned int rcount = 0, laminar = 0;

    for (long unsigned int i = 0; i < nevol; i++){
        map_n(func, xn1, xn, (int)n_return, params);
        d = distance(xn, xn1);
        if (d < threshold && laminar == 0){
            if (xn1[variable] > lbound && xn1[variable] < ubound){
                xreinjected[rcount] = xn1[variable];
                laminar = 1;
                rcount++;
            }
        }
        else if (d > threshold && laminar == 1)
        {
            if(xn1[variable] < lbound || xn1[variable] > ubound){
                laminar = 0;
            }
        }

        xn[0] = xn1[0];
        xn[1] = xn1[1];

        if (rcount >= ntarget) {
            break;
        }
    }

    long double *resizexreinjected = realloc(xreinjected, (long unsigned int)rcount * sizeof(long double));
    xreinjected = resizexreinjected;
    *finalsize = rcount;

    stats_histogram(rpdxi, xi, xreinjected, rcount, nbins);

    free(xreinjected);
    free(xn);
    free(xn1);
}

void rpd_funct_3d(
    void (*func)(
        long double *,
        long double *,
        void *
    ),
    long double *x0,
    long double *xi,
    long double *yi,
    long double **rpdi,
    unsigned int nbins,
    long unsigned int nevol,
    void *params,
    double threshold,
    unsigned int ntarget,
    unsigned int ntransient,
    unsigned int n_return,
    unsigned int *finalsize,
    long double *lbound,
    long double *ubound
){
    long double *xn, *xn1;
    xn = calloc(2, sizeof(long double));
    xn1 = calloc(2, sizeof(long double));
    xn[0] = x0[0];
    xn[1] = x0[1];

    relax_map_n(func, xn1, xn, (int)n_return, ntransient, params);
    xn[0] = xn1[0];
    xn[1] = xn1[1];

    long double d;

    long double *xreinjected = calloc(ntarget, sizeof(long double));
    long double *yreinjected = calloc(ntarget, sizeof(long double));

    unsigned int rcount = 0, laminar = 0;

    time_t prev_time;
    
    double frec = 1;
    
    time(&prev_time);

    for (long unsigned int i = 0; i < nevol; i++){
        map_n(func, xn1, xn, (int)n_return, params);
        d = distance(xn, xn1);
        if (d < threshold && laminar == 0){
            if (xn1[0] > lbound[0] && xn1[0] < ubound[0]){
                if (xn1[1] > lbound[1] && xn1[1] < ubound[1]){
                    xreinjected[rcount] = xn1[0];
                    yreinjected[rcount] = xn1[1];
                    laminar = 1;
                    rcount++;
                }
            }
        }

        else if (d > threshold && laminar == 1){
            if (xn1[0] < lbound[0] || xn1[0] > ubound[0]){
                laminar = 0;
            }
            if (xn1[1] < lbound[1] || xn1[1] > ubound[1]){
                laminar = 0;
            }
        }

        print_progress(rcount, ntarget, &prev_time, frec);

        xn[0] = xn1[0];
        xn[1] = xn1[1];
        
        if (rcount >= ntarget) {
            break;
        }
    }

    long double *resizexreinjected = realloc(xreinjected, (unsigned int)rcount * sizeof(long double));
    long double *resizeyreinjected = realloc(yreinjected, (unsigned int)rcount * sizeof(long double));
    xreinjected = resizexreinjected;
    yreinjected = resizeyreinjected;
    *finalsize = rcount;

    histogram_3d(rpdi, xi, yi, xreinjected, yreinjected, rcount, nbins);
    
    free(xn);
    free(xn1);
}

void rpd_funct_fixedpoints(
    void (*func)(
        long double *,
        long double *,
        void *
    ),
    long double *x0,
    long double *xi,
    long double *rpdxi,
    unsigned int nbins,
    long unsigned int nevol,
    void *params,
    unsigned int ntarget,
    unsigned int ntransient,
    unsigned int n_return,
    unsigned int *finalsize,
    long double *fixed_point,
    long double c,
    unsigned int variable
) {
    long double *xn, *xn1;
    xn = calloc(2, sizeof(long double));
    xn1 = calloc(2, sizeof(long double));
    xn[0] = x0[0];
    xn[1] = x0[1];

    relax_map_n(func, xn1, xn, (int)n_return, ntransient, params);
    xn[0] = xn1[0];
    xn[1] = xn1[1];

    long double *reinjected = calloc(ntarget, sizeof(long double));

    unsigned int rcount = 0;

    time_t prev_time;
    
    double frec = 1;
    
    time(&prev_time);

    for (long unsigned int i = 0; i < nevol; i++){
        map_n(func, xn1, xn, (int)n_return, params);
        
        if (check_reinjection_circular(fixed_point, xn, xn1, 2, c) == 1){
            reinjected[rcount] = xn1[variable];
            rcount++;
        }

        print_progress(rcount, ntarget, &prev_time, frec);

        xn[0] = xn1[0];
        xn[1] = xn1[1];

        if (rcount >= ntarget) {
            break;
        }
    }

    long double *resizereinjected = realloc(reinjected, (long unsigned int)rcount * sizeof(long double));
    reinjected = resizereinjected;
    *finalsize = rcount;

    stats_histogram(rpdxi, xi, reinjected, rcount, nbins);

    free(reinjected);
    free(xn);
    free(xn1);
}

void rpd_funct_3d_fixedpoints(
    void (*func)(
        long double *,
        long double *,
        void *
    ),
    long double *x0,
    long double *xi,
    long double *yi,
    long double **rpdi,
    unsigned int nbins,
    long unsigned int nevol,
    void *params,
    unsigned int ntarget,
    unsigned int ntransient,
    unsigned int n_return,
    unsigned int *finalsize,
    long double *fixed_point,
    long double c
) {
    long double *xn, *xn1;
    xn = calloc(2, sizeof(long double));
    xn1 = calloc(2, sizeof(long double));
    xn[0] = x0[0];
    xn[1] = x0[1];

    relax_map_n(func, xn1, xn, (int)n_return, ntransient, params);
    xn[0] = xn1[0];
    xn[1] = xn1[1];

    long double *xreinjected = calloc(ntarget, sizeof(long double));
    long double *yreinjected = calloc(ntarget, sizeof(long double));

    unsigned int rcount = 0;

    time_t prev_time;
    
    double frec = 1;
    
    time(&prev_time);

    for (long unsigned int i = 0; i < nevol; i++){
        map_n(func, xn1, xn, (int)n_return, params);
        
        if (check_reinjection_circular(fixed_point, xn, xn1, 2, c) == 1){
            xreinjected[rcount] = xn1[0];
            yreinjected[rcount] = xn1[1];
            rcount++;
        }

        print_progress(rcount, ntarget, &prev_time, frec);

        xn[0] = xn1[0];
        xn[1] = xn1[1];

        if (rcount >= ntarget) {
            break;
        }
    }

    long double *resizexreinjected = realloc(xreinjected, (long unsigned int)rcount * sizeof(long double));
    long double *resizeyreinjected = realloc(yreinjected, (long unsigned int)rcount * sizeof(long double));
    xreinjected = resizexreinjected;
    yreinjected = resizeyreinjected;
    *finalsize = rcount;

    histogram_3d(rpdi, xi, yi, xreinjected, yreinjected, rcount, nbins);

    free(xreinjected);
    free(yreinjected);
    free(xn);
    free(xn1);
}