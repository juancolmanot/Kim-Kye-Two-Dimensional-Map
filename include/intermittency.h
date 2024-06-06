#ifndef INTERMITTENCY_H
#define INTERMITTENCY_H

bool check_reinjection_circular(
    long double *fixed_point,
    long double *xn,
    long double *xn1,
    unsigned int size,
    long double c
);

bool check_ejection_circular(
    long double *fixed_point,
    long double *xn,
    long double *xn1,
    unsigned int size,
    long double c
);

bool check_reinjection(
    long double *fixed_point,
    long double *xn,
    long double *xn1,
    unsigned int size,
    long double c
);

bool check_ejection(
    long double *fixed_point,
    long double *xn,
    long double *xn1,
    unsigned int size,
    long double c
);

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
);

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
);

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
);

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
);

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
);

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
);

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
);

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
);

#endif