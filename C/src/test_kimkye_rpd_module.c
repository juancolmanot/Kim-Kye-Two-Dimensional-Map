#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

void progress_status(
    time_t prev_time,
    double frec,
    unsigned int count,
    unsigned int obj
){
    time_t tnow;

    time(&tnow);

    double diff = difftime(tnow, prev_time);

    if (diff >= frec){
        printf("progress: %4.2f %%\n", (double)count / (double)obj);
        prev_time = tnow;
    }
}

long double la_min(long double* x, unsigned int length) {
    long double minimum = x[0];
    for (unsigned int i = 0; i < length; i++) {
        if (x[i] < minimum) {
            minimum = x[i];
        }
    }
    return minimum;
}

long double la_max(long double* x, unsigned int length) {
    long double maximum = 0;
    for (unsigned int i = 0; i < length; i++) {
        if (x[i] > maximum) {
            maximum = x[i];
        }
    }
    return maximum;
}

void histogram_3d(
    long double **z,
    long double *x,
    long double *y,
    long double *datax,
    long double *datay,
    unsigned int data_size,
    unsigned int bin_size
)
{
    long double xmin, xmax, ymin, ymax;
    xmin = la_min(datax, data_size);
    xmax = la_max(datax, data_size);
    ymin = la_min(datay, data_size);
    ymax = la_max(datay, data_size);

    long double dx, dy;
    dx = (xmax - xmin) / (long double)(bin_size - 1);
    dy = (ymax - ymin) / (long double)(bin_size - 1);

    for (unsigned int i = 0; i < bin_size; i++){
        x[i] = xmin + dx * (long double)i;
        y[i] = ymin + dy * (long double)i;
    }

    for (unsigned int k = 0; k < data_size; k++){
        for (unsigned int i = 0; i < bin_size - 1; i++){
            for (unsigned int j = 0; j < bin_size - 1; j++){
                if (datax[k] > (xmin + dx * (long double)i) && datax[k] <= (xmin + dx * (long double)(i + 1))){
                    if (datay[k] > (ymin + (long double)dy * j) && datay[k] <= (ymin + dy * (long double)(j + 1))){
                        z[i][j]++;
                    }
                }
            }
        }
    }

}

void stats_histogram(
    long double *y,
    long double *x,
    long double *data,
    unsigned int data_size,
    unsigned int xsize,
) {
    long double xmax = 0, xmin = 0;
    xmax = la_max(data, data_size);
    xmin = la_min(data, data_size);

    long double dx = (xmax - xmin) / (long double) (xsize - 1);

    for (unsigned int i = 0; i < xsize - 1; i++) {
        for (unsigned int j = 0; j < data_size; j++) {
            if (data[j] > (xmin + dx * i) && data[j] < (xmin + dx * (i + 1))) {
                y[i]++;
            }
        }
        x[i] = xmin + dx * i;
    }
}

typedef struct {
    long double alpha;
    double beta;
} Parameters;

long double distance(long double *xn, long double *xn1) {

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

void map(long double *xn1, long double *xn, void *params) {
    Parameters *p = (Parameters*)params;
    long double alpha = p->alpha;
    long double beta = p->beta;

    xn1[0] = 4 * alpha * xn[0] * (1 - xn[0]) + beta * xn[1] * (1 - xn[0]);
    xn1[1] = 4 * alpha * xn[1] * (1 - xn[1]) + beta * xn[0] * (1 - xn[1]);

}

void map_n(long double *xn1, long double *xn, int n, void *params) {
    long double *x1_aux = calloc(2, sizeof(long double));
    long double *x_aux = calloc(2, sizeof(long double));

    x1_aux[0] = xn1[0];
    x1_aux[1] = xn1[1];
    x_aux[0] = xn[0];
    x_aux[1] = xn[1];

    for (int i = 0; i < n; i++){
        map(x1_aux, x_aux, params);
        x_aux[0] = x1_aux[0];
        x_aux[1] = x1_aux[1];
    }

    xn1[0] = x1_aux[0];
    xn1[1] = x1_aux[1];

    free(x_aux);
    free(x1_aux);
}

void relax_map_n(
    long double *xn1,
    long double *x0,
    int n_return,
    unsigned int transient,
    void *params
) {
    long double *x, *x1;
    x = calloc(2, sizeof(long double));
    x1 = calloc(2, sizeof(long double));
    x[0] = x0[0];
    x[1] = x0[1];

    for (unsigned int i = 0; i < transient; i++){
        map_n(x1, x, n_return, params);
        x[0] = x1[0];
        x[1] = x1[1];
    }

    xn1[0] = x1[0];
    xn1[1] = x1[1];

    free(x);
    free(x1);
}

void rpd_funct(
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
    long double lowerbound,
    long double upperbound,
    unsigned int variable
) {
    long double *xn, *xn1;
    xn = calloc(2, sizeof(long double));
    xn1 = calloc(2, sizeof(long double));
    xn[0] = x0[0];
    xn[1] = x0[1];

    relax_map_n(xn1, xn, (int)n_return, ntransient, params);
    xn[0] = xn1[0];
    xn[1] = xn1[1];

    long double d;

    long double *xreinjected = calloc(ntarget, sizeof(long double));

    unsigned int rcount = 0, laminar = 0;

    time_t prev_time;
    double period = 1;

    time(&prev_time);

    for (long unsigned int i = 0; i < nevol; i++){
        map_n(xn1, xn, (int)n_return, params);
        d = distance(xn, xn1);
        if (d < threshold && laminar == 0){// && xn1[variable] < upperbound && xn1[variable] > lowerbound) {
            xreinjected[rcount] = xn[variable];
            laminar = 1;
            rcount++;
        }
        else if (d > threshold && laminar == 1)//  && (xn1[variable] > upperbound || xn1[variable] < lowerbound))
        {
            laminar = 0;
        }

        progress_status(prev_time, period, rcount, ntarget);

        xn[0] = xn1[0];
        xn[1] = xn1[1];

        if (rcount >= ntarget) {
            break;
        }
    }

    long double *resizexreinjected = realloc(xreinjected, (long unsigned int)rcount * sizeof(long double));
    xreinjected = resizexreinjected;
    *finalsize = rcount;

    long double minxr = la_min(xreinjected, rcount);
    long double maxxr = la_max(xreinjected, rcount);

    stats_histogram(rpdxi, xi, xreinjected, rcount, nbins, minxr, maxxr);

    free(xreinjected);
    free(xn);
    free(xn1);
}

void rpd_funct_3d(
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
    long double *lowerbound,
    long double *upperbound
) {
    long double *xn, *xn1;
    xn = calloc(2, sizeof(long double));
    xn1 = calloc(2, sizeof(long double));
    xn[0] = x0[0];
    xn[1] = x0[1];

    relax_map_n(xn1, xn, (int)n_return, ntransient, params);
    xn[0] = xn1[0];
    xn[1] = xn1[1];

    long double d;

    long double *xreinjected = calloc(ntarget, sizeof(long double));
    long double *yreinjected = calloc(ntarget, sizeof(long double));

    unsigned int rcount = 0, laminar = 0;

    time_t time_prev;
    double period = 1;

    time(&time_prev);

    for (long unsigned int i = 0; i < nevol; i++){
        map_n(xn1, xn, (int)n_return, params);
        d = distance(xn, xn1);
        if (d < threshold && laminar == 0){
            if (xn1[0] < upperbound[0] && xn1[0] > lowerbound[0]) {
                if (xn1[1] < upperbound[1] && xn1[1] > lowerbound[1]){
                    xreinjected[rcount] = xn[0];
                    yreinjected[rcount] = xn[1];
                    laminar = 1;
                    rcount++;
                }
            }
        }
        else if (d > threshold && laminar == 1)
        {
            if (xn1[0] > upperbound[0] || xn1[0] < lowerbound[0]){
                laminar = 0;
            }
            else if (xn1[1] > upperbound[1] || xn1[1] < lowerbound[1]){
                laminar = 0;
            }
        }

        progress_status(time_prev, period, rcount, ntarget);

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


int main() {

    FILE *fp1;
    
    fp1 = fopen("../datafiles/kimkye_rpd_funct_combine_xy_5.dat", "w");
    
    long unsigned int N = 500000000000;
    unsigned int ntarget = 500000;
    double threshold = 1e-5;
    unsigned int transient = 1000;
    unsigned int n_map = 10;
    long double p_alpha = 0.6890015659550103 - expl(-18);
    double p_beta = 0.5;
    long double* xn = calloc(2, sizeof(long double));
    long double* xn1 = calloc(2, sizeof(long double));
    unsigned int finalsize = (unsigned int)N;
    
    const gsl_rng_type *Tr;
    gsl_rng *r;

    gsl_rng_env_setup();

    Tr = gsl_rng_default;
    r = gsl_rng_alloc(Tr);

    xn[0] = (long double)gsl_ran_flat(r, 0.0, 1.0);
    xn[1] = (long double)gsl_ran_flat(r, 0.0, 1.0);
    
    xn1[0] = xn1[1] = 0;

    Parameters *params = malloc(sizeof(Parameters));
    params->alpha = p_alpha;
    params->beta = p_beta;

    long double *xi, *yi;

    unsigned int nbins = 200;

    xi = calloc((size_t)nbins, sizeof(long double));
    yi = calloc((size_t)nbins, sizeof(long double));

    long double *lbound, *ubound;
    lbound = calloc(2, sizeof(long double));
    ubound = calloc(2, sizeof(long double));
    lbound[0] = 0.185;
    ubound[0] = 0.1866;
    lbound[1] = 0.8885;
    ubound[1] = 0.91;

    long double **rpdi = (long double **)calloc(nbins, sizeof(long double));
    for (unsigned int i = 0; i < nbins; i++){
        rpdi[i] = (long double *)calloc(nbins, sizeof(long double));
    }

    // unsigned int variable = 1;
    
    rpd_funct_3d(xn, xi, yi, rpdi, nbins, N, params, threshold, ntarget, transient, n_map, &finalsize, lbound, ubound);

    for (unsigned int i = 0; i < nbins - 1; i++) {
        for (unsigned int j = 0; j < nbins - 1; j++){
            fprintf(fp1, "%5.10Lf %5.10Lf %5.10Lf\n", xi[i], yi[j], rpdi[i][j]);
        }
    }
   
    fclose(fp1);
    free(xn);
    free(xn1);
    free(xi);
    free(yi);
    free(rpdi);
    gsl_rng_free(r);
    free(params);
    return 0;
}

