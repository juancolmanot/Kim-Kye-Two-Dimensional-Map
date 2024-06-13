#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

void read_line_and_parse(
    const char *filename,
    int line_number,
    long double *xrange,
    long double *yrange
){
    FILE *file = fopen(filename, "r");
    if (file == NULL){
        perror("Error opening file\n");
        exit(EXIT_FAILURE);
    }
    
    char line[256];
    int current_line = 0;
    
    while (fgets(line, sizeof(line), file)) {
        current_line++;
        if (current_line == line_number){
            int region;
            sscanf(line, "region%d: x [%Lf - %Lf], y [%Lf - %Lf]", &region, &xrange[0], &yrange[0], &xrange[1], &yrange[1]);
            break;
        }
    }
    
    if (current_line < line_number){
        fprintf(stderr, "Error: Line %d not found in file\n", line_number);
        exit(EXIT_FAILURE);
    }
    
    fclose(file);
}

void print_variables(
    long double *vars,
    char *names[],
    unsigned int n,
    time_t *prev_time,
    double frec
){
    time_t tnow;
    
    time(&tnow);
    
    double diff = difftime(tnow, *prev_time);
    
    if (diff > frec){
        for (unsigned int i = 0; i < n; i++){
            printf("%s: %5.4Lf", names[i], vars[i]);
            if (i < n - 1) {
                printf(" ");
            }
        printf("\n");
        *prev_time = tnow;
        }
    }
}

void print_progress(
    unsigned int count,
    unsigned int target,
    time_t *prev_time,
    double frec
){
    time_t tnow;
    
    time(&tnow);
    
    double diff = difftime(tnow, *prev_time);
    
    if (diff > frec){
        printf("progress: %3.2f %%\n", (double)count * 100 / (double)target);
        *prev_time = tnow;
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
    long double maximum = x[0];
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
        // Find the bin index for datax[k]
        unsigned int i = (unsigned int)((datax[k] - xmin) / dx);
        if (i >= bin_size) i = bin_size - 1;

        // Find the bin index for datay[k]
        unsigned int j = (unsigned int)((datay[k] - ymin) / dy);
        if (j >= bin_size) j = bin_size - 1;

        // Increment the bin count
        z[i][j]++;
    }
}

void stats_histogram(
    long double *y,
    long double *x,
    long double *data,
    unsigned int data_size,
    unsigned int xsize
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

void linear_regression(
    long double *x,
    long double *y,
    unsigned int size,
    long double *m,
    long double *b
){
    long double sumx, sumy, sumxy, sumx2;
    sumx = sumy = sumxy = sumx2 = 0;
    
    for (unsigned int i = 0; i < size; i++){
        sumx += x[i];
        sumy += y[i];
        sumxy += x[i] * y[i];
        sumx2 += x[i] * x[i];
    }
    
    *m = ((long double)size * sumxy - sumx * sumy) / ((long double)size * sumx2 - sumx * sumx);
    *b = (sumy - *m * sumx) / (long double)size;
}

void rotate_state(
    long double *x,
    long double *xrot,
    long double theta
){
    long double theta_rad = theta * M_PI / 180;
    
    xrot[0] = x[0] * cosl(theta_rad) - x[1] * sinl(theta_rad);
    xrot[1] = x[0] * sinl(theta_rad) + x[1] * cosl(theta_rad);
}

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

void reinject_states(
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

    unsigned int rcount = 0, laminar = 0;

    time_t prev_time;
    
    double frec = 1;
    
    time(&prev_time);
    
    for (long unsigned int i = 0; i < nevol; i++){
        map_n(xn1, xn, (int)n_return, params);
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

    unsigned int rcount = 0, laminar = 0;

    time_t prev_time;
    
    double frec = 1;
    
    time(&prev_time);
    
    for (long unsigned int i = 0; i < nevol; i++){
        map_n(xn1, xn, (int)n_return, params);
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

    for (long unsigned int i = 0; i < nevol; i++){
        map_n(xn1, xn, (int)n_return, params);
        d = distance(xn, xn1);
        if (d < threshold && laminar == 0){
            xreinjected[rcount] = xn1[variable];
            laminar = 1;
            rcount++;
        }
        else if (d > threshold && laminar == 1)
        {
            laminar = 0;
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
    unsigned int *finalsize
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

    for (long unsigned int i = 0; i < nevol; i++){
        map_n(xn1, xn, (int)n_return, params);
        d = distance(xn, xn1);
        if (d < threshold && laminar == 0){
            xreinjected[rcount] = xn1[0];
            yreinjected[rcount] = xn1[1];
            laminar = 1;
            rcount++;
        }
        else if (d > threshold && laminar == 1){
            laminar = 0;
        }

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
    
    fp1 = fopen("../datafiles/test_kimkye_rpd_all.dat", "w");

    long unsigned int N = 500000000000;
    unsigned int ntarget = 10000;
    double threshold = 5e-6;
    unsigned int transient = 10000;
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

    unsigned int nbins = 100;
    long double *xi, *yi, **rpdi;

    xi = calloc((size_t)nbins, sizeof(long double));
    yi = calloc((size_t)nbins, sizeof(long double));
    rpdi = (long double **)calloc(nbins, sizeof(long double));
    for (unsigned int i = 0; i < nbins; i++){
        rpdi[i] = (long double *)calloc(nbins, sizeof(long double));
    }
    
    rpd_funct_3d(
        xn,
        xi,
        yi,
        rpdi,
        nbins,
        N,
        params,
        threshold,
        ntarget,
        transient,
        n_map,
        &finalsize
    );

    for (unsigned int i = 0; i < nbins; i++){
        for (unsigned int j = 0; j < nbins; j++){
            fprintf(fp1, "%3.7Lf %3.7Lf %3.7Lf\n", xi[i], yi[j], rpdi[i][j]);
        }
    }

    for (unsigned int i = 0; i < nbins; i++){
        free(rpdi[i]);
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

