#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

typedef struct {
    long double alpha;
    long double beta;
} Parameters;

long double dif(
    long double *xn,
    long double *xn1
){
    long double difference;

    difference = sqrtl(powl(xn1[0] - xn[0], 2) + powl(xn1[1] - xn[1], 2));

    return difference;
}

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

void map(
    long double *xn1,
    long double *xn,
    void *params
) {
    Parameters *p = (Parameters*)params;
    long double alpha = p->alpha;
    long double beta = p->beta;

    xn1[0] = 4 * alpha * xn[0] * (1 - xn[0]) + beta * xn[1] * (1 - xn[0]);
    xn1[1] = 4 * alpha * xn[1] * (1 - xn[1]) + beta * xn[0] * (1 - xn[1]);

}

void map_n(
    long double *xn1,
    long double *xn,
    int n,
    void *params
){
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


void delete_repeated(
    long double *x,
    unsigned int size,
    long double *xuniq,
    unsigned int *uniqsize
){
    unsigned int i, j, k;

    k = 0;
    for (i = 0; i < size; i++) {
        int isUnique = 1;
        for (j = 0; j < k; j++){
            if (x[i] == xuniq[j]){
                isUnique = 0;
                break;
            }
        }
        if (isUnique){
            xuniq[k] = x[i];
            k++;
        }
    }

    *uniqsize = k;

    long double *xureall = realloc(xuniq, *uniqsize * sizeof(long double));
    xuniq = xureall;
}

void find_unique(
    long double *array,
    long double **array_unique,
    unsigned int sizearray,
    unsigned int *sizeunique
){
    
    *array_unique = (long double*)malloc(sizearray * sizeof(long double));
    *sizeunique = 0;
    
    for (int i = 0; i < (int)sizearray;){
        int j = i + 1;
        
        while (j < (int)sizearray && array[i] == array[j]){
            j++;
        }
        
        if (j == i + 1) {
            (*array_unique)[*sizeunique] = array[i];
            (*sizeunique)++;
        }
        i = j;
    }
}

long double find_critical(
    void (*func)(
        long double *,
        long double *,
        int ,
        void *
    ),
    long double *x0,
    unsigned int N,
    unsigned int transient,
    int n_mapa,
    long double *alpha_0,
    long double beta,
    unsigned int max_iter,
    long double dalpha,
    long double tol,
    unsigned int nfixed
){

    long double alpha_prev = 0;
    long double *xn = calloc(2, sizeof(long double));
    long double *xn1 = calloc(2, sizeof(long double));

    long double **x = (long double**)calloc(2, sizeof(long double));
    for (unsigned int i = 0; i < 2; i++){
        x[i] = (long double*)calloc(N, sizeof(long double));
    }
    
    Parameters *p = malloc(sizeof(Parameters));
    p->alpha = *alpha_0;
    p->beta = beta;

    unsigned int jcount = 0;
    
    for (unsigned int j = 0; j < max_iter; j++){
        
        xn[0] = x0[0];
        xn[1] = x0[1];
        
        
        for (unsigned int i = 0; i < transient; i++){
            map_n(xn1, xn, n_mapa, p);
            xn[0] = xn1[0];
            xn[1] = xn1[1];
        }

        for (unsigned int i = 0; i < N; i++){
            map_n(xn1, xn, n_mapa, p);
        
            x[0][i] = xn[0];
            x[1][i] = xn[1];

            xn[0] = xn1[0];
            xn[1] = xn1[1];
        }
        
        long double *xunique, *yunique;
        unsigned int sizexunique, sizeyunique;
        xunique = calloc(N, sizeof(long double));
        yunique = calloc(N, sizeof(long double));
        
        delete_repeated(x[0], N, xunique, &sizexunique);
        delete_repeated(x[1], N, yunique, &sizeyunique);
        
        if (sizexunique == nfixed && sizeyunique == nfixed){
            if (p->alpha - alpha_prev < tol) {
                printf("Critical parameter found: %5.20Lf\n", p->alpha);
                *alpha_0 = p->alpha;
                break;
            }
            dalpha *= 0.95;
            p->alpha = alpha_prev;
        }
        alpha_prev = p->alpha;
        p->alpha = p->alpha + dalpha;
        *alpha_0 = p->alpha;

        jcount++;
    }
    
    if (jcount == max_iter - 1){
        printf("Couldn't find critical parameter\n");
    }
    
    free(xn);
    free(xn1);
    free(x);
    
    return p->alpha;
}

    

int main() {

    FILE *f = fopen("../datafiles/test_kimkye_evol.dat", "w");

    unsigned int N = 1000;
    unsigned int transient = 2000000;
    int n_mapa = 1;
    long double alpha0 = 0.689;
    long double beta = 0.5;
    long double dalpha = expl(-15);
    unsigned int max_iter = 1000;
    
    const gsl_rng_type *Tr;
    gsl_rng *r;
    gsl_rng_env_setup();

    Tr = gsl_rng_default;
    r = gsl_rng_alloc(Tr);

    gsl_rng_set(r, (unsigned int)time(NULL));

    long double *x0 = calloc(2, sizeof(long double));
    x0[0] = (long double)gsl_ran_flat(r, 0.0, 1.0);
    x0[1] = (long double)gsl_ran_flat(r, 0.0, 1.0);

    unsigned int nfixedpoints = 10;
    long double tolalpha = 1e-8;

    find_critical(
        map_n,
        x0,
        N,
        transient,
        n_mapa,
        &alpha0,
        beta,
        max_iter,
        dalpha,
        tolalpha,
        nfixedpoints
    );

    printf("%5.20Lf\n", alpha0);
    
    gsl_rng_free(r);
    free(x0);
    fclose(f);
    return 0;
}

