#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "progress_handle.h"

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