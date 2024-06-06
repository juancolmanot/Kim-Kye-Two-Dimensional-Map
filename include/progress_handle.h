#ifndef PROGRESS_HANDLE_H
#define PROGRESS_HANDLE_H

void print_progress(
    unsigned int count,
    unsigned int target,
    time_t *prev_time,
    double frec
);

void print_variables(
    long double *vars,
    char *names[],
    unsigned int n,
    time_t *prev_time,
    double frec
);

#endif