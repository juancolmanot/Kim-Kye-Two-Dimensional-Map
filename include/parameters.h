#ifndef PARAMETERS1_H
#define PARAMETERS1_H

typedef struct {
    long unsigned int N;
    unsigned int ntarget;
    double threshold;
    unsigned int transient;
    unsigned int n_map;
    long double p_alpha;
    double p_beta;
} Parameters1;

int handler1(
    void *user,
    const char* section,
    const char* name,
    const char* value
);

int handler2(
    void *user,
    const char* section,
    const char* name,
    const char* value
);

void load_parameters_from_file(
    const char* filename,
    void* params,
    int handler(
        void * ,
        const char* ,
        const char* ,
        const char* 
    )
);

#endif
