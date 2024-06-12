#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parameters.h"
#include "ini.h"

int handler1(
    void *user,
    const char* section,
    const char* name,
    const char* value
){
    Parameters1 *p = (Parameters1*)user;

    #define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
    if (MATCH("general", "N")) {
 		p->N = strtoul(value, NULL, 10);
    }
    else if (MATCH("general", "ntarget")) {
    	p->ntarget = (unsigned int)strtoul(value, NULL, 10);
    }
	else if (MATCH("general", "threshold")) {
        p->threshold = strtod(value, NULL);
    } else if (MATCH("general", "transient")) {
        p->transient = (unsigned int)strtoul(value, NULL, 10);
    } else if (MATCH("general", "n_map")) {
        p->n_map = (unsigned int)strtoul(value, NULL, 10);
    } else if (MATCH("general", "p_alpha")) {
        p->p_alpha = strtold(value, NULL);
    } else if (MATCH("general", "p_beta")) {
        p->p_beta = strtod(value, NULL);
    } else {
        return 0;  // unknown section/name, error
    }
    return 1;
}

int handler2(
    void *user,
    const char* section,
    const char* name,
    const char* value
){
    Parameters2 *p = (Parameters2*)user;

    #define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
    if (MATCH("general", "N")) {
        p->N = strtoul(value, NULL, 10);
    }
    else if (MATCH("general", "transient")) {
        p->transient = (unsigned int)strtoul(value, NULL, 10);
    } else if (MATCH("general", "n_map")) {
        p->n_map = (unsigned int)strtoul(value, NULL, 10);
    } else if (MATCH("general", "p_alpha")) {
        p->p_alpha = strtold(value, NULL);
    } else if (MATCH("general", "p_beta")) {
        p->p_beta = strtod(value, NULL);
    } else {
        return 0;  // unknown section/name, error
    }
    return 1;
}

int handler3(
    void *user,
    const char* section,
    const char* name,
    const char* value
){
    Parameters3 *p = (Parameters3*)user;

    #define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
    if (MATCH("general", "n_steps")) {
        p->n_steps = (unsigned int)strtoul(value, NULL, 10);
    }
    else if (MATCH("general", "n")) {
        p->n = (int)strtoul(value, NULL, 10);
    }
    else if (MATCH("general", "n_map")) {
        p->n_map = (int)strtoul(value, NULL, 10);
    }
    else if (MATCH("general", "transient")){
        p->transient = (unsigned int)strtoul(value, NULL, 10);
    }
    else if (MATCH("general", "p_alpha")) {
        p->p_alpha = strtold(value, NULL);
    }
    else if (MATCH("general", "p_beta")) {
        p->p_beta = strtod(value, NULL);
    }
    else if (MATCH("general", "epsilon")) {
        p->epsilon = strtod(value, NULL);
    } else {
        return 0;  // unknown section/name, error
    }
    return 1;
}

int handler4(
    void *user,
    const char* section,
    const char* name,
    const char* value
){
    Parameters4 *p = (Parameters4*)user;

    #define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
    if (MATCH("general", "nevol")) {
        p->nevol = (unsigned int)strtoul(value, NULL, 10);
    }
    else if (MATCH("general", "n")) {
        p->n = (int)strtoul(value, NULL, 10);
    }
    else if (MATCH("general", "n_map")) {
        p->n_map = (int)strtoul(value, NULL, 10);
    }
    else if (MATCH("general", "transient")){
        p->transient = (unsigned int)strtoul(value, NULL, 10);
    }
    else if (MATCH("general", "alpha0")) {
        p->alpha0 = strtold(value, NULL);
    }
    else if (MATCH("general", "alphan")) {
        p->alphan = strtold(value, NULL);
    }    
    else if (MATCH("general", "p_beta")) {
        p->p_beta = strtod(value, NULL);
    }
    else if (MATCH("general", "nalphas")) {
        p->nalphas = (unsigned int)strtoul(value, NULL, 10);
    }
    else if (MATCH("general", "theta")) {
        p->theta = strtold(value, NULL);
    } else {
        return 0;  // unknown section/name, error
    }
    return 1;
}

void load_parameters_from_file(
    const char* filename,
    void* params,
    int handler(
        void * ,
        const char* ,
        const char* ,
        const char* 
    )
){
	if (ini_parse(filename, handler, params) < 0){
		printf("Can't load '%s'\n", filename);
		exit(1);
	}
}
