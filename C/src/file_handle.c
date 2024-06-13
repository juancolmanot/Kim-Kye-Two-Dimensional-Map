#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "file_handle.h"

FILE *open_file(
    const char *filename
){
    FILE *file = fopen(filename, "w");
    if (file == NULL){
        perror("Error opening file for writing");
        exit(EXIT_FAILURE);
    }

    return file;
}

void read_data_file(
    const char *filename,
    long double ***data,
    int *rows,
    int *cols
){
    FILE *file = fopen(filename, "r");
    if (file == NULL){
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }
    
    int col_count = 0;
    int row_count = 0;
    char line[1024];
    if (fgets(line, sizeof(line), file) != NULL){
        char *token = strtok(line, " \t\n");
        while (token != NULL){
            col_count++;
            token = strtok(NULL, " \t\n");
        }
    row_count++;
    }
    
    while (fgets(line, sizeof(line), file) != NULL){
        row_count++;
    }
    rewind(file);

    *data = (long double **)malloc((long unsigned int)row_count * sizeof(long double *));
    for (int i = 0; i < row_count; i++){
        (*data)[i] = (long double *)malloc((long unsigned int)col_count * sizeof(long double));
    }
    
    int row = 0;
    while (fgets(line, sizeof(line), file) != NULL){
        int col = 0;
        char *token = strtok(line, " \t\n");
        while (token != NULL) {
            (*data)[row][col] = strtold(token, NULL);
            col++;
            token = strtok(NULL, " \t\n");
        }
        row++;
    }
    
    *rows = row_count;
    *cols = col_count;
    
    fclose(file);
}

void read_data_file_unsigned(
    const char *filename,
    long double ***data,
    unsigned int *rows,
    unsigned int *cols
){
    FILE *file = fopen(filename, "r");
    if (file == NULL){
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }
    
    unsigned int col_count = 0;
    unsigned int row_count = 0;
    char line[1024];
    if (fgets(line, sizeof(line), file) != NULL){
        char *token = strtok(line, " \t\n");
        while (token != NULL){
            col_count++;
            token = strtok(NULL, " \t\n");
        }
    row_count++;
    }
    
    while (fgets(line, sizeof(line), file) != NULL){
        row_count++;
    }
    rewind(file);

    *data = (long double **)malloc((long unsigned int)row_count * sizeof(long double *));
    for (unsigned int i = 0; i < row_count; i++){
        (*data)[i] = (long double *)malloc((long unsigned int)col_count * sizeof(long double));
    }
    
    unsigned int row = 0;
    while (fgets(line, sizeof(line), file) != NULL){
        unsigned int col = 0;
        char *token = strtok(line, " \t\n");
        while (token != NULL) {
            (*data)[row][col] = strtold(token, NULL);
            col++;
            token = strtok(NULL, " \t\n");
        }
        row++;
    }
    
    *rows = row_count;
    *cols = col_count;
    
    fclose(file);
}