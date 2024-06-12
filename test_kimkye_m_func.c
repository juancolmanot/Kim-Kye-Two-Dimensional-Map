#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

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


void bubble_sort(
    long double *arr,
    int n
) {
    int i, j;
    long double temp;
    for (i = 0; i < n-1; i++) {
        for (j = 0; j < n-i-1; j++) {
            if (arr[j] > arr[j+1]) {
                // Swap arr[j] and arr[j+1]
                temp = arr[j];
                arr[j] = arr[j+1];
                arr[j+1] = temp;
            }
        }
    }
}

void local_slope(
    long double *array_x,
    long double *array_y,
    long double *slopes,
    unsigned int size
){
    
    for (unsigned int i = 1; i < size - 1; i++) {
        slopes[i - 1] = (array_y[i] - array_y[i - 1])/ (array_x[i] - array_x[i - 1]);
    }
}

int main() {
    
    const char *filename = "../datafiles/test_kimkye_reinjected_rotated_region1.dat";
    long double **data;
    int rows, cols;
    
    read_data_file(filename, &data, &rows, &cols);
    
    FILE *fp;
    
    fp = fopen("../datafiles/test_kimkye_m_function_region1.dat", "w");
    
    long double *xreinj, *Mx;
    
    long double Mxac = 0;
    
    xreinj = calloc((size_t)rows, sizeof(long double));
    Mx = calloc((size_t)rows, sizeof(long double));
    
    for (int i = 1; i < rows; i++){
        xreinj[i] = data[i][0];
    }
    
    bubble_sort(xreinj, rows);
    
    for (int i = 1; i < rows; i++){
        Mxac += xreinj[i];
        Mx[i] += Mxac / (long double)(i);
    }
    
    long double *dMdx = calloc((size_t)(rows - 2), sizeof(long double));
    
    local_slope(xreinj, Mx, dMdx, (unsigned int)rows);
    
    for (int i = 1; i < rows - 1; i++){
        if (fabsl(xreinj[i]) > 1e-3){
            fprintf(fp, "%Lf %Lf %Lf\n", xreinj[i], Mx[i], dMdx[i - 1]);
        }
    }
        
    for (int i = 0; i < rows; i++){
        free(data[i]);
    }
    free(data);
        
    fclose(fp);
    free(xreinj);
    return 0;
}

