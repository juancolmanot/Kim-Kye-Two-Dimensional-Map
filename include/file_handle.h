#ifndef FILE_HANDLE_H
#define FILE_HANDLE_H

FILE *open_file(
    const char *filename
);

void read_data_file(
    const char *filename,
    long double ***data,
    int *rows,
    int *cols
);

#endif