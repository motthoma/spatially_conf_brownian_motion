#include <stdlib.h>
#include "array_utils.h"
#include "code_handling.h"

/*
 * Function that allocates memory for a 2 dimensional array of doubles
 */
double **UTILS_calloc_2Ddouble_array(int m, int n)
{
    double **array = malloc(m * sizeof(double *));
    if (!array) return NULL;

    double *data = calloc(m * n, sizeof(double));
    if (!data) {
        free(array);
        return NULL;
    }

    for (int i = 0; i < m; i++)
        array[i] = data + i * n;

    return array;
}

/**
 * Function that allocates memory for a 2 dimensional array of long ints
 */
long int **UTILS_calloc_2Dlint_array(int m, int n) {
    long int **array = malloc(m * sizeof(long int*));
    if (!array) return NULL;

    long int *data = calloc(m * n, sizeof(long int));
    if (!data) { free(array); return NULL; }

    for (int i = 0; i < m; i++) array[i] = data + i * n;

    return array;
}

double UTILS_max_double(double a, double b){
    /* returns max of doubles a and b */
    return a > b ? a : b;
}

void UTILS_copy_code() {
    CODEHAND_copy_file_to_dest("array_utils.c");
    CODEHAND_copy_file_to_dest("array_utils.h");
}
