/*

            EE1103 Assignment - 10 : Gaussian Elimination
            
Developer : EE24B017
Date : 10th October 2024
Purpose : To read N linear equations from a text file and solve the simultaneous system of N linear equations using the Gaussian Elimination method.
Input(s) : filename N
Outputs(s): N values of x
*/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define TOLERANCE 1e-9  //Tolerance for checking zero pivots

//Function prototypes
void Gauss(double** a, double* b, double* x, int n, double tol, int* er);
void Eliminate(double** a, double* s, double* b, int n, double tol, int* er);
void Pivot(double** a, double* b, double* s, int n, int k);
void Substitute(double** a, int n, double* b, double* x);

//Gauss elimination with scaling and pivoting
void Gauss(double** a, double* b, double* x, int n, double tol, int* er) {
    int i;
    double* s = (double*)malloc(n * sizeof(double));
    //Initialize scaling factors
    for (i = 0; i < n; i++) {
        s[i] = fabs(a[i][0]);
        for (int j = 1; j < n; j++) {
            if (fabs(a[i][j]) > s[i]) {
                s[i] = fabs(a[i][j]);
            }
        }
    }
    //Elimination process
    Eliminate(a, s, b, n, tol, er);
    if (*er != -1) {
        Substitute(a, n, b, x);
    }
    free(s);
}

//Elimination phase of Gaussian elimination
void Eliminate(double** a, double* s, double* b, int n, double tol, int* er) {
    int k, i, j;
    double factor;

    for (k = 0; k < n - 1; k++) {
        //Pivot
        Pivot(a, b, s, n, k);

        //Check for zero pivot
        if (fabs(a[k][k]) < tol) {
            *er = -1;
            return;
        }

        // Eliminate
        for (i = k + 1; i < n; i++) {
            factor = a[i][k] / a[k][k];
            for (j = k; j < n; j++) {
                a[i][j] -= factor * a[k][j];
            }
            b[i] -= factor * b[k];
        }
    }

    //Check for zero diagonal element in the last row
    if (fabs(a[n - 1][n - 1]) < tol) {
        *er = -1;
    }
}

//Pivot to find the row with the largest element in the current column
void Pivot(double** a, double* b, double* s, int n, int k) {
    int p = k;
    double big = fabs(a[k][k] / s[k]);
    double dummy;
    for (int i = k + 1; i < n; i++) {
        dummy = fabs(a[i][k] / s[i]);
        if (dummy > big) {
            big = dummy;
            p = i;
        }
    }
    // Swap rows if necessary
    if (p != k) {
        for (int j = 0; j < n; j++) {
            dummy = a[p][j];
            a[p][j] = a[k][j];
            a[k][j] = dummy;
        }

        dummy = b[p];
        b[p] = b[k];
        b[k] = dummy;

        dummy = s[p];
        s[p] = s[k];
        s[k] = dummy;
    }
}

//Back substitution phase of Gaussian elimination
void Substitute(double** a, int n, double* b, double* x) {
    int i, j;
    double sum;

    x[n - 1] = b[n - 1] / a[n - 1][n - 1];

    for (i = n - 2; i >= 0; i--) {
        sum = 0;
        for (j = i + 1; j < n; j++) {
            sum += a[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / a[i][i];
    }
}

//Function to read input from a text file containing the equation
int read_matrix_and_list(const char *filename, int n, double ** matrix, double *list) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("Error opening file.\n");
        return 1;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            //Read the first n elements into the matrix
            if (fscanf(file, "%lf", &matrix[i][j]) != 1) {
                printf("Error reading file.\n");
                fclose(file);
                return 1;
            }
        }
        //Read the (n+1)th element into the list
        if (fscanf(file, "%lf", &list[i]) != 1) {
            printf("Error reading file.\n");
            fclose(file);
            return 1;
        }
    }
    fclose(file);
}


int main(int argc, char **argv) {
    int i, j, er = 0;
    int res = 0;
    double tol = TOLERANCE;

    //Checking for correct number of inputs
    if (argc != 3){
        printf("Usage : filename N\n");
        return 0;
    }

    char *file = argv[1];
    int n = atoi(argv[2]);
    
    // Allocate memory for the matrix and vectors
    double** a = (double**)malloc(n * sizeof(double*));
    for (i = 0; i < n; i++) {
        a[i] = (double*)malloc(n * sizeof(double));
    }
    double* b = (double*)malloc(n * sizeof(double));
    double* x = (double*)malloc(n * sizeof(double));
    
    //Read from the given file
    res = read_matrix_and_list(file, n, a, b);
    if (res == 1){
        return 0;
    }

    // Call the Gauss function
    Gauss(a, b, x, n, tol, &er);

    if (er != -1) {
        for (i = 0; i < n; i++) {
            printf("%f ", x[i]);
        }
        printf("\n");
    } else {
        printf("No unique solution exists.\n");
    }

    //Free allocated memory
    for (i = 0; i < n; i++) {
        free(a[i]);
    }
    free(a);
    free(b);
    free(x);
    return 0;
}
