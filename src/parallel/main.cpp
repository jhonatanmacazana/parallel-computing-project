#include <math.h>

#include <iostream>

#include "utils.hpp"

#define EXPORT
#define NUMBER_SPRINGS  99
#define SPRING_CONSTANT 1.
#define INITIAL_MASS    1.

double** initializeMatrix(int rows, int cols);
double* initializeMass(int size, double initial_m, double initial_n);

void printMatrix(double** matrix, int rows, int cols);

void terminateVector(double* vector);
void terminateMatrix(double** matrix, int size);

/* ---------------------- Main ---------------------- */
int main(int argc, char** argv) {
#ifdef EXPORT
    FILE* fout;
#endif

    const int N = NUMBER_SPRINGS;

    const int rows = N;
    const int cols = N;

    const int n0    = 10;  // Initial mass constant
    const double D  = 1.;  // Spring constant
    const double m0 = 1.;  // Initial mass

    const double T  = 10.;   // time limit
    const double dT = 0.01;  // time slots

    double* e_val;        // vector diagonal
    double* subdiagonal;  // vector subdiagonal
    double* m;            // Vector mass
    double* X;            // Vector position
    double** e_vec;       // Matrix temporal
    double** matrix;      // Matrix to calculate

#ifdef EXPORT
    fout = fopen("results.csv", "w");
    if (fout == NULL) {
        perror("results.csv");
        exit(-1);
    }
#endif

    X      = new double[N];
    m      = initializeMass(N, m0, n0);
    matrix = initializeMatrix(rows, cols);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            matrix[i][j] =
                D * (2 * delta(i, j) - delta(i, j + 1) - delta(i, j - 1)) / sqrt(m[i] * m[j]);
        }
    }

    // z init
    e_vec = initializeMatrix(rows, cols);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            e_vec[i][j] = delta(i, j);
        }
    }


    // diagonal and subdiagonal init
    e_val       = new double[N];
    subdiagonal = new double[N - 1];

    // Pass tred2 algorithm. For evaluation, not necessarily
    tred2(matrix, N, e_val, subdiagonal);


    // for (int i=0; i<N-1; i++)
    //     printf("%f ", e_vec[i][i]);
    // printf("\n");
    // Apply tqli algorithm
    tqli(e_val, subdiagonal, N, e_vec);
    //printMatrix(e_vec, rows, cols);

    // for (int i=0; i<N-1; i++)
    //     printf("%f ", e_vec[i][i]);

    for (double t = 0; t < T; t = t + dT) {  // Replace values in equation of X(t)
        for (int i = 0; i < N; i++) {        // Define X[0]
            X[i] = 10. * double(i);
        }

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                X[i] += e_vec[i][j] * cos(e_val[j] * t) + e_vec[i][j] * sin(e_val[j] * t);
            }

#ifdef EXPORT
            fprintf(fout, "%lf,", X[i]);
#endif
        }
#ifdef EXPORT
        fprintf(fout, "\n");
#endif
    }

    terminateMatrix(e_vec, N);
    return 0;
}

double** initializeMatrix(int rows, int cols) {
    double** output = new double*[rows];
    for (int i = 0; i < rows; ++i) {
        output[i] = new double[cols];
    }
    return output;
}

double* initializeMass(int size, double initial_m, double initial_n) {
    double* output = new double[size];
    for (int i = 0; i < size; ++i) {
        output[i] = initial_m;
        if (i == initial_n) {
            output[i] = 100 * initial_m;
        }
    }
    return output;
}

void terminateVector(double* vector) { free(vector); }

void terminateMatrix(double** matrix, int size) {
    for (int i = 0; i < size; i++) {
        double* tempPointer = matrix[i];
        free(tempPointer);
    }
}

void printMatrix(double** matrix, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%4.2f ", matrix[i][j]);
        }
        printf("\n");
    }
}