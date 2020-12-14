#include <math.h>

#include <iostream>

#include "utils.hpp"

#define EXPORT
#define NUMBER_SPRINGS  99
#define SPRING_CONSTANT 1.
#define INITIAL_MASS    1.

double** initializeMatrix(int rows, int cols);
double* initializeMass(int size, double m0, double n0);
void terminateVector(double* vector);
void terminateMatrix(double** matrix, int size);

/* ---------------------- Main ---------------------- */
int main(int argc, char** argv) {
#ifdef EXPORT
    FILE* fout;
#endif

    double** matrix;
    int N    = NUMBER_SPRINGS;
    int n0   = 10;
    double T = 10., dT = 0.01;

    double D  = SPRING_CONSTANT;
    double m0 = INITIAL_MASS;

    double* m;

    double* diagonal;
    double* X;
    double* subdiagonal;
    double** z;

    int rows = N;
    int cols = N;

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
    z = initializeMatrix(rows, cols);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            z[i][j] = delta(i, j);
        }
    }

    // diagonal and subdiagonal init
    diagonal    = new double[N];
    subdiagonal = new double[N - 1];

    // Pass tred2 algorithm. For evaluation, not necessarily
    tred2(matrix, N, diagonal, subdiagonal);

    // Apply tqli algorithm
    // tqli(diagonal, subdiagonal, N, z);

    for (double t = 0; t < T; t = t + dT) {  // Replace values in equation of X(t)
        for (int i = 0; i < N; i++) {        // Define X[0]
            X[i] = 10. * double(i);
        }

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                X[i] += z[i][j] * cos(diagonal[j] * t) + z[i][j] * sin(diagonal[j] * t);
            }

#ifdef EXPORT
            fprintf(fout, "%lf,", X[i]);
#endif
        }
#ifdef EXPORT
        fprintf(fout, "\n");
#endif
    }

    return 0;
}

double** initializeMatrix(int rows, int cols) {
    double** output = new double*[rows];
    for (int i = 0; i < rows; ++i) {
        output[i] = new double[cols];
    }
    return output;
}

double* initializeMass(int size, double m0, double n0) {
    double* output = new double[size];
    for (int i = 0; i < size; ++i) {
        output[i] = m0;
        if (i == n0) {
            output[i] = 100 * m0;
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