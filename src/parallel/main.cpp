#include <math.h>

#include <iostream>

#include "omp.h"
#include "utils.hpp"

// #define EXPORT  // shows results on CSV
// #define DEBUG   // shows results on screen

#define NUM_ITERATIONS 3

#define NUMBER_SPRINGS  999
#define SPRING_CONSTANT 1.
#define INITIAL_MASS    1.

#define TIME_LIMIT 10.
#define TIME_STEPS 0.1

double** initializeMatrix(int rows, int cols);
double* initializeMass(int size, double initial_m, double initial_n);

void printMatrix(double** matrix, int rows, int cols);
void printVector(double* vector, int size);

void terminateVector(double* vector);
void terminateMatrix(double** matrix, int size);

#if NUMBER_SPRINGS % 2 == 0
#    error NUMBER_SPRINGS must be odd
#endif

/* ---------------------- Main ---------------------- */
int main(int argc, char** argv) {
#ifdef EXPORT
    FILE* fout;
#endif

    // const int N = NUMBER_SPRINGS;

    const int n0    = 10;  // Initial mass constant
    const double D  = 1.;  // Spring constant
    const double m0 = 1.;  // Initial mass

    const double T  = TIME_LIMIT;  // time limit
    const double dT = TIME_STEPS;  // time slots

    double* e_val;        // vector diagonal
    double* subdiagonal;  // vector subdiagonal
    double* m;            // Vector mass
    double* X;            // Vector position
    double** e_vec;       // Matrix temporal
    double** matrix;      // Matrix to calculate

    double t1, t2, t3, t4;

    int i, j;

#ifdef EXPORT
    fout = fopen("results.csv", "w");
    if (fout == NULL) {
        perror("results.csv");
        exit(-1);
    }
#endif

    for (int it = 0; it < NUM_ITERATIONS; it++) {
        for (int numThreads = 4; numThreads <= 16; numThreads += 4) {
            omp_set_num_threads(numThreads);

            for (int N = 99; N < 1000; N += 100) {
                int rows = N;
                int cols = N;

                t1 = omp_get_wtime();

                X      = new double[N];
                m      = initializeMass(N, m0, n0);
                matrix = initializeMatrix(rows, cols);

#pragma omp parallel for default(shared) private(i)
                for (i = 0; i < rows; ++i) {
#pragma omp parallel for default(shared) private(j)
                    for (j = 0; j < cols; ++j) {
                        matrix[i][j] = D * (2 * delta(i, j) - delta(i, j + 1) - delta(i, j - 1)) /
                                       sqrt(m[i] * m[j]);
                    }
                }

                // z init
                e_vec = initializeMatrix(rows, cols);

#pragma omp parallel for default(shared) private(i)
                for (i = 0; i < rows; ++i) {
#pragma omp parallel for default(shared) private(j)
                    for (j = 0; j < cols; ++j) {
                        e_vec[i][j] = delta(i, j);
                    }
                }

                // diagonal and subdiagonal init
                e_val       = new double[N];
                subdiagonal = new double[N - 1];

                t2 = omp_get_wtime();

                // Pass tred2 algorithm. For evaluation, not necessarily
                tred2(matrix, N, e_val, subdiagonal);

                // printMatrix(e_vec, rows, cols);

                // Apply tqli algorithm
                tqli(e_val, subdiagonal, N, e_vec);
                t3 = omp_get_wtime();

#ifdef DEBUG
                printMatrix(e_vec, rows, cols);
#endif

                for (double t = 0; t < T; t = t + dT) {  // Replace values in equation of X(t)
#pragma omp parallel for default(shared) private(i)
                    for (i = 0; i < N; i++) {  // Define X[0]
                        X[i] = 10. * double(i);
                    }

                    //#pragma omp parallel for default(shared)// private(t)
                    for (i = 0; i < rows; i++) {
                        double res = 0;
#pragma omp parallel for default(shared) private(j) reduction(+ : res)
                        for (j = 0; j < cols; j++) {
                            // printf("%d, %d, %f, %d \n", i, j, t, omp_get_thread_num());
                            res +=
                                e_vec[i][j] * cos(e_val[j] * t) + e_vec[i][j] * sin(e_val[j] * t);
                        }
                        X[i] += res;

#ifdef EXPORT
                        fprintf(fout, "%lf,", X[i]);
#endif
                    }
                    // printVector(X, rows);
#ifdef EXPORT
                    fprintf(fout, "\n");
#endif
                }

#ifdef DEBUG
                printMatrix(e_vec, rows, cols);
#endif

                t4 = omp_get_wtime();
                printf("%9.6f, %9.6f, %9.6f, %9.6f\n", (t2 - t1) * 1000, (t3 - t2) * 1000,
                       (t4 - t3) * 1000, (t4 - t1) * 1000);

                terminateMatrix(e_vec, N);
                terminateMatrix(matrix, N);
                terminateVector(e_val);
                terminateVector(subdiagonal);
                terminateVector(m);
                terminateVector(X);
            }
        }
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

void printVector(double* vector, int size) {
    for (int i = 0; i < size; i++) {
        printf("%4.2f ", vector[i]);
    }
    printf("\n");
}