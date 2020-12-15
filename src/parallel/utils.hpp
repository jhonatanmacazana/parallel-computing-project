#ifndef UTILS_H_
#define UTILS_H_

#include <math.h>
#include <omp.h>

#include <iostream>
#include <limits>

/* ---------------------- Definition ---------------------- */
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/* ---------------------- Prototypes ---------------------- */

double pythag(double a, double b);

/*
 *   QL algorithm with implicit shifts, to determine the eigenvalues and eigenvectors of a real,
 *   symmetric, tridiagonal matrix, or of a real, symmetric matrix previously reduced by tred2
 * §11.2.
 *   @inputs:
 *       double* d       diagonal elements of the tridiagonal matrix
 *       double* e       subdiagonal elements of the tridiagonal matrix. e[1] arbitrary
 *       int n           size of square matrix (1..n)
 *       double** z      Z matrix, could be the output of tred2
 *   @outputs
 *       double* d       eigenvalues
 *       double* e       destroyed
 *       double** z      kth column is the normalized eigenvector corresponding to d[k]
 */
void tqli(double* d, double* e, int n, double** z);

/*
 *   delta() function
 *
 *   @inputs:
 *       int a       num1 to test
 *       int b       num2 to test
 *
 *   @outputs
 *       double          1 if a == b. Otherwise, 0
 *
 */
double delta(int a, int b);

/*
 *     Performs a Housholder reduction of a real symmetric matrix
 *     a[][]. On output a[][] is replaced by the orthogonal matrix
 *     effecting the transformation. d[] returns the diagonal elements
 *     of the tri-diagonal matrix, and e[] the off-diagonal elements,
 *     with e[0] = 0.
 *     The function is modified from the version in Numerical recipe.
 * */
void tred2(double** a, int n, double* d, double* e);

void tred2_straight(double** a, int n, double* d, double* e);

/* ---------------------- Functions ---------------------- */

// Impulse function
double delta(int a, int b) { return (a == b) ? 1. : 0.; }

void tred2_straight(double** a, int n, double* d, double* e) {
    for (int i = 0; i < n - 1; i++) {
        d[i] = a[i][i];
        e[i] = a[i + 1][i];
    }
    d[n - 1] = a[n - 1][n - 1];
}

void tred2(double** a, int n, double* d, double* e) {
    int l, k, j, i;
    double scale, hh, h, g, f;

    for (i = n - 1; i > 0; i--) {
        l = i - 1;
        h = scale = 0.0;
        if (l > 0) {
            for (k = 0; k < l + 1; k++) scale += fabs(a[i][k]);
            if (scale == 0.0)
                e[i] = a[i][l];
            else {
                for (k = 0; k < l + 1; k++) {
                    a[i][k] /= scale;
                    h += a[i][k] * a[i][k];
                }
                f    = a[i][l];
                g    = (f >= 0.0 ? -sqrt(h) : sqrt(h));
                e[i] = scale * g;
                h -= f * g;
                a[i][l] = f - g;
                f       = 0.0;
                for (j = 0; j < l + 1; j++) {
                    // Next statement can be omitted if eigenvectors not wanted
                    a[j][i] = a[i][j] / h;
                    g       = 0.0;
                    for (k = 0; k < j + 1; k++) g += a[j][k] * a[i][k];
                    for (k = j + 1; k < l + 1; k++) g += a[k][j] * a[i][k];
                    e[j] = g / h;
                    f += e[j] * a[i][j];
                }
                hh = f / (h + h);
                for (j = 0; j < l + 1; j++) {
                    f    = a[i][j];
                    e[j] = g = e[j] - hh * f;
                    for (k = 0; k < j + 1; k++) {
                        a[j][k] -= (f * e[k] + g * a[i][k]);
                    }
                }
            }
        } else
            e[i] = a[i][l];
        d[i] = h;
    }
    // Next statement can be omitted if eigenvectors not wanted
    d[0] = 0.0;
    e[0] = 0.0;

    // Contents of this loop can be omitted if eigenvectors not
    //	wanted except for statement d[i]=a[i][i];
    for (i = 0; i < n; i++) {
        l = i;
        if (d[i] != 0.0) {
            for (j = 0; j < l; j++) {
                g = 0.0;
                for (k = 0; k < l; k++) g += a[i][k] * a[k][j];
                for (k = 0; k < l; k++) a[k][j] -= g * a[k][i];
            }
        }
        d[i]    = a[i][i];
        a[i][i] = 1.0;

        for (j = 0; j < l; j++) {
            a[j][i] = a[i][j] = 0.0;
        }
    }
}

double pythag(double a, double b) {
    double absa, absb;
    absa = fabs(a);
    absb = fabs(b);
    if (absa > absb)
        return absa * sqrt(1.0 + (absb / absa) * (absb / absa));
    else
        return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + (absa / absb) * (absa / absb)));
}

void tqli(double* d, double* e, int n, double** z) {
    register int m, l, iter, i, k, inner_3_iter;
    double s, r, p, g, f, dd, c, b;
    volatile bool flag_inner_1, flag_inner_2, flag_inner_3;
    int max_iter     = 100;
    const double EPS = std::numeric_limits<double>::epsilon();

    for (i = 1; i < n; i++) {
        e[i - 1] = e[i];
    }

    e[n - 1] = 0.0;

    //#pragma omp parallel for default (shared) private(l, iter, d) reduction(-:d[:n])
    for (l = 0; l < n; l++) {
        iter         = 0;
        flag_inner_3 = true;
        // do {
        //#pragma omp parallel for shared(flag_inner_3, inner_3_iter, iter)
        for (inner_3_iter = 0; inner_3_iter < max_iter; inner_3_iter++)
            if (flag_inner_3) {
                flag_inner_2 = true;
                //#pragma omp parallel for default(shared )
                // for (m = l; m < n - 1; m++) {
                //     if (flag_inner_2){
                //         dd = fabs(d[m]) + fabs(d[m + 1]);
                //         if ((double)(fabs(e[m]) + dd) == dd){
                //             flag_inner_2 = false;
                //         }
                //     }
                // }

                for (m = l; m < n - 1; m++) {
                    dd = abs(d[m]) + abs(d[m + 1]);
                    if (abs(e[m]) <= EPS * dd) break;
                }

                if (m != l) {
                    if (iter++ == 30 || inner_3_iter > max_iter) {
                        // perror("\n\nToo many iterations in tqli.\n");
                        exit(0);
                    }

                    g = (d[l + 1] - d[l]) / (2.0 * e[l]);
                    r = pythag(g, 1.0);
                    g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
                    s = c = 1.0;
                    p     = 0.0;

                    flag_inner_1 = true;

                    //#pragma omp parallel for shared(e, d, flag_inner_1)// reduction(-:d[:n])
                    for (i = m - 1; i >= l; i--) {
                        if (flag_inner_1) {
                            f        = s * e[i];
                            b        = c * e[i];
                            e[i + 1] = (r = pythag(f, g));

                            if (r == 0.0) {
                                d[i + 1] -= p;
                                e[m]         = 0.0;
                                flag_inner_1 = false;
                            }
                            //#pragma omp barrier

                            if (flag_inner_1) {
                                s        = f / r;
                                c        = g / r;
                                g        = d[i + 1] - p;
                                r        = (d[i] - g) * s + 2.0 * c * b;
                                d[i + 1] = g + (p = s * r);
                                g        = c * r - b;
                                // #pragma omp parallel for default(shared) private(k, f)
                                for (k = 0; k < n; k++) {
                                    f           = z[k][i + 1];
                                    z[k][i + 1] = s * z[k][i] + c * f;
                                    z[k][i]     = c * z[k][i] - s * f;
                                } /* end k-loop */
                            }
                        }
                    } /* end i-loop */
                    if (!(r == 0.0 && i >= l)) {
                        d[l] -= p;
                        e[l] = g;
                        e[m] = 0.0;
                    }

                } /* end if-loop for m != 1 */
                if (m != l) flag_inner_3 = false;
            }
        //} while (m != l);
    }
}

#endif  // UTILS_H_