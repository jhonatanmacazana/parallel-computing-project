#ifndef UTILS_H_
#define UTILS_H_

#include <math.h>

#include <iostream>
#include <limits>

/* ---------------------- Definition ---------------------- */
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SQR(a)     ((a) * (a))

/* ---------------------- Prototypes ---------------------- */

double pythag(double a, double b);

/*
 *   Uses the QL algorithm with implicit shifts to determine the eigenvalues and eigenvectors of a
 * real, symmetric, tridiagonal matrix, or of a real, symmetric matrix previously reduced by tred2
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
            for (k = 0; k < i + 1; k++) scale += fabs(a[i][k]);
            if (scale == 0.0)
                e[i] = a[i][l];
            else {
                for (k = 0; k < i; k++) {
                    a[i][k] /= scale;
                    h += a[i][k] * a[i][k];
                }
                f    = a[i][l];
                g    = (f >= 0.0 ? -sqrt(h) : sqrt(h));
                e[i] = scale * g;
                h -= f * g;
                a[i][l] = f - g;
                f       = 0.0;
                for (j = 0; j < i; j++) {
                    // Next statement can be omitted if eigenvectors not wanted
                    a[j][i] = a[i][j] / h;
                    g       = 0.0;
                    for (k = 0; k < j + 1; k++) g += a[j][k] * a[i][k];
                    for (k = j + 1; k < i; k++) g += a[k][j] * a[i][k];
                    e[j] = g / h;
                    f += e[j] * a[i][j];
                }
                hh = f / (h + h);
                for (j = 0; j < i; j++) {
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
            for (j = 0; j < i; j++) {
                g = 0.0;
                for (k = 0; k < i; k++) g += a[i][k] * a[k][j];
                for (k = 0; k < i; k++) a[k][j] -= g * a[k][i];
            }
        }
        d[i]    = a[i][i];
        a[i][i] = 1.0;  // reset A to identity matrix

        for (j = 0; j < i; j++) {
            a[j][i] = a[i][j] = 0.0;
        }
    }
}

double pythag(double a, double b) {
    double absa = fabs(a);
    double absb = fabs(b);

    return (absa > absb ? absa * sqrt(1.0 + SQR(absb / absa))
                        : (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + SQR(absa / absb))));
}

void tqli(double* d, double* e, int n, double** z) {
    int m, l, iter, i, k;
    double s, r, p, g, f, dd, c, b;
    const double EPS = std::numeric_limits<double>::epsilon();

    for (i = 1; i < n; i++) e[i - 1] = e[i];

    e[n - 1] = 0.0;

    for (l = 0; l < n; l++) {
        iter = 0;
        do {
            for (m = l; m < n - 1; m++) {
                dd = abs(d[m]) + abs(d[m + 1]);
                if (abs(e[m]) <= EPS * dd) break;
            }
            if (m != l) {
                if (iter++ == 30) {
                    // perror("\n\nToo many iterations in tqli.\n");
                    exit(1);
                }

                g = (d[l + 1] - d[l]) / (2.0 * e[l]);
                r = pythag(g, 1.0);
                g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
                s = c = 1.0;
                p     = 0.0;
                for (i = m - 1; i >= l; i--) {
                    f        = s * e[i];
                    b        = c * e[i];
                    e[i + 1] = (r = pythag(f, g));

                    if (r == 0.0) {
                        d[i + 1] -= p;
                        e[m] = 0.0;
                        break;
                    }

                    s        = f / r;
                    c        = g / r;
                    g        = d[i + 1] - p;
                    r        = (d[i] - g) * s + 2.0 * c * b;
                    d[i + 1] = g + (p = s * r);
                    g        = c * r - b;

                    for (k = 0; k < n; k++) {
                        f           = z[k][i + 1];
                        z[k][i + 1] = s * z[k][i] + c * f;
                        z[k][i]     = c * z[k][i] - s * f;
                    } /* end k-loop */

                } /* end i-loop */
                if (r == 0.0 && i >= l) continue;
                d[l] -= p;
                e[l] = g;
                // printf("%d ", m);
                e[m] = 0.0;

            } /* end if-loop for m != 1 */
        } while (m != l);
    }
}

#endif  // UTILS_H_