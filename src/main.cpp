#include <math.h>

#include <iostream>

using namespace std;

/* ---------------------- Definition ---------------------- */
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/* ---------------------- Prototypes ---------------------- */
float pythag(float a, float b);
void tqli(float* d, float* e, int n, float** z);

/* ---------------------- Main ---------------------- */
int main(int argc, char** argv) {
    printf("Hello world\n");
    return 0;
}

/* ---------------------- Functions ---------------------- */

/*QL algorithm with implicit shifts, to determine the eigenvalues and eigenvectors of a real,
symmetric, tridiagonal matrix, or of a real, symmetric matrix previously reduced by tred2 §11.2. On
input, d[1..n] contains the diagonal elements of the tridiagonal matrix. On output, it returns
the eigenvalues. The vector e[1..n] inputs the subdiagonal elements of the tridiagonal matrix,
with e[1] arbitrary. On output e is destroyed. When finding only the eigenvalues, several lines
may be omitted, as noted in the comments. If the eigenvectors of a tridiagonal matrix are desired,
the matrix z[1..n][1..n] is input as the identity matrix. If the eigenvectors of a matrix
that has been reduced by tred2 are required, then z is input as the matrix output by tred2.
In either case, the kth column of z returns the normalized eigenvector corresponding to d[k].*/
float pythag(float a, float b) {
    float absa, absb;
    absa = fabs(a);
    absb = fabs(b);
    if (absa > absb)
        return absa * sqrt(1.0 + (absb / absa) * (absb / absa));
    else
        return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + (absa / absb) * (absa / absb)));
}

void tqli(float* d, float* e, int n, float** z) {
    register int m, l, iter, i, k;
    float s, r, p, g, f, dd, c, b;

    for (i = 1; i < n; i++) e[i - 1] = e[i];
    e[n] = 0.0;
    for (l = 0; l < n; l++) {
        iter = 0;
        do {
            for (m = l; m < n - 1; m++) {
                dd = fabs(d[m]) + fabs(d[m + 1]);
                if ((double)(fabs(e[m]) + dd) == dd) break;
            }
            if (m != l) {
                if (iter++ == 30) {
                    cout << "\n\nToo many iterations in tqli.\n";
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
                }     /* end i-loop */
                if (r == 0.0 && i >= l) continue;
                d[l] -= p;
                e[l] = g;
                e[m] = 0.0;
            } /* end if-loop for m != 1 */
        } while (m != l);
    }
}