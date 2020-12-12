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
    // Initialize model: Create tridiagonal matrix 
    double ** matrix;
    int dim = 100;
    
    // Pass tred2 algorithm. For evaluation, not necessarily

    // Read values diagonal and subdiagonal. Create a identity matrix z 
    
    // Apply tqli algorithm

    // Read eigenvalues (d) and eigenvectors (z)

    // Replace values in equation of X(t)

    printf("Hello world\n");
    return 0;
}

/* ---------------------- Functions ---------------------- */
/*
*     ** The function
*     **                tred2()
*     ** perform a Housholder reduction of a real symmetric matrix
*     ** a[][]. On output a[][] is replaced by the orthogonal matrix 
*     ** effecting the transformation. d[] returns the diagonal elements
*     ** of the tri-diagonal matrix, and e[] the off-diagonal elements, 
*     ** with e[0] = 0.
*     ** The function is modified from the version in Numerical recipe.
*     */

void tred2(float **a, int n, float *d, float *e)
{
	   register int    l,k,j,i;
	   float          scale,hh,h,g,f;

	   for(i = n - 1; i > 0; i--) {
	       l = i-1;
	       h = scale= 0.0;
	           if(l > 0) {
		            for(k = 0; k <= l; k++)
		                scale += fabs(a[i][k]);
		                if(scale == 0.0)               // skip transformation
		               e[i] = a[i][l];
			            else {
	               for(k = 0; k <= l; k++) {
	               a[i][k] /= scale;          // used scaled a's for transformation
	                      h       += a[i][k]*a[i][k];
		                  }
		            f       = a[i][l];
	                g       = (f >= 0.0 ? -sqrt(h) : sqrt(h));
		            e[i]    = scale*g;
	                h      -= f * g;
	            a[i][l] = f - g;
	                f       = 0.0;

	            for(j = 0;j <= l;j++) {
	                   a[j][i] = a[i][j]/h;       // can be omitted if eigenvector not wanted
	                  g       = 0.0; 
	                 for(k = 0; k <= j; k++) {
	                   g += a[j][k]*a[i][k];
		                  }
	             for(k = j+1; k <= l; k++)
	                  g += a[k][j]*a[i][k];
		               e[j]=g/h;
	                      f += e[j]*a[i][j];
		                  }
               hh=f/(h+h);
	            for(j = 0; j <= l;j++) {
	                   f = a[i][j];
                 e[j]=g=e[j]-hh*f;
                 for(k = 0; k <= j; k++)
	                   a[j][k] -= (f*e[k]+g*a[i][k]);
		             }
	             }  // end k-loop
	          }  // end if-loop for l > 1
		         else {
	          e[i]=a[i][l];
		        }
		       d[i]=h;
	          }  // end i-loop
	    d[0]  = 0.0;
	    e[0]  = 0.0;

     /* Contents of this loop can be omitted if eigenvectors not
	 * 	 ** wanted except for statement d[i]=a[i][i];
	 * 	          */
          for(i = 0; i < n; i++) {
	        l = i-1;
	      if(d[i]) {
	               for(j = 0; j <= l; j++) {
	                   g= 0.0;
	               for(k = 0; k <= l; k++) {
                      g += a[i][k] * a[k][j];
	                  }
                   for (k = 0; k <= l; k++) {
                  a[k][j] -= g * a[k][i];
	              }
            }
         }
            d[i]    = a[i][i];
          a[i][i] = 1.0;
	        for(j = 0; j <= l; j++)  {
		         a[j][i]=a[i][j] = 0.0;
				       }
			   }
}

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