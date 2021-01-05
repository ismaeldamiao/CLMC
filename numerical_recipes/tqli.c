/*******************************************************************************
Eigenvalue solvers, tred2 and tqli, from "Numerical Recipes in C" (Cambridge
Univ. Press) by W.H. Press, S.A. Teukolsky, W.T. Vetterling, and B.P. Flannery
*******************************************************************************/
#include<math.h>
/*******************************************************************************
Computes (a**2 + b**2)**1/2 without destructive underflow or overflow.
*******************************************************************************/
double pythag(double a, double b){
   double absa = fabs(a), absb = fabs(b), c;
   if(absa > absb){
      c = absb/absa;
      c = absa * sqrt(1.0 + c * c);
   }else{
      c = absa/absb;
      c = (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + c * c));
   }
   return c;
}
#define SIGN(a, b) ((b) < 0.0 ? -fabs(a) : fabs(a))
/*******************************************************************************
Householder reduction of a real, symmetric matrix a[1..n][1..n]. 
On output, a is replaced by the orthogonal matrix Q effecting the
transformation. d[1..n] returns the diagonal elements of the tridiagonal matrix,
and e[1..n] the off-diagonal elements, with e[1]=0. Several statements, as noted
in comments, can be omitted if only eigenvalues are to be found, in which case a
contains no useful information on output. Otherwise they are to be included.
*******************************************************************************/
void tred2(double **a, int n, double d[], double e[]){
   int l, k, j, i;
   double scale, hh, h, g, f;

   for(i = n; i >= 2; i--){
      l = i - 1;
      h = scale = 0.0;
      if(l > 1){
         for(k = 1; k <= l; k++)
         scale += fabs(a[i][k]);
         if(scale == 0.0) /* Skip transformation. */
         e[i] = a[i][l];
         else{
            for(k = 1; k <= l; k++){
               a[i][k] /= scale; /* Use scaled a's for transformation. */
               h += a[i][k] * a[i][k]; /* Form sigma in h. */
            }
            f = a[i][l];
            g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
            e[i] = scale * g;
            h -= f*g; /* Now h is equation (11.2.4). */
            a[i][l] = f - g; /* Store u in the ith row of a. */
            f = 0.0;
            for(j = 1; j <= l; j++) {
               /* Next statement can be omitted if eigenvectors not wanted */
               a[j][i] = a[i][j]/h; /* Store u/H in ith column of a. */
               g = 0.0; /* Form an element of A.u in g. */
               for(k = 1; k <= j; k++)
               g += a[j][k] * a[i][k];
               for(k = j+1; k <= l; k++)
               g += a[k][j] * a[i][k];
               e[j]=g/h;/*Form element of p in temporarily unused element of e*/
               f += e[j] * a[i][j];
            }
            hh = f/(h+h); /* Form K, equation (11.2.11). */
            for(j = 1; j <= l; j++) { /* Form q and store in e overwriting p. */
               f=a[i][j];
               e[j] = g = e[j] - hh*f;
               for(k = 1; k <= j; k++) /* Reduce a, equation (11.2.13). */
               a[j][k] -= (f*e[k]+g*a[i][k]);
            }
         }
      }else
      e[i]=a[i][l];
      d[i]=h;
   }
   /* Next statement can be omitted if eigenvectors not wanted */
   d[1] = 0.0;
   e[1] = 0.0;
   /* Contents of this loop can be omitted if eigenvectors not
   wanted except for statement d[i]=a[i][i]; */
   for(i = 1; i <= n; i++){ /* Begin accumulation of transformation matrices. */
      l = i - 1;
      if(d[i]){ /* This block skipped when i=1. */
         for(j = 1; j <= l; j++){
            g = 0.0;
            for(k = 1; k <= l; k++) /* Use u and u/H stored in a to form P.Q. */
            g += a[i][k]*a[k][j];
            for(k = 1; k <= l; k++)
            a[k][j] -= g * a[k][i];
         }
      }
      d[i] = a[i][i]; /* This statement remains. */
      a[i][i]=1.0;/*Reset row and column of a to identity matrix for next iteration*/
      for(j = 1; j <= l; j++) a[j][i] = a[i][j] = 0.0;
   }
}
/*******************************************************************************
QL algorithm with implicit shifts, to determine the eigenvalues and eigenvectors
of a real, symmetric, tridiagonal matrix, or of a real, symmetric matrix
previously reduced by tred2 sec. 11.2. On input, d[1..n] contains the diagonal
elements of the tridiagonal matrix. On output, it returns the eigenvalues. The
vector e[1..n] inputs the subdiagonal elements of the tridiagonal matrix, with
e[1] arbitrary. On output e is destroyed. When finding only the eigenvalues,
several lines may be omitted, as noted in the comments. If the eigenvectors of
a tridiagonal matrix are desired, the matrix z[1..n][1..n] is input as the
identity matrix. If the eigenvectors of a matrix that has been reduced by tred2
are required, then z is input as the matrix output by tred2. In either case,
the kth column of z returns the normalized eigenvector corresponding to d[k].
*******************************************************************************/
int tqli(double *d, double *e, int n, double **z){
   int m, l, iter, i, k;
   double s, r, p, g, f, dd, c, b;
   /* Convenient to renumber the elements of e. */
   for(i = 2; i <= n; i++) e[i-1] = e[i];
   e[n] = 0.0;
   for(l = 1; l <= n; l++){
      iter = 0;
      do{
			for(m = l; m <= n-1; m++){
				dd = fabs(d[m]) + fabs(d[m+1]);
				if((fabs(e[m]) + dd) == dd) break;
			}
         if(m != l){
            if(iter++ == 30) return 200;
            g = (d[l+1] - d[l]) / (2.0 * e[l]);
            r = sqrt(g*g + 1.0);
            g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
            s = c = 1.0;
            p = 0.0;
            /* A plane rotation as in the original QL, followed by Givens */
            for(i = m-1; i >= l; i--){
               f = s * e[i];/* rotations to restore tridiagonal form. */
               b = c * e[i];
               e[i+1] = r = pythag(f, g);
               if(r == 0.0) { /* Recover from underflow. */
                  d[i+1] -= p;
                  e[m] = 0.0;
                  break;
               }
               s = f / r;
               c = g / r;
               g = d[i+1] - p;
               r = (d[i] - g) * s + 2.0 * c * b;
               d[i+1] = g + (p = s*r);
               g = c * r - b;
               /* Next loop can be omitted if eigenvectors not wanted*/
               for(k = 1; k <= n; k++){ /* Form eigenvectors. */
                  f = z[k][i+1];
                  z[k][i+1] = s * z[k][i] + c * f;
                  z[k][i] = c * z[k][i] - s * f;
               }
            }
            if(r == 0.0 && i >= l) continue;
            d[l] -= p;
            e[l] = g;
            e[m] = 0.0;
         }
      } while(m != l);
   }
   return 0;
}
#undef SIGN
