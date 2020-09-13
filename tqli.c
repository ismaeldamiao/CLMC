/*******************************************************************************
Eigenvalue solvers, tred2 and tqli, from "Numerical Recipes in C" (Cambridge
Univ. Press) by W.H. Press, S.A. Teukolsky, W.T. Vetterling, and B.P. Flannery
*******************************************************************************/
#include<math.h>
/*******************************************************************************
Computes (a**2 + b**2)**1/2 without destructive underflow or overflow.
*******************************************************************************/
double pythag(double a, double b){
   double absa = fabs(a), absb = fabs(b), c, tmp;
   if(absa > absb){
      tmp = absb/absa;
      c = absa * sqrt(1.0 + tmp * tmp);
   }else{
      tmp = absa/absb;
      c = (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + tmp * tmp));
   }
   return c;
}
#define SIGN(a, b) ((b) < 0.0 ? -fabs(a) : fabs(a))
/*******************************************************************************
QL algorithm with implicit shifts, to determine the eigenvalues and eigenvectors
of a real, symmetric, tridiagonal matrix, or of a real, symmetric matrix
previously reduced by tred2 sec. 11.2. On input, d[0..n-1] contains the diagonal
elements of the tridiagonal matrix. On output, it returns the eigenvalues. The
vector e[0..n-1] inputs the subdiagonal elements of the tridiagonal matrix, with
e[0] arbitrary. On output e is destroyed. When finding only the eigenvalues,
several lines may be omitted, as noted in the comments. If the eigenvectors of
a tridiagonal matrix are desired, the matrix z[0..n-1][0..n-1] is input as the
identity matrix. If the eigenvectors of a matrix that has been reduced by tred2
are required, then z is input as the matrix output by tred2. In either case,
the kth column of z returns the normalized eigenvector corresponding to d[k].
*******************************************************************************/
int tqli(double *d, double *e, int n, double **z){
   int m, l, iter, i, k;
   double s, r, p, g, f, dd, c, b;
   /* Convenient to renumber the elements of e. */
   for(i = 1; i < n; i++) e[i-1] = e[i];
   e[n-1] = 0.0;
   for(l = 0; l < n; l++){
      iter = 0;
      do{
			for(m = l; m < n-1; m++){
				dd = fabs(d[m]) + fabs(d[m+1]);
				if((fabs(e[m]) + dd) == dd) break;
			}
         if(m != l){
            if(iter++ == 30) return 1;
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
               for(k = 0; k < n; k++){ /* Form eigenvectors. */
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
