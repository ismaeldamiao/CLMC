/* *****************************************************************************
   Function to compute a set of N correlated numbers using a

   V_i = 

   Where
   *****************************************************************************
   E-mail: ismaellxd@gmail.com
   Site: https://ismaeldamiao.github.io/
   *****************************************************************************
   Copyright (c) 2021 I.F.F. dos SANTOS (Ismael Damiao)

   Permission is hereby granted, free of charge, to any person obtaining a copy 
   of this software and associated documentation files (the “Software”), to 
   deal in the Software without restriction, including without limitation the 
   rights to use, copy, modify, merge, publish, distribute, sublicense, and/or 
   sell copies of the Software, and to permit persons to whom the Software is 
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in 
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS 
   IN THE SOFTWARE.
***************************************************************************** */
#include "./damiao.h"

double *correlated_w_distance(double alpha, int N, int seed){
#if defined(UINT64_MAX)
   uint64_t idum;
#elif defined(UINT32_MAX)
   uint32_t idum;
#else
   unsigned int idum;
#endif
   double *V, *phi, aux0, aux1, aux2, deviation;
   int i, j, k;
   
   phi = (double*)malloc((size_t)((N) * sizeof(double)));
   V = (double*)malloc((size_t)((N) * sizeof(double)));
   
   aux1 = aux2 = 0.0;
   idum = seed;

   for(i = 0; i < N; ++i) phi[i] = 2.0 * random(&idum) - 1.0;

   for(i = 1; i <= N; ++i){
      k = i-1;
      V[k] = 0.0;
      for(j = 1; j <= N; ++j){
         aux0 = (double)abs(i - j) / alpha + 1.0;
         V[k] += phi[j-1] / (aux0*aux0);
      }
      aux1 += V[k];
      aux2 += V[k] * V[k];
   }

   /* Normalize the sequence */
   aux1 /= ((double)N);
   aux2 /= ((double)N);
   deviation = sqrt(aux2 - aux1*aux1);
   for(i = 0; i < N; ++i) V[i] = (V[i] - aux1) / deviation;

   free(phi);

   return V;
}
