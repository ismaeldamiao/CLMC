/* *****************************************************************************
   Function to compute a set of N correlated numbers using a Bernoulli map
   defined as

   X_i =

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

double *correlated_w_bernoulli(double alpha, int N, int seed){
#if defined(UINT64_MAX)
   uint64_t idum;
#elif defined(UINT32_MAX)
   uint32_t idum;
#else
   unsigned int idum;
#endif
   int i;
   double aux, *X, b, caux;

   aux = 0.0;
   b = 1.0e-12;
   caux = pow(2.0, alpha-1.0) * (1.0 - 2.0 * b);
   idum = seed;

   X = (double*)malloc((size_t)((N) * sizeof(double)));

   X[0] = random(&idum);
   for(i = 1; i < N; ++i){
      if((X[i-1] >= 0.0) && (X[i-1] < 0.5))
         aux = pow(X[i-1], alpha);
      else if(X[i-1] >= 0.5)
         aux = -pow(1.0 - X[i-1], alpha);
      X[i] = X[i-1] + caux * aux + b;
   }
   return X;
}
