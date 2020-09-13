/* *****************************************************************************
   Funcao que calcula a dinamica do sistema no tempo t=0
   *****************************************************************************
   E-mail: ismaellxd@gmail.com
   Site: https://ismaeldamiao.github.io/
   *****************************************************************************
   Copyright © 2020 Ismael Damião

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
#include"CLMC.h"
#include<math.h>

double **acoplamento(double eta2, double eta3, double eta4, int N){
   double **eta = dmatriz(N+2, 3);
   
   for(int n = 1; n <= N; ++n){
      eta[n][0] = eta2;
      eta[n][1] = eta3;
      eta[n][2] = eta4;
   }
   return eta;
}

double *DefPosicaoInicial(const int N){
   double *x = dvetor(N+2);
   for(int n = 1; n <= N; ++n) x[n] = 0.0;
   return x;
}

#ifndef DeltaDeKronecker
#define DeltaDeKronecker(x, y) ((x) != (y) ? 0.0 : 1.0)
#endif
double *DefMomentoInicial(int N, double *M, double v0){
   double *P = dvetor(N+2);
   const int N2 = N/2;
   for(int n = 1; n <= N; ++n) P[n] = M[n] * v0 * DeltaDeKronecker(n, N2);
   return P;
}

double CalcEnergiaInicial(int N, double *M, double *x, double *P, double **eta){
   const double c6 = 1.0/6.0;
   double H0 = 0.0, aux[2];
   for(int n = 1, i = 2; n <= N; ++n, ++i){
      aux[0] = x[n+1] - x[n];
      aux[1] = x[n] - x[n-1];
      H0 += 0.5 * pow(P[n], 2.0) / M[n] +\
      (eta[i][0] * pow(aux[0], 2.0) + eta[n][0] * pow(aux[1], 2.0)) * 0.25 +
      (eta[i][1] * pow(aux[0], 3.0) + eta[n][1] * pow(aux[1], 3.0)) * c6 +
      (eta[i][2] * pow(aux[0], 4.0) + eta[n][2] * pow(aux[1], 4.0)) * 0.125;
   }
   return H0;
}
