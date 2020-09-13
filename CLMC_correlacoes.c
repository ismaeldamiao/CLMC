/* *****************************************************************************
   Funcoes de correlacao e funcao que define as massas conforme as correlacoes
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

#ifndef abs
   #define abs(x) ((x) >= 0 ? (x) : -(x))
#endif
/* *****************************************************************************
Geradores de correlacoes
***************************************************************************** */
double *MapaBernoulli(double alpha, long int semente, const int N){
   int i;
   double aux, *X;
   const double b = 1.0e-12, caux = pow(2.0, alpha-1.0) * (1.0 - 2.0 * b);
   
   X = dvetor(N+2); /* Alocar memoria */
   
   X[0] = ran1(&semente);
   for(i = 1; i < N; ++i){
      if((X[i-1] >= 0.0) && (X[i-1] < 0.5)){
         aux = pow(X[i-1], alpha);
      }else{
         if(X[i-1] >= 0.5){
            aux = -pow(1.0 - X[i-1], alpha);
         }
      }
      X[i] = X[i-1] + caux * aux + b;
   }
   return X;
}
double *SemNome1(int N, double alpha, long int semente){
   int Q2 = N/2;
   double *V;
   double phi[Q2], aux1 = 0.0, aux2 = 0.0, var, k, Nr = (double)N, iR;
   
   V = dvetor(N+2); /* Alocar memoria */
   
   for(int i = 0; i < Q2; ++i){
      phi[i] = 2.0 * M_PI * ran1(&semente);
   }
   for(int i = 0; i < N; ++i){
      V[i] = 0;
      iR = (double)i;
      for(int j = 0; j < Q2; ++j){
         k = (double)j + 1.0;
         V[i] += pow(k, -0.5 * alpha) * cos(2.0 * M_PI * iR * k / Nr + phi[j]);
      }
      aux1 += V[i];
      aux2 += pow(V[i], 2.0);
   }
   aux1 /= Nr;
   aux2 /= Nr;
   var = sqrt(aux2 + pow(aux1, 2.0));
   for(int i = 0; i < N; ++i){
      V[i] = (V[i] - aux1) / var;
   }
   return V;
}
double *SemNome2(int N, double alpha, long int semente){
   double E[N], vmedio, vmedio2, aux;
   double *V;
   int i, n;
   
   V = dvetor(N+2); /* Alocar memoria */
   
   vmedio = vmedio2 = 0.0;
   for(i = 1; i <= N; ++i){
      E[i] = 2.0 * ran1(&semente) - 1.0;
   }
   for(n = 1; n <= N; ++n){
      V[n] = 0.0;
      for(i = 1; i <= N; ++i){
         V[n] += E[i] / pow((double)abs(i - n) / alpha + 1.0, 2.0);
      }
      vmedio += V[n];
      vmedio2 += V[n] * V[n];
   }
   vmedio /= N;
   vmedio2 /= N;
   aux = 1.0 / sqrt(vmedio2 - vmedio * vmedio);
   for(n = 1; n <= N; ++n){
      V[n] = (V[n] - vmedio) * aux;
   }
   return V;
}

double *DefMassas(int criterio, double alpha, long int semente,
const int N){
   double *M;

   /* **********************************
   Usando o mapa de Bernoulli para escolher o valor das massas
   *********************************** */
   if(criterio == 1){
      double M0 = 0.5;
      M = MapaBernoulli(alpha, semente, N);
      for(int n = 1; n <= N; ++n) M[n] += M0;
   }

   /* **********************************
   Serie de correlacao usada no meu 1o artigo
   *********************************** */
   if(criterio == 2){
      double R[3] = {1.0, 0.1, 1.2}, MassaTipo[4] = {0.5, 1.0, 1.5, 2.0};
      M = SemNome1(N, alpha, semente);
      for(int n = 1; n <= N; ++n){
         if(M[n] < -R[0]){
            M[n] = MassaTipo[0];
         }else{
            if(M[n] > -R[0] && M[n] < R[1]){
               M[n] = MassaTipo[1];
            }else{
               if(M[n] > R[1] && M[n] < R[2]){
                  M[n] = MassaTipo[2];
               }else{
                  M[n] = MassaTipo[3];
               }
            }
         }
      }
   }

   /* **********************************
   Serie de correlacao que estou usando com o Carlos
   *********************************** */
   if(criterio == 3){
      /* Cadeia com tres tipos de massas */
      double b = 0.5, MassaTipo[3] = {0.5, 1.0, 1.5};
      M = SemNome2(N, alpha, semente);
      for(int n = 1; n <= N; ++n){
         if(M[n] < -b){
            M[n] = MassaTipo[0];
         }else{
            if(M[n] > b){
               M[n] = MassaTipo[2];
            }else{
               M[n] = MassaTipo[1];
            }
         }
      }
   }
   return M;
}
