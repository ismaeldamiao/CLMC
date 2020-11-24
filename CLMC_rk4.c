/* *****************************************************************************
   Resolver a EDO usando Runge-Kutta classico de 4a ordem
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
/* *****************************************************************************
   Bibliotecas
***************************************************************************** */
#include"CLMC.h"
#include<math.h>
/* *****************************************************************************
   Definicoes
***************************************************************************** */
/* Dicretiacao do Runge-Kutta */
#define dt 1.0e-3
#define dt2 0.5e-3
#define pow2(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))
#define pow4(x) ((x)*(x)*(x)*(x))
#define En (0.5 * (P[n]*P[n]) / M[n] +\
      (eta[i][0] * pow2(aux[0]) + eta[n][0] * pow2(aux[1])) * 0.25 +\
      (eta[i][1] * pow3(aux[0]) + eta[n][1] * pow3(aux[1])) * c6 +\
      (eta[i][2] * pow4(aux[0]) + eta[n][2] * pow4(aux[1])) * 0.125)
/* *****************************************************************************
   Funcao do calculo numerico
***************************************************************************** */
int rk4(int N, double *M, double *x, double *P, double **eta, double H0){

   /* Contadores */
   int i, n, contador = 0;
   const int contadorMAX = (int)(2.0/dt);
   /* Medidas de localizacao */
   double xMedio, sigma, *f;
   /* Runge-Kutta */
   double *kP1, *kX1, *kP2, *kX2, *kP3, *kX3, *kP4, *kX4;
   double t;
   const double tf = 0.4 * (double)N, dt6 = dt/6.0;
   /* Tamanho 'dinamico' da cadeia */
   int n0, nf;
   /* Energia e hamiltoniano */
   double H, *E;
   /* *** */
   double aux[2];
   const double c6 = 1.0/6.0;


   /* **********
      A cadeia possui um temanho dinamico da seguinte forma:
      Se a energia nos sitios n0 e nf for pequena, quase nula, entao
      entende-se que nao ha grande interferencia na dinamica do sistema
      por parte desses sitios.
      Se essa energia for rezoavelmente diferente de zero entao entao
      o tamanho da cadeia irah crezcer mudando os valores de n0 e nf para
      valores onde novamente essa energia pode ser considerada nula.
   ********** */
   if(N > 600){
      n0 = N/2 - 300; nf = N/2 + 300;
   }else{
      n0 = 1; nf = N;
   }

   /* ***
   Alocar memoria para os vetores
   *** */
   i = N+2;
   E = dvetor(i);
   f = dvetor(i);
   kX1 = dvetor(i);
   kP1 = dvetor(i);
   kX2 = dvetor(i);
   kP2 = dvetor(i);
   kX3 = dvetor(i);
   kP3 = dvetor(i);
   kX4 = dvetor(i);
   kP4 = dvetor(i);

   /* ***
   Zerar coeficientes que poderiam atrapalhar o calculo
   *** */
   i = N+1;
   kX1[0] = kX2[0] = kX3[0] = kX4[0] = kX1[i] = kX2[i] = kX3[i] = kX4[i] = 0.0;


   for(t = dt; t <= tf; t += dt){
      /* ***********************************************************************
         Rotina do metodo de Runge-kutta de 4a ordem
      *********************************************************************** */
      /* K1 */
      for(n = n0; n <= nf; ++n){
         i = n + 1;
         aux[0] = x[i] - x[n];
         aux[1] = x[n] - x[n-1];
         kP1[n] =
            eta[i][0] * aux[0] +
            eta[i][1] * pow2(aux[0]) +
            eta[i][2] * pow3(aux[0]) -(
            eta[n][0] * aux[1] +
            eta[n][1] * pow2(aux[1]) +
            eta[n][2] * pow3(aux[1]));
         kX1[n] =
            P[n] / M[n];
      }
      /* K2 */
      for(n = n0; n <= nf; ++n){
         i = n + 1;
         aux[0] = x[i] + dt2 * kX1[i] - (x[n] + dt2 * kX1[n]);
         aux[1] = x[n] + dt2 * kX1[n] - (x[n-1] + dt2 * kX1[n-1]);
         kP2[n] =
            eta[i][0] * aux[0] +
            eta[i][1] * pow2(aux[0]) +
            eta[i][2] * pow3(aux[0]) -(
            eta[n][0] * aux[1] +
            eta[n][1] * pow2(aux[1]) +
            eta[n][2] * pow3(aux[1]));
         kX2[n] =
            (P[n] + dt2 * kP1[n]) / M[n];
      }
      /* K3 */
      for(n = n0; n <= nf; ++n){
         i = n + 1;
         aux[0] = x[i] + dt2 * kX2[i] - (x[n] + dt2 * kX2[n]);
         aux[1] = x[n] + dt2 * kX2[n] - (x[n-1] + dt2 * kX2[n-1]);
         kP3[n] =
            eta[i][0] * aux[0] +
            eta[i][1] * pow2(aux[0]) +
            eta[i][2] * pow3(aux[0]) -(
            eta[n][0] * aux[1] +
            eta[n][1] * pow2(aux[1]) +
            eta[n][2] * pow3(aux[1]));
         kX3[n] =
            (P[n] + dt2 * kP2[n]) / M[n];
      }
      /* K4 */
      for(n = n0; n <= nf; ++n){
         i = n + 1;
         aux[0] = x[i] + dt * kX3[i] - (x[n] + dt * kX3[n]);
         aux[1] = x[n] + dt * kX3[n] - (x[n-1] + dt * kX3[n-1]);
         kP4[n] =
            eta[i][0] * aux[0] +
            eta[i][1] * pow2(aux[0]) +
            eta[i][2] * pow3(aux[0]) -(
            eta[n][0] * aux[1] +
            eta[n][1] * pow2(aux[1]) +
            eta[n][2] * pow3(aux[1]));
         kX4[n] =
            (P[n] + dt * kP3[n]) / M[n];
      }
      
      for(n = n0; n <= nf; ++n){
         x[n] += (kX1[n] + 2.0*kX2[n] + 2.0*kX3[n] + kX4[n]) * dt6;
         P[n] += (kP1[n] + 2.0*kP2[n] + 2.0*kP3[n] + kP4[n]) * dt6;
      }

      /* ***********************************************************************
         Rotina para calcular as medidas de localizacao
      *********************************************************************** */
      if(++contador > contadorMAX){
         contador = 0;
         /* ***
         Calculo do hamiltoniano
         *** */
         H = 0.0;
         for(n = 1; n <= N; ++n){
            i = n + 1;
            aux[0] = x[i] - x[n];
            aux[1] = x[n] - x[n-1];
            E[n] = En; /* En expandido pela macro */
            H += E[n];
         }

         /* ***
         Verificar se nao houve problema de numero ilegal ou falta de precisao
         *** */
         if(H != H) return CLMC_ERRO_NAN;
         if(fabs(1.0 - H / H0) > 1.0e-8) return 138;

         /* ***
         Calcular fracao da energia total no sitio n para cada n
         *** */
         xMedio = 0.0;
         for(n = 1; n <= N; ++n){
            f[n] = E[n] / H0;
            xMedio += (double)n * f[n];
         }

         /* ***
         Calcular dispersao da energia na cadeia
         *** */
         sigma = 0.0;
         for(n = 1; n <= N; ++n){
            aux[0] = (double)n - xMedio;
            sigma += aux[0] * aux[0] * f[n];
         }
         sigma = sqrt(sigma);

         /* ***
         Crescer tamanho da cadeia caso energia nas bordas seja relevante
         *** */
         if((N > 600) && ((E[n0] > 1.0e-20) || (E[nf] > 1.0e-20))){
            n0 -= 30; nf += 30;
         }

         /* ***
         Ver arquivo CLMC_arquivos.c
         *** */
         EscreverArquivos(t, N, E, sigma);
      }
   }
   return CLMC_SUCESSO;
}
