/* *****************************************************************************
   Resolver a EDO usando Runge-Kutta de 8a ordem e 11 estagios
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
#define dt 2.5e-2
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
int rk8(int N, double *M, double *x, double *P, double **eta, double H0){

   /* Contadores */
   int i, j, n, contador = 0;
   const int contadorMAX = (int)(2.0/dt);
   /* Medidas de localizacao */
   double xMedio, sigma, *f;
   /* Runge-Kutta */
   const int s = 11; /* Quantidade de estagios */
   double **kP, **kX, **coef;
   double a[s][s], b[s], c[s];
   double t;
   const double tf = 0.4 * (double)N;
   /* Tamanho 'dinamico' da cadeia */
   int n0, nf;
   /* Energia e hamiltoniano */
   double H, *E;
   /* *** */
   double aux[2];
   const double c6 = 1.0/6.0;


   /* **********
      Ler a matriz de Runge-Kutta
   ********** */
   #include"rk8.h"
   for(i = 0; i < s; ++i){
      b[i] *= dt;
      for(j = 0; j < s; ++j) a[i][j] *= dt;
   }


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
   kX = dmatriz(i, s+1);
   kP = dmatriz(i, s+1);
   coef = dmatriz(i, 2);

   /* ***
   Zerar coeficientes que poderiam atrapalhar o calculo
   *** */
   for(n = 0; n <= N; ++n) coef[n][0] = coef[n][1] = 0.0;


   for(t = dt; t <= tf; t += dt){
      /* ***********************************************************************
         Rotina do metodo de Runge-kutta
      *********************************************************************** */
      for(i = 0; i < s; ++i){
      
         for(n = n0; n <= nf; ++n) coef[n][0] = coef[n][1] = 0.0;

         for(n = n0; n <= nf; ++n){
            for(j = 0; j < i; ++j){
               coef[n][0] += a[i][j] * kX[n][j];
               coef[n][1] += a[i][j] * kP[n][j];
            }
         }
         for(n = n0; n <= nf; ++n){
            aux[0] = (x[n+1] + coef[n+1][0]) - (x[n] + coef[n][0]);
            aux[1] = (x[n] + coef[n][0]) - (x[n-1] + coef[n-1][0]);
            kP[n][i] =
               eta[n][0] * aux[0] +
               eta[n][1] * pow2(aux[0]) +
               eta[n][2] * pow3(aux[0]) -(
               eta[n-1][0] * aux[1] +
               eta[n-1][1] * pow2(aux[1]) +
               eta[n-1][2] * pow3(aux[1]));
            kX[n][i] = 
               (P[n] + coef[n][1]) / M[n];
         }
      }
      for(n = n0; n <= nf; ++n){
         for(i = 0; i < s; ++i){
            P[n] += b[i] * kP[n][i];
            x[n] += b[i] * kX[n][i];
         }
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
