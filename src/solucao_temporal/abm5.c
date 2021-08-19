/* *****************************************************************************
   Resolver a EDO usando Adams-Bashforth-Moulton de 10 ordem
   iniciado com Runge-Kutta de 14a ordem e 35 estagios
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
#include "../CLMC.h"
#include "../libdamiao/damiao.h"
//#include <stdio.h>

/* Esta eh a discretizacao */
#define dt 5.0e-3

/* Notacao mais simples para as variaveis. Aqui a notacao descritiva eh deveras
   inconveniente pois deixa os argumentos desnecessariamente longos e feios. */
#define x(n) cadeia[n].posicao
#define P(n) cadeia[n].momento

#define s 35
#define s_abm 5
/* *****************************************************************************
   Funcao do calculo numerico
***************************************************************************** */
int abm5(void){

   int status, i, j, l, n, contador, contadorMAX, n0, nf;
   double t, tf, **kX, **kP, **coef, **kP_abm, **kX_abm, *x_ab, *P_ab;
   long double a[s][s], b[s];
   long double adams[2][s_abm] = {
      // Coeficientes de Adams–Bashforth
      {251.0/720.0, -1274.0/720.0, 2616.0/720.0, -2774.0/720.0, 1901.0/720.0},
      // Coeficientes de Adams–Moulton
      {-19.0/720.0, 106.0/720.0, -264.0/720.0, 646.0/720.0, 251.0/720.0}
   };

   contador = 0;
   contadorMAX = (int)(2.0/dt);


   /* Para economizar tempo e resursos, a cadeia possui limites de calculo,
      as equacoes nao serao resolvidas para particulas fora desse limite. Caso
      o pacote de energia se aproxime desse limite entao ele serah eventualmente
      redefinido para outro limite onde o pacote de energia ainda nao chegou. */
   if(N > 600){
      const int N2 = N/2;
      n0 = N2 - 300;
      nf = N2 + 300;
   }else{
      n0 = 1;
      nf = N;
   }

   /* Alocar memoria para os vetores */
   i = N+2;
   kX = matriz(i, s+1);
   kP = matriz(i, s+1);
   kX_abm = matriz(i, s_abm);
   kP_abm = matriz(i, s_abm);
   coef = matriz(i, 2);
   x_ab = vetor(i);
   P_ab = vetor(i);


   for(i = 0; i < s; ++i){
      b[i] = __rk14_b[i] * dt;
      for(j = 0; j < s; ++j) a[i][j] = __rk14_a[i][j] * dt;
   }

   for(n = 0; n <= N; ++n) coef[n][0] = coef[n][1] = 0.0;

   /* O tempo final de calculo depende das suas necessidades na pesquisa */
   tf = 0.4 * (double)N;


   for(l = 0; l < s_abm; ++l){
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
            kP[n][i] = __forca(n,
            x(n-1) + coef[n-1][0], x(n) + coef[n][0], x(n+1) + coef[n+1][0]);
            kX[n][i] = __velocidade(n, P(n) + coef[n][1]);
         }
      }
      for(n = n0; n <= nf; ++n){
         for(i = 0; i < s; ++i){
            P(n) += b[i] * kP[n][i];
            x(n) += b[i] * kX[n][i];
         }
      }
      for(n = n0; n <= nf; ++n){
         kP_abm[n][l] = __forca(n, x(n-1), x(n), x(n+1));
         kX_abm[n][l] = __velocidade(n, P(n));
      }
   }
   j = s_abm-1;
   for(t = (double)(s_abm)*dt; t <= tf; t += dt){
      /* ***********************************************************************
         Rotina do metodo de Adams-Bashforth-Moulton
      *********************************************************************** */
      for(n = n0; n <= nf; ++n){ //AB
         P_ab[n] = 0.0;
         x_ab[n] = 0.0;
         for(i = 0; i < s_abm; ++i){
            P_ab[n] += adams[0][i] * kP_abm[n][i];
            x_ab[n] += adams[0][i] * kX_abm[n][i];
         }
         P_ab[n] = P(n) + P_ab[n] * dt;
         x_ab[n] = x(n) + x_ab[n] * dt;
      }
      for(n = n0; n <= nf; ++n){
         for(i = 0; i < j; ++i){
            kP_abm[n][i] = kP_abm[n][i+1];
            kX_abm[n][i] = kX_abm[n][i+1];
         }
         kP_abm[n][j] = __forca(n, x_ab[n-1], x_ab[n], x_ab[n+1]);
         kX_abm[n][j] = __velocidade(n, P_ab[n]);
      }
      for(n = n0; n <= nf; ++n){ //AM
         P_ab[n] = 0.0;
         x_ab[n] = 0.0;
         for(i = 0; i < s_abm; ++i){
            P_ab[n] += adams[1][i] * kP_abm[n][i];
            x_ab[n] += adams[1][i] * kX_abm[n][i];
         }
         P(n) += P_ab[n] * dt;
         x(n) += x_ab[n] * dt;
      }

      if(++contador >= contadorMAX){
         contador = 0;

         /* Uma vez que se encontrou a solucao numerica das equacoes em um
            dado tempo t, podemos fazer algumas analises desses dados. Essas
            analises serao feitas pela funcao __escrever_arquivos. */
         status = __escrever_arquivos(t);
         if(status != CLMC_SUCESSO) return status;
         
         //printf("OK\n");
         if(N > 600){
            if(
               (cadeia[n0].energia > 1.0e-20) || (cadeia[nf].energia > 1.0e-20)
            ){
               /* Esta rotina cresce o tamanho dos limites de calculo da cadeia
                  conforme o pacote de energia se aproxima desses limites. */
               if((nf+30) < N){
                  n0 -= 30;
                  nf += 30;
               }else{
                  n0 = 1;
                  nf = N;
               }
            }
         }
      }
   }
   return CLMC_SUCESSO;
}
