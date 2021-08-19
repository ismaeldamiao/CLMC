/* *****************************************************************************
   Esta funcao resolve as equacoes de Hamilton para o sistema
   usando o metodo de Runge-Kutta de 8a ordem e 11 estagios.
   *****************************************************************************
   E-mail: ismaellxd@gmail.com
   Site: https://ismaeldamiao.github.io/
   *****************************************************************************
   Copyright (c) 2020 I.F.F. dos SANTOS (Ismael Damiao)

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

/* Uma caracteristica dos metodos de Runge-Kutta eh que quanto menor for
   o dt, a discretizacao no tempo, entao melhor serah a aproximacao numerica
   para a solucao das equacoes de Hamilton. Essa melhora na aproximacao eh
   conseguida a custo de um maior tempo de calculo. */
#define dt 1.0e-2

/* Notacao mais simples para as variaveis. Aqui a notacao descritiva eh deveras
   inconveniente pois deixa os argumentos desnecessariamente longos e feios. */
#define x(n) cadeia[n].posicao
#define P(n) cadeia[n].momento

/* *****************************************************************************
   Funcao do calculo numerico
***************************************************************************** */
int rk8(void){

   int status, i, j, n, contador, contadorMAX, n0, nf;
   double t, tf, **kX, **kP, **coef;
   const int s = 11;
   long double a[s][s], b[s];

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
   coef = matriz(i, 2);


   for(i = 0; i < s; ++i){
      b[i] = __rk8_b[i] * dt;
      for(j = 0; j < s; ++j) a[i][j] = __rk8_a[i][j] * dt;
   }

   for(n = 0; n <= N; ++n) coef[n][0] = coef[n][1] = 0.0;

   /* O tempo final de calculo depende das suas necessidades na pesquisa */
   tf = 0.4 * (double)N;

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
