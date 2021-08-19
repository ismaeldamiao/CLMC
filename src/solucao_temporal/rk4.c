/* *****************************************************************************
   Esta funcao resolve as equacoes de Hamilton para o sistema
   usando o metodo classico de Runge-Kutta de 4a ordem.
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
//#include <stdio.h>

/* Algumas instrucoes para o pre-processador, caso o compilador seja o GCC */
#if defined(__GNUC__)
#pragma GCC optimize ("O3")
#pragma GCC diagnostic warning "-Wall"
#pragma GCC diagnostic warning "-Wextra"
#pragma GCC diagnostic warning "-Wpedantic"
#endif /* __GNUC__ */

/* Uma caracteristica dos metodos de Runge-Kutta eh que quanto menor for
   o dt, a discretizacao no tempo, entao melhor serah a aproximacao numerica
   para a solucao das equacoes de Hamilton. Essa melhora na aproximacao eh
   conseguida a custo de um maior tempo de calculo. */
#define dt 1.0e-3
#define dt2 0.5e-3 /* Deve ser metade do dt */

/* Notacao mais simples para as variaveis. Aqui a notacao descritiva eh deveras
   inconveniente pois deixa os argumentos desnecessariamente longos e feios. */
#define x(n) cadeia[n].posicao
#define P(n) cadeia[n].momento

/* *****************************************************************************
   Funcao do calculo numerico
***************************************************************************** */
int rk4(void){

   int status, i, j, n, contador, contadorMAX, n0, nf;
   double t, tf, dt6, *kP1, *kX1, *kP2, *kX2, *kP3, *kX3, *kP4, *kX4;

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
   kX1 = vetor(i);
   kP1 = vetor(i);
   kX2 = vetor(i);
   kP2 = vetor(i);
   kX3 = vetor(i);
   kP3 = vetor(i);
   kX4 = vetor(i);
   kP4 = vetor(i);

   /* Zerar coeficientes que poderiam atrapalhar o calculo */
   for(n = 0; n <= N+1; ++n)
   kX1[n] = kX2[n] = kX3[n] = kX4[n] = 0.0;

   /* Alem de dt e dt2, o rk4 precisa tambem do dt6, definido como um sexto
      do dt, usar essas variaveis eh conveniente para evitar calculos repetidos,
      ja que a rotina sera executada entre centenas e trihoes de vezes, ou
      mais. */
   dt6 = dt/6.0;

   /* O tempo final de calculo depende das suas necessidades na pesquisa */
   tf = 0.4 * (double)N;

   for(t = dt; t <= tf; t += dt){
      /* ***********************************************************************
         Rotina do metodo de Runge-kutta de 4a ordem
      *********************************************************************** */
      /* K1 */
      for(n = n0; n <= nf; ++n){
         kP1[n] = __forca(n, x(n-1), x(n), x(n+1));
         kX1[n] = __velocidade(n, P(n));
      }
      /* K2 */
      for(n = n0; n <= nf; ++n){
         i = n + 1;
         j = n - 1;
         kP2[n] = __forca(n,
         x(j) + dt2 * kX1[j], x(n) + dt2 * kX1[n], x(i) + dt2 * kX1[i]);
         kX2[n] = __velocidade(n, P(n) + dt2 * kP1[n]);
      }
      /* K3 */
      for(n = n0; n <= nf; ++n){
         i = n + 1;
         j = n - 1;
         kP3[n] = __forca(n,
         x(j) + dt2 * kX2[j], x(n) + dt2 * kX2[n], x(i) + dt2 * kX2[i]);
         kX3[n] = __velocidade(n, P(n) + dt2 * kP2[n]);
      }
      /* K4 */
      for(n = n0; n <= nf; ++n){
         i = n + 1;
         j = n - 1;
         kP4[n] = __forca(n,
         x(j) + dt * kX3[j], x(n) + dt * kX3[n], x(i) + dt * kX3[i]);
         kX4[n] = __velocidade(n, P(n) + dt * kP3[n]);
      }

      for(n = n0; n <= nf; ++n){
         x(n) += (kX1[n] + 2.0*kX2[n] + 2.0*kX3[n] + kX4[n]) * dt6;
         P(n) += (kP1[n] + 2.0*kP2[n] + 2.0*kP3[n] + kP4[n]) * dt6;
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
