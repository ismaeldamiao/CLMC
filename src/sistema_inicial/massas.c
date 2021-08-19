/* *****************************************************************************
   Esta funcao gera massas correlacionadas para cada particula da cadeia,
   eh possivel usar pelo menos tres tipos de estrategias para gerar correlacao:
   * transformada_de_Fourrier: Usa uma transformada de Fourrier para...
   * mapa_de_Bernoulli:
   * nao_sei_como_chamar_kkk:
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

void __massas(int semente){
   double *V = malloc(sizeof(double*));
   int n, N;

   N = config.quantidade_de_particulas;

   if(config.correlacao_das_massas == transformada_de_Fourrier){
      const double R[3] = {1.0, 0.1, 1.2}, MassaTipo[4] = {0.5, 1.0, 1.5, 2.0};
      V = correlated_w_fourier(config.fator_de_correlacao, N, semente);
      for(n = 0; n < N; ++n){
         if(V[n] < -R[0]){
            V[n] = MassaTipo[0];
         }else if((V[n] > -R[0]) && (V[n] < R[1])){
            V[n] = MassaTipo[1];
         }else if((V[n] > R[1]) && (V[n] < R[2])){
            V[n] = MassaTipo[2];
         }else{
            V[n] = MassaTipo[3];  
         }
      }
   }else if(config.correlacao_das_massas == mapa_de_Bernoulli){
      const double M0 = 0.5;
      V =  correlated_w_bernoulli(config.fator_de_correlacao, N, semente);
      for(n = 0; n < N; ++n)
      V[n] += M0;
   }else if(config.correlacao_das_massas == nao_sei_como_chamar_kkk){
      const double MassaTipo[3] = {0.5, 1.0, 1.5}, r = 0.1;
      if(config.fator_de_correlacao == 0.0) return;
      V = correlated_w_bernoulli(config.fator_de_correlacao, N, semente);
      for(n = 0; n < N; ++n){
         if(V[n] < -r){
            V[n] = MassaTipo[0];
         }else if((V[n] > -r) && (V[n] < r)){
            V[n] = MassaTipo[1];
         }else{
            V[n] = MassaTipo[2];
         }
      }
   }else{
      return;
   }

   for(n = 1; n <= N; ++n)
      cadeia[n].massa = V[n-1];
   cadeia[0].massa = cadeia[N+1].massa = 0.0;
   free(V);

   return;
}
