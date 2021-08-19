/* *****************************************************************************
   Programa para calcular a dinamica dos modos de vibracao em uma cadeia
   linear (unidimencional) e calcular medidas de localizacao da energia na
   cadeia.
   
   As equacoes de hamilton podem ser resolvidas usando os metodos de
   Runge-Kutta de 4a, 8a e 14a ordem.

   $ ./COMPILE && ./clmc
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
#include "./libdamiao/damiao.h"
#include "./CLMC.h"

int main(int argc, char **argv){
   int semente, estado = CLMC_SUCESSO;

   if(argc < 2) semente = 1;
   else semente = atoi(argv[1]);

   config = ler_configuracao();
   /* Verificar possiveis erros de configuracao. */
   if(!(N > 0)){
      fprintf(stderr, "'quantidade_de_particulas' deve ser positiva.\n");
      return CLMC_ERRO_DE_CONFIGURACAO;
   }
   if(!(config.fator_de_correlacao > 0)){
      fprintf(stderr, "'fator_de_correlacao' deve ser positiva.\n");
      return CLMC_ERRO_DE_CONFIGURACAO;
   }

   cadeia = vetorp(N+2);

   /* Em qualquer programa de simulacao de um sistema fisico o primeiro a se
      fazer eh preparar o sistema, de modo similar a como um fisico experimental
      prepara seu experimento. Nesse caso preciso que no instante inicial
      todas as posicoes e todos os momentos tenham valores bem definidos.
      Em muitos sistemas temos, alem de propriedades como posicoes, que
      estao sob o nosso controle, dentro de certo limite, temos tambem
      outras grandezas que dependem do material que estamos trabalhando,
      elas aparecem nas equacoes como parametros que vao controlar os possiveis
      estados finais aos quais pode levar uma condicao inicial, nesse
      sistema esses parametros sao as massas e os termos de acoplamento. */
   __massas(semente);
   __posicoes();
   __momentos();
   __acoplamentos();

   /* */
   __abrir_arquivos(
      config,
      semente,
      "%d_massas-%g_alpha-%g_V0-%g_eta2-%g_eta3-%g_eta4",
      config.quantidade_de_particulas,
      config.fator_de_correlacao,
      config.termo_de_acoplamento_linear,
      config.termo_de_acoplamento_quadratico,
      config.termo_de_acoplamento_cubico
   );

   H0 = 0.0;
   for(int n = 1; n <= N; ++n)
      H0 += __energia(n);

   if(config.metodo_de_solucao == RK4) estado = rk4();
   else if(config.metodo_de_solucao == RK8) estado = rk8();
   else if(config.metodo_de_solucao == RK14) estado = rk14();
   else if(config.metodo_de_solucao == ABM5) estado = abm5();
   else if(config.metodo_de_solucao == ABM10) estado = abm10();

   __fechar_arquivos();

   return estado;
}
