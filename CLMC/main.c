/* *****************************************************************************
   Programa para alguma coisa.
   if ! [ -x COMPILE ]; then chmod 755 COMPILE; fi
   ./COMPILE x; echo $?
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
#include"CLMC_configuracao.h"
/* *****************************************************************************
   Declaracoes globais
***************************************************************************** */   
/* *****************************************************************************
   Funcao principal
***************************************************************************** */
int main(int argc, char *argv[]){

   double **eta, *M, *x, *P, H0;
   int ESTADO;
   long int semente[2] = {-Semente, Semente};

   /* **********
   As massas serao mapeadas/escolhidas com base em series.
   Ver arquivo CLMC_correlacoes.c
   ********** */
   M = DefMassas(criterio, alpha, semente[0], N);

   /* **********
   Os termos de acoplamento, por padrao, sao todos iguais para cada potencial
   eta2 eh o termo de acoplamento do potencial harmonico enquanto que
   eta3 e eta4 sao dos potenciais anarmonicos.
   Ver arquivo CLMC_PVI.c
   ********** */
   eta = acoplamento(eta2, eta3, eta4, N);

   /* **********
   As posicoes e os momentos iniciais sao todos nulos, exceto o momento
   na central da cadeia, ele eh o produto de M[n] com a velocidade v0
   Ver arquivo CLMC_PVI.c
   ********** */
   x = DefPosicaoInicial(N);
   P = DefMomentoInicial(N, M, v0);

   /* **********
   Calcula o valor do hamiltoniano no instante inicial (t=0)
   Ver arquivo CLMC_PVI.c
   ********** */
   H0 = CalcEnergiaInicial(N, M, x, P, eta);

   /* **********
   Ver arquivo CLMC_arquivos.c
   ********** */
   AbrirArquivos("%iMassas_%galpha_%gV0_%geta2_%geta3_%geta4_%02ldsemente",
      N, alpha, v0, eta2, eta3, eta4, semente[1]);

   if(metodo == RK4){
      ESTADO = rk4(N, M, x, P, eta, H0);
   }else if(metodo == RK8){
      ESTADO = rk8(N, M, x, P, eta, H0);
   }else if(metodo == RK14){
      ESTADO = rk14(N, M, x, P, eta, H0);
   }


   if(ESTADO != CLMC_SUCESSO) return ESTADO;
   return CLMC_SUCESSO;
}
