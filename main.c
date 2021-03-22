/* *****************************************************************************
   Programa para calcular a dinamica dos modos de vibracao em uma cadeia
   linear (unidimencional) e calcular medidas de localizacao da energia na
   cadeia.
   
   As equacoes de hamilton podem ser resolvidas usando os metodos de
   Runge-Kutta de 4a, 8a e 14a ordem.

   $ clang main.c -std=c99 -O2 -lm -o clmc && ./clmc
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
double *M, **eta, *x, *P, *E;
double H0, sigma, Z;
/* *****************************************************************************
   Bibliotecas
***************************************************************************** */
#include "CLMC.h"
/* *****************************************************************************
   Declaracoes globais
***************************************************************************** */   
/* *****************************************************************************
   Funcao principal
***************************************************************************** */
int main(int argc, char *argv[]){

   int ESTADO, n;
   long int semente;

   if(argc < 2) semente = 1;
   else semente = atol(argv[1]);
   if(semente < 1) semente = -semente;

   /* **************************************************************************
      Primeiro preparo o sistema fisico
   ************************************************************************** */
   M = __massas__(semente); /* Massas das particulas/atomos */
   eta = __acoplamentos__(); /* Constantes de acoplamento */
   x = __posicoes__(); /* Posicoes iniciais */
   P = __momentos__(); /* Momentos iniciais */
   /* **************************************************************************
      Preparar arquivos de dados a serem escritos
   ************************************************************************** */
   AbrirArquivos(semente,
   "%iMassas_%galpha_%gV0_%geta2_%geta3_%geta4",
   N, __ALPHA__, __V0__, __ETA2__, __ETA3__, __ETA4__);
   /* **************************************************************************
      Resolver a evolucao temporal do sistema fisico
   ************************************************************************** */
   H0 = 0.0;/* Hamiltoniano no tempo t = 0 */
   for(n = 1; n <= N; ++n) H0 += energia(n);

   #if __METODO__ == RK4
      ESTADO = rk4();
   #elif __METODO__ == RK8
      ESTADO = rk8();
   #elif __METODO__ == RK14
      ESTADO = rk14();
   #elif __METODO__ == ABM5
      ESTADO = abm5();
   #elif __METODO__ == ABM10
      ESTADO = abm10();
   #endif

   FecharArquivos();

   return ESTADO;
}
