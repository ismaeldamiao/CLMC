/* *****************************************************************************
   Biblioteca de funcoes CLMC
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
#ifndef CLMC_H
#define CLMC_H 0

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RK4 200
#define RK8 201
#define RK14 202
#define ABM5 203
#define ABM10 204
#include "CLMC_configuracao.h"

#define CLMC_SUCESSO 0
#define CLMC_ERRO_ARQUIVO 3
#define CLMC_ERRO_NAN 4
#define CLMC_ERRO_PRECISAO 5
/* *****************************************************************************
   Formulas, constantes e funcoes matematicas
***************************************************************************** */
#define DeltaDeKronecker(x, y) ((x) != (y) ? 0.0 : 1.0)
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#ifndef abs
   #define abs(x) ((x) >= 0 ? (x) : -(x))
#endif
#ifndef M_PI
   #define M_PI 3.14159265358979323846
#endif

/* *****************************************************************************
   Funcoes do Numerical Recipes
***************************************************************************** */
#include "numerical_recipes/ran1.c"
#include "numerical_recipes/tqli.c"
#include "numerical_recipes/alocar.c"
/* *****************************************************************************
   Funcoes que preparam o sistema fisico
***************************************************************************** */
#include "sistema_inicial/massas.c"
#include "sistema_inicial/acoplamentos.c"
#include "sistema_inicial/posicoes.c"
#include "sistema_inicial/momentos.c"
/* *****************************************************************************
   Funcoes relacionadas com o hamiltoniano ou com suas equacoes
***************************************************************************** */
#include "hamiltoniano/energia.c"
#include "hamiltoniano/forca.c"
#include "hamiltoniano/velocidade.c"
/* *****************************************************************************
   Funcoes para administrar os arquivos de dados
***************************************************************************** */
#include "CLMC_arquivos.c"
/* *****************************************************************************
   Funcoes para resolver a evolucao temporal do sistema fisico
***************************************************************************** */
#if __METODO__ == RK4
   #include "solucao_temporal/rk4.c"
#elif __METODO__ == RK8
   #include "solucao_temporal/rk8.c"
#elif __METODO__ == RK14
   #include "solucao_temporal/rk14.c"
#elif __METODO__ == ABM5
   #include "solucao_temporal/abm5.c"
#elif __METODO__ == ABM10
   #include "solucao_temporal/abm10.c"
#endif

#endif // CLMC_H
