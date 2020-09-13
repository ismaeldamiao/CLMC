/* *****************************************************************************
   Funcoes para alocar memoria para as matrizes
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
#include<stdlib.h>

double **dmatriz(const int lin, const int col){
   int i;
   double **tmp;

   /* *****
   Alocar memoria para as linhas da matriz
   tmp recebe um array do tipo double* de tamanho lin
   ***** */
   tmp = (double**)malloc((size_t)(lin * sizeof(double*)));

   /* *****
   Alocar memoria para as colunas da matriz
   a i-esima posicao de tmp recebe um array do tipo double de tamanho col
   ***** */
   for(i = 0; i < lin; ++i)
      tmp[i] = (double*)malloc((size_t)(col * sizeof(double)));

   return tmp;
}

double *dvetor(const int tam){
   double *tmp;

   /* *****
   Alocar memoria para o vetor
   tmp recebe um array do tipo double de tamanho tam
   ***** */
   tmp = (double*)malloc((size_t)(tam * sizeof(double)));

   return tmp;
}
