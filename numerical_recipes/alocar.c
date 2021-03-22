/* *****************************************************************************
   Macros para alocar memoria para os vetores e matrizes
   *****************************************************************************
   E-mail: ismaellxd@gmail.com
   Site: https://ismaeldamiao.github.io/
***************************************************************************** */
#include<stdlib.h>

/* ***
   void *malloc(size_t)
   Aloca memoria para uma array de tamanho size_t.
   * A variavel array deve ser um poteiro.
   * A variavel type deve ser o tipo de variavel (int, double...).
   * A variavel tam deve ser a dimensao do vetor.
*** */
#define vetor(tam, type, array) {\
array = (type*)malloc((size_t)((tam) * sizeof(type)));}

/* ***
   Para fazer uma matriz a ideia eh primeiro fazer um vetor de ponteiros e,
   para cada componente do vetor, fazer um vetor normal.
   * A variavel array deve ser um poteiro de ponteiro.
   * A variavel type deve ser o tipo de variavel (int, double...).
   * A variavel lin deve ser a quantidade de linhas da matriz.
   * A variavel col deve ser a quantidade de colunas da matriz.
*** */
#define matriz(lin, col, type, array) {\
vetor(lin, type*, array);\
for(int __i_index__ = 0; __i_index__ < (lin); ++__i_index__)\
vetor(col, type, array[__i_index__]);}
