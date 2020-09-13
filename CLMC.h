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

/* ***
   Funcoes do Numerical Recipes
   ran1.c
   tqli.c
   CLMC_memoria.c
*** */
double ran1(long int*);
int tqli(double*, double*, int, double**);
double pythag(double, double);
double **dmatriz(const int, const int);
double *dvetor(const int);

/* ***
   CLMC_correlacoes.c
*** */
double *DefMassas(int,double,long int,int);

/* ***
   CLMC_acoplamento.c
*** */
double **acoplamento(double,double,double,int);

/* ***
   CLMC_PVI.c
*** */
double *DefPosicaoInicial(const int);
double *DefMomentoInicial(int,double*,double);
double CalcEnergiaInicial(int,double*,double*,double*,double**);

/* ***
   CLMC_arquivos.c
*** */
int AbrirArquivos(const char *format, ...);
void EscreverArquivos(double,int,double*,double);
void FecharArquivos(void);

#define RK4 200
int rk4(int,double*,double*,double*,double**,double);

#define RK8 201
int rk8(int,double*,double*,double*,double**,double);

#define RK14 202
int rk14(int,double*,double*,double*,double**,double);

#define CLMC_SUCESSO 0
#define CLMC_ERRO_ARQUIVO 3
#define CLMC_ERRO_NAN 4

#endif // CLMC_H
