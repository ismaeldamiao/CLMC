/* *****************************************************************************
   Rotinas para abrir, escrever e fechar arquivos.
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
#include<stdio.h>
/* *****************************************************************************
   Variaveis globais
***************************************************************************** */
FILE *fsig, *fE;

void AbrirArquivos(const char *format, ...){
   va_list arg;
   char tmp[500];
   char NomeArquivo[200];

   va_start(arg, format);
   vsprintf(tmp, format, arg);

   sprintf(NomeArquivo, "Energia_%s.dat", tmp);
   fE = fopen(NomeArquivo, "w");

   sprintf(NomeArquivo, "Sigma_%s.dat", tmp);
   fsig = fopen(NomeArquivo, "w");

   return;
}

void EscreverArquivos(double t, int N, double *E, double sigma){
   int n;
   for(n = 1; n <= N; ++n){
      fprintf(fE, "%g %d %g\n", t, n, E[n]);
      //fprintf(fX, "%g %d %g\n", t, iM, x[iM]);
   }
   fprintf(fE, "\n");
   //fprintf(fX, "\n");
   fprintf(fsig, "%5.4g %g\n", t, sigma);

   return;
}

void FecharArquivos(void){

   fclose(fE);
   fclose(fsig);

   return;
}
