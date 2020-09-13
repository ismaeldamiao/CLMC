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
#include"CLMC.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<stdarg.h> /* va_list va_start() */

/* *****************************************************************************
   Variaveis globais
***************************************************************************** */
FILE *fX, *fE, *fsig, *fLog;
FILE *gX, *gE, *gsig;
char nX[200], nE[200], nsig[200], NomeArquivo[200];
clock_t tE0, tEf;

int AbrirArquivos(const char *format, ...){
   va_list arg;
   char tmp[500];

   va_start(arg, format);
   vsprintf(tmp, format, arg);

   tE0 = clock();

   sprintf(nE, "Energia_%s.dat", tmp);
   fE = fopen(nE, "w");
   if((void*)fE == NULL) return CLMC_ERRO_ARQUIVO;

   sprintf(nsig, "Sigma_%s.dat", tmp);
   fsig = fopen(nsig, "w");
   if((void*)fsig == NULL) return CLMC_ERRO_ARQUIVO;

   sprintf(NomeArquivo, "Log_%s.log", tmp);
   fLog = fopen(NomeArquivo, "w");
   if((void*)fLog == NULL) return CLMC_ERRO_ARQUIVO;

   return CLMC_SUCESSO;
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

/*void FecharArquivos(double gamma){

   fclose(fX);
   gprintf(gX, "set key off");
   gprintf(gX, "set palette rgbformulae 23,28,3");
   gprintf(gX, "set ticslevel 0");
   gprintf(gX, "set surface");
   gprintf(gX, "set hidden3d");
   gprintf(gX, "set grid");
   gprintf(gX, "set xlabel \"t\"");
   gprintf(gX, "set ylabel \"Massa N\"");
   gprintf(gX, "set zlabel \"Energia\" offset -1,1,0 rotate by 90");
   gprintf(gX, "set view 50.0, 75.0");
   gprintf(gX, "splot \"%s\" with pm3d palette", nX);
   gclose(gX);

   fclose(fE);
   gprintf(gE, "set key off");
   gprintf(gE, "set palette rgbformulae 23,28,3");
   gprintf(gE, "set ticslevel 0");
   gprintf(gE, "set surface");
   gprintf(gE, "set hidden3d");
   gprintf(gE, "set grid");
   gprintf(gE, "set xlabel \"t\"");
   gprintf(gE, "set ylabel \"Massa N\"");
   gprintf(gE, "set zlabel \"Energia\" offset -1,1,0 rotate by 90");
   gprintf(gE, "set view 50.0, 75.0");
   gprintf(gE, "splot \"%s\" with pm3d palette", nE);
   gclose(gE);

   fclose(fsig);
   gprintf(gsig, "set key box left top");
   //gprintf(gsig, "set key off");
   gprintf(gsig, "set xlabel \"t\"");
   gprintf(gsig, "set ylabel \"σ(t)\"");
   gprintf(gsig, "plot \"%s\" w lp lt rgb \"blue\" pt 6 ps 1.2\\", nsig);
   gprintf(gsig, "title \"γ = %g\"", gamma);
   gclose(gsig);

   system(
      "if command -v gnuplot > /dev/null; then env gnuplot *.gnuplot; fi"
   );

   tEf = clock();
   fprintf(fLog, "Tempo de execucao: %g s",
      (double)(tEf - tE0) / CLOCKS_PER_SEC);
   fclose(fLog);

   return;
}*/
