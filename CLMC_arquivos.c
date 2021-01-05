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
   Variaveis globais
***************************************************************************** */
FILE *arquivo_energia, *arquivo_dispersao;

void AbrirArquivos(long int semente, const char *format, ...){
   va_list arg;
   char tmp[512], cmd[512];
   char NomeArquivo[512];

   va_start(arg, format);
   vsprintf(tmp, format, arg);

   sprintf(cmd, "! [ -d %s ] && env mkdir %s", tmp, tmp);
   system(cmd);

   sprintf(cmd,
   "! [ -e info.txt ] && echo \"N=%d\nALPHA=%g\nV0=%g\nETA2=%g\nETA3=%g\nETA4=%g\" >> %s/info.txt",
   N, __ALPHA__, __V0__, __ETA2__, __ETA3__, __ETA4__, tmp);
   system(cmd);

   sprintf(NomeArquivo, "%s/energia_%ld.dat", tmp, semente);
   arquivo_energia = fopen(NomeArquivo, "w");

   sprintf(NomeArquivo, "%s/dispersao_%ld.dat", tmp, semente);
   arquivo_dispersao = fopen(NomeArquivo, "w");

   return;
}

void EscreverArquivos(double t){
   int n;
   for(n = 1; n <= N; ++n)
   fprintf(arquivo_energia, "%g %d %g\n", t, n, E[n]);
   fprintf(arquivo_energia, "\n");
   fprintf(arquivo_dispersao, "%5.4g %g %g\n", t, sigma, Z);
   return;
}

void FecharArquivos(void){
   fclose(arquivo_energia);
   fclose(arquivo_dispersao);
   return;
}
