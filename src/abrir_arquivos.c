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
#include "./libdamiao/damiao.h"
#include "./CLMC.h"

FILE *arquivo_energia, *arquivo_dispersao;

void __abrir_arquivos(
   configuracao config,
   int semente,
   const char *format, ...
){
   va_list arg;
   char tmp[512], cmd[2048];
   char NomeArquivo[1024];

   va_start(arg, format);
   vsprintf(tmp, format, arg);

   sprintf(cmd, "[ -d '%s' ] || env mkdir '%s'", tmp, tmp);
   system(cmd);

   sprintf(cmd,
   "[ -f '%s/info.txt' ] || "
   "echo 'N=%d\nALPHA=%g\nV0=%g\nETA2=%g\nETA3=%g\nETA4=%g' >> '%s/info.txt'",
   tmp, config.quantidade_de_particulas,
   config.fator_de_correlacao,
   config.velocidade_inicial,
   config.termo_de_acoplamento_linear,
   config.termo_de_acoplamento_quadratico,
   config.termo_de_acoplamento_cubico, tmp);
   system(cmd);

   sprintf(NomeArquivo, "%s/energia_%d.dat", tmp, semente);
   arquivo_energia = fopen(NomeArquivo, "w");

   sprintf(NomeArquivo, "%s/dispersao_%d.dat", tmp, semente);
   arquivo_dispersao = fopen(NomeArquivo, "w");

   return;
}

int __escrever_arquivos(double t){
   const int N2 = N/2;
   int n;
   double H, sigma, Z, aux;
   
   /* Calculo do hamiltoniano */
   H = 0.0;
   for(n = 1; n <= N; ++n){
      cadeia[n].energia = __energia(n);
      H += cadeia[n].energia;
   }
   /* Verificar se nao houve problema de numero ilegal ou falta de precisao */
   if(isnan(H)){
      fprintf(stderr,"Algum calculo resultou em nan (=nao eh numero).\n");
      return CLMC_ERRO_NAN;
   }
   if(fabs(1.0 - H / H0) > 1.0e-8){
      fprintf(stderr, "Falta de precisao nos calculos.\n");
      return CLMC_ERRO_PRECISAO;
   }
   /* Calcular fracao da energia total no sitio n para cada n */
   for(n = 1; n <= N; ++n)
   cadeia[n].densidade_energia = cadeia[n].energia / H0;
   /* Calcular dispersao da energia na cadeia */
   sigma = Z = 0.0;
   for(n = 1; n <= N; ++n){
      aux = (double)(n - N2);
      sigma += aux * aux * cadeia[n].densidade_energia;
      Z += cadeia[n].densidade_energia*cadeia[n].densidade_energia;
   }
   sigma = sqrt(sigma);
   for(n = 1; n <= N; ++n)
   fprintf(arquivo_energia, "%g %d %g\n", t, n, cadeia[n].energia);
   fprintf(arquivo_energia, "\n");
   fprintf(arquivo_dispersao, "%5.4g %g %g\n", t, sigma, Z);


   return CLMC_SUCESSO;
}

void __fechar_arquivos(void){
   fclose(arquivo_energia);
   fclose(arquivo_dispersao);
   return;
}
