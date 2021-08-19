/* *****************************************************************************
   Esta funcao ler um arquivo de configuracao que contem algumas informacoes
   sobre o sistema fisico, caso o arquivo nao exista ele serah criado com
   valores padronizados.
   *****************************************************************************
   E-mail: ismaellxd@gmail.com
   Site: https://ismaeldamiao.github.io/
   *****************************************************************************
   Copyright (c) 2020 I.F.F. dos SANTOS (Ismael Damiao)

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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "./CLMC.h"

configuracao ler_configuracao(void){
   const char nome_padrao[] = "CLMC_config";
   char argumento[512], valor[512];
   int status;
   configuracao config;
   FILE *arquivo;

   arquivo = fopen(nome_padrao, "r");

   if(!arquivo){
      arquivo = fopen(nome_padrao, "w");
      fprintf(arquivo,
      "quantidade_de_particulas        200\n"
      "metodos_numericos               RK4\n"
      "correlacao_das_massas           transformada_de_Fourrier\n"
      "fator_de_correlacao             1\n"
      "termo_de_acoplamento_linear     1\n"
      "termo_de_acoplamento_quadratico 1\n"
      "termo_de_acoplamento_cubico     0\n"
      "velocidade_inicial              1");
      fclose(arquivo);
      arquivo = fopen(nome_padrao, "r");
   }
   //printf("OK\n");

   for(;;){
      status = fscanf(arquivo, "%s %s", argumento, valor);
      if(status == EOF) break;

      if(!strcmp(argumento, "quantidade_de_particulas")){
         config.quantidade_de_particulas = atoi(valor);
         N = config.quantidade_de_particulas;
      }else if(!strcmp(argumento, "metodos_numericos")){
         if(!strcmp(valor, "RK4")) config.metodo_de_solucao = RK4;
         else if(!strcmp(valor, "RK8")) config.metodo_de_solucao = RK8;
         else if(!strcmp(valor, "RK14")) config.metodo_de_solucao = RK14;
         else if(!strcmp(valor, "ABM5")) config.metodo_de_solucao = ABM5;
         else if(!strcmp(valor, "ABM10")) config.metodo_de_solucao = ABM10;
         else config.metodo_de_solucao = RK4;
      }else if(!strcmp(argumento, "correlacao_das_massas")){
         if(!strcmp(valor, "transformada_de_Fourrier"))
            config.correlacao_das_massas = transformada_de_Fourrier;
         else if(!strcmp(valor, "mapa_de_Bernoulli"))
            config.correlacao_das_massas = mapa_de_Bernoulli;
         else if(!strcmp(valor, "nao_sei_como_chamar_kkk"))
            config.correlacao_das_massas = nao_sei_como_chamar_kkk;
         else config.correlacao_das_massas = transformada_de_Fourrier;
      }else if(!strcmp(argumento, "fator_de_correlacao")){
         config.fator_de_correlacao = atof(valor);
      }else if(!strcmp(argumento, "termo_de_acoplamento_linear")){
         config.termo_de_acoplamento_linear = atof(valor);
      }else if(!strcmp(argumento, "termo_de_acoplamento_quadratico")){
         config.termo_de_acoplamento_quadratico = atof(valor);
      }else if(!strcmp(argumento, "termo_de_acoplamento_cubico")){
         config.termo_de_acoplamento_cubico = atof(valor);
      }else if(!strcmp(argumento, "velocidade_inicial")){
         config.velocidade_inicial = atof(valor);
      }
   }
   fclose(arquivo);
   return config;
}
