/* *****************************************************************************
   Biblioteca de funcoes CLMC
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
#ifndef CLMC_H
#define CLMC_H 0

/* Essa estrutura guarda informacoes fisicas sobre uma particula, mais adiante
   irei declarar uma cadeia cheia de particulas como sendo um vetor (array)
   do tipo 'particula' (nome que dei para a estrutura). */
typedef struct __particula {
   double massa;
   double posicao;
   double momento;
   double acoplamento_linear;
   double acoplamento_quadratico;
   double acoplamento_cubico;
   double energia;
   double densidade_energia;
} particula;

/* Agora estou usando 'enum' para listar os tipos de erros que, eventualmente,
   podem occorer ao executar o programa. Isso eh uma maneira de depurar o
   programa, isto eh, isso ajuda a identificar possiveis falhas e erros.
   Por padrao o primeiro item (CLMC_SUCESSO) vale 0, o segundo vale 1, etc...
   Isso eh interessante pois em sistemas UNIX como Ubuntu ou Android eh comum
   que um processo retorne um valor diferente de zero sempre que algo dah
   errado e eh possivel verificar qual eh o valor retornado. */
enum status {
   CLMC_SUCESSO,
   CLMC_ERRO_ARQUIVO,
   CLMC_ERRO_DE_CONFIGURACAO,
   CLMC_ERRO_NAN,
   CLMC_ERRO_PRECISAO
};

/* Estou enumerando as diferentes maneiras que o programa consegue usar
   para resolver as equacoes diferenciais alem disso o tipo de variavel
   'metodos_numericos' irah guardar a informacao do metodo a ser utilizado,
   a escolha do metodo serah definida no arquivo de configuracao que eh
   lido ao executar o programa. */
typedef enum __metodos_numericos {
   RK4,
   RK8,
   RK14,
   ABM5,
   ABM10
} metodos_numericos;

/* Este programa pode usar ate tres ataques matematicos diferentes para gerar
   correlacoes, correlacoes podem ser utilizadas para gerar massas aleatorias
   correlacionadas. */
typedef enum __series_de_correlacao {
   transformada_de_Fourrier,
   mapa_de_Bernoulli,
   nao_sei_como_chamar_kkk
} series_de_correlacao;

/* As condicoes iniciais, os valores dos parametros e demais caracteristicas
   do sistema sao informadas atravez de um arquivo de configuracao,
   este struct guarda essas configuracoes. */
typedef struct __configuracao {
   int quantidade_de_particulas;
   metodos_numericos metodo_de_solucao;
   series_de_correlacao correlacao_das_massas;
   double fator_de_correlacao;
   double termo_de_acoplamento_linear;
   double termo_de_acoplamento_quadratico;
   double termo_de_acoplamento_cubico;
   double velocidade_inicial;
} configuracao;

/* *****************************************************************************
   Declaracao das funcoes e variaveis
***************************************************************************** */

#define DeltaDeKronecker(x, y) ((x) != (y) ? 0.0 : 1.0)

/* Cadeia serah uma array externa, uma variavel comum a todos os arquivos do
   programa. Nela serao guardadas as informacoes das particulas em uma
   cadeia classica unidimencional. */
/* Variaveis globais, isto eh, comum a todos os arquivos. */
extern particula *cadeia;
extern configuracao config;
extern double H0;
extern int N;

particula *cadeia;
configuracao config;
double H0;
int N;

/* Funcoes para alocar memoria para as arrays. */
double *vetor(int);
double **matriz(int,int);
particula *vetorp(int);

configuracao ler_configuracao(void);
void __abrir_arquivos(configuracao,int,const char*, ...);
void __fechar_arquivos(void);
int __escrever_arquivos(double);

/* Funcoes que preparam o sistema fisico (veja arquivos) */
void __massas(int);
void __posicoes(void);
void __momentos(void);
void __acoplamentos(void);

/* Funcoes relacionadas com o hamiltoniano. */
double __energia(int);
double __forca(int,double,double,double);
/* Esta funcao calcula a velocidade da n-esima particula da cadeia.
   Para mais detalhes veja a equacao 2 do arquivo CLMC.pdf */
#define __velocidade(n, P) ((P) / cadeia[n].massa)

/* Funcoes que resolvem numericamente as equacoes de Hamilton. */
int rk4(void);
int rk8(void);
int rk14(void);
int abm5(void);
int abm10(void);

#endif // CLMC_H
