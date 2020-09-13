#ifndef CLMC_C
#define CLMC_C 0

/* **********
Tamanho da cadeia, quantidade de massas interagentes
********** */
#define N 800

/* **********
Semente do gerador de numeros aleatorios
********** */
#define Semente 1

/* **********
Fator de correlacao
********** */
#define alpha 0.0

/* **********
Criterio de correlacao
1 - Mapa de Bernoulli
2 - artigo
3 - Carlos
********** */
#define criterio 2

/* **********
Termo de acoplamento
********** */
#define eta2 1.0
#define eta3 1.0
#define eta4 0.0

/* **********
Velocidade inicial
********** */
#define v0 1.0

/* **********
Metodo de resolucao das equacoes de Hamilton, note que se eta3 e eta4 forem
nulos entao o metodo por padrao eh a diagonalizacao, para o caso geral serah
obedecida a seguinte lista:
RK4 - Runge-Kutta classico de 4a ordem
RK8 - Runge-Kutta de 8a ordem
RK14 - Runge-Kutta de 14a ordem
********** */
#define metodo RK14

#endif
