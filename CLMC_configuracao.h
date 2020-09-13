#ifndef CLMC_C
#define CLMC_C 0

/* **********
Tamanho da cadeia, quantidade de massas interagentes (inteiro positivo)
********** */
#define N 800

/* **********
Semente do gerador de numeros aleatorios (inteiro positivo)
********** */
#define Semente 1

/* **********
Criterio de correlacao, define como serao escolhidas as massas.
Escolha um da lista:
1 - Mapa de Bernoulli
2 - Serie que o autor usou em seu 1o artigo
3 - Serie Carlos
********** */
#define criterio 2

/* **********
Fator de correlacao (real positivo)
(No caso do criterio ser 3 alpha deve ser nao nulo)
********** */
#define alpha 0.0

/* **********
Termo de acoplamento (real positivo)
********** */
#define eta2 1.0
#define eta3 1.0
#define eta4 0.0

/* **********
Velocidade inicial (real)
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
