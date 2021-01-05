#ifndef CLMC_C
#define CLMC_C 0

/* **********
Tamanho da cadeia, quantidade de massas interagentes (inteiro positivo)
********** */
#define N 600

/* **********
Criterio de correlacao, define como serao escolhidas as massas.
Escolha um da lista:
1 - Mapa de Bernoulli
2 - Serie que o autor usou em seu 1o artigo
3 - Serie Carlos
********** */
#define __CRITERIO__ 2

/* **********
Fator de correlacao (real positivo)
(No caso do criterio ser 3 alpha deve ser nao nulo)
********** */
#define __ALPHA__ 3.0

/* **********
Termo de acoplamento (real positivo)
********** */
#define __ETA2__ 1.0
#define __ETA3__ 1.0
#define __ETA4__ 0.0

/* **********
Velocidade inicial (real)
********** */
#define __V0__ 1.0

/* **********
Metodo de resolucao das equacoes de Hamilton, note que se eta3 e eta4 forem
nulos entao o metodo por padrao eh a diagonalizacao, para o caso geral serah
obedecida a seguinte lista:
RK4 - Runge-Kutta classico de 4a ordem
RK8 - Runge-Kutta de 8a ordem
RK14 - Runge-Kutta de 14a ordem
********** */
#define __METODO__ RK14

#endif
