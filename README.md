# CLMC - beta 0.2
Cadeia Linear de Massas Correlacionadas

Este programa calcula a dinâmica dos modos de vibração em uma cadeia linear (unidimencional), calcula medidas de localização da energia na cadeia.

É possível resolver o **problema de valor inicial** com diversos metodos iterativos e, em casos particulares, é possível encontrar os modos normais de
vibração e as autofrequências.

Para entender a física no problema veja [CLMC.pdf](CLMC.pdf).

O programa ainda está em versão de de teste (beta) e nem todas as funções estão disponíveis ou são completamente funcionais.

## Download

Para baixar use o comando

```
wget https://github.com/ismaeldamiao/CLMC/archive/master.zip
```

Para descomprimir use o comando

```
unzip master.zip
```

## Compilando

Para compilar entre no diretório CLMC-master e dê permissão de execução para o script [COMPILE](COMPILE)

```
chmod 755 COMPILE
```

E execute ele

```
./COMPILE
```

Se quiser executar o CLMC imediatamente após compilar use

```
./COMPILE x
```

## Configurando

Antes de compilar é possível mudar as condições iniciais do **problema de valor inicial** ou mudar o método como o programa vai resolver o problema.
É possível mudar também a maneira como as massas se autocorrelacionam. O arquivo [CLMC_configuracao.h](CLMC_configuracao.h) está devidamente comentado e é ele
que deve ser utilizado para configurar o programa, basta alterar o valor das macros (#define).

## Estudando

Ao olhar o código fonte sugiro que começe pelo arquivo [main.c](main.c) e tenha o mente o que quer aprender. Por exemplo, se você quiser aprender sobre o
Runge-Kutta de 4ª ordem veja o [main.c](main.c) para entender o que o programa faz e depois veja o arquivo [CLMC_rk4.c](CLMC_rk4.c), não há necessidade
de olhar todos os arquivos do programa.

## Plotando

Use os scripts do diretório gnuplot para plotar usando gnuplot. Por exemplo, se o arquivo de dados da energia gerado pelo CLMC se chamar
`Energia_800Massas_0alpha_1V0_1eta2_1eta3_0eta4_01semente.dat` então você pode usar o script [plot_energia.sh](gnuplot/plot_energia.sh) (lembre de colocar o script no mesmo diretório que o arquivo de dados) para plotar usando o comando

```
./plot_energia.sh Energia_800Massas_0alpha_1V0_1eta2_1eta3_0eta4_01semente.dat
```

A saída será uma imagem em formato `.png`.
