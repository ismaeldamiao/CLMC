# CLMC - beta 0.1
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
Runge-Kutta de 4ª ordem veja o [main.c](main.c) para entender o que o programa faz e depois veja o arquivo [rk4.c](rk4.c), não há necessidade
de olhar todos os arquivos do programa.
