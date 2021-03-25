[![Tip Me via PayPal](https://img.shields.io/badge/PayPal-tip%20me-green.svg?logo=paypal)](https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=D66EM3DGU35EE&source=url)

# CLMC - beta 0.3
Cadeia Linear de Massas Correlacionadas

Este programa estuda como a energia se propaga em uma
cadeia linear (unidimencional) através da dinâmica dos modos de vibração da energia
e de medidas de localização.

É possível resolver o **problema de valor inicial** (PVI) com diversos metodos iterativos,
tais como:
* Runge-Kutta de 4a, 8a e 14a ordem;
* Adams-Bashforth-Moulton de 8a e 10a ordem.

Disponibilizei este programa como material de aprendizagem para quem tem interesse
nos diversos métodos de resolução numérica de equações diferenciais ordinárias
(EDO) associadas a um PVI. Estes métodos podem ser aplicados em equações mais "simples",
entretanto eu resolvi mostrar como aplica-los em um problema mais extenso para
resolver um sistema de `N` equações.

Para entender a física envolvida nesse, bem como para ver mais detalhes
sobre o programa, veja a [wiki](https://ismaeldamiao.github.io/pt/wiki/clmc/intro/)
no meu site.

O programa ainda está em versão de teste (beta) e nem todas as funções estão disponíveis ou são completamente funcionais.

## Download

Para baixar use o comando

```bash
wget https://github.com/ismaeldamiao/CLMC/archive/master.zip
```

Para descomprimir use o comando

```bash
unzip master.zip
```

## Compilando e executando

Antes de compilar se lembre de configurar seu sistema físico no arquivo [CLMC_configuracao.h](CLMC_configuracao.h).

Para compilar primeiro entre no diretório CLMC-master com `cd CLMC-master` e depois compile com um dos dois comandos a seguir.

* Com o clang
```bash
clang main.c -lm -o clmc
```

* Com o gcc
```bash
gcc main.c -lm -o clmc
```

Para executar se lembre de passar uma semente como argumento (`1` no exemplo abaixo).

```bash
./clmc 1
```

Se quiser calcular mais de uma semente para um mesmo sistema utilize o script [loop.sh](scripts/loop.sh) (depois de compilar) para executar,
informando a primeira e a última semente como argumentos (`1` até `10` no exemplo abaixo).

```bash
bash scripts/loop.sh 1 10
```

## Plotando

No diretório `scripts` também há um script para plotar os gráficos usando o GNUplot
(desnecessário informar que, para usá-lo, você deve ter o GNUplot instalado).

A nível de exemplo, suponha que ao executar o programa o diretório de dados que o programa criou se chama
`300Massas_3alpha_1V0_1eta2_1eta3_0eta4`, você deve navegar até esse diretório usando
`cd 300Massas_3alpha_1V0_1eta2_1eta3_0eta4`.

* Para plotar os dados referentes a uma unica semente (por exemplo, que vale `1`) use:
```bash
bash ../scripts/plot.sh 1
```

* Se você rodou várias sementes, antes de plotar tire a média delas:
```bash
wget -q https://github.com/ismaeldamiao/avulsos/raw/master/c/media/media.c
clang media.c -lm -o media.o # Ou use o gcc
./media.o "energia_*.dat"
./media.o "dispersao_*.dat"
bash ../scripts/plot.sh
```

## Configurando

Antes de compilar é possível mudar as condições iniciais do **problema de valor inicial** ou mudar o método como o programa vai resolver o problema.
É possível mudar também a maneira como as massas se autocorrelacionam. O arquivo [CLMC_configuracao.h](CLMC_configuracao.h) está devidamente comentado e é ele
que deve ser utilizado para configurar o programa, basta alterar o valor das macros (`#define`).

## Estudando

Ao olhar o código fonte sugiro que começe pelo arquivo [main.c](main.c) e tenha o mente o que quer aprender.
Por exemplo, se você quiser aprender sobre o
Runge-Kutta de 4ª ordem veja o [main.c](main.c) para entender o que o programa faz e depois veja o arquivo
[rk4.c](solucao_temporal/rk4.c), não há necessidade
de olhar todos os arquivos do programa.
