[![Tip Me via PayPal](https://img.shields.io/badge/PayPal-tip%20me-green.svg?logo=paypal)](https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=D66EM3DGU35EE&source=url)
[![LICENSE](https://img.shields.io/badge/license-MIT-lightgrey.svg)](/LICENSE)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/ismaeldamiao/CLMC)
![GitHub Release Date](https://img.shields.io/github/release-date/ismaeldamiao/CLMC)
![GitHub last commit](https://img.shields.io/github/last-commit/ismaeldamiao/CLMC)

# CLMC
Cadeia Linear de Massas Correlacionadas

Este programa estuda como a energia se propaga em uma
cadeia linear (unidimencional) através da dinâmica dos modos de vibração da energia
e de medidas de localização.

É possível resolver o **problema de valor inicial** (PVI) com diversos metodos iterativos,
tais como:
* Runge-Kutta de 4a, 8a e 14a ordem;
* Adams-Bashforth-Moulton de 5a e 10a ordem.

Disponibilizei este programa como material de aprendizagem para quem tem interesse
nos diversos métodos de resolução numérica de equações diferenciais ordinárias
(EDO) associadas a um PVI. Estes métodos podem ser aplicados em equações mais "simples",
entretanto eu resolvi mostrar como aplica-los em um problema mais extenso para
resolver um sistema de `N` equações.

Para entender a física envolvida nesse, bem como para ver mais detalhes
sobre o programa, veja [CLMC.pdf](tex/CLMC.pdf).

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
