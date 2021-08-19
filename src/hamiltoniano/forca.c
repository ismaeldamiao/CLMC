/* *****************************************************************************
   Esta funcao calcula a forca imprimida ah n-esima particula da cadeia.

   Para mais detalhes veja a equacao 3 do arquivo CLMC.pdf
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
#include "../CLMC.h"

double __forca(int n, double xe, double x, double xd){
   int i;
   double aux[4];

   i = n - 1;
   /* Distancia entre particulas. */
   aux[0] = x - xe;
   aux[1] = xd - x;
   /* Distancia ao quadrado. */
   aux[2] = aux[0]*aux[0];
   aux[3] = aux[1]*aux[1];

   return
   cadeia[n].acoplamento_linear * aux[1] +
   cadeia[n].acoplamento_quadratico * aux[3] +
   cadeia[n].acoplamento_cubico * aux[3]*aux[1] -(
   cadeia[i].acoplamento_linear * aux[0] +
   cadeia[i].acoplamento_quadratico * aux[2] +
   cadeia[i].acoplamento_cubico * aux[2]*aux[0]);
}
