/* *****************************************************************************
   Funcao para calcular a energia na n-esima particula
   *****************************************************************************
   E-mail: ismaellxd@gmail.com
   Site: https://ismaeldamiao.github.io/
   *****************************************************************************
   Copyright © 2020 Ismael Damião

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
const double __um_sexto__ = 1.0/6.0;
double energia(int n){
   int i = n + 1;
   double aux[4] = {x[i] - x[n], x[n] - x[n-1], 0.0, 0.0};
   aux[2] = aux[0]*aux[0];
   aux[3] = aux[1]*aux[1];
   return
   0.5 * P[n] * P[n] / M[n] +
   (eta[i][0] * aux[2] + eta[n][0] * aux[3]) * 0.25 +
   (eta[i][1] * aux[2]*aux[0] + eta[n][1] * aux[3]*aux[1]) * __um_sexto__ +
   (eta[i][2] * aux[2]*aux[2] + eta[n][2] * aux[3]*aux[3]) * 0.125;
}