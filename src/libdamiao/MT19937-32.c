/* *****************************************************************************
   A C- program for MT19937: Real number version
   genrand() generates pseudorandom real number (double)
   which is uniformly distributed on [0, 1] interval, for each
   call.
   ***
   This vesion use words with 32bit (or the size of unsigned int). In the first
   call genrand(seed) set initial values to the working area of 624 words.
   (seed is any 32-bit integer except for 0)
   The function genrand(*y) return the pseudorandom real and
   set y to the pseudorandom integer.
   ***
   The original code has write by Takuji Nishimura.
   *****************************************************************************
   E-mail: ismaellxd@gmail.com
   Site: https://ismaeldamiao.github.io/
   *****************************************************************************
   Copyright (c) 2021 I.F.F. dos SANTOS (Ismael Damiao)

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

#include <stdio.h>
#include <stdint.h> /* Use C99 or latter */

/* Period parameters */
#define WORD_SIZE 32
#define N 624
#define M 397
#define R UINT32_C(0x1f)
#define MATRIX_A UINT32_C(0x9908b0df) /* constant vector a */
#define UPPER_MASK UINT32_C(0x80000000) /* most significant w-r bits */
#define LOWER_MASK UINT32_C(0x7fffffff) /* least significant r bits */
/* Tempering parameters */
#define TEMPERING_MASK_D UINT32_C(0x7fffffff)
#define TEMPERING_MASK_B UINT32_C(0x9d2c5680)
#define TEMPERING_MASK_C UINT32_C(0xefc60000)
#define TEMPERING_SHIFT_U(y) (y >> 11)
#define TEMPERING_SHIFT_S(y) (y << 7)
#define TEMPERING_SHIFT_T(y) (y << 15)
#define TEMPERING_SHIFT_L(y) (y >> 18)

double mt19937(uint32_t *y)
{
   static uint32_t mt[N]; /* the array for the state vector */
   static int mti = N+1; /* mti==N+1 means mt[N] is not initialized */
   static uint32_t mag01[2] = {UINT32_C(0x0), MATRIX_A};
   if(mti >= N){ /* generate N words at one time */
      int kk;
      if(mti == N+1){ /* initializing the array with NONZERO seed */
         if(*y == UINT32_C(0)) *y = UINT32_C(4357); /* a default seed is used */
         /* setting the initial seeds to mt[N] using the generator Line 25 of
         table 1 in [KNUTH 1981, The Art of Computer Programming Vol. 2
         (2nd Ed.), pp102 */
         mt[0] = *y & UINT32_C(0xffffffff);
         for(mti = 1; mti < N; ++mti)
            mt[mti] = (69069 * mt[mti-1]) & UINT32_C(0xffffffff);
      }
      for(kk = 0; kk < N-M; ++kk){
         *y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
         mt[kk] = mt[kk+M] ^ (*y >> 1) ^ mag01[*y & 0x1];
      }
      for(; kk < N-1; ++kk){
         *y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
         mt[kk] = mt[kk + (M-N)] ^ (*y >> 1) ^ mag01[*y & 0x1];
      }
      *y = (mt[N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
      mt[N-1] = mt[M-1] ^ (*y >> 1) ^ mag01[*y & 0x1];
      mti = 0;
   }
   *y = mt[++mti];
   *y ^= TEMPERING_SHIFT_U(*y) & TEMPERING_MASK_D;
   *y ^= TEMPERING_SHIFT_S(*y) & TEMPERING_MASK_B;
   *y ^= TEMPERING_SHIFT_T(*y) & TEMPERING_MASK_C;
   *y ^= TEMPERING_SHIFT_L(*y);
   return ((double)(*y) / (double)UINT32_MAX);
}
