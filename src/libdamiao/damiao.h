/* *****************************************************************************
   This file is a collection of functions and subroutine for use in
   physics, mathematics and engineering.
   This is conformed to the ISO/IEC 9899:1999 (also called C99), but
   contains intructions to compile with some GNU extensions if __GNUC__ is
   defined.
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

#ifndef DAMIAO_H
#define DAMIAO_H 1

/* *****************************************************************************
   C99 standard
***************************************************************************** */

#if !defined(__STDC__) || (__STDC_VERSION__ < 199901L) /* ISO C99 */
#error "Use only ISO C99 and later."
#endif /* ISO C99 */


#include <stdio.h>   /* Many functions for file input and output */
#include <stdlib.h>  /* Numeric conversion functions, memory allocation */
#include <stdbool.h> /* Defines a boolean data type */
#include <stdint.h>  /* Defines exact-width integer types. */
#include <math.h>    /* Many mathematical functions */
#include <complex.h> /* Many complex mathematical functions */
#include <float.h>   /* DBL_MAX */
#include <limits.h>  /* INT_MAX, INT_MIN */
#include <time.h>    /* Defines date- and time-handling functions */
#include <stdarg.h>  /* var-args */
/* There are not a good idea use tgmath.h because in a loop  */
#define _TGMATH_H 1
/* #include <tgmath.h> */


/* Standard pragmas */
#pragma STDC FENV_ACCESS OFF
#pragma STDC FP_CONTRACT ON
#pragma STDC CX_LIMITED_RANGE OFF

/* Some implementations of the C standards use functions in creal and cimag, but
   is best use simples operators (musl implemented like this).
   See
   * ISO/IEC 9899:1999 6.3.1.7 2
   * ISO/IEC 9899:1999 6.2.5 13 */
#ifdef creal
#undef creal
#endif
#define creal(z) ((double)(z))

#ifdef cimag
#undef cimag
#endif
static double _Complex __CIMAG;
#define cimag(z) (__CIMAG = (z), ((double*)&(__CIMAG))[1])




/* *****************************************************************************
   C11 and other newer versions of the C standard
***************************************************************************** */

#if (__STDC_VERSION__ >= 201112L) /* ISO C11 */

#if (__STDC_NO_COMPLEX__ == 1)
#error "Complex arithmetic is not supported."
#endif /* __STDC_NO_COMPLEX__ */

#else

#if defined(_Imaginary_I)
#define CMPLX(x, y) ((double _Complex)((double)(x) + _Imaginary_I *(double)(y)))
#else
#define CMPLX(x, y) ((double _Complex)((double)(x) + _Complex_I * (double)(y)))
#endif /* _Imaginary_I */

#endif /* ISO C11 */




/* *****************************************************************************
   GNU Compiler Collection API and extensions to the C standard
***************************************************************************** */

#if defined(__GNUC__)
/* Always optimize the code */
#pragma GCC optimize ("O3")
/* Always show all warning */
#pragma GCC diagnostic warning "-Wall"
/* Always show extra warning */
#pragma GCC diagnostic warning "-Wextra"
/* To test portabilit enable informes if I use GCC extentions. */
#pragma GCC diagnostic warning "-Wpedantic"
#endif /* __GNUC__ */




/* *****************************************************************************
   POSIX/UNIX extensions to the C standard
***************************************************************************** */

#if !defined(_POSIX_VERSION)
#define M_E             2.7182818284590452354   /* e */
#define M_LOG2E         1.4426950408889634074   /* log_2 e */
#define M_LOG10E        0.43429448190325182765  /* log_10 e */
#define M_LN2           0.69314718055994530942  /* log_e 2 */
#define M_LN10          2.30258509299404568402  /* log_e 10 */
#define M_PI            3.14159265358979323846  /* pi */
#define M_PI_2          1.57079632679489661923  /* pi/2 */
#define M_PI_4          0.78539816339744830962  /* pi/4 */
#define M_1_PI          0.31830988618379067154  /* 1/pi */
#define M_2_PI          0.63661977236758134308  /* 2/pi */
#define M_2_SQRTPI      1.12837916709551257390  /* 2/sqrt(pi) */
#define M_SQRT2         1.41421356237309504880  /* sqrt(2) */
#define M_SQRT1_2       0.70710678118654752440  /* 1/sqrt(2) */
#endif /* _POSIX_VERSION */




/* *****************************************************************************
   Linux kernel API
***************************************************************************** */

#if defined(__linux__)
#endif /* __linux__ */




/* *****************************************************************************
   Android API
***************************************************************************** */

#if defined(__ANDROID__)
#endif /* __ANDROID__ */




/* *****************************************************************************
   Some physical constants
***************************************************************************** */
/* exactly */
#define CONST_SPEED_OF_LIGHT      2.99792458e8              /* m / s */
#define CONST_AVOGADRO            6.02214076e23             /* 1 / mol */
#define CONST_PLANCK              6.62607015e-34            /* kg m^2 / s */
#define CONST_BOLTZMANN           1.380649e-23              /* J / K */
#define CONST_ELEMENTARY_CHARGE   1.602176634e-19           /* A s */
/* approximately */
#define CONST_GRAVITATIONAL       6.674e-11                 /* m^3 / (kg s^2) */
#define CONST_REDUCED_PLANCK      1.054571817646156472e-34  /* kg m^2 / s */
#define CONST_VACUUM_PERMITTIVITY 8.8541878128e-12          /* s^4A^2/(Kg m^3)*/
#define CONST_VACUUM_PERMEABILITY 1.25663706212e-6          /* Kg m^2/(s A)^2 */




/* *****************************************************************************
   ...
***************************************************************************** */

double _Complex atoc(const char *str);

/* Randon numbers. */
#if defined(UINT64_MAX)
double mt19937(uint64_t*);
#elif defined(UINT32_MAX)
double mt19937(uint32_t*);
#else
double mt19937(unsigned int *);
#endif
#define random(idum) mt19937(idum)
double *correlated_w_fourier(double,int,int);
double *correlated_w_bernoulli(double,int,int);
double *correlated_w_distance(double,int,int);

#include "rk8.c"
#include "rk14.c"

#endif /* DAMIAO_H */
