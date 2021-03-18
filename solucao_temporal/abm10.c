/* *****************************************************************************
   Resolver a EDO usando Adams-Bashforth-Moulton de 10 ordem
   iniciado com Runge-Kutta de 14a ordem e 35 estagios
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
/* *****************************************************************************
   Definicoes
***************************************************************************** */
/* Dicretiacao do Runge-Kutta */
#define dt 1.0e-2
#define dt2 0.5e-2
/* *****************************************************************************
   Funcao do calculo numerico
***************************************************************************** */
int abm10(void){

   /* Contadores */
   int i, j, l, n, contador = 0;
   const int contadorMAX = (int)(2.0/dt);
   /* Medidas de localizacao */
   double *f;
   /* Runge-Kutta */
   const int s = 35, s_abm = 10; /* Quantidade de estagios */
   double **kP, **kX, **coef;
   double **kP_abm, **kX_abm, *x_ab, *P_ab;
   double a[s][s], b[s], c[s];
   double t;
   const double tf = 0.4 * (double)N;
   /* Tamanho 'dinamico' da cadeia */
   int n0 = 1, nf = N;
   /* Energia e hamiltoniano */
   double H;
   /* *** */
   const int N2 = N/2;
   double aux;
   double adams[2][s_abm] = 
   {// 10a ordem
      // Coeficientes de Adams–Bashforth
      {-2082753.0/7257600.0, 20884811.0/7257600.0, -94307320.0/7257600.0,
      252618224.0/7257600.0, -444772162.0/7257600.0, 538363838.0/7257600.0,
      -454661776.0/7257600.0, 265932680.0/7257600.0, -104995189.0/7257600.0,
      30277247.0/7257600.0},
      // Coeficientes de Adams–Moulton
      {57281.0/7257600.0, -583435.0/7257600.0, 2687864.0/7257600.0,
      -7394032.0/7257600.0, 13510082.0/7257600.0, -17283646.0/7257600.0,
      16002320.0/7257600.0, -11271304.0/7257600.0, 9449717.0/7257600.0,
      2082753.0/7257600.0}
   };


   /* **********
      Ler a matriz de Runge-Kutta
   ********** */
   #include"rk14.h"
   for(i = 0; i < s; ++i){
      b[i] *= dt;
      for(j = 0; j < s; ++j) a[i][j] *= dt;
   }
   for(i = 0; i < s_abm; ++i){
      adams[0][i] *= dt;
      adams[1][i] *= dt;
   }


   /* **********
      A cadeia possui um temanho dinamico da seguinte forma:
      Se a energia nos sitios n0 e nf for pequena, quase nula, entao
      entende-se que nao ha grande interferencia na dinamica do sistema
      por parte desses sitios.
      Se essa energia for rezoavelmente diferente de zero entao entao
      o tamanho da cadeia irah crezcer mudando os valores de n0 e nf para
      valores onde novamente essa energia pode ser considerada nula.
   ********** */
   #if N > 600
   n0 = N2 - 300;
   nf = N2 + 300;
   #endif


   /* ***
   Alocar memoria para os vetores
   *** */
   i = N+2;
   E = dvetor(i);
   f = dvetor(i);
   kX = dmatriz(i, s+1);
   kP = dmatriz(i, s+1);
   kX_abm = dmatriz(i, s_abm);
   kP_abm = dmatriz(i, s_abm);
   x_ab = dvetor(i);
   P_ab = dvetor(i);
   coef = dmatriz(i, 2);

   /* ***
   Zerar coeficientes que poderiam atrapalhar o calculo
   *** */
   for(n = 0; n <= N; ++n) coef[n][0] = coef[n][1] = 0.0;


   for(l = 0; l < s_abm; ++l){
      /* ***********************************************************************
         Rotina do metodo de Runge-kutta
      *********************************************************************** */
      for(i = 0; i < s; ++i){
      
         for(n = n0; n <= nf; ++n) coef[n][0] = coef[n][1] = 0.0;

         for(n = n0; n <= nf; ++n){
            for(j = 0; j < i; ++j){
               coef[n][0] += a[i][j] * kX[n][j];
               coef[n][1] += a[i][j] * kP[n][j];
            }
         }
         for(n = n0; n <= nf; ++n){
            kP[n][i] = forca(n,
            x[n-1] + coef[n-1][0], x[n] + coef[n][0], x[n+1] + coef[n+1][0]);
            kX[n][i] = velocidade(n, P[n] + coef[n][1]);
         }
      }
      for(n = n0; n <= nf; ++n){
         for(i = 0; i < s; ++i){
            P[n] += b[i] * kP[n][i];
            x[n] += b[i] * kX[n][i];
         }
      }
      for(n = n0; n <= nf; ++n){
         kP_abm[n][l] = forca(n, x[n-1], x[n], x[n+1]);
         kX_abm[n][l] = velocidade(n, P[n]);
      }
   }
   j = s_abm-1;
   for(t = (double)(s_abm)*dt; t <= tf; t += dt){
      /* ***********************************************************************
         Rotina do metodo de Adams-Bashforth-Moulton
      *********************************************************************** */
      for(n = n0; n <= nf; ++n){
         P_ab[n] = P[n];
         x_ab[n] = x[n];
      }
      for(n = n0; n <= nf; ++n){ //AB
         for(i = 0; i < s_abm; ++i){
            P_ab[n] += adams[0][i] * kP_abm[n][i];
            x_ab[n] += adams[0][i] * kX_abm[n][i];
         }
      }
      for(n = n0; n <= nf; ++n){
         for(i = 0; i < j; ++i){
            kP_abm[n][i] = kP_abm[n][i+1];
            kX_abm[n][i] = kX_abm[n][i+1];
         }
         kP_abm[n][j] = forca(n, x_ab[n-1], x_ab[n], x_ab[n+1]);
         kX_abm[n][j] = velocidade(n, P_ab[n]);
      }
      for(n = n0; n <= nf; ++n){ //AM
         for(i = 0; i < s_abm; ++i){
            P[n] += adams[1][i] * kP_abm[n][i];
            x[n] += adams[1][i] * kX_abm[n][i];
         }
      }
      /* ***********************************************************************
         Rotina para calcular as medidas de localizacao
      *********************************************************************** */
      if(++contador > contadorMAX){
         contador = 0;
         /* ***
         Calculo do hamiltoniano
         *** */
         H = 0.0;
         for(n = 1; n <= N; ++n){
            E[n] = energia(n);
            H += E[n];
         }

         /* ***
         Verificar se nao houve problema de numero ilegal ou falta de precisao
         *** */
         if(H != H){
            fprintf(stderr, "Algum calculo resultou em nan (=nao eh numero)\n");
            return CLMC_ERRO_NAN;
         }
         if(fabs(1.0 - H / H0) > 1.0e-8){
            fprintf(stderr, "Falta de precisao nos calculos\n");
            return CLMC_ERRO_PRECISAO;
         }

         /* ***
            Calcular fracao da energia total no sitio n para cada n
         *** */
         for(n = 1; n <= N; ++n) f[n] = E[n] / H0;

         /* ***
            Calcular dispersao da energia na cadeia
         *** */
         sigma = Z = 0.0;
         for(n = 1; n <= N; ++n){
            aux = (double)(n - N2);
            sigma += aux * aux * f[n];
            Z += f[n]*f[n];
         }
         sigma = sqrt(sigma);

         #if N > 600
         /* ***
            Crescer tamanho da cadeia caso energia nas bordas seja relevante
         *** */
         if((E[n0] > 1.0e-20) || (E[nf] > 1.0e-20)){
            n0 -= 30; nf += 30;
         }
         /* Para garantir que a sub-cadeia nao fique maior que a cadeia */
         if(nf > N){
            n0 = 1; nf = N;
         }
         #endif

         /* ***
         Ver arquivo CLMC_arquivos.c
         *** */
         EscreverArquivos(t);
      }
   }
   return CLMC_SUCESSO;
}
