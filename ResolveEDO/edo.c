#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "utils.h"
#include "edo.h"
#include "gaussSeidel_EqDiff.h"

#define MAXIT 50

real_t pp (real_t x);
real_t qq (real_t x);
real_t rr (real_t x);

int main ()
{

  EDo edo = {5, 0, 1, -1, 0, pp, qq, rr}; // definição da EDO do exercício
  real_t *Y; // Resultado da EDO.

  // aplica gauss-seidel para malhas 5, 10, 100 e 1000
  // Mostra resultados
  
  return 0;
}

real_t pp (real_t x)
{
  return x+1;
}

real_t qq (real_t x)
{
  return -2*x;
}

real_t rr (real_t x)
{
  // Escolher apenas um dos retornos abaixo. Qual o melhor?
  // return (1-x*x)*exp(-x);
  // return (1+x)*(1-x)*exp(-x);
}


