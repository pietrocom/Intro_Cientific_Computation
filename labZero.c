#include <stdio.h>
#include <math.h>
#include <float.h>

#include "utils.h"
#include "ZeroFuncao.h"

int main ()
{

  real_t a, b;
  Polinomio pol;

  scanf("%d", &pol.grau);

  for (i=pol.grau; i >=0; --i)
    scanf("%lf", &pol.p[i]);

  scanf("%lf %lf", &a, &b); // intervalo onde está uma das raizes.


  // Restante do programa a partir daqui

  return 0;
}

