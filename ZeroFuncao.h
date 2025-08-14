#ifndef __ZEROFUNCAO_H__
#define __ZEROFUNCAO_H__

#include <float.h>

// Aproximação aceitável como valor zero
#define ZERO DBL_EPSILON

// Parâmetros para teste de convergência
#define MAXIT 500
#define EPS 1.0e-7
#define ULPS 3

typedef struct {
  real_t *p;
  int grau;
} Polinomio;

// Métodos
// Retornam valor do erro quando método finalizou. Este valor depende de tipoErro

real_t newtonRaphson (Polinomio p, real_t x0, int criterioParada, int *it, real_t *raiz);
real_t bisseccao (Polinomio p, real_t a, real_t b, int criterioParada, int *it, real_t *raiz);

// Cálculo de Polinômios
void calcPolinomio_rapido(Polinomio p, real_t x, real_t *px, real_t *dpx );
void calcPolinomio_lento(Polinomio p, real_t x, real_t *px, real_t *dpx );

#endif // __ZEROFUNCAO_H__

