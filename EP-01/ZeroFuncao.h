#ifndef __ZEROFUNCAO_H__
#define __ZEROFUNCAO_H__

#include <float.h>
#include <stdbool.h>

// Aproximação aceitável como valor zero
#define ZERO DBL_EPSILON

// Parâmetros para teste de convergência
#define MAXIT 600
#define EPS 1.0e-7
#define ULPS 3

typedef struct {
  real_t *p;
  int grau;
} Polinomio;

typedef enum {
  CRIT1,
  CRIT2,
  CRIT3
} CriterioParada;

// Métodos
// Retornam valor do erro quando método finalizou. Este valor depende de tipoErro

real_t newtonRaphson (Polinomio pol, real_t x0, int criterioParada, int *it, real_t *raiz);
real_t bisseccao (Polinomio pol, real_t a, real_t b, int criterioParada, int *it, real_t *raiz, bool ehCalcLento);

// Cálculo de Polinômios
void calcPolinomio_rapido(Polinomio pol, real_t x, real_t *px, real_t *dpx );
real_t calcPolinomio_lento(Polinomio pol, real_t x, real_t *dx, real_t *dpx);

#endif // __ZEROFUNCAO_H__

