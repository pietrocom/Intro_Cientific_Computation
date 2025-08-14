#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#include "utils.h"
#include "ZeroFuncao.h"

// Função auxiliar que calcula f(x) do polinômio em questão
real_t func (Polinomio pol, real_t x) {
    real_t soma = 0;
    for (int i = pol.grau; i >= 0; i--) {
        soma += pol.p[i];
    }
    return soma;
}

// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t newtonRaphson (Polinomio pol, real_t x0, int criterioParada, int *it, real_t *raiz)
{

}


// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t bisseccao (Polinomio pol, real_t a, real_t b, CriterioParada wa, int *it, real_t *raiz, bool ehCalcLento)
{
    real_t xm_new, xm_old;

    xm_new = (a + b) / 2;

    real_t fa, dfa, fx, dfx;

    if (ehCalcLento) {
        calcPolinomio_lento(pol, a, &fa, &dfa);
        calcPolinomio_lento(pol, xm_new, &fx, &dfx);
    }
    else {
        calcPolinomio_lento(pol, a, &fa, &dfa);
        calcPolinomio_lento(pol, xm_new, &fx, &dfx);
    }

    switch (wa) {
        case CRIT1:

        case CRIT2:

        case CRIT3:
            
    }

    if (fa * fx < 0) 
        a = xm_new;
    else if (fa * fx > 0)
        b = xm_new;

}


void calcPolinomio_rapido(Polinomio pol, real_t x, real_t *px, real_t *dpx)
{

}


void calcPolinomio_lento(Polinomio pol, real_t x, real_t *px, real_t *dpx)
{
    for (int i = pol.grau; i >= 0; i--) {
        *px += pol.p[i] * x;
        *dpx += *px * i;
    }
}
