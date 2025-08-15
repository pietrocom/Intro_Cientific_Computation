#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#include "utils.h"
#include "ZeroFuncao.h"

// Função auxiliar que decide qual método de calculo de polinômio usar, retornando o resultado do mesmo
real_t calcPolinomio (Polinomio pol, real_t x, real_t *px, real_t *dpx, bool ehCalcLento) {

    return ehCalcLento ? calcPolinomio_lento(pol, x, px, dpx) : calcPolinomio_rapido(pol, x, px, dpx);
}

// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t newtonRaphson (Polinomio pol, real_t x0, int criterioParada, int *it, real_t *raiz)
{
    
}


// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t bisseccao (Polinomio pol, real_t a, real_t b, CriterioParada cp, int *it, real_t *raiz, bool ehCalcLento)
{
    if (it && (*it != 0))
        *it = 0;

    if (a < b) {
        real_t aux = a;
        a = b;
        b = aux;
    }

    real_t xm_new, xm_old;

    xm_new = (a + b) / 2;

    real_t fa, dfa, fx, dfx;

    calcPolinomio(pol, a, &fa, &dfa, ehCalcLento);
    calcPolinomio(pol, xm_new, &fx, &dfx, ehCalcLento);

    if (fa * fx < 0) 
        b = xm_new;
    else if (fa * fx > 0)
        a = xm_new;
    else
        return xm_new;

    bool critParada = false;

    // Laço principal
    while (!critParada && (*it < MAXIT)) {

        xm_old = xm_new;
        xm_new = (a + b) / 2;

        calcPolinomio(pol, a, &fa, &dfa, ehCalcLento);
        calcPolinomio(pol, xm_new, &fx, &dfx, ehCalcLento);

        if (fa * fx < 0) 
            b = xm_new;
        else if (fa * fx > 0)
            a = xm_new;
        else
            return xm_new;

        // Controla qual critério é utilizado para fim do laço
        switch (cp) {
            case CRIT1:
                critParada = (fabs(xm_new - xm_old) / fabs(xm_new)) <= 10e-7;
                break;
            case CRIT2:
                critParada = fabs(fx) <= DBL_EPSILON;
                break;
            case CRIT3:
                critParada = ulp_distance(xm_new, xm_old) <= 3;
                break;
        }

        *it += 1;
    }

}


real_t calcPolinomio_rapido(Polinomio pol, real_t x, real_t *px, real_t *dpx)
{

}


real_t calcPolinomio_lento(Polinomio pol, real_t x, real_t *px, real_t *dpx)
{
    *px = 0.0;
    *dpx = 0.0;
    for (int i = pol.grau; i > 0; i--) {
        *px += pol.p[i] * pow(x, i);
        *dpx += i * pol.p[i] * pow(x, i - 1);
    }
    *px += pol.p[0];

    return *px;
}
