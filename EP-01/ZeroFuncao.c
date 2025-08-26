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
real_t newtonRaphson (Polinomio pol, real_t x0, CriterioParada cp, int *it, real_t *raiz, bool ehCalcLento) {
    if (it) 
        *it = 0;

    real_t x_new = x0;
    real_t x_old;
    real_t fx, dfx;
    real_t erro = INFINITY; // Começa com erro máximo

    bool critParada = false;

    // Função de cálculo de polinômio (lento ou rápido)
    real_t (*calcPol)(Polinomio, real_t, real_t*, real_t*, bool) = calcPolinomio;

    while (!critParada && (*it < MAXIT)) {
        x_old = x_new;

        calcPol(pol, x_old, &fx, &dfx, ehCalcLento);

        // Verifica se a derivada é próxima de zero para evitar divisão por zero
        if (fabs(dfx) < ZERO) {
            fprintf(stderr, "ERRO: Derivada próxima de zero em Newton-Raphson.\n");
            *raiz = x_old;
            return NAN; 
        }

        // Fórmula de Newton-Raphson
        x_new = x_old - (fx / dfx);
        
        *it += 1;

        // Controla qual critério é utilizado para fim do laço
        switch (cp) {
            case CRIT1:
                // Evita divisão por zero se x_new for muito pequeno
                if (fabs(x_new) > ZERO)
                    erro = fabs(x_new - x_old) / fabs(x_new);
                else
                    erro = fabs(x_new - x_old);
                critParada = (erro <= EPS);
                break;
            case CRIT2:
                // O erro é o valor absoluto da função no novo ponto
                calcPol(pol, x_new, &fx, &dfx, ehCalcLento);
                erro = fabs(fx);
                critParada = (erro <= ZERO);
                break;
            case CRIT3:
                erro = ulp_distance(x_new, x_old);
                critParada = (erro <= ULPS);
                break;
        }
    }

    *raiz = x_new;
    return erro;
}


// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t bisseccao (Polinomio pol, real_t a, real_t b, CriterioParada cp, int *it, real_t *raiz, bool ehCalcLento) {
    if (it) 
        *it = 0;

    // Garante que a < b para o intervalo [a, b]
    if (a > b) {
        real_t aux = a;
        a = b;
        b = aux;
    }

    real_t xm_new, xm_old = a; 
    real_t fa, fx, dfx;
    real_t erro = INFINITY;

    // Função de cálculo de polinômio (lento ou rápido)
    real_t (*calcPol)(Polinomio, real_t, real_t*, real_t*, bool) = calcPolinomio;

    calcPol(pol, a, &fa, &dfx, ehCalcLento);

    // Laço principal
    while (*it < MAXIT) {
        xm_new = (a + b) / 2.0;
        calcPol(pol, xm_new, &fx, &dfx, ehCalcLento);
        
        *it += 1;

        // Verifica o critério de parada
        bool critParada = false;
        switch (cp) {
            case CRIT1:
                if (fabs(xm_new) > ZERO)
                    erro = fabs(xm_new - xm_old) / fabs(xm_new);
                else
                    erro = fabs(xm_new - xm_old);
                critParada = (erro <= EPS);
                break;
            case CRIT2:
                erro = fabs(fx);
                critParada = (erro <= ZERO);
                break;
            case CRIT3:
                erro = ulp_distance(xm_new, xm_old);
                critParada = (erro <= ULPS);
                break;
        }

        if (critParada) {
            *raiz = xm_new;
            return erro;
        }
        
        // Atualiza o intervalo
        if (fa * fx < 0) {
            b = xm_new; // A raiz está no intervalo [a, xm]
        } else {
            a = xm_new; // A raiz está no intervalo [xm, b]
            fa = fx;    // Otimização: f(a) para a próxima iteração será o f(x) atual
        }
        
        xm_old = xm_new;
    }

    *raiz = xm_new;
    return erro;
}


real_t calcPolinomio_rapido(Polinomio pol, real_t x, real_t *px, real_t *dpx)
{
    *px = pol.p[pol.grau];
    *dpx = 0.0;
    
    for (int i = pol.grau - 1; i >= 0; i--) {
        *dpx = (*dpx) * x + (*px);
        *px = (*px) * x + pol.p[i];
    }

    return *px;
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
