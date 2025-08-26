// labZero.c
#include <stdio.h>
#include <stdlib.h> // Para malloc
#include <math.h>
#include <float.h>
#include <fenv.h> // Para fesetround

#include "utils.h"
#include "ZeroFuncao.h"

int main ()
{
    // Requisito do enunciado: arredondamento para baixo
    fesetround(FE_DOWNWARD);

    Polinomio pol;
    real_t a, b; // Intervalo para bisseção
    real_t x0;   // Chute inicial para Newton 

    // Lê o grau do polinômio
    scanf("%d", &pol.grau);

    // Aloca memória para os coeficientes
    pol.p = malloc((pol.grau + 1) * sizeof(real_t));
    if (!pol.p) {
        perror("Falha ao alocar memória para o polinômio");
        return 1;
    }

    // Lê os coeficientes do maior para o menor grau
    for (int i = pol.grau; i >= 0; --i) {
        scanf("%lf", &pol.p[i]);
    }

    // Lê o intervalo
    scanf("%lf %lf", &a, &b);
    x0 = (a + b) / 2.0; // Usa o meio do intervalo como chute inicial para Newton

    // Executa para cálculo rápido e lento
    for (int i = 0; i < 2; ++i) {
        bool ehCalcLento = (i == 1);

        printf("%s\n\n", ehCalcLento ? "LENTO" : "RAPIDO");

        // Itera sobre os critérios de parada
        for (CriterioParada cp = CRIT1; cp <= CRIT3; ++cp) {
            int iter_b, iter_n;
            real_t raiz_b, raiz_n;
            real_t erro_b, erro_n;
            double t_inicio, t_fim;

            // Bisseccao 
            t_inicio = timestamp();
            erro_b = bisseccao(pol, a, b, cp, &iter_b, &raiz_b, ehCalcLento);
            t_fim = timestamp();
            printf("bissec  %+.15e %.15e %4d  %.8e\n", raiz_b, erro_b, iter_b, t_fim - t_inicio);
        }

        for (CriterioParada cp = CRIT1; cp <= CRIT3; ++cp) {
            int iter_b, iter_n;
            real_t raiz_b, raiz_n;
            real_t erro_b, erro_n;
            double t_inicio, t_fim;

            // Newton-Raphson
            t_inicio = timestamp();
            erro_n = newtonRaphson(pol, x0, cp, &iter_n, &raiz_n, ehCalcLento);
            t_fim = timestamp();
            printf("newton  %+.15e %.15e %4d  %.8e\n", raiz_n, erro_n, iter_n, t_fim - t_inicio);
        }
        printf("\n");
    }

    free(pol.p);

    return 0;
}