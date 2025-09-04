#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> // Para usar memset

#include "utils.h"
#include "edo.h"
#include "gaussSeidel_EqDiff.h"

#define MAXIT 50

// Protótipos das funções que definem a EDO
real_t pp(real_t x);
real_t qq(real_t x);
real_t rr(real_t x);

int main() {
    // Malhas a serem testadas, conforme especificado no enunciado
    int n_values[] = {5, 10, 100, 1000};
    int num_n_values = sizeof(n_values) / sizeof(int);

    // Loop principal para cada tamanho de malha 'n'
    for (int i = 0; i < num_n_values; ++i) {
        int n = n_values[i];
        
        printf("============================================================\n");
        printf("=> EXECUTANDO PARA MALHA n = %d\n", n);
        printf("============================================================\n\n");

        EDo edo = {n, 0.0, 1.0, -1.0, 0.0, pp, qq, rr};
        real_t *Y = (real_t *)calloc(n, sizeof(real_t));
        if (!Y) {
            perror("Falha ao alocar vetor de solução Y");
            return 1;
        }

        // --- VERSÃO 1: Gauss-Seidel com Vetores ---
        printf("--- VERSÃO 1: Gauss-Seidel com Vetores ---\n");
        Tridiag *sl = genTridiag(&edo);

        // ALTERAÇÃO: Imprime apenas se n for MENOR que 10
        if (n < 10) {
            printf("Matriz Aumentada (Di, D, Ds, B):\n");
            prnTriDiagonal(sl);
        }

        rtime_t tempo1 = gaussSeidel_3Diag(sl, Y, MAXIT);
        real_t norma1 = normaL2_3Diag(sl, Y);

        // ALTERAÇÃO: Imprime apenas se n for MENOR que 10
        if (n < 10) {
            printf("\nSolução Y:\n");
            prnVetor(Y, n);
        }
        
        printf("\nNorma L2 do Resíduo: %e\n", norma1);
        printf("Tempo de execução: %lfs\n\n", tempo1);

        free(sl->D); free(sl->Di); free(sl->Ds); free(sl->B);
        free(sl);

        // --- VERSÃO 2: Gauss-Seidel com Cálculo Direto ---
        printf("--- VERSÃO 2: Gauss-Seidel com Cálculo Direto ---\n");

        memset(Y, 0, n * sizeof(real_t));
        
        rtime_t tempo2 = gaussSeidel_EDO(&edo, Y, MAXIT);
        real_t norma2 = normaL2_EDO(&edo, Y);

        // O vetor da solução da Versão 2 não precisa ser impresso, pois é o mesmo.
        // Se quisesse imprimir, a condição seria if (n < 10) aqui também.
        
        printf("\nNorma L2 do Resíduo: %e\n", norma2);
        printf("Tempo de execução: %lfs\n\n", tempo2);
        
        free(Y);
    }

    return 0;
}

// Implementação das funções p(x), q(x) e r(x) da EDO
real_t pp(real_t x) {
    return x + 1.0;
}

real_t qq(real_t x) {
    return -2.0 * x;
}

real_t rr(real_t x) {
    return (1.0 - x * x) * exp(-x);
}