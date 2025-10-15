// pc24 - GRR20241955

#include <stdio.h>
#include <stdlib.h>
#include <fenv.h> 
#include <likwid.h>
#include "utils.h"
#include "edo.h"

int main() {
    // Define o arredondamento de ponto flutuante para baixo
    fesetround(FE_DOWNWARD);

    EDo edo;

    // Leitura dos parâmetros comuns da família de EDOs
    if (scanf("%d", &edo.n) != 1) return 1;
    if (scanf("%lf %lf", &edo.a, &edo.b) != 2) return 1;
    if (scanf("%lf %lf", &edo.ya, &edo.yb) != 2) return 1;
    if (scanf("%lf %lf", &edo.p, &edo.q) != 2) return 1;

    // Inicializa a biblioteca LIKWID
    LIKWID_MARKER_INIT;

    int edo_count = 0;
    // Laço para ler os coeficientes de r(x) para cada EDO
    while (scanf("%lf %lf %lf %lf", &edo.r1, &edo.r2, &edo.r3, &edo.r4) == 4) {
        
        // Imprime a ordem e a matriz aumentada do sistema linear
        prnEDOsl(&edo);

        // Gera o sistema linear tridiagonal para o solver
        Tridiag *sl = genTridiag(&edo);
        if (!sl) {
            fprintf(stderr, "Erro ao alocar memória para o sistema linear.\n");
            return 1;
        }

        // Aloca vetor para a solução Y
        real_t *Y = (real_t *) malloc(sl->n * sizeof(real_t));
        if (!Y) {
            fprintf(stderr, "Erro ao alocar memória para o vetor solução.\n");
            free(sl);
            return 1;
        }

        real_t residue_norm;
        int iterations;
        rtime_t startTime, endTime;
        
        char marker_name[32];
        sprintf(marker_name, "GaussSeidel_%d", edo_count);

        // Medição de tempo e desempenho com LIKWID
        startTime = timestamp();
        LIKWID_MARKER_START(marker_name);

        iterations = gaussSeidel(sl, Y, &residue_norm);
        
        LIKWID_MARKER_STOP(marker_name);
        endTime = timestamp();

        // Imprime os resultados
        for (int i = 0; i < sl->n; ++i) {
            printf(FORMAT, Y[i]);
        }
        printf("\n");
        printf("%d\n", iterations);
        printf(FORMAT, residue_norm);
        printf("\n");
        printf("%.8e\n", endTime - startTime);

        // Adiciona uma linha em branco entre as saídas de diferentes EDOs, se não for a última
        // Com base no exemplo, um pulo de linha extra é esperado
        printf("\n");

        // Libera a memória alocada para esta iteração
        free(Y);
        free(sl->D);
        free(sl->Di);
        free(sl->Ds);
        free(sl->B);
        free(sl);
        
        edo_count++;
    }

    // Finaliza a biblioteca LIKWID
    LIKWID_MARKER_CLOSE;

    return 0;
}