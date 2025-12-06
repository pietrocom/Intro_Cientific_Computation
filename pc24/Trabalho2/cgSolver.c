/*
Autor: Pietro Comin
GRR:   20241955
Versão: v1 do trabalho 2
*/

#include <stdio.h>
#include <stdlib.h>
#include "sislin.h"
#include "utils.h"

// --- MACROS LIKWID ---
#ifdef LIKWID_PERFMON
#include <likwid.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_SWITCH
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#define LIKWID_MARKER_GET(regionTag, nbEvents, events, time, count)
#endif

#define OUTPUT_VISUALIZATION_MODE 0

int main() {
    LIKWID_MARKER_INIT;

    int n, k, maxit;
    real_t omega, epsilon;
    
    // Uso da nova estrutura otimizada MatrizDiagonal_t
    MatrizDiagonal_t *A = NULL;
    MatrizDiagonal_t *ASP = NULL;
    real_t *B = NULL, *bsp = NULL, *X = NULL;
    
    rtime_t t_pc = 0.0, t_iter = 0.0, t_residuo = 0.0;
    real_t norma_erro = 0.0, norma_residuo = 0.0;
    int total_iterations = 0;
    
    if (scanf("%d %d %lf %d %lf", &n, &k, &omega, &maxit, &epsilon) != 5) return 1;
    
    srandom(20252);

    // Geração Otimizada (DIA)
    criaKDiagonal(n, k, &A, &B);
    genSimetricaPositiva(A, B, n, &ASP, &bsp, &t_pc);

    // ---- OP1: ITERAÇÃO (Medição LIKWID) ----
    LIKWID_MARKER_REGISTER("OP1_ITERACAO");
    LIKWID_MARKER_START("OP1_ITERACAO");
    
    // Solver usa multMatVet otimizado e acesso sequencial
    norma_erro = resolveSL(ASP, bsp, &X, n, maxit, epsilon, omega, &t_iter, &total_iterations);
    
    LIKWID_MARKER_STOP("OP1_ITERACAO");

    // ---- OP2: RESÍDUO (Medição LIKWID) ----
    LIKWID_MARKER_REGISTER("OP2_RESIDUO");
    LIKWID_MARKER_START("OP2_RESIDUO");
    
    norma_residuo = calcResiduoSL(A, B, X, n, &t_residuo);
    
    LIKWID_MARKER_STOP("OP2_RESIDUO");

    if (OUTPUT_VISUALIZATION_MODE == 0) {
        printf("%d\n", n);
        for (int i = 0; i < n; ++i) printf("%.16g ", X[i]);
        printf("\n");
        printf("%.8g\n", norma_erro);
        printf("%.8g\n", norma_residuo);
        printf("%.8g\n", t_pc);
        printf("%.8g\n", t_iter);
        printf("%.8g\n", t_residuo);
    }

    liberaMatrizDiagonal(A);
    if(B) free(B);
    liberaMatrizDiagonal(ASP);
    if(bsp) free(bsp);
    if (X) free(X);

    LIKWID_MARKER_CLOSE;
    return 0;
}