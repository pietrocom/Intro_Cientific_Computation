/*
Autor: Pietro Comin
GRR:   20241955
Versão: v1 (Original Instrumentada com LIKWID)
*/

#include <stdio.h>
#include <stdlib.h>
#include "sislin.h"
#include "utils.h"

// --- MACROS LIKWID (Essenciais para o relatório) ---
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

// Define se imprime a saída detalhada ou apenas o resumo
#define OUTPUT_VISUALIZATION_MODE 0

int main() {
    // 1. Inicializa a biblioteca LIKWID
    LIKWID_MARKER_INIT;

    int n, k, maxit;
    real_t omega, epsilon;
    real_t **A = NULL, **ASP = NULL;
    DLU_matrices_t *dlu_mats = NULL;
    real_t *B = NULL, *bsp = NULL, *X = NULL;
    
    rtime_t t_pc = 0.0, t_iter = 0.0, t_residuo = 0.0;
    real_t norma_erro = 0.0, norma_residuo = 0.0;
    int total_iterations = 0;
    
    // Lê entrada
    if (scanf("%d %d %lf %d %lf", &n, &k, &omega, &maxit, &epsilon) != 5) return 1;
    
    srandom(20252);

    // Gera estruturas (V1 usa matrizes densas real_t**)
    criaKDiagonal(n, k, &A, &B);
    
    // Gera a matriz simétrica. 
    // OBS: O sislin.c desta versão foi otimizado APENAS na geração para não demorar dias,
    // mas a estrutura de dados continua sendo a original do T1.
    genSimetricaPositiva(A, B, n, &ASP, &bsp, &t_pc);

    if (omega >= 1.0 && omega < 2.0) {
        rtime_t t_dlu = 0.0;
        geraDLU(ASP, n, &dlu_mats, &t_dlu); 
        t_pc += t_dlu;
    }

    // ---- OP1: ITERAÇÃO (Medição LIKWID) ----
    // Registra e inicia o marcador
    LIKWID_MARKER_REGISTER("OP1_ITERACAO");
    LIKWID_MARKER_START("OP1_ITERACAO");
    
    norma_erro = resolveSL(ASP, bsp, &X, n, maxit, epsilon, omega, &t_iter, dlu_mats, &total_iterations);
    
    LIKWID_MARKER_STOP("OP1_ITERACAO");

    // ---- OP2: RESÍDUO (Medição LIKWID) ----
    LIKWID_MARKER_REGISTER("OP2_RESIDUO");
    LIKWID_MARKER_START("OP2_RESIDUO");
    
    // Nota: Calcula resíduo usando a matriz original A
    norma_residuo = calcResiduoSL(A, B, X, n, &t_residuo);
    
    LIKWID_MARKER_STOP("OP2_RESIDUO");


    // ---- IMPRESSÃO (Mantém formato para script de comparação) ----
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

    // Limpeza
    destroiKDiagonal(n, A, B);
    destroiKDiagonal(n, ASP, bsp);
    destroiDLU(dlu_mats, n);
    if (X) free(X);

    // Fecha LIKWID
    LIKWID_MARKER_CLOSE;

    return 0;
}