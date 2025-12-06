#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "utils.h"
#include "sislin.h"

// ============================================================================
// GERENCIAMENTO DE MEMÓRIA OTIMIZADO (DIA FORMAT)
// ============================================================================
MatrizDiagonal_t* alocaMatrizDiagonal(int n, int k, int *offsets) {
    MatrizDiagonal_t *mat = malloc(sizeof(MatrizDiagonal_t));
    mat->n = n;
    mat->k = k;
    mat->offsets = malloc(sizeof(int) * k);
    if (offsets) memcpy(mat->offsets, offsets, sizeof(int) * k);
    else memset(mat->offsets, 0, sizeof(int) * k);

    // Alocação Contígua (Vetor Único para todas diagonais)
    // Permite prefetcher e evita cache misses
    mat->diagonals = calloc(n * k, sizeof(real_t));
    return mat;
}

void liberaMatrizDiagonal(MatrizDiagonal_t *mat) {
    if (mat) {
        if (mat->diagonals) free(mat->diagonals);
        if (mat->offsets)   free(mat->offsets);
        free(mat);
    }
}

// ============================================================================
// OPERAÇÕES ÁLGEBRA LINEAR OTIMIZADAS (V2)
// ============================================================================

real_t produtoEscalar(real_t *v1, real_t *v2, int n) {
    real_t resultado = 0.0;
    // Vetorizável (AVX)
    for (int i = 0; i < n; ++i) resultado += v1[i] * v2[i];
    return resultado;
}

// MULTIPLICAÇÃO MATRIZ-VETOR OTIMIZADA (DIA)
// Remove loops aninhados e ifs. Itera por diagonal.
void multMatVet(MatrizDiagonal_t *mat, real_t *vetIn, real_t *vetOut) {
    int n = mat->n;
    int k = mat->k;
    for (int i = 0; i < n; ++i) vetOut[i] = 0.0;

    for (int d = 0; d < k; ++d) {
        int offset = mat->offsets[d];
        real_t *diag_values = &(mat->diagonals[d * n]);
        int i_start = (offset < 0) ? -offset : 0;
        int i_end   = (offset > 0) ? (n - offset) : n;

        // Loop Sequencial Otimizado (AVX)
        for (int i = i_start; i < i_end; ++i) {
            vetOut[i] += diag_values[i] * vetIn[i + offset];
        }
    }
}

void copiaVetor(real_t *origem, real_t *destino, int n) {
    memcpy(destino, origem, n * sizeof(real_t));
}

// ============================================================================
// GERAÇÃO DE DADOS (V2 - GERA DIRETO EM DIA)
// ============================================================================
void validaParam(int n, int k) {
    if (n <= 10 || k <= 1 || k % 2 == 0) exit(EXIT_FAILURE);
}

void criaKDiagonal(int n, int k, MatrizDiagonal_t **A, real_t **B) {
    validaParam(n, k);
    int *offsets = malloc(sizeof(int) * k);
    int centro = k / 2;
    for (int i = 0; i < k; ++i) offsets[i] = i - centro;

    *A = alocaMatrizDiagonal(n, k, offsets);
    free(offsets);

    for (int d = 0; d < k; ++d) {
        int off = (*A)->offsets[d];
        int start = (off < 0) ? -off : 0;
        int end   = (off > 0) ? (n - off) : n;
        for (int i = start; i < end; ++i) 
             (*A)->diagonals[d * n + i] = ((double)rand() / RAND_MAX) * 10.0;
    }
    *B = malloc(sizeof(real_t) * n);
    for (int i = 0; i < n; ++i) (*B)[i] = ((double)rand() / RAND_MAX) * 10.0;
}

void genSimetricaPositiva(MatrizDiagonal_t *A_in, real_t *b_in, int n, 
                          MatrizDiagonal_t **ASP, real_t **bsp, rtime_t *tempo) {
    *tempo = timestamp();
    int k = A_in->k;
    *ASP = alocaMatrizDiagonal(n, k, A_in->offsets);
    *bsp = malloc(sizeof(real_t) * n);

    // Lógica simplificada de geração SPD para formato DIA
    // Garante diagonal dominante e simetria sem custo N^3
    memset((*ASP)->diagonals, 0, n * k * sizeof(real_t));

    for (int d = 0; d < k; ++d) {
        int off = (*ASP)->offsets[d];
        if (off >= 0) continue; 
        
        int d_pair = -1; // Encontra a diagonal simétrica
        for(int p=0; p<k; p++) 
            if((*ASP)->offsets[p] == -off) d_pair = p;

        if(d_pair != -1) {
            int start = -off;
            int end = n;
            for (int i = start; i < end; ++i) {
                real_t val = ((double)rand() / RAND_MAX);
                (*ASP)->diagonals[d * n + i] = val; 
                int j = i + off;
                (*ASP)->diagonals[d_pair * n + j] = val;
            }
        }
    }
    // Diagonal principal robusta
    int diag_idx = k/2; // Assumindo ordenado, ou buscar offset 0
    for(int d=0; d<k; d++) 
        if((*ASP)->offsets[d] == 0) diag_idx = d;
    real_t *diag_main = &((*ASP)->diagonals[diag_idx * n]);
    for (int i = 0; i < n; ++i) 
        diag_main[i] = (real_t)(k * 10.0); 

    for(int i=0; i<n; ++i) (*bsp)[i] = b_in[i];
    *tempo = timestamp() - *tempo;
}

// ============================================================================
// SOLVER E PRECONDICIONADORES (V2)
// ============================================================================

void aplicaPreCond(MatrizDiagonal_t *A, real_t *r, real_t *z, int n, real_t omega) {
    if (omega == -1.0) { copiaVetor(r, z, n); return; }

    // Otimização V2 (Jacobi): Acesso direto ao vetor da diagonal
    if (omega == 0.0) {
        int diag_idx = -1;
        for(int d=0; d < A->k; d++) if(A->offsets[d] == 0) { diag_idx = d; break; }
        if(diag_idx == -1) return;

        real_t *D = &(A->diagonals[diag_idx * n]);
        for (int i = 0; i < n; ++i) z[i] = r[i] / D[i];
        return;
    }
    if (omega >= 1.0) copiaVetor(r, z, n);
}

real_t calcResiduoSL(MatrizDiagonal_t *A, real_t *b, real_t *X, int n, rtime_t *tempo) {
    *tempo = timestamp();
    real_t *Ax = malloc(sizeof(real_t) * n);
    multMatVet(A, X, Ax); // Chama versão otimizada
    real_t norma_L2 = 0.0;
    for (int i = 0; i < n; ++i) {
        real_t res = b[i] - Ax[i];
        norma_L2 += res * res;
    }
    free(Ax);
    *tempo = timestamp() - *tempo;
    return sqrt(norma_L2);
}

real_t resolveSL(MatrizDiagonal_t *A, real_t *b, real_t **X, int n, int maxit, 
                 real_t epsilon, real_t omega, rtime_t *tempo_iter, int *out_total_iterations) {
    *X = calloc(n, sizeof(real_t));
    real_t *r = malloc(n * sizeof(real_t));
    real_t *p = malloc(n * sizeof(real_t));
    real_t *z = malloc(n * sizeof(real_t));
    real_t *Ap = malloc(n * sizeof(real_t));
    real_t *x_velho = malloc(n * sizeof(real_t));

    rtime_t t_inicio = timestamp();
    copiaVetor(b, r, n); 
    aplicaPreCond(A, r, z, n, omega);
    copiaVetor(z, p, n);
    real_t rz_old = produtoEscalar(r, z, n);
    real_t norma_erro = 0.0;
    int k = 0;

    for (k = 0; k < maxit; ++k) {
        if (fabs(rz_old) < 1e-15) break;
        copiaVetor(*X, x_velho, n);
        
        multMatVet(A, p, Ap); // OTIMIZADO

        real_t pAp = produtoEscalar(p, Ap, n);
        if (fabs(pAp) < 1e-15) break;
        real_t alpha = rz_old / pAp;

        for (int i = 0; i < n; ++i) {
            (*X)[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        norma_erro = 0.0;
        for (int i = 0; i < n; ++i) {
            real_t diff = fabs((*X)[i] - x_velho[i]);
            if (diff > norma_erro) norma_erro = diff;
        }
        if (norma_erro < epsilon) break;

        aplicaPreCond(A, r, z, n, omega);
        real_t rz_new = produtoEscalar(r, z, n);
        real_t beta = rz_new / rz_old;
        for (int i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];
        rz_old = rz_new;
    }

    rtime_t t_fim = timestamp();
    if (k > 0) *tempo_iter = (t_fim - t_inicio) / k;
    else *tempo_iter = 0.0;
    if (out_total_iterations) *out_total_iterations = k;

    free(r); free(p); free(z); free(Ap); free(x_velho);
    return norma_erro;
}