#ifndef __SISLIN_H__
#define __SISLIN_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "utils.h"

// Estrutura Otimizada V2 (DIA Format)
typedef struct {
    real_t *diagonals; 
    int *offsets;      
    int n;             
    int k;             
} MatrizDiagonal_t;

real_t produtoEscalar(real_t *v1, real_t *v2, int n);
void multMatVet(MatrizDiagonal_t *mat, real_t *vetIn, real_t *vetOut);
void copiaVetor(real_t *origem, real_t *destino, int n);
void validaParam(int n, int k);

MatrizDiagonal_t* alocaMatrizDiagonal(int n, int k, int *offsets);
void liberaMatrizDiagonal(MatrizDiagonal_t *mat);

void criaKDiagonal(int n, int k, MatrizDiagonal_t **A, real_t **B);
void genSimetricaPositiva(MatrizDiagonal_t *A, real_t *b, int n, MatrizDiagonal_t **ASP, real_t **bsp, rtime_t *tempo);

void aplicaPreCond(MatrizDiagonal_t *A, real_t *r, real_t *z, int n, real_t omega);
real_t calcResiduoSL(MatrizDiagonal_t *A, real_t *b, real_t *X, int n, rtime_t *tempo);
real_t resolveSL(MatrizDiagonal_t *A, real_t *b, real_t **X, int n, int maxit, real_t epsilon, real_t omega, rtime_t *tempo_iter, int *out_total_iterations);

#endif