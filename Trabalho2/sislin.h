#ifndef __SISLIN_H__
#define __SISLIN_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"

// Estrutura de dados para agrupar matrizes para decomposição DLU
typedef struct {
    real_t **D;
    real_t **L;
    real_t **U;
} DLU_matrices_t;

// Funções auxiliares
real_t produtoEscalar(real_t *v1, real_t *v2, int n);
void multMatVet(real_t **mat, real_t *vetIn, real_t *vetOut, int n);
void copiaVetor(real_t *origem, real_t *destino, int n);
void validaParam (int n, int k);

void criaKDiagonal(int n, int k, real_t ***A, real_t **B);
void destroiKDiagonal (int n, real_t **A, real_t *B);

void genSimetricaPositiva(real_t **A, real_t *b, int n, real_t ***ASP, real_t **bsp, rtime_t *tempo);


void geraDLU(real_t **A, int n, DLU_matrices_t **matrices, rtime_t *tempo);
void destroiDLU(DLU_matrices_t *matrices, int n);

void aplicaPreCond(real_t **A, real_t *r, real_t *z, int n, real_t omega, DLU_matrices_t *dlu);
real_t calcResiduoSL(real_t **A, real_t *b, real_t *X, int n, rtime_t *tempo);

#endif // __SISLIN_H__

