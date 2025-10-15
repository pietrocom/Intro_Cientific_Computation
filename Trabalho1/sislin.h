#ifndef __SISLIN_H__
#define __SISLIN_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"

// Funções auxiliares
real_t produtoEscalar(real_t *v1, real_t *v2, int n);
void multMatVet(real_t **mat, real_t *vetIn, real_t *vetOut, int n);
void copiaVetor(real_t *origem, real_t *destino, int n);
void validaParam (int n, int k);

void criaKDiagonal(int n, int k, real_t ***A, real_t **B);
void destroiKDiagonal (int n, real_t **A, real_t *B);

void genSimetricaPositiva(real_t **A, real_t *b, int n, real_t ***ASP, real_t **bsp, rtime_t *tempo);


void geraDLU (double *A, int n, int k, double **D, double **L, double **U, double *tempo);
void geraPreCond(double *D, double *L, double *U, double w, int n, int k, double **M, double *tempo);
double calcResiduoSL (double *A, double *b, double *X, int n, int k, double *tempo);

#endif // __SISLIN_H__

