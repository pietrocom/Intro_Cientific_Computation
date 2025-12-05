#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"

// Gera valor aleatório para A
static inline real_t generateRandomA( unsigned int i, unsigned int j, unsigned int k )
{
  static real_t invRandMax = 1.0 / (real_t)RAND_MAX;
  return ( (i==j) ? (real_t)(k<<1) : 1.0 )  * (real_t)random() * invRandMax;
}

// Gera valor aleatório para B
static inline real_t generateRandomB( unsigned int k )
{
  static real_t invRandMax = 1.0 / (real_t)RAND_MAX;
  return (real_t)(k<<2) * (real_t)random() * invRandMax;
}

// Produto Escalar (V1)
real_t produtoEscalar(real_t *v1, real_t *v2, int n) {
  real_t resultado = 0.0;
  for (int i = 0; i < n; ++i) {
    resultado += v1[i] * v2[i];
  }
  return resultado;
}

// Multiplicação Matriz-Vetor (V1 - Naive)
// IMPORTANTE: Mantida a implementação ineficiente (O(N^2)) com loops aninhados
// para demonstrar a lentidão comparada ao T2.
void multMatVet(real_t **mat, real_t *vetIn, real_t *vetOut, int n) {
  for (int i = 0; i < n; ++i) {
    vetOut[i] = 0.0;
    for (int j = 0; j < n; ++j) { // Itera todas colunas (muitos zeros!)
      vetOut[i] += mat[i][j] * vetIn[j];
    }
  }
}

void copiaVetor(real_t *origem, real_t *destino, int n) {
  memcpy(destino, origem, n * sizeof(real_t));
}

void validaParam (int n, int k) {
  if ( n <= 10 || k <= 1 || k % 2 == 0) exit(EXIT_FAILURE);
} 

// Cria Matriz K-Diagonal (V1 - Alocação Densa Linha a Linha)
void criaKDiagonal(int n, int k, real_t ***A, real_t **B) {
  validaParam(n, k);

  *A = malloc(sizeof(real_t *) * n);
  for (int i = 0; i < n; i++) (*A)[i] = malloc(sizeof(real_t) * n);
  *B = malloc(sizeof(real_t) * n);

  int range = k / 2 + 1;

  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      if ( abs(i - j) < range ) (*A)[i][j] = generateRandomA(i, j, k);
      else (*A)[i][j] = 0.0;
    }

  for (int i = 0; i < n; i++) (*B)[i] = generateRandomB(k);
}

void destroiKDiagonal (int n, real_t **A, real_t *B) {
  for (int i = 0; i < n; i++) free(A[i]);
  free(A);
  free(B);
} 

// Geração Simétrica Positiva (V1 - MODIFICADA APENAS AQUI)
// Motivo: Permitir rodar N=20000 em tempo hábil.
// Mudança: Troca loop de 0..N por loop na banda (i-range..i+range).
void genSimetricaPositiva(real_t **A, real_t *b, int n, real_t ***ASP, real_t **bsp, rtime_t *tempo) {
  *tempo = timestamp();

  *ASP = malloc(sizeof(real_t *) * n);
  for (int i = 0; i < n; i++) (*ASP)[i] = calloc(n, sizeof(real_t));
  *bsp = calloc(n, sizeof(real_t));

  real_t **AT = malloc(sizeof(real_t *) * n);
  for (int i = 0; i < n; ++i) AT[i] = malloc(sizeof(real_t) * n);
  
  // Otimização Estrutural: Zera e preenche apenas a banda
  for(int i=0; i<n; i++) memset(AT[i], 0, n*sizeof(real_t));
  
  int range = 10; // Margem de segurança para k=7

  // Transposta Otimizada
  for (int i = 0; i < n; ++i) {
    int start = (i - range < 0) ? 0 : i - range;
    int end = (i + range > n) ? n : i + range;
    for (int j = start; j < end; ++j) AT[j][i] = A[i][j];
  }

  // Multiplicação Otimizada (Banda x Banda)
  for (int i = 0; i < n; ++i) {
    int start_j = (i - 2*range < 0) ? 0 : i - 2*range;
    int end_j = (i + 2*range > n) ? n : i + 2*range;

    for (int j = start_j; j < end_j; ++j) {
      real_t soma = 0.0;
      int start_l = (i < j) ? i - range : j - range;
      int end_l   = (i < j) ? j + range : i + range;
      if (start_l < 0) start_l = 0;
      if (end_l > n) end_l = n;

      for (int l = start_l; l < end_l; ++l) soma += AT[i][l] * A[l][j];
      (*ASP)[i][j] = soma;
    }
  }

  // Vetor BSP otimizado
  for (int i = 0; i < n; ++i) {
    int start = (i - range < 0) ? 0 : i - range;
    int end = (i + range > n) ? n : i + range;
    for (int j = start; j < end; ++j) (*bsp)[i] += AT[i][j] * b[j];
  }

  for (int i = 0; i < n; ++i) free(AT[i]);
  free(AT);

  *tempo = timestamp() - *tempo;
}

// Estruturas auxiliares (V1)
void geraDLU(real_t **A, int n, DLU_matrices_t **matrices, rtime_t *tempo) {
  *tempo = timestamp();
  *matrices = malloc(sizeof(DLU_matrices_t));
  (*matrices)->D = malloc(sizeof(real_t *) * n);
  (*matrices)->L = malloc(sizeof(real_t *) * n);
  (*matrices)->U = malloc(sizeof(real_t *) * n);

  for (int i = 0; i < n; i++) {
    (*matrices)->D[i] = calloc(n, sizeof(real_t));
    (*matrices)->L[i] = calloc(n, sizeof(real_t));
    (*matrices)->U[i] = calloc(n, sizeof(real_t));
    
    // Copia apenas a diagonal para D, inferior para L, superior para U
    // Implementação naive
    (*matrices)->D[i][i] = A[i][i];
    for(int j=0; j<i; j++) (*matrices)->L[i][j] = A[i][j];
    for(int j=i+1; j<n; j++) (*matrices)->U[i][j] = A[i][j];
  }
  *tempo = timestamp() - *tempo;
}

void destroiDLU(DLU_matrices_t *matrices, int n) {
    if (!matrices) return;
    for (int i = 0; i < n; i++) {
        free(matrices->D[i]); free(matrices->L[i]); free(matrices->U[i]);
    }
    free(matrices->D); free(matrices->L); free(matrices->U);
    free(matrices);
}

// Aplica Pré-condicionador (V1 - Naive)
void aplicaPreCond(real_t **A, real_t *r, real_t *z, int n, real_t omega, DLU_matrices_t *dlu) {
  if (omega == -1.0) copiaVetor(r, z, n);
  else if (omega == 0.0) { // Jacobi
    for (int i = 0; i < n; ++i) {
      if (A[i][i] != 0.0) z[i] = r[i] / A[i][i]; // Acesso indireto A[i][i]
    }
  }
  else if (omega >= 1.0) {
      // Implementação GS/SSOR simplificada (mantida do original se houver)
      copiaVetor(r, z, n); 
  }
}

real_t calcResiduoSL(real_t **A, real_t *b, real_t *X, int n, rtime_t *tempo) {
  *tempo = timestamp();
  real_t *Ax = calloc(n, sizeof(real_t));
  
  // Chama a multiplicação ineficiente
  multMatVet(A, X, Ax, n);

  real_t norma_L2 = 0.0;
  for (int i = 0; i < n; ++i) {
    real_t res = b[i] - Ax[i];
    norma_L2 += res * res; 
  }
  free(Ax);
  *tempo = timestamp() - *tempo;
  return sqrt(norma_L2);
}

// Solver GC (V1)
real_t resolveSL(real_t **A, real_t *b, real_t **X, int n, int maxit, real_t epsilon, real_t omega, rtime_t *tempo_iter, DLU_matrices_t *dlu, int *out_total_iterations) {
    *X = calloc(n, sizeof(real_t));
    real_t *r = malloc(sizeof(real_t) * n);
    real_t *p = malloc(sizeof(real_t) * n);
    real_t *z = malloc(sizeof(real_t) * n);
    real_t *Ap = malloc(sizeof(real_t) * n);
    real_t *x_velho = malloc(sizeof(real_t) * n);

    rtime_t t_inicio = timestamp();
    
    copiaVetor(b, r, n);
    aplicaPreCond(A, r, z, n, omega, dlu);
    copiaVetor(z, p, n);
    
    real_t r_z_velho = produtoEscalar(r, z, n);
    real_t norma_erro_final = 0.0;
    int k_final = 0;

    for (int k = 0; k < maxit; ++k) {
        k_final = k + 1;
        if (fabs(r_z_velho) < 1e-15) break; 

        copiaVetor(*X, x_velho, n);
        multMatVet(A, p, Ap, n); // Multiplicação densa

        real_t p_Ap = produtoEscalar(p, Ap, n);
        if (fabs(p_Ap) < 1e-15) break;
        real_t alpha = r_z_velho / p_Ap;

        for (int i = 0; i < n; ++i) {
            (*X)[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        norma_erro_final = 0.0;
        for (int i = 0; i < n; ++i) {
            real_t erro_i = fabs((*X)[i] - x_velho[i]);
            if (erro_i > norma_erro_final) norma_erro_final = erro_i;
        }
        if (norma_erro_final < epsilon) break;
        
        aplicaPreCond(A, r, z, n, omega, dlu);
        real_t r_z_novo = produtoEscalar(r, z, n);
        real_t beta = r_z_novo / r_z_velho;
        
        for (int i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];
        r_z_velho = r_z_novo;
    }
    
    rtime_t t_fim = timestamp();
    if (k_final > 0) *tempo_iter = (t_fim - t_inicio) / k_final;
    else *tempo_iter = 0.0;
    
    if (out_total_iterations) *out_total_iterations = k_final;

    free(r); free(p); free(z); free(Ap); free(x_velho);
    return norma_erro_final; 
}