// pc24 - GRR20241955
// Implementação das funções para manipulação e solução
// da Equação Diferencial Ordinária (EDO).

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "utils.h"
#include "edo.h"

#define MAX_ITER 100
#define RESIDUE_TOL 1e-5

// Gera a estrutura do sistema linear tridiagonal a partir dos parâmetros da EDO.
Tridiag *genTridiag (EDo *edo)
{
  Tridiag *sl;
  real_t x, rx;
  int n = edo->n;
  
  sl = (Tridiag *) malloc (sizeof(Tridiag));
  if (!sl) return NULL;
  sl->n = edo->n;

  sl->D = (real_t *) calloc(n, sizeof(real_t));
  sl->Di = (real_t *) calloc(n, sizeof(real_t));
  sl->Ds = (real_t *) calloc(n, sizeof(real_t));
  sl->B = (real_t *) calloc(n, sizeof(real_t));
  
  if (!sl->D || !sl->Di || !sl->Ds || !sl->B) return NULL;

  real_t h = (edo->b - edo->a)/(n+1);

  for (int i=0; i < n; ++i) {
    x = edo->a + (i+1)*h;
    rx = edo->r1*x + edo->r2*x*x + edo->r3*cos(x) + edo->r4*exp(x);
    
    sl->B[i] = h*h * rx;
    sl->Di[i] = 1.0 - h * edo->p/2.0;
    sl->D[i] = -2.0 + h*h * edo->q;
    sl->Ds[i] = 1.0 + h * edo->p/2.0;
  }

  sl->B[0] -= edo->ya * (1.0 - h*edo->p/2.0);
  sl->B[n-1] -= edo->yb * (1.0 + h*edo->p/2.0);
  
  return sl;
}

// Imprime a matriz aumentada do sistema linear na saída padrão.
void prnEDOsl (EDo *edoeq)
{
  int n = edoeq->n, i, j;
  real_t x, b, d, di, ds,rx;
  real_t h = (edoeq->b - edoeq->a)/(n+1);

  printf ("%d\n", n);

  for (i=0; i < n; ++i) {
    x = edoeq->a + (i+1)*h;
    rx = edoeq->r1*x + edoeq->r2*x*x + edoeq->r3*cos(x) + edoeq->r4*exp(x);
    
    b = h*h * rx; 
    di = 1.0 - h * edoeq->p/2.0;
    d = -2.0 + h*h * edoeq->q;
    ds = 1.0 + h * edoeq->p/2.0;
      
    for (j=0; j < n; ++j) {
      if (i == j)
        printf (FORMAT,d);
      else if (j == i-1)
        printf (FORMAT,di);
      else if (j == i+1)
        printf (FORMAT,ds);
      else
        printf(FORMAT, 0.0);
    }
      
    if (i == 0)
      b -= edoeq->ya * (1.0 - h*edoeq->p/2.0);
    else if (i == n-1)
      b -= edoeq->yb * (1.0 + h*edoeq->p/2.0);

    printf (FORMAT, b);
      
    printf ("\n");
  }
}

// Soluciona um sistema linear tridiagonal pelo método de Gauss-Seidel.
// Retorna o número de iterações executadas.
int gaussSeidel(Tridiag *SL, real_t *Y, real_t *residue_norm) {
    int n = SL->n;
    
    // Inicializa o vetor solução Y com zeros
    for (int i = 0; i < n; ++i) {
        Y[i] = 0.0;
    }

    int iter;
    for (iter = 1; iter <= MAX_ITER; ++iter) {
        
        // Laço de atualização de Gauss-Seidel
        for (int i = 0; i < n; ++i) {
            real_t sum = SL->B[i];
            if (i > 0) {
                sum -= SL->Di[i] * Y[i-1];
            }
            if (i < n - 1) {
                sum -= SL->Ds[i] * Y[i+1];
            }
            Y[i] = sum / SL->D[i];
        }
        
        // Calcula a norma L2 do resíduo ||b - Ax_k||
        real_t norm_sq = 0.0;
        for (int i = 0; i < n; ++i) {
            real_t Ax_i = SL->D[i] * Y[i];
            if (i > 0) Ax_i += SL->Di[i] * Y[i-1];
            if (i < n-1) Ax_i += SL->Ds[i] * Y[i+1];
            
            real_t residue_i = SL->B[i] - Ax_i;
            norm_sq += residue_i * residue_i;
        }
        *residue_norm = sqrt(norm_sq);

        // Verifica o critério de parada
        if (*residue_norm <= RESIDUE_TOL) {
            return iter;
        }
    }

    // Retorna o número máximo de iterações se não convergiu
    return MAX_ITER;
}