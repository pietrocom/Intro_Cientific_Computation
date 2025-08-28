/* Matriz  'normal' (vetor  de  ponteiros (linhas  matriz) para  vetores
   (colunas da matriz), estilo 'Mazieiro/Prog 2'
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"
#include "eliminacaoGauss.h"

static int encontraMax(SistLinear_t *C, int k, int n)
{
  if (!C)
    return -1;

  real_t max = C->A[k][k];
  int linha = k;
  for (int i = k; i < n; i++) {
    real_t aux = C->A[i][k];
    if (fabs(aux) > fabs(max)) {
      max = aux;
      linha = i;
    }
  }

  return linha;
}

static void trocaLinha (SistLinear_t *C, int k, int p, int n)
{
  if (!C)
    return;

  for (int i = 0; i < n; i++) {
    real_t aux = C->A[k][i];
    C->A[k][i] = C->A[p][i];
    C->A[p][i] = aux;
  }

  real_t aux = C->b[k];
  C->b[k] = C->b[p];
  C->b[p] = aux;
}

/* Seja um S.L. de ordem 'n'
   C = A|B em Ax=B
 */
void triangulariza( SistLinear_t *C )
{
  unsigned int n = C->n;
  for (unsigned int i = 0; i < n - 1; i++) {
    int max = encontraMax(C, i, n);
    trocaLinha(C, i, max, (int)n);

    if (fabs(C->A[i][i]) <= __DBL_EPSILON__) {
      perror("Divisao por 0!\n");
      exit(1);
    }

    for (unsigned int j = i + 1; j < n; j++) {
      real_t m = C->A[j][i] / C->A[i][i];
      C->A[j][i] = 0;
      for (unsigned int k = i + 1; k < n; k++) {
        C->A[j][k] -= (C->A[i][k] * m); 
      }
      C->b[j] -= (C->b[i] * m);
    }
  }
}

// Assume que o sistema já é triangular
void retrosubst( SistLinear_t *C, real_t *X )
{
  int n = C->n;
  for (int i = n - 1; i >= 0; i--) {
    real_t soma = 0.0;
    for (int j = n - 1; j > i; j--) {
      soma += C->A[i][j] * X[j];
    }
    X[i] = (C->b[i] - soma) / C->A[i][i];
  }
}
