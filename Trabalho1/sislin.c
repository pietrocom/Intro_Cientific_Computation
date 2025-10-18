#include <stdio.h>
#include <stdlib.h>    /* for exit e random/srandom */
#include <string.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"

static inline real_t generateRandomA( unsigned int i, unsigned int j, unsigned int k );
static inline real_t generateRandomB( unsigned int k );

/**
 * Função que gera os coeficientes de um sistema linear k-diagonal
 * @param i,j coordenadas do elemento a ser calculado (0<=i,j<n)
 * @param k numero de diagonais da matriz A
 */
static inline real_t generateRandomA( unsigned int i, unsigned int j, unsigned int k )
{
  static real_t invRandMax = 1.0 / (real_t)RAND_MAX;
  return ( (i==j) ? (real_t)(k<<1) : 1.0 )  * (real_t)random() * invRandMax;
}

/**
 * Função que gera os termos independentes de um sistema linear k-diagonal
 * @param k numero de diagonais da matriz A
 */
static inline real_t generateRandomB( unsigned int k )
{
  static real_t invRandMax = 1.0 / (real_t)RAND_MAX;
  return (real_t)(k<<2) * (real_t)random() * invRandMax;
}

// FUNÇÕES AUXILIARES
// Calcula o produto escalar entre dois vetores
real_t produtoEscalar(real_t *v1, real_t *v2, int n) {
  real_t resultado = 0.0;
  for (int i = 0; i < n; ++i) {
    resultado += v1[i] * v2[i];
  }
  return resultado;
}

// Multiplica uma matriz (n x n) por um vetor
void multMatVet(real_t **mat, real_t *vetIn, real_t *vetOut, int n) {
  for (int i = 0; i < n; ++i) {
    vetOut[i] = 0.0;
    for (int j = 0; j < n; ++j) {
      vetOut[i] += mat[i][j] * vetIn[j];
    }
  }
}

// Copia o conteúdo de um vetor para outro
void copiaVetor(real_t *origem, real_t *destino, int n) {
  for (int i = 0; i < n; ++i) {
    destino[i] = origem[i];
  }
}

// Verifica se os parâmetros contêm valores válidos
void validaParam (int n, int k) {
  if ( n <= 10 ) {
    fprintf(stderr, "ERRO: A entrada n eh menor que 10!\n");
    exit(EXIT_FAILURE);
  }

  if ( k <= 1 ) {
    fprintf(stderr, "ERRO: O número de diagonais k deve ser maior que 1!\n");
    exit(EXIT_FAILURE);
  }

  if (k % 2 == 0) {
    fprintf(stderr, "ERRO: O número de diagonais 'k' deve ser ímpar!\n");
    exit(EXIT_FAILURE);
  }
} 

/* 
Cria matriz 'A' k-diagonal e Termos independentes B 
OBS: Para este trabalho específico, a matriz tem dimensão n x n, ou seja, as posições
     fora da diagonal principal contém valor "0". Esta foi uma decisão proposital dado
     o objetivo do trabalho que é comparar dois códigos que fazem a mesma coisa, sendo
     um otimizado e outro não.
     Esta estrutura específica por si só talvez não tenha tanto impacto na performance
     final, sendo mais um problema relacionado com memória. Futuros testes poderão 
     esclarecer esta particularidade.
*/
void criaKDiagonal(int n, int k, real_t ***A, real_t **B) {
  validaParam(n, k);

  *A = malloc(sizeof(real_t *) * n);
  if (!(*A)) {
    fprintf(stderr, "ERRO: Durante alocacao de memoria!\n");
    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < n; i++) {
    (*A)[i] = malloc(sizeof(real_t) * n);
    if (!(*A)[i]) {
      fprintf(stderr, "ERRO: Durante alocacao de memoria!\n");
      exit(EXIT_FAILURE);
    }
  }

  (*B) = malloc(sizeof(real_t) * n);
  if (!(*B)) {
    fprintf(stderr, "ERRO: Durante alocacao de memoria!\n");
    exit(EXIT_FAILURE);
  }

  // Distância da diagonal
  int range = k / 2 + 1;

  // Preenchimento da matriz A
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      if ( abs(i - j) < range )
        (*A)[i][j] = generateRandomA(i, j, k);
      else {
        (*A)[i][j] = 0.0;
      }
    }

  // Preenchimento do vetor B
  for (int i = 0; i < n; i++)
    (*B)[i] = generateRandomB(k);
}


void destroiKDiagonal (int n, real_t **A, real_t *B) {
  if (!A || !B) {
    fprintf(stderr, "ERRO: Durante a desalocacao de memoria!");
    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < n; i++) {
    if (!A[i])
      continue;
    free(A[i]);
  }

  free(A);
  free(B);
} 

/* 
Gera matriz simetrica positiva 
OBS: A função implementada aqui é extremamente custosa, podendo ser otimizada caso
     a implementação da estrutura da matriz k-diagonal também seja ótima.
     Como explicitado anteriormente, a implementação desse trabalho não tem como 
     objetivo ser a mais eficiente, cabendo ao trabalho 2 "apertar os parafusos".
     Com isso em mente, a assinatura original da função foi modificada, pois nessa
     implementação o k não se faz necessário.
*/
void genSimetricaPositiva(real_t **A, real_t *b, int n, real_t ***ASP, real_t **bsp, rtime_t *tempo) {
  *tempo = timestamp();

  // Aloca o ponteiro para as linhas de ASP
  *ASP = malloc(sizeof(real_t *) * n);
  if (!(*ASP)) {
    fprintf(stderr, "ERRO: Falha ao alocar as linhas de ASP");
    exit(EXIT_FAILURE);
  }
  // Aloca cada linha de ASP
  for (int i = 0; i < n; i++) {
    (*ASP)[i] = malloc(sizeof(real_t) * n);
    if (!(*ASP)[i]) {
      fprintf(stderr, "ERRO: Falha ao alocar uma coluna de ASP");
      exit(EXIT_FAILURE);
    }
  }

  // Aloca o vetor bsp
  *bsp = malloc(sizeof(real_t) * n);
  if (!(*bsp)) {
    fprintf(stderr, "ERRO: Falha ao alocar bsp");
    exit(EXIT_FAILURE);
  }

  // Aloca e calcula a Transposta (AT)
  real_t **AT = malloc(sizeof(real_t *) * n);
  for (int i = 0; i < n; ++i) 
    AT[i] = malloc(sizeof(real_t) * n);
  
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      AT[j][i] = A[i][j];
    }
  }

  // Calcula ASP = AT * A (usando *ASP)
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      (*ASP)[i][j] = 0.0;
      for (int l = 0; l < n; ++l) {
        (*ASP)[i][j] += AT[i][l] * A[l][j];
      }
    }
  }

  // Calcula bsp = AT * b (usando *bsp)
  for (int i = 0; i < n; ++i) {
    (*bsp)[i] = 0.0;
    for (int j = 0; j < n; ++j) {
      (*bsp)[i] += AT[i][j] * b[j];
    }
  }

  // Libera memória temporária
  for (int i = 0; i < n; ++i) 
    free(AT[i]);
  free(AT);

  *tempo = timestamp() - *tempo;
}



void geraDLU(real_t **A, int n, DLU_matrices_t **matrices, rtime_t *tempo) {
  if (!A) {
    fprintf(stderr, "ERRO: Matriz de entrada A é nula!\n");
    exit(EXIT_FAILURE);
  }
    
  *tempo = timestamp();

  *matrices = malloc(sizeof(DLU_matrices_t));
  if (!(*matrices)) {
    perror("ERRO: Durante alocação da struct DLU!\n");
    exit(EXIT_FAILURE);
  }

  (*matrices)->D = malloc(sizeof(real_t *) * n);
  (*matrices)->L = malloc(sizeof(real_t *) * n);
  (*matrices)->U = malloc(sizeof(real_t *) * n);
  if (!(*matrices)->D || !(*matrices)->L || !(*matrices)->U) {
    perror("ERRO: Durante alocação da struct DLU!\n");
    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < n; i++) {
    (*matrices)->D[i] = malloc(sizeof(real_t) * n);
    (*matrices)->L[i] = malloc(sizeof(real_t) * n);
    (*matrices)->U[i] = malloc(sizeof(real_t) * n);
    if (!(*matrices)->D[i] || !(*matrices)->L[i] || !(*matrices)->U[i]) {
      perror("ERRO: Durante alocação da struct DLU!\n");
      exit(EXIT_FAILURE);
    }
  }

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i > j) { 
        (*matrices)->L[i][j] = A[i][j];
        (*matrices)->D[i][j] = 0.0;
        (*matrices)->U[i][j] = 0.0;
      } 
      else if (i == j) {
        (*matrices)->L[i][j] = 0.0;
        (*matrices)->D[i][j] = A[i][j];
        (*matrices)->U[i][j] = 0.0;
      }
      else if (i < j) {
        (*matrices)->L[i][j] = 0.0;
        (*matrices)->D[i][j] = 0.0;
        (*matrices)->U[i][j] = A[i][j];
      }
    }
  }

  *tempo = timestamp() - *tempo;
}

void destroiDLU(DLU_matrices_t *matrices, int n) {
    if (!matrices) return;

    for (int i = 0; i < n; i++) {
        free(matrices->D[i]);
        free(matrices->L[i]);
        free(matrices->U[i]);
    }
    free(matrices->D);
    free(matrices->L);
    free(matrices->U);
    free(matrices); // Libera a própria struct no final
}

/*
Houve uma mudança no nome e proposta desta função, que aplica o pré-condicionador para resolver Mz = r.
O condicionador varia a depender do omega.
As matrizes D, L e U podem ser nulas, exceto se for utilizado SSOR.
*/
void aplicaPreCond(real_t **A, real_t *r, real_t *z, int n, real_t omega, DLU_matrices_t *dlu) {
  // Caso 1: Sem pré-condicionador (M = I), então z = r.
  if (omega == -1.0) {
    copiaVetor(r, z, n);
  } 
  // Caso 2: Pré-condicionador de Jacobi (M = D), então zi = ri / Dii.
  else if (omega == 0.0) {
    for (int i = 0; i < n; ++i) {
      if (A[i][i] != 0.0) {
        z[i] = r[i] / A[i][i];
      } 
      else {
        // Divisão por zero
        fprintf(stderr, "ERRO: Divisão por zero no pré-condicionador de Jacobi na posição %d.\n", i);
        exit(EXIT_FAILURE);
      }
    }
  }
  // Caso 3 e 4: Gauss-Seidel (com omega = 1.0) e SSOR (1.0 < omega < 2.0)
  else if (omega >= 1.0 && omega < 2.0) {
    if (!dlu) {
      fprintf(stderr, "ERRO: Matrizes DLU são necessárias para SSOR, mas são nulas.\n");
      exit(EXIT_FAILURE);
    }

    // Aloca um vetor temporário
    real_t *y = malloc(sizeof(real_t) * n);
    if (!y) {
      fprintf(stderr, "ERRO: Durante alocação de vetor temporário!\n");
      exit(EXIT_FAILURE);
    }

    // Substituição para frente
    // Resolve (D/omega + L)y = r
    for (int i = 0; i < n; ++i) {
      real_t soma = 0.0; // Cuida do termo do somatório
      for (int j = 0; j < i; ++j) {
        soma += dlu->L[i][j] * y[j];
      }
      if (dlu->D[i][i] != 0.0) {
        // Fórmula original
        y[i] = (r[i] - omega * soma) / dlu->D[i][i];
      } 
      else {
        fprintf(stderr, "ERRO: Divisão por zero na substituição para frente do SSOR.\n");
        exit(EXIT_FAILURE);
      }
    }

    // Substituição para trás
    // Resolve (D/omega + U)z = (D/omega)*y
    for (int i = n - 1; i >= 0; --i) {
      real_t soma = 0.0;
      for (int j = i + 1; j < n; ++j) {
        soma += dlu->U[i][j] * z[j];
      }
      if (dlu->D[i][i] != 0.0) {
        // Fórmula algébrica simplificada
        z[i] = y[i] - (omega * soma) / dlu->D[i][i];
      } 
      else {
        fprintf(stderr, "ERRO: Divisão por zero na substituição para trás do SSOR.\n");
        exit(EXIT_FAILURE);
      }
    }

    free(y);
  }
  // Omega inválido
  else {
    fprintf(stderr, "ERRO: Valor de omega (%.2f) inválido.\n", omega);
    exit(EXIT_FAILURE);
  }
}


real_t calcResiduoSL(real_t **A, real_t *b, real_t *X, int n, rtime_t *tempo) {
  if (!A || !b || !X) {
    fprintf(stderr, "ERRO: Ponteiros inválidos!");
    exit(EXIT_FAILURE);
  }

  *tempo = timestamp();

  real_t *Ax = calloc(n, sizeof(real_t));
  if (!Ax) {
    perror("ERRO: Durante alocação de memória!");
    exit(EXIT_FAILURE);
  }
  
  multMatVet(A, X, Ax, n);

  // Calcula a norma L2 (b - Ax)
  real_t norma_L2_quadrada = 0.0;
  for (int i = 0; i < n; ++i) {
    real_t residuo_i = b[i] - Ax[i];
    norma_L2_quadrada += residuo_i * residuo_i; 
  }

  free(Ax);

  *tempo = timestamp() - *tempo;

  // Raiz quadrada da soma dos quadrados
  return sqrt(norma_L2_quadrada);
}

