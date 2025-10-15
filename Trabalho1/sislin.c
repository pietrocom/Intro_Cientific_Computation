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



void geraDLU (real_t *A, int n, int k, real_t **D, real_t **L, real_t **U, rtime_t *tempo)
{
  *tempo = timestamp();


  *tempo = timestamp() - *tempo;
}

/**
 * Devolve matriz M⁻¹
 *
 */
void geraPreCond(real_t *D, real_t *L, real_t *U, real_t w, int n, int k, real_t **M, rtime_t *tempo)
{
  *tempo = timestamp();


  *tempo = timestamp() - *tempo;
}


real_t calcResiduoSL(real_t **A, real_t *b, real_t *X, int n, rtime_t *tempo) {
  if (!A || !b || !X) {
    fprintf(stderr, "ERRO: Ponteiros inválidos!");
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


