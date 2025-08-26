#ifndef __SISLIN_H__
#define __SISLIN_H__

// Parâmetros default para teste de convergência
#define MAXIT 50
#define TOL  1.0e-4

// Estrutura para definiçao de um sistema linear qualquer
typedef struct {
  unsigned int n; // tamanho do SL
  real_t **A; // coeficientes
  real_t *b; // termos independentes
} SistLinear_t;

// Alocaçao e desalocação de matrizes
SistLinear_t* alocaSisLin (unsigned int n);
void liberaSisLin (SistLinear_t *SL);

// Leitura e impressão de sistemas lineares
SistLinear_t *lerSisLin ();
SistLinear_t *dupSisLin (SistLinear_t *src);
void prnSisLin (SistLinear_t *SL);
void prnVetor (real_t *vet, unsigned int n);

// Funções auxiliares
real_t normaMax(real_t *X1, real_t *X0, int n);
real_t normaL2(real_t *X, int n);
void residuo(SistLinear_t *SL, real_t *X, real_t *R, int n);

#endif // __SISLIN_H__

