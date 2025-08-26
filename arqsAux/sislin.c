#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "utils.h"
#include "sislin.h"

// Alocaçao de matriz em memória. 
SistLinear_t* alocaSisLin (unsigned int n)
{
  SistLinear_t *SL = (SistLinear_t *) malloc(sizeof(SistLinear_t));
  
  if ( SL ) {
    
    SL->n = n;
    SL->A = (real_t **) calloc (n, sizeof(real_t *));
    SL->A[0] = (real_t *) calloc (n*n, sizeof(real_t));
    
    for (int i=1; i < n; ++i)
      SL->A[i] = SL->A[i-1]+n;
	
    SL->b = (real_t *) calloc (n, sizeof(real_t));

    if (!(SL->A) || !(SL->b)) {
      liberaSisLin(SL);
      return NULL;
    }
  }
  
  return (SL);
}



// Liberacao de memória
void liberaSisLin (SistLinear_t *SL)
{
  if (SL) {
    if (SL->A) {
      if (SL->A[0]) free (SL->A[0]);
      free (SL->A);
    }
    if (SL->b) free(SL->b);
    free(SL);
  }
}


SistLinear_t *lerSisLin ()
{
  unsigned int n;
  SistLinear_t *SL;
  
  scanf("%d",&n);

  SL = alocaSisLin (n);
  
  for(int i=0; i < n; ++i) {
    for(int j=0; j < n; ++j)
      scanf ("%lg", &SL->A[i][j]);
    scanf ("%lg", &SL->b[i]);
  }

#ifdef __DEBUG__
  printf ("\n\n");
  prnSisLin(SL);
#endif /* __DEBUG__ */
  
  return SL;
}


SistLinear_t *dupSisLin (SistLinear_t *src)
{
  int n = src->n;
  
  SistLinear_t *dest = alocaSisLin (n);

  if (dest) {
    for(int i=0; i < n; ++i) {
      for(int j=0; j < n; ++j)
	dest->A[i][j] = src->A[i][j];
      dest->b[i] = src->b[i];
    }
  }
  
  return dest;
}


void prnSisLin (SistLinear_t *SL)
{
  int n=SL->n;

  for(int i=0; i < n; ++i) {
    printf("\n ");
    
    for(int j=0; j < n; ++j)
      printf ("%10.12lf ", SL->A[i][j]);

    printf ("   |   %.12lf", SL->b[i]);
  }
  printf("\n\n");
}

void prnVetor (real_t *v, unsigned int n)
{
  for(int i=0; i < n; ++i)
      printf ("%.12lf ", v[i]);
  printf ("\n");

}


// Calcula a norma máxima do erro absoluto aproximado de 2 vetores
real_t normaMax(real_t *X1, real_t *X0, int n)
{

  real_t k, norma = ABS(X1[0] - X0[0]);

  for (int i=1; i < n; ++i)
    if ((k = ABS(X1[i] - X0[i])) > norma)
      norma = k;

  return norma;
}

// Calcula a norma euclidiana de um vetor
real_t normaL2(real_t *X, int n)
{

  real_t norma = 0.0;

  for (int i=0; i < n; ++i)
    norma += X[i]*X[i];

  return sqrt(norma);
}


/* Calcula o resíduo R de um sistema AX = B
   Vetor *R já deve ter sido alocado previamente.
*/
void residuo(SistLinear_t *SL, real_t *X, real_t *R, int n)
{
  for(int i=0; i < n; ++i) {
    R[i] = SL->b[i];
    for(int j=0; j < n; ++j)
      R[i] -= SL->A[i][j]*X[j];
  }
}

