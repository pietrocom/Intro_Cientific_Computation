// pc24 - GRR20241955

#ifndef __EQDIFF_H__
#define __EQDIFF_H__

typedef double real_t;

// Define o formato de impressão para os valores de ponto flutuante
#define FORMAT " %22.15e"

typedef struct {
  real_t *D, *Di, *Ds, *B; // D: Diagonal Principal, Di: Diagonal Inferior, Ds: Diagonal Superior, B: Termos Independentes
  int n;                   // Ordem do sistema
} Tridiag;

// Estrutura para a Equação Diferencial Ordinária
typedef struct {
  int n;                       // número de pontos internos na malha
  real_t a, b;                 // intervalo [a, b]
  real_t ya, yb;               // condições de contorno y(a) e y(b)
  real_t p, q, r1, r2, r3, r4; // coeficientes da EDO genérica e da função r(x)
} EDo;


// ---- Funções ----

// Gera a estrutura de um sistema linear tridiagonal a partir de uma EDO.
Tridiag *genTridiag (EDo *edoeq);

// Imprime a matriz aumentada do sistema linear na saída padrão.
void prnEDOsl (EDo *edoeq);

// Soluciona um sistema linear tridiagonal pelo método de Gauss-Seidel.
int gaussSeidel(Tridiag *SL, real_t *Y, real_t *residue_norm);

#endif 