#ifndef __EQDIFF_H__
#define __EQDIFF_H__

typedef double real_t;

#define FORMAT "%23.15e"

// Sistema linear Tri-diagonal
typedef struct {
  real_t *D, *Di, *Ds, *B;
  int n;
} Tridiag;

// Equação Diferencial Ordinária
typedef struct {
 int n; // número de pontos internos na malha
 real_t a, b; // intervalo
 real_t ya, yb; // condições contorno
 real_t (* p)(real_t), (* q)(real_t), (* r)(real_t);
} EDo;

#endif // __EQDIFF_H__

