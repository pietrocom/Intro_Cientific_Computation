#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "utils.h"
#include "DoubleType.h"

/*  Retorna tempo em milisegundos

    Forma de uso:
 
    double tempo;
    tempo = timestamp();
    <trecho de programa do qual se deseja medir tempo>
    tempo = timestamp() - tempo;
*/

double timestamp(void)
{
  struct timespec tp;
  clock_gettime(CLOCK_MONOTONIC_RAW, &tp);
  return((double)(tp.tv_sec*1.0e3 + tp.tv_nsec*1.0e-6));
}

long long ulp_distance(double a, double b) {
    // Garante que a e b não são idênticos 
    if (a == b) {
        return 0;
    }

    Double_t u_a, u_b;
    u_a.f = a;
    u_b.f = b;

    if (u_a.i < 0) {
        u_a.i = INT64_MIN - u_a.i;
    }
    if (u_b.i < 0) {
        u_b.i = INT64_MIN - u_b.i;
    }

    return llabs(u_a.i - u_b.i);
}
