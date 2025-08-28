#include <stdlib.h>
#include <stdio.h>
#include <likwid.h>

#include "eliminacaoGauss.h"
#include "gaussSeidel.h"
#include "utils.h"
#include "sislin.h"



int main () {

    SistLinear_t *C = lerSisLin();

    triangulariza(C);

    prnSisLin(C);

    return 0;
    /*
    real_t *X = malloc(sizeof(real_t) * C->n);

    retrosubst(C, X);

    free(X);
    liberaSisLin(C);
    */
}