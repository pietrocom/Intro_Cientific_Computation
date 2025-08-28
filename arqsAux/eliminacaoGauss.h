#ifndef __ELIM_GAUS__
#define __ELIM_GAUS__

#include "utils.h"
#include "sislin.h"

void triangulariza( SistLinear_t *C );
void retrosubst( SistLinear_t *C, real_t *X );

#endif 