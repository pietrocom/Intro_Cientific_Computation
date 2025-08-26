#ifndef __UTILS_H__
#define __UTILS_H__

#include <stdlib.h>
#include <time.h>

typedef double real_t;
typedef double rtime_t;
typedef char * string_t;

double timestamp(void);

long long ulp_distance(double a, double b);

#endif // __UTILS_H__

