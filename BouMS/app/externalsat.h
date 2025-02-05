#ifndef EXTERNALSAT_H
#define EXTERNALSAT_H

#include <BouMS/wcnf.h>
#include <stdbool.h>

typedef enum { UNKNOWN = 0, SAT = 10, UNSAT = 20 } sat_result_t;

sat_result_t preSolveHard(const BouMS_wcnf_t* formula, bool* model, bool* stop);

#endif
