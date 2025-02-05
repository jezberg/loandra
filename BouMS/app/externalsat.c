#include "externalsat.h"

#include <BouMS/common.h>
#include <BouMS/wcnf.h>
#include <stdbool.h>
#include <stddef.h>

#if defined(BOUMS_SATSOLVER_CADICAL)
#  include <ccadical.h>
#  define SATSOLVER_PREFIX ccadical
typedef CCaDiCaL solver_t;
#endif

#define _CONCAT(a, b) a##_##b
#define CONCAT(a, b)  _CONCAT(a, b)
#define IPASIR(f)     CONCAT(SATSOLVER_PREFIX, f)

#define IPASIR_INIT          IPASIR(init)
#define IPASIR_SET_TERMINATE IPASIR(set_terminate)
#define IPASIR_ADD           IPASIR(add)
#define IPASIR_SOLVE         IPASIR(solve)
#define IPASIR_VAL           IPASIR(val)
#define IPASIR_RELEASE       IPASIR(release)

sat_result_t addClauses(solver_t* solver, const BouMS_wcnf_t* formula);

int terminateCb(void* stop);

sat_result_t preSolveHard(const BouMS_wcnf_t* formula, bool* model, bool* stop) {
  sat_result_t ret = UNKNOWN;

  solver_t* solver = IPASIR_INIT();
  if (solver) {
    IPASIR_SET_TERMINATE(solver, stop, terminateCb);

    ret = addClauses(solver, formula);

    if (ret == UNKNOWN) {
      ret = IPASIR_SOLVE(solver);

      if (ret != SAT && ret != UNSAT) {
        ret = UNKNOWN;
      } else if (ret == SAT) {
        for (BouMS_uint_t varIdx = 0; varIdx < formula->numVariables; ++varIdx) {
          const int var = (int)varIdx + 1;
          model[varIdx] = IPASIR_VAL(solver, var) > 0;
        }
      }
    }

    IPASIR_RELEASE(solver);
    solver = NULL;
  }

  return ret;
}

sat_result_t addClauses(solver_t* solver, const BouMS_wcnf_t* formula) {
  sat_result_t ret = UNKNOWN;

  for (BouMS_uint_t clauseIdx = 0; clauseIdx < formula->numClauses; ++clauseIdx) {
    const BouMS_wcnf_clause_t* clause = formula->clauses + clauseIdx;

    if (BouMS_wcnf_isClauseHard(clause)) {
      for (BouMS_uint_t litIdx = 0; litIdx < clause->numLiterals; ++litIdx) {
        const BouMS_wcnf_literal_t* lit = clause->literals + litIdx;
        const int var = (int)BouMS_wcnf_var(lit) + 1;
        IPASIR_ADD(solver, BouMS_wcnf_sign(lit) ? -var : var);
      }
      IPASIR_ADD(solver, 0);
    }
  }

  return ret;
}

int terminateCb(void* stop) {
  return (int)*((bool*)stop);
}
