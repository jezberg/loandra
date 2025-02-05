/**
 * @file common.h
 * @author Ole Lübke (ole.luebke@tuhh.de)
 * @brief Core solving functions
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef COMMON_H
#define COMMON_H

#include <stdbool.h>

#include "BouMS/BouMS.h"
#include "BouMS/common.h"
#include "BouMS/wcnf.h"
#include "private/fixedprec.h"
#include "private/types.h"

/**
 * @brief Initializes variables to current best known assignment, or randomly if no such assignment is known yet
 *
 * @param formula
 * @param result
 * @param initModel
 * @param forceRandom
 */
void initVars(const BouMS_wcnf_t* formula, const BouMS_result_t* result, const bool* initModel, bool forceRandom);

/**
 * @brief Checks the satisfiability of the formula and updates related fields in mem and cost
 *
 * @param formula
 * @param mem
 * @param cost
 */
void initAlgo(const BouMS_wcnf_t* formula, BouMS_memory_t* mem, BouMS_uint_t* cost);

/**
 * @brief Selects a variable to flip; also updates weights if no decreasing variable is found
 *
 * @param formula
 * @param cfg
 * @param mem
 */
BouMS_wcnf_variable_t* selectVariable(const BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem);

/**
 * @brief Flips the value of the given variable and updates satis/falsified status and cost accordingly
 *
 * @param variable
 * @param formula
 * @param mem
 * @param cost
 */
void flipVariable(BouMS_wcnf_variable_t* variable, const BouMS_wcnf_t* formula, BouMS_memory_t* mem, BouMS_uint_t* cost);

/**
 * @brief Updates the scores of variables after the weight of the given clause has been increased
 *
 * @param clauseIdx
 * @param clause
 * @param increase
 * @param mem
 */
void updateScores(BouMS_uint_t clauseIdx, const BouMS_wcnf_clause_t* clause, signed_fixedprec_t increase,
                  BouMS_memory_t* mem);

/**
 * @brief Adds a variable to the decreasing variables if it is not present in the list
 *
 * @param varIdx
 * @param mem
 */
static inline void addDecreasingVar(BouMS_uint_t varIdx, BouMS_memory_t* mem) {
  const score_t score = mem->scores[varIdx];
  if (mem->decreasingVarsIdx[varIdx] == BOUMS_UINT_MAX && score.value > 0) {
    // var wasn't decreasing but is now, append it to the end of the list
    mem->decreasingVars[mem->numDecreasingVars] = varIdx;
    mem->decreasingVarsIdx[varIdx] = mem->numDecreasingVars;
    mem->numDecreasingVars += 1;
  }
}

/**
 * @brief Removes a variable from the decreasing variables if it is present in the list
 *
 * @param varIdx
 * @param mem
 */
static inline void remDecreasingVar(BouMS_uint_t varIdx, BouMS_memory_t* mem) {
  const score_t score = mem->scores[varIdx];
  if (mem->decreasingVarsIdx[varIdx] < BOUMS_UINT_MAX && score.value <= 0) {
    // var was decreasing but isn't anymore
    const BouMS_uint_t replaceIdx = mem->decreasingVarsIdx[varIdx];

    // remove it from the list
    mem->numDecreasingVars -= 1;
    mem->decreasingVarsIdx[varIdx] = BOUMS_UINT_MAX;

    // if it wasn't the last list element, move the last element to fill the gap
    if (replaceIdx < mem->numDecreasingVars) {
      const BouMS_uint_t moveVar = mem->decreasingVars[mem->numDecreasingVars];
      mem->decreasingVars[replaceIdx] = moveVar;
      mem->decreasingVarsIdx[moveVar] = replaceIdx;
    }
  }
}

/**
 * @brief Adds a clause to the specified list of clauses
 *
 * @param clauseIdx
 * @param numClauses
 * @param clauses
 * @param clausesIdx
 */
static inline void addClause(BouMS_uint_t clauseIdx, BouMS_uint_t* numClauses, BouMS_uint_t* clauses,
                             BouMS_uint_t* clausesIdx) {
  // append to the end of the list
  clauses[*numClauses] = clauseIdx;
  clausesIdx[clauseIdx] = *numClauses;
  *numClauses += 1;
}

/**
 * @brief Removes a list from the specified list of clauses
 *
 * @param clauseIdx
 * @param numClauses
 * @param clauses
 * @param clausesIdx
 */
static inline void remClause(BouMS_uint_t clauseIdx, BouMS_uint_t* numClauses, BouMS_uint_t* clauses,
                             BouMS_uint_t* clausesIdx) {
  const BouMS_uint_t replaceIdx = clausesIdx[clauseIdx];

  // remove from the list
  *numClauses -= 1;
  clausesIdx[clauseIdx] = BOUMS_UINT_MAX;

  // if it wasn't the last element, move the last element to fill the gap
  if (replaceIdx < *numClauses) {
    const BouMS_uint_t moveClause = clauses[*numClauses];
    clauses[replaceIdx] = moveClause;
    clausesIdx[moveClause] = replaceIdx;
  }
}

/**
 * @brief Check whether there are any unsatisfied clauses left
 *
 * @param mem
 * @param result
 */
static inline bool done(const BouMS_memory_t* mem, const BouMS_result_t* result) {
  return result->status == BOUMS_OPTIMUM_FOUND || result->status == BOUMS_UNSAT ||
         (mem->numFalsifiedHardClauses == 0 && mem->numFalsifiedSoftClauses == 0);
}

/**
 * @brief Copy a variable assignment
 *
 * @param n The number of variables
 * @param input
 * @param output
 */
static inline void useAssignment(BouMS_uint_t n, const bool* input, BouMS_wcnf_variable_t* output) {
  for (BouMS_uint_t i = 0; i < n; ++i) {
    output[i].value = input[i];
  }
}

/**
 * @brief Copy a variable assignment
 *
 * @param n The number of variables
 * @param input
 * @param output
 */
static inline void saveAssignment(BouMS_uint_t n, const BouMS_wcnf_variable_t* input, bool* output) {
  for (BouMS_uint_t i = 0; i < n; ++i) {
    output[i] = input[i].value;
  }
}

/**
 * @brief Checks whether all hard clauses are satisfied by the current assignment
 * @note Hard clauses must be sorted to the front of the clauses array
 *
 * @param formula
 */
static inline bool isFeasible(const BouMS_wcnf_t* formula) {
  for (BouMS_uint_t clauseIdx = 0; clauseIdx < formula->numHardClauses; ++clauseIdx) {
    const BouMS_wcnf_clause_t* const clause = formula->clauses + clauseIdx;

    bool sat = false;
    for (BouMS_uint_t litIdx = 0; litIdx < clause->numLiterals; ++litIdx) {
      const BouMS_wcnf_literal_t* const literal = clause->literals + litIdx;

      if (BouMS_wcnf_isLiteralSatisfied(formula, literal)) {
        sat = true;
        break;
      }
    }

    if (!sat) {
      return false;
    }
  }
  return true;
}

/**
 * @brief Sets the clause pointers of all literals of the clause to point to this clause
 *
 * @param clause
 */
static inline void fixLitToClausePtrs(const BouMS_wcnf_clause_t* clause) {
  for (BouMS_uint_t litIdx = 0; litIdx < clause->numLiterals; ++litIdx) {
    clause->literals[litIdx].clause = clause;
  }
}

/**
 * @brief Removes a clause by moving it to the end of the array and replacing it with the clause that was there before
 *
 * @param clauseIdx
 * @param clauses
 * @param numClauses
 */
static inline void removeClause(BouMS_uint_t clauseIdx, BouMS_wcnf_clause_t* clauses, BouMS_uint_t numClauses) {
  const BouMS_uint_t lastClauseIdx = numClauses - 1;

  {
    const BouMS_wcnf_clause_t clauseToRemove = clauses[clauseIdx];
    clauses[clauseIdx] = clauses[lastClauseIdx];
    clauses[lastClauseIdx] = clauseToRemove;
  }

  fixLitToClausePtrs(clauses + clauseIdx);
  fixLitToClausePtrs(clauses + lastClauseIdx);
}

#endif
