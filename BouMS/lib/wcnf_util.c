/**
 * @file wcnf_util.c
 * @author Ole Lübke (ole.luebke@tuhh.de)
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "BouMS/wcnf_util.h"

#include <stdlib.h>

#include "BouMS/common.h"
#include "BouMS/wcnf.h"
#include "private/common.h"

/**
 * @brief Ensures that the formula's clauses array is large enough to hold numClauses clauses
 *
 * @note Does not initialize the new clauses with default values
 *
 * @param formula
 * @param numClauses
 * @param realloc
 * @param stop
 * @return true In case of error or stop
 * @return false Otherwise
 */
static inline bool ensureClauses(BouMS_wcnf_t* formula, BouMS_uint_t numClauses, BouMS_wcnf_util_realloc_t realloc,
                                 const bool* stop);

/**
 * @brief Ensures that the formula's variable array is large enough to hold numVariables variables
 *
 * @note Initializes the new variables with default values
 *
 * @param formula
 * @param numVariables
 * @param realloc
 * @param stop
 * @return true In case of error or stop
 * @return false Otherwise
 */
static inline bool ensureVariables(BouMS_wcnf_t* formula, BouMS_uint_t numVariables, BouMS_wcnf_util_realloc_t realloc,
                                   const bool* stop);

/**
 * @brief Ensures that the variable's literal pointer array is large enough to hold numLitPtrs literal pointers
 *
 * @note Does not update the numLiterals field
 *
 * @param variable
 * @param numLitPtrs
 * @param realloc
 */
static inline bool ensureLitPtrs(BouMS_wcnf_variable_t* variable, BouMS_uint_t numLitPtrs,
                                 BouMS_wcnf_util_realloc_t realloc);

/** @brief Counts the number of literals for each variable in numClauses clauses, starting with firstClause
 *
 * @param formula
 * @param firstClause
 * @param numClauses
 * @param realloc
 * @param stop
 * @return Pointer to array with counts (must be freed by caller)
 * @return NULL in case of error or stop
 */
static inline BouMS_uint_t* countLitsPerVar(const BouMS_wcnf_t* formula, BouMS_uint_t firstClause,
                                            BouMS_uint_t numClauses, BouMS_wcnf_util_realloc_t realloc,
                                            BouMS_wcnf_util_free_t free, const bool* stop);

/**
 * @brief Add pointers from variables to their occurences in clauses for newNumClauses, starting at newClausesStart
 *
 * @param formula
 * @param newClausesStart
 * @param numNewClauses
 * @param realloc
 * @param free
 * @param stop
 * @return true In case of error or stop
 * @return false Otherwise
 */
static inline bool addVarToLitPtrs(BouMS_wcnf_t* formula, BouMS_uint_t newClausesStart, BouMS_uint_t numNewClauses,
                                   BouMS_wcnf_util_realloc_t realloc, BouMS_wcnf_util_free_t free, const bool* stop);

void BouMS_wcnf_util_deleteFormula(BouMS_wcnf_t* formula, BouMS_wcnf_util_free_t free, const bool* stop) {
  static const bool STOP_DUMMY = false;
  if (stop == NULL) {
    stop = &STOP_DUMMY;
  }

  if (formula) {
    if (formula->clauses) {
      for (BouMS_uint_t clauseIdx = 0; clauseIdx < formula->numClauses && !*stop; ++clauseIdx) {
        BouMS_wcnf_clause_t* clause = formula->clauses + clauseIdx;
        if (clause->literals) {
          free(clause->literals);
          clause->literals = NULL;
        }
      }

      free(formula->clauses);
      formula->clauses = NULL;
    }

    if (formula->variables) {
      for (BouMS_uint_t varIdx = 0; varIdx < formula->numVariables && !*stop; ++varIdx) {
        BouMS_wcnf_variable_t* variable = formula->variables + varIdx;
        if (variable->literals) {
          free((void*)variable->literals);
          variable->literals = NULL;
        }
      }

      free(formula->variables);
      formula->variables = NULL;
    }
  }
}

bool BouMS_wcnf_util_addClause(BouMS_wcnf_t* formula, BouMS_uint_t weight, const int* literals,
                               BouMS_uint_t numLiterals, BouMS_wcnf_util_realloc_t realloc,
                               BouMS_wcnf_util_free_t free) {
  static const bool STOP_DUMMY = false;

  if (ensureClauses(formula, formula->numClauses + 1, realloc, &STOP_DUMMY)) {
    return true;
  }
  formula->numClauses -= 1;  // only update count when everything went without error

  BouMS_wcnf_clause_t* const clause = formula->clauses + formula->numClauses;

  clause->weight = weight;

  clause->literals = realloc(NULL, numLiterals * sizeof(BouMS_wcnf_literal_t));
  if (!clause->literals) {
    return true;
  }
  clause->numLiterals = numLiterals;

  const BouMS_uint_t oldNumVars = formula->numVariables;

  BouMS_uint_t litIdx;
  for (litIdx = 0; litIdx < numLiterals; ++litIdx) {
    BouMS_wcnf_literal_t* const lit = clause->literals + litIdx;
    lit->clause = clause;
    lit->lit = BouMS_wcnf_util_mkLit(literals[litIdx]);

    const BouMS_uint_t varIdx = BouMS_wcnf_var(lit);
    if (ensureVariables(formula, varIdx + 1, realloc, &STOP_DUMMY)) {
      break;
    }

    BouMS_wcnf_variable_t* const var = formula->variables + varIdx;

    const BouMS_uint_t newNumLits = var->numLiterals + 1;
    if (ensureLitPtrs(var, newNumLits, realloc)) {
      break;
    }
    var->literals[var->numLiterals] = lit;
    var->numLiterals = newNumLits;
  }
  if (litIdx < numLiterals) {  // error, loop was broken
    formula->numVariables = oldNumVars;
    for (BouMS_uint_t undoIdx = 0; undoIdx < litIdx; ++undoIdx) {
      const BouMS_uint_t varIdx = BouMS_wcnf_var(clause->literals + undoIdx);
      if (varIdx < formula->numVariables) {
        formula->variables[varIdx].numLiterals -= 1;
      }
    }

    free(clause->literals);
    clause->literals = NULL;

    return true;
  }

  formula->numClauses += 1;  // cf. above
  if (BouMS_wcnf_isClauseHard(clause)) {
    formula->numHardClauses += 1;
  }

  return false;
}

bool BouMS_wcnf_util_startBatchClauseAdding(BouMS_uint_t expectedNumClauses, BouMS_wcnf_util_realloc_t realloc,
                                            BouMS_wcnf_util_free_t free,
                                            BouMS_wcnf_util_batchClauseAddingState_t* state) {
  *state = (BouMS_wcnf_util_batchClauseAddingState_t){
      .capacity = 0, .size = 0, .clauses = NULL, .maxVar = 0, .realloc = realloc, .free = free};

  if (expectedNumClauses == 0) {
    expectedNumClauses = 2;
  }

  state->clauses = realloc(NULL, expectedNumClauses * sizeof(BouMS_wcnf_clause_t));
  if (!state->clauses) {
    return true;
  }

  state->capacity = expectedNumClauses;

  return false;
}

bool BouMS_wcnf_util_finishBatchClauseAdding(BouMS_wcnf_util_batchClauseAddingState_t* state, BouMS_wcnf_t* formula,
                                             const bool* stop) {
  static const bool STOP_DUMMY = false;
  if (stop == NULL) {
    stop = &STOP_DUMMY;
  }

  if (state->size == 0) {
    state->free(state->clauses);
    state->clauses = NULL;
    state->capacity = 0;
    return false;
  }

  const BouMS_uint_t oldNumClauses = formula->numClauses;
  const BouMS_uint_t oldNumHardClauses = formula->numHardClauses;
  const BouMS_uint_t oldNumVars = formula->numVariables;

  // if formula is empty just reuse clause buffer as formula's clause array
  if (formula->numClauses == 0 && formula->clauses == NULL) {
    formula->numClauses = state->size;
    formula->clauses = state->clauses;
    if (state->capacity > state->size) {  // shrink clause list to size
      BouMS_wcnf_clause_t* const shrunk =
          state->realloc(formula->clauses, formula->numClauses * sizeof(BouMS_wcnf_clause_t));
      if (shrunk) {
        formula->clauses = shrunk;
      }
    }
  }
  // ensure formula's clause array is long enough (no-op if formula was previously empty)
  if (ensureClauses(formula, oldNumClauses + state->size, state->realloc, stop)) {
    formula->numClauses = oldNumClauses;
    return true;
  }
  // if the formula originally wasn't empty, the clauses need to be copied from the clause buffer
  const bool needCopy = formula->clauses != state->clauses;
  // for all new clauses ...
  for (BouMS_uint_t clauseBufIdx = 0; clauseBufIdx < state->size && !*stop; ++clauseBufIdx) {
    BouMS_wcnf_clause_t* const clause = state->clauses + clauseBufIdx;
    // ... copy them to the formula's clause array if necessary
    if (needCopy) {
      formula->clauses[formula->numClauses++] = *clause;
    }
    // .. count them if they're hard
    if (BouMS_wcnf_isClauseHard(clause)) {
      formula->numHardClauses += 1;
    }
  }

  if (*stop) {
    return true;
  }

  // ensure formula's variables array is long enough
  if (ensureVariables(formula, state->maxVar, state->realloc, stop)) {
    formula->numClauses = oldNumClauses;
    formula->numHardClauses = oldNumHardClauses;
    formula->numVariables = oldNumVars;
    return true;
  }

  // set variable occurence pointers
  if (addVarToLitPtrs(formula, oldNumClauses, state->size, state->realloc, state->free, stop)) {
    formula->numClauses = oldNumClauses;
    formula->numHardClauses = oldNumHardClauses;
    formula->numVariables = oldNumVars;
    return true;
  }

  // check if the clause buffer was used without copying above
  if (formula->clauses != state->clauses) {
    state->free(state->clauses);
  }

  // reset clause adding state
  state->clauses = NULL;
  state->capacity = 0;
  state->size = 0;

  return false;
}

void BouMS_wcnf_util_cleanUpBatchClauseAddingAfterError(BouMS_wcnf_util_batchClauseAddingState_t* state) {
  if (state->clauses) {
    state->free(state->clauses);
    state->clauses = NULL;
  }
  state->capacity = 0;
  state->size = 0;
}

static inline bool ensureClauses(BouMS_wcnf_t* formula, BouMS_uint_t numClauses, BouMS_wcnf_util_realloc_t realloc,
                                 const bool* stop) {
  if (numClauses > formula->numClauses) {
    BouMS_wcnf_clause_t* const newClauses = realloc(formula->clauses, numClauses * sizeof(BouMS_wcnf_clause_t));
    if (!newClauses) {
      return true;
    }
    if (formula->clauses != newClauses) {  // memory has been moved
      formula->clauses = newClauses;
      for (BouMS_uint_t clauseIdx = 0; clauseIdx < formula->numClauses && !*stop; ++clauseIdx) {
        fixLitToClausePtrs(formula->clauses + clauseIdx);
      }
    }
    formula->numClauses = numClauses;
  }

  return *stop;
}

static inline bool ensureVariables(BouMS_wcnf_t* formula, BouMS_uint_t numVariables, BouMS_wcnf_util_realloc_t realloc,
                                   const bool* stop) {
  if (numVariables > formula->numVariables) {
    BouMS_wcnf_variable_t* const newVars = realloc(formula->variables, numVariables * sizeof(BouMS_wcnf_variable_t));
    if (!newVars) {
      return true;
    }
    formula->variables = newVars;
    for (BouMS_uint_t newVarIdx = formula->numVariables; newVarIdx < numVariables && !*stop; ++newVarIdx) {
      BouMS_wcnf_variable_t* const var = formula->variables + newVarIdx;
      var->numLiterals = 0;
      var->literals = NULL;
      var->value = false;
    }
    formula->numVariables = numVariables;
  }
  return *stop;
}

static inline bool ensureLitPtrs(BouMS_wcnf_variable_t* variable, BouMS_uint_t numLitPtrs,
                                 BouMS_wcnf_util_realloc_t realloc) {
  if (numLitPtrs > variable->numLiterals) {
    const BouMS_wcnf_literal_t** const newLits =
        (const BouMS_wcnf_literal_t**)realloc((void*)variable->literals, numLitPtrs * sizeof(BouMS_wcnf_literal_t*));
    if (!newLits) {
      return true;
    }
    variable->literals = newLits;
  }
  return false;
}

static inline BouMS_uint_t* countLitsPerVar(const BouMS_wcnf_t* formula, BouMS_uint_t firstClause,
                                            BouMS_uint_t numClauses, BouMS_wcnf_util_realloc_t realloc,
                                            BouMS_wcnf_util_free_t free, const bool* stop) {
  BouMS_uint_t* newLitsPerVar = realloc(NULL, formula->numVariables * sizeof(BouMS_uint_t));
  if (!newLitsPerVar) {
    return NULL;
  }
  for (BouMS_uint_t varIdx = 0; varIdx < formula->numVariables && !*stop; ++varIdx) {
    newLitsPerVar[varIdx] = 0;
  }
  for (BouMS_uint_t clauseIdx = firstClause; clauseIdx < firstClause + numClauses && !*stop; ++clauseIdx) {
    const BouMS_wcnf_clause_t* const clause = formula->clauses + clauseIdx;
    for (BouMS_uint_t litIdx = 0; litIdx < clause->numLiterals; ++litIdx) {
      BouMS_wcnf_literal_t* const lit = clause->literals + litIdx;
      newLitsPerVar[BouMS_wcnf_var(lit)] += 1;
    }
  }

  if (*stop) {
    free(newLitsPerVar);
    newLitsPerVar = NULL;
  }

  return newLitsPerVar;
}

static inline bool addVarToLitPtrs(BouMS_wcnf_t* formula, BouMS_uint_t newClausesStart, BouMS_uint_t numNewClauses,
                                   BouMS_wcnf_util_realloc_t realloc, BouMS_wcnf_util_free_t free, const bool* stop) {
  // count the number of new literals per variable (to reduce need for reallocating)
  {
    BouMS_uint_t* newLitsPerVar = countLitsPerVar(formula, newClausesStart, numNewClauses, realloc, free, stop);
    if (!newLitsPerVar) {
      return true;
    }

    // allocate space for pointers from variables to new literals
    {
      BouMS_uint_t varIdx;
      for (varIdx = 0; varIdx < formula->numVariables && !*stop; ++varIdx) {
        BouMS_wcnf_variable_t* const var = formula->variables + varIdx;
        if (ensureLitPtrs(var, var->numLiterals + newLitsPerVar[varIdx], realloc)) {
          break;
        }
      }
      free(newLitsPerVar);
      newLitsPerVar = NULL;
      if (varIdx < formula->numVariables) {  // loop above was broken
        return true;
      }
    }
  }

  // set pointers from variables to literals
  for (BouMS_uint_t clauseIdx = newClausesStart; clauseIdx < newClausesStart + numNewClauses && !*stop; ++clauseIdx) {
    const BouMS_wcnf_clause_t* const clause = formula->clauses + clauseIdx;
    for (BouMS_uint_t litIdx = 0; litIdx < clause->numLiterals && !*stop; ++litIdx) {
      BouMS_wcnf_literal_t* const lit = clause->literals + litIdx;
      lit->clause = clause;

      const BouMS_uint_t varIdx = BouMS_wcnf_var(lit);
      BouMS_wcnf_variable_t* const var = formula->variables + varIdx;
      var->literals[var->numLiterals++] = lit;
    }
  }

  return *stop;
}
