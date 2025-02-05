/**
 * @file wcnf_util.h
 * @author Ole Lübke (ole.luebke@tuhh.de)
 * @brief Utility functions that rely on dynamic memory management (realloc, free) for creating WCNF formulas
 *
 * @copyright Copyright (c) 2024
 *
 */
#ifndef BOUMS_WCNF_UTIL
#define BOUMS_WCNF_UTIL

#include <stdbool.h>
#include <stdlib.h>

#include "BouMS/common.h"
#include "BouMS/wcnf.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Type of realloc
 */
typedef void* (*BouMS_wcnf_util_realloc_t)(void*, size_t);

/**
 * @brief Type of free
 */
typedef void (*BouMS_wcnf_util_free_t)(void*);

/**
 * @brief Create a literal from MSE style input
 *
 * @see https://maxsat-evaluations.github.io/2024/rules.html#input
 *
 * @param l Input literal, i.e., an integer (negative for negated literal)
 * @return Cf. BouMS_literal_t
 */
static inline BouMS_literal_t BouMS_wcnf_util_mkLit(int l) {
  return BouMS_mkLit(abs(l) - 1, l < 0);
}

/**
 * @brief Creates a new, empty formula
 *
 * @return A new, empty formula
 */
static inline BouMS_wcnf_t BouMS_wcnf_util_newFormula(void) {
  return (BouMS_wcnf_t){.numVariables = 0, .numClauses = 0, .numHardClauses = 0, .variables = NULL, .clauses = NULL};
}

/**
 * @brief Deletes the given formula by free'ing all attached memory and resetting the respective counters
 *
 * @note Does not call free on formula itself
 *
 * @param formula The formula to delete
 * @param free Pointer to free function
 * @param stop NULL or pointer to flag to stop early (i.e., some memory will not be cleaned up)
 */
void BouMS_wcnf_util_deleteFormula(BouMS_wcnf_t* formula, BouMS_wcnf_util_free_t free, const bool* stop);

/**
 * @brief Add a clause from MSE style input
 *
 * @see https://maxsat-evaluations.github.io/2024/rules.html#input
 *
 * @note For adding multiple clauses, better use batch clause adding (below) as it requires fewer reallocations
 *
 * @param formula The formula to add the clause to
 * @param weight The weight of the clause, BOUMS_HARD_CLAUSE_WEIGHT for hard clauses
 * @param literals The input literals in MSE style (i.e., integers)
 * @param numLiterals The number of literals
 * @param realloc Pointer to realloc function for dynamic memory management
 * @param free Pointer to free function for cleanup in case of error
 * @return true In case of error (e.g., out of memory)
 * @return false Otherwise
 */
bool BouMS_wcnf_util_addClause(BouMS_wcnf_t* formula, BouMS_uint_t weight, const int* literals,
                               BouMS_uint_t numLiterals, BouMS_wcnf_util_realloc_t realloc,
                               BouMS_wcnf_util_free_t free);

/**
 * @brief Struct that holds data for batch clause adding
 */
typedef struct {
  BouMS_uint_t capacity;              ///< The size of the clause buffer
  BouMS_uint_t size;                  ///< The number of clauses in the buffer
  BouMS_wcnf_clause_t* clauses;       ///< The clause buffer
  BouMS_uint_t maxVar;                ///< The highest variable number encountered during adding
  BouMS_wcnf_util_realloc_t realloc;  ///< realloc function
  BouMS_wcnf_util_free_t free;        ///< free function
} BouMS_wcnf_util_batchClauseAddingState_t;

/**
 * @brief Start batch clause adding
 *
 * @param expectedNumClauses The expected number of clauses that will be added, 0 if unknown
 * @param realloc realloc function
 * @param free free function
 * @param state Some (unitialized) BouMS_wcnf_util_batchClauseAddingState_t
 * @return true In case of error (e.g., out of memory)
 * @return false Otherwise
 */
bool BouMS_wcnf_util_startBatchClauseAdding(BouMS_uint_t expectedNumClauses, BouMS_wcnf_util_realloc_t realloc,
                                            BouMS_wcnf_util_free_t free,
                                            BouMS_wcnf_util_batchClauseAddingState_t* state);

/**
 * @brief Add all the clauses from the batch to the formula
 *
 * @param state Obtained from BouMS_wcnf_util_startBatchClauseAdding
 * @param formula The formula to add the clauses to
 * @param stop NULL or pointer to flag to stop early (i.e., some clauses will not be added properly)
 * @return true In case of error (e.g., out of memory); the input formula and state are left intact
 * @return true When the stop flag was set; the input formula is not guaranteed to be intact, but the state still is
 * @return false Otherwise (state has been cleaned up)
 */
bool BouMS_wcnf_util_finishBatchClauseAdding(BouMS_wcnf_util_batchClauseAddingState_t* state, BouMS_wcnf_t* formula,
                                             const bool* stop);

/**
 * @brief Clean up BouMS_wcnf_util_batchClauseAddingState_t in case BouMS_wcnf_util_finishBatchClauseAdding failed
 *
 * @param state The state to clean up
 */
void BouMS_wcnf_util_cleanUpBatchClauseAddingAfterError(BouMS_wcnf_util_batchClauseAddingState_t* state);

/**
 * @brief Add a clause from MSE style input to the batch
 *
 * @see https://maxsat-evaluations.github.io/2024/rules.html#input
 *
 * @param state Obtained from BouMS_wcnf_util_startBatchClauseAdding
 * @param weight The weight of the clause, BOUMS_HARD_CLAUSE_WEIGHT for hard clauses
 * @param literals The input literals in MSE style (i.e., integers)
 * @param numLiterals The number of literals
 * @return true In case of error (e.g., out of memory)
 * @return false Otherwise
 */
static inline bool BouMS_wcnf_util_batchAddClause(BouMS_wcnf_util_batchClauseAddingState_t* state, BouMS_uint_t weight,
                                                  const int* literals, BouMS_uint_t numLiterals) {
  if (state->size == state->capacity) {
    const BouMS_uint_t newCapacity = 2 * state->capacity;
    BouMS_wcnf_clause_t* const newClauses =
        (BouMS_wcnf_clause_t*)state->realloc(state->clauses, newCapacity * sizeof(BouMS_wcnf_clause_t));
    if (!newClauses) {
      return true;
    }
    state->clauses = newClauses;
    state->capacity = newCapacity;
  }

  BouMS_wcnf_clause_t* const clause = state->clauses + state->size;

  clause->weight = weight;

  clause->literals = (BouMS_wcnf_literal_t*)state->realloc(NULL, numLiterals * sizeof(BouMS_wcnf_literal_t));
  if (!clause->literals) {
    return true;
  }
  clause->numLiterals = numLiterals;

  for (BouMS_uint_t litIdx = 0; litIdx < numLiterals; ++litIdx) {
    BouMS_wcnf_literal_t* const lit = clause->literals + litIdx;
    lit->lit = BouMS_wcnf_util_mkLit(literals[litIdx]);
    lit->clause = NULL;  // to be set in finishBatchClauseAdding
    const BouMS_uint_t varNum = BouMS_wcnf_var(lit) + 1;
    if (varNum > state->maxVar) {
      state->maxVar = varNum;
    }
  }

  state->size += 1;

  return false;
}

#ifdef __cplusplus
}
#endif

#endif
