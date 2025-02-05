/**
 * @file preprocessing.h
 * @author Ole Lübke
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef BOUMS_PREPROCESSING_H
#define BOUMS_PREPROCESSING_H

#include <stdbool.h>

#include "BouMS/common.h"
#include "BouMS/wcnf.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef bool (*BouMS_clause_predicate_t)(const BouMS_wcnf_clause_t*);

typedef struct {
  BouMS_uint_t unassignedVarsMemReq, unassignedVarsIdxMemReq, hardUnitClausesMemReq, hardUnitClausesIdxMemReq,
      softUnitClausesMemReq, softUnitClausesIdxMemReq, prevDeciStepMemReq, clauseRemovedMemReq, numLiteralsMemReq;
} BouMS_decimation_memoryReq_t;

/**
 * @brief Calculate memory requirements for decimation
 *
 * @param formula The formula the decimation will be run on
 * @param memReq Pointer to a structure to save detailed memory requirement information
 * @return BouMS_uint_t Memory requirements in bytes
 */
BouMS_uint_t BouMS_decimation_calcMemoryRequirements(const BouMS_wcnf_t* formula, BouMS_decimation_memoryReq_t* memReq);

/**
 * @brief Run unit-clause propagation to find an advantageous initial variable assignment
 *
 * @param formula The formula to run decimation on
 * @param memory Memory block of size calculated by BouMS_decimation_calcMemoryRequirements
 * @param memReq Obtained from BouMS_decimation_calcMemoryRequirements
 * @param bmsSize Bucket size for random selection
 * @param reusePrevious Whether to reuse previously generated information to steer away the procedure from the variable
 * assignment that was generated during the previous run
 * @param stop Pointer to Boolean value that, when set to true, terminates the alg. asap
 */
void BouMS_decimation(const BouMS_wcnf_t* formula, void* memory, const BouMS_decimation_memoryReq_t* memReq,
                      BouMS_uint_t bmsSize, bool reusePrevious, const bool* stop);

/**
 * @brief Filters clauses according to the given predicates
 *
 * @param clauses The clauses array to filter
 * @param numClauses The length of the clauses array
 * @param predicates The predicate functions; if a predicate returns true for a clause, the clause is removed
 * @param numPredicates The length of the predicates array
 * @return BouMS_uint_t The number of removed clauses
 */
BouMS_uint_t BouMS_filterClauses(BouMS_wcnf_clause_t* clauses, BouMS_uint_t numClauses,
                                 const BouMS_clause_predicate_t* predicates, BouMS_uint_t numPredicates);

#ifdef __cplusplus
}
#endif

#endif
