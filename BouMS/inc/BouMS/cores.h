/**
 * @file cores.h
 * @author Ole Lübke (ole.luebke@tuhh.de)
 *
 * @copyright Copyright (c) 2025
 */

#ifndef BOUMS_CORES_H
#define BOUMS_CORES_H

#include "BouMS/BouMS.h"
#include "BouMS/dynmem.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  BouMS_uint_t numLiterals;
  BouMS_literal_t* literals;
} BouMS_cores_core_t;

typedef struct {
  BouMS_uint_t numCores;
  const BouMS_cores_core_t** cores;
} BouMS_cores_list_t;

typedef struct {
  BouMS_uint_t numCores;
  const BouMS_cores_core_t* cores;
  const BouMS_uint_t* varToClause;
  BouMS_cores_list_t** varToCores;
  BouMS_uint_t* coreToSatLits;
} BouMS_cores_mem_t;

/**
 * @param formula
 * @param cores
 * @param numCores
 * @param varToClause
 * @param out
 * @param realloc
 * @param free
 * @return true in case of errror
 * @return false otherwise
 */
bool BouMS_cores_init(const BouMS_wcnf_t* formula, const BouMS_cores_core_t* cores, BouMS_uint_t numCores,
                      const BouMS_uint_t* varToClause, BouMS_cores_mem_t* out, BouMS_dynmem_realloc_t realloc,
                      BouMS_dynmem_free_t free);

/**
 * @param formula
 * @param cores
 * @param free
 */
void BouMS_cores_free(const BouMS_wcnf_t* formula, BouMS_cores_mem_t* cores, BouMS_dynmem_free_t free);

/**
 * @brief Tries to solve MaxSAT instance with the help of cores
 *
 * @param formula The formula to try to solve
 * @param cores
 * @param cfg @see BouMS_params_t
 * @param memory A memory block of sufficient size (calculate with BouMS_calcMemoryRequirements)
 * @param memReq @see BouMS_memoryReq_t
 * @param result @see BouMS_result_t
 * @param initModel Pointer to a model from which local search should start or NULL
 * @param map
 * @param maxNonImprovingFlips Stop when no better assignment was found within this number of flips
 * @param stop Pointer to a Boolean value that, when set to true, causes the procedure to terminate asap
 */
void BouMS_cores_solve(BouMS_wcnf_t* formula, BouMS_cores_mem_t* cores, const BouMS_params_t* cfg, void* memory,
                       const BouMS_memoryReq_t* memReq, BouMS_result_t* result, const bool* initModel,
                       BouMS_clauseMap_t* map, BouMS_uint_t maxNonImprovingFlips, const bool* stop);

#ifdef __cplusplus
}
#endif

#endif
