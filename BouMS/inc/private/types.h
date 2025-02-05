/**
 * @file types.h
 * @author Ole Lübke (ole.luebke@tuhh.de)
 * @brief Common (non-public) type definitions
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef TYPES_H
#define TYPES_H

#include "BouMS/BouMS.h"
#include "private/fixedprec.h"

typedef signed_fixedprec_t score_t;

typedef struct {
  unsigned_fixedprec_t* weights;  ///< (updated) weight for each clause
  BouMS_uint_t* numSatLiterals;   ///< the number of sat literals for each clause
  score_t* scores;                ///< the score of each variable
  BouMS_uint_t* satVars;          ///< for each clause, the (only/first) variable which satisfies it

  BouMS_uint_t numDecreasingVars;   ///< the number of decreasing variables
  BouMS_uint_t* decreasingVars;     ///< decreasing variables (with score > 0)
  BouMS_uint_t* decreasingVarsIdx;  ///< for each variable, its index in decrasingVars or BouMS_UINT_MAX

  BouMS_uint_t numFalsifiedHardClauses;   ///< the number of falsified hard clauses
  BouMS_uint_t* falsifiedHardClauses;     ///< falsified hard clauses
  BouMS_uint_t* falsifiedHardClausesIdx;  ///< for each clause, its index in falsifiedHardClauses or BouMS_UINT_MAX

  BouMS_uint_t numFalsifiedSoftClauses;   ///< the number of falsified soft clauses
  BouMS_uint_t* falsifiedSoftClauses;     ///< falsified soft clauses
  BouMS_uint_t* falsifiedSoftClausesIdx;  ///< for each clause, its index in falsifiedSoftClauses or BouMS_UINT_MAX

  BouMS_uint_t numHighWeightSoftClauses;  ///< the number of soft clauses with increased weight
  BouMS_uint_t* highWeightSoftClauses;    ///< soft clauses with increased weight
  BouMS_uint_t*
      highWeightSoftClausesIdx;  ///< for each clause, its index in increasedWeightSoftClauses or BouMS_UINT_MAX

  unsigned_fixedprec_t* tunedWeights;  ///< tuned clause weights (clause weight / avg clause weight)

  void* decimationMem;  ///< memory block for decimation

#ifndef BOUMS_NOFFA
  BouMS_uint_t* varFlipCounts;  ///< counters for how often each variable was flipped
#endif

  BouMS_uint_t numSoftClauses;  ///< the number of soft clauses

  fixedprec_shift_t fixedprecShift;  ///< shift value for fixed precision arithmetic
} BouMS_memory_t;

typedef void (*variable_initializer_t)(BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem,
                                       const BouMS_memoryReq_t* memReq, BouMS_result_t* result, const bool* initModel,
                                       bool* reuseDecimation, bool forceRandomAssignment, const bool* stop);

/**
 * @brief Calculates where pointers in out should point to
 *
 * @param mem
 * @param memReq
 * @param cfg
 * @param formula
 * @param out
 */
void BouMS_initMemory(void* mem, const BouMS_memoryReq_t* memReq, const BouMS_wcnf_t* formula, BouMS_memory_t* out);

#endif
