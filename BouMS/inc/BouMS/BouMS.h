/**
 * @file BouMS.h
 * @author Ole Lübke (ole.luebke@tuhh.de)
 * @brief
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef BOUMS_H
#define BOUMS_H

#include <stdbool.h>

#include "BouMS/common.h"
#include "BouMS/preprocessing.h"
#include "BouMS/wcnf.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Holds the status of incomplete MaxSAT solving
 *
 */
typedef enum {
  BOUMS_OPTIMUM_FOUND = 0,  //!< The hard clauses are satisfiable and the optimum cost has been found
  BOUMS_UNKNOWN,            //!< The hard clauses are satisfiable but the optimum cost is unknown
  BOUMS_UNSAT,              //!< The hard clauses are unsatisfiable
  BOUMS_INFEASIBLE,         //!< The hard clauses may or may not be satisfiable
  BOUMS_STATUS_COUNT        //!< Not an actual status, just the number of possible status
} BouMS_status_t;

/**
 * @brief Holds the result of incomplete MaxSAT solving
 *
 */
typedef struct {
  BouMS_status_t status;  //!< The solution status
  BouMS_uint_t cost;      //!< The solution cost (sum of weights of clauses falsified by the solution)
  bool* assignment;       //!< The solution assignment, memory must be allocated by user
} BouMS_result_t;

/**
 * @brief Holds parameters of the algorithm
 *
 * @ref clause_weighting_params.md
 */
typedef struct {
  BouMS_uint_t maxTries; /*!< The maximum number of restarts (i.e., how often a new initial assignment is generated) */
  BouMS_uint_t maxTriesWOImprovement; /*!< The maximum number of tries without improvement before restart with a new
                                        random assignment is forced */
  BouMS_uint_t maxFlips; /*!< The maximum number of variable flips the algorithm performs within one iteration (total
                             no. of iterations <= maxTries * maxFlips) */
  BouMS_uint_t bmsSize; /*!< The number of decreasing variables selected for "best from multiple selections" strategy */
  BouMS_uint_t hWeightInit;      /*!< The initial weight for hard clauses, only for PMS, WPMS */
  BouMS_uint_t hWeightInc;       /*!< Weight increment for hard clauses during clause reweighting, only for PMS, WPMS */
  BouMS_uint_t sWeightInc;       /*!< Weight increment for soft clauses during clause reweighting */
  BouMS_uint_t sWeightBoundCoef; /*!< Coefficient for calculating soft weight upper bound, only for MS, WMS */
  BouMS_uint_t sWeightBoundOffset; /*!< Offset for calculating soft weight upper bound, only for MS, WMS */
  BouMS_uint_t sSmoothProb;        /*!< Smoothing probability for soft clauses (100000000 == 100%), only for MS, WMS */
  BouMS_uint_t
      fixedPrecSafetyBits; /*!< Number of reserved bits to prevent overflows in fixed precision weight calculations.
                              When NDEBUG is not defined, overflows are asserted. Set to 8*sizeof(BouMS_uint_t) to fall
                              back to usual integer arithmetic. Only fow WMS, WPMS */

  bool isPartial;  /*!< True if there are hard clauses, false otherwise */
  bool isWeighted; /*!< False if all soft clauses have weight 1, true otherwise */
} BouMS_params_t;

/**
 * @brief Memory requirements
 */
typedef struct {
  BouMS_uint_t weightsMemReq, numSatLiteralsMemReq, scoresMemReq, satVarsMemReq, decreasingVarMemReq,
      decreasingVarsIdxMemReq, falsifiedHardClausesMemReq, falsifiedHardClausesIdxMemReq, falsifiedSoftClausesMemReq,
      falsifiedSoftClausesIdxMemReq, highWeightSoftClausesMemReq, hightWeightSoftClausesIdxMemReq, tunedWeightsMemReq,
      varFlipCountMemReq, decimationMemReq;
  BouMS_decimation_memoryReq_t detailedDecimationMemReq;
} BouMS_memoryReq_t;

/**
 * @brief Initializes the library
 *
 */
void BouMS_init(void);

/**
 * @brief Calculates memory requirements for the algorithm
 *
 * @param formula The formula the algorithm is supposed to run on
 * @param memReq
 * @return bsat_uint_t The number of bytes required to run the algorithm on the formula
 */
BouMS_uint_t BouMS_calcMemoryRequirements(const BouMS_wcnf_t* formula, BouMS_memoryReq_t* memReq);

/**
 * @brief Tries to solve MaxSAT instance
 *
 * @param formula The formula to try to solve
 * @param cfg @see BouMS_params_t
 * @param memory A memory block of sufficient size (calculate with BouMS_calcMemoryRequirements)
 * @param memReq @see BouMS_memoryReq_t
 * @param result @see BouMS_result_t
 * @param initModel Pointer to a model from which local search should start or NULL
 * @param maxNonImprovingFlips Stop when no better assignment was found within this number of flips
 * @param stop Pointer to a Boolean value that, when set to true, causes the procedure to terminate asap
 */
void BouMS_solve(BouMS_wcnf_t* formula, const BouMS_params_t* cfg, void* memory, const BouMS_memoryReq_t* memReq,
                 BouMS_result_t* result, const bool* initModel, BouMS_uint_t maxNonImprovingFlips, const bool* stop);

/**
 * @brief Calculates default parameters based on the given formula
 *
 * @param formula
 * @param params
 */
void BouMS_params(const BouMS_wcnf_t* formula, BouMS_params_t* params);

#ifdef __cplusplus
}
#endif

#endif
