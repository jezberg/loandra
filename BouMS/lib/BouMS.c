/**
 * @file BouMS.c
 * @author Ole Lübke (ole.luebke@tuhh.de)
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "BouMS/BouMS.h"

#include <assert.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "BouMS/common.h"
#include "BouMS/preprocessing.h"
#include "BouMS/wcnf.h"
#include "private/align.h"
#include "private/common.h"
#include "private/fixedprec.h"
#include "private/logging.h"
#include "private/nonpartial.h"
#include "private/partial.h"
#include "private/rng.h"
#include "private/types.h"

/**
 * @brief Calculates detailed memory requirements
 *
 * @param formula
 * @param out
 */
static void calcMemoryReq(const BouMS_wcnf_t* formula, BouMS_memoryReq_t* out);

/**
 * @brief Copies weights from formula to weights
 *
 * @param formula
 * @param cfg
 * @param mem
 * @param avgSoftClauseWeight
 */
static void initWeights(const BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem,
                        unsigned_fixedprec_t avgSoftClauseWeight);

/**
 * @brief Updates the algorithm's output
 *
 * @param result
 * @param variables
 * @param numVariables
 * @param infeasible
 * @param cost
 *
 * @return Whether the current assignment is better than the previous best assignment
 */
static bool updateResult(BouMS_result_t* result, const BouMS_wcnf_variable_t* variables, BouMS_uint_t numVariables,
                         bool infeasible, BouMS_uint_t cost);

typedef bool (*clausePredicate_t)(const BouMS_wcnf_clause_t*);

/**
 * @brief In-place moves the clauses fulfilling the predicate to the front of the clauses array
 *
 * @param formula
 * @param clausePredicate
 * @param stop
 * @return The number of clauses fulfilling the predicate
 */
static BouMS_uint_t moveClausesFirst(const BouMS_wcnf_t* formula, clausePredicate_t clausePredicate, const bool* stop);

/**
 * @brief Move hard clauses to the front of the array
 *
 * @param formula
 * @param cfg
 * @param stop
 */
static void moveHardClausesFirst(BouMS_wcnf_t* formula, const BouMS_params_t* cfg, const bool* stop);

/**
 * @brief Preprocesses the formula for solving
 *
 * Checks for empty hard clauses
 * Removes duplicate literals from clauses (those can cause crashes)
 * Removes tautologies
 * Removes 0-weight clauses
 *
 * @param formula
 * @param stop
 * @return true if the formula  has empty hard clauses, false otherwise or stop flag was set
 */
static bool prepare(BouMS_wcnf_t* formula, const bool* stop);

/**
 * @brief Analyzes hard clauses to determine the total weight of empty soft clauses, the average soft clause weight, and
 * the fixed precision shift
 *
 * @param formula
 * @param cfg
 * @param mem
 * @param avgSoftClauseWeight
 * @param totalWeightOfEmptySoftClauses
 * @param numEmptySoftClauses
 * @param stop
 */
static void analyzeSoftClauses(BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem,
                               unsigned_fixedprec_t* avgSoftClauseWeight, BouMS_uint_t* totalWeightOfEmptySoftClauses,
                               BouMS_uint_t* numEmptySoftClauses, const bool* stop);

/**
 * @brief Calculates the cost of the current assignment
 * @note The current assignment must be feasible
 *
 * @param formula
 */
static BouMS_uint_t calcCost(const BouMS_wcnf_t* formula);

/**
 * @brief Ensures that all variables are initialized, either using `initModel` or setting all to false.
 *
 * @param formula
 * @param result
 * @param initModel
 */
static void ensureVarsInitialized(BouMS_wcnf_t* formula, BouMS_result_t* result, const bool* initModel);

void BouMS_init(void) {
  rngInit();
}

BouMS_uint_t BouMS_calcMemoryRequirements(const BouMS_wcnf_t* formula, BouMS_memoryReq_t* memReq) {
  calcMemoryReq(formula, memReq);
  return memReq->weightsMemReq + memReq->numSatLiteralsMemReq + memReq->scoresMemReq + memReq->satVarsMemReq +
         memReq->decreasingVarMemReq + memReq->decreasingVarsIdxMemReq + memReq->falsifiedHardClausesMemReq +
         memReq->falsifiedHardClausesIdxMemReq + memReq->falsifiedSoftClausesMemReq +
         memReq->falsifiedSoftClausesIdxMemReq + memReq->highWeightSoftClausesMemReq +
         memReq->hightWeightSoftClausesIdxMemReq + memReq->tunedWeightsMemReq + memReq->varFlipCountMemReq +
         memReq->decimationMemReq;
}

void BouMS_solve(BouMS_wcnf_t* formula, const BouMS_params_t* cfg, void* memory, const BouMS_memoryReq_t* memReq,
                 BouMS_result_t* result, const bool* initModel, const BouMS_uint_t maxNonImprovingFlips,
                 const bool* stop) {
  INIT_DURATION_MEAS();

  START_DURATION_MEAS();
  if (prepare(formula, stop)) {
    result->status = BOUMS_UNSAT;
    return;
  }
  LOG_VERBOSE_WITH_DURATION("c prepared formula");

  moveHardClausesFirst(formula, cfg, stop);

  BouMS_memory_t mem;
  BouMS_initMemory(memory, memReq, formula, &mem);

  const variable_initializer_t initVars = cfg->isPartial ? p_initVars : np_initVars;

  unsigned_fixedprec_t avgSoftClauseWeight = UFIXEDPREC(0);
  BouMS_uint_t totalWeightOfEmptySoftClauses = 0;
  BouMS_uint_t numEmptySoftClauses = 0;
  analyzeSoftClauses(formula, cfg, &mem, &avgSoftClauseWeight, &totalWeightOfEmptySoftClauses, &numEmptySoftClauses,
                     stop);
  LOG_VERBOSE("c found " BOUMS_UINT_FORMAT " empty soft clauses with total weight " BOUMS_UINT_FORMAT "\n",
              numEmptySoftClauses, totalWeightOfEmptySoftClauses);

  result->status = BOUMS_INFEASIBLE;
  result->cost = BOUMS_UINT_MAX;

  ensureVarsInitialized(formula, result, initModel);

  if (*stop) {
    LOG_WARN("c stop flag was set before solving could start - the output may contain wrong information!\n");
  }

  BouMS_uint_t cost;
  bool reuseDecimation = false;
  bool skipNextRestart = false;
  BouMS_uint_t triesWOImprovement = 0;
  BouMS_uint_t flipsWOImprovement = 0;

  mem.numFalsifiedHardClauses = BOUMS_UINT_MAX;  // so done(...) check below is false on first check
  for (BouMS_uint_t try = 0; try < cfg->maxTries && !done(&mem, result) && !*stop &&
                             (!maxNonImprovingFlips || flipsWOImprovement < maxNonImprovingFlips);
       ++try) {
    LOG_VERBOSE("c try " BOUMS_UINT_FORMAT "\n", try);

    if (!skipNextRestart) {
      const bool forceRandomAssignment = triesWOImprovement >= cfg->maxTriesWOImprovement;
      if (forceRandomAssignment) {
        LOG_VERBOSE("c restarting with random assignment\n");
      }

      START_DURATION_MEAS();
      initVars(formula, cfg, &mem, memReq, result, initModel, &reuseDecimation, forceRandomAssignment, stop);
      LOG_TRACE_WITH_DURATION("c initialized variables");

      START_DURATION_MEAS();
      BouMS_decimation(formula, mem.decimationMem, &memReq->detailedDecimationMemReq, cfg->bmsSize, reuseDecimation,
                       stop);
      reuseDecimation = true;
      LOG_TRACE_WITH_DURATION("c performed decimation");

      START_DURATION_MEAS();
      initWeights(formula, cfg, &mem, avgSoftClauseWeight);
      LOG_TRACE_WITH_DURATION("c initialized weights");

      START_DURATION_MEAS();
      initAlgo(formula, &mem, &cost);
      cost += totalWeightOfEmptySoftClauses;
      LOG_TRACE_WITH_DURATION("c initialized algorithm");
    } else {
      LOG_VERBOSE("c skipping restart\n");
    }

    skipNextRestart = false;

    for (BouMS_uint_t flip = 0; flip < cfg->maxFlips && !done(&mem, result) && !*stop &&
                                (!maxNonImprovingFlips || flipsWOImprovement <= maxNonImprovingFlips);
         ++flip) {
      const bool infeasible = mem.numFalsifiedHardClauses > 0;
      if (updateResult(result, formula->variables, formula->numVariables, infeasible, cost)) {
        flipsWOImprovement = 0;
        skipNextRestart = true;
      }

      // LOG_TRACE("c current cost is " BouMS_UINT_FORMAT " (%s)\n", cost, infeasible ? "infeasible" : "feasible");

      BouMS_wcnf_variable_t* variableToFlip = selectVariable(formula, cfg, &mem);
      flipVariable(variableToFlip, formula, &mem, &cost);
      ++flipsWOImprovement;
    }

    if (updateResult(result, formula->variables, formula->numVariables, mem.numFalsifiedHardClauses > 0, cost)) {
      flipsWOImprovement = 0;
      skipNextRestart = true;
    }

    if (skipNextRestart) {
      triesWOImprovement = 0;
    } else {
      triesWOImprovement += 1;
    }
  }

  formula->numClauses += numEmptySoftClauses;
}

void BouMS_params(const BouMS_wcnf_t* formula, BouMS_params_t* params) {
  params->maxTries = BOUMS_UINT_MAX;
  params->fixedPrecSafetyBits = 32;  // NOLINT(readability-magic-numbers)
  params->isPartial = formula->numHardClauses > 0;
  params->isWeighted = false;

  const BouMS_uint_t numSoftClauses = formula->numClauses - formula->numHardClauses;
  BouMS_uint_t softClauseCnt = 0;
  BouMS_uint_t firstWeight = 0;
  if (numSoftClauses > 0) {
    for (BouMS_uint_t i = 0; i < formula->numClauses; ++i) {
      const BouMS_wcnf_clause_t* const clause = formula->clauses + i;
      if (!BouMS_wcnf_isClauseHard(clause) && clause->weight > 0) {
        if (firstWeight == 0) {
          firstWeight = clause->weight;
        } else if (clause->weight != firstWeight) {
          params->isWeighted = true;
          break;
        }

        if (++softClauseCnt == numSoftClauses) {
          break;
        }
      }
    }
  }

  params->maxTriesWOImprovement = params->isPartial ? 3 : 2;  // NOLINT(readability-magic-numbers)

  /* NuWLS parameter values from
   * Yi Chu, Shaowei Cai, Chuan Luo, NuWLS-c-2023: Solver Description. MaxSAT Evaluation 2023 (source code)
   */
  params->maxFlips = 10000000;  // NOLINT(readability-magic-numbers)
  if (!params->isWeighted) {
    params->sWeightInc = 1;  // NOLINT(readability-magic-numbers)

    if (formula->numHardClauses == 0) {
      params->bmsSize = 94;              // NOLINT(readability-magic-numbers)
      params->hWeightInit = 0;           // NOLINT(readability-magic-numbers)
      params->hWeightInc = 0;            // NOLINT(readability-magic-numbers)
      params->sSmoothProb = 200000;      // NOLINT(readability-magic-numbers)
      params->sWeightBoundCoef = 397;    // NOLINT(readability-magic-numbers)
      params->sWeightBoundOffset = 550;  // NOLINT(readability-magic-numbers)
    } else {
      params->bmsSize = 50;            // NOLINT(readability-magic-numbers)
      params->hWeightInit = 1;         // NOLINT(readability-magic-numbers)
      params->hWeightInc = 1;          // NOLINT(readability-magic-numbers)
      params->sSmoothProb = 0;         // NOLINT(readability-magic-numbers)
      params->sWeightBoundCoef = 0;    // NOLINT(readability-magic-numbers)
      params->sWeightBoundOffset = 0;  // NOLINT(readability-magic-numbers)
    }
  } else {
    if (formula->numHardClauses == 0) {
      params->bmsSize = 22;             // NOLINT(readability-magic-numbers)
      params->hWeightInit = 0;          // NOLINT(readability-magic-numbers)
      params->hWeightInc = 0;           // NOLINT(readability-magic-numbers)
      params->sWeightInc = 3;           // NOLINT(readability-magic-numbers)
      params->sSmoothProb = 100000;     // NOLINT(readability-magic-numbers)
      params->sWeightBoundCoef = 1000;  // NOLINT(readability-magic-numbers)
      params->sWeightBoundOffset = 10;  // NOLINT(readability-magic-numbers)
    } else {
      params->bmsSize = 50;            // NOLINT(readability-magic-numbers)
      params->hWeightInit = 1;         // NOLINT(readability-magic-numbers)
      params->hWeightInc = 5;          // NOLINT(readability-magic-numbers)
      params->sWeightInc = 1;          // NOLINT(readability-magic-numbers)
      params->sSmoothProb = 0;         // NOLINT(readability-magic-numbers)
      params->sWeightBoundCoef = 0;    // NOLINT(readability-magic-numbers)
      params->sWeightBoundOffset = 0;  // NOLINT(readability-magic-numbers)
    }
  }
}

static void calcMemoryReq(const BouMS_wcnf_t* formula, BouMS_memoryReq_t* out) {
  out->weightsMemReq = align(formula->numClauses * sizeof(unsigned_fixedprec_t));
  out->numSatLiteralsMemReq = align(formula->numClauses * sizeof(BouMS_uint_t));
  out->scoresMemReq = align(formula->numVariables * sizeof(score_t));
  out->satVarsMemReq = align(formula->numClauses * sizeof(BouMS_uint_t));

  out->decreasingVarMemReq = align(formula->numVariables * sizeof(BouMS_uint_t));
  out->decreasingVarsIdxMemReq = align(formula->numVariables * sizeof(BouMS_uint_t));

  out->falsifiedHardClausesMemReq = align(formula->numClauses * sizeof(BouMS_uint_t));
  out->falsifiedHardClausesIdxMemReq = align(formula->numClauses * sizeof(BouMS_uint_t));

  out->falsifiedSoftClausesMemReq = align(formula->numClauses * sizeof(BouMS_uint_t));
  out->falsifiedSoftClausesIdxMemReq = align(formula->numClauses * sizeof(BouMS_uint_t));

  out->highWeightSoftClausesMemReq = align(formula->numClauses * sizeof(BouMS_uint_t));
  out->hightWeightSoftClausesIdxMemReq = align(formula->numClauses * sizeof(BouMS_uint_t));

  out->tunedWeightsMemReq = align(formula->numClauses * sizeof(unsigned_fixedprec_t));

  out->varFlipCountMemReq = align(formula->numVariables * sizeof(BouMS_uint_t));

  out->decimationMemReq = align(BouMS_decimation_calcMemoryRequirements(formula, &out->detailedDecimationMemReq));
}

void BouMS_initMemory(void* mem, const BouMS_memoryReq_t* memReq, const BouMS_wcnf_t* formula, BouMS_memory_t* out) {
  out->weights = mem;
  out->numSatLiterals = (BouMS_uint_t*)((uint8_t*)out->weights + memReq->weightsMemReq);
  out->scores = (score_t*)((uint8_t*)out->numSatLiterals + memReq->numSatLiteralsMemReq);
  out->satVars = (BouMS_uint_t*)((uint8_t*)out->scores + memReq->scoresMemReq);

  out->numDecreasingVars = 0;
  out->decreasingVars = (BouMS_uint_t*)((uint8_t*)out->satVars + memReq->satVarsMemReq);
  out->decreasingVarsIdx = (BouMS_uint_t*)((uint8_t*)out->decreasingVars + memReq->decreasingVarMemReq);

  out->numFalsifiedHardClauses = 0;
  out->falsifiedHardClauses = (BouMS_uint_t*)((uint8_t*)out->decreasingVarsIdx + memReq->decreasingVarsIdxMemReq);
  out->falsifiedHardClausesIdx =
      (BouMS_uint_t*)((uint8_t*)out->falsifiedHardClauses + memReq->falsifiedHardClausesMemReq);

  out->numFalsifiedSoftClauses = 0;
  out->falsifiedSoftClauses =
      (BouMS_uint_t*)((uint8_t*)out->falsifiedHardClausesIdx + memReq->falsifiedSoftClausesIdxMemReq);
  out->falsifiedSoftClausesIdx =
      (BouMS_uint_t*)((uint8_t*)out->falsifiedSoftClauses + memReq->falsifiedSoftClausesMemReq);

  out->numHighWeightSoftClauses = 0;
  out->highWeightSoftClauses =
      (BouMS_uint_t*)((uint8_t*)out->falsifiedSoftClausesIdx + memReq->falsifiedSoftClausesIdxMemReq);
  out->highWeightSoftClausesIdx =
      (BouMS_uint_t*)((uint8_t*)out->highWeightSoftClauses + memReq->highWeightSoftClausesMemReq);

  out->tunedWeights =
      (unsigned_fixedprec_t*)((uint8_t*)out->highWeightSoftClausesIdx + memReq->hightWeightSoftClausesIdxMemReq);

  out->decimationMem = (void*)((uint8_t*)out->tunedWeights + memReq->tunedWeightsMemReq);

#ifndef BOUMS_NOFFA
  out->varFlipCounts = (BouMS_uint_t*)((uint8_t*)out->decimationMem + memReq->decimationMemReq);
  // initialize here for global count instead of in initAlgo (try-local count)
  for (BouMS_uint_t varIdx = 0; varIdx < memReq->varFlipCountMemReq / sizeof(BouMS_uint_t); ++varIdx) {
    out->varFlipCounts[varIdx] = 0;
  }
#endif

  out->numSoftClauses = formula->numClauses - formula->numHardClauses;
}

static void initWeights(const BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem,
                        unsigned_fixedprec_t avgSoftClauseWeight) {
  if (cfg->isPartial && cfg->isWeighted) {
    pw_initWeights(formula, cfg, mem, avgSoftClauseWeight);
  } else if (cfg->isPartial) {
    puw_initWeights(formula, cfg, mem);
  } else if (cfg->isWeighted) {
    npw_initWeights(formula, cfg, mem, avgSoftClauseWeight);
  } else {
    npuw_initWeights(formula, cfg, mem);
  }
}

static bool updateResult(BouMS_result_t* result, const BouMS_wcnf_variable_t* variables, BouMS_uint_t numVariables,
                         bool infeasible, BouMS_uint_t cost) {
  if (!infeasible && cost < result->cost) {
    saveAssignment(numVariables, variables, result->assignment);  // save currently best assignment
    result->cost = cost;
    result->status = cost == 0 ? BOUMS_OPTIMUM_FOUND : BOUMS_UNKNOWN;
    LOG_VERBOSE("o " BOUMS_UINT_FORMAT "\n", result->cost);
    return true;
  }
  return false;
}

static BouMS_uint_t moveClausesFirst(const BouMS_wcnf_t* formula, clausePredicate_t clausePredicate, const bool* stop) {
  // in-place sorting clauses to beginning of array
  BouMS_uint_t numClauses = 0;
  BouMS_uint_t lastClauseIdx = 0;

  for (BouMS_uint_t clauseIdx = 0; clauseIdx < formula->numClauses && !*stop; ++clauseIdx) {
    BouMS_wcnf_clause_t* const clause = formula->clauses + clauseIdx;

    if (clausePredicate(clause)) {
      numClauses += 1;
      continue;
    }

    // search next clause
    BouMS_wcnf_clause_t* nextClause = NULL;
    for (BouMS_uint_t nextClauseSearchIdx = (clauseIdx > lastClauseIdx ? clauseIdx : lastClauseIdx) + 1;
         nextClauseSearchIdx < formula->numClauses; ++nextClauseSearchIdx) {
      BouMS_wcnf_clause_t* candidate = formula->clauses + nextClauseSearchIdx;
      if (clausePredicate(candidate)) {
        lastClauseIdx = nextClauseSearchIdx;
        nextClause = candidate;
        break;
      }
    }

    if (!nextClause) {
      break;  // all clauses found
    }

    // swap places
    const BouMS_wcnf_clause_t tmp = *clause;
    *clause = *nextClause;
    *nextClause = tmp;

    // update back pointers
    fixLitToClausePtrs(clause);
    fixLitToClausePtrs(nextClause);

    numClauses += 1;
  }

  return numClauses;
}

static void moveHardClausesFirst(BouMS_wcnf_t* formula, const BouMS_params_t* cfg, const bool* stop) {
  INIT_DURATION_MEAS();

  if (cfg->isPartial) {
    START_DURATION_MEAS();
    const BouMS_uint_t numHardClauses = moveClausesFirst(formula, BouMS_wcnf_isClauseHard, stop);
    LOG_VERBOSE_WITH_DURATION("c separated hard and soft clauses");
    assert(numHardClauses == formula->numHardClauses);
    ((void)numHardClauses);  // suppress -Wunused-variable
  }
}

static bool prepare(BouMS_wcnf_t* formula, const bool* stop) {
  for (BouMS_uint_t clauseIdx = 0; clauseIdx < formula->numClauses && !*stop; ++clauseIdx) {
    BouMS_wcnf_clause_t* const clause = formula->clauses + clauseIdx;
    const bool isHard = BouMS_wcnf_isClauseHard(clause);

    if (BouMS_wcnf_cleanClause(formula, clause) || clause->weight == 0) {
      removeClause(clauseIdx--, formula->clauses, formula->numClauses--);
      if (isHard) {
        formula->numHardClauses -= 1;
      }
    } else if (isHard && BouMS_wcnf_isClauseEmpty(clause)) {
      return true;
    }
  }

  return false;
}

static void analyzeSoftClauses(BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem,
                               unsigned_fixedprec_t* avgSoftClauseWeight, BouMS_uint_t* totalWeightOfEmptySoftClauses,
                               BouMS_uint_t* numEmptySoftClauses, const bool* stop) {
  /* loop over soft clauses to
   *   calculate total weight and determine shift value for fixed precision arithmetic
   *   calculate average soft clause weight for NuWLS 2023 tuned weight calculation (cf. initWeights)
   *   find empty soft clauses
   */
  BouMS_uint_t totalWeight = 0;
  *totalWeightOfEmptySoftClauses = 0;
  *numEmptySoftClauses = 0;

  for (BouMS_uint_t clauseIdx = formula->numHardClauses; clauseIdx < formula->numClauses && !*stop; ++clauseIdx) {
    const BouMS_wcnf_clause_t* const clause = formula->clauses + clauseIdx;
    assert(!BouMS_wcnf_isClauseHard(clause));
    assert(totalWeight <= totalWeight + clause->weight);
    totalWeight += clause->weight;

    // if soft clause is empty, note its weight and ignore it from now on
    if (clause->numLiterals == 0) {
      *totalWeightOfEmptySoftClauses += clause->weight;

      // move clause to end of array ...
      const BouMS_uint_t moveIdx = formula->numClauses - *numEmptySoftClauses - 1;
      const BouMS_wcnf_clause_t moveClause = formula->clauses[moveIdx];
      formula->clauses[moveIdx] = *clause;
      formula->clauses[clauseIdx] = moveClause;

      *numEmptySoftClauses += 1;
      clauseIdx -= 1;
      formula->numClauses -= 1;  // ... and ignore it
    }
  }

  mem->fixedprecShift = cfg->isWeighted ? fixedprec_determineMaxShift(totalWeight, cfg->fixedPrecSafetyBits) : 0;
  LOG_INFO("c fixed precision shift is %u bits\n", mem->fixedprecShift);

  *avgSoftClauseWeight =
      mem->numSoftClauses > 0
          ? fixedprec_udiv(fixedprec_uto(totalWeight, mem->fixedprecShift),
                           fixedprec_uto(mem->numSoftClauses, mem->fixedprecShift), mem->fixedprecShift)
          : UFIXEDPREC(0);
  mem->numSoftClauses -= *numEmptySoftClauses;
#ifndef FIXEDPREC_FLOATING
  LOG_INFO("c average soft clause weight is " BOUMS_UINT_FORMAT "\n",
           fixedprec_ufrom(*avgSoftClauseWeight, mem->fixedprecShift));
#else
  LOG_INFO("c average soft clause weight is %f\n", fixedprec_ufrom(*avgSoftClauseWeight, mem->fixedprecShift));
#endif
}

static BouMS_uint_t calcCost(const BouMS_wcnf_t* formula) {
  BouMS_uint_t cost = 0;

  for (BouMS_uint_t clauseIdx = formula->numHardClauses; clauseIdx < formula->numClauses; ++clauseIdx) {
    const BouMS_wcnf_clause_t* const clause = formula->clauses + clauseIdx;
    assert(!BouMS_wcnf_isClauseHard(clause));

    bool sat = false;
    for (BouMS_uint_t litIdx = 0; litIdx < clause->numLiterals; ++litIdx) {
      const BouMS_wcnf_literal_t* literal = clause->literals + litIdx;

      if (BouMS_wcnf_isLiteralSatisfied(formula, literal)) {
        sat = true;
        break;
      }
    }

    if (!sat) {
      cost += clause->weight;
    }
  }

  return cost;
}

static void ensureVarsInitialized(BouMS_wcnf_t* formula, BouMS_result_t* result, const bool* initModel) {
  if (initModel) {
    useAssignment(formula->numVariables, initModel, formula->variables);
    updateResult(result, formula->variables, formula->numVariables, !isFeasible(formula), calcCost(formula));
  } else {
    // make sure we don't branch on uninitialized variable values by initializing all vars to false
    for (BouMS_uint_t varIdx = 0; varIdx < formula->numVariables; ++varIdx) {
      formula->variables[varIdx].value = false;
    }
  }
}
