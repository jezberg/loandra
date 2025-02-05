/**
 * @file partial.c
 * @author Ole Lübke (ole.luebke@tuhh.de)
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "private/partial.h"

#include <assert.h>
#include <stdbool.h>

#include "BouMS/BouMS.h"
#include "BouMS/common.h"
#include "BouMS/preprocessing.h"
#include "BouMS/wcnf.h"
#include "private/common.h"
#include "private/fixedprec.h"
#include "private/logging.h"
#include "private/types.h"

/**
 * @brief
 *
 * @param formula
 * @param cfg
 * @param mem
 * @param memReq
 * @param result
 * @param initModel
 * @param forceRandom
 * @param reuseDecimation
 * @param stop
 */
static void solveHard(BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem, const BouMS_memoryReq_t* memReq,
                      const BouMS_result_t* result, const bool* initModel, bool forceRandom, bool* reuseDecimation,
                      const bool* stop);

/**
 * @brief
 *
 * @param formula
 * @param cfg
 * @param mem
 */
static void increaseHardWeights(const BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem);

/**
 * @brief
 *
 * @param formula
 * @param cfg
 * @param mem
 */
static void w_increaseSoftWeights(const BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem);

/**
 * @brief
 *
 * @param formula
 * @param cfg
 * @param mem
 */
static void uw_increaseSoftWeights(const BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem);

void p_initVars(BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem, const BouMS_memoryReq_t* memReq,
                BouMS_result_t* result, const bool* initModel, bool* reuseDecimation, bool forceRandomAssignment,
                const bool* stop) {
  if (!isFeasible(formula)) {
    const BouMS_uint_t origNumClauses = formula->numClauses;
    formula->numClauses = formula->numHardClauses;

    BouMS_params_t params = {
        .maxTries = cfg->maxTries,
        .maxFlips = cfg->maxFlips,
        .bmsSize = cfg->bmsSize,
        .isPartial = true,
        .isWeighted = false,
        .hWeightInit = 1,
        .hWeightInc = 1,
        // below params are irrelevant to hard clause solving
        .sWeightInc = 0,
        .sSmoothProb = 0,
        .sWeightBoundCoef = 0,
        .sWeightBoundOffset = 0,
        .maxTriesWOImprovement = 0,
    };

    solveHard(formula, &params, mem, memReq, result, initModel, forceRandomAssignment, reuseDecimation, stop);

    formula->numClauses = origNumClauses;
  }
}

void pw_initWeights(const BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem,
                    unsigned_fixedprec_t avgSoftClauseWeight) {
  const unsigned_fixedprec_t hWeightInit = fixedprec_uto(cfg->hWeightInit, mem->fixedprecShift);

  for (BouMS_uint_t clauseIdx = 0; clauseIdx < formula->numHardClauses; ++clauseIdx) {
    assert(BouMS_wcnf_isClauseHard(formula->clauses + clauseIdx));
    mem->weights[clauseIdx] = hWeightInit;
  }

  for (BouMS_uint_t clauseIdx = formula->numHardClauses; clauseIdx < formula->numClauses; ++clauseIdx) {
    assert(!BouMS_wcnf_isClauseHard(formula->clauses + clauseIdx));
    mem->tunedWeights[clauseIdx] =
        fixedprec_udiv(fixedprec_uto(formula->clauses[clauseIdx].weight, mem->fixedprecShift), avgSoftClauseWeight,
                       mem->fixedprecShift);
    mem->weights[clauseIdx] = UFIXEDPREC(0);
  }
}

void puw_initWeights(const BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem) {
  const unsigned_fixedprec_t hWeightInit = fixedprec_uto(cfg->hWeightInit, mem->fixedprecShift);

  for (BouMS_uint_t clauseIdx = 0; clauseIdx < formula->numHardClauses; ++clauseIdx) {
    mem->weights[clauseIdx] = hWeightInit;
  }

  for (BouMS_uint_t clauseIdx = formula->numHardClauses; clauseIdx < formula->numClauses; ++clauseIdx) {
    mem->weights[clauseIdx] = UFIXEDPREC(0);
  }
}

void p_updateWeights(const BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem) {
  increaseHardWeights(formula, cfg, mem);
  if (mem->numFalsifiedHardClauses == 0) {
    if (cfg->isWeighted) {
      w_increaseSoftWeights(formula, cfg, mem);
    } else {
      uw_increaseSoftWeights(formula, cfg, mem);
    }
  }
}

static void solveHard(BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem, const BouMS_memoryReq_t* memReq,
                      const BouMS_result_t* result, const bool* initModel, bool forceRandom, bool* reuseDecimation,
                      const bool* stop) {
  INIT_DURATION_MEAS();
  START_DURATION_MEAS();

  // basically BouMS_solve
  for (BouMS_uint_t try = 0; try < cfg->maxTries && !done(mem, result) && !*stop; ++try) {
    initVars(formula, result, initModel, forceRandom);

    BouMS_decimation(formula, mem->decimationMem, &memReq->detailedDecimationMemReq, cfg->bmsSize, *reuseDecimation,
                     stop);
    *reuseDecimation = true;

    puw_initWeights(formula, cfg, mem);

    BouMS_uint_t cost;
    initAlgo(formula, mem, &cost);

    for (BouMS_uint_t flip = 0; flip < cfg->maxFlips && !done(mem, result) && !*stop; ++flip) {
      BouMS_wcnf_variable_t* variableToFlip = selectVariable(formula, cfg, mem);
      flipVariable(variableToFlip, formula, mem, &cost);
    }
  }

  // NOLINTNEXTLINE (bugprone-branch-clone)
  if (done(mem, result)) {
    LOG_TRACE_WITH_DURATION("c successfully solved hard clauses");
  } else {
    LOG_TRACE_WITH_DURATION("c couldn't solve hard clauses");
  }
}

static void increaseHardWeights(const BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem) {
  LOG_TRACE("c increasing " BOUMS_UINT_FORMAT " hard clause weights\n", mem->numFalsifiedHardClauses);

  const unsigned_fixedprec_t increase = fixedprec_uto(cfg->hWeightInc, mem->fixedprecShift);

  for (BouMS_uint_t i = 0; i < mem->numFalsifiedHardClauses; ++i) {
    // for each falsified hard clause
    const BouMS_uint_t clauseIdx = mem->falsifiedHardClauses[i];
    const BouMS_wcnf_clause_t* const clause = formula->clauses + clauseIdx;

    // increment the weight
    mem->weights[clauseIdx] = fixedprec_uadd(mem->weights[clauseIdx], increase);

    for (BouMS_uint_t litIdx = 0; litIdx < clause->numLiterals; ++litIdx) {
      // for all variables of the clause
      const BouMS_wcnf_literal_t* const literal = clause->literals + litIdx;
      // increase the score by the increment
      const BouMS_uint_t varIdx = BouMS_wcnf_var(literal);
      mem->scores[varIdx] = fixedprec_sadd(mem->scores[varIdx], fixedprec_utos(increase));
      addDecreasingVar(varIdx, mem);
    }
  }
}

static void w_increaseSoftWeights(const BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem) {
  LOG_TRACE("c increasing " BOUMS_UINT_FORMAT " soft clause weights\n", mem->numSoftClauses);

  // for each soft clause ...
  for (BouMS_uint_t clauseIdx = formula->numHardClauses; clauseIdx < formula->numClauses; ++clauseIdx) {
    const BouMS_wcnf_clause_t* const clause = formula->clauses + clauseIdx;
    assert(!BouMS_wcnf_isClauseHard(clause));

    const unsigned_fixedprec_t sWeightInc = fixedprec_uto(cfg->sWeightInc, mem->fixedprecShift);
    const unsigned_fixedprec_t increase = fixedprec_umul(sWeightInc, mem->tunedWeights[clauseIdx], mem->fixedprecShift);
    // LOG_TRACE("c increasing weight of soft clause " BouMS_UINT_FORMAT " by " BouMS_UINT_FORMAT "\n", clauseIdx,
    //           fixedprec_ufrom(mem->tunedWeights[clauseIdx], mem->fixedprecShift));

    // ... increase the weight
    mem->weights[clauseIdx] = fixedprec_uadd(mem->weights[clauseIdx], increase);

    updateScores(clauseIdx, clause, fixedprec_utos(increase), mem);
  }
}

static void uw_increaseSoftWeights(const BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem) {
  LOG_TRACE("c increasing " BOUMS_UINT_FORMAT " soft clause weights\n", mem->numSoftClauses);

  const unsigned_fixedprec_t sWeightInc = fixedprec_uto(cfg->sWeightInc, mem->fixedprecShift);

  // for each soft clause ...
  for (BouMS_uint_t clauseIdx = formula->numHardClauses; clauseIdx < formula->numClauses; ++clauseIdx) {
    const BouMS_wcnf_clause_t* const clause = formula->clauses + clauseIdx;
    assert(!BouMS_wcnf_isClauseHard(clause));

    const unsigned_fixedprec_t increase = sWeightInc;

    // ... increase the weight
    mem->weights[clauseIdx] = fixedprec_uadd(mem->weights[clauseIdx], increase);

    updateScores(clauseIdx, clause, fixedprec_utos(increase), mem);
  }
}
