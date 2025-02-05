/**
 * @file nonpartial.c
 * @author Ole Lübke (ole.luebke@tuhh.de)
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "private/nonpartial.h"

#include <assert.h>
#include <stdbool.h>

#include "BouMS/BouMS.h"
#include "BouMS/common.h"
#include "BouMS/wcnf.h"
#include "private/common.h"
#include "private/fixedprec.h"
#include "private/logging.h"
#include "private/rng.h"
#include "private/types.h"

/**
 * @brief
 *
 * @param cfg
 * @param mem
 */
static void smoothSoftWeights(const BouMS_params_t* cfg, BouMS_memory_t* mem);

/**
 * @brief
 *
 * @param formula
 * @param cfg
 * @param mem
 */
static void increaseSoftWeights(const BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem);

/**
 * @param formula
 * @param mem
 */
static void resetHighWeightSoftClauses(const BouMS_wcnf_t* formula, BouMS_memory_t* mem);

// NOLINTBEGIN(readability-non-const-parameter)
void np_initVars(BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem, const BouMS_memoryReq_t* memReq,
                 BouMS_result_t* result, const bool* initModel, bool* reuseDecimation, bool forceRandomAssignment,
                 const bool* stop)
// NOLINTEND(readability-non-const-parameter)
{
  ((void)cfg);
  ((void)mem);
  ((void)memReq);
  ((void)reuseDecimation);
  ((void)stop);

  initVars(formula, result, initModel, forceRandomAssignment);
}

void npw_initWeights(const BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem,
                     unsigned_fixedprec_t avgSoftClauseWeight) {
  const unsigned_fixedprec_t sWeightBoundCoef = fixedprec_uto(cfg->sWeightBoundCoef, mem->fixedprecShift);
  const unsigned_fixedprec_t sWeightInc = fixedprec_uto(cfg->sWeightInc, mem->fixedprecShift);

  resetHighWeightSoftClauses(formula, mem);

  for (BouMS_uint_t i = 0; i < formula->numClauses; ++i) {
    mem->tunedWeights[i] = fixedprec_umul(sWeightBoundCoef,
                                          fixedprec_udiv(fixedprec_uto(formula->clauses[i].weight, mem->fixedprecShift),
                                                         avgSoftClauseWeight, mem->fixedprecShift),
                                          mem->fixedprecShift);
    mem->weights[i] = mem->tunedWeights[i];
    if (mem->weights[i].value > sWeightInc.value) {
      addClause(i, &mem->numHighWeightSoftClauses, mem->highWeightSoftClauses, mem->highWeightSoftClausesIdx);
    }
  }
}

void npuw_initWeights(const BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem) {
  const unsigned_fixedprec_t sWeightBoundCoef = fixedprec_uto(cfg->sWeightBoundCoef, mem->fixedprecShift);
  const unsigned_fixedprec_t sWeightInc = fixedprec_uto(cfg->sWeightInc, mem->fixedprecShift);

  resetHighWeightSoftClauses(formula, mem);

  for (BouMS_uint_t i = 0; i < formula->numClauses; ++i) {
    mem->tunedWeights[i] = sWeightBoundCoef;
    mem->weights[i] = sWeightBoundCoef;
    if (mem->weights[i].value > sWeightInc.value) {
      addClause(i, &mem->numHighWeightSoftClauses, mem->highWeightSoftClauses, mem->highWeightSoftClausesIdx);
    }
  }
}

void np_updateWeights(const BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem) {
  if (flipBiasedCoin(cfg->sSmoothProb)) {
    smoothSoftWeights(cfg, mem);
  } else {
    increaseSoftWeights(formula, cfg, mem);
  }
}

static void smoothSoftWeights(const BouMS_params_t* cfg, BouMS_memory_t* mem) {
  LOG_TRACE("c smoothing soft clause weights\n");

  const unsigned_fixedprec_t decrease = fixedprec_uto(cfg->sWeightInc, mem->fixedprecShift);

  // for each soft clause with increased weight ...
  for (long long highWeightClauseIdx = 0; highWeightClauseIdx < (long long)mem->numHighWeightSoftClauses;
       ++highWeightClauseIdx) {
    const BouMS_uint_t clauseIdx = mem->highWeightSoftClauses[highWeightClauseIdx];

    // ... if it is satisfied ...
    if (mem->falsifiedSoftClausesIdx[clauseIdx] == BOUMS_UINT_MAX) {
      // ... decrease the weight
      mem->weights[clauseIdx] = fixedprec_usub(mem->weights[clauseIdx], decrease);

      // if now its weight isn't high anymore ...
      if (mem->weights[clauseIdx].value < decrease.value) {
        // ... remove it from the respective list
        remClause(clauseIdx, &mem->numHighWeightSoftClauses, mem->highWeightSoftClauses, mem->highWeightSoftClausesIdx);
        highWeightClauseIdx -= 1;
      }

      // if the clause is only satisfied by a single literal ...
      if (mem->numSatLiterals[clauseIdx] == 1) {
        // ... increase the score of the literal's variable
        const BouMS_uint_t varIdx = mem->satVars[clauseIdx];
        mem->scores[varIdx] = fixedprec_sadd(mem->scores[varIdx], fixedprec_utos(decrease));
        addDecreasingVar(varIdx, mem);
      }
    }
  }
}

static void increaseSoftWeights(const BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem) {
  LOG_TRACE("c increasing " BOUMS_UINT_FORMAT " soft clause weights\n", mem->numFalsifiedSoftClauses);

  const unsigned_fixedprec_t sWeightInc = fixedprec_uto(cfg->sWeightInc, mem->fixedprecShift);
  const unsigned_fixedprec_t sWeightBoundOffset = fixedprec_uto(cfg->sWeightBoundOffset, mem->fixedprecShift);

  // for each falsified soft clause ...
  for (BouMS_uint_t i = 0; i < mem->numFalsifiedSoftClauses; ++i) {
    const BouMS_uint_t clauseIdx = mem->falsifiedSoftClauses[i];
    const BouMS_wcnf_clause_t* const clause = formula->clauses + clauseIdx;
    assert(!BouMS_wcnf_isClauseHard(clause));

    const unsigned_fixedprec_t increase = sWeightInc;

    // ... if the upper weight bound hasn't been reached ...
    if (mem->weights[clauseIdx].value >= fixedprec_uadd(mem->tunedWeights[clauseIdx], sWeightBoundOffset).value) {
      continue;
    }

    // ... increase the weight
    mem->weights[clauseIdx] = fixedprec_uadd(mem->weights[clauseIdx], increase);
    // if the clause is not yet in the high weight soft clauses list ...
    if (mem->highWeightSoftClausesIdx[clauseIdx] == BOUMS_UINT_MAX) {
      // ... add it
      addClause(clauseIdx, &mem->numHighWeightSoftClauses, mem->highWeightSoftClauses, mem->highWeightSoftClausesIdx);
    }

    updateScores(clauseIdx, clause, fixedprec_utos(increase), mem);
  }
}

static void resetHighWeightSoftClauses(const BouMS_wcnf_t* formula, BouMS_memory_t* mem) {
  mem->numHighWeightSoftClauses = 0;
  for (BouMS_uint_t clauseIdx = 0; clauseIdx < formula->numClauses; ++clauseIdx) {
    mem->highWeightSoftClausesIdx[clauseIdx] = BOUMS_UINT_MAX;
  }
}
