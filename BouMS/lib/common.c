/**
 * @file common.c
 * @author Ole Lübke (ole.luebke@tuhh.de)
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "private/common.h"

#include <assert.h>
#include <stdbool.h>

#include "BouMS/BouMS.h"
#include "BouMS/common.h"
#include "BouMS/wcnf.h"
#include "private/fixedprec.h"
#include "private/logging.h"
#include "private/nonpartial.h"
#include "private/partial.h"
#include "private/rng.h"
#include "private/types.h"

/**
 * @brief
 *
 * @param numSatLiterals
 * @param varIdx
 * @param clause
 * @param clauseIdx
 * @param clauseWeight
 * @param cost
 * @param mem
 */
static void onLiteralSatisfied(BouMS_uint_t numSatLiterals, BouMS_uint_t varIdx, const BouMS_wcnf_clause_t* clause,
                               BouMS_uint_t clauseIdx, unsigned_fixedprec_t clauseWeight, BouMS_uint_t* cost,
                               BouMS_memory_t* mem);

/**
 * @brief
 *
 * @param numSatLiterals
 * @param formula
 * @param clause
 * @param clauseIdx
 * @param clauseWeight
 * @param cost
 * @param mem
 */
static void onLiteralFalsified(BouMS_uint_t numSatLiterals, const BouMS_wcnf_t* formula,
                               const BouMS_wcnf_clause_t* clause, BouMS_uint_t clauseIdx,
                               unsigned_fixedprec_t clauseWeight, BouMS_uint_t* cost, BouMS_memory_t* mem);

void initVars(const BouMS_wcnf_t* formula, const BouMS_result_t* result, const bool* initModel, bool forceRandom) {
  const BouMS_uint_t nVars = formula->numVariables;
  BouMS_wcnf_variable_t* const vars = formula->variables;

  if (!forceRandom && result->status != BOUMS_INFEASIBLE) {
    useAssignment(nVars, result->assignment, vars);  // best known assignment
  } else {
    if (!forceRandom && initModel) {
      useAssignment(nVars, initModel, vars);  // user-provided assignment
    } else {
      for (BouMS_uint_t varIdx = 0; varIdx < nVars; ++varIdx) {  // random assignment
        formula->variables[varIdx].value = flipCoin();
      }
    }
  }
}

void initAlgo(const BouMS_wcnf_t* formula, BouMS_memory_t* mem, BouMS_uint_t* cost) {
  // zero all memory that we write to (by adding) later
  *cost = 0;
  mem->numDecreasingVars = 0;
  mem->numFalsifiedHardClauses = 0;
  mem->numFalsifiedSoftClauses = 0;

  for (BouMS_uint_t varIdx = 0; varIdx < formula->numVariables; ++varIdx) {
    mem->scores[varIdx] = SFIXEDPREC(0);
    mem->decreasingVarsIdx[varIdx] = BOUMS_UINT_MAX;
    // mem->varFlipCounts[varIdx] = 0; see initMemory
  }

  for (BouMS_uint_t clauseIdx = 0; clauseIdx < formula->numClauses; ++clauseIdx) {
    mem->falsifiedHardClausesIdx[clauseIdx] = BOUMS_UINT_MAX;
    mem->falsifiedSoftClausesIdx[clauseIdx] = BOUMS_UINT_MAX;

    const BouMS_wcnf_clause_t* const clause = formula->clauses + clauseIdx;
    const unsigned_fixedprec_t clauseWeight = mem->weights[clauseIdx];

    // count satisfied literals per clause
    BouMS_uint_t* numSatLiterals = mem->numSatLiterals + clauseIdx;
    *numSatLiterals = 0;
    for (BouMS_uint_t litIdx = 0; litIdx < clause->numLiterals; ++litIdx) {
      const BouMS_wcnf_literal_t* const literal = clause->literals + litIdx;

      if (BouMS_wcnf_isLiteralSatisfied(formula, literal)) {
        *numSatLiterals += 1;

        if (*numSatLiterals == 1) {
          mem->satVars[clauseIdx] = BouMS_wcnf_var(literal);
        }
      }
    }

    const score_t clauseScore = fixedprec_utos(clauseWeight);
    if (*numSatLiterals == 1) {
      // penalize score of only satVar b/c if it is flipped, total cost increases
      mem->scores[mem->satVars[clauseIdx]] = fixedprec_ssub(mem->scores[mem->satVars[clauseIdx]], clauseScore);
    } else if (*numSatLiterals == 0) {
      // increase scores of all variables of the unsat clause
      for (BouMS_uint_t litIdx = 0; litIdx < clause->numLiterals; ++litIdx) {
        const BouMS_uint_t varIdx = BouMS_wcnf_var(clause->literals + litIdx);
        mem->scores[varIdx] = fixedprec_sadd(mem->scores[varIdx], clauseScore);
      }

      // add unsat clause to respective list
      if (BouMS_wcnf_isClauseHard(clause)) {
        addClause(clauseIdx, &mem->numFalsifiedHardClauses, mem->falsifiedHardClauses, mem->falsifiedHardClausesIdx);
      } else {
        addClause(clauseIdx, &mem->numFalsifiedSoftClauses, mem->falsifiedSoftClauses, mem->falsifiedSoftClausesIdx);
        // add weight of unsat clause to total cost
        *cost += clause->weight;
      }
    }
  }

  // now that initial scores for all variables are calculated, filter those that are decreasing
  for (BouMS_uint_t varIdx = 0; varIdx < formula->numVariables; ++varIdx) {
    addDecreasingVar(varIdx, mem);
  }
}

BouMS_wcnf_variable_t* selectVariable(const BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem) {
  {
    BouMS_uint_t bestVar;
    score_t bestScore = SFIXEDPREC(0);

    if (mem->numDecreasingVars > cfg->bmsSize) {
      // draw some decreasing vars randomly (w/ replacement) and select the one with highest score
      bestVar = mem->decreasingVars[randIdx(mem->numDecreasingVars)];
      bestScore = mem->scores[bestVar];

      for (BouMS_uint_t i = 1; i < cfg->bmsSize; ++i) {
        BouMS_uint_t curVar = mem->decreasingVars[randIdx(mem->numDecreasingVars)];
        score_t curScore = mem->scores[curVar];
        if (curScore.value > bestScore.value) {
          bestVar = curVar;
          bestScore = curScore;
        }
      }
      // LOG_TRACE("c picking var " BouMS_UINT_FORMAT " after BMS\n", bestVar);
      return formula->variables + bestVar;
    }

    if (mem->numDecreasingVars > 0) {
      // there are few decreasing vars, check all and select the one with the highest score
      bestVar = mem->decreasingVars[0];
      bestScore = mem->scores[bestVar];

      for (BouMS_uint_t i = 1; i < mem->numDecreasingVars; ++i) {
        BouMS_uint_t curVar = mem->decreasingVars[i];
        score_t curScore = mem->scores[curVar];
        if (curScore.value > bestScore.value) {
          bestVar = curVar;
          bestScore = curScore;
        }
      }

      // LOG_TRACE("c picking var " BouMS_UINT_FORMAT " b/c it's the best we have\n", bestVar);
      return formula->variables + bestVar;
    }
  }

  // if there are no decreasing vars, update weights
  LOG_TRACE("c reached local optimum\n");
  INIT_DURATION_MEAS();
  START_DURATION_MEAS();
  if (cfg->isPartial) {
    p_updateWeights(formula, cfg, mem);
  } else {
    np_updateWeights(formula, cfg, mem);
  }
  LOG_TRACE_WITH_DURATION("c updated weights");

  // select a clause randomly (prefer unsatisfied hard clauses)
  BouMS_uint_t numClausesToChooseFrom;
  const BouMS_uint_t* clausesToChooseFrom;
  if (mem->numFalsifiedHardClauses > 0) {
    // LOG_TRACE("c picking a var from an unsat hard clause\n");
    numClausesToChooseFrom = mem->numFalsifiedHardClauses;
    clausesToChooseFrom = mem->falsifiedHardClauses;
  } else {
    // LOG_TRACE("c picking a var from an unsat soft clause\n");
    numClausesToChooseFrom = mem->numFalsifiedSoftClauses;
    clausesToChooseFrom = mem->falsifiedSoftClauses;
  }

  const BouMS_wcnf_clause_t* const selectedClause =
      formula->clauses + clausesToChooseFrom[randIdx(numClausesToChooseFrom)];

  // select variable from selected clause
  assert(selectedClause->numLiterals > 0);
  BouMS_uint_t selectedVarIdx = BouMS_wcnf_var(selectedClause->literals);
  for (BouMS_uint_t litIdx = 1; litIdx < selectedClause->numLiterals; ++litIdx) {
    const BouMS_wcnf_literal_t* const literal = selectedClause->literals + litIdx;
    const BouMS_uint_t varIdx = BouMS_wcnf_var(literal);
#ifndef BOUMS_NOFFA
    if (mem->varFlipCounts[varIdx] < mem->varFlipCounts[selectedVarIdx])
#else
    if (mem->scores[varIdx].value > mem->scores[selectedVarIdx].value)
#endif
    {
      selectedVarIdx = varIdx;
    }
  }

  return formula->variables + selectedVarIdx;
}

void flipVariable(BouMS_wcnf_variable_t* variable, const BouMS_wcnf_t* formula, BouMS_memory_t* mem, BouMS_uint_t* cost) {
  // first things first, flip the variable
  variable->value = !variable->value;
  // and update its flip count
  const BouMS_uint_t varIdx = variable - formula->variables;

#ifndef BOUMS_NOFFA
  assert(mem->varFlipCounts[varIdx] < mem->varFlipCounts[varIdx] + 1);
  mem->varFlipCounts[varIdx] += 1;
#endif

  for (BouMS_uint_t litIdx = 0; litIdx < variable->numLiterals; ++litIdx) {
    // for each literal (clause) the flipped variable occurs in
    const BouMS_wcnf_literal_t* const literal = variable->literals[litIdx];
    const BouMS_wcnf_clause_t* const clause = literal->clause;
    const BouMS_uint_t clauseIdx = clause - formula->clauses;

    if (clauseIdx >= formula->numClauses) {
      continue;  // the clause is out of range, may happen when solving hard clauses only
    }

    const unsigned_fixedprec_t clauseWeight = mem->weights[clauseIdx];

    BouMS_uint_t* const numSatLiterals = mem->numSatLiterals + clauseIdx;

    // update the number of satisfied literals
    if (BouMS_wcnf_isLiteralSatisfied(formula, literal)) {
      *numSatLiterals += 1;
      onLiteralSatisfied(*numSatLiterals, varIdx, clause, clauseIdx, clauseWeight, cost, mem);
    } else {
      *numSatLiterals -= 1;
      onLiteralFalsified(*numSatLiterals, formula, clause, clauseIdx, clauseWeight, cost, mem);
    }
  }
}

void updateScores(BouMS_uint_t clauseIdx, const BouMS_wcnf_clause_t* clause, signed_fixedprec_t increase,
                  BouMS_memory_t* mem) {
  // if the clause is falsified ...
  if (mem->falsifiedSoftClausesIdx[clauseIdx] < BOUMS_UINT_MAX) {
    // ... increase the score of each variable in the clause
    for (BouMS_uint_t litIdx = 0; litIdx < clause->numLiterals; ++litIdx) {
      const BouMS_uint_t varIdx = BouMS_wcnf_var(clause->literals + litIdx);
      mem->scores[varIdx] = fixedprec_sadd(mem->scores[varIdx], increase);
      addDecreasingVar(varIdx, mem);
    }
  }
  // if the clause is satisfied by only one literal ...
  else if (mem->numSatLiterals[clauseIdx] == 1) {
    const BouMS_uint_t varIdx = mem->satVars[clauseIdx];
    // ... decrease the corresponding variable's score
    mem->scores[varIdx] = fixedprec_ssub(mem->scores[varIdx], increase);
    remDecreasingVar(varIdx, mem);
  }
}

static void onLiteralSatisfied(BouMS_uint_t numSatLiterals, BouMS_uint_t varIdx, const BouMS_wcnf_clause_t* clause,
                               BouMS_uint_t clauseIdx, unsigned_fixedprec_t clauseWeight, BouMS_uint_t* cost,
                               BouMS_memory_t* mem) {
  const score_t clauseScore = fixedprec_utos(clauseWeight);

  if (numSatLiterals == 1) {
    // if the clause is satisfied now
    // update the scores of all involved variables
    mem->scores[varIdx] = fixedprec_ssub(mem->scores[varIdx], clauseScore);
    for (BouMS_uint_t lIdx = 0; lIdx < clause->numLiterals; ++lIdx) {
      const BouMS_uint_t vIdx = BouMS_wcnf_var(clause->literals + lIdx);
      mem->scores[vIdx] = fixedprec_ssub(mem->scores[vIdx], clauseScore);
      remDecreasingVar(vIdx, mem);
    }

    // save the only satisfying variable
    mem->satVars[clauseIdx] = varIdx;

    // remove the clause from the falsified clauses
    if (BouMS_wcnf_isClauseHard(clause)) {
      remClause(clauseIdx, &mem->numFalsifiedHardClauses, mem->falsifiedHardClauses, mem->falsifiedHardClausesIdx);
    } else {
      // decrease total cost by the clause's weight
      assert(*cost >= *cost - clause->weight);
      *cost -= clause->weight;
      remClause(clauseIdx, &mem->numFalsifiedSoftClauses, mem->falsifiedSoftClauses, mem->falsifiedSoftClausesIdx);
    }
  } else if (numSatLiterals == 2) {
    // if the clause was satisfied already, but is now satisfied by 2 literals
    // update the score of the previously only satisfying variable
    const BouMS_uint_t vIdx = mem->satVars[clauseIdx];
    mem->scores[vIdx] = fixedprec_sadd(mem->scores[vIdx], clauseScore);
    addDecreasingVar(vIdx, mem);
  }
}

static void onLiteralFalsified(BouMS_uint_t numSatLiterals, const BouMS_wcnf_t* formula,
                               const BouMS_wcnf_clause_t* clause, BouMS_uint_t clauseIdx,
                               unsigned_fixedprec_t clauseWeight, BouMS_uint_t* cost, BouMS_memory_t* mem) {
  const score_t clauseScore = fixedprec_utos(clauseWeight);

  if (numSatLiterals == 1) {
    // if the clause is only satisfied by a single variable now
    // find that variable...
    BouMS_uint_t litIdx;
    for (litIdx = 0; litIdx < clause->numLiterals; ++litIdx) {
      const BouMS_wcnf_literal_t* const literal = clause->literals + litIdx;
      if (BouMS_wcnf_isLiteralSatisfied(formula, literal)) {
        // ...update its score..
        const BouMS_uint_t varIdx = BouMS_wcnf_var(literal);
        mem->scores[varIdx] = fixedprec_ssub(mem->scores[varIdx], clauseScore);
        remDecreasingVar(varIdx, mem);
        // ...and save it
        mem->satVars[clauseIdx] = varIdx;
        break;
      }
    }
    assert(litIdx < clause->numLiterals);
  } else if (numSatLiterals == 0) {
    // if the clause is falsified now
    // update the scores of all involved variables
    mem->scores[mem->satVars[clauseIdx]] = fixedprec_sadd(mem->scores[mem->satVars[clauseIdx]], clauseScore);
    for (BouMS_uint_t litIdx = 0; litIdx < clause->numLiterals; ++litIdx) {
      const BouMS_uint_t varIdx = BouMS_wcnf_var(clause->literals + litIdx);
      mem->scores[varIdx] = fixedprec_sadd(mem->scores[varIdx], clauseScore);
      addDecreasingVar(varIdx, mem);
    }

    // unsave the only satisfying variable
    mem->satVars[clauseIdx] = BOUMS_UINT_MAX;

    // add the clause to the falsified clauses
    if (BouMS_wcnf_isClauseHard(clause)) {
      addClause(clauseIdx, &mem->numFalsifiedHardClauses, mem->falsifiedHardClauses, mem->falsifiedHardClausesIdx);
    } else {
      // increase total costs by the clause's weight
      assert(*cost <= *cost + clause->weight);
      *cost += clause->weight;
      addClause(clauseIdx, &mem->numFalsifiedSoftClauses, mem->falsifiedSoftClauses, mem->falsifiedSoftClausesIdx);
    }
  }
}
