/**
 * @file decimation.c
 * @author Ole Lübke (ole.luebke@tuhh.de)
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "BouMS/common.h"
#include "BouMS/preprocessing.h"
#include "BouMS/wcnf.h"
#include "private/align.h"
#include "private/rng.h"

typedef struct {
  BouMS_uint_t numUnassignedVars;
  BouMS_wcnf_variable_t** unassignedVars;
  BouMS_uint_t* unassignedVarsIdx;

  BouMS_uint_t numHardUnitClauses;
  const BouMS_wcnf_literal_t** hardUnitClauses;
  BouMS_uint_t* hardUnitClausesIdx;

  BouMS_uint_t numSoftUnitClauses;
  const BouMS_wcnf_literal_t** softUnitClauses;
  BouMS_uint_t* softUnitClausesIdx;

  BouMS_uint_t* prevDeciStep;
  bool* clauseRemoved;
  BouMS_uint_t* numLiterals;
} memory_t;

/**
 * @brief Calculates detailed memory requirements
 *
 * @param formula
 * @param out
 */
static void calcMemoryReq(const BouMS_wcnf_t* formula, BouMS_decimation_memoryReq_t* out);

/**
 * @brief Calculates where pointers in out should point to
 *
 * @param mem
 * @param memReq
 * @param out
 */
static void initMemory(void* mem, const BouMS_decimation_memoryReq_t* memReq, memory_t* out);

/**
 * @brief
 *
 * @param formula
 * @param mem
 * @param reusePrevious
 */
static void initAlgorithm(const BouMS_wcnf_t* formula, memory_t* mem, bool reusePrevious);

/**
 * @brief
 *
 * @param mem
 * @return BouMS_wcnf_variable_t*
 */
static BouMS_wcnf_variable_t* pickVariable(const memory_t* mem, const BouMS_wcnf_t* formula, BouMS_uint_t bmsSize);

/**
 * @brief
 *
 * @param mem
 * @param bmsSize
 * @param hard
 * @return const BouMS_wcnf_literal_t*
 */
static const BouMS_wcnf_literal_t* pickUnitClause(const memory_t* mem, BouMS_uint_t bmsSize, bool hard);

/**
 * @brief
 *
 * @param formula
 * @param mem
 * @param variable
 */
static void simplify(const BouMS_wcnf_t* formula, const BouMS_wcnf_variable_t* variable, memory_t* mem);

/**
 * @brief Returns the first satisfied literal in the given clause
 *
 * @param clause
 * @param mem
 * @return const BouMS_wcnf_literal_t*
 */
static const BouMS_wcnf_literal_t* findUnitLiteral(const BouMS_wcnf_clause_t* clause, const memory_t* mem);

/**
 * @brief
 *
 * @param varIdx
 * @param formula
 * @param mem
 */
static void remUnassignedVar(BouMS_uint_t varIdx, const BouMS_wcnf_t* formula, memory_t* mem);

/**
 * @brief
 *
 * @param literal
 * @param formula
 * @param mem
 */
static void addUnitClause(const BouMS_wcnf_literal_t* literal, const BouMS_wcnf_t* formula, memory_t* mem);

/**
 * @brief
 * @param literal
 * @param formula
 * @param mem
 */
static void remUnitClause(const BouMS_wcnf_literal_t* literal, const BouMS_wcnf_t* formula, memory_t* mem);

BouMS_uint_t BouMS_decimation_calcMemoryRequirements(const BouMS_wcnf_t* formula,
                                                     BouMS_decimation_memoryReq_t* memReq) {
  calcMemoryReq(formula, memReq);
  return memReq->unassignedVarsMemReq + memReq->unassignedVarsIdxMemReq + memReq->hardUnitClausesMemReq +
         memReq->hardUnitClausesIdxMemReq + memReq->softUnitClausesMemReq + memReq->softUnitClausesIdxMemReq +
         memReq->prevDeciStepMemReq + memReq->clauseRemovedMemReq + memReq->numLiteralsMemReq;
}

void BouMS_decimation(const BouMS_wcnf_t* formula, void* memory, const BouMS_decimation_memoryReq_t* memReq,
                      BouMS_uint_t bmsSize, bool reusePrevious, const bool* stop) {
  memory_t mem;
  initMemory(memory, memReq, &mem);

  initAlgorithm(formula, &mem, reusePrevious);

  for (BouMS_uint_t iteration = 0; iteration < formula->numVariables && !*stop; ++iteration) {
    BouMS_wcnf_variable_t* pickedVar = NULL;

    const bool hardUnitClausesExist = mem.numHardUnitClauses > 0;
    const bool softUnitClausesExist = mem.numSoftUnitClauses > 0;
    if (hardUnitClausesExist || softUnitClausesExist) {
      const BouMS_wcnf_literal_t* const pickedLit = pickUnitClause(&mem, bmsSize, hardUnitClausesExist);
      pickedVar = formula->variables + BouMS_wcnf_var(pickedLit);
      pickedVar->value = !BouMS_wcnf_sign(pickedLit);
    } else {
      pickedVar = pickVariable(&mem, formula, bmsSize);
    }

    const BouMS_uint_t pickedVarIdx = pickedVar - formula->variables;
    remUnassignedVar(pickedVarIdx, formula, &mem);
    mem.prevDeciStep[pickedVarIdx] = iteration;
    simplify(formula, pickedVar, &mem);
  }
}

static void calcMemoryReq(const BouMS_wcnf_t* formula, BouMS_decimation_memoryReq_t* out) {
  out->unassignedVarsMemReq = align(formula->numVariables * sizeof(BouMS_wcnf_variable_t*));
  out->unassignedVarsIdxMemReq = align(formula->numVariables * sizeof(BouMS_uint_t));

  out->hardUnitClausesMemReq = align(formula->numClauses * sizeof(BouMS_wcnf_literal_t*));
  out->hardUnitClausesIdxMemReq = align(formula->numClauses * sizeof(BouMS_uint_t));

  out->softUnitClausesMemReq = align(formula->numClauses * sizeof(BouMS_wcnf_literal_t*));
  out->softUnitClausesIdxMemReq = align(formula->numClauses * sizeof(BouMS_uint_t));

  out->prevDeciStepMemReq = align(formula->numVariables * sizeof(BouMS_uint_t));
  out->clauseRemovedMemReq = align(formula->numClauses * sizeof(bool));
  out->numLiteralsMemReq = align(formula->numClauses * sizeof(BouMS_uint_t));
}

static void initMemory(void* mem, const BouMS_decimation_memoryReq_t* memReq, memory_t* out) {
  out->numUnassignedVars = 0;
  out->unassignedVars = (BouMS_wcnf_variable_t**)mem;
  out->unassignedVarsIdx = (BouMS_uint_t*)((uint8_t*)out->unassignedVars + memReq->unassignedVarsMemReq);

  out->numHardUnitClauses = 0;
  out->hardUnitClauses =
      (const BouMS_wcnf_literal_t**)((uint8_t*)out->unassignedVarsIdx + memReq->unassignedVarsIdxMemReq);
  out->hardUnitClausesIdx = (BouMS_uint_t*)((uint8_t*)out->hardUnitClauses + memReq->hardUnitClausesMemReq);

  out->numSoftUnitClauses = 0;
  out->softUnitClauses =
      (const BouMS_wcnf_literal_t**)((uint8_t*)out->hardUnitClausesIdx + memReq->hardUnitClausesIdxMemReq);
  out->softUnitClausesIdx = (BouMS_uint_t*)((uint8_t*)out->softUnitClauses + memReq->softUnitClausesMemReq);

  out->prevDeciStep = (BouMS_uint_t*)((uint8_t*)out->softUnitClausesIdx + memReq->softUnitClausesIdxMemReq);
  out->clauseRemoved = (bool*)((uint8_t*)out->prevDeciStep + memReq->prevDeciStepMemReq);
  out->numLiterals = (BouMS_uint_t*)((uint8_t*)out->clauseRemoved + memReq->clauseRemovedMemReq);
}

static void initAlgorithm(const BouMS_wcnf_t* formula, memory_t* mem, bool reusePrevious) {
  mem->numUnassignedVars = formula->numVariables;
  for (BouMS_uint_t varIdx = 0; varIdx < formula->numVariables; ++varIdx) {
    BouMS_wcnf_variable_t* const variable = formula->variables + varIdx;
    mem->unassignedVars[varIdx] = variable;
    mem->unassignedVarsIdx[varIdx] = varIdx;
    if (!reusePrevious) {
      mem->prevDeciStep[varIdx] = BOUMS_UINT_MAX;
    }
  }

  for (BouMS_uint_t clauseIdx = 0; clauseIdx < formula->numClauses; ++clauseIdx) {
    BouMS_wcnf_clause_t* const clause = formula->clauses + clauseIdx;

    mem->hardUnitClausesIdx[clauseIdx] = BOUMS_UINT_MAX;
    mem->softUnitClausesIdx[clauseIdx] = BOUMS_UINT_MAX;
    mem->clauseRemoved[clauseIdx] = false;
    mem->numLiterals[clauseIdx] = clause->numLiterals;

    if (clause->numLiterals == 1) {
      addUnitClause(clause->literals, formula, mem);
    }
  }
}

static BouMS_wcnf_variable_t* pickVariable(const memory_t* mem, const BouMS_wcnf_t* formula, BouMS_uint_t bmsSize) {
  BouMS_wcnf_variable_t* variable = mem->unassignedVars[randIdx(mem->numUnassignedVars)];
  BouMS_uint_t variablePrevDeciStep = mem->prevDeciStep[variable - formula->variables];
  for (BouMS_uint_t i = 1; i < bmsSize; ++i) {
    BouMS_wcnf_variable_t* candidateVariable = mem->unassignedVars[randIdx(mem->numUnassignedVars)];
    BouMS_uint_t candidatePrevDeciStep = mem->prevDeciStep[candidateVariable - formula->variables];
    if (candidatePrevDeciStep > variablePrevDeciStep) {
      variablePrevDeciStep = candidatePrevDeciStep;
      variable = candidateVariable;
    }
  }

  return variable;
}

static const BouMS_wcnf_literal_t* pickUnitClause(const memory_t* mem, BouMS_uint_t bmsSize, bool hard) {
  BouMS_uint_t numUnitClauses;
  const BouMS_wcnf_literal_t** unitClauses;

  if (hard) {
    numUnitClauses = mem->numHardUnitClauses;
    unitClauses = mem->hardUnitClauses;
  } else {
    numUnitClauses = mem->numSoftUnitClauses;
    unitClauses = mem->softUnitClauses;
  }

  const BouMS_wcnf_literal_t* literal = unitClauses[randIdx(numUnitClauses)];
  BouMS_uint_t litVarPrevDeciStep = mem->prevDeciStep[BouMS_wcnf_var(literal)];
  for (BouMS_uint_t i = 1; i < bmsSize; ++i) {
    const BouMS_wcnf_literal_t* candidateLiteral = unitClauses[randIdx(numUnitClauses)];
    const BouMS_uint_t candiatePrevDeciStep = mem->prevDeciStep[BouMS_wcnf_var(candidateLiteral)];
    if (candiatePrevDeciStep > litVarPrevDeciStep) {
      litVarPrevDeciStep = candiatePrevDeciStep;
      literal = candidateLiteral;
    }
  }

  return literal;
}

static void simplify(const BouMS_wcnf_t* formula, const BouMS_wcnf_variable_t* variable, memory_t* mem) {
  for (BouMS_uint_t litIdx = 0; litIdx < variable->numLiterals; ++litIdx) {
    const BouMS_wcnf_literal_t* const literal = variable->literals[litIdx];
    const BouMS_wcnf_clause_t* const clause = literal->clause;
    const BouMS_uint_t clauseIdx = clause - formula->clauses;

    if (clauseIdx >= formula->numClauses) {
      continue;  // the clause is out of range, may happen when solving hard clauses only
    }

    if (!mem->clauseRemoved[clauseIdx]) {
      if (BouMS_wcnf_isLiteralSatisfied(formula, literal)) {
        mem->clauseRemoved[clauseIdx] = true;

        if (mem->numLiterals[clauseIdx] == 1) {
          remUnitClause(literal, formula, mem);
        }
      } else {
        mem->numLiterals[clauseIdx] -= 1;

        if (mem->numLiterals[clauseIdx] == 1) {
          const BouMS_wcnf_literal_t* unitLiteral = findUnitLiteral(clause, mem);
          // NULL check b/c finding the lit. with unass. var. may fail when var. occurs in clause multiple times
          if (unitLiteral) {
            addUnitClause(unitLiteral, formula, mem);
          } else {
            mem->clauseRemoved[clauseIdx] = true;
          }
        } else if (mem->numLiterals[clauseIdx] == 0) {
          mem->clauseRemoved[clauseIdx] = true;
          remUnitClause(literal, formula, mem);
        }
      }
    }
  }
}

static const BouMS_wcnf_literal_t* findUnitLiteral(const BouMS_wcnf_clause_t* clause, const memory_t* mem) {
  for (BouMS_uint_t litIdx = 0; litIdx < clause->numLiterals; ++litIdx) {
    if (mem->unassignedVarsIdx[BouMS_wcnf_var(clause->literals + litIdx)] < BOUMS_UINT_MAX) {
      return clause->literals + litIdx;
    }
  }
  return NULL;
}

static void remUnassignedVar(BouMS_uint_t varIdx, const BouMS_wcnf_t* formula, memory_t* mem) {
  const BouMS_uint_t replaceIdx = mem->unassignedVarsIdx[varIdx];

  mem->unassignedVarsIdx[varIdx] = BOUMS_UINT_MAX;
  mem->numUnassignedVars -= 1;

  if (mem->numUnassignedVars > 0 && replaceIdx < mem->numUnassignedVars) {
    BouMS_wcnf_variable_t* const moveVar = mem->unassignedVars[mem->numUnassignedVars];
    mem->unassignedVars[replaceIdx] = moveVar;
    mem->unassignedVarsIdx[moveVar - formula->variables] = replaceIdx;
  }
}

static void addUnitClause(const BouMS_wcnf_literal_t* literal, const BouMS_wcnf_t* formula, memory_t* mem) {
  BouMS_uint_t* numUnitClauses;
  const BouMS_wcnf_literal_t** unitClauses;
  BouMS_uint_t* unitClausesIdx;

  if (BouMS_wcnf_isClauseHard(literal->clause)) {
    numUnitClauses = &mem->numHardUnitClauses;
    unitClauses = mem->hardUnitClauses;
    unitClausesIdx = mem->hardUnitClausesIdx;
  } else {
    numUnitClauses = &mem->numSoftUnitClauses;
    unitClauses = mem->softUnitClauses;
    unitClausesIdx = mem->softUnitClausesIdx;
  }

  unitClauses[*numUnitClauses] = literal;
  unitClausesIdx[literal->clause - formula->clauses] = *numUnitClauses;
  *numUnitClauses += 1;
}

static void remUnitClause(const BouMS_wcnf_literal_t* literal, const BouMS_wcnf_t* formula, memory_t* mem) {
  const BouMS_uint_t clauseIdx = literal->clause - formula->clauses;

  BouMS_uint_t* numUnitClauses;
  const BouMS_wcnf_literal_t** unitClauses;
  BouMS_uint_t* unitClausesIdx;

  if (BouMS_wcnf_isClauseHard(literal->clause)) {
    numUnitClauses = &mem->numHardUnitClauses;
    unitClauses = mem->hardUnitClauses;
    unitClausesIdx = mem->hardUnitClausesIdx;
  } else {
    numUnitClauses = &mem->numSoftUnitClauses;
    unitClauses = mem->softUnitClauses;
    unitClausesIdx = mem->softUnitClausesIdx;
  }

  const BouMS_uint_t replaceIdx = unitClausesIdx[clauseIdx];

  unitClausesIdx[clauseIdx] = BOUMS_UINT_MAX;
  *numUnitClauses -= 1;

  if (*numUnitClauses > 0 && replaceIdx < *numUnitClauses) {
    const BouMS_wcnf_literal_t* const moveLit = unitClauses[*numUnitClauses];
    unitClauses[replaceIdx] = moveLit;
    unitClausesIdx[moveLit->clause - formula->clauses] = replaceIdx;
  }
}
