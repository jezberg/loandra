/**
 * @file wcnf.h
 * @author Ole Lübke (ole.luebke@tuhh.de)
 * @brief Types for representing a MaxSAT instance in WCNF and associated functions
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef BOUMS_FORMULAS_WCNF_H
#define BOUMS_FORMULAS_WCNF_H

#include <assert.h>
#include <stdbool.h>

#include "BouMS/common.h"

#ifdef __cplusplus
extern "C" {
#endif

#define BOUMS_HARD_CLAUSE_WEIGHT BOUMS_UINT_MAX

struct BouMS_wcnf_literal_s;
struct BouMS_wcnf_clause_s;

/**
 * @brief A variable in a BouMS_wcnf_t struct
 *
 */
typedef struct {
  bool value;                                    ///> The truth value of the variable
  BouMS_uint_t numLiterals;                      ///> The number of literals the variable occurs in
  const struct BouMS_wcnf_literal_s** literals;  ///> Pointers to the literals the variable occurs in
} BouMS_wcnf_variable_t;

/** The bottom bit is the sign, the remaining bits are the 0-based variable index */
typedef BouMS_uint_t BouMS_literal_t;

/**
 * @brief Helper function to construct the `lit` field of an BouMS_wcnf_literal_t
 *
 * @param varIdx The index of the literal's variable
 * @param sign Whether the literal is negated
 * @return Cf. BouMS_literal_t
 */
static inline BouMS_literal_t BouMS_mkLit(BouMS_uint_t varIdx, bool sign) {
  return (varIdx << 1) | sign;
}

/**
 * @brief A literal, i.e., a possibly negated variable, in a BouMS_wcnf_clause_t struct
 *
 */
typedef struct BouMS_wcnf_literal_s {
  BouMS_literal_t lit;                       ///> Cf. BouMS_literal_t
  const struct BouMS_wcnf_clause_s* clause;  ///> Pointer to the clause the literal is part of
} BouMS_wcnf_literal_t;

/**
 * @brief Get the variable index of a literal
 *
 * @param lit
 */
static inline BouMS_uint_t BouMS_wcnf_var(const BouMS_wcnf_literal_t* lit) {
  return lit->lit >> 1;
}

/**
 * @brief Get the sign of a literal
 *
 * @param lit
 */
static inline bool BouMS_wcnf_sign(const BouMS_wcnf_literal_t* lit) {
  return lit->lit & 1;
}

/**
 * @brief A weighted CNF clause in a BouMS_wcnf_t struct
 *
 */
typedef struct BouMS_wcnf_clause_s {
  BouMS_uint_t weight;             ///< The weight of the clause, 0 for hard clauses, >= 1 for soft clauses
  BouMS_uint_t numLiterals;        ///< The number of literals in the clause
  BouMS_wcnf_literal_t* literals;  ///< Array which hold the literals
} BouMS_wcnf_clause_t;

/**
 * @brief Checks whether a clause is hard
 *
 * @param clause The clause to check
 * @return true When the clause is hard
 * @return false Otherwise
 */
static inline bool BouMS_wcnf_isClauseHard(const BouMS_wcnf_clause_t* clause) {
  return clause->weight == BOUMS_HARD_CLAUSE_WEIGHT;
}

/**
 * @brief Checks whether a clause is empty
 *
 * @param clause The clause to check
 * @return true When the clause is empty
 * @return false Otherwise
 */
static inline bool BouMS_wcnf_isClauseEmpty(const BouMS_wcnf_clause_t* clause) {
  return clause->numLiterals == 0;
}

/**
 * @brief A formula in WCNF, i.e., a conjunction of weighted clauses
 *
 */
typedef struct {
  BouMS_uint_t numVariables;         ///< The total number of variables of the formula
  BouMS_uint_t numClauses;           ///< The total number of clauses of the formula
  BouMS_uint_t numHardClauses;       ///< The number of hard clauses in the formula
  BouMS_wcnf_variable_t* variables;  ///< Array which hold the variables
  BouMS_wcnf_clause_t* clauses;      ///< Array which holds the clauses
} BouMS_wcnf_t;

/**
 * @brief Checks whether a literal is satisfied
 *
 * @param formula The associated formula
 * @param literal The literal to check
 * @return true When the literal is satisfied
 * @return false Otherwise
 */
static inline bool BouMS_wcnf_isLiteralSatisfied(const BouMS_wcnf_t* formula, const BouMS_wcnf_literal_t* literal) {
  return formula->variables[BouMS_wcnf_var(literal)].value != BouMS_wcnf_sign(literal);
}

/**
 * @brief Removes the literal at the given index
 *
 * @param formula The formula the clause is part of
 * @param clause The clause to remove the literal from
 * @param litIdx The index of the literal in the clause
 */
static inline void BouMS_wcnf_removeLiteralAt(BouMS_wcnf_t* formula, BouMS_wcnf_clause_t* clause, BouMS_uint_t litIdx) {
  assert(litIdx < clause->numLiterals);
  assert((BouMS_uint_t)(clause - formula->clauses) < formula->numClauses);

  const BouMS_wcnf_literal_t* const litToRemove = clause->literals + litIdx;
  const BouMS_wcnf_literal_t* const litToMove = clause->literals + clause->numLiterals - 1;

  {  // remove pointer to literal from var
    const BouMS_uint_t varIdx = BouMS_wcnf_var(litToRemove);
    BouMS_wcnf_variable_t* const var = formula->variables + varIdx;
    {
      BouMS_uint_t varLitIdx;
      for (varLitIdx = 0; varLitIdx < var->numLiterals; ++varLitIdx) {
        const BouMS_wcnf_literal_t* const varLit = var->literals[varLitIdx];
        if (varLit == litToRemove) {
          var->literals[varLitIdx] = var->literals[var->numLiterals - 1];
          break;
        }
      }
      assert(varLitIdx < var->numLiterals);
      var->numLiterals -= 1;
    }
  }

  if (litToRemove != litToMove) {  // update pointer to moved literal from var
    const BouMS_uint_t varIdx = BouMS_wcnf_var(litToMove);
    BouMS_wcnf_variable_t* const var = formula->variables + varIdx;
    {
      BouMS_uint_t varLitIdx;
      for (varLitIdx = 0; varLitIdx < var->numLiterals; ++varLitIdx) {
        const BouMS_wcnf_literal_t* const varLit = var->literals[varLitIdx];
        if (varLit == litToMove) {
          var->literals[varLitIdx] = litToRemove;
          break;
        }
      }
      assert(varLitIdx < var->numLiterals);
    }
  }

  clause->literals[litIdx] = *litToMove;
  clause->numLiterals -= 1;
}

/**
 * @brief Removes duplicate literals from a clause and checks whether it is a tautology
 *
 * @param formula The formula the clause is part of
 * @param clause The clause to clean
 * @return true When the clause is a tautology
 * @return false Otherwise
 */
static inline bool BouMS_wcnf_cleanClause(BouMS_wcnf_t* formula, BouMS_wcnf_clause_t* clause) {
  assert((BouMS_uint_t)(clause - formula->clauses) < formula->numClauses);

  for (BouMS_uint_t litIdx = 0; litIdx < clause->numLiterals; ++litIdx) {
    const BouMS_wcnf_literal_t* lit = clause->literals + litIdx;
    for (BouMS_uint_t cmpIdx = litIdx + 1; cmpIdx < clause->numLiterals; ++cmpIdx) {
      const BouMS_wcnf_literal_t* cmp = clause->literals + cmpIdx;
      if (BouMS_wcnf_var(lit) == BouMS_wcnf_var(cmp)) {
        if (BouMS_wcnf_sign(lit) != BouMS_wcnf_sign(cmp)) {
          return true;
        }
        BouMS_wcnf_removeLiteralAt(formula, clause, cmpIdx--);
      }
    }
  }
  return false;
}

/**
 * @brief Count the number of satisfied literals in a clause
 *
 * @param formula The associated formula
 * @param clause The clause to check
 * @return BouMS_uint_t The number of satisfied literals in the clause
 */
static inline BouMS_uint_t BouMS_wcnf_countSatLiterals(const BouMS_wcnf_t* formula, const BouMS_wcnf_clause_t* clause) {
  BouMS_uint_t count = 0;
  for (BouMS_uint_t l = 0; l < clause->numLiterals; ++l) {
    const BouMS_wcnf_literal_t* const literal = clause->literals + l;
    if (BouMS_wcnf_isLiteralSatisfied(formula, literal)) {
      count += 1;
    }
  }
  return count;
}

#ifdef __cplusplus
}
#endif

#endif
