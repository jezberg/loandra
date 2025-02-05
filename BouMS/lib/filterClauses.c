/**
 * @file filterClauses.c
 * @author Ole Lübke (ole.luebke@tuhh.de)
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "BouMS/common.h"
#include "BouMS/preprocessing.h"
#include "BouMS/wcnf.h"
#include "private/common.h"

BouMS_uint_t BouMS_filterClauses(BouMS_wcnf_clause_t* clauses, BouMS_uint_t numClauses,
                                 const BouMS_clause_predicate_t* predicates, BouMS_uint_t numPredicates) {
  BouMS_uint_t numRemovedClauses = 0;

  for (BouMS_uint_t clauseIdx = 0; clauseIdx < numClauses; ++clauseIdx) {
    for (BouMS_uint_t predIdx = 0; predIdx < numPredicates; ++predIdx) {
      if (predicates[predIdx](clauses + clauseIdx)) {
        removeClause(clauseIdx, clauses, numClauses);
        numClauses -= 1;
        numRemovedClauses += 1;
        break;
      }
    }
  }

  return numRemovedClauses;
}
