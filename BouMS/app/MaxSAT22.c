/**
 * @file MaxSAT22.c
 * @author Ole Lübke (ole.luebke@tuhh.de)
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "MaxSAT22.h"

#include <BouMS/common.h>
#include <BouMS/wcnf.h>
#include <BouMS/wcnf_util.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"

#define INITIAL_LIT_BUF_SZ 256

// only for testing!
// #define NO_BATCH_CLAUSE_ADDING

typedef enum { PARSE_RESULT_OK = 0, PARSE_RESULT_ERROR, PARSE_RESULT_OOM } parseResult_t;

/**
 * @brief Advances pointer until next non-whitespace character
 *
 * @param str
 */
void skipWhitespace(char** str);

/**
 * @brief Counts the number of clauses in the stream
 *
 * @param stream
 * @return BouMS_uint_t
 */
BouMS_uint_t findNumClauses(FILE* stream);

/**
 * @brief Reads the next unsigned integer from string and advances str to behind the number
 *
 * @param str
 * @param num
 * @return true In case of error
 * @return false Otherwise
 */
bool readWeight(char** str, BouMS_uint_t* num);

/**
 * @brief Count the number of literals in a clause line
 *
 * @param str Needs to point AFTER the weight
 * @return BouMS_uint_t
 */
BouMS_uint_t findNumLiterals(char* str);

/**
 * @brief Parses a clause line
 *
 * @param line
 * @param clauseAddingState
 * @param litBufSz
 * @param litBufP
 * @return true
 * @return false
 */
parseResult_t parseClause(char* line,
#ifndef NO_BATCH_CLAUSE_ADDING
                          BouMS_wcnf_util_batchClauseAddingState_t* clauseAddingState,
#else
                          BouMS_wcnf_t* formula,
#endif
                          size_t* litBufSz, int** litBufP);

BouMS_wcnf_t* wcnfFromMaxSAT22(FILE* stream, BouMS_uint_t* errLine, const bool* stop) {
  BouMS_wcnf_t* formula = calloc(1, sizeof(BouMS_wcnf_t));
  if (!formula) {
    return NULL;
  }

  BouMS_uint_t lineNo = 0;

#ifndef NO_BATCH_CLAUSE_ADDING
  BouMS_wcnf_util_batchClauseAddingState_t clauseAddingState;
  if (!BouMS_wcnf_util_startBatchClauseAdding(findNumClauses(stream), realloc, free, &clauseAddingState)) {
#endif
    size_t litBufSz = INITIAL_LIT_BUF_SZ;
    int* litBuf = malloc(litBufSz * sizeof(int));
    if (litBuf) {
      bool error = false;
      char* line = NULL;
      size_t numChars;
      while (getline(&line, &numChars, stream) != -1 && !error && !*stop) {
        lineNo += 1;

#ifndef NO_BATCH_CLAUSE_ADDING
        const parseResult_t parseResult = parseClause(line, &clauseAddingState, &litBufSz, &litBuf);
#else
      const parseResult_t parseResult = parseClause(line, formula, &litBufSz, &litBuf);
#endif
        if (parseResult != PARSE_RESULT_OK) {
          error = true;
          if (parseResult == PARSE_RESULT_OOM) {
            lineNo = 0;
          }
          break;
        }
      }

      free(line);
      line = NULL;
      free(litBuf);
      litBuf = NULL;

      if (*stop) {
        lineNo = 0;
      } else if (!error) {
#ifndef NO_BATCH_CLAUSE_ADDING
        if (!BouMS_wcnf_util_finishBatchClauseAdding(&clauseAddingState, formula, stop)) {
#endif
          return formula;
#ifndef NO_BATCH_CLAUSE_ADDING
        }
        lineNo = 0;
#endif
      }
    }
#ifndef NO_BATCH_CLAUSE_ADDING
  }
  BouMS_wcnf_util_cleanUpBatchClauseAddingAfterError(&clauseAddingState);
#endif

  *errLine = lineNo;
  BouMS_wcnf_util_deleteFormula(formula, free, stop);
  formula = NULL;

  return NULL;
}

bool wcnfToMaxSAT22(FILE* stream, const BouMS_wcnf_t* formula) {
  for (BouMS_uint_t clauseIdx = 0; clauseIdx < formula->numClauses; ++clauseIdx) {
    const BouMS_wcnf_clause_t* const clause = (formula->clauses + clauseIdx);
    if (BouMS_wcnf_isClauseHard(clause)) {
      if (fputs("h ", stream) == EOF) {
        return true;
      }
    } else {
      // NOLINTNEXTLINE(clang-analyzer-security.insecureAPI.DeprecatedOrUnsafeBufferHandling)
      if (fprintf(stream, BOUMS_UINT_FORMAT " ", clause->weight) < 0) {
        return true;
      }
    }

    for (BouMS_uint_t litIdx = 0; litIdx < clause->numLiterals; ++litIdx) {
      const BouMS_wcnf_literal_t* const literal = clause->literals + litIdx;
      if (BouMS_wcnf_sign(literal)) {
        if (fputc('-', stream) == EOF) {
          return true;
        }
      }
      // NOLINTNEXTLINE(clang-analyzer-security.insecureAPI.DeprecatedOrUnsafeBufferHandling)
      if (fprintf(stream, BOUMS_UINT_FORMAT " ", BouMS_wcnf_var(literal) + 1) < 0) {
        return true;
      }
    }

    if (fputs("0\n", stream) == EOF) {
      return true;
    }
  }

  return false;
}

bool wcnfVariablesToMaxSAT22(FILE* stream, const bool* variables, BouMS_uint_t numVariables) {
  if (fputs("v ", stream) == EOF) {
    return true;
  }

  for (BouMS_uint_t varIdx = 0; varIdx < numVariables; ++varIdx) {
    const char value = variables[varIdx] ? '1' : '0';
    if (fputc(value, stream) == EOF) {
      return true;
    }
  }

  if (fputc('\n', stream) == EOF) {
    return true;
  }

  return false;
}

void skipWhitespace(char** str) {
  while (isspace(**str) && **str != 0) {
    ++*str;
  }
}

BouMS_uint_t findNumClauses(FILE* stream) {
  char* line = NULL;
  size_t numChars = 0;
  BouMS_uint_t numClauses = 0;

  while (getline(&line, &numChars, stream) != -1) {
    if (line[0] != 'c' && strlen(line) > 0) {
      numClauses += 1;
    }
  }

  free(line);
  rewind(stream);

  return numClauses;
}

bool readWeight(char** str, BouMS_uint_t* num) {
  if (*str[0] == 'h') {
    *num = BOUMS_HARD_CLAUSE_WEIGHT;
    *str += 1;
  } else {
    *num = strtoull(*str, str, BOUMS_DECIMAL_SYSTEM_BASE);
    if (*num == ULLONG_MAX) {
      return true;
    }
  }
  return false;
}

BouMS_uint_t findNumLiterals(char* str) {
  // count literals by counting spaces in line
  BouMS_uint_t numLits = 0;
  for (; *str != '\n'; ++str) {
    if (isspace(*str)) {
      numLits += 1;
      skipWhitespace(&str);
    }
  }
  return numLits;
}

parseResult_t parseClause(char* line,
#ifndef NO_BATCH_CLAUSE_ADDING
                          BouMS_wcnf_util_batchClauseAddingState_t* clauseAddingState,
#else
                          BouMS_wcnf_t* formula,
#endif
                          size_t* litBufSz, int** litBufP) {
  if (line[0] == 'c' || strlen(line) == 0) {  // comment line
    return PARSE_RESULT_OK;
  }

  BouMS_uint_t weight = 0;
  if (readWeight(&line, &weight)) {
    return PARSE_RESULT_ERROR;
  }

  // ignore 0-weight soft clauses
  if (weight == 0) {
    return PARSE_RESULT_OK;
  }

  skipWhitespace(&line);

  const BouMS_uint_t numLiterals = findNumLiterals(line);
  if (numLiterals == 0) {
    return PARSE_RESULT_OK;
  }

  if (numLiterals > *litBufSz) {
    free(*litBufP);
    *litBufSz = 0;
    const size_t newLitBufSz = 2 * numLiterals;
    *litBufP = malloc(newLitBufSz * sizeof(int));
    if (!litBufP) {
      return PARSE_RESULT_OOM;
    }
    *litBufSz = newLitBufSz;
  }

  // read variables
  BouMS_uint_t litIdx = 0;
  while (*line != '0') {
    const int lit = (int)strtol(line, &line, BOUMS_DECIMAL_SYSTEM_BASE);
    if (errno == ERANGE || lit == 0) {
      return PARSE_RESULT_ERROR;
    }

    (*litBufP)[litIdx++] = lit;

    skipWhitespace(&line);
  }
  if (litIdx != numLiterals) {
    return PARSE_RESULT_ERROR;
  }

  if (
#ifndef NO_BATCH_CLAUSE_ADDING
      BouMS_wcnf_util_batchAddClause(clauseAddingState, weight, *litBufP, numLiterals)
#else
      BouMS_wcnf_util_addClause(formula, weight, *litBufP, numLiterals, realloc, free)
#endif
  ) {
    return PARSE_RESULT_OOM;
  }

  return PARSE_RESULT_OK;
}
