/**
 * @file main.c
 * @author Ole Lübke (ole.luebke@tuhh.de)
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <BouMS/BouMS.h>
#include <BouMS/common.h>
#include <BouMS/wcnf.h>
#include <BouMS/wcnf_util.h>
#include <argp.h>
#include <assert.h>
#include <errno.h>
#include <private/logging.h>
#include <signal.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "MaxSAT22.h"
#include "common.h"

#if BOUMS_USE_SATSOLVER
#  include "externalsat.h"
#endif

typedef struct {
  const char* file;
  BouMS_params_t params;
} arguments_t;

static void signalHandler(int signal);

static void registerSignalHandlers(void);

// NOLINTNEXTLINE(readability-non-const-parameter)
static error_t parseOpt(int key, char* arg, struct argp_state* state);

static void parseUIntOpt(BouMS_uint_t* val, const char* name, const char* arg, struct argp_state* state);

static void setParameters(arguments_t* args, const BouMS_wcnf_t* formula);

static void setParamIfNotUserProvided(BouMS_uint_t* param, BouMS_uint_t value);

static void printResult(const BouMS_result_t* result, BouMS_uint_t numVariables);

static int run(BouMS_wcnf_t* formula, const arguments_t* args);

const char* argp_program_version = "BouMS";
const char* argp_program_bug_address = "<ole.luebke@tuhh.de>";

static bool stop = false;

int main(int argc, char** argv) {
  INIT_DURATION_MEAS();

  char doc[] = "Try to solve MaxSAT problem read from FILE using BouMS";
  char argsDoc[] = "FILE";
  struct argp_option options[] = {
      {"maxTries", 't', "NUM", 0, "Maximum number of retries (outer loop).\nDefault: UINT64_MAX", 0},
      {"maxFlips", 'l', "NUM", 0,
       "Maximum number of flips (inner loop).\nIf not specified, maxFlipsFactor is used to determine this.", 0},
      {"maxRetriesWOImprovement", 'r', "NUM", 0,
       "Maximum number of retries without improvement, before restart with new random assignment is forced", 0},
      {"bmsSize", 'b', "NUM", 0, "Size for Best of Multiple Selection", 0},
      {"hWeightInit", 'i', "NUM", 0, "Hard clause initial weight, only for PMS, WPMS", 0},
      {"hWeightInc", 'h', "NUM", 0, "Hard clause weight increment, only for PMS, WPMS", 0},
      {"sWeightInc", 's', "NUM", 0, "Soft clause weight increment", 0},
      {"sSmoothProb", 'p', "NUM", 0,
       "Smoothing probability for soft clause weights, only for MS, WMS\n100000000 = 100%", 0},
      {"sWeightBoundCoef", 'c', "NUM", 0, "Coefficient for calculating soft clause weight bounds, only for MS, WMS", 0},
      {"sWeightBoundOffset", 'o', "NUM", 0, "Offset for calculating soft clause weight bounds, only for MS, WMS", 0},
      {0}};
  struct argp argp = {options, parseOpt, argsDoc, doc, NULL, NULL, NULL};
  arguments_t args = {
      .file = NULL,
      .params = {.maxTries = BOUMS_UINT_MAX,
                 .maxFlips = BOUMS_UINT_MAX,
                 .maxTriesWOImprovement = BOUMS_UINT_MAX,
                 .bmsSize = BOUMS_UINT_MAX,
                 .hWeightInit = BOUMS_UINT_MAX,
                 .hWeightInc = BOUMS_UINT_MAX,
                 .sWeightInc = BOUMS_UINT_MAX,
                 .sSmoothProb = BOUMS_UINT_MAX,
                 .sWeightBoundCoef = BOUMS_UINT_MAX,
                 .sWeightBoundOffset = BOUMS_UINT_MAX,
                 .fixedPrecSafetyBits = BOUMS_UINT_MAX},
  };

  registerSignalHandlers();

  argp_parse(&argp, argc, argv, 0, NULL, &args);

  int exitCode = 0;

  FILE* inputStream = fopen(args.file, "r");
  if (inputStream) {
    LOG_INFO("c starting to read problem\n");
    START_DURATION_MEAS();
    BouMS_uint_t errLine = 0;
    BouMS_wcnf_t* formula = wcnfFromMaxSAT22(inputStream, &errLine, &stop);
    LOG_INFO_WITH_DURATION("c read problem");
    // NOLINTNEXTLINE
    fclose(inputStream);

    if (formula) {
      setParameters(&args, formula);

      LOG_INFO("c Parameters:\n");
      LOG_INFO("c   maxTries:                " BOUMS_UINT_FORMAT "\n", args.params.maxTries);
      LOG_INFO("c   maxRetriesWOImprovement: " BOUMS_UINT_FORMAT "\n", args.params.maxTriesWOImprovement);
      LOG_INFO("c   maxFlips:                " BOUMS_UINT_FORMAT "\n", args.params.maxFlips);
      LOG_INFO("c   bmsSize:                 " BOUMS_UINT_FORMAT "\n", args.params.bmsSize);
      LOG_INFO("c   hWeightInit:             " BOUMS_UINT_FORMAT "\n", args.params.hWeightInit);
      LOG_INFO("c   hWeightInc:              " BOUMS_UINT_FORMAT "\n", args.params.hWeightInc);
      LOG_INFO("c   sWeightInc:              " BOUMS_UINT_FORMAT "\n", args.params.sWeightInc);
      LOG_INFO("c   sSmoothProb:             " BOUMS_UINT_FORMAT "\n", args.params.sSmoothProb);
      LOG_INFO("c   sWeightBoundCoef:        " BOUMS_UINT_FORMAT "\n", args.params.sWeightBoundCoef);
      LOG_INFO("c   sWeightBoundOffset:      " BOUMS_UINT_FORMAT "\n", args.params.sWeightBoundOffset);

      LOG_INFO("c running on a formula with " BOUMS_UINT_FORMAT " variables and " BOUMS_UINT_FORMAT
               " clauses (" BOUMS_UINT_FORMAT " hard + " BOUMS_UINT_FORMAT " soft)\n",
               formula->numVariables, formula->numClauses, formula->numHardClauses,
               formula->numClauses - formula->numHardClauses);

      exitCode = run(formula, &args);

      BouMS_wcnf_util_deleteFormula(formula, free, &stop);
      free(formula);
      formula = NULL;
    } else {
      if (errLine) {
        // NOLINTNEXTLINE
        fprintf(stderr, "c Syntax error on line " BOUMS_UINT_FORMAT ".\n", errLine);
        exitCode = EXIT_FAILURE;
      } else {
        // NOLINTNEXTLINE
        fprintf(stderr, "c Ran out of memory while parsing or stop flag was set.\n");
        printf("s UNKNOWN\n");
      }
    }
  } else {
    // NOLINTNEXTLINE
    fprintf(stderr, "c Could not open file \"%s\".\n", args.file);
    exitCode = EXIT_FAILURE;
  }

  return exitCode;
}

static void signalHandler(int signal) {
  switch (signal) {
    case SIGINT:  // fall-through
    case SIGTERM:
      /* printResult(&result, numVariables);
      exit(EXIT_SUCCESS); */
      stop = true;
      break;
    default:
      break;
  }
}

static void registerSignalHandlers(void) {
  // NOLINTNEXTLINE(performance-no-int-to-ptr)
  if (signal(SIGTERM, signalHandler) == SIG_ERR) {
    // NOLINTNEXTLINE
    fprintf(stderr, "c Failed to register SIGTERM handler\n");
    exit(EXIT_FAILURE);
  }
  // NOLINTNEXTLINE(performance-no-int-to-ptr)
  if (signal(SIGINT, signalHandler) == SIG_ERR) {
    // NOLINTNEXTLINE
    fprintf(stderr, "c Failed to register SIGINT handler\n");
    exit(EXIT_FAILURE);
  }
}

static error_t parseOpt(int key, char* arg, struct argp_state* state) {
  arguments_t* const arguments = state->input;

  switch (key) {
    case 't':
      parseUIntOpt(&arguments->params.maxTries, "maxTries", arg, state);
      break;
    case 'l':
      parseUIntOpt(&arguments->params.maxFlips, "maxFlips", arg, state);
      break;
    case 'r':
      parseUIntOpt(&arguments->params.maxTriesWOImprovement, "maxRetriesWOImprovement", arg, state);
      break;
    case 'b':
      parseUIntOpt(&arguments->params.bmsSize, "bmsSize", arg, state);
      if (arguments->params.bmsSize < 1) {
        argp_error(state, "bmsSize should be >= 1");
      }
      break;
    case 'i':
      parseUIntOpt(&arguments->params.hWeightInit, "hWeightInit", arg, state);
      break;
    case 'h':
      parseUIntOpt(&arguments->params.hWeightInc, "hWeightInc", arg, state);
      break;
    case 's':
      parseUIntOpt(&arguments->params.sWeightInc, "sWeightInc", arg, state);
      break;
    case 'p':
      parseUIntOpt(&arguments->params.sSmoothProb, "sSmoothProb", arg, state);
      if (arguments->params.sSmoothProb > 100000000) {  // NOLINT(readability-magic-numbers)
        argp_error(state, "0 <= sSmoothProb <= 100000000");
      }
      break;
    case 'c':
      parseUIntOpt(&arguments->params.sWeightBoundCoef, "sWeightBoundCoef", arg, state);
      break;
    case 'o':
      parseUIntOpt(&arguments->params.sWeightBoundOffset, "sWeightBoundOffset", arg, state);
      break;
    case ARGP_KEY_ARG:
      if (state->arg_num >= 1) {
        argp_usage(state);
      }
      arguments->file = arg;
      break;
    case ARGP_KEY_END:
      if (!arguments->file) {
        argp_usage(state);
      }
      break;
    default:
      return ARGP_ERR_UNKNOWN;
  }

  return EXIT_SUCCESS;
}

static void parseUIntOpt(BouMS_uint_t* val, const char* name, const char* arg, struct argp_state* state) {
  if (arg == NULL) {
    argp_error(state, "Missing value for %s", name);
    return;
  }
  *val = strtoull(arg, NULL, BOUMS_DECIMAL_SYSTEM_BASE);
  if (errno == ERANGE) {
    argp_error(state, "Could not parse %s: %s", name, arg);
  }
}

static void setParamIfNotUserProvided(BouMS_uint_t* param, BouMS_uint_t value) {
  if (*param == BOUMS_UINT_MAX) {
    *param = value;
  }
}

static void setParameters(arguments_t* args, const BouMS_wcnf_t* formula) {
  BouMS_params_t defaultParams;
  BouMS_params(formula, &defaultParams);

  args->params.isPartial = defaultParams.isPartial;
  args->params.isWeighted = defaultParams.isWeighted;

  setParamIfNotUserProvided(&args->params.maxTries, defaultParams.maxTries);
  setParamIfNotUserProvided(&args->params.maxTriesWOImprovement, defaultParams.maxTriesWOImprovement);
  setParamIfNotUserProvided(&args->params.maxFlips, defaultParams.maxFlips);
  setParamIfNotUserProvided(&args->params.bmsSize, defaultParams.bmsSize);
  setParamIfNotUserProvided(&args->params.hWeightInit, defaultParams.hWeightInit);
  setParamIfNotUserProvided(&args->params.hWeightInc, defaultParams.hWeightInc);
  setParamIfNotUserProvided(&args->params.sWeightInc, defaultParams.sWeightInc);
  setParamIfNotUserProvided(&args->params.sSmoothProb, defaultParams.sSmoothProb);
  setParamIfNotUserProvided(&args->params.sWeightBoundCoef, defaultParams.sWeightBoundCoef);
  setParamIfNotUserProvided(&args->params.sWeightBoundOffset, defaultParams.sWeightBoundOffset);
  setParamIfNotUserProvided(&args->params.fixedPrecSafetyBits, defaultParams.fixedPrecSafetyBits);
}

static void printResult(const BouMS_result_t* result, BouMS_uint_t numVariables) {
  static const char* const msgs[BOUMS_STATUS_COUNT] = {"OPTIMUM FOUND", "SATISFIABLE", "UNSATISFIABLE", "UNKNOWN"};

  assert(result->status >= 0 && result->status < BOUMS_STATUS_COUNT);
  const char* const msg = msgs[result->status];

  printf("s %s\n", msg);
  if (result->status != BOUMS_INFEASIBLE && result->status != BOUMS_UNSAT) {
    printf("o " BOUMS_UINT_FORMAT "\n", result->cost);

    if (wcnfVariablesToMaxSAT22(stdout, result->assignment, numVariables)) {
      // NOLINTNEXTLINE
      fputs("c Failed to write result.\n", stderr);
    }
  }
}

static int run(BouMS_wcnf_t* formula, const arguments_t* args) {
  static const int ret[BOUMS_STATUS_COUNT] = {30, 10, 20, 0};

  BouMS_result_t result = {.status = BOUMS_INFEASIBLE, .cost = BOUMS_UINT_MAX, .assignment = NULL};

  bool* initModel = NULL;

#ifdef BOUMS_USE_SATSOLVER
  if (formula->numHardClauses > 0) {
    initModel = malloc(formula->numVariables * sizeof(bool));
    if (initModel) {
      LOG_INFO("c pre-solving hard clauses\n");
      INIT_DURATION_MEAS();
      START_DURATION_MEAS();
      const sat_result_t satResult = preSolveHard(formula, initModel, &stop);
      LOG_INFO_WITH_DURATION("c pre-solving finished");

      if (satResult == UNKNOWN || satResult == UNSAT) {
        free(initModel);
        initModel = NULL;

        if (satResult == UNKNOWN) {
          result.status = BOUMS_INFEASIBLE;
          // NOLINTNEXTLINE
          fputs("c failed to pre-solve hard clauses\n", stderr);
        } else {
          result.status = BOUMS_UNSAT;
        }
        printResult(&result, 0);

        return ret[result.status];
      }
    } else {
      // NOLINTNEXTLINE
      fputs("c Out of memory\n", stderr);
      return ret[result.status];
    }

    result.status = BOUMS_UNKNOWN;
    LOG_INFO("c hard clauses are SAT\n");
  }
#endif

  BouMS_init();

  BouMS_memoryReq_t detailMemReq;
  const BouMS_uint_t memReq = BouMS_calcMemoryRequirements(formula, &detailMemReq);
  void* memory = malloc(memReq);

  result.assignment = malloc(formula->numVariables * sizeof(bool));

  if (memory && result.assignment) {
    BouMS_solve(formula, &args->params, memory, &detailMemReq, &result, initModel, 0, &stop);
    printResult(&result, formula->numVariables);
  } else {
    // NOLINTNEXTLINE
    fputs("c Out of memory\n", stderr);
  }

  if (initModel) {
    free(initModel);
    initModel = NULL;
  }

  if (memory) {
    free(memory);
    memory = NULL;
  }
  if (result.assignment) {
    free(result.assignment);
    result.assignment = NULL;
  }

  assert(result.status >= 0 && result.status < BOUMS_STATUS_COUNT);
  return ret[result.status];
}
