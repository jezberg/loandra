/**
 * @file cores.c
 * @author Ole Lübke (ole.luebke@tuhh.de)
 *
 * @copyright Copyright (c) 2025
 */

#include "BouMS/cores.h"

#include <stdlib.h>

#include "BouMS/common.h"
#include "BouMS/dynmem.h"
#include "BouMS/wcnf.h"

bool BouMS_cores_init(const BouMS_wcnf_t* formula, const BouMS_cores_core_t* cores, BouMS_uint_t numCores,
                      const BouMS_uint_t* varToClause, BouMS_cores_mem_t* out, BouMS_dynmem_realloc_t realloc,
                      BouMS_dynmem_free_t free) {
  BouMS_cores_free(formula, out, free);

  out->numCores = numCores;
  out->cores = cores;
  out->varToClause = varToClause;

  out->coreToSatLits = (BouMS_uint_t*)realloc(NULL, numCores * sizeof(BouMS_uint_t));
  if (!out->coreToSatLits) {
    return true;
  }

  out->varToCores = (BouMS_cores_list_t**)realloc(NULL, formula->numVariables * sizeof(BouMS_cores_list_t*));
  if (!out->varToCores) {
    free(out->coreToSatLits);
    out->coreToSatLits = NULL;

    return true;
  }

  for (BouMS_uint_t vIdx = 0; vIdx < formula->numVariables; ++vIdx) {
    out->varToCores[vIdx] = NULL;
  }

  bool oom = false;
  for (BouMS_uint_t cIdx = 0; cIdx < numCores && !oom; ++cIdx) {
    const BouMS_cores_core_t* const core = cores + cIdx;

    for (BouMS_uint_t lIdx = 0; lIdx < core->numLiterals && !oom; ++lIdx) {
      const BouMS_literal_t lit = core->literals[lIdx];
      const BouMS_uint_t var = BouMS_var(lit);
      BouMS_cores_list_t** const coreListPtr = out->varToCores + var;

      if (!*coreListPtr) {
        *coreListPtr = (BouMS_cores_list_t*)realloc(NULL, sizeof(BouMS_cores_list_t));
        if (!*coreListPtr) {
          oom = true;
          break;
        }

        BouMS_cores_list_t* const coreList = *coreListPtr;
        coreList->numCores = 1;
        coreList->cores = (const BouMS_cores_core_t**)realloc(NULL, sizeof(BouMS_cores_core_t*));
        if (!coreList->cores) {
          oom = true;
          break;
        }

        coreList->cores[0] = core;
      } else {
        BouMS_cores_list_t* const coreList = *coreListPtr;
        const BouMS_uint_t newNumCores = coreList->numCores + 1;
        const BouMS_cores_core_t** const newCoreList =
            (const BouMS_cores_core_t**)realloc((void*)coreList->cores, newNumCores * sizeof(BouMS_cores_core_t*));
        if (!newCoreList) {
          oom = true;
          break;
        }

        coreList->cores = newCoreList;
        coreList->cores[coreList->numCores] = core;
        coreList->numCores = newNumCores;
      }
    }
  }

  if (oom) {
    BouMS_cores_free(formula, out, free);
    return true;
  }

  return false;
}

void BouMS_cores_free(const BouMS_wcnf_t* formula, BouMS_cores_mem_t* cores, BouMS_dynmem_free_t free) {
  cores->numCores = 0;
  cores->cores = NULL;

  if (cores->coreToSatLits) {
    free(cores->coreToSatLits);
    cores->coreToSatLits = NULL;
  }

  if (cores->varToCores) {
    for (BouMS_uint_t vIdx = 0; vIdx < formula->numVariables; ++vIdx) {
      BouMS_cores_list_t** const coreListPtr = cores->varToCores + vIdx;
      if (*coreListPtr) {
        BouMS_cores_list_t* coreList = *coreListPtr;

        if (coreList->cores) {
          free((void*)coreList->cores);
          coreList->cores = NULL;
          coreList->numCores = 0;
        }

        free(coreList);
        *coreListPtr = NULL;
      }
    }
  }
}
