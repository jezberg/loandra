/**
 * @file partial.h
 * @author Ole Lübke (ole.luebke@tuhh.de)
 * @brief Functions for solving partial instances
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef PARTIAL_H
#define PARTIAL_H

#include <stdbool.h>

#include "BouMS/BouMS.h"
#include "BouMS/wcnf.h"
#include "private/fixedprec.h"
#include "private/types.h"

void p_initVars(BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem, const BouMS_memoryReq_t* memReq,
                BouMS_result_t* result, const bool* initModel, bool* reuseDecimation, bool forceRandomAssignment,
                const bool* stop);

void pw_initWeights(const BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem,
                    unsigned_fixedprec_t avgSoftClauseWeight);
void puw_initWeights(const BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem);

void p_updateWeights(const BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem);

#endif
