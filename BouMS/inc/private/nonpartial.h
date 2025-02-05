/**
 * @file nonpartial.h
 * @author Ole Lübke (ole.luebke@tuhh.de)
 * @brief Functions for solving non-partial instances
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef NONPARTIAL_H
#define NONPARTIAL_H

#include <stdbool.h>

#include "BouMS/BouMS.h"
#include "BouMS/wcnf.h"
#include "private/fixedprec.h"
#include "private/types.h"

void np_initVars(BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem, const BouMS_memoryReq_t* memReq,
                 BouMS_result_t* result, const bool* initModel, bool* reuseDecimation, bool forceRandomAssignment,
                 const bool* stop);

void npw_initWeights(const BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem,
                     unsigned_fixedprec_t avgSoftClauseWeight);
void npuw_initWeights(const BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem);

void np_updateWeights(const BouMS_wcnf_t* formula, const BouMS_params_t* cfg, BouMS_memory_t* mem);

#endif
