/**
 * @file map.h
 * @author Ole Lübke (ole.luebke@tuhh.de)
 *
 * @copyright Copyright (c) 2025
 */

#ifndef BOUMS_MAP_H
#define BOUMS_MAP_H

#include "BouMS/common.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  BouMS_uint_t* ex2In;
  BouMS_uint_t* in2Ex;
} BouMS_clauseMap_t;

#ifdef __cplusplus
}
#endif

#endif
