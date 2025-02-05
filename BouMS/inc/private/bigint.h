/**
 * @file bigint.h
 * @author Ole Lübke (ole.luebke@tuhh.de)
 * @brief Big integer type
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef BOUMS_PRIVATE_BIGINT_H
#define BOUMS_PRIVATE_BIGINT_H

#include <stdbool.h>

#include "BouMS/common.h"

typedef struct {
  bool negative;
  BouMS_uint_t value;
} bigint_t;

static inline void bigint_add(bigint_t* sum, BouMS_uint_t summand) {
  if (sum->negative) {
    if (sum->value > summand) {
      sum->value -= summand;
    } else {
      sum->negative = false;
      sum->value = (summand - sum->value);
    }
  } else {
    sum->value += summand;
  }
}

static inline void bigint_sub(bigint_t* diff, BouMS_uint_t subtrahend) {
  if (!diff->negative) {
    if (diff->value >= subtrahend) {
      diff->value -= subtrahend;
    } else {
      diff->negative = true;
      diff->value = (subtrahend - diff->value);
    }
  } else {
    diff->value += subtrahend;
  }
}

static inline bool bigint_gt(const bigint_t* a, const bigint_t* b) {
  if (a->negative && b->negative) {
    return a->value < b->value;
  } else if (!a->negative && !b->negative) {
    return a->value > b->value;
  } else {
    return b->negative;
  }
}

#endif
