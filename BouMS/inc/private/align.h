/**
 * @file align.h
 * @author Ole Lübke (ole.luebke@tuhh.de)
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef BOUMS_PRIVATE_ALIGN_H
#define BOUMS_PRIVATE_ALIGN_H

#include "BouMS/common.h"

/**
 * @brief Make sure x is divisible by the pointer size of the system
 *
 * @param x
 * @return BouMS_uint_t x or the next number after x that is divisible by the pointer size of the system
 */
static inline BouMS_uint_t align(BouMS_uint_t x) {
  static const BouMS_uint_t ptrSize = sizeof(void*);
  const BouMS_uint_t mod = x % ptrSize;
  return mod != 0 ? x + ptrSize - mod : x;
}

#endif
