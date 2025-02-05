/**
 * @file fixedprec.h
 * @author Ole Lübke (ole.luebke@tuhh.de)
 * @brief Fixed precision math functions
 *
 * @copyright Copyright (c) 2024
 *
 */

#ifndef FIXEDPREC_H
#define FIXEDPREC_H

// #define FIXEDPREC_FLOATING

#ifndef FIXEDPREC_FLOATING
#  include <assert.h>

#  include "BouMS/common.h"
#endif

#define UFIXEDPREC(x) ((unsigned_fixedprec_t){.value = (x)})
#define SFIXEDPREC(x) ((signed_fixedprec_t){.value = (x)})

typedef unsigned int fixedprec_shift_t;

#ifndef FIXEDPREC_FLOATING
typedef BouMS_uint_t fixedprec_orig_t;
#else
typedef double fixedprec_orig_t;
#endif

typedef struct {
  fixedprec_orig_t value;
} unsigned_fixedprec_t;

#ifndef FIXEDPREC_FLOATING
typedef long long fixedprec_sorig_t;
typedef struct {
  fixedprec_sorig_t value;
} signed_fixedprec_t;
#else
typedef double fixedprec_sorig_t;
typedef unsigned_fixedprec_t signed_fixedprec_t;
#endif

fixedprec_shift_t fixedprec_determineMaxShift(fixedprec_orig_t maxUnshifted, fixedprec_shift_t safetyBits);

static inline unsigned_fixedprec_t fixedprec_uto(fixedprec_orig_t val, fixedprec_shift_t shift) {
#ifndef FIXEDPREC_FLOATING
  return UFIXEDPREC(val << shift);
#else
  ((void)shift);
  return UFIXEDPREC(val);
#endif
}

static inline fixedprec_orig_t fixedprec_ufrom(unsigned_fixedprec_t val, fixedprec_shift_t shift) {
#ifndef FIXEDPREC_FLOATING
  return val.value >> shift;
#else
  ((void)shift);
  return val.value;
#endif
}

static inline signed_fixedprec_t fixedprec_utos(unsigned_fixedprec_t val) {
  return SFIXEDPREC(val.value);
}

static inline unsigned_fixedprec_t fixedprec_stou(signed_fixedprec_t val) {
  return UFIXEDPREC(val.value);
}

static inline unsigned_fixedprec_t fixedprec_uadd(unsigned_fixedprec_t summand1, unsigned_fixedprec_t summand2) {
  const fixedprec_orig_t sum = summand1.value + summand2.value;
#ifndef FIXEDPREC_FLOATING
  assert(sum >= summand1.value && sum >= summand2.value);
#endif
  return UFIXEDPREC(sum);
}

static inline unsigned_fixedprec_t fixedprec_usub(unsigned_fixedprec_t minuend, unsigned_fixedprec_t subtrahend) {
  const fixedprec_orig_t diff = minuend.value - subtrahend.value;
#ifndef FIXEDPREC_FLOATING
  assert(diff <= minuend.value);
#endif
  return UFIXEDPREC(diff);
}

static inline unsigned_fixedprec_t fixedprec_umul(unsigned_fixedprec_t factor1, unsigned_fixedprec_t factor2,
                                                  fixedprec_shift_t shift) {
#ifndef FIXEDPREC_FLOATING
  if (factor1.value >= factor2.value) {
    factor1.value >>= shift;
  } else {
    factor2.value >>= shift;
  }
#else
  ((void)shift);
#endif
  const fixedprec_orig_t prod = factor1.value * factor2.value;
#ifndef FIXEDPREC_FLOATING
  assert(((factor1.value == 0 || factor2.value == 0) && prod == 0) || (prod >= factor1.value && prod >= factor2.value));
#endif
  return UFIXEDPREC(prod);
}

static inline unsigned_fixedprec_t fixedprec_udiv(unsigned_fixedprec_t numerator, unsigned_fixedprec_t denominator,
                                                  fixedprec_shift_t shift) {
#ifndef FIXEDPREC_FLOATING
  const fixedprec_orig_t quot = numerator.value / (denominator.value >> shift);
  assert(quot <= numerator.value);
  return UFIXEDPREC(quot);
#else
  ((void)shift);
  return UFIXEDPREC(numerator.value / denominator.value);
#endif
}

static inline signed_fixedprec_t fixedprec_sadd(signed_fixedprec_t summand1, signed_fixedprec_t summand2) {
  const fixedprec_sorig_t sum = summand1.value + summand2.value;
#ifndef FIXEDPREC_FLOATING
  assert((summand1.value >= 0 && summand2.value >= 0 && sum >= summand1.value && sum >= summand2.value) ||
         (summand1.value >= 0 && summand2.value <= 0 && sum <= summand1.value && sum >= summand2.value) ||
         (summand1.value <= 0 && summand2.value >= 0 && sum >= summand1.value && sum <= summand2.value) ||
         (summand1.value <= 0 && summand2.value <= 0 && sum <= summand1.value && sum <= summand2.value));
#endif
  return SFIXEDPREC(sum);
}

static inline signed_fixedprec_t fixedprec_ssub(signed_fixedprec_t minuend, signed_fixedprec_t subtrahend) {
  const fixedprec_sorig_t diff = minuend.value - subtrahend.value;
#ifndef FIXEDPREC_FLOATING
  assert((minuend.value >= 0 && subtrahend.value >= 0 && diff <= minuend.value) ||
         (minuend.value >= 0 && subtrahend.value <= 0 && diff >= minuend.value) ||
         (minuend.value <= 0 && subtrahend.value >= 0 && diff <= minuend.value) ||
         (minuend.value <= 0 && subtrahend.value <= 0 && diff >= minuend.value));
#endif
  return SFIXEDPREC(diff);
}

#endif
