/**
 * @file fixedprec.c
 * @author Ole Lübke (ole.luebke@tuhh.de)
 *
 * @copyright Copyright (c) 2024
 *
 */

#include "private/fixedprec.h"

#ifndef FIXEDPREC_FLOATING
#  include "BouMS/common.h"
#endif

fixedprec_shift_t fixedprec_determineMaxShift(fixedprec_orig_t maxUnshifted, fixedprec_shift_t safetyBits) {
#ifndef FIXEDPREC_FLOATING
  const unsigned int TOTAL_BITS = 8 * sizeof(BouMS_uint_t);

  fixedprec_shift_t msb = safetyBits;
  while (maxUnshifted >>= 1) {
    msb += 1;
  }

  int ret = ((int)TOTAL_BITS - 1) - (int)msb;
  return ret > 0 ? ret : 0;
#else
  ((void)maxUnshifted);
  ((void)safetyBits);
  return 0;
#endif
}
