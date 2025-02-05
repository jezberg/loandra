/**
 * @file rng.h
 * @author Ole Lübke (ole.luebke@tuhh.de)
 * @brief Random number generation
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef RNG_H
#define RNG_H

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <time.h>

#include "BouMS/common.h"

#define RNG_MAXIMUM_PROBABILITY   100000000
#define RNG_FAIR_COIN_PROBABILITY (RNG_MAXIMUM_PROBABILITY / 2)

/**
 * @brief Initializes the RNG module
 */
static inline void rngInit(void) {
  // TODO don't use standard library for RNG?
  srand(time(NULL));
}

/**
 * @brief Generates a random index < max
 *
 * @param max
 * @return BouMS_uint_t
 */
static inline BouMS_uint_t randIdx(BouMS_uint_t max) {
  // TODO don't use standard library?
  // TODO truly uniform RNG for [0, max)?
  return rand() % max;
}

/**
 * @brief Returns true with (bias / RNG_MAXIMUM_PROBABILITY / 100)% probability, otherwise returns false
 *
 * @param bias [0, RNG_MAXIMUM_PROBABILITY]
 * @return true
 * @return false
 */
static inline bool flipBiasedCoin(BouMS_uint_t bias) {
  return randIdx(RNG_MAXIMUM_PROBABILITY) < bias;
}

/**
 * @brief Returns true with 50% probability, otherwise returns false
 *
 * @return true
 * @return false
 */
static inline bool flipCoin(void) {
  return flipBiasedCoin(RNG_FAIR_COIN_PROBABILITY);
}
#endif
