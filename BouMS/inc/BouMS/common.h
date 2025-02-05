/**
 * @file common.h
 * @author Ole Lübke (ole.luebke@tuhh.de)
 * @brief Common definitions used across other files of the project.
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef BOUMS_COMMON_H
#define BOUMS_COMMON_H

#include <limits.h>

#ifdef __cplusplus
extern "C" {
#endif

#define BOUMS_UINT_FORMAT "%llu"
#define BOUMS_UINT_MAX    ULLONG_MAX
typedef unsigned long long BouMS_uint_t;

#ifdef __cplusplus
}
#endif

#endif
