/**
 * @file dynmem.h
 * @author Ole Lübke (ole.luebke@tuhh.de)
 *
 * @copyright Copyright (c) 2025
 */

#ifndef BOUMS_DYNMEM_H
#define BOUMS_DYNMEM_H

#include <stdlib.h>

/**
 * @brief Type of realloc
 */
typedef void* (*BouMS_dynmem_realloc_t)(void*, size_t);

/**
 * @brief Type of free
 */
typedef void (*BouMS_dynmem_free_t)(void*);

#endif
