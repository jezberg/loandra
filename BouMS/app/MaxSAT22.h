/**
 * @file MaxSAT22.h
 * @author Ole Lübke <ole.luebke@tuhh.de)
 * @brief MaxSAT Evaluation 2022 WCNF IO functions
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef MAXSAT22_H
#define MAXSAT22_H

#include <BouMS/common.h>
#include <BouMS/wcnf.h>
#include <stdbool.h>
#include <stdio.h>

/**
 * @brief Reads a MaxSAT problem in
 * [MaxSAT Evaluation 2022 format](https://maxsat-evaluations.github.io/2022/rules.html#input).
 *
 * @param stream The input stream to read from, e.g., a file or STDIN.
 * @param errLine In case of a parse error, `errLine` is set to the line number
 * on which the error occurred.
 * When an error occurred which was not a parse error, `errLine` is set to 0.
 * `errLine` is an optional argument, i.e., it may be NULL.
 * @param stop Flag to signal the function to exit early (may cause memory leak!)
 * @return BouMS_wcnf_t* Pointer to the resulting formula, should be freed
 * with `deleteWCNF()`.
 */
BouMS_wcnf_t* wcnfFromMaxSAT22(FILE* stream, BouMS_uint_t* errLine, const bool* stop);

/**
 * @brief Writes a given formula to a given stream in MaxSAT Evaluation 2022 format.
 *
 * @param stream The stream to write to.
 * @param formula The formula to write.
 * @return true In case of error.
 * @return false When no error occurred.
 */
bool wcnfToMaxSAT22(FILE* stream, const BouMS_wcnf_t* formula);

/**
 * @brief Writes the variables of a formula to the given stream in MaxSAT Evaluation 2022 format.
 *
 * @param stream The stream to write to.
 * @param variables
 * @param numVariables
 * @return true In case of error.
 * @return false When no error occured.
 */
bool wcnfVariablesToMaxSAT22(FILE* stream, const bool* variables, BouMS_uint_t numVariables);

#endif
