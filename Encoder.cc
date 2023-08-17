/*!
 * \author Ruben Martins - ruben@sat.inesc-id.pt
 *
 * @section LICENSE
 *
 * Open-WBO, Copyright (c) 2013-2018, Ruben Martins, Vasco Manquinho, Ines Lynce
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */

#include "Encoder.h"

using namespace openwbo;





/************************************************************************************************
 //
 // Encoding of cardinality constraints
 //
 ************************************************************************************************/
//
// Manages the encoding of cardinality encodings.
void Encoder::encodeCardinality(CaDiCaL::Solver * SC, vec<Lit> &lits, int64_t rhs) {

  vec<Lit> lits_copy;
  lits.copyTo(lits_copy);

  switch (cardinality_encoding) {
  case _CARD_TOTALIZER_:
    totalizer.build(SC, lits_copy, rhs);
    if (totalizer.hasCreatedEncoding())
      totalizer.update(SC, rhs);
    break;

  default:
    printf("c Error: Invalid cardinality encoding.\n");
    printf("s UNKNOWN\n");
    exit(_ERROR_);
  }
}


void Encoder::addCardinality(CaDiCaL::Solver * SC, Encoder &enc, int64_t rhs) {
  if (cardinality_encoding == _CARD_TOTALIZER_ &&
      enc.cardinality_encoding == _CARD_TOTALIZER_) {
    totalizer.add(SC, enc.totalizer, rhs);
  } else {
    printf("c Error: Cardinality encoding does not support incrementality.\n");
    printf("s UNKNOWN\n");
    exit(_ERROR_);
  }
}

// Manages the update of cardinality constraints.
void Encoder::updateCardinality(CaDiCaL::Solver * SC, int64_t rhs) {

  switch (cardinality_encoding) {
  case _CARD_TOTALIZER_:
    totalizer.update(SC, rhs);
    break;


  default:
    printf("c Error: Invalid cardinality encoding.\n");
    printf("s UNKNOWN\n");
    exit(_ERROR_);
  }
}


// Incremental methods for cardinality encodings:
//
// Manages the building of cardinality encodings.
// Currently is only used for incremental solving.
void Encoder::buildCardinality(CaDiCaL::Solver * SC, vec<Lit> &lits, int64_t rhs) {
  assert(incremental_strategy != _INCREMENTAL_NONE_);

  vec<Lit> lits_copy;
  lits.copyTo(lits_copy);

  switch (cardinality_encoding) {
  case _CARD_TOTALIZER_:
    totalizer.build(SC, lits_copy, rhs);
    break;

  default:
    printf("c Error: Cardinality encoding does not support incrementality.\n");
    printf("s UNKNOWN\n");
    exit(_ERROR_);
  }
}

// Manages the incremental update of cardinality constraints.
void Encoder::incUpdateCardinality(CaDiCaL::Solver * SC, vec<Lit> &join, vec<Lit> &lits,
                                   int64_t rhs, vec<Lit> &assumptions) {
  assert(incremental_strategy == _INCREMENTAL_ITERATIVE_ ||
         incremental_strategy == _INCREMENTAL_WEAKENING_);

  vec<Lit> join_copy;
  join.copyTo(join_copy);
  vec<Lit> lits_copy;
  lits.copyTo(lits_copy);
  // Note: the assumption vector will be updated in this procedure

  switch (cardinality_encoding) {
  case _CARD_TOTALIZER_:
    if (join.size() > 0)
      totalizer.join(SC, join_copy, rhs);

    assert(lits.size() > 0);
    totalizer.update(SC, rhs, lits_copy, assumptions);
    break;

  default:
    printf("c Error: Cardinality encoding does not support incrementality.\n");
    printf("s UNKNOWN\n");
    exit(_ERROR_);
  }
}

void Encoder::joinEncoding(CaDiCaL::Solver * SC, vec<Lit> &lits, int64_t rhs) {

  switch (cardinality_encoding) {
  case _CARD_TOTALIZER_:
    totalizer.join(SC, lits, rhs);
    break;

  default:
    printf("c Error: Cardinality encoding does not support incrementality.\n");
    printf("s UNKNOWN\n");
    exit(_ERROR_);
  }
}

/************************************************************************************************
 //
 // Encoding of pseudo-Boolean constraints
 //
 ************************************************************************************************/
//
// Manages the encoding of PB encodings.
void Encoder::encodePB(CaDiCaL::Solver * SC, vec<Lit> &lits, vec<uint64_t> &coeffs,
                       uint64_t rhs) {

  vec<Lit> lits_copy;
  lits.copyTo(lits_copy);
  vec<uint64_t> coeffs_copy;
  coeffs.copyTo(coeffs_copy);

  switch (pb_encoding) {

  case _PB_GTE_:
    gte.encode(SC, lits_copy, coeffs_copy, rhs);
    break;

  case _PB_ADDER_:
    adder.encode(SC, lits_copy, coeffs_copy, rhs);
    break;

  default:
    printf("c Error: Invalid PB encoding.\n");
    printf("s UNKNOWN\n");
    exit(_ERROR_);
  }
}

int Encoder::predictPB(Solver *S, vec<Lit> &lits, vec<uint64_t> &coeffs,
                       uint64_t rhs) {

  vec<Lit> lits_copy;
  lits.copyTo(lits_copy);
  vec<uint64_t> coeffs_copy;
  coeffs.copyTo(coeffs_copy);

  switch (pb_encoding) {
  case _PB_SWC_:
    return -1;
    break;

  case _PB_GTE_:
    return gte.predict(S, lits_copy, coeffs_copy, rhs);
    break;

  case _PB_ADDER_:
    return -1;
    break;

  default:
    printf("Error: Invalid PB encoding.\n");
    printf("s UNKNOWN\n");
    exit(_ERROR_);
  }
}


// Manages the update of PB encodings.
void Encoder::updatePB(CaDiCaL::Solver * SC, uint64_t rhs) {

  switch (pb_encoding) {

  case _PB_GTE_:
    gte.update(SC, rhs);
    break;

  case _PB_ADDER_:
    adder.update(SC, rhs);
    break;

  default:
    printf("Error: Invalid PB encoding.\n");
    printf("s UNKNOWN\n");
    exit(_ERROR_);
  }
}

vec<Lit> &Encoder::lits() {
  assert(cardinality_encoding == _CARD_TOTALIZER_ &&
         incremental_strategy == _INCREMENTAL_ITERATIVE_);

  return totalizer.lits();
}

vec<Lit> &Encoder::outputs() {
  assert(cardinality_encoding == _CARD_TOTALIZER_ &&
         incremental_strategy == _INCREMENTAL_ITERATIVE_);

  return totalizer.outputs();
}

/************************************************************************************************
 //
 // Other
 //
 ************************************************************************************************/
// Returns true if the cardinality encoding was built, false otherwise.
bool Encoder::hasCardEncoding() {

  if (cardinality_encoding == _CARD_TOTALIZER_)
    return totalizer.hasCreatedEncoding();

  return false;
}

// Returns true if the PB encoding was built, false otherwise.
bool Encoder::hasPBEncoding() {
  if (pb_encoding == _PB_GTE_)
    return gte.hasCreatedEncoding();
  else if (pb_encoding == _PB_ADDER_)
    return adder.hasCreatedEncoding();
  
   printf("Error: should not be here.\n");
  return false;
}
