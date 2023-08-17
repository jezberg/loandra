/*!
 * \author Ruben Martins - ruben@sat.inesc-id.pt
 *
 * @section LICENSE
 *
 * Open-WBO, Copyright (c) 2013-2018, Ruben Martins, Vasco Manquinho, Ines Lynce
 * PBLib,    Copyright (c) 2012-2013  Peter Steinke
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

#include "Enc_Adder.h"
#include <algorithm>
#include <numeric>

using namespace openwbo;

void Adder::FA_extra ( Solver *S, CaDiCaL::Solver * SC, Lit xc, Lit xs, Lit a, Lit b, Lit c )
{
  
  clause.clear();
  addTernaryClause(S, SC, ~xc, ~xs, a);
  addTernaryClause(S, SC, ~xc, ~xs, b);
  addTernaryClause(S, SC, ~xc, ~xs, c);

  addTernaryClause(S, SC, xc, xs, ~a);
  addTernaryClause(S, SC, xc, xs, ~b);
  addTernaryClause(S, SC, xc, xs, ~c);
}


Lit Adder::FA_carry ( Solver *S, CaDiCaL::Solver * SC, Lit a, Lit b, Lit c ) {
  
  Lit x = mkLit(S->newVar(), false);

  addTernaryClause(S, SC, b, c, ~x);
  addTernaryClause(S, SC, a, c, ~x);
  addTernaryClause(S, SC, a, b, ~x);

  addTernaryClause(S, SC, ~b, ~c, x);
  addTernaryClause(S, SC, ~a, ~c, x);
  addTernaryClause(S, SC, ~a, ~b, x);

  return x;
}

Lit Adder::FA_sum ( Solver *S, CaDiCaL::Solver * SC, Lit a, Lit b, Lit c )
{
    Lit x = mkLit(S->newVar(), false);

    addQuaternaryClause(S, SC, a, b, c, ~x);
    addQuaternaryClause(S, SC, a, ~b, ~c, ~x);
    addQuaternaryClause(S, SC, ~a, b, ~c, ~x);
    addQuaternaryClause(S, SC, ~a, ~b, c, ~x);

    addQuaternaryClause(S, SC, ~a, ~b, ~c, x);
    addQuaternaryClause(S, SC, ~a, b, c, x);
    addQuaternaryClause(S, SC, a, ~b, c, x);
    addQuaternaryClause(S, SC, a, b, ~c, x);

    return x;
}

Lit Adder::HA_carry ( Solver *S, CaDiCaL::Solver * SC, Lit a, Lit b) // a AND b
{  
  Lit x = mkLit(S->newVar(), false);

  addBinaryClause(S, SC, a, ~x);
  addBinaryClause(S, SC, b, ~x);
  addTernaryClause(S, SC, ~a, ~b, x);

  return x;
}

Lit Adder::HA_sum ( Solver *S, CaDiCaL::Solver * SC, Lit a, Lit b ) // a XOR b
{
  Lit x = mkLit(S->newVar(), false);

  addTernaryClause(S, SC, ~a, ~b, ~x);
  addTernaryClause(S, SC, a, b, ~x);
  
  addTernaryClause(S, SC, ~a, b, x);
  addTernaryClause(S, SC,  a, ~b, x);

  return x;
}


void Adder::adderTree (Solver *S, CaDiCaL::Solver * SC, std::vector< std::queue< Lit > > & buckets, vec< Lit >& result ) {
  Lit x,y,z;
  Lit u = lit_Undef;

  for ( int i = 0; i < buckets.size(); i++ ) {
      if ( buckets[i].size() == 0 )
    continue;

      if ( i == buckets.size() - 1 && buckets[i].size() >= 2 ) {
    buckets.push_back ( std::queue<Lit>() );
    result.push ( u );
    }

      while ( buckets[i].size() >= 3 ) {
    x = buckets[i].front();
    buckets[i].pop();
    y = buckets[i].front();
    buckets[i].pop();
    z = buckets[i].front();
    buckets[i].pop();
    Lit xs = FA_sum ( S, SC, x,y,z );
    Lit xc = FA_carry ( S, SC, x,y,z );
    buckets[i  ].push ( xs );
    buckets[i+1].push ( xc );
    FA_extra(S, SC, xc, xs, x, y, z);
    }

      if ( buckets[i].size() == 2 ) {
    x = buckets[i].front();
    buckets[i].pop();
    y = buckets[i].front();
    buckets[i].pop();
    buckets[i  ].push ( HA_sum ( S, SC,  x,y ) );
    buckets[i+1].push ( HA_carry ( S, SC, x,y ) );
    }


      result[i] = buckets[i].front();
      buckets[i].pop();
      }

  }

  // Generates clauses for “xs <= ys”, assuming ys has only constant signals (0 or 1).
// xs and ys must have the same size

void Adder::lessThanOrEqual (Solver *S, CaDiCaL::Solver * SC, vec< Lit > & xs, std::vector< uint64_t > & ys) {
  assert ( xs.size() == ys.size() );
  vec<Lit> clause;
  bool skip;
  for ( int i = 0; i < xs.size(); ++i ) {
      if ( ys[i] == 1 || xs[i] == lit_Undef )
    continue;
      
      clause.clear();

      skip = false;

      for ( int j = i + 1; j < xs.size(); ++j )
      {
    if ( ys[j] == 1 )
    {
        if ( xs[j] == lit_Undef )
        {
      skip = true;
      break;
        }

        clause.push ( ~xs[j] );
    }
    else
    {
        assert ( ys[j] == 0 );

        if ( xs[j] == lit_Undef )
      continue;

        clause.push ( xs[j] );
    }
      }

      if ( skip )
    continue;

      clause.push ( ~xs[i] );

      //formula.addClause( clause );
      S->addClause(clause);
      ICadical::addClause(SC, clause);
      }

}

void Adder::lessThanOrEqualInc (Solver *S, CaDiCaL::Solver * SC, vec< Lit > & xs, std::vector< uint64_t > & ys, vec<Lit>& assumptions) {
  assert ( xs.size() == ys.size() );
  vec<Lit> clause;
  bool skip;
  for ( int i = 0; i < xs.size(); ++i ) {
      if ( ys[i] == 1 || xs[i] == lit_Undef )
    continue;
      
      clause.clear();

      skip = false;

      for ( int j = i + 1; j < xs.size(); ++j )
      {
    if ( ys[j] == 1 )
    {
        if ( xs[j] == lit_Undef )
        {
      skip = true;
      break;
        }

        clause.push ( ~xs[j] );
    }
    else
    {
        assert ( ys[j] == 0 );

        if ( xs[j] == lit_Undef )
      continue;

        clause.push ( xs[j] );
    }
      }

      if ( skip )
    continue;

      clause.push ( ~xs[i] );

      //formula.addClause( clause );
      Lit t = mkLit(S->newVar(), false);
      clause.push(t);
      assumptions.push(~t);
      S->addClause(clause);
      ICadical::addClause(SC, clause);
      }

}

void Adder::numToBits ( std::vector<uint64_t> & bits, uint64_t n, uint64_t number ) {
    bits.clear();

  
  for ( int64_t i = n - 1; i >= 0; --i ) {
      int64_t tmp = ((int64_t)1) << i;
      if ( number < tmp ) {
    bits.push_back ( 0 );
    }
      else {
    bits.push_back ( 1 );
    number -= tmp;
    }
      }

  reverse ( bits.begin(), bits.end() );
}

void Adder::encode(Solver *S, CaDiCaL::Solver * SC, vec<Lit> &lits, vec<uint64_t> &coeffs, uint64_t rhs){

    _output.clear();

    uint64_t nb = ld64(rhs); // number of bits
    Lit u = lit_Undef;

    for ( int iBit = 0; iBit < nb; ++iBit ) {
        _buckets.push_back ( std::queue<Lit>() );
        _output.push ( u );
        for ( int iVar = 0; iVar < lits.size(); ++iVar ) {
            if ( ( ( ((int64_t)1) << iBit ) & coeffs[iVar] ) != 0 )
                _buckets.back().push ( lits[iVar] );
            }
        }

    std::vector<uint64_t> kBits;

    adderTree (S, SC, _buckets, _output);
  
    numToBits (kBits, _buckets.size(), rhs );

    lessThanOrEqual (S, SC, _output, kBits);
    hasEncoding = true;
}

void Adder::encodeInc(Solver *S, CaDiCaL::Solver * SC, vec<Lit> &lits, vec<uint64_t> &coeffs, uint64_t rhs, vec<Lit> &assumptions){
    _output.clear();

    uint64_t nb = ld64(rhs); // number of bits
    Lit u = lit_Undef;

    for ( int iBit = 0; iBit < nb; ++iBit ) {
        _buckets.push_back ( std::queue<Lit>() );
        _output.push ( u );
        for ( int iVar = 0; iVar < lits.size(); ++iVar ) {
            if ( ( ( ((int64_t)1) << iBit ) & coeffs[iVar] ) != 0 )
                _buckets.back().push ( lits[iVar] );
            }
        }

    std::vector<uint64_t> kBits;

    adderTree (S, SC, _buckets, _output);
    numToBits (kBits, _buckets.size(), rhs );

    lessThanOrEqualInc (S, SC, _output, kBits, assumptions);
    hasEncoding = true;

}

void Adder::updateInc(Solver *S, CaDiCaL::Solver * SC, uint64_t rhs, vec<Lit>& assumptions){
      
      std::vector<uint64_t> kBits;
      numToBits (kBits, _buckets.size(), rhs );
      lessThanOrEqualInc (S, SC,  _output, kBits, assumptions);
}

void Adder::update(Solver *S,CaDiCaL::Solver * SC, uint64_t rhs){
      
      std::vector<uint64_t> kBits;
      numToBits (kBits, _buckets.size(), rhs );
      lessThanOrEqual (S, SC,  _output, kBits);
}


uint64_t Adder::ld64(const uint64_t x)
{
  return (sizeof(uint64_t) << 3) - __builtin_clzll (x);
//   cout << "x " << x << endl;
//   int ldretutn = 0;
//   for (int i = 0; i < 63; ++i)
//   {
//     if ((x & (1 << i)) > 0)
//     {
//       cout << "ldretutn " << ldretutn << endl;
//       ldretutn = i + 1;
//     }
//   }
//   
//   return ldretutn;
}












