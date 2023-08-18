/*!
 * \author Jeremias Berg - jeremiasberg@gmail.com
 * 
 * @section LICENSE
 *  Loandra, Copyright (c) 2018 Jeremias Berg
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

#ifndef ICadical_h
#define ICadical_h

#ifdef SIMP
#include "simp/SimpSolver.h"
#else
#include "core/Solver.h"
#endif


#include "MaxTypes.h"
#include "utils/System.h"
#include "cadical/src/cadical.hpp"

using NSPACE::vec;
using NSPACE::Lit;
using NSPACE::lit_Undef;
using NSPACE::mkLit;
using NSPACE::lbool;

// TODO: refactor this class to maintain the cadical solver.... 

namespace openwbo {

class ICadical {
    public:
        static lbool searchSATSolver(CaDiCaL::Solver * solver, vec<Lit> & assumptions);
        static void addClause(CaDiCaL::Solver * solver, vec<Lit> & clause);
        static void getCore(CaDiCaL::Solver * solver, vec<Lit> & assumptions, vec<Lit> & core_out);
        static void getModel(CaDiCaL::Solver * solver, vec<lbool> & model_out);
        static CaDiCaL::Solver * newSATSolver();
          

    protected:
        static int lit2Int(Lit l);
        static Lit int2Lit(int l);


};
} // namespace openwbo

#endif
