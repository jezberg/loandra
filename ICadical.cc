/*!
 * \author Jeremias Berg - jeremiasberg@hmail.com
 *
 * @section LICENSE
 *   Loandra, Copyright (c) 2018, Jeremias Berg, Emir Demirovic, Peter Stuckey
 *  Open-WBO, Copyright (c) 2013-2017, Ruben Martins, Vasco Manquinho, Ines Lynce
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

#include "ICadical.h"

using namespace openwbo;

lbool ICadical::searchSATSolver(CaDiCaL::Solver * solver, vec<Lit> & assumptions) {
    for (int i = 0; i < assumptions.size(); i++) {
        Lit l = assumptions[i];
        solver->assume(lit2Int(l));
    }
    int res = solver->solve();
    if (res == 20) {
        return l_False;
    }
    else if (res == 10) {
        return l_True;
    }
    else {
        assert(res == 0);
        return l_Undef;
    }
}

void ICadical::addClause(CaDiCaL::Solver * solver, vec<Lit> & clause) {
    for (int i = 0; i < clause.size(); i++) {
        Lit l = clause[i];
        solver->add(lit2Int(l));
    }
    solver->add(0);
}

void ICadical::getCore(CaDiCaL::Solver * solver, vec<Lit> & assumptions, vec<Lit>  & core_out) {
    assert(core_out.size() == 0);
    for (int i = 0; i < assumptions.size(); i ++) {
        Lit l = assumptions[i];
        if (solver->failed(lit2Int(l))) {
            core_out.push(~l);
        }
    }
}

void ICadical::getModel(CaDiCaL::Solver * solver, vec<lbool> & model_out) {
    assert(model_out.size() == 0);
    for (int i = 1; i <= solver->vars(); i++) {
        int v = solver->val(i);
        if (v == i) model_out.push(l_True);
        else if (v == (-1)*i) model_out.push(l_False);
        else model_out.push(l_Undef);
    }
}

CaDiCaL::Solver* ICadical::newSATSolver() {
    return new CaDiCaL::Solver;
}



int ICadical::lit2Int(Lit l) {
	if (sign(l)) {
		return  -(var(l) + 1);
	}
	else {
		return var(l) + 1; 
	}
}

Lit ICadical::int2Lit(int l) {
	int var = abs(l) - 1;
	bool sign = l > 0;
	return sign ? mkLit(var) : ~mkLit(var);
}
