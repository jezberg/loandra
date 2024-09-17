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

#include "SISPropagator.h"

using namespace openwbo;

bool SISPropagator::cb_check_found_model (const std::vector<int> &model) {
  //TODO compute the reduced cost, enforce a new bound based on it. 
  //TODO check best solution cost. 
  // update sol
  //update next bound rustSAT

  check_sol_global_cost(model);
  uint64_t reduced_cost = compute_cost_sol(precision, model);
  assert ( reduced_cost < precision_ub);
  if (reduced_cost == 0) {
    return true;
  }
  uint64_t coarse_ub = RustSAT::dpw_coarse_ub(dpw, reduced_cost);
  int  num_vars = -1; // TODO figure out how to get the current number of variables

  RustSAT::dpw_encode_ub(dpw, coarse_ub, coarse_ub, &num_vars, &dpw_clause_collector, static_cast<void*>(&clauses_to_add));
  std::vector<int> assumps;
  RustSAT::MaybeError ret = RustSAT::dpw_enforce_ub(dpw, coarse_ub, &dpw_assumps, static_cast<void*>(&assumps));
  if (ret == RustSAT::MaybeError::NotEncoded) {
      ///TODO handle error 
  }
  assert(ret == RustSAT::MaybeError::Ok);
  assert(assumps.size() > 0);
  clauses_to_add.push_front(assumps[0]);
  clauses_to_add.push_front(0);

  return false; 
}


int SISPropagator::cb_add_external_clause_lit () {
  int lit = clauses_to_add.back();
  clauses_to_add.pop_back();
  return lit;
};

void SISPropagator::build_objective_function(const vec<Lit> & objFunction, const vec<uint64_t> & coeffs) {
  assert(objFunction.size() == coeffs.size());
  for (int i = 0; i < objFunction.size(); i++) {
    int lit = lit2Int(objFunction[i]);
    objective_coefficients[abs(lit)] = coeffs[i]; 
    objective_lit_incurs_cost[abs(lit)] = lit > 0; 
  }
}

bool SISPropagator::literal_incurs_cost_in_objective(int lit) {
  return (lit > 0) == objective_lit_incurs_cost[abs(lit)];
}

void SISPropagator::check_sol_global_cost(const  std::vector<int> &model) {
  uint64_t cost = compute_cost_sol(1, model);
  if (cost < ub) {
    ub = cost; 
    best_model.clear();
    for (int i = 0; i < model.size(); i++) {
      best_model.push_back(model[i]);
    }
  }
  return; 
}

uint64_t SISPropagator::compute_cost_sol(uint64_t precision, const  std::vector<int> &model) {
  uint64_t cost = 0;
  for (int i = 0; i < model.size(); i++) {
    int lit = model[i];
    if (literal_incurs_cost_in_objective(lit)) {
        cost += objective_coefficients[abs(lit)] / precision;
    }

  }
  return cost; 
}

void SISPropagator::dpw_clause_collector(int lit, void *ptr) {
  std::deque<int>* queue = static_cast<std::deque<int>*>(ptr);
  queue->push_front(lit);
  return;
}

void SISPropagator::dpw_assumps(int lit, void *assumps) {
  std::vector<int> * asssums_vec = static_cast<std::vector<int>*>(assumps);
  asssums_vec->push_back(lit);
}


 










