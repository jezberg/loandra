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
  if (verbose) {
    std::cout << "c cb_check_found_model size " << model.size() << std::endl;
    std::cout << "c model:";
    for (int i : model) {
      std::cout << " " << i; 
    } 
    std::cout << endl;
  }
  check_sol_global_cost(model);
  uint64_t reduced_cost = compute_cost_sol(precision, model);

  if (verbose) {std::cout << "c reduced_cost " << reduced_cost << endl;}
  
  if (reduced_cost == 0) {
    return true;
  }
  if (!has_encoded ) {
    if (precision_ub < reduced_cost) {
      reduced_cost = precision_ub;
    }
    if (reduced_cost > sum_coeffs) {
      reduced_cost = sum_coeffs;
    }
    if (verbose) {cout << "c extra call to encode with bound " << reduced_cost -1 << endl;}
    RustSAT::dpw_encode_ub(dpw, reduced_cost-1, reduced_cost-1, &num_vars, &dpw_clause_collector, static_cast<void*>(&clauses_to_add));
  }

  if (verbose) {cout << "c call to coarse " << reduced_cost -1 << endl;}
  uint64_t coarse_ub = RustSAT::dpw_coarse_ub(dpw, reduced_cost - 1);
  if (verbose) {cout << "c coarse ub " << coarse_ub << endl;}
  RustSAT::dpw_encode_ub(dpw, coarse_ub, coarse_ub, &num_vars, &dpw_clause_collector, static_cast<void*>(&clauses_to_add));
  std::vector<int> assumps;
  assumps.clear();
  if (verbose) {cout << "c call to dpw_enforce_ub: coarse " << coarse_ub << endl;}
  RustSAT::MaybeError ret = RustSAT::dpw_enforce_ub(dpw, coarse_ub, &dpw_assumps, static_cast<void*>(&assumps));
  if (ret == RustSAT::MaybeError::NotEncoded) {
      std::cout << "c Rustsat returned error " << endl;
      exit(1);
  }
  assert(ret == RustSAT::MaybeError::Ok);
  if (assumps.size() > 1) {
    if (verbose) {std::cout << "c got more than 1 assumption, done " << endl;}
    return true;
  }
  if (verbose) {std::cout << "c adding the unit " << assumps[0] << " size of vector " << assumps.size() << endl;}
  clauses_to_add.push_front(assumps[0]);
  clauses_to_add.push_front(0);
  has_encoded = true;
  return false; 
}


int SISPropagator::cb_add_external_clause_lit () {
  int lit = clauses_to_add.back();
  clauses_to_add.pop_back();
  return lit;
};

void SISPropagator::resetPropagator(uint64_t _precision, uint64_t _precision_ub, int vars) {
  if (verbose) {cout << "c resetting propagator to precision " << _precision << " ub " << _precision_ub << endl;}
  has_encoded = false; 
  precision = _precision;
  num_vars = vars;
  precision_ub = _precision_ub;
  if (dpw != NULL) {
      RustSAT::dpw_drop(dpw);
      dpw = NULL;
  }
  dpw = RustSAT::dpw_new();
  sum_coeffs = 0;
  for (int i = 0; i < objective_vars.size(); i++) {
        int lit = objective_vars[i];
        if ((objective_coefficients[abs(lit)] / precision) > 0) {
          RustSAT::dpw_add(dpw, lit, objective_coefficients[abs(lit)] / precision);
          sum_coeffs += objective_coefficients[lit] / precision;
          if (verbose) {std::cout << "c Adding to  RustSAT: " << lit << " weight " <<  objective_coefficients[abs(lit)] / precision << endl;}
        }
  }
}

void SISPropagator::build_objective_function( const vec<Lit> & objFunction, const vec<uint64_t> & coeffs,
                                              const vec<Lit> & orig_objFunction, const vec<uint64_t> & orig_coeffs) {
  assert(objFunction.size() == coeffs.size());
  for (int i = 0; i < objFunction.size(); i++) {
    int lit = lit2Int(objFunction[i]);
    objective_vars.push_back(lit);

    while (objective_coefficients.size() < abs(lit) +1) objective_coefficients.push_back(0);
    while (objective_lit_incurs_cost.size() < abs(lit) +1) objective_lit_incurs_cost.push_back(false);
    objective_coefficients[abs(lit)] = coeffs[i];
    objective_lit_incurs_cost[abs(lit)] = lit > 0;
  }

  assert(orig_objFunction.size() == orig_coeffs.size());
  for (int i = 0; i < orig_objFunction.size(); i++) {
    int lit = lit2Int(orig_objFunction[i]);
    while (orig_objective_coefficients.size() < abs(lit) +1) orig_objective_coefficients.push_back(0);
    while (orig_objective_lit_incurs_cost.size() < abs(lit) +1) orig_objective_lit_incurs_cost.push_back(false);
    orig_objective_coefficients[abs(lit)] = orig_coeffs[i];
    orig_objective_lit_incurs_cost[abs(lit)] = lit > 0;
  }
  if (verbose) {
    cout << "c end of build objective_vars:";
    for (int i = 0; i < objective_vars.size(); i++) {
      cout << " " << objective_vars[i];
    } 
    cout << endl;

  }
}

bool SISPropagator::literal_incurs_cost_in_objective(int lit) {
  return (lit > 0) == objective_lit_incurs_cost[abs(lit)];
}

void SISPropagator::check_sol_global_cost(const  std::vector<int> &model) {
  uint64_t cost = 0;
  for (int i = 0; i < model.size(); i++) {
    int lit = model[i];
    if ((lit > 0) == orig_objective_lit_incurs_cost[abs(lit)]) {
        uint64_t cost_to_add = abs(lit) < orig_objective_coefficients.size() ? orig_objective_coefficients[abs(lit)] : 0;
        cost += cost_to_add;
    }
  }

  if (verbose) {std::cout << "c cost of sol " << cost << " ub " << ub << std::endl;}
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
  if (verbose) {std::cout << "c computing cost with precision " << precision << endl;}
  uint64_t cost = 0;
  for (int i = 0; i < model.size(); i++) {
    int lit = model[i];
    if (literal_incurs_cost_in_objective(lit)) {
      uint64_t cost_to_add = abs(lit) < objective_coefficients.size() ? objective_coefficients[abs(lit)] / precision : 0;
      cost += cost_to_add;
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


 










