/*!
 * \author Jeremias Berg - jeremiasberg@gmail.com
 *  This implementation is heavily based on Open-WBO, thanks to the authors! 
 * 
 * @section LICENSE
 *  Loandra, Copyright (c) 2018 Jeremias Berg
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

#ifndef SIS_Propagate_h
#define SIS_Propagate_h

#include "utils/System.h"
#include "cadical/src/cadical.hpp"
#include "rustsat/capi/rustsat.h"
#include "MaxSAT.h"


namespace openwbo {

class SISPropagator : public CaDiCaL::ExternalPropagator {
//TODO the propagator should be allowed to compute costs of the actual old objective 
public:
  // NOTE: currently the encoding is not set as an input parameter.
  SISPropagator(  uint64_t global_ub, uint64_t _precision_ub, int _num_vars, vec<Lit> & objFunction, const vec<uint64_t> & coeffs, 
                  uint64_t _precision, vec<Lit> & orig_objFunction, const vec<uint64_t> & orig_coeffs, bool _verbosity) {
    precision = _precision;
    verbose = _verbosity;
    ub = global_ub;
    precision_ub = _precision_ub;
    build_objective_function(objFunction, coeffs, orig_objFunction, orig_coeffs);
    clauses_to_add.clear();
    num_vars = _num_vars;
    has_encoded = false;
    dpw = RustSAT::dpw_new();
    sum_coeffs = 0;
    for (int i = 0; i < objective_vars.size(); i++) {
        int lit = objective_vars[i];
        if ((objective_coefficients[abs(lit)] / precision) > 0) {
          RustSAT::dpw_add(dpw, lit, objective_coefficients[abs(lit)] / precision);
          if (verbose) {std::cout << "c adding to  RustSAT: " << lit << " weight " <<  objective_coefficients[abs(lit)] / precision << endl;}
          sum_coeffs +=objective_coefficients[lit] / precision;
        }
    }
    if (verbose) {
      std::cout << "c end of constructor objective_coefficients, size: " << objective_coefficients.size() << " elements:"; 
      for (int i = 0; i < objective_coefficients.size(); i++) {
        std::cout  << " " << i << "/" << objective_coefficients[i]; 
      }
      std::cout << endl;
      std::cout << "c end of constructor orig_objective_coefficients, size: " << orig_objective_coefficients.size() << " elements:"; 
      for (int i = 0; i < orig_objective_coefficients.size(); i++) {
        std::cout  << " " << i << "/" << orig_objective_coefficients[i]; 
      }
      std::cout << endl;
    }   
  }

  ~SISPropagator() {
    clauses_to_add.clear();
    objective_coefficients.clear();
    objective_lit_incurs_cost.clear();
    best_model.clear();
    if (dpw != NULL) {
      RustSAT::dpw_drop(dpw);
      dpw = NULL;
    }
  }

  bool is_lazy = true; // lazy propagator only checks complete assignments
  bool are_reasons_forgettable = false; // Reason external clauses can be deleted
  bool is_tainting = true; 

  void notify_assignment (const std::vector<int>& lits)  {};
  void notify_new_decision_level () {};
  void notify_backtrack (size_t new_level) {};

  // Check by the external propagator the found complete solution (after
  // solution reconstruction). If it returns false, the propagator must
  // provide an external clause during the next callback.
  //
  bool cb_check_found_model (const std::vector<int> &model);

  int cb_decide () { return 0; };

  // TODO: think about if we could propagate something extra here, probably not though... 
  //
  int cb_propagate () { return 0; };

  int cb_add_reason_clause_lit (int propagated_lit) {
    return 0;
  };

  bool cb_has_external_clause (bool& is_forgettable) {is_forgettable = false; return clauses_to_add.size() > 0;};

  // The actual function called to add the external clause.
  //
  int cb_add_external_clause_lit ();

  //this will store the best assignment to the objective variables in terms of the full cost. 
  std::vector<int> best_model;
  uint64_t get_best_cost() {return ub;};
  void resetPropagator(uint64_t _precision, uint64_t _precision_ub, int vars);


protected:
  RustSAT::DynamicPolyWatchdog *dpw;
  bool has_encoded; 
  uint64_t sum_coeffs;
  
  void build_objective_function(const vec<Lit> & objFunction, const vec<uint64_t> & coeffs,
                                const vec<Lit> & orig_objFunction, const vec<uint64_t> & orig_coeffs);
  
  std::vector<int> objective_vars;
  std::vector<uint64_t> objective_coefficients;
  // the entry in the in the i:th field here indicates the polarity that incurs cost 
  std::vector<bool> objective_lit_incurs_cost;
  bool literal_incurs_cost_in_objective(int lit);

  std::vector<uint64_t> orig_objective_coefficients;
  // the entry in the in the i:th field here indicates the polarity that incurs cost 
  std::vector<bool> orig_objective_lit_incurs_cost;

  void check_sol_global_cost(const  std::vector<int> &model);
  uint64_t compute_cost_sol(uint64_t precision, const  std::vector<int> &model) ;
  
  static void dpw_clause_collector(int lit, void *ptr);
  static void dpw_assumps(int lit, void *assumps);

  std::deque<int> clauses_to_add; 

  uint64_t precision;
  uint64_t ub; 
  uint64_t precision_ub;

  int num_vars; 
  bool verbose;

  void store_best_model(); 

  int lit2Int(Lit l) {
	  if (sign(l)) {
		  return  -(var(l) + 1);
	  }
	  else {
		  return var(l) + 1; 
	  }
  }

};
} // namespace openwbo

#endif
