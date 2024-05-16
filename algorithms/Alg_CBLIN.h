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

#ifndef Alg_CBLIN_h
#define Alg_CBLIN_h

#ifdef SIMP
#include "simp/SimpSolver.h"
#else
#include "core/Solver.h"
#endif

#include "../Encoder.h"
#include "../MaxSAT.h"
#include "../MaxTypes.h"
#include "Alg_NUWLS.h"
#include "../rustsat/capi/rustsat.h"
#include "utils/System.h"
#include <map>
#include <set>
#include <utility>
#include <iostream>
#include <stdio.h>


namespace openwbo {

class CBLIN : public MaxSAT {

public:
  // NOTE: currently the encoding is not set as an input parameter.
  CBLIN(int verb = _VERBOSITY_MINIMAL_, int weight = _WEIGHT_NORMAL_, 
        int linear = 0, bool delsol = false, 
        int gcLim = -1, bool r2strat = false, bool incrementalV = false, 
        bool reconstruct_sol_ = false, bool minimize_sol_ = true, int m_strat = 0, bool use_dpw = false, 
        bool dpw_coarse_ = false, bool dpw_inc_ = false, bool extend_models_ = true, bool local_s = false, uint64_t _non_inc_precision = 10 , 
        bool _harden_in_SIS = false) {
    
    solver = NULL;
    verbosity = verb;

    nbCurrentSoft = 0;
    weightStrategy = weight;
    num_hardened = 0;
    clauses_added = 0;
    softs_added = 0;
    vars_added = 0;
    maxw_nothardened = 0;
    
    lins = linear;
    enc = NULL;
    did_harden = false;
    known_gap = UINT64_MAX;
    timeLimitCores = gcLim;
    relaxBeforeStrat = r2strat;
    reconstruct_sol = reconstruct_sol_;
    reconstruct_iter = false;
    incrementalVarres = incrementalV;

    inLinSearch = false;
    wrong_eval_cg = 0;
    wrong_eval_lin = 0;
    init_rhs = 0;
    delete_before_lin = delsol;

    minimize_sol = minimize_sol_;
    minimize_strat = m_strat;

    extend_models = extend_models_;
    max_weight_after_cg = 0;
    dpw = NULL;
    use_DPW = use_dpw;
    dpw_coarse = dpw_coarse_;
    dpw_fine_convergence_after = false;
    incremental_DPW = dpw_inc_;
    have_encoded_precision = false;
    weight_map_setup = false;
    if (incremental_DPW) {
      assert(use_DPW);
    }

    use_local_search = local_s;
    skip_local_search = false;
    harden_in_SIS = _harden_in_SIS;

    non_inc_precision = _non_inc_precision;


  }

  ~CBLIN() {
    //Formula is deleted in MaxSAT.cc
    if (solver != NULL)
      delete solver;
  }

  StatusCode search(); // WBO search.

protected:
  // Rebuild MaxSAT solver
  //
  // Rebuild MaxSAT solver with weight-based strategy.
  Solver *updateSolver();
  int clauses_added;
  int softs_added;
  int vars_added;
  int lins;
  bool delete_before_lin;
  int timeLimitCores; //-1 = noBudget;
  bool relaxBeforeStrat;

  void softsSatisfied();
  void updateCurrentWeight(int strategy); // Updates 'currentWeight'.
  uint64_t
  findNextWeight(uint64_t weight); // Finds the next weight for 'currentWeight'.
  uint64_t
  findNextWeightDiversity(uint64_t weight); // Finds the next weight for
                                            // 'currentWeight' using diversify
                                            // heuristic.

  
  // Utils for core management
  //
  void encodeMaxRes(vec<Lit> &core, uint64_t weightCore); // Encodes exactly one constraint.
  void relaxCore(vec<Lit> &conflict, uint64_t weightCore);            // Relaxes a core.
  uint64_t computeCostCore(const vec<Lit> &conflict); // Computes the cost of a core.
  void setAssumptions(vec<Lit> &assumps);
  int num_hardened;

  Solver * hardenClauses();
  bool harden_in_SIS; 
  void hardenClausesSIS(uint64_t reduced_cost, vec<lbool> &currentModel);
  Solver * resetSolver();
  uint64_t maxw_nothardened;
  uint64_t max_coeff_nothardened_sis;

  uint64_t known_gap;

  uint64_t init_rhs;

  void checkGap();
  bool inLinSearch;

  //Varying Resolutio
  bool weight_map_setup;
  bool incrementalVarres;
  uint64_t get_Maximum_Weight();
  void update_SIS_precision();
  int  moreThanWeight(uint64_t w);
  void init_SIS_precision();
  void initializePBConstraint(uint64_t rhs);

  void updateBoundLinSearch (uint64_t newBound);

  bool checkModel(bool from_local_search = false);
  uint64_t computeCostReducedWeights(vec<lbool> &currentModel);

  void setPBencodings();
  Encoder * enc;
  
  ///DPW
  struct SolverWithBuffer {
    Solver *solver_b;
    vec<Lit> buffer;
    int clauses_added;
    int verbosity;
};

  RustSAT::DynamicPolyWatchdog *dpw;
  bool use_DPW;
  bool have_encoded_precision;
  bool dpw_coarse;
  bool dpw_fine_convergence_after;
  uint64_t fine_bound;
  uint64_t dpw_next_precision();
  void dpw_encode_and_enforce(uint64_t rhs);
  static void dpw_assumps(int lit, void *assumps);
  static void dpw_clause_collector(int lit, void *ptr);

  bool incremental_DPW;

  ///DPW
  
  
  vec<lbool> bestModel;
  void flipValueinBest(Lit l);

  void extendBestModel();
  void setCardVars(bool prepro);
  
  vec<bool> isSoft; 

  std::string print_timeSinceStart();
  time_t timeSinceStart();
  time_t timeSincePrepro();

  uint64_t precision_factors();
  uint64_t non_inc_precision;

  // Core guided division
  std::vector<uint64_t> reducedWeight;
  StatusCode weightDisjointCoresDivision();


  StatusCode unsatSearch();  // Search using only hard clauses.
  StatusCode weightSearch(); // Search using weight-based methods.
 
  StatusCode setup(); // unsat search and other setups
  StatusCode coreGuidedLinearSearch();
  uint64_t max_weight_after_cg;
  int exponent(uint64_t weight);
  uint64_t raise_to(int exponent);
  std::vector<int> coeff_counter;
  void set_up_objective_counter(uint64_t init);
  //These are subroutines in other searches and should not be 
  StatusCode linearSearch();
  StatusCode weightDisjointCores(); // LB phase
  void build_objective_func_and_coeffs();
  vec<Lit> objFunction; // Literals to be used in the constraint that excludes
                        // models.
  vec<uint64_t> coeffs; // Coefficients of the literals that are used in the
                        // constraint that excludes models.
 
  //DEBUGGING
  vec<Lit> objFunction_;
  vec<uint64_t> coeffs_;
  uint64_t rhs_;
  uint64_t num_literals_;



  vec<Lit> minimisable_lits;
 
  void savePhase();
  time_t time_start;
  time_t time_prepro;
	time_t time_best_solution;

  // Other
  // Initializes assumptions and core extraction.
  void initAssumptions();

  void printProgress();

  
  void addSoftClauseAndAssumptionVar(uint64_t weight, vec<Lit> &clause);
  uint64_t computeCostOfModel(vec<lbool> &currentModel);

  int nRealSoft();
  bool shouldUpdate();
  bool did_harden;

  // SAT solver
  Solver *solver;  // SAT solver used as a black box.

  // Variables used  in 'weightSearch'
  //
  int nbCurrentSoft;  // Current number of soft clauses used by the MaxSAT
                      // solver.
  int weightStrategy; // Weight strategy to be used in 'weightSearch'.

  // Core extraction
  //
  std::map<Lit, int> coreMapping; // Maps the assumption literal to the number
                                  // of the soft clause.
  vec<Lit> assumptions; // Stores the assumptions to be used in the extraction
                        // of the core.

  StatusCode getModelAfterCG();

 int wrong_eval_cg;
 int wrong_eval_lin;
  
 bool reconstruct_sol; 
 bool reconstruct_iter;
 bool minimize_sol;
 int  minimize_strat;
 void minimizelinearsolution( vec<lbool> & sol);
 bool use_local_search;
 bool skip_local_search;
 void localsearch(vec<lbool> & sol);

  bool extend_models;

 void test_pb_enc();

};
} // namespace openwbo

#endif
