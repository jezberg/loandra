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
#include "../SISPropagator.h"
#include "../MaxTypes.h"
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
  CBLIN(int verb = _VERBOSITY_MINIMAL_, int weight = _WEIGHT_DIVERSIFY_, 
        int linear = 0, bool delsol = false, 
        int gcLim = -1, bool r2strat = false, bool incrementalV = false, 
        bool reconstruct_sol_ = false, bool minimize_sol_ = true, int m_strat = 0,
        bool dpw_coarse_ = false, bool dpw_inc_ = false, bool extend_models_ = true, bool local_s = false, uint64_t _non_inc_precision = 10 , 
        bool _harden_in_SIS = false, bool opt_phase_save = false, bool _sis_in_propagator = false) {
    
    use_propagator = _sis_in_propagator; 
    solverCad = NULL;
    timer = new Timer(gcLim);
    has_flipped = false;
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
    relaxBeforeStrat = r2strat;
    reconstruct_sol = reconstruct_sol_;
    reconstruct_iter = false;

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
    dpw_coarse = dpw_coarse_;
    dpw_fine_convergence_after = false;
    incremental_DPW = dpw_inc_;
    have_encoded_precision = false;
    weight_map_setup = false;


    use_local_search = local_s;
    if (use_local_search)
      minimize_sol = true;

    skip_local_search = false;
    harden_in_SIS = _harden_in_SIS;

    non_inc_precision = _non_inc_precision;
    optimistic = opt_phase_save;

  }

  ~CBLIN() {
    //Formula is deleted in MaxSAT.cc
    if (solverCad != NULL)
      delete solverCad;
  }

  StatusCode search(); // WBO search.

protected:
  // Rebuild MaxSAT solver
  //
  // Rebuild MaxSAT solver with weight-based strategy.
  void updateSolver();
  int clauses_added;
  int softs_added;
  int vars_added;
  int lins;
  bool delete_before_lin; //-1 = noBudget;
  bool relaxBeforeStrat;

  void softsSatisfied();
  void clearFixingsonSoft();
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

  void hardenClauses();
  bool harden_in_SIS; 
  void hardenClausesSIS(uint64_t reduced_cost);
  void resetSolver();
  uint64_t maxw_nothardened;
  uint64_t max_coeff_nothardened_sis;

  uint64_t known_gap;

  uint64_t init_rhs;

  void checkGap();
  bool inLinSearch;
  void update_current_soft(uint64_t precision);
  //Varying Resolution
  bool weight_map_setup;
  bool incrementalVarres;
  uint64_t get_Maximum_Weight();
  void update_SIS_precision();
  int  moreThanWeight(uint64_t w);
  uint64_t compute_first_precision();
  uint64_t compute_next_SIS_precision(uint64_t current_precision);
  void init_SIS_precision();
  void harden_incremental();
  void initializePBConstraint(uint64_t rhs);

  void updateBoundLinSearch (uint64_t newBound);

  bool checkModel(bool from_local_search = false, bool improve_better = false);

  //auto litValFromCadical() { return [this](Lit l){ return solverCad->val(MaxSAT::lit2Int(l)); }; }
  //auto litValFromModel(vec<Lit> & model) {return [&model, this](Lit l){ return literalTrueInModel(l, model);};  }

  template <typename LitVal>
  uint64_t computeCostReducedWeights(LitVal* lit_true) {
      return computeCostReducedWeights_prec(lit_true, maxsat_formula->getMaximumWeight());
  }

  template <typename LitVal>
  uint64_t computeCostReducedWeights_prec (LitVal* lit_true, uint64_t precision) {
    logPrint("Computing cost of reduced precision");
  
      uint64_t tot_reducedCost = 0;

      for (int i = 0; i < maxsat_formula->nSoft(); i++) {
        assert(maxsat_formula->getSoftClause(i).clause.size() == 1);
        Lit l = maxsat_formula->getSoftClause(i).clause[0];
        if (!(*lit_true)(l)) {
          tot_reducedCost += (maxsat_formula->getSoftClause(i).weight / precision);
        }

      }
      logPrint("reduced cost " , tot_reducedCost, " gap ", known_gap / precision);
      return tot_reducedCost;
  }

  void setPBencodings();
  Encoder * enc;
  


  RustSAT::DynamicPolyWatchdog *dpw;
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

  bool use_propagator;
  StatusCode linearSearch_propagator();
  void set_observed_vars(vec<Lit> & orig_obj);
  uint64_t compute_ub_red_cost(uint64_t precision);


  StatusCode weightDisjointCores(); // LB phase
  void build_objective_func_and_coeffs();
  void build_objective_func_and_coeffs_prop();
  vec<Lit> objFunction; // Literals to be used in the constraint that excludes
                        // models.
  vec<uint64_t> coeffs; // Coefficients of the literals that are used in the
                        // constraint that excludes models.
 
 

  vec<Lit> minimisable_lits;
  bool optimistic;
  void savePhase();
  time_t time_start;
  time_t time_prepro;
	time_t time_best_solution;

  // Other
  // Initializes assumptions and core extraction.
  void initAssumptions();

  void printProgress();

  
  void addSoftClauseAndAssumptionVar(uint64_t weight, vec<Lit> &clause);
  bool literal_sat_in_cadical(Lit l);

  template <typename LitVal>
  uint64_t computeCostOfModel(LitVal* lit_true) { 
    logPrint("Compute cost ");
    if (!do_preprocess) {
        return computeCostOriginalClauses(lit_true);
    }
    if (reconstruct_sol && reconstruct_iter) {
      vec<lbool> model;
      assert(bestModel.size() > 0);
      for (int i = 1; i <= bestModel.size(); i++) {
        Lit l = mkLit(i, true);
        if ((*lit_true)(l)) model.push(l_True);
        else if (!(*lit_true)(l)) model.push(l_False);
        else model.push(l_Undef);
      }
      vec<lbool> reconstructed;
      reconstruct_model_prepro(model, reconstructed); 
      auto lambda = [this, &reconstructed](Lit l){return literalTrueInModel(l, reconstructed);};
      return computeCostOriginalClauses(&lambda);
    }
    else {
      return computeCostObjective(lit_true);
    }
  }
  
  
  bool did_harden;
  int nRealSoft();
  bool shouldUpdate();
  void flipValueinBestModel(Lit l);

  // SAT solver
  Encoder encoder; // Interface for the encoder of constraints to CNF.

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
 void minimizelinearsolution(vec<lbool> & sol);
 
 //TODO refactor these into some other class. 
 bool has_flipped;
 bool flipLiterals();
 void freezeObjective();
 
 //Cadical Related
 CaDiCaL:: Solver* solverCad;
 
 // Timer
 Timer* timer;

 bool use_local_search;
 bool skip_local_search;
 void localsearch(vec<lbool> & sol);

  bool extend_models;


};
} // namespace openwbo

#endif
