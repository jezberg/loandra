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
#include "../rustsat/capi/rustsat.h"
#include "utils/System.h"
#include "BouMS/BouMS.h"
#include "BouMS/cores.h"
#include "BouMS/wcnf.h"
#include "BouMS/wcnf_util.h"
#include <map>
#include <set>
#include <utility>
#include <iostream>
#include <stdio.h>
#include <functional>


namespace openwbo {

class CBLIN : public MaxSAT {

public:
  // NOTE: currently the encoding is not set as an input parameter.
  CBLIN(int verb = _VERBOSITY_MINIMAL_, int weight = _WEIGHT_DIVERSIFY_, 
        int linear = 0, bool delsol = false, 
        int gcLim = -1, bool r2strat = false, bool incrementalV = false, 
        bool reconstruct_sol_ = false, bool minimize_sol_ = true, int m_strat = 0, bool use_dpw = false, 
        bool dpw_coarse_ = false, bool dpw_inc_ = false, bool extend_models_ = true, bool local_s = false, uint64_t _non_inc_precision = 10 , 
        bool _harden_in_SIS = false, int ls_init_level_ = 0, bool ls_dyn_prec_ = false, bool ls_sis_ = false,
        bool ls_merge_assign_ = false, bool ls_min_ = false, bool ls_cores_ = false, double zero_weight_core_fact_ = 3)
  {
    
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
    if (use_local_search) {
      minimize_sol = true;
      ls_init_level = ls_init_level_;
      ls_dyn_prec = ls_dyn_prec_;
      ls_sis = ls_sis_;
      ls_merge_assign = ls_merge_assign_;
      ls_min = ls_min_;
      ls_cores = ls_cores_;
      zero_weight_core_fact = zero_weight_core_fact_;
    }

    skip_local_search = false;
    harden_in_SIS = _harden_in_SIS;

    non_inc_precision = _non_inc_precision;


  }

  ~CBLIN() {
    //Formula is deleted in MaxSAT.cc
    if (solver != NULL)
      delete solver;

    if (boums_cores) {
      BouMS_cores_free(&boums_inst, boums_cores, free);
      delete boums_cores;
      boums_cores = NULL;
    }
    if (core_var_to_clause) {
      delete[] core_var_to_clause;
      core_var_to_clause = NULL;
    }
    BouMS_wcnf_util_deleteFormula(&boums_inst, free, NULL);
    if (boums_mem) {
      free(boums_mem);
      boums_mem = NULL;
    }
    if (boums_assignment) {
      delete[] boums_assignment;
      boums_assignment = NULL;
    }
    if (boums_clause_map.ex2In) {
      delete[] boums_clause_map.ex2In;
      boums_clause_map.ex2In = NULL;
    }
    if (boums_clause_map.in2Ex) {
      delete[] boums_clause_map.in2Ex;
      boums_clause_map.in2Ex = NULL;
    }
    if (orig_maxsat_formula) {
      delete orig_maxsat_formula;
      orig_maxsat_formula = NULL;
    }
    if (init_ls_ub_assign) {
      delete[] init_ls_ub_assign;
      init_ls_ub_assign = NULL;
    }
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
  void harden_incremental();
  void initializePBConstraint(uint64_t rhs);

  void updateBoundLinSearch (uint64_t newBound);

  bool checkModel(bool from_local_search = false, bool improve_better = false);

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
 bool localsearch(vec<lbool> & sol);

  bool extend_models;

  // BEGIN LS w/ BouMS
  int ls_init_level = 0; // 0=disabled, 1=only on preprocessed, 2=on preprocessed then on original
  bool ls_dyn_prec = false; // run LS in dynamic resolution, more precisely in initializePBconstraint
  bool ls_sis = false; // run LS in SIS
  bool ls_merge_assign = false; // use assignment merging for LS in dyn prec and sis
  bool ls_min = false; // run LS for solution minimization
  bool ls_cores = false; // use cores in LS
  double zero_weight_core_fact = 3;
  MaxSATFormula* orig_maxsat_formula = NULL;
  uint64_t init_ls_ub = UINT64_MAX;
  bool* init_ls_ub_assign = NULL;
  BouMS_wcnf_t boums_inst = BouMS_wcnf_util_newFormula();
  BouMS_memoryReq_t boums_mem_req;
  BouMS_uint_t boums_bytes = 0;
  void* boums_mem = NULL;
  BouMS_params_t boums_params;
  bool* boums_assignment = NULL;
  BouMS_clauseMap_t boums_clause_map = { .ex2In = NULL, .in2Ex = NULL };
  bool boums_broken = false;
  vec<lbool> ls_merged_assign;
  vec<lbool>* ls_usual_init_assign = &bestModel;
  uint64_t old_sis_precision;
  std::vector<BouMS_cores_core_t> cores;
  BouMS_cores_mem_t* boums_cores = NULL;
  BouMS_uint_t* core_var_to_clause = NULL;
  void updateBouMSInstance(); // return true in case of error
  template<typename V, typename v> unsigned int mergeAssignments(vec<lbool>& dst, const V& src, const std::function<lbool(const v)>& tolbool) {
    unsigned int num_disagreements = 0;
    for (unsigned int vIdx = 0; vIdx < dst.size(); ++vIdx) {
      if (dst[vIdx] != tolbool(src[vIdx])) {
        dst[vIdx] = rand() % 2 ? l_True : l_False;
        ++num_disagreements;
      }
    }
    return num_disagreements;
  }
  static inline lbool tolbool(bool b) {
    return b ? l_True : l_False;
  }
  virtual void loadFormula(MaxSATFormula *maxsat) override;
  virtual void setup_formula() override;
  virtual void printAnswer(int) override;
  // END LS w/ BouMS

};
} // namespace openwbo

#endif
