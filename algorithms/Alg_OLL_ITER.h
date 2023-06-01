/*!
 * \author Ruben Martins - ruben@sat.inesc-id.pt
 *
 * @section LICENSE
 *
 * Open-WBO, Copyright (c) 2013-2017, Ruben Martins, Vasco Manquinho, Ines Lynce
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

#ifndef Alg_OLL_ITER_h
#define Alg_OLL_ITER_h

#ifdef SIMP
#include "simp/SimpSolver.h"
#else
#include "core/Solver.h"
#endif

#include "../Encoder.h"
#include "../MaxSAT.h"
#include <map>
#include <set>
//#include <unorderedset>
#include <iostream>


//preprocessing
#include "../maxpre2/src/preprocessorinterface.hpp"



namespace openwbo {

//=================================================================================================
class OLL_ITER : public MaxSAT {

public:
  OLL_ITER( int verb = _VERBOSITY_MINIMAL_, int enc = _CARD_TOTALIZER_, bool pre = false, bool fix = false) 
  {
    solver = NULL;
    verbosity = verb;
    incremental_strategy = _INCREMENTAL_ITERATIVE_;
    encoding = enc;
    min_weight = 1;
    num_hardened = 0;
    nonreformulatedCores = 0;
    num_card_dropped = 0;    
    do_preprocess = pre;
    weightRemoved = 0;
    use_reconstruct = fix;
    num_hardened_me = 0;
    nOrigVars  = 0;
    notHardened_ind.clear();
  }
  ~OLL_ITER() {
    if (solver != NULL)
      delete solver;
  }

  StatusCode search();

  // Print solver configuration.
  void printConfiguration() {

    if(!print) return;

    printf("c ==========================================[ Solver Settings "
           "]============================================\n");
    printf("c |                                                                "
           "                                       |\n");
    printf("c |  Algorithm: %23s                                             "
           "                      |\n",
           "OLL-ITERATIVE");
    print_Card_configuration(encoding);
    printf("c |                                                                "
           "                                       |\n");
  }

protected:
  // Rebuild MaxSAT solver
  //
  Solver *rebuildSolver(); // Rebuild MaxSAT solver.

  // Other
  void initAssumptions(); // Relaxes soft clauses.

  StatusCode unweighted();
  StatusCode weighted();

  Solver *solver;  // SAT Solver used as a black box.

  // Controls the incremental strategy used by OLL algorithms.
  int incremental_strategy;
  // Controls the cardinality encoding used by OLL algorithms.
  int encoding;

  std::map<Lit, int> coreMapping; // Mapping between the assumption literal and
                                  // the respective soft clause.

  // lit -> <ID, bound, weight>
  std::map<Lit, std::pair<std::pair<int, uint64_t>, uint64_t>> boundMapping;

  // Soft clauses that are currently in the MaxSAT formula.

  uint64_t findNextWeightDiversity(uint64_t weight);
  uint64_t findNextWeight(uint64_t weight);

  bool literalTrueInModel(Lit l, const vec<lbool> & curModel);
  uint64_t min_weight;

  void logPrint(std::string s, int verb = _VERBOSITY_SOME_);

  //////////////////
  int setAssumptions(vec<Lit> & assumptions_out);

  //RETURNS UNSAT IF NO SOLUTIONS
  StatusCode setup();
  MaxSATFormula * standardizeMaxSATFormula();

  time_t time_start;
	time_t time_best_solution;
  time_t time_best_lb;

  std::string print_timeSinceStart();
  time_t timeSinceStart();
  bool checkModel();
  int num_hardened;
  int num_card_dropped; 
  int nRealSoft();
  uint64_t computeCostCore(const vec<Lit> &conflict);
  void decreaseWeights(const vec<Lit> &core, const uint64_t min_core,  vec<Encoder *> & soft_cardinality);
  void reformulateCore(const vec<Lit> &core, const uint64_t min_core,  vec<Encoder *> & soft_cardinality);
  void increaseBound(Lit cardAssump, const uint64_t min_core, vec<Encoder *> & soft_cardinality);
  void processCore(const vec<Lit> & orig_core, vec<Lit> &processed_core );
  vec<uint64_t> origWeights;
  uint64_t computeCostFromLabels(vec<lbool> &currentModel);
  void printProgress();

  ///WCE
  int nonreformulatedCores; 
  vec< vec<Lit> > cores_found;
  vec < uint64_t > core_weights; 

  bool reformulateDelayed(vec<Encoder *> & soft_cardinality);
  
  vec<lbool> bestModel;

  //Preprocessor 
  maxPreprocessor::PreprocessorInterface * pif;
  bool do_preprocess;
  MaxSATFormula * preprocess();
  int lit2Int(Lit l);
  Lit int2Lit(int l);
  void ppClause2SolClause(vec<Lit>  & solClause_out, const std::vector<int> & ppClause);
  void solClause2ppClause(const vec<Lit>  & solClause,  std::vector<int> & ppClause_out);
  uint64_t weightRemoved;

  /// HARDENING
  uint64_t maxw_nothardened;
  uint64_t ubLabelCost;
  vec<lbool>  hardeningModel; 
  
  int num_hardened_me;
  Solver* hardenClauses(vec<Encoder *> & soft_cardinality);

  void findCardinality(Lit p, int64_t & cur_bound_out, uint64_t & bound_w_out, Encoder * & bound_enc_out, int & soft_card_id, vec<Encoder *> & soft_cardinality); 
  void extendCardEnc(Lit p, vec<Encoder *> & soft_cardinality); 

  void resetActivities(); 

  MaxSATFormula* cost_formula;
  void extendModel();
  uint64_t computeCostFromClauses(vec<lbool> &currentModel);
  void reconstruct(vec<lbool> &currentModel, vec<lbool> &reconstructed_out);
  bool use_reconstruct;

  int nOrigVars; 
  


  struct LitEqual {
    public:
      bool operator()(const Lit & l1, const Lit & l2) const {
      int val1 = sign(l1) ? -(var(l1) + 1) : var(l1) + 1;
      int val2 = sign(l2) ? -(var(l2) + 1) : var(l2) + 1;
      return val1 == val2;
    }
  };

  struct LitHash {
    public:
	    size_t operator()(const Lit & l) const {
		      int size = sign(l) ? -(var(l) + 1) : var(l) + 1;
		      return std::hash<int>()(size);
	    }
};
std::unordered_set<Lit, LitHash, LitEqual> notHardened_ind;
std::unordered_set<Lit, LitHash, LitEqual> cardinality_assumptions;




};
} // namespace openwbo

#endif
