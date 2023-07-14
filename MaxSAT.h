/*!
 * \author Ruben Martins - ruben@sat.inesc-id.pt
 *
 * @section LICENSE
 *
 * MiniSat,  Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
 *           Copyright (c) 2007-2010, Niklas Sorensson
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

#ifndef MaxSAT_h
#define MaxSAT_h

#ifdef SIMP
#include "simp/SimpSolver.h"
#else
#include "core/Solver.h"
#endif

#include "MaxSATFormula.h"
#include "MaxTypes.h"
#include "utils/System.h"
#include "maxpre2/src/preprocessorinterface.hpp"
#include <algorithm>
#include <map>
#include <set>
#include <utility>
#include <iostream>
#include <vector>

using NSPACE::vec;
using NSPACE::Lit;
using NSPACE::lit_Undef;
using NSPACE::mkLit;
using NSPACE::lbool;
using NSPACE::Solver;
using NSPACE::cpuTime;

namespace openwbo {

class MaxSAT {

public:

  MaxSAT(MaxSATFormula *mx) {
    maxsat_formula = mx;
    original_labels = NULL; 
    full_original_scla = NULL;

    searchStatus = _UNKNOWN_;

    // 'ubCost' will be set to the sum of the weights of soft clauses
    //  during the parsing of the MaxSAT formula.
    ubCost = 0;
    lbCost = 0;

    off_set = 0;

    // Statistics
    nbSymmetryClauses = 0;
    nbCores = 0;
    nbSatisfiable = 0;
    sumSizeCores = 0;

    print_model = false;
    print_soft = false;
    print = false;
    unsat_soft_file = NULL;

    pif = NULL;
    prepro_techs = "";
    preprocess_time_limit = -1;
    prepro_verb = 0;
    do_preprocess = false;
    ub_prepro = UINT64_MAX;
    gate_extraction = false; 
    label_matching = true;
    skip_technique = 20; 

  }

  MaxSAT() {
    maxsat_formula = NULL;
    original_labels = NULL; 
    full_original_scla = NULL;

    searchStatus = _UNKNOWN_;

    // 'ubCost' will be set to the sum of the weights of soft clauses
    //  during the parsing of the MaxSAT formula.
    ubCost = 0;
    lbCost = 0;

    off_set = 0;

    // Statistics
    nbSymmetryClauses = 0;
    nbCores = 0;
    nbSatisfiable = 0;
    sumSizeCores = 0;

    print_model = false;
    print_soft = false;
    print = false;
    unsat_soft_file = NULL;
    
    pif = NULL;
    prepro_techs = "";
    preprocess_time_limit = -1;
    prepro_verb = 0;
    do_preprocess = false;
    ub_prepro = UINT64_MAX;
    
    gate_extraction = false; 
    label_matching = true;
    skip_technique = 20; 
  }

  virtual ~MaxSAT() {
    if (maxsat_formula != NULL)
      delete maxsat_formula;
    if (original_labels != NULL)
      delete original_labels;
    if (full_original_scla != NULL)
      delete full_original_scla;
  }

  void setInitialTime(double initial); // Set initial time.

  // Print configuration of the MaxSAT solver.
  // virtual void printConfiguration();
  void printConfiguration();

  // Encoding information.
  void print_AMO_configuration(int encoding);
  void print_PB_configuration(int encoding);
  void print_Card_configuration(int encoding);

  // Incremental information.
  void print_Incremental_configuration(int incremental);

  virtual StatusCode search();      // MaxSAT search.
  virtual void printAnswer(int type); // Print the answer.

  bool hasSolution() {return model.size() > 0;};

  // Tests if a MaxSAT formula has a lexicographical optimization criterion.
  bool isBMO(bool cache = true);

  void loadFormula(MaxSATFormula *maxsat) {
    maxsat_formula = maxsat;
    maxsat_formula->setInitialVars(maxsat_formula->nVars());

    if (maxsat_formula->getObjFunction() != NULL) {
      off_set = maxsat_formula->getObjFunction()->_const;
      maxsat_formula->convertPBtoMaxSAT();
    }

    ubCost = maxsat_formula->getSumWeights();
  }

  void blockModel(Solver *solver);

  Soft &getSoftClause(int i) { return maxsat_formula->getSoftClause(i); }
  Hard &getHardClause(int i) { return maxsat_formula->getHardClause(i); }
  Lit getAssumptionLit(int soft) {
    return maxsat_formula->getSoftClause(soft).assumption_var;
  }
  Lit getRelaxationLit(int soft, int i = 0) {
    return maxsat_formula->getSoftClause(soft).relaxation_vars[i];
  }

  int64_t getOffSet() { return off_set; }

  MaxSATFormula *getMaxSATFormula() { return maxsat_formula; }

  void setPrintModel(bool model) { print_model = model; }
  bool getPrintModel() { return print_model; }

  void setPrint(bool doPrint) { print = doPrint; }
  bool getPrint() { return print; }

  void setPrintSoft(const char* file) { 
    if (file != NULL){
      unsat_soft_file = (char*)malloc(sizeof(char) * (sizeof(file)));
      strcpy(unsat_soft_file,file) ; 
      print_soft = true;
    }
  }
  bool isPrintSoft() { return print_soft; }
  char * getPrintSoftFilename() { return unsat_soft_file; }

  /** return status of current search
   *
   *  This method helps to extract the status in case the solver is used as a
   *  library object without printing solutions.
   */
  StatusCode getStatus() { return searchStatus; }

  /** return truth values for variables
   *
   *  This method returns the truth value for a variable in the internal
   *  format. However, the return value reflects the external format, e.g.
   *  getValue(0) will return 1 or -1, depending on the sign of the variable
   *  in the model.
   */
  int getValue(const NSPACE::Var v)
  {
    if (v > model.size()) return 0;
    if (model[v] == l_True) return v+1;
    return -(int)v - 1;
  }

  /** Preprocessing methods
   * 
   * 
   * 
  */
  void set_preprocessing_parameters
        (double preprocess_time_limit = 30, std::string pre_techs = "[u]#[uvsrgVGc]",
        bool gate_extraction = false, bool label_matching = true, int skip_technique = 20); 
  void setup_formula();
    //replaces the internal maxsat_formula with the preprocessed one 
 void print_statistics();


 protected:
  //Preprocessing relatied 
  MaxSATFormula* preprocessed_formula();
  MaxSATFormula* standardized_formula();


  void reconstruct_model_prepro(vec<lbool> &currentModel, vec<lbool> &reconstructed_out);
  
  maxPreprocessor::PreprocessorInterface * pif;
  std::string prepro_techs;
  double preprocess_time_limit;
  int prepro_verb;
  bool do_preprocess;
  bool gate_extraction;
  bool label_matching;
  int skip_technique; 
  uint64_t ub_prepro;
  
  
  int lit2Int(Lit l);
  Lit int2Lit(int l);
  void ppClause2SolClause(vec<Lit>  & solClause_out, const std::vector<int> & ppClause);
  void solClause2ppClause(const vec<Lit>  & solClause,  std::vector<int> & ppClause_out);

  



  // Interface with the SAT solver
  //
  Solver *newSATSolver(); // Creates a SAT solver.
  // Solves the formula that is currently loaded in the SAT solver.
  lbool searchSATSolver(Solver *S, vec<Lit> &assumptions, bool pre = false);
  lbool searchSATSolver(Solver *S, bool pre = false);

  void newSATVariable(Solver *S); // Creates a new variable in the SAT solver.

  void reserveSATVariables(Solver *S, unsigned maxVariable); // Reserve space for multiple variables in the SAT solver.

  // Properties of the MaxSAT formula
  //
  vec<lbool> model; // Stores the best satisfying model.
  vec<lbool> model_of_original; // stores the best known solution of the non-åreårocessed isntance, obtained by calling reconstruct_model_prepro
  StatusCode searchStatus; // Stores the current state of the formula

  // Statistics
  //
  int nbCores;           // Number of cores.
  int nbSymmetryClauses; // Number of symmetry clauses.
  uint64_t sumSizeCores; // Sum of the sizes of cores.
  int nbSatisfiable;     // Number of satisfiable calls.

  // Bound values
  //
  uint64_t ubCost; // Upper bound value.
  uint64_t lbCost; // Lower bound value.
  uint64_t cost_removed_preprocessing; 
  int64_t off_set; // Offset of the objective function for PB solving.

  MaxSATFormula *maxsat_formula;
  MaxSATFormula *original_labels;
  MaxSATFormula *full_original_scla;

  // Others
  // int currentWeight;  // Initialized to the maximum weight of soft clauses.
  double initialTime; // Initial time.
  int verbosity;      // Controls the verbosity of the solver.
  bool print_model;   // Controls if the model is printed at the end.
  bool print;         // Controls if data should be printed at all
  bool print_soft;    // Controls if the unsatified soft clauses are printed at the end.
  char * unsat_soft_file;  // Name of the file where the unsatisfied soft clauses will be printed.

  // Different weights that corresponds to each function in the BMO algorithm.
  std::vector<uint64_t> orderWeights;

  // Utils for model management
  //
  bool literalTrueInModel(Lit l, vec<lbool> &model);
  void saveModel(vec<lbool> &currentModel); // Saves a Model.
  // Compute the cost of a model.
  uint64_t computeCostOriginalClauses(vec<lbool> &reconstructed_model); // checks the soft clauses of the input formula reconstructing solutions if needed
  uint64_t computeCostObjective(vec<lbool> &premodel); // computes cost from the standardized formula or preprocessed formuals objective funciton.

  // Utils for printing
  //
  void printBound(int64_t bound); // Print the current bound.
  void printModel(); // Print the best satisfying model.
  void printStats(); // Print search statistics.
  std::string printSoftClause(int id); // Prints a soft clause.
  void printUnsatisfiedSoftClauses(); // Prints unsatisfied soft clauses.
  void logPrint(std::string s); // Prints string if verbosity > 0

  // Greater than comparator.
  bool static greaterThan(uint64_t i, uint64_t j) { return (i > j); }
};
} // namespace openwbo

#endif
