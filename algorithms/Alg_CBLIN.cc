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

#include "Alg_CBLIN.h"
#include "BouMS/BouMS.h"
#include "BouMS/cores.h"
#include "BouMS/wcnf.h"

using namespace openwbo;

/************************************************************************************************
 //
 // Rebuild MaxSAT solver
 //
 ************************************************************************************************/

/*_________________________________________________________________________________________________
  |
  |  rebuildWeightSolver : (strategy : int)  ->  [Solver *]
  |
  |   Description:
  |
  |    Rebuilds a SAT solver with the current MaxSAT formula using a
  |    weight-based strategy.
  |    Only soft clauses with weight greater or equal to 'currentWeight' are
  |    considered in the working MaxSAT formula.
  |
  |   For further details see:
  |     * Ruben Martins, Vasco Manquinho, Inês Lynce: On Partitioning for
  |       Maximum Satisfiability. ECAI 2012: 913-914
  |	    * Carlos Ansótegui, Maria Luisa Bonet, Joel Gabàs, Jordi Levy:
  |       Improving SAT-Based Weighted MaxSAT Solvers. CP 2012: 86-101
  |
  |   Pre-conditions:
  |     * Assumes that 'currentWeight' has been previously updated.
  |     * Assumes that the weight strategy is either '_WEIGHT_NORMAL_' or
  |       '_WEIGHT_DIVERSIFY_'.
  |
  |   Post-conditions:
  |    
  |
  |________________________________________________________________________________________________@*/
void CBLIN::updateSolver() {

  solverCad->reserve(maxsat_formula->nVars());
  vars_added = maxsat_formula->nVars();

  for (int i = clauses_added; i < maxsat_formula->nHard(); i++) {
    ICadical::addClause(solverCad, maxsat_formula->getHardClause(i).clause);
  }
  
  clauses_added = maxsat_formula->nHard();
  softs_added = maxsat_formula->nSoft();


  //We do not support these
  assert(maxsat_formula->nPB() == 0); 
  //We do not support these
  assert(maxsat_formula->nCard() == 0);
}

/*_________________________________________________________________________________________________
  |
  |  updateCurrentWeight : (strategy : int)  ->  [void]
  |
  |  Description:
  |
  |    Updates the value of 'currentWeight' with a predefined strategy.
  |
  |  Pre-conditions:
  |    * Assumes that the weight strategy is either '_WEIGHT_NORMAL_' or
  |      '_WEIGHT_DIVERSIFY_'.
  |
  |  Post-conditions:
  |    * 'currentWeight' is updated by this method.
  |
  |________________________________________________________________________________________________@*/
void CBLIN::updateCurrentWeight(int strategy) {

  assert(strategy == _WEIGHT_NORMAL_ || strategy == _WEIGHT_DIVERSIFY_);

  if (strategy == _WEIGHT_NORMAL_)
    maxsat_formula->setMaximumWeight(findNextWeight(maxsat_formula->getMaximumWeight()));
  else if (strategy == _WEIGHT_DIVERSIFY_) {
    maxsat_formula->setMaximumWeight(findNextWeightDiversity(maxsat_formula->getMaximumWeight()));
  }
  
  logPrint("next_strat_weight=" ,maxsat_formula->getMaximumWeight(), " time " , timeSinceStart());
}


/*_________________________________________________________________________________________________
  |
  |  findNextWeight : (weight : uint64_t)  ->  [uint64_t]
  |
  |  Description:
  |
  |    Finds the greatest weight that is smaller than the 'weight'.
  |
  |  For further details see:
  |    * Ruben Martins, Vasco Manquinho, Inês Lynce: On Partitioning for Maximum
  |      Satisfiability. ECAI 2012: 913-914
  |
  |________________________________________________________________________________________________@*/
uint64_t CBLIN::findNextWeight(uint64_t weight) {

  uint64_t nextWeight = 1;
  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    if (maxsat_formula->getSoftClause(i).weight > nextWeight &&
        maxsat_formula->getSoftClause(i).weight < weight)
      nextWeight = maxsat_formula->getSoftClause(i).weight;
  }

  return nextWeight;
}

/*_________________________________________________________________________________________________
  |
  |  findNextWeightDiversity : (weight : uint64_t)  ->  [uint64_t]
  |
  |  Description:
  |
  |  Finds the greatest weight that is smaller than the 'currentWeight' and
  |  respects a given ratio
  |  between the number of different weights and the number of soft clauses.
  |
  |  Pre-conditions:
  |    * Assumes that the weight strategy is '_WEIGHT_DIVERSIFY_'.
  |    * Assumes that 'unsatSearch' was call before (this implies that
  |      nbSatisfable > 0).
  |
  |  For further details see:
  |    * Carlos Ansótegui, Maria Luisa Bonet, Joel Gabàs, Jordi Levy: Improving
  |      SAT-Based Weighted MaxSAT Solvers. CP 2012: 86-101
  |
  |________________________________________________________________________________________________@*/
uint64_t CBLIN::findNextWeightDiversity(uint64_t weight) {

  assert(weightStrategy == _WEIGHT_DIVERSIFY_);
  assert(nbSatisfiable > 0); // Assumes that unsatSearch was done before.

  uint64_t nextWeight = weight;
  int nbClauses = 0;
  std::set<uint64_t> nbWeights;
  float alpha = 1.25;

  bool findNext = false;

  for (;;) {
    if (nbSatisfiable > 1 || findNext)
      nextWeight = findNextWeight(nextWeight);

    nbClauses = 0;
    nbWeights.clear();
    for (int i = 0; i < maxsat_formula->nSoft(); i++) {
      if (maxsat_formula->getSoftClause(i).weight >= nextWeight) {
        nbClauses++;
        nbWeights.insert(maxsat_formula->getSoftClause(i).weight);
      }
    }

    if ( ( (float)nbClauses / nbWeights.size() > alpha && nbClauses > nbCurrentSoft) ||
        nbClauses == nRealSoft())
      break;

    if (nbSatisfiable == 1 && !findNext)
      findNext = true;
  }

  return nextWeight;
}
/*_________________________________________________________________________________________________
  |
  |  hardenclauses : void
  |
  |  Description:
  |
  |  Hardens soft clauses that are heavier than the current known gap between the upper and lower bound
  |   
  |  Pre-conditions:
  |    
  |
  |  For further details see:
  |     Antonio Morgado, Federico Heras, and Joao Marques-Silva. 2012. Improvements to core-guided binary search for MaxSAT. 
  |     In Proceedings of the 15th international conference on Theory and Applications of Satisfiability Testing (SAT'12), 
  |     Alessandro Cimatti and Roberto Sebastiani (Eds.). Springer-Verlag, Berlin, Heidelberg, 284-297. 
  |
  |
  |________________________________________________________________________________________________@*/

  void CBLIN::hardenClauses() { 
	   uint64_t gap = ubCost - lbCost;
     logPrint("hardening at gap: " , gap);
     int num_hardened_round = 0;
	   maxw_nothardened = 0;
     vec<Lit> toAdd;
	   for (int i = 0; i < softs_added; i++)
		  {
			bool satisfied = false;
      Lit l =  maxsat_formula->getSoftClause(i).clause[0];
      satisfied = literalTrueInModel(l, bestModel);  
			if (maxsat_formula->getSoftClause(i).weight > gap || (maxsat_formula->getSoftClause(i).weight == gap && satisfied) ) {  // 
				toAdd.push(l);
				maxsat_formula->getSoftClause(i).weight = 0;
        maxsat_formula->getSoftClause(i).assumption_var = lit_Undef;
				num_hardened++;
				num_hardened_round++;
			}
		else if (maxsat_formula->getSoftClause(i).weight > maxw_nothardened) {
				maxw_nothardened = maxsat_formula->getSoftClause(i).weight;
		} 	
		}

    for (int i = 0; i < toAdd.size(); i++) {
      Lit l = toAdd[i];
      vec<Lit> clause;
      clause.push(l);
      ICadical::addClause(solverCad, clause);
      maxsat_formula->addHardClause(clause); 
    }
		logPrint("Hardened in total: " + std::to_string(num_hardened_round) + " clauses");
    logPrint("Hardening again at gap " + std::to_string(maxw_nothardened));
   }



   
void CBLIN::flipValueinBest(Lit l) {
  assert(var(l) < bestModel.size());
  assert(bestModel[var(l)] != l_Undef);
  if (bestModel[var(l)] == l_False) {
    bestModel[var(l)] = l_True;
  }
  else {
    bestModel[var(l)] = l_False;
  }
}

  void CBLIN::hardenClausesSIS(uint64_t reduced_cost) { 
      if (incremental_DPW) {
        logPrint("can not harden based on reduced cost with incremental DPW");
        return; 
      }
      if (solverCad->status() != 10) {
        logPrint("Solver not in satisfiable state");
        return;
      }
     int num_hardened_round = 0;
	   max_coeff_nothardened_sis = 0;
     uint64_t precision = maxsat_formula->getMaximumWeight();

     logPrint("hardening in SIS, reduced cost ", reduced_cost, " precision ", precision);
     vec<Lit> toHarden; 
     toHarden.clear();
	   for (int i = 0; i < maxsat_formula->nSoft(); i++)
		  {
      if (maxsat_formula->getSoftClause(i).weight <= precision) {
        continue;
      }
      Lit l =  maxsat_formula->getSoftClause(i).clause[0];
      uint64_t red_weight = maxsat_formula->getSoftClause(i).weight / precision;
      assert(l != lit_Undef);
			if (red_weight > reduced_cost || (red_weight == reduced_cost && literal_sat_in_cadical(l)) ) {  //
        toHarden.push(l);
				num_hardened_round++;
			}
			else if (red_weight > maxw_nothardened) {
				max_coeff_nothardened_sis = red_weight;
			} 
		}
    for (int i = 0; i < toHarden.size(); i++) {
        Lit l = toHarden[i];
        vec<Lit> clause;
				clause.clear();
				clause.push(l);
				ICadical::addClause(solverCad, clause);
    }
		logPrint("hardened in total: " , num_hardened_round, " literals");
    logPrint("Hardening again at red-cost ", max_coeff_nothardened_sis);
   }
/************************************************************************************************
 //
 // Varying resolution 
 //
 ************************************************************************************************/  
  uint64_t CBLIN::precision_factors() {
    return incremental_DPW ? 2 : non_inc_precision;
  }
  
  uint64_t CBLIN::get_Maximum_Weight() {
    uint64_t maxW = 1;
    for (int i = 0; i < maxsat_formula->nSoft(); i++) {
      if (maxsat_formula->getSoftClause(i).weight > maxW) {
        maxW = maxsat_formula->getSoftClause(i).weight;
      }
    }
    return maxW;
  }

  uint64_t CBLIN::compute_first_precision() {
    uint64_t weightCand;
    if (incremental_DPW) {
        have_encoded_precision = true; 
        weightCand = dpw_next_precision();
    }
    else {      
        uint64_t varresFactor = precision_factors();
        uint64_t maxW =  get_Maximum_Weight();
        int counter = 0; //64 bits, isokay  
        while (maxW > 0) {
          counter++;
          maxW = maxW / varresFactor; 
      }
      weightCand = pow(varresFactor, counter - 1);   
    }
    max_weight_after_cg = weightCand;
    return weightCand;
  }

  void CBLIN::init_SIS_precision() {
    maxsat_formula->setMaximumWeight(compute_first_precision());
    logPrint("first precision for SIS " , maxsat_formula->getMaximumWeight());
  }

  uint64_t CBLIN::compute_next_SIS_precision(uint64_t current_precision) {
    uint64_t precision_factor = precision_factors();
    uint64_t nextFactor;

    if (incremental_DPW && have_encoded_precision) {
      nextFactor = dpw_next_precision();
    }
    else {
     nextFactor =  current_precision / precision_factor;
    }
     
    while (moreThanWeight(nextFactor) == nbCurrentSoft && nextFactor > 1 ) {
      nextFactor /= precision_factor; 
    }
    return nextFactor;
  }

  void CBLIN::update_SIS_precision() { 
    maxsat_formula->setMaximumWeight(compute_next_SIS_precision(maxsat_formula->getMaximumWeight()));
    logPrint("new precision for SIS " , maxsat_formula->getMaximumWeight());
  }

  void CBLIN::set_up_objective_counter(uint64_t init) {
      logPrint("building objective counters");

      int maxExponent = exponent(init);
      for (int i = 0; i <= maxExponent; i++) coeff_counter.push_back(0);

      for (int i = 0; i < maxsat_formula->nSoft(); i++) {
        uint64_t w = maxsat_formula->getSoftClause(i).weight;
        for (int i = 0; i <= maxExponent; i++) {
          if (w >= raise_to(i) ) {
            coeff_counter[i]++;
          }
          else {
            break;
          }
        }  
      }

      logPrint("building done");
      weight_map_setup = true;
  }

   uint64_t CBLIN::raise_to(int exponent) {
      uint64_t precision_factor = precision_factors(); 
      if (exponent == 0) {
        return 1;
      }
      if (exponent == 1) {
        return precision_factor;
      }
      else if (exponent % 2 == 0) {
        return raise_to(exponent / 2) * raise_to(exponent / 2);
      }
      else {
        return precision_factor * raise_to(exponent - 1);
      }
   }

  int CBLIN::exponent(uint64_t weight) {
    assert(weight > 0);
    uint64_t precision_factor = precision_factors(); 
    if (weight == 1) {
      return 0; 
    }
    int exponent = 0; 
    if (precision_factor == 2) {
      while (weight >>= 1) ++exponent;
    }
    else {
      while (weight) {
        exponent++;
        weight /= precision_factor;
      }
    }
    return exponent;
  }

  int CBLIN::moreThanWeight(uint64_t weightCand) {
    if (!weight_map_setup) {
          set_up_objective_counter(max_weight_after_cg);
    }
    
    return coeff_counter[exponent(weightCand)];
  }


/************************************************************************************************
 //
 // Utils for core management
 //
 ************************************************************************************************/

/*_________________________________________________________________________________________________
  |
  |  PMRES : (lits : vec<Lit>&)  ->  [void]
  |
  |  Description:
  |
  |    The PMRES transformation
  |
  |  For further details see:
  |    * Nina Narodytska and Fahiem Bacchus. 2014. Maximum satisfiability using core-guided MAXSAT resolution. 
  |      In Proceedings of the Twenty-Eighth AAAI Conference on Artificial Intelligence (AAAI'14). AAAI Press 2717-2723.
  |
  |  Pre-conditions:
  |    * Assumes that 'core' is not empty.
  |
  |  Post-conditions:
  |    * 'hardClauses' are updated with the clauses that encode the PMRES
  |      constraint. Soft clauses are immediately added as hard and onyl the assumption is used in soft
  |
  |________________________________________________________________________________________________@*/
void CBLIN::encodeMaxRes(vec<Lit> &core, uint64_t weightCore)
{
	assert(core.size() != 0); 

  int n = core.size(); 
  vec<Lit> dVars;
	vec<Lit> clause; 

  for (int i = 0; i < n - 1; i++)
    {
      Lit p = maxsat_formula->newLiteral();
      dVars.push(p);
    }
  
  // NEW HARD CLAUSES
  // lins == 0 -> only run PMRES 
  if (lins == 0) {
    clause.clear();
    core.copyTo(clause);
    maxsat_formula->addHardClause(clause);
  }

  if (n > 2) {
		for (int i = 0; i < n-2; i++) {
			// d_i -> (b_{i+1} v d_{i+1})
			// clause = { ~dVars[i], dVars[i + 1], core[i + 1] };
      // Not needed for completeness 
			
      if (lins == 0) {
      	clause.clear();
				clause.push(~dVars[i]);
				clause.push(dVars[i + 1]);
				clause.push(core[i + 1]);
				maxsat_formula->addHardClause(clause);
      }
		
		
			// (b_{i+1} v d_{i+1}) -> d_i
      
			// d_i v -b_{i+1}
			// clause = { dVars[i], ~core[i + 1] };
			clause.clear();
			clause.push(dVars[i]);
			clause.push(~core[i + 1]);
			maxsat_formula->addHardClause(clause);
			
			 // d_i v -d_{i+1}
			 // clause = { dVars[i], ~dVars[i + 1] };
			clause.clear();
			clause.push(dVars[i]);
			clause.push(~dVars[i + 1]);
			maxsat_formula->addHardClause(clause);	
		}
	}
    
    if (n > 1) {
		 // handle i = p - 1 case
		 // clause = { dVars[p - 2], ~core[p - 1] };
		 clause.clear();
		 clause.push(dVars[n - 2]);
		 clause.push(~core[n - 1]);
		 maxsat_formula->addHardClause(clause);
		 
		 // clause = { ~dVars[p - 2], core[p - 1] };
		 clause.clear();
		 clause.push(~dVars[n - 2]);
		 clause.push(core[n - 1]);
		 maxsat_formula->addHardClause(clause);
	}
	
	// NEW SOFT CLAUSES
    for (int i = 0; i < n-1; i++) {
		//clause = { ~b_i, ~d_i };
		clause.clear();
		clause.push(~core[i]);
		clause.push(~dVars[i]);
		addSoftClauseAndAssumptionVar(weightCore, clause);
	}
}

/*_________________________________________________________________________________________________
  |
  |  relaxCore : (conflict : vec<Lit>&) (weightCore : int) 
  |              ->  [void]
  |
  |  Description:
  |
  |    Relaxes the core as described in the original WBO paper.
  |
  |  For further details see:
  |    * Vasco Manquinho, Joao Marques-Silva, Jordi Planes: Algorithms for
  |      Weighted Boolean Optimization. SAT 2009: 495-508
  |
  |  Pre-conditions:
  |    * Assumes that the core ('conflict') is not empty.
  |    * Assumes that the weight of the core is not 0 (should always be greater
  |      than or equal to 1).
  |
  |  Post-conditions:
  | 
  |    * If the weight of the soft clause is not the same as the weight of the
  |      core:
  |      - 'softClauses[indexSoft].weight' is decreased by the weight of the
  |        core.
  |    * 'sumSizeCores' is updated.
  |
  |________________________________________________________________________________________________@*/
void CBLIN::relaxCore(vec<Lit> &core, uint64_t weightCore) {

  assert(core.size() > 0);
  assert(weightCore > 0);


  for (int i = 0; i < core.size(); i++) {
    int indexSoft = coreMapping[core[i]];
    assert(maxsat_formula->getSoftClause(indexSoft).weight >= weightCore);
    maxsat_formula->getSoftClause(indexSoft).weight -= weightCore;

    if(maxsat_formula->getSoftClause(indexSoft).weight == 0) {
      maxsat_formula->getSoftClause(indexSoft).assumption_var = lit_Undef;
      num_hardened++;
    }
  }
  encodeMaxRes(core, weightCore);
  sumSizeCores += core.size();

  if (ls_cores) {
    BouMS_cores_core_t c;
    c.numLiterals = core.size();
    c.literals = new BouMS_literal_t[c.numLiterals];
    if (!c.literals) {
      logPrint("Failed to allocate memory for a BouMS core!");
    } else {
      for (unsigned int lIdx = 0; lIdx < c.numLiterals; ++lIdx) {
        const auto lit = core[lIdx];
        c.literals[lIdx] = BouMS_mkLit(var(lit), sign(lit));
      }
      cores.push_back(c);
      logPrint("Added core of length ", c.numLiterals, " to BouMS's core list");
    }
  }
}

/*_________________________________________________________________________________________________
  |
  |  computeCostCore : (conflict : vec<Lit>&)  ->  [int]
  |
  |    Description:
  |
  |      Computes the cost of the core. The cost of a core is the minimum coefficient
  |      of the objective literals that appear in that core.
  |
  |    Pre-conditions:
  |      * Assumes that 'conflict' is not empty.
  |
  |________________________________________________________________________________________________@*/
uint64_t CBLIN::computeCostCore(const vec<Lit> &core) {

  assert(core.size() != 0);

  if (maxsat_formula->getProblemType() == _UNWEIGHTED_) {
    return 1;
  }

  uint64_t coreCost = UINT64_MAX;
  for (int i = 0; i < core.size(); i++) {
    int indexSoft = coreMapping[core[i]];
    if (maxsat_formula->getSoftClause(indexSoft).weight < coreCost)
      coreCost = maxsat_formula->getSoftClause(indexSoft).weight;
  }

  return coreCost;
}

/************************************************************************************************
 //
 // SEARCHES
 //
 ************************************************************************************************/

/*_________________________________________________________________________________________________
  |
  |  unsatSearch : [void] ->  [void]
  |
  |  Description:
  |
  |    Calls a SAT solver only on the hard clauses of the MaxSAT formula.
  |    If the hard clauses are unsatisfiable then the MaxSAT solver terminates
  |    and returns 'UNSATISFIABLE'.
  |    Otherwise, a model has been found and it is stored. Without this call,
  |    the termination of the MaxSAT solver is not guaranteed.
  |
  |  For further details see:
  |    * Carlos Ansótegui, Maria Luisa Bonet, Jordi Levy: SAT-based MaxSAT
  |      algorithms. Artif. Intell. 196: 77-105 (2013)
  |
  |  Post-conditions:
  |   * If the hard clauses are satisfiable then 'ubCost' is updated to the cost
  |     of the model.
  |   * If the working formula is satisfiable, then 'nbSatisfiable' is increased
  |     by 1. Otherwise, 'nbCores' is increased by 1.
  |
  |________________________________________________________________________________________________@*/
StatusCode CBLIN::unsatSearch() {

  assert(assumptions.size() == 0);

  updateSolver();
  freezeObjective();
  softsSatisfied();
  //logPrint("In Unsat in solver " + std::to_string(solver->nVars()) + " vars amd " + std::to_string(solver->nClauses()) + " clauses" );
  lbool res = ICadical::searchSATSolver(solverCad, assumptions);
  has_flipped = false;

  clearFixingsonSoft();
  

  if (res == l_False) {
    nbCores++;
    printAnswer(_UNSATISFIABLE_);
    return _UNSATISFIABLE_;
  } else if (res == l_True) {
    //flipLiterals();
    nbSatisfiable++;
    uint64_t beforecheck = ubCost;
    checkModel(false, true);
    
    uint64_t aftercheck = ubCost;
    assert(beforecheck >= aftercheck);    
  }

  return _SATISFIABLE_;
}
/*_________________________________________________________________________________________________
  |
  |  weightDisjointCores : [void] ->  [void]
  |
  |  Description:
  |
  |    Runs the sat solver repeadetly on the current settings, extracting and relaxing cores
  |    does not rebuild the SAT-solver, alter strat weight or harden clauses. 
  |
  |  For further details see:
  |    * Berg, J., & Järvisalo, M. (2017). Weight-Aware Core Extraction in SAT-Based MaxSAT Solving. CP.
  |
  |  Pre-conditions:
  |    * Assumes the setup method has been called
  |
  |  Post-conditions:
  |     * The hard clauses in formula reflect the found and relaxed cores. 
  |     * LB is updated.
  |________________________________________________________________________________________________@*/

  StatusCode CBLIN::weightDisjointCores() {
    
    for (;;) {
      if(timer->terminate()) {
        return _UNKNOWN_;
      }
      setAssumptions(assumptions);
     
     lbool res = ICadical::searchSATSolver(solverCad, assumptions);
     has_flipped = false;

      if (res == l_Undef) {
        //Interrupted
        solverCad->disconnect_terminator();
        return _UNKNOWN_;
      }

      if (res == l_False) {
      
        nbCores++;
        vec<Lit> cad_core;
        ICadical::getCore(solverCad, assumptions, cad_core);
        assert(cad_core.size() > 0);

        //logPrint("Core " + core_2_string(cad_core));

        uint64_t coreCost = computeCostCore(cad_core);
        lbCost += coreCost;
        checkGap();
        logPrint("LB ", lbCost, " core size ", cad_core.size(), " core-min-cost " , coreCost); 
        relaxCore(cad_core, coreCost);
      }

      if (res == l_True) {
        return _SATISFIABLE_; 
      }
      if (lbCost > ubCost) {
        logPrint("LB bigger than UB, something fishy is going on....");
        return _ERROR_;
      }

    }
  }

/*_________________________________________________________________________________________________
  |
  |  setup : [void] ->  [void]
  |
  |  Description:
  |
  |    Makes SAT solver and checks that solutions exist. Most other search methods assume this has been run-  
  |    RETURNS unsat if no solutions exist
  |
  |  Post-conditions:
  |     * SAT solver exists 
  |     * Solutions exists
  |________________________________________________________________________________________________@*/
  StatusCode CBLIN::setup() {

      if (maxsat_formula->nHard() == 0) {
        if (!do_preprocess && maxsat_formula->nInitialVars() > 0) {
          vec<lbool> currentModel;
          for (int i = 0; i < maxsat_formula->nInitialVars(); i++ ) currentModel.push(l_False);
          for (int i = 0; i < maxsat_formula->nSoft(); i++) {
            assert( maxsat_formula->getSoftClause(i).clause.size() == 1);
            Lit l = maxsat_formula->getSoftClause(i).clause[0];
            if (!literalTrueInModel(l, currentModel)) {
              currentModel[var(l)] = l_True;
            }
          }
          saveModel(currentModel);
        }
        return _OPTIMUM_;
      }

      while (isSoft.size() < maxsat_formula->nVars()) isSoft.push(false);
      maxw_nothardened = 1; 
      for (int i = 0; i < maxsat_formula->nSoft(); i++)  {
          assert( maxsat_formula->getSoftClause(i).clause.size() == 1);
          Lit l = maxsat_formula->getSoftClause(i).clause[0];
          if (maxsat_formula->getSoftClause(i).weight > maxw_nothardened) {
            maxw_nothardened = maxsat_formula->getSoftClause(i).weight;
          }
          assert(var(l) < isSoft.size());
          isSoft[var(l)] = true; 
          if ( maxsat_formula->getSoftClause(i).weight > maxw_nothardened) {
            maxw_nothardened = maxsat_formula->getSoftClause(i).weight;
          }
      }
      
      initAssumptions();  
      
      solverCad = ICadical::newSATSolver();

      if (dynamic_cadical) {
        ICadical::sat_mode(solverCad);
      }

      StatusCode rs = unsatSearch();
      if (rs == _UNSATISFIABLE_) return rs;
      
      //Here we know that the formula is SAT
      if (maxsat_formula->nSoft() == 0 || ubCost == lbCost) {
          return _OPTIMUM_; //Solved by preprocessing
      }    
      
      if (dynamic_cadical) {
       resetSolver(0); // 0 default 1 unsat 2 sat
      }

      updateCurrentWeight(weightStrategy);
      
      return rs;
  }




/*_________________________________________________________________________________________________
  |
  |  weightSearch : [void] ->  [void]
  |
  |  Description:
  |
  |    MaxSAT weight-based search. Considers the weights of soft clauses to find
  |    cores with larger weights first.
  |
  |  For further details see:
  |    * Ruben Martins, Vasco Manquinho, Inês Lynce: On Partitioning for Maximum
  |      Satisfiability. ECAI 2012: 913-914
  |    * Carlos Ansótegui, Maria Luisa Bonet, Joel Gabàs, Jordi Levy: Improving
  |      SAT-Based Weighted MaxSAT Solvers. CP 2012: 86-101
  |
  |  Pre-conditions:
  |    * Assumes 'weightStrategy' to be '_WEIGHT_NORMAL_' or
  |      '_WEIGHT_DIVERSIFY_'.
  |
  |  Post-conditions:
  |    * 'lbCost' is updated.
  |    * 'ubCost' is updated.
  |    * 'nbSatisfiable' is updated.
  |    * 'nbCores' is updated.
  |________________________________________________________________________________________________@*/
StatusCode CBLIN::weightSearch() {

  assert(weightStrategy == _WEIGHT_NORMAL_ ||
         weightStrategy == _WEIGHT_DIVERSIFY_);
  inLinSearch = false;
  
  for (;;) {
    StatusCode us = weightDisjointCores(); 

    //LB phase proves optimality, current model is not for the current formula. 
    if (us == _OPTIMUM_) {
        logPrint("LB = UB");
        return getModelAfterCG();
    }

    //At this point solver returned true and as such has a model
   
    nbSatisfiable++;

    checkModel(false, false);

    if (lbCost == ubCost) {
      if (verbosity > 0)
        logPrint("LB = UB");          
        printAnswer(_OPTIMUM_);
        return _OPTIMUM_;
    }
    if (nbCurrentSoft == nRealSoft()) {
      assert(ubCost == lbCost);
      printBound(lbCost);
      printAnswer(_OPTIMUM_);
      return _OPTIMUM_;
    } 
    if (ubCost - lbCost < maxw_nothardened) {
      hardenClauses();
    }
    if (shouldUpdate()) {
      updateSolver();
    } 
    else {
      updateCurrentWeight(weightStrategy);
    }
    
  }
}



/*_________________________________________________________________________________________________
  |
  |  weightPMRES+Linear: [void] ->  [void]
  |
  |  Description:
  |
  |    New idea, run weight aware disjoint phase and then do the rest by Linear search 
  |
  | 
  |________________________________________________________________________________________________@*/
StatusCode CBLIN::coreGuidedLinearSearch() {
  inLinSearch = false;
  timer->start_timer();
  solverCad->connect_terminator(timer);
  if (ls_learn_clauses) {
    solverCad->connect_learner(&small_clause_learner);
  }
  for (;;) {
    StatusCode us = weightDisjointCores(); 
    if (us == _OPTIMUM_) {
        logPrint("LB = UB");
        return getModelAfterCG();
    }

    if (us == _UNKNOWN_ ) {
        logPrint("interrupted core guided phase");
        if(shouldUpdate()) {
          logPrint("Updating solver at ", timeSinceStart());
          updateSolver();
        }
        return linearSearch();
    }
    

    //At this point solver returned true and as such has a model
    assert(us == _SATISFIABLE_ );

    logPrint("SAT during core guided phase at " , timeSinceStart());
    nbSatisfiable++;
    checkModel();

    if (lbCost == ubCost) {
      if (verbosity > 0)
        logPrint("LB = UB");
        printAnswer(_OPTIMUM_);
        return _OPTIMUM_;
    }

   
     if (nbCurrentSoft == nRealSoft()) {
      checkModel();
      if (lbCost < ubCost) {
        ubCost = lbCost;
        vec<lbool> model_out; 
        ICadical::getModel(solverCad, model_out);
        saveModel(model_out);
        printBound(lbCost);
      }
      printAnswer(_OPTIMUM_);
      return _OPTIMUM_;
    } 


   //if code gets here algorithm cant terminate yet   
   if (ubCost - lbCost < maxw_nothardened) {
        hardenClauses();
      }

   if(relaxBeforeStrat) {
      logPrint("Relax 2 Strat");
      if(shouldUpdate()) {
          logPrint("Updating solver at ", timeSinceStart());
          updateSolver();
      }
      else if (maxsat_formula->getMaximumWeight() > 1) {
              logPrint("weight update at " , timeSinceStart());
              updateCurrentWeight(weightStrategy); 
              if (maxsat_formula->getMaximumWeight() == 1) {
                logPrint("Weight = 1 -> Done with cores at ", timeSinceStart());
                return linearSearch();
              }
      }
      else {
        return _ERROR_; // Should not get here
      }
      
   }
  else {
      logPrint("Strat 2 Relax");
      if (maxsat_formula->getMaximumWeight() > 1) {
              logPrint("weight update at " , timeSinceStart());
              updateCurrentWeight(weightStrategy); 
                      
      }
      if (maxsat_formula->getMaximumWeight() == 1 && nbCores > 0) {
        if(shouldUpdate()) {
          logPrint("Updating solver at ", timeSinceStart());
          updateSolver();
        }
        return linearSearch();
        
      }
    }
  }
  //Code never gets here 
  return _ERROR_;
}

/*
  only used if the cg phase proves optimality. 
 */
StatusCode CBLIN::getModelAfterCG() {
  if (!shouldUpdate()) {
    logPrint("ERROR: CG phase proves UNSAT without finding new cores");
  }
  updateSolver();
  setAssumptions(assumptions);
  lbool res = ICadical::searchSATSolver(solverCad, assumptions); 
  has_flipped = false;

  assert(res == l_True);
  checkModel();
  assert(ubCost == lbCost);
  printAnswer(_OPTIMUM_);
  return _OPTIMUM_;
}

void CBLIN::set_observed_vars(vec<Lit> & original_objective) {
  assert(objFunction.size() == coeffs.size());
  if (verbosity > 1) {
    std::cout << "c Observing";
  }
  for (int i = 0; i < objFunction.size(); i ++) {
      int lit = lit2Int(objFunction[i]);
      if (verbosity > 1) {
      std::cout << " " << abs(lit);
      }
      solverCad->add_observed_var(abs(lit));    
  }
  for (int i = 0; i < original_objective.size(); i ++) {
      int lit = lit2Int(original_objective[i]);
      if (verbosity > 1) {
      std::cout << " " << abs(lit);
      }
      solverCad->add_observed_var(abs(lit));    
  }
  if (verbosity > 1) {
      std::cout << endl;
      }
}

uint64_t CBLIN::compute_ub_red_cost(uint64_t precision) {

  if (bestModel.size() < maxsat_formula->nVars() || solverCad->status() != 10) {
      logPrint("Extending best model to full formula");
      extendBestModel();
  }
  auto lambda = [this](Lit l){return literalTrueInModel(l, bestModel);};
  uint64_t reduced_cost = computeCostReducedWeights_prec(&lambda, precision);   
  uint64_t red_gap = known_gap / maxsat_formula->getMaximumWeight();

  if (red_gap < reduced_cost) {
      logPrint("Setting reduced_cost to reduced gap " + std::to_string(red_gap));
      reduced_cost = red_gap;
  } 

  if (ls_min) {
    //TODO maybe with reduced objective?
    localsearch(bestModel);
    auto lambda = [this](Lit l){return literalTrueInModel(l, bestModel);};
    uint64_t min_cost = computeCostReducedWeights_prec(&lambda, precision);
    if (min_cost < reduced_cost) {
      reduced_cost = min_cost;
    }
  }
   
  // if the bound is obtained from preprocessing, we can not set variables in encoding according to a model. 
  bool bound_set_by_prepro = false;
  if (do_preprocess) {
    uint64_t red_p_gap = (ub_prepro - lbCost) / precision;
    if (reduced_cost > red_p_gap) {
        logPrint("reduced cost from preprocessor gap: " ,red_p_gap, " better than best model " ,reduced_cost);
        reduced_cost = red_p_gap;
        bound_set_by_prepro = true;
    }
  }     
  setCardVars(bound_set_by_prepro);
  return reduced_cost;
}

void CBLIN::update_current_soft(uint64_t precision) {
  nbCurrentSoft = 0;
  for (int i = 0; i < objFunction.size(); i++) {
        if ((coeffs[i] / precision) > 0) {
          nbCurrentSoft++;
        }
  }
  logPrint("There are ", nbCurrentSoft, " variables in this precision");
}

StatusCode CBLIN::linearSearch_propagator() {
  logPrint( "Starting lin search with a propagator: LB: ",lbCost, " UB: " , ubCost,
            " gap: " , known_gap, " time " , timeSinceStart() );
  
  inLinSearch = true;
  assumptions.clear();

  assert(bestModel.size() > 0);
  savePhase();

  uint64_t precision = compute_first_precision();
  build_objective_func_and_coeffs_prop();

  vec<Lit> original_objective;
  vec<uint64_t> original_coeffs; 

  for (int i = 0; i < original_labels->nSoft(); i++) {
    Lit l = original_labels->getSoftClause(i).clause[0];
    uint64_t coeff = original_labels->getSoftClause(i).weight;
    original_objective.push(l);
    original_coeffs.push(coeff);
  }
  SISPropagator* prop = NULL;
  while (true) {
    logPrint("Solving precision ", precision);
    update_current_soft(precision);
    resetSolver();
    uint64_t red_cost_ub = compute_ub_red_cost(precision);
    if (red_cost_ub > 0) {
      // Call the SAT solver 
      logPrint("New propagator: ubCost ", ubCost, " red_cost_ub ", red_cost_ub, " vars ", solverCad->vars(), " obj size ", objFunction.size(), " precision ", precision);
      if (prop == NULL) {
        prop = new SISPropagator(ubCost, red_cost_ub, solverCad->vars(), objFunction, coeffs, precision, original_objective, original_coeffs,  verbosity > 1);
      }
      else {
        prop->resetPropagator(precision, red_cost_ub, solverCad->vars());
      }
      solverCad->connect_external_propagator(prop);
      set_observed_vars(original_objective);
      int res = solverCad->solve();
      solverCad->disconnect_external_propagator();
     // assert(res == 20);
//      assert(prop->best_model.size() == objFunction.size());

      //get the best model 
      resetSolver();
      freezeObjective();
      if (verbosity > 1) {
        cout << "c best model from propagator:" ;
      }
      for (int i = 0; i < (prop->best_model).size(); i++ ) {
        if (verbosity > 1) {cout << " " << (prop->best_model)[i];}
        solverCad->assume((prop->best_model)[i]);
      }
      if (verbosity > 1) { cout << endl; } 
      res = solverCad->solve();
      assert(res == 10);
      flipLiterals();
      checkModel(false, false);
    }
    if (precision > 1) {
      if (ubCost == lbCost) {
        logPrint("Terminating on bounds");
        printAnswer(_OPTIMUM_);
        return _OPTIMUM_;
      }
      else { 
        logPrint("New precision");
        precision = compute_next_SIS_precision(precision);
      }
    }
    else {
        //TODO so far incomplete, need to od more here ot get 
        logPrint("Precision 1, stopping");
        printAnswer(_SATISFIABLE_);
        return _SATISFIABLE_;
    }
  }
  return _ERROR_;
}


StatusCode CBLIN::linearSearch() {
  if (use_propagator) {
    return linearSearch_propagator();
  }

  logPrint( "Starting lin search with: LB: ",lbCost, " UB: " ,ubCost,
            " gap: " ,ubCost-lbCost, " time " , timeSinceStart() );

  inLinSearch = true;
  assumptions.clear();
  

  assert(bestModel.size() > 0);
  savePhase();

  
  if(delete_before_lin) {
    resetSolver(2); 
  }
   
  if (incremental_DPW) {
    // add all literals into the encoding, these steps do not yet encode anything. 
    assert(dpw == NULL); 
    dpw = RustSAT::dpw_new();

    if (verbosity > 1) {
            cout << "c Adding to RustSAT: ";
          }
    for (int i = 0; i < maxsat_formula->nSoft(); i++) {
      if (maxsat_formula->getSoftClause(i).weight > 0) {
          Lit l = maxsat_formula->getSoftClause(i).assumption_var; 
          assert (l != lit_Undef);
          if (verbosity > 1) {
            cout << " " << lit2Int(l) << "/" << maxsat_formula->getSoftClause(i).weight ;
          }
          RustSAT::dpw_add(dpw, lit2Int(l),  maxsat_formula->getSoftClause(i).weight);
      }
    }
    if (verbosity > 1) {
            cout <<  endl;
    }
  }

  if (ls_learn_clauses) {
    solverCad->disconnect_learner();
    learned_clauses = small_clause_learner.extract_clauses();
    if (learned_clauses.size() == 0) {
      logPrint("Did not learn any clauses for LS?!");
    }  
  }

  updateBouMSInstance();

  if (ls_cores > 1 && !boums_broken && cores.size() > 0) {
    bool deleteCores = false;

    core_var_to_clause = new BouMS_uint_t[boums_inst.numVariables];
    if (core_var_to_clause) {
      for (BouMS_uint_t coreIdx = 0; coreIdx < cores.size(); ++coreIdx) {
        const auto& core = cores[coreIdx];
        for (BouMS_uint_t litIdx = 0; litIdx < core.numLiterals; ++litIdx) {
          const auto lit = core.literals[litIdx];
          core_var_to_clause[BouMS_var(lit)] =
            maxsat_formula->nHard() + coreMapping[mkLit(BouMS_var(lit), BouMS_sign(lit))];
        }
      }

      boums_cores = new BouMS_cores_mem_t;
      if (boums_cores) {
        boums_cores->varToCores = NULL;
        boums_cores->coreToSatLits = NULL;
        if (BouMS_cores_init(&boums_inst, cores.data(), cores.size(), core_var_to_clause, boums_cores, realloc, free)) {
          deleteCores = true;
        }
      } else {
        deleteCores = true;
      }

      if (deleteCores) {
        logPrint("Could not allocate memory for BouMS cores!");

        if (boums_cores) {
          delete boums_cores;
        }

        for (auto& c : cores) {
          delete[] c.literals;
        }
        cores.clear();
      } else {
        logPrint("Added ", cores.size(), " cores to BouMS!");
      }
    } else {
      logPrint("Failed to allocate memory for BouMS' core variable to clause mapping");
    }
  }

  if (ls_replace_sis) {
    boums_params.maxTries = 1;

    localsearch(bestModel);

    if (ubCost == lbCost) {
      logPrint("LB = UB");
      printAnswer(_OPTIMUM_);
      return _OPTIMUM_;
    }

    return _SATISFIABLE_;
  }

  if (bestModel.size() < maxsat_formula->nVars() || solverCad->status() != 10) {
      logPrint("Extending best model to full formula");
      extendBestModel();
  }

  if (!boums_broken && (ls_dyn_prec || ls_sis) && ls_merge_assign) {
    bestModel.copyTo(ls_merged_assign);
    ls_usual_init_assign = &ls_merged_assign;
  }

  init_SIS_precision();
  setPBencodings();
  
  lbool res = l_True;
  bool minimize_iteration = true;
  reconstruct_iter = true;

  // int file_name_counter = 0;

  while (res == l_True) {

    logPrint("SAT Call at " , timeSinceStart(), " # assumptions " , assumptions.size(), " clauses in SAT solver " ,solverCad->irredundant());

    if (verbosity > 1) {
      cout << "c assumptions:";
      for (int i = 0; i < assumptions.size(); i++) {
        cout << " " << lit2Int(assumptions[i]);
      }
      cout << endl;
    }

    res = ICadical::searchSATSolver(solverCad, assumptions);  
    has_flipped = false;
    
    if (res == l_True) {
      nbSatisfiable++;
      
      assert(solverCad->status() == 10);
      flipLiterals();
      
      auto lambda = [this](Lit l){ return literal_sat_in_cadical(l); };
      uint64_t new_reduced_cost = computeCostReducedWeights(&lambda);
      bool better = checkModel(false, false);
      
      vec<lbool> cadModel; 
      ICadical::getModel(solverCad, cadModel);

      if (!incremental_DPW && assumptions.size() > 0) {     
        if (harden_in_SIS && !incremental_DPW && new_reduced_cost  < max_coeff_nothardened_sis) {
          logPrint("harden in sis");
          hardenClausesSIS(new_reduced_cost);
        }
        vec<Lit> clause; 
        clause.push(assumptions[0]);
        ICadical::addClause(solverCad, clause);

      }
      if (better && incremental_DPW) {
        harden_incremental();
      }

      if (ls_sis && !skip_local_search && new_reduced_cost > 0) {
        if (ls_merge_assign) {
          static const std::function<lbool(const lbool&)> lboolid = [](const lbool& b) { return b; };
          const auto num_disagree = mergeAssignments(ls_merged_assign, cadModel, lboolid);
          logPrint("Merged SIS SAT solver assignment, agreeing variables: ", boums_inst.numVariables - num_disagree,
                   ", disagreeing variables: ", num_disagree);
        }

        if (localsearch(*ls_usual_init_assign)) {
          const auto lambda = [this](Lit l) { return boums_assignment[var(l)] != sign(l); };
          const auto ls_reduced_cost = computeCostReducedWeights(&lambda);
          if (ls_reduced_cost < new_reduced_cost) {
            new_reduced_cost = ls_reduced_cost;
            logPrint("LS found better reduced cost");
            if (ls_merge_assign) {
              const auto num_disagree = mergeAssignments<bool*, bool>(ls_merged_assign, boums_assignment, tolbool);
              logPrint("Merged SIS LS assignment, agreeing variables: ", boums_inst.numVariables - num_disagree,
                       ", disagreeing variables: ", num_disagree);
            }
          }
        }
      }
      

      if (minimize_sol && new_reduced_cost > 0 && minimize_iteration && minimize_strat > 0) {
        uint64_t t = new_reduced_cost;
        minimizelinearsolution(cadModel);
        if (minimize_strat == 2) {
          minimize_iteration = false;
        }

        assert(solverCad->status() == 10);
        //cadModel.clear();
        //ICadical::getModel(solverCad, cadModel);
        auto lambda = [this](Lit l){ return literal_sat_in_cadical(l); };
        new_reduced_cost = computeCostReducedWeights(&lambda); 
        if ( t != new_reduced_cost )
          logPrint("cost minimized, before: ",  t , " after " , new_reduced_cost);
        assert(t >= new_reduced_cost);
      }

      if (reconstruct_iter && minimize_strat == 2) reconstruct_iter = false;
      

      if (ubCost == lbCost) {
        logPrint("LB = UB");
        printAnswer(_OPTIMUM_);
        return _OPTIMUM_;
      }

      if (new_reduced_cost > 0) {
        updateBoundLinSearch(new_reduced_cost - 1);
      }
      else {
        bool incremental_done = RustSAT::dpw_is_max_precision(dpw) && incremental_DPW;
        if (maxsat_formula->getMaximumWeight() == 1 || incremental_done) {
            logPrint("new reduced cost " , new_reduced_cost, " at precision 1, stopping.");
            // No need to check for fine convergence because here we have a model whose cost matches the lb proven by core-guided search
            printAnswer(_OPTIMUM_);
            return _OPTIMUM_;
        }
        else {
          logPrint("rebuilding after SAT");
          minimize_iteration = true;
          reconstruct_iter = true;
          if (!(incrementalVarres || incremental_DPW)) {  
            resetSolver(2); 
          }
          update_SIS_precision();
          setPBencodings();
        }
      }

    } 
    else { //res = false
         bool incremental_done = RustSAT::dpw_is_max_precision(dpw) && incremental_DPW;
       if (maxsat_formula->getMaximumWeight() == 1 || incremental_done) {
          if (dpw_fine_convergence_after) {
            logPrint("stopping coarse convergence");
            dpw_coarse = false;
            dpw_fine_convergence_after = false;
            updateBoundLinSearch(fine_bound); 
            res = l_True;
          }
          else {
            logPrint("UNSAT at precision 1, stopping.");
            printAnswer(_OPTIMUM_);
            return _OPTIMUM_;
          }
        }
        else {
          logPrint("Rebuilding after UNSAT");
          if (!(incrementalVarres || incremental_DPW)) {
            resetSolver(2);
          }
          minimize_iteration = true;
          reconstruct_iter = true;
          update_SIS_precision();
          setPBencodings();
          res = l_True;
        }
      
    }
  }

  return _ERROR_;
}

void CBLIN::harden_incremental() {
  uint64_t global_ub_dpw = ubCost - lbCost;
  logPrint("hardening in incremental DPW");
  RustSAT::dpw_limit_range(dpw, 0, global_ub_dpw, &dpw_clause_collector, static_cast<void*>(solverCad));
  logPrint("hardening incremental DPW bound: " , global_ub_dpw) ;
}


uint64_t CBLIN::dpw_next_precision() {
  assert(have_encoded_precision);
  uint64_t next_prec = RustSAT::dpw_next_precision(dpw);
  have_encoded_precision = false; 
  return next_prec;
}

void CBLIN::dpw_encode_and_enforce(uint64_t rhs) {
    int num_vars = solverCad->vars();
    RustSAT::dpw_encode_ub(dpw, rhs, rhs, &num_vars, &dpw_clause_collector, static_cast<void*>(solverCad));
    logPrint("rhs in encode and enforce ", rhs) ;
    assumptions.clear();
    RustSAT::MaybeError ret = RustSAT::dpw_enforce_ub(dpw, rhs, &dpw_assumps, &assumptions);
    if (ret == RustSAT::MaybeError::NotEncoded) {
      logPrint("rustsat returned not encoded");
    }
    assert(ret == RustSAT::MaybeError::Ok);
    have_encoded_precision = true;
}


void CBLIN::dpw_assumps(int lit, void *assumps) {
  ((vec<Lit> *)assumps)->push(MaxSAT::int2Lit(lit));
}

/// TDOD: with cadical this could be much more direct... 
void CBLIN::dpw_clause_collector(int lit, void *ptr) {
  CaDiCaL::Solver * solverC = static_cast<CaDiCaL::Solver *>(ptr);
  solverC->add(lit);
}


void CBLIN::updateBoundLinSearch (uint64_t newBound) {  
  logPrint("new bound to enforce: " , newBound, " at ", timeSinceStart());
  
    if (dpw_coarse) {
      uint64_t coarse_b = RustSAT::dpw_coarse_ub(dpw, newBound);
      dpw_fine_convergence_after = (coarse_b != newBound);
      fine_bound = newBound;
      newBound = coarse_b;
      logPrint("Coarse convergence bound: " , coarse_b);
    } 
    dpw_encode_and_enforce(newBound);
  
}


// Sets according to current maxweight
void CBLIN::setPBencodings() {
  
  if (bestModel.size() < maxsat_formula->nVars() || solverCad->status() != 10) {
      logPrint("Extending best model to full formula");
      extendBestModel();
  }
  nbCurrentSoft = 0; 
  max_coeff_nothardened_sis = 0;
  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    uint64_t reducedWeight = maxsat_formula->getSoftClause(i).weight / maxsat_formula->getMaximumWeight();
    if (reducedWeight > 0) { //i.e. if it wasnt hardened in PMRES step OR left out by varres. 
          nbCurrentSoft++;
          if (reducedWeight > max_coeff_nothardened_sis) {
            max_coeff_nothardened_sis = reducedWeight;
          }
      }
    }
  // nRealSoft() = maxsat_formula->nSoft() - num_hardened
  logPrint("there are " , nbCurrentSoft, " of ", nRealSoft(),  " objective lits on this precision with maxcoeff " , max_coeff_nothardened_sis);
  
  auto lambda = [this](Lit l){return literalTrueInModel(l, bestModel);};

  uint64_t reduced_cost = computeCostReducedWeights(&lambda); 
  if (reduced_cost == 0 && maxsat_formula->getMaximumWeight() > 1) {
      update_SIS_precision();
      setPBencodings(); 
      return; 
  }
  logPrint("building new PB");
  initializePBConstraint(reduced_cost); 
}

void CBLIN::initializePBConstraint(uint64_t rhs) {
  build_objective_func_and_coeffs();

  uint64_t red_gap = known_gap / maxsat_formula->getMaximumWeight();

  bool ls_feasible = false;
  if (ls_dyn_prec) {
    for (unsigned int scIdx = 0; scIdx < maxsat_formula->nSoft(); ++scIdx) {
      const auto boumsIdx = boums_clause_map.ex2In[maxsat_formula->nHard() + scIdx]; // hard clauses were added first
      assert(!BouMS_wcnf_isClauseHard(boums_inst.clauses + boumsIdx));
      boums_inst.clauses[boumsIdx].weight =
        maxsat_formula->getSoftClause(scIdx).weight / maxsat_formula->getMaximumWeight();
    }

    ls_feasible = localsearch(*ls_usual_init_assign);
  }

  bool ls_improved = false;
  const auto lambda = [this](Lit l){return literalTrueInModel(l, bestModel);};
  uint64_t min_cost = computeCostReducedWeights(&lambda);
  if (min_cost < rhs) {
    if (ls_feasible && skip_local_search) {
      ls_improved = true;
      logPrint("LS found better global UB and RHS for PB, old RHS: ", rhs, ", new RHS: ", min_cost);
    }
    rhs = min_cost;
  } else if (ls_feasible && ls_dyn_prec) {
    const auto lambda = [this](Lit l) { return boums_assignment[var(l)] != sign(l); };
    min_cost = computeCostReducedWeights(&lambda);
    if (min_cost < rhs) {
      ls_improved = true;
      logPrint("LS found better RHS for PB, old RHS: ", rhs, ", new RHS: ", min_cost);
      rhs = min_cost;
    }
  }
  if (ls_improved && ls_merge_assign) {
    const auto num_disagree = mergeAssignments<bool*, bool>(ls_merged_assign, boums_assignment, tolbool);
    logPrint("Merged dyn. prec. LS assignment, agreeing variables: ", boums_inst.numVariables - num_disagree,
             ", disagreeing variables: ", num_disagree);
  }
  
  if (red_gap < rhs) {
      logPrint("Setting rhs to reduced gap " + std::to_string(red_gap));
      rhs = red_gap;
  }    
  
  
  // if the bound is obtained from preprocessing, we can not set variables in encoding according to a model. 
  bool bound_set_by_prepro = false;
  if (do_preprocess) {
    uint64_t red_p_gap = (ub_prepro - lbCost) / maxsat_formula->getMaximumWeight();
    if (rhs > red_p_gap) {
        logPrint("reduced cost from preprocessor gap: " ,red_p_gap, " better than best model " ,rhs);
        rhs = red_p_gap;
        bound_set_by_prepro = true;
    }
  }

  if (rhs == 0 && maxsat_formula->getMaximumWeight() > 1) {
      update_SIS_precision();
      setPBencodings(); 
      return; 
  }

  logPrint("encoding PB with UB: " ,rhs, " obj size: " ,nbCurrentSoft, " precision: " ,maxsat_formula->getMaximumWeight());


  if (incremental_DPW) {
      assert(dpw != NULL);
      RustSAT::dpw_set_precision(dpw, maxsat_formula->getMaximumWeight());
  }
  else {
    if (dpw != NULL) {
      RustSAT::dpw_drop(dpw);
      dpw = NULL;
    }
    dpw = RustSAT::dpw_new();
    for (int i = 0; i < objFunction.size(); i++) {
        RustSAT::dpw_add(dpw, lit2Int(objFunction[i]), coeffs[i]);
      }
  }
  dpw_encode_and_enforce(rhs);
  logPrint("Encoding Done");        
  setCardVars(bound_set_by_prepro);
}

void CBLIN::build_objective_func_and_coeffs_prop() {
  objFunction.clear();
  coeffs.clear();
  if (verbosity > 1) {
          cout << "c objective and weights;"; 
        }
  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    uint64_t weight = maxsat_formula->getSoftClause(i).weight;
    if (weight > 0) { //i.e. if it wasnt hardened in PMRES step OR left out by varres. 
        Lit l = maxsat_formula->getSoftClause(i).clause[0]; 
        assert (l != lit_Undef);
        objFunction.push(~l);
        coeffs.push(weight);
        if (verbosity > 1) {
          cout << " l:" << lit2Int(~l) << ",w:" << weight; 
        }
    }
    
  }
  if (verbosity > 1) {
          cout << endl; 
        }
}


void CBLIN::build_objective_func_and_coeffs() {
  if (incremental_DPW) {
    return; // in incremental mode, all objective literals are collected and added in the beginning of SIS.
  }
  objFunction.clear();
  coeffs.clear();

  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    uint64_t reducedWeight = maxsat_formula->getSoftClause(i).weight / maxsat_formula->getMaximumWeight();

    if (reducedWeight > 0) { //i.e. if it wasnt hardened in PMRES step OR left out by varres. 
        Lit l = maxsat_formula->getSoftClause(i).assumption_var; 
          assert (l != lit_Undef);
          objFunction.push(l);
          coeffs.push(reducedWeight);
      }
    }
  maxsat_formula->setProblemType(_WEIGHTED_);
  
}

void CBLIN::setCardVars(bool prepro_bound) {
    if (!extend_models) {
      return;
    }
    logPrint("Setting Card Vars currently: " + std::to_string(solverCad->vars()) + " / orig " + std::to_string(isSoft.size()));
    vec<Lit> cardAssumps;

    if (!prepro_bound) {
      assert(isSoft.size() <= bestModel.size());
      for (int i = 0; i < isSoft.size(); i++ ) {
        if (isSoft[i]) continue; 
        Lit l = mkLit(i, false);
        if (literalTrueInModel(l, bestModel)) {
          cardAssumps.push(l);
        }
        else {
          cardAssumps.push(~l);
        }
        
      }
    }
    logPrint("assumps - done");

    lbool res = ICadical::searchSATSolver(solverCad, cardAssumps);
    
    logPrint("Cadical call done");
    has_flipped = false;

    if (res == l_False) {
      logPrint("Warning: UNSAT in card setting");
      return;
    }

    assert(res == l_True);
    checkModel(false, true);
    savePhase();
    logPrint("CardVars DONE  ");
}

void CBLIN::extendBestModel() {

   logPrint("Debug: extending, current UB: " + std::to_string(ubCost) + " size of best model " + std::to_string(bestModel.size()));
   logPrint("Debug: Variables in formula: " + std::to_string(maxsat_formula->nVars()) + " variables in cadical " + std::to_string(solverCad->vars()));

    vec<Lit> modelAssumps;

    for (int i = 0; i < isSoft.size(); i++ ) {
      if (!isSoft[i]) continue;
      Lit l = mkLit(i, true); 
      if (literalTrueInModel(l, bestModel)) {
        logPrint("Pushing: " + std::to_string(lit2Int(l)));
        modelAssumps.push(l);
      }     
      else {
        logPrint("Pushing: " + std::to_string(lit2Int(~l)));
        modelAssumps.push(~l);
      }
      //if (isSoft[i]) continue;
      //modelAssumps.push(mkLit(i,  bestModel[i] == l_False));
    }
    lbool res =  ICadical::searchSATSolver(solverCad, modelAssumps);
    has_flipped = false;
    assert(res == l_True);
    assert(solverCad->status() == 10);
    checkModel(false, true);
    if (ls_extend) {
      localsearch(bestModel);
    }
    logPrint("Debug: after extending, current UB: " + std::to_string(ubCost) + " size of best model " + std::to_string(bestModel.size()));
}

bool CBLIN::localsearch(vec<lbool> & sol) {
    if (boums_broken) {
      logPrint("BouMS is currently unusable due to a previous error");
      return false;
    }

    BouMS_result_t boums_result;
    boums_result.assignment = boums_assignment;
  
    // BouMS: set initial assignment
    for (int i = 0; i < maxsat_formula->nVars(); ++i) {
      boums_assignment[i] = sol[i] == l_True ? true : false;
    }

    {
      // BouMS: solve
      const auto num_clauses = boums_inst.numClauses;
      const auto num_hard_clauses = boums_inst.numHardClauses;
      const bool boums_stop_dummy = false;
      BouMS_cores_solve(&boums_inst, boums_cores, &boums_params, boums_mem, &boums_mem_req, &boums_result,
                        boums_assignment, &boums_clause_map, boums_params.maxFlips, &boums_stop_dummy);
      // make sure we don't lose clauses e.g., when their weights are set to 0
      boums_inst.numClauses = num_clauses;
      boums_inst.numHardClauses = num_hard_clauses;
    }

    if (boums_result.status == BOUMS_UNKNOWN || boums_result.status == BOUMS_OPTIMUM_FOUND) {
      vec<lbool> local_search_best(maxsat_formula->nVars());
      for (int i = 0; i < maxsat_formula->nVars(); ++i) {
        local_search_best[i] = boums_result.assignment[i] ? l_True : l_False;
      }

      auto lambda = [&local_search_best, this](Lit l){return literalTrueInModel(l, local_search_best);};
      uint64_t local_search_cost =  computeCostOfModel(&lambda);
      if (local_search_cost < ubCost) {
        vec<Lit> local_search_model;
        for (int i = 0; i < maxsat_formula->nVars(); i++ ) {
          Lit l = mkLit(i, true); 
          if (literalTrueInModel(l, local_search_best)) {
            local_search_model.push(l);
          }     
          else {
            local_search_model.push(~l);
          }
        }
        lbool res = ICadical::searchSATSolver(solverCad, local_search_model); 
        assert(res == l_True);
        checkModel(true, true);
      }
      return true;
    }
    else {
      logPrint("Local search found no solution");
      return false;
    }
}

void CBLIN::flipValueinBestModel(Lit l) {
  assert(var(l) < bestModel.size());
  if (bestModel[var(l)] == l_True) bestModel[var(l)] = l_False;
  else if (bestModel[var(l)] == l_False) bestModel[var(l)] = l_True;
}


void CBLIN::minimizelinearsolution(vec<lbool> & model) {

  if (ls_min) {
    if (!skip_local_search) {
       localsearch(model);
    }
    return;
  }

  if (objFunction.size() == 0) {
    return;
  }

  vec<Lit> minimizable; 
  vec<bool> skip; 
  vec<Lit> fixed_assumptions; 
  time_t rec = time(NULL);
  for (int i = 0; i < isSoft.size(); i ++) {
    if (isSoft[i]) continue; 
    Lit l = mkLit(i, true); 
    if (literalTrueInModel(l, model)) {
      fixed_assumptions.push(l);
    }
    else  {
      fixed_assumptions.push(~l);
    }
  }

  for (int i = 0; i < objFunction.size(); i++) {
    Lit l = objFunction[i]; 
    assert(var(l) >= isSoft.size() || isSoft[var(l)] );
    if (literalTrueInModel(l, model)) {
      // induces cost
      minimizable.push(l);
      skip.push(false);
    }
    else {
      fixed_assumptions.push(~l);
    }
  }

  vec<Lit> assumps; 
  lbool res = l_False; 

  for (int i = 0; i < minimizable.size() ; i ++) {
    if (skip[i]) continue;
    Lit l = minimizable[i];
    assumps.clear();
    fixed_assumptions.copyTo(assumps);
    assumps.push(~l);
    res =  ICadical::searchSATSolver(solverCad, assumps);
    has_flipped = false;
    if (res == l_True) {
      fixed_assumptions.push(~l);
      for (int j = i+1; j < minimizable.size(); j++) {
        if (skip[j]) continue;
        Lit n = minimizable[j];
        if (literal_sat_in_cadical(~n)) {
          skip[j] = true; 
          fixed_assumptions.push(~n);
        }
      }
    } else if (res == l_False) {
      fixed_assumptions.push(l);
    } else {
      logPrint("undef in model minimisation");
      exit(_ERROR_);
    }
  }
  
  if (res == l_False) {
    res = ICadical::searchSATSolver(solverCad, fixed_assumptions);
    has_flipped = false;
    assert(res == l_True);
  }
  checkModel(false, true);
  time_t done = time(NULL);
  logPrint("minimization time " , done - rec, " init minsize " , minimizable.size());

} 

// Set assumptions on what soft clauses to consider, saves rebuilding the solver; 
void CBLIN::setAssumptions(vec<Lit> &assumps) {
    nbCurrentSoft = 0;
    assumps.clear();
    for (int i = 0; i < softs_added ; i++) {
      assert(maxsat_formula->getSoftClause(i).clause.size() == 1);
      
      bool shouldAdd = 
            (maxsat_formula->getSoftClause(i).weight >= maxsat_formula->getMaximumWeight() ) ||
            (maxsat_formula->getSoftClause(i).weight / maxsat_formula->getMaximumWeight()  > 0 ); 

      if (shouldAdd) {
        Lit l = maxsat_formula->getSoftClause(i).assumption_var;
        assert(l != lit_Undef); 
        assumps.push(~l);
        nbCurrentSoft++;     
      }
    }
  }



// Public search method
StatusCode CBLIN::search() {
  if (weightStrategy == _WEIGHT_NONE_) {
    logPrint("changing weight strategy to normal");
    weightStrategy = _WEIGHT_NORMAL_;
  }
  logPrint("parameters");
  logPrint("linear_strat=", lins);
  logPrint("relax_before_strat=", relaxBeforeStrat);
  logPrint("incremental_varying_res_GTE=", incrementalVarres);
  logPrint("precision_varres=" , non_inc_precision);
  logPrint("incremental_DPW=" ,incremental_DPW);
  logPrint("dpw_coarse=" , dpw_coarse);
  logPrint("minimize_sol=" , minimize_sol);
  logPrint("minimize_strat=" , minimize_strat);
  logPrint("ls_init_level=", ls_init_level);
  logPrint("ls_min=", ls_min);
  logPrint("ls_dyn_prec=", ls_dyn_prec);
  logPrint("ls_sis=", ls_sis);
  logPrint("ls_merge_assign=", ls_merge_assign);
  logPrint("ls_cores=", ls_cores);
  logPrint("zero_weight_core_fact=", zero_weight_core_fact);
  logPrint("ls_learn_clauses=", ls_learn_clauses);
  logPrint("ls_learn_clauses_fact=", ls_learn_clauses_fact);
  logPrint("ls_extend=", ls_extend);
  logPrint("ls_replace_sis=", ls_replace_sis);

  logPrint("Before search: UB ", ubCost, " LB ", lbCost, " off_set ", off_set, 
            " standardization_removed ", standardization_removed, " preprocessing_removed ", 
                cost_removed_preprocessing);

  time_start = time(NULL);
	time_best_solution = time_start;

  StatusCode r = setup(); 

  if (r == _UNSATISFIABLE_) {
         logPrint("clauses unsat, no solutions");
         return _UNSATISFIABLE_;
  }
  if (r == _OPTIMUM_) {
    ubCost = lbCost;
    printBound(ubCost);
    printAnswer(_OPTIMUM_);
    return _OPTIMUM_;
  }

  if (ls_init_level > 0) {
    logPrint("Running LS on preprocessed instance");
    localsearch(bestModel);

    // use current model as starting point for ls on original
    vec<lbool>* m;
    if (do_preprocess) {
      model_of_original.clear();
      reconstruct_model_prepro(model, model_of_original);
      m = &model_of_original;
    } else {
      m = &model;
    }

    if (ls_init_level > 1) {
      // run LS on original input instance to potentially improve UB
      auto tmp_maxsat_formula = maxsat_formula;
      maxsat_formula = orig_maxsat_formula;
      updateBouMSInstance();

      if (!boums_broken) {
        for (BouMS_uint_t vIdx = 0; vIdx < boums_inst.numVariables; ++vIdx) {
          boums_assignment[vIdx] = (*m)[vIdx] == l_True;
        }

        BouMS_result_t res;
        res.assignment = boums_assignment;

        logPrint("Running LS on original instance");
        const bool dummy = false;
        BouMS_solve(&boums_inst, &boums_params, boums_mem, &boums_mem_req, &res, boums_assignment, &boums_clause_map,
                    boums_params.maxFlips, &dummy);
        assert(res.status == BOUMS_OPTIMUM_FOUND || res.status == BOUMS_UNKNOWN);
        if (res.cost < ubCost) {
          logPrint("LS on original found better UB, old: ", ubCost, ", new: ", res.cost);
          init_ls_ub_assign = new bool[boums_inst.numVariables];
          if (init_ls_ub_assign) {
            memcpy(init_ls_ub_assign, res.assignment, boums_inst.numVariables * sizeof(bool));
            init_ls_ub = res.cost;
            printBound(init_ls_ub);
          } else {
            logPrint("Could not save LS model, OOM");
          }
        }
      }
      maxsat_formula = tmp_maxsat_formula;
      updateBouMSInstance();
    }
  }

  switch (lins) {
    case 0:
      return weightSearch();
      break;
    
    case 1:
      return coreGuidedLinearSearch();
      break;

    case 2: 
      return linearSearch();
      break;

    default:
      logPrint("Error: Invalid variation value.");
      cout << "s UNKNOWN" << std::endl;
      exit(_ERROR_);
    }
}

/************************************************************************************************
 //
 // Other protected methods
 //
 ************************************************************************************************/

/*_________________________________________________________________________________________________
  |
  |  initAssumptions :  [void]
  |
  |  Description:
  |
  |    Defines the new assumption literal for each soft clause  Assumptions are used to
  |    extract cores. Assumes all soft clauses are of length 1.
  |  Post-conditions:
  |    * Map the literal in each soft clause as the assumption. I.e. set 'softClauses[i].assumptionVar' to equal the negation of the 
        literal in the clause.
  |    * 'coreMapping' is updated by mapping each assumption literal with the
  |      corresponding index of each soft clause.
  |    * original weights tracks the initial weights of the formula
  |
  |________________________________________________________________________________________________@*/
void CBLIN::initAssumptions() {
  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    assert(maxsat_formula->getSoftClause(i).clause.size() == 1);
    Lit l = maxsat_formula->getSoftClause(i).clause[0];
    maxsat_formula->getSoftClause(i).assumption_var = ~l;
    coreMapping[~l] = i;    
  }
}



void CBLIN::printProgress() {
  logPrint(inLinSearch ? "LIN " : "CG ", "UB " , ubCost, " LB " , lbCost, " time " , time_best_solution - time_start );  
}

time_t CBLIN::timeSinceStart() {
  time_t cur = time(NULL);
  return cur - time_start;
}

time_t CBLIN::timeSincePrepro() {
  time_t cur = time(NULL);
  return cur - time_prepro;
}

void CBLIN::addSoftClauseAndAssumptionVar(uint64_t weight, vec<Lit> &clause) {
    Lit l = maxsat_formula->newLiteral();
    clause.push(l);
    maxsat_formula->addHardClause(clause);

    clause.clear();
    clause.push(~l);
    maxsat_formula->addSoftClause(weight, clause);

    maxsat_formula->getSoftClause(maxsat_formula->nSoft() - 1).assumption_var = l;
    coreMapping[l] = maxsat_formula->nSoft() - 1;  // Map the new soft clause to its assumption literal.
}

int CBLIN::nRealSoft() {
  return maxsat_formula->nSoft() - num_hardened;
}  


void CBLIN::resetSolver(int mode) { //  0 default 1 unsat 2 sat
    logPrint("Deleting solver");
    delete solverCad;

    solverCad = ICadical::newSATSolver();

    if (dynamic_cadical) {
      ICadical::set_mode(mode, solverCad);
    }

    clauses_added = 0;
    softs_added = 0;
    vars_added = 0;
    return updateSolver();
} 

bool CBLIN::shouldUpdate() {
  return clauses_added < maxsat_formula->nHard();
}

// save polarity from last model 
 void CBLIN::savePhase() {
    logPrint("save phase");
		for (int i = 0; i < bestModel.size(); i++){
      Lit l = mkLit(i, false);
      if (isSoft[var(l)] && optimistic) continue;

      if (literalTrueInModel(l, bestModel)) {
        solverCad->phase(lit2Int(l));
      }
      else {
        solverCad->phase(lit2Int(~l));
      }
		}
    if (optimistic) {
      softsSatisfied();
    }
 }

 void CBLIN::softsSatisfied() {
    for (int i = 0; i < original_labels->nSoft(); i++) {
      assert(original_labels->getSoftClause(i).clause.size() == 1); 
      Lit l = original_labels->getSoftClause(i).clause[0];
      int indexSoft = coreMapping[l];

      if (maxsat_formula->getSoftClause(indexSoft).weight > 0) {
        solverCad->phase(MaxSAT::lit2Int(l));
      }
    } 
 }

 void CBLIN::clearFixingsonSoft() {
    for (int i = 0; i < maxsat_formula->nSoft(); i++) {
        assert(maxsat_formula->getSoftClause(i).clause.size() == 1 );
        Lit l =  maxsat_formula->getSoftClause(i).clause[0];
        solverCad->unphase(lit2Int(l));
    }

 }

  void CBLIN::freezeObjective() {
    logPrint("freeze literals");
    for (int i = 0; i < original_labels->nSoft(); i++) {
      Lit l = original_labels->getSoftClause(i).clause[0];
      solverCad->freeze(MaxSAT::lit2Int(l));
    } 
  }

  bool CBLIN::flipLiterals() {
    assert(solverCad->status() == 10 );
    if (has_flipped) {
      return false;
    }
//    return false;
    uint64_t flips = 0;
    uint64_t failed_flips = 0;
    for (int i = 0; i < original_labels->nSoft(); i++) {
      Lit l = original_labels->getSoftClause(i).clause[0];
      if(!literal_sat_in_cadical(l)) {
        if (solverCad->flip(lit2Int(l))) {
          flips++;
        }
        else {
          failed_flips++;
        }
      }
    } 
    logPrint("Flipped " + std::to_string(flips) + " failed flips " + std::to_string(failed_flips));
    has_flipped = true;
    return (flips > 0);
  }

//TODO parametrize on the model... 
bool CBLIN::checkModel(bool from_local_search, bool improve_better) {
  flipLiterals();
  auto lambda = [this](Lit l){ return literal_sat_in_cadical(l) ;};

  uint64_t modelCost = computeCostOfModel(&lambda);
  logPrint("cost of new ", modelCost, " ub ", ubCost);
  
  bool isBetter = modelCost < ubCost;
  if (isBetter) {
        vec<lbool> cadModel; 
        ICadical::getModel(solverCad, cadModel);
        ubCost = modelCost;
        time_best_solution = time(NULL);
        printProgress();
        saveModel(cadModel);
        bestModel.clear();
        cadModel.copyTo(bestModel);
        printBound(ubCost);
        checkGap();
        skip_local_search = from_local_search;
    }
  else if (improve_better && (modelCost == ubCost) && solverCad->vars() >= bestModel.size()) {
      vec<lbool> cadModel; 
      ICadical::getModel(solverCad, cadModel);
      logPrint("Found same cost model covering more or the same amount of variables");
      saveModel(cadModel);
      bestModel.clear();
      cadModel.copyTo(bestModel);
      isBetter = true;
  }
  if (isBetter && inLinSearch) {
        savePhase();
  }
  return isBetter;
 }

 void CBLIN::checkGap() {
   uint64_t currentGap = ubCost - lbCost;
   if (currentGap < known_gap) {
     known_gap = currentGap;
     if (inLinSearch)
        logPrint("LIN gap ", known_gap , " at " , timeSinceStart());
     else 
        logPrint("CG gap ", known_gap, " at ",  timeSinceStart());
   }
 }

 bool  CBLIN::literal_sat_in_cadical(Lit l) {
        assert(solverCad->status() == 10);
        int lit = MaxSAT::lit2Int(l); 
        if (lit > 0) {
          return solverCad->val(lit) == lit;
        }
        else {
          return solverCad->val(lit) != lit;
        } 
 }

void CBLIN::updateBouMSInstance() {
  boums_broken = false;

  const auto oldNumVars = boums_inst.numVariables;
  const auto oldNumClauses = boums_inst.numClauses;
  BouMS_wcnf_util_deleteFormula(&boums_inst, free, NULL);
  boums_inst.numClauses = 0;
  boums_inst.numHardClauses = 0;
  boums_inst.numVariables = 0;

  BouMS_wcnf_util_batchClauseAddingState_t boums_clause_adder;
  if (BouMS_wcnf_util_startBatchClauseAdding(maxsat_formula->nHard() + maxsat_formula->nSoft(), realloc, free,
                                           &boums_clause_adder)) {
    logPrint("Error initializing BouMS clause adding");
    boums_broken = true;
    return;
  }

  const auto convert_clause = [](const vec<Lit>& clause, vec<int>& converted) {
    converted.growTo(clause.size());
    for (int litIdx = 0; litIdx < clause.size(); ++litIdx) {
      const auto& lit = clause[litIdx];
      converted[litIdx] = (sign(lit) ? -1 : 1) * (var(lit) + 1);
    }
  };
  vec<int> converted_clause;

  bool oom = false;
  for (int clauseIdx = 0; clauseIdx < maxsat_formula->nHard() && !oom; ++clauseIdx) {
    const auto& clause = maxsat_formula->getHardClause(clauseIdx).clause;
    convert_clause(clause, converted_clause);
    if (BouMS_wcnf_util_batchAddClause(&boums_clause_adder, BOUMS_HARD_CLAUSE_WEIGHT, &converted_clause[0],
                                       clause.size())) {
      oom = true;
    }
  }

  for (int clauseIdx = 0; clauseIdx < maxsat_formula->nSoft() && !oom; ++clauseIdx) {
    const auto& clause = maxsat_formula->getSoftClause(clauseIdx);
    convert_clause(clause.clause, converted_clause);
    if (BouMS_wcnf_util_batchAddClause(&boums_clause_adder, clause.weight, &converted_clause[0],
                                       clause.clause.size())) {
      oom = true;
    }
  }

  const auto numCores = cores.size();
  if (numCores > 0 && ls_cores == 1) {
    for (int coreIdx = 0; coreIdx < numCores && !oom; ++coreIdx) {
      const auto& core = cores.at(coreIdx);
      const auto coreSize = core.numLiterals;
      std::vector<int> copy;
      copy.reserve(coreSize);
      for (int lIdx = 0; lIdx < coreSize; ++lIdx) {
        const auto lit = core.literals[lIdx];
        const int var = BouMS_var(lit) + 1;
        const bool sign = BouMS_sign(lit);
        copy.push_back(sign ? -var : var);
      }
      if (BouMS_wcnf_util_batchAddClause(&boums_clause_adder, BOUMS_HARD_CLAUSE_WEIGHT, copy.data(), copy.size())) {
        oom = true;
      }
    }
    if (!oom) {
      logPrint("Added ", numCores, " core clauses to BouMS");
    }
  }

  const auto numLearnedClauses = learned_clauses.size();
  if (numLearnedClauses > 0) {
    size_t maxClauseSize = 0;
    for (int clauseIdx = 0; clauseIdx < numLearnedClauses && !oom; ++ clauseIdx) {
      const auto& clause = learned_clauses.at(clauseIdx);
      const auto clauseSize = clause.size();
      if (clauseSize > maxClauseSize) {
        maxClauseSize = clauseSize;
      }
      if (BouMS_wcnf_util_batchAddClause(&boums_clause_adder, BOUMS_HARD_CLAUSE_WEIGHT, clause.data(), clauseSize)) {
        oom = true;
      }
    }
    if (!oom) {
      logPrint("Added ", numLearnedClauses, " learned clauses to BouMS (max ", small_clause_learner.max_num_clauses,
               "), max clause size = ", maxClauseSize);
    }
  }

  if (oom) {
    boums_broken = true;
    BouMS_wcnf_util_cleanUpBatchClauseAddingAfterError(&boums_clause_adder);
    logPrint("Error adding clauses to BouMS");
    return;
  } else if (BouMS_wcnf_util_finishBatchClauseAdding(&boums_clause_adder, &boums_inst, NULL)) {
    boums_broken = true;
    BouMS_wcnf_util_cleanUpBatchClauseAddingAfterError(&boums_clause_adder);
    BouMS_wcnf_util_deleteFormula(&boums_inst, free, NULL);
    logPrint("Error adding clauses to BouMS");
    return;
  }

  const auto new_boums_bytes = BouMS_calcMemoryRequirements(&boums_inst, &boums_mem_req);
  if (new_boums_bytes > boums_bytes || boums_mem == NULL) {
    boums_bytes = new_boums_bytes;
    free(boums_mem);
    boums_mem = malloc(boums_bytes);
    if (!boums_mem) {
      boums_broken = true;
      BouMS_wcnf_util_deleteFormula(&boums_inst, free, NULL);
      logPrint("Error allocating memory for BouMS");
      return;
    }
  }

  if (oldNumVars < boums_inst.numVariables || boums_assignment == NULL) {
    delete[] boums_assignment;
    boums_assignment = new bool[boums_inst.numVariables];
    if (!boums_assignment) {
      boums_broken = true;
      free(boums_mem);
      BouMS_wcnf_util_deleteFormula(&boums_inst, free, NULL);
      logPrint("Error allocating memory for BouMS' assignment");
      return;
    }
  }

  if (oldNumClauses < boums_inst.numClauses ||  boums_clause_map.ex2In == NULL || boums_clause_map.in2Ex == NULL) {
    if (boums_clause_map.ex2In != NULL) {
      delete[] boums_clause_map.ex2In;
      boums_clause_map.ex2In = NULL;
    }
    if (boums_clause_map.in2Ex != NULL) {
      delete[] boums_clause_map.in2Ex;
      boums_clause_map.in2Ex = NULL;
    }

    boums_clause_map.ex2In = new BouMS_uint_t[boums_inst.numClauses];
    boums_clause_map.in2Ex = new BouMS_uint_t[boums_inst.numClauses];

    if (boums_clause_map.ex2In == NULL || boums_clause_map.in2Ex == NULL) {
      boums_broken = true;
      if (boums_clause_map.ex2In != NULL) {
        delete[] boums_clause_map.ex2In;
        boums_clause_map.ex2In = NULL;
      }
      if (boums_clause_map.in2Ex != NULL) {
        delete[] boums_clause_map.in2Ex;
        boums_clause_map.in2Ex = NULL;
      }
      delete[] boums_assignment;
      free(boums_mem);
      BouMS_wcnf_util_deleteFormula(&boums_inst, free, NULL);
      logPrint("Error allocating memory for BouMS' clause map");
      return;
    }

    for (unsigned int cIdx = 0; cIdx < boums_inst.numClauses; ++cIdx) {
      boums_clause_map.ex2In[cIdx] = cIdx;
      boums_clause_map.in2Ex[cIdx] = cIdx;
    }
  }

}

void CBLIN::loadFormula(MaxSATFormula *maxsat) {
  MaxSAT::loadFormula(maxsat);
  orig_maxsat_formula = maxsat_formula->copyMaxSATFormula();
  small_clause_learner.max_num_clauses = ls_learn_clauses_fact * orig_maxsat_formula->nHard();
  if (small_clause_learner.max_num_clauses == 0) {
    small_clause_learner.max_num_clauses = 1;
  }
}

void CBLIN::setup_formula() {
  MaxSAT::setup_formula();
  updateBouMSInstance();
  BouMS_params(&boums_inst, &boums_params);
  // we always provide an initial assignment that we don't want to replace by a random one ever (?)
  boums_params.maxTriesWOImprovement = BOUMS_UINT_MAX;
  boums_params.zeroWeightCoreFactor = zero_weight_core_fact;
  if (ls_cores) {
    if (ls_cores == 2) {
      boums_params.coreWeightingMode = 0;
    } else if (ls_cores == 3) {
      boums_params.coreWeightingMode = 1;
    }
  }
}

void CBLIN::printAnswer(int type) {
  if (init_ls_ub < ubCost) {
    ubCost = init_ls_ub;
    std::cout << "s SATISFIABLE" << std::endl;
    printBound(ubCost);
    std::stringstream s;
    s << "v ";
    for (unsigned int vIdx = 0; vIdx < orig_maxsat_formula->nVars(); ++vIdx) {
      s << (init_ls_ub_assign[vIdx] ? "1" : "0");
    }
    s << std::endl;
    std::cout << s.str();
  } else {
    MaxSAT::printAnswer(type);
  }
}

CBLIN::SmallClauseLearner::SmallClauseLearner(size_t max_num_clauses)
  : max_num_clauses(max_num_clauses)
  , clauses()
  , current_clause()
{}

bool CBLIN::SmallClauseLearner::learning(int size) {
  return clauses.size() < max_num_clauses || size <= clauses.top().size();
}

void CBLIN::SmallClauseLearner::learn(int lit) {
  if (lit != 0) {
    current_clause.push_back(lit);
    return;
  }

  if (current_clause.size() > 0) {
    if (clauses.size() == max_num_clauses) {
      clauses.pop();
    }
    clauses.push(std::move(current_clause));
    current_clause.clear();
  }
}

CBLIN::SmallClauseLearner::container_type CBLIN::SmallClauseLearner::extract_clauses() {
  container_type c;
  c.reserve(clauses.size());
  while (!clauses.empty()) {
    c.push_back(std::move(clauses.top()));
    clauses.pop();
  }
  return c;
}

size_t CBLIN::SmallClauseLearner::size() const {
  return clauses.size();
}
