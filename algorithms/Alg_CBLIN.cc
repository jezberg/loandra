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
Solver *CBLIN::updateSolver() {

  reserveSATVariables(solver, maxsat_formula->nVars());

  for (int i = vars_added; i < maxsat_formula->nVars(); i++)
    newSATVariable(solver);

  vars_added = maxsat_formula->nVars();

  for (int i = clauses_added; i < maxsat_formula->nHard(); i++)
    solver->addClause(maxsat_formula->getHardClause(i).clause);
  
  clauses_added = maxsat_formula->nHard();

  softs_added = maxsat_formula->nSoft();


  //We do not support these
  assert(maxsat_formula->nPB() == 0); 
  //We do not support these
  assert(maxsat_formula->nCard() == 0);
  return solver;
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

  Solver *CBLIN::hardenClauses() { 
	   uint64_t gap = ubCost - lbCost;
     logPrint("hardening at gap: " , gap);
     int num_hardened_round = 0;
	   maxw_nothardened = 0;
	   for (int i = 0; i < softs_added; i++)
		  {
      Lit l =  maxsat_formula->getSoftClause(i).clause[0];
      assert(l != lit_Undef);
			if (maxsat_formula->getSoftClause(i).weight > gap || (maxsat_formula->getSoftClause(i).weight == gap && literalTrueInModel(l, bestModel)) ) {  // 
        
        assert(var(l) < solver->nVars());
        if (!literalTrueInModel(l, bestModel)) {
          flipValueinBest(l);
        }
      	vec<Lit> clause;
				clause.clear();
				clause.push(l);
				solver->addClause(clause);
        maxsat_formula->addHardClause(clause); 
        
				maxsat_formula->getSoftClause(i).weight = 0;
        maxsat_formula->getSoftClause(i).assumption_var = lit_Undef;

				num_hardened++;
				num_hardened_round++;
        did_harden = true;
			}
			else if (maxsat_formula->getSoftClause(i).weight > maxw_nothardened) {
				maxw_nothardened = maxsat_formula->getSoftClause(i).weight;
			} 

			
		}

    if (num_hardened_round) {
      updateBouMSInstance();
    }

		logPrint("hardened in total: " , num_hardened_round, " literals");
    logPrint("hardening again at gap ", maxw_nothardened);
		return solver;
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

  void CBLIN::hardenClausesSIS(uint64_t reduced_cost, vec<lbool> &currentModel) { 
      if (incremental_DPW) {
        logPrint("can not harden based on reduced cost with incremental DPW");
        return; 
      }
     int num_hardened_round = 0;
	   max_coeff_nothardened_sis = 0;
     uint64_t precision = maxsat_formula->getMaximumWeight();

     logPrint("hardening in SIS, reduced cost ", reduced_cost, " precision ", precision);

	   for (int i = 0; i < maxsat_formula->nSoft(); i++)
		  {
      if (maxsat_formula->getSoftClause(i).weight <= precision) {
        continue;
      }
      Lit l =  maxsat_formula->getSoftClause(i).clause[0];
      uint64_t red_weight = maxsat_formula->getSoftClause(i).weight / precision;
      assert(l != lit_Undef);
			if (red_weight > reduced_cost || (red_weight == reduced_cost && literalTrueInModel(l, currentModel)) ) {  // 
        assert(var(l) < solver->nVars());
      	vec<Lit> clause;
				clause.clear();
				clause.push(l);
				solver->addClause(clause);
				num_hardened_round++;
			}
			else if (red_weight > maxw_nothardened) {
				max_coeff_nothardened_sis = red_weight;
			} 
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

  void CBLIN::init_SIS_precision() {
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
    maxsat_formula->setMaximumWeight(weightCand);
    max_weight_after_cg = weightCand;
    logPrint("first precision for SIS " , maxsat_formula->getMaximumWeight());
  }

  void CBLIN::update_SIS_precision() {
    uint64_t precision_factor = precision_factors();
    uint64_t nextFactor;

    if (incremental_DPW && have_encoded_precision) {
      nextFactor = dpw_next_precision();
    }
    else {
     nextFactor = maxsat_formula->getMaximumWeight() / precision_factor;
    }
     
    while (moreThanWeight(nextFactor) == nbCurrentSoft && nextFactor > 1 ) {
      nextFactor /= precision_factor; 
    }
    maxsat_formula->setMaximumWeight(nextFactor);
    logPrint("new precision for SIS " , nextFactor);
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


  solver = updateSolver();

  softsSatisfied();
  lbool res = searchSATSolver(solver, assumptions);
  solver->resetFixes();

  if (res == l_False) {
    nbCores++;
    printAnswer(_UNSATISFIABLE_);
    return _UNSATISFIABLE_;
  } else if (res == l_True) {
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
      if(timeLimitCores > 0 && (time_t)timeLimitCores- timeSinceStart() <= 0 ) {
        return _UNKNOWN_;
      }
      if(timeLimitCores > 0) {
        logPrint("cg time remaining " , (time_t)timeLimitCores - timeSinceStart());
        solver->setTimeBudget(timeLimitCores- timeSinceStart());
      }
      setAssumptions(assumptions);
      lbool res; 
      res = searchSATSolver(solver, assumptions);

      if (res == l_Undef) {
        //Interrupted
        return _UNKNOWN_;
      }

      if (res == l_False) {
      
        nbCores++;
        assert(solver->conflict.size() > 0);
        uint64_t coreCost = computeCostCore(solver->conflict);
        lbCost += coreCost;
        checkGap();
        logPrint("LB ", lbCost, " core size ", solver->conflict.size(), " core-min-cost " , coreCost); 
        relaxCore(solver->conflict, coreCost);
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

      maxw_nothardened = 0;

      for (int i = 0; i < maxsat_formula->nSoft(); i++)  {
          assert( maxsat_formula->getSoftClause(i).clause.size() == 1);
          Lit l = maxsat_formula->getSoftClause(i).clause[0];
          assert(var(l) < isSoft.size());
          isSoft[var(l)] = true; 
          if ( maxsat_formula->getSoftClause(i).weight > maxw_nothardened) {
            maxw_nothardened = maxsat_formula->getSoftClause(i).weight;
          }
      }
      
      initAssumptions();  
      solver = newSATSolver();
      solver->setSolutionBasedPhaseSaving(false);
      StatusCode rs = unsatSearch();
      if (rs == _UNSATISFIABLE_) return rs;
      
      //Here we know that the formula is SAT
      if (maxsat_formula->nSoft() == 0 || ubCost == lbCost) {
          return _OPTIMUM_; //Solved by preprocessing
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
  assert(timeLimitCores < 0);
  
  for (;;) {
    StatusCode us = weightDisjointCores(); 

    //LB phase proves optimality, current model is not for the current formula. 
    if (us == _OPTIMUM_) {
        logPrint("LB = UB");
        return getModelAfterCG();
    }

    //At this point solver returned true and as such has a model
   
    nbSatisfiable++;
    auto lambda = [this](Lit l){ return literalTrueInModel(l, solver->model); };
    uint64_t modelCost = computeCostOfModel(&lambda);
    if (modelCost < ubCost) {
        ubCost = modelCost;
        saveModel(solver->model);
        printBound(ubCost);
    }
    if (lbCost == ubCost) {
      if (verbosity > 0)
        logPrint("LB = UB");          
        printAnswer(_OPTIMUM_);
        return _OPTIMUM_;
    }
    if (nbCurrentSoft == nRealSoft()) {
      assert(modelCost == lbCost);
      if (lbCost < ubCost) {
        ubCost = lbCost;
        saveModel(solver->model);
        printBound(lbCost);
      }
        printAnswer(_OPTIMUM_);
        return _OPTIMUM_;
    } 
    if (ubCost - lbCost < maxw_nothardened) {
      solver = hardenClauses();
    }
    if (shouldUpdate()) {
      solver = updateSolver();
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
  for (;;) {
    StatusCode us = weightDisjointCores(); 
    if (us == _OPTIMUM_) {
        logPrint("LB = UB");
        return getModelAfterCG();
    }

    if (us == _UNKNOWN_ ) {
        logPrint("interrupted core guided phase");
        if(shouldUpdate()) {
          logPrint("updating solver at ",  timeSinceStart());
          solver = updateSolver();
        }
        return linearSearch();
    }
    

    //At this point solver returned true and as such has a model
    assert(us == _SATISFIABLE_ );

    logPrint("SAT-During core guided phase at " , timeSinceStart());
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
        saveModel(solver->model);
        printBound(lbCost);
      }
        printAnswer(_OPTIMUM_);
        return _OPTIMUM_;
    } 


   //if code gets here algorithm cant terminate yet   
   if (ubCost - lbCost < maxw_nothardened) {
        solver = hardenClauses();
      }

   if(relaxBeforeStrat) {
      logPrint("Relax 2 Strat");
      if(shouldUpdate()) {
          logPrint("updating solver at ", timeSinceStart());
          solver = updateSolver();
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
          logPrint("updating solver at " , timeSinceStart());
          solver = updateSolver();
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
  solver = updateSolver();
  setAssumptions(assumptions);
  lbool res; 
  res = searchSATSolver(solver, assumptions);
  assert(res == l_True);

  auto lambda = [this](Lit l){ return literalTrueInModel(l, solver->model); };
  uint64_t modelCost = computeCostOfModel(&lambda);
  assert(modelCost == lbCost);
  if (lbCost < ubCost) {
    ubCost = lbCost;
    saveModel(solver->model);
  }
  printAnswer(_OPTIMUM_);
  return _OPTIMUM_;
}




StatusCode CBLIN::linearSearch() {
  logPrint( "Starting lin search with: LB: ",lbCost, " UB: " ,ubCost,
            " UB - LB: " ,ubCost-lbCost, " time " , timeSinceStart() );

  inLinSearch = true;
  solver->budgetOff();
  assumptions.clear();
  

  assert(bestModel.size() > 0);
  
  if(delete_before_lin) {
    solver = resetSolver();
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

  if (bestModel.size() < maxsat_formula->nVars() || !solver->okay() ) {
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

    if (!(incrementalVarres || use_DPW)) {
      assumptions.clear();
    }   
    logPrint("SAT Call at " , timeSinceStart(), " # assumptions " , assumptions.size(), " clauses in SAT solver " ,solver->nClauses());  
    if (verbosity > 1) {
      cout << "c assumptions:";
      for (int i = 0; i < assumptions.size(); i++) {
        cout << " " << lit2Int(assumptions[i]);
      }
      cout << endl;
    }

    res = searchSATSolver(solver, assumptions);

    if (res == l_True) {
      nbSatisfiable++;
      
      auto lambda = [this](Lit l){ return literalTrueInModel(l, solver->model); };
      uint64_t new_reduced_cost = computeCostReducedWeights(&lambda);
      bool better = checkModel(false, false);
      
      if (use_DPW && !incremental_DPW && assumptions.size() > 0) {
        vec<Lit> clause; 
        clause.push(assumptions[0]);
        solver->addClause(clause);

        if (harden_in_SIS && !incremental_DPW && new_reduced_cost  < max_coeff_nothardened_sis) {
          hardenClausesSIS(new_reduced_cost, solver->model);
        }

      }
      if (better && incremental_DPW) {
        harden_incremental();
      }

      if (ls_sis && !skip_local_search && new_reduced_cost > 0) {
        if (ls_merge_assign) {
          static const std::function<lbool(const lbool&)> lboolid = [](const lbool& b) { return b; };
          const auto num_disagree = mergeAssignments(ls_merged_assign, solver->model, lboolid);
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
        vec<lbool> temp; 
        solver->model.copyTo(temp);
        minimizelinearsolution(temp);
        if (minimize_strat == 2) {
          minimize_iteration = false;
        }
        auto lambda = [this, &temp](Lit l){ return literalTrueInModel(l, temp); };
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
            solver = resetSolver();
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
          logPrint("rebuilding after UNSAT");
          if (!(incrementalVarres || incremental_DPW)) {
            solver = resetSolver();
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
  SolverWithBuffer solver_with_buffer{.solver_b = solver, .clauses_added = 0, .verbosity = verbosity};
  int num_vars = solver->nVars();
  RustSAT::dpw_limit_range(dpw, 0, global_ub_dpw, &dpw_clause_collector, &solver_with_buffer);
  logPrint("hardening incremental DPW bound: " , global_ub_dpw, " clauses added " , solver_with_buffer.clauses_added) ;
}


uint64_t CBLIN::dpw_next_precision() {
  assert(have_encoded_precision);
  uint64_t next_prec = RustSAT::dpw_next_precision(dpw);
  have_encoded_precision = false; 
  return next_prec;
}

void CBLIN::dpw_encode_and_enforce(uint64_t rhs) {
    SolverWithBuffer solver_with_buffer{.solver_b = solver, .clauses_added = 0, .verbosity = verbosity};
    int num_vars = solver->nVars();
    RustSAT::dpw_encode_ub(dpw, rhs, rhs, &num_vars, &dpw_clause_collector, &solver_with_buffer);
    logPrint("clauses added in encode and enforce " , solver_with_buffer.clauses_added, " rhs " , rhs) ;
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

void CBLIN::dpw_clause_collector(int lit, void *ptr) {
  SolverWithBuffer *solver_with_buffer = (SolverWithBuffer *)ptr;
  if (lit) {
    while (solver_with_buffer->solver_b->nVars() < abs(lit))  solver_with_buffer->solver_b->newVar();
    solver_with_buffer->buffer.push(MaxSAT::int2Lit(lit));
    return;
  }
  if (solver_with_buffer->verbosity > 1) {
    cout << "c RUSTSAT clause:";
    for (int i = 0; i < solver_with_buffer->buffer.size(); i++) {
      cout << " " << MaxSAT::lit2Int(solver_with_buffer->buffer[i]); 
    }
    cout << endl;
  }
  solver_with_buffer->clauses_added += 1;
  solver_with_buffer->solver_b->addClause(solver_with_buffer->buffer);
  solver_with_buffer->buffer.clear();
}


void CBLIN::updateBoundLinSearch (uint64_t newBound) {  
  logPrint("new bound to enforce: " , newBound, " at ", timeSinceStart());
  
  if (use_DPW) {
    if (dpw_coarse) {
      uint64_t coarse_b = RustSAT::dpw_coarse_ub(dpw, newBound);
      dpw_fine_convergence_after = (coarse_b != newBound);
      fine_bound = newBound;
      newBound = coarse_b;
      logPrint("Coarse convergence bound: " , coarse_b);
    } 
    dpw_encode_and_enforce(newBound);
  }
  else{
    if (enc->hasPBEncoding()) {
      if(!incrementalVarres) {
        if (maxsat_formula->getProblemType() == _WEIGHTED_) {
          enc->updatePB(solver, newBound);
        } else {
          enc->updateCardinality(solver, newBound);
        }
      }
      else {
        assert(maxsat_formula->getProblemType() == _WEIGHTED_ );
        assumptions.clear();
        enc->updatePBA(assumptions, newBound);
      }
    }
    else {
      logPrint("no encoding");
      int added = 0;
      for (int i = 0 ; i < objFunction.size(); i ++) {
        if (coeffs[i] > newBound && coeffs[i] <= init_rhs) { //the second condition is here because literals that have coefficients higher than init_rhs are fixed to dfalse in the encoder
            if (!incrementalVarres) {
              solver->addClause({~objFunction[i]});
              added++;
            }
            else {
              assumptions.clear();
              assumptions.push(~objFunction[i]);
            }  
        }
      }
      assert(added > 0);
    }
  }
} 


// Sets according to current maxweight
void CBLIN::setPBencodings() {
  
  if (bestModel.size() < maxsat_formula->nVars()) {
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

  auto lambda = [this](Lit l){ return literalTrueInModel(l, bestModel); };
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

  if (use_DPW) {
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
  }
  else {
    if (enc != NULL)
      delete enc;
    enc = new Encoder(_INCREMENTAL_NONE_, _CARD_MTOTALIZER_,
                               _AMO_LADDER_, _PB_GTE_);
    assert(!enc->hasPBEncoding());
    enc->encodePB(solver, objFunction, coeffs, rhs);
    init_rhs = rhs; 
  }
  
  logPrint("Encoding done #assumptions " , assumptions.size());        
  setCardVars(bound_set_by_prepro);
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
    logPrint("setting Card Vars currently: " , solver->nVars(),  " / orig ", isSoft.size());
    solver->setSolutionBasedPhaseSaving(false);
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
    lbool res = searchSATSolver(solver, cardAssumps);
    if (res == l_False) {
      logPrint("Warning: UNSAT in card setting");
      //DEBUG
      //test_pb_enc();
      return;
    }
    assert(res == l_True);
    checkModel(false, true);
    solver->setSolutionBasedPhaseSaving(true);
}

/*
  After this method, solver->model() is a model of all of the variables in the SAT solver that matches 
  the best known model in terms of the original objective.
*/
void CBLIN::extendBestModel() {
    logPrint("extending best model to full formula");
    vec<Lit> modelAssumps;

    for (int i = 0; i < isSoft.size(); i++ ) {
      if (!isSoft[i]) continue;
      Lit l = mkLit(i, true); 
      if (literalTrueInModel(l, bestModel)) {
        modelAssumps.push(l);
      }     
      else {
        modelAssumps.push(~l);
      }
      //if (isSoft[i]) continue;
      //modelAssumps.push(mkLit(i,  bestModel[i] == l_False));
    }

    solver->setSolutionBasedPhaseSaving(false);
    lbool res = searchSATSolver(solver, modelAssumps);
    assert(res == l_True);
    solver->setSolutionBasedPhaseSaving(true);  
    checkModel(false, true);  
    if (ls_extend) {
      localsearch(bestModel);
    }
  //  logPrint("Debug: after extending, current UB: " + std::to_string(ubCost) + " size of best model " + std::to_string(bestModel.size()));
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

      solver->setSolutionBasedPhaseSaving(false);
      lbool res = searchSATSolver(solver, local_search_model);
      assert(res == l_True);
      solver->setSolutionBasedPhaseSaving(true);
      checkModel(true, true);
      }
      return true;
    }
    else {
      logPrint("Local search found no solution");
      return false;
    }
}

void CBLIN::minimizelinearsolution(vec<lbool> & sol) {
  if (ls_min) {
    if (!skip_local_search) {
       localsearch(sol);
    }
    return;
  }
  

  if (nbCurrentSoft == 0) {
    logPrint("No softs");
    return;
  }

  //TODO, check the model..... 
  vec<Lit> minimizable; 
  vec<bool> skip; 
  vec<Lit> fixed_assumptions; 
  time_t rec = time(NULL);
  for (int i = 0; i < isSoft.size(); i ++) {
    if (isSoft[i]) continue; 
    Lit l = mkLit(i, true); 
    if (literalTrueInModel(l, sol)) {
      fixed_assumptions.push(l);
    }
    else {
      fixed_assumptions.push(~l);
    }
  }

  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    uint64_t reducedWeight = maxsat_formula->getSoftClause(i).weight / maxsat_formula->getMaximumWeight();

    if (reducedWeight > 0) { 
            Lit l = maxsat_formula->getSoftClause(i).assumption_var; 
            assert (l != lit_Undef);
            assert(var(l) >= isSoft.size() || isSoft[var(l)] );
              if (literalTrueInModel(l, sol)) {
                // induces cost
                minimizable.push(l);
                skip.push(false);
              }
              else {
                fixed_assumptions.push(~l);
              }
      }
    }




  vec<Lit> assumps; 
  lbool res = l_False; 

  int skipped = 0;

  for (int i = 0; i < minimizable.size() ; i ++) {
    if (skip[i]) continue;
    Lit l = minimizable[i];
    assumps.clear();
    fixed_assumptions.copyTo(assumps);
    assumps.push(~l);
    res = searchSATSolver(solver, assumps);
    if (res == l_True) {
      fixed_assumptions.push(~l);
      
      for (int j = i+1; j < minimizable.size(); j++) {
        if (skip[j]) continue;
        Lit n = minimizable[j];
        if (literalTrueInModel(~n, solver->model)) {
          skip[j] = true; 
          fixed_assumptions.push(~n);
          skipped++;
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
    res = searchSATSolver(solver, fixed_assumptions);
    assert(res == l_True);
  }
  checkModel(false, true);
  time_t done = time(NULL);
  logPrint("minimization time " , done - rec, " init minsize " , minimizable.size(), " skipped ", skipped);

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
  logPrint("use_DPW=" , use_DPW);
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
  logPrint("ls_extend=", ls_extend);

  logPrint("Before search: UB ", ubCost, " LB ", lbCost, " off_set ", off_set, 
            " standarddization_removed ", standardization_removed, " preprocessing_removed ", 
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
  timeLimitCores += (time(NULL) - time_start);

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
      timeLimitCores = -1;
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

Solver * CBLIN::resetSolver() {
    logPrint("deleting solver");
    delete solver; 
    solver = newSATSolver();
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
    //TODO: think about how to do this better with the DPW, new variables are mostl ikely 1. 
    logPrint("save phase");
    solver->_user_phase_saving.clear();
		for (int i = 0; i < bestModel.size(); i++){
			solver->_user_phase_saving.push(bestModel[i]);		
		}
 }

 void CBLIN::softsSatisfied() {
     for (int i = 0; i < maxsat_formula->nSoft(); i++) {
        assert(maxsat_formula->getSoftClause(i).clause.size() == 1 );
        Lit l =  maxsat_formula->getSoftClause(i).clause[0];
        solver->setPolarity(var(l), sign(l) ? true : false);
    }
 }

//TODO parametrize on the model... 
 bool CBLIN::checkModel(bool from_local_search, bool improve_better) {
   logPrint("checkingModel size_of_model " , solver->model.size());

   auto lambda = [this](Lit l){ return literalTrueInModel(l, solver->model); };
   uint64_t modelCost = computeCostOfModel(&lambda);
   bool isBetter = modelCost < ubCost;
   if (isBetter) {
        ubCost = modelCost;
        time_best_solution = time(NULL);
        printProgress();
        saveModel(solver->model);
        bestModel.clear();
        solver->model.copyTo(bestModel);
        printBound(ubCost);
        checkGap();
        skip_local_search = from_local_search;
    }
  else if (improve_better && (modelCost == ubCost) && solver->nVars() >= bestModel.size()) {
      logPrint("Found same cost model covering more or the same amount of variables");
      saveModel(solver->model);
      bestModel.clear();
      solver->model.copyTo(bestModel);
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
