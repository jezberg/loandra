/*!
 * \author Jeremias Berg jeremiasberg@gmail.com, based on code by Ruben Martins - ruben@sat.inesc-id.pt
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

#include "Alg_OLL_ITER.h"

using namespace openwbo;

uint64_t OLL_ITER::findNextWeight(uint64_t weight) {

  uint64_t nextWeight = 1;

  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    if (maxsat_formula->getSoftClause(i).weight > nextWeight &&
        maxsat_formula->getSoftClause(i).weight < weight)
      nextWeight = maxsat_formula->getSoftClause(i).weight;
  }

  for (auto it = cardinality_assumptions.begin();
       it != cardinality_assumptions.end(); ++it) {
    assert(boundMapping.find(*it) != boundMapping.end());
    std::pair<std::pair<int, uint64_t>, uint64_t> soft_id = boundMapping[*it];
    if (soft_id.second > nextWeight && soft_id.second < weight)
      nextWeight = soft_id.second;
  }

  return nextWeight;
}

uint64_t OLL_ITER::findNextWeightDiversity(uint64_t weight) {

  assert(nbSatisfiable > 0); // Assumes that unsatSearch was done before.
  uint64_t nextWeight = weight;
  int nbClauses = 0;
  std::set<uint64_t> nbWeights;
  float alpha = 1.25;
  bool findNext = false;
  int safetyiteration = 0; 

  for (;;) {
    safetyiteration++;
    if (safetyiteration > 1000) {
      logPrint("Check weights on this instance ");
      break;
    }

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

    for (auto it = cardinality_assumptions.begin();
         it != cardinality_assumptions.end(); ++it) {
      assert(boundMapping.find(*it) != boundMapping.end());
      std::pair<std::pair<int, uint64_t>, uint64_t> soft_id = boundMapping[*it];
      if (soft_id.second >= nextWeight) {
        nbClauses++;
        nbWeights.insert(soft_id.second);
      }
    }

    if (nextWeight == 1 || (float)nbClauses / nbWeights.size() > alpha ||
        (unsigned)nbClauses >=
            (unsigned)nRealSoft() + cardinality_assumptions.size() )
      break;

    if (nbSatisfiable == 1 && !findNext)
      findNext = true;
  }
  return nextWeight;
}

StatusCode OLL_ITER::weighted() {
  lbool res = l_True;
  vec<Lit> assumptions;
  vec<Encoder *> soft_cardinality;

  assert(nbSatisfiable == 1);
  assert(solver->vars() == maxsat_formula->nVars());
  assert(maxsat_formula->getFormat() != _FORMAT_PB_);

  min_weight = findNextWeightDiversity(min_weight);
  
  for (;;) {

    if (lbCost > ubCost) {
      logPrint("LB > UB ERROR"); 
      return _ERROR_;
    }

    int not_considered = setAssumptions(assumptions);
    logPrint("Doing a SAT call at " + print_timeSinceStart());
    res = ICadical::searchSATSolver(solver, assumptions);

    if (res == l_True) {
      nbSatisfiable++;

      checkModel();
      
      
      if (lbCost == ubCost) {
        logPrint("LB == UB");
        logPrint("Terminating with " + std::to_string(nonreformulatedCores) + " unrelaxed out of " +  std::to_string(nbCores), 0);
        printBound(ubCost);
        logPrint("HARDENED: " + std::to_string(num_hardened_me));
        printAnswer(_OPTIMUM_);
        return _OPTIMUM_;
      }
     
      /* 
      if (lbCost == 803) {
        logPrint("Terminating with " + std::to_string(nonreformulatedCores) + " unrelaxed out of " +  std::to_string(nbCores), 0);
        printBound(ubCost);
        logPrint("HARDENED: " + std::to_string(num_hardened_me));
        printAnswer(_OPTIMUM_);
        return _OPTIMUM_;
      }
      */
      
     

      
      solver = hardenClauses(soft_cardinality);
      
      
      if (!reformulateDelayed(soft_cardinality) ) {
        if (not_considered == 0) {
          logPrint("WEIRDNESS IS GOING ON");
          logPrint("Cur W: " + std::to_string(min_weight));
          return _ERROR_;
        }
        else {
          logPrint("Updating weight: " + std::to_string(min_weight));
          min_weight = findNextWeightDiversity(min_weight);
        }
         //   logPrint("New weight: " + std::to_string(min_weight));
      }
     
    }
    else  {
      assert(res == l_False);
      vec<Lit> core;
      vec<Lit> conflict;
      ICadical::getCore(solver, conflict, assumptions);
      processCore(conflict, core);

      uint64_t min_core = computeCostCore(core);

      lbCost += min_core;
      time_best_lb = time(NULL);
      logPrint("LB : " + std::to_string(lbCost) +  " CS : " + std::to_string(core.size()) + " W: " + std::to_string(min_core)); 

      nbCores++;
      
      if (lbCost == ubCost) {
        assert(nbSatisfiable > 0);
        logPrint("LB == UB");
        logPrint("Terminating with " + std::to_string(nonreformulatedCores) + " unrelaxed out of " + std::to_string(nbCores), 0  );
        logPrint("HARDENED: " + std::to_string(num_hardened_me));
        printBound(ubCost);
        printAnswer(_OPTIMUM_);
        return _OPTIMUM_;
      }

      sumSizeCores += core.size();
      decreaseWeights(core, min_core, soft_cardinality);
      /*
      if ( ( ubCost - lbCost < maxw_nothardened ) ) {
        solver = hardenClauses(cardinality_assumptions, soft_cardinality);
      }
      */

      cores_found.push();
      new (&cores_found[cores_found.size() - 1]) vec<Lit>();
      core.copyTo( cores_found[cores_found.size() - 1]  );
      core_weights.push(min_core);
      nonreformulatedCores++;
      
      
    }
  }
}



StatusCode OLL_ITER::search() {

  if (encoding != _CARD_TOTALIZER_) {
    if(print) {
      printf("Error: Currently algorithm OLL with iterative encoding only "
             "supports the totalizer encoding.\n");
      printf("s UNKNOWN\n");
    }
    throw MaxSATException(__FILE__, __LINE__, "OLL only supports totalizer");
    return _UNKNOWN_;
  }

  printf("c WARNING: the OLL implementation has not been tested properly");

  printConfiguration();

  logPrint("OLL ALGORITHM ");
  logPrint("OLL PREPROCESS = " + std::to_string(do_preprocess));

  
  StatusCode set = setup();
  if (set == _UNSATISFIABLE_) {
    logPrint("No solutions to instance");
    return _UNSATISFIABLE_;
  }

  if (set == _OPTIMUM_) {
    assert(pif != NULL);
    logPrint("Solved by preprocessing");
    ubCost =  cost_removed_preprocessing; ;
    printBound(ubCost);
    printAnswer(_OPTIMUM_);
    return _OPTIMUM_;
  }
  return weighted();
}
/************************************************************************************************
 //
 // HELPERS
 //
 ************************************************************************************************/
  StatusCode OLL_ITER::setup() {
    time_start = time(NULL);
	  time_best_solution = time_start;
    time_best_lb = time_start;
    bestModel.clear();
  
    initAssumptions();
    solver = rebuildSolver();

    nOrigVars = solver->vars(); 

    min_weight = maxsat_formula->getMaximumWeight();
    maxw_nothardened = maxsat_formula->getMaximumWeight();
    logPrint("Init hard w: " + std::to_string(maxw_nothardened));

    vec<Lit> dummy;
    lbool res; 
    

    res = ICadical::searchSATSolver(solver, dummy);
    if (res != l_True) {
      return _UNSATISFIABLE_;
    }
    checkModel();

    nbSatisfiable++;

    if (maxsat_formula->nSoft() == 0) return _OPTIMUM_;
    else return _UNKNOWN_;
  }


  bool OLL_ITER::literalTrueInModel(Lit l, const vec<lbool> & curModel) {
    assert(var(l) < curModel.size());
    return ((sign(l) && curModel[var(l)] == l_False) || (!sign(l) && curModel[var(l)] == l_True));
  }

  //Returns number assumptions missed;
  int OLL_ITER::setAssumptions(vec<Lit> & assumptions_out) {
      assumptions_out.clear();
      int not_considered = 0;
      for (int i = 0; i < maxsat_formula->nSoft(); i++) {
        if (maxsat_formula->getSoftClause(i).weight >= min_weight) {
          assumptions_out.push(~maxsat_formula->getSoftClause(i).assumption_var);
        }
        else if (maxsat_formula->getSoftClause(i).weight > 0) {
          not_considered++;
        }
      }

      for (auto it = cardinality_assumptions.begin(); it != cardinality_assumptions.end(); ++it) {
        assert(boundMapping.find(*it) != boundMapping.end());
        std::pair<std::pair<int, uint64_t>, uint64_t> soft_id = boundMapping[*it];
        if (soft_id.second  >=   min_weight ) {
              assumptions_out.push(~(*it));
        }
        else if (soft_id.second > 0 ) {
          not_considered++;
        }
      }
      return not_considered; 
  }

uint64_t OLL_ITER::computeCostCore(const vec<Lit> &conflict) {

  assert(conflict.size() != 0);

  if (maxsat_formula->getProblemType() == _UNWEIGHTED_) {
    return 1;
  }

  uint64_t min_core = UINT64_MAX;
  for (int i = 0; i < conflict.size(); i++) {
    Lit p = conflict[i];
    if (coreMapping.find(p) != coreMapping.end()) {
      int indexCore = coreMapping[p];
      assert(maxsat_formula->getSoftClause(indexCore).weight > 0);
      assert(maxsat_formula->getSoftClause(indexCore).assumption_var == p);
      // NEW SOFT CLAUSE 
      if (maxsat_formula->getSoftClause(indexCore).weight < min_core) {
            min_core = maxsat_formula->getSoftClause(indexCore).weight;
      }
    }
      //EXISTING SOFT CLAUSE
    if (boundMapping.find(p) != boundMapping.end()) {
      std::pair<std::pair<int, uint64_t>, uint64_t> soft_id = boundMapping[ p ];
      if (soft_id.second < min_core) {
        min_core = soft_id.second;
      }
    }
  }
  return min_core;
}

// TODO if core monoimization can be added. 
void OLL_ITER::processCore(const vec<Lit> & orig_core, vec<Lit> &processed_core ) {
  orig_core.copyTo(processed_core);
}

void OLL_ITER::decreaseWeights(const vec<Lit> &core, const uint64_t min_core, vec<Encoder *> & soft_cardinality) {
     for (int i = 0; i < core.size(); i++) {
        Lit p = core[i];
        bool isSoftClause = coreMapping.find(p) != coreMapping.end();
        if (isSoftClause) {
            int indexSoft = coreMapping[p];
            assert(maxsat_formula->getSoftClause(indexSoft).weight - min_core >= 0);
            assert(maxsat_formula->getSoftClause(indexSoft).assumption_var == p);
          
            maxsat_formula->getSoftClause(indexSoft).weight -= min_core;
            if (maxsat_formula->getSoftClause(indexSoft).weight == 0) {
              num_hardened++;
              notHardened_ind.erase(p);

              maxsat_formula->getSoftClause(indexSoft).assumption_var = lit_Undef; //this field should never be touched again...
            }
        }
        else {
          assert(cardinality_assumptions.find(p) != cardinality_assumptions.end());
          assert(boundMapping.find(p) != boundMapping.end());
          
          //<<ID, bound>, weight>
          std::pair<std::pair<int, uint64_t>, uint64_t> soft_id = boundMapping[p];
          int64_t cur_bound = soft_id.first.second;
          uint64_t cur_w = soft_id.second;
          int soft_card_id = soft_id.first.first;

          assert(cur_w - min_core >= 0); 
          boundMapping[p] = std::make_pair(std::make_pair(soft_card_id, cur_bound), cur_w - min_core);
          
          extendCardEnc(p, soft_cardinality);

          if (cur_w - min_core == 0) { 
            num_card_dropped++;
            num_hardened++;
            notHardened_ind.erase(p);
          }
        }
     }
}

void OLL_ITER::reformulateCore(const vec<Lit> &core, const uint64_t min_core, vec<Encoder *> & soft_cardinality) {
  vec<Lit> newCore;
  int soft_relax = 0;
  int card_relax = 0; 
  for (int i = 0; i < core.size(); i++) {
        Lit p = core[i];
        bool isSoftClause = coreMapping.find(p) != coreMapping.end();
        if (isSoftClause) {
          newCore.push(p);
          soft_relax++;
        }
        else {
          assert(cardinality_assumptions.find(p) != cardinality_assumptions.end());
          increaseBound(p, min_core, soft_cardinality);
          
          newCore.push(p);
          card_relax++;
        }
  }

  if (soft_relax + card_relax > 1) {

    Encoder *e = new Encoder();
    e->setIncremental(_INCREMENTAL_ITERATIVE_);
    //Terrible hack need to integrate cadical to make it work
    e->buildCardinality( solver, newCore, 1);

    soft_cardinality.push(e);
    assert(e->outputs().size() > 1);

    Lit out = e->outputs()[1];
    boundMapping[out] = std::make_pair(std::make_pair(soft_cardinality.size() - 1, 1), min_core);
    cardinality_assumptions.insert(out);  
    notHardened_ind.insert(out);   
  }
  nonreformulatedCores--;
}

void OLL_ITER::findCardinality(Lit p, int64_t & cur_bound_out, uint64_t & bound_w_out, 
                              Encoder * & bound_enc_out, int & soft_card_id,  vec<Encoder *> & soft_cardinality) {
  assert(boundMapping.find(p) != boundMapping.end());
  //<<ID, bound>, weight>
  std::pair<std::pair<int, uint64_t>, uint64_t> soft_id = boundMapping[p];
  soft_card_id = soft_id.first.first;
  
  cur_bound_out = soft_id.first.second;
  bound_enc_out = soft_cardinality[soft_card_id];
  bound_w_out = soft_id.second;
}

void OLL_ITER::extendCardEnc(Lit p, vec<Encoder *> & soft_cardinality) {
  int64_t cur_bound = 0; 
  Encoder * cur = NULL;
  uint64_t dummy = 0; 
  int soft_card_id = 0;
  findCardinality(p, cur_bound, dummy, cur, soft_card_id, soft_cardinality);
  int64_t newBound = cur_bound + 1;

  vec<Lit> joinObjFunction;
  vec<Lit> encodingAssumptions;
  joinObjFunction.clear();
  encodingAssumptions.clear();

  cur->incUpdateCardinality(  solver, joinObjFunction, cur->lits(), newBound, encodingAssumptions);

} 


void OLL_ITER::increaseBound(Lit p, const uint64_t min_core,  vec<Encoder *> & soft_cardinality) {
  int64_t cur_bound = 0; 
  Encoder * cur = NULL;
  uint64_t dummy = 0; 
  int soft_card_id = 0;
  findCardinality(p, cur_bound, dummy, cur, soft_card_id, soft_cardinality);
  int64_t newBound = cur_bound + 1;

  // only enforce if new bound is less than total number of literals 
  if (newBound < (unsigned)cur->outputs().size()) {
      Lit out = cur->outputs()[newBound];
              
      auto it_c = cardinality_assumptions.find(out);
     
      bool new_assump = it_c == cardinality_assumptions.end();

      if (new_assump) {
        boundMapping[out] = std::make_pair(std::make_pair(soft_card_id, newBound), min_core);
        cardinality_assumptions.insert(out);
      }
        
       else {
        Encoder * dummy = NULL; 
        int64_t old_bound = 0; 
        uint64_t old_w = 0; 
        int old_id = 0;

        findCardinality(out, old_bound, old_w, dummy, old_id, soft_cardinality);
        assert(old_bound == newBound);
        boundMapping[out] = std::make_pair(std::make_pair(soft_card_id, newBound), min_core + old_w);
      }
      
  }          
}

 CaDiCaL::Solver *OLL_ITER::hardenClauses(vec<Encoder *> & soft_cardinality) { 
    uint64_t bound;

    bool use_new = false && do_preprocess && use_reconstruct;
    
    if (use_new) {
      bound = ubLabelCost - lbCost;
    }
    else {
      bound = ubCost - lbCost;
    }
    if ( bound >= maxw_nothardened  ) {
      return solver;
    }


    vec<lbool> & curModel = use_new ? hardeningModel : bestModel;
     logPrint("Hardening with gap: " + std::to_string(bound));
     

     int num_hardened_round = 0;
	   maxw_nothardened = 0;
     bool sat; 
    
    auto it = notHardened_ind.begin();
    while (it != notHardened_ind.end()) {
        auto current = it++;
        
        Lit card = *current;

        bool isSoft = coreMapping.find(card) != coreMapping.end();
        uint64_t cur_weight;


        Encoder * cur; 
        int soft_card_id = 0;
        int64_t cur_bound;


        if (isSoft) {
          cur_weight = maxsat_formula->getSoftClause(coreMapping[card]).weight;
        }
        else {
          findCardinality(card, cur_bound, cur_weight, cur, soft_card_id, soft_cardinality);
        }
        assert(cur_weight > 0);

        if (var(card) >= bestModel.size()) {
          sat = false; 
        }
        else {
          sat = literalTrueInModel(~card, curModel);
        }

        bool shouldHarden = cur_weight > bound   || ( (cur_weight == bound)  && sat) ;

        if (shouldHarden) {
          vec<Lit> clause;
				  clause.clear();

          clause.push(~card);
				  ICadical::addClause(solver, clause);
          
          if (isSoft) {
            int i = coreMapping[card];
            maxsat_formula->getSoftClause(i).weight = 0;
            maxsat_formula->getSoftClause(i).assumption_var = lit_Undef;
          }
          else {				    
            boundMapping[card] = std::make_pair(std::make_pair(soft_card_id, cur_bound), 0);
				    num_card_dropped++;
          }


          num_hardened++;
          num_hardened_me++;
				  num_hardened_round++;
          notHardened_ind.erase(current);

        }
        else if (cur_weight > maxw_nothardened) {
			  	maxw_nothardened = cur_weight;
			  }  
        
        
      
			  
    }

    logPrint("Hardened in total: " + std::to_string(num_hardened_round) + " clauses");
    logPrint("Hardening again at gap " + std::to_string(maxw_nothardened));
		return solver;		
   }


/************************************************************************************************
 //
 // Rebuild MaxSAT solver
 //
 ************************************************************************************************/

/*_________________________________________________________________________________________________
  |
  |  rebuildSolver : [void]  ->  [Solver *]
  |
  |  Description:
  |
  |    Rebuilds a SAT solver with the current MaxSAT formula.
  |
  |________________________________________________________________________________________________@*/
CaDiCaL::Solver *OLL_ITER::rebuildSolver() {

  CaDiCaL::Solver *S = ICadical::newSATSolver();
  S->reserve(maxsat_formula->nVars());

  for (int i = 0; i < maxsat_formula->nHard(); i++)
    ICadical::addClause(S, maxsat_formula->getHardClause(i).clause);

  assert(maxsat_formula->nPB() == 0);
  assert(maxsat_formula->nCard() == 0);
  
  return S;
}


/************************************************************************************************
 //
 // Other protected methods
 //
 ************************************************************************************************/

/*_________________________________________________________________________________________________
  |
  |  initRelaxation : [void] ->  [void]
  |
  |  Description:
  |
  |    Initializes the assumptions and the core mapping
  |
  | 
  |________________________________________________________________________________________________@*/
void OLL_ITER::initAssumptions() {

  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    assert(maxsat_formula->getSoftClause(i).clause.size() == 1);
    Lit l = maxsat_formula->getSoftClause(i).clause[0];
    maxsat_formula->getSoftClause(i).assumption_var = ~l;
    coreMapping[~l] = i;
    origWeights.push(maxsat_formula->getSoftClause(i).weight);
    notHardened_ind.insert(~l);
  }

}

void OLL_ITER::logPrint(std::string s, int verb_bound) {
  if (verbosity >= verb_bound) {
    std::cout << "c " << s << std::endl;
  }
}

void OLL_ITER::printProgress() {
  logPrint("Progress at: " + print_timeSinceStart(), 0);
  logPrint("UB: " + std::to_string(ubCost) + " at " + std::to_string(time_best_solution - time_start), 0);
  logPrint("LB: " + std::to_string(lbCost) + " at " + std::to_string(time_best_lb - time_start), 0 );
  time_t latestProgress = time_best_solution > time_best_lb ? time_best_solution: time_best_lb;
  logPrint("GAP: " + std::to_string(ubCost- lbCost) + " at " + std::to_string(latestProgress - time_start), 0 );
  logPrint("HARDENED: " + std::to_string(num_hardened_me));

}

  int OLL_ITER::nRealSoft() {
    return maxsat_formula->nSoft() - num_hardened;
  }  

std::string OLL_ITER::print_timeSinceStart() {
  return std::to_string(timeSinceStart());
}

time_t OLL_ITER::timeSinceStart() {
  time_t cur = time(NULL);
  return cur - time_start;
}

bool OLL_ITER::reformulateDelayed(vec<Encoder *> & soft_cardinality) {
  int allCores = nonreformulatedCores;
  assert(nonreformulatedCores ==  cores_found.size());
  assert(nonreformulatedCores ==  core_weights.size());
  

  while (cores_found.size()) {
    int cur_ind = cores_found.size() - 1; 
    reformulateCore(cores_found[cur_ind], core_weights[cur_ind], soft_cardinality);
    cores_found.pop(); 
    core_weights.pop();
  }
  logPrint("Relaxed: " + std::to_string(allCores) +  " cores") ;
  return allCores  > 0;
}

bool OLL_ITER::checkModel() {
   uint64_t modelCost;
   uint64_t clausecost;
   uint64_t labelCost;


   auto lambda = [ this](Lit l){return solver->val(MaxSAT::lit2Int(l) > 0);};
   labelCost = computeCostObjective(&lambda);
   if (do_preprocess && use_reconstruct) {
      vec<lbool> model;
      ICadical::getModel(solver, model);
      vec<lbool> reconstruct; 
      reconstruct_model_prepro(model, reconstruct);
      auto lambda_oirg = [this, &reconstruct](Lit l){return literalTrueInModel(l, reconstruct);};
      clausecost = computeCostOriginalClauses(&lambda_oirg);
   } 
   else {
      if (!do_preprocess) clausecost = computeCostOriginalClauses(&lambda);
      else clausecost = labelCost;
   }


   if (labelCost != clausecost) { 
      assert(clausecost <= labelCost);
      modelCost = clausecost;
   }
   else {
     modelCost = labelCost;
   }  

   if (do_preprocess && use_reconstruct) {
     bool better_hardening = labelCost < ubLabelCost;
     if (better_hardening) {
        logPrint("New hardening model");
        ubLabelCost = labelCost;
        hardeningModel.clear();
        model.copyTo(hardeningModel);
     }
   }

   bool isBetter = (modelCost < ubCost) || (nbSatisfiable == 1);
   if (isBetter) {
        ubCost = modelCost;
        time_best_solution = time(NULL);
        saveModel(model);
        bestModel.clear();
        model.copyTo(bestModel);
        printBound(ubCost);
        printProgress();
    }
    if (modelCost == ubCost && bestModel.size() == 0) {
      bestModel.clear();
      model.copyTo(bestModel);
      saveModel(model);
    }
    /*
    if (bestModel.size() < solver->nVars()) {
      logPrint("Extending model");
      extendModel();
      assert(bestModel.size() == solver->nVars());
    }
    */
    return isBetter;
 }

 void OLL_ITER::extendModel() {
   ///Extend 

    vec<Lit> dummy;
    lbool res; 
    
    for (int i = 0; i < nOrigVars ; i ++) {
      dummy.push(mkLit(i, bestModel[i] == l_False));
    }

    res = ICadical::searchSATSolver(solver, dummy);

    assert( res == l_True);
   
    bestModel.clear();
    vec<lbool> model;
    ICadical::getModel(solver, model);
    model.copyTo(bestModel);
    
 }  






  
