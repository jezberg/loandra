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

#include "MaxSAT.h"

#include <sstream>

using namespace openwbo;

/************************************************************************************************
 //
 // Public methods
 //
 ************************************************************************************************/

StatusCode MaxSAT::search() {
  if(print) printf("Error: Invalid MaxSAT algoritm.\n");
  throw MaxSATException(__FILE__, __LINE__, "Did not implement MaxSAT search");
  return _ERROR_;
}

void MaxSAT::setInitialTime(double initial) {
  initialTime = initial;
} // Sets the initial time.



/************************************************************************************************
 //
 // Utils for model management
 //
 ************************************************************************************************/

/*_________________________________________________________________________________________________
  |
  |  saveModel : (currentModel : vec<lbool>&)  ->  [void]
  |
  |  Description:
  |
  |    Saves the current model found by the SAT solver.
  |
  |  Pre-conditions:
  |    * Assumes that 'nbInitialVariables' has been initialized.
  |    * Assumes that 'currentModel' is not empty.
  |
  |  Post-conditions:
  |    * 'model' is updated to the current model.
  |
  |________________________________________________________________________________________________@*/
void MaxSAT::saveModel(vec<lbool> &currentModel) {
  assert(maxsat_formula->nInitialVars() != 0);
  assert(currentModel.size() != 0);

  model.clear();
  // Only store the value of the variables that belong to the
  // original MaxSAT formula.
  for (int i = 0; i < maxsat_formula->nInitialVars(); i++)
    model.push(currentModel[i]);
}

/*_________________________________________________________________________________________________
  |
  |  computeCostOriginalClauses : (reconstructed_model : vec<lbool>&)
  |
  |  Description:
  |
  |    Computes the cost of 'reconstructed_model'. The cost of a model is the sum of
  |    the weights of the unsatisfied soft clauses.
  |    
  |    if preprocessing is employed, the model should be reconstructed before this mehtod is invoked. 
  |
  |  Pre-conditions:
  |    * Assumes that 'currentModel' is not empty.
  |
  |________________________________________________________________________________________________@*/
uint64_t MaxSAT::computeCostOriginalClauses_legacy(vec<lbool> &reconstructed_model) {
  assert(reconstructed_model.size() != 0 || full_original_scla->nSoft() == 0);
  uint64_t currentCost = 0;

  for (int i = 0; i < full_original_scla->nSoft(); i++) {
    bool unsatisfied = true;
    for (int j = 0; j < full_original_scla->getSoftClause(i).clause.size(); j++) {
      Lit l = full_original_scla->getSoftClause(i).clause[j]; 

      assert(var(l) < reconstructed_model.size());
      if (literalTrueInModel(l, reconstructed_model)) {
        unsatisfied = false;
        break;
      }
    }

    if (unsatisfied) {
      currentCost += full_original_scla->getSoftClause(i).weight;
    }
  }
  return currentCost;
}

uint64_t MaxSAT::computeCostOriginalClauses(CaDiCaL::Solver* solver) {
  uint64_t currentCost = 0;
  assert(solver->vars() >= full_original_scla->nVars());

  for (int i = 0; i < full_original_scla->nSoft(); i++) {
    bool unsatisfied = true;
    for (int j = 0; j < full_original_scla->getSoftClause(i).clause.size(); j++) {
      Lit l = full_original_scla->getSoftClause(i).clause[j]; 
      int sat = solver->val(lit2Int(l));
      if (sat > 0) {
        unsatisfied = false;
        break;
      }
    }

    if (unsatisfied) {
      currentCost += full_original_scla->getSoftClause(i).weight;
    }
  }
  return currentCost;
}

uint64_t MaxSAT::computeCostObjective(CaDiCaL::Solver* solver) {
  uint64_t currentCost = standardization_removed; //this is 0 if preprocessing

  for (int i = 0; i < original_labels->nSoft(); i++) {
    assert(original_labels->getSoftClause(i).clause.size() == 1); 
    Lit l = original_labels->getSoftClause(i).clause[0];
    int res = solver->val(lit2Int(l));

    //// satisfied soft
    if (res > 0) {
      continue;
    }
    currentCost += original_labels->getSoftClause(i).weight;  
  } 
  if (do_preprocess) {
    currentCost += cost_removed_preprocessing;
  }
  return currentCost;
}

uint64_t MaxSAT::computeCostObjective_legacy(vec<lbool> &model) {
  assert(model.size() != 0);
  uint64_t currentCost = standardization_removed; //this is 0 if no preprocessing

  for (int i = 0; i < original_labels->nSoft(); i++) {
    assert(original_labels->getSoftClause(i).clause.size() == 1); 
    Lit l = original_labels->getSoftClause(i).clause[0];
    
    //// satisfied soft
    if (literalTrueInModel(l, model)) {
      continue;
    }
    currentCost += original_labels->getSoftClause(i).weight;  
  } 
  if (do_preprocess) {
    currentCost += cost_removed_preprocessing;
  }
  return currentCost;
}


/*_________________________________________________________________________________________________
  |
  |  isBMO : (cache : bool)  ->  [void]
  |
  |  Description:
  |
  |    Tests if the MaxSAT formula has lexicographical optimization criterion.
  |
  |  For further details see:
  |    * Joao Marques-Silva, Josep Argelich, Ana Graça, Inês Lynce: Boolean
  |      lexicographic optimization: algorithms & applications. Ann. Math.
  |      Artif. Intell. 62(3-4): 317-343 (2011)
  |
  |  Post-conditions:
  |    * 'orderWeights' is updated with the weights in lexicographical order if
  |      'cache' is true.
  |
  |________________________________________________________________________________________________@*/
bool MaxSAT::isBMO(bool cache) {
  assert(orderWeights.size() == 0);
  bool bmo = true;
  std::set<uint64_t> partitionWeights;
  std::map<uint64_t, uint64_t> nbPartitionWeights;

  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    partitionWeights.insert(maxsat_formula->getSoftClause(i).weight);
    nbPartitionWeights[maxsat_formula->getSoftClause(i).weight]++;
  }

  for (std::set<uint64_t>::iterator iter = partitionWeights.begin();
       iter != partitionWeights.end(); ++iter) {
    orderWeights.push_back(*iter);
  }

  std::sort(orderWeights.begin(), orderWeights.end(), greaterThan);

  uint64_t totalWeights = 0;
  for (int i = 0; i < (int)orderWeights.size(); i++)
    totalWeights += orderWeights[i] * nbPartitionWeights[orderWeights[i]];

  for (int i = 0; i < (int)orderWeights.size(); i++) {
    totalWeights -= orderWeights[i] * nbPartitionWeights[orderWeights[i]];
    if (orderWeights[i] < totalWeights) {
      bmo = false;
      break;
    }
  }

  if (!cache)
    orderWeights.clear();

  return bmo;
}

/************************************************************************************************
 //
 // Utils for printing
 //
 ************************************************************************************************/

// Prints information regarding the AMO encoding.
void MaxSAT::print_AMO_configuration(int encoding) {
  switch (encoding) {
  case _AMO_LADDER_:
    printf("c |  AMO Encoding:         %12s                      "
           "                                             |\n",
           "Ladder");
    break;

  default:
    printf("c Error: Invalid AMO encoding.\n");
    printf("s UNKNOWN\n");
    break;
  }
}

// Prints information regarding the PB encoding.
void MaxSAT::print_PB_configuration(int encoding) {
  switch (encoding) {
  case _PB_SWC_:
    printf("c |  PB Encoding:         %13s                        "
           "                                           |\n",
           "SWC");
    break;

  case _PB_GTE_:
    printf("c |  PB Encoding:         %13s                        "
           "                                           |\n",
           "GTE");
    break;

  default:
    printf("c Error: Invalid PB encoding.\n");
    printf("s UNKNOWN\n");
    break;
  }
}

// Prints information regarding the cardinality encoding.
void MaxSAT::print_Card_configuration(int encoding) {
  switch (encoding) {
  case _CARD_CNETWORKS_:
    printf("c |  Cardinality Encoding: %12s                                "
           "                                   |\n",
           "CNetworks");
    break;

  case _CARD_TOTALIZER_:
    printf("c |  Cardinality Encoding: %12s                                "
           "                                   |\n",
           "Totalizer");
    break;

  case _CARD_MTOTALIZER_:
    printf("c |  Cardinality Encoding:    %19s                             "
           "                            |\n",
           "Modulo Totalizer");
    break;

  default:
    printf("c Error: Invalid cardinality encoding.\n");
    printf("s UNKNOWN\n");
    break;
  }
}

void MaxSAT::blockModel(Solver *solver) {
  assert(model.size() != 0);

  vec<Lit> blocking;

  printf("v ");
  for (int i = 0; i < model.size(); i++) {
    indexMap::const_iterator iter = maxsat_formula->getIndexToName().find(i);
    if (iter != maxsat_formula->getIndexToName().end()) {
      if (model[i] == l_False)
        printf("-");
      printf("%s ", iter->second.c_str());
    }
  }
  printf("\n");

  for (int i = 0; i < model.size(); i++) {
    blocking.push((model[i] == l_True) ? ~mkLit(i) : mkLit(i));
  }

  solver->addClause(blocking);
}

void MaxSAT::logPrint(std::string s) {
  if (verbosity > 0) {
    std::cout << "c " << s << std::endl;
  }
}

void MaxSAT::printBound(uint64_t bound)
{
  if(!print) return;

  // print bound only, if its below the hard weight
  // FIXME: possible issue for PB instances when bound is negative; in MaxSAT bound is always positive
  if( bound <= maxsat_formula->getHardWeight() ) printf("o %" PRIu64 "\n", bound);
}

// Prints the best satisfying model if it found one.
void MaxSAT::printModel() {
  
  if (model.size() == 0)
      return;

  assert(maxsat_formula->getFormat() != _FORMAT_PB_); //disabled for now

  std::stringstream s;
  s << "v ";
  if (do_preprocess) {
    assert(model_of_original.size() > 0);
    for(int i = 0; i < model_of_original.size(); i++){
      //for (int i = 0; i < model.size(); i++) {
      if (model_of_original[i] == l_True) {
        s << "1";
      }
      else {
        s << "0";
      }
    }
  }
  else {
    for(int i = 0; i < maxsat_formula->nInitialVars(); i++){
      if (model[i] == l_True) {
        s << "1";
      }
      else {
        s << "0";
      }
    }
  }
  s << std::endl;
  std::cout << s.str(); 
}

std::string MaxSAT::printSoftClause(int id){
  assert (maxsat_formula->getFormat() == _FORMAT_MAXSAT_);
  assert (id < maxsat_formula->nSoft());

  std::stringstream ss;
  ss << maxsat_formula->getSoftClause(id).weight << " ";

  for (int j = 0; j < maxsat_formula->getSoftClause(id).clause.size(); j++){
    if (sign(maxsat_formula->getSoftClause(id).clause[j]))
      ss << "-";
    ss << (var(maxsat_formula->getSoftClause(id).clause[j])+1) << " ";
  }
  ss << "0\n";
  return ss.str();
}

void MaxSAT::printUnsatisfiedSoftClauses() {
  assert (model.size() != 0);

  std::stringstream s;
  int soft_size = 0;
  
  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    bool unsatisfied = true;
    for (int j = 0; j < maxsat_formula->getSoftClause(i).clause.size(); j++) {

      assert(var(maxsat_formula->getSoftClause(i).clause[j]) <
             model.size());
      if ((sign(maxsat_formula->getSoftClause(i).clause[j]) &&
           model[var(maxsat_formula->getSoftClause(i).clause[j])] ==
               l_False) ||
          (!sign(maxsat_formula->getSoftClause(i).clause[j]) &&
           model[var(maxsat_formula->getSoftClause(i).clause[j])] ==
               l_True)) {
        unsatisfied = false;
        break;
      }
    }

    if (unsatisfied) {
      s << printSoftClause(i);
      soft_size++;
    }
  }
  FILE * file = fopen (getPrintSoftFilename(),"w");
  fprintf(file,"p cnf %d %d\n",maxsat_formula->nInitialVars(),soft_size);
  fprintf(file,"%s", s.str().c_str());
}

// Prints search statistics.
void MaxSAT::printStats() {
  double totalTime = cpuTime();
  float avgCoreSize = 0;
  if (nbCores != 0)
    avgCoreSize = (float)sumSizeCores / nbCores;

  printf("c\n");
  if (model.size() == 0)
    printf("c  Best solution:          %12s\n", "-");
  else
    printf("c  Best solution:          %12" PRIu64 "\n", ubCost);
  printf("c  Total time:             %12.2f s\n", totalTime - initialTime);
  printf("c  Nb SAT calls:           %12d\n", nbSatisfiable);
  printf("c  Nb UNSAT calls:         %12d\n", nbCores);
  printf("c  Average core size:      %12.2f\n", avgCoreSize);
  printf("c  Nb symmetry clauses:    %12d\n", nbSymmetryClauses);
  printf("c\n");
}

// Prints the corresponding answer.
void MaxSAT::printAnswer(int type) {
  if (verbosity > 0 && print)
    printStats();

  if (type == _UNKNOWN_ && model.size() > 0)
    type = _SATISFIABLE_;
  
  if (type == _UNSATISFIABLE_) {
     printf("s UNSATISFIABLE\n");
     return;
  }

  if (do_preprocess) {
    assert(maxsat_formula->nHard() == 0 || model.size() > 0);
    model_of_original.clear(); 

    reconstruct_model_prepro(model, model_of_original);
    uint64_t newCost = computeCostOriginalClauses_legacy(model_of_original);

    if (newCost < ubCost) {
      logPrint("cost improved after reconstruction");
    }
    if (newCost > ubCost) {
      logPrint("warning: cost worsened after final reconstruct");
    }
    ubCost = newCost; 
  }

  printBound(ubCost);

  // store type in member variable
  searchStatus = (StatusCode)type;
  if(!print) return;

  switch (type) {
  case _SATISFIABLE_:
    printf("s SATISFIABLE\n");
    if (print_model)
      printModel();
    if (print_soft)
      printUnsatisfiedSoftClauses();
    break;
  case _OPTIMUM_:
    printf("s OPTIMUM FOUND\n");
    if (print_model)
      printModel();
    if (print_soft)
      printUnsatisfiedSoftClauses();
    break;
  case _UNKNOWN_:
    printf("s UNKNOWN\n");
    break;
  default:
    printf("c Error: Invalid answer type.\n");
  }
}

void MaxSAT::set_preprocessing_parameters
    (double pre_timeLimit, std::string pre_techs, bool gate_extraction_, bool label_matching_,
    int skip_technique_){
    //These could be supported by first encoding them into CNF. 
    assert(maxsat_formula->nCard() == 0);
    assert(maxsat_formula->nPB() == 0);

    do_preprocess = true;
    assert(pif == NULL);
    
    preprocess_time_limit = pre_timeLimit; 
	  prepro_verb = 0;
    prepro_techs = pre_techs;

    gate_extraction = gate_extraction_;
    label_matching = label_matching_;
    skip_technique = skip_technique_;

    logPrint("Preprocessing parameters:"); 
    logPrint("preprocess_time_limit: " + std::to_string(pre_timeLimit) );
    logPrint("prepro_verb: " + std::to_string(prepro_verb) );
    logPrint("prepro_techs: " + pre_techs );
    logPrint("gate_extraction: " + std::to_string(gate_extraction) );
    logPrint("label_matching: " + std::to_string(label_matching) );
    logPrint("skip_technique: " + std::to_string(skip_technique) );
  
}

void MaxSAT::setup_formula() {
  MaxSATFormula* temp_formula = NULL;
  full_original_scla = maxsat_formula->copySoftsFromFormula();

  if (do_preprocess) {
    temp_formula = preprocessed_formula();  
  }
  else {   
    temp_formula = standardized_formula();
  }
  delete maxsat_formula;
  maxsat_formula = temp_formula; 

  original_labels = maxsat_formula->copySoftsFromFormula();

}

MaxSATFormula* MaxSAT::standardized_formula() {
  MaxSATFormula *copymx = new MaxSATFormula();
  copymx->setInitialVars(maxsat_formula->nVars());

  for (int i = 0; i < maxsat_formula->nVars(); i++)
    copymx->newVar();

  for (int i = 0; i < maxsat_formula->nHard(); i++)
    copymx->addHardClause(maxsat_formula->getHardClause(i).clause);

  vec<Lit> clause; 
  std::map<int, uint64_t> existing_units;
  for (int i = 0; i < maxsat_formula->nSoft(); i++) {
    clause.clear();
    maxsat_formula->getSoftClause(i).clause.copyTo(clause);
    if (clause.size() != 1){
      Lit l = copymx->newLiteral();
      clause.push(l);
      copymx->addHardClause(clause);
    
      clause.clear();
      clause.push(~l);
      copymx->addSoftClause(maxsat_formula->getSoftClause(i).weight, clause);
    }  
    else {
      Lit l = clause[0];
      uint64_t w_mew = maxsat_formula->getSoftClause(i).weight;
      if (existing_units.find(lit2Int(l)) == existing_units.end()) {
        existing_units[lit2Int(l)] = w_mew;
      }
      else {
        existing_units[lit2Int(l)] += w_mew;
      }
    }
  }
  std::set<int> to_be_removed;
  for (auto iter = existing_units.begin(); iter != existing_units.end(); iter++ ) {
    int lit = iter->first;
    int negation = lit * (-1);
    uint64_t weight = iter->second;
    if (to_be_removed.find(lit) != to_be_removed.end() ) {
      continue;
    }
    bool contradicting_literal = existing_units.find(negation) != existing_units.end();
    clause.clear();
    if (contradicting_literal) {
      to_be_removed.insert(negation);
      uint64_t contr_w = existing_units[negation];
      if (contr_w > weight) {
        clause.push(int2Lit(negation));
        copymx->addSoftClause(contr_w - weight, clause);
        lbCost += weight; 
        standardization_removed += weight;
      }
      else if (contr_w < weight) {
        clause.push(int2Lit(lit));
        copymx->addSoftClause(weight - contr_w, clause);
        lbCost += contr_w;
        standardization_removed += contr_w;
      }
      else { // contradicting labels have equal weight
        lbCost += contr_w;
        standardization_removed += contr_w;
        continue;
      }
    }
    else{
      clause.push(int2Lit(lit));
      copymx->addSoftClause(weight, clause);
    }

  }
       


  //TODO, add the literals in the map... map is stored "key is hte exact polarity it appears in the soft clause"

  copymx->setProblemType(maxsat_formula->getProblemType());
  copymx->updateSumWeights(maxsat_formula->getSumWeights());
  copymx->setMaximumWeight(maxsat_formula->getMaximumWeight());
  copymx->setHardWeight(maxsat_formula->getHardWeight());
  return copymx;
}


MaxSATFormula* MaxSAT::preprocessed_formula() {

    int cla_before = maxsat_formula->nSoft() + maxsat_formula->nHard();
    int softs_before = maxsat_formula->nSoft(); 

    std::vector<std::vector<int> > clauses_out;
		std::vector<uint64_t> weights_out;
    std::vector<int> ppClause; 

    uint64_t top_orig = maxsat_formula->getSumWeights() + 1;

    for (int i = 0; i < maxsat_formula->nHard(); i++) {
        solClause2ppClause(maxsat_formula->getHardClause(i).clause, ppClause);
        assert(maxsat_formula->getHardClause(i).clause.size() == ppClause.size());
        clauses_out.push_back(ppClause);
        weights_out.push_back(top_orig);
    }

    for (int i = 0; i < maxsat_formula->nSoft(); i++) {
        bool isHardened = maxsat_formula->getSoftClause(i).weight == 0;
        if (!isHardened) {
          solClause2ppClause(maxsat_formula->getSoftClause(i).clause, ppClause);
          assert(maxsat_formula->getSoftClause(i).clause.size() == ppClause.size());
          clauses_out.push_back(ppClause);
          weights_out.push_back(maxsat_formula->getSoftClause(i).weight);
        } 
    }

    pif = new maxPreprocessor::PreprocessorInterface (clauses_out, weights_out, top_orig, false);
		
    pif->setBVEGateExtraction(gate_extraction);	
	  pif->setLabelMatching(label_matching);
	  pif->setSkipTechnique(skip_technique);

    pif->preprocess(prepro_techs, prepro_verb, preprocess_time_limit);
    ub_prepro = pif->getUpperBound();
    uint64_t lb = pif->getRemovedWeight();
    

    cost_removed_preprocessing += pif->getRemovedWeight();
    lbCost += cost_removed_preprocessing;

    //COLLECT NEW
    std::vector<std::vector<int> > pre_Clauses; 
		std::vector<uint64_t> pre_Weights; 
		std::vector<int> pre_Labels; //will not be used	
		pif->getInstance(pre_Clauses, pre_Weights, pre_Labels);
		uint64_t top = pif->getTopWeight();

    MaxSATFormula *copymx = new MaxSATFormula();
    copymx->setProblemType(maxsat_formula->getProblemType());
    copymx->setHardWeight(top);

    uint64_t init_vars = 0;
    uint64_t sum_of_weights = 0;
    uint64_t max_weight = 0;
    
    std::set<int> forced_hards;

    assert(pre_Weights.size() == pre_Clauses.size());
    for (int i = 0; i < pre_Weights.size(); i++) {
        uint64_t cur = pre_Weights[i];
        if (cur < top) {
          if (cur > max_weight) {
            max_weight = cur;
          }
          sum_of_weights += cur;
        }
        else {
          if(pre_Clauses[i].size() == 1) {
            forced_hards.insert(abs(pre_Clauses[i][0]));
          }
        }
        for (int j = 0; j < pre_Clauses[i].size(); j++) {
          int var = pre_Clauses[i][j];
          
          if (abs(var) > init_vars) {
            init_vars = abs(var);
          }
        }

    }

    copymx->setInitialVars(init_vars);
    copymx->updateSumWeights(sum_of_weights);
    copymx->setMaximumWeight(max_weight);
    
    for (int i = 0; i < init_vars; i++) {
      copymx->newVar();
    }

    vec<Lit> sol_cla;		
    int num_skipped = 0;
		for (int i = 0; i < pre_Clauses.size(); i++) {
			sol_cla.clear();				
			ppClause2SolClause(sol_cla, pre_Clauses[i]);
			assert(sol_cla.size() == pre_Clauses[i].size());
						
			uint64_t weight = pre_Weights[i];
			if (weight < top) {
				//SOFT 
				assert(sol_cla.size() == 1);
				assert(weight > 0);

        if (weight == cost_removed_preprocessing && forced_hards.find(abs(pre_Clauses[i][0])) != forced_hards.end() && num_skipped == 0) {
           num_skipped++; 
        }
        else{
          copymx->addSoftClause(weight, sol_cla);
        }
			}
			else {
				copymx->addHardClause(sol_cla);
			}			
		}
    //logPrint("Preprocess removed weight: "  + std::to_string(cost_removed_preprocessing) + " num_skipped " + std::to_string(num_skipped))  ;
    assert(cost_removed_preprocessing == 0 || num_skipped == 1);

      //logPrint("Preprocess time: " + print_timeSinceStart() + " removed weight: "  + std::to_string(cost_removed_preprocessing)) ;
    logPrint("Preprocessing left " + std::to_string(copymx->nHard()) + " clauses and " + std::to_string(copymx->nSoft()) + " softs");
    logPrint("Preprocessing removed " + std::to_string(cla_before - copymx->nHard()) + " clauses and " + std::to_string(softs_before - copymx->nSoft()) + " softs");
    logPrint("Preprocessing obtained ub: " + std::to_string(ub_prepro) + " lb " +std::to_string(lb));
  return copymx;
}

void MaxSAT::solClause2ppClause(const vec<Lit>  & solClause,  std::vector<int> & ppClause_out) {
	ppClause_out.clear();
	for (int i = 0; i < solClause.size(); i++) {
    Lit l = solClause[i];
    assert( int2Lit ( lit2Int( l ) ) == l ); 
		ppClause_out.push_back( lit2Int( l ));
	}
}

void MaxSAT::ppClause2SolClause(vec<Lit>  & solClause_out, const std::vector<int> & ppClause) {
	solClause_out.clear();
	for (int i = 0; i < ppClause.size(); i++) {
    int int_var = ppClause[i];
    assert( int_var == lit2Int ( int2Lit(int_var) )  ) ;
		solClause_out.push( int2Lit( int_var ));
	}
}

int MaxSAT::lit2Int(Lit l) {
	if (sign(l)) {
		return  -(var(l) + 1);
	}
	else {
		return var(l) + 1; 
	}
}

Lit MaxSAT::int2Lit(int l) {
	int var = abs(l) - 1;
	bool sign = l > 0;
	return sign ? mkLit(var) : ~mkLit(var);
}

void MaxSAT::reconstruct_model_prepro(vec<lbool> &currentModel, vec<lbool> &reconstructed_out) {
    time_t rec = time(NULL);
    std::vector<int> trueLiterals;
    for (int i = 0 ; i < currentModel.size() ; i++) {
      Lit l = mkLit(i, true); 
      if (literalTrueInModel(l, currentModel)) {
          trueLiterals.push_back(lit2Int(l));
      }
      else {
        assert(literalTrueInModel(~l, currentModel));
        trueLiterals.push_back(lit2Int(~l));
      }
    }
    std::vector<int> true_model = pif->reconstruct(trueLiterals);

    assert(true_model.size() == full_original_scla->nVars());
    for (int i = 0; i < true_model.size() ; i ++) {
      if (true_model[i] > 0) 
        reconstructed_out.push(l_True);
      else {
        reconstructed_out.push(l_False);
      }
    }
    time_t done = time(NULL);
    logPrint("Reconstruction time: " + std::to_string( done - rec ));
}

bool MaxSAT::literalTrueInModel(Lit l, vec<lbool> &model) {

  if (var(l) >= model.size()) {
    logPrint("Error, asking for truthness of literal beyond model size");
    assert(var(l) < model.size());
  }

  if (model[var(l)] == l_Undef) {
    logPrint("Warning, undef literal treated as false");
    return false;
  }

  return (sign(l) && model[var(l)] == l_False) || (!sign(l) && model[var(l)] == l_True);
}

void MaxSAT::print_statistics() {
   printf("c |                                                                "
           "                                       |\n");
    printf("c ========================================[ Problem Statistics "
           "]===========================================\n");
    printf("c |                                                                "
           "                                       |\n");

    if (maxsat_formula->getFormat() == _FORMAT_MAXSAT_)
      printf(
          "c |  Problem Format:  %17s                                         "
          "                          |\n",
          "MaxSAT");
    else
      printf(
          "c |  Problem Format:  %17s                                         "
          "                          |\n",
          "PB");

    if (maxsat_formula->getProblemType() == _UNWEIGHTED_)
      printf("c |  Problem Type:  %19s                                         "
             "                          |\n",
             "Unweighted");
    else
      printf("c |  Problem Type:  %19s                                         "
             "                          |\n",
             "Weighted");

    printf("c |  Number of variables:  %12d                                    "
           "                               |\n",
           maxsat_formula->nVars());
    printf("c |  Number of hard clauses:    %7d                                "
           "                                   |\n",
           maxsat_formula->nHard());
    printf("c |  Number of soft clauses:    %7d                                "
           "                                   |\n",
           maxsat_formula->nSoft());
    printf("c |  Number of cardinality:     %7d                                "
           "                                   |\n",
           maxsat_formula->nCard());
    printf("c |  Number of PB :             %7d                                "
           "                                   |\n",
           maxsat_formula->nPB());
}

std::string MaxSAT::core_2_string(vec<Lit> & core) {
  std::string s = "";
  for (int i = 0; i < core.size() ; i++) {
    s += " " + std::to_string(lit2Int(core[i]));
  }
  return s;
 }


std::string MaxSAT::model_2_string(vec<lbool> & model) {
  std::string s = "";
  for (int i = 0; i < model.size() ; i++) {
    if (model[i] == l_True) {
      s += "T";
    }
    else if (model[i] == l_True) {
      s += "F";
    }
    else {
      s += "U";
    }
    
  }
  return s;
}


