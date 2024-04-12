/*!
 * \author Jeremias Berg - jeremias.berg@helsinki.fi
 *
 * @section LICENSE
 *  
 * MiniSat,  Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
 *           Copyright (c) 2007-2010, Niklas Sorensson
 * Open-WBO, Copyright (c) 2013-2017, Ruben Martins, Vasco Manquinho, Ines Lynce
 * NuWLS -- Copyright (c) 2021-2022, Yi Chu, Xiang He
 * Loandra    Copyright (c) 2018-2024, Jeremias Berg, Emir Demirovic, Peter Stuckey
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

#include "utils/Options.h"
#include "utils/ParseUtils.h"
#include "utils/System.h"
#include <errno.h>
#include <signal.h>
#include <zlib.h>


#include <fstream>
#include <iostream>
#include <map>
#include <stdlib.h>
#include <string>
#include <vector>

#ifdef SIMP
#include "simp/SimpSolver.h"
#else
#include "core/Solver.h"
#endif

#include "MaxSAT.h"
#include "MaxTypes.h"
#include "ParserMaxSAT.h"
#include "ParserPB.h"

// Algorithms
#include "algorithms/Alg_OLL_ITER.h"
#include "algorithms/Alg_CBLIN.h"

#define VER1_(x) #x
#define VER_(x) VER1_(x)
#define SATVER VER_(SOLVERNAME)
#define VER VER_(VERSION)

using NSPACE::cpuTime;
using NSPACE::OutOfMemoryException;
using NSPACE::IntOption;
using NSPACE::BoolOption;
using NSPACE::StringOption;
using NSPACE::IntRange;
using NSPACE::parseOptions;
using namespace openwbo;

//=================================================================================================

static MaxSAT *mxsolver;

static void SIGINT_exit(int signum) {
  mxsolver->printAnswer(_UNKNOWN_);
  exit(_UNKNOWN_);
}



//=================================================================================================
// Main:


int main(int argc, char **argv) {
  printf("c\nc LOANDRA:\t an extension of Open-WBO to core-boosted linear search.\n");
  printf("c Version:\t July 2019 2018 -- Release: 1.3\n");
  printf("c Authors:\t Jeremias Berg\n");
  printf("c Contributors:\t Emir Demirovic, Peter Stuckey, Christoph Jabs, Marcus Leivo, Matti Järvisalo\n");
  printf("c We thank the developers of Open-WBO and NuWLS for their work\n");
  printf(
      "c\nc Open-WBO:\t a Modular MaxSAT Solver -- based on %s (%s version)\n",
      SATVER, VER);
  printf("c Version:\t September 2018 -- Release: 2.1\n");
  printf("c Authors:\t Ruben Martins, Vasco Manquinho, Ines Lynce\n");
  printf("c Contributors:\t Miguel Neves, Saurabh Joshi, Norbert Manthey, Mikolas Janota\n");
  printf("c Author of NuWLS:\t Yi Chu, Xiang He\n");
  char* _emergencyMemory = new char[16384];
  try {
    NSPACE::setUsageHelp("c USAGE: %s [options] <input-file>\n\n");

#if defined(__linux__)
    fpu_control_t oldcw, newcw;
    _FPU_GETCW(oldcw);
    newcw = (oldcw & ~_FPU_EXTENDED) | _FPU_DOUBLE;
    _FPU_SETCW(newcw);
    printf(
        "c WARNING: for repeatability, setting FPU to use double precision\n");
#endif

  BoolOption printmodel("Open-WBO", "print-model", "Print the solution found.\n", false);
    
  BoolOption oldformat("Open-WBO", "old-format", "Parse WCNF files in the pre 2022 format.\n", false);

  StringOption printsoft("Open-WBO", "print-unsat-soft", "Print unsatisfied soft claueses in the solution found.\n", NULL);

  IntOption verbosity("Open-WBO", "verbosity",
                        "Verbosity level (0=minimal, 1=more).\n", 0,
                        IntRange(0, 1));

  IntOption algorithm("Open-WBO", "algorithm",
                        "MaxSAT algorithm "
                        "(0=core-boosted linear search (default),2=oll)."
                        "\n",
                        0, IntRange(0, 1));

  IntOption cardinality("Encodings", "cardinality",
                          "Cardinality encoding (0=cardinality networks, "
                          "1=totalizer, 2=modulo totalizer).\n",
                          1, IntRange(0, 2));

  IntOption formula("Open-WBO", "formula",
                      "Type of formula (0=WCNF, 1=OPB).\n", 0, IntRange(0, 1));

  IntOption weight(
        "CBLIN", "weight-strategy",
        "Weight strategy (0=none, 1=weight-based, 2=diversity-based).\n", 2,
        IntRange(0, 2));

  IntOption pmreslin("CBLIN", "cb", "Run sat-unsat search in conjunction with core-guided search (i.e. core-boosted search): "
                                            "(0=not at all, 1=first cores then sat-unsat search 2=only sat-unsat) .\n", 1,
                  IntRange(0, 3));

   BoolOption pmreslin_delsol("CBLIN", "cb-del", "Reinitialise the SAT solver between core guided and linear phase.\n", true);
   BoolOption pmreslin_varres("CBLIN", "cb-varres", "Do varying resolution.\n", true);
   BoolOption pmreslin_relax2strat("CBLIN", "cb-r-2-s", "Relax cores before lowering the stratification bound.\n", false);
   BoolOption pmreslin_varresCG("CBLIN", "cb-varCG", "Do varying resolution during core-guided search.\n", false);
   BoolOption pmreslin_incvarres("CBLIN", "cb-i-varres", "Do varying resolution incrementally, without reinitialising the SAT solver.\n", false);
   IntOption pmreslin_cgLim("CBLIN", "cb-cglim", "Time limit for core guided phase (s): "
                                            "(-1=unlimited) .\n", 30,
                  IntRange(-1, INT_MAX));
   
  BoolOption pmreslin_dpw("CBLIN", "cb-DPW", "Use the dynamic polynomial watchdog (default=false in which case the generalized totalizer is used).\n", false);
  BoolOption pmreslin_dpw_coarse("CBLIN", "cb-DPW-coarse", "Only do coarse-convergence with the DPW for resolutions higher than 1.\n", false);
  BoolOption pmreslin_local_search("CBLIN", "cb-local-search", "Use NuWLS for solution minimization.\n", false);
  BoolOption extend("CBLIN", "extend-models", "Extend models to the variables in cardinality constraints.\n", true);



  BoolOption prepro_rec("PREPROCESS", "pr-rec", "Reconstruct solutions before computing their costs (only applicable when preprocessing).\n", false);
  BoolOption prepro_min("PREPROCESS", "pr-min", "Minimize solutions locally after preprocessing.\n", true);
  IntOption prepro_min_strat("PREPROCESS", "pr-min-strat", "Strategy for solution minimization: 1=agressive (all solutions), 2=only the two first in each resolution: "
                                            "(0=only the best after each resolution) .\n", 0,
                  IntRange(0, 2));
  StringOption prT("PREPROCESS", "pr-tech", "Preprocess techniques used (see MaxPRE documentation for more details).\n", "[u]#[uvsrgVGc]");
  BoolOption preprocess("PREPROCESS", "preprocess", "Preprocess the instance prior to search.\n", true);

  

    parseOptions(argc, argv, true);
    std::string preTechs(prT);

    double initial_time = cpuTime();
    MaxSAT *S = NULL;

    switch ((int)algorithm) {
    
    case _ALGORITHM_CBLIN_:
      S = new CBLIN(verbosity, weight, pmreslin, pmreslin_delsol, pmreslin_varres, pmreslin_varresCG, 
                    pmreslin_cgLim, pmreslin_relax2strat, pmreslin_incvarres, prepro_rec, 
                    prepro_min,prepro_min_strat, pmreslin_dpw, pmreslin_dpw_coarse, extend, pmreslin_local_search);
      break;
    
    case _ALGORITHM_OLLITER_:
      S = new OLL_ITER(verbosity, cardinality, prepro_rec);
      break;

    default:
      printf("c Error: Invalid MaxSAT algorithm.\n");
      printf("s UNKNOWN\n");
      exit(_ERROR_);
    }

    signal(SIGXCPU, SIGINT_exit);
    signal(SIGTERM, SIGINT_exit);

    if (argc == 1) {
      printf("c Error: no filename.\n");
      printf("s UNKNOWN\n");
      exit(_ERROR_);
    }


    gzFile in = (argc == 1) ? gzdopen(0, "rb") : gzopen(argv[1], "rb");
    if (in == NULL)
      printf("c ERROR! Could not open file: %s\n",
             argc == 1 ? "<stdin>" : argv[1]),
          printf("s UNKNOWN\n"), exit(_ERROR_);

    MaxSATFormula *maxsat_formula = new MaxSATFormula();

    if ((int)formula == _FORMAT_MAXSAT_) {
      parseMaxSATFormula(in, oldformat, maxsat_formula);
      maxsat_formula->setFormat(_FORMAT_MAXSAT_);
    } else {
      ParserPB *parser_pb = new ParserPB();
      parser_pb->parsePBFormula(argv[1], maxsat_formula);
      maxsat_formula->setFormat(_FORMAT_PB_);
    }
    gzclose(in);


  
   

    if (S->getMaxSATFormula() == NULL)
      S->loadFormula(maxsat_formula);
    
    printf("c Before Setup: \n");
    S->print_statistics();
    double parsed_time = cpuTime();
    printf("c |  Parse time:           %12.2f s                                "
            "                                 |\n",
            parsed_time - initial_time);
    printf("c |                                                                "
            "                                       |\n");
    
    
    S->setPrintModel(printmodel);
    S->setPrintSoft((const char *)printsoft);
    S->setInitialTime(initial_time);
    
    mxsolver = S;
    mxsolver->setPrint(true);
    if (preprocess) {
      //TODO abstract all of these into parameters: timelimit, techs, gate extraction, label matching, skiptechnique.
      mxsolver->set_preprocessing_parameters(30, preTechs, false, true, 20);
    }
    mxsolver->setup_formula();
    printf("c After Setup: \n");
    mxsolver->print_statistics(); 

    int ret = (int)mxsolver->search();
    delete mxsolver; // S
    return ret;
  } catch (OutOfMemoryException &) {
    sleep(1);
    delete[] _emergencyMemory;
    //TODO print the solution here.
    std::cout << "c Error: Out of memory." << std::endl;
    if (mxsolver->hasSolution()) {
      mxsolver->printAnswer(_UNKNOWN_);
      exit(_UNKNOWN_);
    }
    else {
      std::cout << "s UNKNOWN" << std::endl;
      exit(_ERROR_);
    }
  } catch(MaxSATException &e) {
    sleep(1);
    delete[] _emergencyMemory;
    //TODO print the solution
    std::cout << "c Error: MaxSAT Exception: %s" << std::endl;
    std::cout << e.getMsg() << std::endl;
    std::cout <<  "s UNKNOWN" << std::endl;
    exit(_ERROR_);
  }
}
