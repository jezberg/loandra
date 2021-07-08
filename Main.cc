/*!
 * \author Jeremias Berg - ruben@sat.inesc-id.pt
 *
 * @section LICENSE
 *  
 * MiniSat,  Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
 *           Copyright (c) 2007-2010, Niklas Sorensson
 * Open-WBO, Copyright (c) 2013-2017, Ruben Martins, Vasco Manquinho, Ines Lynce
 * Loandra    Copyright (c) 2018-2019, Jeremias Berg, Emir Demirovic, Peter Stuckey
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
#include "algorithms/Alg_LinearSU.h"
#include "algorithms/Alg_MSU3.h"
#include "algorithms/Alg_OLL.h"
#include "algorithms/Alg_OLL_ITER.h"
#include "algorithms/Alg_PartMSU3.h"
#include "algorithms/Alg_WBO.h"
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
  printf("c Version:\t July 2019 2018 -- Release: 1.0\n");
  printf("c Authors:\t Jeremias Berg, Emir Demirovic, Peter Stuckey\n");
  printf("c We thank the developers of Open-WBO for their work\n");
  printf(
      "c\nc Open-WBO:\t a Modular MaxSAT Solver -- based on %s (%s version)\n",
      SATVER, VER);
  printf("c Version:\t September 2018 -- Release: 2.1\n");
  printf("c Authors:\t Ruben Martins, Vasco Manquinho, Ines Lynce\n");
  printf("c Contributors:\t Miguel Neves, Saurabh Joshi, Norbert Manthey, Mikolas Janota\n");
  printf("c Contact:\t open-wbo@sat.inesc-id.pt -- "
         "http://sat.inesc-id.pt/open-wbo/\nc\n");
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

    BoolOption printmodel("Open-WBO", "print-model", "Print model.\n", true);

    StringOption printsoft("Open-WBO", "print-unsat-soft", "Print unsatisfied soft claues in the optimal assignment.\n", NULL);

    IntOption verbosity("Open-WBO", "verbosity",
                        "Verbosity level (0=minimal, 1=more).\n", 0,
                        IntRange(0, 1));

    IntOption algorithm("Open-WBO", "algorithm",
                        "Search algorithm "
                        "(0=wbo,1=CBLIN,2=linear-su,3=msu3,4=part-msu3,5=oll,6=oll_iter, 7=best)."
                        "\n",
                        1, IntRange(0, 6));

    IntOption partition_strategy("PartMSU3", "partition-strategy",
                                 "Partition strategy (0=sequential, "
                                 "1=sequential-sorted, 2=binary)"
                                 "(only for unsat-based partition algorithms).",
                                 2, IntRange(0, 2));

    IntOption graph_type("PartMSU3", "graph-type",
                         "Graph type (0=vig, 1=cvig, 2=res) (only for unsat-"
                         "based partition algorithms).",
                         2, IntRange(0, 2));

    BoolOption bmo("Open-WBO", "bmo", "BMO search.\n", true);

    IntOption cardinality("Encodings", "cardinality",
                          "Cardinality encoding (0=cardinality networks, "
                          "1=totalizer, 2=modulo totalizer).\n",
                          1, IntRange(0, 2));

    IntOption amo("Encodings", "amo", "AMO encoding (0=Ladder).\n", 0,
                  IntRange(0, 0));

    IntOption pb("Encodings", "pb", "PB encoding (0=SWC,1=GTE,2=Adder).\n", 1,
                 IntRange(0, 2));

    IntOption formula("Open-WBO", "formula",
                      "Type of formula (0=WCNF, 1=OPB).\n", 0, IntRange(0, 1));

    IntOption weight(
        "WBO", "weight-strategy",
        "Weight strategy (0=none, 1=weight-based, 2=diversity-based).\n", 2,
        IntRange(0, 2));

    BoolOption symmetry("WBO", "symmetry", "Symmetry breaking.\n", true);
    
        IntOption symmetry_lim(
        "WBO", "symmetry-limit",
        "Limit on the number of symmetry breaking clauses.\n", 500000,
        IntRange(0, INT32_MAX));




    IntOption pmreslin("PMRES", "pmreslin", "Run linear search in conjunction with PMRES: "
                                            "(0=not att all, 1=first cores then lin 2=only lins) .\n", 1,
                  IntRange(0, 3));

   BoolOption pmreslin_delsol("PMRES", "pmreslin-del", "Delete Solver between core guided and linear search.\n", true);
   BoolOption pmreslin_varres("PMRES", "pmreslin-varres", "Do varying resolution.\n", true);
   BoolOption pmreslin_relax2strat("PMRES", "pmreslin-r-b-s", "RelaxCores before strat.\n", false);
   BoolOption pmreslin_varresCG("PMRES", "pmreslin-varresCG", "Do varying resolution for CG.\n", false);
   BoolOption pmreslin_incvarres("PMRES", "pmreslin-v-inc", "Do varying resolution incrementally.\n", false);
   IntOption pmreslin_cgLim("PMRES", "pmreslin-cglim", "Time limit for core guided phase (s): "
                                            "(-1=not at all) .\n", 30,
                  IntRange(-1, INT_MAX));
    
  BoolOption prepro_rec("PMRES", "preprocess-rec", "Reconstruct solutions during preprocessing.\n", false);
  BoolOption prepro_min("PMRES", "preprocess-min", "Minimize solutions locally.\n", true);
  IntOption prepro_min_strat("PMRES", "preprocess-min-strat", "1=agressive (all solutions), 2=only the two first in each resolution: "
                                            "(0=only the best after each resolution) .\n", 0,
                  IntRange(0, 2));


   BoolOption preprocess("PMRES", "preprocess", "Preprocess the instance.\n", true);

   
    parseOptions(argc, argv, true);

    double initial_time = cpuTime();
    MaxSAT *S = NULL;

    switch ((int)algorithm) {
    case _ALGORITHM_WBO_:
      S = new WBO(verbosity, weight, symmetry, symmetry_lim);
      break;
    
    case _ALGORITHM_PMRES_:
      S = new CBLIN(verbosity, weight, pmreslin, pmreslin_delsol, pmreslin_varres, pmreslin_varresCG, 
                    pmreslin_cgLim, pmreslin_relax2strat, pmreslin_incvarres, preprocess, prepro_rec, prepro_min,prepro_min_strat );
      break;

    case _ALGORITHM_LINEAR_SU_:
      S = new LinearSU(verbosity, bmo, cardinality, pb);
      break;

    case _ALGORITHM_PART_MSU3_:
      S = new PartMSU3(verbosity, partition_strategy, graph_type, cardinality);
      break;

    case _ALGORITHM_MSU3_:
      S = new MSU3(verbosity);
      break;

    case _ALGORITHM_OLL_:
      S = new OLL(verbosity, cardinality);
      break;
    
    case _ALGORITHM_OLLITER_:
      S = new OLL_ITER(verbosity, cardinality, preprocess, prepro_rec);
      break;

    case _ALGORITHM_BEST_:
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
      parseMaxSATFormula(in, maxsat_formula);
      maxsat_formula->setFormat(_FORMAT_MAXSAT_);
    } else {
      ParserPB *parser_pb = new ParserPB();
      parser_pb->parsePBFormula(argv[1], maxsat_formula);
      maxsat_formula->setFormat(_FORMAT_PB_);
    }
    gzclose(in);

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
    double parsed_time = cpuTime();

    printf("c |  Parse time:           %12.2f s                                "
           "                                 |\n",
           parsed_time - initial_time);
    printf("c |                                                                "
           "                                       |\n");

    if (algorithm == _ALGORITHM_BEST_) {
      assert(S == NULL);

      if (maxsat_formula->getProblemType() == _UNWEIGHTED_) {
        // Unweighted
        S = new PartMSU3(_VERBOSITY_MINIMAL_, _PART_BINARY_, RES_GRAPH,
                         cardinality);
        S->loadFormula(maxsat_formula);

        if (((PartMSU3 *)S)->chooseAlgorithm() == _ALGORITHM_MSU3_) {
          // FIXME: possible memory leak
          S = new MSU3(_VERBOSITY_MINIMAL_);
        }

      } else {
        // Weighted
        S = new OLL(_VERBOSITY_MINIMAL_, cardinality);
      }
    }

    if (S->getMaxSATFormula() == NULL)
      S->loadFormula(maxsat_formula);
    S->setPrintModel(printmodel);
    S->setPrintSoft((const char *)printsoft);
    S->setInitialTime(initial_time);
    mxsolver = S;
    mxsolver->setPrint(true);

    int ret = (int)mxsolver->search();
    delete S;
    return ret;
  } catch (OutOfMemoryException &) {
    sleep(1);
    printf("c Error: Out of memory.\n");
    printf("s UNKNOWN\n");
    exit(_ERROR_);
  } catch(MaxSATException &e) {
    sleep(1);
    printf("c Error: MaxSAT Exception: %s\n", e.getMsg());
    printf("s UNKNOWN\n");
    exit(_ERROR_);
  }
}
