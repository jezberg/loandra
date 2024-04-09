/*!
 * \author Jeremias Berg - jeremiasberg@gmail.com
 *
 * @section LICENSE
  * NuWLS -- Copyright (c) 2021-2022, Yi Chu, Xiang He
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

#include "Alg_NUWLS.h"

using namespace openwbo;

void NUWLS::build_nuwls_clause_structure(openwbo::MaxSATFormula* maxsat_formula)
{

  nuwls_nvars = maxsat_formula->nVars(); // - maxsat_formula->nSoft();
  nuwls_nclauses = maxsat_formula->nHard() + maxsat_formula->nSoft();
  nuwls_topclauseweight = maxsat_formula->getHardWeight();

  int nuwls_num_hclauses = maxsat_formula->nHard();
  nuwls_clause_lit = new clauselit *[nuwls_nclauses + 10];
  nuwls_clause_lit_count = new int[nuwls_nclauses + 10];
  nuwls_clause_weight = new long long[nuwls_nclauses + 10];

  int *redunt_test = new int[nuwls_nvars + 10];
  memset(redunt_test, 0, sizeof(int) * (nuwls_nvars + 10));

  int tem_v, tem_sense, tem_lit_count;
  bool clause_reduent;
  int c = 0;
  for (int i = 0; i < maxsat_formula->nHard(); ++i)
  {
    nuwls_clause_lit_count[c] = maxsat_formula->getHardClause(i).clause.size();
    nuwls_clause_lit[c] = new clauselit[nuwls_clause_lit_count[c] + 1];
    clause_reduent = false;
    tem_lit_count = 0;
    for (int j = 0; j < nuwls_clause_lit_count[c]; ++j)
    {
      tem_v = var(maxsat_formula->getHardClause(i).clause[j]) + 1;
      tem_sense = 1 - sign(maxsat_formula->getHardClause(i).clause[j]);
      if (redunt_test[tem_v] == 0)
      {
        redunt_test[tem_v] = tem_sense - 2;

        // nuwls_clause_lit[c][tem_lit_count].clause_num = c;
        nuwls_clause_lit[c][tem_lit_count].var_num = tem_v;
        nuwls_clause_lit[c][tem_lit_count].sense = tem_sense;

        tem_lit_count++;
      }
      else if (redunt_test[tem_v] == tem_sense - 2)
      {
        continue;
      }
      else
      {
        clause_reduent = true;
        break;
      }
    }
    for (int j = 0; j < nuwls_clause_lit_count[c]; ++j)
      redunt_test[var(maxsat_formula->getHardClause(i).clause[j]) + 1] = 0;
    if (clause_reduent == false)
    {
      nuwls_clause_weight[c] = nuwls_topclauseweight;
      nuwls_clause_lit[c][tem_lit_count].var_num = 0;
      // nuwls_clause_lit[c][tem_lit_count].clause_num = -1;
      nuwls_clause_lit_count[c] = tem_lit_count;
      c++;
    }
    else
    {
      delete nuwls_clause_lit[c];
    }
  }
  for (int i = nuwls_num_hclauses; i < nuwls_nclauses; ++i)
  {
    nuwls_clause_lit_count[c] = maxsat_formula->getSoftClause(i - nuwls_num_hclauses).clause.size();
    nuwls_clause_lit[c] = new clauselit[nuwls_clause_lit_count[c] + 1];
    clause_reduent = false;
    tem_lit_count = 0;
    for (int j = 0; j < nuwls_clause_lit_count[c]; ++j)
    {
      tem_v = var(maxsat_formula->getSoftClause(i - nuwls_num_hclauses).clause[j]) + 1;
      tem_sense = 1 - sign(maxsat_formula->getSoftClause(i - nuwls_num_hclauses).clause[j]);
      if (redunt_test[tem_v] == 0)
      {
        redunt_test[tem_v] = tem_sense - 2;

        // nuwls_clause_lit[c][tem_lit_count].clause_num = c;
        nuwls_clause_lit[c][tem_lit_count].var_num = tem_v;
        nuwls_clause_lit[c][tem_lit_count].sense = tem_sense;

        tem_lit_count++;
      }
      else if (redunt_test[tem_v] == tem_sense - 2)
      {
        continue;
      }
      else
      {
        clause_reduent = true;
        break;
      }
    }
    for (int j = 0; j < nuwls_clause_lit_count[c]; ++j)
      redunt_test[var(maxsat_formula->getSoftClause(i - nuwls_num_hclauses).clause[j]) + 1] = 0;
    if (clause_reduent == false)
    {
      nuwls_clause_weight[c] = maxsat_formula->getSoftClause(i - nuwls_num_hclauses).weight;
      nuwls_clause_lit[c][tem_lit_count].var_num = 0;
      // nuwls_clause_lit[c][tem_lit_count].clause_num = -1;
      nuwls_clause_lit_count[c] = tem_lit_count;
      c++;
    }
    else
    {
      delete nuwls_clause_lit[c];
    }
  }
  nuwls_nclauses = c;
}

void NUWLS::local_search() {
    start_timing();
    int time_limit_for_ls = NUWLS_TIME_LIMIT;
    if (if_using_neighbor)
    {
      for (int step = 1; step < max_flips; ++step)
      {
        if (hard_unsat_nb == 0)
        {
          local_soln_feasible = 1;
          if (soft_unsat_weight < opt_unsat_weight)
          {
            max_flips = step + max_non_improve_flip;
            time_limit_for_ls = get_runtime() + NUWLS_TIME_LIMIT;

            best_soln_feasible = 1;
            opt_unsat_weight = soft_unsat_weight;
            for (int v = 1; v <= num_vars; ++v) {
              best_soln[v] = cur_soln[v];
            }
            if (opt_unsat_weight == 0)
              break;
          }
        }
        int flipvar = pick_var();
        flip2(flipvar);
        time_stamp[flipvar] = step;

        if (step % 1000 == 0)
        {
          if (get_runtime() > time_limit_for_ls)
            break;
        }
      }
    }
    else
    {
      for (int step = 1; step < max_flips; ++step)
      {
        if (hard_unsat_nb == 0)
        {
          local_soln_feasible = 1;
          if (soft_unsat_weight < opt_unsat_weight)
          {
            max_flips = step + max_non_improve_flip;
            time_limit_for_ls = get_runtime() + NUWLS_TIME_LIMIT;

            best_soln_feasible = 1;
            opt_unsat_weight = soft_unsat_weight;
            for (int v = 1; v <= num_vars; ++v)
            {
              best_soln[v] = cur_soln[v];
            }

            if (opt_unsat_weight == 0)
              break;
          }
        }
        int flipvar = pick_var();
        flip(flipvar);
        time_stamp[flipvar] = step;

        if (step % 1000 == 0)
        {
          if (get_runtime() > time_limit_for_ls)
            break;
        }
      }
    }
    cout << "c nuwls search done!" << endl;
  }








