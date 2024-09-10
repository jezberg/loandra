# Loandra Anytime MaxSAT solver
## Version 2.2 June 2024

The master branch matches the 2024 version of Loandra. 
The main additions are the dynamic polynomial watchdog encoding. 

### Citing
If you use Loandra in your work please cite: 

- Berg, J., Demirović, E. and Stuckey, P.J., 2019. Core-boosted linear search for incomplete MaxSAT. In Integration of Constraint Programming, Artificial Intelligence, and Operations Research: 16th International Conference, CPAIOR 2019, Thessaloniki, Greece, June 4–7, 2019, Proceedings 16 (pp. 39-56). Springer International Publishing.

### Building
First ensure you have [Rust installed](https://www.rust-lang.org/tools/install). 
Clone this repo with submodules. This can be done using `git clone --recurse-submodules https://github.com/jezberg/loandra.git`.
Alternatively, make sure to go to the `maxpre2` and rustsat folders, and run the `git submodule init` and `git submodule update`.


Afterwards, run make in the base folder.

### Most significant command-line arguments:
Run ```./loandra --help or --help-verb``` for more information. 
Also see the command-line arguments for [open-wbo](https://github.com/sat-group/open-wbo) and [maxpre2](https://bitbucket.org/coreo-group/maxpre2/src/master/)

- **print-model, no-print-model** (default: off)
  - Print the final solution found after search. 

- **old-format, no-print-model** (default: off)
  - Parse WCNF files in the pre 2022 format.

- **cb-DPW, no-cb-DPW** (default: on)
  - Use the dynamic polynomial watchdog encoding during the solution improving phase. Otherwise, use the generalized totalizer. 

- **cb-r-2-s, no-cb-r-2-s**                 (default: off)
  - Relax the cores found on each stratification level before lowering the stratification bound. 
  - The bound is only lowered after no more cores can be found. 

- **cb-del, no-cb-del**                     (default: on)
  - Reinitialize the SAT solver between core guided and solution improving search.

- **cb-cglim**     = <int32>  [  -1 .. imax] (default: 30)
  - Time limit for core guided phase (s): (-1=unlimited) .

- **cb**           = <int32>  [   0 ..    3] (default: 1)
  - Strategy for running linear search in conjunction with PMRES: (0=not at all, 1=first core-guided then linear search.  2=only linnear search) .

- **preprocess, no-preprocess**             (default: on)
  - Preprocess the instance with MaxPRE before search.

- **pr-rec, no-pr-rec**                     (default: off)
  - Reconstruct solutions prior to computing their costs. (Only applicable when preprocessing)

- **pr-tech**    = <string>
  - Which techniques should maxpre run? Please see MaxPREs documentation for details. 

- **pr-min, no-pr-min**                     (default: on)
  - Attempt to find trivial inprovements to preprocessed instances by flipping objective literals to false. 

- **cb-local-search, no-cb-local-searchn**                     (default: on)
  - Call the NuWLS local search solver prior to SIS search and between each resolution. 

- **pr-min-strat** = <int32>  [   0 ..    2] (default: 0)
  - Which solutions should be minimized: 0=only the best after each resolution, 1=all solutions, 2=only the two first in each resolution.


### Additional references

- Ihalainen, H., Berg, J. and Järvisalo, M., 2022, August. Clause Redundancy and Preprocessing in Maximum Satisfiability. In Automated Reasoning: 11th International Joint Conference, IJCAR 2022, Haifa, Israel, August 8–10, 2022, Proceedings (pp. 75-94). Cham: Springer International Publishing.

- Demirović, E. and Stuckey, P.J., 2019. Techniques inspired by local search for incomplete maxsat and the linear algorithm: Varying resolution and solution-guided search. In Principles and Practice of Constraint Programming: 25th International Conference, CP 2019, Stamford, CT, USA, September 30–October 4, 2019, Proceedings 25 (pp. 177-194). Springer International Publishing.


# Open-WBO MaxSAT Solver
## Version 2.1 - September 2018

Open-WBO is an extensible and modular open-source MaxSAT Solver.
Open-WBO was one of the best solvers in the partial MaxSAT categories at 
MaxSAT Evaluations 2014, 2015, 2016 and 2017 and in the decision and 
optimization for SMALLINT categories at PB Evaluation 2016.

## MaxSAT Evaluation 2018
The default algorithms used by Open-WBO in the complete track are: 
* unweighted: Part-MSU3
* weighted: OLL

Usage of the solver:
./loandra [options] <input-file>

The following options are available in Open-WBO:

## Global Options
### Formula type (0=MaxSAT, 1=PB)
```-formula      = <int32>  [   0 ..    1] (default: 0)```

### Print model
```-print-model, -no-print-model (default on)```

### Write unsatisfied soft clauses to file
```-print-unsat-soft = <output-file>```

### Verbosity level (0=minimal, 1=more)
```-verbosity    = <int32>  [   0 ..    1] (default: 1)```

### Search algorithm (0=wbo,1=linear-su,2=msu3,3=part-msu3,4=oll,5=best)
```-algorithm    = <int32>  [   0 ..    1] (default: 5)```

### BMO search 
```-bmo,-no-bmo (default on)```

### Pseudo-Boolean encodings (0=SWC,1=GTE, 2=Adder)
```-pb           = <int32>  [   0 ..    1] (default: 1)```

### At-most-one encodings (0=ladder)
```-amo          = <int32>  [   0 ..    0] (default: 0)```

### Cardinality encodings (0=cardinality networks, 1=totalizer, 2=modulo totalizer)
```-cardinality  = <int32>  [   0 ..    2] (default: 1)```

       
## WBO Options (algorithm=0, unsatisfiability-based algorithm)
### Weight strategy (0=none, 1=weight-based, 2=diversity-based)
```-weight-strategy = <int32>  [   0 ..    2] (default: 2)```

### Symmetry breaking
```-symmetry, -no-symmetry (default on)```

### Limit on the number of symmetry breaking clauses
```-symmetry-limit = <int32>  [   0 .. imax] (default: 500000)```

## PartMSU3 OPTIONS (algorithm=3, partition-based algorithm)
### Graph type (0=vig, 1=cvig, 2=res)
```-graph-type   = <int32>  [   0 ..    2] (default: 2)```

### Partition strategy (0=sequential, 1=sequential-sorted, 2=binary)
```-partition-strategy = <int32>  [   0 ..    2] (default: 2)```

## Output of solver
Open-WBO follows the standard output of MaxSAT solvers:
* Comments ("c " lines) 
* Solution Status ("s " line):
  * s OPTIMUM FOUND : an optimum solution was found
  * s UNSATISFIABLE : the hard clauses are unsatisfiable
  * s SATISFIABLE   : a solution was found but optimality was not proven
* Solution Cost Line ("o " lines):
  * This represents the cost of the best solution found by the solver. The cost 
  of a solution is given by the sum of the weights of the unsatisfied soft clause.
* Solution Values (Truth Assignment) ("v " lines): 
  * This represents the truth assignment (true/false) assigned to each variable. 
  A literal is denoted by an integer that identifies the variable and the negation 
  of a literal is denoted by a minus sign immediately followed by the integer of 
  the variable.

> Authors: Ruben Martins, Vasco Manquinho, Ines Lynce

> Contributors: Miguel Neves, Norbert Manthey, Saurabh Joshi, Mikolas Janota

> To contact the authors please send an email to:  open-wbo@sat.inesc-id.pt
