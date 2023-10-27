# Loandra Anytime MaxSAT solver
## Version 2.0 June 2023

The version of the Master branch was finally updated to match the 2022 version of Loandra. 
In addition to a new version of MaxPre, there's been a numerous quality of life improvements and better naming of options. 

### Citing
If you use Loandra in your work please cite: 

- Berg, J., Demirović, E. and Stuckey, P.J., 2019. Core-boosted linear search for incomplete MaxSAT. In Integration of Constraint Programming, Artificial Intelligence, and Operations Research: 16th International Conference, CPAIOR 2019, Thessaloniki, Greece, June 4–7, 2019, Proceedings 16 (pp. 39-56). Springer International Publishing.

### Building
Clone this repo with submodules. This can be done using `git clone --recurse-submodules git@github.com:jezberg/loandra.git`.
Alternatively, make sure to go to the `maxpre2` folder, and run the `git submodule init` and `git submodule update`.

Afterwards, run make in the base folder.

### Most significant command-line arguments:
Please also see the command-line arguments for open-wbo and maxpre 2

- **cb-i-varres, no-cb-i-varres**           (default: off)
  - Do varying resolution incrementally, i.e. do not delete the SAT solver between resolutions. 
  - In practice incremental varying resolution by adding an extra assumption variable to all PB encodings, which is then used to delete clauses afterwards. 

- **cb-varCG, no-cb-varCG**                 (default: off)
  - Do a varying resolution inspired stratifictation in the core-guided phase. 
  - In practice this is a more aggressive stratification heuristic that decreases the weight more between all interations. 

- **cb-r-2-s, no-cb-r-2-s**                 (default: off)
  - Relax the cores found on each stratification level before lowering the stratification bound. 
  - The bound is only lowered after no more cores can be found. 

- **cb-varres, no-cb-varres**               (default: on)
  - Do varying resolution. Note that if not, then the full PB encoding is immediately built, which memory intensive. 


- **cb-del, no-cb-del**                     (default: on)
  - Reinitialize the SAT solver between core guided and linear search.

- **cb-cglim**     = <int32>  [  -1 .. imax] (default: 30)
  - Time limit for core guided phase (s): (-1=unlimited) .

- **cb**           = <int32>  [   0 ..    3] (default: 1)
  - Strategy for running linear search in conjunction with PMRES: (0=not at all, 1=first core-guided then linear search.  2=only linnear search) .

- **preprocess, no-preprocess**             (default: on)
  - Preprocess the instance with MaxPRE before search.

- **pr-rec, no-pr-rec**                     (default: off)
  - Reconstruct solutions prior to computing their costs. (Only applicable when preprocessing)


- **pr-min, no-pr-min**                     (default: on)
  - Attempt to find trivial inprovements to preprocessed instances by flipping objective literals to false. 

- **pr-min-stra**t** = <int32>  [   0 ..    2] (default: 0)
  - Which solutions should be minimized: 0=only the best after each resolution, 1=all solutions, 2=only the two first in each resolution.

- **pr-tech**    = <string>
  - Which techniques should maxpre run? Please see MaxPREs documentation for details. 

### Additional references

- Ihalainen, H., Berg, J. and Järvisalo, M., 2022, August. Clause Redundancy and Preprocessing in Maximum Satisfiability. In Automated Reasoning: 11th International Joint Conference, IJCAR 2022, Haifa, Israel, August 8–10, 2022, Proceedings (pp. 75-94). Cham: Springer International Publishing.

- Demirović, E. and Stuckey, P.J., 2019. Techniques inspired by local search for incomplete maxsat and the linear algorithm: Varying resolution and solution-guided search. In Principles and Practice of Constraint Programming: 25th International Conference, CP 2019, Stamford, CT, USA, September 30–October 4, 2019, Proceedings 25 (pp. 177-194). Springer International Publishing.



## Version 1.0 July 2019

Loandra is an extension of the OpenWBO MaxSAT solver to core-boosted linear search. Loandra was one of the best performing solvers in the incomlete category of the 2019 MaxSAT Evaluation. 

This project includes the 2.1 version of Open-WBO as well as the algorithm CBLIN that implements core-boosted linear search as describedi in [1].

## Command line arguments 
Loandra accepts all of the same arguments as Open wbo + the following:


### core guided, core boosted or linear search
```-pmreslin = <int32> [0 (core-guided), 1 (core-boosted), 2 (only linear search)] (default 1)``` 


### set time limit (seconds) for core guided phase during core-boosted search 
```-pmreslin-cglim" = <int32> [-1 (unlimited, run intil stratification bound gets to 1), ...] (default 30) ```
   
   
### delete SAT solver between core-guided and linear phases
```-pmreslin-del, -no-pmreslin-del (default on) ```

### Do varying resolution during linear search: 
```-pmreslin-varres, -no-pmreslin-varres (default on) ```

### Encode maxres before lowering strat bound: 
```-pmreslin-relax2strat, -no-pmreslin-relax2strat (default off) ```

###  Keep satsolver between resolutions: 
```-pmreslin-incvarres, -no-pmreslin-incvarres (default off) ```

We thank the developers of Open-WBO!

[1] Berg, Jeremias & Demirović, Emir & Stuckey, Peter. (2019). Core-Boosted Linear Search for Incomplete MaxSAT. 10.1007/978-3-030-19212-9_3. 



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
