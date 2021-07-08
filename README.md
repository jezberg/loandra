# Loandra Incomplete MaxSAT solver
## Version 1.0 July 2019

Loandra is an extension of the OpenWBO MaxSAT solver to core-boosted linear search. Loandra was one of the best performing solvers in the incomlete category of the 2019 MaxSAT Evaluation and the 2021 MaxSAT Evaluation-

This project includes the 2.1 version of Open-WBO as well as the algorithm CBLIN that implements core-boosted linear search as describedi in [1].
THis branch also includes the MaxPRE [2] preprocessor.

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
```-pmreslin-r-b-s, -no-pmreslin-r-b-s (default off) ```

###  Keep satsolver between resolutions: 
```-pmreslin-v-inc, -no-pmreslin-v-inc (default off) ```

## preprocess the instance with MaxPre
```-preprocess, -no-preprocess (default on) ```

## reconstruct all solutions found during search 
```-preprocess-rec, -no-preprocess-rec (default on) ```

## locally minimize solutions
```-preprocess-min, -no-preprocess-min (default on) ```

We thank the developers of Open-WBO!

## Citing Loandra
If you use Loandra in your research, we kindly ask you to cite the follwoing publication 

[1] Berg, Jeremias & Demirović, Emir & Stuckey, Peter. (2019). Core-Boosted Linear Search for Incomplete MaxSAT. 10.1007/978-3-030-19212-9_3. 

[2] Korhonen, T., Berg, J., Saikko, P. and Järvisalo, M., (2017), MaxPre: an extended MaxSAT preprocessor. 10.1007/978-3-319-66263-3\_28. 

### Additional references
[3] Berg, J. and Järvisalo, M., (2017), August. Weight-aware core extraction in SAT-based MaxSAT solving. 10.1007/978-3-319-66158-2\_42

[4] Demirović, Emir & Stuckey, Peter. (2019). Techniques Inspired by Local Search for Incomplete MaxSAT and the Linear Algorithm: Varying Resolution and Solution-Guided Search. 10.1007/978-3-030-30048-7_11

[5] Leivo, M., Berg, J. and Järvisalo, M., (2020). Preprocessing in incomplete maxsat solving, 10.3233/FAIA200112

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
