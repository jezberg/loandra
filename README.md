# SLS-Enhanced Core-Boosted Linear Search for Anytime Maximum Satisfiability -- Implementation

Implementation of the techniques introduced in

*Ole Lübke and Jeremias Berg. SLS-Enhanced Core-Boosted Linear Search for Anytime Maximum Satisfiability.
31st International Conference on Principles and Practice of Constraint Programming (CP 2025).*

## Build Instructions

1. Ensure you have [Rust](https://www.rust-lang.org/tools/install), `make`, and a C++ compiler installed.
2. `./configure`
3. `make rs`
4. The compiled solver is called `loandra_static`.

## Running the Solver Configurations from the Paper

Each solver configuration from the paper corresponds to a certain set of command line arguments.
Almost all configurations use `-no-preprocess`.
Indeed, this is the only flag required for the Base solver.
The following flags correspond to the techniques from the paper:

- Init: `-ls-init-level=1`
- ExtInit: `-ls-init-level=2`
- AssignExt: `ls-extend`
- IncPre: `-ls-dyn-prec`
- CoreClauses: `-ls-cores=1`
- CoreScores: `-ls-cores=2 -ls-cores-factor=1`

The Preproc configuration is achieved by omitting `-no-preprocess`.
Combinations of techniques are achieved by combining the arguments, e.g., for ExtInit+AssignExt+CoreClauses:
`-no-preprocess -ls-init-level=2 -ls-extend -ls-cores=1`.

## More Information on Loandra

Please consult [Loandra's original README](README-Loandra.md), [it's GitHub page](https://github.com/jezberg/loandra),
and/or `loandra_static --help`.
