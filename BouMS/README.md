# BouMS

Bounded MaxSAT solving library, where bounded refers to

- time: only `for` loops with clear upper bounds are used to facilitate timing analysis.
- memory: given the number of clauses and variables, the library can predetermine its memory consumption.

The MaxSAT solving algorithm falls into the incomplete/anytime category, i.e., it

- can find feasible solutions quickly
- improves the solution quality over time

but in general it cannot determine whether a solution

- is optimal
- exists.

The *app* directory includes a solver application that uses the library.

## Build Instructions

Building is only tested on GNU/Linux.
The main reason is that BouMS requires [GNU Argp](https://www.gnu.org/software/libc/manual/html_node/Argp.html).
Building on macOS could be possible with [argp-standalone](https://formulae.brew.sh/formula/argp-standalone).
For Windows users, we recommend setting up the
[Windows subsytem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install).

BouMS uses [CMake](https://cmake.org/) as its build system
and provides [CMake presets](https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html) for different build
configurations.
Available presets are listed in *CMakePresets.json*.

1. Decide which preset you want to build (in the following, `release` is used as an example).
2. `cmake -S . --preset=release`
   This generates the build files in a new directory `build/release`.
3. `cmake --build build/release`
   This starts the actual build.

Alternatively, there is a makefile that automates these steps.
The default build target is `release-loginfo`,
but there are many more predefined configurations that can be found in the makefile.
To build the default target using the makefile, just run `make`.
To build another target, e.g., `debug-logverbose`, run `make debug-logverbose`.
To remove any and all build artifacts (i.e., the build directory), run `make clean`.

The build presets control the build process by setting compile-time definitions.


### Build Definitions for both Library and Application

- `BOUMS_LOG_LEVEL` takes a value between 0 (no logging) and 4 (trace logging). The log level is specified at compile
  time to allow for removing all logging code and thus the dependence on `stdio.h`.
- `BOUMS_PROFILING` adds the `-pg` flag for the compiler and linker for profiling with
  [GNU gprof](https://ftp.gnu.org/old-gnu/Manuals/gprof-2.9.1/html_mono/gprof.html)

### Build Definitions for the Library

- `BOUMS_NOFFA` makes the solver pick the variable with highest score (from some unsatisfied clause) instead of the
  least frequently flipped one in local optima. In benchmark results seen so far this may slightly increase the
  performance on unweighted instances, but significantly decrease it for weighted ones.

### Build Definitions for the Application

- `BOUMS_SATSOLVER_CADICAL` makes the solver use [CaDiCaL](https://github.com/arminbiere/cadical) once to find an
  initial, feasible assignment.

## Application Programming Interface

The header files of BouMS can be found in the `inc` directory.
In the following application programming interface (API) description, paths are relative to this directory.
The core API is defined in `BouMS/BouMS.h` and `BouMS/wcnf.h`;
headers in the `private` directory are intended for internal use only.
The `CMakeLists.txt` specifies a build target named `docs` that creates [Doxygen](https://www.doxygen.nl/)
documentation in HTML and LaTeX (to be found in the build directory's subdirectories `html` and `latex`).
The BouMS application serves as an example on how to use the API, its sources are located in the `app` directory.

### Constructing an Instance

Before solving a MaxSAT instance, it needs to be constructed in memory.
BouMS uses the type `BouMS_wcnf_t`, defined in `BouMS/wcnf.h`, for that.
Its main purpose is to hold arrays of the instance's variables and clauses.

A clause, represented by `BouMS_wcnf_clause_t`, holds a weight and the clause's literals.
For hard clauses, set the weight to `BOUMS_HARD_CLAUSE_WEIGHT`.

A literal, represented by `BouMS_wcnf_literal_t`,
holds the literal itself and a pointer to the clause it appears in.
The literal is represented by `BouMS_literal_t`, and can be constructed using the function `BouMS_mkLit`.
Note, that variable numbers are also used for indexing in BouMS and are 0-based.

Variables, represented by `BouMS_wcnf_variable_t`,
hold their value as well as an array of pointers to all the literals they occur in.

The `BouMS/wcnf.h` file also contains functions to work with the types its specifies.
For more details, please refer to the file or the Doxygen documentation.

Since it can be tedious and error-prone to set up an instance manually,
the `BouMS/wcnf_util.h` header contains convenience functions that you can use if you have access to functions
(that behave like) `realloc` and `free`.
As a first step, initialize a new instance with `BouMS_wcnf_util_newFormula` and set up batch clause adding with
`BouMS_wcnf_util_startBatchClauseAdding`.
Clauses can also be added using `BouMS_wcnf_util_addClause`, but batch clause adding requires fewer reallocations
and should consequently be faster if you intend to add more than one clause.
Clauses are added to the batch using `BouMS_wcnf_util_batchAddClause`,
and finally, the batch is added to the instance with `BouMS_wcnf_util_finishBatchClauseAdding`.
Note, that for these convenience functions,
literals are represented by signed integers and variable numbers are 1-based,
as in the [MaxSAT Evaluation's](https://maxsat-evaluations.github.io/) input file format.
After solving, memory should be freed with `BouMS_wcnf_util_deleteFormula`.
Please refer to the header or the Doxygen documentation for more details.

### Solving an Instance

The solver API is contained in `BouMS/BouMS.h`.
As a first step, call `BouMS_init` to initialize the library.
*Currently, this only calls `srand(time(NULL))`, but this is subject to change.
If your application already calls `srand`, you can and should skip this step.*

Next, having constructed an instance already,
call `BouMS_calcMemoryRequirements` to calculate the amount of memory BouMS will require for solving.
You can use this number to allocate a sufficiently-sized chunk of memory,
or check that your statically allocated memory is indeed large enough.
In any way, calling this function is required as it also fills a `BouMS_memoryReq_t` structure that the solving
routine requires to further subdivide the memory chunk it is given.

Then, specify the parameters for solving in a `BouMS_params_t` structure.
For detailed information on the available parameters, please consult the header file or Doxygen documentation,
and `BouMS/clause_weighting_params.md`.
The function `BouMS_params` sets reasonable default parameters.

Finally, you can call `BouMS_solve` with the instance, memory information, and parameters.
You can also pass an initial variable assignment that the solver should use as a starting point
(e.g., obtained from a SAT solver), and a pointer to a stop flag that is regularly checked.
When the stop flag is true (e.g., set from a signal handler, interrupt or another thread),
the solving routine exits as soon as possible.

### Accessing Solver Internals

It is possible to inspect the solver's internal memory using functionality from the `private/types.h` header.
There, `BouMS_memory_t` is defined as a structure that holds all the information the solver relies on during solving.
Using `BouMS_initMemory` on the pointer to the memory chunk you passed to `BouMS_solve` rebuilds this data structure
(by merely setting pointers accordingly).
For more detailed information, pleaser refer to the file.
