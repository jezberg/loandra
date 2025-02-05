# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Released versions can be found on [zenodo](https://zenodo.org/records/14191233).

## [Unreleased]

### Added

- Build presets for building only the library.
- C++ compatibility (`extern "C"`).

### Changed

- `BouMS_params` ignores zero-weight soft clauses.
- `BouMS_params` condition for weighted instances: an instance is weighted when it has at least
  - old: one soft clause weight not equal to one.
  - new: two different soft clause weights.

### Fixed

- Division-by-0 when there are zero-weight soft clauses so that the average soft clause weight is less than 1.
  The solution is to ignore those clauses.

## [2.0.0](https://doi.org/10.5281/zenodo.14744890) - 2025-01-27

### Added

- Early termination criterion `maxNonImprovingFlips`.
- Changelog

### Changed

- Move `initMemory` function to *private/types.h* to allow users to inspect internal data structures.
- Number of safety bits for fixed-precision scheme is now user-controlled parameter instead of compile-time constant.
- Improve README.md

### Fixed

- Replacing fixed-precision scheme with `double` through `FIXEDPREC_FLOATING`.

## [1.0.0](https://doi.org/10.5281/zenodo.14191233) - 2024-11-20

Initial release as BouMS, formerly known as [noSAT-MaxSAT](https://doi.org/10.15480/882.13564).
Changes listed here refer to the MSE2024 version (noSAT-MaxSATv3).

### Added

- Information on clause weighting scheme and parameters in *inc/BouMS/clause_weighting_params.md*.
- Convenience functions for adding clauses when the system has `realloc` and `free`.
- Compile-time switch to disable FFA-inspired variable selection, `BOUMS_NOFFA`.
- Compile-time switch to pre-solve hard clauses with [CaDiCaL](https://github.com/arminbiere/cadical),
  `BOUMS_SATSOLVER_CADICAL`.

### Changed

- Check `stop` flag in more places (input reading, preprocessing).
- Weight initialization for PMS (parameterized).
- Propagate user parameters to hard clause solving as well.
- Exit with `EXIT_FAILURE` in case of parse error or input file not found.

