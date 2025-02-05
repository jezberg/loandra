cmake_minimum_required(VERSION 3.21)

project(BouMS C)

# options that apply both to the library and the application
option(BOUMS_LOG_LEVEL "Log level (0 = off, 4 = trace)" 2)
option(BOUMS_PROFILING "Build with BOUMS_PROFILING information" OFF)

# options that only apply to the library
option(BOUMS_NOFFA "Select variable with best score instead of the least frequently flipped one in local optima" OFF)

# options that only apply to the application
option(BOUMS_SATSOLVER_CADICAL "Use CaDiCaL to obtain a feasible initial assignment" OFF)

if(BOUMS_SATSOLVER_CADICAL)
  include(ExternalProject)
  include(FetchContent)
  set(BOUMS_USE_SATSOLVER ON)
endif()

# Flags

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_compile_options(-Wall -Wextra -Wpedantic -Werror)

# Library Documentation
find_package(Doxygen OPTIONAL_COMPONENTS dot)
if(DOXYGEN_FOUND)
  set(DOXYGEN_GENERATE_LATEX true)
  set(DOXYGEN_USE_MDFILE_AS_MAINPAGE README.md)
  doxygen_add_docs(docs
    README.md
    inc/BouMS/common.h
    inc/BouMS/wcnf.h
    inc/BouMS/BouMS.h
    inc/BouMS/preprocessing.h
    inc/BouMS/wcnf_util.h
    ALL
    USE_STAMP_FILE
)
endif(DOXYGEN_FOUND)

# Solver description
set(SOLVER_DESC_TEX ${CMAKE_CURRENT_SOURCE_DIR}/doc/mse/noSAT-MaxSAT.tex)
if(EXISTS ${SOLVER_DESC_TEX})
  find_program(LATEXMK latexmk)
  if(LATEXMK)
    set(SOLVER_DESC_PDF noSAT-MaxSAT.pdf)
    add_custom_command(OUTPUT ${SOLVER_DESC_PDF}
                      COMMAND latexmk -pdf -cd ${SOLVER_DESC_TEX} -outdir=${CMAKE_BINARY_DIR}
                      MAIN_DEPENDENCY ${SOLVER_DESC_TEX}
    )
    add_custom_target(SolverDescription DEPENDS ${SOLVER_DESC_PDF})
  endif(LATEXMK)
endif(EXISTS ${SOLVER_DESC_TEX})

# Profiling
if(BOUMS_PROFILING)
  add_compile_options(-pg)
  add_link_options(-pg)
endif(BOUMS_PROFILING)

# Static Analysis
find_program(CLANG_TIDY clang-tidy)
if(CLANG_TIDY)
  set(CMAKE_C_CLANG_TIDY clang-tidy -checks=bugprone-*,clang-analyzer-*,misc-*,performance-*,portability-*,readability-*,cert-*,misc-*,-bugprone-easily-swappable-parameters)
endif(CLANG_TIDY)

find_program(IWYU_PATH NAMES include-what-you-use iwyu)

# Library
add_library(BouMS
  lib/BouMS.c
  lib/common.c
  lib/decimation.c
  lib/filterClauses.c
  lib/fixedprec.c
  lib/nonpartial.c
  lib/partial.c
  lib/wcnf_util.c
)
target_include_directories(BouMS
  PUBLIC inc
)
if(BOUMS_LOG_LEVEL)
  target_compile_definitions(BouMS PRIVATE BOUMS_LOG_LEVEL=${BOUMS_LOG_LEVEL})
endif(BOUMS_LOG_LEVEL)
if(BOUMS_NOFFA)
  target_compile_definitions(BouMS PRIVATE BOUMS_NOFFA)
endif(BOUMS_NOFFA)
if(IWYU_PATH)
  set_property(TARGET BouMS PROPERTY C_INCLUDE_WHAT_YOU_USE ${IWYU_PATH})
endif(IWYU_PATH)

# Executable
add_executable(app
  app/main.c
  app/MaxSAT22.c
)
target_include_directories(app PRIVATE app)
set_target_properties(app
  PROPERTIES OUTPUT_NAME BouMS
)
target_link_libraries(app
  PUBLIC BouMS
)
if(BOUMS_LOG_LEVEL)
  target_compile_definitions(app PRIVATE BOUMS_LOG_LEVEL=${BOUMS_LOG_LEVEL})
endif(BOUMS_LOG_LEVEL)
if(IWYU_PATH)
  set_property(TARGET app PROPERTY C_INCLUDE_WHAT_YOU_USE ${IWYU_PATH})
endif(IWYU_PATH)

# Optional SAT Solvers
if(BOUMS_SATSOLVER_CADICAL)
  set(SATSOLVER_NAME CaDiCaL)
  set(SATSOLVER_GITURL https://github.com/arminbiere/cadical.git)
  set(SATSOLVER_GITTAG rel-2.1.0)
  set(SATSOLVER_LIBPATH build/libcadical.a)
endif()

if(BOUMS_USE_SATSOLVER)
  string(TOLOWER ${SATSOLVER_NAME} SATSOLVER_NAME_LC)
  string(CONCAT SATSOLVER_SOURCEDIR ${SATSOLVER_NAME_LC} _SOURCE_DIR)
  string(CONCAT SATSOLVER_LIBNAME Lib ${SATSOLVER_NAME})

  FetchContent_Declare(${SATSOLVER_NAME}
                       GIT_REPOSITORY ${SATSOLVER_GITURL}
                       GIT_TAG ${SATSOLVER_GITTAG}
                       GIT_PROGRESS ON
  )
  FetchContent_MakeAvailable(${SATSOLVER_NAME})

  ExternalProject_Add(${SATSOLVER_NAME}
                      SOURCE_DIR ${${SATSOLVER_SOURCEDIR}}
                      CONFIGURE_COMMAND cd <SOURCE_DIR> && ./configure
                      BUILD_COMMAND cd <SOURCE_DIR> && make -j${CMAKE_BUILD_PARALLEL_LEVEL} ${SATSOLVER_NAME_LC}
                      INSTALL_COMMAND ""
                      BUILD_BYPRODUCTS <SOURCE_DIR>/${SATSOLVER_LIBPATH}
  )
  add_library(${SATSOLVER_LIBNAME} STATIC IMPORTED)
  add_dependencies(${SATSOLVER_LIBNAME} ${SATSOLVER_NAME})
  set_target_properties(${SATSOLVER_LIBNAME} PROPERTIES
                        IMPORTED_LOCATION ${${SATSOLVER_SOURCEDIR}}/${SATSOLVER_LIBPATH}
                        INTERFACE_INCLUDE_DIRECTORIES ${${SATSOLVER_SOURCEDIR}}/src
  )

  target_sources(app PRIVATE app/externalsat.c)
  target_link_libraries(app
                        PRIVATE ${SATSOLVER_LIBNAME}
                        PRIVATE m
  )
  target_compile_definitions(app PRIVATE BOUMS_USE_SATSOLVER)
  if(BOUMS_SATSOLVER_CADICAL)
    target_compile_definitions(app PRIVATE BOUMS_SATSOLVER_CADICAL)
    target_link_libraries(app PRIVATE stdc++)
  endif()
endif(BOUMS_USE_SATSOLVER)
