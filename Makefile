# Define which solver to use as backend, this can be a name of a file in the
# solvers directory.
SOLVER     ?= glucose4.1
#
# The following values should be defined in the included file:
# VERSION    = core or simp 
# SOLVERNAME = name of the SAT solver
# SOLVERDIR  = subdirectory of the SAT solver
# NSPACE     = namespace of the SAT solver
#
include $(PWD)/solvers/$(SOLVER).mk

# THE REMAINING OF THE MAKEFILE SHOULD BE LEFT UNCHANGED
EXEC       = loandra
PREPRO_DIR = maxpre2
DPW_DIR = rustsat/target/release
DEPDIR     += mtl utils core
#LIBD 	   = -L $(PWD)/maxpre/lib -l preprocessor
##INCLUDES   = $(PWD)/maxpre2/
DEPDIR     +=  ../../encodings ../../algorithms ../../graph ../../classifier
MROOT      ?= $(PWD)/solvers/$(SOLVERDIR)
LFLAGS     += -lgmpxx -lgmp
CFLAGS     += -Wall -Wno-parentheses -std=c++11 -DNSPACE=$(NSPACE) -DSOLVERNAME=$(SOLVERNAME) -DVERSION=$(VERSION)
ifeq ($(VERSION),simp)
DEPDIR     += simp
CFLAGS     += -DSIMP=1 
ifeq ($(SOLVERDIR),glucored)
LFLAGS     += -pthread
CFLAGS     += -DGLUCORED
DEPDIR     += reducer glucored
endif
endif



# Some solvers do not have a template.mk file any more
# E.g.: Minisat or Riss
ifeq ($(SOLVERDIR),$(filter $(SOLVERDIR),minisat riss))
include $(PWD)/mtl/template.mk
else
include $(MROOT)/mtl/template.mk
endif
