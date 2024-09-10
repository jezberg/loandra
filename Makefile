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
CADICAL_DIR = cadical
DPW_DIR 	= rustsat
DEPDIR     += mtl utils core
DEPDIR     +=  ../../encodings ../../algorithms ../../graph ../../classifier
## potential bug, depdir has  ../../maxpre2
MROOT      ?= $(PWD)/solvers/$(SOLVERDIR)
LFLAGS     += -lgmpxx -lgmp -pthread -ldl
CFLAGS     += -Wall -Wc++17-extensions -Wno-parentheses -std=c++11 -DNSPACE=$(NSPACE) -DSOLVERNAME=$(SOLVERNAME) -DVERSION=$(VERSION)
ifeq ($(VERSION),simp)
DEPDIR     += simp
CFLAGS     += -DSIMP=1 
ifeq ($(SOLVERDIR),glucored)
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
