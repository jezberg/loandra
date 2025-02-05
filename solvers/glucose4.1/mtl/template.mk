##
##  Template makefile for Standard, Profile, Debug, Release, and Release-static versions
##
##    eg: "make rs" for a statically linked release version.
##        "make d"  for a debug version (no optimizations).
##        "make"    for the standard version (optimized, but with debug information and assertions active)

PWD        = $(shell pwd)
EXEC      ?= $(notdir $(PWD))

CSRCS      = $(wildcard $(PWD)/*.cc) 
DSRCS      = $(foreach dir, $(DEPDIR), $(filter-out $(MROOT)/$(dir)/Main.cc, $(wildcard $(MROOT)/$(dir)/*.cc)))
CHDRS      = $(wildcard $(PWD)/*.h)
COBJS      = $(CSRCS:.cc=.o) $(DSRCS:.cc=.o)
PREOBJ	   = $(wildcard $(PREPRO_DIR)/src/lib/*.a) 
CADOBJ	   = $(wildcard $(CADICAL_DIR)/build/*.a) 

DPWOBJ	   = $(wildcard $(DPW_DIR)/target/release/*.a) 
BOUMSRELOBJ = $(BOUMS_DIR)/build/release-nologging/libBouMS.a
BOUMSDBGOBJ = $(BOUMS_DIR)/build/debug-logverbose/libBouMS.a

PCOBJS     = $(addsuffix p,  $(COBJS))
DCOBJS     = $(addsuffix d,  $(COBJS))
RCOBJS     = $(addsuffix r,  $(COBJS))

CXX       ?= g++
CFLAGS    ?= -Wall -Wno-parentheses -std=c++11
LFLAGS    ?= -Wall -pthread 

COPTIMIZE ?= -O3

CFLAGS    += -I$(MROOT) -D __STDC_LIMIT_MACROS -D __STDC_FORMAT_MACROS
LFLAGS    += -lz

.PHONY : s p d r rs clean 

s:  BOUMSOBJ=$(BOUMSRELOBJ) 
p:  BOUMSOBJ=$(BOUMSRELOBJ) 
d:  BOUMSOBJ=$(BOUMSDBGOBJ) 
r:  BOUMSOBJ=$(BOUMSRELOBJ) 
rs:	BOUMSOBJ=$(BOUMSRELOBJ) 

s:  builddeps $(EXEC) 
p:  builddeps $(EXEC)_profile
d:  builddeps $(EXEC)_debug 
r:  builddeps $(EXEC)_release 
rs: builddeps $(EXEC)_static

libs:	lib$(LIB)_standard.a
libp:	lib$(LIB)_profile.a
libd:	lib$(LIB)_debug.a
libr:	lib$(LIB)_release.a

## Compile options
%.o:			CFLAGS +=$(COPTIMIZE) -g -D DEBUG
%.op:			CFLAGS +=$(COPTIMIZE) -pg -g -D NDEBUG
%.od:			CFLAGS +=-O0 -g -D DEBUG
%.or:			CFLAGS +=$(COPTIMIZE) -g -D NDEBUG

## Link options
$(EXEC):		LFLAGS += -g
$(EXEC)_profile:	LFLAGS += -g -pg
$(EXEC)_debug:		LFLAGS += -g
#$(EXEC)_release:	LFLAGS += ...
$(EXEC)_static:		LFLAGS += --static

## Dependencies
$(EXEC):		$(COBJS)
$(EXEC)_profile:	$(PCOBJS)
$(EXEC)_debug:		$(DCOBJS)
$(EXEC)_release:	$(RCOBJS)
$(EXEC)_static:		$(RCOBJS)

lib$(LIB)_standard.a:	$(filter-out */Main.o,  $(COBJS))
lib$(LIB)_profile.a:	$(filter-out */Main.op, $(PCOBJS))
lib$(LIB)_debug.a:	$(filter-out */Main.od, $(DCOBJS))
lib$(LIB)_release.a:	$(filter-out */Main.or, $(RCOBJS))


## Build rule
%.o %.op %.od %.or:	%.cc
	@echo Compiling: $(subst $(MROOT)/,,$@)
	@$(CXX) $(CFLAGS) -c -o $@ $<

## Linking rules (standard/profile/debug/release)
$(EXEC) $(EXEC)_profile $(EXEC)_debug $(EXEC)_release $(EXEC)_static: 
	@echo Linking: "$@ ( $(foreach f,$^,$(subst $(MROOT)/,,$f)) )"
	@echo preprocessor and DPW library: $(DPWOBJ)  $(PREOBJ)
	@echo CaDiCaL: $(CADOBJ)
	@echo BouMS library: $(BOUMSOBJ)
	@echo @$(CXX) $^ $(DPWOBJ) $(PREOBJ) $(BOUMSOBJ) $(LFLAGS) -o $@  
	@$(CXX) $^ $(DPWOBJ) $(PREOBJ) $(CADOBJ) $(BOUMSOBJ) $(LFLAGS) -o $@  

## Library rules (standard/profile/debug/release)
lib$(LIB)_standard.a lib$(LIB)_profile.a lib$(LIB)_release.a lib$(LIB)_debug.a:
	@echo Making library: "$@ ( $(foreach f,$^,$(subst $(MROOT)/,,$f)) )"
	@$(AR) -rcsv $@ $^

## Library Soft Link rule:
libs libp libd libr:
	@echo "Making Soft Link: $^ -> lib$(LIB).a"
	@ln -sf $^ lib$(LIB).a

## Clean rule
allclean: clean
	@rm -f ../simp/*.o ../simp/*.or ../simp/*.od  ../core/*.o ../core/*.or ../core/*.od
	$(MAKE) -C $(PREPRO_DIR) clean
	$(MAKE) -C $(CADICAL_DIR) clean
	cd $(DPW_DIR)/capi && cargo clean

clean:
	rm -f $(EXEC) $(EXEC)_profile $(EXEC)_debug $(EXEC)_release $(EXEC)_static \
	  $(COBJS) $(PCOBJS) $(DCOBJS) $(RCOBJS) *.core depend.mk
	$(MAKE) -C $(PREPRO_DIR) clean
	cd $(DPW_DIR)/capi && cargo clean
	$(MAKE) -C $(BOUMS_DIR) clean

builddeps:
	@echo Making MaxPre
	$(MAKE) -C $(PREPRO_DIR) lib with_zlib=false
	@echo Making RustSAT
	cd $(DPW_DIR)/capi && cargo build --release
	@echo Making cadical
	$(MAKE) -C $(CADICAL_DIR)
	if [ $(BOUMSOBJ) = $(BOUMSRELOBJ) ]; then \
		@echo "Making BouMS (release-nologging)"; \
		$(MAKE) -C $(BOUMS_DIR) libonly-release-nologging; \
	else \
		@echo "Making BouMS (debug-logverbose)"; \
		$(MAKE) -C $(BOUMS_DIR) libonly-debug-logverbose; \
	fi

## Make dependencies
depend.mk: $(CSRCS) $(CHDRS)
	@echo Making dependencies
	@$(CXX) $(CFLAGS) -I$(MROOT) \
		$(CSRCS) -MM | sed 's|\(.*\):|$(PWD)/\1 $(PWD)/\1r $(PWD)/\1d $(PWD)/\1p:|' > depend.mk
	@for dir in $(DEPDIR); do \
			if [ -r $(MROOT)/$${dir}/depend.mk ]; then \
			echo Depends on: $${dir}; \
			cat $(MROOT)/$${dir}/depend.mk >> depend.mk; \
			fi; \
	  done

-include $(MROOT)/mtl/config.mk
-include depend.mk
