################################################################
################################################################
####
#### This makefile is part of the source of 
#### "Parallel Computing for Science and Engineering"
#### by Victor Eijkhout, copyright 2013-2024
####
#### makefile for diff2d in OpenMP/C++ usig range algorithms
####
################################################################
################################################################

info ::
	@echo "================================================================"
	@echo "    OpenMP codes for diff2d in C++ using range algorithms"
	@echo "================================================================"
	@echo && echo "Available rules:" && echo

APPS = diff2d
WRONGS = 

MODE = omp
LANGUAGE = CXX
CATEGORY = range
MODULES = cxxopts mdspan

include ../Make.depend

##
## compilation
##
include ../../../makefiles/Make.compile
ifeq "${INTEL23}" "23"
  CPPSTANDARD = 2b
else
  CPPSTANDARD = 23
endif
EXTRA_INCLUDES := -I${TACC_MDSPAN_INC} \
    $$( pkg-config --cflags tbb ) \
    $$( pkg-config --cflags cxxopts ) 
EXTRA_LIBS = \
    $$( pkg-config --libs tbb ) 

##
## Runs
##

RUNS = 
RUNS += run_default
include ../Make.runs
include ../../../makefiles/Make.omp_run

.PHONY: ${RUNS}

##
## Clean
##

include  ../../../makefiles/Make.clean
MAKEFILESDIR := ../../${MAKEFILESDIR}
MAKETRACESTRING := "${MAKETRACESTRING} .. .. .."
clean :: codeclean
	@/bin/rm -rf ${APPS} ${WRONGS}
