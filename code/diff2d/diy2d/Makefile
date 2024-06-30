################################################################
################################################################
####
#### This makefile is part of the source of 
#### "Parallel Computing for Science and Engineering"
#### by Victor Eijkhout, copyright 2013-2024
####
#### makefile for OpenMP in C++ using mdspan and cartesian
####
################################################################
################################################################

info ::
	@echo "================================================================"
	@echo "    OpenMP codes for diff2d in C++"
	@echo "================================================================"
	@echo && echo "Available rules:" && echo

APPS = diff2d
WRONGS = 

##
## use of c++20 range algorithms in OpenMP
## note, we need the ranges library regardless
## for the cartesian product and such
##
RANGE = 0

MODE = omp
LANGUAGE = CXX
CATEGORY = diy2d
MODULES := cxxopts mdspan 

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
EXTRA_INCLUDES = \
    -I../lib \
    -I${TACC_MDSPAN_INC} \
    $$( pkg-config --cflags cxxopts ) 

##
## profiling
##
include ../Make.prof

##
## Runs
##

RUNS = 
include ../../../makefiles/Make.omp_run
include ../Make.runs
.PHONY: ${RUNS}

##
## Clean
##

include  ../../../makefiles/Make.clean
MAKEFILESDIR := ../../${MAKEFILESDIR}
MAKETRACESTRING := "${MAKETRACESTRING} .. .. .."
clean :: codeclean
	@/bin/rm -rf ${APPS} ${WRONGS}
