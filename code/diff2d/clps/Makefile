################################################################
################################################################
####
#### This makefile is part of the source of 
#### "Parallel Computing for Science and Engineering"
#### by Victor Eijkhout, copyright 2013-2024
####
#### makefile for diff2d in OpenMP/C++, omp collapse
####
################################################################
################################################################

info ::
	@echo "================================================================"
	@echo "    OpenMP codes for diff2d in C++"
	@echo "================================================================"
	@echo && echo "Available rules:" && echo

APPS = diff2d

MODE = omp
LANGUAGE = CXX
CATEGORY = clps
MODULES = cxxopts mdspan

##
## compilation
##
include ../Make.depend
include ../../../makefiles/Make.compile
CPPSTANDARD = 23

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
