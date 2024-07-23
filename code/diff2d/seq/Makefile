################################################################
################################################################
####
#### This makefile is part of the source of 
#### "Parallel Computing for Science and Engineering"
#### by Victor Eijkhout, copyright 2013-2024
####
#### makefile for sequential C++ programs
####
################################################################
################################################################

info ::
	@echo "================================================================"
	@echo "    Sequential codes for diff2d in C++"
	@echo "================================================================"
	@echo && echo "Available rules:" && echo


APPS = diff2d

MODE = seq
LANGUAGE = CXX
CATEGORY = seq
MODULES = cxxopts mdspan 

##
## compilation
##
include ../Make.depend
include ../../../makefiles/Make.compile
ifeq "${TACC_FAMILY_COMPILER}" "intel"
  CPPSTANDARD = 2b
else
  CPPSTANDARD = 23
endif

##
## Runs
##

RUNS = 
RUNS += run_default
include ../Make.runs
include ../../../makefiles/Make.omp_run

.PHONY: ${RUNS}

##
## cmake
##

CMAKEFILESDIR = ../../../makefiles
CMAKE_CXX_FLAGS = \
include ../../../makefiles/Make.cmake

##
## Clean
##

include  ../../../makefiles/Make.clean
MAKEFILESDIR := ../../${MAKEFILESDIR}
MAKETRACESTRING := "${MAKETRACESTRING} .. .. .."
clean :: codeclean
	@/bin/rm -rf ${APPS} ${WRONGS}
