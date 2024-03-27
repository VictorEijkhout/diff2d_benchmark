################################################################
################################################################
####
#### This makefile is part of the source of 
#### "Parallel Computing for Science and Engineering"
#### by Victor Eijkhout, copyright 2013-2023
####
#### makefile for diff2d in OpenMP/C++, omp with tile parameter
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
CATEGORY = tile
MODULES = cxxopts fmtlib mdspan

##
## dependencies
##
EXTRA_OBJS = linalg.o
diff2d.o : ../main.cpp ../options.cpp ../linalg.hpp
linalg.o : ../linalg.hpp linalg.cpp
## reuse the omp main unchanged
diff2d.cpp : ../omp/diff2d.cpp
	@echo " .. copying $<"
	@cp $< .
## we use the omp implementation, adding a `tile' parameter
linalg.cpp : ../omp/linalg.cpp
	@echo " .. copying $<"
	@cat $< \
	| sed -e '/parallel for/s/$$/ tile(30,30)/' \
	> $@
clean ::
	@rm -f diff2d.cpp linalg.cpp

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
    $$( pkg-config --cflags cxxopts ) \
    $$( pkg-config --cflags fmt )
EXTRA_LIBS = \
    $$( pkg-config --libs fmt )
RPATH = ${TACC_FMTLIB_LIB}

##
## profiling
##
include ../Make.prof

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
