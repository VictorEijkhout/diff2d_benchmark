################################################################
################################################################
####
#### This makefile is part of the source of 
#### "Parallel Computing for Science and Engineering"
#### by Victor Eijkhout, copyright 2013-2024
####
#### This is the makefile for the diff2d using MPL
####
################################################################
################################################################

info ::
	@echo "================================================================"
	@echo "    MPL/MPI codes for diff2d in C++"
	@echo "================================================================"
include ../../../makefiles/Make.info

APPS = diff2d
PROGRAM = diff2d
WRONGS = 

MODE = mpl
LANGUAGE = CXX
CATEGORY = dist
MODULES = cxxopts mdspan mpl

## NEEDS TOO MUCH CUSTOMMIZATION include ../Make.depend
##
## dependencies
##
EXTRA_OBJS = ../lib/${CATEGORY}.o ../lib/seq.o ../lib/base.o ../lib/options.o
EXTRA_OBJS_diff2d = ../lib/${CATEGORY}.o ../lib/seq.o ../lib/base.o ../lib/options.o
diff2d.o : ../main.cpp ../lib/options.hpp ../lib/${CATEGORY}.hpp ../lib/base.hpp
diff2d.o : modules ${EXTRA_OBJS}
../lib/${CATEGORY}.o : ../lib/${CATEGORY}.cpp ../lib/${CATEGORY}.hpp
	@cd ../lib \
	 && make --no-print-directory ECHO=${ECHO} CATEGORY=${CATEGORY} \
	      EXTRA_INCLUDES="${EXTRA_INCLUDES}" \
	      EXTRA_OPTIONS="${EXTRA_OPTIONS}" \
	      ${CATEGORY}.o
../lib/seq.o : ../lib/seq.cpp ../lib/seq.hpp
	@cd ../lib \
	 && make --no-print-directory ECHO=${ECHO} CATEGORY=seq \
	      EXTRA_INCLUDES="${EXTRA_INCLUDES}" \
	      EXTRA_OPTIONS="${EXTRA_OPTIONS}" \
	      seq.o
../lib/base.o : ../lib/base.cpp ../lib/base.hpp
	@cd ../lib \
	 && make --no-print-directory ECHO=${ECHO} CATEGORY=${CATEGORY} \
	      EXTRA_INCLUDES="${EXTRA_INCLUDES}" \
	      EXTRA_OPTIONS="${EXTRA_OPTIONS}" \
	      base.o
../lib/options.o : ../lib/options.cpp ../lib/options.hpp
	@cd ../lib \
	 && make --no-print-directory ECHO=${ECHO} CATEGORY=${CATEGORY} \
	      EXTRA_INCLUDES="${EXTRA_INCLUDES}" \
	      EXTRA_OPTIONS="${EXTRA_OPTIONS}" \
	      options.o
##
## compilation
##
EXTRA_INCLUDES = -I${TACC_MDSPAN_INC} \
    $$( pkg-config --cflags cxxopts ) \
    -I${TACC_MPL_INC}

include ../../../makefiles/Make.compile
CPPSTANDARD = 23

RUNS =
info ::
	@echo "make run_XXX (in ${RUNS})"
include ../Make.runs
include ../../../makefiles/Make.mpi_run

include ../../../makefiles/Make.clean
CATCHER = ../../../teaching/course-support/catchout

clean ::
	@rm -f *.dat
