################################################################
################################################################
####
#### This makefile is part of the source of 
#### "Parallel Computing for Science and Engineering"
#### by Victor Eijkhout, copyright 2013-2024
####
#### makefile for diff2d in OpenMP/C++ using SYCL
####
################################################################
################################################################

info ::
	@echo "================================================================"
	@echo "    OpenMP codes for diff2d in C++ using SYCL"
	@echo "================================================================"
	@echo && echo "Available rules:" && echo

APPS = diff2d
WRONGS = 

MODE = sycl
LANGUAGE = CXX
CATEGORY = sycl
MODULES = cxxopts 

include ../Make.depend

##
## compilation
##
include ../../../makefiles/Make.compile
CPPSTANDARD = 23
#
# sycl build for explicit target:
# https://www.intel.com/content/www/us/en/docs/dpcpp-cpp-compiler/developer-guide-reference/2024-2/ahead-of-time-compilation.html
#
EXTRA_OPTIONS = -fsycl-targets=spir64_x86_64
EXTRA_LINK_FLAGS = -fsycl-targets=spir64_x86_64 -Xsycl-target-backend  "-march=avx2"

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
