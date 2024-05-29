################################################################
################################################################
####
#### This makefile is part of the source of 
#### "Parallel Computing for Science and Engineering"
#### by Victor Eijkhout, copyright 2013-2024
####
################################################################
################################################################

.PHONY: info
info ::

# list cxx last because that's sequential 
# and we're not terribly interested
VARIANTS = oned kokkos range span sycl \
    seq
# clps tile tbb 
.PHONY: compile
info ::
	@echo "make compile : [ VARIANTS=... (default: ${VARIANTS}) ]"
compile :
	@for v in ${VARIANTS} ; do \
	    ( echo "Making variant $${v}" \
	       && cd $${v} \
	       && make --no-print-directory diff2d \
	    ) \
	; done 2>&1 | tee compile.log && echo "See: compile.log"
.PHONY: run_scaling
info ::
	@echo "run_scaling : big run across all variants"
	@echo "    [ NSIZE= (default: ${NSIZE}) ]"
	@echo "    [ GITADD= (1 to rename as runout and add to repo)"
NSIZE = 40000
run_scaling :
	@if [ -z "${TACC_PROC_CODE}" ] ; then \
	    echo "Set variable: TACC_PROC_CODE" && exit 1 ; fi 
	@for v in ${VARIANTS} ; do \
	  ( cd $${v} \
	     && make --no-print-directory run_scaling \
	          TACC_SYSTEM=${TACC_PROC_CODE} \
	          NDIV=10 NSIZE=${NSIZE} ECHO=1 GITADD=${GITADD} \
	  ) \
	 ; done 2>&1 | tee run.log && echo "See: run.log"

include ../../makefiles/Make.clean
clean :: recursive_clean
MAKEFILESDIR := ../${MAKEFILESDIR}
MAKETRACESTRING := "${MAKETRACESTRING} .. .."
