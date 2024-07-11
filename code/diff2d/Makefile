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
VARIANTS = oned clps kokkos span iota sycl diy2d diy2e \
    seq
.PHONY: variants
variants: 
	@echo ${VARIANTS}
# clps tile tbb 
.PHONY: compile
info ::
	@echo "================ Program building"
	@echo "make compile : [ VARIANTS=... (default: ${VARIANTS}) ]"
compile :
	@variants="$$( echo ${VARIANTS} | tr ',' ' ' )" \
	 && for v in $$variants ; do \
	    ( echo "Making variant $${v}" \
	       && cd $${v} \
	       && make --no-print-directory diff2d \
	    ) \
	    ; done 2>&1 | tee compile.log && echo "See: compile.log" \
	 && ls -l $$( for v in $${variants} ; do echo $$v/diff2d ; done )

.PHONY: cmake
info ::
	@echo "make cmake [ VARIANTS=.... ]"
cmake :
	@for v in ${VARIANTS} ; do \
	  if [ -d "$${v}" ] ; then \
	    rm -rf build $${v}/CMakeLists.txt \
	     && mkdir -p build && ln cmake/CMakeLists.txt $${v}/ \
	     && ( cd build  \
	         && cmake \
	          -D CMAKE_CXX_COMPILER=${TACC_CXX} \
	          -D CMAKE_VERBOSE_MAKEFILE=ON \
	          -D VARIANT=$${v} \
	          -D CMAKE_INSTALL_PREFIX=../codes \
	          ../$${v} \
	         && make V=1 && make install \
	      ) \
	  ; else echo "No such variant: <<$$v>>" && exit 1 \
	  ; fi \
	; done

##
## Scaling
## here we loop over models,
## do scaling on each
##

.PHONY: run_scaling
info ::
	@echo "================ Running"
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

##
## GPROF profiling
## here we loop over models,
## do profiling on each
##
info ::
	echo "================ Profiling"
	echo "run_gprof : gprof on variants"
	echo "    [ NSIZE=123456 (default: ${NSIZE}) ]"
.PHONY: run_gprof
run_gprof :
	@for v in ${VARIANTS} ; do \
	  ( cd $${v} \
	     && rm -f diff2d && make --no-print-directory diff2d PROF=1 ECHO=${ECHO} \
	     && ./diff2d -m ${NSIZE} -n ${NSIZE} -i 5 -t \
	     && make --no-print-directory profile \
	     && echo "see $${v}/diff2d.prof" \
	  ) \
	; done

##
## PERF profiling
## here we loop over models,
## do profiling on each
##
info ::
	echo "================ Profiling"
	echo "run_perf : perf on variants"
	echo "    [ NSIZE=123456 (default: ${NSIZE}) ]"
.PHONY: run_perf
run_perf :
	@for v in ${VARIANTS} ; do \
	  ( cd $${v} \
	     && rm -f diff2d && make --no-print-directory diff2d ECHO=${ECHO} \
	     && perf record \
	            ./diff2d -m ${NSIZE} -n ${NSIZE} -i 5 -t \
	     && make --no-print-directory perf \
	  ) \
	; done
##           --call-graph fp --freq=max

##
## Clean
##

include ../../makefiles/Make.clean
clean :: recursive_clean
	@echo " .. local clean"
	@rm -rf build
	@rm -rf r[0-9]*hs *.log
	@for v in ${VARIANTS} ; do rm -f $${v}/CMakeLists.txt ; done
noclean :
	@echo " .. code clean" && cd codes && make --no-print-directory clean 
MAKEFILESDIR := ../${MAKEFILESDIR}
MAKETRACESTRING := "${MAKETRACESTRING} .. .."
