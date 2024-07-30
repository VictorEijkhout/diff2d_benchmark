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
	@if [ ! -z "${VARIANT}" ] ; then \
	    echo "The option is VARIANTS, not VARIANT" && exit 1 ; fi
	@variants="$$( echo ${VARIANTS} | tr ',' ' ' )" \
	 && for v in $$variants ; do \
	    ( echo "Making variant $${v}" \
	       && cd $${v} \
	       && make --no-print-directory diff2d \
	       && cp diff2d ../bin/$${v} \
	    ) \
	    ; done 2>&1 | tee compile.log && echo "See: compile.log" \
	 && ls -l $$( for v in $${variants} ; do echo bin/$$v ; done )

.PHONY: cmake
info ::
	@echo "make cmake [ VARIANTS=.... ]"
cmake :
	@if [ ! -z "${VARIANT}" ] ; then \
	    echo "The option is VARIANTS, not VARIANT" && exit 1 ; fi
	@for v in ${VARIANTS} ; do \
	  if [ -d "$${v}" ] ; then \
	    rm -rf build $${v}/CMakeLists.txt \
	     && mkdir -p build && ln cmake/CMakeLists.txt $${v}/CMakeLists.txt \
	     && ( cd build  \
	         && cmake \
	          $$( if [ ! -z "${TACC_CXX}" ] ; then echo "-D CMAKE_CXX_COMPILER=${TACC_CXX}" ; fi ) \
	          -D CMAKE_VERBOSE_MAKEFILE=ON \
	          -D VARIANT=$${v} \
	          -D CMAKE_INSTALL_PREFIX=../bin \
	          ../$${v} \
	         && make V=1 && make install \
	      ) \
	  ; else echo "No such variant: <<$$v>>" && exit 1 \
	  ; fi \
	 ; done \
	 && echo "Generated:" \
	 && for v in ${VARIANTS} ; do ls -l bin/$$v ; done

##
## Do scaling runs over all VARIANTS
## This calls "make run_scaling" in each directory;
## that rule is in Make.runs
##
info ::
	@echo "================ Running"
	@echo "run_scaling : big run across all variants"
	@echo "    [ NSIZE= (default: ${NSIZE}) ]"
	@echo "    [ GITADD= (1 to rename as runout and add to repo)"
NSIZE = 40000
.PHONY: run_scaling
run_scaling :
	@if [ -z "${TACC_PROC_CODE}" ] ; then \
	    echo "Set variable: TACC_PROC_CODE" && exit 1 ; fi 
	@for v in ${VARIANTS} ; do \
	  ( cd $${v} \
	     && make --no-print-directory run_scaling PROGRAM=../bin/$${v} \
	          TACC_SYSTEM=${TACC_PROC_CODE} \
	          NDIV=10 NSIZE=${NSIZE} ECHO=1 GITADD=${GITADD} \
	  ) \
	 ; done 2>&1 | tee run.log && echo "See: run.log"

include Make.prof

##
## GPROF profiling
## here we loop over models,
## do profiling on each
##
info ::
	@echo "================ Profiling"
	@echo "run_gprof : gprof on variants"
	@echo "    [ NSIZE=123456 (default: ${NSIZE}) ]"
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
	@echo "================ Profiling"
	@echo "run_perf : perf on variants"
	@echo "    [ NSIZE=123456 (default: ${NSIZE}) ]"
	@echo "    [ NORUN= (nonzero for redisplay of perf output) ]"
NORUN =
.PHONY: run_perf
run_perf :
	@for v in ${VARIANTS} ; do \
	  ( cd $${v} \
	     && perfdata=perf_$${v}_${TACC_FAMILY_COMPILER}.out \
	     && if [ -z "${NORUN}" ] ; then \
	          rm -f diff2d && make --no-print-directory diff2d ECHO=${ECHO} \
	           && cmdline="perf record --output=$${perfdata} \
	                 --call-graph fp \
	                ./diff2d -m ${NSIZE} -n ${NSIZE} -i 5 -t" \
	           && echo "perf cmdline: $$cmdline" && eval $$cmdline \
	        ; fi \
	     && make --no-print-directory perf PERDATA=$${perfdata} \
	  ) \
	; done
clean ::
	@rm -f */perf_*.out */perf.data*


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
