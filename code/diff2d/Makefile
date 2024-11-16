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

VARIANTS = oned clps kokkos2d span iota sycl diy2d diy2e \
    seq f77
.PHONY: variants
variants: 
	@echo ${VARIANTS}
# tile tbb 
.PHONY: compile
info ::
	@echo "================ Program building"
	@echo "make compile : [ VARIANTS=... (default: ${VARIANTS}) ]"
compile :
	@compiler=$$( make --no-print-directory set_compiler ) \
	 && bindir=bin_$${compiler} \
	 && mkdir -p $${bindir} \
	 && variants="$$( make --no-print-directory variants_list )" \
	 && makelog="$$( make --no-print-directory make_log COMPILER=$${compiler} )" \
	 \
	 && for v in $$variants ; do \
	    ( echo && echo "Making variant $${v}" \
	       && cd $${v} \
	       && make --no-print-directory diff2d \
	            LIBDIR=$${libdir} ECHO=1 \
	            EXTRA_LIBS=${EXTRA_LIBS} \
	          2>&1 | tee ../$${bindir}/$${v}_build.txt \
	       && cp diff2d ../$${bindir}/$${v} \
	    ) \
	    ; done 2>&1 | tee $${makelog} \
	 && make --no-print-directory build_summary \
	      ECHO=${ECHO} MAKELOG=$${makelog}
.PHONY: set_compiler variants_list make_log build_summary
set_compiler :
	@if [ -z "${TACC_FAMILY_COMPILER}" ] ; then \
	    echo "IMPORTANT set TACC_FAMILY_COMPILER to your compiler name" && exit 1 \
	 ; fi \
	 && echo ${TACC_FAMILY_COMPILER}
variants_list :
	@if [ ! -z "${VARIANT}" ] ; then \
	      variants="${VARIANT}" \
	    ; else \
	      variants="${VARIANTS}" \
	    ; fi \
	 && echo $${variants} | tr ',' ' ' 
make_log :
	@mkdir -p logs \
	 && compiler=$${COMPILER} \
	 && echo logs/build-$${compiler:-compiler}.log
clean ::
	@rm -rf logs
build_summary :
	@compiler=${TACC_FAMILY_COMPILER} \
	 && bindir=bin_$${compiler} \
	 && echo "Generated:" \
	 && ls -l $$( for v in ${VARIANTS} ; do echo $${bindir}/$$v ; done ) || /bin/true \
	 && if [ "${ECHO}"  = "1" ] ; then echo "(see also ${MAKELOG})" ; fi

.PHONY: cmake
info ::
	@echo "make cmake [ VARIANTS=.... ] [ ECHO=1 for verbose compilation ]"
cmake :
	@compiler=$$( make --no-print-directory set_compiler ) \
	 && bindir=bin_$${compiler} \
	 && mkdir -p $${bindir} \
	 && variants="$$( make --no-print-directory variants_list )" \
	 && makelog="$$( make --no-print-directory make_log COMPILER=$${compiler} )" \
	 \
	 && for v in $${variants} ; do \
	      if [ -d "$${v}" ] ; then \
	        buildlog=bin/build_$${v}_$${compiler}.log \
	         && rm -rf build $${v}/CMakeLists.txt $${buildlog} \
	         && mkdir -p build && ln cmake/CMakeLists.txt $${v}/CMakeLists.txt \
	         && ( cd build  \
	             && cmake \
	              -D CMAKE_INSTALL_PREFIX=../$${bindir} \
	              $$( if [ ! -z "${TACC_CXX}" ] ; then echo "-D CMAKE_CXX_COMPILER=${TACC_CXX}" ; fi ) \
	              $$( if [ "$${v}" = "f77" ] ; then echo "-D USEFORTRAN=ON" ; fi ) \
	              -D VARIANT=$${v} \
	              -D CMAKE_VERBOSE_MAKEFILE=ON \
	              ../$${v} \
	             && make --no-print-directory V=1 \
	                  2>&1 | tee ../$${bindir}/$${v}_build.txt \
	             && make --no-print-directory install \
	          ) 2>&1 | tee $${buildlog} \
	      ; else echo "No such variant: <<$$v>>" && exit 1 \
	      ; fi \
	     ; done 2>&1 | tee $${makelog} \
	 && make --no-print-directory build_summary \
	      ECHO=${ECHO} MAKELOG=$${makelog}

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
.PHONY: clean_installation
info ::
	@echo "make clean_installation [VARIANT=] : clean for compiler [variant only]"
clean_installation : clean 
	@compiler=${TACC_FAMILY_COMPILER} \
	 && bindir=bin_$${compiler} \
	 && if [ -z "${VARIANT}" ] ; then \
	      echo " .. removing bin dir <<$${bindir}>>" \
	       && rm -rf $${bindir} \
	    ; else \
	      echo " .. removing variant <<$${bindir}/${VARIANT}>>" \
	       &&  rm -f $${bindir}/${VARIANT} \
	    ; fi
total_clean :: 
	@rm -rf bin_*
clean :: recursive_clean
	@echo " .. local clean"
	@rm -rf build
	@rm -rf r[0-9]*hs *.log
	@for v in ${VARIANTS} ; do rm -f $${v}/CMakeLists.txt ; done
noclean :
	@echo " .. code clean" && cd codes && make --no-print-directory clean 
MAKEFILESDIR := ../${MAKEFILESDIR}
MAKETRACESTRING := "${MAKETRACESTRING} .. .."
