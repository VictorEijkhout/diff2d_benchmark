#!/bin/bash
################################################################
####
#### Compare programming paradigms
####
################################################################

function usage () {
    echo "Usage: $0 [ -c c1,c2,c3 (default: ${codes}) ] [ -t (trace) ]"
}

queue=spr
trace=
codes=span,kokkos,sycl
while [ $# -gt 0 ] ; do
    if [ $1 = "-h" ] ; then 
	usage && exit 1
    elif [ $1 = "-c" ] ; then 
	shift && codes=$1 && shift
    elif [ $1 = "-t" ] then
	 trace=1 && shift
    else 
	echo "Undefined option: $1" && exit 1
    fi
done

for code in $( echo $codes | tr ',' ' ' ) ; do
    echo && echo "================ submit to queue=$queue, proc code=$code" && echo
    export TACC_SYSTEM=spr
    ( cd ${code} \
       && srun -p $queue -t 0:30:0 -N 1 -n 1 -A A-ccsc \
	    make run_scaling NSIZE=20000 GITADD=1 \
	      THREADSYSTEM=$( \
	        if [ "${code}" = "sycl" ] ; then echo dpcpp ; else echo omp ; fi ) \
	      $( if [ ! -z "${trace}" ] ; then echo "ECHO=1 D2D_OPTIONS=--trace" ; fi ) \
	)
done
