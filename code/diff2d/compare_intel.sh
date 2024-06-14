#!/bin/bash
################################################################
####
#### Compare programming paradigms on one chip
####
################################################################

function usage () {
    echo "Usage: $0 [ -c c1,c2,c3 (default: ${codes}) ] }"
    echo "    [ -g (git add ) ] [ -s (srun; default: ${srun}) ]"
    echo "    [ -q q1,q2,q3 (default: ${queues}) ] [ -t (trace) ] "
}

gitadd=0
srun=1
trace=
codes=oned
queues=skx,icx,spr,systest-grx
logfile=intel.log
while [ $# -gt 0 ] ; do
    if [ $1 = "-h" ] ; then 
	usage && exit 1
    elif [ $1 = "-c" ] ; then 
	shift && codes=$1 && shift
    elif [ $1 = "-g" ] ; then 
	gitadd=1 && shift
    elif [ $1 = "-q" ] ; then
	shift && queues=$1 && shift
    elif [ $1 = "-s" ] ; then 
	srun=1 && shift
    elif [ $1 = "-t" ] ; then
	 trace=1 && shift
    else
	echo "Undefined option: $1" && exit 1
    fi
done

for code in $( echo $codes | tr ',' ' ' ) ; do
    for queue in $( echo $queues | tr ',' ' ' ) ; do
	cpus=112
	mask=$( python3 ../utils/maskgen.py ${cpus} 1 )
	echo && echo "================ submit to queue=$queue, proc code=$code" && echo
	( cd ${code} \
	      && if [ "${srun}" = "1" ] ; then \
	      cmdline="srun -p $queue -t 0:30:0 -N 1 -n 1 -A A-ccsc \
               --cpu-bind=verbose,mask_cpu=${mask}" \
	      ; else \
	      cmdline="" \
	      ; fi \
	      && cmdline="$cmdline \
	    make run_scaling NSIZE=25000 GITADD=${gitadd} \
	      TACC_SYSTEM=$queue \
	      THREADSYSTEM=$( \
	        if [ \"${code}\" = \"sycl\" ] ; then echo dpcpp ; else echo omp ; fi ) \
	      $( if [ ! -z \"${trace}\" ] ; then echo "ECHO=1 D2D_OPTIONS=--trace" ; fi ) \
	      " \
	      && echo $cmdline \
	      && eval $cmdline \
	    )
    done
done \
    2>&1 | tee ${logfile} \
    && echo && echo "See: ${logfile}" && echo
