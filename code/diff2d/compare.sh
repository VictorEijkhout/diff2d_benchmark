#!/bin/bash

function usage () {
    echo "Usage: $0 [ -c c1,c2,c3 ]"
}

codes=skx,icx,spr
while [ $# -gt 0 ] ; do
    if [ $1 = "-h" ] ; then 
	usage && exit 1
    elif [ $1 = "-c" ] ; then 
	shift && codes=$1 && shift
    else 
	echo "Undefined option: $1" && exit 1
    fi
done

for code in $( echo $codes | tr ',' ' ' ) ; do
    if [ $code = "csx" ] ; then 
	queue=small ; else queue=$code ; fi
    echo && echo "================ submit to queue=$queue, proc code=$code" && echo
    export TACC_SYSTEM=$code
    srun -p $queue -t 0:30:0 -N 1 -n 1 -A A-ccsc \
	 make run_scaling NSIZE=20000 GITADD=1
done
