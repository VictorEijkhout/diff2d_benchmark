#!/bin/bash
################################################################
####
#### Compare programming paradigms on one chip
####
#### default values are very TACC centric, mostly stampede3 queues
####
################################################################

commandline=$*

function usage () {
    echo && echo "srun a code on multiple queues" && echo
    echo "Usage: $0 "
    echo "    [ -c c1,c2,c3 (default: ${codes}) ] }"
    echo "    [ -q q1,q2,q3 (default: ${queues}) ]"
    echo "    [ -n 12345 (default: ${nsize} ]"
    echo "    [ -g (git add ) ] [ -t (trace) ]"
}

codes=oned
queues=skx,icx,spr,grx
cpus_skx=48
cpus_icx=80
cpus_spr=112
cpus_grx=240
logfile=intel.log

nsize=25000
gitadd=0
##srun=1
trace=

while [ $# -gt 0 ] ; do
    if [ $1 = "-h" ] ; then 
	usage && exit 1
    elif [ $1 = "-c" ] ; then 
	shift && codes=$1 && shift
    elif [ $1 = "-g" ] ; then 
	gitadd=1 && shift
    elif [ $1 = "-n" ] ; then
	shift && nsize=$1 && shift
    elif [ $1 = "-q" ] ; then
	shift && queues=$1 && shift
    # elif [ $1 = "-s" ] ; then 
    # 	srun=1 && shift
    elif [ $1 = "-t" ] ; then
	 trace=1 && shift
    else
	echo "Undefined option: $1" && exit 1
    fi
done

echo
echo "Runnning codes: <<${codes}>> on queues: <<${queues}>>"
echo

( echo && echo "Test: $commandline" && echo ) | tee ${logfile}
for code in $( echo $codes | tr ',' ' ' ) ; do
    for queue in $( echo $queues | tr ',' ' ' ) ; do
	if [ -d "${code}" ] ; then
	    eval cpus=\${cpus_${queue}}
	    if [ ${queue} = "grx" ] ; then queue=systest-grx; fi
	    mask=$( python3 ../utils/maskgen.py ${cpus} 1 )
	    echo && echo "================ submit to queue=$queue, code=$code cores=$cpus" && echo
	    pushd "${code}"
	    cmdline="srun -p $queue -t 0:30:0 -N 1 -n 1 -A A-ccsc \
               --cpu-bind=verbose,mask_cpu=${mask}"
	    cmdline="$cmdline \
	        make run_scaling PROGRAM=../bin/${code} NSIZE=${nsize} \
	          TACC_SYSTEM=$queue \
	          THREADSYSTEM=$( \
	            if [ \"${code}\" = \"sycl\" ] ; then echo dpcpp ; else echo omp ; fi ) \
	 	  $( if [ ! -z \"${trace}\" ] ; then echo "ECHO=1 D2D_OPTIONS=--trace" ; fi ) \
		  GITADD=${gitadd} \
	          "
	    echo $cmdline
	    eval $cmdline
	    popd
	else
	    echo && echo "ERROR Unknown code: <<${code}>" && echo && exit 1
	fi
    done
done \
    2>&1 | tee -a ${logfile}
echo
echo "See: ${logfile}"
echo

