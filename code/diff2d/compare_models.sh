#!/bin/bash
################################################################
####
#### Compare programming paradigms on one chip
####
################################################################

function usage () {
    echo
    echo "This generates commandlines:" && echo
    echo "[ srun -p QUEUE -t ${runtime} -N 1 -n 1 -A ACCOUNT "
    echo "      --cpu-bind=verbose,mask_cpu=MASK ]"
    echo "make run_scaling ## complicated script!"
    echo "    PROGRAM=../bin/VARIANT NSIZE=\${NSIZE} TACC_SYSTEM=\${CPU} CATEGORY=\${VARIANT}"
    echo "    [ GITADD=1 ] "
    echo && echo " where parameters can be set:" && echo
    echo "Usage: $0 "
    echo "    [ -c cpu ]"
    echo "    [ -n 123456 (default ${nsize}) ] [ -p 1123 (default=${cores}) ]"
    echo "    [ -g (git add ) ] [ -s (prepend srun) ]"
    echo "    [ -q queue (default: ${queue}) ] [ -t (trace) ] "
    echo "    c1,c2,c3 (from: ${allcodes} or \"all\" for all)"
}

cpu=cpu
nsize=25000
procs=112
gitadd=0
srun=
trace=
runtime=0:30:0
## TACC specific:
queue=spr
## Victor specific:
account=A-ccsc
while [ $# -gt 0 ] ; do
    if [ $1 = "-h" ] ; then 
	usage && exit 1
    elif [ $1 = "-c" ] ; then 
	shift && cpu=$1 && shift
    elif [ $1 = "-g" ] ; then 
	gitadd=1 && shift
    elif [ $1 = "-n" ] ; then
	shift && nsize=$1 && shift
    elif [ $1 = "-p" ] ; then
	shift && procs=$1 && shift
    elif [ $1 = "-q" ] ; then
	shift && queue=$1 && shift
    elif [ $1 = "-s" ] ; then 
	srun=1 && shift
    elif [ $1 = "-t" ] ; then
	 trace=1 && shift
    else
	break
    fi
done
codes=$1
if [ ${codes} = "all" ] ; then
    codes=oned,clps,span,iota,kokkos,sycl,
fi

echo "================ Testing codes: ${codes}"
echo " problem size $nsize"
echo " using $procs cores"
echo " cpu designation: $cpu"
echo " queue $queue"

for code in $( echo $codes | tr ',' ' ' ) ; do
    mask=$( python3 ../utils/maskgen.py ${procs} 1 )
    echo && echo "================ submit diff2d code=$code" && echo
    if [ ! -f bin/${code} ] ; then continue ; fi
    ( cd ${code} \
       && if [ "${srun}" = "1" ] ; then \
	     cmdline="srun -p $queue -t ${runtime} -N 1 -n 1 -A ${account} \
               --cpu-bind=verbose,mask_cpu=${mask}" \
	  ; else \
	      cmdline="" \
	  ; fi \
       && cmdline="$cmdline \
	    make run_scaling PROGRAM=../bin/${code} NSIZE=${nsize} GITADD=${gitadd} \
	      TACC_SYSTEM=${cpu} \
	      CATEGORY=${code} \
	      THREADSYSTEM=$( \
	        if [ \"${code}\" = \"sycl\" ] ; then echo dpcpp ; else echo omp ; fi ) \
	      $( if [ ! -z \"${trace}\" ] ; then echo "ECHO=1 D2D_OPTIONS=--trace" ; fi ) \
	      " \
       && echo $cmdline \
       && eval $cmdline \
	)
done
for code in $( echo $codes | tr ',' ' ' ) ; do
    ls -l ${code}/diff2d-scaling-*-${cpu}-${TACC_FAMILY_COMPILER}-${nsize}.*out
done

