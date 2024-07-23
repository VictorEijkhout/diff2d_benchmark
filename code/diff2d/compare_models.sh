#!/bin/bash
################################################################
####
#### Compare programming paradigms on one chip
####
################################################################

function usage () {
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
queue=spr
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
    codes=oned,clps,kokkos,span,sycl,diy2d
fi

echo "================ Testing codes: ${codes}"
echo " problem size $nsize"
echo " using $procs cores"
echo " cpu designation: $cpu"
echo " queue $queue"

for code in $( echo $codes | tr ',' ' ' ) ; do
    mask=$( python3 ../utils/maskgen.py ${procs} 1 )
    echo && echo "================ submit diff2d code=$code" && echo
    ( cd ${code} \
       && if [ "${srun}" = "1" ] ; then \
	     cmdline="srun -p $queue -t 0:30:0 -N 1 -n 1 -A A-ccsc \
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

