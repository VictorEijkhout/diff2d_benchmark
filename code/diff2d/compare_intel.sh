#!/bin/bash
################################################################
####
#### Compare programming paradigms on one chip
####
#### default values are very TACC centric, mostly stampede3 procs
####
################################################################

commandline=$*

function usage () {
    echo && echo "srun a code on multiple procs" && echo
    echo "Usage: $0 "
    echo "    [ -c c1,c2,c3 (default: ${codes}) ] }"
    echo "    [ -p p1,p2,p3 (default: ${procs}) ]"
    echo "        ( known procs: ${allprocs} )"
    echo "    [ -n 12345 (default: ${nsize} ]"
    echo "    [ -g (git add ) ] [ -t (trace) ]"
}

codes=oned
procs=skx,icx,spr,grx
allprocs=clx,${procs}

# frontera
cpus_clx=56

# stampede3
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
    elif [ $1 = "-p" ] ; then
	shift && procs=$1 && shift
    elif [ $1 = "-q" ] ; then
	echo ">>>> OBSOLETE OPTION: use -p for processors <<<<" && exit 1
    elif [ $1 = "-t" ] ; then
	 trace=1 && shift
    else
	echo "Undefined option: $1" && exit 1
    fi
done

echo
echo "Runnning codes: <<${codes}>> on procs: <<${procs}>>"
echo

( echo && echo "Test: $commandline" && echo ) | tee ${logfile}
for code in $( echo $codes | tr ',' ' ' ) ; do
    for proc in $( echo $procs | tr ',' ' ' ) ; do
	if [ -d "${code}" ] ; then
	    eval cpus=\${cpus_${proc}}
	    case ${proc} in \
		( "grx" ) queue=systest-grx;; \
		( "clx" ) queue=small;; \
		( * ) queue=${proc} ;;
		esac
	    mask=$( python3 ../utils/maskgen.py ${cpus} 1 )
	    echo && echo "================ submit to queue=$queue, code=$code cores=$cpus" && echo
	    pushd "${code}"
	    cmdline="srun -p $queue -t 0:30:0 -N 1 -n 1 -A A-ccsc \
               --cpu-bind=verbose,mask_cpu=${mask}"
	    cmdline="$cmdline \
	        make run_scaling PROGRAM=../bin/${code} NSIZE=${nsize} \
	          TACC_SYSTEM=$proc \
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

if [ "${gitadd}" = "1" ] ; then
    for code in ${codes} ; do
	ls -l ${code}/*${nsize}*.runout
    done
fi \
    2>&1 | tee -a ${logfile}

echo
echo "See: ${logfile}"
echo

