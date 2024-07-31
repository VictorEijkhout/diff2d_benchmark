#!/bin/bash

function usage () {
    echo "$0 [ -h ]"
    echo "    [ -n 123456 (default: ${nsize}) ]"
    echo "    [ -c cpu (default: ${cpu}) ]"
    echo "    [ -t (for test) ] [ -v (for verbose) ]"
}

nsize=25000
cpu=spr
test=0
verbose=0
while [ $# -gt 0 ] ; do
    if [ $1 = "-h" ] ; then
	usage && exit 0
    elif [ $1 = "-n" ] ; then
	shift && nsize=$1 && shift
    elif [ $1 = "-p" ] ; then
	shift && cpu=$1 && shift
    elif [ $1 = "-t" ] ; then
	test=1 && shift
    elif [ $1 = "-v" ] ; then
	verbose=1 && shift
    else
	echo "Unknown option <<$1>>" && exit 1
    fi
done

plotdir=../../writing/plots
if [ ! -d "${plotdir}" ] ; then
    echo "ERROR can not find plot dir: ${plotdir}" && exit 1
fi
for m in oned clps span kokkos ; do \
    ifile=${m}/diff2d-scaling-${m}-${cpu}-intel-${nsize}.runout
    gfile=${m}/diff2d-scaling-${m}-${cpu}-gcc-${nsize}.runout
    if [ ! -f ${ifile} -o ! -f ${gfile} ] ; then
	echo "Not both of ${ifile} / ${gfile}"
	continue
    fi
    python3 ../../scripts/ratio_graphs_extract.py \
	    --path ${plotdir} \
	    --name ${m}-ratios-iccgcc \
	    $( if [ "${verbose}" = "1" ] ; then echo --verbose ; fi ) \
	    $( if [ "${test}" = "1" ] ; then echo --test ; fi ) \
	    ${ifile}:intel ${gfile}:gnu \
	    Execute:4 Time:1
done
