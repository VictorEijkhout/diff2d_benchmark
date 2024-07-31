#!/bin/bash

function usage () {
    echo "$0 [ -h ]"
    echo "    [ -n 123456 (default: ${nsize}) ]"
    echo "    [ -c cpu (default: ${cpu}) ] [ -g (gcc, otherwise intel) ]"
    echo "    [ -t (for test) ] [ -v (for verbose) ]"
}

nsize=25000
compiler=intel
cpu=spr
test=0
verbose=0
while [ $# -gt 0 ] ; do
    if [ $1 = "-h" ] ; then
	usage && exit 0
    elif [ $1 = "-g" ] ; then
	compiler=gcc && shift
    elif [ $1 = "-n" ] ; then
	shift && nsize=$1 && shift
    elif [ $1 = "-p" ] ; then
	shift && cpu=$1 && shift
    elif [ $1 = "-y" ] ; then
	sycl=1 && shift
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
python3 ../../scripts/ratio_graphs_extract.py \
	--path ${plotdir} \
	--name ${cpu}-ratios-${compiler} \
	$( if [ "${verbose}" = "1" ] ; then echo --verbose ; fi ) \
	$( if [ "${test}" = "1" ] ; then echo --test ; fi ) \
	$( for m in oned clps span iota kokkos sycl ; do \
	       file=${m}/diff2d-scaling-${m}-${cpu}-${compiler}-${nsize}.runout \
		   && if [ -f "${file}" ] ; then echo ${file}:${m} ; fi \
	       ; done ) \
	Execute:4 Time:1
