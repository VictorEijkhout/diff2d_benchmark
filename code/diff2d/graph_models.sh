#!/bin/bash

function usage () {
    echo
    echo "This makes for all models on a given cpu"
    echo
    echo "$0 [ -h ] [ -n 123456 ] "
    echo "    [ -c cpu (default: ${cpu}) ] [ -g (gcc, otherwise intel) ]"
}

nsize=25000
compiler=intel
cpu=spr
test=0
while [ $# -gt 0 ] ; do
    if [ $1 = "-h" ] ; then
	usage && exit 0
    elif [ $1 = "-c" ] ; then
	shift && cpu=$1 && shift
    elif [ $1 = "-g" ] ; then
	compiler=gcc && shift
    elif [ $1 = "-n" ] ; then
	shift && nsize=$1 && shift
    elif [ $1 = "-t" ] ; then
	test=1 && shift
    elif [ $1 = "-y" ] ; then
	sycl=1 && shift
    else
	echo "Unknown option <<$1>>" && exit 1
    fi
done

plotdir=../../writing/plots
if [ ! -d "${plotdir}" ] ; then
    echo "ERROR can not find plot dir: ${plotdir}" && exit 1
fi
python3 ../../scripts/multi_graphs_extract.py \
	-p ${plotdir} \
	-n ${cpu}-models-${compiler}-${nsize} \
	$( if [ "${test}" = "1" ] ; then echo --test ; fi ) \
	$( for m in oned clps span iota kokkos2d sycl ; do \
	       file=${m}/diff2d-scaling-${m}-${cpu}-${compiler}-${nsize}.runout \
		   && if [ -f "${file}" ] ; then echo ${file}:${m} ; fi \
	       ; done ) \
	Execute:4 Time:1
## lines such as:
## 	clps/diff2d-scaling-clps-${cpu}-${compiler}-${nsize}.runout:clps
