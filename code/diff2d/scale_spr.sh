#!/bin/bash

function usage () {
    echo "$0 [ -h ] [ -g (gcc, otherwise intel) ] [ -n 123456 ] [ -y (include sycl) ]"
}

nsize=25000
compiler=intel
sycl=0
while [ $# -gt 0 ] ; do
    if [ $1 = "-h" ] ; then
	usage && exit 0
    elif [ $1 = "-g" ] ; then
	compiler=gcc && shift
    elif [ $1 = "-n" ] ; then
	shift && nsize=$1 && shift
    elif [ $1 = "-y" ] ; then
	sycl=1 && shift
    else
	echo "Unknown option <<$1>>" && exit 1
    fi
done

python3 ../../scripts/multi_graphs_extract.py \
	-p ../../plots -n spr-models-${compiler} \
	oned/diff2d-scaling-oned-spr-${compiler}-${nsize}.runout:oned \
	clps/diff2d-scaling-clps-spr-${compiler}-${nsize}.runout:clps \
	span/diff2d-scaling-span-spr-${compiler}-${nsize}.runout:span \
	kokkos/diff2d-scaling-kokkos-spr-${compiler}-${nsize}.runout:kokkos \
	$( yfile=sycl/diff2d-scaling-sycl-spr-${compiler}-${nsize}.runout && if [ "${sycl}" = "1" ] ; then echo ${yfile}:sycl ; fi ) \
	Execute:4 Time:1
