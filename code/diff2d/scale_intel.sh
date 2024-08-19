#!/bin/bash

function usage () {
    echo "Usage: $0 "
    echo "    [ -g (gnu, otherwise intel) ]"
    echo "    [ -n 12345 (default: ${nsize}) ] "
    echo "    [ -v variant (default: ${variant}) ]"
    echo "    [ -x extension (default: ${extension}) ]"
}

nsize=20000
variant=oned
compiler=intel
extension=sp
while [ $# -gt 0 ] ; do
    if [ $1 = "-h" ] ; then
	usage && exit 0
    elif [ $1 = "-g" ] ; then 
	shift && compiler=gcc
    elif [ $1 = "-n" ] ; then 
	shift && nsize=$1 && shift
    elif [ $1 = "-v" ] ; then 
	shift && variant=$1 && shift
    elif [ $1 = "-x" ] ; then 
	shift && extension=$1 && shift
    else
	echo "Unknown option: <<$1>>" && exit 1
    fi
done

for queue in skx clx icx spr grx ; do 
    file=${variant}/diff2d-scaling-${variant}-${queue}-${compiler}-${nsize}.runout
    if [ ! -f ${file} ] ; then
	echo "ERROR could not find file ${file}" && exit 1
    fi
    if [ "${extension}" = "bw" ] ; then
	python3 ../../scripts/multi_graphs_extract.py \
	    -p ../../writing/plots -n ${compiler}-${variant}-${queue}-${nsize}-bw \
	    ${file}:${queue} \
	    Execute:4 omp-BW:3
    else
	python3 ../../scripts/multi_graphs_extract.py \
	    -p ../../writing/plots -n ${compiler}-${variant}-${queue}-${nsize} \
	    ${file}:${queue} \
	    Execute:4 Time:1
    fi
done
