#!/bin/bash

function usage () {
    echo "Usage: $0 [ -n 12345 (default: ${nsize}) ] [ -v variant (default: ${variant}) ]"
}

nsize=20000
variant=oned
while [ $# -gt 0 ] ; do
    if [ $1 = "-h" ] ; then
	usage && exit 0
    elif [ $1 = "-n" ] ; then 
	shift && nsize=$1 && shift
    elif [ $1 = "-v" ] ; then 
	shift && variant=$1 && shift
    else
	echo "Unknown option: <<$1>>" && exit 1
    fi
done

for queue in skx icx spr grx ; do 
    file=${variant}/diff2d-scaling-${variant}-${queue}-intel-${nsize}.runout
    if [ ! -f ${file} ] ; then
	echo "ERROR could not find file ${file}" && exit 1
    fi
    python3 ../../scripts/multi_graphs_extract.py \
	    -p ../../plots -n intel-${variant}-${queue}-${nsize} \
	    ${file}:${queue} \
	    Execute:4 Time:1
done
