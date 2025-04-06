#!/bin/bash

if [ "${TACC_FAMILY_COMPILER}" != "gcc" ] ; then
    echo "ERROR This script is for gcc"
    exit 1
fi

function usage () {
    echo "Usage: $0 [ -h ] [ --make : recompile ]"
    echo " .. recompiles and runs models on SPR"
    echo " .. can be used in/out of idev"
}

compile=0
while [ $# -gt 0 ] ; do
    if [ "$1" = "-h" ] ; then
	usage
	exit 0
    elif [ "$1" == "--make" ] ; then
	compile=1
    else
	echo "Unrecognized option <<$1>>"
	exit 1
    fi
done

if [ "$compile" = "1" ] ; then
    make clean
    rm -rf bin_gcc
    make cmake VARIANTS=oned,clps,span,iota,kokkos2d
fi

./compare_models.sh \
    -c spr -n 26000 -g -t \
    $( if [ -z "${SLURM_NPROCS}" ] ; then echo "-s" ; fi ) \
    all
