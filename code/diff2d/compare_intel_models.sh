#!/bin/bash

if [ "${TACC_FAMILY_COMPILER" != "intel" ] ; then
    echo "ERROR This script is for intel"
    exit 1
fi

make clean
make cmake VARIANTS=oned,clps,span,iota,kokkos2d,sycl
./compare_models.sh -c spr -n 26000 -g -t all
