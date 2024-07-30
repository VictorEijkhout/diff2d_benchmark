#!/bin/bash

nsize=30000
if [ $# -gt 0 ] ; then
    nsize=$1
fi
python3 ../../scripts/multi_graphs_extract.py \
	omp/diff2d-scaling-omp-frontera${nsize}.runout:omp \
	clps/diff2d-scaling-clps-frontera${nsize}.runout:clps \
	rng/diff2d-scaling-range-frontera${nsize}.runout:rng \
	span/diff2d-scaling-span-frontera${nsize}.runout:span \
	Execute:4 Time:1
