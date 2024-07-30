#!/bin/bash

nsize=30000
if [ $# -gt 0 ] ; then
    nsize=$1
fi
python3 ../../scripts/multi_graphs_extract.py \
	mpl/diff2d-scaling-frontera-N4-40000.runout:mpl \
	Procs:1 Time:1
