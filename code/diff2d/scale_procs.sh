#!/bin/bash

nsize=40000
if [ $# -gt 0 ] ; then
    nsize=$1
fi
for arch in skx clx icx ; do
    infile=omp/diff2d-scaling-oned-${arch}-intel-40000.runout
    ## time
    basename=${infile}
    basename=${basename%%.runout}
    basename=${basename##*/}
    outfile=../../plots/${basename}.csv
    python3 ../../scripts/graphs_extract.py \
	    ${infile} \
	    ${outfile} \
	    Execute:4 Time:1
    git add $outfile 
    ## bandwidth
    basename=${infile}
    basename=${basename%%.runout}
    basename=${basename##*/}
    outfile=../../plots/${basename}-bw.csv
    python3 ../../scripts/graphs_extract.py \
	    ${infile} \
	    ${outfile} \
	    Execute:4 omp-BW:3
    git add $outfile 
done
