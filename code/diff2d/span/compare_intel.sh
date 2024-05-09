#!/bin/bash

for q in skx icx spr ; do
    export TACC_SYSTEM=$q
    srun -p $q -t 0:10:0 \
	 make run_scaling NSIZE=40000
done
