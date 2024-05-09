#!/bin/bash

for q in skx icx spr ; do
    echo && echo "================ submit to: $q" && echo
    export TACC_SYSTEM=$q
    srun -p $q -t 0:30:0 -N 1 -n 1 -A A-ccsc \
	 make run_scaling NSIZE=20000 GITADD=1
done
