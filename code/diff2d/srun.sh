#!/bin/bash

function usage () {
    echo "$0 [ -h ]"
    echo "    [ -q queue (default: $queue) ]"
    echo "    [ -t time (default: $time) ]"
    echo "    [ -v vtune options ] [ -r vtune_out_dir ]"
    echo "    program options"
    exit 0
}

queue=spr
time=0:20:0
resultdir=
vtune=

while [ $# -gt 0 ] ; do
    if [ $1 = "-h" ] ; then
	usage && exit 0
    elif [ $1 = "-q" ] ; then
	shift && queue=$1 && shift
    elif [ $1 = "-t" ] ; then
	shift && time=$1 && shift
    elif [ $1 = "-v" ] ; then
	shift && vtune="$1" && shift
    elif [ $1 = "-r" ] ; then
	shift && resultdir=$1 && shift
    else
	break
    fi
done
if [ ! -z "${vtune}" ] ; then
    vtune="vtune ${vtune}"
    if [ ! -z "${resultdir}" ] ; then
	vtune="${vtune} -result-dir ${resultdir}"
    fi
    vtune="${vtune} --"
fi


cmdline="srun -A A-ccsc -p ${queue} -t ${time} -N 1 -n 1 ${vtune} $*"
echo "cmdline=${cmdline}"
eval ${cmdline}

