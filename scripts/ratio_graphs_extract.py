# -*- python -*-
################################################################
################################################################
####
#### This text file is part of the source of 
#### `Parallel Programming in MPI and OpenMP'
#### by Victor Eijkhout, copyright 2012-2024
####
#### ratio_graphs_extract.py : extract graphs from multiple files
####
#### 1. generate absolute time graph
#### 2. generate ratios to the first graph given
####
################################################################
################################################################

import os
import re
import sys

import argparse

parser = argparse.ArgumentParser(description="Extract graphs from multiple files")

# Optional verbose flag
parser.add_argument('-v', '--verbose', action='store_true', help='set verbose mode')
# Optional test flag
parser.add_argument('-t', '--test', action='store_true', help='test only')

# Destination path
parser.add_argument('-n', '--name', type=str, default='graphs', help='Base name of the csv file')
parser.add_argument('-p', '--path', type=str, default='.', help='Path to the directory')

# Positional arguments
parser.add_argument('colonlist', nargs='*', help='List of items')

# Parse arguments
parsed_args = parser.parse_args()

test = parsed_args.test
verbose = parsed_args.verbose
destpath = parsed_args.path
destname = parsed_args.name

if parsed_args.colonlist:
    args = parsed_args.colonlist
else:
    print("No positional arguments were provided."); sys.exit(1)

##
## first bunch of arguments is filename:key
## parse until something isn't a filename anymore
##
print("================")
graphs = {}; files = {}; firstfile = False; filekeys = []
while len(args)>0:
    fk = args[0]; fk = fk.split(":"); f = fk[0]; k = fk[1]
    if verbose:
        print( f"Testing <<{f}>> : <<{k}>>" )
    if os.path.isfile(f):
        print(f"File: {f}")
        if not firstfile: firstfile = f
        filekeys.append(k)
        ## store file handle for `k'
        files[k] = open(f,"r")
        ## create an empty graph for `k'
        graphs[k] = []
        ## parse next argument
        args = args[1:]
    else: break
if len(files)==0:
    print("No files found"); sys.exit(1)
else:
    print(f"Processing files: <<{files.keys()}>>")

##
## Remaining arguments are searchkey:column
## First key will the x axis, others are multiple graphs
##
keys = args
fields = {}; columns = []
xaxis = None
for k in keys:
    fc = k.split(":"); f = fc[0]; c = fc[1]
    print(f"Retrieving column <<{c}>> at field <<{f}>>")
    if not xaxis: xaxis = f
    columns.append(f)
    fields[f] = int(c)
    graphs[f] = []

##
## Go through the files
##
for fk,fh in files.items():
    xaxis = False
    for line in fh:
        line = line.strip()
        for fc in fields.items():
            f = fc[0]; c = fc[1]
            if key_line:=re.match(f'{f}.*$',line):
                ##
                ## find field 'f' starting the current line
                ## if found, extract column 'c'
                ## and append to 'graphs[fk]
                ##
                key_line = key_line.group().split()
                #print(f"Found <<{f}>>=<<{key_line[c]}>> in line <<{line}>>")
                if not xaxis:
                    xaxis = int( key_line[c] )
                else:
                    graphs[fk].append( [xaxis,key_line[c]] )
                    xaxis = False
for f in files.keys():
    print(f"Graph <<{f}>>: <<{graphs[f]}>>")

if test: sys.exit(0)

##
## write all data to csv file
## this is probably timing 
##
# firstfile = re.sub(r'.*/','',firstfile)
# csvfile = re.sub('runout','csv',firstfile)
csvfile = f"{destpath}/{destname}.csv"
print(f"Writing to csv file <<{csvfile}>>>")
with open(csvfile,"w") as csv:
    ## column names
    for c in ["cores"]+[ k for k in files.keys() ]:
        csv.write(c+", ")
    csv.write("\n")
    for kvs in zip( *[ graphs[k] for k in files.keys() ] ):
        pvalue = False
        values_for_p = [] # time values for specific p value and all files
        #print(kvs)
        for k,v in kvs:
            #print(f"{k}:{v}")
            if not pvalue: pvalue = k
            values_for_p.append(v)
        for v in [ pvalue ] + values_for_p:
            csv.write( f"{v}, " )
        csv.write("\n")
print(f"Written plottable data: <<{csvfile}>>")

##
## assuming above plots times,
## plot rat to first method
##
csvfile = re.sub('.csv','-rat.csv',csvfile)
with open(csvfile,"w") as csv:
    ## column names
    for c in ["cores"]+[ k for k in files.keys() ]:
        csv.write(c+", ")
    ## values, first x axis, then all y values
    csv.write("\n")
    for kvs in zip( *[ graphs[k] for k in files.keys() ] ):
        # go through a zipped tuple of proc/time pairs
        # each zipped tuple correspond to one proc value
        # so we gather the time values for that p value
        # the first time value is used as reference
        reference_value = False; proc_value = False; values_for_p = []
        for k,v in kvs:
            # we note the proc value. should really test for consistency
            if not proc_value: proc_value = k
            # the first graph is what we measure against
            if not reference_value: reference_value = float(v)
            # record the ratio of this time to the reference time
            values_for_p.append(float(v)/reference_value)
        csv.write( f"{proc_value}, " )
        for v in values_for_p:
            csv.write(f"{v}, ")
        csv.write("\n")
print(f"Written ratio data: <<{csvfile}>>")
