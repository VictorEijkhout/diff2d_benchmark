# -*- python -*-
################################################################
################################################################
####
#### This text file is part of the source of 
#### `Parallel Programming in MPI and OpenMP'
#### by Victor Eijkhout, copyright 2012-2024
####
#### graphs_extract.py : extract multiple graphs from a single file.
####
################################################################
################################################################

import re
import sys

import argparse

parser = argparse.ArgumentParser(description="Extract single graph from multiple files")

# Optional test and verbose flags
parser.add_argument('-t', '--test', action='store_true', help='set test mode')
parser.add_argument('-v', '--verbose', action='store_true', help='set verbose mode')

# Destination path
parser.add_argument('-n', '--name', type=str, default='graphs', help='Base name of the csv file')
parser.add_argument('-p', '--path', type=str, default='.', help='Path to the directory')

# Positional arguments
parser.add_argument('colonlist', nargs='*', help='List of items')

# Parse arguments
parsed_args = parser.parse_args()

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
## Remaning arguments are searchkey:column
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

with open(file,"r") as data:
    for line in data:
        line = line.strip()
        ## print(line)
        for c in fields:
            f = fields[c]
            if key_line:=re.match(f'{c}.*$',line):
                #print(f"Found <<{c}>> in line <<{line}>>")
                key_line = key_line.group().split()
                val = int( float( key_line[f] ) )
                graphs[c].append(val)

    for c in fields:
        print(f"Graph <<{c}>>: <<{graphs[c]}>>")

## this should really test any extension
# csvfile = re.sub('runout','csv',file)
# if csvfile==file:
#     print(f"Could not generate csv name from: <<{file}>>")
#     sys.exit(2)
with open(csvfile,"w") as csv:
    ## column names
    for c in columns:
        csv.write(c+", ")
    csv.write("\n")
    for values in zip( *[ graphs[c] for c in columns ] ):
        #
        # values = [ xaxis, yval1, yval2, ...]
        #
        csv.write( f"{values[0]}, " ); values = values[1:]
        for v in values:
            csv.write(f"{v}, ")
        csv.write("\n")
print(f"Written plottable data: <<{csvfile}>>")

sys.exit(0)

csvfile = re.sub('.csv','-eff.csv',csvfile)
if csvfile==file:
    print(f"Could not generate eff.csv name from: <<{file}>>")
    sys.exit(2)
with open(csvfile,"w") as csv:
    ## column names
    for c in columns:
        csv.write(c+", ")
    ## values, first x axis, then all y values
    csv.write("\n")
    reference_values = False
    for values in zip( *[ graphs[c] for c in columns ] ):
        #
        # values = [ xaxis, yval1, yval2, ...]
        #
        xvalue = values[0]
        csv.write( f"{xvalue}, " ); values = values[1:]
        if not reference_values: reference_values = values
        for col,yvalue in enumerate(values):
            csv.write(f"{yvalue*xvalue/reference_values[col]}, ")
        csv.write("\n")
print(f"Written plottable data: <<{csvfile}>>")
