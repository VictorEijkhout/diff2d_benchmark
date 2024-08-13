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
args = sys.argv
script = args[0]
args = args[1:]
if len( sys.argv )<2:
    print("Usage: {script} infile outfile")
    sys.exit(1)

file = args[0]
csvfile = args[1]
print(f"Processing file: <<{file}>> into <<{csvfile}>>")
keys = args[2:]
graphs = {}; fields = {}; columns = []
xaxis = None
for k in keys:
    cf = k.split(":"); c = cf[0]; f = cf[1]
    print(f"Retrieving column <<{c}>> at field <<{f}>>")
    if not xaxis: xaxis = c
    columns.append(c)
    fields[c] = int(f)
    graphs[c] = []

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
