#!/usr/bin/env python
################################################################
####
#### This source is part of
#### Parallel Programming for Science and Engineering
#### Victor Eijkhout eijkhout@tacc.utexas.edu
#### copyright 2010-2024
####
#### python program to strip initial spaces from source snippets
####
################################################################

import re
import sys

if len(sys.argv)<2:
    print("Usage: cutspace.py fn"); sys.exit(1)

fn = sys.argv[1]
f = open(fn,"r")
firstline = True
nwhite = None
lines = []
for line in f.readlines():
    lines.append(line)
    if firstline and ( re.match(r' *//',line) or re.match(r' *##',line) ):
        continue # file name, don't measure
    white = re.search("^([ ]*)([^ ])",line).groups()[0]
    white = len(white)
    if nwhite is None or (white<nwhite and white>0):
        nwhite = white
    firstline = False
f.close()
f = open(fn,"w")
firstline = True
for line in lines:
    if firstline and  ( re.match(r' *//',line) or re.match(r' *##',line) ): # file name, align
        f.write( re.sub("^ *","",line) )
    else: # strip minimum number of spaces
        f.write( re.sub("^"+nwhite*" ","",line) )
    firstline = False
f.flush()
f.close()
