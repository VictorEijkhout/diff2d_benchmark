# -*- makefile -*-
################################################################
################################################################
####
#### This make file is part of the source of 
#### `Parallel Programing in Science and Engineering'
#### by Victor Eijkhout, copyright 2023-2024
####
#### Root makefile for the diff2d writeup
####
################################################################
################################################################

info ::
	@echo "make pdf [ FILE= (from ${FILES}, default=${FILE}) ]"

FILE = diff2d
FILES= diff2d cppcon2024
pdf : 
	@echo rsync ../../plots/*.csv ./plots/
	@TEXINPUTS=../../macros:../../graphics:../../snippets:../../chapters:${TEXINPUTS} \
	    pdflatex -shell-escape ${FILE}
	cp ${FILE}.pdf ${HOME}/Dropbox/Scicomp/
	cp ${FILE}.pdf ${HOME}/Current/istc/diff2d_benchmark/reports/

include ../../makefiles/Make.clean
MAKEFILESDIR := ../../${MAKEFILESDIR}
MAKETRACESTRING := "${MAKETRACESTRING} .."
