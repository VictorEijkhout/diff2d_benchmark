# -*- makefile -*-
################################################################
################################################################
####
#### This make file is part of the source of 
#### the diff2d_benchmark
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
	@if [ ! -d plots ] ; then ln -s ../plots . ; fi 
	@export TEXINPUTS=../macros:../graphics:../snippets:../:${TEXINPUTS} \
	 && pdflatex -shell-escape ${FILE}

include ../../makefiles/Make.clean
MAKEFILESDIR := ../../${MAKEFILESDIR}
MAKETRACESTRING := "${MAKETRACESTRING} .."
