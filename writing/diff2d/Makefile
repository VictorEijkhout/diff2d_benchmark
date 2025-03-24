# -*- makefile -*-
################################################################
################################################################
####
#### This make file is part of the source of 
#### the diff2d_benchmark
#### by Victor Eijkhout, copyright 2023-2025
####
#### Root makefile for the diff2d writeup
####
################################################################
################################################################

info ::
	@echo "make pdf [ FILE= (from ${FILES}, default=${FILE}) ]"

FILE = diff2d
FILES= diff2d cppcon2024proposal cppcon2024slides ixpug24eijkhout
LATEX = pdflatex
##xelatex
.PHONY: pdf onepass upload
pdf : getfiles
	@echo rsync ../../plots/*.csv ./plots/
	@if [ ! -d plots ] ; then ln -s ../plots . ; fi 
	@make --no-print-directory onepass
	@bibtex "${FILE}"
	@make --no-print-directory onepass
	@make --no-print-directory upload FILE="${FILE}"
upload :
	@if [ "$$(whoami)" = "eijkhout" ] ; then \
	    echo " .. copy to dropbox" \
	     && cp "${FILE}.pdf" ${HOME}/Dropbox/Scicomp/diff2d/ \
	 ; fi
onepass :
	@export TEXINPUTS=../macros:../graphics:../beamer:../graphs:${SNIPPETSDIR}:../:${TEXINPUTS} \
	 && ${LATEX} -shell-escape ${FILE}

.PHONY: getfiles
ISP  = ../../../introduction-scientific-programming-private
PCSE = ../../../parallel-programming-private
STC  = ../../../scientific-computing-private
getfiles :
	@if [ "$$(whoami)" = "eijkhout" ] ; then \
	    if [ -f ../../../scientific-computing-private/scicompbook.aux ] ; then \
	      rsync ../../../scientific-computing-private/scicompbook.aux . \
	     ; fi \
	     && rsync \
	            ${ISP}/macros/ispidxmacs.tex \
	            ${STC}/macros/{book,common,idx,listing,snippet}macs.tex \
	            ${PCSE}/macros/pcse{idx,snippet}macs.tex \
	            ../macros/ \
	     && rsync \
	            ${STC}/makefiles/Make.snippets \
	            ../../makefiles \
	 ; fi
clean ::
	@rm -f plots/*.{dep,dpth,idx,log,mw}

SNIPPETSOURCESROOT = ../../code/diff2d/
SCRIPTSDIR = ../../scripts
SNIPPETSSOURCES = . clps diy2d diy2e f77 iota kokkos2d lib mpl oned range seq span sycl
include ../../makefiles/Make.snippets
# SNIPPETSDIR = ../snippets
# CODEDIR = ../../code
# snippets-dir :
# 	@if [ ! -d ${SNIPPETSDIR} ] ; then \
# 	    mkdir -p ${SNIPPETSDIR} \
# 	    && echo "## Automatically generated files only here" > ${SNIPPETSDIR}/README \
# 	    ; fi
# snippets : snippets-dir
# 	@echo "Extracting snippets"
# 	@( export pwd=`pwd` \
# 	  && cd ${SNIPPETSDIR} \
# 	  && echo "C/C++ snippets" \
# 	  && for x in \
# 	        $${pwd}/${CODEDIR}/diff2d/*.{cpp,hpp} \
# 	        $${pwd}/${CODEDIR}/diff2d/*/*.{cpp,hpp} \
# 	        ; do \
# 	      echo "  file <<$$x>>" \
# 	      && cat $$x \
# 	            | sed -e 's/CHKERRQ(ierr);//' \
# 	            | awk -v fn=$${x##*/} \
# 	                'BEGIN {f=0; ff=0; e=0; o=0} \
# 	                 /snippet end/ {f=0; e=1} \
# 	                 /omitend/     {o=0} \
# 	                 f==1 && e==1  {print "   /* code omitted */ " > file; e=0 } \
# 			 f==1 && ff==0 {print "// " fn > file; ff=1} \
# 	                 f==1 && !/solution/ && !/solved/ && !/snipthis/ && !/snippetskip/ && o==0 && !/omit/ \
# 	                               {print $$0 > file} \
# 	                 /snippet /     {f=1; file=$$2 } \
# 	                 /omitbegin/    {o=1} \
# 	                ' \
# 	      ; \
# 	     done \
# 	  && for s in `ls` ; do \
# 	       python $${pwd}/../../scripts/cutspace.py $$s ; \
# 	     done \
# 	 ) >snippets.log 2>&1
# 	@echo "(see snippets.log)"

include ../../makefiles/Make.clean
MAKEFILESDIR := ../../${MAKEFILESDIR}
MAKETRACESTRING := "${MAKETRACESTRING} .."
