################################################################
################################################################
####
#### This makefile is part of the source of 
#### "Parallel Computing for Science and Engineering"
#### by Victor Eijkhout, copyright 2013-2024
####
#### compile rules for diff2d library routines
####
################################################################
################################################################

##
## Compilation rules
##
include ../Make.depend
LANGUAGE = CXX
MODE = omp
include ../../../makefiles/Make.compile
CPPSTANDARD = 23

##
## Clean
##

include  ../../../makefiles/Make.clean
MAKEFILESDIR := ../../${MAKEFILESDIR}
MAKETRACESTRING := "${MAKETRACESTRING} .. .. .."
clean :: codeclean
	@/bin/rm -rf ${APPS} ${WRONGS}
