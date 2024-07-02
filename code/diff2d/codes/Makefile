include ../../../../makefiles/Make.clean
recursive_clean ::
clean :: local_clean
	@rm -f $$( make --no-print-directory variants -f ../Makefile )
MAKEFILESDIR := ../../../${MAKEFILESDIR}
MAKETRACESTRING := "${MAKETRACESTRING} .. .. .. "
