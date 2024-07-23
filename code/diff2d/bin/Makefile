include ../../../../makefiles/Make.clean
recursive_clean :: localclean
localclean :: 
	@variants="$$( cd .. && make --no-print-directory variants )" \
	 && echo " .. variants: $${variants}" >/dev/null \
	 && rm -f $${variants}
MAKEFILESDIR := ../../../${MAKEFILESDIR}
MAKETRACESTRING := "${MAKETRACESTRING} .. .. .. "
