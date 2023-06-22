################################################################
####
#### Makefile for JULIA installation
####
#### build instructions:
#### https://github.com/JuliaLang/julia/blob/master/doc/src/devdocs/build/build.md
####
################################################################

PACKAGE = JULIA
URL = https://julialang.org/
PACKAGEVERSION = git
MODE = seq

include ${MAKEINCLUDES}/Make.info

## JULIA_INSTALL = ${STOCKYARD}/julia/installation-${TACC_SYSTEM}-${TACC_FAMILY_COMPILER}
.PHONY: configure build
configure ::
	@echo "there is no configure rule"
build ::
	@source ${MAKEINCLUDES}/names.sh \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	 && ( \
	    echo "in-source build in $$srcdir" \
	     && cd $$srcdir \
	     && echo prefix=$$installdir > Make.user \
	     && make && make install \
	    ) 2>&1 | tee $$installlog
include ${MAKEINCLUDES}/Make.public
build ::
	@make public

GITREPO = https://github.com/JuliaLang/julia.git
include ${MAKEINCLUDES}/Make.git
