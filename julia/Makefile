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
	    echo "In-source build in $$srcdir" \
	     && rm -rf $$builddir && mkdir -p $$builddir \
	     && cd $$srcdir \
	     && echo " .. git reset and clean" \
	     && git reset --hard && git clean -x -f -d \
	     && echo prefix=$$installdir > Make.user \
	     && echo " .. Make.user : " && cat Make.user \
	     && echo " .. starting make" \
	     && make --no-print-directory \
	     && make --no-print-directory install \
	    ) 2>&1 | tee $$installlog
##  BUILDROOT=$$builddir
include ${MAKEINCLUDES}/Make.public
build ::
	@make public

GITREPO = https://github.com/JuliaLang/julia.git
include ${MAKEINCLUDES}/Make.git
include ${MAKEINCLUDES}/Make.clean

