################################################################
####
#### Makefile for LibCeed installation
####
################################################################

PACKAGE = LIBCEED
URL = https://github.com/ceed/libceed
PACKAGEVERSION = git
CHMODR = 1
MODULES = petsc
MODE = mpi

include ${MAKEINCLUDES}/Make.info

JCOUNT = 12
.PHONY: configure build
PKGCONFIG = lib/pkgconfig
configure ::
	@echo "there is no separate configure rule"
build :: modules
	@\
	source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" "${VARIANT}" "${MODULENAME}" "${MODE}" \
	 \
	 && ( \
	    source ${MAKEINCLUDES}/compilers.sh \
	     && reportcompilers \
	     && requirenonzero srcdir \
	     && pushd $$srcdir \
	     && cmdline="make --no-print-directory \
	            -B OPT=-g -j ${JCOUNT} V=1 \
	            BUILDDIR=$$builddir/ \
	            prefix=$$installdir \
	            install "\
	     && echo "making: $$cmdline" && eval $$cmdline \
	     && popd && make module \
	     ) 2>&1 | tee $$installlog \
	 && echo "See $$installlog"
include ${MAKEINCLUDES}/Make.vars
default_install :
	@make configure build PACKAGEVERSION=git JCOUNT=${JCOUNT}

# ####
# #### Examples
# ####

# info :: 
# 	@echo "make examples : copy and build examples"
# .PHONY: examples
# examples : 
# 	@make --no-print-directory modules MODULES=libceed
# 	@examples=examples-${PETSC_ARCH} \
# 	 && if [ "${TAU}" = "1" ] ; then \
# 	     examples=$${examples}-tau \
# 	    ; fi \
# 	 && ( \
# 	    source ${MAKEINCLUDES}/names.sh \
# 	     && setnames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" \
# 	     && setmodulenames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" \
# 	     && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
# 	     && reportnames \
# 	     && source ${MAKEINCLUDES}/compilers.sh \
# 	     && setcompilers && reportcompilers \
# 	     \
# 	     && ( cp -r $${srcdir}/examples $${examples} \
# 	         && cd $${examples} \
# 	         && builddir=$${builddir}-$${PETSC_ARCH} \
# 	         && if [ "${TAU}" = "1" ] ; then \
# 	                module load tau && export TAU_OPTIONS="$${TAU_OPTIONS} -optRevert" \
# 	                 && CCEXTRA="CC=tau_cc.sh" \
# 	                 && builddir=$${builddir}-tau \
# 	            ; fi \
# 	         && export CEED_DIR=$${PETSC_DIR}/$${PETSC_ARCH} \
# 	         && export LIBCEED_DIR=$${installdir} \
# 	         && echo "Make examples in `pwd`" \
# 	         && make --no-print-directory $${CCEXTRA} \
# 	        ) \
# 	    ) 2>&1 | tee $${examples}.log \
# 	 && echo && echo "Built examples in: $${examples}" \
# 	 && echo "see also $${examples}.log"

####
#### Clean
####

include ${MAKEINCLUDES}/Make.clean

####
#### Download stuff
####

include ${MAKEINCLUDES}/Make.git
GITREPO = https://github.com/CEED/libCEED
