################################################################
####
#### Makefile for LibCeed installation
####
################################################################

PACKAGE = LIBCEED
PACKAGEVERSION = git
CHMODR = 1
MODULES = petsc

include ${MAKEINCLUDES}/Make.info

JCOUNT = 12
.PHONY: configure build
configure ::
	@echo "there is no separate configure rule"
build :: modules
	@( \
	source ${MAKEINCLUDES}/names.sh \
	 && setnames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" \
	 && setmodulenames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	 \
	 && source ${MAKEINCLUDES}/compilers.sh \
	 && reportcompilers \
	 \
	 && requirenonzero srcdir \
	 && cd $$srcdir \
	 && make --no-print-directory \
	        -B OPT=-g -j ${JCOUNT} V=1 \
	        BUILDDIR=$$builddir/ \
	        prefix=$$installdir \
	        all install \
	 ) 2>&1 | tee install.log \
	 && echo "See install.log"

info ::
	@echo "make examples : copy and build examples"
.PHONY: examples
examples : petsc-installation
	@examples=examples-${PETSC_ARCH} \
	 && if [ "${TAU}" = "1" ] ; then \
	     examples=$${examples}-tau \
	    ; fi \
	 && rm -rf $${examples} \
	 && ( \
	    source ${MAKEINCLUDES}/names.sh \
	     && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} \
	     && ( cp -r $${srcdir}/examples $${examples} \
	         && cd $${examples} \
	         && builddir=$${builddir}-$${PETSC_ARCH} \
	         && if [ "${TAU}" = "1" ] ; then \
	                module load tau && export TAU_OPTIONS="$${TAU_OPTIONS} -optRevert" \
	                 && CCEXTRA="CC=tau_cc.sh" \
	                 && builddir=$${builddir}-tau \
	            ; fi \
	         && pcloc=$${installdir}/lib/pkgconfig \
	             && if [ -f $${pcloc}/libceed.pc ] ; then \
	                  echo "using $${pcloc}/libceed.pc:" \
	                   && cat $${pcloc}/libceed.pc \
	                ; else \
	                  echo "Can not find libceed.pc at: $${pcloc}" && exit 1 \
	                ; fi \
	         && export PKG_CONFIG_PATH=$${PKG_CONFIG_PATH}:$$pcloc \
	         && export CEED_DIR=$${PETSC_DIR}/$${PETSC_ARCH} \
	         && export LIBCEED_DIR=$${installdir} \
	         && echo "Make examples in `pwd`" \
	         && make --no-print-directory $${CCEXTRA} \
	        ) \
	    ) 2>&1 | tee $${examples}.log \
	 && echo && echo "Built examples in: $${examples}" \
	 && echo "see also $${examples}.log"

####
#### Clean
####

include ${MAKEINCLUDES}/Make.clean

####
#### Download stuff
####

include ${MAKEINCLUDES}/Make.git
GITREPO = https://github.com/CEED/libCEED
