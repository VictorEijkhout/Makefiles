################################################################
####
#### Makefile for CORRFUNC installation
####
################################################################

PACKAGE = CORRFUNC
ABOUT = Blazing fast correlation functions on the CPU.
PACKAGEVERSION = 2.5.3
URL = https://github.com/manodeep/Corrfunc
MODULES = gsl
MODE = seq

include ${MAKEINCLUDES}/Make.info

##
## Notes
## we have our own common.mk to set some compiler flags
## https://github.com/manodeep/Corrfunc/issues/322
##
.PHONY: build
build :
	@source ${MAKEINCLUDES}/names.sh \
	 && source ${MAKEINCLUDES}/compilers.sh \
	 && export scriptdir=$$(pwd) \
	 && rm -rf /tmp/pip-build-env* \
	 && setcompilers \
	 \
	 && setnames \
	    "${PACKAGE}"  "${PACKAGEVERSION}"  "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	    "${INSTALLPATH}"  "${INSTALLROOT}"  "${INSTALLEXT}"  "${INSTALLVARIANT}" \
	    "${HOMEDIR}"  "${BUILDDIRROOT}"  "${MODE}" \
	    "${PREFIXOPTION}"  "${PREFIXEXTRA}" \
	     \
	 && cp common.mk $${srcdir} \
	 && ( \
	    cd $${srcdir} \
	     && echo && echo "==== Regular compilation" && echo \
	     && mkdir -p $${prefixdir} \
	     && make realclean && make && make install \
	     && echo && echo "==== Python pip" && echo \
	     && python3 -m pip install . --user || /bin/true \
	    ) 2>&1 | tee $${installlog}
include ${MAKEINCLUDES}/Make.vars

.PHONY: default_install
info ::
	@echo "make default_install : build only"
default_install : build varsmodule public

DOWNLOADURL = https://github.com/manodeep/Corrfunc/archive/refs/tags/${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

