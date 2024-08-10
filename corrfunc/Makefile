################################################################
####
#### Makefile for CORRFUNC installation
####
################################################################

PACKAGE = CORRFUNC
PACKAGEVERSION = 2.5.2
URL = https://github.com/manodeep/Corrfunc
MODE = seq

include ${MAKEINCLUDES}/Make.info

## include ${MAKEINCLUDES}/Make.configure
## include ${MAKEINCLUDES}/Make.install
install :
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
	     && make realclean && make && make install \
	     && echo && echo "==== Python pip" && echo \
	     && python3 -m pip install . --user \
	    ) 2>&1 | tee $${installlog}

.PHONY: default_install
info ::
	@echo "make default_install : configure build"
default_install : configure build public

TGZURL = https://github.com/manodeep/Corrfunc/archive/refs/tags/${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

