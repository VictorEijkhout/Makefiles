################################################################
####
#### Makefile for Python installation
####
################################################################

PACKAGE = python
PACKAGEVERSION = 3.12.4
URL = https://www.python.org/
MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info

info ::
	@echo "================ installation"
	@echo "make default_install"
.PHONY: default_install
default_install :
	@source ${MAKEINCLUDES}/names.sh \
	 && setnames \
	    "${PACKAGE}"  "${PACKAGEVERSION}"  "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	    "${INSTALLPATH}"  "${INSTALLROOT}"  "${INSTALLEXT}"  "${INSTALLVARIANT}" \
	    "${HOMEDIR}"  "${BUILDDIRROOT}"  "${MODE}" \
	    "${PREFIXOPTION}"  "${PREFIXEXTRA}" \
	 && source ${MAKEINCLUDES}/compilers.sh && setcompilers \
	 && ( \
	    requirenonzero prefixdir && mkdir -p $$prefixdir && rm -rf $$prefixdir/* \
	     && $${scriptdir}/install.sh -v ${PACKAGEVERSION} \
	          --prefix $${prefixdir} \
	    ) 2>&1 | tee $$installlog \
	 && make --no-print-directory varsmodule \
	      CMAKEPREFIXPATHSET=1

include ${MAKEINCLUDES}/Make.vars
include ${MAKEINCLUDES}/Make.public

TGZURL = https://www.python.org/ftp/python/${pythonver}/Python-${pythonver}.tgz
include ${MAKEINCLUDES}/Make.download

include ${MAKEINCLUDES}/Make.clean

