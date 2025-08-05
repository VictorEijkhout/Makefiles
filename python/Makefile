################################################################
####
#### Makefile for Python installation
####
################################################################

PACKAGE = python
ABOUT = programming language
PACKAGEVERSION = 3.12.5
## problems: 13.0
URL = https://www.python.org/
MODE = mpi
MODULES = hdf5 mkl sqlite

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
	    echo "Calling install.sh script" \
	     && requirenonzero srcdir \
	     && requirenonzero builddir && mkdir -p $$builddir && rm -rf $$builddir/* \
	     && requirenonzero prefixdir && mkdir -p $$prefixdir && rm -rf $$prefixdir/* \
	     && $${scriptdir}/install.sh \
	          -v ${PACKAGEVERSION} \
	          --build $${builddir} \
	          --srcdir $${srcdir} \
	          --prefix $${prefixdir} \
	    ) 2>&1 | tee $$installlog \
	 && make --no-print-directory varsmodule

config-help :
	@source ${MAKEINCLUDES}/names.sh \
	 && setnames \
	    "${PACKAGE}"  "${PACKAGEVERSION}"  "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	    "${INSTALLPATH}"  "${INSTALLROOT}"  "${INSTALLEXT}"  "${INSTALLVARIANT}" \
	    "${HOMEDIR}"  "${BUILDDIRROOT}"  "${MODE}" \
	    "${PREFIXOPTION}"  "${PREFIXEXTRA}" \
	 && cd $${srcdir} && ./configure --help 

CMAKEPREFIXPATHSET = 1
HASBIN = 1

include ${MAKEINCLUDES}/Make.vars
include ${MAKEINCLUDES}/Make.public

TGZURL = https://www.python.org/ftp/python/${PACKAGEVERSION}/Python-${PACKAGEVERSION}.tgz
include ${MAKEINCLUDES}/Make.download

include ${MAKEINCLUDES}/Make.clean

