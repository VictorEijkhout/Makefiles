################################################################
####
#### Makefile for Vtkhdf installation
####
#### for ticket
#### https://consult.tacc.utexas.edu/Ticket/Display.html?id=102066
####
################################################################

PACKAGE = vtkhdf
URL = https://github.com/johodges/fds.git
PACKAGEVERSION = git

MODE = mpi
MODULES = phdf5 

include ${MAKEINCLUDES}/Make.info

.PHONY: default_install
default_install : 
	@source ${MAKEINCLUDES}/names.sh \
	 && setnames \
	    "${PACKAGE}" "${PACKAGEVERSION}" "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	    "${INSTALLPATH}" "${INSTALLROOT}" "${INSTALLEXT}" "${INSTALLVARIANT}" \
	    "${HOMEDIR}" "${BUILDDIRROOT}"  "${MODE}" \
	    "${PREFIXOPTION}" "${PREFIXEXTRA}" \
	 && source ${MAKEINCLUDES}/compilers.sh && setcompilers \
	 && ( cd $$srcdir \
	     && cd Build && make clean \
	     && export HDF5_HOME=${TACC_HDF5_DIR} \
	     && cd impi_intel_linux \
	     && ./make_fds.sh \
	    ) 2>&1 | tee $$installlog \

GITREPO = https://github.com/johodges/fds.git
BRANCH= vtkhdf
include ${MAKEINCLUDES}/Make.git

include ${MAKEINCLUDES}/Make.clean

