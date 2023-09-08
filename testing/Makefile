################################################################
####
#### Makefile for infrastructure testing
####
################################################################

PACKAGE = test 
PACKAGEVERSION = 0.0.1
MODE = seq

include ${MAKEINCLUDES}/Make.info

info ::
	@echo "make test"
.PHONY: test
test :
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} && export INSTALLROOT=${INSTALLROOT} \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	 && source ${MAKEINCLUDES}/compilers.sh \
	 && if [ "${MODE}" = "mpi" ] ; then \
	      setmpicompilers ; else setcompilers ; fi \
	 && cmake --version | head -n 1 \
	 && reportcompilers && echo
