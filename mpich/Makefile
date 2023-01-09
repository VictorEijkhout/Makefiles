################################################################
####
#### Makefile for Mpich installation
####
#### https://www.mpich.org/2022/01/12/mpich-4-0rc3-released/
####
################################################################

PACKAGE = MPICH
PACKAGEROOT = ${HOME}/Installation
PACKAGEVERSION = 4.1a1
MAKEINCLUDES = ${HOME}/Software/makefiles

include ${MAKEINCLUDES}/Make.info

MPICH_HOME=${HOME}/Installation/mpich/
MPICH_SRC=${MPICH_HOME}/mpich-${PACKAGEVERSION}
MPICH_BUILD=${MPICH_HOME}/build-${PACKAGEVERSION}
MPICH_INSTALL=${MPICH_HOME}/installation-${PACKAGEVERSION}

MPICH_BIN=${MPICH_INSTALL}/bin

info ::
	@echo "make configure install"
configure ::
	@/bin/rm -rf ${MPICH_BUILD} ${MPICH_INSTALL}
	@( source ${MAKEINCLUDES}/names.sh \
	 && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} \
	 && cd $$srcdir \
	 && rm -rf $$installdir && mkdir -p $$installdir \
	 && export CC=clang && export CXX=clang++ && export FC=gfortran \
	 && if [ "${PACKAGEVERSION}" = "git" ] ; then \
	        ./autogen.sh \
	    ; fi \
	 && ./configure --prefix=$$installdir\
	) 2>&1 | tee configure.log

include ${MAKEINCLUDES}/Make.install

TGZURL = http://www.mpich.org/static/downloads/${PACKAGEVERSION}/mpich-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

SUBMODULE = submodule
## include ${MAKEINCLUDES}/Make.git
