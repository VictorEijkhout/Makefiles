################################################################
####
#### Makefile for Blis installation
####
################################################################

PACKAGE = blis
PACKAGEROOT = ${STOCKYARD}
PACKAGEVERSION = git
MODE = seq

include ${MAKEINCLUDES}/Make.info

CONFIGUREFLAGS = generic
include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install

GITREPO = https://github.com/flame/blis.git
include ${MAKEINCLUDES}/Make.git

old_install ::
	@/bin/rm -rf ${BLIS_BUILD} ${BLIS_INSTALL}
	@( \
	cd ${BLIS_SRC} \
	 && mkdir -p ${BLIS_BUILD} ${BLIS_INSTALL} \
	 && ./configure --prefix=${BLIS_INSTALL} \
	        `if [ ! -z "${MT}" ] ; then echo "--enable-threading=openmp" ; fi` \
	        CC=clang generic \
	 && make \
	 && make install \
	) 2>&1 | tee install.log
