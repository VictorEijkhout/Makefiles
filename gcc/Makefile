################################################################
####
#### Makefile for Gcc installation
####
################################################################

PACKAGE = gcc
PACKAGEVERSION = 14.2
URL = 
MODE = core
MODULES = 

include ${MAKEINCLUDES}/Make.info

info ::
	@echo "================ bootstrap install"
	@echo "make build"
TOOLSET = ${TACC_FAMILY_COMPILER}

#
# Instructions:
# https://gcc.gnu.org/install/configure.html
#
include ${MAKEINCLUDES}/Make.configure
CONFIGINBUILDDIR = 1
CONFIGUREFLAGS =  --disable-multilib \
    --with-isl-include=${TACC_GCC_INC} --with-isl-lib=${TACC_GCC_LIB}
##--with-mpc=${TACC_GCC_INC}

#
# https://gcc.gnu.org/install/build.html
# https://gcc.gnu.org/install/finalinstall.html
#
build ::
	@source ${MAKEINCLUDES}/names.sh \
	 && setnames \
	    "${PACKAGE}" "${PACKAGEVERSION}" "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	    "${INSTALLPATH}" "${INSTALLROOT}" "${INSTALLEXT}" "${INSTALLVARIANT}" \
	    "${HOMEDIR}" "${BUILDDIRROOT}" "${MODE}" \
	    "${PREFIXOPTION}" "${PREFIXEXTRA}" \
	 && source ${MAKEINCLUDES}/compilers.sh && setcompilers \
	 && ( \
	    requirenonzero builddir \
	     && cd $${builddir} \
	     && make && make install \
	    ) 2>&1 | tee $${installlog}
# module 
# 	 && make --no-print-directory varsmodule \
# 	      CMAKEPREFIXPATHSET=1
include ${MAKEINCLUDES}/Make.vars
HASBIN = 1
include ${MAKEINCLUDES}/Make.public

GITREPO = https://gcc.gnu.org/git/gcc.git
BRANCH = releases/gcc-14
include ${MAKEINCLUDES}/Make.git
