################################################################
####
#### Makefile for OpenFoam-extended
####
################################################################

PACKAGE = foamx
PACKAGEVERSION = git
MODE = mpi
MODULES = cmake qt5 \
    swr metis
#ospray paraview

include ${MAKEINCLUDES}/Make.info

####
#### Totally custom configure/build
####

.PHONY: configure build
info ::
	@echo "make configure build"
configure : modules
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} && export INSTALLROOT=${INSTALLROOT} \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	 \
	 && if [ ! -f "bashrc_${TACC_SYSTEM}" ] ; then \
	        echo "No bashrc for system ${TACC_SYSTEM}" && exit 1 ; fi \
	 && cp bashrc_${TACC_SYSTEM}   $${srcdir}/etc/bashrc \
	 && if [ ! -f "prefs.sh.${TACC_SYSTEM}" ] ; then \
	        echo "No prefs.sh for system ${TACC_SYSTEM}" && exit 1 ; fi \
	 && cp prefs.sh.${TACC_SYSTEM} $${srcdir}/etc/prefs.sh \
	 && echo " .. set up files copied .. "

build : modules
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} && export INSTALLROOT=${INSTALLROOT} \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	 \
	 && export WM_COMPILER=Icc \
	 && export FOAM_INST_DIR=$${homedir} \
	 && export WM_PROJECT_DIR=$${srcdir} \
	 && ( \
	    cd $${srcdir} \
	     && set -x \
	     && source etc/bashrc \
	     && ./Allwmake.firstInstall \
	    ) 2>&1 | tee $${scriptdir}/install.log

####
####
####

GITREPO = git://git.code.sf.net/p/foam-extend/foam-extend-3.2
include ${MAKEINCLUDES}/Make.git

