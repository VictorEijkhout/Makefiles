################################################################
####
#### Makefile for OpenFoam-extended
####
################################################################

PACKAGE = foamx
PACKAGEVERSION = git
MODE = mpi
MODULES = cmake qt5 \
    swr ospray paraview

include ${MAKEINCLUDES}/Make.info

####
#### Totally custom configure/build
####

configure :
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} && export INSTALLROOT=${INSTALLROOT} \
	 && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} \
	 && setmodulenames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	 \
	 && cp bashrc_${TACC_SYSTEM}   $${srcdir}/etc/bashrc \
	 && cp prefs.sh.${TACC_SYSTEM} $${srcdir}/etc/prefs.sh \
	 && echo " .. set up files copied .. "

build :
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} && export INSTALLROOT=${INSTALLROOT} \
	 && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} \
	 && setmodulenames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" \
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

