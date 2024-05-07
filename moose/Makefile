################################################################
####
#### Makefile for Moose installation
####
################################################################

PACKAGE = moose
URL = https://mooseframework.inl.gov/
PACKAGEVERSION = git

MODE = mpi
MODULES = 

include ${MAKEINCLUDES}/Make.info

.PHONY: configure install default_install
default_install : configure install
include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install
notconfigure : modules
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "" "" "" "${MODE}" \
	 \
	 && source ${MAKEINCLUDES}/compilers.sh \
	 && setmpicompilers \
	    \
	 && requirenonzero builddir \
	 && rm -rf $$builddir && mkdir -p $$builddir \
	 && ( \
	    cd $${srcdir} \
	     && echo "not cd $${builddir}" \
	     && reportnames && echo \
	     && reportcompilers && echo \
	     && ( cd scripts \
	         && export MOOSE_JOBS=6 METHODS=opt \
	         && ./update_and_rebuild_petsc.sh \
	         && ./update_and_rebuild_libmesh.sh \
	         && ./update_and_rebuild_wasp.sh \
	        ) \
	     && ( cd test \
	         && make -j ${JCOUNT} \
	        ) \
	    ) 2>&1 | tee $$configurelog
	@echo && echo "CMake configuration ended: $$( date )" && echo 

# tagged versions fail, seem to think they are in a repository?
# TGZURL = https://github.com/idaholab/moose/archive/refs/tags/${PACKAGEVERSION}.tar.gz
# include ${MAKEINCLUDES}/Make.download
GITREPO = https://github.com/idaholab/moose.git
SUBMODULE = 1
BRANCH = master
include ${MAKEINCLUDES}/Make.git

include ${MAKEINCLUDES}/Make.clean

