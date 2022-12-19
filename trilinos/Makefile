################################################################
####
#### Makefile for Trilinos installation
####
################################################################

PACKAGE = trilinos
PACKAGEROOT = ${STOCKYARD}
PACKAGEVERSION = 13-4-1
MAKEINCLUDES = ${STOCKYARD}/makefiles
MODE = mpi
MODULES =
#boost swig

include ${MAKEINCLUDES}/Make.info

TRILINOS_HOME=${STOCKYARD}/trilinos
TRILINOS_SRC=${TRILINOS_HOME}/trilinos-${PACKAGEVERSION}
TRILINOS_INSTALL=${TRILINOS_HOME}/installation-${TACC_FAMILY_COMPILER}

info ::
	@echo "================ cmake configure:"
	@echo "make configure install"
	@echo "    [ PACKAGEVERSION=... (default=${PACKAGEVERSION}) ]"
.PHONY: configure install config-help
configure :: 
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} \
	 && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} "" ${PACKAGEBASENAME} \
	 && rm -rf $$builddir && mkdir -p $$builddir \
	 && ( \
	    cd $${builddir} \
	    \
	     && export HAS_PYTHON=0 \
	     && export PYTHON_VERSION=${TACC_FAMILY_PYTHON_VERSION} \
	     && export PYTHON_VERSION=$${PYTHON_VERSION%.*} \
	     && export PYTHON_LOAD_FLAG=${TACC_PYTHON_LIB}/libpython$${PYTHON_VERSION}.so \
	    \
	     && export HAS_HDF5=0 && export HAS_NETCDF=0 && export HAS_MUMPS=0 \
	     && export HAS_SEACAS=0 && export HAS_STK=0 \
	    \
	     && export HAS_YAMLCPP=0 \
	     && if [ "$${HAS_YAMLCPP}" = "1" ] ; then \
	            source $${scriptdir}/../yaml-cpp/yamlcpp_vars.sh ; fi \
	    \
	     && export TRILINOS_LOCATION=${TRILINOS_SCRIPTDIR} \
	     && export COPTFLAGS="-O2" && export MKLFLAG="-mkl" \
	     && source $${scriptdir}/trilinos_cmake.sh \
	     && ( \
	          echo "export TACC_TRILINOS_INC=$${INSTALL_PREFIX}/include" \
	           && echo "export TACC_TRILINOS_LIB=$${INSTALL_PREFIX}/lib64" \
	        ) > $${varfile} \
	     && echo "Set environment by sourcing: trilinos_vars.sh" \
	    ) 2>&1 | tee configure.log \
	 && make --no-print-directory varsfile VARSFILE=$$varfile \
	        PACKAGE=${PACKAGE} \
	        INSTALLDIR="$${installdir}" \
	        LIBDIR="$${installdir}/lib" \
	        INCDIR="$${installdir}/include" \
	        PKGCONFIGSET="${PKGCONFIGSET}" 

include ${MAKEINCLUDES}/Make.vars
include ${MAKEINCLUDES}/Make.cbuild

TGZURL = https://github.com/trilinos/Trilinos/archive/refs/tags/trilinos-release-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = git@github.com:trilinos/Trilinos.git
include ${MAKEINCLUDES}/Make.git
