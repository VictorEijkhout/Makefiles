################################################################
####
#### Makefile for Trilinos installation
####
#### https://www.docs.trilinos.org/files/TrilinosBuildReference.html
####
################################################################

PACKAGE = trilinos
PACKAGEVERSION = 13.4.1
PACKAGEVERSION_DOWNLOAD = 13-4-1
MODE = mpi
MODULES = boost
#boost swig

include ${MAKEINCLUDES}/Make.info

info ::
	@echo "================ cmake configure:"
	@echo "make configure"
	@echo "    [ PACKAGEVERSION=... (default=${PACKAGEVERSION}) ]"
	@echo "    [ PYTHON=ON/OFF (default: ${PYTHON} ) ]"
.PHONY: configure
LIBDIR = lib64
CMAKE_PREFIXPATH_SET = lib/cmake
PYTHON = ON
configure :: 
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} \
	 && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} "" ${PACKAGEBASENAME} \
	 && setmodulenames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	 \
	 && source ${MAKEINCLUDES}/compilers.sh \
	 && setmpicompilers \
	    \
	 && requirenonzero builddir \
	 && rm -rf $$builddir && mkdir -p $$builddir \
	 && ( \
	    cd $${builddir} \
	     && reportnames && echo \
	     && reportcompilers && echo \
	     \
	     && export HAS_HDF5=0 && export HAS_NETCDF=0 && export HAS_MUMPS=0 \
	     && export HAS_SEACAS=0 && export HAS_STK=0 \
	     \
	     && export HAS_YAMLCPP=0 \
	     && if [ "$${HAS_YAMLCPP}" = "1" ] ; then \
	            source $${scriptdir}/../yaml-cpp/yamlcpp_vars.sh ; fi \
	     \
	     && export TRILINOS_LOCATION=${TRILINOS_SCRIPTDIR} \
	     && if [ "${LMOD_FAMILY_COMPILER}" = "intel" ] ; then \
	          export MKLFLAG="-mkl" ; fi \
	     && export LMOD_MKL_INC=${MKLROOT}/include \
	     && export LMOD_MKL_LIB=${MKLROOT}/lib/intel64 \
	     \
	     && export HAS_PYTHON=${PYTHON} \
	     \
	     && set -x \
	     && if [ -f $${scriptdir}/trilinos_cmake_${TACC_SYSTEM}.sh ] ; then \
	            source $${scriptdir}/trilinos_cmake_${TACC_SYSTEM}.sh \
	        ; else \
	            source $${scriptdir}/trilinos_cmake.sh \
	        ; fi \
	     \
	     && export PYTHON_VERSION=3 \
	     && export PYTHON_VERSION=$${PYTHON_VERSION%.*} \
	     && export PYTHON_LOAD_FLAG=${TACC_PYTHON_LIB}/libpython$${PYTHON_VERSION}.so \
	    ) 2>&1 | tee $$configurelog \
	 && make --no-print-directory varsmodule

include ${MAKEINCLUDES}/Make.cbuild

TGZURL = https://github.com/trilinos/Trilinos/archive/refs/tags/trilinos-release-${PACKAGEVERSION_DOWNLOAD}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = git@github.com:trilinos/Trilinos.git
include ${MAKEINCLUDES}/Make.git
