################################################################
####
#### Makefile for Trilinos installation
####
#### https://www.docs.trilinos.org/files/TrilinosBuildReference.html
####
#### ExodusII / Seacas issue: https://github.com/trilinos/Trilinos/issues/12022
####
################################################################

PACKAGE = trilinos
URL = https://github.com/trilinos/Trilinos \
    https://www.docs.trilinos.org/files/TrilinosBuildReference.html
PACKAGEVERSION = 16.1.0
PACKAGEVERSION_DOWNLOAD = 16-1-0
#15-1-0
MODE = hybrid
MODULES = boost swig phdf5 pnetcdf mkl
MODULE_MINIMUM_VERSION_intel = 20

include ${MAKEINCLUDES}/Make.info

info ::
	@echo "================ cmake configure:"
	@echo "make configure"
	@echo "    [ PACKAGEVERSION=... (default=${PACKAGEVERSION}) ]"
	@echo "    [ PYTHON=ON/OFF (default: ${PYTHON} ) ]"
	@echo "    [ TEUCHOS=ON/OFF (default: ${TEUCHOS} ) ]"
PYTHON = OFF

CMAKEPREFIXPATHSET = 1
LANGUAGE = CXX
CPPSTANDARD = 20
# ThyraEpetraAdapters don't work in 13
HAS_ADAPTERS = ON

# dealii expects Amesos Epetra AztecOO
# the epetra flag triggers the other two
HAS_EPETRA = ON

HAS_BOOST = ON
HAS_EXODUS = OFF # doesn't work anymore in 16
HAS_FEI = ${HAS_ADAPTERS}
HAS_KOKKOS = ON
HAS_NETCDF = OFF
HAS_PNETCDF = ON
HAS_HDF5 = ON
HAS_OPENMP = ON
HAS_Teka = ${HAS_ADAPTERS}
.PHONY: configure config-help
info ::
	@echo "make config-help"
config-help :
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "" "" "" "${MODE}" \
	 \
	 && cd $${srcdir} && rm -rf config-help-dir && mkdir config-help-dir \
	     && cd config-help-dir && cmake .. 2>&1 | grep "Final set of .*enabled packages"
configure :: modules
	@source ${MAKEINCLUDES}/names.sh \
	 && export scriptdir=$$(pwd) \
	 \
	 && setnames \
	    "${PACKAGE}"  "${PACKAGEVERSION}"  "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	    "${INSTALLPATH}"  "${INSTALLROOT}"  "${INSTALLEXT}"  "${INSTALLVARIANT}" \
	    "${HOMEDIR}"  "${BUILDDIRROOT}"  "${MODE}" \
	    "${PREFIXOPTION}"  "${PREFIXEXTRA}" \
	     \
	 && source ${MAKEINCLUDES}/compilers.sh \
	 && setmpicompilers \
	     \
	 && ( \
	    reportnames	&& reportcompilers \
	     && echo \
	     && requirenonzero builddir \
	     && rm -rf $$builddir && mkdir -p $$builddir \
	     && cd $${builddir} \
	     \
	     && export HAS_MUMPS=0 \
	     && export HAS_STK=0 \
	     \
	     && export HAS_YAMLCPP=0 \
	     && if [ "$${HAS_YAMLCPP}" = "1" ] ; then \
	            source $${scriptdir}/../yaml-cpp/yamlcpp_vars.sh ; fi \
	     \
	     && export TRILINOS_LOCATION=$${scriptdir} \
	     && if [ [ "a" = "b" ] ; then \
	          exit 1 \
	        ; elif [ "${TACC_SYSTEM}" = "vista" ] ; then \
	          export MKLFLAG="-Mnvpl -Mscalapack" \
	           && export TACC_MKL_INC=${TACC_NVPL_INC} \
	           && export TACC_MKL_LIB=${TACC_NVPL_LIB} \
	           && export TRILINOS_BLAS_LIBS='nvpl_blas_lp64_seq;nvpl_blas_core' \
	           && export TRILINOS_LAPACK_LIBS='nvpl_lapack_lp64_seq;nvpl_lapack_core' \
	        ; elif [ "${TACC_FAMILY_COMPILER}" = "intel" ] ; then \
	          export MKLFLAG="-mkl" \
	           && export TACC_MKL_INC=${MKLROOT}/include \
	           && export TACC_MKL_LIB=${MKLROOT}/lib/intel64 \
	           && export TRILINOS_BLAS_LIBS="mkl_intel_lp64;mkl_core;pthread" \
	           && export TRILINOS_LAPACK_LIBS="mkl_intel_lp64;mkl_sequential;mkl_core;pthread" \
	           && echo "removed iomp5" >/dev/null \
	        ; elif [ "${TACC_FAMILY_COMPILER}" = "gcc" ] ; then \
	          export TACC_MKL_INC=${MKLROOT}/include \
	           && export TACC_MKL_LIB=${MKLROOT}/lib/intel64 \
	           && export TRILINOS_BLAS_LIBS='mkl_intel_lp64;mkl_core;pthread' \
	           && export TRILINOS_LAPACK_LIBS='mkl_intel_lp64;mkl_sequential;mkl_core;pthread' \
	        ; fi \
	     \
	     && export HAS_PYTHON=${PYTHON} \
	     && export HAS_BOOST=${HAS_BOOST} \
	     && export HAS_HDF5=${HAS_HDF5} \
	     && export HAS_NETCDF=${HAS_NETCDF} \
	     && export HAS_OPENMP=${HAS_OPENMP} \
	     \
	     && export DEPRECATED=OFF \
	     && export HAS_TEUCHOS=ON \
	     && export HAS_EPETRA=${HAS_EPETRA} \
	     \
	     && set -x \
	     && if [ -f $${scriptdir}/trilinos_cmake_${TACC_SYSTEM}.sh ] ; then \
	            source $${scriptdir}/trilinos_cmake_${TACC_SYSTEM}.sh \
	        ; else \
	            BUILDDIR=$$builddir SRCDIR=$$srcdir \
	              source $${scriptdir}/trilinos_cmake.sh \
	        ; fi \
	     \
	     && export PYTHON_VERSION=3 \
	     && export PYTHON_VERSION=$${PYTHON_VERSION%.*} \
	     && export PYTHON_LOAD_FLAG=${TACC_PYTHON_LIB}/libpython$${PYTHON_VERSION}.so \
	    ) 2>&1 | tee $$configurelog
	@echo && echo "CMake configuration ended: $$( date )" && echo 
##	 && make --no-print-directory varsmodule

include ${MAKEINCLUDES}/Make.cbuild
.PHONY: netcdfcopy notnetcdfcopy
build :: notnetcdfcopy
notnetcdfcopy :
	@echo 
	@echo "Warning: not copying netcdf files"
	@echo 
netcdfcopy ::
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "" "" "" "${MODE}" \
	 && requirenonzero TACC_PNETCDF_INC \
	 && cmdline="cp ${TACC_NETCDF_INC}/netcdf*.h $${installdir}/include" \
	 && echo "copy netcdf.h: " && echo $$cmdline \
	 && eval $$cmdline \
	 && if [ $$? -gt 0 ] ; then error "ERROR copy failed" ; fi

.PHONY: default_install
info ::
	@echo "make default_install : configure build"
default_install : configure build

TGZURL = https://github.com/trilinos/Trilinos/archive/refs/tags/trilinos-release-${PACKAGEVERSION_DOWNLOAD}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = git@github.com:trilinos/Trilinos.git
include ${MAKEINCLUDES}/Make.git

# Which netcdf?
# https://github.com/trilinos/Trilinos/issues/12854

# Build fail on Chaco
# https://github.com/trilinos/Trilinos/issues/12857
