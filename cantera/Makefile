################################################################
####
#### Makefile for CANTERA installation
####
################################################################

PACKAGE = CANTERA
ABOUT = 
PACKAGEVERSION = 3.1.0
URL = https://github.com/Cantera/cantera
DOCURL = https://cantera.org/stable/install/index.html
MODE = seq
MODULES = boost eigen hdf5 highfive sundials mkl

include ${MAKEINCLUDES}/Make.info

##
## CMake
##

# something missing
# https://github.com/Cantera/cantera/issues/1971
CMAKESUBDIR = interfaces/python_sdist

# Blas / Lapack
ifeq "${TACC_SYSTEM}" "vista"
  BLASLIB=${TACC_NVPL_LIB}/libnvpl_blas_lp64_seq.so 
  LAPACKLIB=${TACC_NVPL_LIB}/libnvpl_lapack_lp64_seq.so
else
  # mkl_intel_lp64;mkl_sequential;mkl_core;iomp5;pthread"
  BLASLIB=${TACC_MKL_LIB}/libmkl_sequential.so ${TACC_MKL_LIB}/libmkl_core.so 
  LAPACKLIB=${TACC_MKL_LIB}/libmkl_intel_lp64.so ${TACC_MKL_LIB}/libmkl_sequential.so 
endif

CMAKEFLAGS = -Wno-dev \
    -D LAPACK_LIBRARIES='${LAPACKLIB}' \
    -D BLAS_LIBRARIES='${BLASLIB}'

include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild

info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build public
cons_install :
	@source ${MAKEINCLUDES}/names.sh \
	 && setnames \
	    "${PACKAGE}" "${PACKAGEVERSION}" "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	    "${INSTALLPATH}" "${INSTALLROOT}" "${INSTALLEXT}" "${INSTALLVARIANT}" \
	    "${HOMEDIR}" "${BUILDDIRROOT}" "${MODE}" \
	    "${PREFIXOPTION}" "${PREFIXEXTRA}" \
	 && requirenonzero prefixdir && mkdir -p $${prefixdir} && rm -rf $${prefixdir}/* \
	 && requirenonzero builddir && rm -rf $$builddir && mkdir -p $$builddir \
	 && requirenonzero logfilesdir \
	 && echo "logfiles dir: $${logfilesdir}." \
	 && mkdir -p $${logfilesdir} \
	 && source ${MAKEINCLUDES}/compilers.sh && setcompilers \
	 && cp boost-build-nvcc/nvcc.jam $${srcdir}/tools/build/src/tools \
	 && ( \
	    cd $${srcdir} \
	    scons build env_vars=all CC=icx CXX=icpx FORTRAN=ifx blas_lapack_libs=mkl_rt blas_lapack_dir=$(MKLROOT)/lib/intel64 \
	    ) 2>&1 | tee $$installlog \
	 && requirenonzero logfilesdir \
	 && cp $${installlog} $${logfilesdir}/ \
	 && make --no-print-directory varsmodule \
	      CMAKEPREFIXPATHSET=1
# configure build public

####
#### Download
####
TGZURL = https://github.com/Cantera/cantera/releases/download/v${PACKAGEVERSION}/cantera-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = https://github.com/Cantera/cantera.git
include ${MAKEINCLUDES}/Make.git
