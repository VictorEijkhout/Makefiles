################################################################
####
#### Makefile for CANTERA installation
####
#### user group:
#### https://groups.google.com/g/cantera-users?pli=1
####
################################################################

PACKAGE = CANTERA
ABOUT = Chemical kinetics, thermodynamics, and transport tool suite
PACKAGEVERSION = 3.1.0
URL = https://github.com/Cantera/cantera
DOCURL = https://cantera.org/stable/develop/index.html#sec-compiling
MODE = seq
MODULES = boost eigen fmtlib hdf5 highfive ssundials yamlcpp \
    scons mkl

include ${MAKEINCLUDES}/Make.info

##
## CMake
##

# CMake build seems not supoprted
# https://github.com/Cantera/cantera/issues/1971

# Blas / Lapack
ifeq "${TACC_SYSTEM}" "vista"
  BLASLIB=${TACC_NVPL_LIB}/libnvpl_blas_lp64_seq.so 
  LAPACKLIB=${TACC_NVPL_LIB}/libnvpl_lapack_lp64_seq.so
else
  # mkl_intel_lp64;mkl_sequential;mkl_core;iomp5;pthread"
  BLASLIB=${TACC_MKL_LIB}/libmkl_sequential.so ${TACC_MKL_LIB}/libmkl_core.so 
  LAPACKLIB=${TACC_MKL_LIB}/libmkl_intel_lp64.so ${TACC_MKL_LIB}/libmkl_sequential.so 
endif

# include ${MAKEINCLUDES}/Make.cmake
# include ${MAKEINCLUDES}/Make.cbuild

info ::
	@echo "make default_install : configure build"
.PHONY: default_install
SCONS_EXTRA_INCS = ${TACC_EIGEN_INC}:${TACC_FMTLIB_INC}:${TACC_HIGHFIVE_INC}:${TACC_SUNDIALS_INC}:${TACC_YAMLCPP_INC}
SCONS_EXTRA_LIBS = ${TACC_EIGEN_LIB}:${TACC_FMTLIB_LIB}:${TACC_HIGHFIVE_LIB}:${TACC_SUNDIALS_LIB}:${TACC_YAMLCPP_LIB}
default_install : modules scons_config scons_build
.PHONY: scons_config scons_build
scons_config:
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
	 && ( \
	    cd $${srcdir} \
	     && echo "Cons build" \
	     && scons build env_vars=all \
	          CC=${TACC_CC} CXX=${TACC_CXX} FORTRAN=${TACC_FC} \
	          prefix=$${prefixdir} \
	          system_eigen='y' system_fmt='y' system_highfive='y' \
	          system_sundials='y' system_yamlcpp='y' \
	          extra_inc_dirs=${SCONS_EXTRA_INCS} \
	          extra_lib_dirs=${SCONS_EXTRA_LIBS} \
	          boost_inc_dir=${TACC_BOOST_INC} \
	          \
	          hdf_support='y' hdf_include=${TACC_HDF5_INC} hdf_libdir=${TACC_HDF5_LIB} \
	          googletest='none' \
	          blas_lapack_libs=mkl_rt blas_lapack_dir=$(MKLROOT)/lib/intel64 \
	    ) 2>&1 | tee $$installlog \
	 && requirenonzero logfilesdir \
	 && cp $${installlog} $${logfilesdir}/
scons_build :
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
	 && ( \
	    cd $${srcdir} \
	     && echo "Cons install" \
	     && scons install \
	    ) \
	 && make --no-print-directory varsmodule
include ${MAKEINCLUDES}/Make.vars
HASBIN = 1
PKGCONFIGLIB = pkgconfiglib
DEPENDSON = fmtlib hdf5 ssundials yamlcpp
EXTRAINSTALLPATHS = PYTHONPATH=lib/python3.9/site-packages
.PHONY: configlog
configlog :
	@source ${MAKEINCLUDES}/names.sh \
	 && setnames \
	    "${PACKAGE}" "${PACKAGEVERSION}" "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	    "${INSTALLPATH}" "${INSTALLROOT}" "${INSTALLEXT}" "${INSTALLVARIANT}" \
	    "${HOMEDIR}" "${BUILDDIRROOT}" "${MODE}" \
	    "${PREFIXOPTION}" "${PREFIXEXTRA}" \
	 && ls $${srcdir}/config.log


####
#### Download
####
DOWNLOADURL = https://github.com/Cantera/cantera/releases/download/v${PACKAGEVERSION}/cantera-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = https://github.com/Cantera/cantera.git
include ${MAKEINCLUDES}/Make.git
