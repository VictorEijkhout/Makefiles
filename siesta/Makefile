################################################################
####
#### Makefile for Siesta installation
####
#### Fortran debug flags: https://gitlab.com/siesta-project/siesta/-/issues/397
#### libfdf cloning in F jail: https://gitlab.com/siesta-project/siesta/-/issues/417
#### scalapack fortran interface: https://gitlab.com/siesta-project/siesta/-/issues/438
####
#### Vista:
#### https://gitlab.com/siesta-project/siesta/-/issues/441
####
################################################################

PACKAGE = siesta
URL = https://siesta-project.org/siesta/index.html \
    https://gitlab.com/siesta-project/siesta/-/releases
## see below if this ever reaches 5.4
PACKAGEVERSION = 5.4.0

MODE = mpi
MODULES = fftw3 parallelnetcdf \
    mkl 

include ${MAKEINCLUDES}/Make.info

PKGCONFIGLIB = pkgconfig

##
## cmake version
##
# https://docs.siesta-project.org/projects/siesta/en/stable/installation/build-issues.html#siesta-does-not-recognize-intel-mkl-libraries
## CMAKEPREP = find . -name \*make -exec sed -i -e 's/\-xHost//' -e 's/\-check//' {} \;
## CMAKEBUILDDEBUG = 1
# for nvidia add flag: -Mextend ?
# --log-level=debug --log-context
CMAKEFLAGS =  \
    -D SIESTA_WITH_FLOOK=FALSE \
    -D Fortran_FLAGS_DEBUG=-O0 \
    -D CMAKE_Fortran_FLAGS=-fPIE \
    \
    -D SIESTA_WITH_WANNIER90=ON \
    -D NetCDF_ROOT=${TACC_PARALLELNETCDF_DIR} \
    \
    -D SIESTA_WITH_FFTW=TRUE

##
## Where do we find scalapack?
##
ifdef TACC_SCALAPACK_DIR
  CMAKEFLAGS += \
    -D SCALAPACK_LIBRARY_DIR=${TACC_SCALAPACK_LIB} \
    -D SCALAPACK_LIBRARY='-L${TACC_SCALAPACK_LIB} -lscalapack -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl' \
    -D LAPACK_LIBRARY='-L${TACC_MKL_LIB} -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl'

else ifdef TACC_MKL_DIR
  CMAKEFLAGS += \
    -D LAPACK_LIBRARY='-L${TACC_MKL_LIB} -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl' \
    -D SCALAPACK_LIBRARY='-L${TACC_MKL_LIB} -lmkl_scalapack_lp64 -Wl,--start-group -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -Wl,--end-group -lpthread -lm -ldl'

else ifdef TACC_NVPL_DIR
  CMAKEFLAGS += \
    -D BLAS_LIBRARY_DIR=${TACC_NVPL_LIB} \
    -D BLAS_LIBRARY='-L${TACC_NVPL_LIB} -lnvpl_blas_lp64_seq -lnvpl_blas_core' \
    -D LAPACK_LIBRARY_DIR=${TACC_NVPL_LIB} \
    -D LAPACK_LIBRARY='-L${TACC_NVPL_LIB} -lnvpl_lapack_lp64_seq -lnvpl_lapack_core -lnvpl_blas_lp64_seq -lnvpl_blas_core' \
    -D SCALAPACK_LIBRARY_DIR=${TACC_NVPL_LIB} \
    -D SCALAPACK_LIBRARY='-L${TACC_NVPL_LIB} -lnvpl_scalapack_lp64 -lnvpl_blacs_lp64_openmpi5'
endif

##
## Intel no xHost FFS
## https://gitlab.com/siesta-project/siesta/-/issues/501
## will be fixed in 5.4
##
# ifeq "${TACC_FAMILY_COMPILER}" "intel"
#   BOOTSTRAP = find . -type f -exec sed -e 's/-xHost//' -e 's/-xHOST//' -i {} \;
# #sed -e 's/-xHost//' -e 's/-xHOST//' -i ./Config/cmake/toolchains/*.cmake
# endif

HASBIN = 1
CPTOINSTALLDIR = Examples
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build


TGZURL = https://gitlab.com/siesta-project/siesta/-/releases/${PACKAGEVERSION}/downloads/siesta-${PACKAGEVERSION}.tar.gz
# do not use "assets": https://gitlab.com/siesta-project/siesta/-/archive/${PACKAGEVERSION}/siesta-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://gitlab.com/siesta-project/siesta.git
SUBMODULE = 1
include ${MAKEINCLUDES}/Make.git

# Diagnose scalapack linking problem by 
# editing SiestaCheckLinalg.cmake
# to read:
# check_fortran_source_compiles(
# "
# external :: blacs_gridinit
# call blacs_gridinit()
# end
# "
# scalapack_has_blacs_gridinit SRC_EXT F90 OUTPUT_VARIABLE scalapack_compile_output )

# Nvidia / Openmpi doesn't like something:
# [  2%] Building Fortran object _deps/libgridxc-build/CMakeFiles/libgridxc-lib.dir/src/mesh3d.F90.o
# cd /work/00434/eijkhout/siesta/build-5.2.1-vista-nvidia24.9-openmpi5.0.5_nvc249/_deps/libgridxc-build && /home1/apps/nvidia/Linux_aarch64/24.9/compilers/bin/nvfortran -DHAVE_MPI -Dlibgridxc_lib_EXPORTS -I/work/00434/eijkhout/siesta/build-5.2.1-vista-nvidia24.9-openmpi5.0.5_nvc249/_deps/libgridxc-build/include -isystem /opt/apps/nvidia24/openmpi/5.0.5_nvc249/include -isystem /opt/apps/nvidia24/openmpi/5.0.5_nvc249/lib -fPIE  -O2 -gopt -module include -fPIC -c /work/00434/eijkhout/siesta/build-5.2.1-vista-nvidia24.9-openmpi5.0.5_nvc249/_deps/libgridxc-src/src/mesh3d.F90 -o CMakeFiles/libgridxc-lib.dir/src/mesh3d.F90.o
# NVFORTRAN-W-0189-Argument number 1 to mpi_allreduce: association of scalar actual argument to array dummy argument (/work/00434/eijkhout/siesta/build-5.2.1-vista-nvidia24.9-openmpi5.0.5_nvc249/_deps/libgridxc-src/src/mesh3d.F90: 2016)
# NVFORTRAN-W-0189-Argument number 2 to mpi_allreduce: association of scalar actual argument to array dummy argument (/work/00434/eijkhout/siesta/build-5.2.1-vista-nvidia24.9-openmpi5.0.5_nvc249/_deps/libgridxc-src/src/mesh3d.F90: 2016)
#   0 inform,   2 warnings,   0 severes, 0 fatal for optimizetransferorder
# nvfortran-Fatal-/home1/apps/nvidia/Linux_aarch64/24.9/compilers/bin/tools/fort1 TERMINATED by signal 11
