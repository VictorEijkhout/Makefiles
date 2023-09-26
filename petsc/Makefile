################################################################
####
#### Makefile for PETSc local installation
####
################################################################

PACKAGE = PETSC
PACKAGEVERSION = 3.19.5
MODE = mpi

include ${MAKEINCLUDES}/Make.info
list_installations ::
	@source ${MAKEINCLUDES}/names.sh \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" "${VARIANT}" "${MODULENAME}" "${MODE}" \
	 && echo "installations in $${srcdir}: $$( cd $$srcdir && for d in * ; do if [ -d $${d} ] ; then echo $${d} ; fi ; done )"
modules ::
	@if [ "${TACC_FAMILY_COMPILER}" = "gcc" ] ; then \
	    make --no-print-directory moduletest MODULE=mkl ; fi

info ::
	@echo "make build "
	@echo "     [CUSTOMEXT=customextension]"

JCOUNT = 4

PETSC_SRC_DIR = ${PACKAGEROOT}/petsc/petsc-${PACKAGEVERSION}
OPTFLAGS=-O2 -fPIC
LOPTFLAGS=-O1 -fPIC

MPI_FLAGS = --with-mpi-compilers=1
ifeq "${TACC_SYSTEM}" "nvarm"
  BLAS_FLAGS = --download-fblaslapack=yes
  MPI_FLAGS = --with-mpi-compilers=1 --with-cc=`which mpicc` --with-cxx=`which mpicxx` 
else ifeq "${TACC_SYSTEM}" "macbookair"
  BLAS_FLAGS=
##    --with-cc=`which mpicc` --with-cxx=`which mpicxx` 
else
  BLAS_FLAGS=--with-blas-lapack-dir=${MKLROOT}
endif
# ifeq "${FORTRAN}" "1"
#   MPI_FLAGS += --with-fc=`which mpif90`
# endif

info ::
	@echo "make installed : list installations"
.PHONY: installed
installed :
	@ls -d ${PETSC_SRC_DIR}/${PETSCARCH}*

CUDA = 0
info ::
	@echo "     [CUDA=0/1 (default: ${CUDA})]"
CUDA_OPTIONS = --with-cuda=1 --with-cuda-dir=${TACC_CUDA_DIR} --with-cuda-gencodearch=70

TIME =

##
## Petsc Architecture
##

.PHONY: petscarch petscshortarch
petscshortarch :
	@export PETSC_EXT="" \
	  \
	 && if [ "${SCALAR}" = "complex" ] ; then \
	    export PETSC_EXT=$${PETSC_EXT}complex ; fi \
	  \
	 && if [ "${PRECISION}" = "single" ] ; then \
	    export PETSC_EXT=$${PETSC_EXT}single ; fi \
	  \
	 && if [ "${PRECISION}" = "quad" ] ; then \
	    export PETSC_EXT=$${PETSC_EXT}single ; fi \
	  \
	 && if [ "${INT}" = "64" ] ; then \
	    export PETSC_EXT=$${PETSC_EXT}i64 ; fi \
	 \
	 && if [ ! -z "${ALIGN}" ] ; then \
	    export PETSC_EXT=$${PETSC_EXT}a${ALIGN} ; fi \
	  \
	 && if [ "${CUDA}" = "1" ] ; then \
	    export PETSC_EXT=$${PETSC_EXT}cuda ; fi \
	  \
	 && if [ "${DEBUG}" = "1" ] ; then \
	    export PETSC_EXT=$${PETSC_EXT}debug ; fi \
	  \
	 && if [ ! -z "${CUSTOMEXT}" ] ; then \
	    export PETSC_EXT=${CUSTOMEXT} ; fi \
	 && echo $${PETSC_EXT}
petscarch :
	@source ${MAKEINCLUDES}/names.sh \
	 && systemnames && compilernames \
	 && setnames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" "${VARIANT}" "${MODULENAME}" "${MODE}" \
	 && requirenonzero compilercode && requirenonzero mpicode \
	 && petscarch=$${taccsystemcode}-$${compilercode}$${compilershortversion}-$${mpicode} \
	 && petscext=$$( make --no-print-directory petscshortarch ) \
	 && if [ -z $${petscext} ] ; then \
	      echo $${petscarch} ; \
	    else \
	      echo $${petscarch}-$${petscext} ; \
	    fi

####
#### Packages
####

include Makefile.packages

X = 0
info ::
	@echo "    [ X=0/1, default=${X}]"
WITH_X = "--with-x=1"

configure :
	@echo && echo "there is no configure rule" && echo
include ${MAKEINCLUDES}/Make.vars
NOINSTALL = 1
build :: modules
	@if [ ${CUDA} -eq 1 -a -z "${TACC_CUDA_DIR}" ] ; then \
	    echo && echo "Please load cuda module" && echo && exit 1 \
	    ; fi
	@source ${MAKEINCLUDES}/names.sh \
	 && export INSTALLROOT=${INSTALLROOT} \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" "${VARIANT}" "${MODULENAME}" "${MODE}" \
	 \
	 && export petscarch=`make --no-print-directory \
	                SCALAR=${SCALAR} PRECISION=${PRECISION} INT=${INT} \
	                ALIGN=${ALIGN} \
	                CUDA=${CUDA} DEBUG=${DEBUG} CUSTOMEXT=${CUSTOMEXT} \
	                petscarch` \
	 && export package_options="$$( make --no-print-directory package_options \
	        CHACO=${CHACO} EIGEN=${EIGEN} FFTW=${FFTW} HDF5=${HDF5} \
	        HYPRE=${HYPRE} KOKKOS=${KOKKOS} LIBCEED=${LIBCEED} MPE=${MPE} \
	        MUMPS=${MUMPS} PETSC4PY=${PETSC4PY} PTSCOTCH=${PTSCOTCH} SLEPC=${SLEPC} \
	        SUPERLU=${SUPERLU} TRIANGLE=${TRIANGLE} ZOLTAN=${ZOLTAN} \
	        )" \
	 && export param_settings="$$( make --no-print-directory param_settings \
	        )" \
	 && ( \
	    if [ ! -d $$srcdir ] ; then \
	        echo && echo "SRCDIR does not exist: $$srcdir" && echo && exit 1 ; fi \
	     && reportnames \
	     && cd $$srcdir \
	     && echo "Configuring version=${PACKAGEVERSION} with arch=$${petscarch} in `pwd`" \
	     && echo "Installing with packages: $${package_options}" \
	     && echo "Installing with parameters: $${param_settings}" \
	     && echo \
	     && rm -rf $${petscarch} \
	     && export PETSC_DIR=`pwd` \
	     && export I_MPI_FABRICS=shm:tmi \
	     && export PETSC_ARCH=$${petscarch} \
	     && cmdline="${TIME} python3 ./configure --prefix=$${installdir} \
	        ${MPI_FLAGS} \
	        --useThreads=0 \
	        --with-np=${JCOUNT} \
	        CFLAGS=\"${OPTFLAGS}\" CXXFLAGS=\"${OPTFLAGS}\" \
	        `if [ \"${FORTRAN}\" = 1 ] ; then echo FFLAGS=\"${OPTFLAGS}\" ; fi `\
	        `if [ ${X} -eq 1 ]    ; then echo ${WITH_X} ; fi ` \
	        `if [ ${CUDA} -eq 1 ] ; then echo ${CUDA_OPTIONS} ; fi ` \
	        "${BLAS_FLAGS}" \
	        $${package_options} $${param_settings} \
	        " \
	     && echo cmdline=$$cmdline && eval $$cmdline \
	     && ${TIME} make MAKE_NP=${JCOUNT} \
	     && make install \
	  ) 2>&1 | tee $${installlog}
info ::
	@echo "make default_install : with 0 & 1 debug"
.PHONY: default_install
default_install :
	@for d in 0 1 ; do \
	    make configure build DEBUG=$d JCOUNT=${JCOUNT} \
	; done 
build :: makemodule
.PHONY: makemodule
info ::
	@echo "make makemodule"
makemodule :
	@source ${MAKEINCLUDES}/names.sh \
	 && export INSTALLROOT=${INSTALLROOT} \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" "${VARIANT}" "${MODULENAME}" "${MODE}" \
	 \
	 && export petscarch=`make --no-print-directory \
	                SCALAR=${SCALAR} PRECISION=${PRECISION} INT=${INT} \
	                ALIGN=${ALIGN} \
	                CUDA=${CUDA} DEBUG=${DEBUG} CUSTOMEXT=${CUSTOMEXT} \
	                petscarch` \
	 && shortversion=${PACKAGEVERSION} && shortversion=$${shortversion%.*} \
	 && petscshortarch=$$( make --no-print-directory petscshortarch ) \
	 && echo "making module with version=$${shortversion} extra=<<$${petscshortarch}>>" \
	 && ( make --no-print-directory varsmodule \
	        PACKAGEVERSION=${PACKAGEVERSION} \
	        PACKAGEMODULEVERSION=$${shortversion} \
	        MODULEVERSIONEXTRA=$${petscshortarch} \
	        PKGCONFIGLIB=pkgconfig \
	        EXTRAVARS="PETSC_DIR=$${installdir} PETSC_ARCH=$${petscarch}" \
	  ) 2>&1 | tee -a $${installlog}

info ::
	@echo "make biginstall : lots of packages"
.PHONY: biginstall allinstall
allinstall :
	@for ext in \
	        "" debug i64 i64-debug single single-debug \
	        complex complex-debug complexi64 complexi64-debug \
	        nohdf5 i64nohdf5 \
	; do \
	    echo "Installing for extension <<$${ext}>>" \
	     && make --no-print-directory biginstall EXT=$${ext} \
	; done
biginstall :
	@ext=${EXT} \
	   && case $${ext} in ( *debug* ) DEBUG=1 ;; ( * ) DEBUG= ;; esac \
	   && case $${ext} in ( *i64* ) INT=64 ;; ( * ) INT=32 ;; esac \
	   && case $${ext} in ( *single* ) PRECISION=single ;; ( * ) PRECISION=double ;; esac \
	   && case $${ext} in ( *complex* ) SCALAR=complex ;; ( * ) SCALAR=real ;; esac \
	   \
	   && case $${ext} in ( *nohdf5* ) HDF= ;; ( *single* ) HDF5= ;; ( * ) HDF5=1 ;; esac \
	   \
	   && case ${TACC_FAMILY_COMPILER} in ( gcc ) FORTRAN= ;; ( * ) FORTRAN=1 ;; esac \
	   \
	   && cmdline="make --no-print-directory configure build public \
	      PACKAGEVERSION=${PACKAGEVERSION} JCOUNT=${JCOUNT} \
	      DEBUG=$${DEBUG} INT=$${INT} PRECISION=$${PRECISION} SCALAR=$${SCALAR} \
	      FORTRAN=$${FORTRAN} \
	      HDF5=$${HDF5} \
	      " \
	   && echo "Big install: $$cmdline" && eval $$cmdline \

public ::
	@echo "public rule disabled"

info ::
	@echo "make config-help config-view"
	@echo "     [PACKAGEVERSION= (default=${PACKAGEVERSION})]"
config-help :
	@cd ${PETSC_SRC_DIR} \
	 && ./configure -help
config-view :
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} && export INSTALLROOT=${INSTALLROOT} \
	 && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} \
	 \
	 && export petscarch=`make --no-print-directory \
	                SCALAR=${SCALAR} INT=${INT} PRECISION=${PRECISION} ALIGN=${ALIGN} \
	                CUDA=${CUDA} DEBUG=${DEBUG} CUSTOMEXT=${CUSTOMEXT} \
	                petscarch` \
	 && cat $$srcdir/configure.log
log-view :
	@cat ${PACKAGEROOT}/petsc/petsc-${PACKAGEVERSION}/configure.log

include ${MAKEINCLUDES}/Make.git
GITREPO = https://gitlab.com/petsc/petsc.git
BRANCH = main

include ${MAKEINCLUDES}/Make.download
TGZURL = https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-${PACKAGEVERSION}.tar.gz

clean ::
	@/bin/rm -f *~
