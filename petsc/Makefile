################################################################
####
#### Makefile for PETSc local installation
####
################################################################

PACKAGE = PETSC
PACKAGEVERSION = 3.20.2
MODE = mpi

include ${MAKEINCLUDES}/Make.info
list_installations ::
	@source ${MAKEINCLUDES}/names.sh \
	 && setdirlognames \
	      "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	      "${PACKAGEBASENAME}" "${VARIANT}" "${MODULENAME}" "${MODE}" \
	 && echo "installations in $${srcdir}: $$( cd $$srcdir && for d in * ; do if [ -d $${d} ] ; then echo $${d} ; fi ; done )"
modules ::
	@if [ "${TACC_FAMILY_COMPILER}" = "gcc" -a "${TACC_SYSTEM}" != "macbookair" ] ; then \
	    make --no-print-directory moduletest MODULE=mkl ; fi

info ::
	@echo "make build "
	@echo "     [CUSTOMEXT=customextension]"

JCOUNT = 4

OPTFLAGS=-O2 -fPIC
LOPTFLAGS=-O1 -fPIC

ifeq "${TACC_SYSTEM}" "nvarm"
  BLAS_FLAGS = --download-fblaslapack=yes
else ifeq "${TACC_SYSTEM}" "macbookair"
  BLAS_FLAGS= --download-fblaslapack=yes
else ifeq "${TACC_SYSTEM}" "grace"
  BLAS_FLAGS = --download-fblaslapack=yes
else
  BLAS_FLAGS=--with-blas-lapack-dir=${MKLROOT}
endif

CUDA = 0
info ::
	@echo "     [CUDA=0/1 (default: ${CUDA})]"
CUDA_OPTIONS = --with-cuda=1 --with-cuda-dir=${TACC_CUDA_DIR} \
	--with-cudac=nvcc --with-cuda-gencodearch=70

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
	 && if [ "${PRECISION}" = "single" ] ; then \
	    export PETSC_EXT=$${PETSC_EXT}single ; fi \
	 && if [ "${PRECISION}" = "quad" ] ; then \
	    export PETSC_EXT=$${PETSC_EXT}single ; fi \
	 && if [ "${INT}" = "64" ] ; then \
	    export PETSC_EXT=$${PETSC_EXT}i64 ; fi \
	 && if [ ! -z "${ALIGN}" ] ; then \
	    export PETSC_EXT=$${PETSC_EXT}a${ALIGN} ; fi \
	 && if [ "${CUDA}" = "1" ] ; then \
	    export PETSC_EXT=$${PETSC_EXT}cuda ; fi \
	 && if [ "${DEBUG}" = "1" ] ; then \
	    export PETSC_EXT=$${PETSC_EXT}debug ; fi \
	  \
	 && if [ ! -z "${CUSTOMEXT}" ] ; then \
	    export PETSC_EXT=${CUSTOMEXT} ; fi \
	 && echo $${PETSC_EXT}
petscarch :
	@source ${MAKEINCLUDES}/names.sh \
	 && systemnames && compilernames \
	 && notetoself "TODO dir names are no needed?" \
	 && setdirlognames \
	      "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	      "${PACKAGEBASENAME}" "${VARIANT}" "${MODULENAME}" "${MODE}" \
	 && requirenonzero compilercode && requirenonzero mpicode \
	 && petscarch=$${taccsystemcode}-$${compilercode}$${compilershortversion}-$${mpicode} \
	 && petscext=$$( make --no-print-directory petscshortarch ) \
	 && if [ -z $${petscext} ] ; then \
	      echo ${PACKAGEVERSION}; \
	    else \
	      echo ${PACKAGEVERSION}-$${petscext} ; \
	    fi

####
#### Packages
####

include Makefile.packages

configure :
	@echo && echo "there is no configure rule" && echo
build :: modules
	@if [ ${CUDA} -eq 1 -a -z "${TACC_CUDA_DIR}" ] ; then \
	    echo && echo "Please load cuda module" && echo && exit 1 \
	    ; fi
	@source ${MAKEINCLUDES}/names.sh \
	 && export INSTALLROOT=${INSTALLROOT} \
	 && export petscarch=`make --no-print-directory \
	                SCALAR=${SCALAR} PRECISION=${PRECISION} INT=${INT} \
	                ALIGN=${ALIGN} \
	                CUDA=${CUDA} DEBUG=${DEBUG} CUSTOMEXT=${CUSTOMEXT} \
	                petscarch` \
	 && echo && echo "Deduce petc arch <<$$petscarch>>" && echo \
	 && notetoself "VARIANT is petsc arch" \
	 && notetoself "INSTALLEXT zero: ext=petsc-ls6-intel24-impi22-arch" \
	 && notetoself "used for build/install-ext" \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "" \
	    "${PACKAGEBASENAME}" "$${petscarch}" "${MODULENAME}" "${MODE}" \
	 && echo "Deduce install dir <<$$installdir>>" && echo \
	 \
	 && export package_options="$$( make --no-print-directory package_options \
	        AMGX=${AMGX} CHACO=${CHACO} EIGEN=${EIGEN} FFTW=${FFTW} HDF5=${HDF5} \
	        HYPRE=${HYPRE} KOKKOS=${KOKKOS} LIBCEED=${LIBCEED} MPE=${MPE} \
	        MUMPS=${MUMPS} PTSCOTCH=${PTSCOTCH} SLEPC=${SLEPC} \
	        SUPERLU=${SUPERLU} TRIANGLE=${TRIANGLE} ZOLTAN=${ZOLTAN} \
	        \
	        PETSC4PY=${PETSC4PY} SLEPC4PY=${SLEPC4PY} \
	        )" \
	 && export param_settings="$$( make --no-print-directory param_settings )" \
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
	     && echo "Kludge: adding <<$${installdir}/lib>> to pythonpath" \
	     && export PYTHONPATH=${PYTHONPATH}:$${installdir}/lib \
	     && echo "PYTHONPATH=$${PYTHONPATH}" \
	     && cmdline="${TIME} python3 ./configure --prefix=$${installdir} \
	        CFLAGS=\"${OPTFLAGS}\" CXXFLAGS=\"${OPTFLAGS}\" \
	        `if [ \"${FORTRAN}\" = 1 ] ; then echo FFLAGS=\"${OPTFLAGS}\" ; fi `\
	        `if [ ${CUDA} -eq 1 ] ; then echo ${CUDA_OPTIONS} ; else echo --with-cudac=0 ; fi ` \
	        "${BLAS_FLAGS}" \
	        $${package_options} $${param_settings} \
	        " \
	     && echo "Cmdline=$$cmdline" && eval $$cmdline \
	     && echo "PYTHONPATH=$${PYTHONPATH}" \
	     && make PETSC_DIR=`pwd` PETSC_ARCH=$${PETSC_ARCH} MAKE_NP=${JCOUNT} all \
	     && make PETSC_DIR=`pwd` PETSC_ARCH=$${PETSC_ARCH} MAKE_NP=${JCOUNT} install \
	  ) 2>&1 | tee $${installlog}
include ${MAKEINCLUDES}/Make.vars
include ${MAKEINCLUDES}/Make.public
## build :: installation_public

info ::
	@echo "make installation_public : open installation (also done as part of build)"
.PHONY: installation_public
## build :: installation_public
## I CAN NOT FIGURE OUT THE RIGHT WAY TO GENERATE THE NAME OF THE INSTSALL DIRECTORY
installation_public :
	@ export petscarch=`make --no-print-directory \
	         SCALAR=${SCALAR} PRECISION=${PRECISION} INT=${INT} \
	         ALIGN=${ALIGN} \
	         CUDA=${CUDA} DEBUG=${DEBUG} CUSTOMEXT=${CUSTOMEXT} \
	         petscarch` \
	 && make public INSTALLEXT="" \
	    $$( if [ ! -z ${CUSTOMEXT} ] ; then \
	                echo VARIANT=${CUSTOMEXT} \
	         ; else echo VARIANT=$${petscarch} ; fi )

info ::
	@echo "make makemodule : custom rule for module"
build :: makemodule
.PHONY: makemodule
makemodule :
	@source ${MAKEINCLUDES}/names.sh \
	 && export INSTALLROOT=${INSTALLROOT} \
	 && export petscarch=`make --no-print-directory \
	                SCALAR=${SCALAR} PRECISION=${PRECISION} INT=${INT} \
	                ALIGN=${ALIGN} \
	                CUDA=${CUDA} DEBUG=${DEBUG} CUSTOMEXT=${CUSTOMEXT} \
	                petscarch` \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" \
	    "" "${PACKAGEBASENAME}" "$${petscarch}" "${MODULENAME}" "${MODE}" \
	 && echo "makemodule: deduce install dir <<$$installdir>>" \
	 \
	 && shortversion=${PACKAGEVERSION} && shortversion=$${shortversion%.*} \
	 && petscshortarch=$$( make --no-print-directory petscshortarch ) \
	 && echo "making module with version=$${shortversion} extra=<<$${petscshortarch}>>" \
	 && ( make --no-print-directory varsmodule  \
	        INSTALLDIR="$${installdir}" \
	        PACKAGEVERSION=${PACKAGEVERSION} \
	        INSTALLEXT= VARIANT=$${petscarch} \
	        PACKAGEMODULEVERSION=$${shortversion} \
	        MODULEVERSIONEXTRA=$${petscshortarch} \
	        PKGCONFIGLIB=pkgconfig \
	        EXTRAVARS="PETSC_DIR=$${installdir} PETSC_ARCH=$${petscarch}" \
	        EXTRAINSTALLPATHS="PYTHONPATH=lib" \
	    ) 2>&1 | tee -a $${installlog}

##
## installation modes
##

.PHONY: default_install
default_install : configure build

info ::
	@echo "make biginstall [ EXT=i64debug and such ]"
	@echo "    build specific configuration"
	@echo "    packages specified by HDF5=1 and such"
.PHONY: biginstall allinstall
EXTENSIONS = \
    "" debug i64 i64-debug single single-debug \
    complex complex-debug complexi64 complexi64-debug \
    nohdf5 i64nohdf5
info ::
	@echo "make allinstall [ JCOUNT= PACKAGEVERSION= ]"
	@echo "    biginstall for all extensions except cuda"
allinstall :
	@log=$$( pwd )/allinstall.log \
	 && for ext in ${EXTENSIONS} ; do \
	      echo "================ Installing for extension <<$${ext}>>" \
	        | tee -a $$log \
	       && make --no-print-directory biginstall EXT=$${ext} \
	            JCOUNT=${JCOUNT} PACKAGEVERSION=${PACKAGEVERSION} \
	            OUTERLOG=$${log} \
	    ; done \
	 && echo && echo "See: $${log}" && echo
biginstall :
	@ext=${EXT} \
	   && echo "Selecting basic options" \
	   && case $${ext} in ( *debug* ) DEBUG=1 ;; ( * ) DEBUG= ;; esac \
	   && case $${ext} in ( *i64* ) INT=64 ;; ( * ) INT=32 ;; esac \
	   && case $${ext} in ( *single* ) PRECISION=single ;; ( * ) PRECISION=double ;; esac \
	   && case $${ext} in ( *complex* ) SCALAR=complex ;; ( * ) SCALAR=real ;; esac \
	   \
	   && echo "Selecting packages" \
	   && if [ "${CUDA}" != "1" -a "${AMGX}" = "1" ] ; then \
	        echo " .. disable AMGX because no CUDA" && AMGX=0 \
	      ; else AMGX=${AMGX} ; fi \
	   && case $${ext} in ( *nohdf5* ) HDF= ;; ( *single* ) HDF5= ;; ( * ) HDF5=${HDF5} ;; \
	        esac \
	   && if [ "${PARMETIS}" = "1" ] ; then MUMPS=1 ; else MUMPS=${MUMPS}; fi \
	   \
	   && echo "Compiler issues" \
	   && case ${TACC_FAMILY_COMPILER} in \
	        ( gcc ) echo " .. no FORTRAN for gcc" && FORTRAN= ;; \
	        ( * ) FORTRAN=${FORTRAN} ;; esac \
	   \
	   && cmdline="make --no-print-directory configure build  \
	      PACKAGEVERSION=${PACKAGEVERSION} JCOUNT=${JCOUNT} \
	      DEBUG=$${DEBUG} INT=$${INT} PRECISION=$${PRECISION} SCALAR=$${SCALAR} \
	      $$( if [ ! -z ${CUSTOMEXT} ] ; then echo CUSTOMEXT=${CUSTOMEXT} ; fi ) \
	      \
	      AMGX=$${AMGX} EIGEN=${EIGEN} FFTW3=${FFTW3} HDF5=$${HDF5} HYPRE=${HYPRE} \
	      MUMPS=$${MUMPS} PARMETIS=${PARMETIS} \
	      \
	      CUDA=${CUDA} FORTRAN=$${FORTRAN} \
	      PETSC4PY=${PETSC4PY} SLEPC4PY=${SLEPC4PY} \
	      " \
	   && if [ ! -z "${OUTERLOG}" ] ; then \
	        echo "Big install: $$cmdline" | tee -a "${OUTERLOG}" \
	      ; else \
	        echo "Big install: $$cmdline"\
	      ; fi \
	   && eval $$cmdline 2>&1 | tee big.log \
	   && echo && echo "See: big.log" && echo

info ::
	@echo "make config-help config-view"
	@echo "     [PACKAGEVERSION= (default=${PACKAGEVERSION})]"
config-help :
	@cd ${PACKAGEROOT}/petsc/petsc-${PACKAGEVERSION} \
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
TGZURL = https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-${PACKAGEVERSION}.tar.gz

clean ::
	@/bin/rm -f *~
