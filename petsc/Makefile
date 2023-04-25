################################################################
####
#### Makefile for PETSc local installation
####
################################################################

PACKAGE = PETSC
PACKAGEVERSION = 3.18.3
MODE = mpi

include ${MAKEINCLUDES}/Make.info
list_installations ::
	@source ${MAKEINCLUDES}/names.sh \
	 && setnames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	 && setmodulenames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	 && echo "installations in $${srcdir}: $$( cd $$srcdir && for d in * ; do if [ -d $${d} ] ; then echo $${d} ; fi ; done )"

info ::
	@echo "make build "
	@echo "     [CUSTOMEXT=customextension]"

JCOUNT = 4

PETSC_SRC_DIR = ${PACKAGEROOT}/petsc/petsc-${PACKAGEVERSION}
OPTFLAGS=-O2 -fPIC
LOPTFLAGS=-O1 -fPIC

ifeq "${TACC_SYSTEM}" "nvarm"
  BLAS_FLAGS = --download-fblaslapack=yes
  MPI_FLAGS = --with-mpi-compilers=1 --with-cc=`which mpicc` --with-cxx=`which mpicxx` 
else ifeq "${TACC_SYSTEM}" "macbookair"
  BLAS_FLAGS=
  MPI_FLAGS = --with-mpi-compilers=1 --with-cc=`which mpicc` --with-cxx=`which mpicxx` 
else
  BLAS_FLAGS=--with-blas-lapack-dir=${MKLROOT}
  MPI_FLAGS = --with-mpi-compilers=1 --with-cc=`which mpicc` --with-cxx=`which mpicxx` 
endif

info ::
	@echo "make installed : list installations"
.PHONY: installed
installed :
	@ls -d ${PETSC_SRC_DIR}/${PETSCARCH}*

PRECISION = double
info ::
	@echo "    [ PRECISION=single/double/quad (default: ${PRECISION})]"
.PHONY: precstring
precstring :
	@if [ "${PRECISION}" = "quad" ] ; then \
	    echo "__float128 --download-f2cblaslapack" \
	; else echo ${PRECISION} ; fi

SCALAR = real
info ::
	@echo "    [ SCALAR=real/complex (default: ${SCALAR})]"

INT = 32
info ::
	@echo "    [ INT=32/64 (default: ${INT}) ]"

ALIGN = 16
info ::
	@echo "    [ ALIGN=nn (default: ${ALIGN})]"

CUDA = 0
info ::
	@echo "     [CUDA=0/1 (default: ${CUDA})]"
CUDA_OPTIONS = --with-cuda=1 --with-cuda-dir=${LMOD_CUDA_DIR} --with-cuda-gencodearch=70

DEBUG = 0
info ::
	@echo "     [DEBUG=0/1, default=${DEBUG}]"
TIME =

TRIANGLE = 0
info ::
	@echo "    [ TRIANGLE=0/1 (default=${TRIANGLE}) ]"

##
## Petsc Architecture
##

.PHONY: petscarch
petscarch :
	@source ${MAKEINCLUDES}/names.sh \
	 && MODE=mpi setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} ${INSTALLEXT} \
	 && export PETSC_EXT="" \
	  \
	 && if [ "${SCALAR}" = "complex" ] ; then \
	    export PETSC_EXT=$${PETSC_EXT}complex ; fi \
	  \
	 && if [ ${PRECISION} = "single" ] ; then \
	    export PETSC_EXT=$${PETSC_EXT}single ; fi \
	  \
	 && if [ ${PRECISION} = "quad" ] ; then \
	    export PETSC_EXT=$${PETSC_EXT}single ; fi \
	  \
	 && if [ ${INT} = "64" ] ; then \
	    export PETSC_EXT=$${PETSC_EXT}i64 ; fi \
	 \
	 && if [ ${ALIGN} -ne 16 ] ; then \
	    export PETSC_EXT=$${PETSC_EXT}a${ALIGN} ; fi \
	  \
	 && if [ "${CUDA}" -eq 1 ] ; then \
	    export PETSC_EXT=$${PETSC_EXT}cuda ; fi \
	  \
	 && if [ "${DEBUG}" -eq 1 ] ; then \
	    export PETSC_EXT=$${PETSC_EXT}debug ; fi \
	  \
	 && if [ ! -z "${CUSTOMEXT}" ] ; then \
	    export PETSC_EXT=${CUSTOMEXT} ; fi \
	  \
	 && petscarch=$${taccsystemcode}-${LMOD_FAMILY_COMPILER}-${LMOD_FAMILY_MPI} \
	  && if [ -z $${PETSC_EXT} ] ; then \
	      echo $${petscarch} ; \
	    else \
	      echo $${petscarch}-$${PETSC_EXT} ; \
	    fi

####
#### Packages
####

PACKAGES = 
PACKAGEINFO =
include Makefile.packages

X = 0
info ::
	@echo "    [ X=0/1, default=${X}]"
WITH_X = "--with-x=1"

configure :
	@echo && echo "there is no configure rule" && echo
include ${MAKEINCLUDES}/Make.vars
NOINSTALL = 1
build ::
	@if [ ${CUDA} -eq 1 -a -z "${LMOD_CUDA_DIR}" ] ; then \
	    echo && echo "Please load cuda module" && echo && exit 1 \
	    ; fi
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} && export INSTALLROOT=${INSTALLROOT} \
	 && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} \
	 && setmodulenames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	 \
	 && export petscarch=`make --no-print-directory \
	                SCALAR=${SCALAR} PRECISION=${PRECISION} INT=${INT} \
	                ALIGN=${ALIGN} \
	                CUDA=${CUDA} DEBUG=${DEBUG} CUSTOMEXT=${CUSTOMEXT} \
	                petscarch` \
	 && ( \
	    pushd $$srcdir \
	     && echo "Configuring version=${PACKAGEVERSION} with arch=$${petscarch} in `pwd`" \
	     && echo \
	     && rm -rf $${petscarch} \
	     && export PETSC_DIR=`pwd` \
	     && export I_MPI_FABRICS=shm:tmi \
	     && date && pwd \
	     && echo "Installing with packages: ${PACKAGEINFO}" \
	     && export PETSC_ARCH=$${petscarch} \
	     && ${TIME} python3 ./configure \
	        ${MPI_FLAGS} \
	        --useThreads=0 --download-triangle=${TRIANGLE} \
	        --with-scalar-type=${SCALAR} \
	        $$( if [ "${INT}" = "64" ] ; then echo --with-64-bit-indices=1 ; fi ) \
	        --with-precision=$$( cd $$scriptdir \
	                && make --no-print-directory precstring PRECISION=${PRECISION} ) \
	        --with-memalign=${ALIGN} \
	        --with-np=${JCOUNT} \
	        CFLAGS="${OPTFLAGS}" CXXFLAGS="${OPTFLAGS}" FFLAGS="${OPTFLAGS}" \
	        `if [ ${X} -eq 1 ]    ; then echo ${WITH_X} ; fi ` \
	        `if [ ${CUDA} -eq 1 ] ; then echo ${CUDA_OPTIONS} ; fi ` \
	        --with-debugging=${DEBUG} \
	        "${BLAS_FLAGS}" \
	        ${PACKAGES} \
	     && ${TIME} make MAKE_NP=${JCOUNT} \
	     && popd \
	         && make --no-print-directory varsmodule \
	            MODULEVERSIONEXTRA=$${petscarch} \
	            PKGCONFIGSET=$${petscarch}/lib/pkgconfig \
	            EXTRAVARS="PETSC_DIR=$${srcdir} PETSC_ARCH=$${petscarch}" \
	  ) 2>&1 | tee $${scriptdir}/install-$${petscarch}.log
## $${installlog}

public ::
	@echo "public rule disabled"
nopublic :
	@source ${MAKEINCLUDES}/names.sh \
	 && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} ${INSTALLEXT} \
	 && setmodulenames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	 && export MODE=${MODE} && export MODULEROOT=${MODULEROOT} \
	 && export INSTALLPATH=${PACKAGEROOT}/petsc/petsc-${PACKAGEVERSION}/$$petscarch \
	 && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} ${INSTALLEXT} \
	 && chmod -R g+rX,o+rX $$installdir
info ::
	@echo "make petsc_varsfile"
.PHONY: petsc_varsfile
petsc_varsfile :
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} && export MODULEROOT=${MODULEROOT} \
	 && export INSTALLPATH=${PACKAGEROOT}/petsc/petsc-${PACKAGEVERSION}/$$petscarch \
	 && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} ${INSTALLEXT} \
	 && setmodulenames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" "${PACKAGEBASENAME}" \
	 && setdirlognames "${PACKAGEROOT}" "${PACKAGE}" "${PACKAGEVERSION}" "${INSTALLEXT}" \
	 && make --no-print-directory varsmodule \
	 && modulefile=$$( make --no-print-directory varsmodulename ) \
	 && echo " .. adding petsc vars to: $${modulefile}" \
	 && ( \
	    echo "setenv( \"PETSC_DIR\",  \"${PACKAGEROOT}/petsc/petsc-${PACKAGEVERSION}\" )" ; \
	    echo "setenv( \"PETSC_ARCH\", \"$${petscarch}\" )" \
	    ) >>$${modulefile}

info ::
	@echo "make config-help config-view"
	@echo "     [PACKAGEVERSION= (default=${PACKAGEVERSION})]"
config-help :
	cd ${PETSC_SRC_DIR} ; \
	  ./configure -help
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
	/bin/rm -f *~
