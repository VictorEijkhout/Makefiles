################################################################
####
#### Makefile for PETSc local installation
####
################################################################

include ${STOCKYARD}/makefiles/Make.info

info ::
	@echo "make install "
	@echo "     [PACKAGEVERSION= (default=${PACKAGEVERSION})]"
	@echo "     [PRECISION=single/double (default=${PRECISION})]"
	@echo "     [CUSTOMEXT=customextension]"

PACKAGE = PETSC
PACKAGEROOT = ${STOCKYARD}
PACKAGEVERSION = 3.18.2
MAKEINCLUDES = ${STOCKYARD}/makefiles
MODULEROOT = ${HOME}/modulefiles
MODE = mpi

PETSC_SRC_DIR = ${PACKAGEROOT}/petsc/petsc-${PACKAGEVERSION}
info ::
	@echo "     [TACC_SYSTEM=system (default: ${TACC_SYSTEM})]"
OPTFLAGS=-O2 -fPIC
LOPTFLAGS=-O1 -fPIC

ifeq "${TACC_SYSTEM}" "nvarm"
  BLAS_FLAGS = --download-fblaslapack=yes
  MPI_FLAGS = --with-mpi-compilers=1 --with-cc=`which mpicc` --with-cxx=`which mpicxx` 
else ifeq "${TACC_SYSTEM}" "oneapi"
  BLAS_FLAGS=--with-blas-lapack-dir=${MKLROOT}
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

SCALAR = real
info ::
	@echo "    [SCALAR=real/complex (default: ${SCALAR})]"

ALIGN = 16
info ::
	@echo "     [ALIGN=nn (default: ${ALIGN})]"

CUDA = 0
info ::
	@echo "     [CUDA=0/1 (default: ${CUDA})]"
CUDA_OPTIONS = --with-cuda=1 --with-cuda-dir=${LMOD_CUDA_DIR} --with-cuda-gencodearch=70

DEBUG = 0
info ::
	@echo "     [DEBUG=0/1, default=${DEBUG}]"
TIME =

# ifeq "${SCALAR}${PRECISION}${CUDA}" "realdouble0"
#   PACKAGES += --download-sundials2=1
# endif

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
install ::
	@if [ ${CUDA} -eq 1 -a -z "${LMOD_CUDA_DIR}" ] ; then \
	    echo && echo "Please load cuda module" && echo && exit 1 \
	    ; fi \
	 && export scriptdir=`pwd` \
	 && export petscarch=`make --no-print-directory \
	                SCALAR=${SCALAR} PRECISION=${PRECISION} ALIGN=${ALIGN} \
	                CUDA=${CUDA} DEBUG=${DEBUG} CUSTOMEXT=${CUSTOMEXT} \
	                petscarch` \
	 && export logfile=petsc-$${petscarch}.log \
	 && pushd ${PETSC_SRC_DIR} \
	 && echo \
	 && echo "Configuring version=${PACKAGEVERSION} with arch=$${petscarch} in `pwd`" \
	 && echo \
	 && rm -rf $${petscarch} \
	 && export PETSC_DIR=`pwd` \
	 && export I_MPI_FABRICS=shm:tmi \
	 && ( set -x && date && pwd \
	    && echo "Installing with packages: ${PACKAGEINFO}" \
	    && export PETSC_ARCH=$${petscarch} \
	    && ${TIME} python3 ./configure \
	        ${MPI_FLAGS} \
	        --useThreads=0 --download-triangle \
	        --with-scalar-type=${SCALAR} \
	        --with-precision=${PRECISION} \
	        --with-memalign=${ALIGN} \
	        --with-np=4 \
	        CFLAGS="${OPTFLAGS}" CXXFLAGS="${OPTFLAGS}" FFLAGS="${OPTFLAGS}" \
	        `if [ ${X} -eq 1 ]    ; then echo ${WITH_X} ; fi ` \
	        `if [ ${CUDA} -eq 1 ] ; then echo ${CUDA_OPTIONS} ; fi ` \
	        --with-debugging=${DEBUG} \
	        "${BLAS_FLAGS}" \
	        ${PACKAGES} \
	    && ${TIME} make MAKE_NP=4 \
	  ) 2>&1 | tee $${scriptdir}/$${logfile} \
	 && popd \
	 && make --no-print-directory petsc_varsfile
public ::
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} && export MODULEROOT=${MODULEROOT} \
	 && export INSTALLPATH=${PACKAGEROOT}/petsc/petsc-${PACKAGEVERSION}/$$petscarch \
	 && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} ${INSTALLEXT} \
	 && chmod -R g+rX,o+rX ${INSTALLPATH}
info ::
	@echo "make petsc_varsfile"
.PHONY: petsc_varsfile
petsc_varsfile :
	@source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} && export MODULEROOT=${MODULEROOT} \
	 && export INSTALLPATH=${PACKAGEROOT}/petsc/petsc-${PACKAGEVERSION}/$$petscarch \
	 && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} ${INSTALLEXT} \
	 && make --no-print-directory varsmodule \
	 && modulefile=$$( make --no-print-directory varsmodulename ) \
	 && echo " .. adding petsc vars to: $${modulefile}" \
	 && ( \
	    echo "setenv( \"PETSC_DIR\",  \"${PACKAGEROOT}/petsc/petsc-${PACKAGEVERSION}\" )" ; \
	    echo "setenv( \"PETSC_ARCH\", \"$${petscarch}\" )" \
	    ) >>$${modulefile}
old_vars_rule :
	@export petscarch=$$( make --no-print-directory \
	                SCALAR=${SCALAR} PRECISION=${PRECISION} \
	                CUDA=${CUDA} DEBUG=${DEBUG} CUSTEXT=${CUSTOMEXT} \
	                petscarch ) \
	 && varsfile=vars-$${petscarch}.sh \
	 && installdir=${PACKAGEROOT}/petsc/petsc-${PACKAGEVERSION} \
	 && make --no-print-directory varsfile VARSFILE=$${varsfile} \
	        INSTALLDIR=$$installdir \
	        LIBDIR=$${installdir}/$$petscarch/lib \
	        INCDIR=$${installdir}/$$petscarch/include \
	        PKGCONFIGSET=$${petscarch}/lib/pkgconfig \
	 && echo "export PETSC_DIR=$$installdir" >>$${varsfile} \
	 && echo "export PETSC_ARCH=$$petscarch"                           >>$${varsfile} \
	 && echo && echo "Variable settings in: $${varsfile}" && echo

info ::
	@echo "make config-help config-view"
	@echo "     [PACKAGEVERSION= (default=${PACKAGEVERSION})]"
config-help :
	cd ${PETSC_SRC_DIR} ; \
	  ./configure -help
config-view :
	cat ${PETSC_SRC_DIR}/configure.log
#	less ${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf/configure.log
log-view :
	@cat ${PACKAGEROOT}/petsc/petsc-${PACKAGEVERSION}/configure.log

include ${STOCKYARD}/makefiles/Make.git
GITREPO = https://gitlab.com/petsc/petsc.git
BRANCH = main

include ${STOCKYARD}/makefiles/Make.download
TGZURL = https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-${PACKAGEVERSION}.tar.gz

clean ::
	/bin/rm -f *~
