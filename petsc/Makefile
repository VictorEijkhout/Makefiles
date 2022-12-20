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
PACKAGEVERSION = 3.18.1
MAKEINCLUDES = ${STOCKYARD}/makefiles
MODULEROOT = ${HOME}/modulefiles

PETSC_SRC_DIR = ${PACKAGEROOT}/petsc/petsc-${PACKAGEVERSION}
info ::
	@echo "     [TACC_SYSTEM=system (default: ${TACC_SYSTEM})]"
ifeq "${TACC_SYSTEM}" "stampede2"
  TACC_SYSTEM=skx
else ifeq "${TACC_SYSTEM}" "frontera"
  TACC_SYSTEM=clx
else ifeq "${TACC_SYSTEM}" "frontera"
  TACC_SYSTEM=milan
endif
PETSCARCH = ${TACC_SYSTEM}-${TACC_FAMILY_COMPILER}-${TACC_FAMILY_MPI}
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
  MPI_FLAGS = --download-mpich=yes
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
CUDA_OPTIONS = --with-cuda=1 --with-cuda-dir=${TACC_CUDA_DIR} --with-cuda-gencodearch=70
  # --with-cuda=<bool>
  #      Indicate if you wish to test for cuda  current: 0
  # --with-cuda-dir=<dir>
  #      Indicate the root directory of the cuda installation
  # --with-cuda-pkg-config=<dir>
  #      Look for cuda using pkg-config utility optional directory to look in
  # --with-cuda-include=<dirs>
  #      Indicate the directory of the cuda include files
  # --with-cuda-lib=<libraries: e.g. [/Users/..../libcuda.a,...]>
  #      Indicate the cuda libraries
  # --with-cuda-gencodearch
  #      Cuda architecture for code generation, for example 70, (this may be used by external packages), use all to build a fat binary for distribution

nocuda :
	        echo "--with-cuda=yes --with-cudac=nvcc --with-cuda-include=${TACC_CUDA_DIR}/targets/x86_64-linux/include/ --with-cuda-lib=[${TACC_CUDA_DIR}/targets/x86_64-linux/lib/libcufft.so,libcusparse.so,libcusolver.so,libcublas.so,libcuda.so,libcudart.so] --with-precision=single --download-cusp=yes CUDAFLAGS=-arch=sm_70 " 

PACKAGES = 
PACKAGEINFO =

FORTRAN = 1
info ::
	@echo "    [FORTRAN=0/1 (default: ${FORTRAN})]"
ifeq "${FORTRAN}" "1"
  PACKAGES += --with-fortran-bindings=1 --with-fc=`which mpif90`
else 
  PACKAGES += --with-fortran-bindings=0 --with-fc=0
endif

FFTW = 0
info ::
	@echo "    [FFTW=0/1 (default: ${FFTW})]"
ifeq "${FFTW}${CUDA}" "10"
  PACKAGES += --with-fftw=1 --download-fftw=1
endif

MPE = 0
info ::
	@echo "    [MPE=0/1 (default: ${MPE})]"
ifeq "${MPE}" "1"
  PACKAGES += --with-mpe=1 --download-mpe=1
endif

PTSCOTCH = 0
info ::
	@echo "    [PTSCOTCH=0/1 (default: ${PTSCOTCH})]"
ifeq "${PTSCOTCH}${CUDA}" "10"
  PACKAGES += --with-ptscotch=1 --download-ptscotch=1
endif

SLEPC = 0
info ::
	@echo "    [SLEPC=0/1 (default: ${SLEPC})]"
ifeq "${SLEPC}" "1"
  PACKAGES += --with-slepc=1 --download-slepc=1
endif

SUPERLU = 0
info ::
	@echo "    [SUPERLU=0/1 (default: ${SUPERLU})]"
ifeq "${SUPERLU}" "1"
  PACKAGES +=  --with-superlu=1 --download-superlu=1 --with-superlu_dist=1 --download-superlu_dist=1 --with-cxx-dialect=C++11 --with-openmp=1
endif

ZOLTAN = 0
info ::
	@echo "    [ZOLTAN=0/1 (default: ${ZOLTAN})]"
ifeq "${ZOLTAN}${CUDA}" "10"
  PACKAGES += --with-zoltan=1 --download-zoltan=1
endif

DEBUG = 0
info ::
	@echo "     [DEBUG=0/1, default=${DEBUG}]"
TIME =

# ifeq "${SCALAR}${PRECISION}${CUDA}" "realdouble0"
#   PACKAGES += --download-sundials2=1
# endif

.PHONY: petscarch
petscarch :
	@export PETSC_EXT="" \
	  \
	  && if [ "${SCALAR}" = "complex" ] ; then \
	    export PETSC_EXT=$${PETSC_EXT}complex ; fi \
	  \
	  && if [ ${PRECISION} = "single" ] ; then \
	    export PETSC_EXIT=$${PETSC_EXT}single ; fi \
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
	  && if [ -z $${PETSC_EXT} ] ; then \
	      echo ${PETSCARCH} ; \
	    else \
	      echo ${PETSCARCH}-$${PETSC_EXT} ; \
	    fi

CEED = 0
info ::
	@echo "    [ CEED=0/1 (default: ${CEED})]"
ifeq "${CEED}"  "1"
  PACKAGES += --download-libceed=yes --download-libceed-commit=origin/main
endif

CHACO = 0
info ::
	@echo "    [ CHACO=0/1 (default: ${CHACO})]"
ifeq "${CHACO}"  "1"
  PACKAGES += --download-chaco=yes
endif

EIGEN = 0
info ::
	@echo "    [ EIGEN=0/1 (default: ${EIGEN})]"
ifeq "${EIGEN}" "1"
  PACKAGES += --download-eigen=yes
endif

HDF5 = 0
info ::
	@echo "    [ HDF5=0/1 (default: ${HDF5})]"
ifeq "${HDF5}" "1"
  PACKAGES += --download-hdf5=yes
endif

HYPRE = 0
FEI = 0
info ::
	@echo "    [ HYPRE=0/1 (default: ${HYPRE})]"
	@echo "    [FEI=0/1 (default: ${FEI})]"
ifeq "${HYPRE}" "1"
  PACKAGES += --download-hypre=yes
endif
ifeq "${FEI}" "1"
  PACKAGES += --download-hypre-configure-arguments=--with-fei
endif

KOKKOS = 0
info ::
	@echo "    [KOKKOS=0/1 (default: ${KOKKOS})]"
PACKAGEINFO += kokkos:${KOKKOS},
ifeq "${KOKKOS}" "1"
  PACKAGES += --with-kokkos=1 --download-kokkos=1
endif

MUMPS = 0
info ::
	@echo "    [MUMPS=0/1 (default: ${MUMPS})]"
ifeq "${MUMPS}" "1"
  PACKAGES += --with-mumps=1 --download-mumps=yes --download-scalapack=yes --download-blacs=yes --download-parmetis=yes --download-metis=yes
endif

P4P = 0
info ::
	@echo "    [ P4P=0/1 (default: ${P4P})]"
ifeq "${P4P}" "1"
  PACKAGES += --download-petsc4py=yes --with-python=1 --with-python-exec=python3
endif

X = 0
info ::
	@echo "    [ X=0/1, default=${X}]"
WITH_X = "--with-x=1"

configure :
	@echo && echo "there is no configure rule" && echo
include ${MAKEINCLUDES}/Make.install
NOINSTALL = 1
install ::
	@if [ ${CUDA} -eq 1 -a -z "${TACC_CUDA_DIR}" ] ; then \
	    echo && echo "Please load cuda module" && echo && exit 1 \
	    ; fi
	export scriptdir=`pwd` \
	&& export petscarch=`make --no-print-directory \
	                SCALAR=${SCALAR} PRECISION=${PRECISION} ALIGN=${ALIGN} \
	                CUDA=${CUDA} DEBUG=${DEBUG} CUSTEXT=${CUSTOMEXT} \
	                petscarch` \
	 && export logfile=petsc-$${petscarch}.log \
	 && pushd ${PETSC_SRC_DIR} \
	  && echo \
	  && echo "Configuring version=${PACKAGEVERSION} with arch=$${petscarch} in `pwd`" ; echo \
	  && rm -rf $${petscarch} \
	  && export PETSC_DIR=`pwd` \
	  && export I_MPI_FABRICS=shm:tmi \
	  && ( set -x && date && pwd \
	    && echo "Installing with packages: ${PACKAGEINFO}" \
	    && PETSC_ARCH=$$petscarch ${TIME} python3 ./configure \
	        ${MPI_FLAGS} \
	        --useThreads=0 --download-triangle \
	        --with-precision=${PRECISION} \
	        --with-scalar-type=${SCALAR} \
	        --with-np=4 \
	        --with-memalign=${ALIGN} \
	        CFLAGS="${OPTFLAGS}" CXXFLAGS="${OPTFLAGS}" FFLAGS="${OPTFLAGS}" \
	        `if [ ${X} -eq 1 ]    ; then echo ${WITH_X} ; fi ` \
	        `if [ ${CUDA} -eq 1 ] ; then echo ${CUDA_OPTIONS} ; fi ` \
	        --with-debugging=${DEBUG} \
	        "${BLAS_FLAGS}" \
	        ${PACKAGES} \
	    && ${TIME} make MAKE_NP=4 \
	  ) 2>&1 | tee $${scriptdir}/$${logfile} \
	 && popd \
	 && varsfile=vars-$${petscarch}.sh \
	 && make --no-print-directory petsc_varsfile VARSFILE=$$varsfile \
	 && chmod -R g+rX,o+rX ${PETSC_SRC_DIR}/$${petscarch}
info ::
	@echo "make petsc_varsfile"
.PHONY: petsc_varsfile
petsc_varsfile :
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
	 && echo && echo "Variable settings in: $${varsfile}" && echo \
	\
	 && source ${MAKEINCLUDES}/names.sh \
	 && export MODE=${MODE} && export MODULEROOT=${MODULEROOT} \
	 && export INSTALLPATH=${PACKAGEROOT}/petsc/petsc-${PACKAGEVERSION}/$$petscarch \
	 && setnames ${PACKAGEROOT} ${PACKAGE} ${PACKAGEVERSION} ${INSTALLEXT} \
	 && make --no-print-directory varsmodule

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
