################################################################
####
#### Makefile for PETSc local installation
####
################################################################

PACKAGE = PETSC
URL = https://petsc.org/release/install/download/
PACKAGEVERSION = 3.23.4
MODE = mpi
UPSTREAM = chyps dealii libmesh mfem p4est precice

include ${MAKEINCLUDES}/Make.info
list_installations ::
	@source ${MAKEINCLUDES}/names.sh \
	 && echo "installations in $${srcdir}: $$( cd $$srcdir && for d in * ; do if [ -d $${d} ] ; then echo $${d} ; fi ; done )"
modules ::
	@\
	  case "${TACC_SYSTEM}" in \
	    ( "vista" )      make --no-print-directory moduletest MODULE=nvpl ;; \
	    ( * ) \
	      if [ "${TACC_FAMILY_COMPILER}" = "gcc" ] ; then \
	        make --no-print-directory moduletest MODULE=mkl \
	      ; fi ;; \
	  esac

info ::
	@echo "make build "
	@echo "     [CUSTOMEXT=customextension]"

JCOUNT = 4

OPTFLAGS=-O2 -fPIC
LOPTFLAGS=-O1 -fPIC
ifeq "${TACC_FAMILY_COMPILER}" "gcc"
  OPTFLAGS +=  -fabi-version=16
endif

ifeq "${TACC_SYSTEM}" "vista"
  BLAS_FLAGS = "\
    --with-lapack-lib=[${TACC_NVPL_LIB}/libnvpl_lapack_lp64_seq.so,libnvpl_lapack_core.so] \
    --with-blas-lib=[${TACC_NVPL_LIB}/libnvpl_blas_lp64_seq.so,libnvpl_blas_core.so] \
    "
else ifeq "${TACC_SYSTEM}" "macbookair"
  BLAS_FLAGS= --download-fblaslapack=yes
else
  BLAS_FLAGS=--with-blas-lapack-dir=${TACC_MKL_DIR}
endif

DOWNLOADCMAKE =
## --download-cmake


CUDA = 0
info ::
	@echo "     [CUDA=0/1 (default: ${CUDA})]"

TIME =

##
## Petsc Architecture
##

info ::
	@echo "================ PETSc specific options"

.PHONY: petscarch petscshortarch
petscshortarch :
	@export PETSC_EXT="" \
	 \
	 && if [ "${NOHDF5}" = "1" ] ; then \
	    export PETSC_EXT=$${PETSC_EXT}nohdf5 ; fi \
	 && if [ "${SCALAR}" = "complex" ] ; then \
	    export PETSC_EXT=$${PETSC_EXT}complex ; fi \
	 && if [ "${PRECISION}" = "single" ] ; then \
	    export PETSC_EXT=$${PETSC_EXT}single ; fi \
	 && if [ "${PRECISION}" = "quad" ] ; then \
	    export PETSC_EXT=$${PETSC_EXT}single ; fi \
	 && if [ "${FORTRAN}" = "08" ] ; then \
	    export PETSC_EXT=$${PETSC_EXT}f08 ; fi \
	 && if [ "${INT}" = "64" ] ; then \
	    export PETSC_EXT=$${PETSC_EXT}i64 ; fi \
	 && if [ ! -z "${ALIGN}" ] ; then \
	    export PETSC_EXT=$${PETSC_EXT}a${ALIGN} ; fi \
	 && if [ "${CUDA}" = "1" ] ; then \
	    export PETSC_EXT=$${PETSC_EXT}cuda ; fi \
	 && if [ "${DEBUGGING}" = "1" ] ; then \
	    export PETSC_EXT=$${PETSC_EXT}debug ; fi \
	  \
	 && if [ ! -z "${CUSTOMEXT}" ] ; then \
	    export PETSC_EXT=${CUSTOMEXT} ; fi \
	 && echo $${PETSC_EXT}
petscarch :
	@source ${MAKEINCLUDES}/names.sh \
	 && systemnames && compilernames \
	 && notetoself "TODO dir names are no needed?" \
	 && compilernames && requirenonzero compilercode && requirenonzero mpicode \
	 && systemnames \
	 && petscarch=$${taccsystemcode}-$${compilercode}$${compilershortversion}-$${mpicode} \
	 && petscext=$$( make --no-print-directory petscshortarch ) \
	 && if [ -z "$${petscext}" -o "$${petscext}" = "vanilla" ] ; then \
	      echo ${PACKAGEVERSION}; \
	    else \
	      echo ${PACKAGEVERSION}-$${petscext} ; \
	    fi
.PHONY: configurelog
info ::
	@echo "make configurelog : display current configure.log"
configurelog :
	@source ${MAKEINCLUDES}/names.sh \
	 && systemnames && compilernames \
	 && notetoself "TODO dir names are no needed?" \
	 && setnames \
	    "${PACKAGE}"  "${PACKAGEVERSION}"  "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	 && echo "$${srcdir}/configure.log"
oldloglog :
	 && cd $${srcdir} \
	 && ( ls -l $${srcdir}/configure.log && echo ">>>>>>>>>>>>>>>>" && cat ./configure.log )

####
#### Packages
####

include Makefile.packages

info ::
	@echo "make default_install (there are no configure/install rules)"
.PHONY: default_install
default_install :: modules
	@if [ "${CUDA}" -eq 1 -a -z "${TACC_CUDA_DIR}" ] ; then \
	    echo && echo "Please load cuda module" && echo && exit 1 \
	 ; else \
	    taccver=${TACC_CUDA_VER} \
	     && if [ ! -z "$${taccver}" ] ; then \
	          vermin=$${taccver##*.} \
	           && if [ $$vermin -lt 8 ] ; then \
	                echo "CUDA version needs to be at least 12.8" \
	                 && exit 1 \
	              ; fi \
	        ; fi \
	 ; fi
	@source ${MAKEINCLUDES}/names.sh \
	 && export INSTALLROOT=${INSTALLROOT} \
	 && export petscarch=`make --no-print-directory \
	                SCALAR=${SCALAR} PRECISION=${PRECISION} INT=${INT} \
	                ALIGN=${ALIGN} \
	                CUDA=${CUDA} DEBUGGING=${DEBUGGING} CUSTOMEXT=${CUSTOMEXT} \
	                petscarch` \
	 && echo && echo "Deduce petc arch <<$$petscarch>>" && echo \
	 && setnames \
	    "${PACKAGE}"  "${PACKAGEVERSION}"  "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	    "${INSTALLPATH}"  "${INSTALLROOT}"  "${INSTALLEXT}"  "$${petscarch}" \
	    "${HOMEDIR}"  "${BUILDDIRROOT}"  "${MODE}" \
	    "${PREFIXOPTION}"  "${PREFIXEXTRA}" \
	 && echo "Deduce install dir <<$$prefixdir>>" && echo \
	 && source ${MAKEINCLUDES}/compilers.sh && setcompilers && reportcompilers \
	 \
	 && echo "Finding options and settings" \
	 && export package_options="$$( make --no-print-directory package_options \
	      AMGX=${AMGX} CHACO=${CHACO} EIGEN=${EIGEN} FFTW=${FFTW} HDF5=${HDF5} \
	      HYPRE=${HYPRE} KOKKOS=${KOKKOS} LIBCEED=${LIBCEED} MPE=${MPE} \
	      MUMPS=${MUMPS} PTSCOTCH=${PTSCOTCH} SLEPC=${SLEPC} \
	      STRUMPACK=${STRUMPACK} SUPERLU=${SUPERLU} SUPERLUDIST=${SUPERLUDIST} \
	      TRIANGLE=${TRIANGLE} \
	      \
	      PETSC4PY=${PETSC4PY} SLEPC4PY=${SLEPC4PY} \
	      )" \
	 && export param_settings="$$( make --no-print-directory param_settings \
	        DEBUGGING=${DEBUGGING} INT=${INT} PRECISION=${PRECISION} SCALAR=${SCALAR} \
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
	     && echo "Kludge: adding <<$${prefixdir}/lib>> to pythonpath" \
	     && export PYTHONPATH=${PYTHONPATH}:$${prefixdir}/lib \
	     && echo "PYTHONPATH=$${PYTHONPATH}" \
	     && echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VLE configure" \
	     && if [ ${CUDA} -eq 1 ] ; then \
	          export NVCC_APPEND_FLAGS=-allow-unsupported-compiler \
	           && echo "Cuda environment:" && env | grep -i cuda \
	        ; fi \
	     && cmdline="${TIME} python3 ./configure --prefix=$${prefixdir} \
	        ${DOWNLOADCMAKE} \
	        COPTFLAGS=\"${OPTFLAGS}\" \
	        COPTXXFLAGS=\"${OPTFLAGS}\" \
	        FOPTFLAGS=\"${OPTFLAGS}\" \
	        `if [ ${CUDA} -eq 1 ] ; then \
	                echo \"--with-cuda=1 --with-cuda-dir=${TACC_CUDA_DIR} \
	                  --with-cudac=nvcc \
	                  --with-cuda-gencodearch=70\" \
	          ; else echo --with-cudac=0 ; fi ` \
	        "${BLAS_FLAGS}" \
i	        $${package_options} $${param_settings} \
	        " \
	     && echo "Cmdline=$$cmdline" && retcode=0 && eval $$cmdline || retcode=$$? \
	    \
	     && if [ $${retcode} -gt 0 ] ; then \
	        echo "Build failed for $${PETSC_ARCH}" && exit 1 ; fi \
	     && echo "................................ VLE configure" \
	     \
	     && echo "PYTHONPATH=$${PYTHONPATH}" \
	     && echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VLE make" \
	     && make PETSC_DIR=`pwd` PETSC_ARCH=$${PETSC_ARCH} MAKE_NP=${JCOUNT} all \
	     && make PETSC_DIR=`pwd` PETSC_ARCH=$${PETSC_ARCH} MAKE_NP=${JCOUNT} install \
	     && echo " .. open up prefix dir" \
	     && chmod -R g+rX,o+rX $${prefixdir} \
	     && echo " .. sys error disable" \
	     && sed -e 's/#error/#warning/' -i $${prefixdir}/include/petscsys.h \
	     && echo "................................ VLE make" \
	  ) 2>&1 | tee $${installlog}
include ${MAKEINCLUDES}/Make.vars
include ${MAKEINCLUDES}/Make.public

info ::
	@echo "make installation_public : open installation (also done as part of build)"
.PHONY: installation_public
## build :: installation_public
## I CAN NOT FIGURE OUT THE RIGHT WAY TO GENERATE THE NAME OF THE INSTSALL DIRECTORY
installation_public :
	@ export petscarch=`make --no-print-directory \
	         SCALAR=${SCALAR} PRECISION=${PRECISION} INT=${INT} \
	         ALIGN=${ALIGN} \
	         CUDA=${CUDA} DEBUGGING=${DEBUGGING} CUSTOMEXT=${CUSTOMEXT} \
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
	@echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VLE make module" \
	 && source ${MAKEINCLUDES}/names.sh \
	 && export INSTALLROOT=${INSTALLROOT} \
	 && export petscarch=`make --no-print-directory \
	                SCALAR=${SCALAR} PRECISION=${PRECISION} INT=${INT} \
	                ALIGN=${ALIGN} \
	                CUDA=${CUDA} DEBUGGING=${DEBUGGING} CUSTOMEXT=${CUSTOMEXT} \
	                petscarch` \
	 && setnames \
	    "${PACKAGE}"  "${PACKAGEVERSION}"  "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	    "${INSTALLPATH}"  "${INSTALLROOT}"  "${INSTALLEXT}"  "$${petscarch}" \
	    "${HOMEDIR}"  "${BUILDDIRROOT}"  "${MODE}" \
	    "${PREFIXOPTION}"  "${PREFIXEXTRA}" \
	 && echo "makemodule: deduce install dir <<$$prefixdir>>" \
	 \
	 && shortversion=${PACKAGEVERSION} && shortversion=$${shortversion%.*} \
	 && petscshortarch=$$( make --no-print-directory petscshortarch ) \
	 && echo "making module with version=$${shortversion} extra=<<$${petscshortarch}>>" \
	 && ( make --no-print-directory varsmodule  \
	        MODULEROOT=${MODULEROOT} MODULEPATH=${MODULEPATH} \
	        INSTALLPATH="$${prefixdir}" \
	        PACKAGEMODULEVERSION=$${shortversion} \
	        MODULEVERSIONEXTRA=$${petscshortarch} \
	        PKGCONFIGLIB=pkgconfig \
	        EXTRAVARS="PETSC_DIR=$${prefixdir} PETSC_ARCH=$${petscarch}" \
	        EXTRAINSTALLPATHS="PYTHONPATH=lib" \
	    ) 2>&1 | tee -a $${installlog}

##
## installation modes
##

info ::
	@echo "make biginstall [ EXT=i64debug and such ]"
	@echo "    build specific configuration"
	@echo "    packages specified by HDF5=1 and such"
.PHONY: listextensions biginstall allinstall
EXTENSIONS = \
    "" debug i64 i64-debug single single-debug \
    complex complex-debug complexi64 complexi64-debug \
    nohdf5 i64nohdf5
info ::
	@echo "make allinstall [ JCOUNT= PACKAGEVERSION= ]"
	@echo "    biginstall for all extensions except cuda"
listextensions :
	@echo "vanilly ${EXTENSIONS}"
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
	   && echo "Compiler issues" \
	   && FORTRAN=${FORTRAN} && SLEPC=${SLEPC} \
	   && if [ -z "${TACC_FAMILY_COMPILER}}" ] ; then \
	        echo "ERROR compiler does not define the usual variables" && exit 1 ; fi \
	   && case ${TACC_FAMILY_COMPILER} in \
	        ( gcc ) echo " .. no FORTRAN for gcc" && echo FORTRAN=0 ;; \
	        esac \
	   && echo "Selecting basic options for extension=$${ext}:" \
	   && case $${ext} in ( *debug* ) DEBUGGING=1 ;; ( * ) DEBUGGING=0 ;; esac \
	   && case $${ext} in ( *i64* ) INT=64 ;; ( * ) INT=32 ;; esac \
	   && case $${ext} in ( *single* ) PRECISION=single ;; ( * ) PRECISION=double ;; esac \
	   && case $${ext} in ( *complex* ) SCALAR=complex ;; ( * ) SCALAR=real ;; esac \
	   \
	   && echo "Selecting packages" \
	   && AMGX=${AMGX} && CHACO=${CHACO} && CUDA=${CUDA} && EIGEN=${EIGEN} && FFTW=${FFTW} && HDF5=${HDF5} && HYPRE=${HYPRE} && METIS=${METIS} \
	   && MUMPS=${MUMPS} && PARMETIS=${PARMETIS} && PTSCOTCH=${PTSCOTCH} \
	   && STRUMPACK=${STRUMPACK} \
	   && SUPERLU=${SUPERLU} && SUPERLUDIST=${SUPERLUDIST} \
	   && SCALAPACK=${SCALAPACK} && ZOLTAN=${ZOLTAN} \
	   && compiler_major_version=${TACC_FAMILY_COMPILER_VERSION} \
	       && compiler_major_version=$${compiler_major_version%%.*} \
	   && echo "Incompatibilities mode/package:" \
	   && if [ "$${CUDA}" = "1" ] ; then \
	        echo " .. disable STRUMPACK: don't work with cuda" && STRUMPACK=0 \
	         && echo " .. disable PARMETIS, PTSCOTCH: something about 64 bit blas/lapack fftw" \
	         && PARMETIS=0 && PTSCOTCH=0 && METIS=0 && EIGEN=0 && CHACO=0 \
	         && TRIANGLE=0 \
	         && echo " .. disable AMGX: don't work with cuda, but should; investigate!" && AMGX=0 \
	         && echo " .. more cuda disables" && HDF5=0 \
	      ; elif [ "${TACC_FAMILY_COMPILER}" = "gcc" -a "$${compiler_major_version}" = "15" ] ; then \
	        echo " .. disable gcc15 stuff: CHACO HYPRE SUPERLUDIST" \
	         && CHACO=0 && HYPRE=0 && SUPERLUDIST=0 \
	      ; else \
	        echo " .. disable AMGX because no cuda" && AMGX=0 \
	      ; fi \
	   && if [ "$${FORTRAN}" = "0" ] ; then \
	        echo " .. disable MUMPS, STrumpack, Scalapack because no Fortran" \
	         && MUMPS=0 && SCALAPACK=0 && STRUMPACK=0 \
	      ; elif [ "$${PARMETIS}" = "1" ] ; then \
	          MUMPS=1 \
	      ; fi \
	   && AMGX=0 \
	   && echo "Incompatibilities mode/package" \
	   && case $${ext} in \
	        ( *nohdf5* ) HDF5=0 ; PHDF5=0 ;; esac \
	   && case $${ext} in \
	        ( *single* ) echo " .. disable fftw hdf5, Hypre, Mumps because single" \
	                && FFTW=0 && HDF5=0 && PHDF5=0 && HYPRE=0 && MUMPS=0 ;; esac \
	   && case $${ext} in \
	        ( *complexsingle* ) echo " .. disable superlu_dist because complex single" \
	                && SUPERLUDIST=0 ;; esac \
	   && case $${ext} in \
	        ( *i64* ) CHACO=0 && PARMETIS=0 && SUPERLU=0 SUPERLUDIST=0 ;; esac \
	   && case $${ext} in \
	        ( *complex* ) \
	                echo " .. disable Hypre, Metis, Mumps because complex" \
	                 && HYPRE=0 \
	      ;; esac \
	   \
	   && if [ "$${METIS}" = "0" ] ; then PARMETIS=0 && PTSCOTCH=0 ; fi \
	   && if [ "$${MUMPS}" = "0" ] ; then SCALAPACK=0 ; fi \
	   \
	   && echo "Start configure/build at $$( date )" \
	   && cmdline="make --no-print-directory default_install makemodule \
	      PACKAGEVERSION=${PACKAGEVERSION} JCOUNT=${JCOUNT} \
	      DEBUGGING=$${DEBUGGING} INT=$${INT} PRECISION=$${PRECISION} SCALAR=$${SCALAR} \
	      $$( if [ ! -z ${CUSTOMEXT} ] ; then echo CUSTOMEXT=${CUSTOMEXT} ; fi ) \
	      \
	      AMGX=$${AMGX} CHACO=$${CHACO} EIGEN=$${EIGEN} FFTW=$${FFTW} \
	      HDF5=$${HDF5} HYPRE=$${HYPRE} METIS=$${METIS} \
	      MUMPS=$${MUMPS} PARMETIS=$${PARMETIS} PTSCOTCH=$${PTSCOTCH} \
	      STRUMPACK=$${STRUMPACK} \
	      SUPERLU=$${SUPERLU} SUPERLUDIST=$${SUPERLUDIST} \
	      TRIANGLE=$${TRIANGLE} \
	      SCALAPACK=$${SCALAPACK} ZOLTAN=$${ZOLTAN} \
	      \
	      CUDA=${CUDA} FORTRAN=$${FORTRAN} SLEPC=$${SLEPC} \
	      PETSC4PY=${PETSC4PY} SLEPC4PY=${SLEPC4PY} \
	      " \
	   && if [ ! -z "${OUTERLOG}" ] ; then \
	        echo "Big install: $$cmdline" | tee -a "${OUTERLOG}" \
	      ; else \
	        echo "Big install: $$cmdline"\
	      ; fi \
	   && installlog=install.log \
	   && set -e && set -o pipefail && eval $${cmdline} 2>&1 | tee $${installlog} \
	   && echo && echo "See: $${installlog}" && echo

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
	                CUDA=${CUDA} DEBUGGING=${DEBUGGING} CUSTOMEXT=${CUSTOMEXT} \
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
	@/bin/rm -f *job.o[0-9]* archs-* *setnames petscinstalljob* 
