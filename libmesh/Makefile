################################################################
####
#### Makefile for Libmesh installation
####
#### even the release seems to want to submodule
####
################################################################

PACKAGE = libmesh
URL = http://libmesh.github.io/ \
    https://github.com/libMesh/libmesh
PACKAGEVERSION = 1.8.0
PACKAGEVERSIONDOWNLOAD = 1.8.0-rc1
# 1.7.6
## git

MODE = mpi
# note: petsc installation needs to have hypre
MODULES = boost phdf5 petsc

include ${MAKEINCLUDES}/Make.info

info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build

##
## configure version
##
include ${MAKEINCLUDES}/Make.configure
CONFIGUREFLAGS=\
      --enable-silent-rules \
      --enable-unique-id \
      --disable-maintainer-mode \
      --disable-warnings \
    \
      --enable-hdf5 \
      --enable-petsc-hypre-required \
    --enable-dirichlet --enable-nodeconstraint \
    --enable-metaphysicl-required \
    --enable-parmesh \
    --enable-xdr-required \
    \
    --enable-fortran \
      --with-thread-model=openmp \
      --with-cxx-std-min=2014 \
    --with--boost=${TACC_BOOST_DIR} \
      --without-gdb-command \
      --with-methods=opt 
HASBIN = 1
PKGCONFIGLIB = pkgconfig
DEPENDSONCURRENTMINOR = phdf5 petsc

include ${MAKEINCLUDES}/Make.install

TGZURL = https://github.com/libMesh/libmesh/releases/download/v${PACKAGEVERSIONDOWNLOAD}/libmesh-${PACKAGEVERSIONDOWNLOAD}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/libMesh/libmesh.git
SUBMODULE = 1
include ${MAKEINCLUDES}/Make.git

include ${MAKEINCLUDES}/Make.clean

