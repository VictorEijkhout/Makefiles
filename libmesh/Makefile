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
PACKAGEVERSION = 1.8.0-rc1
# 1.7.6
## git

MODE = mpi
MODULES = phdf5 petsc

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
      --disable-warnings \
      --with-thread-model=openmp \
      --disable-maintainer-mode \
      --enable-hdf5 \
      --enable-petsc-hypre-required \
      --enable-metaphysicl-required \
      --enable-xdr-required \
      --with-cxx-std-min=2014 \
      --without-gdb-command \
      --with-methods=opt 

include ${MAKEINCLUDES}/Make.install

TGZURL = https://github.com/libMesh/libmesh/releases/download/v${PACKAGEVERSION}/libmesh-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/libMesh/libmesh.git
SUBMODULE = 1
include ${MAKEINCLUDES}/Make.git

include ${MAKEINCLUDES}/Make.clean

