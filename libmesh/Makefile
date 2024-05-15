################################################################
####
#### Makefile for Libmesh installation
####
#### even the release seems to want to submodule
####
################################################################

PACKAGE = libmesh
URL = 
PACKAGEVERSION = git
## 1.7.1

MODE = mpi
MODULES = phdf5

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

TGZURL = https://github.com/libMesh/libmesh/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = https://github.com/libMesh/libmesh.git
SUBMODULE = 1
include ${MAKEINCLUDES}/Make.git

include ${MAKEINCLUDES}/Make.clean

