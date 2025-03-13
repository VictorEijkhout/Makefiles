################################################################
####
#### Makefile for CCCL installation
####
#### https://github.com/NVIDIA/cccl/issues/4101
####
################################################################

PACKAGE = CCCL
PACKAGEVERSION = 2.8.0
URL = https://github.com/NVIDIA/cccl
MODE = seq
MODULES = cuda

include ${MAKEINCLUDES}/Make.info

CMAKEFLAGS = \
    -D ARG_ALLOW_EXTERNAL=OFF
include ${MAKEINCLUDES}/Make.cmake
PKGCONFIGLIB = pkgconfig
CMAKEPREFIXPATHSET = 1
## NOLIB = 1
include ${MAKEINCLUDES}/Make.cbuild
.PHONY: default_install
info ::
	@echo "make default_install : configure build"
default_install : configure build public

TGZURL = https://github.com/NVIDIA/cccl/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = https://github.com/NVIDIA/cccl.git

include ${MAKEINCLUDES}/Make.clean

