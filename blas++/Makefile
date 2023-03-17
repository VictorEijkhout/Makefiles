################################################################
####
#### Makefile for Blas++
####
################################################################

PACKAGE = BLASPP
PACKAGEVERSION = git
MODE = seq

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
CMAKEFLAGS = -Dblas=mkl
LIBDIR = lib64
# .cmake files
CMAKE_PREFIXPATH_SET = lib64/blaspp

include ${MAKEINCLUDES}/Make.cbuild

GITREPO= git@github.com:icl-utk-edu/blaspp.git
include ${MAKEINCLUDES}/Make.git
