################################################################
####
#### Makefile for Blas++
####
################################################################

PACKAGE = BLASPP
PACKAGEROOT = ${STOCKYARD}
PACKAGEVERSION = git
MAKEINCLUDES = ${STOCKYARD}/makefiles
MODE = seq

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.cmake
CMAKEFLAGS = -Dblas=mkl
LIBDIR = lib64
# .cmake files are on the CMAKE_MODULE_PATH
CMAKE_MODULEPATH_ADD = lib64/blaspp

include ${MAKEINCLUDES}/Make.cbuild

GITREPO= git@github.com:icl-utk-edu/blaspp.git
include ${MAKEINCLUDES}/Make.git
