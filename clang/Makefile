################################################################
####
#### Makefile for LLVM/CLANG installation
####
################################################################

PACKAGE = clang
URL = https://clang.llvm.org/get_started.html
PACKAGEVERSION = git
MODE = seq

include ${MAKEINCLUDES}/Make.info

NOLIB = 1
HASBIN = 1
CMAKESOURCE = llvm
CMAKEFLAGS = \
    -D LLVM_ENABLE_PROJECTS=clang \
    -D CMAKE_BUILD_TYPE=Release
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
default_install : configure build public

GITREPO = https://github.com/llvm/llvm-project.git
CLONEARGS = --depth=1
include ${MAKEINCLUDES}/Make.git
