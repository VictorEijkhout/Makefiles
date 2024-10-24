################################################################
####
#### Makefile for BISON installation
####
################################################################

PACKAGE = BISON
URL = https://www.gnu.org/software/bison/ \
    https://ftp.gnu.org/gnu/bison/
PACKAGEVERSION = 3.8
MODE = core
MODULES = 

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.configure
NOLIB = 1
HASBIN = 1
include ${MAKEINCLUDES}/Make.install

.PHONY: default_install
default_install : configure build 

TXZURL = https://ftp.gnu.org/gnu/bison/bison-${PACKAGEVERSION}.tar.xz
include ${MAKEINCLUDES}/Make.download

include ${MAKEINCLUDES}/Make.clean
