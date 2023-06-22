################################################################
####
#### Makefile for SQLITE installation
#### https://www.sqlite.org/download.html
####
################################################################

PACKAGE = SQLITE
## version could be newer:
## I don't know how to download a specific version
PACKAGEVERSION = 3.43.0
MODE = seq
MODULES = 

include ${MAKEINCLUDES}/Make.info

# https://www.sqlite.org/howtocompile.html
CONFIGUREFLAGS = --disable-tcl
include ${MAKEINCLUDES}/Make.configure

EXTRABUILDTARGETS = all sqlite3.pc
include ${MAKEINCLUDES}/Make.install

TGZURL = https://www.sqlite.org/src/tarball/sqlite.tar.gz
## ?r=release
include ${MAKEINCLUDES}/Make.download
