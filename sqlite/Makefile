################################################################
####
#### Makefile for SQLITE installation
####
################################################################

PACKAGE = SQLITE
ABOUT = SQL database engine
URL = https://www.sqlite.org https://github.com/sqlite/sqlite
## version could be newer:
## I don't know how to download a specific version
PACKAGEVERSION = 3.49.2
# 3.51.0 is out
MODE = seq
MODULES = 

include ${MAKEINCLUDES}/Make.info

# https://www.sqlite.org/howtocompile.html
CONFIGUREFLAGS = \
    --all \
    --disable-tcl
HASBIN = 1
include ${MAKEINCLUDES}/Make.configure

EXTRABUILDTARGETS = all sqlite3.pc
PKGCONFIGLIB = pkgconfig
include ${MAKEINCLUDES}/Make.install

.PHONY: default_install
info ::
	@echo "make default_install"
default_install : configure build public 

DOWNLOADURL = https://github.com/sqlite/sqlite/archive/refs/tags/version-${PACKAGEVERSION}.tar.gz
# https://www.sqlite.org/src/tarball/sqlite.tar.gz
## ?r=release
include ${MAKEINCLUDES}/Make.download
