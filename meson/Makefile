################################################################
####
#### Makefile for Meson installation
####
################################################################

PACKAGE = meson
ABOUT = The Meson Build System
PACKAGEVERSION = 1.4.1
# 1.10.0 is out
URL = https://github.com/mesonbuild/meson
MODE = core
MODULE = 

include ${MAKEINCLUDES}/Make.info

PKGCONFIGLIB = pkgconfig
CMAKEPREFIXPATHSET = 1
NOLIB = 1
HASBIN = 1
include ${MAKEINCLUDES}/Make.vars
.PHONY: default_install
info ::
	@echo "make default_install : pip in the source directory"
default_install : 
	@source ${MAKEINCLUDES}/names.sh \
	 && setnames \
	    "${PACKAGE}"  "${PACKAGEVERSION}"  "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	    "${INSTALLPATH}"  "${INSTALLROOT}"  "${INSTALLEXT}"  "${INSTALLVARIANT}" \
	    "${HOMEDIR}"  "${BUILDDIRROOT}"  "${MODE}" \
	    "${PREFIXOPTION}"  "${PREFIXEXTRA}" \
	 && requirenonzero srcdir \
	 && requirenonzero prefixdir \
	 && ( cd $${srcdir} && pip3 install --target=$${prefixdir} . ) \
	 && make --no-print-directory varsmodule \
	        EXTRAINSTALLPATHS=PYTHONPATH=

DOWNLOADURL = https://github.com/mesonbuild/meson/archive/refs/tags/${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = https://github.com/mesonbuild/meson.git

