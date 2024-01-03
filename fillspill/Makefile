################################################################
####
#### Makefile for FILLSPILLMERGE installation
####
################################################################

PACKAGE = FILLSPILLMERGE
PACKAGEVERSION = git
# 1.0.0
URL = https://github.com/r-barnes/Barnes2020-FillSpillMerge
MODE = seq
MODULES = gdal jsonc proj

include ${MAKEINCLUDES}/Make.info

PKGCONFIGLIB = pkgconfig
CMAKEFLAGS = \
    -D USE_GDAL=ON
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
HASBIN = 1
NOLIB = 1
#NOINSTALL = 1
.PHONY: default_install
default_install : configure build public

TGZURL = https://github.com/r-barnes/Barnes2020-FillSpillMerge/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = git@github.com:r-barnes/Barnes2020-FillSpillMerge.git
SUBMODULE = submodules
include ${MAKEINCLUDES}/Make.git
