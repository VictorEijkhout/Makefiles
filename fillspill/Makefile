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

include ${MAKEINCLUDES}/Make.info

## PKGCONFIGLIB = pkgconfig
## NOLIB = 1
include ${MAKEINCLUDES}/Make.cmake
include ${MAKEINCLUDES}/Make.cbuild
default_install : configure build public

TGZURL = https://github.com/r-barnes/Barnes2020-FillSpillMerge/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

GITREPO = git@github.com:r-barnes/Barnes2020-FillSpillMerge.git
SUBMODULE = 1
include ${MAKEINCLUDES}/Make.git
