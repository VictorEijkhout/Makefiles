################################################################
####
#### Makefile for CDO installation
####
################################################################

PACKAGE = CDO
URL = https://code.mpimet.mpg.de/projects/cdo/wiki/Cdo#Download-Compile-Install
PACKAGEVERSION = 2.5.0
## note: version is not enough for the download url
MODE = seq
MODULES = netcdff netcdf

include ${MAKEINCLUDES}/Make.info

HASBIN = 1
include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install
NOINC =1
HASBIN = 1
NOLIB = 1
.PHONY: default_install
info ::
	@echo "make default_install: configure build"
default_install : configure build public

TGZURL = https://code.mpimet.mpg.de/attachments/download/29786/cdo-${PACKAGEVERSION}.tar.gz
## 2.3.0: TGZURL = https://code.mpimet.mpg.de/attachments/download/29019/cdo-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
