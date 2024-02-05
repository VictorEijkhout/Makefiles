################################################################
####
#### Makefile for CDO installation
####
################################################################

PACKAGE = CDO
URL = https://code.mpimet.mpg.de/projects/cdo/wiki/Cdo#Download-Compile-Install
PACKAGEVERSION = 2.3.0
MODE = seq
MODULES = netcdff netcdf

include ${MAKEINCLUDES}/Make.info

HASBIN = 1
include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install
default_install : configure build public

TGZURL = https://code.mpimet.mpg.de/attachments/download/29019/cdo-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
