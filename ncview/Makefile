################################################################
####
#### Makefile for NCVIEW installation
####
################################################################

PACKAGE = NCVIEW
URL = https://cirrus.ucsd.edu/ncview/
PACKAGEVERSION = 2.1.9
MODE = seq
MODULES = netcdf
 
include ${MAKEINCLUDES}/Make.info

HASBIN = 1
include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install
NOLIB = 1
CONFIGUREFLAGS = \
    --with-nc-config=${TACC_NETCDF_BIN}/nc-config
default_install : 
	@export TACC_CC=$$( which ${TACC_CC} ) \
	 && make --no-print-directory configure build public

TGZURL = https://cirrus.ucsd.edu/~pierce/ncview/ncview-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
