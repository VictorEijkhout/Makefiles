################################################################
####
#### Makefile for NCVIEW installation
####
#### install on vista requires:
#### export cc=$( which ${TACC_CC} )
#### contemplate doing that in general?
####
################################################################

PACKAGE = NCVIEW
ABOUT = NetCDF viewer
URL = https://cirrus.ucsd.edu/ncview/
PACKAGEVERSION = 2.1.9
MODE = core
MODULES = netcdf udunits

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install
HASBIN = 1
NOLIB = 1
CONFIGUREFLAGS = \
    --with-nc-config=${TACC_NETCDF_BIN}/nc-config \
    -with-udunits2_incdir=${TACC_UDUNITS_INC} \
    -with-udunits2_libdir=${TACC_UDUNITS_LIB}
EXPLICITCOMPILERS=1

.PHONY: default_install
info ::
	@echo "make default_install : configure build"
default_install : 
	@export EXPLICITCOMPILERS=1 \
	 && make --no-print-directory configure build public \
	      PACKAGEVERSION=$(PACKAGEVERSION) JCOUNT=${JCOUNT}

DOWNLOADURL = https://cirrus.ucsd.edu/~pierce/ncview/ncview-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
