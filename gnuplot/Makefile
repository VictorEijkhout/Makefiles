################################################################
####
#### Makefile for Gnuplot installation
####
################################################################

PACKAGE = GNUPLOT
URL = https://gnuplot.sourceforge.net \
    https://sourceforge.net/p/gnuplot/gnuplot-main/ci/master/tree/INSTALL
PACKAGEVERSION = 6.0.2
MODE = seq
MODULES = qt5

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.configure
CONFIGUREFLAGS = \
    --with-qt=qt5

HASBIN = 1
NOLIB = 1
NOINC = 1
CPTOINSTALLDIR = demo
include ${MAKEINCLUDES}/Make.install

.PHONY: default_install
info ::
	@echo "make default_install : configure build"
default_install : configure build

TGZURL = https://sourceforge.net/projects/gnuplot/files/gnuplot/${PACKAGEVERSION}/gnuplot-${PACKAGEVERSION}.tar.gz/download
include ${MAKEINCLUDES}/Make.download

