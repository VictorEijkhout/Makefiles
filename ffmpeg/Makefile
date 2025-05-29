################################################################
####
#### Makefile for Ffmpeg installation
####
################################################################

PACKAGE = FFMPEG
URL = https://ffmpeg.org/download.html
PACKAGEVERSION = git

MODE = core
MODULES = 

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.configure
CONFIGUREFLAGS = \
    --disable-x86asm
include ${MAKEINCLUDES}/Make.install
HASBIN=1
NOINC = 1
CPTOINSTALLDIR = doc fftools libav*

info ::
	@echo "make default_install : configure build"
.PHONY: default_install
default_install : configure build public

##
## download
##
BZ2URL = https://ffmpeg.org/releases/ffmpeg-snapshot-git.tar.bz2
include ${MAKEINCLUDES}/Make.download

##
## git clone
##
GITREPO = https://git.ffmpeg.org/ffmpeg.git
include ${MAKEINCLUDES}/Make.git

##
## Clean
##
include ${MAKEINCLUDES}/Make.clean
