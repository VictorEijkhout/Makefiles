################################################################
####
#### Makefile for TMUX installation
####
################################################################

PACKAGE = TMUX
ABOUT = Terminal Multi pleXer
URL = https://github.com/tmux/tmux
PACKAGEVERSION = 3.5a
MODE = core
MODULES = 

include ${MAKEINCLUDES}/Make.info

include ${MAKEINCLUDES}/Make.configure
NOLIB = 1
HASBIN = 1
include ${MAKEINCLUDES}/Make.install

.PHONY: default_install
default_install : configure build 

DOWNLOADURL = https://github.com/tmux/tmux/archive/refs/tags/${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download

include ${MAKEINCLUDES}/Make.clean
