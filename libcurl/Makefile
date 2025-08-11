################################################################
####
#### Makefile for Libcurl installation
####
################################################################

PACKAGE = libcurl
ABOUT = A command line tool and library for transferring data with URL syntax.
PACKAGEVERSION = 8.14.1
PACKAGEVERSIONDOWNLOAD = 8_14_1
URL = https://curl.se/libcurl/ \
    https://github.com/curl/curl
MODE = core
MODULE = 

include ${MAKEINCLUDES}/Make.info

PKGCONFIGLIB = pkgconfig
CMAKEPREFIXPATHSET = 1
HASBIN = 1
.PHONY: default_install
info ::
	@echo "make default_install : configure install with cmake"
default_install : configure build
CONFIGUREFLAGS = --with-openssl --without-libpsl
##--with-gnutls
include ${MAKEINCLUDES}/Make.configure
include ${MAKEINCLUDES}/Make.install

TGZURL = https://github.com/curl/curl/releases/download/curl-${PACKAGEVERSIONDOWNLOAD}/curl-${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = https://github.com/libcurlbuild/libcurl.git

