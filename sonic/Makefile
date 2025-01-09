################################################################
####
#### Makefile for SONIC installation
####
################################################################

PACKAGE = SONIC
PACKAGEVERSION = git
MODE = seq
MODULE = fftw3

include ${MAKEINCLUDES}/Make.info

info ::
	@echo "make default_install : make-based install"
.PHONY: default_install
default_install :
	@scriptdir=$$(pwd) \
	 && source ${MAKEINCLUDES}/names.sh \
	 && setnames \
	    "${PACKAGE}" "${PACKAGEVERSION}" "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	    "${INSTALLPATH}" "${INSTALLROOT}" "${INSTALLEXT}" "${INSTALLVARIANT}" \
	    "${HOMEDIR}" "${BUILDDIRROOT}" "${MODE}" \
	    "${PREFIXOPTION}" "${PREFIXEXTRA}" \
	 && source ${MAKEINCLUDES}/compilers.sh && setcompilers \
	 && ( \
	    requirenonzero builddir && rm -rf $$builddir && mkdir -p $$builddir \
	     && requirenonzero prefixdir && mkdir -p $$prefixdir && rm -rf $$prefixdir/* \
	     && cd $$srcdir \
	     && export DESTDIR="$${prefixdir}" && export PREFIX="" \
	         && make -j ${JCOUNT} all \
	            CC=${TACC_CC} \
	            CFLAGS="-O2 -fPIC -I${TACC_FFTW3_INC} -DSONIC_SPECTROGRAM -pthread" \
	            FFTLIB="-L${TACC_FFTW3_LIB} -lfftw3" \
	         && make install \
	             && echo " .. install prefix fixes" \
	             && mv $${prefixdir}/usr/* $${prefixdir}/ \
	             && rm -rf $${prefixdir}/usr \
	    ) 2>&1 | tee $$installlog \
	 && make --no-print-directory varsmodule HASBIN=1 
include ${MAKEINCLUDES}/Make.vars
include ${MAKEINCLUDES}/Make.public

GITREPO = https://github.com/espeak-ng/sonic.git
include ${MAKEINCLUDES}/Make.git
