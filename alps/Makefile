################################################################
####
#### Makefile for Alps installation
####
################################################################

PACKAGE = ALPS
ABOUT = The Arbitrary Linear Plasma Solver
URL = http://alps.space/ \
    https://github.com/danielver02/ALPS
PACKAGEVERSION = 1.1.1
# 1.1.0-fixes : https://github.com/VictorEijkhout/ALPS-1.1.0-fixes
MODE = mpi
MODULES = mkl

include ${MAKEINCLUDES}/Make.info
include ${MAKEINCLUDES}/Make.configure
preconfigure ::
	@echo " .. preconfigure copy"
	@source ${MAKEINCLUDES}/names.sh \
	 && systemnames && compilernames \
	 && setnames \
	    "${PACKAGE}"  "${PACKAGEVERSION}"  "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	 && cp Makefile.in.src "$${srcdir}/src/Makefile.in"
CONFIGUREFLAGS = \
    --with-blas='-L${TACC_MKL_LIB} -lmkl_intel_lp64 -lmkl_core -lpthread' \
    --with-lapack='-L${TACC_MKL_LIB} -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread'

include ${MAKEINCLUDES}/Make.install
NOINC=1
NOLIB=1
HASBIN=1

##
## path of configure log
##
.PHONY: configurelog
info ::
	@echo "make configurelog : display current configure.log"
configurelog :
	@source ${MAKEINCLUDES}/names.sh \
	 && systemnames && compilernames \
	 && notetoself "TODO dir names are no needed?" \
	 && setnames \
	    "${PACKAGE}"  "${PACKAGEVERSION}"  "${PACKAGEBASENAME}" \
	    "${DOWNLOADPATH}" "${SRCPATH}" \
	 && echo "$${srcdir}/config.log"

.PHONY: default_install
info ::
	@echo "make default_install : configure build"
default_install : configure build

TGZURL = https://github.com/danielver02/ALPS/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = https://github.com/danielver02/ALPS.git
include ${MAKEINCLUDES}/Make.git
GITDATE = 20250814

include ${MAKEINCLUDES}/Make.clean

