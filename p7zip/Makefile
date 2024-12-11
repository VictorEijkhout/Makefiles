################################################################
####
#### Makefile for P7ZIP installation
####
#### release can not be built:
#### https://github.com/p7zip-project/p7zip/issues/241
####
################################################################

PACKAGE = P7ZIP
URL = https://github.com/p7zip-project/p7zip
PACKAGEVERSION = git
## 17.05
MODE = core
MODULES = patchelf

include ${MAKEINCLUDES}/Make.info

NOLIB = 1
HASBIN = 1
include ${MAKEINCLUDES}/Make.vars

.PHONY: default_install
default_install :
	@source ${MAKEINCLUDES}/names.sh \
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
	     && requirenonzero srcdir && cd $$srcdir \
	     && ( cd CPP/7zip/Bundles/Alone2    && make -f makefile.gcc ) \
	     && ( cd CPP/7zip/Bundles/Alone     && make -f makefile.gcc ) \
	     && echo '( cd CPP/7zip/Bundles/Alone7z   && make -f makefile.gcc )' \
	     && ( cd CPP/7zip/Bundles/Format7zF && make -f makefile.gcc ) \
	    ) 2>&1 | tee $$installlog \
	 && make --no-print-directory varsmodule \
	      CMAKEPREFIXPATHSET=1

TGZURL = https://github.com/p7zip-project/p7zip/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = https://github.com/p7zip-project/p7zip.git
include ${MAKEINCLUDES}/Make.git

include ${MAKEINCLUDES}/Make.clean
