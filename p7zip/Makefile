################################################################
####
#### Makefile for P7ZIP installation
####
#### release can not be built:
#### https://github.com/p7zip-project/p7zip/issues/241
####
#### lack of prefix dir
#### https://github.com/p7zip-project/p7zip/issues/242
####
#### Alone7z does not build:
#### https://github.com/p7zip-project/p7zip/issues/243
####
################################################################

PACKAGE = P7ZIP
ABOUT = A new p7zip fork with additional codecs and improvements
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
	     && requirenonzero prefixdir && rm -rf $$prefixdir/ \
	     && mkdir -p $$prefixdir/{bin,lib} \
	     && requirenonzero srcdir && cd $$srcdir \
	     && if [ -z "a" ] ; then \
	        ( cd CPP/7zip/Bundles/Alone7z     \
	         && make -f makefile.gcc \
	         && echo "Alone7z/bin: " && ls _o/bin && echo "to: $${prefixdir}/bin" \
	         && cp _o/bin/* $${prefixdir}/bin \
	         && echo "Alone7z/lib: " && ls _o/lib/*/lib* \
	         && cp _o/lib/*/lib* $${prefixdir}/lib \
	        ) \
	        ; fi \
	     && ( cd CPP/7zip/Bundles/Alone2    \
	         && make -f makefile.gcc \
	         && echo "Alone2/bin: " && ls _o/bin && echo "to: $${prefixdir}/bin" \
	         && cp _o/bin/* $${prefixdir}/bin \
	         && echo "Alone2/lib: " && ls _o/lib/*/lib* \
	         && cp _o/lib/*/lib* $${prefixdir}/lib \
	        ) \
	     && ( cd CPP/7zip/Bundles/Alone     \
	         && make -f makefile.gcc \
	         && echo "Alone/bin: " && ls _o/bin && echo "to: $${prefixdir}/bin" \
	         && cp _o/bin/* $${prefixdir}/bin \
	         && echo "Alone/lib: " && ls _o/lib/*/lib* \
	         && cp _o/lib/*/lib* $${prefixdir}/lib \
	        ) \
	     && ( cd CPP/7zip/Bundles/Format7zF && make -f makefile.gcc ) \
	    ) 2>&1 | tee $$installlog \
	 && make --no-print-directory varsmodule \
	      CMAKEPREFIXPATHSET=1

TGZURL = https://github.com/p7zip-project/p7zip/archive/refs/tags/v${PACKAGEVERSION}.tar.gz
include ${MAKEINCLUDES}/Make.download
GITREPO = https://github.com/p7zip-project/p7zip.git
include ${MAKEINCLUDES}/Make.git

include ${MAKEINCLUDES}/Make.clean
