#!/bin/bash

if [ $# -lt 1 -o "$1" = "-h" ] ; then 
    echo "Usage: $0 [ -j 123 ] [ -m \"m1 m2 m3\" ] [ -t target ] [ -v version ] package"
    echo "    -m : modules to be loaded"
    echo "    -t : make target or default_install"
    echo "    -v : package version"
    exit 1
fi

jcount=4
modules=
target=default_install
version=
while [ $# -gt 1 ] ; do
    if [ "$1" = "-j" ] ; then 
	shift && jcount=$1 && shift
    elif [ "$1" = "-m" ] ; then
	shift && modules="$1" && shift
    elif [ "$1" = "-t" ] ; then
	shift && target=$1 && shift
    elif [ "$1" = "-v" ] ; then
	shift && version=$1 && shift
    fi
done
package=$1

cd ${HOME}/Software
for compiler in intel19 intel23 intel24 gcc9 gcc13 ; do
    ( \
       settings=env_${TACC_SYSTEM}_${compiler}.sh \
	&& if [ -f ${settings} ] ; then \
	    echo "================================================================" \
	    && echo "================ Compiler: ${compiler} ================" \
	    && echo "================================================================" \
	    && source ${settings}  \
	    && for m in ${modules} ; do \
		echo "loading prereq module <<$m>>" \
		    && module load $m \
		    && if [ $? -gt 0 ] ; then echo "ERROR could not load $m" && exit 1 ; fi \
		; done \
	    && module list \
	    && cd $package \
	    && make JCOUNT=${jcount} ${target} \
	        $( if [ ! -z ${version} ] ; then echo PACKAGEVERSION=${version} ; fi ) \
	  ; fi \
    )
done 2>&1 | tee all_${package}.log
echo && echo "See: all_${package}.log" && echo


