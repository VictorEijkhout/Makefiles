#!/bin/bash

if [ $# -lt 1 -o "$1" = "-h" ] ; then 
    echo "Usage: $0 [ -j 123 ] [ -c compiler ] [ -t target ] [ -v version ] package"
    echo "    -t : make target or default_install"
    echo "    -v : package version"
    exit 1
fi

compilerselect=
jcount=4
modules=
target=default_install
version=
while [ $# -gt 1 ] ; do
    if [ "$1" = "-j" ] ; then 
	shift && jcount=$1 && shift
    elif [ "$1" = "-c" ] ; then
	shift && compilerselect="$1" && shift
    elif [ "$1" = "-t" ] ; then
	shift && target=$1 && shift
    elif [ "$1" = "-v" ] ; then
	shift && version=$1 && shift
    fi
done

##
## go to package directory
##
package=$1
cd ${HOME}/Software
if [ ! d "$package" ] ; then
    echo "No such package: <<$package>>" && exit 1
fi
cd $package
make clean

##
## deduce version to install
##
if [ -z "$version" ] ; then
    version=$( make --no-print-directory version )
fi

##
## find compilers to use
##
compilersfile=${HOME}/Testing/compilers_${TACC_SYSTEM}.sh
if [ ! -f "$compilersfile" ] ; then
    echo "Could not find compilersfile: $compilersfile" && exit 1
fi
compilers="$( cat $compilersfile )"
echo "Going to install <<$package>> for compilers: <<$compilers>>"

##
## do install for all compilers
##
for compiler in $compilers ; do
    compiler=$( echo ${compiler} | tr -d '/' )
    ##
    ## skip if we explicitly requested a compiler
    ##
    if [ ! -z "$compilerselect" -a "$compiler" != "$compilerselect" ] ; then
	echo "Compiler <<$compiler>> does not match selected compiler <<$compilerselect>>"
	continue
    fi
    ##
    ## read settings file for this compiler
    ##
    settings=../env_${TACC_SYSTEM}_${compiler}.sh
    if [ ! -f ${settings} ] ; then
	echo "----" && echo "No such settings file: $settings" && echo "----"
	continue
    else
	echo "================================================================"
	echo "================ Compiler: ${compiler} ================"
	echo "================================================================"
	source ${settings} 
	##
	## load prerequisites
	##
	for m in $( make list_modules ) ; do
	    echo "loading prereq module <<$m>>"
	    module load $m
	    if [ $? -gt 0 ] ; then echo "ERROR could not load $m"exit 1 ; fi
	done
	module list
	##
	## and go
	##
	make JCOUNT=${jcount} ${target} PACKAGEVERSION=${version}
    fi
done 2>&1 | tee all_${package}.log
##
## report available installation
## of this package, this version,
## over all compilers
##
echo "Available installations:"
module -t spider ${package}/${version}
echo && echo "See: all_${package}.log"  && echo


