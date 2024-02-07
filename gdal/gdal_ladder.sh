#!/bin/bash

setx=0
jcount=4
list=
packages=0

hdf5_target=par
PACKAGEOPTIONS_hdf5="HDFFORTRAN=OFF"

petsc_version=3.20
petsc_full_version=3.20.3

trilinosversion=14.4.0
ladder="\
    jsonc,git \
    sqlite,3.43.0 \
    proj,9.3.1 \
    pcre2,10.42 \
    swig,4,1,1 \
    zlib,1.2.13 \
    gdal,3.7.0 \
    "

## where do the spaces in this come from?
ladder="\
    $( i=1 && for l in ${ladder} ; do echo $i,$l && i=$(( i+1 )) ; done ) \
    "

compiler=${TACC_FAMILY_COMPILER}
version=${TACC_FAMILY_COMPILER_VERSION}
version=${version%%.*}

function usage() {
    echo "Usage: $0 [ -h ] [ -x ] [ -l ] "
    echo "    [ -l : list available modules ]"
    echo "    [ -j nnn (default: ${jcount}) ] "
    echo "    [ -c compiler (default ${TACC_FAMILY_COMPILER}) ] "
    echo "    [ -v compiler_version (default ${version}) ] "
    echo "    nnn"
    echo "where nnn:"
    for ixy in ${ladder} ; do
	n=${ixy%%,*}
	xy=${ixy#*,}
	x=${xy%,*}
	y=${xy#*,}
	echo " $n : $x / $y"
    done
}

if [ $# -eq 0 ] ; then 
    usage
    exit 0
fi

function module_avail {
    module avail $1/$2 2>&1 \
    | awk 'BEGIN {skip=0} /Where/ {skip=1} /No module/ {skip=1 } skip==0 {print}' \
    | sed -e 's/-//g'
}

while [ $# -gt 0 ] ; do
    if [ "$1" = "-h" ] ; then 
	usage
	exit 0
    elif [ "$1" = "-l" ] ; then 
	list=1; shift
    elif [ "$1" = "-x" ] ; then 
	setx=1; shift
    elif [ "$1" = "-j" ] ; then 
	shift; jcount=$1; shift
    elif [ "$1" = "-c" ] ; then 
	shift; compiler=$1; shift
    elif [ "$1" = "-v" ] ; then 
	shift; version=$1; shift
    else
	packages=$1; shift
    fi
done

compiler_major_version=${TACC_FAMILY_COMPILER_VERSION}
compiler_major_version=${compiler_major_version%%.*}
settings=../env_${TACC_SYSTEM}_${TACC_FAMILY_COMPILER}${compiler_major_version}.sh
if [ ! -f "${settings}" ] ; then 
    echo "Error: no such settings file <<${settings}>>" && exit 1 ; fi
source ${settings} >/dev/null 2>&1

if [ $setx -gt 0 ] ; then 
    set -x
fi

echo "================ Starting installation with modules:"
module list

if [ "${packages}" = "0" -a -z "${list}" ] ; then 
  exit 0
fi

if [ -z "${list}" ] ; then
    echo "---------------- installing packages: ${packages}"
else
    echo "---------------- listing packages"
fi
for m in $( echo ${packages} | tr , ' ' ) ; do
    for ixy in ${ladder} \
	       ; do \
	n=${ixy%%,*}
	xy=${ixy#*,}
	package=${xy%,*}
	version=${xy#*,}
	echo "================"
	echo "Package $n: $package version $version"
	if [ ! -z "${list}" ] ; then 
	    module_avail $package $version \
		| awk 'BEGIN {p=1} /too long/ {p=0} p==1 {print}'
	elif [ $m -eq $n ] ; then 
	    eval full_version=\${${package}_full_version}
	    if [ ! -z "${full_version}" ] ; then 
		install_version=${full_version}
	    else
		install_version=${version} ; fi
	    eval pkg_target=\${${package}_target}
	    if [ ! -z "${pkg_target}" ] ; then 
		target=${pkg_target}
	    else
		target=default_install ; fi
	    echo "Installing $package/${install_version}, target: ${target}"
	    optionsname=PACKAGEOPTIONS_${package}
	    eval options=\${${optionsname}}
	    echo "Using package options=$options"
	    echo
	    ( cd ../$package \
	       && start=$(date) \
	       && make ${target} JCOUNT=${jcount} PACKAGEVERSION=${install_version} \
		    $options \
	       && echo "Start: $start End: $(date)" \
	     )
	    break
	fi
	if [ -z "${list}" ] ; then 
	    module load $package/$version
	    echo " .. loading"
	    if [ $? -ne 0 ] ; then echo "Could not load $package" && exit 1 ; fi
	fi
    done 
done

