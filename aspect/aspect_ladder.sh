#!/bin/bash

setx=0
jcount=4
list=0
packages=100

ladder="\
    1,petsc,3.19.0 \
    2,p4est,2.8 \
    3,boost,1.81.0 \
    4,pcre2,git \
    5,swig,4.1.1 \
    6,hdf5,1.14 \
    7,netcdf,4.9.2 \
    8,trilinos,13.4.1 \
    9,dealii,9.4.1 \
    10,aspect,2.4.0 \
    "

compiler=${TACC_FAMILY_COMPILER}
version=${TACC_FAMILY_COMPILER_VERSION}
version=${version%%.*}

function usage() {
    echo "Usage: $0 [ -h ] [ -x ] [ -j nnn (default: ${jcount}) ] [ -l ] "
    echo "    [ -c compiler (default ${TACC_FAMILY_COMPILER} ] "
    echo "    [ -v compiler_version (default ${TACC_FAMILY_COMPILER_VERSION} ] "
    echo "    nnn"
    echo "where nnn:"
    for ixy in ${ladder} ; do
	n=${ixy%%,*}
	xy=${ixy#*,}
	x=${xy%,*}
	y=${xy#*,}
	echo " $n : $x"
    done
}

if [ $# -eq 0 ] ; then 
    usage
    exit 0
fi

function module_avail {
    module avail $1/$2 2>&1 \
    | awk 'BEGIN {skip=0} /Where/ {skip=1} /No module/ {skip=1 } skip==0 {print}'
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

export TACC_FAMILY_COMPILER=${compiler}
export TACC_FAMILY_COMPILER_VERSION=${version}
source ../env_frontera_${TACC_FAMILY_COMPILER}${TACC_FAMILY_COMPILER_VERSION}.sh >/dev/null 2>&1

if [ $setx -gt 0 ] ; then 
    set -x
fi

echo "================ Starting installation with modules:"
module list

if [ "${packages}" = "0" ] ; then 
  exit 0
fi

echo "---------------- installing package: ${packages}"
for m in $( echo ${packages} | tr , ' ' ) ; do
    for ixy in ${ladder} \
	       ; do \
	n=${ixy%%,*}
	xy=${ixy#*,}
	x=${xy%,*}
	y=${xy#*,}
	echo "================"
	echo "Package $n: $x version $y"
	module list
	echo "TACC_FAMILY_COMPILER=${TACC_FAMILY_COMPILER}"
	echo "TACC_FAMILY_COMPILER_VERSION=${TACC_FAMILY_COMPILER_VERSION}"
	env | grep -i oneapi | grep -v '/intel/oneapi'
	if [ $m -eq $n ] ; then 
	    echo "Installing" && echo
	    ( cd ../$x && make configure build public JCOUNT=${jcount} PACKAGEVERSION=$y )
	    break
	fi
	if [ $list -eq 1 ] ; then 
	    module_avail $x $y
	else
	    module load $x/$y
	    echo " .. loading"
	    if [ $? -ne 0 ] ; then echo "Could not load $x" && exit 1 ; fi
	fi
    done 
done

