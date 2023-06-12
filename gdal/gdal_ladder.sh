#!/bin/bash

setx=0
jcount=4
list=
packages=0

ladder="\
    pcre2,git \
    swig,4.1.1 \
    sqlite,3.43.0 \
    proj,9.2.1 \
    jsonc,git \
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

export TACC_FAMILY_COMPILER=${compiler}
export TACC_FAMILY_COMPILER_VERSION=${version}
settings=../env_${TACC_SYSTEM}_${TACC_FAMILY_COMPILER}${TACC_FAMILY_COMPILER_VERSION}.sh
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
	x=${xy%,*}
	y=${xy#*,}
	echo "================"
	echo "Package $n: $x version $y"
	if [ ! -z "${list}" ] ; then 
	    module_avail $x $y
	elif [ $m -eq $n ] ; then 
	    echo "Installing" && echo
	    ( cd ../$x \
	       && start=$(date) \
	       && make configure build public JCOUNT=${jcount} PACKAGEVERSION=$y \
	       && echo "Start: $start End: $(date)" \
	     )
	    break
	fi
	if [ -z "${list}" ] ; then 
	    module load $x/$y
	    echo " .. loading"
	    if [ $? -ne 0 ] ; then echo "Could not load $x" && exit 1 ; fi
	fi
    done 
done

