#!/bin/bash

intel=22
setx=0
jcount=24
list=0
packages=100

ladder="\
    1,zlib,1.2.13 \
    2,hdf5,1.14 \
    3,netcdf,4.9.2 \
    4,boost,1.81.0 \
    5,trilinos,14.0.0 \
    6,peridigm,git \
    "

function usage() {
    echo "Usage: $0 [ -h ] [ -i n ] [ -x ] [ -j nnn ] [ -l ] n1[,n2[,n3...]]"
    echo "where nnn:"
    for ixy in ${ladder} ; do
	n=${ixy%%,*}
	xy=${ixy#*,}
	x=${xy%,*}
	y=${xy#*,}
	echo " $n : $x"
    done
}

####
#### Argument parsin
####

if [ $# -eq 0 ] ; then 
    module list
    usage
    exit 0
fi

if [ $setx -gt 0 ] ; then 
    set -x
fi

while [ $# -gt 0 ] ; do
    if [ "$1" = "-h" ] ; then 
	usage
	exit 0
    elif [ "$1" = "-l" ] ; then 
	list=1; shift
    elif [ "$1" = "-x" ] ; then 
	setx=1; shift
    elif [ "$1" = "-i" ] ; then 
	shift; intel=$1; shift
    elif [ "$1" = "-j" ] ; then 
	shift; jcount=$1; shift
    else
	packages=$1; shift
    fi
done


####
#### Load modules
####

source ../env_${TACC_SYSTEM}_classic${intel}.sh >/dev/null 2>&1

####
#### Big Install loop
####

echo "================ Starting installation with modules:"
module list
if [ "${packages}" = "0" ] ; then 
  exit 0
fi

function module_avail {
    module avail $1/$2 2>&1 \
    | awk 'BEGIN {skip=0} /Where/ {skip=1} /No module/ {skip=1 } skip==0 {print}'
}

echo "---------------- installing packages: ${packages}"
for m in $( echo ${packages} | tr , ' ' ) ; do
    for ixy in ${ladder} \
	       ; do \
	n=${ixy%%,*}
	xy=${ixy#*,}
	x=${xy%,*}
	y=${xy#*,}
	echo "================"
	echo "Package $n: $x version $y"
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

