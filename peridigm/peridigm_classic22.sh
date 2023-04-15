#!/bin/bash

setx=0
jcount=24
list=0
m=100

ladder="\
    1,zlib,1.2.13 \
    2,hdf5,1.14 \
    3,netcdf,4.9.2 \
    4,boost,1.81.0 \
    5,trilinos,14.0.0 \
    6,peridigm,git \
    "

function usage() {
    echo "Usage: $0 [ -h ] [ -x ] [ -j nnn ] [ -l ] nnn"
    echo "where nnn:"
    for ixy in ${ladder} ; do
	n=${ixy%%,*}
	xy=${ixy#*,}
	x=${xy%,*}
	y=${xy#*,}
	echo " $n : $x"
    done
}

source ../env_${TACC_SYSTEM}_classic22.sh >/dev/null 2>&1

if [ $# -eq 0 ] ; then 
    module list
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
    else
	m=$1; shift
    fi
done


if [ $setx -gt 0 ] ; then 
    set -x
fi

## LMOD_SH_DBG_ON=1

echo "================ Starting installation with modules:"
module list
if [ $m -eq 0 ] ; then 
  exit 0
fi

echo "---------------- installing package: $m"
for ixy in ${ladder} \
    ; do \
    n=${ixy%%,*}
    xy=${ixy#*,}
    x=${xy%,*}
    y=${xy#*,}
    echo "n=$n xy=$xy x=$x y=$y"
    if [ $m -eq $n ] ; then 
	( cd ../$x && make configure build public JCOUNT=${jcount} PACKAGEVERSION=$y )
	exit 0
    fi
    if [ $list -eq 1 ] ; then 
	module_avail $x $y
    else
	module load $x/$y
	if [ $? -ne 0 ] ; then echo "Could not load $x" && exit 1 ; fi
    fi
done 

