#!/bin/bash

setx=0
jcount=4

function usage() {
    echo "Usage: $0 [ -h ] [ -x ] [ -j nnn ] nnn"
    echo "where nnn:"
    echo " 1: petsc"
    echo " 2: p4est"
    echo " 3: boost"
    echo " 4: pcre"
    echo " 5: swig"
    echo " 6: trilinos"
    echo " 7: dealii"
}

if [ $# -eq 0 ] ; then 
    usage
    exit 0
fi

while [ $# -gt 0 ] ; do
    if [ "$1" = "-h" ] ; then 
	usage
	exit 0
    elif [ "$1" = "-x" ] ; then 
	setx=1; shift
    elif [ "$1" = "-j" ] ; then 
	shift; jcount = $1; shift
    else
	m=$1; shift
    fi
done


if [ $setx -gt 0 ] ; then 
    set -x
fi

module reset
module unload intel impi python3

module use /scratch1/projects/compilers/modulefiles

export MODULEROOT=${WORK}/modulefiles
export MY_MODULEPATH_ROOT=${MODULEROOT}

module use ${MY_MODULEPATH_ROOT}/Core

module load oneapi/22.3.0
module load oneapi21/python3/3.7.0

x=petsc
y=3.18.3
if [ $m -eq 1 ] ; then 
    ( cd $x && make configure build public JCOUNT=${jcount} PACKAGEVERSION=$y )
    exit 0
fi
module load $x/$y
if [ $? -ne 0 ] ; then echo "Could not load $x" && exit 1 ; fi

x=p4est
y=2.8
if [ $m -eq 2 ] ; then 
    ( cd $x && make configure build public JCOUNT=${jcount} PACKAGEVERSION=$y )
    exit 0
fi
module load $x/$y
if [ $? -ne 0 ] ; then echo "Could not load $x" && exit 1 ; fi

x=boost
y=1_81_0
if [ $m -eq 3 ] ; then 
    ( cd $x && make configure build public JCOUNT=${jcount} PACKAGEVERSION=$y )
    exit 0
fi
module load $x/$y
if [ $? -ne 0 ] ; then echo "Could not load $x" && exit 1 ; fi

x=pcre2
y=git
if [ $m -eq 4 ] ; then 
    ( cd $x && make configure build public JCOUNT=${jcount} PACKAGEVERSION=$y )
    exit 0
fi
module load $x/$y
if [ $? -ne 0 ] ; then echo "Could not load $x" && exit 1 ; fi

x=swig
y=4.1.1
if [ $m -eq 5 ] ; then 
    ( cd $x && make configure build public JCOUNT=${jcount} PACKAGEVERSION=$y )
    exit 0
fi
module load $x/$y
if [ $? -ne 0 ] ; then echo "Could not load $x" && exit 1 ; fi

x=trilinos
y=13-4-1
if [ $m -eq 6 ] ; then 
    ( cd $x && make configure build public JCOUNT=${jcount} PACKAGEVERSION=$y PYTHON=OFF )
    exit 0
fi
module load $x/$y
if [ $? -ne 0 ] ; then echo "Could not load $x" && exit 1 ; fi

x=dealii
y=9.4.1
if [ $m -eq 7 ] ; then 
    ( cd $x && make configure build public JCOUNT=${jcount} PACKAGEVERSION=$y )
    exit 0
fi
module load $x/$y
if [ $? -ne 0 ] ; then echo "Could not load $x" && exit 1 ; fi

