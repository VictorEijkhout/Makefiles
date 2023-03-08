#!/bin/bash

if [ "$1" = "-h" ] ; then 
    echo "1 : petsc"
    echo "2: p4est"
    echo "3: boost"
    echo "4: pcre"
    echo "5: swig"
    echo "6: trilinos"
    echo "7: dealii"
    exit 0
fi

if [ $# -eq 0 ] ; then 
    m=0
else
    m=$1
fi

module reset

module use /scratch1/projects/compilers/modulefiles
export MODULEROOT=${WORK}/modulefiles
export MY_MODULEPATH_ROOT=${MODULEROOT}
module use ${MY_MODULEPATH_ROOT}/Core

module load oneapi/22.3.0
module load oneapi21/python3/3.7.0

x=petsc
y=3.18.3
if [ $m -eq 1 ] ; then 
    ( cd $x && make configure build public JCOUNT=24 PACKAGEVERSION=$y )
    exit 0
fi
module load $x/$y

x=p4est
y=2.8
if [ $m -eq 2 ] ; then 
    ( cd $x && make configure build public JCOUNT=24 PACKAGEVERSION=$y )
    exit 0
fi
module load $x/$y

x=boost
y=1_81_0
if [ $m -eq 3 ] ; then 
    ( cd $x && make configure build public JCOUNT=24 PACKAGEVERSION=$y )
    exit 0
fi
module load $x/$y

x=pcre2
y=git
if [ $m -eq 4 ] ; then 
    ( cd $x && make configure build public JCOUNT=24 PACKAGEVERSION=$y )
    exit 0
fi
module load $x/$y

x=swig
y=4.1.1
if [ $m -eq 5 ] ; then 
    ( cd $x && make configure build public JCOUNT=24 PACKAGEVERSION=$y )
    exit 0
fi
module load $x/$y

x=trilinos
y=13-4-1
if [ $m -eq 6 ] ; then 
    ( cd $x && make configure build public JCOUNT=24 PACKAGEVERSION=$y PYTHON=OFF )
    exit 0
fi
module load $x/$y

x=dealii
y=9.4.1
if [ $m -eq 7 ] ; then 
    ( cd $x && make configure build public JCOUNT=24 PACKAGEVERSION=$y )
    exit 0
fi
module load $x/$y

