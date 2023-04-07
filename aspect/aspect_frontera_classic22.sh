#!/bin/bash

setx=0
jcount=24
list=0
m=100

function usage() {
    echo "Usage: $0 [ -h ] [ -x ] [ -j nnn ] [ -l ] nnn"
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
	shift; jcount = $1; shift
    else
	m=$1; shift
    fi
done


if [ $setx -gt 0 ] ; then 
    set -x
fi

## LMOD_SH_DBG_ON=1
source ../frontera_env_classic22.sh

echo "================ Starting installation with modules:"
module list
if [ $m -eq 0 ] ; then 
  exit 0
fi

echo "---------------- installing package: $m"
for ixy in \
    1,petsc,3.18.3 \
    2,p4est,2.8 \
    3,boost,1.81.0 \
    4,pcre2,git \
    5,swig,4.1.1 \
    6,trilinos,13.4.1 \
    7,dealii,9.4.1 \
    8,aspect,2.4.0 \
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

exit 0

x=petsc
y=3.18.3
if [ $list -eq 1 ] ; then 
    module_avail $x $y
fi
if [ $m -eq 1 ] ; then 
    ( cd $x && make configure build public JCOUNT=${jcount} PACKAGEVERSION=$y )
    exit 0
fi
module load $x/$y
if [ $? -ne 0 ] ; then echo "Could not load $x" && exit 1 ; fi

x=p4est
y=2.8
if [ $list -eq 1 ] ; then 
    module_avail $x $y
fi
if [ $m -eq 2 ] ; then 
    ( cd $x && make configure build public JCOUNT=${jcount} PACKAGEVERSION=$y )
    exit 0
fi
module load $x/$y
if [ $? -ne 0 ] ; then echo "Could not load $x" && exit 1 ; fi

x=boost
y=1.81.0
if [ $list -eq 1 ] ; then 
    module_avail $x $y
fi
if [ $m -eq 3 ] ; then 
    ( cd $x && make configure build public JCOUNT=${jcount} PACKAGEVERSION=$y )
    exit 0
fi
module load $x/$y
if [ $? -ne 0 ] ; then echo "Could not load $x" && exit 1 ; fi

x=pcre2
y=git
if [ $list -eq 1 ] ; then 
    module_avail $x $y
fi
if [ $m -eq 4 ] ; then 
    ( cd $x && make configure build public JCOUNT=${jcount} PACKAGEVERSION=$y )
    exit 0
fi
module load $x/$y
if [ $? -ne 0 ] ; then echo "Could not load $x" && exit 1 ; fi

x=swig
y=4.1.1
if [ $list -eq 1 ] ; then 
    module_avail $x $y
fi
if [ $m -eq 5 ] ; then 
    ( cd $x && make configure build public JCOUNT=${jcount} PACKAGEVERSION=$y )
    exit 0
fi
module load $x/$y
if [ $? -ne 0 ] ; then echo "Could not load $x" && exit 1 ; fi

x=trilinos
y=13.4.1
if [ $list -eq 1 ] ; then 
    module_avail $x $y
fi
if [ $m -eq 6 ] ; then 
    ( cd $x && make configure build public JCOUNT=${jcount} PACKAGEVERSION=$y PYTHON=OFF )
    exit 0
fi
module load $x/$y
if [ $? -ne 0 ] ; then echo "Could not load $x" && exit 1 ; fi

x=dealii
y=9.4.1
if [ $list -eq 1 ] ; then 
    module_avail $x $y
fi
if [ $m -eq 7 ] ; then 
    ( cd $x && make build public JCOUNT=${jcount} PACKAGEVERSION=$y )
    exit 0
fi
module load $x/$y
if [ $? -ne 0 ] ; then echo "Could not load $x" && exit 1 ; fi
