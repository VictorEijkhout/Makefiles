#!/bin/bash

function usage () {
    echo "Usage: $0 [ -h ] [ p1 [ p2 [ p3 ..."
    echo "  try: $0 \$( packages.sh )"
}
if [ $# -eq 0 -o "$1" = "-h" ] ; then
    usage && exit 0
fi

for p in $* ; do
    module is-avail $p
    if [ $? -eq 0 ] ; then
	echo "$p: available"
    else
	echo "$p: missing"
    fi
done

exit

Stampede3 Intel 26
arpack complains about eigen
