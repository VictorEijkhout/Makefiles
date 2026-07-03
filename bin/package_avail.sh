#!/bin/bash

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
