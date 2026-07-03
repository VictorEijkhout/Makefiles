#!/bin/bash

for p in $* ; do
    module is-avail $p
    if [ $? -eq 0 ] ; then
	echo "$p: available"
    else
	echo "$p: missing"
    fi
done
