#!/bin/bash

rm -rf build
mkdir build
cd build

cmake -D PACKAGE=$1 \
      ../find_targets


