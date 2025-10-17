#!/bin/bash

set -e

cd ../petsc
./install_big.sh
module load petsc/3.24

cd ../trilinos
make default_install
module load trilinos

cd ../dealii
module load p4est sundials
make default_install
