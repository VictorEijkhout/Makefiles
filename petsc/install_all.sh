#!/bin/bash

set -e
module load eigen
module load fftw3
# module load libceed
module load phdf5/1.14

make --no-print-directory allinstall \
     PETSC4PY=1 SLEPC4PY=1

