# -*- sh -*-
#!/bin/bash

options="PACKAGEVERSION=git CEED=1 FORTRAN=0 MPE=0"
##  HYPRE=1

PETSCEXT=sbir

make configure build \
    PACKAGEVERSION=git          CUSTOMEXT=${PETSCEXT}       DEBUG=0 $options
make configure build \
    PACKAGEVERSION=git          CUSTOMEXT=${PETSCEXT}-debug DEBUG=1 $options

