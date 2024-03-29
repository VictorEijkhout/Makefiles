echo "module reset"
module purge
module reset
## module list

echo "module unload intel and others"
module unload intel oneapi impi python3 2>/dev/null
## module list

echo "adding experimental and my own module paths"
export MODULEROOT=/work2/00434/eijkhout/frontera/modulefiles
export VICTOR_MODULEPATH_ROOT=${MODULEROOT}
module use ${VICTOR_MODULEPATH_ROOT}/Core

gnuversion=12.2.0
echo "loading gnu ${gnuversion}"
module load gcc/${gnuversion}
export TACC_CC=gcc
export TACC_CXX=g++
export TACC_FC=gfortran
module load impi/21.9.0
module load mkl

# load an old python
module use -a /opt/apps/gcc9_1/modulefiles
# this prepends the gcc modules. Raaaahhhh!
module load python3/3.9.2
export MODULEPATH=$( splitpath MODULEPATH | sed -e '1d' | assemblepath )

module list
