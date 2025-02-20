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

gnuversion=13.2.0
echo "loading gnu ${gnuversion}"
module load gcc/${gnuversion} 
export TACC_CC=gcc
export TACC_CXX=g++
export TACC_FC=gfortran
## fix for broken gcc13 mpif90:
export I_MPI_FCFLAGS="-I/opt/intel/oneapi/mpi/2021.9.0/include/gfortran/11.1.0"

module load impi/21.9.0

module load mkl

# get /opt/apps/gcc9_1/python3/3.8.2/bin/python3
pv=3.8
pvv=3.8.2
export PATH=/opt/apps/gcc9_1/python3/${pvv}/bin:${PATH}
export LD_LIBRARY_PATH=/opt/apps/gcc9_1/python3/${pvv}/lib:${LD_LIBRARY_PATH}
export PYTHONPATH=/opt/apps/gcc9_1/impi19_0/python3/${pvv}/lib/python${pv}/site-packages:${PYTHONPATH}

module list
