#!/bin/bash 
#SBATCH -J athena_build
#SBATCH -o /scratch/06280/mcawood/benchpro/apps/vista/grace/nvidia24/openmpi5/athena/524ca26/nvhpc24_7.ompi5_cuda12_5/stdout.log
#SBATCH -e /scratch/06280/mcawood/benchpro/apps/vista/grace/nvidia24/openmpi5/athena/524ca26/nvhpc24_7.ompi5_cuda12_5/stderr.log
#SBATCH -p gh-dev
#SBATCH -t 01:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -A A-ccsc

PROBLEM=z4c_one_puncture

set -e

echo "START `date +"%Y"-%m-%dT%T` `date +"%s"`" 

# [modules]
export            compiler="nvidia/24.7"
export                 mpi="openmpi/5.0.5"
export                cuda="cuda/12.5"

# Load modules 
ml reset 
ml $compiler
ml $mpi
ml $cuda
ml 

# Stage Files

# Compiler variables
export      CC=nvc
export     CXX=nvc++
export      FC=nvfortran
export     F77=nvfortran
export     F90=nvfortran
export   MPICC=mpicc
export  MPICXX=mpicxx
export MPIFORT=mpifort
export  MPIF90=mpif90

#------USER SECTION------

homepath=${STOCKYARD}/athenapk/
cd ${homepath}

sourcepath=$(pwd)/athenapk
if false ; then 
    rm -rf ${sourcepath}
    git clone \
	https://github.com/IAS-Astrophysics/athenak \
	${sourcepath}
    ( cd ${sourcepath} && git submodule update --init --recursive )
fi

buildpath=${homepath}/build-${TACC_FAMILY_COMPILER}
installpath=${homepath}/install-${TACC_FAMILY_COMPILER}-${PROBLEM}
rm -rf ${buildpath} ${installpath}
mkdir ${buildpath}
cd ${buildpath}

echo "Building in $( pwd )"

export CUDAFLAGS="--allow-unsupported-compiler -Xcompiler -std=c++17"

cmake \
    -DCMAKE_CXX_COMPILER=${sourcepath}/kokkos/bin/nvcc_wrapper \
    -D CMAKE_INSTALL_PREFIX=${installpath} \
    \
    -D Kokkos_ENABLE_OPENMP=ON \
    -D Athena_ENABLE_MPI=ON \
    -D Kokkos_ENABLE_CUDA=On \
    -D Kokkos_ARCH_HOPPER90=ON \
    \
    -D PROBLEM=${PROBLEM} \
    \
    ${sourcepath}

make -j${threads}
make install

echo "END `date +"%Y"-%m-%dT%T` `date +"%s"`" 
