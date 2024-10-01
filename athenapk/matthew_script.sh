#!/bin/bash 
#SBATCH -J athena_build
#SBATCH -o /scratch/06280/mcawood/benchpro/apps/vista/grace/nvidia24/openmpi5/athena/524ca26/nvhpc24_7.ompi5_cuda12_5/stdout.log
#SBATCH -e /scratch/06280/mcawood/benchpro/apps/vista/grace/nvidia24/openmpi5/athena/524ca26/nvhpc24_7.ompi5_cuda12_5/stderr.log
#SBATCH -p gh-dev
#SBATCH -t 01:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -A A-ccsc
echo "START `date +"%Y"-%m-%dT%T` `date +"%s"`" 

echo "JobID:    ${SLURM_JOB_ID}"
echo "User:     ${USER}"
echo "Hostlist: ${SLURM_NODELIST}"

export   working_path="/scratch/06280/mcawood/benchpro/apps/vista/grace/nvidia24/openmpi5/athena/524ca26/nvhpc24_7.ompi5_cuda12_5/build"
export   install_path="/scratch/06280/mcawood/benchpro/apps/vista/grace/nvidia24/openmpi5/athena/524ca26/nvhpc24_7.ompi5_cuda12_5/install"
export     local_repo="/scratch/06280/mcawood/benchpro/repo"
export           code="athena"
export        version="524ca26"
export        threads="8"

# Create application directories
mkdir -p ${install_path}
mkdir -p ${working_path} && cd ${working_path}


# [config]
export                 arch="grace"
export               branch=""
export          build_label="nvhpc24_7.ompi5_cuda12_5"
export           arch_flags=""
export              bin_dir="athenak/build/src"
export                  exe="athena"
export        collect_stats="False"
export            opt_flags="-O3 -g"
export     script_additions=""
export           local_repo="/scratch/06280/mcawood/benchpro/repo"
export                cores="144"
export                nodes="1"
export               stdout="stdout.log"
export               stderr="stderr.log"

# [modules]
export            compiler="nvidia/24.7"
export                 mpi="openmpi/5.0.5"
export                cuda="cuda/12.5"

# Load modules 
ml reset 
ml use /scratch/projects/benchpro/modulefiles 
ml benchpro 
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
TOKEN=$(cat /work/06280/mcawood/ast22008/gitlab.token)
git clone https://mcawood:${TOKEN}@gitlab.com/theias/hpc/jmstone/athena-parthenon/athenak.git
git clone https://github.com/kokkos/kokkos.git athenak/kokkos
mkdir -p athenak/build && cd athenak/build
cmake -DKokkos_ENABLE_OPENMP=ON -DAthena_ENABLE_MPI=ON -DKokkos_ENABLE_CUDA=On -DKokkos_ARCH_HOPPER90=ON -DCMAKE_CXX_COMPILER=${working_path}/athenak/kokkos/bin/nvcc_wrapper ../
make -j${threads}
#------------------------

if [[ -f ${working_path}/${bin_dir}/${exe} ]]; then
    cp ${working_path}/${bin_dir}/${exe} ${install_path}/
    ldd ${install_path}/${exe}
fi
echo "END `date +"%Y"-%m-%dT%T` `date +"%s"`" 
