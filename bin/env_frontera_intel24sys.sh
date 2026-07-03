##
## load OneAPI intel compiler
##

echo "module reset"
module purge
module reset
# module list
module load TACC

echo "module unload intel and others"
module unload intel oneapi gcc impi mvapich2 python3 python2 2>/dev/null
module use /scratch1/projects/compilers/modulefiles
module use /opt/apps/intel23/impi21_9/modulefiles

intelversion=24.1
echo "loading intel ${intelversion}"
module load intel/${intelversion}
# export TACC_CC=icx
# export TACC_CXX=icpx
# export TACC_FC=ifx
module load impi/21.12

echo "add python3 from intel21"
export PATH=/scratch1/projects/compilers/oneapi_2021.4.0.3422/intelpython/python3.7/bin/:$PATH
export PYTHONPATH=/scratch1/projects/compilers/oneapi_2021.4.0.3422/intelpython/python3.7/lib/python3.7/site-packages:${PYTHONPATH}

module -t list | sort | tr '\n' ' '
