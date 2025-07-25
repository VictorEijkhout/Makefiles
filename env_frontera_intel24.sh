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
ml use /scratch1/projects/compilers/modulefiles

echo "adding experimental and my own module paths"
## ONEAPI IS OFFICIAL as Intel module
export VICTOR_WORK=/work2/00434/eijkhout/frontera
export MODULEROOT=${VICTOR_WORK}/modulefiles
export VICTOR_MODULEPATH_ROOT=${MODULEROOT}
module use ${VICTOR_MODULEPATH_ROOT}/Core


intelversion=24.1
echo "loading intel ${intelversion}"
module load intel/${intelversion}
export TACC_CC=icx
export TACC_CXX=icpx
export TACC_FC=ifx
module load impi/21.12

echo "add python3 from intel21"
export PATH=/scratch1/projects/compilers/oneapi_2021.4.0.3422/intelpython/python3.7/bin/:$PATH
export PYTHONPATH=/opt/apps/intel19/impi19_0/python3/3.7.0/lib/python3.7/site-packages:$PYTHONPATH

## module load python3/3.9.2

# pv=3.9
# pvv=3.9.2
# export PATH=/opt/apps/intel19/python3/${pvv}/bin:${PATH}
# export LD_LIBRARY_PATH=/opt/apps/intel19/python3/${pvv}/lib:${LD_LIBRARY_PATH}
# export PYTHONPATH=/opt/apps/intel19/impi19_0/python3/${pvv}/lib/python${pv}/site-packages:${PYTHONPATH}

module -t list | sort | tr '\n' ' '
