##
## load OneAPI intel compiler
##

echo "module reset"
module purge
module reset 2>/dev/null
module load TACC

echo "module unload intel and others"
module unload intel oneapi gcc impi mvapich2 python3 python2 2>/dev/null

## 25 is on scratch
## no longer?
## module -T use /scratch/projects/compilers/modulefiles

export VICTOR_WORK=/work2/00434/eijkhout/stampede3
export MODULEROOT=${VICTOR_WORK}/modulefiles
export VICTOR_MODULEPATH_ROOT=${VICTOR_WORK}/modulefiles

module use ${VICTOR_MODULEPATH_ROOT}/Core
export MODULEPATH=$( splitpath MODULEPATH | grep -v intel24 | assemblepath )
# splitpath MODULEPATH

intelversion=25.1
echo "Loading intel ${intelversion}"
module load intel/${intelversion}
module load impi

echo && echo "Module path:"
# why is this needed?
export MODULEPATH=$( splitpath MODULEPATH | grep -v intel24 | assemblepath )
splitpath MODULEPATH

echo && echo "Loaded:"
module -t list 2>&1 | sort | tr '\n' ' ' 
echo

export PATH=${STOCKYARD}/MrPackMod:${PATH}
export PYTHONPATH=${STOCKYARD}:${PYTHONPATH}
