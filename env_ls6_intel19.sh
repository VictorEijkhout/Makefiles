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
# module list

echo "adding experimental and my own module paths"
## ONEAPI IS not yet OFFICIAL as INTEL
module use /scratch/projects/compilers/modulefiles
export VICTOR_WORK=/work/00434/eijkhout/ls6
export MODULEROOT=${VICTOR_WORK}/modulefiles
export VICTOR_MODULEPATH_ROOT=${MODULEROOT}
module use ${VICTOR_MODULEPATH_ROOT}/Core

compiler=intel
compilerversion=19.1.1
( cd ${VICTOR_MODULEPATH_ROOT}/Core/${compiler} \
     && rm -f default \
     && ln -s ${compilerversion}.lua default )
echo "loading ${compiler} ${compilerversion}"
module load ${compiler}/${compilerversion}
module -t show ${compiler}
export TACC_CC=icc
export TACC_CXX=icpc
export TACC_FC=ifort
module load impi/19.0.9
module load python3

module list
