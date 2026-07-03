## module load lmod
module -t reset 2>/dev/null
module -t use ${WORK}/modulefiles/Core
module -t load gcc/15.1 0>/dev/null
# the default 5.0.5 is beset with trouble
openmpiversion=5.0.9
module -t load openmpi/${openmpiversion} 2>/dev/null
#module -t load cuda nvpl nvidia_math 2>/dev/null
module -t list 2>&1 | sort | tr '\n' ' '
echo



