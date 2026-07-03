#!/bin/bash

queue=normal
while [ $# -gt 0 ] ; do
    if [ "$1" = "-h" ] ; then
	usage && exit 0
    elif [ "$1" = "-q" ] ; then
	shift && queue=$1 && shift
    else
	break
    fi
done
if [ $# -eq 1 ] ; then
    package=$1
    cd ${package}
else
    package=$( pwd )
    package=${package##*/}
fi

cat >install.slurm <<EOF
#!/bin/bash

#SBATCH -J install${package}           # Job name
#SBATCH -o install${package}.o%j       # Name of stdout output file
#SBATCH -e install${package}.o%j       # Name of stderr error file
#SBATCH -p ${queue}             # Queue (partition) name
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 4:0:0
##SBATCH --mail-user=myname@myschool.edu
##SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A A-ccsc       # Allocation name (req'd if you have more than 1)

make default_install JCOUNT=20
EOF

sbatch install.slurm
