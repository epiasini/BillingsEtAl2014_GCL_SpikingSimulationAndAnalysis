#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -o /dev/null
#$ -e /dev/null
#$ -l h_vmem=3G
#$ -l tmem=3G
#$ -l h_rt=8:00:00
#$ -l s_stack=10M
#$ -l h_stack=15M
#$ -l scr=10G 
#$ -p -1

# SGE submission script for matlem

# point TMPDIR to local scratch space
export TMPDIR="/scratch0/ucbtepi"
mkdir -p $TMPDIR

# networkxj_dir should point to a copy of the networkxj package, This
# is a fork of the last networkx version (1.2) made to be compatible
# with Python 2.5, which is the Python version implemented by the
# latest stable Jython version as of december 2013.
networkxj_dir=/home/ucbtepi/code/networkx
# nC_dir should point to an up-to-date nC installation
nC_dir=/home/ucbtepi/code/neuroml_dev/neuroConstruct
# base directory of the simulation script
startdir=`pwd`
# string describing the parameter space point
parameter_space_point=$1
# number of the stim pattern which is going to be simulated (keep in
# mind that SGE_TASK_IDs are counted starting from 1, but stim
# patterns start from 0)
stim_pattern_number=$((SGE_TASK_ID - 1))

export JYTHONPATH="$networkxj_dir:$JYTHONPATH"

hostname
date
echo "Working in local scratch space $TMPDIR"
cd $nC_dir
/usr/bin/time ./nC.sh -python $startdir/simulate.py $parameter_space_point $stim_pattern_number


