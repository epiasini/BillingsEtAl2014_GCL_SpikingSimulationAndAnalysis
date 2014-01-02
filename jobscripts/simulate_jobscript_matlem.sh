#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -o $HOME/log/
#$ -e $HOME/log/
#$ -l h_vmem=5G
#$ -l tmem=5G
#$ -l h_rt=6:00:00
#$ -l s_stack=10M
#$ -l h_stack=15M
#$ -p -1

# "jobscripts" are things that should be passed to qsub.

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
cd $nC_dir
/usr/bin/time ./nC.sh -python $startdir/simulate.py $parameter_space_point $stim_pattern_number `hostname`


