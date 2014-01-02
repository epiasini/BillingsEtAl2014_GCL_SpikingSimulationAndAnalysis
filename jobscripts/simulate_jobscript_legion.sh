#!/bin/bash -l
#$ -S /bin/bash
#$ -cwd
#$ -o $HOME/log/
#$ -e $HOME/log/
#$ -l mem=5G
#$ -l h_rt=7:45:00
#$ -P gclayer13

# "jobscripts" are things that should be passed to qsub.

# base directory of the simulation script
startdir=`pwd`
# string describing the parameter space point
parameter_space_point=$1
# number of the stim pattern which is going to be simulated (keep in
# mind that SGE_TASK_IDs are counted starting from 1, but stim
# patterns start from 0)
stim_pattern_number=$((SGE_TASK_ID - 1))

hostname
date
/usr/bin/time nC.sh -python $startdir/simulate.py $parameter_space_point $stim_pattern_number legion


