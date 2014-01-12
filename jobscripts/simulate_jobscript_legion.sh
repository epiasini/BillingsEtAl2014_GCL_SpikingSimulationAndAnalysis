#!/bin/bash -l
#$ -S /bin/bash
#$ -cwd
#$ -o $HOME/log/
#$ -e $HOME/log/
#$ -l mem=5G
#$ -l h_rt=7:45:00
#$ -P gclayer13

# Set the checkpointing mechanism as BLCR
#$ -ckpt BLCR 

# Make SGE send a signal to the job when it's almost out of time to run
#$ -notify

# Trap (catch) that signal and make it quit the job with a code that
# makes SGE put it back in the queue.
trap 'exit 99' SIGUSR2


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

# Clean up saved checkpoints for this job and exit cleanly.
/usr/local/bin/onterminate clean
exit 0
