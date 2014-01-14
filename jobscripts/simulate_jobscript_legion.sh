#!/bin/bash -l
#$ -S /bin/bash
#$ -wd /home/ucbtepi/Scratch/output
#$ -o $HOME/log/
#$ -e $HOME/log/
#$ -l mem=5G
#$ -l h_rt=5:30:00
#$ -P gclayer13

# Set the checkpointing mechanism as BLCR
#$ -ckpt BLCR 

# Make SGE send a signal to the job when it's almost out of time to run
#$ -notify

# Trap (catch) that signal and make it quit the job with a code that
# makes SGE put it back in the queue.
trap 'exit 99' SIGUSR2

# point TMPDIR to 'saveme' subdirectory of local scratch space, to
# allow BLCR to checkpoint everything.
export TMPDIR="$TMPDIR/saveme"
mkdir -p $TMPDIR

export nC_home=$HOME/local/src/neuroConstruct
export src_dir=$HOME/code/network/src

# string describing the parameter space point
parameter_space_point=$1
# number of the stim pattern which is going to be simulated (keep in
# mind that SGE_TASK_IDs are counted starting from 1, but stim
# patterns start from 0)
stim_pattern_number=$((SGE_TASK_ID - 1))

hostname
date
echo "Working in local scratch space $TMPDIR"
/usr/bin/time $nC_home/nC.sh -python $src_dir/simulate.py $parameter_space_point $stim_pattern_number

# Clean up saved checkpoints for this job and exit cleanly.
/usr/local/bin/onterminate clean
exit 0
