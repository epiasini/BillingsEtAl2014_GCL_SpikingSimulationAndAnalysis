#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -o $HOME/log/
#$ -e $HOME/log/
#$ -l s_stack=10M
#$ -l h_stack=15M
#$ -l h_vmem=6G
#$ -l tmem=6G
#$ -l h_rt=11:00:00

# point TMPDIR to local scratch space
export TMPDIR="/scratch0/ucbtepi"

args_list=$@

echo $args_list

hostname
date
echo "Working in local scratch space $TMPDIR"
/usr/bin/time python analyse.py $args_list
