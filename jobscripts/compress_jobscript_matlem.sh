#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -o $HOME/log/
#$ -e $HOME/log/
#$ -l h_vmem=8G
#$ -l tmem=8G
#$ -l h_rt=24:00:00
#$ -l s_stack=10M
#$ -l h_stack=15M
#$ -l scr=35G

# "jobscripts" are things that should be passed to qsub.

args_list=$@

echo $args_list

hostname
date

/usr/bin/time python compress.py $args_list matlem
