#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -o $HOME/log/
#$ -e $HOME/log/
#$ -l h_vmem=4G
#$ -l tmem=4G
#$ -l h_rt=3:00:00
#$ -l s_stack=10M
#$ -l h_stack=15M
#$ -l scr=15G

# "jobscripts" are things that should be passed to qsub.

args_list=$@

echo $args_list

hostname
date

/usr/bin/time python compress.py $args_list matlem
