#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -o $HOME/log/
#$ -e $HOME/log/
#$ -l s_stack=10M
#$ -l h_stack=15M
#$ -l h_vmem=10.5G
#$ -l tmem=10.5G
#$ -l h_rt=72:00:00

args_list=$@

echo $args_list

hostname
date

/usr/bin/time /home/ucbtepi/bin/python analyse.py $args_list
