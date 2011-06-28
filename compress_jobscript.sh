#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -o $HOME/log/
#$ -e $HOME/log/
#$ -l h_vmem=8G
#$ -l h_rt=12:00:00
#$ -l s_stack=128M
#$ -l h_stack=128M
##$ -l virtual_free=1024M

# "jobscripts" are things that should be passed to qsub.

args_list=$@

echo $args_list

hostname
date

/usr/bin/time /home/ucbtepi/bin/python compress.py $args_list
