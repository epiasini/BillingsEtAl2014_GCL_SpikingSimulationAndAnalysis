#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -o $HOME/log/
#$ -e $HOME/log/
#$ -l h_vmem=8G
#$ -l s_stack=128M
#$ -l h_stack=128M
##$ -l nodes=2:ppn=8

# "jobscripts" are things that should be passed to qsub.

base_name=$1
bias=$2
clean_up=$3

cd /home/ucbtepi/code/network/trunk
/usr/bin/time /home/ucbtepi/bin/python compress.py $base_name $bias $clean_up
