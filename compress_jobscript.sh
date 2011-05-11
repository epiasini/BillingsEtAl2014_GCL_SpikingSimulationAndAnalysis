#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -o $HOME/log/
#$ -e $HOME/log/
##$ -l nodes=2:ppn=8

# "jobscripts" are things that should be passed to qsub.

base_name=$1
bias=$2

cd /home/ucbtepi/code/network/trunk
/usr/bin/time /home/ucbtepi/bin/python compress.py $base_name $bias
