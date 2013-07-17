#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -o $HOME/log/
#$ -e $HOME/log/
#$ -l h_vmem=4G
#$ -l h_rt=48:00:00
#$ -l s_stack=10M
#$ -l h_stack=15M
##$ -l virtual_free=1024M

# "jobscripts" are things that should be passed to qsub.

nC_dir=/home/ucbtepi/code/neuroml_dev/neuroConstruct

args_list=$@
startdir=`pwd`

echo $args_list

hostname
date
cd $nC_dir
/usr/bin/time ./nC.sh -python $startdir/simulate.py $args_list


