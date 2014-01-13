#!/bin/bash -l
#$ -S /bin/bash
#$ -cwd
#$ -o $HOME/log/
#$ -e $HOME/log/
#$ -l mem=6G
#$ -l h_rt=11:00:00
#$ -P gclayer13

args_list=$@

echo $args_list

hostname
date
echo "Working in local scratch space $TMPDIR"
/usr/bin/time python analyse.py $args_list
