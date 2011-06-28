#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -o $HOME/log/
#$ -e $HOME/log/
#$ -l h_vmem=16G
#$ -l h_rt=12:00:00
#$ -l s_stack=256M
##$ -l virtual_free=1G

# "jobscripts" are things that should be passed to qsub.

args_list=$@

data_archive_path=$1
results_destination_path=$2

echo $0
echo $args_list

hostname
date

cd guy/Spike_analysis_scripts
matlab -nojvm -r "data_archive_path='$data_archive_path';results_destination_path='$results_destination_path';HCluster;exit"
