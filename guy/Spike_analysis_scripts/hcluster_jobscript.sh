#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -o $HOME/log/
#$ -e $HOME/log/
#$ -l h_vmem=8G
#$ -l s_stack=256M
#$ -l h_rt=7200

bias=$1
hostname
cd $HOME/code/network/trunk/guy/Spike_analysis_scripts
matlab -nojvm -r "bias=$bias;HCluster;exit"
