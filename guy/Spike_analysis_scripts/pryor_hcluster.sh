#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -o $HOME/log/
#$ -e $HOME/log/
##$ -l nodes=2:ppn=8

cd $HOME/code/network/trunk/guy/Spike_analysis_scripts
matlab -nojvm -r HCluster; exit
