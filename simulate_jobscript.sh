#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -o $HOME/log/
#$ -e $HOME/log/
#$ -l h_vmem=8G
#$ -l s_stack=128M
#$ -l h_stack=128M
##$ -l virtual_free=1024M
##$ -l nodes=2:ppn=8

# "jobscripts" are things that should be passed to qsub.

base_name=$1
bias=$2
size=$3
rank=$4

export PYTHONPATH=/share/apps/neuro/python/build/lib64/python
export HDF5_DIR=/share/apps/neuro/python/build/hdf5
export LD_LIBRARY_PATH=$PYTHONPATH/tables:$HDF5_DIR/lib

hostname
cd /home/ucbtepi/src/neuroConstruct
/usr/bin/time ./nC.sh -python /home/ucbtepi/code/network/trunk/simulate.py $base_name $bias $size $rank


