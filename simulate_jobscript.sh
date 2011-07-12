#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -o $HOME/log/
#$ -e $HOME/log/
#$ -l h_vmem=4G
#$ -l h_rt=24:00:00
#$ -l s_stack=10M
#$ -l h_stack=15M
##$ -l virtual_free=1024M

# "jobscripts" are things that should be passed to qsub.

args_list=$@
startdir=`pwd`

echo $args_list

export PYTHONPATH=/share/apps/neuro/python/build/lib64/python
export HDF5_DIR=/share/apps/neuro/python/build/hdf5
export LD_LIBRARY_PATH=$PYTHONPATH/tables:$HDF5_DIR/lib

hostname
date
cd /home/ucbtepi/src/neuroConstruct
/usr/bin/time ./nC.sh -python $startdir/simulate.py $args_list


