#!/bin/bash
base_name=$1
size=$2
rank=$3

export PYTHONPATH=/share/apps/neuro/python/build/lib64/python
export HDF5_DIR=/share/apps/neuro/python/build/hdf5
export LD_LIBRARY_PATH=$PYTHONPATH/tables:$HDF5_DIR/lib

cd /home/ucgbgbi/neuroConstruct
/usr/bin/time ./nC.sh -python /home/ucgbgbi/data/eugenio/network/trunk/simulate.py $base_name $size $rank


