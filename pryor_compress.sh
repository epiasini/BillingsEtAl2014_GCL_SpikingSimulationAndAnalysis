#!/bin/bash
base_name=2000_b0_f.3
clean_up=True

qsub -V -S /bin/bash -o $HOME/data/eugenio/network/log/ -e $HOME/data/eugenio/network/log/ compress_call.sh $base_name $clean_up
