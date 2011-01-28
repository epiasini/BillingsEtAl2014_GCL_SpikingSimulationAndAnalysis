#!/bin/bash
base_name=example
size=2
qsub -V -S /bin/bash -o $HOME/data/eugenio/network/log/ -e $HOME/data/eugenio/network/log/ -t 1-$size simulate_call.sh $base_name $size
