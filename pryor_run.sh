#!/bin/bash
qsub -V -S /bin/bash -o $HOME/data/eugenio/network/log/ -e $HOME/data/eugenio/network/log/ simulate_call.sh
