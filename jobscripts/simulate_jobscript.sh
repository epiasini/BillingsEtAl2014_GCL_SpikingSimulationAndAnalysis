#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -o $HOME/log/
#$ -e $HOME/log/
#$ -l h_vmem=4G
#$ -l tmem=4G
#$ -l h_rt=48:00:00
#$ -l s_stack=10M
#$ -l h_stack=15M
##$ -l virtual_free=1024M

# "jobscripts" are things that should be passed to qsub.

# networkxj_dir should point to a copy of the networkxj package, This
# is a fork of the last networkx version (1.2) made to be compatible
# with Python 2.5, which is the Python version implemented by the
# latest stable Jython version as of december 2013.
networkxj_dir=/home/ucbtepi/code/networkx
# nC_dir should point to an up-to-date nC installation
nC_dir=/home/ucbtepi/code/neuroml_dev/neuroConstruct

args_list=$@
startdir=`pwd`

echo $args_list
export JYTHONPATH="$networkxj_dir:$JYTHONPATH"

hostname
date
cd $nC_dir
/usr/bin/time ./nC.sh -python $startdir/simulate.py $args_list


