#!/bin/bash
base_name=$1
clean_up=$2
cd /home/ucgbgbi/data/eugenio/network/trunk
/usr/bin/time /home/ucgbgbi/data/eugenio/bin/python compress_data.py $base_name $clean_up


