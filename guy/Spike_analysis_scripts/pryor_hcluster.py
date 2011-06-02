#! /usr/bin/env python
# -*- coding: utf-8 -*-
import subprocess
import numpy as np

base_name = "20_f.5_s1.00"

conf_path = '/home/ucbtepi/code/network/trunk/'

# read the configuration file and extract the variables that will be used
conf_file = open(conf_path+base_name+'.conf.txt')
conf = eval(conf_file.read())
conf_file.close()

bias_settings = conf['bias_settings']
network_scale = conf['network_scale']

bias_values = np.unique(np.rint(np.linspace(bias_settings['start'],bias_settings['stop'],bias_settings['num'])).astype(np.int)) #this somewhat complex operation is required to ensure that the bias values are unique integers (usin integers is somewhat of an arbitrary constraint, and the need for a uniqueness check is a consequence of generating the values as floats and then rounding and casting them to integers).

for bias in bias_values:
    subprocess.call(['qsub', '/home/ucbtepi/code/network/trunk/guy/Spike_analysis_scripts/hcluster_jobscript.sh', str(network_scale), str(bias)])
