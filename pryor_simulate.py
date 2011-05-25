#! /usr/bin/env python
# -*- coding: utf-8 -*-
import subprocess
import random
import numpy as np

from itertools import combinations

from utils import conn_pattern_filename, stim_pattern_filename

base_name = "20_f.5_s1.33"
size = 10

conf_path = '/home/ucbtepi/code/network/trunk/'

# read the configuration file and extract the variables that will be used
conf_file = open(conf_path+base_name+'.conf.txt')
conf = eval(conf_file.read())
conf_file.close()

network_scale = conf['network_scale']
grc_mf_ratio = conf['grc_mf_ratio']
min_mf_number = conf['min_mf_number']
n_grc_dend = conf['n_grc_dend']
active_mf_fraction = conf['active_mf_fraction']
n_stim_patterns = conf['n_stim_patterns']
sim_path = conf['sim_path']
bias_settings = conf['bias_settings']

bias_values = np.unique(np.rint(np.linspace(bias_settings['start'],bias_settings['stop'],bias_settings['num'])).astype(np.int)) #this somewhat complex operation is required to ensure that the bias values are unique integers (using integers is somewhat of an arbitrary constraint, and the need for a uniqueness check is a consequence of generating the values as floats and then rounding and casting them to integers).
n_mf = int(round(min_mf_number * network_scale))
n_gr = int(round(n_mf * grc_mf_ratio))

# create connection pattern and save it in a text file
conn_pattern = [random.sample(range(n_mf), n_grc_dend) for each in range(n_gr)]
conn_pattern_file = open(sim_path+conn_pattern_filename(base_name),"w")
for gr in range(n_gr):
    for mf in conn_pattern[gr]:
        conn_pattern_file.write(str(mf) + " ")
    conn_pattern_file.write("\n")
conn_pattern_file.close()

# generate random stimulation patterns and save them in a text file
active_mf_number = int(round(n_mf*active_mf_fraction))
stim_patterns = []
stim_pattern_file = open(sim_path + stim_pattern_filename(base_name), "w")
for spn in range(n_stim_patterns):
    while True:
        sp = sorted(random.sample(range(n_mf), active_mf_number))
        if sp not in stim_patterns:
            break
    stim_patterns.append(sp)
    for mf in sp:
        stim_pattern_file.write(str(mf) + " ")
    stim_pattern_file.write("\n")

for bias in bias_values:
    for rank in range(size):
        print 'pryor_run.py', base_name, str(bias), str(size), str(rank)
        subprocess.call(['qsub', '/home/ucbtepi/code/network/trunk/simulate_jobscript.sh', base_name, str(bias), str(size), str(rank)])

