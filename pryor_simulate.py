#! /usr/bin/env python
# -*- coding: utf-8 -*-
import subprocess
import random
import numpy as np

from utils import conn_pattern_filename

base_name = "10_f.3"
size = 10

bias_values = np.unique(np.rint(np.linspace(50,-10,3)).astype(np.int)) #this somewhat complex operation is required to ensure that the bias values are unique integers (using integers is somewhat of an arbitrary constraint, and the need for a uniqueness check is a consequence of generating the values as floats and then rounding and casting them to integers).

# create connection pattern
n_mf = 230
n_gr = 690
n_grc_dend = 4
sim_path = "/home/ucbtepi/nC_projects/if_gl/simulations/"
random.seed(1234)
conn_pattern = [random.sample(range(n_mf), n_grc_dend) for each in range(n_gr)]
random.seed()

# record the connection pattern in a text file
conn_pattern_file=open(sim_path+conn_pattern_filename(base_name),"w")
for gr in range(n_gr):
    for mf in conn_pattern[gr]:
        conn_pattern_file.write(str(mf) + " ")
    conn_pattern_file.write("\n")
conn_pattern_file.close()

for bias in bias_values:
    for rank in range(size):
        print 'pryor_run.py', base_name, str(bias), str(size), str(rank)
        subprocess.call(['qsub', '/home/ucbtepi/code/network/trunk/simulate_jobscript.sh', base_name, str(bias), str(size), str(rank)])

