#! /usr/bin/env python
# -*- coding: utf-8 -*-
import subprocess
import numpy as np

base_name = "50_f.3"
size = 25

bias_values = np.unique(np.rint(np.linspace(0,-2,3)).astype(np.int)) #this somewhat complex operation is required to ensure that the bias values are unique integers (usin integers is somewhat of an arbitrary constraint, and the need for a uniqueness check is a consequence of generating the values as floats and then rounding and casting them to integers).

for bias in bias_values:
    for rank in range(size):
        print 'pryor_run.py', base_name, str(bias), str(size), str(rank)
        subprocess.call(['qsub', '/home/ucbtepi/code/network/trunk/simulate_jobscript.sh', base_name, str(bias), str(size), str(rank)])
