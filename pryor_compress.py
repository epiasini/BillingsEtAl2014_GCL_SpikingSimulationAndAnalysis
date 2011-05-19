#! /usr/bin/env python
# -*- coding: utf-8 -*-
import subprocess
import numpy as np

base_name = "10_f.3"

bias_values = np.unique(np.rint(np.linspace(50,-10,3)).astype(np.int)) #this somewhat complex operation is required to ensure that the bias values are unique integers (usin integers is somewhat of an arbitrary constraint, and the need for a uniqueness check is a consequence of generating the values as floats and then rounding and casting them to integers).

for bias in bias_values:
    print 'pryor_compress.py', base_name, str(bias)
    subprocess.call(['qsub', '/home/ucbtepi/code/network/trunk/compress_jobscript.sh', base_name, str(bias)])
