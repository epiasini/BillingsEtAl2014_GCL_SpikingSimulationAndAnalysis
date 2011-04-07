#! /usr/bin/env python
# -*- coding: utf-8 -*-
import subprocess

base_name = "50_b0_f.3"
size = 25

for rank in range(size):
    subprocess.call(['qsub', '-V', '-S', '/bin/bash', '-o', '$HOME/log/', '-e', '$HOME/log/', '/home/ucbtepi/code/network/trunk/simulate_call.sh', base_name, str(size), str(rank)])
