#! /usr/bin/env python
# -*- coding: utf-8 -*-
import os
import time
import random
import itertools
from math import factorial
from subprocess import Popen, PIPE

from utils import data_folder_path_ctor, data_archive_path_ctor, conn_pattern_filename, stim_pattern_filename

class QueueError(Exception):
    pass

class NetworkSizeError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def update_job_sets(initial_jobs_set):
    stat_running = Popen(['qstat', '-s', 'r'], stdout=PIPE, stderr=PIPE).communicate()[0]
    stat_waiting = Popen(['qstat', '-s', 'p'], stdout=PIPE, stderr=PIPE).communicate()[0]
    active_jobs_set = set([j for j in initial_jobs_set if str(j) in stat_running])
    waiting_jobs_set = set([j for j in initial_jobs_set if str(j) in stat_waiting])
    return active_jobs_set, waiting_jobs_set

#+++++fixed parameters+++++++
project_path = '/home/ucbtepi/nC_projects/if_gl' # hardcoded in simulate.py
sim_path = '/home/ucbtepi/nC_projects/if_gl/simulations' # hardcoded in simulate.py
project_filename = 'if_gl.ncx' # hardcoded in simulate.py
sim_config_name = 'Default Simulation Configuration' # hardcoded in simulate.py
nC_seed = 1234 # hardcoded in simulate.py
sim_duration = 300.0 # hardcoded in simulate.py
n_stim_patterns = 20
n_trials = 50
min_mf_number = 6
grc_mf_ratio = 3.
#+++++parameter ranges+++++++++++++
n_grc_dend_range = [3]
network_scale_range = [1.00, 1.33, 1.66]
active_mf_fraction_range = [0.4, 0.5]
bias_range = [20., 10., 0., -10., -20.]
#++++++++++++++++++++++++++

###testing
n_stim_patterns = 10
n_trials = 15

ranges = [n_grc_dend_range, network_scale_range, active_mf_fraction_range, bias_range]
parameter_space = itertools.product(*ranges)

#----parallel settings----
size_per_simulation = 10

#----parameter consistency control
min_nmf = round(min(network_scale_range)*min_mf_number)
min_active_mf = round(min_nmf*min(active_mf_fraction_range))
max_active_mf = round(min_nmf*max(active_mf_fraction_range))
min_coding_inputs = min(min_active_mf, min_nmf - max_active_mf)
max_patterns_encoded = float(factorial(min_nmf))/(factorial(min_coding_inputs)*factorial(min_nmf - min_coding_inputs))
if n_stim_patterns > max_patterns_encoded:
    raise NetworkSizeError("Network size inferior limit too small for chosen number of patterns! %d MFs can't represent %d patterns for at least one of the requested sparsity values." % (min_nmf, n_stim_patterns))

master_list = [{'params': parameter_point} for parameter_point in parameter_space]
initial_jobs = set()

for sim_dict in master_list:
    n_grc_dend = sim_dict['params'][0]
    scale = sim_dict['params'][1]
    active_mf_fraction = sim_dict['params'][2]
    bias = sim_dict['params'][3]
    # prepare directory tree
    try:
        os.makedirs(data_folder_path_ctor(grc_mf_ratio, n_grc_dend, scale, active_mf_fraction, bias))
    except OSError:
        # this means that the directory is already there
        pass
    
    # create connection pattern and save it in a text file
    conn_pattern_path = conn_pattern_filename(grc_mf_ratio, n_grc_dend, scale)
    n_mf = int(round(min_mf_number * scale))
    n_gr = int(round(n_mf * grc_mf_ratio))
    if not os.path.exists(conn_pattern_path):
        conn_pattern = [random.sample(range(n_mf), n_grc_dend) for each in range(n_gr)]
        conn_pattern_file = open(conn_pattern_path, "w")
        for gr in range(n_gr):
            for mf in conn_pattern[gr]:
                conn_pattern_file.write(str(mf) + " ")
            conn_pattern_file.write("\n")
        conn_pattern_file.close()

    # generate random stimulation patterns and save them in a text file
    stim_pattern_path = stim_pattern_filename(grc_mf_ratio, n_grc_dend, scale, active_mf_fraction, n_stim_patterns)
    if not os.path.exists(stim_pattern_path):
        active_mf_number = int(round(n_mf*active_mf_fraction))
        stim_patterns = []
        stim_pattern_file = open(stim_pattern_path, "w")
        for spn in range(n_stim_patterns):
            while True:
                sp = sorted(random.sample(range(n_mf), active_mf_number))
                if sp not in stim_patterns:
                    break
            stim_patterns.append(sp)
            for mf in sp:
                stim_pattern_file.write(str(mf) + " ")
            stim_pattern_file.write("\n")

    # submit simulations to the queue
    data_archive_path = data_archive_path_ctor(grc_mf_ratio, n_grc_dend, scale, active_mf_fraction, bias, n_stim_patterns, n_trials)
    if not os.path.exists(data_archive_path):
        sim_dict['jids'] = []
        sim_dict['sim_qsub_handles'] = []
        for rank in range(size_per_simulation):
            popen_command = itertools.chain(['qsub', 'simulate_jobscript.sh'], [str(min_mf_number), str(grc_mf_ratio)], [str(x) for x in sim_dict['params']], [str(n_stim_patterns), str(n_trials)], [str(size_per_simulation), str(rank)])
            handle = Popen(popen_command, stdin=PIPE, stdout=PIPE, stderr=PIPE)
            sim_dict['sim_qsub_handles'].append(handle)
            jid = int((handle.communicate()[0]).split(' ')[2])
            sim_dict['jids'].append(jid)
            initial_jobs.add(jid)

active_jobs = set()
waiting_jobs = set(initial_jobs)

#print('Job IDs are: %s' % str(jids))

if any(h.returncode!=0 for s in master_list for h in s['sim_qsub_handles']):
    raise QueueError()

while active_jobs or waiting_jobs:
    #print([(s['params'], [j for j in s['jids'] if j in active_jobs]) for s in master_list if any(j in active_jobs for j in s['jids'])])
    time.sleep(10)
    active_jobs, waiting_jobs = update_job_sets(initial_jobs)
    print(active_jobs)
    print(waiting_jobs)

print("Jobs seem to be completed.")
    
    
