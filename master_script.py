#! /usr/bin/env python
# -*- coding: utf-8 -*-
import os
import time
import random
import itertools
from math import factorial
from subprocess import Popen, PIPE, call

from utils.paths import data_folder_path_ctor, data_archive_path_ctor, conn_pattern_filename, stim_pattern_filename

class QueueError(Exception):
    pass

class NetworkSizeError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class ProcessManager(object):
    def __init__(self):
        self.__managed_jobs = set() # all the jobs in the queue that are being managed by this manager
        self.running_jobs = set() # all (managed) running jobs
        self.waiting_jobs = set() # all (managed) waiting jobs
        self.other_jobs = set() # all the other managed jobs (e.g. jobs that are being transferred)
        self.__qsub_commands = dict() # historical record of the qsub commands that generated each managed job
    def __managed_jids_from_qstat(self, *qstat_argument_list):
        stat_lines = Popen(itertools.chain(['qstat'], qstat_argument_list), stdout=PIPE, stderr=PIPE).communicate()[0].split('\n')
        return set(int(l.split()[0]) for l in stat_lines[2:] if len(l)>0).intersection(self.__managed_jobs)
        
    def submit_job(self, qsub_argument_list):
        popen_command = itertools.chain(['qsub'], qsub_argument_list)
        handle = Popen(popen_command, stdin=PIPE, stdout=PIPE, stderr=PIPE)
#        if handle.returncode!=0:
#            raise QueueError()
        jid = int((handle.communicate()[0]).split(' ')[2])
        self.__qsub_commands[jid] = qsub_argument_list
        self.__managed_jobs.add(jid)
        self.update_job_sets()
        return jid
    def delete_job(self, jid):
        call(['qdel', str(jid)])
        self.update_job_sets()
    def update_job_sets(self):
        self.running_jobs = self.__managed_jids_from_qstat('-s', 'r')
        self.waiting_jobs = self.__managed_jids_from_qstat('-s', 'p')
        self.other_jobs = set(jid for jid in self.__managed_jids_from_qstat() if jid not in self.running_jobs.union(self.waiting_jobs))
    def update_and_check_for_CME(self):
        self.update_job_sets()
        for jid in self.running_jobs:
            try:
                with open('/home/ucbtepi/log/simulate_jobscript.sh.e{jid}'.format(jid=jid), 'r') as f:
                    if 'ConcurrentModificationException' in f.read():
                        self.delete_job(jid)
                        new_jid = self.submit_job(self.__qsub_commands[jid])
                        print('Found ConcurrentModificationException in job {old_jid}. Deleting job and resubmitted as {new_jid}.'.format(old_jid=jid, new_jid=new_jid))
            except IOError as (errno, strerror):
                print('I/O Error {errno}:{strerror}'.format(errno=errno, strerror=strerror))

#++++general controls+++++
simulate = True
compress = True
analyse = False

#+++++fixed parameters+++++++
project_path = '/home/ucbtepi/nC_projects/if_gl' # hardcoded in simulate.py
sim_path = '/home/ucbtepi/nC_projects/if_gl/simulations' # hardcoded in simulate.py
project_filename = 'if_gl.ncx' # hardcoded in simulate.py
sim_config_name = 'Default Simulation Configuration' # hardcoded in simulate.py
nC_seed = 1234 # hardcoded in simulate.py
sim_duration = 300.0 # hardcoded in simulate.py
n_stim_patterns = 20
n_trials = 500
min_mf_number = 6
grc_mf_ratio = 3.
#+++++parameter ranges+++++++++++++
n_grc_dend_range = [4]
network_scale_range = [1.66]
active_mf_fraction_range = [0.5]
bias_range = [0., -10., -25., -45.]
#++++++++++++++++++++++++++

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
process_manager = ProcessManager()
############################
##====SIMULATION STAGE====##
############################
initial_jobs = set()

if simulate:
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
            stim_pattern_file.close()

        # submit simulations to the queue
        data_archive_path = data_archive_path_ctor(grc_mf_ratio, n_grc_dend, scale, active_mf_fraction, bias, n_stim_patterns, n_trials)
        if not os.path.exists(data_archive_path):
            for rank in range(size_per_simulation):
                qsub_argument_list = itertools.chain(['simulate_jobscript.sh', str(min_mf_number), str(grc_mf_ratio)], [str(x) for x in sim_dict['params']], [str(n_stim_patterns), str(n_trials), str(size_per_simulation), str(rank)])
                process_manager.submit_job(qsub_argument_list)
                
    while process_manager.running_jobs or process_manager.waiting_jobs or process_manager.other_jobs:
        # check for jobs that may have stalled due to the dreaded ConcurrentModificationException by parsing the error log files, since I don't seem to be able to catch that exception at the jython level.
        process_manager.update_and_check_for_CME()
        print('{rj} running, {wj} waiting, {oj} other jobs'.format(rj=len(process_manager.running_jobs), wj=len(process_manager.waiting_jobs), oj=len(process_manager.other_jobs)))
        time.sleep(60)
    print("Simulation stage complete.")

#############################
##====COMPRESSION STAGE====##
#############################
if compress:
    print("Entering compression stage.")
    clean_up_after_compress = True
    initial_compr_jobs = set()
    
    for sim_dict in master_list:
        n_grc_dend = sim_dict['params'][0]
        scale = sim_dict['params'][1]
        active_mf_fraction = sim_dict['params'][2]
        bias = sim_dict['params'][3]
        data_archive_path = data_archive_path_ctor(grc_mf_ratio, n_grc_dend, scale, active_mf_fraction, bias, n_stim_patterns, n_trials)
        if not os.path.exists(data_archive_path):
            qsub_argument_list = itertools.chain(['compress_jobscript.sh', str(min_mf_number), str(grc_mf_ratio)], [str(x) for x in sim_dict['params']], [str(n_stim_patterns), str(n_trials), str(int(clean_up_after_compress))])
            process_manager.submit_job(qsub_argument_list)            
            #popen_command = itertools.chain(['qsub', 'compress_jobscript.sh'], [str(min_mf_number), str(grc_mf_ratio)], [str(x) for x in sim_dict['params']], [str(n_stim_patterns), str(n_trials)], [str(int(clean_up_after_compress))])
    
    while process_manager.running_jobs or process_manager.waiting_jobs or process_manager.other_jobs:
        process_manager.update_job_sets()
        print('{rj} running, {wj} waiting, {oj} other_jobs'.format(rj=len(process_manager.running_jobs), wj=len(process_manager.waiting_jobs), oj=len(process_manager.other_jobs)))
        time.sleep(60)
    print("Compression stage seems to be complete.")

##########################
##====ANALYSIS STAGE====##
##########################
if analyse:
    print("Entering analysis stage.")
    initial_an_jobs = set()
    for sim_dict in master_list:
        n_grc_dend = sim_dict['params'][0]
        scale = sim_dict['params'][1]
        active_mf_fraction = sim_dict['params'][2]
        bias = sim_dict['params'][3]
        data_archive_path = data_archive_path_ctor(grc_mf_ratio, n_grc_dend, scale, active_mf_fraction, bias, n_stim_patterns, n_trials)
        an_result_path = an_result_path_ctor(grc_mf_ratio, n_grc_dend, scale, active_mf_fraction, bias, n_stim_patterns, n_trials)
        popen_command = itertools.chain(['qsub', 'hcluster_jobscript.sh'], [data_archive_path, an_result_path])

        if not os.path.exists(an_result_path):
            handle = Popen(popen_command, stdin=PIPE, stdout=PIPE, stderr=PIPE)
            sim_dict['an_qsub_handle'] = handle
            jid = int((handle.communicate()[0]).split(' ')[2])
            sim_dict['an_jid'] = jid
            initial_an_jobs.add(jid)


    active_an_jobs = set()
    waiting_an_jobs = set(initial_an_jobs)
    other_an_jobs = set()

    if any(s['an_qsub_handle'].returncode!=0 for s in master_list if 'an_qsub_handle' in s.keys()):
        raise QueueError()
    
    while active_an_jobs or waiting_an_jobs or other_an_jobs:
        time.sleep(10)
        active_an_jobs, waiting_an_jobs, other_an_jobs = update_job_sets(initial_an_jobs)
        print(active_an_jobs)
        print(waiting_an_jobs)

print("Master script execution terminated.")
