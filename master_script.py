#! /usr/bin/env python
# -*- coding: utf-8 -*-
import os
import time
import random
import itertools
from math import factorial

from utils.paths import data_folder_path_ctor, data_archive_path_ctor, conn_pattern_filename, stim_pattern_filename
from utils.queue import ProcessManager

#+++++Exceptions+++++
class NetworkSizeError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class TrainingSetSizeError(Exception):
    pass

#++++general controls+++++
simulate = True
compress = True
analyse = True


#+++++fixed parameters+++++++
project_path = '/home/ucbtepi/nC_projects/if_gl' # hardcoded in simulate.py
sim_path = '/home/ucbtepi/nC_projects/if_gl/simulations' # hardcoded in simulate.py
project_filename = 'if_gl.ncx' # hardcoded in simulate.py
sim_config_name = 'Default Simulation Configuration' # hardcoded in simulate.py
nC_seed = 1234 # hardcoded in simulate.py
sim_duration = 300.0 # hardcoded in simulate.py
n_stim_patterns = 20
n_trials = 200
min_mf_number = 6
grc_mf_ratio = 2.
#+++++parameter ranges+++++++++++++
n_grc_dend_range = [4]
network_scale_range = [5]
active_mf_fraction_range = [.1, .2, .3, .4, .5, .6, .7, .8, .9]
bias_range = [0., -5., -10., -15., -20., -25., -30.]
stim_rate_mu_range = [120]
stim_rate_sigma_range = [30]
noise_rate_mu_range = [40, 60, 70]
noise_rate_sigma_range = [10]
#++++++++++++++++++++++++++
###analysis
an_tau = 5
an_dt = 2
multineuron_metric_mixing_range = [0.]
training_size_range = [40]
linkage_methods_range = ['ward']



ranges = [n_grc_dend_range, network_scale_range, active_mf_fraction_range, bias_range, stim_rate_mu_range, stim_rate_sigma_range, noise_rate_mu_range, noise_rate_sigma_range]
parameter_space = itertools.product(*ranges)

#----parallel settings----
size_per_simulation = 20

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
if simulate:
    for sim_dict in master_list:
        n_grc_dend = sim_dict['params'][0]
        scale = sim_dict['params'][1]
        active_mf_fraction = sim_dict['params'][2]
        bias = sim_dict['params'][3]
        stim_rate_mu = sim_dict['params'][4]
        stim_rate_sigma = sim_dict['params'][5]
        noise_rate_mu = sim_dict['params'][6]
        noise_rate_sigma = sim_dict['params'][7]
        # prepare directory tree
        try:
            os.makedirs(data_folder_path_ctor(grc_mf_ratio, n_grc_dend, scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma))
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
        data_archive_path = data_archive_path_ctor(grc_mf_ratio, n_grc_dend, scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma, n_stim_patterns, n_trials)
        if not os.path.exists(data_archive_path):
            for rank in range(size_per_simulation):
                qsub_argument_list = itertools.chain(['simulate_jobscript.sh', str(min_mf_number), str(grc_mf_ratio)], [str(x) for x in sim_dict['params']], [str(n_stim_patterns), str(n_trials), str(size_per_simulation), str(rank)])
                process_manager.submit_job(qsub_argument_list)
                
    while process_manager.queue_is_not_empty():
        # check for jobs that may have stalled due to the dreaded ConcurrentModificationException by parsing the error log files, since I don't seem to be able to catch that exception at the jython level.
        process_manager.update_jobs_and_check_for_CME()
        process_manager.update_prequeue()
        print('{rj} running, {wj} waiting, {oj} other jobs, {pqj} in the pre-queue'.format(rj=len(process_manager.running_jobs), wj=len(process_manager.waiting_jobs), oj=len(process_manager.other_jobs), pqj=process_manager.get_prequeue_length()))
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
        data_archive_path = data_archive_path_ctor(grc_mf_ratio, n_grc_dend, scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma, n_stim_patterns, n_trials)
        if not os.path.exists(data_archive_path):
            qsub_argument_list = itertools.chain(['compress_jobscript.sh', str(min_mf_number), str(grc_mf_ratio)], [str(x) for x in sim_dict['params']], [str(n_stim_patterns), str(n_trials), str(int(clean_up_after_compress))])
            process_manager.submit_job(qsub_argument_list)            
    
    while process_manager.queue_is_not_empty():
        process_manager.update_job_sets()
        process_manager.update_prequeue()        
        print('{rj} running, {wj} waiting, {oj} other jobs, {pqj} in the pre-queue'.format(rj=len(process_manager.running_jobs), wj=len(process_manager.waiting_jobs), oj=len(process_manager.other_jobs), pqj=process_manager.get_prequeue_length()))
        time.sleep(60)
    print("Compression stage seems to be complete.")

##########################
##====ANALYSIS STAGE====##
##########################
if analyse:
    print("Entering analysis stage.")
    #----parameter consistency control
    if any([s >= n_trials for s in training_size_range]):
        raise TrainingSetSizeError()
    
    an_ranges = [[min_mf_number], [grc_mf_ratio], n_grc_dend_range, network_scale_range, active_mf_fraction_range, bias_range, stim_rate_mu_range, stim_rate_sigma_range, noise_rate_mu_range, noise_rate_sigma_range, [n_stim_patterns], [n_trials], [sim_duration], [an_tau], [an_dt], multineuron_metric_mixing_range, training_size_range, linkage_methods_range]
    an_parameter_space = [[str(value) for value in combination] for combination in itertools.product(*an_ranges)]

    for k,par_combination in enumerate(an_parameter_space):
        qsub_argument_list = itertools.chain(['analyse_jobscript.sh'], par_combination)
        process_manager.submit_job(qsub_argument_list)
        
    while process_manager.queue_is_not_empty():
        process_manager.update_job_sets()
        process_manager.update_prequeue()        
        print('{rj} running, {wj} waiting, {oj} other jobs, {pqj} in the pre-queue'.format(rj=len(process_manager.running_jobs), wj=len(process_manager.waiting_jobs), oj=len(process_manager.other_jobs), pqj=process_manager.get_prequeue_length()))
        time.sleep(60)
    print("Analysis stage complete.")

print("Master script execution terminated.")
