#! /usr/bin/env python
# -*- coding: utf-8 -*-
import sys

from utils.analysis import analyse_single_configuration

#+++++Exceptions+++++
class TrainingSetSizeError(Exception):
    pass

min_mf_number = int(sys.argv[1])
grc_mf_ratio = float(sys.argv[2])
n_grc_dend = int(sys.argv[3])
network_scale = float(sys.argv[4])
active_mf_fraction = float(sys.argv[5])
bias = float(sys.argv[6])
n_stim_patterns = int(sys.argv[7])
n_trials = int(sys.argv[8])
sim_duration = 300.0 # hardcoded in simulate.py
tau = float(sys.argv[9])
dt = float(sys.argv[10])
multineuron_metric_mixing = float(sys.argv[11])
training_size = int(sys.argv[12])
linkage_method = sys.argv[13]

#----parameter consistency control
if any([s >= min(n_trials_range) for s in training_size_range]):
    raise TrainingSetSizeError()


analyse_single_configuration(min_mf_number, grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, n_stim_patterns, n_trials, sim_duration, tau, dt, multineuron_metric_mixing, training_size, linkage_method)






    
