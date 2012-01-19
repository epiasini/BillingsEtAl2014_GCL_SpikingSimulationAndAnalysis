#! /usr/bin/env python
# -*- coding: utf-8 -*-
import sys

from utils.analysis import analyse_single_configuration, analyse_synchrony

min_mf_number = int(sys.argv[1])
grc_mf_ratio = float(sys.argv[2])
n_grc_dend = int(sys.argv[3])
network_scale = float(sys.argv[4])
active_mf_fraction = float(sys.argv[5])
bias = float(sys.argv[6])
stim_rate_mu = float(sys.argv[7])
stim_rate_sigma = float(sys.argv[8])
noise_rate_mu = float(sys.argv[9])
noise_rate_sigma = float(sys.argv[10])
n_stim_patterns = int(sys.argv[11])
n_trials = int(sys.argv[12])
sim_duration = float(sys.argv[13])
tau = float(sys.argv[14])
dt = float(sys.argv[15])
multineuron_metric_mixing = float(sys.argv[16])
training_size = int(sys.argv[17])
linkage_method = sys.argv[18]

analyse_single_configuration(min_mf_number, grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma, n_stim_patterns, n_trials, sim_duration, tau, dt, multineuron_metric_mixing, training_size, linkage_method)

#analyse_synchrony(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma, n_stim_patterns, n_trials, sim_duration, tau, dt)





    
