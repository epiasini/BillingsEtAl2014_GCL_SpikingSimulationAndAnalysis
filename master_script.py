#! /usr/bin/env python
# -*- coding: utf-8 -*-
import socket

from utils.queue import BatchManager
from utils.parameters import ParameterSpace
from utils.parameters import PSlice as psl

#++++general controls+++++
force_rerun_simulations = False
clean_up_simulation_files = True

simulate = True
compress = True
analyse = False

#+++++parameter ranges+++++++++++++
n_grc_dend = psl(1,21,1)
connectivity_rule = psl(0) # 0: tissue model, 1: random bipartite graph
input_spatial_correlation_scale = psl(0) # 0: uncorrelated
active_mf_fraction = psl(.1,1.,.1)
gaba_scale = psl(1)
dta = psl(0.3)
modulation_frequency = psl(0)
stim_rate_mu = psl(80)
stim_rate_sigma = psl(0)
noise_rate_mu = psl(10)
noise_rate_sigma = psl(0)
n_stim_patterns = psl(128)
n_trials = psl(1) # typically, 1 when simulating and 60 when analysing
sim_duration = psl(3750.) # typically, 3750 when simulating and 180 when analysing
ana_duration = psl(30.0) # when analysing, must be < min(sim_duration)
training_size = psl(30) # when analysing, must be < min(n_trials)
multineuron_metric_mixing = psl(0.)
linkage_method = psl(1) # 0: ward, 1: kmeans
tau = psl(5)
dt = psl(2)

#----parameter consistency check
if analyse and training_size.realstop > n_trials.start:
    raise Exception("Decoder training set size must always be smaller than the number of trials!")
if analyse and ana_duration.realstop > sim_duration.start:
    raise Exception("Simulation length must always be greater than analysis time window!")

#----determine which cluster we're runnning on
hostname = socket.gethostname()
if 'pryor' in hostname:
    system = 'matlem'
elif 'login' in hostname:
    system = 'legion'
else:
    raise Exception("This script is designed to be used only on matlem or legion!")

#---parameter space creation
parameter_space = ParameterSpace(n_grc_dend,
                                 connectivity_rule,
                                 input_spatial_correlation_scale,
                                 active_mf_fraction,
                                 gaba_scale,
                                 dta,
                                 modulation_frequency,
                                 stim_rate_mu,
                                 stim_rate_sigma,
                                 noise_rate_mu,
                                 noise_rate_sigma,
                                 n_stim_patterns,
                                 n_trials,
                                 sim_duration,
                                 ana_duration,
                                 training_size,
                                 multineuron_metric_mixing,
                                 linkage_method,
                                 tau,
                                 dt)

batch_manager = BatchManager(parameter_space, system=system)

############################################
##====SIMULATION AND COMPRESSION STAGE====##
############################################
if simulate and compress:
    print("Submitting simulation and compression jobs")
    batch_manager.start_simulation_and_compression(force=force_rerun_simulations,
                                                   clean_up=clean_up_simulation_files)
else:
    if simulate:
        print("Submitting simulation jobs")
        batch_manager.start_simulation(force=force_rerun_simulations)
    if compress:
        print("Submitting compression jobs")
        batch_manager.start_compression(clean_up=clean_up_simulation_files)
##########################
##====ANALYSIS STAGE====##
##########################
if analyse:
    print("Submitting analysis jobs.")
    batch_manager.start_analysis()

print("All jobs submitted.")
