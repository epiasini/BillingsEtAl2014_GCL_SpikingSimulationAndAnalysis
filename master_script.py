#! /usr/bin/env python
# -*- coding: utf-8 -*-
import time
from math import factorial

from utils.queue import BatchManager
from utils.parameters import ParameterSpace
from utils.parameters import PSlice as psl

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

#+++++parameter ranges+++++++++++++
sim_duration = psl(300.0)
min_mf_number = psl(6)
grc_mf_ratio = psl(2.)
n_grc_dend = psl(4)
network_scale = psl(5.)
active_mf_fraction = psl(.1, 1., .1)
bias = psl(-30., 5., 5.)
stim_rate_mu = psl(120)
stim_rate_sigma = psl(30)
noise_rate_mu = psl(10, 80, 10)
noise_rate_sigma = psl(10)
n_stim_patterns = psl(20) # must be > SimpleParameterPoint.SIZE_PER_SIMULATION
n_trials = psl(200)
training_size = psl(40) # must be < min(n_trials)
multineuron_metric_mixing = psl(0.)
linkage_method = psl(0) # 0: ward
tau = psl(5)
dt = psl(2)

#----parameter consistency check
min_nmf = round((network_scale.start)*(min_mf_number.start))
min_active_mf = round(min_nmf*(active_mf_fraction.start))
max_active_mf = round(min_nmf*(active_mf_fraction.realstop))
min_coding_inputs = min(min_active_mf, min_nmf - max_active_mf)
max_patterns_encoded = float(factorial(min_nmf))/(factorial(min_coding_inputs)*factorial(min_nmf - min_coding_inputs))
if n_stim_patterns.realstop > max_patterns_encoded:
    raise NetworkSizeError("Network size inferior limit too small for chosen number of patterns! %d MFs can't represent %d patterns for at least one of the requested sparsity values." % (min_nmf, n_stim_patterns.realstop))


#---parameter space creation
parameter_space = ParameterSpace(sim_duration,
                                 min_mf_number,
                                 grc_mf_ratio,
                                 n_grc_dend,
                                 network_scale,
                                 active_mf_fraction,
                                 bias,
                                 stim_rate_mu,
                                 stim_rate_sigma,
                                 noise_rate_mu,
                                 noise_rate_sigma,
                                 n_stim_patterns,
                                 n_trials,
                                 training_size,
                                 multineuron_metric_mixing,
                                 linkage_method,
                                 tau,
                                 dt)

batch_manager = BatchManager(parameter_space)

############################
##====SIMULATION STAGE====##
############################
print("Entering simulation stage.")
batch_manager.start_simulation(force=False)
while batch_manager.simulation.queue_is_not_empty():
    time.sleep(60)
    batch_manager.update_status()
print("Simulation stage complete. Entering compression stage.")
#############################
##====COMPRESSION STAGE====##
#############################
batch_manager.start_compression()
while batch_manager.compression.queue_is_not_empty():
    time.sleep(60)
    batch_manager.update_status()
print("Compression stage complete. Entering analysis stage.")
##########################
##====ANALYSIS STAGE====##
##########################
batch_manager.start_analysis()
while batch_manager.analysis.queue_is_not_empty():
    time.sleep(60)
    batch_manager.update_status()
print("Analysis stage complete. Master script execution terminated.")
