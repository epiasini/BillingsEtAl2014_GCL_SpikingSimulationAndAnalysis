#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import h5py
import itertools
from matplotlib import pyplot as plt

from utils.analysis import analyse_single_configuration

#+++++debugging stuff+++++
import pdb
#np.seterr(all='raise') # to convert warnings to exceptions (this allow to track the offending line)
np.seterr(divide='ignore') # to suppress 'divide by zero' warnings

#+++++Exceptions+++++
class TrainingSetSizeError(Exception):
    pass

#+++++fixed parameters+++++++
sim_duration = 300.0 # hardcoded in simulate.py
n_stim_patterns = 20
min_mf_number = 6
grc_mf_ratio = 3.
tau = 5.
dt = 2.
plotting_mode = 'precision'
#+++++parameter ranges+++++++++++++
n_grc_dend_range = [4]
network_scale_range = [1.00]
active_mf_fraction_range = [0.5]
bias_range = [20.]
n_trials_range = [200]
training_size_range = [40]
multineuron_metric_mixing_range = [0.]
linkage_methods_range = ['ward']
#++++++++++++++++++++++++++

#----parameter consistency control
if any([s >= min(n_trials_range) for s in training_size_range]):
    raise TrainingSetSizeError()

ranges = [n_grc_dend_range, network_scale_range, active_mf_fraction_range, bias_range, n_trials_range, training_size_range, multineuron_metric_mixing_range, linkage_methods_range]
parameter_space = itertools.product(*ranges)




colors = 'bgrcmyk'
fig = plt.figure()
ax = fig.add_subplot(111)

for k,par_combination in enumerate(parameter_space):
    # load parameter combination
    n_grc_dend = par_combination[0]
    network_scale = par_combination[1]
    active_mf_fraction = par_combination[2]
    bias = par_combination[3]
    n_trials = par_combination[4]
    training_size = par_combination[5]
    multineuron_metric_mixing = par_combination[6]
    linkage_method = par_combination[7]
    # analyse data
    tr_direct_mi, ts_decoded_mi_plugin, ts_decoded_mi_qe, decoder_precision = analyse_single_configuration(min_mf_number, grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, n_stim_patterns, n_trials, sim_duration, tau, dt, multineuron_metric_mixing, training_size, linkage_method)
    # plot
    color = colors[k%len(colors)]
    if plotting_mode == 'precision':
        ax.semilogx(decoder_precision, tr_direct_mi, label='method: %s, direct' % (linkage_method), linestyle='-.', color=color)
        ax.semilogx(decoder_precision, ts_decoded_mi_plugin, label='method: %s, plugin' % (linkage_method), color=color)
        #ax.plot(decoder_precision, ts_decoded_mi_qe, label='n_trials%d, train%d, qe' % (n_trials, training_size), color=color)
        ax.plot([decoder_precision[n_stim_patterns], decoder_precision[n_stim_patterns]], ax.get_ylim(), linestyle='--', color=color)
    elif plotting_mode == 'alphabet_size':
        ax.plot(tr_direct_mi, label='mixing%.2f, direct' % (multineuron_metric_mixing), linestyle='-.', color=color)
        ax.plot(ts_decoded_mi_qe, label='mixing%.2f, qe' % (multineuron_metric_mixing), color=color)
        ax.plot([n_stim_patterns, n_stim_patterns], ax.get_ylim(), linestyle='--', color=color)

#ax.plot([n_stim_patterns, n_stim_patterns], ax.get_ylim(), linestyle='--', color='black')
ax.legend(loc='upper left')
if plotting_mode == 'precision':
    ax.set_xlabel('decoder precision (1/cluster separation)')
elif plotting_mode == 'alphabet_size':
    ax.set_xlabel('alphabet size (clusters in the decoder)')
ax.set_ylabel('MI (bits)')
#if not prec_thresholds[0] in ax.get_xticks():
#    ax.set_xticks(list(ax.get_xticks()) + [prec_thresholds[0]])
#plt.show()
    
