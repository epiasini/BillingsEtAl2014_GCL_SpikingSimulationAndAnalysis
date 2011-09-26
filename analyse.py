#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import h5py
import itertools
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram

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
network_scale_range = [1.66]
active_mf_fraction_range = [0.5]
bias_range = [0.]
n_trials_range = [200]
training_size_range = [40]
multineuron_metric_mixing_range = [0.]
linkage_methods_range = ['ward']
#++++++++++++++++++++++++++

#----parameter consistency control
if any([s >= min(n_trials_range) for s in training_size_range]):
    raise TrainingSetSizeError()

ranges = [n_grc_dend_range, network_scale_range, active_mf_fraction_range, bias_range, n_trials_range, training_size_range, multineuron_metric_mixing_range, linkage_methods_range]
parameter_space = list(itertools.product(*ranges))




colors = 'bgrcmyk'
fig = plt.figure()
ax = fig.add_subplot(111)
#fig2 = plt.figure()
#ax2 = fig2.add_subplot(111)
dend_fig = plt.figure()
dend_axes = []

info_at_npatterns = [[None for each in bias_range] for each in network_scale_range]

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
    tr_direct_mi, ts_decoded_mi_plugin, ts_decoded_mi_qe, decoder_precision, tr_tree = analyse_single_configuration(min_mf_number, grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, n_stim_patterns, n_trials, sim_duration, tau, dt, multineuron_metric_mixing, training_size, linkage_method)
    info_at_npatterns[network_scale_range.index(network_scale)][bias_range.index(bias)] = ts_decoded_mi_qe[n_stim_patterns]
    # plot
    color = colors[k%len(colors)]
    if plotting_mode == 'precision':
        ax.semilogx(decoder_precision, ts_decoded_mi_plugin, label='test set size: %d, decoded (plugin)' % (n_trials-training_size), linestyle='-.', color=color)
        ax.semilogx(decoder_precision, ts_decoded_mi_qe, label='test set size: %d, decoded (qe)' % (n_trials-training_size), color=color)
        ax.plot([decoder_precision[n_stim_patterns], decoder_precision[n_stim_patterns]], ax.get_ylim(), linestyle='--', color=color)
        #ax2.semilogy(tr_direct_mi/np.log2(n_stim_patterns), decoder_precision, linestyle='', marker='+')
    elif plotting_mode == 'alphabet_size':
        ax.plot(tr_direct_mi, label='mixing%.2f, direct' % (multineuron_metric_mixing), linestyle='-.', color=color)
        ax.plot(ts_decoded_mi_qe, label='mixing%.2f, qe' % (multineuron_metric_mixing), color=color)
        ax.plot([n_stim_patterns, n_stim_patterns], ax.get_ylim(), linestyle='--', color=color)
    d_ax = dend_fig.add_subplot(1,len(parameter_space),k)
    d_ax.set_title('Network scale parameter: %.2f' % network_scale)
    d_ax.set_xlabel('Datapoint index')
    d_ax.set_ylabel('Inter-node multiunit V.R. distance (a.u.)')
    dend_axes.append(d_ax)
    dendrogram(tr_tree, color_threshold=tr_tree[-n_stim_patterns+1,2], no_labels=True)




#for k,l in enumerate(info_at_npatterns):
#    ax.plot(network_scale_range, l, label='b:%d' % (bias_range[k]))

# i_at_npatt_plot = ax2.imshow(info_at_npatterns, interpolation='nearest', cmap='hot')
# i_at_npatt_cbar = fig2.colorbar(i_at_npatt_plot)
# i_at_npatt_cbar.set_label('MI at $|\mathcal{A}_{out}| = |\mathcal{A}_{in}|$ (bits)')
# ax2.set_xticks(np.arange(len(bias_range)))
# ax2.set_xticklabels([str(x) for x in bias_range])
# ax2.set_xlabel('Threshold current (pA)')
# ax2.set_yticks(np.arange(len(network_scale_range)))
# ax2.set_yticklabels([str(y) for y in network_scale_range])
# ax2.set_ylabel('Network scale factor')
# ax2.set_title('Effect of network size and threshold current')

#ax.plot([n_stim_patterns, n_stim_patterns], ax.get_ylim(), linestyle='--', color='black')
ax.legend(loc='upper left')
if plotting_mode == 'precision':
    ax.set_xlabel('decoder precision (1/cluster separation)')
elif plotting_mode == 'alphabet_size':
    ax.set_xlabel('alphabet size (clusters in the decoder)')
ax.set_ylabel('MI (bits)')
#if not prec_thresholds[0] in ax.get_xticks():
#    ax.set_xticks(list(ax.get_xticks()) + [prec_thresholds[0]])

dend_axes_ylim = (0, max([ax.get_ylim()[1] for ax in dend_axes]))
for ax in dend_axes:
    ax.set_ylim(dend_axes_ylim)

plt.show()
    
