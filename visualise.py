#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import h5py
import itertools
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram

from utils.graphic import plot_data_vector
from utils.analysis import analyse_single_configuration, cluster_centroids, kl_divergence_from_flat_p

#+++++debugging stuff+++++
import pdb
#np.seterr(all='raise') # to convert warnings to exceptions (this allow to track the offending line)
np.seterr(divide='ignore') # to suppress 'divide by zero' warnings

#+++++Exceptions+++++
class TrainingSetSizeError(Exception):
    pass

plot_mi_detail = False
plot_dendrograms = False
plot_mutual_information = True
plot_kl_divergence = False
plot_barcodes = False

#+++++fixed parameters+++++++
sim_duration = 300.0 # hardcoded in simulate.py
n_stim_patterns = 20
min_mf_number = 6
grc_mf_ratio = 2.
tau = 5.
dt = 2.
plotting_mode = 'precision'
#+++++parameter ranges+++++++++++++
n_grc_dend_range = [4]
network_scale_range = [5.]
active_mf_fraction_range = list(np.arange(.1, 1, .1))
bias_range = list(np.arange(0., -50., -5.))
n_trials_range = [500]
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
if plot_mi_detail:
    mi_det_fig = plt.figure()
    mi_det_ax = mi_det_fig.add_subplot(111)
if plot_dendrograms:
    dend_fig = plt.figure()
    dend_axes = []
if plot_mutual_information:
    info_at_npatterns = [[None for each in bias_range] for each in active_mf_fraction_range]
if plot_kl_divergence:
    kl_div_values = [[None for each in bias_range] for each in active_mf_fraction_range]
if plot_barcodes:
    centroid_figs = []
    centroid_sets = []

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
    tr_direct_mi, ts_decoded_mi_plugin, ts_decoded_mi_qe, decoder_precision, tr_tree, px_at_same_size_point = analyse_single_configuration(min_mf_number, grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, n_stim_patterns, n_trials, sim_duration, tau, dt, multineuron_metric_mixing, training_size, linkage_method)
    if plot_mutual_information:
        info_at_npatterns[active_mf_fraction_range.index(active_mf_fraction)][bias_range.index(bias)] = ts_decoded_mi_qe[n_stim_patterns]
    if plot_kl_divergence:
        kl_div_values[active_mf_fraction_range.index(active_mf_fraction)][bias_range.index(bias)] = kl_divergence_from_flat_p(q=px_at_same_size_point)
    # plot
    color = colors[k%len(colors)]
    if plot_mi_detail and plotting_mode == 'precision':
        mi_det_ax.semilogx(decoder_precision, tr_direct_mi, linestyle='-.', color=color)
        mi_det_ax.semilogx(decoder_precision, ts_decoded_mi_qe, label=r'mf act. prob: %.1f' % (active_mf_fraction), color=color)
        mi_det_ax.plot([decoder_precision[n_stim_patterns], decoder_precision[n_stim_patterns]], mi_det_ax.get_ylim(), linestyle='--', color=color)
        #mi_det_ax2.semilogy(tr_direct_mi/np.log2(n_stim_patterns), decoder_precision, linestyle='', marker='+')
    elif plotting_mode == 'alphabet_size':
        mi_det_ax.plot(np.arange(1, n_stim_patterns*training_size), tr_direct_mi, label=r'$\vartheta = %.2f^{\circ}$, direct' % (np.arccos(multineuron_metric_mixing)/np.pi*180.), linestyle='-.', color=color)
        mi_det_ax.plot(np.arange(1, n_stim_patterns*training_size), ts_decoded_mi_qe, label=r'$\vartheta = %.2f^{\circ}$, decoded' % (np.arccos(multineuron_metric_mixing)/np.pi*180.), color=color)
        if k==len(parameter_space)-1:
            mi_det_ax.plot([n_stim_patterns, n_stim_patterns], mi_det_ax.get_ylim(), linestyle='--', color='k')
    if plot_dendrograms:
        d_ax = dend_fig.add_subplot(1,len(parameter_space),k)
        d_ax.set_title('Network scale parameter: %.2f' % (network_scale))
        d_ax.set_xlabel('Datapoint index')
        d_ax.set_ylabel('Inter-node multiunit V.R. distance (a.u.)')
        dend_axes.append(d_ax)
        dendrogram(tr_tree, color_threshold=tr_tree[-n_stim_patterns+1,2], no_labels=True)

    if plot_barcodes:
        # build and plot a representation of the "typical" centroids at the |A_in|=|A_out| point.
        clust_idxs, centroids = cluster_centroids(min_mf_number, grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, n_stim_patterns, n_trials, sim_duration, tau, dt, multineuron_metric_mixing, training_size, linkage_method, n_clusts=n_stim_patterns)
        centroid_sets.append(centroids)
        centroids_fig = plt.figure()
        centroid_figs.append(centroids_fig)
        conv_plots = []
        for j,c in enumerate(centroids):
           c_ax = centroids_fig.add_subplot(max(int(np.ceil(len(centroids)/3.)), 1), 3, j)
           conv_plots.append(plot_data_vector(c_ax, c))
        clims = [plot.get_clim()[1] for plot in conv_plots]
        max_clim = max(clims)
        for plot in conv_plots:
           plot.set_clim(vmin=None, vmax=max_clim)
        centroids_fig.suptitle('Centroids at $|\mathcal{A}_{out}| = |\mathcal{A}_{in}|$.\nnet scale %.2f, act. prob. %.1f, bias %.2f, mixing %.2f ,MI at $|\mathcal{A}_{out}| = |\mathcal{A}_{in}|$ %.2f' % (network_scale, active_mf_fraction, bias, multineuron_metric_mixing, ts_decoded_mi_qe[n_stim_patterns]))


if plot_mutual_information:
    mi_fig = plt.figure()
    mi_ax = mi_fig.add_subplot(111)
    i_at_npatt_plot = mi_ax.imshow(info_at_npatterns, interpolation='none', cmap='coolwarm', origin='lower')
    i_at_npatt_cbar = mi_fig.colorbar(i_at_npatt_plot, use_gridspec=True)
    i_at_npatt_cbar.set_label('MI at $|\mathcal{A}_{out}| = |\mathcal{A}_{in}|$ (bits)')
    mi_ax.set_xticks(np.arange(len(bias_range)))
    mi_ax.set_xticklabels([str(x) for x in bias_range])
    mi_ax.set_xlabel('Threshold current (pA)')
    mi_ax.set_yticks(np.arange(len(active_mf_fraction_range)))
    mi_ax.set_yticklabels([str(y) for y in active_mf_fraction_range])
    mi_ax.set_ylabel('MF activation probability')
    mi_ax.set_title('Effect of sparseness and inhibition\n$n_{dend}$=%d, $n_{MF}$=%d,  $|\mathcal{A}_{in}|$=%d, $N_{trials}$=%d' % (n_grc_dend, int(round(network_scale_range[0]*min_mf_number)), n_stim_patterns, n_trials))

if plot_kl_divergence:
    kl_fig = plt.figure()
    kl_ax = kl_fig.add_subplot(111)
    kl_plot = kl_ax.imshow(kl_div_values, interpolation='none', cmap='coolwarm', origin='lower')
    kl_cbar = kl_fig.colorbar(kl_plot, use_gridspec=True)
    kl_cbar.set_label('KL divergence (bits)')
    kl_ax.set_xticks(np.arange(len(bias_range)))
    kl_ax.set_xticklabels([str(x) for x in bias_range])
    kl_ax.set_xlabel('Threshold current (pA)')
    kl_ax.set_yticks(np.arange(len(active_mf_fraction_range)))
    kl_ax.set_yticklabels([str(y) for y in active_mf_fraction_range])
    kl_ax.set_ylabel('MF activation probability')
    kl_ax.set_title('Effect of sparseness and inhibition\n$n_{dend}$=%d, $n_{MF}$=%d,  $|\mathcal{A}_{in}|$=%d, $N_{trials}$=%d' % (n_grc_dend, int(round(network_scale_range[0]*min_mf_number)), n_stim_patterns, n_trials))

if plot_mi_detail:
    mi_det_ax.legend(loc='upper left')
    if plotting_mode == 'precision':
        mi_det_ax.set_xlabel('decoder precision (1/cluster separation)')
    elif plotting_mode == 'alphabet_size':
        mi_det_ax.set_xlabel('alphabet size (clusters in the decoder)')
    mi_det_ax.set_ylabel('MI (bits)')

if plot_dendrograms:
    dend_axes_ylim = (0, max([ax.get_ylim()[1] for ax in dend_axes]))
    for ax in dend_axes:
        ax.set_ylim(dend_axes_ylim)

plt.show()
    
