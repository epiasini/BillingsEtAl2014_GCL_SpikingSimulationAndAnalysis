#! /usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import h5py
import sys
import random
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
import functools
import itertools
import pyentropy as pe
from matplotlib import pyplot as plt

from utils.paths import data_archive_path_ctor, mi_archive_path_ctor
from utils.analysis import loadspikes, convolve, multineuron_distance

import pdb

class TrainingSetSizeError(Exception):
    pass

#+++++fixed parameters+++++++
sim_duration = 300.0 # hardcoded in simulate.py
n_stim_patterns = 20
min_mf_number = 6
grc_mf_ratio = 3.
tau = 5.
dt = 2.
plotting_mode = 'alphabet_size'
#+++++parameter ranges+++++++++++++
n_grc_dend_range = [4]
network_scale_range = [1.00]
active_mf_fraction_range = [0.5]
bias_range = [20.]
n_trials_range = [200]
training_size_range = [50]
multineuron_metric_mixing_range = [0., 0.05]
#++++++++++++++++++++++++++

#----parameter consistency control
if any([s >= min(n_trials_range) for s in training_size_range]):
    raise TrainingSetSizeError()

ranges = [n_grc_dend_range, network_scale_range, active_mf_fraction_range, bias_range, n_trials_range, training_size_range, multineuron_metric_mixing_range]
parameter_space = itertools.product(*ranges)


def analyse_single_configuration(min_mf_number, grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, n_stim_patterns, n_trials, sim_duration, tau, dt, multineuron_metric_mixing, training_size):
    # prepare hdf5 archive
    mi_archive = h5py.File(mi_archive_path_ctor(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias))
    nspg = mi_archive.require_group('sp%d' % n_stim_patterns)
    ntrg = nspg.require_group('t%d' % n_trials)
    mixg = ntrg.require_group('mix%.2f' % multineuron_metric_mixing)
    target_group = mixg.require_group('train%d' % training_size)
    if all([ds in target_group.keys() for ds in ['tr_indexes', 'tr_linkage', 'tr_direct_mi', 'ts_decoded_mi_plugin', 'ts_decoded_mi_qe']]):
        print('found previously computed results in hdf5 archive.')
        decoder_precision = (1./np.array(target_group['tr_linkage'])[:,2])[::-1]
        tr_direct_mi = np.array(target_group['tr_direct_mi'])
        ts_decoded_mi_plugin = np.array(target_group['ts_decoded_mi_plugin'])
        ts_decoded_mi_qe = np.array(target_group['ts_decoded_mi_qe'])
    else:
        print('no previous results found. computing from the simulation data.')
        cell_type = 'grc'
        n_obs = n_stim_patterns * n_trials

        # load data
        spikes = loadspikes(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, n_stim_patterns, n_trials, cell_type)
        n_cells = spikes.shape[1]

        # choose training and testing set: trials are picked at random, but every stim pattern is represented equally (i.e., get the same number of trials) in both sets. Trials are ordered with respect to their stim pattern.
        n_tr_obs_per_sp = training_size
        n_ts_obs_per_sp = n_trials - n_tr_obs_per_sp
        train_idxs = list(itertools.chain(*([x+n_trials*sp for x in random.sample(range(n_trials), n_tr_obs_per_sp)] for sp in range(n_stim_patterns))))
        test_idxs = [x for x in range(n_obs) if x not in train_idxs]
        n_tr_obs = len(train_idxs)
        n_ts_obs = len(test_idxs)
        tr_input_signals = [int(x/n_trials) for x in train_idxs]
        tr_spikes = spikes[train_idxs]
        ts_spikes = spikes[test_idxs]

        # convolve spike train to generate vector fields
        tr_fields = convolve(tr_spikes, sim_duration, tau, dt)
        ts_fields = convolve(ts_spikes, sim_duration, tau, dt)
        n_timepoints = tr_fields.shape[2]

        # prepare multineuron distance function by partial application and calculate distances
        print('computing distances between training observations')
        fixed_c_multineuron_distance = functools.partial(multineuron_distance, c=multineuron_metric_mixing)
        tr_distances = []
        for h in range(n_tr_obs):
            for k in range(h+1, n_tr_obs):
                tr_distances.append(fixed_c_multineuron_distance(tr_fields[h], tr_fields[k]))

        # cluster training data
        print('clustering training data')
        tr_tree = linkage(tr_distances, method='weighted')

        # create an array describing the precision of each clustering step
        decoder_precision = (1./tr_tree[:,2])[::-1]

        # compute mutual information by using direct clustering on training data
        # --note: fcluster doesn't work in the border case with n_clusts=n_obs, as it never returns the trivial clustering. Cluster number 0 is never present in a clustering.
        tr_direct_mi = np.zeros(n_tr_obs-1)
        Xn = 1 # the output is effectively one-dimensional
        Ym = n_stim_patterns
        Ny = np.array([n_tr_obs_per_sp for each in range(n_stim_patterns)])
        for n_clusts in range(1,n_tr_obs):
            Xm = n_clusts
            X = fcluster(tr_tree, t=n_clusts, criterion='maxclust') - 1 # cluster number 0 is never present in fcluster's output, so the elements of X live in [1,n_clusts] 
            X_dims = (Xn, Xm)
            s = pe.SortedDiscreteSystem(X, X_dims, Ym, Ny)
            s.calculate_entropies(method='plugin', sampling='naive', calc=['HX', 'HXY'])
            tr_direct_mi[n_clusts-1] = s.I()

        # train the decoder and use it to calculate mi on the testing dataset
        print("training the decoder and using it to calculate mi on test data")
        relevant_tr_clusts = np.zeros(n_tr_obs+n_tr_obs-1)
        # --first step: prepare the 'square part' of the distance matrix. special case for n_clusts==n_tr_obs.
        relevant_tr_clusts = range(n_tr_obs)
        tr_clustering = fcluster(tr_tree, t=n_tr_obs, criterion='maxclust')
        out_alphabet = np.zeros(shape=(n_tr_obs+n_tr_obs-1, n_cells, n_timepoints))
        out_alphabet[0:n_tr_obs] = tr_fields
        distances = np.zeros(shape=(n_ts_obs, n_tr_obs+n_tr_obs-1))
        for n, observation in enumerate(ts_fields):
            for m, symbol in enumerate(out_alphabet[relevant_tr_clusts]):
                distances[n,m] = fixed_c_multineuron_distance(observation, symbol)
        # --now iterate over the number of clusters and, step by step, train the decoder and use it to calculate mi
        ts_decoded_mi_plugin = np.zeros(n_tr_obs-1)
        ts_decoded_mi_qe = np.zeros(n_tr_obs-1)
        Ny = np.array([n_ts_obs_per_sp for each in range(n_stim_patterns)])
        for n_clusts in range(n_tr_obs-1,0,-1):
            clust_idx = n_tr_obs + n_tr_obs - n_clusts - 1 # n_tr_obs, n_tr_obs+1, ..., n_tr_obs+n_tr_obs-2
            joined = tr_tree[clust_idx-n_tr_obs, 0:2]
            [relevant_tr_clusts.remove(c) for c in joined] # from now on, ingore clusters that have just been merged..
            relevant_tr_clusts.append(clust_idx) # ..but include the newly formed cluster in the next computations
            tr_clustering[tr_clustering==joined[0]] = clust_idx # this is to avoid calling fcluster again
            tr_clustering[tr_clustering==joined[1]] = clust_idx
            # compute new symbol as the weighted average of the joined clusters' centroids (that is, the centroid of the new cluster)
            # prepare weights for weighted average
            if joined[0] <= n_tr_obs:
                left_weight = 1
            else:
                left_weight = tr_tree[joined[0]-n_tr_obs][3]
            if joined[1] <= n_tr_obs:
                right_weight = 1
            else:
                right_weight = tr_tree[joined[1]-n_tr_obs][3]
            out_alphabet[clust_idx] = (out_alphabet[joined[0]] *  left_weight+ out_alphabet[joined[1]] * right_weight)/(left_weight + right_weight)
            # fill in the column in the distance matrix corresponding to the newly created symbol
            for n, observation in enumerate(ts_fields):
                distances[n, clust_idx] = fixed_c_multineuron_distance(observation, out_alphabet[clust_idx])
            # decode test data with the updated alphabet
            decoded_output = distances[:, relevant_tr_clusts].argmin(axis=1)
            # compute MI on the decoded data
            Xm = n_clusts
            X_dims = (Xn, Xm)
            X = decoded_output
            s = pe.SortedDiscreteSystem(X, X_dims, Ym, Ny)
            s.calculate_entropies(method='plugin', sampling='naive', calc=['HX', 'HXY'])
            ts_decoded_mi_plugin[n_clusts-1] = s.I()
            s.calculate_entropies(method='qe', sampling='naive', calc=['HX', 'HXY'], qe_method='plugin')
            ts_decoded_mi_qe[n_clusts-1] = s.I()    

        # save analysis results in a separate file
        datasets_to_be_deleted = [ds for ds in ['tr_indexes', 'tr_linkage', 'tr_direct_mi', 'ts_decoded_mi_plugin', 'ts_decoded_mi_qe'] if ds in target_group.keys()]
        for ds in datasets_to_be_deleted:
            del target_group[ds]
        target_group.create_dataset('tr_indexes', data=np.array(train_idxs))
        target_group.create_dataset('tr_linkage', data=tr_tree)
        target_group.create_dataset('tr_direct_mi', data=tr_direct_mi)
        target_group.create_dataset('ts_decoded_mi_plugin', data=ts_decoded_mi_plugin)
        target_group.create_dataset('ts_decoded_mi_qe', data=ts_decoded_mi_qe)    
    # close archive and return results
    mi_archive.close()
    return tr_direct_mi, ts_decoded_mi_plugin, ts_decoded_mi_qe, decoder_precision


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
    # analyse data
    tr_direct_mi, ts_decoded_mi_plugin, ts_decoded_mi_qe, decoder_precision = analyse_single_configuration(min_mf_number, grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, n_stim_patterns, n_trials, sim_duration, tau, dt, multineuron_metric_mixing, training_size)
    # plot
    color = colors[k%len(colors)]
    if plotting_mode == 'precision':
        ax.plot(decoder_precision, tr_direct_mi, label='train%d, direct' % (training_size), linestyle='-.', color=color)
        #ax.plot(ts_decoded_mi_plugin, label='bias%d, trials%d, plugin' % (bias, n_trials))
        ax.plot(decoder_precision, ts_decoded_mi_qe, label='train%d, qe' % (training_size), color=color)
        ax.plot([decoder_precision[n_stim_patterns], decoder_precision[n_stim_patterns]], ax.get_ylim(), linestyle='--', color=color)
    elif plotting_mode == 'alphabet_size':
        ax.plot(tr_direct_mi, label='mixing%.2f, direct' % (multineuron_metric_mixing), linestyle='-.', color=color)
        ax.plot(ts_decoded_mi_qe, label='mixing%.2f, qe' % (multineuron_metric_mixing), color=color)
        ax.plot([n_stim_patterns, n_stim_patterns], ax.get_ylim(), linestyle='--', color=color)

#ax.plot([n_stim_patterns, n_stim_patterns], ax.get_ylim(), linestyle='--', color='black')
ax.legend(loc='upper right')
if plotting_mode == 'precision':
    ax.set_xlabel('decoder precision (1/cluster separation)')
elif plotting_mode == 'alphabet_size':
    ax.set_xlabel('alphabet size (clusters in the decoder)')
ax.set_ylabel('MI (bits)')
#if not prec_thresholds[0] in ax.get_xticks():
#    ax.set_xticks(list(ax.get_xticks()) + [prec_thresholds[0]])
#plt.show()
    
