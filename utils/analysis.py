import numpy as np
import h5py
import random
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
import functools
import itertools
import fcntl

import pyentropy as pe

from .paths import data_archive_path_ctor, mi_archive_path_ctor

def loadspikes(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, n_stim_patterns, n_trials, cell_type):
    
    archive_filename = data_archive_path_ctor(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, n_stim_patterns, n_trials)
    archive = h5py.File(archive_filename)
    
    n_cells = archive['000']['00']['{0}_spiketimes'.format(cell_type)].shape[1]

    n_expcted_obs = n_stim_patterns * n_trials
    
    observation_list = [x[1]['{0}_spiketimes'.format(cell_type)] for s in archive.items() if isinstance(s[1], h5py.highlevel.Group) for x in s[1].items() if isinstance(x[1], h5py.highlevel.Group)]
    
    max_n_spikes = max([o.shape[0] for o in observation_list])

    spikes = -1 * np.ones(shape=(n_expcted_obs, n_cells, max_n_spikes))
    
    for k, o in enumerate(observation_list):
        spikes[k][:,0:o.shape[0]] = np.array(o).transpose()

    archive.close()
    return spikes

def convolve(obs_array, sim_length, tau, dt):
    """Convolve with exponential kernel."""
    dt = float(dt)
    n_obs, n_cells, max_n_spikes = obs_array.shape
    kernel = np.exp(-np.arange(0, 10*tau, dt)/tau)
    kernel_length = 10*tau/dt
    n_bins = sim_length/dt
    conv = np.zeros(shape=(n_obs, n_cells, n_bins+kernel_length))

    for o, obs in enumerate(obs_array):
        for c, cell in enumerate(obs):
            for spike_index in [spike_time/dt for spike_time in cell[cell > 0]]:
                # here we could optimise this by writing it as a list comprehension, by using the .__add__ method
                conv[o, c, spike_index:spike_index+kernel_length] += kernel
    return conv

def multineuron_distance(p,q, theta):
    delta = p-q
    return np.sqrt(np.einsum('mt,mt', np.einsum('nt,nm->mt', delta, theta), delta))

def multineuron_distance_labeled_line(p,q):
    delta = p-q
    return np.sqrt(np.einsum('nt,nt', delta, delta))

def open_archive(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, n_stim_patterns, n_trials, multineuron_metric_mixing, training_size, linkage_method):
    filename = mi_archive_path_ctor(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias)
    archive_lock = open(filename, 'a')
    fcntl.lockf(archive_lock, fcntl.LOCK_EX)
    mi_archive = h5py.File(filename)
    nspg = mi_archive.require_group('sp%d' % n_stim_patterns)
    ntrg = nspg.require_group('t%d' % n_trials)
    mixg = ntrg.require_group('mix%.2f' % multineuron_metric_mixing)
    trsg = mixg.require_group('train%d' % training_size)
    target_group = trsg.require_group('method_%s' % linkage_method)
    return mi_archive, target_group, archive_lock

def analyse_single_configuration(min_mf_number, grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, n_stim_patterns, n_trials, sim_duration, tau, dt, multineuron_metric_mixing, training_size, linkage_method):
    # prepare hdf5 archive
    mi_archive, target_group, archive_lock = open_archive(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, n_stim_patterns, n_trials, multineuron_metric_mixing, training_size, linkage_method)
    tr_tree = None        
    if all([ds in target_group.keys() for ds in ['tr_indexes', 'tr_linkage', 'tr_direct_mi', 'ts_decoded_mi_plugin', 'ts_decoded_mi_qe', 'px_at_same_size_point']]):
        print('found previously computed results in hdf5 archive.')
        decoder_precision = (1./np.array(target_group['tr_linkage'])[:,2])[::-1]
        tr_direct_mi = np.array(target_group['tr_direct_mi'])
        ts_decoded_mi_plugin = np.array(target_group['ts_decoded_mi_plugin'])
        ts_decoded_mi_qe = np.array(target_group['ts_decoded_mi_qe'])
        tr_tree = np.array(target_group['tr_linkage'])
        px_at_same_size_point = np.array(target_group['px_at_same_size_point'])
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
        if multineuron_metric_mixing!=0:
            theta = (np.eye(n_cells) + multineuron_metric_mixing*(np.ones((n_cells, n_cells)) - np.eye(n_cells)))
            fixed_c_multineuron_distance = functools.partial(multineuron_distance, theta=theta)
        else:
            fixed_c_multineuron_distance = multineuron_distance_labeled_line
        tr_distances = []
        for h in range(n_tr_obs):
            for k in range(h+1, n_tr_obs):
                tr_distances.append(fixed_c_multineuron_distance(tr_fields[h], tr_fields[k]))

        # cluster training data
        print('clustering training data')
        tr_tree = linkage(tr_distances, method=linkage_method)

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
            if n_clusts == n_stim_patterns:
                px_at_same_size_point = s.PX
            

        # save analysis results in a separate file
        datasets_to_be_deleted = [ds for ds in ['tr_indexes', 'tr_linkage', 'tr_direct_mi', 'ts_decoded_mi_plugin', 'ts_decoded_mi_qe', 'px_at_same_size_point'] if ds in target_group.keys()]
        for ds in datasets_to_be_deleted:
            del target_group[ds]
        target_group.create_dataset('tr_indexes', data=np.array(train_idxs))
        target_group.create_dataset('tr_linkage', data=tr_tree)
        target_group.create_dataset('tr_direct_mi', data=tr_direct_mi)
        target_group.create_dataset('ts_decoded_mi_plugin', data=ts_decoded_mi_plugin)
        target_group.create_dataset('ts_decoded_mi_qe', data=ts_decoded_mi_qe)
        target_group.create_dataset('px_at_same_size_point', data=px_at_same_size_point)
    # close archive and return results
    mi_archive.close()
    fcntl.lockf(archive_lock, fcntl.LOCK_UN)
    return tr_direct_mi, ts_decoded_mi_plugin, ts_decoded_mi_qe, decoder_precision, tr_tree, px_at_same_size_point


def cluster_centroids(min_mf_number, grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, n_stim_patterns, n_trials, sim_duration, tau, dt, multineuron_metric_mixing, training_size, linkage_method, n_clusts, cell_type='grc'):
    '''Returns cluster centroids for the given analysis and number of clusters. Useful to build representations of the "typical" network activities at a particular resolution.'''
    mi_archive, target_group, archive_lock = open_archive(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, n_stim_patterns, n_trials, multineuron_metric_mixing, training_size, linkage_method)
    tr_spikes = loadspikes(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, n_stim_patterns, n_trials, cell_type)[target_group['tr_indexes']]
    tr_vectors = convolve(tr_spikes, sim_duration, tau, dt)
    flat_clustering = fcluster(target_group['tr_linkage'], n_clusts, criterion='maxclust')
    clust_idxs = sorted(list(set(flat_clustering)))
    centroids = np.array([np.mean(tr_vectors[flat_clustering==c], axis=0) for c in clust_idxs])
    mi_archive.close()
    fcntl.lockf(archive_lock, fcntl.LOCK_UN)
    return clust_idxs, centroids
    
def kl_divergence(p,q):
    return (p * np.log2(p/q)).sum()

def kl_divergence_from_flat_p(q):
    p = np.ones_like(q, dtype=np.float)/len(q)
    return kl_divergence(p,q)
