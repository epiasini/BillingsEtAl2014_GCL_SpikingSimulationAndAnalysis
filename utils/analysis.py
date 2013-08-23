import numpy as np
import h5py
import random
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist, squareform
import functools
from math import floor

def convolve(obs_array, sim_length, tau, dt):
    """Convolve with exponential kernel."""
    dt = float(dt)
    n_obs, n_cells, max_n_spikes = obs_array.shape
    kernel = np.exp(-np.arange(0, 10*tau, dt)/tau)
    kernel_length = len(kernel)#10*tau/dt
    n_bins = int(round(1.1*sim_length/dt))
    conv = np.zeros(shape=(n_obs, n_cells, n_bins))

    for o, obs in enumerate(obs_array):
        for c, cell in enumerate(obs):
            for spike_index in [int(floor(spike_time/dt)) for spike_time in cell]:
                # here we could optimise this by writing it as a list comprehension, by using the .__add__ method
                available_bins = n_bins - spike_index
                if available_bins > kernel_length:
                    conv[o, c, spike_index:spike_index+kernel_length] += kernel
                else:
                    conv[o, c, spike_index:] += kernel[:available_bins]
    return conv

def multineuron_distance(p,q, theta):
    delta = p-q
    return np.sqrt(np.einsum('mt,mt', np.einsum('nt,nm->mt', delta, theta), delta))

def multineuron_distance_labeled_line(p,q):
    delta = p-q
    return np.sqrt(np.einsum('nt,nt', delta, delta))

def output_level(spike_array):
    """if spike_array is a n_spiketrains*n_cells*max_spikes array, return the number of spikes per trial per cell."""
    o_level_array = (spike_array>0).sum(axis=2)
    o_level_hist_values, o_level_hist_edges = np.histogram(o_level_array, bins=10)
    return o_level_array, o_level_hist_values, o_level_hist_edges

def population_sparseness(level_array):
    """if level_array is a n_stimuli*n_cells array of spike
    numbers/firing rates, return the Treves-Rolls population
    sparseness measure. Any stimulus that doesn't evoke any response
    is excluded from the final averaging."""
    square_of_average_by_stimulus = np.square(np.mean(level_array, axis=1)) # n_responses * 1
    average_of_square_by_stimulus = np.mean(np.square(level_array), axis=1)
    # now we average across all stimuli to get the average population sparseness for the dataset. Note that if the sparseness is not well defined if, for at least one stimulus, the average number of spikes is 0 across all neurons.
    return 1 - np.mean((square_of_average_by_stimulus/average_of_square_by_stimulus)[average_of_square_by_stimulus > 0])

def cluster_centroids(min_mf_number, grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma, n_stim_patterns, n_trials, sim_duration, tau, dt, multineuron_metric_mixing, training_size, linkage_method, n_clusts, cell_type='grc'):
    '''Returns cluster centroids for the given analysis and number of clusters. Useful to build representations of the "typical" network activities at a particular resolution.'''
    mi_archive, target_group, archive_lock = open_mi_archive(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma, n_stim_patterns, n_trials, multineuron_metric_mixing, training_size, linkage_method)
    tr_spikes = loadspikes(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma, n_stim_patterns, n_trials, cell_type)[target_group['tr_indexes']]
    tr_vectors = convolve(tr_spikes, sim_duration, tau, dt)
    flat_clustering = fcluster(target_group['tr_linkage'], n_clusts, criterion='maxclust')
    clust_idxs = sorted(list(set(flat_clustering)))
    clust_sizes = [(flat_clustering==c).sum() for c in clust_idxs]
    centroids = np.array([np.mean(tr_vectors[flat_clustering==c], axis=0) for c in clust_idxs])
    close_archive(mi_archive, archive_lock)
    return clust_idxs, centroids, clust_sizes

def kl_divergence(p,q):
    return (p * np.log2(p/q)).sum()

def kl_divergence_from_flat_p(q):
    p = np.ones_like(q, dtype=np.float)/len(q)
    return kl_divergence(p,q)

def entropy(p):
    return -(p[p>0]*np.log2(p[p>0])).sum()

def output_sparsity(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma, n_stim_patterns, n_trials):
    out_spike_array = loadspikes(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma, n_stim_patterns, n_trials, cell_type='grc')
    grc_act_prob = np.mean([(ob[:,0]>-1).sum()/float(ob.shape[0]) for ob in out_spike_array])
    return grc_act_prob

def output_level_array(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma, n_stim_patterns, n_trials):
    out_spike_array = loadspikes(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, n_stim_patterns, n_trials, cell_type='grc')
    # TODO: need convolution?
    out_spike_array[out_spike_array==-1] = 0
    activity_level_out_array = np.array(out_spike_array, dtype=np.bool).sum(axis=2)
    return activity_level_out_array

def synchrony(p):
    """Population-wide average of Schreiber's (2003) bivariate correlation-based measure of spike-timing reliability. Meant to be used with already-filtered spike trains, so in this context they are probably going to be exponential- and not gaussian-filtered."""
    n = p.shape[0]
    sync = 0
    for i, cell1 in enumerate(p):
        for j, cell2 in enumerate(p[i:]):
            sync += np.nan_to_num(np.correlate(cell1, cell2)/np.sqrt((np.square(cell1).sum() * np.square(cell2).sum())))
    return (2/float(n*(n-1)))*sync

def mean_synchrony(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma, n_stim_patterns, n_trials, sim_duration, tau, dt):
    out_spikes = loadspikes(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma, n_stim_patterns, n_trials, cell_type='grc')
    out_spikes = out_spikes[random.sample(range(out_spikes.shape[0]), 100)]
    out_vectors = convolve(out_spikes, sim_duration, tau, dt)
    return np.mean([synchrony(p) for p in out_vectors])

def open_sync_archive(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma, n_stim_patterns, n_trials):
    filename = mi_archive_path_ctor(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma)
    archive_lock = open(filename, 'a')
    fcntl.lockf(archive_lock, fcntl.LOCK_EX)
    archive = h5py.File(filename)
    nspg = archive.require_group('sp%d' % n_stim_patterns)
    ntrg = nspg.require_group('t%d' % n_trials)
    target_group = ntrg
    return archive, target_group, archive_lock

def analyse_synchrony(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma, n_stim_patterns, n_trials, sim_duration, tau, dt):
    archive, target_group, archive_lock = open_sync_archive(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma, n_stim_patterns, n_trials)
    if 'mean_synchrony' in target_group.keys():
        print('Found synchrony value in %s' % archive.filename)
        synchrony = np.array(target_group['mean_synchrony'])
    else:
        print('Analysing synchrony from %s' % archive.filename)
        close_archive(archive, archive_lock)
        synchrony = mean_synchrony(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma, n_stim_patterns, n_trials, sim_duration, tau, dt)
        archive, target_group, archive_lock = open_sync_archive(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma, n_stim_patterns, n_trials)
        target_group.create_dataset('mean_synchrony', data=synchrony)
    close_archive(archive, archive_lock)
    return synchrony

def distance_matrix(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma, n_stim_patterns, n_trials, sim_duration, tau, dt):
    out_spikes = loadspikes(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma, n_stim_patterns, n_trials, cell_type='grc')
    out_vectors = convolve(out_spikes, sim_duration, tau, dt)
    distances = []
    for h, vec in enumerate(out_vectors):
        for k in range(h+1, out_vectors.shape[0]):
            distances.append(multineuron_distance_labeled_line(out_vectors[h], out_vectors[k]))
    distances = squareform(np.array(distances))
    return distances




