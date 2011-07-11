import numpy as np
import h5py
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

def multineuron_distance(p, q, c=0):
    delta = p-q
    E = np.dot(delta, delta.transpose())
    weighted_distances = E * (np.eye(E.shape[0]) + c*(np.ones_like(E) - np.eye(E.shape[0])))
    d = np.sqrt(weighted_distances.sum())
    if np.isnan(d):
        d = 0
    return d
