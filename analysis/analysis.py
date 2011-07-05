import numpy as np
import h5py
import sys
import random
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist
import functools

sys.path.append("../")
from utils import data_archive_path_ctor

import pdb

def loadspikes(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, n_stim_patterns, n_trials, cell_type):
    
    archive_filename = data_archive_path_ctor(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, n_stim_patterns, n_trials)
    archive = h5py.File(archive_filename)
    
    #n_stim_patterns = len([x for x in archive.items() if isinstance(x[1], h5py.highlevel.Group)])
    #n_trials = len([x for x in archive['000'].items() if isinstance(x[1], h5py.highlevel.Group)])
    n_cells = archive['000']['00']['{0}_spiketimes'.format(cell_type)].shape[1]

    n_expcted_obs = n_stim_patterns * n_trials
    
    observation_list = [x[1]['{0}_spiketimes'.format(cell_type)] for s in archive.items() if isinstance(s[1], h5py.highlevel.Group) for x in s[1].items() if isinstance(x[1], h5py.highlevel.Group)]
    
    max_n_spikes = max([o.shape[0] for o in observation_list])

    spikes = -1 * np.ones(shape=(n_expcted_obs, n_cells, max_n_spikes))
    
    for k, o in enumerate(observation_list):
        spikes[k][:,0:o.shape[0]] = np.array(o).transpose()

    archive.close()
    return spikes

def convolve(obs_array, tau, dt):
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

TRAINING_RELATIVE_SIZE = .8

## min_mf_number = int(sys.argv[1])
## grc_mf_ratio = float(sys.argv[2])
## n_grc_dend = int(sys.argv[3])
## network_scale = float(sys.argv[4])
## active_mf_fraction = float(sys.argv[5])
## bias = float(sys.argv[6])
## n_stim_patterns = int(sys.argv[7])
## n_trials = int(sys.argv[8])

min_mf_number = 6
grc_mf_ratio = 3
n_grc_dend = 4
network_scale = 1.00
active_mf_fraction = 0.5
bias = 20
n_stim_patterns = 20
n_trials = 20

sim_length = 300.
tau = 5.
dt = 2.
multineuron_metric_mixing = 0
cell_type = 'grc'

n_obs = n_stim_patterns * n_trials


# load data
spikes = loadspikes(grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, n_stim_patterns, n_trials, cell_type)

# choose training and testing set
n_train_points = int(TRAINING_RELATIVE_SIZE * n_obs)
train_idxs = random.sample(xrange(n_obs), n_train_points)
test_idxs = [x for x in xrange(n_obs) if x not in train_idxs]

tr_spikes = spikes[train_idxs]
ts_spikes = spikes[test_idxs]

# convolve spike train to generate vector fields
tr_fields = convolve(tr_spikes, tau, dt)
ts_fields = convolve(ts_spikes, tau, dt)

# prepare multineuron distance function by partial application
fixed_c_multineuron_distance = functools.partial(multineuron_distance, c=multineuron_metric_mixing)

tr_distances = []
for h in range(n_train_points):
    for k in range(h+1, n_train_points):
        tr_distances.append(fixed_c_multineuron_distance(tr_fields[h], tr_fields[k]))

