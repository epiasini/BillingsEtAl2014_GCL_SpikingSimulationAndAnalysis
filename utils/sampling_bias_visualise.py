from matplotlib import pyplot as plt
import numpy as np
import pdb

from utils.analysis import analyse_single_configuration

#+++++fixed parameters+++++++
sim_duration = 300.0 # hardcoded in simulate.py
n_stim_patterns = 20
min_mf_number = 6
grc_mf_ratio = 3.
tau = 5.
dt = 2.
plotting_mode = 'alphabet_size'
n_grc_dend = 4
network_scale = 1.00
active_mf_fraction = 0.5
bias = 20.
training_size = 40
multineuron_metric_mixing = 0.
#+++++parameter ranges+++++++++++++
n_trials_range = [50, 100, 200]
#++++++++++++++++++++++++++

ideal_limit_n_trials = 500

fig = plt.figure()
ax = fig.add_subplot(111)

ideal_tr_direct_mi, ideal_ts_decoded_mi_plugin, ideal_ts_decoded_mi_qe, ideal_decoder_precision = analyse_single_configuration(min_mf_number, grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, n_stim_patterns, ideal_limit_n_trials, sim_duration, tau, dt, multineuron_metric_mixing, training_size)

sampling_bias = np.zeros(shape=(training_size*n_stim_patterns-1, len(n_trials_range)))

for k, n_trials in enumerate(n_trials_range):
    tr_direct_mi, ts_decoded_mi_plugin, ts_decoded_mi_qe, decoder_precision = analyse_single_configuration(min_mf_number, grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, n_stim_patterns, n_trials, sim_duration, tau, dt, multineuron_metric_mixing, training_size)
    sampling_bias[:, k] = ts_decoded_mi_plugin - ideal_ts_decoded_mi_qe

alphabet_size_values = [20,40,60,80,100]
for alphabet_size in alphabet_size_values:
    ax.plot(n_trials_range, sampling_bias[alphabet_size, :], label='%d symbols' % alphabet_size, marker='o')

ax.set_xlabel('total number of observations per input pattern')
ax.set_ylabel('MI bias (bits)')
ax.legend(loc='lower left')
