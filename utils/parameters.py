import numpy as np

class Param_space_point(object):
    def __init__(self, sim_duration, min_mf_number, grc_mf_ratio, n_grc_dend, network_scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma, n_stim_patterns, n_trials):
        self.sim_duration = sim_duration
        self.min_mf_number = min_mf_number
        self.grc_mf_ratio = grc_mf_ratio
        self.n_grc_dend = n_grc_dend
        self.network_scale = network_scale
        self.active_mf_fraction = active_mf_fraction
        self.bias = bias
        self.stim_rate_mu = stim_rate_mu
        self.stim_rate_sigma = stim_rate_sigma
        self.noise_rate_mu = noise_rate_mu
        self.noise_rate_sigma = noise_rate_sigma
        self.n_stim_patterns = n_stim_patterns
        self.n_trials = n_trials
    def __repr__(self):
        return "sim_dur: {0} | min_mf_number: {1} | gmr: {2} | gd: {3} | s: {4} | mf: {5} | b: {6} | sr_mu: {7} | sr_s: {8} | nr_mu: {9} | nr_s: {10} | np: {11} | t: {12}".format(self.sim_duration, self.min_mf_number, self.grc_mf_ratio, self.n_grc_dend, self.network_scale, self.active_mf_fraction, self.bias, self.stim_rate_mu, self.stim_rate_sigma, self.noise_rate_mu, self.noise_rate_sigma, self.n_stim_patterns, self.n_trials)

Param_space_mesh = np.vectorize(Param_space_point)

class Param_space(np.ndarray):
    def __new__(cls, sim_duration_slice, min_mf_number_slice, grc_mf_ratio_slice, n_grc_dend_slice, network_scale_slice, active_mf_fraction_slice, bias_slice, stim_rate_mu_slice, stim_rate_sigma_slice, noise_rate_mu_slice, noise_rate_sigma_slice, n_stim_patterns_slice, n_trials_slice):
        m = np.mgrid[sim_duration_slice, min_mf_number_slice, grc_mf_ratio_slice, n_grc_dend_slice, network_scale_slice, active_mf_fraction_slice, bias_slice, stim_rate_mu_slice, stim_rate_sigma_slice, noise_rate_mu_slice, noise_rate_sigma_slice, n_stim_patterns_slice, n_trials_slice]
        obj = Param_space_mesh(*m).view(cls)
        return obj
    def __array_finalize__(self, obj):
        self.DIDXS = {
            'sim_duration': 0,
            'min_mf_number': 1,
            'grc_mf_ratio': 2,
            'n_grc_dend': 3,
            'network_scale': 4,
            'active_mf_fraction': 5,
            'bias': 6,
            'stim_rate_mu': 7,
            'stim_rate_sigma': 8,
            'noise_rate_mu': 9,
            'noise_rate_sigma': 10,
            'n_stim_patterns': 11,
            'n_trials': 12
            }
    def get_range(self, parameter):
        # TODO: improving by vectorising getattr
        return np.unique(np.array([getattr(x, parameter, None) for x in self.flatten()]))
    def get_idx(self, parameter, value):
        return np.searchsorted(self.get_range(parameter), value)
    def embedded_subspace(self, parameter, value):
        idx = self.get_idx(parameter, value)
        return np.split(self, self.shape[self.DIDXS[parameter]], self.DIDXS[parameter])[idx]
    def subspace(self, parameter, value):
        # WARNING: not consistent yet (doesn't update the DIDXS)
        return np.squeeze(self.embedded_subspace(parameter, value))

