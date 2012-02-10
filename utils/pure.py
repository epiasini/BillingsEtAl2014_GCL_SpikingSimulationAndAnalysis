BASE_DIR = "/home/ucbtepi/code/network/data"

class ParamSpacePoint(object):
    def __init__(self,
                 sim_duration,
                 min_mf_number,
                 grc_mf_ratio,
                 n_grc_dend,
                 network_scale,
                 active_mf_fraction,
                 bias,
                 stim_rate_mu,
                 stim_rate_sigma,
                 noise_rate_mu,
                 noise_rate_sigma,
                 n_stim_patterns,
                 n_trials):
        #--parameter space coordinates
        self.sim_duration = sim_duration
        self.min_mf_number = int(round(min_mf_number))
        self.grc_mf_ratio = grc_mf_ratio
        self.n_grc_dend = int(round(n_grc_dend))
        self.network_scale = network_scale
        self.active_mf_fraction = active_mf_fraction
        self.bias = bias
        self.stim_rate_mu = stim_rate_mu
        self.stim_rate_sigma = stim_rate_sigma
        self.noise_rate_mu = noise_rate_mu
        self.noise_rate_sigma = noise_rate_sigma
        self.n_stim_patterns = int(round(n_stim_patterns))
        self.n_trials = int(round(n_trials))
        #--relevant filenames
        self.net_structure_folder_path = "%s/gmr%.02f/gd%d/s%.02f" % (BASE_DIR,
                                                                      self.grc_mf_ratio,
                                                                      self.n_grc_dend,
                                                                      self.network_scale)
        self.conn_pattern_filename = "%s/gmr%.02f_gd%d_s%.02f_conn.txt" % (self.net_structure_folder_path,
                                                                           self.grc_mf_ratio,
                                                                           self.n_grc_dend,
                                                                           self.network_scale)
        self.stim_pattern_filename = "%s/f%.02f/gmr%.02f_gd%d_s%.02f_sp%d_stim.txt" % (self.net_structure_folder_path,
                                                                                       self.active_mf_fraction,
                                                                                       self.grc_mf_ratio,
                                                                                       self.n_grc_dend,
                                                                                       self.network_scale,
                                                                                       self.n_stim_patterns)
        self.data_folder_path = "%s/f%.02f/b%02d/sm%d/ss%d/nm%d/ns%d" % (self.net_structure_folder_path,
                                                                         self.active_mf_fraction,
                                                                         self.bias,
                                                                         self.stim_rate_mu,
                                                                         self.stim_rate_sigma,
                                                                         self.noise_rate_mu,
                                                                         self.noise_rate_sigma)
    def __repr__(self):
        return "sim_dur: {0} | min_mf_number: {1} | gmr: {2} | gd: {3} | s: {4} | mf: {5} | b: {6} | sr_mu: {7} | sr_s: {8} | nr_mu: {9} | nr_s: {10} | np: {11} | t: {12}".format(self.sim_duration, self.min_mf_number, self.grc_mf_ratio, self.n_grc_dend, self.network_scale, self.active_mf_fraction, self.bias, self.stim_rate_mu, self.stim_rate_sigma, self.noise_rate_mu, self.noise_rate_sigma, self.n_stim_patterns, self.n_trials)
