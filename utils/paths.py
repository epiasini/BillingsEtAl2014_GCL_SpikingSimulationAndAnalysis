BASE_DIR = "/home/ucbtepi/code/network/data"

def ref_ctor(n_stim_patterns, n_trials, stimulus_pattern_index, trial):
    """Returns the simulation reference (the name of the relative subdirectory) of a given element in a batch of simulations."""
    sim_ref = "sp%d_t%d_spn%d_tn%d" % (n_stim_patterns, n_trials, stimulus_pattern_index, trial)
    return sim_ref

def conn_pattern_filename(grc_mf_ratio, n_grc_dend, scale):
    return "%s/gmr%.02f_gd%d_s%.02f_conn.txt" % (net_structure_folder_path_ctor(grc_mf_ratio, n_grc_dend, scale), grc_mf_ratio, n_grc_dend, scale)

def stim_pattern_filename(grc_mf_ratio, n_grc_dend, scale, active_mf_fraction, n_stim_patterns):
    return "%s/f%.02f/gmr%.02f_gd%d_s%.02f_sp%d_stim.txt" % (net_structure_folder_path_ctor(grc_mf_ratio, n_grc_dend, scale), active_mf_fraction, grc_mf_ratio, n_grc_dend, scale, n_stim_patterns)    

def net_structure_folder_path_ctor(grc_mf_ratio, n_grc_dend, scale):
    return "%s/gmr%.02f/gd%d/s%.02f" % (BASE_DIR, grc_mf_ratio, n_grc_dend, scale)

def data_folder_path_ctor(grc_mf_ratio, n_grc_dend, scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma):
    return "%s/f%.02f/b%02d/sm%d/ss%d/nm%d/ns%d" % (net_structure_folder_path_ctor(grc_mf_ratio, n_grc_dend, scale), active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma)

def datafiles_base_name(grc_mf_ratio, n_grc_dend, scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma, n_stim_patterns, n_trials):
    return "%s/sp%d_t%d" % (data_folder_path_ctor(grc_mf_ratio, n_grc_dend, scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma), n_stim_patterns, n_trials)

def data_archive_path_ctor(grc_mf_ratio, n_grc_dend, scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma, n_stim_patterns, n_trials):
    return "%s.hdf5" % datafiles_base_name(grc_mf_ratio, n_grc_dend, scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma, n_stim_patterns, n_trials)

def mi_archive_path_ctor(grc_mf_ratio, n_grc_dend, scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma):
    return "%s/mi.hdf5" % data_folder_path_ctor(grc_mf_ratio, n_grc_dend, scale, active_mf_fraction, bias, stim_rate_mu, stim_rate_sigma, noise_rate_mu, noise_rate_sigma)
