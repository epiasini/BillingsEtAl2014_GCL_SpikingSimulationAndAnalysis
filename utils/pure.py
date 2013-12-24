import networkx as nx

BASE_DIR = "/home/ucbtepi/code/network/data"

class SimpleParameterSpacePoint(object):
    """Used in the simulation script and as a base class for ParameterSpacePoint"""
    #--class constants
    BASE_DIR = BASE_DIR
    def __init__(self,
                 n_grc_dend,
                 connectivity_rule,
                 input_spatial_correlation_scale,
                 active_mf_fraction,
                 extra_tonic_inhibition,
                 stim_rate_mu,
                 stim_rate_sigma,
                 noise_rate_mu,
                 noise_rate_sigma,
                 n_stim_patterns,
                 n_trials,
                 sim_duration):
        #--parameter space coordinates
        self.n_grc_dend = int(round(n_grc_dend))
        self.connectivity_rule = int(round(connectivity_rule ))
        self.input_spatial_correlation_scale = input_spatial_correlation_scale
        self.active_mf_fraction = active_mf_fraction
        self.bias = bias
        self.stim_rate_mu = stim_rate_mu
        self.stim_rate_sigma = stim_rate_sigma
        self.noise_rate_mu = noise_rate_mu
        self.noise_rate_sigma = noise_rate_sigma
        self.n_stim_patterns = int(round(n_stim_patterns))
        self.SIZE_PER_SIMULATION = self.n_stim_patterns # TODO: rename this to something more sensible as "workers_per_parameter_space_point", or remove it altogether.
        self.n_trials = int(round(n_trials))
        self.sim_duration = sim_duration
        #--relevant filenames
        self.net_structure_folder_path = "%s/gd%d/cr%d" % (self.BASE_DIR,
                                                           self.n_grc_dend,
                                                           self.connectivity_rule)
        self.graphml_network_filename = "%s/gd%d_cr%d.graphml" % (self.net_structure_folder_path,
                                                                  self.n_grc_dend,
                                                                  self.connectivity_rule)
        self.stim_pattern_folder_path = "%s/iscs%.02f/f%.02f" % (self.net_structure_folder_path,
                                                                 self.input_spatial_correlation_scale,
                                                                 self.active_mf_fraction)
        self.stim_pattern_filename = "%s/sp%d_stim.txt" % (self.stim_pattern_folder_path,
                                                           self.n_stim_patterns)
        self.data_folder_path = "%s/b%04d/sm%d/ss%d/nm%d/ns%d" % (self.stim_pattern_folder_path,
                                                                  self.extra_tonic_inhibition,
                                                                  self.stim_rate_mu,
                                                                  self.stim_rate_sigma,
                                                                  self.noise_rate_mu,
                                                                  self.noise_rate_sigma)
        #--useful quantities
        self.network_graph = nx.read_graphml(self.graphml_network_filename,
                                             node_type=int)
        self.graph_mf_nodes = [n for n,d in self.network_graph.nodes(data=True) if d['bipartite']==0]
        self.graph_grc_nodes = [n for n,d in self.network_graph.nodes(data=True) if d['bipartite']==1]
        self.n_mf = len(self.graph_mf_nodes)
        self.n_grc = len(self.graph_grc_nodes)
        assert self.n_mf + self.n_grc == self.network_graph.number_of_nodes()
    def __repr__(self):
        # MUST NOT HAVE SPACES (see how simulations are submitted)
        return "SimpleParameterSpacePoint(%d,%f,%d,%d,%d,%d,%d,%d,%d,%d)" % (self.n_grc_dend, self.active_mf_fraction, self.extra_tonic_inibition, self.stim_rate_mu, self.stim_rate_sigma, self.noise_rate_mu, self.noise_rate_sigma, self.n_stim_patterns, self.n_trials, self.sim_duration)
    def __str__(self):
        return "gd: %d | mf: %.1f | b: %d | sr_mu: %d | sr_s: %d | nr_mu: %d | nr_s: %d | nsp: %d | t: %d | sdur: %d" % (self.n_grc_dend, self.active_mf_fraction, self.extra_tonic_inhibition, self.stim_rate_mu, self.stim_rate_sigma, self.noise_rate_mu, self.noise_rate_sigma, self.n_stim_patterns, self.n_trials, self.sim_duration)
    def get_simulation_reference(self, stimulus_pattern_index, trial):
        """
        Return the simulation reference (the name of the relative subdirectory) of a given element in a batch of simulations.
        """
        return "sp%d_t%d_spn%d_tn%d" % (self.n_stim_patterns, self.n_trials, stimulus_pattern_index, trial)
    def nC_cell_index_from_graph_node(self, node):
        """
        given the name of a graph node, return the nC CellGroup and the
        index of the cell in the group.
        """
        if node <= self.n_mf:
            return node - 1, 'MFs'
        else:
            return node - (self.n_mf + 1), 'GrCs'
