import glob
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
        self.extra_tonic_inhibition = extra_tonic_inhibition
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
        # the spike archive the point gets associated with can be an
        # archive for a larger set of simulations. For example, if
        # this point has n_stim_patterns=128, n_trials=50 and
        # sim_duration=150 but an archive is available for a point
        # which only differs from this one in having
        # n_stim_patterns=1024, n_trials=100 and sim_duration=200, the
        # larger archive will be reused to avoid rerunning the same
        # simulations.
        self.spike_archive_path = self.get_existing_spike_archive_path()
        if not self.spike_archive_path:
            self.spike_archive_path = "{0}/sp{1}_t{2}_sdur{3}.hdf5".format(self.data_folder_path,
                                                                           self.n_stim_patterns,
                                                                           self.n_trials,
                                                                           self.sim_duration)
        #--useful quantities
        self.network_graph = nx.read_graphml(self.graphml_network_filename,
                                             node_type=int)
        self.graph_mf_nodes = [n for n,d in self.network_graph.nodes(data=True) if d['bipartite']==0]
        self.graph_grc_nodes = [n for n,d in self.network_graph.nodes(data=True) if d['bipartite']==1]
        self.n_mf = len(self.graph_mf_nodes)
        self.n_grc = len(self.graph_grc_nodes)
        assert self.n_mf + self.n_grc == self.network_graph.number_of_nodes()
    def representation(self):
        # MUST NOT HAVE SPACES (see how simulations are submitted)
        return "SimpleParameterSpacePoint(%d,%d,%f,%f,%d,%d,%d,%d,%d,%d,%d,%d)" % (self.n_grc_dend, self.connectivity_rule, self.input_spatial_correlation_scale, self.active_mf_fraction, self.extra_tonic_inhibition, self.stim_rate_mu, self.stim_rate_sigma, self.noise_rate_mu, self.noise_rate_sigma, self.n_stim_patterns, self.n_trials, self.sim_duration)
    def representation_without_commas(self):
        # sanitised version of the Point representation, with commas
        # replaced by | signs. This is needed because of a known bug
        # in Legion's version of JSV which freaks out when script
        # arguments contain commas.
        return self.representation().replace(',', '+')
    def __repr__(self):
        return self.representation()
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
    def get_tar_simulation_archive_path(stimulus_pattern_index):
        # temporary tar archive for storing simulation data for a
        # pattern while the simulations are running. The
        # information contained in these archives will be
        # restructured later in the compression step.
        return "{0}/sp{1}_t{2}_sdur{3}_spn{4}.tar".format(self.data_folder_path,
                                                          self.n_stim_patterns,
                                                          self.n_trials,
                                                          self.sim_duration,
                                                          stimulus_pattern_index)
    def get_existing_spike_archive_path(self):
        """If it exists, return the path of a spike archive for this
        point. This means having all coordinates equal to the point
        except the sim duration which can be longer and the number of
        trials which can be larger. If no such archive exists, return
        None.

        """
        candidate_archive_paths = glob.glob('%s/sp*_t*_sdur*.hdf5' % (self.data_folder_path))
        compatible_archive_paths = sorted([path for path in candidate_archive_paths if float(path.rstrip('.hdf5').partition('sdur')[2])>=self.sim_duration and float(path.rpartition('_t')[2].partition('_sdur')[0])>=self.n_trials and float(path.rpartition('_t')[0].rpartition('sp')[2])>=self.n_stim_patterns])
        if compatible_archive_paths:
            return compatible_archive_paths[0]
        else:
            return None
