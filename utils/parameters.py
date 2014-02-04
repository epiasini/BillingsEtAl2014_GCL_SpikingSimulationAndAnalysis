import numpy as np
import itertools
import functools
import random
import os.path
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.cluster import KMeans
import pyentropy as pe

from pure import SimpleParameterSpacePoint
from archival import SpikesArchive, ResultsArchive
from analysis import convolve, multineuron_distance, multineuron_distance_labeled_line, hoyer_sparseness, activity_sparseness

class PSlice(object):
    """
    TO BE USED ONLY AS AN 'ASCENDING' SLICE
    Conceptually, a subclass of the builtin slice class. Defaults to a single-point slice when just one argument is given.
    """
    def __init__(self, start, stop=None, step=1):
        self.start = start
        if stop==None:
            self.stop = start + 1
        else:
            self.stop = stop
        self.step = step
        self.realstop = self.stop - self.step

class ParameterSpacePoint(SimpleParameterSpacePoint):
    def __init__(self,
                 n_grc_dend,
                 connectivity_rule,
                 input_spatial_correlation_scale,
                 active_mf_fraction,
                 gaba_scale,
                 dta,
                 exc_cond_scaling,
                 modulation_frequency,
                 stim_rate_mu,
                 stim_rate_sigma,
                 noise_rate_mu,
                 noise_rate_sigma,
                 n_stim_patterns,
                 n_trials,
                 sim_duration,
                 ana_duration,
                 training_size,
                 multineuron_metric_mixing,
                 linkage_method,
                 tau,
                 dt):
        #--analysis-specific coordinates
        self.ana_duration = ana_duration
        self.training_size = int(round(training_size))
        self.multineuron_metric_mixing = multineuron_metric_mixing
        self.linkage_method = int(round(linkage_method))
        self.tau = tau
        self.dt = dt
        self.linkage_method_string = ['ward', 'kmeans'][self.linkage_method]
        super(ParameterSpacePoint, self).__init__(n_grc_dend,
                                                  connectivity_rule,
                                                  input_spatial_correlation_scale,
                                                  active_mf_fraction,
                                                  gaba_scale,
                                                  dta,
                                                  exc_cond_scaling,
                                                  modulation_frequency,
                                                  stim_rate_mu,
                                                  stim_rate_sigma,
                                                  noise_rate_mu,
                                                  noise_rate_sigma,
                                                  n_stim_patterns,
                                                  n_trials,
                                                  sim_duration)
        #--useful quantities
        self.sim_transient_time = self.sim_duration - self.ana_duration
        #--archive objects
        self.spikes_arch = SpikesArchive(self)
        self.results_arch = ResultsArchive(self)
    def __repr__(self):
        simple = self.simple_representation()
        simple_args = simple.split('(')[1].split(')')[0]
        return('ParameterSpacePoint({0},{1},{2},{3},{4},{5},{6})'.format(simple_args, self.ana_duration, self.training_size, self.multineuron_metric_mixing, self.linkage_method, self.tau, self.dt))
    def __str__(self):
        analysis_specific_repr = " |@| adur: {0} | train: {1} | mix: {2} | link: {3} | tau: {4} | dt: {5}".format(self.ana_duration, self.training_size, self.multineuron_metric_mixing, self.linkage_method_string[self.linkage_method], self.tau, self.dt)
        return super(ParameterSpacePoint, self).__str__() + analysis_specific_repr
    def simple_representation(self):
        """Describe the point as if it were a SimpleParameterSpacePoint."""
        return super(ParameterSpacePoint, self).__repr__()
    def simple_representation_without_commas(self):
        return super(ParameterSpacePoint, self).representation_without_commas()
    def representation_without_commas(self):
        # sanitised version of the Point representation, with commas
        # replaced by | signs. This is needed because of a known bug
        # in Legion's version of JSV which freaks out when script
        # arguments contain commas.
        return self.__repr__().replace(',', '+')

    def is_spike_archive_compatible(self, path):
        """Check if the archive at the given path, if present, is suitable
        for providing the simulation data necessary to perform the
        analysis specified by this data point. Simulation duration and
        number of stimulus patterns need to be greater in the archive
        than in the analysis settings, while the number of trials that
        can be extracted from the archive can depend (via time
        slicing) on the length of the original simulations compared to
        the length of the analysis duration and the transient time (ie
        sim_duration-ana_duration) we are asking for.

        """
        path_sdur = float(path.rstrip('.hdf5').partition('sdur')[2])
        path_n_trials = float(path.rpartition('_t')[2].partition('_sdur')[0]) * (1 + max(0, (path_sdur - self.sim_duration)//(self.ana_duration + self.SIM_DECORRELATION_TIME)))
        path_spn = float(path.rpartition('_t')[0].rpartition('sp')[2])
        sdur_c = path_sdur >= self.sim_duration
        n_trials_c = path_n_trials >= self.n_trials
        spn_c = path_spn >= self.n_stim_patterns
        return all([sdur_c, n_trials_c, spn_c])
    #-------------------
    # Simulation methods
    #-------------------
    
    #-------------------
    # Compression methods
    #-------------------
    def run_compression(self):
        pass
    #-------------------
    # Analysis methods
    #-------------------
    def run_analysis(self):
        if self.results_arch.load():
            # we have the results already (loaded in memory or on the disk)
            pass
        else:
            # check if the spikes archive to analyse is actually present on disk
            if not os.path.isfile(self.spike_archive_path):
                raise Exception("Spike archive {} not found! aborting analysis.".format(self.spike_archive_path))
            # we actually need to calculate them
            print("Analysing for: {0} from spike archive: {1}".format(self, self.spike_archive_path))
            n_obs = self.n_stim_patterns * self.n_trials
            # load data
            min_clusts_analysed = int(round(self.n_stim_patterns * 1.0))
            max_clusts_analysed = int(round(self.n_stim_patterns * 1.0))
            clusts_step = max(int(round(self.n_stim_patterns * 0.05)), 1)
            # choose training and testing set: trials are picked at random, but every stim pattern is represented equally (i.e., get the same number of trials) in both sets. Trials are ordered with respect to their stim pattern.
            n_tr_obs_per_sp = self.training_size
            n_ts_obs_per_sp = self.n_trials - n_tr_obs_per_sp
            train_idxs = list(itertools.chain(*([x+self.n_trials*sp for x in random.sample(range(self.n_trials), n_tr_obs_per_sp)] for sp in range(self.n_stim_patterns))))
            test_idxs = [x for x in range(n_obs) if x not in train_idxs]
            n_tr_obs = len(train_idxs)
            n_ts_obs = len(test_idxs)
            Ym = self.n_stim_patterns
            Ny = np.array([n_ts_obs_per_sp for each in range(self.n_stim_patterns)])
            Xn = 1 # the output is effectively one-dimensional
            # initialize data structures for storage of results
            ts_decoded_mi_plugin = np.zeros(n_obs)
            ts_decoded_mi_qe = np.zeros(n_obs)
            ts_decoded_mi_pt = np.zeros(n_obs)
            ts_decoded_mi_nsb = np.zeros(n_obs)

            # compute mutual information by using direct clustering on training data (REMOVED)
            # --note: fcluster doesn't work in the border case with n_clusts=n_obs, as it never returns the trivial clustering. Cluster number 0 is never present in a clustering.
            print('counting spikes in output spike trains')
            i_level_array = self.spikes_arch.get_spike_counts(cell_type='mf')
            o_level_array = self.spikes_arch.get_spike_counts(cell_type='grc')
            print('computing input and output sparsity')
            i_sparseness_hoyer = hoyer_sparseness(i_level_array)
            i_sparseness_activity = activity_sparseness(i_level_array)
            o_sparseness_hoyer = hoyer_sparseness(o_level_array)
            o_sparseness_activity = activity_sparseness(o_level_array)
            print('input sparseness: hoyer {:.2f}, activity {:.2f}'.format(i_sparseness_hoyer, i_sparseness_activity))
            print('output sparseness: hoyer {:.2f}, activity {:.2f}'.format(o_sparseness_hoyer, o_sparseness_activity))
            if self.linkage_method_string == 'kmeans':
                spike_counts = o_level_array
                # divide spike count data in training and testing set
                tr_spike_counts = np.array([spike_counts[o] for o in train_idxs])
                ts_spike_counts = np.array([spike_counts[o] for o in test_idxs])
                for n_clusts in range(min_clusts_analysed, max_clusts_analysed+1, clusts_step):
                    clustering = KMeans(n_clusters=n_clusts)
                    print('performing k-means clustering on training set (training the decoder) for k='+str(n_clusts))
                    clustering.fit(tr_spike_counts)
                    print('using the decoder trained with k-means clustering to classify data points in testing set')
                    decoded_output = clustering.predict(ts_spike_counts)
                    # calculate MI
                    print('calculating MI')
                    Xm = n_clusts
                    X_dims = (Xn, Xm)
                    X = decoded_output
                    s = pe.SortedDiscreteSystem(X, X_dims, Ym, Ny)
                    s.calculate_entropies(method='plugin', calc=['HX', 'HXY'])
                    ts_decoded_mi_plugin[n_clusts-1] = s.I()
                    s.calculate_entropies(method='qe', sampling='naive', calc=['HX', 'HXY'], qe_method='plugin')
                    ts_decoded_mi_qe[n_clusts-1] = s.I()
                    s.calculate_entropies(method='pt', sampling='naive', calc=['HX', 'HXY'])
                    ts_decoded_mi_pt[n_clusts-1] = s.I()
                    s.calculate_entropies(method='nsb', sampling='naive', calc=['HX', 'HXY'])
                    ts_decoded_mi_nsb[n_clusts-1] = s.I()            
            else:
                tr_tree = np.zeros(shape=(n_tr_obs-1, 3))
                import pymuvr
                spikes = self.spikes_arch.get_spikes(cell_type='grc')
                self.spikes_arch.load_attrs()
                tr_spikes = [spikes[o] for o in train_idxs]
                ts_spikes = [spikes[o] for o in test_idxs]

                # compute multineuron distance between each pair of training observations
                print('calculating distances between training observations')
                tr_distances = pymuvr.square_distance_matrix(tr_spikes,
                                                             self.multineuron_metric_mixing,
                                                             self.tau)
                # cluster training data
                print('clustering training data')
                tr_tree = linkage(tr_distances, method=self.linkage_method_string)

                # train the decoder and use it to calculate mi on the testing dataset
                print("training the decoder and using it to calculate mi on test data")

                tr_distances_square = np.square(tr_distances)

                for n_clusts in range(min_clusts_analysed, max_clusts_analysed+1):
                    # iterate over the number of clusters and, step by
                    # step, train the decoder and use it to calculate mi
                    tr_clustering = fcluster(tr_tree, t=n_clusts, criterion='maxclust')
                    out_alphabet = []
                    for c in range(1,n_clusts+1):
                        # every cluster is represented in the output
                        # alphabet by the element which minimizes the sum
                        # of intra-cluster square distances
                        obs_in_c = [ob for ob in range(n_tr_obs) if tr_clustering[ob]==c]
                        sum_of_intracluster_square_distances = tr_distances_square[obs_in_c,:][:,obs_in_c].sum(axis=1)
                        out_alphabet.append(tr_spikes[np.argmin(sum_of_intracluster_square_distances)])
                    distances = pymuvr.distance_matrix(ts_spikes,
                                                       out_alphabet,
                                                       self.multineuron_metric_mixing,
                                                       self.tau)
                    # each observation in the testing set is decoded by
                    # assigning it to the cluster whose representative
                    # element it's closest to
                    decoded_output = distances.argmin(axis=1)
                    # calculate MI
                    Xm = n_clusts
                    X_dims = (Xn, Xm)
                    X = decoded_output
                    s = pe.SortedDiscreteSystem(X, X_dims, Ym, Ny)
                    s.calculate_entropies(method='qe', sampling='naive', calc=['HX', 'HXY'], qe_method='plugin')
                    ts_decoded_mi_qe[n_clusts-1] = s.I()
                    s.calculate_entropies(method='pt', sampling='naive', calc=['HX', 'HXY'])
                    ts_decoded_mi_pt[n_clusts-1] = s.I()
                    s.calculate_entropies(method='nsb', sampling='naive', calc=['HX', 'HXY'])
                    ts_decoded_mi_nsb[n_clusts-1] = s.I()            
                    if n_clusts == self.n_stim_patterns:
                        px_at_same_size_point = s.PX 
                # save linkage tree to results archive (only if
                # performing hierarchical clustering)
                self.results_arch.update_result('tr_linkage', data=tr_tree)

            # save analysis results in the archive
            print('updating results archive')
            self.results_arch.update_result('ts_decoded_mi_plugin', data=ts_decoded_mi_plugin)
            self.results_arch.update_result('ts_decoded_mi_qe', data=ts_decoded_mi_qe)
            self.results_arch.update_result('ts_decoded_mi_pt', data=ts_decoded_mi_pt)
            self.results_arch.update_result('ts_decoded_mi_nsb', data=ts_decoded_mi_nsb)

            self.results_arch.update_result('i_sparseness_hoyer', data=i_sparseness_hoyer)
            self.results_arch.update_result('i_sparseness_activity', data=i_sparseness_activity)
            self.results_arch.update_result('o_sparseness_hoyer', data=o_sparseness_hoyer)
            self.results_arch.update_result('o_sparseness_activity', data=o_sparseness_activity)
            # update attributes
            self.results_arch.load()


# A numpy ndarray with object dtype, and composed of (ParameterSpacePoint)s.
ParameterSpaceMesh = np.vectorize(ParameterSpacePoint)

class ParameterSpace(np.ndarray):
    ABSOLUTE_DIDXS = [
        'n_grc_dend',
        'connectivity_rule',
        'input_spatial_correlation_scale',
        'active_mf_fraction',
        'gaba_scale',
        'dta',
        'exc_cond_scaling',
        'modulation_frequency',
        'stim_rate_mu',
        'stim_rate_sigma',
        'noise_rate_mu',
        'noise_rate_sigma',
        'n_stim_patterns',
        'n_trials',
        'sim_duration',
        'ana_duration',
        'training_size',
        'multineuron_metric_mixing',
        'linkage_method',
        'tau',
        'dt']
    def __new__(cls,
                n_grc_dend_slice,
                connectivity_rule_slice,
                input_spatial_correlation_scale_slice,
                active_mf_fraction_slice,
                gaba_scale_slice,
                dta_slice,
                exc_cond_scaling_slice,
                modulation_frequency_slice,
                stim_rate_mu_slice,
                stim_rate_sigma_slice,
                noise_rate_mu_slice,
                noise_rate_sigma_slice,
                n_stim_patterns_slice,
                n_trials_slice,
                sim_duration_slice,
                ana_duration_slice,
                trainig_size_slice,
                multineuron_metric_mixing_slice,
                linkage_method_slice,
                tau_slice,
                dt_slice):
        # Create the ndarray instance of our type, given the usual
        # ndarray input arguments.  This will call the standard
        # ndarray constructor, but return an object of our type.
        # It also triggers a call to ParamSpace.__array_finalize__
        m = np.mgrid[n_grc_dend_slice,
                     connectivity_rule_slice,
                     input_spatial_correlation_scale_slice,
                     active_mf_fraction_slice,
                     gaba_scale_slice,
                     dta_slice,
                     exc_cond_scaling_slice,
                     modulation_frequency_slice,
                     stim_rate_mu_slice,
                     stim_rate_sigma_slice,
                     noise_rate_mu_slice,
                     noise_rate_sigma_slice,
                     n_stim_patterns_slice,
                     n_trials_slice,
                     sim_duration_slice,
                     ana_duration_slice,
                     trainig_size_slice,
                     multineuron_metric_mixing_slice,
                     linkage_method_slice,
                     tau_slice,
                     dt_slice]
        obj = ParameterSpaceMesh(*m).view(cls)
        # Dimensional InDeXS. Dictionary or list? I like the 'semantic' nature of a dictionary, but lists have a natural ordering and a natural way of updating when a dimension is added or removed.
        obj.DIDXS = list(ParameterSpace.ABSOLUTE_DIDXS)
        obj.ABBREVIATIONS = {
            'n_grc_dend': 'gd',
            'connectivity_rule': 'cr',
            'input_spatial_correlation_scale': 'iscs',
            'active_mf_fraction': 'mf',
            'gaba_scale': 'b',
            'dta': 'dta',
            'exc_cond_scaling': 'ecs',
            'modulation_frequency': 'mod',
            'stim_rate_mu': 'sm',
            'stim_rate_sigma': 'ss',
            'noise_rate_mu': 'nm',
            'noise_rate_sigma': 'ns',
            'n_stim_patterns': 'sp',
            'n_trials': 't',
            'sim_duration_slice': 'sdur',
            'ana_duration_slice': 'adur',
            'training_size': 'tr'}
        # Finally, we must return the newly created object:
        return obj
    def __array_finalize__(self, obj):
        # ``self`` is a new object resulting from
        # ndarray.__new__(Param_space, ...), therefore it only has
        # attributes that the ndarray.__new__ constructor gave it -
        # i.e. those of a standard ndarray.
        #
        # We could have got to the ndarray.__new__ call in 3 ways:
        # From an explicit constructor - e.g. ParamSpace():
        #    obj is None
        #    (we're in the middle of the ParamSpace.__new__
        #    constructor, and self.DIDXS will be set when we return to
        #    ParamSpace.__new__)
        if obj is None: return
        # From view casting - e.g arr.view(ParamSpace):
        #    obj is arr
        #    (type(obj) can be ParamSpace)
        # From new-from-template - e.g paramspace[:3]
        #    type(obj) is ParamSpace
        #
        # Note that it is here, rather than in the __new__ method,
        # that we set the default value for DIDXS, because this
        # method sees all creation of default objects - with the
        # ParamSpace.__new__ constructor, but also with
        # arr.view(ParamSpace).
        self.DIDXS = list(getattr(obj, 'DIDXS', []))
        self.ABBREVIATIONS = dict(getattr(obj, 'ABBREVIATIONS', {}))
        self.fixed_parameters = dict(getattr(obj, 'fixed_parameters', {}))
        self.fixed_parameters_short = dict(getattr(obj, 'fixed_parameters_short', {}))
        # We do not need to return anything
    #-----------------------------------------------
    # Basic space and coordinate manipulation methods
    #-----------------------------------------------
    def _didx(self, parameter):
        return self.DIDXS.index(parameter)
    def _param(self, didx):
        return self.DIDXS[didx]
    def absolute_didx(self, parameter):
        return ParameterSpace.ABSOLUTE_DIDXS.index(parameter)
    def get_sorted_fixed_params(self):
        sorted_fixed_params = sorted(self.fixed_parameters, key=self.absolute_didx)
        return 
    def get_range(self, parameter):
        # TODO: improving by vectorising getattr
        return np.unique(np.array([getattr(x, parameter, None) for x in self.flat]))
    def _get_attribute_array(self, attribute):
        """
        Return an array containing the value of the requested attribute
        for all the points in self.  The shape of the returned array
        is of course self.shape. Default to np.nan when a point
        doesn't have the requested attribute: this is to increase
        robustness against missing datapoints (typical situation:
        trying to plot a 2d heatmap where we lack one of the values).

        """
        return np.array([getattr(x, attribute, np.nan) for x in self.flat]).reshape(self.shape)
    def _get_idx_from_value(self, parameter, value):
        param_range = self.get_range(parameter)
        if value not in param_range:
            # Trying to find the index for a parameter (i.e.,
            # coordinate) value that's not on the mesh raises an
            # exception.
            raise ValueError('Parameter value ({0}, {1}) not present on the mesh!'.format(value, parameter))
        return np.searchsorted(param_range, value)
    def _fix_dimension(self, parameter, value):
        self.fixed_parameters[parameter] = value
        del self.DIDXS[self._didx(parameter)]
        try:
            self.fixed_parameters_short[self.ABBREVIATIONS[parameter]] = value
        except KeyError:
            pass
    def _get_embedded_hyperplane(self, parameter, value):
        """
        Return a new Parameter_space object made only of those points that
        satisfy the condition parameter==value.  Don't change the
        dimensionality of the space (keep the hyperplane 'embedded' in
        the higher-dimensional space).

        """
        idx = self._get_idx_from_value(parameter, value)
        return np.split(self, self.shape[self._didx(parameter)], self._didx(parameter))[idx]
    def _get_hyperplane(self, parameter, value):
        hp = self._get_embedded_hyperplane(parameter, value)
        hp = hp.squeeze(hp._didx(parameter))
        hp._fix_dimension(parameter, value)
        return hp
    def get_subspace(self, *parameter_value_pairs):
        """
        Return a new Parameter_space object where one or more parameters
        (dimensions) have been fixed by setting parameter==value for
        every parameter, value pair in the parameter_value_pairs
        iterable.  Reduce the dimensionality of the space, forgetting
        about the dimensions that have been fixed.

        """
        temp = self._get_hyperplane(*parameter_value_pairs[0])
        for pair in parameter_value_pairs[1:]:
            temp = temp._get_hyperplane(*pair)
        return temp
    def get_nontrivial_subspace(self, *parameter_value_pairs):
        """
        Return a new Parameter_space object where one or more parameters
        (dimensions) have been fixed by setting parameter==value for
        every (parameter, value) pair that is given as an argument.
        Reduce the dimensionality of the space, forgetting about the
        dimensions that have been fixed.  Finally, reduce further the
        dimensionality by fogetting about the _trivial_ dimensions,
        i.e. the ones where our points don't have any variability.

        """
        temp = self.get_subspace(*parameter_value_pairs)
        for parameter in [temp._param(k) for k,s in enumerate(temp.shape) if s==1]:
            value = temp.get_range(parameter)[0]
            temp = temp.squeeze(temp._didx(parameter))
            temp._fix_dimension(parameter, value)
        return temp
    #-------------------
    # Simulation methods
    #-------------------

    #-------------------
    # Compression methods
    #-------------------

    #-------------------
    # Analysis methods
    #-------------------
    def run_analysis(self):
        for p in self.flat:
            p.run_analysis()
    def load_analysis_results(self):
        return [p.results_arch.load() for p in self.flat]
    #-------------------
    # Visualisation methods
    #-------------------
    def plot_2d_heatmap(self, heat_dim, fig_title=''):
        """Plot a bidimensional heatmap for the heat_dim quantity"""        
        if len(self.shape) > 2:
            raise Exception()
        x_param = self._param(1)
        y_param = self._param(0)        
        fig, ax = plt.subplots()
        plot = ax.imshow(self._get_attribute_array(heat_dim), interpolation='none', cmap='coolwarm', origin='lower')
        cbar = fig.colorbar(plot, use_gridspec=True)
        cbar.set_label(heat_dim)
        ax.set_xticks(np.arange(self.shape[1]))
        ax.set_xticklabels([str(x) for x in self.get_range(x_param)])
        ax.set_xlabel(x_param)
        ax.set_yticks(np.arange(self.shape[0]))
        ax.set_yticklabels([str(y) for y in self.get_range(y_param)])
        ax.set_ylabel(y_param)
        ax.set_title(fig_title)
        return fig, ax

