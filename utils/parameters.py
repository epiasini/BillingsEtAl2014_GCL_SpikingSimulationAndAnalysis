import numpy as np
import itertools
import functools
import random
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.cluster import MiniBatchKMeans
import pyentropy as pe

import pymuvr
from pure import SimpleParameterSpacePoint
from archival import SpikesArchive, ResultsArchive
from analysis import convolve, multineuron_distance, multineuron_distance_labeled_line, output_level, synchrony

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
                 n_trials,
                 training_size,
                 multineuron_metric_mixing,
                 linkage_method,
                 tau,
                 dt):
        super(ParameterSpacePoint, self).__init__(sim_duration,
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
                                                  n_trials)        
        #--analysis-specific coordinates
        self.training_size = int(round(training_size))
        self.multineuron_metric_mixing = multineuron_metric_mixing
        self.linkage_method = int(round(linkage_method))
        self.tau = tau
        self.dt = dt
        self.linkage_method_string = ['ward', 'kmeans'][self.linkage_method]
        #--archive objects
        self.spikes_arch = SpikesArchive(self)
        self.results_arch = ResultsArchive(self)
    def __repr__(self):
        simple = self.simple_representation()
        simple_args = simple.split('(')[1].split(')')[0]
        return('ParameterSpacePoint({0},{1},{2},{3},{4},{5})'.format(simple_args, self.training_size, self.multineuron_metric_mixing, self.linkage_method, self.tau, self.dt))
    def __str__(self):
        analysis_specific_repr = " |@| train: {0} | mix: {1} | link: {2} | tau: {3} | dt: {4}".format(self.training_size, self.multineuron_metric_mixing, self.linkage_method_string[self.linkage_method], self.tau, self.dt)
        return super(ParameterSpacePoint, self).__str__() + analysis_specific_repr
    def simple_representation(self):
        """Describe the point as if it were a SimpleParameterSpacePoint."""
        return super(ParameterSpacePoint, self).__repr__()
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
            # we actually need to calculate them
            print("Analysing for: {0}".format(self))
            n_obs = self.n_stim_patterns * self.n_trials
            # load data
            min_clusts_analysed = int(round(self.n_stim_patterns * 0.4))
            max_clusts_analysed = int(round(self.n_stim_patterns * 1.6))
            clusts_step = int(round(self.n_stim_patterns * 0.05))
            # choose training and testing set: trials are picked at random, but every stim pattern is represented equally (i.e., get the same number of trials) in both sets. Trials are ordered with respect to their stim pattern.
            n_tr_obs_per_sp = self.training_size
            n_ts_obs_per_sp = self.n_trials - n_tr_obs_per_sp
            train_idxs = list(itertools.chain(*([x+self.n_trials*sp for x in random.sample(range(self.n_trials), n_tr_obs_per_sp)] for sp in range(self.n_stim_patterns))))
            test_idxs = [x for x in range(n_obs) if x not in train_idxs]
            n_tr_obs = len(train_idxs)
            n_ts_obs = len(test_idxs)
            # initialize data structures for storage of results
            ts_decoded_mi_plugin = np.zeros(n_obs)
            ts_decoded_mi_bootstrap = np.zeros(n_obs)
            ts_decoded_mi_qe = np.zeros(n_obs)
            ts_decoded_mi_pt = np.zeros(n_obs)
            ts_decoded_mi_nsb = np.zeros(n_obs)
            tr_tree = np.zeros(shape=(n_tr_obs-1, 3))

            # compute mutual information by using direct clustering on training data (REMOVED)
            # --note: fcluster doesn't work in the border case with n_clusts=n_obs, as it never returns the trivial clustering. Cluster number 0 is never present in a clustering.
            tr_direct_mi = np.zeros(n_tr_obs-1)

            if self.linkage_method_string == 'kmeans':
                print('counting spikes in output spike trains')
                spike_counts = self.spikes_arch.get_spike_counts()
                Ym = self.n_stim_patterns
                Ny = np.array([self.n_trials for each in range(self.n_stim_patterns)])
                Xn = 1 # the output is effectively one-dimensional
                for n_clusts in range(min_clusts_analysed, max_clusts_analysed+1, clusts_step):
                    clustering = MiniBatchKMeans(n_clusters=n_clusts,
                                                 batch_size=self.n_stim_patterns*10)
                    print('performing k-means clustering for k='+str(n_clusts))
                    clustering.fit(spike_counts)
                    print('using k-means clustering to classify data points')
                    decoded_output = clustering.predict(spike_counts)
                    # calculate MI
                    print('calculating MI')
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
                    
            else:
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
                Ym = self.n_stim_patterns
                Ny = np.array([n_ts_obs_per_sp for each in range(self.n_stim_patterns)])
                Xn = 1 # the output is effectively one-dimensional

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

            # compute number of spikes fired by output cells
            print("Calculating output levels")
            o_level_array, o_level_hist_values, o_level_hist_edges = np.array([0]),np.array([0]),np.array([0]) #output_level(spikes)
            # compute output layer synchronisation (pairwise average of Schreiber's reliability measure)
            #random_observations_subset = spikes[random.sample(range(spikes.shape[0]), 1000)]
            #random_fields_subset = convolve(random_observations_subset, self.sim_duration, self.tau, self.dt)
            o_synchrony = np.array([0]) #np.mean([synchrony(p) for p in random_fields_subset])
            # save analysis results in the archive
            self.results_arch.update_result('tr_indexes', data=np.array(train_idxs))
            self.results_arch.update_result('tr_linkage', data=tr_tree)
            self.results_arch.update_result('tr_direct_mi', data=tr_direct_mi)
            self.results_arch.update_result('ts_decoded_mi_plugin', data=ts_decoded_mi_plugin)
            self.results_arch.update_result('ts_decoded_mi_bootstrap', data=ts_decoded_mi_bootstrap)
            self.results_arch.update_result('ts_decoded_mi_qe', data=ts_decoded_mi_qe)
            self.results_arch.update_result('ts_decoded_mi_pt', data=ts_decoded_mi_pt)
            self.results_arch.update_result('ts_decoded_mi_nsb', data=ts_decoded_mi_nsb)
            self.results_arch.update_result('px_at_same_size_point', data=px_at_same_size_point)
            self.results_arch.update_result('o_level_array', data=o_level_array)
            self.results_arch.update_result('o_level_hist_values', data=o_level_hist_values)
            self.results_arch.update_result('o_level_hist_edges', data=o_level_hist_edges)
            self.results_arch.update_result('o_synchrony', data=o_synchrony)
            # update attributes
            self.results_arch.load()


# A numpy ndarray with object dtype, and composed of (ParameterSpacePoint)s.
ParameterSpaceMesh = np.vectorize(ParameterSpacePoint)

class ParameterSpace(np.ndarray):
    ABSOLUTE_DIDXS = [
            'sim_duration',
            'min_mf_number',
            'grc_mf_ratio',
            'n_grc_dend',
            'network_scale',
            'active_mf_fraction',
            'bias',
            'stim_rate_mu',
            'stim_rate_sigma',
            'noise_rate_mu',
            'noise_rate_sigma',
            'n_stim_patterns',
            'n_trials',
            'training_size',
            'multineuron_metric_mixing',
            'linkage_method',
            'tau',
            'dt']
    def __new__(cls,
                sim_duration_slice,
                min_mf_number_slice,
                grc_mf_ratio_slice,
                n_grc_dend_slice,
                network_scale_slice,
                active_mf_fraction_slice,
                bias_slice,
                stim_rate_mu_slice,
                stim_rate_sigma_slice,
                noise_rate_mu_slice,
                noise_rate_sigma_slice,
                n_stim_patterns_slice,
                n_trials_slice,
                trainig_size_slice,
                multineuron_metric_mixing_slice,
                linkage_method_slice,
                tau_slice,
                dt_slice):
        # Create the ndarray instance of our type, given the usual
        # ndarray input arguments.  This will call the standard
        # ndarray constructor, but return an object of our type.
        # It also triggers a call to ParamSpace.__array_finalize__
        m = np.mgrid[sim_duration_slice,
                     min_mf_number_slice,
                     grc_mf_ratio_slice,
                     n_grc_dend_slice,
                     network_scale_slice,
                     active_mf_fraction_slice,
                     bias_slice,
                     stim_rate_mu_slice,
                     stim_rate_sigma_slice,
                     noise_rate_mu_slice,
                     noise_rate_sigma_slice,
                     n_stim_patterns_slice,
                     n_trials_slice,
                     trainig_size_slice,
                     multineuron_metric_mixing_slice,
                     linkage_method_slice,
                     tau_slice,
                     dt_slice]
        obj = ParameterSpaceMesh(*m).view(cls)
        # Dimensional InDeXS. Dictionary or list? I like the 'semantic' nature of a dictionary, but lists have a natural ordering and a natural way of updating when a dimension is added or removed.
        obj.DIDXS = list(ParameterSpace.ABSOLUTE_DIDXS)
        obj.ABBREVIATIONS = {
            'grc_mf_ratio': 'gmr',
            'n_grc_dend': 'gd',
            'network_scale': 's',
            'active_mf_fraction': 'mf',
            'bias': 'b',
            'stim_rate_mu': 'sm',
            'stim_rate_sigma': 'ss',
            'noise_rate_mu': 'nm',
            'noise_rate_sigma': 'ns',
            'n_stim_patterns': 'sp',
            'n_trials': 't',
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
        Return an array containing the value of the requested attribute for all the points in self.
        The shape of the returned array is of course self.shape. Default to np.nan when a point doesn't
        have the requested attribute: this is to increase robustness against missing datapoints
        (typical situation: trying to plot a 2d heatmap where we lack one of the values).
        """
        return np.array([getattr(x, attribute, np.nan) for x in self.flat]).reshape(self.shape)
    # def get_attribute_array(self, attribute, **kwargs):
    #     attr_list = []
    #     for x in self.flat:
    #         attr = getattr(x, attribute, None)
    #         if callable(attr):
    #             attr = attr(**kwargs)
    #         attr_list.append(attr)
    #     return np.array(attr_list).reshape(self.shape)
    def _get_idx_from_value(self, parameter, value):
        param_range = self.get_range(parameter)
        if value not in param_range:
            # Trying to find the index for a parameter (i.e., coordinate) value that's not on the mesh
            #   raises an exception.
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
        Return a new Parameter_space object made only of those points that satisfy the condition parameter==value.
        Don't change the dimensionality of the space (keep the hyperplane 'embedded' in the higher-dimensional space).
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
        Return a new Parameter_space object where one or more parameters (dimensions) have been fixed
        by setting parameter==value for every parameter, value pair in the parameter_value_pairs iterable.
        Reduce the dimensionality of the space, forgetting about the dimensions that have been fixed.
        """
        temp = self._get_hyperplane(*parameter_value_pairs[0])
        for pair in parameter_value_pairs[1:]:
            temp = temp._get_hyperplane(*pair)
        return temp
    def get_nontrivial_subspace(self, *parameter_value_pairs):
        """
        Return a new Parameter_space object where one or more parameters (dimensions) have been fixed
        by setting parameter==value for every (parameter, value) pair that is given as an argument.
        Reduce the dimensionality of the space, forgetting about the dimensions that have been fixed.
        Finally, reduce further the dimensionality by fogetting about the _trivial_ dimensions,
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

