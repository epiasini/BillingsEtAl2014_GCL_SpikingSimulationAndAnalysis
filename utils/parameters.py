import numpy as np
from matplotlib import pyplot as plt
from pure import ParamSpacePoint

class PSlice(object):
    """Conceptually, a subclass of the builtin slice class. Defaults to a single-point slice when just one argument is given."""
    def __init__(self, start, stop=None, step=1):
        self.start = start
        if stop==None:
            self.stop = start + 1
        else:
            self.stop = stop
        self.step = step

# A numpy ndarray with object dtype, and composed of (ParamSpacePoint)s.
ParamSpaceMesh = np.vectorize(ParamSpacePoint)

class ParamSpace(np.ndarray):
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
        obj = ParamSpaceMesh(*m).view(cls)
        # Dimensional InDeXS. Dictionary or list? I like the 'semantic' nature of a dictionary, but lists have a natural ordering and a natural way of updating when a dimension is added or removed.
        obj.DIDXS = [
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
            'mutineuron_metric_mixing',
            'linkage_method',
            'tau',
            'dt']
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
        self.DIDXS = getattr(obj, 'DIDXS', [])[:]
        # We do not need to return anything
    #-----------------------------------------------
    # Basic space and coordinate manipulation methods
    #-----------------------------------------------
    def _didx(self, parameter):
        return self.DIDXS.index(parameter)
    def _param(self, didx):
        return self.DIDXS[didx]
    def _get_range(self, parameter):
        # TODO: improving by vectorising getattr
        return np.unique(np.array([getattr(x, parameter, None) for x in self.flat]))
    def _get_attribute_array(self, attribute):
        return np.array([getattr(x, attribute, None) for x in self.flat]).reshape(self.shape)
    # def get_attribute_array(self, attribute, **kwargs):
    #     attr_list = []
    #     for x in self.flat:
    #         attr = getattr(x, attribute, None)
    #         if callable(attr):
    #             attr = attr(**kwargs)
    #         attr_list.append(attr)
    #     return np.array(attr_list).reshape(self.shape)
    def _get_idx_from_value(self, parameter, value):
        return np.searchsorted(self._get_range(parameter), value)
    def _remove_dimensional_index(self, parameter):
        del self.DIDXS[self._didx(parameter)]
    def _get_embedded_hyperplane(self, parameter, value):
        """
        Return a new Parameter_space object made only of those points that satisfy the condition parameter==value. Don't change the dimensionality of the space (keep the hyperplane 'embedded' in the higher-dimensional space).
        """
        idx = self._get_idx_from_value(parameter, value)
        return np.split(self, self.shape[self._didx(parameter)], self._didx(parameter))[idx]
    def _get_hyperplane(self, parameter, value):
        hp = self._get_embedded_hyperplane(parameter, value)
        hp = hp.squeeze(hp._didx(parameter))
        hp._remove_dimensional_index(parameter)
        return hp
    def get_subspace(self, *parameter_value_pairs):
        """
        Return a new Parameter_space object where one or more parameters (dimensions) have been fixed by setting parameter==value for every parameter, value pair in the parameter_value_pairs iterable. Reduce the dimensionality of the space, forgetting about the dimensions that have been fixed.
        """
        temp = self._get_hyperplane(*parameter_value_pairs[0])
        for pair in parameter_value_pairs[1:]:
            temp = temp._get_hyperplane(*pair)
        return temp
    def get_nontrivial_subspace(self, *parameter_value_pairs):
        """
        Return a new Parameter_space object where one or more parameters (dimensions) have been fixed by setting parameter==value for every (parameter, value) pair that is given as an argument. Reduce the dimensionality of the space, forgetting about the dimensions that have been fixed. Finally, reduce further the dimensionality by fogetting about the _trivial_ dimensions, i.e. the ones where our points don't have any variability.
        """
        temp = self.get_subspace(*parameter_value_pairs)
        for parameter in [temp._param(k) for k,s in enumerate(temp.shape) if s==1]:
            temp = temp.squeeze(temp._didx(parameter))
            temp._remove_dimensional_index(parameter)
        return temp
    #------------------
    # Analysis methods
    #------------------
    def run_analysis(self):
        for p in self.flat:
            p.run_analysis()
    def load_analysis_results(self):
        return [p.load_analysis_results() for p in self.flat]
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
        ax.set_xticklabels([str(x) for x in self._get_range(x_param)])
        ax.set_xlabel(x_param)
        ax.set_yticks(np.arange(self.shape[0]))
        ax.set_yticklabels([str(y) for y in self._get_range(y_param)])
        ax.set_ylabel(y_param)
        ax.set_title(fig_title)
        return fig, ax

