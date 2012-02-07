import numpy as np
from matplotlib import pyplot as plt

class ParamSpacePoint(object):
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
                n_trials_slice):
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
                     n_trials_slice]
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
            'n_trials'
            ]
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
    def didx(self, parameter):
        return self.DIDXS.index(parameter)
    def param(self, didx):
        return self.DIDXS[didx]
    def get_range(self, parameter):
        # TODO: improving by vectorising getattr
        return np.unique(np.array([getattr(x, parameter, None) for x in self.flat]))
    def get_attribute_array(self, attribute):
        return np.array([getattr(x, attribute, None) for x in self.flat]).reshape(self.shape)
    # def get_attribute_array(self, attribute, **kwargs):
    #     attr_list = []
    #     for x in self.flat:
    #         attr = getattr(x, attribute, None)
    #         if callable(attr):
    #             attr = attr(**kwargs)
    #         attr_list.append(attr)
    #     return np.array(attr_list).reshape(self.shape)
    def get_idx_from_value(self, parameter, value):
        return np.searchsorted(self.get_range(parameter), value)
    def _remove_dimensional_index(self, parameter):
        # to be called only by get_hyperplane
        del self.DIDXS[self.didx(parameter)]
    def get_embedded_hyperplane(self, parameter, value):
        """Returns a new Parameter_space object made only of those points that satisfy the condition parameter==value. The dimensionality of the space is unchanged."""
        idx = self.get_idx_from_value(parameter, value)
        return np.split(self, self.shape[self.didx(parameter)], self.didx(parameter))[idx]
    def get_hyperplane(self, parameter, value):
        hp = self.get_embedded_hyperplane(parameter, value)
        hp = hp.squeeze(hp.didx(parameter))
        hp._remove_dimensional_index(parameter)
        return hp
    def get_subspace(self, *parameter_value_pairs):
        """Returns a new Parameter_space object where one or more parameters (dimensions) have been fixed, as an intersection of embedded hyperplanes. The dimensionality of the space is reduced to account for this (this is why the subspace is not 'embedded')."""
        temp = self.get_hyperplane(*parameter_value_pairs[0])
        for pair in parameter_value_pairs[1:]:
            temp = temp.get_hyperplane(*pair)
        return temp
    def get_nontrivial_subspace(self, *parameter_value_pairs):
        temp = self.get_subspace(*parameter_value_pairs)
        for parameter in [temp.param(k) for k,s in enumerate(temp.shape) if s==1]:
            temp = temp.squeeze(temp.didx(parameter))
            temp._remove_dimensional_index(parameter)
        return temp
    def plot_2d_heatmap(self, heat_dim, fig_title=''):
        if len(self.shape) > 2:
            raise Exception()
        x_param = self.param(1)
        y_param = self.param(0)        
        fig, ax = plt.subplots()
        plot = ax.imshow(self.get_attribute_array(heat_dim), interpolation='none', cmap='coolwarm', origin='lower')
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
        
