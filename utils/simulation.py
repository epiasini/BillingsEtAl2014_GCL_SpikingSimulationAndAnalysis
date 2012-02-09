from pure import ParamSpacePoint

def plast_correction_factor(f, syn_type):
    if syn_type=='AMPA':
        return 0.92 - 0.004*f + 1e-5 * f**2
    elif syn_type=='NMDA':
        return 0.94 - 0.0032*f + 8.5e-6 * f**2

class SimulationPoint(ParamSpacePoint):
    def get_ref(self, stim_pattern_index, trial):
        """Returns the simulation reference (the name of the relative subdirectory) of a given element in a batch of simulations."""
        sim_ref = "sp%d_t%d_spn%d_tn%d" % (self.n_stim_patterns, self.n_trials, stim_pattern_index, trial)
        return sim_ref
