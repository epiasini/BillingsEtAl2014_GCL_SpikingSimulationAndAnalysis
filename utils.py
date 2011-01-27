def ref_constructor(base_name, connection_pattern_index, stimulus_pattern_index, bias, trial):
    """Returns the simulation reference (the name of the relative subdirectory) of a given element in a batch of simulations."""
    sim_ref = base_name + "_cp" + str(connection_pattern_index) + "_sp" + str(stimulus_pattern_index) + "_b" + str(bias) + "_t" + str(trial)
    return sim_ref